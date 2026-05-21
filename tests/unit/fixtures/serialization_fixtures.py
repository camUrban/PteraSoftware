"""This module contains functions to create objects for use in serialization tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware._panel import Panel

# noinspection PyProtectedMember
from pterasoftware.geometry.airfoil import Airfoil
from pterasoftware.geometry.airplane import Airplane
from pterasoftware.geometry.wing import Wing
from pterasoftware.geometry.wing_cross_section import WingCrossSection
from pterasoftware.movements.airplane_movement import AirplaneMovement
from pterasoftware.movements.movement import Movement
from pterasoftware.movements.operating_point_movement import OperatingPointMovement
from pterasoftware.movements.wing_cross_section_movement import (
    WingCrossSectionMovement,
)
from pterasoftware.movements.wing_movement import WingMovement
from pterasoftware.operating_point import OperatingPoint
from pterasoftware.problems import SteadyProblem, UnsteadyProblem


def make_basic_panel_fixture() -> Panel:
    """This method makes a fixture that is a Panel with basic configuration for
    serialization testing.

    :return basic_panel_fixture: Panel
        This is a unit square Panel that is a leading edge Panel.
    """
    return Panel(
        Frpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Flpp_G_Cg=np.array([1.0, 1.0, 0.0]),
        Blpp_G_Cg=np.array([0.0, 1.0, 0.0]),
        Brpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        is_leading_edge=True,
        is_trailing_edge=False,
    )


def make_meshed_airplane_fixture() -> Airplane:
    """This method makes a fixture that is a meshed Airplane for serialization testing.

    :return meshed_airplane_fixture: Airplane
        This is a meshed Airplane with one Wing containing 2 x 2 Panels.
    """
    wing = Wing(
        wing_cross_sections=[
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=2,
                chord=1.0,
            ),
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=1.0,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
            ),
        ],
        num_chordwise_panels=2,
        chordwise_spacing="uniform",
    )
    airplane = Airplane(wings=[wing])
    operating_point = OperatingPoint()
    SteadyProblem(airplanes=[airplane], operating_point=operating_point)
    return airplane


def make_steady_problem_fixture() -> SteadyProblem:
    """This method makes a fixture that is a SteadyProblem for serialization testing.

    :return steady_problem_fixture: SteadyProblem
        This is a SteadyProblem with one meshed Airplane and 2 x 2 Panels.
    """
    wing = Wing(
        wing_cross_sections=[
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=2,
                chord=1.0,
            ),
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=1.0,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
            ),
        ],
        num_chordwise_panels=2,
        chordwise_spacing="uniform",
    )
    airplane = Airplane(wings=[wing])
    operating_point = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
    return SteadyProblem(airplanes=[airplane], operating_point=operating_point)


def _make_two_airplane_wing() -> Wing:
    """This method makes a Wing for use in multi Airplane fixture functions.

    :return two_airplane_wing: Wing
        This is an unmeshed Wing with 2 x 2 Panels.
    """
    return Wing(
        wing_cross_sections=[
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=2,
                chord=1.0,
            ),
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=1.0,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
            ),
        ],
        num_chordwise_panels=2,
        chordwise_spacing="uniform",
    )


def make_formation_steady_problem_fixture() -> SteadyProblem:
    """This method makes a fixture that is a SteadyProblem with two Airplanes for
    formation flight serialization testing.

    :return formation_steady_problem_fixture: SteadyProblem
        This is a SteadyProblem with two meshed Airplanes.
    """
    airplane1 = Airplane(wings=[_make_two_airplane_wing()])
    airplane2 = Airplane(
        wings=[_make_two_airplane_wing()],
        Cg_GP1_CgP1=(0.0, 10.0, 0.0),
    )
    operating_point = OperatingPoint()
    return SteadyProblem(
        airplanes=[airplane1, airplane2], operating_point=operating_point
    )


def make_unsteady_problem_fixture() -> UnsteadyProblem:
    """This method makes a fixture that is an UnsteadyProblem for serialization testing.

    :return unsteady_problem_fixture: UnsteadyProblem
        This is an UnsteadyProblem with a static Movement and 2 x 2 Panels.
    """
    wing_cross_section_movement_root = WingCrossSectionMovement(
        base_wing_cross_section=WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=2,
            chord=1.0,
        ),
    )
    wing_cross_section_movement_tip = WingCrossSectionMovement(
        base_wing_cross_section=WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=None,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
        ),
    )
    wing_movement = WingMovement(
        base_wing=Wing(
            wing_cross_sections=[
                wing_cross_section_movement_root.base_wing_cross_section,
                wing_cross_section_movement_tip.base_wing_cross_section,
            ],
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
        ),
        wing_cross_section_movements=[
            wing_cross_section_movement_root,
            wing_cross_section_movement_tip,
        ],
    )
    airplane_movement = AirplaneMovement(
        base_airplane=Airplane(wings=[wing_movement.base_wing]),
        wing_movements=[wing_movement],
    )
    operating_point_movement = OperatingPointMovement(
        base_operating_point=OperatingPoint(),
    )
    movement = Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_steps=3,
    )
    return UnsteadyProblem(movement=movement)


def make_formation_unsteady_problem_fixture() -> UnsteadyProblem:
    """This method makes a fixture that is an UnsteadyProblem with two Airplanes for
    formation flight serialization testing.

    :return formation_unsteady_problem_fixture: UnsteadyProblem
        This is an UnsteadyProblem with two Airplanes, 3 time steps, and 2 x 2 Panels
        per Wing.
    """
    wing_cross_section_movements = [
        WingCrossSectionMovement(
            base_wing_cross_section=WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=2,
                chord=1.0,
            ),
        ),
        WingCrossSectionMovement(
            base_wing_cross_section=WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=1.0,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
            ),
        ),
    ]
    wing_movement1 = WingMovement(
        base_wing=Wing(
            wing_cross_sections=[
                wing_cross_section_movement.base_wing_cross_section
                for wing_cross_section_movement in wing_cross_section_movements
            ],
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
        ),
        wing_cross_section_movements=wing_cross_section_movements,
    )
    airplane_movement1 = AirplaneMovement(
        base_airplane=Airplane(wings=[wing_movement1.base_wing]),
        wing_movements=[wing_movement1],
    )

    wing_cross_section_movements2 = [
        WingCrossSectionMovement(
            base_wing_cross_section=WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=2,
                chord=1.0,
            ),
        ),
        WingCrossSectionMovement(
            base_wing_cross_section=WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=1.0,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
            ),
        ),
    ]
    wing_movement2 = WingMovement(
        base_wing=Wing(
            wing_cross_sections=[
                wing_cross_section_movement.base_wing_cross_section
                for wing_cross_section_movement in wing_cross_section_movements2
            ],
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
        ),
        wing_cross_section_movements=wing_cross_section_movements2,
    )
    airplane_movement2 = AirplaneMovement(
        base_airplane=Airplane(
            wings=[wing_movement2.base_wing],
            Cg_GP1_CgP1=(0.0, 10.0, 0.0),
        ),
        wing_movements=[wing_movement2],
    )

    movement = Movement(
        airplane_movements=[airplane_movement1, airplane_movement2],
        operating_point_movement=OperatingPointMovement(
            base_operating_point=OperatingPoint(),
        ),
        num_steps=3,
    )
    return UnsteadyProblem(movement=movement)
