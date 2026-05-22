"""This module contains classes to test the FreeFlightUnsteadyProblem class."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import movement_fixtures


def _make_minimal_free_flight_movement():
    """Build a minimal FreeFlightMovement for testing FreeFlightUnsteadyProblem.

    Uses a single static type-1 Wing with no prescribed motion so that symmetry
    types are guaranteed to remain constant across all time steps.

    :return: A FreeFlightMovement suitable for use as a test fixture.
    """
    airfoil = ps.geometry.airfoil.Airfoil(name="naca0012")

    root_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=airfoil,
        num_spanwise_panels=1,
        chord=1.0,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spanwise_spacing="uniform",
    )
    tip_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=airfoil,
        num_spanwise_panels=None,
        chord=0.5,
        Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
        spanwise_spacing=None,
    )

    wing = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
        symmetric=False,
        mirror_only=False,
        num_chordwise_panels=2,
        chordwise_spacing="uniform",
    )

    airplane = ps.geometry.airplane.Airplane(
        wings=[wing],
        name="Free Flight Test Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )

    root_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
        base_wing_cross_section=airplane.wings[0].wing_cross_sections[0],
    )
    tip_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
        base_wing_cross_section=airplane.wings[0].wing_cross_sections[1],
    )

    wing_movement = ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=[
            root_wing_cross_section_movement,
            tip_wing_cross_section_movement,
        ],
    )

    airplane_movement = (
        ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement(
            base_airplane=airplane,
            wing_movements=[wing_movement],
        )
    )

    operating_point = ps.operating_point.OperatingPoint(
        rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0, externalFX_W=0.0, nu=15.06e-6
    )

    operating_point_movement = ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
        base_operating_point=operating_point,
    )

    return ps.movements.free_flight_movement.FreeFlightMovement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        delta_time=0.1,
        prescribed_num_steps=1,
        free_num_steps=1,
    )


class TestFreeFlightUnsteadyProblem(unittest.TestCase):
    """Tests for the FreeFlightUnsteadyProblem class."""

    def setUp(self) -> None:
        """Create a FreeFlightUnsteadyProblem fixture.

        :return: None
        """
        self.movement = _make_minimal_free_flight_movement()
        self.problem = ps.problems.FreeFlightUnsteadyProblem(movement=self.movement)

    def test_initialization_rejects_non_free_flight_movement(self) -> None:
        """FreeFlightUnsteadyProblem must raise TypeError when passed a Movement
        that is not a FreeFlightMovement.

        :return: None
        """
        basic_movement = movement_fixtures.make_basic_movement_fixture()
        with self.assertRaises(TypeError):
            ps.problems.FreeFlightUnsteadyProblem(movement=basic_movement)

    def test_initialization_accepts_free_flight_movement(self) -> None:
        """FreeFlightUnsteadyProblem must accept a valid FreeFlightMovement.

        :return: None
        """
        self.assertIsInstance(self.problem, ps.problems.FreeFlightUnsteadyProblem)

    def test_num_steps_equals_movement_num_steps(self) -> None:
        """num_steps must equal the FreeFlightMovement's num_steps (prescribed +
        free).

        :return: None
        """
        self.assertEqual(self.problem.num_steps, self.movement.num_steps)

    def test_initial_steady_problems_count_is_one(self) -> None:
        """Exactly one SteadyProblem must exist immediately after initialization
        (the problem at time step zero).

        :return: None
        """
        self.assertEqual(len(self.problem.steady_problems), 1)

    def test_initialize_next_problem_raises_not_implemented(self) -> None:
        """initialize_next_problem must raise NotImplementedError because the free
        flight solver is not yet implemented.

        :return: None
        """
        with self.assertRaises(NotImplementedError):
            self.problem.initialize_next_problem(solver=None)  # type: ignore[arg-type]
