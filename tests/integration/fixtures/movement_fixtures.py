"""This module creates movement objects to be used as fixtures."""

import pterasoftware as ps
from tests.integration.fixtures import airplane_fixtures
from tests.integration.fixtures import operating_point_fixtures


def make_static_validation_movement():
    """This function creates a Movement with static geometry to be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a Movement with static geometry to be used as a fixture.
    """
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    unsteady_validation_root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0]
        )
    )

    unsteady_validation_tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
        )
    )

    unsteady_validation_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_root_wing_cross_section_movement,
            unsteady_validation_tip_wing_cross_section_movement,
        ],
    )

    unsteady_validation_airplane_movement = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=unsteady_validation_airplane,
            wing_movements=[unsteady_validation_wing_movement],
        )
    )

    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    unsteady_validation_movement = ps.movements.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_chords=6,
    )

    return unsteady_validation_movement


def make_variable_validation_movement():
    """This function creates a Movement with variable geometry to be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a Movement with variable geometry to be used as a fixture.
    """
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    unsteady_validation_root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0],
        )
    )

    unsteady_validation_tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_izyx=(30.0, 30.0, 30.0),
            periodAngles_Wcsp_to_Wcs_izyx=(1.0, 0.5, 0.5),
            spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        )
    )

    unsteady_validation_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_root_wing_cross_section_movement,
            unsteady_validation_tip_wing_cross_section_movement,
        ],
    )

    unsteady_validation_airplane_movement = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=unsteady_validation_airplane,
            wing_movements=[
                unsteady_validation_wing_movement,
            ],
        )
    )

    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    unsteady_validation_movement = ps.movements.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
    )

    return unsteady_validation_movement


def make_multiple_wing_static_validation_movement():
    """This function creates a Movement with static, multi-wing geometry to be used
    as a fixture.

    :return unsteady_validation_movement: Movement
        This is a Movement with variable geometry to be used as a fixture.
    """
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_multiple_wing_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    unsteady_validation_main_wing_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0]
        )
    )

    unsteady_validation_main_wing_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
        )
    )

    unsteady_validation_hstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[0]
        )
    )

    unsteady_validation_hstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[1],
        )
    )

    unsteady_validation_vstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[0]
        )
    )

    unsteady_validation_vstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[1],
        )
    )

    unsteady_validation_main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_main_wing_root_cross_section_movement,
            unsteady_validation_main_wing_tip_cross_section_movement,
        ],
    )

    unsteady_validation_hstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[1],
        wing_cross_section_movements=[
            unsteady_validation_hstab_root_cross_section_movement,
            unsteady_validation_hstab_tip_cross_section_movement,
        ],
    )

    unsteady_validation_vstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[2],
        wing_cross_section_movements=[
            unsteady_validation_vstab_root_cross_section_movement,
            unsteady_validation_vstab_tip_cross_section_movement,
        ],
    )

    unsteady_validation_airplane_movement = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=unsteady_validation_airplane,
            wing_movements=[
                unsteady_validation_main_wing_movement,
                unsteady_validation_hstab_movement,
                unsteady_validation_vstab_movement,
            ],
        )
    )

    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    unsteady_validation_movement = ps.movements.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=8,
        delta_time=1 / 8 / 10,
    )

    return unsteady_validation_movement


def make_multiple_wing_variable_validation_movement():
    """This function creates a Movement with variable, multi-wing geometry to
    be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a Movement with variable geometry to be used as a fixture.
    """
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_multiple_wing_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    unsteady_validation_main_wing_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0],
        )
    )

    unsteady_validation_main_wing_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_izyx=(30.0, 15.0, 15.0),
            periodAngles_Wcsp_to_Wcs_izyx=(1.0, 1.0, 0.5),
            spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        )
    )

    unsteady_validation_hstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[0]
        )
    )

    unsteady_validation_hstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[1],
        )
    )

    unsteady_validation_vstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[0]
        )
    )

    unsteady_validation_vstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[1],
        )
    )

    unsteady_validation_main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_main_wing_root_cross_section_movement,
            unsteady_validation_main_wing_tip_cross_section_movement,
        ],
    )

    unsteady_validation_hstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[1],
        wing_cross_section_movements=[
            unsteady_validation_hstab_root_cross_section_movement,
            unsteady_validation_hstab_tip_cross_section_movement,
        ],
    )

    unsteady_validation_vstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[2],
        wing_cross_section_movements=[
            unsteady_validation_vstab_root_cross_section_movement,
            unsteady_validation_vstab_tip_cross_section_movement,
        ],
    )

    unsteady_validation_airplane_movement = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=unsteady_validation_airplane,
            wing_movements=[
                unsteady_validation_main_wing_movement,
                unsteady_validation_hstab_movement,
                unsteady_validation_vstab_movement,
            ],
        )
    )

    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    unsteady_validation_movement = ps.movements.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=20,
        delta_time=1 / 8 / 10,
    )

    return unsteady_validation_movement
