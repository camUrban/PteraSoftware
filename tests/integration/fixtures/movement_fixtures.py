# NOTE: I haven't yet started refactoring this module.
"""This module creates movement objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_static_validation_movement: This function creates a movement object with
    static geometry to be used as a fixture.

    make_variable_validation_movement: This function creates a movement object with
    variable geometry to be used as a fixture.

    make_multiple_wing_static_validation_movement: This function creates a movement
    object with static, multi-wing geometry to be used as a fixture.

    make_multiple_wing_variable_validation_movement: This function creates a movement
    object with variable, multi-wing geometry to be used as a fixture.
"""

import pterasoftware as ps
from tests.integration.fixtures import airplane_fixtures
from tests.integration.fixtures import operating_point_fixtures


def make_static_validation_movement():
    """This function creates a movement with static geometry to be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a movement object with static geometry to be used as a fixture.
    """

    # Construct an airplane object and an operating point object.
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    # Create a wing cross section movement object associated with this airplane's
    # root wing cross section.
    unsteady_validation_root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0]
        )
    )

    # Create a wing cross section movement object associated with this airplane's tip
    # wing cross section.
    unsteady_validation_tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
        )
    )

    # Create a wing movement object associated with this airplane's wing.
    unsteady_validation_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_root_wing_cross_section_movement,
            unsteady_validation_tip_wing_cross_section_movement,
        ],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_root_wing_cross_section_movement
    del unsteady_validation_tip_wing_cross_section_movement

    # Create an airplane movement object associated with this airplane.
    unsteady_validation_airplane_movement = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=unsteady_validation_airplane,
            wing_movements=[unsteady_validation_wing_movement],
        )
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane
    del unsteady_validation_wing_movement

    # Create an operating point movement object associated with this operating point.
    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_operating_point

    # Create a movement object associated with this airplane and operating point.
    unsteady_validation_movement = ps.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_chords=6,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement

    # Return the movement fixture.
    return unsteady_validation_movement


def make_variable_validation_movement():
    """This function creates a movement with variable geometry to be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a movement object with variable geometry to be used as a fixture.
    """

    # Construct an airplane object and an operating point object.
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    # Create a wing cross section movement object associated with this airplane's
    # root wing cross section.
    unsteady_validation_root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0],
        )
    )

    # Create a wing cross section movement object associated with this airplane's tip
    # wing cross section.
    unsteady_validation_tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
            sweeping_amplitude=30.0,
            sweeping_period=1.0,
            sweeping_spacing="sine",
            pitching_amplitude=30.0,
            pitching_period=0.5,
            pitching_spacing="sine",
            heaving_amplitude=30.0,
            heaving_period=0.5,
            heaving_spacing="sine",
        )
    )

    # Create a wing movement object associated with this airplane's wing.
    unsteady_validation_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_root_wing_cross_section_movement,
            unsteady_validation_tip_wing_cross_section_movement,
        ],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_root_wing_cross_section_movement
    del unsteady_validation_tip_wing_cross_section_movement

    # Create an airplane movement object associated with this airplane.
    unsteady_validation_airplane_movement = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=unsteady_validation_airplane,
            wing_movements=[
                unsteady_validation_wing_movement,
            ],
        )
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane
    del unsteady_validation_wing_movement

    # Create an operating point movement object associated with this operating point.
    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_operating_point

    # Create a movement object associated with this airplane and operating point.
    unsteady_validation_movement = ps.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement

    # Return the movement fixture.
    return unsteady_validation_movement


def make_multiple_wing_static_validation_movement():
    """This function creates a movement object with static, multi-wing geometry to
    be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a movement object with variable geometry to be used as a fixture.
    """

    # Construct an airplane object and an operating point object.
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_multiple_wing_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    # Create a wing cross section movement object associated with this airplane's
    # main wing's root wing cross section.
    unsteady_validation_main_wing_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0]
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # main wing's tip wing cross section.
    unsteady_validation_main_wing_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # horizontal stabilizer's root wing cross section.
    unsteady_validation_hstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[0]
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # horizontal stabilizer's tip wing cross section.
    unsteady_validation_hstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[1],
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # vertical stabilizer's root wing cross section.
    unsteady_validation_vstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[0]
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # vertical stabilizer's tip wing cross section.
    unsteady_validation_vstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[1],
        )
    )

    # Create a wing movement object associated with this airplane's main wing.
    unsteady_validation_main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_main_wing_root_cross_section_movement,
            unsteady_validation_main_wing_tip_cross_section_movement,
        ],
    )

    # Create a wing movement object associated with this airplane's horizontal
    # stabilizer.
    unsteady_validation_hstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[1],
        wing_cross_section_movements=[
            unsteady_validation_hstab_root_cross_section_movement,
            unsteady_validation_hstab_tip_cross_section_movement,
        ],
    )

    # Create a wing movement object associated with this airplane's vertical stabilizer.
    unsteady_validation_vstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[2],
        wing_cross_section_movements=[
            unsteady_validation_vstab_root_cross_section_movement,
            unsteady_validation_vstab_tip_cross_section_movement,
        ],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_main_wing_root_cross_section_movement
    del unsteady_validation_main_wing_tip_cross_section_movement
    del unsteady_validation_hstab_root_cross_section_movement
    del unsteady_validation_hstab_tip_cross_section_movement
    del unsteady_validation_vstab_root_cross_section_movement
    del unsteady_validation_vstab_tip_cross_section_movement

    # Create an airplane movement object associated with this airplane.
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

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane
    del unsteady_validation_main_wing_movement
    del unsteady_validation_hstab_movement
    del unsteady_validation_vstab_movement

    # Create an operating point movement object associated with this operating point.
    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_operating_point

    # Create a movement object associated with this airplane and operating point.
    unsteady_validation_movement = ps.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=8,
        delta_time=1 / 8 / 10,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement

    # Return the movement fixture.
    return unsteady_validation_movement


def make_multiple_wing_variable_validation_movement():
    """This function creates a movement object with variable, multi-wing geometry to
    be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a movement object with variable geometry to be used as a fixture.
    """

    # Construct an airplane object and an operating point object.
    unsteady_validation_airplane = (
        airplane_fixtures.make_symmetric_multiple_wing_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    # Create a wing cross section movement object associated with this airplane's
    # main wing's root wing cross section.
    unsteady_validation_main_wing_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[0],
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # main wing's tip wing cross section.
    unsteady_validation_main_wing_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                0
            ].wing_cross_sections[1],
            sweeping_amplitude=30.0,
            sweeping_period=1.0,
            sweeping_spacing="sine",
            heaving_amplitude=15.0,
            heaving_period=1.0,
            heaving_spacing="sine",
            pitching_amplitude=15.0,
            pitching_period=0.5,
            pitching_spacing="sine",
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # horizontal stabilizer's root wing cross section.
    unsteady_validation_hstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[0]
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # horizontal stabilizer's tip wing cross section.
    unsteady_validation_hstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                1
            ].wing_cross_sections[1],
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # vertical stabilizer's root wing cross section.
    unsteady_validation_vstab_root_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[0]
        )
    )

    # Create a wing cross section movement object associated with this airplane's
    # vertical stabilizer's tip wing cross section.
    unsteady_validation_vstab_tip_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=unsteady_validation_airplane.wings[
                2
            ].wing_cross_sections[1],
        )
    )

    # Create a wing movement object associated with this airplane's main wing.
    unsteady_validation_main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_section_movements=[
            unsteady_validation_main_wing_root_cross_section_movement,
            unsteady_validation_main_wing_tip_cross_section_movement,
        ],
    )

    # Create a wing movement object associated with this airplane's horizontal
    # stabilizer.
    unsteady_validation_hstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[1],
        wing_cross_section_movements=[
            unsteady_validation_hstab_root_cross_section_movement,
            unsteady_validation_hstab_tip_cross_section_movement,
        ],
    )

    # Create a wing movement object associated with this airplane's vertical stabilizer.
    unsteady_validation_vstab_movement = ps.movements.wing_movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[2],
        wing_cross_section_movements=[
            unsteady_validation_vstab_root_cross_section_movement,
            unsteady_validation_vstab_tip_cross_section_movement,
        ],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_main_wing_root_cross_section_movement
    del unsteady_validation_main_wing_tip_cross_section_movement
    del unsteady_validation_hstab_root_cross_section_movement
    del unsteady_validation_hstab_tip_cross_section_movement
    del unsteady_validation_vstab_root_cross_section_movement
    del unsteady_validation_vstab_tip_cross_section_movement

    # Create an airplane movement object associated with this airplane.
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

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane
    del unsteady_validation_main_wing_movement
    del unsteady_validation_hstab_movement
    del unsteady_validation_vstab_movement

    # Create an operating point movement object associated with this operating point.
    unsteady_validation_operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=unsteady_validation_operating_point
        )
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_operating_point

    # Create a movement object associated with this airplane and operating point.
    unsteady_validation_movement = ps.movement.Movement(
        airplane_movements=[unsteady_validation_airplane_movement],
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=20,
        delta_time=1 / 8 / 10,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement

    # Return the movement fixture.
    return unsteady_validation_movement
