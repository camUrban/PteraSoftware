""" This module creates movement objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_static_validation_movement: This function creates a movement object with static geometry to be used as a
                                     fixture.
    make_variable_validation_movement: This function creates a movement object with variable geometry to be used as a
                                       fixture.
"""

import aviansoftwareminimumviableproduct as asmvp
import tests.integration.fixtures.airplane_fixtures
import tests.integration.fixtures.operating_point_fixtures


def make_static_validation_movement():
    """ This function creates a movement with static geometry to be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a movement object with static geometry to be used as a fixture.
    """

    # Construct an airplane object and an operating point object.
    unsteady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_asymmetric_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )

    # Create a wing cross section movement object associated with this airplane's root wing cross section.
    unsteady_validation_root_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
        base_wing_cross_section=unsteady_validation_airplane.wings[
            0
        ].wing_cross_sections[0]
    )

    # Create a wing cross section movement object associated with this airplane's tip wing cross section.
    unsteady_validation_tip_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
        base_wing_cross_section=unsteady_validation_airplane.wings[
            0
        ].wing_cross_sections[1],
    )

    # Create a wing movement object associated with this airplane's wing.
    unsteady_validation_wing_movement = asmvp.movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_sections_movements=[
            unsteady_validation_root_wing_cross_section_movement,
            unsteady_validation_tip_wing_cross_section_movement,
        ],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_root_wing_cross_section_movement
    del unsteady_validation_tip_wing_cross_section_movement

    # Create an airplane movement object associated with this airplane.
    unsteady_validation_airplane_movement = asmvp.movement.AirplaneMovement(
        base_airplane=unsteady_validation_airplane,
        wing_movements=[unsteady_validation_wing_movement],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane
    del unsteady_validation_wing_movement

    # Create an operating point movement object associated with this operating point.
    unsteady_validation_operating_point_movement = asmvp.movement.OperatingPointMovement(
        base_operating_point=unsteady_validation_operating_point
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_operating_point

    # Create a movement object associated with this airplane and operating point.
    unsteady_validation_movement = asmvp.movement.Movement(
        airplane_movement=unsteady_validation_airplane_movement,
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=60,
        delta_time=1 / 6 / 10,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement

    # Return the movement fixture.
    return unsteady_validation_movement


def make_variable_validation_movement():
    """ This function creates a movement with variable geometry to be used as a fixture.

    :return unsteady_validation_movement: Movement
        This is a movement object with variable geometry to be used as a fixture.
    """

    # Construct an airplane object and an operating point object.
    unsteady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_symmetric_unsteady_validation_airplane()
    )
    unsteady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )

    # Create a wing cross section movement object associated with this airplane's root wing cross section.
    unsteady_validation_root_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
        base_wing_cross_section=unsteady_validation_airplane.wings[
            0
        ].wing_cross_sections[0],
    )

    # Create a wing cross section movement object associated with this airplane's tip wing cross section.
    unsteady_validation_tip_wing_cross_section_movement = asmvp.movement.WingCrossSectionMovement(
        base_wing_cross_section=unsteady_validation_airplane.wings[
            0
        ].wing_cross_sections[1],
        z_le_amplitude=2.0,
        z_le_period=1.0,
        z_le_spacing="sine",
    )

    # Create a wing movement object associated with this airplane's wing.
    unsteady_validation_wing_movement = asmvp.movement.WingMovement(
        base_wing=unsteady_validation_airplane.wings[0],
        wing_cross_sections_movements=[
            unsteady_validation_root_wing_cross_section_movement,
            unsteady_validation_tip_wing_cross_section_movement,
        ],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_root_wing_cross_section_movement
    del unsteady_validation_tip_wing_cross_section_movement

    # Create an airplane movement object associated with this airplane.
    unsteady_validation_airplane_movement = asmvp.movement.AirplaneMovement(
        base_airplane=unsteady_validation_airplane,
        wing_movements=[unsteady_validation_wing_movement],
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane
    del unsteady_validation_wing_movement

    # Create an operating point movement object associated with this operating point.
    unsteady_validation_operating_point_movement = asmvp.movement.OperatingPointMovement(
        base_operating_point=unsteady_validation_operating_point,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_operating_point

    # Create a movement object associated with this airplane and operating point.
    unsteady_validation_movement = asmvp.movement.Movement(
        airplane_movement=unsteady_validation_airplane_movement,
        operating_point_movement=unsteady_validation_operating_point_movement,
        num_steps=60,
        delta_time=1 / 6 / 10,
    )

    # Delete the now extraneous constructing fixtures.
    del unsteady_validation_airplane_movement
    del unsteady_validation_operating_point_movement

    # Return the movement fixture.
    return unsteady_validation_movement
