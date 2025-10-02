"""This module contains functions to create Movements for use in tests."""

import pterasoftware as ps

from . import airplane_movement_fixtures
from . import operating_point_fixtures


def make_static_movement_fixture():
    """This method makes a fixture that is a Movement with all static components.

    :return static_movement_fixture: Movement
        This is the Movement with no motion.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_static_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the static Movement.
    static_movement_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_chords=3,
    )

    # Delete the constructing fixtures.
    del airplane_movements
    del operating_point_movement

    # Return the Movement fixture.
    return static_movement_fixture


def make_basic_movement_fixture():
    """This method makes a fixture that is a Movement with general-purpose values.

    :return basic_movement_fixture: Movement
        This is the Movement with general-purpose values for testing.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the basic Movement.
    basic_movement_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_cycles=1,
    )

    # Delete the constructing fixtures.
    del airplane_movements
    del operating_point_movement

    # Return the Movement fixture.
    return basic_movement_fixture


def make_static_movement_with_explicit_num_steps_fixture():
    """This method makes a fixture that is a Movement with static motion and
    explicitly set num_steps.

    :return static_movement_with_explicit_num_steps_fixture: Movement
        This is the Movement with static motion and explicit num_steps.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_static_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with explicit num_steps.
    static_movement_with_explicit_num_steps_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_steps=5,
    )

    # Delete the constructing fixtures.
    del airplane_movements
    del operating_point_movement

    # Return the Movement fixture.
    return static_movement_with_explicit_num_steps_fixture


def make_non_static_movement_with_explicit_num_steps_fixture():
    """This method makes a fixture that is a Movement with non-static motion and
    explicitly set num_steps.

    :return non_static_movement_with_explicit_num_steps_fixture: Movement
        This is the Movement with non-static motion and explicit num_steps.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with explicit num_steps.
    non_static_movement_with_explicit_num_steps_fixture = (
        ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=10,
        )
    )

    # Delete the constructing fixtures.
    del airplane_movements
    del operating_point_movement

    # Return the Movement fixture.
    return non_static_movement_with_explicit_num_steps_fixture


def make_movement_with_custom_delta_time_fixture():
    """This method makes a fixture that is a Movement with custom delta_time.

    :return movement_with_custom_delta_time_fixture: Movement
        This is the Movement with custom delta_time parameter.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with custom delta_time.
    movement_with_custom_delta_time_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        delta_time=0.05,
        num_cycles=1,
    )

    # Delete the constructing fixtures.
    del airplane_movements
    del operating_point_movement

    # Return the Movement fixture.
    return movement_with_custom_delta_time_fixture


def make_movement_with_multiple_airplanes_fixture():
    """This method makes a fixture that is a Movement with multiple AirplaneMovements.

    :return movement_with_multiple_airplanes_fixture: Movement
        This is the Movement with multiple AirplaneMovements.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with multiple AirplaneMovements.
    movement_with_multiple_airplanes_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_cycles=1,
    )

    # Delete the constructing fixtures.
    del airplane_movements
    del operating_point_movement

    # Return the Movement fixture.
    return movement_with_multiple_airplanes_fixture
