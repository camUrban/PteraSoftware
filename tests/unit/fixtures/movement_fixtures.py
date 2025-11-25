"""This module contains functions to create Movements and CoupledMovements for use
in tests."""

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

    # Return the Movement fixture.
    return movement_with_custom_delta_time_fixture


def make_cyclic_movement_fixture():
    """This method makes a fixture that is a Movement with cyclic motion.

    :return cyclic_movement_fixture: Movement
        This is the Movement with cyclic motion.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the cyclic Movement.
    cyclic_movement_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_cycles=2,
    )

    # Return the Movement fixture.
    return cyclic_movement_fixture


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

    # Return the Movement fixture.
    return movement_with_multiple_airplanes_fixture


def make_basic_coupled_movement_fixture():
    """This method makes a fixture that is a CoupledMovement with general-purpose
    values for testing.

    :return basic_coupled_movement_fixture: CoupledMovement
        This is the CoupledMovement with general-purpose values for testing.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_basic_coupled_operating_point_fixture()
    )

    # Create the basic CoupledMovement.
    basic_coupled_movement_fixture = ps.movements.movement.CoupledMovement(
        airplane_movement=airplane_movement,
        initial_coupled_operating_point=initial_coupled_operating_point,
        delta_time=0.1,
        prescribed_num_steps=5,
        free_num_steps=5,
    )

    # Return the CoupledMovement fixture.
    return basic_coupled_movement_fixture


def make_static_coupled_movement_fixture():
    """This method makes a fixture that is a CoupledMovement with all static
    components.

    :return static_coupled_movement_fixture: CoupledMovement
        This is the CoupledMovement with no motion.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_static_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_basic_coupled_operating_point_fixture()
    )

    # Create the static CoupledMovement.
    static_coupled_movement_fixture = ps.movements.movement.CoupledMovement(
        airplane_movement=airplane_movement,
        initial_coupled_operating_point=initial_coupled_operating_point,
        delta_time=0.1,
        prescribed_num_steps=3,
        free_num_steps=3,
    )

    # Return the CoupledMovement fixture.
    return static_coupled_movement_fixture


def make_coupled_movement_with_angular_speed_fixture():
    """This method makes a fixture that is a CoupledMovement with non zero angular
    speed in its initial CoupledOperatingPoint.

    :return coupled_movement_with_angular_speed_fixture: CoupledMovement
        This is the CoupledMovement with non zero initial angular speed.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_with_angular_speed_coupled_operating_point_fixture()
    )

    # Create the CoupledMovement with angular speed.
    coupled_movement_with_angular_speed_fixture = ps.movements.movement.CoupledMovement(
        airplane_movement=airplane_movement,
        initial_coupled_operating_point=initial_coupled_operating_point,
        delta_time=0.1,
        prescribed_num_steps=5,
        free_num_steps=5,
    )

    # Return the CoupledMovement fixture.
    return coupled_movement_with_angular_speed_fixture


def make_coupled_movement_with_attitude_angles_fixture():
    """This method makes a fixture that is a CoupledMovement with non zero attitude
    angles in its initial CoupledOperatingPoint.

    :return coupled_movement_with_attitude_angles_fixture: CoupledMovement
        This is the CoupledMovement with non zero initial attitude angles.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_with_attitude_angles_coupled_operating_point_fixture()
    )

    # Create the CoupledMovement with attitude angles.
    coupled_movement_with_attitude_angles_fixture = (
        ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.1,
            prescribed_num_steps=5,
            free_num_steps=5,
        )
    )

    # Return the CoupledMovement fixture.
    return coupled_movement_with_attitude_angles_fixture


def make_full_coupled_movement_fixture():
    """This method makes a fixture that is a CoupledMovement with all parameters set
    to non zero values for comprehensive testing.

    :return full_coupled_movement_fixture: CoupledMovement
        This is the CoupledMovement with all parameters set to non zero values.
    """
    # Initialize the constructing fixtures.
    airplane_movement = (
        airplane_movement_fixtures.make_multiple_periods_airplane_movement_fixture()
    )
    initial_coupled_operating_point = (
        operating_point_fixtures.make_full_coupled_operating_point_fixture()
    )

    # Create the full CoupledMovement.
    full_coupled_movement_fixture = ps.movements.movement.CoupledMovement(
        airplane_movement=airplane_movement,
        initial_coupled_operating_point=initial_coupled_operating_point,
        delta_time=0.05,
        prescribed_num_steps=10,
        free_num_steps=15,
    )

    # Return the CoupledMovement fixture.
    return full_coupled_movement_fixture


def make_coupled_movement_with_custom_delta_time_fixture():
    """This method makes a fixture that is a CoupledMovement with a custom delta time.

    :return coupled_movement_with_custom_delta_time_fixture: CoupledMovement
        This is the CoupledMovement with custom delta time.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_basic_coupled_operating_point_fixture()
    )

    # Create the CoupledMovement with custom delta time.
    coupled_movement_with_custom_delta_time_fixture = (
        ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.02,
            prescribed_num_steps=5,
            free_num_steps=5,
        )
    )

    # Return the CoupledMovement fixture.
    return coupled_movement_with_custom_delta_time_fixture


def make_coupled_movement_with_more_prescribed_steps_fixture():
    """This method makes a fixture that is a CoupledMovement with more prescribed
    steps than free steps.

    :return coupled_movement_with_more_prescribed_steps_fixture: CoupledMovement
        This is the CoupledMovement with more prescribed steps.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_basic_coupled_operating_point_fixture()
    )

    # Create the CoupledMovement with more prescribed steps.
    coupled_movement_with_more_prescribed_steps_fixture = (
        ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.1,
            prescribed_num_steps=15,
            free_num_steps=3,
        )
    )

    # Return the CoupledMovement fixture.
    return coupled_movement_with_more_prescribed_steps_fixture


def make_coupled_movement_with_more_free_steps_fixture():
    """This method makes a fixture that is a CoupledMovement with more free steps
    than prescribed steps.

    :return coupled_movement_with_more_free_steps_fixture: CoupledMovement
        This is the CoupledMovement with more free steps.
    """
    # Initialize the constructing fixtures.
    airplane_movement = airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    initial_coupled_operating_point = (
        operating_point_fixtures.make_basic_coupled_operating_point_fixture()
    )

    # Create the CoupledMovement with more free steps.
    coupled_movement_with_more_free_steps_fixture = (
        ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.1,
            prescribed_num_steps=3,
            free_num_steps=15,
        )
    )

    # Return the CoupledMovement fixture.
    return coupled_movement_with_more_free_steps_fixture
