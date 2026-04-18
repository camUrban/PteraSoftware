"""This module contains functions to create CoreMovements for use in tests."""

# noinspection PyProtectedMember
from pterasoftware._core import CoreMovement

from . import core_airplane_movement_fixtures, core_operating_point_movement_fixtures


def make_static_core_movement_fixture():
    """This method makes a fixture that is a CoreMovement with all static components.

    :return static_core_movement_fixture: CoreMovement
        This is the CoreMovement with no motion.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        core_airplane_movement_fixtures.make_static_core_airplane_movement_fixture()
    ]
    operating_point_movement = (
        core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
    )

    # Create the static CoreMovement.
    static_core_movement_fixture = CoreMovement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        delta_time=0.01,
        num_steps=5,
    )

    # Return the CoreMovement fixture.
    return static_core_movement_fixture


def make_basic_core_movement_fixture():
    """This method makes a fixture that is a CoreMovement with general-purpose values.

    :return basic_core_movement_fixture: CoreMovement
        This is the CoreMovement with general-purpose values for testing.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        core_airplane_movement_fixtures.make_basic_core_airplane_movement_fixture()
    ]
    operating_point_movement = (
        core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
    )

    # Create the basic CoreMovement.
    basic_core_movement_fixture = CoreMovement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        delta_time=0.01,
        num_steps=50,
    )

    # Return the CoreMovement fixture.
    return basic_core_movement_fixture
