"""This module contains functions to create FreeFlightOperatingPointMovements for use in
tests."""

import pterasoftware as ps

from . import operating_point_fixtures


def make_basic_free_flight_operating_point_movement_fixture() -> (
    ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement
):
    """This method makes a fixture that is a FreeFlightOperatingPointMovement for
    general testing.

    :return basic_free_flight_operating_point_movement_fixture:
        FreeFlightOperatingPointMovement This is the FreeFlightOperatingPointMovement
        configured for general testing.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the basic FreeFlightOperatingPointMovement.
    basic_free_flight_operating_point_movement_fixture = ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
        base_operating_point=base_operating_point,
    )

    # Return the FreeFlightOperatingPointMovement fixture.
    return basic_free_flight_operating_point_movement_fixture
