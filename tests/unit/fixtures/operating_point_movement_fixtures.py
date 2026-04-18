"""This module contains functions to create OperatingPointMovements for use in
tests."""

import pterasoftware as ps

from . import operating_point_fixtures


def make_static_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with all
    parameters zero (no movement).

    :return static_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with no movement.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the static OperatingPointMovement.
    static_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=0.0,
            periodVCg__E=0.0,
            spacingVCg__E="sine",
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return static_operating_point_movement_fixture


def make_sine_spacing_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with sine
    spacing for vCg__E.

    :return sine_spacing_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with sine spacing for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the sine spacing OperatingPointMovement.
    sine_spacing_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=10.0,
            periodVCg__E=1.0,
            spacingVCg__E="sine",
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return sine_spacing_operating_point_movement_fixture
