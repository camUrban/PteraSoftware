"""This module contains functions to create AeroelasticOperatingPointMovements for use
in tests."""

import pterasoftware as ps

from . import operating_point_fixtures


def make_sine_spacing_aeroelastic_operating_point_movement_fixture():
    """This method makes a fixture that is an AeroelasticOperatingPointMovement with
    sine spacing for vCg__E.

    :return sine_spacing_aeroelastic_operating_point_movement_fixture:
        AeroelasticOperatingPointMovement
        This is the AeroelasticOperatingPointMovement with sine spacing for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the sine spacing AeroelasticOperatingPointMovement.
    sine_spacing_aeroelastic_operating_point_movement_fixture = ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=10.0,
        periodVCg__E=1.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )

    # Return the AeroelasticOperatingPointMovement fixture.
    return sine_spacing_aeroelastic_operating_point_movement_fixture
