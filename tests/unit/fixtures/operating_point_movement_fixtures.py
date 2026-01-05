"""This module contains functions to create OperatingPointMovements for use in
tests."""

import numpy as np

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


def make_uniform_spacing_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with uniform
    spacing for vCg__E.

    :return uniform_spacing_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with uniform spacing for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the uniform spacing OperatingPointMovement.
    uniform_spacing_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=10.0,
            periodVCg__E=1.0,
            spacingVCg__E="uniform",
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return uniform_spacing_operating_point_movement_fixture


def make_phase_offset_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with
    non-zero phase offset for vCg__E.

    :return phase_offset_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with phase offset for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the phase offset OperatingPointMovement.
    phase_offset_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=20.0,
            periodVCg__E=2.0,
            spacingVCg__E="sine",
            phaseVCg__E=90.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return phase_offset_operating_point_movement_fixture


def make_custom_spacing_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with a
    custom spacing function for vCg__E.

    :return custom_spacing_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with custom spacing for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Define a custom harmonic spacing function.
    def custom_harmonic(x):
        """Custom harmonic spacing function: normalized combination of harmonics.

        This function satisfies all requirements: starts at 0, returns to 0 at
        2*pi, has zero mean, has amplitude of 1, and is periodic.

        :param x: (N,) ndarray of floats
            The input angles in radians.

        :return: (N,) ndarray of floats
            The output values.
        """
        return (3.0 / (2.0 * np.sqrt(2.0))) * (
            np.sin(x) + (1.0 / 3.0) * np.sin(3.0 * x)
        )

    # Create the custom spacing OperatingPointMovement.
    custom_spacing_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=15.0,
            periodVCg__E=1.5,
            spacingVCg__E=custom_harmonic,
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return custom_spacing_operating_point_movement_fixture


def make_basic_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with
    general-purpose moderate values.

    :return basic_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with general-purpose values.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the basic OperatingPointMovement.
    basic_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=5.0,
            periodVCg__E=2.0,
            spacingVCg__E="sine",
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return basic_operating_point_movement_fixture


def make_large_amplitude_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with large
    amplitude relative to base speed.

    :return large_amplitude_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with large amplitude for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the large amplitude OperatingPointMovement.
    large_amplitude_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=50.0,
            periodVCg__E=1.0,
            spacingVCg__E="sine",
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return large_amplitude_operating_point_movement_fixture


def make_long_period_operating_point_movement_fixture():
    """This method makes a fixture that is an OperatingPointMovement with a long
    period.

    :return long_period_operating_point_movement_fixture: OperatingPointMovement
        This is the OperatingPointMovement with a long period for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the long period OperatingPointMovement.
    long_period_operating_point_movement_fixture = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=3.0,
            periodVCg__E=10.0,
            spacingVCg__E="sine",
            phaseVCg__E=0.0,
        )
    )

    # Return the OperatingPointMovement fixture.
    return long_period_operating_point_movement_fixture
