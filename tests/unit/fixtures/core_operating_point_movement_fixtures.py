"""This module contains functions to create CoreOperatingPointMovements for use in
tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware._core import CoreOperatingPointMovement

from . import operating_point_fixtures


def make_static_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with all
    parameters zero (no movement).

    :return static_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with no movement.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the static CoreOperatingPointMovement.
    static_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return static_core_operating_point_movement_fixture


def make_sine_spacing_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with sine
    spacing for vCg__E.

    :return sine_spacing_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with sine spacing for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the sine spacing CoreOperatingPointMovement.
    sine_spacing_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=10.0,
        periodVCg__E=1.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return sine_spacing_core_operating_point_movement_fixture


def make_uniform_spacing_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with uniform
    spacing for vCg__E.

    :return uniform_spacing_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with uniform spacing for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the uniform spacing CoreOperatingPointMovement.
    uniform_spacing_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=10.0,
        periodVCg__E=1.0,
        spacingVCg__E="uniform",
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return uniform_spacing_core_operating_point_movement_fixture


def make_phase_offset_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with
    non-zero phase offset for vCg__E.

    :return phase_offset_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with phase offset for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the phase offset CoreOperatingPointMovement.
    phase_offset_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=20.0,
        periodVCg__E=2.0,
        spacingVCg__E="sine",
        phaseVCg__E=90.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return phase_offset_core_operating_point_movement_fixture


def make_custom_spacing_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with a
    custom spacing function for vCg__E.

    :return custom_spacing_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with custom spacing for vCg__E.
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

    # Create the custom spacing CoreOperatingPointMovement.
    custom_spacing_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=15.0,
        periodVCg__E=1.5,
        spacingVCg__E=custom_harmonic,
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return custom_spacing_core_operating_point_movement_fixture


def make_basic_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with
    general-purpose moderate values.

    :return basic_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with general-purpose values.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the basic CoreOperatingPointMovement.
    basic_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=5.0,
        periodVCg__E=2.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return basic_core_operating_point_movement_fixture


def make_large_amplitude_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with large
    amplitude relative to base speed.

    :return large_amplitude_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with large amplitude for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = (
        operating_point_fixtures.make_high_speed_operating_point_fixture()
    )

    # Create the large amplitude CoreOperatingPointMovement.
    large_amplitude_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=50.0,
        periodVCg__E=1.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return large_amplitude_core_operating_point_movement_fixture


def make_long_period_core_operating_point_movement_fixture():
    """This method makes a fixture that is an CoreOperatingPointMovement with a long
    period.

    :return long_period_core_operating_point_movement_fixture: CoreOperatingPointMovement
        This is the CoreOperatingPointMovement with a long period for vCg__E.
    """
    # Initialize the constructing fixture.
    base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the long period CoreOperatingPointMovement.
    long_period_core_operating_point_movement_fixture = CoreOperatingPointMovement(
        base_operating_point=base_operating_point,
        ampVCg__E=3.0,
        periodVCg__E=10.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )

    # Return the CoreOperatingPointMovement fixture.
    return long_period_core_operating_point_movement_fixture
