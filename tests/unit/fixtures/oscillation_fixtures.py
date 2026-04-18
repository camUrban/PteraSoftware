"""This module contains functions to create fixtures for oscillation tests."""

import numpy as np


def make_static_parameters_fixture():
    """This method makes a fixture with static parameters for testing oscillating
    functions.

    :return tuple of scalars
        This returns a tuple containing scalar values representing static motion where
        amp and period are both 0.0.
    """
    amp = 0.0
    period = 0.0
    phase = 0.0
    base = 2.0

    return amp, period, phase, base


def make_phase_offset_parameters_fixture():
    """This method makes a fixture with a non zero phase offset parameter for testing
    oscillating functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with a non zero phase.
    """
    amp = 1.0
    period = 1.0
    phase = 90.0
    base = 0.0

    return amp, period, phase, base


def make_max_phase_parameters_fixture():
    """This method makes a fixture with maximum phase parameters for testing
    oscillating functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with phase at the maximum
        allowed value.
    """
    amp = 1.0
    period = 1.0
    phase = 180.0
    base = 0.0

    return amp, period, phase, base


def make_time_fixture():
    """This method makes a fixture for a time parameter.

    :return float
        This returns a time.
    """
    time = 5.32

    return time


def make_valid_custom_sine_function_fixture():
    """This method makes a fixture that is a valid custom sine function for testing
    oscillating_customspaces.

    :return callable
        This returns a valid custom function that satisfies all requirements for
        oscillating_customspaces.
    """

    def custom_sine(x):
        return np.sin(x)

    return custom_sine
