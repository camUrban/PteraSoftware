"""This module contains functions to create fixtures for movements functions tests."""

import numpy as np


def make_scalar_parameters_fixture():
    """This method makes a fixture with scalar parameters for testing oscillating
    functions.

    :return tuple of scalars
        This returns a tuple containing scalar values for amps, periods, phases, and
        bases.
    """
    amps = 1.0
    periods = 1.0
    phases = 0.0
    bases = 0.0

    return amps, periods, phases, bases


def make_array_parameters_fixture():
    """This method makes a fixture with array parameters for testing oscillating
    functions.

    :return tuple of (3,) ndarrays of floats
        This returns a tuple containing array values for amps, periods, phases, and
        bases.
    """
    amps = np.array([1.0, 2.0, 0.5], dtype=float)
    periods = np.array([1.0, 2.0, 0.5], dtype=float)
    phases = np.array([0.0, 90.0, -45.0], dtype=float)
    bases = np.array([0.0, 1.0, -0.5], dtype=float)

    return amps, periods, phases, bases


def make_static_parameters_fixture():
    """This method makes a fixture with static parameters for testing oscillating
    functions.

    :return tuple of scalars
        This returns a tuple containing scalar values representing static motion where
        amps and periods are both 0.0.
    """
    amps = 0.0
    periods = 0.0
    phases = 0.0
    bases = 2.0

    return amps, periods, phases, bases


def make_mixed_static_parameters_fixture():
    """This method makes a fixture with mixed static and dynamic parameters for testing
    oscillating functions.

    :return tuple of (3,) ndarrays of floats
        This returns a tuple containing array values where some elements are static
        and some are dynamic.
    """
    amps = np.array([1.0, 0.0, 0.5], dtype=float)
    periods = np.array([1.0, 0.0, 0.5], dtype=float)
    phases = np.array([0.0, 0.0, -45.0], dtype=float)
    bases = np.array([0.0, 1.0, -0.5], dtype=float)

    return amps, periods, phases, bases


def make_phase_offset_parameters_fixture():
    """This method makes a fixture with phase-offset parameters for testing oscillating
    functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with a non-zero phase.
    """
    amps = 1.0
    periods = 1.0
    phases = 90.0
    bases = 0.0

    return amps, periods, phases, bases


def make_large_amplitude_parameters_fixture():
    """This method makes a fixture with large amplitude parameters for testing
    oscillating functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with a large amplitude.
    """
    amps = 10.0
    periods = 1.0
    phases = 0.0
    bases = 5.0

    return amps, periods, phases, bases


def make_small_period_parameters_fixture():
    """This method makes a fixture with small period parameters for testing oscillating
    functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with a small period.
    """
    amps = 1.0
    periods = 0.1
    phases = 0.0
    bases = 0.0

    return amps, periods, phases, bases


def make_negative_phase_parameters_fixture():
    """This method makes a fixture with negative phase parameters for testing
    oscillating functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with a negative phase.
    """
    amps = 1.0
    periods = 1.0
    phases = -90.0
    bases = 0.0

    return amps, periods, phases, bases


def make_max_phase_parameters_fixture():
    """This method makes a fixture with maximum phase parameters for testing
    oscillating functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with phase at the maximum
        allowed value.
    """
    amps = 1.0
    periods = 1.0
    phases = 180.0
    bases = 0.0

    return amps, periods, phases, bases


def make_min_phase_parameters_fixture():
    """This method makes a fixture with minimum phase parameters for testing
    oscillating functions.

    :return tuple of scalars
        This returns a tuple containing scalar values with phase just above the
        minimum allowed value.
    """
    amps = 1.0
    periods = 1.0
    phases = -179.9
    bases = 0.0

    return amps, periods, phases, bases


def make_num_steps_and_delta_time_fixture():
    """This method makes a fixture for num_steps and delta_time parameters.

    :return tuple of int and float
        This returns a tuple containing num_steps as an int and delta_time as a
        float.
    """
    num_steps = 100
    delta_time = 0.01

    return num_steps, delta_time


def make_small_num_steps_fixture():
    """This method makes a fixture for small num_steps parameter.

    :return int
        This returns num_steps as a small positive int.
    """
    num_steps = 10

    return num_steps


def make_large_num_steps_fixture():
    """This method makes a fixture for large num_steps parameter.

    :return int
        This returns num_steps as a large positive int.
    """
    num_steps = 1000

    return num_steps


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


def make_valid_custom_triangle_function_fixture():
    """This method makes a fixture that is a valid custom triangle function for testing
    oscillating_customspaces.

    :return callable
        This returns a valid custom triangle function that satisfies all requirements
        for oscillating_customspaces.
    """

    def custom_triangle(x):
        return 2.0 / np.pi * np.arcsin(np.sin(x))

    return custom_triangle


def make_valid_custom_harmonic_function_fixture():
    """This method makes a fixture that is a valid custom harmonic function for testing
    oscillating_customspaces.

    :return callable
        This returns a valid custom harmonic function that satisfies all requirements
        for oscillating_customspaces.
    """

    def custom_harmonic(x):
        # Create harmonic with fundamental and second harmonic.
        raw = np.sin(x) + 0.5 * np.sin(2 * x)
        # Normalize to have amplitude of 1.
        raw_min = float(np.min(raw))
        raw_max = float(np.max(raw))
        amplitude = (raw_max - raw_min) / 2.0
        return raw / amplitude

    return custom_harmonic


def make_invalid_custom_function_wrong_start_fixture():
    """This method makes a fixture that is an invalid custom function that doesn't
    start at 0.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_start(x):
        return np.sin(x) + 0.5

    return invalid_start


def make_invalid_custom_function_wrong_end_fixture():
    """This method makes a fixture that is an invalid custom function that doesn't
    return to 0 after one period.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_end(x):
        return np.sin(x * 1.01)

    return invalid_end


def make_invalid_custom_function_wrong_mean_fixture():
    """This method makes a fixture that is an invalid custom function with non-zero
    mean.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_mean(x):
        return np.sin(x) + 0.1 * abs(np.sin(2 * x))

    return invalid_mean


def make_invalid_custom_function_wrong_amplitude_fixture():
    """This method makes a fixture that is an invalid custom function with amplitude
    not equal to 1.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_amplitude(x):
        return 2.0 * np.sin(x)

    return invalid_amplitude


def make_invalid_custom_function_not_periodic_fixture():
    """This method makes a fixture that is an invalid custom function that is not
    periodic.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_periodic(x):
        return np.sin(x) * (1.0 + 0.01 * x)

    return invalid_periodic


def make_invalid_custom_function_returns_nan_fixture():
    """This method makes a fixture that is an invalid custom function that returns NaN.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_nan(x):
        result = np.sin(x)
        result[len(result) // 2] = np.nan
        return result

    return invalid_nan


def make_invalid_custom_function_returns_inf_fixture():
    """This method makes a fixture that is an invalid custom function that returns Inf.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_inf(x):
        result = np.sin(x)
        result[len(result) // 2] = np.inf
        return result

    return invalid_inf


def make_invalid_custom_function_wrong_shape_fixture():
    """This method makes a fixture that is an invalid custom function that returns the
    wrong shape.

    :return callable
        This returns an invalid custom function for testing validation.
    """

    def invalid_shape(x):
        return np.sin(x)[:-1]

    return invalid_shape
