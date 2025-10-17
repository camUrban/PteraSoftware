"""This module contains useful functions for the movement classes."""

import numpy as np
from scipy import signal

from .. import _parameter_validation


def oscillating_sinspaces(amps, periods, phases, bases, num_steps, delta_time):
    """This function returns a (..., num_steps) ndarray of floats that are calculated
    by inputting a vector of linearly spaced time steps into a sine function defined
    with the parameters given by the scalars or array-like objects amp, period,
    phase, and base.

    :param amps: number or array-like of numbers

        The amplitudes of the fluctuation. It must be a non-negative number (int or
        float) or an array-like object of non-negative numbers. All elements will be
        converted to floats internally. If any of its elements are 0.0, then the
        corresponding periods element must also be 0.0, and the corresponding results
        will have no fluctuations. Its units can be anything so long as they
        correspond with the units of base.

    :param periods: number or array-like of numbers

        The periods of the fluctuations. It must be a non-negative number (int or
        float) or an array-like object of non-negative numbers. All elements will be
        converted to floats internally. If any of its elements are 0.0, then the
        corresponding amps element must also be 0.0, and the corresponding results
        will have no fluctuations. If an array-like object, its shape must match that
        of amps. Its units are in seconds.

    :param phases: number or array-like of numbers

        The phase offsets of the fluctuation. It must be a number (int or float),
        or an array-like object of numbers, in the range (-180.0, 180.0], and will be
        converted to a float internally. Positive values correspond to phase lead. If
        a given result has no fluctuations (corresponding elements in amps and
        periods are 0.0), the corresponding element in phases must be 0.0. If an
        array-like object, its shape must match that of amps. Its units are in degrees.

    :param bases: number or array-like of numbers

        The mean values about which the fluctuation occurs. It must be a number (int
        or float), or array-like object of numbers. All elements will be converted to
        floats internally. If an array-like object, its shape must match that of
        amps. Its units can be anything so long as they correspond with the units of
        amps.

    :param num_steps: int

        The number of time steps to iterate through. It must be a positive int.

    :param delta_time: number

        The change in time between each time step. It must be a positive number (int
        or float), and will be converted to a float internally. Its units are in
        seconds.

    :return: (num_steps,) or (..., num_steps) ndarray of floats

        The resulting ndarray of sinusoidally varying values. It will be a ndarray of
        floats with shape (num_steps,) (for scalar parameters) or (S, num_steps) (for
        array-like parameters of shape S). Its units will match those of amp and base.
    """
    amps, periods, phases, bases, num_steps, delta_time, mask_static = (
        _validate_oscillating_function_parameters(
            amps, periods, phases, bases, num_steps, delta_time
        )
    )

    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristic ndarrays into ndarrays of classic wave
    # function constants. This also adds a trailing dimension of length 1 to each
    # ndarray, so that broadcasting can occur when calculating the results.
    a = amps[..., None]

    b = np.zeros_like(periods, dtype=float)
    b[~mask_static] = 2 * np.pi / periods[~mask_static]
    b = b[..., None]

    h = np.deg2rad(phases)[..., None]
    k = bases[..., None]

    return a * np.sin(b * times + h) + k


def oscillating_linspaces(amps, periods, phases, bases, num_steps, delta_time):
    """This function returns a (..., num_steps) ndarray of floats that are calculated
    by inputting a vector of linearly spaced time steps into a triangular wave
    function defined with the parameters given by the scalars or array-like objects
    amp, period, phase, and base.

    :param amps: number or array-like of numbers

        The amplitudes of the fluctuation. It must be a non-negative number (int or
        float) or an array-like object of non-negative numbers. All elements will be
        converted to floats internally. If any of its elements are 0.0, then the
        corresponding periods element must also be 0.0, and the corresponding results
        will have no fluctuations. Its units can be anything so long as they
        correspond with the units of base.

    :param periods: number or array-like of numbers

        The periods of the fluctuations. It must be a non-negative number (int or
        float) or an array-like object of non-negative numbers. All elements will be
        converted to floats internally. If any of its elements are 0.0, then the
        corresponding amps element must also be 0.0, and the corresponding results
        will have no fluctuations. If an array-like object, its shape must match that
        of amps. Its units are in seconds.

    :param phases: number or array-like of numbers

        The phase offsets of the fluctuation. It must be a number (int or float),
        or an array-like object of numbers, in the range (-180.0, 180.0], and will be
        converted to a float internally. Positive values correspond to phase lead. If
        a given result has no fluctuations (corresponding elements in amps and
        periods are 0.0), the corresponding element in phases must be 0.0. If an
        array-like object, its shape must match that of amps. Its units are in degrees.

    :param bases: number or array-like of numbers

        The mean values about which the fluctuation occurs. It must be a number (int
        or float), or array-like object of numbers. All elements will be converted to
        floats internally. If an array-like object, its shape must match that of
        amps. Its units can be anything so long as they correspond with the units of
        amps.

    :param num_steps: int

        The number of time steps to iterate through. It must be a positive int.

    :param delta_time: number

        The change in time between each time step. It must be a positive number (int
        or float), and will be converted to a float internally. Its units are in
        seconds.

    :return: (num_steps,) or (..., num_steps) ndarray of floats

        The resulting ndarray of varying values. It will be a ndarray of floats with
        shape (num_steps,) (for scalar parameters) or (S, num_steps) (for array-like
        parameters of shape S). Its units will match those of amp and base.
    """
    amps, periods, phases, bases, num_steps, delta_time, mask_static = (
        _validate_oscillating_function_parameters(
            amps, periods, phases, bases, num_steps, delta_time
        )
    )

    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristic ndarrays into ndarrays of classic wave
    # function constants. This also adds a trailing dimension of length 1 to each
    # ndarray, so that broadcasting can occur when calculating the results.
    a = amps[..., None]

    b = np.zeros_like(periods, dtype=float)
    b[~mask_static] = 2 * np.pi / periods[~mask_static]
    b = b[..., None]

    h = (np.pi / 2) + np.deg2rad(phases)[..., None]
    k = bases[..., None]

    # Calculate and return the values.
    return a * signal.sawtooth((b * times + h), 0.5) + k


def oscillating_customspaces(
    amps, periods, phases, bases, num_steps, delta_time, custom_function
):
    """This function returns a (..., num_steps) ndarray of floats that are calculated
    by inputting a vector of linearly spaced time steps into a custom oscillating
    function defined with the parameters given by the scalars or array-like objects
    amp, period, phase, and base.

    Note: This function is intended for advanced users. The custom function is
    validated to ensure it meets requirements, but users should thoroughly test their
    functions before use in simulations.

    Custom Function Requirements:
        The function must start at 0 with f(0) = 0, and must return to 0 after one
        period with f(2*pi) = 0. The function must have amplitude of 1, meaning (max -
        min) / 2 = 1.0. The function must be periodic with period 2*pi such that f(x) =
        f(x + 2*pi). The function must return finite values only with no NaN or Inf.
        The function must accept a ndarray as input and return a ndarray of the same
        shape.

        Note: Functions with non-zero mean are allowed but will shift the effective 
        center of oscillation away from the base value. This can be useful for 
        creating asymmetric motion (e.g., faster upstroke than downstroke in flapping).

    Parameter Interaction:
        The custom function is transformed by the amps, periods, phases, and bases
        parameters. The output is calculated as amps * custom_function(2*pi * time /
        periods + deg2rad(phases)) + bases. The amps parameter scales the vertical
        amplitude of the custom function. The periods parameter scales the horizontal
        period of the custom function. The phases parameter shifts the function
        horizontally in degrees. The bases parameter shifts the function vertically.

    :param amps: number or array-like of numbers

        The amplitudes of the fluctuation. It must be a non-negative number (int or
        float) or an array-like object of non-negative numbers. All elements will be
        converted to floats internally. If any of its elements are 0.0, then the
        corresponding periods element must also be 0.0, and the corresponding results
        will have no fluctuations. Its units can be anything so long as they
        correspond with the units of base.

    :param periods: number or array-like of numbers

        The periods of the fluctuations. It must be a non-negative number (int or
        float) or an array-like object of non-negative numbers. All elements will be
        converted to floats internally. If any of its elements are 0.0, then the
        corresponding amps element must also be 0.0, and the corresponding results
        will have no fluctuations. If an array-like object, its shape must match that
        of amps. Its units are in seconds.

    :param phases: number or array-like of numbers

        The phase offsets of the fluctuation. It must be a number (int or float),
        or an array-like object of numbers, in the range (-180.0, 180.0], and will be
        converted to a float internally. Positive values correspond to phase lead. If
        a given result has no fluctuations (corresponding elements in amps and
        periods are 0.0), the corresponding element in phases must be 0.0. If an
        array-like object, its shape must match that of amps. Its units are in degrees.

    :param bases: number or array-like of numbers

        The mean values about which the fluctuation occurs. It must be a number (int
        or float), or array-like object of numbers. All elements will be converted to
        floats internally. If an array-like object, its shape must match that of
        amps. Its units can be anything so long as they correspond with the units of
        amps.

    :param num_steps: int

        The number of time steps to iterate through. It must be a positive int.

    :param delta_time: number

        The change in time between each time step. It must be a positive number (int
        or float), and will be converted to a float internally. Its units are in
        seconds.

    :param custom_function: callable

        A custom oscillating function that defines the waveform shape. The function
        must meet all requirements listed above. It must accept a ndarray as input
        and return a ndarray of the same shape. The function will be scaled and
        shifted by the amps, periods, phases, and bases parameters. Example valid
        functions, assuming numpy is imported as np, include np.sin for a standard
        sine wave, lambda x: 2 * np.sin(x) - np.sin(2 * x) for a custom harmonic,
        or lambda x: np.where(x < np.pi, x / np.pi, 2 - x / np.pi) for a triangle
        wave. Custom functions are validated before use, and if validation fails,
        a detailed error message will indicate which requirement was not met.

    :return: (num_steps,) or (..., num_steps) ndarray of floats

        The resulting ndarray of varying values. It will be a ndarray of floats with
        shape (num_steps,) (for scalar parameters) or (S, num_steps) (for array-like
        parameters of shape S). Its units will match those of amp and base.
    """
    amps, periods, phases, bases, num_steps, delta_time, mask_static = (
        _validate_oscillating_function_parameters(
            amps, periods, phases, bases, num_steps, delta_time
        )
    )

    # Validate the custom function before using it.
    _validate_custom_spacing_function(custom_function)

    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristic ndarrays into ndarrays of classic wave
    # function constants. This also adds a trailing dimension of length 1 to each
    # ndarray, so that broadcasting can occur when calculating the results.
    a = amps[..., None]

    b = np.zeros_like(periods, dtype=float)
    b[~mask_static] = 2 * np.pi / periods[~mask_static]
    b = b[..., None]

    h = np.deg2rad(phases)[..., None]
    k = bases[..., None]

    # Calculate the output or raise an exception if custom_functions throws.
    try:
        output = np.asarray(a * custom_function(b * times + h) + k)
    except Exception as e:
        raise ValueError(
            f"Calling your custom_function on the inputs resulted in the following exception:\n{e}"
        )

    output_shape = output.shape
    expected_shape = amps.shape + (num_steps,)

    if output_shape != expected_shape:
        raise ValueError(
            f"Calling custom_function on your inputs resulted in an ndarray of shape {output_shape}, but the expected shape is {expected_shape}."
        )
    return output


def _validate_oscillating_function_parameters(
    amps, periods, phases, bases, num_steps, delta_time
):
    """Validates and returns the conditioned parameters for the oscillating_*
    functions. See their docstrings for details on the requirements for the
    parameters. It also returns the array mask for identifying static cases.
    """
    amps = _parameter_validation.arrayLike_of_numbers_in_range_return_float(
        amps, "amps", 0.0, True, None, None
    )
    periods = _parameter_validation.arrayLike_of_numbers_in_range_return_float(
        periods, "periods", 0.0, True, None, None
    )
    phases = _parameter_validation.arrayLike_of_numbers_in_range_return_float(
        phases, "phases", -180.0, False, 180.0, True
    )
    bases = _parameter_validation.arrayLike_of_numbers_in_range_return_float(
        bases, "bases", None, None, None, None
    )

    expected_shape = amps.shape
    values_to_check = (periods, phases, bases)
    value_names = ("periods", "phases", "bases")
    for value_id in range(len(values_to_check)):
        value = values_to_check[value_id]
        value_name = value_names[value_id]

        value_shape = value.shape
        if value_shape != expected_shape:
            raise ValueError(
                f"After conversion to a ndarray, {value_name} must have the same shape as amps, which is {expected_shape}, but its shape is {value_shape}."
            )

    mask_invalid = (amps == 0.0) ^ (periods == 0.0)
    if np.any(mask_invalid):
        raise ValueError(
            "If an element in amps is 0.0, the corresponding element in periods must also be 0.0."
        )
    mask_static = amps == 0.0
    if np.any(mask_static & (phases != 0.0)):
        raise ValueError(
            "If the elements at a given location in amps and periods are 0.0, the corresponding element in phases must also be 0.0."
        )

    num_steps = _parameter_validation.positive_int_return_int(num_steps, "num_steps")
    delta_time = _parameter_validation.positive_number_return_float(
        delta_time, "delta_time"
    )

    return [amps, periods, phases, bases, num_steps, delta_time, mask_static]


def _validate_custom_spacing_function(custom_function):
    """Validates that a custom spacing function meets requirements for use in
    oscillating_customspaces.

    The function must start at 0 with f(0) approximately equal to 0, and return to
    0 after one period with f(2*pi) approximately equal to 0. The function must have
    amplitude of 1 with (max - min) / 2 approximately equal to 1.0. The function must
    be periodic such that f(x) is approximately equal to f(x + 2*pi). The function
    must return finite values only.

    :param custom_function: callable
        The custom spacing function to validate.

    :raises ValueError: If the function doesn't meet the requirements.
    """
    # Test the function over two full periods. Us an odd number of points so that one
    # lies exactly on 2*pi.
    test_input = np.linspace(0, 4 * np.pi, 201)

    try:
        test_output = custom_function(test_input)
    except Exception as e:
        raise ValueError(
            f"Custom spacing function failed when called with test input: {e}"
        )

    # Convert to ndarray and check shape.
    test_output = np.asarray(test_output)
    if test_output.shape != test_input.shape:
        raise ValueError(
            f"Custom spacing function must return a ndarray of the same shape as its input. "
            f"Input shape: {test_input.shape}, output shape: {test_output.shape}."
        )

    # Check for finite values.
    if not np.isfinite(test_output).all():
        raise ValueError(
            "Custom spacing function must return finite values only (no NaN or Inf)."
        )

    # Extract one period of data for validation (first period).
    first_period_indices = test_input < 2 * np.pi
    first_period_output = test_output[first_period_indices]

    tolerance = 0.05

    # Check that function starts at 0.
    start_value = test_output[0]
    if not np.isclose(start_value, 0.0, atol=tolerance):
        raise ValueError(
            f"Custom spacing function must start at 0. f(0) = {start_value:.4f}, "
            f"but should be within {tolerance} of 0."
        )

    # Check that function returns to 0 after one period.
    # Find the index closest to 2*pi.
    end_period_idx = np.argmin(np.abs(test_input - 2 * np.pi))
    end_value = test_output[end_period_idx]
    if not np.isclose(end_value, 0.0, atol=tolerance):
        raise ValueError(
            f"Custom spacing function must return to 0 after one period. "
            f"f(2*pi) = {end_value:.4f}, but should be within {tolerance} of 0."
        )

    # Check amplitude = 1.
    max_value = float(np.max(first_period_output))
    min_value = float(np.min(first_period_output))
    amplitude = (max_value - min_value) / 2.0
    if not np.isclose(amplitude, 1.0, atol=tolerance):
        raise ValueError(
            f"Custom spacing function must have amplitude of 1. "
            f"Amplitude = {amplitude:.4f}, but should be within {tolerance} of 1."
        )

    # Check periodicity by comparing first and second periods.
    second_period_indices = (test_input >= 2 * np.pi) & (test_input < 4 * np.pi)
    second_period_output = test_output[second_period_indices]

    # They should have the same length if properly sampled.
    if len(first_period_output) == len(second_period_output):
        if not np.allclose(first_period_output, second_period_output, atol=tolerance):
            max_diff = np.max(np.abs(first_period_output - second_period_output))
            raise ValueError(
                f"Custom spacing function must be periodic with period 2*pi. "
                f"Maximum difference between first and second period: {max_diff:.4f}, "
                f"but should be within {tolerance}."
            )
