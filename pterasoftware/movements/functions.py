# TODO: Consider making this module private.
"""This module contains useful functions for the movement classes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:

    oscillating_sinspaces: This function returns a (..., num_steps) ndarray of floats
    that are calculated by inputting a vector of linearly spaced time steps into one
    or more sine functions defined with the parameters given by the scalars or
    array-like objects amp, period, phase, and base.

    oscillating_linspace: This function returns a (..., num_steps) ndarray of floats
    that are calculated by inputting a vector of linearly spaced time steps into one
    or more triangular wave functions defined with the parameters given by the
    scalars or array-like objects amp, period, phase, and base.

    oscillating_customspace: This function returns a (..., num_steps) ndarray of
    floats that are calculated by inputting a vector of linearly spaced time steps
    into a custom oscillating function defined with the parameters given by the
    scalars or array-like objects amp, period, phase, and base."""

import numpy as np
from scipy import signal

from .. import parameter_validation


# TODO: Create unit tests for this function.
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


# TODO: Create unit tests for this function.
def oscillating_linspace(amps, periods, phases, bases, num_steps, delta_time):
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


# TODO: Create unit tests for this function.
def oscillating_customspace(
    amps, periods, bases, phases, num_steps, delta_time, custom_function
):
    """This function returns a (..., num_steps) ndarray of floats that are calculated
    by inputting a vector of linearly spaced time steps into a custom oscillating
    function defined with the parameters given by the scalars or array-like objects
    amp, period, phase, and base.

    Note: This function performs vary basic validation on custom_function,
    but it could fail in unexpected ways. It is intended for advanced users.

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

    :param custom_function: function

        This is a custom oscillating function used to return the values. For example,
        it could be np.cos or np.sinh (assuming numpy had previously been imported as
        np). It will be horizontally scaled by periods and vertically scaled by amps.
        For example, say the function has an internal amplitude of 2 units,
        an internal period of 3 units, a particular location in amps is set to 4
        units and in periods to 5 units. The result will have a net amplitude of 8
        units and a net period of 15 units. It will also be shifted vertically by the
        bases and horizontally by phases. The function must take a ndarray as an
        input, operate elementwise, and return a ndarray of the same shape (or a
        scalar if the inputs are scalars/0D ndarrays).

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
    amps = parameter_validation.arrayLike_of_numbers_in_range_return_float(
        amps, "amps", 0.0, True, None, None
    )
    periods = parameter_validation.arrayLike_of_numbers_in_range_return_float(
        periods, "periods", 0.0, True, None, None
    )
    phases = parameter_validation.arrayLike_of_numbers_in_range_return_float(
        phases, "phases", -180.0, False, 180.0, True
    )
    bases = parameter_validation.arrayLike_of_numbers_in_range_return_float(
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

    num_steps = parameter_validation.positive_int_return_int(num_steps, "num_steps")
    delta_time = parameter_validation.positive_number_return_float(
        delta_time, "delta_time"
    )

    return [amps, periods, phases, bases, num_steps, delta_time, mask_static]
