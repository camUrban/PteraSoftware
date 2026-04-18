"""Contains useful functions for the movement classes."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import scipy.signal as sp_sig


def oscillating_sin_at_time(
    amp: float,
    period: float,
    phase: float,
    base: float,
    time: float,
) -> float:
    """Returns the result of a customizable sine function evaluated at a time.

    **Note:**

    This function doesn't perform any input validation on its parameters. The
    requirements for the parameters must be validated before passing them in.

    :param amp: The amplitude of the sine function. It must be non negative and can have
        any units as long as they correspond with the units of base.
    :param period: The period of the sine function. It must be non negative. If 0.0, amp
        must also be 0.0, and the function will return 0.0. Its units are in seconds.
    :param phase: The phase offset of the sine function. It must be in the range
        (-180.0, 180.0]. Positive values correspond to phase lead. If both amp and
        period are 0.0, phase must also be 0.0. Its units are in degrees.
    :param base: The mean value about which the sine function oscillates. It can have
        any units as long as they correspond with the units of amp.
    :param time: The time at which to evaluate the sine function. It must be non
        negative. Its units are in seconds.
    :return: The value of the sine function evaluated at given time. Its units will
        match those of amp and base.
    """
    # Convert the function characteristics into classic wave function constants.
    a = amp
    b = 0.0
    if amp != 0.0:
        b = 2.0 * np.pi / period
    h = np.deg2rad(phase)
    k = base

    return float(a * np.sin(b * time + h) + k)


def oscillating_lin_at_time(
    amp: float,
    period: float,
    phase: float,
    base: float,
    time: float,
) -> float:
    """Returns the result of a customizable triangular wave function evaluated at a
    time.

    **Note:**

    This function doesn't perform any input validation on its parameters. The
    requirements for the parameters must be validated before passing them in.

    :param amp: The amplitude of the triangular wave function. It must be non negative,
        and can have any units as long as they correspond with the units of base.
    :param period: The period of the triangular wave function. It must be non negative.
        If 0.0, amp must also be 0.0, and the function will return 0.0. Its units are in
        seconds.
    :param phase: The phase offset of the triangular wave function. It must be in the
        range (-180.0, 180.0]. Positive values correspond to phase lead. If both amp and
        period are 0.0, phase must also be 0.0. Its units are in degrees.
    :param base: The mean value about which the triangular wave function oscillates. It
        can have any units as long as they correspond with the units of amp.
    :param time: The time at which to evaluate the triangular wave function. It must be
        non negative. Its units are in seconds.
    :return: The value of the triangular wave function evaluated at the given time. Its
        units will match those of amp and base.
    """
    # Convert the function characteristics into classic wave function constants.
    a = amp
    b = 0.0
    if amp != 0.0:
        b = 2.0 * np.pi / period
    h = (np.pi / 2.0) + np.deg2rad(phase)
    k = base

    return float(a * sp_sig.sawtooth((b * time + h), 0.5) + k)


def oscillating_custom_at_time(
    amp: float,
    period: float,
    phase: float,
    base: float,
    time: float,
    custom_function: Callable[[float], float],
) -> float:
    """Returns the result of a custom oscillating function evaluated at a time.

    This function is intended for advanced users. The custom function is validated to
    ensure it meets requirements, but users should thoroughly test their functions
    before use in simulations.

    **Note:**

    This function only performs input validation on the custom_function parameter. The
    requirements for the other parameters must be validated before passing them in.

    **Custom Function Requirements:**

    Must start at 0.0 with f(0.0) = 0.0.

    Must return to 0.0 after one period with f(2.0 * pi) = 0.0.

    Must have amplitude of 1.0, meaning (max - min) / 2.0 = 1.0

    Must be periodic with period 2.0 * pi such that f(x) = f(x + 2.0 * pi)

    Must return finite values only with no NaN or Inf

    Must accept a float as input and return a float

    Functions with non zero mean are allowed but will shift the effective center of
    oscillation away from the base value. This can be useful for creating asymmetric
    motion (e.g., faster upstroke than downstroke in flapping).

    **Parameter Interaction:**

    The custom function is transformed by the amps, periods, phases, and bases
    parameters. The output is calculated as amps * custom_function(2.0 * pi * time /
    periods + deg2rad(phases)) + bases. The amps parameter scales the vertical amplitude
    of the custom function. The periods parameter scales the horizontal period of the
    custom function. The phases parameter shifts the function horizontally in degrees.
    The bases parameter shifts the function vertically.

    :param amp: The amplitude of the custom function. It must be non negative, and can
        have any units as long as they correspond with the units of base.
    :param period: The period of the custom function. It must be non negative. If 0.0,
        amp must also be 0.0, and the function will return 0.0. Its units are in
        seconds.
    :param phase: The phase offset of the custom function. It must be in the range
        (-180.0, 180.0]. Positive values correspond to phase lead. If both amp and
        period are 0.0, phase must also be 0.0. Its units are in degrees.
    :param base: The mean value about which the custom function oscillates. It can have
        any units as long as they correspond with the units of amp.
    :param time: The time at which to evaluate the custom function. It must be non
        negative. Its units are in seconds.
    :param custom_function: A custom oscillating function that defines the waveform
        shape. The function must meet all requirements listed above. It must accept a
        float as input and return a float. The function will be scaled and shifted by
        the amps, periods, phases, and bases parameters. Example valid functions,
        assuming numpy is imported as np, include np.sin for a standard sine wave,
        lambda x: 2.0 * np.sin(x) - np.sin(2.0 * x) for a custom harmonic, or lambda x:
        np.where(x < np.pi, x / np.pi, 2.0 - x / np.pi) for a triangle wave. Custom
        functions are validated before use, and if validation fails, a detailed error
        message will indicate which requirement was not met.
    :return: The value of the custom function evaluated at the given time. Its units
        will match those of amp and base.
    """
    # Validate the custom function before using it.
    _validate_custom_spacing_function(custom_function)

    # Convert the function characteristics into classic wave function constants.
    a = amp
    b = 0.0
    if amp != 0.0:
        b = 2.0 * np.pi / period
    h = np.deg2rad(phase)
    k = base

    # Calculate the output or raise an exception if custom_functions throws.
    try:
        return float(a * custom_function(b * time + h) + k)
    except Exception as e:  # pragma: no cover
        raise ValueError(
            f"Calling your custom_function on the inputs resulted in the following "
            f"exception:\n{e}"
        )


def _validate_custom_spacing_function(
    custom_function: Callable[[float], float],
) -> None:
    """Validates that a custom spacing function meets requirements for use in
    oscillating_custom_at_time.

    See the oscillating_custom_at_time docstring for the exact requirements for the
    custom function.

    :param custom_function: The custom spacing function to validate.
    :return: None
    """
    # Test the function over two full periods. Use an odd number of points so that one
    # point lies exactly on 2.0 * pi.
    test_times = np.linspace(0.0, 4.0 * np.pi, 201, dtype=float)

    test_output = np.zeros_like(test_times)
    try:
        for this_test_time_id, this_test_time in enumerate(test_times):
            test_output[this_test_time_id] = custom_function(this_test_time)
    except Exception as e:
        raise ValueError(
            f"The custom spacing function failed when called with test input: {e}."
        )

    for this_output_id, this_output in enumerate(test_output):
        if not isinstance(this_output, float):
            raise ValueError(
                f"The custom spacing function must return a float for a float input, "
                f"but it returned {type(this_output)} for the input "
                f"{test_times[this_output_id]}."
            )

    # Check for finite values.
    if not np.isfinite(test_output).all():
        raise ValueError(
            "Custom spacing function must return finite values only (no NaN or Inf)."
        )

    # Extract one period of data for validation (first period).
    first_period_indices = test_times < 2.0 * np.pi
    first_period_output = test_output[first_period_indices]

    tolerance = 0.05

    # Check that the function starts at 0.0.
    start_value = test_output[0]
    if not np.isclose(start_value, 0.0, atol=tolerance):
        raise ValueError(
            f"Custom spacing function must start at 0.0. f(0.0) = {start_value:.4f}, "
            f"but should be within {tolerance} of 0.0."
        )

    # Check that the function returns to 0.0 after one period.
    # Find the index closest to 2.0 * pi.
    end_period_idx = np.argmin(np.abs(test_times - 2.0 * np.pi))
    end_value = test_output[end_period_idx]
    if not np.isclose(end_value, 0.0, atol=tolerance):
        raise ValueError(
            f"Custom spacing function must return to 0.0 after one period. "
            f"f(2.0 * pi) = {end_value:.4f}, but should be within {tolerance} of 0.0."
        )

    # Check that the amplitude = 1.0.
    max_value = float(np.max(first_period_output))
    min_value = float(np.min(first_period_output))
    amplitude = (max_value - min_value) / 2.0
    if not np.isclose(amplitude, 1.0, atol=tolerance):
        raise ValueError(
            f"Custom spacing function must have an amplitude of 1.0. "
            f"Amplitude = {amplitude:.4f}, but should be within {tolerance} of 1.0."
        )

    # Check periodicity by comparing the first and second periods.
    second_period_indices = (test_times >= 2.0 * np.pi) & (test_times < 4.0 * np.pi)
    second_period_output = test_output[second_period_indices]

    # They should have the same length if properly sampled.
    if len(first_period_output) == len(second_period_output):
        if not np.allclose(first_period_output, second_period_output, atol=tolerance):
            max_diff = np.max(np.abs(first_period_output - second_period_output))
            raise ValueError(
                f"Custom spacing function must be periodic with period 2.0 * pi. "
                f"Maximum difference between first and second period: {max_diff:.4f}, "
                f"but should be within {tolerance}."
            )
