# NOTE: I haven't yet started refactoring this module.
"""This module contains useful functions for the movement classes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
from scipy import signal


# NOTE: I haven't yet started refactoring this function.
def oscillating_sinspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D array of values that are calculated by inputting a
    vector of linearly spaced time steps into a sine function.

    :param amplitude: float
        This is the amplitude of the value fluctuation.
    :param period: float
        This is the period of the value fluctuation.
    :param base_value: float
        This is the starting value.
    :param num_steps: int
        This is the number of time steps to iterate through.
    :param delta_time: float
        This is the change in time between each time step.
    :return: 1D array of floats
        This is the resulting vector of sinusoidally spaced values
    """
    # If either the amplitude or the period are 0, return a vector with length equal
    # to the number of steps, and all the values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = 0
    k = base_value

    # Calculate and return the values.
    return a * np.sin(b * (times - h)) + k


# NOTE: I haven't yet started refactoring this function.
def oscillating_linspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D array of values that are calculated by inputting a
    vector of linearly spaced time steps into a triangle function.

    :param amplitude: float
        This is the amplitude of the value fluctuation.
    :param period: float
        This is the period of the value fluctuation.
    :param base_value: float
        This is the starting value.
    :param num_steps: int
        This is the number of time steps to iterate through.
    :param delta_time: float
        This is the change in time between each time step.
    :return: 1D array of floats
        This is the resulting vector of uniformly spaced values
    """
    # If either the amplitude or the period are 0, return a vector with length equal
    # to the number of steps, and all the values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = np.pi / 2
    k = base_value

    # Calculate and return the values.
    return a * signal.sawtooth((b * times + h), 0.5) + k


# NOTE: I haven't yet started refactoring this function.
def oscillating_customspace(
    amplitude, period, base_value, num_steps, delta_time, custom_function
):
    """This function returns a 1D array of values that are calculated by inputting a
    vector of linearly spaced time steps into a custom function.

    :param amplitude: float
        This is the amplitude of the value fluctuation.
    :param period: float
        This is the period of the value fluctuation.
    :param base_value: float
        This is the starting value.
    :param num_steps: int
        This is the number of time steps to iterate through.
    :param delta_time: float
        This is the change in time between each time step.
    :param custom_function: function
        This is a custom function used to return the values. For example, it could be
        np.cos or np.sinh (assuming numpy had previously been imported as np). It
        will be horizontally scaled by the period, vertically scaled by the
        amplitude. For example, say the function has an internal amplitude of 2
        units, an internal period of 3 units, amplitude is set to 4 units and period
        is set to 5 units. The result will have a net amplitude of 8 units and a net
        period of 15 units.
    :return: 1D array of floats
        This is the resulting vector of custom spaced values
    """
    # If either the amplitude or the period are 0, return a vector with length equal
    # to the number of steps, and all the values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = 0
    k = base_value

    # Calculate and return the values.
    return a * custom_function(b * (times - h)) + k
