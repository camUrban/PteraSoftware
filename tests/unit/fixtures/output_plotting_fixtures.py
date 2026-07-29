"""This module contains functions to create fixtures for the output plotting tests."""

import numpy as np


def make_times_fixture() -> np.ndarray:
    """Makes a fixture that is the time at each of a short time history's five steps.

    :return: A (5,) ndarray of floats representing the time, in seconds, at each time
        step.
    """
    return np.array([0.0, 0.1, 0.2, 0.3, 0.4], dtype=float)


def make_single_series_fixture() -> list[np.ndarray]:
    """Makes a fixture that is one time history series.

    The values rise by one each time step, so a test can tell one time step's value from
    another's.

    :return: A list of one (5,) ndarray of floats, one value per time step.
    """
    return [np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=float)]


def make_three_series_fixture() -> list[np.ndarray]:
    """Makes a fixture that is three time history series.

    The three series rise, fall, and hold constant, so a test can tell one series from
    another and one time step's value from another's.

    :return: A list of three (5,) ndarrays of floats, one value per time step.
    """
    return [
        np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=float),
        np.array([0.0, -0.5, -1.0, -1.5, -2.0], dtype=float),
        np.array([1.0, 1.0, 1.0, 1.0, 1.0], dtype=float),
    ]


def make_three_labels_fixture() -> list[str]:
    """Makes a fixture that is the legend labels for three series.

    :return: A list of three legend labels, one per series.
    """
    return ["Induced Drag", "Side Force", "Lift"]


def make_three_colors_fixture() -> list[str]:
    """Makes a fixture that is the line colors for three series.

    :return: A list of three line colors, one per series.
    """
    return ["tab:blue", "tab:orange", "tab:green"]


def make_three_headers_fixture() -> list[str]:
    """Makes a fixture that is the CSV column headers for three columns.

    :return: A list of three column headers, one per column, not counting time.
    """
    return [
        "Induced Drag in Wind Axes (N)",
        "Side Force in Wind Axes (N)",
        "Lift in Wind Axes (N)",
    ]


def make_figure_size_fixture() -> tuple[float, float]:
    """Makes a fixture that is a figure's width and height.

    :return: A tuple of two floats representing the figure's width and height in inches.
    """
    return 6.0, 4.0
