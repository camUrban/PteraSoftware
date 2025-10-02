"""This module contains functions to create fixtures for aerodynamics function tests."""

import numpy as np


def make_single_point_fixture():
    """This method makes a fixture that is a single evaluation point for testing
    velocity calculation functions.

    :return single_point_fixture: (1, 3) ndarray of floats
        This is a single evaluation point at the origin (in geometry axes,
        relative to the CG).
    """
    single_point_fixture = np.array([[0.0, 0.0, 0.0]], dtype=float)

    return single_point_fixture


def make_grid_of_points_fixture():
    """This method makes a fixture that is a grid of evaluation points for testing
    velocity calculation functions.

    :return grid_of_points_fixture: (25, 3) ndarray of floats
        This is a 5x5 grid of evaluation points in the xz-plane (in geometry axes,
        relative to the CG), spanning from -2 to 2 in both x and z directions.
    """
    x = np.linspace(-2.0, 2.0, 5, dtype=float)
    z = np.linspace(-2.0, 2.0, 5, dtype=float)
    xx, zz = np.meshgrid(x, z)
    grid_of_points_fixture = np.column_stack(
        [xx.ravel(), np.zeros(25, dtype=float), zz.ravel()]
    )

    return grid_of_points_fixture


def make_line_of_points_fixture():
    """This method makes a fixture that is a line of evaluation points for testing
    velocity calculation functions.

    :return line_of_points_fixture: (10, 3) ndarray of floats
        This is a line of 10 evaluation points along the x-axis (in geometry axes,
        relative to the CG), from x=-5.0 to x=5.0.
    """
    x = np.linspace(-5.0, 5.0, 10, dtype=float)
    line_of_points_fixture = np.column_stack(
        [x, np.zeros(10, dtype=float), np.zeros(10, dtype=float)]
    )

    return line_of_points_fixture


def make_random_points_fixture():
    """This method makes a fixture that is a set of random evaluation points for
    testing velocity calculation functions.

    :return random_points_fixture: (20, 3) ndarray of floats
        This is a set of 20 random evaluation points (in geometry axes, relative
        to the CG) within a cube of side length 10 centered at the origin.
    """
    np.random.seed(42)
    random_points_fixture = np.random.uniform(-5.0, 5.0, (20, 3)).astype(float)

    return random_points_fixture


def make_simple_ring_vortex_arrays_fixture():
    """This method makes a fixture containing arrays describing a single RingVortex
    for testing velocity calculation functions.

    :return tuple of ndarrays
        This returns a tuple containing:
        - stackBrrvp_G_Cg: (1, 3) ndarray of floats
        - stackFrrvp_G_Cg: (1, 3) ndarray of floats
        - stackFlrvp_G_Cg: (1, 3) ndarray of floats
        - stackBlrvp_G_Cg: (1, 3) ndarray of floats
        - strengths: (1,) ndarray of floats
    """
    stackBrrvp_G_Cg = np.array([[1.0, 0.5, 0.0]], dtype=float)
    stackFrrvp_G_Cg = np.array([[0.0, 0.5, 0.0]], dtype=float)
    stackFlrvp_G_Cg = np.array([[0.0, -0.5, 0.0]], dtype=float)
    stackBlrvp_G_Cg = np.array([[1.0, -0.5, 0.0]], dtype=float)
    strengths = np.array([1.0], dtype=float)

    return stackBrrvp_G_Cg, stackFrrvp_G_Cg, stackFlrvp_G_Cg, stackBlrvp_G_Cg, strengths


def make_multiple_ring_vortex_arrays_fixture():
    """This method makes a fixture containing arrays describing multiple RingVortices
    for testing velocity calculation functions.

    :return tuple of ndarrays
        This returns a tuple containing arrays for 3 RingVortices:
        - stackBrrvp_G_Cg: (3, 3) ndarray of floats
        - stackFrrvp_G_Cg: (3, 3) ndarray of floats
        - stackFlrvp_G_Cg: (3, 3) ndarray of floats
        - stackBlrvp_G_Cg: (3, 3) ndarray of floats
        - strengths: (3,) ndarray of floats
    """
    stackBrrvp_G_Cg = np.array(
        [[1.0, 0.5, 0.0], [2.0, 0.5, 0.0], [3.0, 0.5, 0.0]], dtype=float
    )
    stackFrrvp_G_Cg = np.array(
        [[0.0, 0.5, 0.0], [1.0, 0.5, 0.0], [2.0, 0.5, 0.0]], dtype=float
    )
    stackFlrvp_G_Cg = np.array(
        [[0.0, -0.5, 0.0], [1.0, -0.5, 0.0], [2.0, -0.5, 0.0]], dtype=float
    )
    stackBlrvp_G_Cg = np.array(
        [[1.0, -0.5, 0.0], [2.0, -0.5, 0.0], [3.0, -0.5, 0.0]], dtype=float
    )
    strengths = np.array([1.0, 1.5, 2.0], dtype=float)

    return stackBrrvp_G_Cg, stackFrrvp_G_Cg, stackFlrvp_G_Cg, stackBlrvp_G_Cg, strengths


def make_simple_horseshoe_vortex_arrays_fixture():
    """This method makes a fixture containing arrays describing a single
    HorseshoeVortex for testing velocity calculation functions.

    :return tuple of ndarrays
        This returns a tuple containing:
        - stackBrhvp_G_Cg: (1, 3) ndarray of floats
        - stackFrhvp_G_Cg: (1, 3) ndarray of floats
        - stackFlhvp_G_Cg: (1, 3) ndarray of floats
        - stackBlhvp_G_Cg: (1, 3) ndarray of floats
        - strengths: (1,) ndarray of floats
    """
    stackBrhvp_G_Cg = np.array([[20.0, 0.5, 0.0]], dtype=float)
    stackFrhvp_G_Cg = np.array([[0.0, 0.5, 0.0]], dtype=float)
    stackFlhvp_G_Cg = np.array([[0.0, -0.5, 0.0]], dtype=float)
    stackBlhvp_G_Cg = np.array([[20.0, -0.5, 0.0]], dtype=float)
    strengths = np.array([1.0], dtype=float)

    return stackBrhvp_G_Cg, stackFrhvp_G_Cg, stackFlhvp_G_Cg, stackBlhvp_G_Cg, strengths


def make_multiple_horseshoe_vortex_arrays_fixture():
    """This method makes a fixture containing arrays describing multiple
    HorseshoeVortices for testing velocity calculation functions.

    :return tuple of ndarrays
        This returns a tuple containing arrays for 2 HorseshoeVortices:
        - stackBrhvp_G_Cg: (2, 3) ndarray of floats
        - stackFrhvp_G_Cg: (2, 3) ndarray of floats
        - stackFlhvp_G_Cg: (2, 3) ndarray of floats
        - stackBlhvp_G_Cg: (2, 3) ndarray of floats
        - strengths: (2,) ndarray of floats
    """
    stackBrhvp_G_Cg = np.array([[20.0, 0.5, 0.0], [20.0, 1.5, 0.0]], dtype=float)
    stackFrhvp_G_Cg = np.array([[0.0, 0.5, 0.0], [0.0, 1.5, 0.0]], dtype=float)
    stackFlhvp_G_Cg = np.array([[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]], dtype=float)
    stackBlhvp_G_Cg = np.array([[20.0, -0.5, 0.0], [20.0, 0.5, 0.0]], dtype=float)
    strengths = np.array([1.0, -0.5], dtype=float)

    return stackBrhvp_G_Cg, stackFrhvp_G_Cg, stackFlhvp_G_Cg, stackBlhvp_G_Cg, strengths


def make_ages_fixture():
    """This method makes a fixture that is an ndarray of ages for wake vortices.

    :return ages_fixture: (3,) ndarray of floats
        This is an ndarray of ages in seconds for 3 vortices.
    """
    ages_fixture = np.array([0.1, 0.5, 1.0], dtype=float)

    return ages_fixture


def make_zero_ages_fixture():
    """This method makes a fixture that is an ndarray of zero ages for bound vortices.

    :return zero_ages_fixture: (3,) ndarray of floats
        This is an ndarray of zero ages for 3 vortices (simulating bound vortices
        with no core radius).
    """
    zero_ages_fixture = np.array([0.0, 0.0, 0.0], dtype=float)

    return zero_ages_fixture


def make_kinematic_viscosity_fixture():
    """This method makes a fixture that is a kinematic viscosity value for air.

    :return kinematic_viscosity_fixture: float
        This is the kinematic viscosity of air at standard conditions, approximately
        1.5e-5 meters squared per second.
    """
    kinematic_viscosity_fixture = 1.5e-5

    return kinematic_viscosity_fixture
