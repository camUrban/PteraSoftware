"""This module contains functions to analyze the trim conditions of steady and
unsteady solvers.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    analyze_steady_trim: This function attempts to calculate a trim condition of a
    steady solver by varying the angles of attack and sideslip until the steady
    pitching and yawing moments are zero. If a trim condition can be found,
    it returns the angles of attack and sideslip. Otherwise, it returns NaN angles
    and logs the failure.

    analyze_unsteady_trim: This function attempts to calculate a cycle-averaged trim
    condition of an unsteady solver by varying the angles of attack and sideslip
    until the unsteady, cycle-averaged pitching and yawing moments are zero. If a
    trim condition can be found, it returns the angles of attack and sideslip.
    Otherwise, it returns NaN angles and logs the failure. """
import logging

import numpy as np

from . import steady_horseshoe_vortex_lattice_method


# ToDo: Document this function.
def analyze_steady_trim(
    problem,
    velocity_offsets=(-5, 5),
    alpha_offsets=(-5, 5),
    beta_offsets=(-5, 5),
    num_samples=5,
    weight=100,
    thrust=25,
):
    """This function attempts to calculate a trim condition of a steady solver by
    varying the angles of attack and sideslip until the steady pitching and yawing
    moments are zero. If a trim condition can be found, it returns the angles of
    attack and sideslip. Otherwise, it returns NaN angles and logs the failure.

    :return:
    """
    if len(problem.airplanes) != 1:
        logging.critical(
            "The problem objects for trim analyses must have only one airplane."
        )

    base_velocity = problem.operating_point.velocity
    base_alpha = problem.operating_point.alpha
    base_beta = problem.operating_point.beta

    neg_velocity_offset, pos_velocity_offset = velocity_offsets
    neg_alpha_offset, pos_alpha_offset = alpha_offsets
    neg_beta_offset, pos_beta_offset = beta_offsets

    min_velocity = max(0, base_velocity + neg_velocity_offset)
    max_velocity = base_velocity + pos_velocity_offset
    min_alpha = max(-90, base_alpha + neg_alpha_offset)
    max_alpha = min(90, base_alpha + pos_alpha_offset)
    min_beta = max(-90, base_beta + neg_beta_offset)
    max_beta = min(90, base_beta + pos_beta_offset)

    velocities = np.linspace(min_velocity, max_velocity, num_samples)
    alphas = np.linspace(min_alpha, max_alpha, num_samples)
    betas = np.linspace(min_beta, max_beta, num_samples)

    net_forces = np.zeros((num_samples, num_samples, num_samples))
    net_moments = np.zeros((num_samples, num_samples, num_samples))

    ext_force = np.array([thrust, 0, weight])

    for i, velocity in enumerate(velocities):
        for j, alpha in enumerate(alphas):
            for k, beta in enumerate(betas):

                problem.operating_point.velocity = velocity
                problem.operating_point.alpha = alpha
                problem.operating_point.beta = beta

                solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                    steady_problem=problem
                )

                solver.run()

                airplane = solver.airplanes[0]

                net_force = np.linalg.norm(
                    airplane.total_near_field_force_wind_axes + ext_force
                )
                net_moment = np.linalg.norm(airplane.total_near_field_moment_wind_axes)

                net_forces[i, j, k] = net_force
                net_moments[i, j, k] = net_moment


# ToDo: Document this function.
def analyze_unsteady_trim():
    """This function attempts to calculate a cycle-averaged trim condition of an
    unsteady solver by varying the angles of attack and sideslip until the unsteady,
    cycle-averaged pitching and yawing moments are zero. If a trim condition can be
    found, it returns the angles of attack and sideslip. Otherwise, it returns NaN
    angles and logs the failure.

    :return:
    """
    pass
