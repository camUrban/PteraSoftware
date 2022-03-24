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
import scipy.optimize

from . import steady_horseshoe_vortex_lattice_method


# ToDo: Document this function.
def analyze_steady_trim(
    problem,
    weight=250,
    base_thrust=10,
    num_iter=100,
):
    """

    :return:
    """
    if len(problem.airplanes) != 1:
        logging.critical(
            "The problem objects for trim analyses must have only one airplane."
        )

    obj_cut_off = ((weight + base_thrust) / 2) / 100

    base_velocity = problem.operating_point.velocity
    base_alpha = problem.operating_point.alpha
    base_beta = problem.operating_point.beta

    # ToDo: Document this function.
    def obj(arguments):
        """

        :param arguments:
        :return:
        """
        velocity, alpha, beta, thrust = arguments

        problem.operating_point.velocity = velocity
        problem.operating_point.alpha = alpha
        problem.operating_point.beta = beta
        ext_force = np.array([thrust, 0, weight])

        solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
            steady_problem=problem
        )

        solver.run()

        airplane = solver.airplanes[0]

        net_force = np.linalg.norm(
            airplane.total_near_field_force_wind_axes + ext_force
        )
        net_moment = np.linalg.norm(airplane.total_near_field_moment_wind_axes)
        return abs(net_force) + abs(net_moment)

    logging.info("Starting local optimization.")
    result_local = scipy.optimize.minimize(
        fun=obj,
        x0=(base_velocity, base_alpha, base_beta, base_thrust),
    )
    logging.info("Local optimization function executed.")

    if result_local.fun < obj_cut_off:
        logging.info("Acceptable local minima found.")
        return result_local.x

    logging.warning("No acceptable local minima found. Starting global search.")
    result_global = scipy.optimize.dual_annealing(
        func=obj,
        bounds=[(5, 15), (-10, 10), (-10, 10), (0, 10)],
        x0=(base_velocity, base_alpha, base_beta, base_thrust),
        maxfun=num_iter,
    )

    if result_global.fun < obj_cut_off:
        logging.info("Acceptable global minima found.")
        return result_global.x

    logging.error("No trim condition found.")


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
