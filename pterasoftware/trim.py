"""This module contains functions to analyze the trim conditions of steady and
unsteady solvers.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    analyze_steady_trim: This function attempts to calculate a trim condition of a
    steady solver by varying the operating point's velocity, angle of attack,
    angle of sideslip, and external thrust until the net force and net moment on the
    aircraft are sufficient low. If a trim condition can be found, it returns the
    trimmed operating point values. Otherwise, it logs an error.

    analyze_unsteady_trim: This function attempts to calculate a trim condition of an
    unsteady solver by varying the operating point's velocity, angle of attack,
    angle of sideslip, and external thrust until the net cycle-averaged force and net
    cycle-averaged moment on the aircraft are sufficient low. If a trim condition can
    be found, it returns the trimmed operating point values. Otherwise, it logs an
    error. """
import logging

import numpy as np
import scipy.optimize

from . import steady_horseshoe_vortex_lattice_method
from . import unsteady_ring_vortex_lattice_method
from . import movement
from . import problems


trim_logger = logging.getLogger("trim")
trim_logger.setLevel(logging.DEBUG)


# ToDo: Document this function.
def analyze_steady_trim(
    problem,
    velocity_bounds,
    alpha_bounds,
    beta_bounds,
    external_thrust_bounds,
    objective_cut_off=0.01,
    num_calls=50,
):
    """This function attempts to calculate a trim condition of a steady solver by
    varying the operating point's velocity, angle of attack, angle of sideslip,
    and external thrust until the net force and net moment on the aircraft are
    sufficient low. If a trim condition can be found, it returns the trimmed
    operating point values. Otherwise, it logs an error.

    :return:
    """
    if len(problem.airplanes) != 1:
        trim_logger.error(
            "The problem objects for trim analyses must have only one airplane."
        )

    weight = problem.airplanes[0].weight
    base_velocity = problem.operating_point.velocity
    base_alpha = problem.operating_point.alpha
    base_beta = problem.operating_point.beta
    base_external_thrust = problem.operating_point.external_thrust

    if base_velocity < velocity_bounds[0] or base_velocity > velocity_bounds[1]:
        trim_logger.error(
            "The operating point's velocity must be within the specified velocity "
            "bounds."
        )
    if base_alpha < alpha_bounds[0] or base_alpha > alpha_bounds[1]:
        trim_logger.error(
            "The operating point's alpha must be within the specified alpha bounds."
        )
    if base_beta < beta_bounds[0] or base_beta > beta_bounds[1]:
        trim_logger.error(
            "The operating point's beta must be within the specified beta bounds."
        )
    if (
        base_external_thrust < external_thrust_bounds[0]
        or base_external_thrust > external_thrust_bounds[1]
    ):
        trim_logger.error(
            "The operating point's external thrust must be within the specified "
            "external thrust bounds."
        )

    # ToDo: Document this function.
    def objective_function(arguments):
        """

        :param arguments:
        :return:
        """
        velocity, alpha, beta, external_thrust = arguments

        problem.operating_point.velocity = velocity
        problem.operating_point.alpha = alpha
        problem.operating_point.beta = beta
        dynamic_pressure = problem.operating_point.calculate_dynamic_pressure()
        s_ref = problem.airplanes[0].s_ref

        external_forces = np.array([external_thrust, 0, weight])
        external_force_coefficients = external_forces / dynamic_pressure / s_ref

        solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
            steady_problem=problem
        )

        solver.run()

        airplane = solver.airplanes[0]

        net_force_coefficient = np.linalg.norm(
            airplane.total_near_field_force_coefficients_wind_axes
            - external_force_coefficients
        )
        net_moment_coefficient = np.linalg.norm(
            airplane.total_near_field_moment_coefficients_wind_axes
        )

        objective = (abs(net_force_coefficient) + abs(net_moment_coefficient)) / 2

        return objective

    initial_guess = np.array(
        [base_velocity, base_alpha, base_beta, base_external_thrust]
    )
    bounds = (velocity_bounds, alpha_bounds, beta_bounds, external_thrust_bounds)

    trim_logger.info("Starting local optimization.")
    result_local = scipy.optimize.minimize(
        fun=objective_function,
        x0=initial_guess,
        bounds=bounds,
        options={"maxiter": num_calls},
    )
    trim_logger.info("Local optimization function executed.")

    if result_local.fun < objective_cut_off:
        trim_logger.info("Acceptable local minima found.")
        return result_local.x

    trim_logger.warning("No acceptable local minima found. Starting global search.")
    result_global = scipy.optimize.dual_annealing(
        func=objective_function,
        bounds=bounds,
        x0=initial_guess,
        maxfun=num_calls,
    )

    if result_global.fun < objective_cut_off:
        trim_logger.info("Acceptable global minima found.")
        return result_global.x

    trim_logger.error(
        "No trim condition found. Try increasing the bounds and the maximum number of "
        "iterations."
    )


# ToDo: Document this function.
def analyze_unsteady_trim(
    airplane_movement,
    operating_point,
    velocity_bounds,
    alpha_bounds,
    beta_bounds,
    objective_cut_off=0.01,
    num_calls=50,
):
    """analyze_unsteady_trim: This function attempts to calculate a trim condition of
    an unsteady solver by varying the operating point's velocity, angle of attack,
    angle of sideslip, and external thrust until the net cycle-averaged force and net
    cycle-averaged moment on the aircraft are sufficient low. If a trim condition can
    be found, it returns the trimmed operating point values. Otherwise, it logs an
    error.

    :return:
    """
    weight = airplane_movement.base_airplane.weight
    base_velocity = operating_point.velocity
    base_alpha = operating_point.alpha
    base_beta = operating_point.beta

    if base_velocity < velocity_bounds[0] or base_velocity > velocity_bounds[1]:
        trim_logger.error(
            "The operating point's velocity must be within the specified velocity "
            "bounds."
        )
    if base_alpha < alpha_bounds[0] or base_alpha > alpha_bounds[1]:
        trim_logger.error(
            "The operating point's alpha must be within the specified alpha bounds."
        )
    if base_beta < beta_bounds[0] or base_beta > beta_bounds[1]:
        trim_logger.error(
            "The operating point's beta must be within the specified beta bounds."
        )

    # ToDo: Document this function.
    def objective_function(arguments):
        """

        :param arguments:
        :return:
        """
        velocity, alpha, beta = arguments

        operating_point.velocity = velocity
        operating_point.alpha = alpha
        operating_point.beta = beta
        external_forces = np.array([0, 0, weight])
        dynamic_pressure = operating_point.calculate_dynamic_pressure()
        s_ref = airplane_movement.base_airplane.s_ref
        external_force_coefficients = external_forces / dynamic_pressure / s_ref

        operating_point_movement = movement.OperatingPointMovement(
            base_operating_point=operating_point
        )

        this_movement = movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
        )

        this_problem = problems.UnsteadyProblem(
            movement=this_movement, only_final_results=True
        )

        this_solver = (
            unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                unsteady_problem=this_problem
            )
        )

        this_solver.run(logging_level="Critical")

        force_coefficients = (
            this_solver.unsteady_problem.final_total_near_field_force_coefficients_wind_axes
        )
        moment_coefficients = (
            this_solver.unsteady_problem.final_total_near_field_moment_coefficients_wind_axes
        )

        net_force_coefficients = np.linalg.norm(
            force_coefficients - external_force_coefficients
        )
        net_moment_coefficients = np.linalg.norm(moment_coefficients)

        objective = (abs(net_force_coefficients) + abs(net_moment_coefficients)) / 2

        v_str = str(round(velocity, 3))
        a_str = str(round(alpha, 3))
        b_str = str(round(beta, 3))
        o_str = str(round(objective, 3))

        state_msg = "State: velocity=" + v_str + ", alpha=" + a_str + ", beta=" + b_str
        obj_msg = "Objective: " + o_str

        trim_logger.info(state_msg)
        trim_logger.info(obj_msg)
        return objective

    initial_guess = np.array([base_velocity, base_alpha, base_beta])
    bounds = (velocity_bounds, alpha_bounds, beta_bounds)

    trim_logger.info("Starting local optimization.")
    result_local = scipy.optimize.minimize(
        fun=objective_function,
        x0=initial_guess,
        bounds=bounds,
        options={"maxiter": num_calls},
    )
    trim_logger.info("Local optimization function executed.")

    if result_local.fun < objective_cut_off:
        trim_logger.info("Acceptable local minima found.")
        return result_local.x

    trim_logger.info("No acceptable local minima found. Starting global search.")
    result_global = scipy.optimize.dual_annealing(
        func=objective_function,
        bounds=bounds,
        x0=initial_guess,
        maxfun=num_calls,
    )

    if result_global.fun < objective_cut_off:
        trim_logger.info("Acceptable global minima found.")
        return result_global.x

    trim_logger.critical(
        "No trim condition found. Try increasing the bounds and the maximum number of "
        "iterations."
    )
