"""Contains functions to analyze the trim conditions of SteadyProblems and
UnsteadyProblems.

This module contains the following classes:
    None

This module contains the following functions:
    analyze_steady_trim: This function attempts to calculate a trim condition of a
    SteadyProblem by varying the OperatingPoint's vCg__E, alpha, beta,
    and externalFX_W parameters until the net loads are sufficient low. If a trim
    condition can be found, it returns the trimmed OperatingPoint parameters.

    analyze_unsteady_trim: This function attempts to calculate a trim condition of an
    UnsteadyProblem by varying the OperatingPointMovement's base OperatingPoint's
    vCg__E, alpha, beta, and externalFX_W parameters until the net final loads are
    sufficient low. If a trim condition can be found, it returns the trimmed base
    OperatingPoint parameters.
"""

import logging
from typing import Sequence, Any

import numpy as np
import scipy.optimize

from . import _parameter_validation
from . import movements
from . import problems
from . import steady_horseshoe_vortex_lattice_method
from . import steady_ring_vortex_lattice_method
from . import unsteady_ring_vortex_lattice_method

trim_logger = logging.getLogger("trim")
trim_logger.setLevel(logging.DEBUG)
logging.basicConfig()

# Set a seed for reproducibility in the dual annealing optimizer.
_seed = 42


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
# TODO: It would be awesome if we could incorporate control surface deflections into
#  the steady trim analysis.
def analyze_steady_trim(
    problem: problems.SteadyProblem,
    solver_type: str,
    boundsVCg__E: tuple[float, float],
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    boundsExternalFX_W: tuple[float, float],
    objective_cut_off: float = 0.01,
    num_calls: int = 100,
) -> tuple[float, float, float, float] | tuple[None, None, None, None]:
    """This function attempts to calculate a trim condition of a SteadyProblem by
    varying the OperatingPoint's vCg__E, alpha, beta, and externalFX_W parameters
    until the net loads are sufficient low. If a trim condition can be found,
    it returns the trimmed OperatingPoint parameters.

    Trim is found by minimizing an objective function that combines the net load
    coefficients' magnitudes. The function uses a two-stage optimization approach:
    first attempting a local search using L-BFGS-B, and if that fails to converge,
    it performs a global search using dual annealing.

    The optimization process varies four parameters (vCg__E, alpha, beta,
    and externalFX_W) within their specified bounds. At each iteration, a steady
    solver is run to compute aerodynamic loads. These loads are combined with the
    external loads to compute a net load coefficient.

    The search terminates early if the objective function falls below the specified
    cutoff value. If no trim condition is found within the maximum number of
    function calls, the function returns None values and logs a critical error.

    :param problem: SteadyProblem

        This is the SteadyProblem whose trim condition will be found. The
        SteadyProblem must contain exactly one Airplane. The problem's OperatingPoint
        will be modified during the trim search.

    :param solver_type: str

        This parameter determines what type of steady solver will be used to analyze
        the SteadyProblem. The options are "steady horseshoe vortex lattice method"
        and "steady ring vortex lattice method".

    :param boundsVCg__E: tuple of two positive numbers

        This parameter is a tuple of two positive numbers (ints or floats) in
        ascending order. It determines the range of speeds of the Airplane's CG (in
        the Earth frame) to search, in units of meters per second. The
        SteadyProblem's OperatingPoint's initial vCg__E must be within these bounds.

    :param alpha_bounds: tuple of two numbers

        This parameter is a tuple of two numbers (ints or floats) in ascending order.
        It determines the range of angles of attack to search, in units of degrees.
        The SteadyProblem's OperatingPoint's initial alpha must be within these bounds.

    :param beta_bounds: tuple of two numbers

        This parameter is a tuple of two numbers (ints or floats) in ascending order.
        It determines the range of sideslip angles to search, in units of degrees.
        The SteadyProblem's OperatingPoint's initial beta must be within these bounds.

    :param boundsExternalFX_W: tuple of two numbers

        This parameter is a tuple of two numbers (ints or floats) in ascending order.
        It determines the range of external forces (in the wind axes' x-direction) to
        search, in units of Newtons. The SteadyProblem's OperatingPoint's initial
        externalFX_W must be within these bounds.

    :param objective_cut_off: positive number, optional

        This parameter must be a positive number (int or float), and it determines
        the convergence threshold for the trim search. When the objective function
        falls below this value, the search terminates successfully. Lower values
        result in tighter trim conditions but may require more iterations. The
        default value is 0.01.

    :param num_calls: positive int, optional

        This parameter must be a positive int, and it determines the maximum number
        of objective function evaluations allowed during optimization. This limit
        applies separately to the local search and global search stages. Higher
        values allow for more thorough searching but increase computation time. The
        default value is 100.

    :return: tuple of four floats or tuple of four Nones

        This function returns a tuple of four floats representing the trimmed
        parameters found. In order, they are: vCg__E (in meters per second),
        alpha (in degrees), beta (in degrees), and externalFX_W (in Newtons). If no
        trim condition was found, it will instead return a tuple of four Nones.
    """
    # Validate the problem parameter.
    if not isinstance(problem, problems.SteadyProblem):
        raise TypeError("problem must be a SteadyProblem.")
    if len(problem.airplanes) != 1:
        raise ValueError(
            "The SteadyProblem must contain exactly one Airplane for trim analysis."
        )

    # Validate the solver_type parameter.
    if solver_type not in (
        "steady horseshoe vortex lattice method",
        "steady ring vortex lattice method",
    ):
        raise ValueError(
            'solver_type must be either "steady horseshoe vortex lattice method" or '
            '"steady ring vortex lattice method".'
        )

    # Validate the boundsVCg__E parameter.
    if not (isinstance(boundsVCg__E, tuple) and len(boundsVCg__E) == 2):
        raise TypeError("boundsVCg__E must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in boundsVCg__E):
        raise TypeError("Both values in boundsVCg__E must be numbers.")
    if boundsVCg__E[0] > boundsVCg__E[1]:
        raise ValueError(
            "The first value in boundsVCg__E must be less than or equal to the second "
            "value."
        )
    if boundsVCg__E[0] <= 0:
        raise ValueError("Both values in boundsVCg__E must be positive.")

    # Validate the alpha_bounds parameter.
    if not (isinstance(alpha_bounds, tuple) and len(alpha_bounds) == 2):
        raise TypeError("alpha_bounds must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in alpha_bounds):
        raise TypeError("Both values in alpha_bounds must be numbers.")
    if alpha_bounds[0] > alpha_bounds[1]:
        raise ValueError(
            "The first value in alpha_bounds must be less than or equal to the second "
            "value."
        )

    # Validate the beta_bounds parameter.
    if not (isinstance(beta_bounds, tuple) and len(beta_bounds) == 2):
        raise TypeError("beta_bounds must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in beta_bounds):
        raise TypeError("Both values in beta_bounds must be numbers.")
    if beta_bounds[0] > beta_bounds[1]:
        raise ValueError(
            "The first value in beta_bounds must be less than or equal to the second "
            "value."
        )

    # Validate the boundsExternalFX_W parameter.
    if not (isinstance(boundsExternalFX_W, tuple) and len(boundsExternalFX_W) == 2):
        raise TypeError("boundsExternalFX_W must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in boundsExternalFX_W):
        raise TypeError("Both values in boundsExternalFX_W must be numbers.")
    if boundsExternalFX_W[0] > boundsExternalFX_W[1]:
        raise ValueError(
            "The first value in boundsExternalFX_W must be less than or equal to the "
            "second value."
        )

    # Validate the objective_cut_off parameter.
    objective_cut_off = _parameter_validation.positive_number_return_float(
        objective_cut_off, "objective_cut_off"
    )

    # Validate the num_calls parameter.
    num_calls = _parameter_validation.positive_int_return_int(num_calls, "num_calls")

    # Get the SteadyProblem's OperatingPoint's initial parameter values and check
    # that the ones that will vary to find a trim condition are within the specified
    # bounds.
    weight = problem.airplanes[0].weight
    baseVCg__E = problem.operating_point.vCg__E
    base_alpha = problem.operating_point.alpha
    base_beta = problem.operating_point.beta
    baseExternalFX_W = problem.operating_point.externalFX_W

    if baseVCg__E < boundsVCg__E[0] or baseVCg__E > boundsVCg__E[1]:
        raise ValueError(
            "The OperatingPoint's vCg__E must be within the specified vCg__E bounds."
        )
    if base_alpha < alpha_bounds[0] or base_alpha > alpha_bounds[1]:
        raise ValueError(
            "The OperatingPoint's alpha must be within the specified alpha bounds."
        )
    if base_beta < beta_bounds[0] or base_beta > beta_bounds[1]:
        raise ValueError(
            "The OperatingPoint's beta must be within the specified beta bounds."
        )
    if (
        baseExternalFX_W < boundsExternalFX_W[0]
        or baseExternalFX_W > boundsExternalFX_W[1]
    ):
        raise ValueError(
            "The OperatingPoint's externalFX_W must be within the specified "
            "externalFX_W bounds."
        )

    current_arguments = [np.nan, np.nan, np.nan, np.nan]

    def objective_function(arguments: tuple[float, float, float, float]) -> float:
        """Computes the trim objective function for a given set of OperatingPoint
        parameters.

        This function evaluates the trim quality by running a steady solver and
        computing the average magnitude of net load coefficients. It updates the
        outer scope's current_arguments list to communicate the last evaluated state
        back to the optimizer.

        If the objective falls below the cutoff threshold, this function raises
        StopIteration to signal successful trim convergence.

        :param arguments: array-like of four numbers

            The OperatingPoint parameters to evaluate: [vCg__E, alpha, beta,
            externalFX_W].

        :return: float

            The trim objective value, computed as the average of the magnitude of the
            net load coefficients.
        """
        vCg__E, alpha, beta, externalFX_W = arguments

        current_arguments.clear()
        current_arguments.extend([vCg__E, alpha, beta, externalFX_W])

        problem.operating_point.vCg__E = vCg__E
        problem.operating_point.alpha = alpha
        problem.operating_point.beta = beta

        qInf__E = problem.operating_point.qInf__E
        s_ref = problem.airplanes[0].s_ref

        # To my knowledge, there isn't a standard way to define "external" force
        # coefficients. However, simply checking trim against forces and moments in
        # Newtons and Newton meters isn't a good solution because, for example,
        # an imbalance of 0.01 N may be trivial to a bird-scale UAV but critical to a
        # flying microrobot. I'm choosing to non-dimensionalize with reference area,
        # as that is what is used for aerodynamic force coefficients. If we later
        # allow users to apply external moments we may need to come up with a better
        # approach, as moment coefficients non-dimensionalize using different
        # dimensions.
        externalForces_W = np.array([externalFX_W, 0, weight])
        externalForceCoefficients_W = externalForces_W / qInf__E / s_ref

        if solver_type == "steady horseshoe vortex lattice method":
            solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                steady_problem=problem
            )
        else:
            solver = (
                steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
                    steady_problem=problem
                )
            )

        solver.run()

        airplane = solver.airplanes[0]

        netForceCoefficient_W = np.linalg.norm(
            airplane.forceCoefficients_W + externalForceCoefficients_W
        )
        netMomentCoefficient_W_CgP1 = np.linalg.norm(airplane.momentCoefficients_W_CgP1)

        objective = (netForceCoefficient_W + netMomentCoefficient_W_CgP1) / 2

        v_str = str(round(vCg__E, 2))
        a_str = str(round(alpha, 2))
        b_str = str(round(beta, 2))
        f_str = str(round(externalFX_W, 2))
        o_str = str(round(objective, 3))

        state_msg = (
            "\tState: vCg__E="
            + v_str
            + ", alpha="
            + a_str
            + ", beta="
            + b_str
            + ", externalFX_W="
            + f_str
        )
        obj_msg = "\t\tObjective: " + o_str

        trim_logger.info(state_msg)
        trim_logger.info(obj_msg)

        if objective < objective_cut_off:
            raise StopIteration

        current_arguments.clear()
        current_arguments.extend([np.nan, np.nan, np.nan, np.nan])

        return objective

    initial_guess = np.array([baseVCg__E, base_alpha, base_beta, baseExternalFX_W])
    bounds: Sequence[tuple[float, float]] = [
        (boundsVCg__E[0], boundsVCg__E[1]),
        (alpha_bounds[0], alpha_bounds[1]),
        (beta_bounds[0], beta_bounds[1]),
        (boundsExternalFX_W[0], boundsExternalFX_W[1]),
    ]

    trim_logger.info("Starting local search.")
    try:
        options: Any = {"maxfun": num_calls, "eps": 0.01}
        scipy.optimize.minimize(
            fun=objective_function,
            x0=initial_guess,
            bounds=bounds,
            method="L-BFGS-B",
            options=options,
        )
    except StopIteration:
        trim_logger.info("Acceptable value reached with local search.")
        return (
            current_arguments[0],
            current_arguments[1],
            current_arguments[2],
            current_arguments[3],
        )

    trim_logger.warning(
        "No acceptable value reached with local search. Starting global search."
    )
    try:
        options: Any = {"maxfun": num_calls, "eps": 0.01}
        minimizer_kwargs: Any = {
            "method": "L-BFGS-B",
            "options": options,
        }
        scipy.optimize.dual_annealing(
            func=objective_function,
            bounds=bounds,
            x0=initial_guess,
            maxfun=num_calls,
            minimizer_kwargs=minimizer_kwargs,
            seed=_seed,
        )
    except StopIteration:
        trim_logger.info("Acceptable global minima found.")
        return (
            current_arguments[0],
            current_arguments[1],
            current_arguments[2],
            current_arguments[3],
        )

    trim_logger.critical(
        "No trim condition found. Try increasing the bounds and the maximum number of "
        "iterations."
    )
    return None, None, None, None


# TEST: Consider adding unit tests for this function.
# TEST: Consider adding integration tests for this function.
# TODO: Add the ability to specify running the solver with a prescribed or free wake.
# TODO: It would be awesome if we could incorporate particular amplitude and period
#  parameters of the AirplaneMovement's WingMovements and WingCrossSectionMovements
#  into the unsteady trim analysis. Incorporating control surface deflection would
#  also be great but less important.
def analyze_unsteady_trim(
    problem: problems.UnsteadyProblem,
    boundsVCg__E: tuple[float, float],
    alpha_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    boundsExternalFX_W: tuple[float, float],
    objective_cut_off: float = 0.01,
    num_calls: int = 100,
) -> tuple[float, float, float, float] | tuple[None, None, None, None]:
    """This function attempts to calculate a trim condition of an UnsteadyProblem by
    varying the OperatingPointMovement's base OperatingPoint's vCg__E, alpha, beta,
    and externalFX_W parameters until the net final loads are sufficient low. If a
    trim condition can be found, it returns the trimmed base OperatingPoint parameters.

    Trim is found by minimizing an objective function that combines the net final
    load coefficients' magnitudes. For variable geometry cases, the final load
    coefficients are the RMS values from the final motion cycle. For static geometry
    cases, they are the load coefficients at the final time step. The function uses a
    two-stage optimization approach: first attempting a local search using L-BFGS-B,
    and if that fails to converge, it performs a global search using dual annealing.

    The optimization process varies four of the base OperatingPoint's parameters (
    vCg__E, alpha, beta, externalFX_W) within their specified bounds. At each
    iteration, an UnsteadyRingVortexLatticeMethodSolver is run to compute aerodynamic
    loads. These loads are combined with the external loads to compute a net load
    coefficient.

    The search terminates early if the objective function falls below the specified
    cutoff value. If no trim condition is found within the maximum number of function
    calls, the function returns None values and logs a critical error.

    :param problem: UnsteadyProblem

        This is the UnsteadyProblem whose trim condition will be found. The
        UnsteadyProblem's Movement must contain exactly one AirplaneMovement. The
        problem's OperatingPointMovement's base OperatingPoint will be modified
        during the trim search.

    :param boundsVCg__E: tuple of two positive numbers

        This parameter is a tuple of two positive numbers (ints or floats) in
        ascending order. It determines the range of base speeds of the Airplane's CG
        (in the Earth frame) to search, in units of meters per second. The base
        OperatingPoint's initial vCg__E must be within these bounds.

    :param alpha_bounds: tuple of two numbers

        This parameter is a tuple of two numbers (ints or floats) in ascending order.
        It determines the range of angles of attack to search, in units of degrees.
        The base OperatingPoint's initial alpha must be within these bounds.

    :param beta_bounds: tuple of two numbers

        This parameter is a tuple of two numbers (ints or floats) in ascending order.
        It determines the range of sideslip angles to search, in units of degrees.
        The base OperatingPoint's initial beta must be within these bounds.

    :param boundsExternalFX_W: tuple of two numbers

        This parameter is a tuple of two numbers (ints or floats) in ascending order.
        It determines the range of external forces (in the wind axes' x-direction) to
        search, in units of Newtons. The base OperatingPoint's initial
        externalFX_W must be within these bounds.

    :param objective_cut_off: positive number, optional

        This parameter must be a positive number (int or float), and it determines
        the convergence threshold for the trim search. When the objective function
        falls below this value, the search terminates successfully. Lower values
        result in tighter trim conditions but may require more iterations. The
        default value is 0.01.

    :param num_calls: positive int, optional

        This parameter must be a positive int, and it determines the maximum number
        of objective function evaluations allowed during optimization. This limit
        applies separately to the local search and global search stages. Higher
        values allow for more thorough searching but increase computation time. The
        default value is 100.

    :return: tuple of four floats or tuple of four Nones

        This function returns a tuple of four floats representing the trimmed
        parameters found. In order, they are: vCg__E (in meters per second),
        alpha (in degrees), beta (in degrees), and externalFX_W (in Newtons). If no
        trim condition was found, it will instead return a tuple of four Nones.
    """
    # Validate the problem parameter.
    if not isinstance(problem, problems.UnsteadyProblem):
        raise TypeError("problem must be an UnsteadyProblem.")
    if len(problem.movement.airplane_movements) != 1:
        raise ValueError(
            "The UnsteadyProblem's Movement must contain exactly one AirplaneMovement "
            "for trim analysis."
        )

    # Validate the boundsVCg__E parameter.
    if not (isinstance(boundsVCg__E, tuple) and len(boundsVCg__E) == 2):
        raise TypeError("boundsVCg__E must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in boundsVCg__E):
        raise TypeError("Both values in boundsVCg__E must be numbers.")
    if boundsVCg__E[0] > boundsVCg__E[1]:
        raise ValueError(
            "The first value in boundsVCg__E must be less than or equal to the second "
            "value."
        )
    if boundsVCg__E[0] <= 0:
        raise ValueError("Both values in boundsVCg__E must be positive.")

    # Validate the alpha_bounds parameter.
    if not (isinstance(alpha_bounds, tuple) and len(alpha_bounds) == 2):
        raise TypeError("alpha_bounds must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in alpha_bounds):
        raise TypeError("Both values in alpha_bounds must be numbers.")
    if alpha_bounds[0] > alpha_bounds[1]:
        raise ValueError(
            "The first value in alpha_bounds must be less than or equal to the second "
            "value."
        )

    # Validate the beta_bounds parameter.
    if not (isinstance(beta_bounds, tuple) and len(beta_bounds) == 2):
        raise TypeError("beta_bounds must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in beta_bounds):
        raise TypeError("Both values in beta_bounds must be numbers.")
    if beta_bounds[0] > beta_bounds[1]:
        raise ValueError(
            "The first value in beta_bounds must be less than or equal to the second "
            "value."
        )

    # Validate the boundsExternalFX_W parameter.
    if not (isinstance(boundsExternalFX_W, tuple) and len(boundsExternalFX_W) == 2):
        raise TypeError("boundsExternalFX_W must be a tuple with length 2.")
    if not all(isinstance(bound, (int, float)) for bound in boundsExternalFX_W):
        raise TypeError("Both values in boundsExternalFX_W must be numbers.")
    if boundsExternalFX_W[0] > boundsExternalFX_W[1]:
        raise ValueError(
            "The first value in boundsExternalFX_W must be less than or equal to the "
            "second value."
        )

    # Validate the objective_cut_off parameter.
    objective_cut_off = _parameter_validation.positive_number_return_float(
        objective_cut_off, "objective_cut_off"
    )

    # Validate the num_calls parameter.
    num_calls = _parameter_validation.positive_int_return_int(num_calls, "num_calls")

    base_operating_point = (
        problem.movement.operating_point_movement.base_operating_point
    )

    # Get the base OperatingPoint's initial parameter values and check
    # that the ones that will vary to find a trim condition are within the specified
    # bounds.
    weight = problem.movement.airplane_movements[0].base_airplane.weight
    baseVCg__E = base_operating_point.vCg__E
    base_alpha = base_operating_point.alpha
    base_beta = base_operating_point.beta
    baseExternalFX_W = base_operating_point.externalFX_W

    if baseVCg__E < boundsVCg__E[0] or baseVCg__E > boundsVCg__E[1]:
        raise ValueError(
            "The base OperatingPoint's vCg__E must be within the specified vCg__E "
            "bounds."
        )
    if base_alpha < alpha_bounds[0] or base_alpha > alpha_bounds[1]:
        raise ValueError(
            "The base OperatingPoint's alpha must be within the specified alpha bounds."
        )
    if base_beta < beta_bounds[0] or base_beta > beta_bounds[1]:
        raise ValueError(
            "The base OperatingPoint's beta must be within the specified beta bounds."
        )
    if (
        baseExternalFX_W < boundsExternalFX_W[0]
        or baseExternalFX_W > boundsExternalFX_W[1]
    ):
        raise ValueError(
            "The base OperatingPoint's externalFX_W must be within the specified "
            "externalFX_W bounds."
        )

    current_arguments = [np.nan, np.nan, np.nan, np.nan]

    def objective_function(arguments: tuple[float, float, float, float]) -> float:
        """Computes the trim objective function for a given set of OperatingPoint
        parameters.

        This function evaluates the trim quality by running an
        UnsteadyRingVortexLatticeMethodSolver and computing the average magnitude of
        net load coefficients. It updates the outer scope's current_arguments list to
        communicate the last evaluated state back to the optimizer.

        If the objective falls below the cutoff threshold, this function raises
        StopIteration to signal successful trim convergence.

        :param arguments: array-like of four numbers

            The base OperatingPoint parameters to evaluate: [vCg__E, alpha, beta,
            externalFX_W].

        :return: float

            The trim objective value, computed as the average of the magnitude of the
            net load coefficients.
        """
        vCg__E, alpha, beta, externalFX_W = arguments

        current_arguments.clear()
        current_arguments.extend([vCg__E, alpha, beta, externalFX_W])

        base_operating_point.vCg__E = vCg__E
        base_operating_point.alpha = alpha
        base_operating_point.beta = beta

        qInf__E = base_operating_point.qInf__E
        s_ref = problem.movement.airplane_movements[0].base_airplane.s_ref

        # To my knowledge, there isn't a standard way to define "external" force
        # coefficients. However, simply checking trim against forces and moments in
        # Newtons and Newton meters isn't a good solution because, for example,
        # an imbalance of 0.01 N may be trivial to a bird-scale UAV but critical to a
        # flying microrobot. I'm choosing to non-dimensionalize with reference area,
        # as that is what is used for aerodynamic force coefficients. If we later
        # allow users to apply external moments we may need to come up with a better
        # approach, as moment coefficients non-dimensionalize using different
        # dimensions.
        externalForces_W = np.array([externalFX_W, 0, weight])
        externalForceCoefficients_W = externalForces_W / qInf__E / s_ref

        this_operating_point_movement = (
            movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point
            )
        )

        this_movement = movements.movement.Movement(
            airplane_movements=[problem.movement.airplane_movements[0]],
            operating_point_movement=this_operating_point_movement,
            num_steps=problem.movement.num_steps,
        )

        this_problem = problems.UnsteadyProblem(
            movement=this_movement, only_final_results=True
        )

        this_solver = (
            unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                unsteady_problem=this_problem
            )
        )

        this_solver.run(logging_level="Critical", prescribed_wake=True)

        finalForceCoefficients_W = this_solver.unsteady_problem.finalForceCoefficients_W
        finalMomentCoefficients_W_Cg = (
            this_solver.unsteady_problem.finalMomentCoefficients_W_CgP1
        )

        netForceCoefficients_W = np.linalg.norm(
            finalForceCoefficients_W + externalForceCoefficients_W
        )
        netMomentCoefficients_W_Cg = np.linalg.norm(finalMomentCoefficients_W_Cg)

        objective = (netForceCoefficients_W + netMomentCoefficients_W_Cg) / 2

        v_str = str(round(vCg__E, 2))
        a_str = str(round(alpha, 2))
        b_str = str(round(beta, 2))
        f_str = str(round(externalFX_W, 2))
        o_str = str(round(objective, 3))

        state_msg = (
            "\tState: vCg__E="
            + v_str
            + ", alpha="
            + a_str
            + ", beta="
            + b_str
            + ", externalFX_W="
            + f_str
        )
        obj_msg = "\t\tObjective: " + o_str

        trim_logger.info(state_msg)
        trim_logger.info(obj_msg)

        if objective < objective_cut_off:
            raise StopIteration

        current_arguments.clear()
        current_arguments.extend([np.nan, np.nan, np.nan, np.nan])

        return objective

    initial_guess = np.array([baseVCg__E, base_alpha, base_beta, baseExternalFX_W])
    bounds: Sequence[tuple[float, float]] = [
        (boundsVCg__E[0], boundsVCg__E[1]),
        (alpha_bounds[0], alpha_bounds[1]),
        (beta_bounds[0], beta_bounds[1]),
        (boundsExternalFX_W[0], boundsExternalFX_W[1]),
    ]

    trim_logger.info("Starting local search.")
    try:
        options: Any = {"maxfun": num_calls, "eps": 0.01}
        scipy.optimize.minimize(
            fun=objective_function,
            x0=initial_guess,
            bounds=bounds,
            method="L-BFGS-B",
            options=options,
        )
    except StopIteration:
        trim_logger.info("Acceptable value reached with local search.")
        return (
            current_arguments[0],
            current_arguments[1],
            current_arguments[2],
            current_arguments[3],
        )

    trim_logger.warning(
        "No acceptable value reached with local search. Starting global search."
    )
    try:
        options: Any = {"maxfun": num_calls, "eps": 0.01}
        minimizer_kwargs: Any = {
            "method": "L-BFGS-B",
            "options": options,
        }
        scipy.optimize.dual_annealing(
            func=objective_function,
            bounds=bounds,
            x0=initial_guess,
            maxfun=num_calls,
            minimizer_kwargs=minimizer_kwargs,
            seed=_seed,
        )
    except StopIteration:
        trim_logger.info("Acceptable global minima found.")
        return (
            current_arguments[0],
            current_arguments[1],
            current_arguments[2],
            current_arguments[3],
        )

    trim_logger.critical(
        "No trim condition found. Try increasing the bounds and the maximum number of "
        "iterations."
    )
    return None, None, None, None
