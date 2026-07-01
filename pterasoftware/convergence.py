"""Contains functions for analyzing the convergence of SteadyProblems and
UnsteadyProblems.

**Contains the following classes:**

None

**Contains the following functions:**

analyze_steady_convergence: Finds the converged parameters of a SteadyProblem solved
using a given steady solver.

analyze_unsteady_convergence: Finds the converged parameters of an UnsteadyProblem
solved using the UnsteadyRingVortexLatticeMethodSolver.
"""

from __future__ import annotations

import time
from collections.abc import Callable

import numpy as np

from . import (
    _logging,
    _parameter_validation,
    geometry,
    movements,
    problems,
    steady_horseshoe_vortex_lattice_method,
    steady_ring_vortex_lattice_method,
    unsteady_ring_vortex_lattice_method,
)

_logger = _logging.get_logger("convergence")

# The labels of the six load coefficients that convergence is checked against, in the
# order they are stored: the three force coefficients followed by the three moment
# coefficients.
_COEFFICIENT_LABELS = ("cFX", "cFY", "cFZ", "cMX", "cMY", "cMZ")


# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
def analyze_steady_convergence(
    ref_problem: problems.SteadyProblem,
    solver_type: str,
    panel_aspect_ratio_bounds: tuple[int, int] = (4, 1),
    num_chordwise_panels_bounds: tuple[int, int] = (3, 12),
    rtol: float | int = 0.05,
    atol: float | int = 0.001,
    coefficient_mask: tuple[bool, bool, bool, bool, bool, bool] | None = None,
    resolve_converged_solver: bool | np.bool_ = False,
) -> (
    tuple[
        int,
        int,
        steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
        | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
        | None,
    ]
    | tuple[None, None, None]
):
    """Finds the converged parameters of a SteadyProblem solved using a given steady
    solver.

    **Procedure:**

    Convergence is found by varying the SteadyProblem's Airplanes' Panels' aspect ratios
    their Wings' numbers of chordwise Panels. These values are iterated over via two
    nested for loops (with the number of chordwise Panels as the inner loop).

    With each new combination of these values, the SteadyProblem is solved, and each
    Airplane's six load coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ) are stored. Then,
    convergence is checked per coefficient between this iteration and the iterations
    with incrementally coarser meshes in the two parameter directions (Panel aspect
    ratio and number of chordwise Panels). A coefficient is converged in a parameter
    direction when its absolute change from the coarser iteration is at most atol + rtol
    * max(abs(this), abs(coarser)). A parameter direction is converged when every
    unmasked coefficient of every Airplane is converged in that direction, so a multi-
    Airplane study converges only once all its Airplanes have.

    If an iteration is converged in both parameter directions, then we exit the nested
    for loops and return the converged parameters. However, the converged parameters are
    actually the values incrementally coarser than the final values, because refining
    from the coarser values to the final ones changed every unmasked coefficient by
    within the tolerance.

    **Notes:**

    There are two edge cases to this function. The first is if the user inputs equal
    values for the coarsest and finest values of either the Panel aspect ratio or the
    number of chordwise Panels (e.g. panel_aspect_ratio_bounds=(2, 2)). Then, this
    parameter will not be iterated over, and convergence will only be checked for the
    other parameter.

    The second edge case happens if the Panel aspect ratio has not converged at a value
    of 1. This is the gold standard value for Panel aspect ratio, so this function will
    return 1 for the converged value of Panel aspect ratio. In the code below, this
    state is referred to as a "saturated" Panel aspect ratio case.

    Each Wing is refined according to how its spanwise mesh was defined. A non-edge-
    defined Wing (one built from WingCrossSections) is refined by sweeping its number of
    spanwise Panels to hit the target Panel aspect ratio, holding its WingCrossSections
    fixed. An edge-defined Wing (one built from edge curves with Wing.from_edge_points)
    is refined by resampling its stored edge curves into the number of WingCrossSections
    that hits the target Panel aspect ratio, preserving its planform shape. An Airplane
    may hold both kinds of Wing at once. A Wing that has been exploded into single-panel
    strips cannot be refined and is rejected.

    :param ref_problem: The SteadyProblem whose converged parameters will be found.
    :param solver_type: Determines what type of steady solver will be used to analyze
        the SteadyProblem. The options are "steady horseshoe vortex lattice method" and
        "steady ring vortex lattice method".
    :param panel_aspect_ratio_bounds: A tuple of two ints, in descending order, that
        determines the range of Panel aspect ratios to consider, from largest to
        smallest. This value dictates the Panels' average y component length (in wing
        cross section parent axes) divided their average x component width (in wing
        cross section parent axes). Historically, these values range between 5 and 1.
        Values above 5 can be used for a coarser mesh, but the minimum value cannot be
        less than 1. The default is (4, 1).
    :param num_chordwise_panels_bounds: A tuple of two ints, in ascending order, that
        determines the range of values to use for the Wings' numbers of chordwise
        panels. The default is (3, 12).
    :param rtol: A positive number (int or float) giving the relative tolerance for the
        per-coefficient convergence check. A coefficient is converged in a parameter
        direction when its absolute change from the incrementally coarser iteration is
        at most atol + rtol * max(abs(this), abs(coarser)). Set this smaller for a
        stricter convergence. Values are converted to floats internally. The default is
        0.05.
    :param atol: A positive number (int or float) giving the absolute tolerance floor
        for the per-coefficient convergence check. It keeps coefficients that sit near
        zero from being held to an unreachable relative tolerance. Values are converted
        to floats internally. The default is 0.001.
    :param coefficient_mask: A tuple of six bools that determines which of the six load
        coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ, in that order) must converge, or
        None to require all six. At least one element must be True. Use this to ignore
        coefficients that are physically irrelevant to the analysis. The default is
        None.
    :param resolve_converged_solver: A bool for whether to recreate and run the solver
        at the converged parameters and return it. Because finding convergence is
        expensive, this defaults to False, in which case the returned solver is None.
        When True, the solver is rebuilt at the converged parameters (which are
        frequently coarser than the last iteration run) and run with streamlines
        calculated, so the returned solver is ready to use. The default is False.
    :return: A tuple of two ints and a solver, or a tuple of three Nones. In order, the
        first two elements are the converged Panel aspect ratio and the converged number
        of chordwise Panels. The third element is the converged solver if
        resolve_converged_solver is True, otherwise None. If the function could not find
        a set of converged parameters, it returns (None, None, None).
    """
    # Validate the ref_problem parameter.
    if not isinstance(ref_problem, problems.SteadyProblem):
        raise TypeError("ref_problem must be a SteadyProblem.")

    # Validate the solver_type parameter.
    if solver_type not in (
        "steady horseshoe vortex lattice method",
        "steady ring vortex lattice method",
    ):
        raise ValueError(
            'solver_type must be either "steady horseshoe vortex lattice method" or '
            '"steady ring vortex lattice method".'
        )

    # Validate the panel_aspect_ratio_bounds parameter.
    _validate_panel_aspect_ratio_bounds(panel_aspect_ratio_bounds)

    # Validate the num_chordwise_panels_bounds parameter.
    _validate_num_chordwise_panels_bounds(num_chordwise_panels_bounds)

    # Validate the rtol parameter.
    rtol = _parameter_validation.number_in_range_return_float(
        rtol, "rtol", min_val=0.0, min_inclusive=False
    )

    # Validate the atol parameter.
    atol = _parameter_validation.number_in_range_return_float(
        atol, "atol", min_val=0.0, min_inclusive=False
    )

    # Validate the coefficient_mask parameter.
    coefficient_mask_array = _validate_coefficient_mask(coefficient_mask)

    # Validate the resolve_converged_solver parameter.
    resolve_converged_solver = _parameter_validation.boolLike_return_bool(
        resolve_converged_solver, "resolve_converged_solver"
    )

    run_start_time = time.time()
    _logger.info("Beginning convergence analysis...")

    ref_airplanes = ref_problem.airplanes

    # Reject any Wing that cannot be refined (any spanwise mesh other than trapezoidal
    # or edge-defined).
    _reject_unrefinable_wings(ref_airplanes, "analyze_steady_convergence")

    # Create lists containing each Panel aspect ratio and each number of chordwise
    # Panels to test.
    panel_aspect_ratios_list = list(
        range(panel_aspect_ratio_bounds[0], panel_aspect_ratio_bounds[1] - 1, -1)
    )
    num_chordwise_panels_list = list(
        range(num_chordwise_panels_bounds[0], num_chordwise_panels_bounds[1] + 1)
    )

    # Initialize some empty ndarrays to hold variables related to each iteration.
    # Going forward, an "iteration" refers to a SteadyProblem containing one of the
    # combinations of Panel aspect ratio and number of chordwise Panels.
    iter_times = np.zeros(
        (len(panel_aspect_ratios_list), len(num_chordwise_panels_list)), dtype=float
    )

    # Each iteration stores the six load coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ) of
    # each Airplane, in that order, along the last axis.
    coefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplanes),
            6,
        ),
        dtype=float,
    )

    iteration = 0
    num_iterations = len(panel_aspect_ratios_list) * len(num_chordwise_panels_list)

    # This is a cache to store previously calculated numbers of spanwise Panels for
    # specific combinations of parameters to avoid redundant calculations. The key is
    # a tuple of 5 ints: ar_id, chord_id, ref_airplane_id, ref_wing_id,
    # ref_wing_cross_section_id,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int] = {}

    # This is the analogous cache for edge-defined Wings, which are refined by their
    # number of WingCrossSections rather than by their number of spanwise Panels. The key
    # is a tuple of 4 ints: ar_id, chord_id, ref_airplane_id, ref_wing_id.
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int] = {}

    # Begin iterating through the outer loop of Panel aspect ratios.
    for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):
        _logger.info("\tPanel aspect ratio: " + str(panel_aspect_ratio))

        # Begin iterating through the inner loop of number of chordwise Panels.
        for chord_id, num_chordwise_panels in enumerate(num_chordwise_panels_list):
            _logger.info("\t\tChordwise Panels: " + str(num_chordwise_panels))

            iteration += 1
            _logger.info(
                "\t\t\tIteration number: " + str(iteration) + "/" + str(num_iterations)
            )

            # Build this iteration's SteadyProblem with the current Panel aspect ratio
            # and number of chordwise Panels.
            this_problem = _build_steady_problem(
                ar_id,
                chord_id,
                panel_aspect_ratio,
                num_chordwise_panels,
                ref_problem,
                num_spanwise_panels_cache,
                num_wing_cross_sections_cache,
            )
            these_airplanes = this_problem.airplanes

            # Create this iteration's steady solver based on the type specified.
            this_solver: (
                steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
                | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
            )
            if solver_type == "steady horseshoe vortex lattice method":
                this_solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                    steady_problem=this_problem,
                )
            else:
                this_solver = steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
                    steady_problem=this_problem,
                )

            _logger.info("\t\t\tStarting simulation...")

            # Run the steady solver and time how long it takes to execute. Skip the
            # streamline trace since it does not affect convergence metrics.
            iter_start = time.time()
            this_solver.run(calculate_streamlines=False)
            iter_stop = time.time()
            this_iter_time = iter_stop - iter_start

            # Create and fill an ndarray with each of this iteration's Airplanes' six
            # load coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ).
            theseCoefficients = np.zeros((len(these_airplanes), 6), dtype=float)

            for airplane_id, airplane in enumerate(these_airplanes):
                _forceCoefficients_W = airplane.forceCoefficients_W
                assert _forceCoefficients_W is not None

                _momentCoefficients_W_CgP1 = airplane.momentCoefficients_W_CgP1
                assert _momentCoefficients_W_CgP1 is not None

                theseCoefficients[airplane_id, 0:3] = _forceCoefficients_W
                theseCoefficients[airplane_id, 3:6] = _momentCoefficients_W_CgP1

            # Populate the ndarray that stores information from all the iterations with
            # the data from this iteration.
            coefficients[ar_id, chord_id, :, :] = theseCoefficients
            iter_times[ar_id, chord_id] = this_iter_time

            _logger.info(
                "\t\t\tSimulation completed in " + str(round(this_iter_time, 3)) + " s"
            )

            # Check per-coefficient convergence in each parameter direction against the
            # incrementally coarser iteration. A direction cannot be checked on its first
            # value, where there is no coarser iteration to compare against, so it starts
            # as not converged.
            ar_converged = False
            chord_converged = False

            # If this isn't the first Panel aspect ratio, check convergence in the Panel
            # aspect ratio direction.
            if ar_id > 0:
                ar_converged, ar_metric, ar_limiting_id = (
                    _check_coefficient_convergence(
                        theseCoefficients,
                        coefficients[ar_id - 1, chord_id, :, :],
                        rtol,
                        atol,
                        coefficient_mask_array,
                    )
                )

                _logger.info(
                    "\t\t\tPanel aspect ratio convergence metric: "
                    + str(round(ar_metric, 2))
                    + "% (limiting coefficient: "
                    + _COEFFICIENT_LABELS[ar_limiting_id]
                    + ")"
                )
            else:
                _logger.info(
                    "\t\t\tPanel aspect ratio convergence metric: not yet checked"
                )

            # If this isn't the first number of chordwise Panels, check convergence in
            # the number of chordwise Panels direction.
            if chord_id > 0:
                chord_converged, chord_metric, chord_limiting_id = (
                    _check_coefficient_convergence(
                        theseCoefficients,
                        coefficients[ar_id, chord_id - 1, :, :],
                        rtol,
                        atol,
                        coefficient_mask_array,
                    )
                )

                _logger.info(
                    "\t\t\tNumber of chordwise Panels convergence metric: "
                    + str(round(chord_metric, 2))
                    + "% (limiting coefficient: "
                    + _COEFFICIENT_LABELS[chord_limiting_id]
                    + ")"
                )
            else:
                _logger.info(
                    "\t\t\tNumber of chordwise Panels convergence metric: not yet "
                    "checked"
                )

            # Consider the Panel aspect ratio value to be saturated if it is equal to
            # 1. This is because a Panel aspect ratio of 1 is considered the maximum
            # degree of fineness.
            ar_saturated = panel_aspect_ratio == 1

            # Check if only one value for either the Panel aspect ratio or the number
            # of chordwise Panels were specified.
            single_ar = len(panel_aspect_ratios_list) == 1
            single_chord = len(num_chordwise_panels_list) == 1

            # Consider each convergence parameter to have "passed" if it is
            # converged, single, or saturated.
            ar_passed = ar_converged or single_ar or ar_saturated
            chord_passed = chord_converged or single_chord

            # If both convergence parameters have passed, then a converged or
            # semi-converged combination of parameters has been found and will be
            # returned.
            if ar_passed and chord_passed:
                converged_ar_id = _converged_parameter_id(
                    ar_id, single_ar, ar_converged
                )
                converged_chord_id = _converged_parameter_id(
                    chord_id, single_chord, chord_converged
                )

                converged_aspect_ratio = panel_aspect_ratios_list[converged_ar_id]
                converged_chordwise_panels = num_chordwise_panels_list[
                    converged_chord_id
                ]
                converged_iter_time = float(
                    iter_times[converged_ar_id, converged_chord_id]
                )

                if single_ar or single_chord:
                    _logger.info("The analysis found a semi-converged case:")
                    if single_ar:
                        _logger.warning(
                            "Panel aspect ratio convergence was not checked"
                        )
                    if single_chord:
                        _logger.warning("Chordwise panels convergence was not checked")
                else:
                    _logger.info("The analysis found a converged case:")

                _logger.info("\tPanel aspect ratio: " + str(converged_aspect_ratio))
                _logger.info("\tChordwise Panels: " + str(converged_chordwise_panels))
                _logger.info(
                    "\tSimulation time: " + str(round(converged_iter_time, 3)) + " s"
                )
                _logger.info("\tSpanwise Panels:")
                for airplane_id, airplane in enumerate(ref_airplanes):
                    _logger.info("\t\t" + airplane.name + ":")
                    for wing_id, wing in enumerate(airplane.wings):
                        _logger.info("\t\t\t" + wing.name + ":")

                        # An edge-defined Wing is refined by its number of
                        # WingCrossSections rather than its number of spanwise Panels, so
                        # report that instead.
                        if wing.spanwise_mesh == "edge_defined":
                            num_wing_cross_sections_key = (
                                converged_ar_id,
                                converged_chord_id,
                                airplane_id,
                                wing_id,
                            )
                            _logger.info(
                                "\t\t\t\tWingCrossSections: "
                                + str(
                                    num_wing_cross_sections_cache[
                                        num_wing_cross_sections_key
                                    ]
                                )
                            )
                            continue

                        for wing_cross_section_id, wing_cross_section in enumerate(
                            wing.wing_cross_sections
                        ):
                            if (
                                wing_cross_section_id
                                < len(wing.wing_cross_sections) - 1
                            ):
                                # Not the last WingCrossSection, retrieve from cache.
                                num_spanwise_panels_key = (
                                    converged_ar_id,
                                    converged_chord_id,
                                    airplane_id,
                                    wing_id,
                                    wing_cross_section_id,
                                )
                                num_spanwise_panels = num_spanwise_panels_cache[
                                    num_spanwise_panels_key
                                ]
                            else:
                                # Last WingCrossSection.
                                num_spanwise_panels = None
                            _logger.info(
                                "\t\t\t\tWingCrossSection "
                                + str(wing_cross_section_id + 1)
                                + ": "
                                + str(num_spanwise_panels)
                            )

                # If requested, recreate and run the solver at the converged
                # parameters so it can be returned. The converged parameters are
                # frequently coarser than this iteration's, so the converged solver is
                # rebuilt rather than reusing this iteration's finer solver.
                converged_solver: (
                    steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
                    | steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
                    | None
                ) = None
                if resolve_converged_solver:
                    _logger.info("Recreating and running the converged solver...")
                    converged_problem = _build_steady_problem(
                        converged_ar_id,
                        converged_chord_id,
                        converged_aspect_ratio,
                        converged_chordwise_panels,
                        ref_problem,
                        num_spanwise_panels_cache,
                        num_wing_cross_sections_cache,
                    )
                    if solver_type == "steady horseshoe vortex lattice method":
                        converged_solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                            steady_problem=converged_problem,
                        )
                    else:
                        converged_solver = steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
                            steady_problem=converged_problem,
                        )
                    converged_solver.run(calculate_streamlines=True)

                _logger.info(
                    "Convergence analysis completed in "
                    + str(round(time.time() - run_start_time, 3))
                    + " s"
                )

                return (
                    converged_aspect_ratio,
                    converged_chordwise_panels,
                    converged_solver,
                )

    # If all iterations have been checked and none of them resulted in both
    # convergence parameters passing, then indicate that no converged case was found
    # and return values of None for the converged parameters.
    _logger.info("The analysis did not find a converged case within the given bounds")
    _logger.info(
        "Convergence analysis completed in "
        + str(round(time.time() - run_start_time, 3))
        + " s"
    )
    return None, None, None


# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
def analyze_unsteady_convergence(
    ref_problem: problems.UnsteadyProblem,
    prescribed_wake: bool | np.bool_ = True,
    free_wake: bool | np.bool_ = True,
    num_cycles_bounds: tuple[int, int] | None = None,
    num_chords_bounds: tuple[int, int] | None = None,
    panel_aspect_ratio_bounds: tuple[int, int] = (4, 1),
    num_chordwise_panels_bounds: tuple[int, int] = (3, 12),
    rtol: float | int = 0.05,
    atol: float | int = 0.001,
    coefficient_mask: tuple[bool, bool, bool, bool, bool, bool] | None = None,
    show_solver_progress: bool | np.bool_ = True,
    resolve_converged_solver: bool | np.bool_ = False,
) -> (
    tuple[
        bool,
        int,
        int,
        int,
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
        | None,
    ]
    | tuple[None, None, None, None, None]
):
    """Finds the converged parameters of an UnsteadyProblem solved using the
    UnsteadyRingVortexLatticeMethodSolver.

    **Procedure:**

    Convergence is found by varying the UnsteadyRingVortexLatticeMethodSolver's wake
    state (prescribed or free), the final length of the UnsteadyProblem's wake (in
    number of chord lengths for static geometry or number of maximum-period motion
    cycles for variable geometry), the Airplanes' Wings' Panels' aspect ratios, and the
    Airplanes' Wings' numbers of chordwise Panels. These values are iterated over via
    four nested loops. The outermost loop is the wake state. The next loop is the wake
    length. The loop after that is the Panel aspect ratios, and the innermost loop is
    the number of chordwise Panels.

    With each new combination of these values, the UnsteadyProblem is solved, and each
    Airplane's six final load coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ) are stored. As
    this function deals with UnsteadyProblems, it considers the final load coefficients
    to be the final cycle's mean load coefficients (the signed time-average over the
    final cycle) for UnsteadyProblems with variable geometry, and the final time step's
    load coefficients for static geometry cases. Then, convergence is checked per
    coefficient between this iteration and the iterations with incrementally coarser
    meshes in all four parameter directions (wake state, wake length, Panel aspect
    ratio, and number of chordwise Panels). A coefficient is converged in a parameter
    direction when its absolute change from the coarser iteration is at most atol + rtol
    * max(abs(this), abs(coarser)). A parameter direction is converged when every
    unmasked coefficient of every Airplane is converged in that direction, so a multi-
    Airplane study converges only once all its Airplanes have.

    If an iteration is converged in all four parameter directions, then we exit the
    nested for loops and return the converged parameters. However, the converged
    parameters are actually the values incrementally coarser than the final values,
    because refining from the coarser values to the final ones changed every unmasked
    coefficient by within the tolerance.

    **Notes:**

    There are two edge cases to this function. The first occurs when the user indicates
    that they only want check a single value for any of the four parameters (e.g.
    panel_aspect_ratio_bounds=(2, 2), both prescribed_wake=True and free_wake=False,
    etc.). Then, this parameter will not be iterated over, and convergence will only be
    checked for the other parameters.

    The second edge case happens if the Panel aspect ratio has not converged at a value
    of 1 or if the wake state hasn't converged once it's set to a free wake. These
    conditions are the gold standards for Panel aspect ratio and wake state, so this
    function will return 1 for the converged value of Panel aspect ratio and a free wake
    for the converged wake state. In the code below, this state is referred to as a
    "saturated" Panel aspect ratio or wake state.

    Each Wing is refined according to how its spanwise mesh was defined. A non-edge-
    defined Wing (one built from WingCrossSections) is refined by sweeping its number of
    spanwise Panels to hit the target Panel aspect ratio, holding its WingCrossSections
    fixed. An edge-defined Wing (one built from edge curves with Wing.from_edge_points)
    is refined by resampling its stored edge curves into the number of WingCrossSections
    that hits the target Panel aspect ratio, preserving its planform shape. An Airplane
    may hold both kinds of Wing at once. A Wing that has been exploded into single-panel
    strips cannot be refined and is rejected. Because resampling an edge-defined Wing
    changes its number of WingCrossSections, its WingCrossSectionMovements must all be
    static (the WingMovement's own flapping and sweeping motion is preserved); a non-
    static WingCrossSectionMovement on an edge-defined Wing is rejected.

    The time step is not a convergence parameter. Each iteration instead optimizes its
    own delta_time for its mesh with Movement's iterative optimizer, so every solve uses
    the time step that is correct for its geometry and motion. The optimum depends only
    on the mesh, so it is computed once per Panel aspect ratio and number of chordwise
    Panels and reused across every wake state and wake length at that mesh.

    :param ref_problem: The UnsteadyProblem whose converged parameters will be found.
        This must be a standard UnsteadyProblem, not a FreeFlightUnsteadyProblem or an
        AeroelasticUnsteadyProblem, neither of which is supported.
    :param prescribed_wake: Determines if a prescribed wake state should be analyzed. If
        this parameter is False, then the ``free_wake`` parameter must be set to True.
        Can be a bool or a numpy bool and will be converted to a bool internally. The
        default is True.
    :param free_wake: Determines if a free wake state should be analyzed. If this
        parameter is False, then the ``prescribed_wake`` parameter must be set to True.
        Can be a bool or a numpy bool and will be converted to a bool internally. The
        default is True.
    :param num_cycles_bounds: For problems with non static geometry, determines the
        range of wake lengths (measured in number of maximum-period motion cycles) to
        simulate. For problems with static geometry, this must be None, and the
        ``num_chords_bounds`` parameter will control the range of wake lengths instead.
        Otherwise, it must be a tuple of two positive ints with the first value less
        than or equal to the second value. Reasonable values range from 1 to 10,
        depending strongly on the Strouhal number. The default is None.
    :param num_chords_bounds: For problems with static geometry, determines the range of
        wake lengths (measured in number of reference chords) to simulate. For problems
        with non static geometry, it must be None, and the ``num_cycles_bounds``
        parameter will control the wake length instead. Otherwise, it must be a tuple of
        two positive ints with the first value less than or equal to the second value.
        Reasonable values range from 3 to 20. The default is None.
    :param panel_aspect_ratio_bounds: A tuple of two ints, in descending order, that
        determines the range of Panel aspect ratios to consider, from largest to
        smallest. This value dictates the Panels' average y component length (in wing
        cross section parent axes) divided their average x component width (in wing
        cross section parent axes). Historically, these values range between 5 and 1.
        Values above 5 can be used for a coarser mesh, but the minimum value cannot be
        less than 1. The default is (4, 1).
    :param num_chordwise_panels_bounds: A tuple of two ints, in ascending order, that
        determines the range of values to use for the Wings' numbers of chordwise
        panels. The default is (3, 12).
    :param rtol: A positive number (int or float) giving the relative tolerance for the
        per-coefficient convergence check. A coefficient is converged in a parameter
        direction when its absolute change from the incrementally coarser iteration is
        at most atol + rtol * max(abs(this), abs(coarser)). Set this smaller for a
        stricter convergence. Values are converted to floats internally. The default is
        0.05.
    :param atol: A positive number (int or float) giving the absolute tolerance floor
        for the per-coefficient convergence check. It keeps coefficients that sit near
        zero from being held to an unreachable relative tolerance. Values are converted
        to floats internally. The default is 0.001.
    :param coefficient_mask: A tuple of six bools that determines which of the six load
        coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ, in that order) must converge, or
        None to require all six. At least one element must be True. Use this to ignore
        coefficients that are physically irrelevant to the analysis. The default is
        None.
    :param show_solver_progress: Set this to True to show the TQDM progress bar during
        each run of the unsteady solver. For showing progress bars and displaying log
        statements, set up logging using the setup_logging function. It can be a bool or
        a numpy bool and will be converted internally to a bool. The default is True.
    :param resolve_converged_solver: A bool for whether to recreate and run the solver
        at the converged parameters and return it. Because finding convergence is
        expensive, this defaults to False, in which case the returned solver is None.
        When True, the solver is rebuilt at the converged parameters (which are
        frequently coarser than the last iteration run) and run with streamlines
        calculated, so the returned solver is ready to use. The default is False.
    :return: A tuple of one bool, three ints, and a solver, or a tuple of five Nones. In
        order, the first four elements are the converged wake state (prescribed=True and
        free=False), the converged wake length (in number of cycles for non static
        geometries and number of chords for static geometries), the converged Panel
        aspect ratio, and the converged number of chordwise Panels. The fifth element is
        the converged solver if resolve_converged_solver is True, otherwise None. If the
        function could not find a set of converged parameters, it returns (None, None,
        None, None, None).
    """
    # Validate the ref_problem parameter.
    if not isinstance(ref_problem, problems.UnsteadyProblem):
        raise TypeError(
            "ref_problem must be a standard UnsteadyProblem, not a "
            "FreeFlightUnsteadyProblem or an AeroelasticUnsteadyProblem."
        )

    # Validate the wake type parameters.
    prescribed_wake = _parameter_validation.boolLike_return_bool(
        prescribed_wake, "prescribed_wake"
    )
    free_wake = _parameter_validation.boolLike_return_bool(free_wake, "free_wake")
    if not (prescribed_wake or free_wake):
        raise ValueError("At least one of prescribed_wake or free_wake must be True.")

    # Validate the wake length bounds parameters.
    ref_movement: movements.movement.Movement = ref_problem.movement
    static = ref_movement.static
    if static:
        if not num_cycles_bounds is None:
            raise ValueError(
                "num_cycles_bounds must be None for UnsteadyProblems "
                "with static geometry."
            )
        if not (isinstance(num_chords_bounds, tuple) and len(num_chords_bounds) == 2):
            raise TypeError("num_chords_bounds must be a tuple with length 2.")
        if not all(isinstance(bound, int) for bound in num_chords_bounds):
            raise TypeError("Both values in num_chords_bounds must be ints.")
        if num_chords_bounds[1] < num_chords_bounds[0]:
            raise ValueError(
                "The second value in num_chords_bounds must be greater than or equal "
                "to the first value."
            )
        if num_chords_bounds[1] <= 0:
            raise ValueError("Both values in num_chords_bounds must be positive.")
    else:
        if not num_chords_bounds is None:
            raise ValueError(
                "num_chords_bounds must be None for UnsteadyProblems "
                "with variable geometry."
            )
        if not (isinstance(num_cycles_bounds, tuple) and len(num_cycles_bounds) == 2):
            raise TypeError("num_cycles_bounds must be a tuple with length 2.")
        if not all(isinstance(bound, int) for bound in num_cycles_bounds):
            raise TypeError("Both values in num_cycles_bounds must be ints.")
        if num_cycles_bounds[1] < num_cycles_bounds[0]:
            raise ValueError(
                "The second value in num_cycles_bounds must be greater than or equal "
                "to the first value."
            )
        if num_cycles_bounds[1] <= 0:
            raise ValueError("Both values in num_cycles_bounds must be positive.")

    # Validate the panel_aspect_ratio_bounds parameter.
    _validate_panel_aspect_ratio_bounds(panel_aspect_ratio_bounds)

    # Validate the num_chordwise_panels_bounds parameter.
    _validate_num_chordwise_panels_bounds(num_chordwise_panels_bounds)

    # Validate the rtol parameter.
    rtol = _parameter_validation.number_in_range_return_float(
        rtol, "rtol", min_val=0.0, min_inclusive=False
    )

    # Validate the atol parameter.
    atol = _parameter_validation.number_in_range_return_float(
        atol, "atol", min_val=0.0, min_inclusive=False
    )

    # Validate the coefficient_mask parameter.
    coefficient_mask_array = _validate_coefficient_mask(coefficient_mask)

    # Validate the show_solver_progress parameter.
    show_solver_progress = _parameter_validation.boolLike_return_bool(
        show_solver_progress, "show_solver_progress"
    )

    # Validate the resolve_converged_solver parameter.
    resolve_converged_solver = _parameter_validation.boolLike_return_bool(
        resolve_converged_solver, "resolve_converged_solver"
    )

    run_start_time = time.time()
    _logger.info("Beginning convergence analysis...")

    ref_airplane_movements = ref_movement.airplane_movements

    # Reject any Wing that cannot be refined (any spanwise mesh other than trapezoidal or
    # edge-defined).
    _reject_unrefinable_wings(
        tuple(
            ref_airplane_movement.base_airplane
            for ref_airplane_movement in ref_airplane_movements
        ),
        "analyze_unsteady_convergence",
    )

    # Create the list of wake states to iterate over.
    wake_list = []
    if prescribed_wake:
        wake_list.append(True)
    if free_wake:
        wake_list.append(False)

    # Create the list of wake lengths to iterate over.
    if static:
        assert num_chords_bounds is not None
        wake_lengths_list = list(range(num_chords_bounds[0], num_chords_bounds[1] + 1))
    else:
        assert num_cycles_bounds is not None
        wake_lengths_list = list(range(num_cycles_bounds[0], num_cycles_bounds[1] + 1))

    # Create the lists of Panel aspect ratios and number of chordwise Panels to
    # iterate over.
    panel_aspect_ratios_list = list(
        range(panel_aspect_ratio_bounds[0], panel_aspect_ratio_bounds[1] - 1, -1)
    )
    num_chordwise_panels_list = list(
        range(num_chordwise_panels_bounds[0], num_chordwise_panels_bounds[1] + 1)
    )

    # Initialize some empty ndarrays to hold variables regarding each iteration.
    # Going forward, an "iteration" refers to an UnsteadyProblem containing one of
    # the combinations of the wake state, wake length, Panel aspect ratio, and number
    # of chordwise Panels.
    iter_times = np.zeros(
        (
            len(wake_list),
            len(wake_lengths_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
        ),
        dtype=float,
    )
    # Each iteration stores the six final load coefficients (cFX, cFY, cFZ, cMX, cMY,
    # cMZ) of each Airplane, in that order, along the last axis.
    finalCoefficients = np.zeros(
        (
            len(wake_list),
            len(wake_lengths_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplane_movements),
            6,
        ),
        dtype=float,
    )

    iteration = 0
    num_iterations = (
        len(wake_list)
        * len(wake_lengths_list)
        * len(panel_aspect_ratios_list)
        * len(num_chordwise_panels_list)
    )

    # This is a cache to store previously calculated numbers of spanwise Panels for
    # specific combinations of parameters to avoid redundant calculations. The key is
    # a tuple of 5 ints: ar_id, chord_id, ref_base_airplane_id, ref_base_wing_id,
    # ref_base_wing_cross_section_id,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int] = {}

    # This is the analogous cache for edge-defined Wings, which are refined by their
    # number of WingCrossSections rather than by their number of spanwise Panels. The key
    # is a tuple of 4 ints: ar_id, chord_id, ref_base_airplane_id, ref_base_wing_id.
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int] = {}

    # This caches the optimized delta_time for each mesh, keyed on a tuple of 2 ints:
    # ar_id, chord_id. The optimum depends only on the mesh (the geometry and motion),
    # not on the wake state or wake length, so it is reused across every wake-state and
    # wake-length combination at that mesh.
    delta_time_cache: dict[tuple[int, int], float] = {}

    # Begin iterating through the outermost loop of wake states.
    for wake_id, wake in enumerate(wake_list):
        if wake:
            _logger.info("\tWake type: prescribed")
        else:
            _logger.info("\tWake type: free")

        # Begin iterating through the second loop of wake lengths.
        for length_id, wake_length in enumerate(wake_lengths_list):
            if static:
                _logger.info("\t\tChord lengths: " + str(wake_length))
            else:
                _logger.info("\t\tCycles: " + str(wake_length))

            # Begin iterating through the third loop of Panel aspect ratios.
            for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):
                _logger.info("\t\t\tPanel aspect ratio: " + str(panel_aspect_ratio))

                # Begin iterating through the innermost loop of number of chordwise
                # Panels.
                for chord_id, num_chordwise_panels in enumerate(
                    num_chordwise_panels_list
                ):
                    _logger.info(
                        "\t\t\t\tChordwise Panels: " + str(num_chordwise_panels)
                    )

                    iteration += 1
                    _logger.info(
                        "\t\t\t\t\tIteration number: "
                        + str(iteration)
                        + "/"
                        + str(num_iterations)
                    )

                    # Build this iteration's UnsteadyProblem with the current wake
                    # length, Panel aspect ratio, and number of chordwise Panels.
                    this_problem = _build_unsteady_problem(
                        ar_id,
                        chord_id,
                        panel_aspect_ratio,
                        num_chordwise_panels,
                        wake_length,
                        ref_problem,
                        num_spanwise_panels_cache,
                        num_wing_cross_sections_cache,
                        delta_time_cache,
                    )

                    # Create and run this iteration's
                    # UnsteadyRingVortexLatticeMethodSolver and time how long it
                    # takes to execute.
                    this_solver = unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        unsteady_problem=this_problem
                    )

                    _logger.info("\t\t\t\t\tStarting simulation...")

                    iter_start = time.time()
                    this_solver.run(
                        prescribed_wake=wake,
                        calculate_streamlines=False,
                        show_progress=show_solver_progress,
                    )
                    iter_stop = time.time()
                    this_iter_time = iter_stop - iter_start

                    # Create and fill an ndarray with each of this iteration's Airplanes'
                    # six final load coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ).
                    theseFinalCoefficients = np.zeros(
                        (len(ref_airplane_movements), 6), dtype=float
                    )

                    for airplane_id in range(len(ref_airplane_movements)):
                        # For static geometry, the final load coefficients are the final
                        # time step's. For variable geometry, they are the mean over the
                        # final cycle (the signed time-average).
                        if static:
                            theseFinalCoefficients[airplane_id, 0:3] = (
                                this_problem.finalForceCoefficients_W[airplane_id]
                            )
                            theseFinalCoefficients[airplane_id, 3:6] = (
                                this_problem.finalMomentCoefficients_W_CgP1[airplane_id]
                            )
                        else:
                            theseFinalCoefficients[airplane_id, 0:3] = (
                                this_problem.finalMeanForceCoefficients_W[airplane_id]
                            )
                            theseFinalCoefficients[airplane_id, 3:6] = (
                                this_problem.finalMeanMomentCoefficients_W_CgP1[
                                    airplane_id
                                ]
                            )

                    # Populate the ndarray that stores information from all the
                    # iterations with the data from this iteration.
                    finalCoefficients[wake_id, length_id, ar_id, chord_id, :, :] = (
                        theseFinalCoefficients
                    )
                    iter_times[wake_id, length_id, ar_id, chord_id] = this_iter_time

                    _logger.info(
                        "\t\t\t\t\tSimulation completed in "
                        + str(round(this_iter_time, 3))
                        + " s"
                    )

                    # Check per-coefficient convergence in each parameter direction
                    # against the incrementally coarser iteration. A direction cannot be
                    # checked on its first value, where there is no coarser iteration to
                    # compare against, so it starts as not converged.
                    wake_converged = False
                    length_converged = False
                    ar_converged = False
                    chord_converged = False

                    # If this isn't the first wake state, check convergence in the wake
                    # state direction.
                    if wake_id > 0:
                        wake_converged, wake_metric, wake_limiting_id = (
                            _check_coefficient_convergence(
                                theseFinalCoefficients,
                                finalCoefficients[
                                    wake_id - 1, length_id, ar_id, chord_id, :, :
                                ],
                                rtol,
                                atol,
                                coefficient_mask_array,
                            )
                        )

                        _logger.info(
                            "\t\t\t\t\tWake type convergence metric: "
                            + str(round(wake_metric, 2))
                            + "% (limiting coefficient: "
                            + _COEFFICIENT_LABELS[wake_limiting_id]
                            + ")"
                        )
                    else:
                        _logger.info(
                            "\t\t\t\t\tWake type convergence metric: not yet checked"
                        )

                    # If this isn't the first wake length, check convergence in the wake
                    # length direction.
                    if length_id > 0:
                        length_converged, length_metric, length_limiting_id = (
                            _check_coefficient_convergence(
                                theseFinalCoefficients,
                                finalCoefficients[
                                    wake_id, length_id - 1, ar_id, chord_id, :, :
                                ],
                                rtol,
                                atol,
                                coefficient_mask_array,
                            )
                        )

                        _logger.info(
                            "\t\t\t\t\tWake length convergence metric: "
                            + str(round(length_metric, 2))
                            + "% (limiting coefficient: "
                            + _COEFFICIENT_LABELS[length_limiting_id]
                            + ")"
                        )
                    else:
                        _logger.info(
                            "\t\t\t\t\tWake length convergence metric: not yet checked"
                        )

                    # If this isn't the first Panel aspect ratio, check convergence in
                    # the Panel aspect ratio direction.
                    if ar_id > 0:
                        ar_converged, ar_metric, ar_limiting_id = (
                            _check_coefficient_convergence(
                                theseFinalCoefficients,
                                finalCoefficients[
                                    wake_id, length_id, ar_id - 1, chord_id, :, :
                                ],
                                rtol,
                                atol,
                                coefficient_mask_array,
                            )
                        )

                        _logger.info(
                            "\t\t\t\t\tPanel aspect ratio convergence metric: "
                            + str(round(ar_metric, 2))
                            + "% (limiting coefficient: "
                            + _COEFFICIENT_LABELS[ar_limiting_id]
                            + ")"
                        )
                    else:
                        _logger.info(
                            "\t\t\t\t\tPanel aspect ratio convergence metric: not yet "
                            "checked"
                        )

                    # If this isn't the first number of chordwise Panels, check
                    # convergence in the number of chordwise Panels direction.
                    if chord_id > 0:
                        chord_converged, chord_metric, chord_limiting_id = (
                            _check_coefficient_convergence(
                                theseFinalCoefficients,
                                finalCoefficients[
                                    wake_id, length_id, ar_id, chord_id - 1, :, :
                                ],
                                rtol,
                                atol,
                                coefficient_mask_array,
                            )
                        )

                        _logger.info(
                            "\t\t\t\t\tNumber of chordwise Panels convergence metric: "
                            + str(round(chord_metric, 2))
                            + "% (limiting coefficient: "
                            + _COEFFICIENT_LABELS[chord_limiting_id]
                            + ")"
                        )
                    else:
                        _logger.info(
                            "\t\t\t\t\tNumber of chordwise Panels convergence metric: "
                            "not yet checked"
                        )

                    # Consider the Panel aspect ratio value to be saturated if it is
                    # equal to 1. This is because a Panel aspect ratio of 1 is
                    # considered the maximum degree of fineness. Consider the wake
                    # state to be saturated if it False (which corresponds to a free
                    # wake), as this is considered to be the most accurate wake state.
                    wake_saturated = not wake
                    ar_saturated = panel_aspect_ratio == 1

                    # Check if only one value was specified for any of the four
                    # convergence parameters.
                    single_wake = len(wake_list) == 1
                    single_length = len(wake_lengths_list) == 1
                    single_ar = len(panel_aspect_ratios_list) == 1
                    single_chord = len(num_chordwise_panels_list) == 1

                    # Consider each convergence parameter to have "passed" if it is
                    # converged, single, or saturated.
                    wake_passed = wake_converged or single_wake or wake_saturated
                    length_passed = length_converged or single_length
                    ar_passed = ar_converged or single_ar or ar_saturated
                    chord_passed = chord_converged or single_chord

                    # If all four convergence parameters have passed, then a
                    # converged or semi-converged combination of parameters has been
                    # found and will be returned.
                    if wake_passed and length_passed and ar_passed and chord_passed:
                        converged_wake_id = _converged_parameter_id(
                            wake_id, single_wake, wake_converged
                        )
                        converged_length_id = _converged_parameter_id(
                            length_id, single_length, length_converged
                        )
                        converged_ar_id = _converged_parameter_id(
                            ar_id, single_ar, ar_converged
                        )
                        converged_chord_id = _converged_parameter_id(
                            chord_id, single_chord, chord_converged
                        )

                        converged_wake = wake_list[converged_wake_id]
                        converged_wake_length = wake_lengths_list[converged_length_id]
                        converged_chordwise_panels = num_chordwise_panels_list[
                            converged_chord_id
                        ]
                        converged_aspect_ratio = panel_aspect_ratios_list[
                            converged_ar_id
                        ]
                        converged_iter_time = float(
                            iter_times[
                                converged_wake_id,
                                converged_length_id,
                                converged_ar_id,
                                converged_chord_id,
                            ]
                        )

                        if single_wake or single_length or single_ar or single_chord:
                            _logger.info("The analysis found a semi-converged case:")
                            if single_wake:
                                _logger.warning("Wake type convergence not checked")
                            if single_length:
                                _logger.warning("Wake length convergence not checked")
                            if single_ar:
                                _logger.warning(
                                    "Panel aspect ratio convergence not checked"
                                )
                            if single_chord:
                                _logger.warning(
                                    "Chordwise Panels convergence not checked"
                                )
                        else:
                            _logger.info("The analysis found a converged case:")

                        if converged_wake:
                            _logger.info("\tWake type: prescribed")
                        else:
                            _logger.info("\tWake type: free")

                        if static:
                            _logger.info(
                                "\tChord lengths: " + str(converged_wake_length)
                            )
                        else:
                            _logger.info("\tCycles: " + str(converged_wake_length))

                        _logger.info(
                            "\tPanel aspect ratio: " + str(converged_aspect_ratio)
                        )
                        _logger.info(
                            "\tChordwise Panels: " + str(converged_chordwise_panels)
                        )
                        _logger.info(
                            "\tSimulation completed in "
                            + str(round(converged_iter_time, 3))
                            + " s"
                        )
                        _logger.info("\tSpanwise Panels:")
                        for airplane_movement_id, airplane_movement in enumerate(
                            ref_airplane_movements
                        ):
                            base_airplane = airplane_movement.base_airplane
                            _logger.info("\t\t" + base_airplane.name + ":")
                            for wing_movement_id, wing_movement in enumerate(
                                airplane_movement.wing_movements
                            ):
                                base_wing = wing_movement.base_wing
                                _logger.info("\t\t\t" + base_wing.name + ":")

                                # An edge-defined Wing is refined by its number of
                                # WingCrossSections rather than its number of spanwise
                                # Panels, so report that instead.
                                if base_wing.spanwise_mesh == "edge_defined":
                                    num_wing_cross_sections_key = (
                                        converged_ar_id,
                                        converged_chord_id,
                                        airplane_movement_id,
                                        wing_movement_id,
                                    )
                                    _logger.info(
                                        "\t\t\t\tWingCrossSections: "
                                        + str(
                                            num_wing_cross_sections_cache[
                                                num_wing_cross_sections_key
                                            ]
                                        )
                                    )
                                    continue

                                for (
                                    wing_cross_section_movement_id,
                                    wing_cross_section_movement,
                                ) in enumerate(
                                    wing_movement.wing_cross_section_movements
                                ):
                                    if (
                                        wing_cross_section_movement_id
                                        < len(
                                            wing_movement.wing_cross_section_movements
                                        )
                                        - 1
                                    ):
                                        # Not the last WingCrossSection, retrieve
                                        # from cache.
                                        num_spanwise_panels_key = (
                                            converged_ar_id,
                                            converged_chord_id,
                                            airplane_movement_id,
                                            wing_movement_id,
                                            wing_cross_section_movement_id,
                                        )
                                        num_spanwise_panels = num_spanwise_panels_cache[
                                            num_spanwise_panels_key
                                        ]
                                    else:
                                        # Last WingCrossSection.
                                        num_spanwise_panels = None
                                    _logger.info(
                                        "\t\t\t\tWingCrossSection "
                                        + str(wing_cross_section_movement_id + 1)
                                        + ": "
                                        + str(num_spanwise_panels)
                                    )

                        # If requested, recreate and run the solver at the converged
                        # parameters so it can be returned. The converged parameters are
                        # frequently coarser than this iteration's, so the converged
                        # solver is rebuilt rather than reusing this iteration's finer
                        # solver.
                        converged_solver: (
                            unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
                            | None
                        ) = None
                        if resolve_converged_solver:
                            _logger.info(
                                "Recreating and running the converged solver..."
                            )
                            converged_problem = _build_unsteady_problem(
                                converged_ar_id,
                                converged_chord_id,
                                converged_aspect_ratio,
                                converged_chordwise_panels,
                                converged_wake_length,
                                ref_problem,
                                num_spanwise_panels_cache,
                                num_wing_cross_sections_cache,
                                delta_time_cache,
                            )
                            converged_solver = unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                                unsteady_problem=converged_problem
                            )
                            converged_solver.run(
                                prescribed_wake=converged_wake,
                                calculate_streamlines=True,
                                show_progress=show_solver_progress,
                            )

                        _logger.info(
                            "Convergence analysis completed in "
                            + str(round(time.time() - run_start_time, 3))
                            + " s"
                        )

                        return (
                            converged_wake,
                            converged_wake_length,
                            converged_aspect_ratio,
                            converged_chordwise_panels,
                            converged_solver,
                        )

    # If all iterations have been checked and none of them resulted in all
    # convergence parameters passing, then indicate that no converged solution was
    # found and return values of None for the converged parameters.
    _logger.info("The analysis did not find a converged case within the given bounds")
    _logger.info(
        "Convergence analysis completed in "
        + str(round(time.time() - run_start_time, 3))
        + " s"
    )
    return None, None, None, None, None


def _build_steady_problem(
    ar_id: int,
    chord_id: int,
    panel_aspect_ratio: int,
    num_chordwise_panels: int,
    ref_problem: problems.SteadyProblem,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
) -> problems.SteadyProblem:
    """Builds the SteadyProblem for one convergence iteration.

    Each of the reference SteadyProblem's Airplanes is copied with the given Panel
    aspect ratio and number of chordwise Panels. Each Wing is refined by its spanwise
    mesh: a trapezoidal Wing has each non-tip WingCrossSection's number of spanwise
    Panels resolved (and cached) to achieve the desired Panel aspect ratio, while an
    edge-defined Wing has its stored edge curves resampled into the number of
    WingCrossSections (also resolved and cached) that achieves it. The copied Airplanes
    are wrapped in a new SteadyProblem with the reference OperatingPoint.

    :param ar_id: The index of the current Panel aspect ratio within the list of Panel
        aspect ratios being tested.
    :param chord_id: The index of the current number of chordwise Panels within the list
        of numbers of chordwise Panels being tested.
    :param panel_aspect_ratio: The Panel aspect ratio to target for this iteration.
    :param num_chordwise_panels: The number of chordwise Panels to use for this
        iteration.
    :param ref_problem: The reference SteadyProblem whose Airplanes and OperatingPoint
        are copied.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels, used for
        trapezoidal Wings. It is read and updated in place.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections, used for edge-defined Wings.
        It is read and updated in place.
    :return: The SteadyProblem for this iteration.
    """
    ref_airplanes = ref_problem.airplanes
    ref_operating_point = ref_problem.operating_point

    # Initialize an empty list to hold this iteration's Airplanes. Then, fill the list by
    # making new copies of each of the Airplanes with modified values for Panel aspect
    # ratio and number of chordwise Panels.
    these_airplanes = []
    for ref_airplane_id, ref_airplane in enumerate(ref_airplanes):
        ref_wings = ref_airplane.wings
        these_wings = []

        for ref_wing_id, ref_wing in enumerate(ref_wings):

            # An edge-defined Wing is refined by resampling its stored edge curves into
            # the number of WingCrossSections that achieves the desired Panel aspect
            # ratio, rather than by sweeping the number of spanwise Panels.
            if ref_wing.spanwise_mesh == "edge_defined":
                this_num_wing_cross_sections = _resolve_num_wing_cross_sections(
                    ar_id,
                    chord_id,
                    ref_airplane_id,
                    ref_wing_id,
                    ref_airplane.name,
                    ref_wing.name,
                    "\t\t\t\t",
                    num_wing_cross_sections_cache,
                    lambda start: _get_num_wing_cross_sections_for_panel_ar(
                        panel_aspect_ratio,
                        num_chordwise_panels,
                        ref_wing,
                        start,
                    ),
                )
                these_wings.append(
                    _build_edge_defined_wing(
                        ref_wing, num_chordwise_panels, this_num_wing_cross_sections
                    )
                )
                continue

            ref_wing_cross_sections = ref_wing.wing_cross_sections
            these_wing_cross_sections = []

            for ref_wing_cross_section_id, ref_wing_cross_section in enumerate(
                ref_wing_cross_sections
            ):

                # If this is not the last WingCrossSection, find the number of spanwise
                # Panels to use for this section of the Wing, based on the desired Panel
                # aspect ratio and number of chordwise Panels.
                if ref_wing_cross_section_id < (len(ref_wing_cross_sections) - 1):
                    this_num_spanwise_panels = _resolve_num_spanwise_panels(
                        ar_id,
                        chord_id,
                        ref_airplane_id,
                        ref_wing_id,
                        ref_wing_cross_section_id,
                        ref_airplane.name,
                        ref_wing.name,
                        "\t\t\t\t",
                        num_spanwise_panels_cache,
                        lambda start: _get_wing_section_num_spanwise_panels(
                            panel_aspect_ratio,
                            num_chordwise_panels,
                            ref_wing.chordwise_spacing,
                            ref_wing_cross_section,
                            ref_wing_cross_sections[ref_wing_cross_section_id + 1],
                            start,
                        ),
                    )
                else:
                    this_num_spanwise_panels = None

                these_wing_cross_sections.append(
                    geometry.wing_cross_section.WingCrossSection(
                        # These values are copied from the reference WingCrossSection.
                        chord=ref_wing_cross_section.chord,
                        Lp_Wcsp_Lpp=ref_wing_cross_section.Lp_Wcsp_Lpp,
                        angles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                        control_surface_symmetry_type=ref_wing_cross_section.control_surface_symmetry_type,
                        control_surface_hinge_point=ref_wing_cross_section.control_surface_hinge_point,
                        control_surface_deflection=ref_wing_cross_section.control_surface_deflection,
                        spanwise_spacing=ref_wing_cross_section.spanwise_spacing,
                        # These values change.
                        num_spanwise_panels=this_num_spanwise_panels,
                        airfoil=geometry.airfoil.Airfoil(
                            name=ref_wing_cross_section.airfoil.name,
                            outline_A_lp=ref_wing_cross_section.airfoil.outline_A_lp,
                            resample=ref_wing_cross_section.airfoil.resample,
                            n_points_per_side=ref_wing_cross_section.airfoil.n_points_per_side,
                        ),
                    )
                )

            these_wings.append(
                geometry.wing.Wing(
                    # These values are copied from the reference Wing.
                    name=ref_wing.name,
                    Ler_Gs_Cgs=ref_wing.Ler_Gs_Cgs,
                    angles_Gs_to_Wn_ixyz=ref_wing.angles_Gs_to_Wn_ixyz,
                    symmetric=ref_wing.symmetric,
                    mirror_only=ref_wing.mirror_only,
                    symmetryNormal_G=ref_wing.symmetryNormal_G,
                    symmetryPoint_G_Cg=ref_wing.symmetryPoint_G_Cg,
                    chordwise_spacing=ref_wing.chordwise_spacing,
                    # These values change.
                    wing_cross_sections=these_wing_cross_sections,
                    num_chordwise_panels=num_chordwise_panels,
                )
            )

        these_airplanes.append(
            geometry.airplane.Airplane(
                # These values are copied from the reference Airplane.
                name=ref_airplane.name,
                Cg_GP1_CgP1=ref_airplane.Cg_GP1_CgP1,
                weight=ref_airplane.weight,
                # These values change.
                wings=these_wings,
                s_ref=None,
                c_ref=None,
                b_ref=None,
            )
        )

    # Create a new SteadyProblem for this iteration.
    return problems.SteadyProblem(
        airplanes=these_airplanes, operating_point=ref_operating_point
    )


def _build_unsteady_problem(
    ar_id: int,
    chord_id: int,
    panel_aspect_ratio: int,
    num_chordwise_panels: int,
    wake_length: int,
    ref_problem: problems.UnsteadyProblem,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
    delta_time_cache: dict[tuple[int, int], float],
) -> problems.UnsteadyProblem:
    """Builds the UnsteadyProblem for one convergence iteration.

    Each of the reference Movement's AirplaneMovements (and its nested WingMovements and
    WingCrossSectionMovements) is copied with the given Panel aspect ratio and number of
    chordwise Panels, preserving every motion parameter. Each Wing is refined by its
    spanwise mesh: a trapezoidal Wing has each non-tip WingCrossSection's number of
    spanwise Panels resolved (and cached) to achieve the desired Panel aspect ratio,
    while an edge-defined Wing has its stored edge curves resampled into the number of
    WingCrossSections (also resolved and cached) that achieves it. An edge-defined
    Wing's WingCrossSectionMovements must all be static, because resampling changes the
    WingCrossSection count and so cannot preserve per-WingCrossSection motion; the
    refined WingCrossSections are wrapped in motion free WingCrossSectionMovements and
    the WingMovement's own motion parameters are copied. The copied AirplaneMovements
    are wrapped in a new Movement (with the given wake length) and UnsteadyProblem.

    The Movement's time step is optimized for this iteration's mesh with Movement's
    iterative delta_time optimizer. The optimum depends only on the mesh, so it is
    resolved once per (Panel aspect ratio index, number of chordwise Panels index) and
    cached, then reused across every wake-state and wake-length combination at that
    mesh.

    :param ar_id: The index of the current Panel aspect ratio within the list of Panel
        aspect ratios being tested.
    :param chord_id: The index of the current number of chordwise Panels within the list
        of numbers of chordwise Panels being tested.
    :param panel_aspect_ratio: The Panel aspect ratio to target for this iteration.
    :param num_chordwise_panels: The number of chordwise Panels to use for this
        iteration.
    :param wake_length: The wake length (number of chords if static, number of cycles if
        variable) to use for this iteration.
    :param ref_problem: The reference UnsteadyProblem whose Movement is copied.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels, used for
        trapezoidal Wings. It is read and updated in place.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections, used for edge-defined Wings.
        It is read and updated in place.
    :param delta_time_cache: The cache mapping a (Panel aspect ratio index, number of
        chordwise Panels index) tuple to a previously optimized delta_time for that
        mesh. It is read and updated in place.
    :return: The UnsteadyProblem for this iteration.
    """
    ref_movement = ref_problem.movement
    static = ref_movement.static
    ref_airplane_movements = ref_movement.airplane_movements
    ref_operating_point_movement = ref_movement.operating_point_movement

    # Initialize an empty list for this iteration's base AirplaneMovements.
    these_base_airplanes = []

    # Create an empty list for the AirplaneMovement copies.
    these_airplane_movements = []

    # Now, we will begin iterating through this iteration's reference AirplaneMovements,
    # WingMovements, and WingCrossSectionMovements, and creating copies of them. These
    # copies will have identical parameters to their respective reference movements
    # except for the number of spanwise Panels (which is based on the Panel aspect
    # ratio), and the number of chordwise Panels.
    #
    # To do this, we iterate over the AirplaneMovements and perform a several step
    # procedure:
    # 1. Reference the AirplaneMovement's base Airplane.
    # 2. Reference the AirplaneMovement's list of WingMovements.
    # 3. Create an empty list for the WingMovements' base Wing copies.
    # 4. Create an empty list for the WingMovement copies.
    # 5. Iterate over the WingMovements.
    #   5.1. Reference the WingMovement's base Wing.
    #   5.2. Reference the WingMovement's list of WingCrossSectionMovements.
    #   5.3. Create an empty list for the WingCrossSectionMovements' base WingCrossSection
    #        copies.
    #   5.4. Create an empty list for the WingCrossSectionMovement copies.
    #   5.5. Iterate over the WingCrossSectionMovements.
    #     5.5.1. Reference the WingCrossSectionMovement's base WingCrossSection.
    #     5.5.2. Calculate the number of spanwise Panels that corresponds to the desired
    #            combination of Panel aspect ratio and number of chordwise Panels.
    #     5.5.3. Create a copy of the base WingCrossSection.
    #     5.5.4. Create a copy of the WingCrossSectionMovement.
    #     5.5.5. Append the base WingCrossSection copy to the list of base WingCrossSection
    #            copies.
    #     5.5.6. Append the WingCrossSectionMovement copy to the list of
    #            WingCrossSectionMovement copies.
    #   5.6. Create a copy of the base Wing.
    #   5.7. Create a copy of the WingMovement.
    #   5.8. Append the base Wing copy to the list  of base Wing copies.
    #   5.9. Append the WingMovement copy to the list of WingMovement copies.
    # 6. Create a copy of the base Airplane.
    # 7. Create a copy of the AirplaneMovement.
    # 8. Append the base Airplane copy to the list of base Airplane copies.
    # 9. Append the AirplaneMovement copy to the list of AirplaneMovement copies.
    ref_airplane_movement: movements.airplane_movement.AirplaneMovement
    for ref_airplane_movement_id, ref_airplane_movement in enumerate(
        ref_airplane_movements
    ):
        # 1. Reference the AirplaneMovement's base Airplane.
        ref_base_airplane = ref_airplane_movement.base_airplane

        # 2. Reference the AirplaneMovement's list of WingMovements.
        ref_wing_movements = ref_airplane_movement.wing_movements

        # 3. Create an empty list for the WingMovements' base Wing copies.
        these_base_wings = []

        # 4. Create an empty list for the WingMovement copies.
        these_wing_movements = []

        # 5. Iterate over the WingMovements.
        for ref_wing_movement_id, ref_wing_movement in enumerate(ref_wing_movements):
            # 5.1. Reference the WingMovement's base Wing.
            ref_base_wing = ref_wing_movement.base_wing

            # An edge-defined Wing is refined by resampling its stored edge curves into
            # the number of WingCrossSections that achieves the desired Panel aspect
            # ratio. Resampling changes the WingCrossSection count, so it cannot preserve
            # per-WingCrossSection motion. Every WingCrossSectionMovement must therefore
            # be static; the refined WingCrossSections are wrapped in motion free
            # WingCrossSectionMovements and only the WingMovement's own motion is carried
            # over.
            if ref_base_wing.spanwise_mesh == "edge_defined":
                for (
                    ref_wing_cross_section_movement_id,
                    ref_wing_cross_section_movement,
                ) in enumerate(ref_wing_movement.wing_cross_section_movements):
                    if np.any(ref_wing_cross_section_movement.ampLp_Wcsp_Lpp) or np.any(
                        ref_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz
                    ):
                        raise ValueError(
                            "analyze_unsteady_convergence cannot refine an edge-defined "
                            "Wing whose WingCrossSectionMovements are not all static, "
                            "because resampling the Wing changes its number of "
                            "WingCrossSections. The WingCrossSectionMovement at index "
                            f"{ref_wing_cross_section_movement_id} of the WingMovement "
                            f'for the Wing named "{ref_base_wing.name}" is not static.'
                        )

                this_num_wing_cross_sections = _resolve_num_wing_cross_sections(
                    ar_id,
                    chord_id,
                    ref_airplane_movement_id,
                    ref_wing_movement_id,
                    ref_base_airplane.name,
                    ref_base_wing.name,
                    "\t\t\t\t\t\t",
                    num_wing_cross_sections_cache,
                    lambda start: _get_num_wing_cross_sections_for_panel_ar(
                        panel_aspect_ratio,
                        num_chordwise_panels,
                        ref_base_wing,
                        start,
                    ),
                )

                this_base_wing = _build_edge_defined_wing(
                    ref_base_wing, num_chordwise_panels, this_num_wing_cross_sections
                )

                these_wing_cross_section_movements = [
                    movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=this_base_wing_cross_section
                    )
                    for this_base_wing_cross_section in (
                        this_base_wing.wing_cross_sections
                    )
                ]

                this_wing_movement = movements.wing_movement.WingMovement(
                    # These values are copied from the reference WingMovement.
                    ampLer_Gs_Cgs=ref_wing_movement.ampLer_Gs_Cgs,
                    periodLer_Gs_Cgs=ref_wing_movement.periodLer_Gs_Cgs,
                    spacingLer_Gs_Cgs=ref_wing_movement.spacingLer_Gs_Cgs,
                    phaseLer_Gs_Cgs=ref_wing_movement.phaseLer_Gs_Cgs,
                    ampAngles_Gs_to_Wn_ixyz=ref_wing_movement.ampAngles_Gs_to_Wn_ixyz,
                    periodAngles_Gs_to_Wn_ixyz=(
                        ref_wing_movement.periodAngles_Gs_to_Wn_ixyz
                    ),
                    spacingAngles_Gs_to_Wn_ixyz=(
                        ref_wing_movement.spacingAngles_Gs_to_Wn_ixyz
                    ),
                    phaseAngles_Gs_to_Wn_ixyz=(
                        ref_wing_movement.phaseAngles_Gs_to_Wn_ixyz
                    ),
                    # These values change.
                    base_wing=this_base_wing,
                    wing_cross_section_movements=these_wing_cross_section_movements,
                )

                these_base_wings.append(this_base_wing)
                these_wing_movements.append(this_wing_movement)
                continue

            # 5.2. Reference the WingMovement's list of WingCrossSectionMovements.
            ref_wing_cross_section_movements = (
                ref_wing_movement.wing_cross_section_movements
            )

            # 5.3. Create an empty list for the WingCrossSectionMovements' base
            # WingCrossSection copies.
            these_base_wing_cross_sections = []

            # 5.4. Create an empty list for the WingCrossSectionMovement copies.
            these_wing_cross_section_movements = []

            # 5.5. Iterate over the WingCrossSectionMovements.
            for (
                ref_wing_cross_section_movement_id,
                ref_wing_cross_section_movement,
            ) in enumerate(ref_wing_cross_section_movements):
                # 5.5.1. Reference the WingCrossSectionMovement's base WingCrossSection.
                ref_base_wing_cross_section = (
                    ref_wing_cross_section_movement.base_wing_cross_section
                )

                # 5.5.2. Calculate the number of spanwise Panels that corresponds to the
                # desired combination of Panel aspect ratio and number of chordwise
                # Panels.
                if ref_wing_cross_section_movement_id < (
                    len(ref_wing_cross_section_movements) - 1
                ):
                    this_num_spanwise_panels = _resolve_num_spanwise_panels(
                        ar_id,
                        chord_id,
                        ref_airplane_movement_id,
                        ref_wing_movement_id,
                        ref_wing_cross_section_movement_id,
                        ref_base_airplane.name,
                        ref_base_wing.name,
                        "\t\t\t\t\t\t",
                        num_spanwise_panels_cache,
                        lambda start: _get_wing_section_movement_num_spanwise_panels(
                            panel_aspect_ratio,
                            num_chordwise_panels,
                            ref_base_wing.chordwise_spacing,
                            ref_movement.airplanes[ref_airplane_movement_id],
                            ref_wing_movement_id,
                            ref_wing_cross_section_movement_id,
                            ref_wing_cross_section_movement_id + 1,
                            start,
                            ref_problem.first_averaging_step,
                        ),
                    )
                else:
                    this_num_spanwise_panels = None

                # 5.5.3. Create a copy of the base WingCrossSection.
                this_base_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                    # These values are copied from the reference base WingCrossSection.
                    chord=ref_base_wing_cross_section.chord,
                    Lp_Wcsp_Lpp=ref_base_wing_cross_section.Lp_Wcsp_Lpp,
                    angles_Wcsp_to_Wcs_ixyz=ref_base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                    control_surface_symmetry_type=ref_base_wing_cross_section.control_surface_symmetry_type,
                    control_surface_hinge_point=ref_base_wing_cross_section.control_surface_hinge_point,
                    control_surface_deflection=ref_base_wing_cross_section.control_surface_deflection,
                    spanwise_spacing=ref_base_wing_cross_section.spanwise_spacing,
                    # These values change.
                    num_spanwise_panels=this_num_spanwise_panels,
                    airfoil=geometry.airfoil.Airfoil(
                        name=ref_base_wing_cross_section.airfoil.name,
                        outline_A_lp=ref_base_wing_cross_section.airfoil.outline_A_lp,
                        resample=ref_base_wing_cross_section.airfoil.resample,
                        n_points_per_side=ref_base_wing_cross_section.airfoil.n_points_per_side,
                    ),
                )

                # 5.5.4. Create a copy of the WingCrossSectionMovement.
                this_wing_cross_section_movement = movements.wing_cross_section_movement.WingCrossSectionMovement(
                    # These values are copied from the reference WingCrossSectionMovement.
                    ampLp_Wcsp_Lpp=ref_wing_cross_section_movement.ampLp_Wcsp_Lpp,
                    periodLp_Wcsp_Lpp=ref_wing_cross_section_movement.periodLp_Wcsp_Lpp,
                    spacingLp_Wcsp_Lpp=ref_wing_cross_section_movement.spacingLp_Wcsp_Lpp,
                    phaseLp_Wcsp_Lpp=ref_wing_cross_section_movement.phaseLp_Wcsp_Lpp,
                    ampAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz,
                    periodAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.periodAngles_Wcsp_to_Wcs_ixyz,
                    spacingAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.spacingAngles_Wcsp_to_Wcs_ixyz,
                    phaseAngles_Wcsp_to_Wcs_ixyz=ref_wing_cross_section_movement.phaseAngles_Wcsp_to_Wcs_ixyz,
                    # These values change.
                    base_wing_cross_section=this_base_wing_cross_section,
                )

                # 5.5.5. Append the base WingCrossSection copy to the list of base
                # WingCrossSection copies.
                these_base_wing_cross_sections.append(this_base_wing_cross_section)

                # 5.5.6. Append the WingCrossSectionMovement copy to the list of
                # WingCrossSectionMovement copies.
                these_wing_cross_section_movements.append(
                    this_wing_cross_section_movement
                )

            # 5.6. Create a copy of base Wing.
            this_base_wing = geometry.wing.Wing(
                # These values are copied from the reference Wing.
                name=ref_base_wing.name,
                Ler_Gs_Cgs=ref_base_wing.Ler_Gs_Cgs,
                angles_Gs_to_Wn_ixyz=ref_base_wing.angles_Gs_to_Wn_ixyz,
                symmetric=ref_base_wing.symmetric,
                mirror_only=ref_base_wing.mirror_only,
                symmetryNormal_G=ref_base_wing.symmetryNormal_G,
                symmetryPoint_G_Cg=ref_base_wing.symmetryPoint_G_Cg,
                chordwise_spacing=ref_base_wing.chordwise_spacing,
                # These values change.
                wing_cross_sections=these_base_wing_cross_sections,
                num_chordwise_panels=num_chordwise_panels,
            )

            # 5.7. Create a copy of the WingMovement.
            this_wing_movement = movements.wing_movement.WingMovement(
                # These values are copied from the reference WingMovement.
                ampLer_Gs_Cgs=ref_wing_movement.ampLer_Gs_Cgs,
                periodLer_Gs_Cgs=ref_wing_movement.periodLer_Gs_Cgs,
                spacingLer_Gs_Cgs=ref_wing_movement.spacingLer_Gs_Cgs,
                phaseLer_Gs_Cgs=ref_wing_movement.phaseLer_Gs_Cgs,
                ampAngles_Gs_to_Wn_ixyz=ref_wing_movement.ampAngles_Gs_to_Wn_ixyz,
                periodAngles_Gs_to_Wn_ixyz=ref_wing_movement.periodAngles_Gs_to_Wn_ixyz,
                spacingAngles_Gs_to_Wn_ixyz=ref_wing_movement.spacingAngles_Gs_to_Wn_ixyz,
                phaseAngles_Gs_to_Wn_ixyz=ref_wing_movement.phaseAngles_Gs_to_Wn_ixyz,
                # These values change.
                base_wing=this_base_wing,
                wing_cross_section_movements=these_wing_cross_section_movements,
            )

            # 5.8. Append the base Wing copy to the list of base Wing copies.
            these_base_wings.append(this_base_wing)

            # 5.9. Append the WingMovement copy to the list of WingMovement copies.
            these_wing_movements.append(this_wing_movement)

        # 6. Create a copy of the base Airplane.
        this_base_airplane = geometry.airplane.Airplane(
            # These values are copied from the reference Airplane.
            name=ref_base_airplane.name,
            Cg_GP1_CgP1=ref_base_airplane.Cg_GP1_CgP1,
            weight=ref_base_airplane.weight,
            # These values change.
            wings=these_base_wings,
            s_ref=None,
            c_ref=None,
            b_ref=None,
        )

        # 7. Create a copy of the AirplaneMovement.
        this_airplane_movement = movements.airplane_movement.AirplaneMovement(
            # These values are copied from the reference AirplaneMovement.
            ampCg_GP1_CgP1=ref_airplane_movement.ampCg_GP1_CgP1,
            periodCg_GP1_CgP1=ref_airplane_movement.periodCg_GP1_CgP1,
            spacingCg_GP1_CgP1=ref_airplane_movement.spacingCg_GP1_CgP1,
            phaseCg_GP1_CgP1=ref_airplane_movement.phaseCg_GP1_CgP1,
            # These values change.
            base_airplane=this_base_airplane,
            wing_movements=these_wing_movements,
        )

        # 8. Append the base Airplane copy to the list of base Airplane copies.
        these_base_airplanes.append(this_base_airplane)

        # 9. Append the AirplaneMovement copy to the list of AirplaneMovement copies.
        these_airplane_movements.append(this_airplane_movement)

    # Resolve the time step for this mesh. The optimum depends only on the mesh, so on a
    # cache miss the iterative optimizer finds it (delta_time="optimize"), and on a hit
    # the cached value is reused. The Movement resolves the string to a float, which is
    # read back and cached below.
    delta_time_key = (ar_id, chord_id)
    this_delta_time: str | float = delta_time_cache.get(delta_time_key, "optimize")

    # Create a new Movement for this iteration.
    if static:
        this_movement = movements.movement.Movement(
            airplane_movements=these_airplane_movements,
            operating_point_movement=ref_operating_point_movement,
            num_chords=wake_length,
            delta_time=this_delta_time,
        )
    else:
        this_movement = movements.movement.Movement(
            airplane_movements=these_airplane_movements,
            operating_point_movement=ref_operating_point_movement,
            num_cycles=wake_length,
            delta_time=this_delta_time,
        )

    # Cache the optimized delta_time on a miss so it is reused across the other
    # wake-state and wake-length combinations at this mesh.
    if delta_time_key not in delta_time_cache:
        delta_time_cache[delta_time_key] = this_movement.delta_time
        _logger.info(
            "\t\t\t\t\tOptimized delta_time: "
            + str(round(this_movement.delta_time, 6))
            + " s"
        )
    else:
        _logger.info(
            "\t\t\t\t\tCached delta_time: "
            + str(round(this_movement.delta_time, 6))
            + " s"
        )

    # Create a new UnsteadyProblem for this iteration.
    return problems.UnsteadyProblem(
        movement=this_movement,
        only_final_results=True,
    )


def _get_wing_section_movement_num_spanwise_panels(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_airplanes: tuple[geometry.airplane.Airplane, ...],
    ref_wing_id: int,
    ref_root_wing_cross_section_id: int,
    ref_tip_wing_cross_section_id: int,
    start_val: int,
    first_applicable_time_step_id: int,
) -> int:
    """Calculates the number of spanwise Panels to use for the wing section of a
    WingMovement based on a desired average Panel aspect ratio.

    :param desired_average_panel_aspect_ratio: The target average Panel aspect ratio to
        achieve. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes). It must be a positive int.
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param chordwise_spacing: The type of spacing between the chordwise Panels. Can be
        "cosine" or "uniform".
    :param ref_airplanes: A tuple of the Airplanes at each time step.
    :param ref_wing_id: The index of the Wing within each Airplane in ref_airplanes. It
        must be a non negative int.
    :param ref_root_wing_cross_section_id: The index of the root WingCrossSection of the
        wing section within the Wing. It must be a non negative int.
    :param ref_tip_wing_cross_section_id: The index of the tip WingCrossSection of the
        wing section within the Wing. It must be a non negative int and greater than
        ``ref_root_wing_cross_section_id``.
    :param start_val: The initial number of spanwise Panels to start the search from. It
        must be a positive int. Using a higher value can speed up the search if a lower
        bound is already known.
    :param first_applicable_time_step_id: The index within ref_airplanes of the first
        time step to consider. It must be a non negative int. All Airplanes from this
        index onward will be analyzed.
    :return: The maximum number of spanwise Panels needed across all applicable time
        steps to achieve the desired average Panel aspect ratio.
    """
    # Slice the tuple of Airplanes to only the applicable ones. For cases with static
    # geometry, this is just the last time step's Airplane. For cases with variable
    # geometry, this is the last max-period cycle's time steps' Airplanes
    ref_airplanes = ref_airplanes[first_applicable_time_step_id:]

    num_time_steps = len(ref_airplanes)
    these_num_spanwise_panels = np.zeros_like(ref_airplanes, dtype=int)

    for time_step_id, ref_airplane_at_time_step in enumerate(ref_airplanes):
        ref_wing_at_time_step = ref_airplane_at_time_step.wings[ref_wing_id]
        ref_root_wing_cross_section_at_time_step = (
            ref_wing_at_time_step.wing_cross_sections[ref_root_wing_cross_section_id]
        )
        ref_tip_wing_cross_section_at_time_step = (
            ref_wing_at_time_step.wing_cross_sections[ref_tip_wing_cross_section_id]
        )

        _logger.debug(
            f"\t\t\t\t\t\t\tCalculating the number of spanwise Panels for time step "
            f"{time_step_id+1}/{num_time_steps}..."
        )

        num_spanwise_panels_at_step = _get_wing_section_num_spanwise_panels(
            desired_average_panel_aspect_ratio=desired_average_panel_aspect_ratio,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
            ref_root_wing_cross_section=ref_root_wing_cross_section_at_time_step,
            ref_tip_wing_cross_section=ref_tip_wing_cross_section_at_time_step,
            start_val=start_val,
        )

        these_num_spanwise_panels[time_step_id] = num_spanwise_panels_at_step

        _logger.debug(
            f"\t\t\t\t\t\t\tNumber of spanwise Panels: "
            f"{num_spanwise_panels_at_step}"
        )

    return int(max(these_num_spanwise_panels))


def _get_wing_section_num_spanwise_panels(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_root_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    ref_tip_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    start_val: int,
) -> int:
    """Calculates the number of spanwise Panels to use for the wing section of a Wing
    based on a desired average Panel aspect ratio.

    :param desired_average_panel_aspect_ratio: The target average Panel aspect ratio to
        achieve. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes). It must be a positive int.
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param chordwise_spacing: The type of spacing between the chordwise Panels. Can be
        "cosine" or "uniform".
    :param ref_root_wing_cross_section: The root WingCrossSection of the wing section.
    :param ref_tip_wing_cross_section: The tip WingCrossSection of the wing section.
    :param start_val: The initial number of spanwise Panels to start the search from. It
        must be a positive int. Using a higher value can speed up the search if a lower
        bound is already known.
    :return: The number of spanwise Panels that results in an average Panel aspect ratio
        closest to the desired value.
    """
    this_num_spanwise_panels = start_val
    average_panel_aspect_ratios = []

    while True:
        this_average_panel_aspect_ratio = _get_wing_section_average_panel_aspect_ratio(
            num_chordwise_panels,
            chordwise_spacing,
            ref_root_wing_cross_section,
            ref_tip_wing_cross_section,
            num_spanwise_panels=this_num_spanwise_panels,
        )
        average_panel_aspect_ratios.append(this_average_panel_aspect_ratio)

        if this_average_panel_aspect_ratio <= desired_average_panel_aspect_ratio:
            break

        this_num_spanwise_panels += 1

    if len(average_panel_aspect_ratios) < 2:
        return this_num_spanwise_panels

    this_aspect_ratio_difference = abs(
        average_panel_aspect_ratios[-1] - desired_average_panel_aspect_ratio
    )
    last_aspect_ratio_difference = abs(
        average_panel_aspect_ratios[-2] - desired_average_panel_aspect_ratio
    )

    if last_aspect_ratio_difference < this_aspect_ratio_difference:
        this_num_spanwise_panels -= 1
    return this_num_spanwise_panels


def _get_wing_section_average_panel_aspect_ratio(
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_root_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    ref_tip_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
    num_spanwise_panels: int,
) -> float:
    """Calculates the average aspect ratio of Panels in a wing section with a particular
    number of chordwise and spanwise Panels.

    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param chordwise_spacing: The type of spacing between the chordwise Panels. Can be
        "cosine" or "uniform".
    :param ref_root_wing_cross_section: The root WingCrossSection of the wing section.
    :param ref_tip_wing_cross_section: The tip WingCrossSection of the wing section.
    :param num_spanwise_panels: The number of spanwise Panels to use. It must be a
        positive int.
    :return: The average Panel aspect ratio for the wing section with the given Panel
        counts. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes).
    """
    this_airplane = geometry.airplane.Airplane(
        wings=[
            geometry.wing.Wing(
                wing_cross_sections=[
                    geometry.wing_cross_section.WingCrossSection(
                        airfoil=geometry.airfoil.Airfoil(
                            name=ref_root_wing_cross_section.airfoil.name,
                            outline_A_lp=ref_root_wing_cross_section.airfoil.outline_A_lp,
                            resample=ref_root_wing_cross_section.airfoil.resample,
                            n_points_per_side=ref_root_wing_cross_section.airfoil.n_points_per_side,
                        ),
                        num_spanwise_panels=num_spanwise_panels,
                        chord=ref_root_wing_cross_section.chord,
                        spanwise_spacing=ref_root_wing_cross_section.spanwise_spacing,
                    ),
                    geometry.wing_cross_section.WingCrossSection(
                        airfoil=geometry.airfoil.Airfoil(
                            name=ref_tip_wing_cross_section.airfoil.name,
                            outline_A_lp=ref_tip_wing_cross_section.airfoil.outline_A_lp,
                            resample=ref_tip_wing_cross_section.airfoil.resample,
                            n_points_per_side=ref_tip_wing_cross_section.airfoil.n_points_per_side,
                        ),
                        num_spanwise_panels=None,
                        chord=ref_tip_wing_cross_section.chord,
                        Lp_Wcsp_Lpp=ref_tip_wing_cross_section.Lp_Wcsp_Lpp,
                        angles_Wcsp_to_Wcs_ixyz=ref_tip_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                    ),
                ],
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing,
            )
        ]
    )

    _average_panel_aspect_ratio = this_airplane.wings[0].average_panel_aspect_ratio
    assert _average_panel_aspect_ratio is not None

    return _average_panel_aspect_ratio


def _build_edge_defined_wing(
    ref_wing: geometry.wing.Wing,
    num_chordwise_panels: int,
    num_wing_cross_sections: int,
) -> geometry.wing.Wing:
    """Rebuilds an edge-defined Wing from its stored edge curves at a new resolution.

    The reference Wing's stored leading edge and trailing edge curves are resampled into
    the requested number of WingCrossSections with Wing.from_edge_points, reusing every
    other geometric parameter (position, orientation, symmetry, chordwise spacing, and
    tip trim) from the reference Wing. The single Airfoil is read from the reference
    Wing's first WingCrossSection and shared across the rebuilt WingCrossSections, which
    is safe because Airfoils are immutable.

    :param ref_wing: The reference edge-defined Wing whose stored edge curves and
        parameters are reused. Its spanwise_mesh must be "edge_defined".
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param num_wing_cross_sections: The number of WingCrossSections to resample the edge
        curves into. It must be an int of at least 2.
    :return: The rebuilt edge-defined Wing.
    """
    # An edge-defined Wing always stores its edge curves and tip trim fraction, so these
    # are never None here. The assertions narrow the accessors' optional types.
    leadingEdgePoints_Wn_Ler = ref_wing.leadingEdgePoints_Wn_Ler
    trailingEdgePoints_Wn_Ler = ref_wing.trailingEdgePoints_Wn_Ler
    tip_trim_fraction = ref_wing.tip_trim_fraction
    assert leadingEdgePoints_Wn_Ler is not None
    assert trailingEdgePoints_Wn_Ler is not None
    assert tip_trim_fraction is not None

    return geometry.wing.Wing.from_edge_points(
        leadingEdgePoints_Wn_Ler=leadingEdgePoints_Wn_Ler,
        trailingEdgePoints_Wn_Ler=trailingEdgePoints_Wn_Ler,
        num_wing_cross_sections=num_wing_cross_sections,
        airfoil=ref_wing.wing_cross_sections[0].airfoil,
        name=ref_wing.name,
        Ler_Gs_Cgs=ref_wing.Ler_Gs_Cgs,
        angles_Gs_to_Wn_ixyz=ref_wing.angles_Gs_to_Wn_ixyz,
        symmetric=ref_wing.symmetric,
        mirror_only=ref_wing.mirror_only,
        symmetryNormal_G=ref_wing.symmetryNormal_G,
        symmetryPoint_G_Cg=ref_wing.symmetryPoint_G_Cg,
        num_chordwise_panels=num_chordwise_panels,
        chordwise_spacing=ref_wing.chordwise_spacing,
        tip_trim_fraction=tip_trim_fraction,
    )


def _get_num_wing_cross_sections_for_panel_ar(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    ref_wing: geometry.wing.Wing,
    start_val: int,
) -> int:
    """Calculates the number of WingCrossSections that gives an edge-defined Wing a
    target average Panel aspect ratio.

    An edge-defined Wing is refined by resampling its stored edge curves into more
    WingCrossSections, not by adding spanwise Panels to fixed WingCrossSections, so the
    number of WingCrossSections is the refinement knob. Increasing it monotonically
    lowers the meshed average Panel aspect ratio.

    Rather than estimate the count from a closed form, this searches for it by meshing.
    A closed form built from the mean chord would undershoot the count on a tapered
    planform, because average_panel_aspect_ratio is the mean of each Panel's individual
    aspect ratio and is dominated by the small-chord tip Panels. Instead, the reference
    Wing is rebuilt at trial WingCrossSection counts (with _build_edge_defined_wing),
    each trial is meshed, and its average Panel aspect ratio is read. The search
    proportionally jumps to bracket the target, then bisects to the crossing, and
    returns whichever of the two counts straddling the target gives the closer average
    Panel aspect ratio. This matches how the trapezoidal
    _get_wing_section_num_spanwise_panels searches, so a target Panel aspect ratio means
    the same physical thing for an edge-defined Wing as for a trapezoidal one. Only
    meshing is performed here, never solving.

    :param desired_average_panel_aspect_ratio: The target average Panel aspect ratio to
        achieve. The Panel aspect ratio is the Panels' average y component length (in
        wing cross section parent axes) divided by their average x component width (in
        wing cross section parent axes). It must be a positive int.
    :param num_chordwise_panels: The number of chordwise Panels to use. It must be a
        positive int.
    :param ref_wing: The reference edge-defined Wing whose stored edge curves are
        resampled. Its spanwise_mesh must be "edge_defined".
    :param start_val: The initial number of WingCrossSections to start the search from.
        It must be a positive int and a lower bound on the result. Using a higher value
        can speed up the search if a lower bound is already known.
    :return: The number of WingCrossSections that gives an average Panel aspect ratio
        closest to the desired value.
    """
    meshed_average_panel_aspect_ratios: dict[int, float] = {}

    def average_panel_aspect_ratio_at(num_wing_cross_sections: int) -> float:
        if num_wing_cross_sections not in meshed_average_panel_aspect_ratios:
            refined_wing = _build_edge_defined_wing(
                ref_wing, num_chordwise_panels, num_wing_cross_sections
            )
            refined_airplane = geometry.airplane.Airplane(wings=[refined_wing])
            this_average_panel_aspect_ratio = refined_airplane.wings[
                0
            ].average_panel_aspect_ratio
            assert this_average_panel_aspect_ratio is not None
            meshed_average_panel_aspect_ratios[num_wing_cross_sections] = (
                this_average_panel_aspect_ratio
            )
        return meshed_average_panel_aspect_ratios[num_wing_cross_sections]

    # from_edge_points requires at least two WingCrossSections.
    lower_num = max(int(start_val), 2)
    lower_average_panel_aspect_ratio = average_panel_aspect_ratio_at(lower_num)

    # If the starting count already meets the target, it is the smallest count that does,
    # so return it.
    if lower_average_panel_aspect_ratio <= desired_average_panel_aspect_ratio:
        return lower_num

    # Proportionally jump to bracket the target from above. The average Panel aspect ratio
    # scales roughly as 1 / (num_wing_cross_sections - 1), so this estimate lands near the
    # crossing in a few steps. Each jump strictly increases the count, and the aspect
    # ratio falls toward zero as the count grows, so the loop terminates.
    upper_num = lower_num
    upper_average_panel_aspect_ratio = lower_average_panel_aspect_ratio
    while upper_average_panel_aspect_ratio > desired_average_panel_aspect_ratio:
        projected_num = 1 + round(
            upper_average_panel_aspect_ratio
            * (upper_num - 1)
            / desired_average_panel_aspect_ratio
        )
        upper_num = max(projected_num, upper_num + 1)
        upper_average_panel_aspect_ratio = average_panel_aspect_ratio_at(upper_num)

    # Bisect between the count above the target (lower_num) and the count at or below it
    # (upper_num) for the integer crossing.
    while upper_num - lower_num > 1:
        middle_num = (lower_num + upper_num) // 2
        if (
            average_panel_aspect_ratio_at(middle_num)
            > desired_average_panel_aspect_ratio
        ):
            lower_num = middle_num
        else:
            upper_num = middle_num

    # Return whichever of the two straddling counts gives the closer average Panel aspect
    # ratio. On a tie, prefer the finer mesh (upper_num), matching the trapezoidal search.
    lower_difference = abs(
        average_panel_aspect_ratio_at(lower_num) - desired_average_panel_aspect_ratio
    )
    upper_difference = abs(
        average_panel_aspect_ratio_at(upper_num) - desired_average_panel_aspect_ratio
    )
    if lower_difference < upper_difference:
        return lower_num
    return upper_num


def _resolve_num_wing_cross_sections(
    panel_aspect_ratio_id: int,
    num_chordwise_panels_id: int,
    airplane_id: int,
    wing_id: int,
    airplane_name: str,
    wing_name: str,
    log_indent: str,
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
    compute_num_wing_cross_sections: Callable[[int], int],
) -> int:
    """Resolves the number of WingCrossSections for one edge-defined Wing, using and
    updating the shared cache.

    The result for a given (Panel aspect ratio, number of chordwise Panels, Airplane,
    Wing) combination is returned from the cache if present. Otherwise, the search
    starts from a conservative lower bound (the smallest number of WingCrossSections
    already found for this Wing at an incrementally coarser mesh, since the current
    finer mesh must need at least that many), ``compute_num_wing_cross_sections`` is
    called with that starting value to find the count, and the result is cached.

    :param panel_aspect_ratio_id: The index of the current Panel aspect ratio within the
        list of Panel aspect ratios being tested.
    :param num_chordwise_panels_id: The index of the current number of chordwise Panels
        within the list of numbers of chordwise Panels being tested.
    :param airplane_id: The index of the current Airplane.
    :param wing_id: The index of the current Wing within the Airplane.
    :param airplane_name: The name of the current Airplane, used in the log messages.
    :param wing_name: The name of the current Wing, used in the log messages.
    :param log_indent: The leading whitespace prepended to the log messages so they nest
        under the calling function's other log output.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a
        previously calculated number of WingCrossSections. It is read and updated in
        place.
    :param compute_num_wing_cross_sections: A callable that takes a starting number of
        WingCrossSections and returns the number of WingCrossSections that achieves the
        desired Panel aspect ratio. It is called only on a cache miss.
    :return: The number of WingCrossSections for the Wing.
    """
    num_wing_cross_sections_key = (
        panel_aspect_ratio_id,
        num_chordwise_panels_id,
        airplane_id,
        wing_id,
    )

    if num_wing_cross_sections_key in num_wing_cross_sections_cache:
        _logger.debug(
            f"{log_indent}Getting the cached number of WingCrossSections calculated for "
            f"{airplane_name}'s {wing_name}..."
        )
        this_num_wing_cross_sections = num_wing_cross_sections_cache[
            num_wing_cross_sections_key
        ]
    else:
        # Start the search from a conservative lower bound: the smallest number of
        # WingCrossSections already found for this Wing at an incrementally coarser mesh
        # (in Panel aspect ratio, number of chordwise Panels, or both), since the current
        # finer mesh must use at least that many. The floor is two, the fewest
        # WingCrossSections an edge-defined Wing can have.
        starting_num_wing_cross_sections = 2
        last_ar_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id,
            airplane_id,
            wing_id,
        )
        last_chord_key = (
            panel_aspect_ratio_id,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
        )
        last_ar_and_chord_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
        )
        last_cache_val = min(
            num_wing_cross_sections_cache.get(last_ar_key, np.inf),
            num_wing_cross_sections_cache.get(last_chord_key, np.inf),
            num_wing_cross_sections_cache.get(last_ar_and_chord_key, np.inf),
        )
        if last_cache_val != np.inf:
            starting_num_wing_cross_sections = int(last_cache_val)

        _logger.debug(
            f"{log_indent}Calculating the number of WingCrossSections for "
            f"{airplane_name}'s {wing_name}, with a starting value of "
            f"{starting_num_wing_cross_sections}..."
        )

        this_num_wing_cross_sections = compute_num_wing_cross_sections(
            starting_num_wing_cross_sections
        )

        num_wing_cross_sections_cache[num_wing_cross_sections_key] = (
            this_num_wing_cross_sections
        )

    _logger.debug(
        f"{log_indent}Number of WingCrossSections: {this_num_wing_cross_sections}"
    )

    return this_num_wing_cross_sections


def _validate_panel_aspect_ratio_bounds(
    panel_aspect_ratio_bounds: tuple[int, int],
) -> None:
    """Validates the panel_aspect_ratio_bounds parameter shared by the convergence
    analysis functions.

    :param panel_aspect_ratio_bounds: A tuple of two ints, in descending order, giving
        the range of Panel aspect ratios to consider, from largest to smallest. The
        second value must be positive.
    :return: None
    """
    if not (
        isinstance(panel_aspect_ratio_bounds, tuple)
        and len(panel_aspect_ratio_bounds) == 2
    ):
        raise TypeError("panel_aspect_ratio_bounds must be a tuple with length 2.")
    if not all(isinstance(bound, int) for bound in panel_aspect_ratio_bounds):
        raise TypeError("Both values in panel_aspect_ratio_bounds must be ints.")
    if panel_aspect_ratio_bounds[0] < panel_aspect_ratio_bounds[1]:
        raise ValueError(
            "The first value in panel_aspect_ratio_bounds must be greater than or "
            "equal to the second value."
        )
    if panel_aspect_ratio_bounds[1] <= 0:
        raise ValueError("Both values in panel_aspect_ratio_bounds must be positive.")


def _validate_num_chordwise_panels_bounds(
    num_chordwise_panels_bounds: tuple[int, int],
) -> None:
    """Validates the num_chordwise_panels_bounds parameter shared by the convergence
    analysis functions.

    :param num_chordwise_panels_bounds: A tuple of two ints, in ascending order, giving
        the range of numbers of chordwise Panels to consider. The first value must be
        positive.
    :return: None
    """
    if not (
        isinstance(num_chordwise_panels_bounds, tuple)
        and len(num_chordwise_panels_bounds) == 2
    ):
        raise TypeError("num_chordwise_panels_bounds must be a tuple with length 2.")
    if not all(isinstance(bound, int) for bound in num_chordwise_panels_bounds):
        raise TypeError("Both values in num_chordwise_panels_bounds must be ints.")
    if num_chordwise_panels_bounds[1] < num_chordwise_panels_bounds[0]:
        raise ValueError(
            "The first value in num_chordwise_panels_bounds must be less than or "
            "equal to the second value."
        )
    if num_chordwise_panels_bounds[0] <= 0:
        raise ValueError("Both values in num_chordwise_panels_bounds must be positive.")


def _reject_unrefinable_wings(
    ref_airplanes: tuple[geometry.airplane.Airplane, ...],
    analyze_function_name: str,
) -> None:
    """Raises if any Wing in the reference Airplanes cannot be refined.

    The convergence analysis functions refine a Wing by its spanwise mesh: a trapezoidal
    Wing by sweeping the number of spanwise Panels while holding its WingCrossSections
    fixed, and an edge-defined Wing (one built with Wing.from_edge_points) by resampling
    its stored edge curves into more WingCrossSections. Both are supported. An exploded
    Wing (one built with explode_into_strips=True) is already in single-panel-strip form
    but carries no edge curves, so it can be neither resampled nor meaningfully swept
    for Panel density, and it is rejected.

    :param ref_airplanes: A tuple of the reference Airplanes whose Wings are checked.
    :param analyze_function_name: The name of the calling analysis function, used in the
        error message.
    :return: None
    """
    for ref_airplane in ref_airplanes:
        for ref_wing in ref_airplane.wings:
            if ref_wing.spanwise_mesh not in ("trapezoidal", "edge_defined"):
                raise ValueError(
                    f"{analyze_function_name} cannot refine the Wing named "
                    f'"{ref_wing.name}", whose spanwise mesh is '
                    f'"{ref_wing.spanwise_mesh}". Convergence analysis can refine an '
                    "edge-defined Wing (one built from edge curves) or a non-edge-defined "
                    "Wing built from WingCrossSections, but not a Wing already exploded "
                    "into single-panel strips."
                )


def _validate_coefficient_mask(
    coefficient_mask: tuple[bool, bool, bool, bool, bool, bool] | None,
) -> np.ndarray:
    """Validates the coefficient_mask parameter shared by the convergence analysis
    functions and returns it as an ndarray.

    :param coefficient_mask: A tuple of six bools that determines which of the six load
        coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ) must converge, or None to require
        all six. At least one element must be True.
    :return: A (6,) ndarray of bools representing the validated mask.
    """
    if coefficient_mask is None:
        return np.ones(6, dtype=bool)
    if not isinstance(coefficient_mask, tuple):
        raise TypeError("coefficient_mask must be a tuple or None.")
    if len(coefficient_mask) != 6:
        raise ValueError("coefficient_mask must have exactly six elements.")
    if not all(isinstance(element, bool) for element in coefficient_mask):
        raise TypeError("Every element in coefficient_mask must be a bool.")
    if not any(coefficient_mask):
        raise ValueError("At least one element in coefficient_mask must be True.")
    return np.array(coefficient_mask, dtype=bool)


def _check_coefficient_convergence(
    these_coefficients: np.ndarray,
    coarser_coefficients: np.ndarray,
    rtol: float,
    atol: float,
    coefficient_mask: np.ndarray,
) -> tuple[bool, float, int]:
    """Checks per-coefficient convergence between an iteration's load coefficients and
    those of an incrementally coarser iteration.

    For each of the six load coefficients of each Airplane, the error is the absolute
    difference between this iteration's coefficient and the coarser iteration's, and the
    tolerance is atol + rtol * max(abs(this), abs(coarser)). A coefficient is converged
    when its error is less than or equal to its tolerance. The iteration is converged
    only when every unmasked coefficient of every Airplane has converged, so a multi-
    Airplane study converges only once all its Airplanes have.

    A convergence metric accompanies the result for logging. It is a percentage that is
    100.0 when a coefficient has converged and falls toward 0.0 as its error grows past
    its tolerance, so a lower metric means a coefficient is further from converging. The
    minimum metric across the unmasked coefficients of all Airplanes and the index of
    the coefficient that attains it are returned to identify the limiting coefficient.

    :param these_coefficients: A (M,6) ndarray of floats, where M is the number of
        Airplanes, of this iteration's six load coefficients (cFX, cFY, cFZ, cMX, cMY,
        cMZ) for each Airplane.
    :param coarser_coefficients: A (M,6) ndarray of floats, of the same shape, of the
        incrementally coarser iteration's six load coefficients for each Airplane.
    :param rtol: The relative tolerance. It must be a positive float.
    :param atol: The absolute tolerance floor for coefficients near zero. It must be a
        positive float.
    :param coefficient_mask: A (6,) ndarray of bools that determines which of the six
        load coefficients must converge. At least one element is True.
    :return: A tuple of a bool, a float, and an int. In order, whether every unmasked
        coefficient of every Airplane has converged, the minimum convergence metric
        (percentage) across the unmasked coefficients of all Airplanes, and the index
        within the six coefficients of the limiting coefficient (the one attaining that
        minimum metric).
    """
    errors = np.abs(these_coefficients - coarser_coefficients)
    tolerances = atol + rtol * np.maximum(
        np.abs(these_coefficients), np.abs(coarser_coefficients)
    )
    converged = errors <= tolerances

    all_converged = bool(np.all(converged[:, coefficient_mask]))

    # The metric is 100.0 for a converged coefficient (including one with zero error) and
    # falls toward 0.0 as the error grows past its tolerance.
    with np.errstate(divide="ignore", invalid="ignore"):
        metrics = np.where(
            errors == 0.0,
            100.0,
            100.0 * np.minimum(1.0, tolerances / errors),
        )

    # Exclude the masked-out coefficients from the limiting-coefficient search by setting
    # their metrics to infinity. The mask has shape (6,) and broadcasts across Airplanes.
    masked_metrics = np.where(coefficient_mask, metrics, np.inf)
    min_metric = float(np.min(masked_metrics))
    limiting_coefficient_id = int(np.argmin(masked_metrics) % 6)

    return all_converged, min_metric, limiting_coefficient_id


def _resolve_num_spanwise_panels(
    panel_aspect_ratio_id: int,
    num_chordwise_panels_id: int,
    airplane_id: int,
    wing_id: int,
    wing_cross_section_id: int,
    airplane_name: str,
    wing_name: str,
    log_indent: str,
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    compute_num_spanwise_panels: Callable[[int], int],
) -> int:
    """Resolves the number of spanwise Panels for one Wing section, using and updating
    the shared cache.

    The result for a given (Panel aspect ratio, number of chordwise Panels, Airplane,
    Wing, WingCrossSection) combination is returned from the cache if present.
    Otherwise, the search starts from a conservative lower bound (the smallest number of
    spanwise Panels already found for this Wing section at an incrementally coarser
    mesh, since the current finer mesh must need at least that many),
    ``compute_num_spanwise_panels`` is called with that starting value to find the
    count, and the result is cached.

    :param panel_aspect_ratio_id: The index of the current Panel aspect ratio within the
        list of Panel aspect ratios being tested.
    :param num_chordwise_panels_id: The index of the current number of chordwise Panels
        within the list of numbers of chordwise Panels being tested.
    :param airplane_id: The index of the current Airplane.
    :param wing_id: The index of the current Wing within the Airplane.
    :param wing_cross_section_id: The index of the current WingCrossSection within the
        Wing.
    :param airplane_name: The name of the current Airplane, used in the log messages.
    :param wing_name: The name of the current Wing, used in the log messages.
    :param log_indent: The leading whitespace prepended to the log messages so they nest
        under the calling function's other log output.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a previously calculated number of spanwise Panels. It is read
        and updated in place.
    :param compute_num_spanwise_panels: A callable that takes a starting number of
        spanwise Panels and returns the number of spanwise Panels that achieves the
        desired Panel aspect ratio. It is called only on a cache miss.
    :return: The number of spanwise Panels for the Wing section.
    """
    num_spanwise_panels_key = (
        panel_aspect_ratio_id,
        num_chordwise_panels_id,
        airplane_id,
        wing_id,
        wing_cross_section_id,
    )

    if num_spanwise_panels_key in num_spanwise_panels_cache:
        _logger.debug(
            f"{log_indent}Getting the cached number of spanwise Panels calculated for "
            f"the #{wing_cross_section_id + 1} WingCrossSection of {airplane_name}'s "
            f"{wing_name}..."
        )
        this_num_spanwise_panels = num_spanwise_panels_cache[num_spanwise_panels_key]
    else:
        # Start the search from a conservative lower bound: the smallest number of
        # spanwise Panels already found for this Wing section at an incrementally
        # coarser mesh (in Panel aspect ratio, number of chordwise Panels, or both),
        # since the current finer mesh must use at least that many.
        starting_num_spanwise_panels = 1
        last_ar_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id,
            airplane_id,
            wing_id,
            wing_cross_section_id,
        )
        last_chord_key = (
            panel_aspect_ratio_id,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
            wing_cross_section_id,
        )
        last_ar_and_chord_key = (
            panel_aspect_ratio_id - 1,
            num_chordwise_panels_id - 1,
            airplane_id,
            wing_id,
            wing_cross_section_id,
        )
        last_cache_val = min(
            num_spanwise_panels_cache.get(last_ar_key, np.inf),
            num_spanwise_panels_cache.get(last_chord_key, np.inf),
            num_spanwise_panels_cache.get(last_ar_and_chord_key, np.inf),
        )
        if last_cache_val != np.inf:
            starting_num_spanwise_panels = int(last_cache_val)

        _logger.debug(
            f"{log_indent}Calculating the number of spanwise Panels for the "
            f"#{wing_cross_section_id + 1} WingCrossSection of {airplane_name}'s "
            f"{wing_name}, with a starting value of {starting_num_spanwise_panels}..."
        )

        this_num_spanwise_panels = compute_num_spanwise_panels(
            starting_num_spanwise_panels
        )

        num_spanwise_panels_cache[num_spanwise_panels_key] = this_num_spanwise_panels

    _logger.debug(f"{log_indent}Number of spanwise Panels: {this_num_spanwise_panels}")

    return this_num_spanwise_panels


def _converged_parameter_id(
    this_id: int,
    single: bool,
    converged: bool,
) -> int:
    """Selects the index of the converged value for one convergence parameter.

    This is only meaningful once the parameter has passed, that is once it is converged,
    single, or saturated. When only a single value of the parameter was tested, this
    iteration's index is returned, since there is no coarser value to compare against.
    Otherwise, when this iteration is converged against the incrementally coarser one,
    the coarser index is returned, because refining from the coarser value to this one
    changed the result by less than the convergence criteria, so the coarser value is
    the converged one. When this iteration passed without converging (the parameter is
    saturated at its finest setting), this iteration's own index is returned.

    :param this_id: The index of this iteration's value within the list of values tested
        for the parameter.
    :param single: Whether only a single value of the parameter was tested.
    :param converged: Whether this iteration is converged against the incrementally
        coarser iteration for the parameter.
    :return: The index of the converged value within the list of values tested.
    """
    if single:
        return this_id
    if converged:
        return this_id - 1
    return this_id
