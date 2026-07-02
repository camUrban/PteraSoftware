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
from pathlib import Path

import numpy as np

from . import (
    _convergence_cache,
    _convergence_meshing,
    _logging,
    _parameter_validation,
    _serialization,
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


def analyze_steady_convergence(
    ref_problem: problems.SteadyProblem,
    solver_type: str,
    panel_aspect_ratio_bounds: tuple[int, int] = (4, 1),
    num_chordwise_panels_bounds: tuple[int, int] = (3, 12),
    rtol: float | int = 0.05,
    atol: float | int = 0.001,
    coefficient_mask: tuple[bool, bool, bool, bool, bool, bool] | None = None,
    resolve_converged_solver: bool | np.bool_ = False,
    cache_path: str | Path | None = None,
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
    :param cache_path: An optional path (a str or Path, which must end with ".json") to
        a JSON file that caches each iteration's solved load coefficients, keyed on the
        reference problem and the mesh parameters, together with the mesh counts it
        resolved (the spanwise Panel and WingCrossSection counts). When given, an
        iteration already in the cache reuses those stored coefficients and counts
        instead of re-running the solver and re-resolving the mesh, and each new
        iteration is written through to the file, so an interrupted or repeated study
        reuses the work it has already done. The mesh counts are keyed on the absolute
        mesh rather than a sweep index, so a later run over different bounds still
        reuses any iteration it shares. When None, no cache is read or written. The
        default is None.
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

    # Normalize the cache_path parameter to a Path, or leave it as None to disable
    # caching. Require a ".json" suffix.
    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.suffix != ".json":
            raise ValueError(
                f"cache_path must end with '.json', got '{cache_path.name}'."
            )
        if cache_path.is_dir():
            raise ValueError(
                f"cache_path must be a file path, got directory '{cache_path}'."
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

    # Compute the reference problem's content hash once. It anchors each iteration's
    # solve-cache and memo keys to this specific geometry, operating point, and motion.
    ref_problem_hash = _serialization.hash_object(ref_problem)

    # Load the JSON solve cache, or start empty when caching is disabled or the file does
    # not yet exist. Solved coefficients are looked up from and written through to it.
    solve_cache = _convergence_cache.load_solve_cache(cache_path)

    # Pre-populate the in-loop memo caches from the cache file's memo section, translating
    # each stored memo from its absolute mesh back onto this run's sweep indices. Memos
    # outside this run's bounds or written for another reference problem are dropped. A
    # cache with a populated memo section lets a warm run reuse these counts instead of
    # re-resolving them. The number of spanwise Panels cache is keyed on a tuple of 5
    # ints (ar_id, chord_id, ref_airplane_id, ref_wing_id, ref_wing_cross_section_id) for
    # trapezoidal Wings, and the analogous number of WingCrossSections cache is keyed on a
    # tuple of 4 ints (ar_id, chord_id, ref_airplane_id, ref_wing_id) for edge-defined
    # Wings. A steady analysis has no delta_time, so that cache is discarded.
    num_spanwise_panels_cache, num_wing_cross_sections_cache, _ = (
        _convergence_cache.memos_from_disk(
            _convergence_cache.load_memo_cache(cache_path),
            ref_problem_hash,
            panel_aspect_ratios_list,
            num_chordwise_panels_list,
        )
    )

    # Persist the whole cache (solved coefficients and memos) to disk on each solve miss.
    # The memo caches are translated back to their absolute mesh keys at write time. This
    # is None when caching is disabled, in which case nothing is written.
    persist_cache: Callable[[], None] | None = None
    if cache_path is not None:

        def persist_cache() -> None:
            assert cache_path is not None
            _convergence_cache.write_cache(
                cache_path,
                solve_cache,
                _convergence_cache.memos_to_disk(
                    ref_problem_hash,
                    panel_aspect_ratios_list,
                    num_chordwise_panels_list,
                    num_spanwise_panels_cache,
                    num_wing_cross_sections_cache,
                    {},
                ),
            )

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
            # and number of chordwise Panels. The mesh is always built so its spanwise
            # Panel and WingCrossSection counts are recorded, even when the solve itself
            # is served from the cache.
            this_problem = _convergence_meshing.build_steady_problem(
                ar_id,
                chord_id,
                panel_aspect_ratio,
                num_chordwise_panels,
                ref_problem,
                num_spanwise_panels_cache,
                num_wing_cross_sections_cache,
            )
            these_airplanes = this_problem.airplanes

            # Build this mesh's solve-cache key from the reference problem, the solver
            # type, and the mesh parameters that determine the solve.
            solve_cache_key = _convergence_cache.solve_cache_key(
                ref_problem_hash, solver_type, panel_aspect_ratio, num_chordwise_panels
            )
            cached = solve_cache_key in solve_cache

            def solve_this_mesh() -> np.ndarray:
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

                # Run the steady solver. Skip the streamline trace since it does not
                # affect convergence metrics.
                this_solver.run(calculate_streamlines=False)

                # Create and fill an ndarray with each of this iteration's Airplanes' six
                # load coefficients (cFX, cFY, cFZ, cMX, cMY, cMZ).
                solved_coefficients = np.zeros((len(these_airplanes), 6), dtype=float)

                for this_airplane_id, this_airplane in enumerate(these_airplanes):
                    _forceCoefficients_W = this_airplane.forceCoefficients_W
                    assert _forceCoefficients_W is not None

                    _momentCoefficients_W_CgP1 = this_airplane.momentCoefficients_W_CgP1
                    assert _momentCoefficients_W_CgP1 is not None

                    solved_coefficients[this_airplane_id, 0:3] = _forceCoefficients_W
                    solved_coefficients[this_airplane_id, 3:6] = (
                        _momentCoefficients_W_CgP1
                    )

                return solved_coefficients

            # Reuse the cached coefficients on a hit, otherwise solve this mesh and store
            # them. Time the whole step: it is the solve time on a miss and negligible on
            # a hit.
            iter_start = time.time()
            theseCoefficients = _convergence_cache.cached_solve(
                solve_cache, solve_cache_key, solve_this_mesh, persist_cache
            )
            this_iter_time = time.time() - iter_start

            # Populate the ndarray that stores information from all the iterations with
            # the data from this iteration.
            coefficients[ar_id, chord_id, :, :] = theseCoefficients
            iter_times[ar_id, chord_id] = this_iter_time

            if cached:
                _logger.info("\t\t\tLoaded cached solution.")
            else:
                _logger.info(
                    "\t\t\tSimulation completed in "
                    + str(round(this_iter_time, 3))
                    + " s"
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
                    + f"{ar_metric:#.4G}"
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
                    + f"{chord_metric:#.4G}"
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
                    converged_problem = _convergence_meshing.build_steady_problem(
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
    cache_path: str | Path | None = None,
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
    :param cache_path: An optional path (a str or Path, which must end with ".json") to
        a JSON file that caches each iteration's solved load coefficients, keyed on the
        reference problem and the wake state, wake length, and mesh parameters, together
        with the mesh values it resolved (the spanwise Panel and WingCrossSection counts
        and the optimized delta_time). When given, an iteration already in the cache
        reuses those stored coefficients and values instead of re-running the solver and
        re-resolving the mesh, so a warm run also skips the delta_time optimizer, and
        each new iteration is written through to the file, so an interrupted or repeated
        study reuses the work it has already done. The mesh values are keyed on the
        absolute mesh rather than a sweep index, so a later run over different bounds
        still reuses any iteration it shares. When None, no cache is read or written.
        The default is None.
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

    # Normalize the cache_path parameter to a Path, or leave it as None to disable
    # caching. Require a ".json" suffix.
    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.suffix != ".json":
            raise ValueError(
                f"cache_path must end with '.json', got '{cache_path.name}'."
            )
        if cache_path.is_dir():
            raise ValueError(
                f"cache_path must be a file path, got directory '{cache_path}'."
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

    # Compute the reference problem's content hash once. It anchors each iteration's
    # solve-cache and memo keys to this specific geometry, operating point, and motion.
    ref_problem_hash = _serialization.hash_object(ref_problem)

    # Load the JSON solve cache, or start empty when caching is disabled or the file does
    # not yet exist. Solved coefficients are looked up from and written through to it.
    solve_cache = _convergence_cache.load_solve_cache(cache_path)

    # Pre-populate the in-loop memo caches from the cache file's memo section, translating
    # each stored memo from its absolute mesh back onto this run's sweep indices. Memos
    # outside this run's bounds or written for another reference problem are dropped. A
    # cache with a populated memo section lets a warm run reuse these values instead of
    # re-resolving them, which for the delta_time cache skips the iterative optimizer. The
    # number of spanwise Panels cache is keyed on a tuple of 5 ints (ar_id, chord_id,
    # ref_base_airplane_id, ref_base_wing_id, ref_base_wing_cross_section_id) for
    # trapezoidal Wings, the number of WingCrossSections cache is keyed on a tuple of 4
    # ints (ar_id, chord_id, ref_base_airplane_id, ref_base_wing_id) for edge-defined
    # Wings, and the delta_time cache is keyed on a tuple of 2 ints (ar_id, chord_id),
    # since the optimum depends only on the mesh and is reused across every wake state and
    # wake length at that mesh.
    (
        num_spanwise_panels_cache,
        num_wing_cross_sections_cache,
        delta_time_cache,
    ) = _convergence_cache.memos_from_disk(
        _convergence_cache.load_memo_cache(cache_path),
        ref_problem_hash,
        panel_aspect_ratios_list,
        num_chordwise_panels_list,
    )

    # Persist the whole cache (solved coefficients and memos) to disk on each solve miss.
    # The memo caches are translated back to their absolute mesh keys at write time. This
    # is None when caching is disabled, in which case nothing is written.
    persist_cache: Callable[[], None] | None = None
    if cache_path is not None:

        def persist_cache() -> None:
            assert cache_path is not None
            _convergence_cache.write_cache(
                cache_path,
                solve_cache,
                _convergence_cache.memos_to_disk(
                    ref_problem_hash,
                    panel_aspect_ratios_list,
                    num_chordwise_panels_list,
                    num_spanwise_panels_cache,
                    num_wing_cross_sections_cache,
                    delta_time_cache,
                ),
            )

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
                    # length, Panel aspect ratio, and number of chordwise Panels. The
                    # mesh is always built so its spanwise Panel and WingCrossSection
                    # counts and its optimized delta_time are recorded, even when the
                    # solve itself is served from the cache.
                    this_problem = _convergence_meshing.build_unsteady_problem(
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

                    # Build this mesh's solve-cache key from the reference problem and
                    # the wake state, wake length, and mesh parameters that determine the
                    # solve.
                    solve_cache_key = _convergence_cache.solve_cache_key(
                        ref_problem_hash,
                        wake,
                        wake_length,
                        panel_aspect_ratio,
                        num_chordwise_panels,
                    )
                    cached = solve_cache_key in solve_cache

                    def solve_this_mesh() -> np.ndarray:
                        # Create and run this iteration's
                        # UnsteadyRingVortexLatticeMethodSolver.
                        this_solver = unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                            unsteady_problem=this_problem
                        )

                        _logger.info("\t\t\t\t\tStarting simulation...")

                        this_solver.run(
                            prescribed_wake=wake,
                            calculate_streamlines=False,
                            show_progress=show_solver_progress,
                        )

                        # Create and fill an ndarray with each of this iteration's
                        # Airplanes' six final load coefficients (cFX, cFY, cFZ, cMX, cMY,
                        # cMZ).
                        solved_coefficients = np.zeros(
                            (len(ref_airplane_movements), 6), dtype=float
                        )

                        for airplane_id in range(len(ref_airplane_movements)):
                            # For static geometry, the final load coefficients are the
                            # final time step's. For variable geometry, they are the mean
                            # over the final cycle (the signed time-average).
                            if static:
                                solved_coefficients[airplane_id, 0:3] = (
                                    this_problem.finalForceCoefficients_W[airplane_id]
                                )
                                solved_coefficients[airplane_id, 3:6] = (
                                    this_problem.finalMomentCoefficients_W_CgP1[
                                        airplane_id
                                    ]
                                )
                            else:
                                solved_coefficients[airplane_id, 0:3] = (
                                    this_problem.finalMeanForceCoefficients_W[
                                        airplane_id
                                    ]
                                )
                                solved_coefficients[airplane_id, 3:6] = (
                                    this_problem.finalMeanMomentCoefficients_W_CgP1[
                                        airplane_id
                                    ]
                                )

                        return solved_coefficients

                    # Reuse the cached coefficients on a hit, otherwise solve this mesh
                    # and store them. Time the whole step: it is the solve time on a miss
                    # and negligible on a hit.
                    iter_start = time.time()
                    theseFinalCoefficients = _convergence_cache.cached_solve(
                        solve_cache, solve_cache_key, solve_this_mesh, persist_cache
                    )
                    this_iter_time = time.time() - iter_start

                    # Populate the ndarray that stores information from all the
                    # iterations with the data from this iteration.
                    finalCoefficients[wake_id, length_id, ar_id, chord_id, :, :] = (
                        theseFinalCoefficients
                    )
                    iter_times[wake_id, length_id, ar_id, chord_id] = this_iter_time

                    if cached:
                        _logger.info("\t\t\t\t\tLoaded cached solution.")
                    else:
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
                            + f"{wake_metric:#.4G}"
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
                            + f"{length_metric:#.4G}"
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
                            + f"{ar_metric:#.4G}"
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
                            + f"{chord_metric:#.4G}"
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
                            converged_problem = (
                                _convergence_meshing.build_unsteady_problem(
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
