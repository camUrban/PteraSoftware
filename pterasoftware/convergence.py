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

convergence_logger = _logging.get_logger("convergence")


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
# TODO: If a converged mesh was found, consider also returning the converged solver.
def analyze_steady_convergence(
    ref_problem: problems.SteadyProblem,
    solver_type: str,
    panel_aspect_ratio_bounds: tuple[int, int] = (4, 1),
    num_chordwise_panels_bounds: tuple[int, int] = (3, 12),
    convergence_criteria: float | int = 5.0,
) -> tuple[int, int] | tuple[None, None]:
    """Finds the converged parameters of a SteadyProblem solved using a given steady
    solver.

    **Procedure:**

    Convergence is found by varying the SteadyProblem's Airplanes' Panels' aspect ratios
    their Wings' numbers of chordwise Panels. These values are iterated over via two
    nested for loops (with the number of chordwise Panels as the inner loop).

    With each new combination of these values, the SteadyProblem is solved, and its
    resultant load coefficients are stored. Each Airplanes' three force coefficients are
    combined by taking their root-sum-square to find the resultant force coefficient.
    Then, the absolute percent change (APE) of each Airplanes' resultant force
    coefficient is found between this iteration, and the iterations with incrementally
    coarser meshes in the two parameter directions (panel aspect ratio and number of
    chordwise panels). These two steps are repeated for the three moment coefficients.

    The maximums of the resultant force coefficient APEs and resultant moment
    coefficient APEs are found. This leaves us with two maximum APEs, one for each
    parameter direction, per Airplane. Next, we take the maximum of each parameter
    directions' APEs across all Airplanes, leaving us with two maximum APEs total. If
    either of the parameter direction APEs is below the convergence criteria, then this
    iteration has found a converged solution for that parameter direction.

    If an iteration's two APEs are both below the converged criteria, then we exit the
    nested for loops and return the converged parameters. However, the converged
    parameters are actually the values incrementally coarser than the final values
    (because the incrementally coarser values were found to be within the convergence
    criteria percent difference from the final values).

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
    :param convergence_criteria: A positive number (int or float) that determines the
        point at which the function considers the simulation to have converged.
        Specifically, it is the maximum absolute percent change in the combined load
        coefficients. Therefore, it is in units of percent. Refer to the description in
        this function's docstring for more details on how it affects the solver. In
        short, set this value to 5.0 for a lenient convergence, and 1.0 for a strict
        convergence. Values are converted to floats internally. The default is 5.0.
    :return: A tuple of two ints or a tuple of two Nones. In order, they are the
        converged of Panel aspect ratio and the converged number of chordwise Panels. If
        the function could not find a set of converged parameters, it returns (None,
        None).
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

    # Validate the num_chordwise_panels_bounds parameter.
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

    # Validate the convergence_criteria parameter.
    convergence_criteria = _parameter_validation.number_in_range_return_float(
        convergence_criteria, "convergence_criteria", min_val=0.0, min_inclusive=False
    )

    convergence_logger.info("Beginning convergence analysis...")

    ref_operating_point = ref_problem.operating_point
    ref_airplanes = ref_problem.airplanes

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
    combinedForceCoefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplanes),
        ),
        dtype=float,
    )
    combinedMomentCoefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplanes),
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

    # Begin iterating through the outer loop of Panel aspect ratios.
    for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):
        convergence_logger.info("\tPanel aspect ratio: " + str(panel_aspect_ratio))

        # Begin iterating through the inner loop of number of chordwise Panels.
        for chord_id, num_chordwise_panels in enumerate(num_chordwise_panels_list):
            convergence_logger.info(
                "\t\tChordwise Panels: " + str(num_chordwise_panels)
            )

            iteration += 1
            convergence_logger.info(
                "\t\t\tIteration number: " + str(iteration) + "/" + str(num_iterations)
            )

            # Initialize an empty list to hold this iteration's Airplanes. Then,
            # fill the list by making new copies of each of the Airplanes with
            # modified values for Panel aspect ratio and number of chordwise Panels.
            these_airplanes = []
            for ref_airplane_id, ref_airplane in enumerate(ref_airplanes):
                ref_wings = ref_airplane.wings
                these_wings = []

                for ref_wing_id, ref_wing in enumerate(ref_wings):
                    ref_wing_cross_sections = ref_wing.wing_cross_sections
                    these_wing_cross_sections = []

                    for (
                        ref_wing_cross_section_id,
                        ref_wing_cross_section,
                    ) in enumerate(ref_wing_cross_sections):

                        # If this is not the last WingCrossSection, find the number
                        # of spanwise Panels to use for this section of the Wing,
                        # based on the desired Panel aspect ratio and number of
                        # chordwise Panels.
                        if ref_wing_cross_section_id < (
                            len(ref_wing_cross_sections) - 1
                        ):
                            # Check if we've already calculated the number of
                            # spanwise Panels for this case/combination of parameters.
                            num_spanwise_panels_key = (
                                ar_id,
                                chord_id,
                                ref_airplane_id,
                                ref_wing_id,
                                ref_wing_cross_section_id,
                            )
                            if num_spanwise_panels_key in num_spanwise_panels_cache:
                                convergence_logger.debug(
                                    f"\t\t\t\tGetting the cached number of spanwise "
                                    f"Panels calculated for the #"
                                    f"{ref_wing_cross_section_id + 1} "
                                    f"WingCrossSection of {ref_airplane.name}'s "
                                    f"{ref_wing.name}..."
                                )

                                this_num_spanwise_panels = num_spanwise_panels_cache[
                                    num_spanwise_panels_key
                                ]
                            else:
                                # The way we calculate the correct number of spanwise
                                # Panels is to make skeleton Airplanes containing
                                # only one Wing with only the two WingCrossSections
                                # that make up the current Wing section. During
                                # initialization, the Airplane meshes its Wing,
                                # and we can then access the Wing's
                                # average_panel_aspect_ratio property. We repeat this
                                # process with increasing numbers of spanwise Panels,
                                # until we find the value that results in
                                # average_panel_aspect_ratio most closely matches the
                                # desired Panel aspect ratio. Initially, the first
                                # skeleton Airplane uses num_spanwise_panels=1.
                                # However, if we've already calculated a number of
                                # spanwise Panels for this Wing section with a
                                # coarser mesh (either in Panel aspect ratio,
                                # number of chordwise Panels, or both), then we know
                                # the current mesh must use at least this many
                                # spanwise Panels. Therefore, we can start the
                                # iterations with a higher number of spanwise Panels.
                                starting_num_spanwise_panels = 1

                                # Get the keys for the three coarser cases.
                                last_ar_key = (
                                    ar_id - 1,
                                    chord_id,
                                    ref_airplane_id,
                                    ref_wing_id,
                                    ref_wing_cross_section_id,
                                )
                                last_chord_key = (
                                    ar_id,
                                    chord_id - 1,
                                    ref_airplane_id,
                                    ref_wing_id,
                                    ref_wing_cross_section_id,
                                )
                                last_ar_and_chord_key = (
                                    ar_id - 1,
                                    chord_id - 1,
                                    ref_airplane_id,
                                    ref_wing_id,
                                    ref_wing_cross_section_id,
                                )

                                # Initialize the three coarser cases number of
                                # spanwise to be infinity, and update them if they
                                # exist in the cache.
                                last_ar_cache_val = np.inf
                                if last_ar_key in num_spanwise_panels_cache:
                                    last_ar_cache_val = num_spanwise_panels_cache[
                                        last_ar_key
                                    ]
                                last_chord_cache_val = np.inf
                                if last_chord_key in num_spanwise_panels_cache:
                                    last_chord_cache_val = num_spanwise_panels_cache[
                                        last_chord_key
                                    ]
                                last_ar_and_chord_cache_val = np.inf
                                if last_ar_and_chord_key in num_spanwise_panels_cache:
                                    last_ar_and_chord_cache_val = (
                                        num_spanwise_panels_cache[last_ar_and_chord_key]
                                    )

                                # To be conservative, take the minimum
                                # num_spanwise_panels of the three coarser cases. If
                                # at least one of the three cases has already been
                                # calculated, use that num_spanwise_panels as the
                                # starting value instead of 1.
                                last_cache_val = min(
                                    last_ar_cache_val,
                                    last_chord_cache_val,
                                    last_ar_and_chord_cache_val,
                                )
                                if last_cache_val != np.inf:
                                    starting_num_spanwise_panels = int(last_cache_val)

                                next_ref_wing_cross_section = ref_wing_cross_sections[
                                    ref_wing_cross_section_id + 1
                                ]

                                convergence_logger.debug(
                                    f"\t\t\t\tCalculating the number of spanwise "
                                    f"Panels for the #{ref_wing_cross_section_id + 1} "
                                    f"WingCrossSection of {ref_airplane.name}'s "
                                    f"{ref_wing.name}, with a starting value of "
                                    f"{starting_num_spanwise_panels}..."
                                )

                                # Iteratively find the correct number of spanwise
                                # Panels.
                                this_num_spanwise_panels = (
                                    _get_wing_section_num_spanwise_panels(
                                        panel_aspect_ratio,
                                        num_chordwise_panels,
                                        ref_wing.chordwise_spacing,
                                        ref_wing_cross_section,
                                        next_ref_wing_cross_section,
                                        starting_num_spanwise_panels,
                                    )
                                )

                                # Cache the calculated number of spanwise Panels for
                                # future use.
                                num_spanwise_panels_cache[num_spanwise_panels_key] = (
                                    this_num_spanwise_panels
                                )

                            convergence_logger.debug(
                                f"\t\t\t\tNumber of spanwise Panels: "
                                f"{this_num_spanwise_panels}"
                            )
                        else:
                            this_num_spanwise_panels = None

                        these_wing_cross_sections.append(
                            geometry.wing_cross_section.WingCrossSection(
                                # These values are copied from the reference
                                # WingCrossSection.
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
            this_problem = problems.SteadyProblem(
                airplanes=these_airplanes, operating_point=ref_operating_point
            )

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

            convergence_logger.info("\t\t\tStarting simulation...")

            # Run the steady solver and time how long it takes to execute.
            iter_start = time.time()
            this_solver.run()
            iter_stop = time.time()
            this_iter_time = iter_stop - iter_start

            # Create and fill ndarrays with each of this iteration's Airplanes'
            # combined load coefficients.
            theseCombinedForceCoefficients = np.zeros(len(these_airplanes), dtype=float)
            theseCombinedMomentCoefficients = np.zeros(
                len(these_airplanes), dtype=float
            )

            for airplane_id, airplane in enumerate(these_airplanes):
                _forceCoefficients_W = airplane.forceCoefficients_W
                assert _forceCoefficients_W is not None

                theseCombinedForceCoefficients[airplane_id] = np.linalg.norm(
                    _forceCoefficients_W
                )

                _momentCoefficients_W_CgP1 = airplane.momentCoefficients_W_CgP1
                assert _momentCoefficients_W_CgP1 is not None

                theseCombinedMomentCoefficients[airplane_id] = np.linalg.norm(
                    _momentCoefficients_W_CgP1
                )

            # Populate the ndarrays that store information from all the iterations with
            # the data from this iteration.
            combinedForceCoefficients[ar_id, chord_id, :] = (
                theseCombinedForceCoefficients
            )
            combinedMomentCoefficients[ar_id, chord_id, :] = (
                theseCombinedMomentCoefficients
            )
            iter_times[ar_id, chord_id] = this_iter_time

            convergence_logger.info(
                "\t\t\tSimulation completed in " + str(round(this_iter_time, 3)) + " s"
            )

            max_ar_pc = np.inf
            max_chord_pc = np.inf

            # If this isn't the first Panel aspect ratio, calculate the Panel aspect
            # ratio APE.
            if ar_id > 0:
                lastArCombinedForceCoefficients = combinedForceCoefficients[
                    ar_id - 1, chord_id, :
                ]
                lastArCombinedMomentCoefficients = combinedMomentCoefficients[
                    ar_id - 1, chord_id, :
                ]
                max_ar_force_pc = max(
                    100
                    * np.abs(
                        (
                            theseCombinedForceCoefficients
                            - lastArCombinedForceCoefficients
                        )
                        / lastArCombinedForceCoefficients
                    )
                )
                max_ar_moment_pc = max(
                    100
                    * np.abs(
                        (
                            theseCombinedMomentCoefficients
                            - lastArCombinedMomentCoefficients
                        )
                        / lastArCombinedMomentCoefficients
                    )
                )
                max_ar_pc = max(max_ar_force_pc, max_ar_moment_pc)

                convergence_logger.info(
                    "\t\t\tMaximum combined coefficient change from Panel aspect "
                    "ratio: " + str(round(max_ar_pc, 2)) + "%"
                )
            else:
                convergence_logger.info(
                    "\t\t\tMaximum combined coefficient change from Panel aspect "
                    "ratio: " + str(max_ar_pc)
                )

            # If this isn't the first number of chordwise Panels, calculate the
            # number of chordwise Panels APE.
            if chord_id > 0:
                lastChordCombinedForceCoefficients = combinedForceCoefficients[
                    ar_id, chord_id - 1, :
                ]
                lastChordCombinedMomentCoefficients = combinedMomentCoefficients[
                    ar_id, chord_id - 1, :
                ]
                max_chord_force_pc = max(
                    100
                    * np.abs(
                        (
                            theseCombinedForceCoefficients
                            - lastChordCombinedForceCoefficients
                        )
                        / lastChordCombinedForceCoefficients
                    )
                )
                max_chord_moment_pc = max(
                    100
                    * np.abs(
                        (
                            theseCombinedMomentCoefficients
                            - lastChordCombinedMomentCoefficients
                        )
                        / lastChordCombinedMomentCoefficients
                    )
                )
                max_chord_pc = max(max_chord_force_pc, max_chord_moment_pc)

                convergence_logger.info(
                    "\t\t\tMaximum combined coefficient change from number of "
                    "chordwise Panels: " + str(round(max_chord_pc, 2)) + "%"
                )
            else:
                convergence_logger.info(
                    "\t\t\tMaximum combined coefficient change from number of "
                    "chordwise Panels: " + str(max_chord_pc)
                )

            # Consider the Panel aspect ratio value to be saturated if it is equal to
            # 1. This is because a Panel aspect ratio of 1 is considered the maximum
            # degree of fineness.
            ar_saturated = panel_aspect_ratio == 1

            # Check if only one value for either the Panel aspect ratio or the number
            # of chordwise Panels were specified.
            single_ar = len(panel_aspect_ratios_list) == 1
            single_chord = len(num_chordwise_panels_list) == 1

            # Check if this iteration is converged with respect to the Panel aspect
            # ratio and/or the number of chordwise Panels.
            ar_converged = max_ar_pc < convergence_criteria
            chord_converged = max_chord_pc < convergence_criteria

            # Consider each convergence parameter to have "passed" if it is
            # converged, single, or saturated.
            ar_passed = ar_converged or single_ar or ar_saturated
            chord_passed = chord_converged or single_chord

            # If both convergence parameters have passed, then a converged or
            # semi-converged combination of parameters has been found and will be
            # returned.
            if ar_passed and chord_passed:
                if single_ar:
                    converged_ar_id = ar_id
                else:
                    # More than one Panel aspect ratio was tested.
                    if ar_converged:
                        # There is no big difference between this Panel aspect ratio
                        # and the last (coarser) Panel aspect ratio. Therefore,
                        # the last (coarser) Panel aspect ratio is converged.
                        converged_ar_id = ar_id - 1
                    else:
                        # There is a big difference between this Panel aspect ratio
                        # and the last (coarser) Panel aspect ratio. However,
                        # the Panel aspect ratio is one, so it's saturated.
                        # Therefore, this Panel aspect ratio is converged.
                        converged_ar_id = ar_id

                if single_chord:
                    converged_chord_id = chord_id
                else:
                    converged_chord_id = chord_id - 1

                converged_aspect_ratio = panel_aspect_ratios_list[converged_ar_id]
                converged_chordwise_panels = num_chordwise_panels_list[
                    converged_chord_id
                ]
                converged_iter_time = float(
                    iter_times[converged_ar_id, converged_chord_id]
                )

                if single_ar or single_chord:
                    convergence_logger.info("The analysis found a semi-converged case:")
                    if single_ar:
                        convergence_logger.warning(
                            "Panel aspect ratio convergence was not checked"
                        )
                    if single_chord:
                        convergence_logger.warning(
                            "Chordwise panels convergence was not checked"
                        )
                else:
                    convergence_logger.info("The analysis found a converged case:")

                convergence_logger.info(
                    "\tPanel aspect ratio: " + str(converged_aspect_ratio)
                )
                convergence_logger.info(
                    "\tChordwise Panels: " + str(converged_chordwise_panels)
                )
                convergence_logger.info(
                    "\tSimulation time: " + str(round(converged_iter_time, 3)) + " s"
                )
                convergence_logger.info("\tSpanwise Panels:")
                for airplane_id, airplane in enumerate(ref_airplanes):
                    convergence_logger.info("\t\t" + airplane.name + ":")
                    for wing_id, wing in enumerate(airplane.wings):
                        convergence_logger.info("\t\t\t" + wing.name + ":")
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
                            convergence_logger.info(
                                "\t\t\t\tWingCrossSection "
                                + str(wing_cross_section_id + 1)
                                + ": "
                                + str(num_spanwise_panels)
                            )

                return (
                    converged_aspect_ratio,
                    converged_chordwise_panels,
                )

    # If all iterations have been checked and none of them resulted in both
    # convergence parameters passing, then indicate that no converged case was found
    # and return values of None for the converged parameters.
    convergence_logger.info(
        "The analysis did not find a converged case within the given bounds"
    )
    return None, None


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
# TODO: If a converged mesh was found, consider also returning the converged solver.
def analyze_unsteady_convergence(
    ref_problem: problems.UnsteadyProblem,
    prescribed_wake: bool | np.bool_ = True,
    free_wake: bool | np.bool_ = True,
    num_cycles_bounds: tuple[int, int] | None = None,
    num_chords_bounds: tuple[int, int] | None = None,
    panel_aspect_ratio_bounds: tuple[int, int] = (4, 1),
    num_chordwise_panels_bounds: tuple[int, int] = (3, 12),
    convergence_criteria: float | int = 5.0,
    show_solver_progress: bool | np.bool_ = True,
) -> tuple[bool, int, int, int] | tuple[None, None, None, None]:
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
    Airplanes' final load coefficients are stored. As this function deals with
    UnsteadyProblems, it considers the final load coefficients to be the final-cycle's
    RMS load coefficients for UnsteadyProblems with variable geometry, and the final
    time step's load coefficients for static geometry cases. For each Airplane, the
    three final force coefficients are combined by taking their root-sum-square to find
    the resultant final force coefficient for that Airplane. Then, the absolute percent
    change (APE) of this Airplane's resultant final force coefficient is found between
    this iteration and the iterations with incrementally coarser meshes in all four
    parameter directions (wake state, wake length, Panel aspect ratio, and number of
    chordwise Panels). These steps are repeated for the three final moment coefficients.

    For each Airplane, the maximums of the resultant final force coefficient APEs and
    resultant final moment coefficient APEs are found. This leaves us with four maximum
    APEs per Airplane, one for each parameter direction. Next, we find the maximum
    across all the Airplanes for each parameter directions' maximum APE. Now, we are
    left with four maximum APEs total. If any of the four parameter direction APEs is
    below the convergence criteria, then this iteration has found a converged solution
    for that parameter direction.

    If an iteration's four APEs are all below the converged criteria, then we exit the
    nested for loops and return the converged parameters. However, the converged
    parameters are actually the values incrementally coarser than the final values (
    because the incrementally coarser values were found to be within the convergence
    criteria percent difference from the final values).

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

    :param ref_problem: The UnsteadyProblem whose converged parameters will be found.
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
    :param convergence_criteria: A positive number (int or float) that determines the
        point at which the function considers the simulation to have converged.
        Specifically, it is the maximum absolute percent change in the combined load
        coefficients. Therefore, it is in units of percent. Refer to the description in
        this function's docstring for more details on how it affects the solver. In
        short, set this value to 5.0 for a lenient convergence, and 1.0 for a strict
        convergence. Values are converted to floats internally. The default is 5.0.
    :param show_solver_progress: Set this to True to show the TQDM progress bar during
        each run of the unsteady solver. For showing progress bars and displaying log
        statements, set up logging using the setup_logging function. It can be a bool or
        a numpy bool and will be converted internally to a bool. The default is True.
    :return: A tuple of one bool and three ints. In order, they are the converged wake
        state (prescribed=True and free=False), the converged wake length (in number of
        cycles for non static geometries and number of chords for static geometries),
        the converged Panel aspect ratio, and the converged number of chordwise Panels.
        If the function could not find a set of converged parameters, it returns (None,
        None, None, None).
    """
    # Validate the ref_problem parameter.
    if not isinstance(ref_problem, problems.UnsteadyProblem):
        raise TypeError("ref_problem must be an UnsteadyProblem.")

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

    # Validate the num_chordwise_panels_bounds parameter.
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

    # Validate the convergence_criteria parameter.
    convergence_criteria = _parameter_validation.number_in_range_return_float(
        convergence_criteria, "convergence_criteria", min_val=0.0, min_inclusive=False
    )

    # Validate the show_solver_progress parameter.
    show_solver_progress = _parameter_validation.boolLike_return_bool(
        show_solver_progress, "show_solver_progress"
    )

    convergence_logger.info("Beginning convergence analysis...")

    ref_airplane_movements = ref_movement.airplane_movements
    ref_operating_point_movement = ref_movement.operating_point_movement

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
    combinedFinalLoadCoefficients = np.zeros(
        (
            len(wake_list),
            len(wake_lengths_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplane_movements),
            2,
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

    # Begin iterating through the outermost loop of wake states.
    for wake_id, wake in enumerate(wake_list):
        if wake:
            convergence_logger.info("\tWake type: prescribed")
        else:
            convergence_logger.info("\tWake type: free")

        # Begin iterating through the second loop of wake lengths.
        for length_id, wake_length in enumerate(wake_lengths_list):
            if static:
                convergence_logger.info("\t\tChord lengths: " + str(wake_length))
            else:
                convergence_logger.info("\t\tCycles: " + str(wake_length))

            # Begin iterating through the third loop of Panel aspect ratios.
            for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):
                convergence_logger.info(
                    "\t\t\tPanel aspect ratio: " + str(panel_aspect_ratio)
                )

                # Begin iterating through the innermost loop of number of chordwise
                # Panels.
                for chord_id, num_chordwise_panels in enumerate(
                    num_chordwise_panels_list
                ):
                    convergence_logger.info(
                        "\t\t\t\tChordwise Panels: " + str(num_chordwise_panels)
                    )

                    iteration += 1
                    convergence_logger.info(
                        "\t\t\t\t\tIteration number: "
                        + str(iteration)
                        + "/"
                        + str(num_iterations)
                    )

                    # Initialize an empty list for this iteration's base
                    # AirplaneMovements.
                    these_base_airplanes = []

                    # Create an empty list for the AirplaneMovement copies.
                    these_airplane_movements = []

                    # Now, we will begin iterating through this iteration's reference
                    # AirplaneMovements, WingMovements,
                    # and WingCrossSectionMovements, and creating copies of them.
                    # These copies will have identical parameters to their respective
                    # reference movements except for the number of spanwise Panels (
                    # which is based on the Panel aspect ratio), and the number of
                    # chordwise Panels.
                    #
                    # To do this, we iterate over the AirplaneMovements and perform a
                    # several step procedure:
                    # 1:    Reference the AirplaneMovement's base Airplane.
                    # 2:    Reference the AirplaneMovement's list of WingMovements.
                    # 3:    Create an empty list for the WingMovements' base Wing
                    #       copies.
                    # 4:    Create an empty list for the WingMovement copies.
                    # 5:    Iterate over the WingMovements.
                    #       5.1:    Reference the WingMovement's base Wing.
                    #       5.2:    Reference the WingMovement's list of
                    #               WingCrossSectionMovements.
                    #       5.3:    Create an empty list for the
                    #               WingCrossSectionMovements' base WingCrossSection
                    #               copies.
                    #       5.4:    Create an empty list for the
                    #               WingCrossSectionMovement copies.
                    #       5.5:    Iterate over the WingCrossSectionMovements.
                    #               5.5.1:  Reference the WingCrossSectionMovement's
                    #                       base WingCrossSection.
                    #               5.5.2:  Calculate the number of spanwise Panels
                    #                       that corresponds to the desired combination
                    #                       of Panel aspect ratio and number of
                    #                       chordwise Panels.
                    #               5.5.3:  Create a copy of the base WingCrossSection.
                    #               5.5.4:  Create a copy of the
                    #                       WingCrossSectionMovement.
                    #               5.5.5:  Append the base WingCrossSection copy to
                    #                       the list of base WingCrossSection copies.
                    #               5.5.6:  Append the WingCrossSectionMovement copy to
                    #                       the list of WingCrossSectionMovement
                    #                       copies.
                    #       5.6:    Create a copy of the base Wing.
                    #       5.7:    Create a copy of the WingMovement.
                    #       5.8:    Append the base Wing copy to the list  of base Wing
                    #               copies.
                    #       5.9:    Append the WingMovement copy to the list of
                    #               WingMovement copies.
                    # 6:    Create a copy of the base Airplane.
                    # 7:    Create a copy of the AirplaneMovement.
                    # 8:    Append the base Airplane copy to the list of base Airplane
                    #       copies.
                    # 9:    Append the AirplaneMovement copy to the list of
                    #       AirplaneMovement copies.
                    ref_airplane_movement: movements.airplane_movement.AirplaneMovement
                    for ref_airplane_movement_id, ref_airplane_movement in enumerate(
                        ref_airplane_movements
                    ):
                        # 1. Reference the AirplaneMovement's base Airplane.
                        ref_base_airplane = ref_airplane_movement.base_airplane

                        # 2. Reference the AirplaneMovement's list of WingMovements.
                        ref_wing_movements = ref_airplane_movement.wing_movements

                        # 3. Create an empty list for the WingMovements' base Wing
                        # copies.
                        these_base_wings = []

                        # 4: Create an empty list for the WingMovement copies.
                        these_wing_movements = []

                        # 5: Iterate over the WingMovements.
                        for ref_wing_movement_id, ref_wing_movement in enumerate(
                            ref_wing_movements
                        ):
                            # 5.1: Reference the WingMovement's base Wing.
                            ref_base_wing = ref_wing_movement.base_wing

                            # 5.2: Reference the WingMovement's list of
                            # WingCrossSectionMovements.
                            ref_wing_cross_section_movements = (
                                ref_wing_movement.wing_cross_section_movements
                            )

                            # 5.3: Create an empty list for the
                            # WingCrossSectionMovements' base WingCrossSection copies.
                            these_base_wing_cross_sections = []

                            # 5.4: Create an empty list for the
                            # WingCrossSectionMovement copies.
                            these_wing_cross_section_movements = []

                            # 5.5: Iterate over the WingCrossSectionMovements.
                            for (
                                ref_wing_cross_section_movement_id,
                                ref_wing_cross_section_movement,
                            ) in enumerate(ref_wing_cross_section_movements):
                                # 5.5.1: Reference the WingCrossSectionMovement's
                                # base WingCrossSection.
                                ref_base_wing_cross_section = (
                                    ref_wing_cross_section_movement.base_wing_cross_section
                                )

                                # 5.5.2: Calculate the number of spanwise Panels that
                                # corresponds to the desired combination of Panel
                                # aspect ratio and number of chordwise Panels.
                                if ref_wing_cross_section_movement_id < (
                                    len(ref_wing_cross_section_movements) - 1
                                ):
                                    # Check if we've already calculated the number of
                                    # spanwise Panels for this case/combination of
                                    # parameters.
                                    num_spanwise_panels_key = (
                                        ar_id,
                                        chord_id,
                                        ref_airplane_movement_id,
                                        ref_wing_movement_id,
                                        ref_wing_cross_section_movement_id,
                                    )
                                    if (
                                        num_spanwise_panels_key
                                        in num_spanwise_panels_cache
                                    ):
                                        convergence_logger.debug(
                                            f"\t\t\t\t\t\tGetting the cached number of "
                                            f"spanwise Panels calculated for the #"
                                            f"{ref_wing_cross_section_movement_id + 1} "
                                            f"WingCrossSection of "
                                            f"{ref_base_airplane.name}'s "
                                            f"{ref_base_wing.name}..."
                                        )

                                        this_num_spanwise_panels = (
                                            num_spanwise_panels_cache[
                                                num_spanwise_panels_key
                                            ]
                                        )
                                    else:
                                        # The way we calculate the correct number of
                                        # spanwise Panels is to make skeleton
                                        # Airplanes containing only one Wing with only
                                        # the two WingCrossSections that make up the
                                        # current Wing section. During
                                        # initialization, the Airplane meshes its
                                        # Wing, and we can then access the Wing's
                                        # average_panel_aspect_ratio property. We
                                        # repeat this process with increasing numbers
                                        # of spanwise Panels, until we find the value
                                        # that results in average_panel_aspect_ratio
                                        # most closely matches the desired Panel
                                        # aspect ratio. Initially, the first skeleton
                                        # Airplane uses num_spanwise_panels=1.
                                        # However, if we've already calculated a
                                        # number of spanwise Panels for this Wing
                                        # section with a coarser mesh (either in
                                        # Panel aspect ratio, number of chordwise
                                        # Panels, or both), then we know the current
                                        # mesh must use at least this many spanwise
                                        # Panels. Therefore, we can start the
                                        # iterations with a higher number of spanwise
                                        # Panels.
                                        starting_num_spanwise_panels = 1

                                        # Get the keys for the three coarser cases.
                                        last_ar_key = (
                                            ar_id - 1,
                                            chord_id,
                                            ref_airplane_movement_id,
                                            ref_wing_movement_id,
                                            ref_wing_cross_section_movement_id,
                                        )
                                        last_chord_key = (
                                            ar_id,
                                            chord_id - 1,
                                            ref_airplane_movement_id,
                                            ref_wing_movement_id,
                                            ref_wing_cross_section_movement_id,
                                        )
                                        last_ar_and_chord_key = (
                                            ar_id - 1,
                                            chord_id - 1,
                                            ref_airplane_movement_id,
                                            ref_wing_movement_id,
                                            ref_wing_cross_section_movement_id,
                                        )

                                        # Initialize the three coarser cases number
                                        # of spanwise to be infinity, and update them
                                        # if they exist in the cache.
                                        last_ar_cache_val = np.inf
                                        if last_ar_key in num_spanwise_panels_cache:
                                            last_ar_cache_val = (
                                                num_spanwise_panels_cache[last_ar_key]
                                            )
                                        last_chord_cache_val = np.inf
                                        if last_chord_key in num_spanwise_panels_cache:
                                            last_chord_cache_val = (
                                                num_spanwise_panels_cache[
                                                    last_chord_key
                                                ]
                                            )
                                        last_ar_and_chord_cache_val = np.inf
                                        if (
                                            last_ar_and_chord_key
                                            in num_spanwise_panels_cache
                                        ):
                                            last_ar_and_chord_cache_val = (
                                                num_spanwise_panels_cache[
                                                    last_ar_and_chord_key
                                                ]
                                            )

                                        # To be conservative, take the minimum
                                        # num_spanwise_panels of the three coarser
                                        # cases. If at least one of the three cases
                                        # has already been calculated, use that
                                        # num_spanwise_panels as the starting value
                                        # instead of 1.
                                        last_cache_val = min(
                                            last_ar_cache_val,
                                            last_chord_cache_val,
                                            last_ar_and_chord_cache_val,
                                        )
                                        if last_cache_val != np.inf:
                                            starting_num_spanwise_panels = int(
                                                last_cache_val
                                            )

                                        convergence_logger.debug(
                                            f"\t\t\t\t\t\tCalculating the number of "
                                            f"spanwise Panels for the #"
                                            f"{ref_wing_cross_section_movement_id + 1} "
                                            f"WingCrossSection of "
                                            f"{ref_base_airplane.name}'s "
                                            f"{ref_base_wing.name}, with a starting "
                                            f"value of "
                                            f"{starting_num_spanwise_panels}..."
                                        )

                                        # Iteratively find the correct number of
                                        # spanwise Panels.
                                        this_num_spanwise_panels = _get_wing_section_movement_num_spanwise_panels(
                                            panel_aspect_ratio,
                                            num_chordwise_panels,
                                            ref_base_wing.chordwise_spacing,
                                            ref_movement.airplanes[
                                                ref_airplane_movement_id
                                            ],
                                            ref_wing_movement_id,
                                            ref_wing_cross_section_movement_id,
                                            ref_wing_cross_section_movement_id + 1,
                                            starting_num_spanwise_panels,
                                            ref_problem.first_averaging_step,
                                        )

                                        # Cache the calculated number of spanwise
                                        # Panels for future use.
                                        num_spanwise_panels_cache[
                                            num_spanwise_panels_key
                                        ] = this_num_spanwise_panels

                                    convergence_logger.debug(
                                        f"\t\t\t\t\t\tNumber of spanwise Panels: "
                                        f"{this_num_spanwise_panels}"
                                    )
                                else:
                                    this_num_spanwise_panels = None

                                # 5.5.3: Create a copy of the base WingCrossSection.
                                this_base_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                                    # These values are copied from the reference base
                                    # WingCrossSection.
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

                                # 5.5.4: Create a copy of the WingCrossSectionMovement.
                                this_wing_cross_section_movement = movements.wing_cross_section_movement.WingCrossSectionMovement(
                                    # These values are copied from the reference
                                    # WingCrossSectionMovement.
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

                                # 5.5.5: Append the base WingCrossSection copy to the
                                # list of base WingCrossSection copies.
                                these_base_wing_cross_sections.append(
                                    this_base_wing_cross_section
                                )

                                # 5.5.6: Append the WingCrossSectionMovement copy to
                                # the list of WingCrossSectionMovement copies.
                                these_wing_cross_section_movements.append(
                                    this_wing_cross_section_movement
                                )

                            # 5.6: Create a copy of base Wing.
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

                            # 5.7: Create a copy of the WingMovement.
                            this_wing_movement = movements.wing_movement.WingMovement(
                                # These values are copied from the reference
                                # WingMovement.
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

                            # 5.8: Append the base Wing copy to the list of base Wing
                            # copies.
                            these_base_wings.append(this_base_wing)

                            # 5.9: Append the WingMovement copy to the list of
                            # WingMovement copies.
                            these_wing_movements.append(this_wing_movement)

                        # 6: Create a copy of the base Airplane.
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
                            # These values are copied from the reference
                            # AirplaneMovement.
                            ampCg_GP1_CgP1=ref_airplane_movement.ampCg_GP1_CgP1,
                            periodCg_GP1_CgP1=ref_airplane_movement.periodCg_GP1_CgP1,
                            spacingCg_GP1_CgP1=ref_airplane_movement.spacingCg_GP1_CgP1,
                            phaseCg_GP1_CgP1=ref_airplane_movement.phaseCg_GP1_CgP1,
                            # These values change.
                            base_airplane=this_base_airplane,
                            wing_movements=these_wing_movements,
                        )

                        # 8. Append the base Airplane copy to the list of base
                        # Airplane copies.
                        these_base_airplanes.append(this_base_airplane)

                        # 9. Append the AirplaneMovement copy to the list of
                        # AirplaneMovement copies.
                        these_airplane_movements.append(this_airplane_movement)

                    # Create a new Movement for this iteration.
                    if static:
                        this_movement = movements.movement.Movement(
                            airplane_movements=these_airplane_movements,
                            operating_point_movement=ref_operating_point_movement,
                            num_chords=wake_length,
                        )
                    else:
                        this_movement = movements.movement.Movement(
                            airplane_movements=these_airplane_movements,
                            operating_point_movement=ref_operating_point_movement,
                            num_cycles=wake_length,
                        )

                    # Create a new UnsteadyProblem for this iteration.
                    this_problem = problems.UnsteadyProblem(
                        movement=this_movement,
                        only_final_results=True,
                    )

                    # Create and run this iteration's
                    # UnsteadyRingVortexLatticeMethodSolver and time how long it
                    # takes to execute.
                    this_solver = unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        unsteady_problem=this_problem
                    )

                    convergence_logger.info("\t\t\t\t\tStarting simulation...")

                    iter_start = time.time()
                    this_solver.run(
                        prescribed_wake=wake,
                        calculate_streamlines=False,
                        show_progress=show_solver_progress,
                    )
                    iter_stop = time.time()
                    this_iter_time = iter_stop - iter_start

                    # Create and fill ndarrays with each of this iteration's Airplanes'
                    # combined final load coefficients.
                    theseCombinedFinalLoadCoefficients = np.zeros(
                        (len(these_airplane_movements), 2), dtype=float
                    )

                    for airplane_id, airplane in enumerate(these_airplane_movements):
                        # If this UnsteadyProblem is static, then get it's combined
                        # final load coefficients. If it's variable, get the combined
                        # final RMS load coefficients.
                        if static:
                            combinedFinalForceCoefficient = np.linalg.norm(
                                this_problem.finalForceCoefficients_W[airplane_id]
                            )
                            combinedFinalMomentCoefficient = np.linalg.norm(
                                this_problem.finalMomentCoefficients_W_CgP1[airplane_id]
                            )
                        else:
                            combinedFinalForceCoefficient = np.linalg.norm(
                                this_problem.finalRmsForceCoefficients_W[airplane_id]
                            )
                            combinedFinalMomentCoefficient = np.linalg.norm(
                                this_problem.finalRmsMomentCoefficients_W_CgP1[
                                    airplane_id
                                ]
                            )

                        theseCombinedFinalLoadCoefficients[airplane_id, 0] = (
                            combinedFinalForceCoefficient
                        )
                        theseCombinedFinalLoadCoefficients[airplane_id, 1] = (
                            combinedFinalMomentCoefficient
                        )

                    # Populate the ndarrays that store information from all the
                    # iterations with the data from this iteration.
                    combinedFinalLoadCoefficients[
                        wake_id, length_id, ar_id, chord_id, :, :
                    ] = theseCombinedFinalLoadCoefficients
                    iter_times[wake_id, length_id, ar_id, chord_id] = this_iter_time

                    convergence_logger.info(
                        "\t\t\t\t\tSimulation completed in "
                        + str(round(this_iter_time, 3))
                        + " s"
                    )

                    max_wake_pc = np.inf
                    max_length_pc = np.inf
                    max_ar_pc = np.inf
                    max_chord_pc = np.inf

                    # If this isn't the first wake state, calculate the wake state APE.
                    if wake_id > 0:
                        lastWakeCombinedFinalLoadCoefficients = (
                            combinedFinalLoadCoefficients[
                                wake_id - 1, length_id, ar_id, chord_id, :, :
                            ]
                        )
                        max_wake_pc = np.max(
                            100
                            * np.abs(
                                (
                                    theseCombinedFinalLoadCoefficients
                                    - lastWakeCombinedFinalLoadCoefficients
                                )
                                / lastWakeCombinedFinalLoadCoefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from wake "
                            "type: " + str(round(max_wake_pc, 2)) + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from wake "
                            "type: " + str(max_wake_pc)
                        )

                    # If this isn't the first wake length, calculate the wake length
                    # APE.
                    if length_id > 0:
                        lastLengthCombinedFinalLoadCoefficients = (
                            combinedFinalLoadCoefficients[
                                wake_id, length_id - 1, ar_id, chord_id, :, :
                            ]
                        )
                        max_length_pc = np.max(
                            100
                            * np.abs(
                                (
                                    theseCombinedFinalLoadCoefficients
                                    - lastLengthCombinedFinalLoadCoefficients
                                )
                                / lastLengthCombinedFinalLoadCoefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from wake "
                            "length: " + str(round(max_length_pc, 2)) + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from wake "
                            "length: " + str(max_length_pc)
                        )

                    # If this isn't the first Panel aspect ratio, calculate the Panel
                    # aspect ratio APE.
                    if ar_id > 0:
                        lastArCombinedFinalLoadCoefficients = (
                            combinedFinalLoadCoefficients[
                                wake_id, length_id, ar_id - 1, chord_id, :, :
                            ]
                        )
                        max_ar_pc = np.max(
                            100
                            * np.abs(
                                (
                                    theseCombinedFinalLoadCoefficients
                                    - lastArCombinedFinalLoadCoefficients
                                )
                                / lastArCombinedFinalLoadCoefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from Panel "
                            "aspect ratio: " + str(round(max_ar_pc, 2)) + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from Panel "
                            "aspect ratio: " + str(max_ar_pc)
                        )

                    # If this isn't the first number of chordwise Panels, calculate
                    # the number of chordwise Panels APE.
                    if chord_id > 0:
                        lastChordCombinedFinalLoadCoefficients = (
                            combinedFinalLoadCoefficients[
                                wake_id, length_id, ar_id, chord_id - 1, :, :
                            ]
                        )
                        max_chord_pc = np.max(
                            100
                            * np.abs(
                                (
                                    theseCombinedFinalLoadCoefficients
                                    - lastChordCombinedFinalLoadCoefficients
                                )
                                / lastChordCombinedFinalLoadCoefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from "
                            "chordwise Panels: " + str(round(max_chord_pc, 2)) + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum combined coefficient change from "
                            "chordwise Panels: " + str(max_chord_pc)
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

                    # Check if the iteration is converged with respect to any of the
                    # four convergence parameters.
                    wake_converged = max_wake_pc < convergence_criteria
                    length_converged = max_length_pc < convergence_criteria
                    ar_converged = max_ar_pc < convergence_criteria
                    chord_converged = max_chord_pc < convergence_criteria

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
                        if single_wake:
                            converged_wake_id = wake_id
                        else:
                            # We've tested both prescribed and free wakes.
                            if wake_converged:
                                # There isn't a big difference between the prescribed
                                # wake and free wake, so the prescribed wake is
                                # converged.
                                converged_wake_id = wake_id - 1
                            else:
                                # There is a big different difference between the
                                # prescribed wake and free wake, so the free wake is
                                # converged.
                                converged_wake_id = wake_id

                        if single_length:
                            converged_length_id = length_id
                        else:
                            converged_length_id = length_id - 1

                        if single_ar:
                            converged_ar_id = ar_id
                        else:
                            # We've tested more than one Panel aspect ratio.
                            if ar_converged:
                                # There is no big difference between this Panel aspect
                                # ratio and the last (coarser) Panel aspect ratio.
                                # Therefore, the last (coarser) Panel aspect ratio is
                                # converged.
                                converged_ar_id = ar_id - 1
                            else:
                                # There is a big difference between this Panel aspect
                                # ratio and the last (coarser) Panel aspect ratio.
                                # However, the Panel aspect ratio is one, so it's
                                # saturated. Therefore, this Panel aspect ratio is
                                # converged.
                                converged_ar_id = ar_id

                        if single_chord:
                            converged_chord_id = chord_id
                        else:
                            converged_chord_id = chord_id - 1

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
                            convergence_logger.info(
                                "The analysis found a semi-converged case:"
                            )
                            if single_wake:
                                convergence_logger.warning(
                                    "Wake type convergence not checked"
                                )
                            if single_length:
                                convergence_logger.warning(
                                    "Wake length convergence not checked"
                                )
                            if single_ar:
                                convergence_logger.warning(
                                    "Panel aspect ratio convergence not checked"
                                )
                            if single_chord:
                                convergence_logger.warning(
                                    "Chordwise Panels convergence not checked"
                                )
                        else:
                            convergence_logger.info(
                                "The analysis found a converged case:"
                            )

                        if converged_wake:
                            convergence_logger.info("\tWake type: prescribed")
                        else:
                            convergence_logger.info("\tWake type: free")

                        if static:
                            convergence_logger.info(
                                "\tChord lengths: " + str(converged_wake_length)
                            )
                        else:
                            convergence_logger.info(
                                "\tCycles: " + str(converged_wake_length)
                            )

                        convergence_logger.info(
                            "\tPanel aspect ratio: " + str(converged_aspect_ratio)
                        )
                        convergence_logger.info(
                            "\tChordwise Panels: " + str(converged_chordwise_panels)
                        )
                        convergence_logger.info(
                            "\tSimulation completed in "
                            + str(round(converged_iter_time, 3))
                            + " s"
                        )
                        convergence_logger.info("\tSpanwise Panels:")
                        for airplane_movement_id, airplane_movement in enumerate(
                            ref_airplane_movements
                        ):
                            base_airplane = airplane_movement.base_airplane
                            convergence_logger.info("\t\t" + base_airplane.name + ":")
                            for wing_movement_id, wing_movement in enumerate(
                                airplane_movement.wing_movements
                            ):
                                base_wing = wing_movement.base_wing
                                convergence_logger.info("\t\t\t" + base_wing.name + ":")
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
                                    convergence_logger.info(
                                        "\t\t\t\tWingCrossSection "
                                        + str(wing_cross_section_movement_id + 1)
                                        + ": "
                                        + str(num_spanwise_panels)
                                    )

                        return (
                            converged_wake,
                            converged_wake_length,
                            converged_aspect_ratio,
                            converged_chordwise_panels,
                        )

    # If all iterations have been checked and none of them resulted in all
    # convergence parameters passing, then indicate that no converged solution was
    # found and return values of None for the converged parameters.
    convergence_logger.info(
        "The analysis did not find a converged case within the given bounds"
    )
    return None, None, None, None


# TEST: Consider adding unit tests for this function.
def _get_wing_section_movement_num_spanwise_panels(
    desired_average_panel_aspect_ratio: int,
    num_chordwise_panels: int,
    chordwise_spacing: str,
    ref_airplanes: list[geometry.airplane.Airplane],
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
    :param ref_airplanes: A list of the Airplanes at each time step.
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
    # Slice the list of Airplanes to only the applicable ones. For cases with static
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

        convergence_logger.debug(
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

        convergence_logger.debug(
            f"\t\t\t\t\t\t\tNumber of spanwise Panels: "
            f"{num_spanwise_panels_at_step}"
        )

    return int(max(these_num_spanwise_panels))


# TEST: Consider adding unit tests for this function.
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


# TEST: Consider adding unit tests for this function.
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
