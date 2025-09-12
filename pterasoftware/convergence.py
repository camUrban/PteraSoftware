"""This module contains functions for analyzing the convergence of steady and
unsteady problems.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    analyze_steady_convergence: This function finds the converged parameters of a
    steady problem.

    analyze_unsteady_convergence: This function finds the converged parameters of an
    unsteady problem."""

import logging
import math
import time

import numpy as np

from . import geometry
from . import problems
from . import movement

from . import unsteady_ring_vortex_lattice_method
from . import steady_horseshoe_vortex_lattice_method
from . import steady_ring_vortex_lattice_method

convergence_logger = logging.getLogger("convergence")
convergence_logger.setLevel(logging.INFO)
logging.basicConfig()


def analyze_steady_convergence(
    ref_problem,
    solver_type,
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 12),
    convergence_criteria=5.0,
):
    """This function finds the converged parameters of a steady problem.

    Convergence is found by varying the problem's aircraft objects' panel aspect
    ratios and numbers of chordwise panels. These values are iterated over via two
    nested for loops (with the number of chordwise panels as the inner loop).

    With each new combination of these values, the problem is solved, and its
    resultant force and moment coefficients are stored. The force coefficients are
    combined by taking the vector norm. This is repeated for the moment coefficients.
    Then, absolute percent change (APE) of the resultant force coefficient is found
    between this iteration, and the iterations with incrementally coarser meshes
    in both the number of chordwise panels and panel aspect ratio. The process is
    repeated for to find the resultant moment coefficient APE.

    The maximums of the resultant force coefficient APEs and resultant moment
    coefficient APEs are found. This leaves us with two maximum APEs, one for the
    difference in the number of chordwise panels, and one for the difference in panel
    aspect ratio. If the chordwise panels APE is below the convergence criteria,
    this iteration has found a converged number of chordwise panels. The same is true
    for the panel aspect ratio APE.

    If an iteration's chordwise panels and the panel aspect ratio are both converged,
    then the solver will exit the loops and return the converged number of chordwise
    panels and panel aspect ratio. However, the converged parameters are actually
    the values incrementally coarser than the final values (because the
    incrementally coarser values were found to be within the convergence criteria
    percent difference from the final values).

    There are two edge cases to this function. The first is if the user inputs
    equal values for the coarsest and finest values of either the panel aspect ratio
    or the number of chordwise panels (i.e. panel_aspect_ratio_bounds=(2, 2)). Then,
    this parameter will not be iterated over, and convergence will only be checked
    for the other parameter.

    The second edge case happens if the panel aspect ratio has not converged at a
    value of 1. This is the gold standard value for panel aspect ratio, so the solver
    will return 1 for the converged value of panel aspect ratio. In the code below,
    this state is referred to as a "saturated" panel aspect ratio case.

    :param ref_problem: SteadyProblem
        This is the SteadyProblem object whose convergence will be analyzed.
    :param solver_type: str
        This parameter determines what type of steady solver will be used to analyze
        the problem. The options are "steady horseshoe vortex lattice method" and
        "steady ring vortex lattice method".
    :param panel_aspect_ratio_bounds: tuple, optional
        This parameter determines the range of panel aspect ratio sizes, from largest
        to smallest. For a given wing section, this value dictates the average panel
        body-frame-y length divided by the average body-frame-x width. Historically,
        these values range between 5 and 1. Values above 5 can be uses for a coarser
        mesh, but the minimum value should not be less than 1. The first value must
        be greater than or equal to the second value. The default value is (4, 1).
    :param num_chordwise_panels_bounds: tuple, optional
        This parameter determines the range of each wing section's number of
        chordwise panels from smallest to largest. The first value must be less than
        or equal to the second value. The default value is (3, 12).
    :param convergence_criteria: float, optional
        This parameter determines at what point the function continues the problem
        converged. Specifically, it is the absolute percent change in the resultant
        force coefficient or moment coefficient (whichever is higher). Therefore,
        it is in units of percent. Refer to the description above for more details on
        how it affects the solver. In short, set this value to 5.0 for a lenient
        convergence, and 1.0 for a strict convergence. The default value is 5.0.
    :return: list
        This function returns a list of two ints. In order, they are the converged of
        panel aspect ratio and the converged number of chordwise panels. If the
        function could not find a set of converged parameters, it returns values of
        None for all items in the list.
    """
    convergence_logger.info("Beginning convergence analysis.")

    ref_operating_point = ref_problem.operating_point
    ref_airplanes = ref_problem.airplanes

    panel_aspect_ratios_list = list(
        range(panel_aspect_ratio_bounds[0], panel_aspect_ratio_bounds[1] - 1, -1)
    )
    num_chordwise_panels_list = list(
        range(num_chordwise_panels_bounds[0], num_chordwise_panels_bounds[1] + 1)
    )

    # Initialize some empty arrays to hold attributes regarding each iteration. Going
    # forward, an "iteration" refers to a problem containing one of the combinations
    # of panel aspect ratio and number of chordwise panels.
    iter_times = np.zeros(
        (len(panel_aspect_ratios_list), len(num_chordwise_panels_list))
    )
    force_coefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplanes),
        )
    )
    moment_coefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplanes),
        )
    )

    iteration = 0
    num_iterations = len(panel_aspect_ratios_list) * len(num_chordwise_panels_list)

    # Begin iterating through the outer loop of panel aspect ratios.
    for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):
        convergence_logger.info("\tPanel aspect ratio: " + str(panel_aspect_ratio))

        # Begin iterating through the inner loop of number of chordwise panels.
        for chord_id, num_chordwise_panels in enumerate(num_chordwise_panels_list):
            convergence_logger.info(
                "\t\tChordwise panels: " + str(num_chordwise_panels)
            )

            iteration += 1
            convergence_logger.info(
                "\t\t\tIteration Number: " + str(iteration) + "/" + str(num_iterations)
            )

            # Initialize an empty list to hold this iteration's problem's airplanes.
            # Then, fill the list by making new copies of each of the reference
            # problem's airplanes with modified values for panel aspect ratio and
            # number of chordwise panels.
            these_airplanes = []
            for ref_airplane in ref_airplanes:

                ref_wings = ref_airplane.wings
                these_wings = []
                for ref_wing in ref_wings:

                    ref_wing_cross_sections = ref_wing.wing_cross_sections
                    these_wing_cross_sections = []
                    for (
                        ref_wing_cross_section_id,
                        ref_wing_cross_section,
                    ) in enumerate(ref_wing_cross_sections):

                        if ref_wing_cross_section_id < (
                            len(ref_wing_cross_sections) - 1
                        ):
                            next_ref_wing_cross_section = ref_wing_cross_sections[
                                ref_wing_cross_section_id + 1
                            ]
                            # ToDo: Modify this to allow for new geometry with custom
                            #  planes for the wing cross sections. As of now,
                            #  it fails for vertical wings.
                            section_length = (
                                next_ref_wing_cross_section.y_le
                                - ref_wing_cross_section.y_le
                            )
                            root_chord = ref_wing_cross_section.chord
                            tip_chord = next_ref_wing_cross_section.chord
                            section_area = section_length * (root_chord + tip_chord) / 2
                            section_standard_mean_chord = section_area / section_length

                            # As we can't directly specify the panel aspect ratio,
                            # calculate the number of spanwise panels that
                            # corresponds to the desired panel aspect ratio.
                            this_num_spanwise_panels = round(
                                (section_length * num_chordwise_panels)
                                / (section_standard_mean_chord * panel_aspect_ratio)
                            )

                            this_num_spanwise_panels = math.ceil(
                                this_num_spanwise_panels
                            )
                        else:
                            this_num_spanwise_panels = 0

                        these_wing_cross_sections.append(
                            geometry.wing_cross_section.WingCrossSection(
                                # These values are copied from the reference wing
                                # cross section.
                                x_le=ref_wing_cross_section.x_le,
                                y_le=ref_wing_cross_section.y_le,
                                z_le=ref_wing_cross_section.z_le,
                                chord=ref_wing_cross_section.chord,
                                twist=ref_wing_cross_section.twist,
                                control_surface_type=ref_wing_cross_section.control_surface_type,
                                control_surface_hinge_point=ref_wing_cross_section.control_surface_hinge_point,
                                control_surface_deflection=ref_wing_cross_section.control_surface_deflection,
                                spanwise_spacing=ref_wing_cross_section.spanwise_spacing,
                                # These values change.
                                num_spanwise_panels=this_num_spanwise_panels,
                                airfoil=geometry.airfoil.Airfoil(
                                    name=ref_wing_cross_section.airfoil.name,
                                    coordinates=ref_wing_cross_section.airfoil.coordinates,
                                    repanel=ref_wing_cross_section.airfoil.repanel,
                                    n_points_per_side=ref_wing_cross_section.airfoil.n_points_per_side,
                                ),
                            )
                        )

                    these_wings.append(
                        geometry.wing.Wing(
                            # These values are copied from this reference wing.
                            name=ref_wing.name,
                            x_le=ref_wing.x_le,
                            y_le=ref_wing.y_le,
                            z_le=ref_wing.z_le,
                            symmetric=ref_wing.symmetric,
                            chordwise_spacing=ref_wing.chordwise_spacing,
                            # These values change.
                            num_chordwise_panels=num_chordwise_panels,
                            wing_cross_sections=these_wing_cross_sections,
                        )
                    )

                these_airplanes.append(
                    geometry.airplane.Airplane(
                        # These values are copied from the reference airplane.
                        name=ref_airplane.name,
                        x_ref=ref_airplane.x_ref,
                        y_ref=ref_airplane.y_ref,
                        z_ref=ref_airplane.z_ref,
                        weight=ref_airplane.weight,
                        # These are kept as None so that they are recalculated with
                        # this airplane's mesh.
                        s_ref=None,
                        c_ref=None,
                        b_ref=None,  # This value changes.
                        wings=these_wings,
                    )
                )

            # Create a new problem for this iteration.
            this_problem = problems.SteadyProblem(
                airplanes=these_airplanes, operating_point=ref_operating_point
            )

            # Create this iteration's solver based on the type specified.
            if solver_type == "steady horseshoe vortex lattice method":
                this_solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                    steady_problem=this_problem,
                )
            elif solver_type == "steady ring vortex lattice method":
                this_solver = steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
                    steady_problem=this_problem,
                )
            else:
                raise Exception("You entered an invalid type of solver.")

            del this_problem

            # Run the solver and time how long it takes to execute.
            iter_start = time.time()
            this_solver.run(logging_level="Critical")
            iter_stop = time.time()
            this_iter_time = iter_stop - iter_start

            # Create and fill arrays with each of this iteration's airplane's
            # resultant force and moment coefficients.
            these_force_coefficients = np.zeros(len(these_airplanes))
            these_moment_coefficients = np.zeros(len(these_airplanes))
            for airplane_id, airplane in enumerate(these_airplanes):
                these_force_coefficients[airplane_id] = np.linalg.norm(
                    airplane.total_near_field_force_coefficients_wind_axes
                )
                these_moment_coefficients[airplane_id] = np.linalg.norm(
                    airplane.total_near_field_moment_coefficients_wind_axes
                )

            # Populate the arrays that store information of all the iterations with
            # the data from this iteration.
            force_coefficients[ar_id, chord_id, :] = these_force_coefficients
            moment_coefficients[ar_id, chord_id, :] = these_moment_coefficients
            iter_times[ar_id, chord_id] = this_iter_time

            convergence_logger.info(
                "\t\t\tIteration Time: " + str(round(this_iter_time, 3)) + " s"
            )

            max_ar_pc = np.inf
            max_chord_pc = np.inf

            # If this isn't the first panel aspect ratio, calculate the panel aspect
            # ratio APE.
            if ar_id > 0:
                last_ar_force_coefficients = force_coefficients[ar_id - 1, chord_id, :]
                last_ar_moment_coefficients = moment_coefficients[
                    ar_id - 1, chord_id, :
                ]
                max_ar_force_pc = max(
                    100
                    * np.abs(
                        (these_force_coefficients - last_ar_force_coefficients)
                        / last_ar_force_coefficients
                    )
                )
                max_ar_moment_pc = max(
                    100
                    * np.abs(
                        (these_moment_coefficients - last_ar_moment_coefficients)
                        / last_ar_moment_coefficients
                    )
                )
                max_ar_pc = max(max_ar_force_pc, max_ar_moment_pc)

                convergence_logger.info(
                    "\t\t\tMaximum coefficient change from panel aspect ratio: "
                    + str(round(max_ar_pc, 2))
                    + "%"
                )
            else:
                convergence_logger.info(
                    "\t\t\tMaximum coefficient change from panel aspect ratio: "
                    + str(max_ar_pc)
                )

            # If this isn't the first number of chordwise panels, calculate the
            # number of chordwise panels APE.
            if chord_id > 0:
                last_chord_force_coefficients = force_coefficients[
                    ar_id, chord_id - 1, :
                ]
                last_chord_moment_coefficients = moment_coefficients[
                    ar_id, chord_id - 1, :
                ]
                max_chord_force_pc = max(
                    100
                    * np.abs(
                        (these_force_coefficients - last_chord_force_coefficients)
                        / last_chord_force_coefficients
                    )
                )
                max_chord_moment_pc = max(
                    100
                    * np.abs(
                        (these_moment_coefficients - last_chord_moment_coefficients)
                        / last_chord_moment_coefficients
                    )
                )
                max_chord_pc = max(max_chord_force_pc, max_chord_moment_pc)

                convergence_logger.info(
                    "\t\t\tMaximum coefficient change from chordwise panels: "
                    + str(round(max_chord_pc, 2))
                    + "%"
                )
            else:
                convergence_logger.info(
                    "\t\t\tMaximum coefficient change from chordwise panels: "
                    + str(max_chord_pc)
                )

            # Consider the panel aspect ratio value to be saturated if it is equal to
            # 1. This is because a panel aspect ratio of 1 is considered the maximum
            # degree of fineness.
            ar_saturated = panel_aspect_ratio == 1

            # Check if the user only specified one value for either the panel aspect
            # ratio or the number of chordwise panels.
            single_ar = len(panel_aspect_ratios_list) == 1
            single_chord = len(num_chordwise_panels_list) == 1

            # Check if the iteration calculated that it is converged with respect to
            # the panel aspect ratio and or the number of chordwise panels.
            ar_converged = max_ar_pc < convergence_criteria
            chord_converged = max_chord_pc < convergence_criteria

            # Consider each convergence parameter to have passed it is converged,
            # single, or saturated.
            ar_passed = ar_converged or single_ar or ar_saturated
            chord_passed = chord_converged or single_chord

            # If both convergence parameters have passed, then the solver has found a
            # converged or semi-converged value and will return the converged
            # parameters.
            if ar_passed and chord_passed:
                if single_ar:
                    converged_ar_id = ar_id
                else:
                    # We've tested more than one panel aspect ratio.
                    if ar_converged:
                        # There is no big difference between this panel aspect ratio
                        # and the last (coarser) panel aspect ratio. Therefore,
                        # the last (coarser) panel aspect ratio is converged.
                        converged_ar_id = ar_id - 1
                    else:
                        # There is a big difference between this panel aspect ratio
                        # and the last (coarser) panel aspect ratio. However,
                        # the panel aspect ratio is one, so it's saturated.
                        # Therefore, this panel aspect ratio is converged.
                        converged_ar_id = ar_id

                if single_chord:
                    converged_chord_id = chord_id
                else:
                    converged_chord_id = chord_id - 1

                converged_chordwise_panels = num_chordwise_panels_list[
                    converged_chord_id
                ]
                converged_aspect_ratio = panel_aspect_ratios_list[converged_ar_id]
                converged_iter_time = iter_times[converged_ar_id, converged_chord_id]

                if single_ar or single_chord:
                    convergence_logger.info("The analysis found a semi-converged mesh:")
                    if single_ar:
                        convergence_logger.warning(
                            "Panel aspect ratio convergence not checked."
                        )
                    if single_chord:
                        convergence_logger.warning(
                            "Chordwise panels convergence not checked."
                        )
                else:
                    convergence_logger.info("The analysis found a converged mesh:")

                convergence_logger.info(
                    "\tPanel aspect ratio: " + str(converged_aspect_ratio)
                )
                convergence_logger.info(
                    "\tChordwise panels: " + str(converged_chordwise_panels)
                )
                convergence_logger.info(
                    "\tIteration time: " + str(round(converged_iter_time, 3)) + " s"
                )

                return [
                    converged_aspect_ratio,
                    converged_chordwise_panels,
                ]

    # If all iterations have been checked and none of them resulted in both
    # convergence parameters passing, then indicate that no converged solution was
    # found and return values of None for the converged parameters.
    convergence_logger.info("The analysis did not find a converged mesh.")
    return [None, None]


# ToDo: Add the new parameters to the documentation.
def analyze_unsteady_convergence(
    ref_problem,
    prescribed_wake=True,
    free_wake=True,
    num_cycles_bounds=(2, 4),
    num_chords_bounds=(3, 7),
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 12),
    coefficient_mask=None,
    convergence_criteria=5.0,
):
    """This function finds the converged parameters of an unsteady problem.

    Convergence is found by varying if the solver's wake state (prescribed or free),
    the final length of the problem's wake (in number of chord lengths for static
    geometry and number of flap cycles for variable geometry), the airplanes' panel
    aspect ratios, and the airplanes' numbers of chordwise panels. These values are
    iterated over via four nested loops. The outermost loop is the wake state. The
    next loop is the wake length. The loop after that is the panel aspect ratios,
    and the final loop is the number of chordwise panels.

    With each new combination of these values, the problem is solved, and its
    resultant force and moment coefficients are stored. The force coefficients are
    combined by taking the vector norm. This is repeated for the moment coefficients.
    Then, absolute percent change (APE) of the resultant force coefficient is found
    between this iteration, and the iterations with the incrementally coarser meshes
    in all four parameters (wake state, wake length, panel aspect ratio,
    and number of chordwise panels). The process is repeated for to find the
    resultant moment coefficient APE.

    The maximums of the resultant force coefficient APEs and resultant moment
    coefficient APEs are found. This leaves us with four maximum APEs, one for each
    parameter. If any of the parameter's APE is below the iteration has found a
    converged solution for that parameter.

    If an iteration's four APEs are all below the converged criteria, then the solver
    will exit the loops and return the converged parameters. However, the converged
    parameters are actually the values incrementally coarser than the final values (
    because the incrementally coarser values were found to be within the convergence
    criteria percent difference from the final values).

    There are two edge cases to this function. The first occurs when the user
    indicates that they only want check a single value for any of the four parameters
    (i.e. panel_aspect_ratio_bounds=(2, 2) or (prescribed_wake=True and
    free_wake=False)). Then, this parameter will not be iterated over,
    and convergence will only be checked for the other parameters.

    The second edge case happens if the panel aspect ratio has not converged at a
    value of 1 or if the wake state hasn't converged once it is equal to a free wake.
    These conditions are the gold standards for panel aspect ratio and wake state,
    so the solver will return 1 for the converged value of panel aspect ratio and a
    free wake for the converged wake state. In the code below, this state is referred
    to as a "saturated" panel aspect ratio or wake state.

    :param ref_problem: UnsteadyProblem
        This is the UnsteadyProblem object whose convergence will be analyzed.
    :param prescribed_wake: bool, optional
        This parameter determines if a prescribed wake type should be analyzed. If
        this parameter is false, then the free_wake parameter must be set to True.
        The default value is True.
    :param free_wake: bool, optional
        This parameter determines if a free wake type should be analyzed.  If this
        parameter is false, then the prescribed_wake parameter must be set to
        True. The default value is True.
    :param num_cycles_bounds: tuple, optional
        This parameter determines the range of wake lengths, measured in number of
        cycles to simulate. If the problem has static geometry, it will be ignored,
        and the num_chords_bounds parameter will control the wake lengths instead.
        Reasonable values range from 1 to 10, depending strongly on the Strouhal
        number. The first value must be less than or equal to the second value. The
        default value is (1, 4).
    :param num_chords_bounds: tuple, optional
        This parameter determines the range of wake lengths, measured in number of
        reference chords the wake should extend to. If the problem has variable
        geometry, it will be ignored, and the num_cycles_bounds parameter will
        control the wake lengths instead. Reasonable values range from 3 to 20. The
        first value must be less than or equal to the second value. The default value
        is (3, 7).
    :param panel_aspect_ratio_bounds: tuple, optional
        This parameter determines the range of panel aspect ratios, from largest to
        smallest. For a given wing section, this value dictates the average panel
        body-frame-y length divided by the average body-frame-x width. Historically,
        these values range between 5 and 1. Values above 5 can be uses for a coarser
        mesh, but the minimum value should not be less than 1. The first value must
        be greater than or equal to the second value. The default value is ( , 1).
    :param num_chordwise_panels_bounds: tuple, optional
        This parameter determines the range of each wing section's number of
        chordwise panels from smallest to largest. The first value must be less than
        or equal to the second value. The default value is (3, 12).
    :param convergence_criteria: float, optional
        This parameter determines at what point the function continues the problem
        converged. Specifically, it is the absolute percent change in the resultant
        force coefficient or moment coefficient (whichever is higher). Therefore,
        it is in units of percent. Refer to the description above for more details on
        how it affects the solver. In short, set this value to 5.0 for a lenient
        convergence, and 1.0 for a strict convergence. The default value is 5.0.
    :return: list
        This function returns a list of four ints. In order, they are the converged
        wake state, the converged wake length, the converged of panel aspect ratio
        and the converged number of chordwise panels. If the function could not find
        a set of converged parameters, it returns values of None for all items in the
        list.
    """
    convergence_logger.info("Beginning convergence analysis.")

    ref_movement = ref_problem.movement
    is_static = ref_movement.get_max_period() == 0

    ref_airplane_movements = ref_movement.airplane_movements
    ref_operating_point_movement = ref_movement.operating_point_movement

    if coefficient_mask is None:
        coefficient_mask = [True, True, True, True, True, True]

    num_coefficients = coefficient_mask.count(True)

    wake_list = []
    if prescribed_wake:
        wake_list.append(True)
    if free_wake:
        wake_list.append(False)

    # If this problem has static geometry, base the wake length on the number of
    # chords parameter. Otherwise, base it on the number of cycles parameter.
    if is_static:
        wake_lengths_list = list(range(num_chords_bounds[0], num_chords_bounds[1] + 1))
    else:
        wake_lengths_list = list(range(num_cycles_bounds[0], num_cycles_bounds[1] + 1))

    panel_aspect_ratios_list = list(
        range(panel_aspect_ratio_bounds[0], panel_aspect_ratio_bounds[1] - 1, -1)
    )
    num_chordwise_panels_list = list(
        range(num_chordwise_panels_bounds[0], num_chordwise_panels_bounds[1] + 1)
    )

    # Initialize some empty arrays to hold attributes regarding each iteration. Going
    # forward, an "iteration" refers to a problem containing one of the combinations
    # of the wake state, wake length, panel aspect ratio, and number of chordwise
    # panels parameters.
    iter_times = np.zeros(
        (
            len(wake_list),
            len(wake_lengths_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
        )
    )
    coefficients = np.zeros(
        (
            len(wake_list),
            len(wake_lengths_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplane_movements),
            num_coefficients,
        )
    )

    iteration = 0
    num_iterations = (
        len(wake_list)
        * len(wake_lengths_list)
        * len(panel_aspect_ratios_list)
        * len(num_chordwise_panels_list)
    )

    # Begin iterating through the first loop of wake states.
    for wake_id, wake in enumerate(wake_list):
        if wake:
            convergence_logger.info("\tWake type: prescribed")
        else:
            convergence_logger.info("\tWake type: free")

        # Begin iterating through the second loop of wake lengths.
        for length_id, wake_length in enumerate(wake_lengths_list):
            if is_static:
                convergence_logger.info("\t\tChord lengths: " + str(wake_length))
            else:
                convergence_logger.info("\t\tCycles: " + str(wake_length))

            # Begin iterating through the third loop of panel aspect ratios.
            for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):
                convergence_logger.info(
                    "\t\t\tPanel aspect ratio: " + str(panel_aspect_ratio)
                )

                # Begin iterating through the fourth and innermost loop of number of
                # chordwise panels.
                for chord_id, num_chordwise_panels in enumerate(
                    num_chordwise_panels_list
                ):
                    convergence_logger.info(
                        "\t\t\t\tChordwise panels: " + str(num_chordwise_panels)
                    )

                    iteration += 1
                    convergence_logger.info(
                        "\t\t\t\t\tIteration Number: "
                        + str(iteration)
                        + "/"
                        + str(num_iterations)
                    )

                    # Create an empty list for the airplane movement base objects.
                    these_base_airplanes = []

                    # Create an empty list for the airplane movement copies.
                    these_airplane_movements = []

                    # Now, we will begin iterating through the reference movement's
                    # sub-movements, and making copies. These copies will have
                    # identical parameters to their respective reference
                    # (sub-)movements except for the number of spanwise panels (which
                    # is based on the panel aspect ratio), and the number of
                    # chordwise panels. For each (sub-)movement, there is a 9-step
                    # process:
                    # 1. Reference this (sub-)movement's base object.
                    # 2. Reference this (sub-)movement's list of (sub-)sub-movements.
                    # 3. Create an empty list for the (sub-)sub-movement base objects.
                    # 4: Create an empty list for the (sub-)sub-movement copies.
                    # 5: Iterate over the (sub-)sub-movements.
                    # 6: Create a copy of the base object.
                    # 7. Create a copy of the new (sub-)movement.
                    # 8. Append the new base object to the list of new base objects.
                    # 9. Append the new (sub-)movement to the list of new (
                    # sub-)movements.
                    for ref_airplane_movement in ref_airplane_movements:
                        # 1. Reference this (sub-)movement's base object.
                        ref_base_airplane = ref_airplane_movement.base_airplane

                        # 2. Reference this (sub-)movement's list of (
                        # sub-)sub-movements.
                        ref_wing_movements = ref_airplane_movement.wing_movements

                        # 3. Create an empty list for the (sub-)sub-movement base
                        # objects.
                        these_base_wings = []

                        # 4: Create an empty list for the (sub-)sub-movement copies.
                        these_wing_movements = []

                        # 5: Iterate over the (sub-)sub-movements.
                        for ref_wing_movement in ref_wing_movements:
                            # 1. Reference this (sub-)movement's base object.
                            ref_base_wing = ref_wing_movement.base_wing

                            # 2. Reference this (sub-)movement's list of (
                            # sub-)sub-movements.
                            ref_wing_cross_section_movements = (
                                ref_wing_movement.wing_cross_section_movements
                            )

                            # 3. Create an empty list for the (sub-)sub-movement base
                            # objects.
                            these_base_wing_cross_sections = []

                            # 4: Create an empty list for the sub-movement copies.
                            these_wing_cross_section_movements = []

                            # 5: Iterate over the (sub-)sub-movements.
                            for (
                                ref_wing_cross_section_movement_id,
                                ref_wing_cross_section_movement,
                            ) in enumerate(ref_wing_cross_section_movements):
                                # 1. Reference this (sub-)movement's base object.
                                ref_base_wing_cross_section = (
                                    ref_wing_cross_section_movement.base_wing_cross_section
                                )

                                # As we can't directly specify the panel aspect
                                # ratio, calculate the number of spanwise panels that
                                # corresponds to the desired panel aspect ratio.
                                if ref_wing_cross_section_movement_id < (
                                    len(ref_wing_cross_section_movements) - 1
                                ):
                                    next_ref_base_wing_cross_section = (
                                        ref_wing_cross_section_movements[
                                            ref_wing_cross_section_movement_id + 1
                                        ].base_wing_cross_section
                                    )
                                    # ToDo: Modify this to allow for new geometry
                                    #  with custom planes for the wing cross
                                    #  sections. As of now, it fails for vertical wings.
                                    section_length = (
                                        next_ref_base_wing_cross_section.y_le
                                        - ref_base_wing_cross_section.y_le
                                    )
                                    root_chord = ref_base_wing_cross_section.chord
                                    tip_chord = next_ref_base_wing_cross_section.chord
                                    section_area = (
                                        section_length * (root_chord + tip_chord) / 2
                                    )
                                    section_standard_mean_chord = (
                                        section_area / section_length
                                    )

                                    this_num_spanwise_panels = round(
                                        (section_length * num_chordwise_panels)
                                        / (
                                            section_standard_mean_chord
                                            * panel_aspect_ratio
                                        )
                                    )

                                    this_num_spanwise_panels = math.ceil(
                                        this_num_spanwise_panels
                                    )
                                else:
                                    this_num_spanwise_panels = 0

                                # 2. Reference this (sub-)movement's list of (
                                # sub-)sub-movements.
                                # N/A

                                # 3. Create an empty list for the (sub-)sub-movement
                                # base objects.
                                # N/A

                                # 4: Create an empty list for the sub-movement copies.
                                # N/A

                                # 5: Iterate over the (sub-)sub-movements.
                                # N/A

                                # 6: Create a copy of the base object.
                                this_base_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                                    # These values are copied from the reference base
                                    # wing cross section.
                                    x_le=ref_base_wing_cross_section.x_le,
                                    y_le=ref_base_wing_cross_section.y_le,
                                    z_le=ref_base_wing_cross_section.z_le,
                                    chord=ref_base_wing_cross_section.chord,
                                    twist=ref_base_wing_cross_section.twist,
                                    control_surface_type=ref_base_wing_cross_section.control_surface_type,
                                    control_surface_hinge_point=ref_base_wing_cross_section.control_surface_hinge_point,
                                    control_surface_deflection=ref_base_wing_cross_section.control_surface_deflection,
                                    spanwise_spacing=ref_base_wing_cross_section.spanwise_spacing,
                                    # These values change.
                                    num_spanwise_panels=this_num_spanwise_panels,
                                    airfoil=geometry.airfoil.Airfoil(
                                        name=ref_base_wing_cross_section.airfoil.name,
                                        coordinates=ref_base_wing_cross_section.airfoil.coordinates,
                                        repanel=ref_base_wing_cross_section.airfoil.repanel,
                                        n_points_per_side=ref_base_wing_cross_section.airfoil.n_points_per_side,
                                    ),
                                )

                                # 7. Create a copy of the new (sub-)movement.
                                this_wing_cross_section_movement = movement.WingCrossSectionMovement(
                                    # These values are copied from the reference
                                    # object.
                                    sweeping_amplitude=ref_wing_cross_section_movement.sweeping_amplitude,
                                    sweeping_period=ref_wing_cross_section_movement.sweeping_period,
                                    sweeping_spacing=ref_wing_cross_section_movement.sweeping_spacing,
                                    custom_sweep_function=ref_wing_cross_section_movement.custom_sweep_function,
                                    pitching_amplitude=ref_wing_cross_section_movement.pitching_amplitude,
                                    pitching_period=ref_wing_cross_section_movement.pitching_period,
                                    pitching_spacing=ref_wing_cross_section_movement.pitching_spacing,
                                    custom_pitch_function=ref_wing_cross_section_movement.custom_pitch_function,
                                    heaving_amplitude=ref_wing_cross_section_movement.heaving_amplitude,
                                    heaving_period=ref_wing_cross_section_movement.heaving_period,
                                    heaving_spacing=ref_wing_cross_section_movement.heaving_spacing,
                                    custom_heave_function=ref_wing_cross_section_movement.custom_heave_function,
                                    # This value is new.
                                    base_wing_cross_section=this_base_wing_cross_section,
                                )

                                # 8. Append the new base object to the list of new
                                # base objects.
                                these_base_wing_cross_sections.append(
                                    this_base_wing_cross_section
                                )

                                # 9. Append the new (sub-)movement to the list of new
                                # (sub-)movements.
                                these_wing_cross_section_movements.append(
                                    this_wing_cross_section_movement
                                )

                            # 6: Create a copy of the base object.
                            this_base_wing = geometry.wing.Wing(
                                # These values are copied from this reference base
                                # wing.
                                name=ref_base_wing.name,
                                x_le=ref_base_wing.x_le,
                                y_le=ref_base_wing.y_le,
                                z_le=ref_base_wing.z_le,
                                symmetric=ref_base_wing.symmetric,
                                chordwise_spacing=ref_base_wing.chordwise_spacing,
                                # These values change.
                                num_chordwise_panels=num_chordwise_panels,
                                wing_cross_sections=these_base_wing_cross_sections,
                            )

                            # 7. Create a copy of the new (sub-)movement.
                            this_wing_movement = movement.WingMovement(
                                # These values are copied from this reference wing
                                # movement.
                                x_le_amplitude=ref_wing_movement.x_le_amplitude,
                                x_le_period=ref_wing_movement.x_le_period,
                                x_le_spacing=ref_wing_movement.x_le_spacing,
                                y_le_amplitude=ref_wing_movement.y_le_amplitude,
                                y_le_period=ref_wing_movement.y_le_period,
                                y_le_spacing=ref_wing_movement.y_le_spacing,
                                z_le_amplitude=ref_wing_movement.z_le_amplitude,
                                z_le_period=ref_wing_movement.z_le_period,
                                z_le_spacing=ref_wing_movement.z_le_spacing,
                                # These values change.
                                base_wing=this_base_wing,
                                wing_cross_sections_movements=these_wing_cross_section_movements,
                            )

                            # 8. Append the new base object to the list of new base
                            # objects.
                            these_base_wings.append(this_base_wing)

                            # 9. Append the new (sub-)movement to the list of new (
                            # sub-)movements.
                            these_wing_movements.append(this_wing_movement)

                        # 6: Create a copy of the base object.
                        this_base_airplane = geometry.airplane.Airplane(
                            # These values are copied from the reference base airplane.
                            name=ref_base_airplane.name,
                            x_ref=ref_base_airplane.x_ref,
                            y_ref=ref_base_airplane.y_ref,
                            z_ref=ref_base_airplane.z_ref,
                            weight=ref_base_airplane.weight,
                            # These are kept as None so that they are recalculated
                            # with this airplane's mesh.
                            s_ref=None,
                            c_ref=None,
                            b_ref=None,  # This value changes.
                            wings=these_base_wings,
                        )

                        # 7. Create a copy of the new (sub-)movement.
                        this_airplane_movement = movement.AirplaneMovement(
                            # These values are copied from this reference airplane
                            # movement.
                            x_ref_amplitude=ref_airplane_movement.x_ref_amplitude,
                            x_ref_period=ref_airplane_movement.x_ref_period,
                            x_ref_spacing=ref_airplane_movement.x_ref_spacing,
                            y_ref_amplitude=ref_airplane_movement.y_ref_amplitude,
                            y_ref_period=ref_airplane_movement.y_ref_period,
                            y_ref_spacing=ref_airplane_movement.y_ref_spacing,
                            z_ref_amplitude=ref_airplane_movement.z_ref_amplitude,
                            z_ref_period=ref_airplane_movement.z_ref_period,
                            z_ref_spacing=ref_airplane_movement.z_ref_spacing,
                            # These values change.
                            base_airplane=this_base_airplane,
                            wing_movements=these_wing_movements,
                        )

                        # 8. Append the new base object to the list of new base
                        # objects.
                        these_base_airplanes.append(this_base_airplane)

                        # 9. Append the new (sub-)movement to the list of new (
                        # sub-)movements.
                        these_airplane_movements.append(this_airplane_movement)

                    # Create a new movement for this iteration.
                    if is_static:
                        this_movement = movement.Movement(
                            airplane_movements=these_airplane_movements,
                            operating_point_movement=ref_operating_point_movement,
                            num_chords=wake_length,
                        )
                    else:
                        this_movement = movement.Movement(
                            airplane_movements=these_airplane_movements,
                            operating_point_movement=ref_operating_point_movement,
                            num_cycles=wake_length,
                        )

                    # Create a new problem for this iteration.
                    this_problem = problems.UnsteadyProblem(
                        movement=this_movement,
                        only_final_results=True,
                    )

                    # Create and run this iteration's solver and time how long it
                    # takes to execute.
                    this_solver = unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        unsteady_problem=this_problem
                    )
                    iter_start = time.time()
                    this_solver.run(
                        logging_level="Critical",
                        prescribed_wake=wake,
                        calculate_streamlines=False,
                    )
                    iter_stop = time.time()
                    this_iter_time = iter_stop - iter_start

                    # Create and fill arrays with each of this iteration's airplane's
                    # resultant force and moment coefficients.
                    these_coefficients = np.zeros(
                        (len(these_airplane_movements), num_coefficients)
                    )
                    for airplane_id, airplane in enumerate(these_airplane_movements):

                        # If this problem is static, then get it's final load
                        # coefficients. If it's variable, get the final RMS load
                        # coefficients.
                        if is_static:
                            all_force_coefficients = this_problem.final_near_field_force_coefficients_wind_axes[
                                airplane_id
                            ]
                            all_moment_coefficients = this_problem.final_near_field_moment_coefficients_wind_axes[
                                airplane_id
                            ]
                        else:
                            all_force_coefficients = this_problem.final_rms_near_field_force_coefficients_wind_axes[
                                airplane_id
                            ]
                            all_moment_coefficients = this_problem.final_rms_near_field_moment_coefficients_wind_axes[
                                airplane_id
                            ]

                        all_coefficients = np.concatenate(
                            [all_force_coefficients, all_moment_coefficients]
                        )

                        these_coefficients[airplane_id] = all_coefficients[
                            coefficient_mask
                        ]

                    # Populate the arrays that store information of all the
                    # iterations with the data from this iteration.
                    coefficients[wake_id, length_id, ar_id, chord_id, :, :] = (
                        these_coefficients
                    )
                    iter_times[wake_id, length_id, ar_id, chord_id] = this_iter_time

                    convergence_logger.info(
                        "\t\t\t\t\tIteration Time: "
                        + str(round(this_iter_time, 3))
                        + " s"
                    )

                    max_wake_pc = np.inf
                    max_length_pc = np.inf
                    max_ar_pc = np.inf
                    max_chord_pc = np.inf

                    # If this isn't the first wake state, calculate the wake state APE.
                    if wake_id > 0:
                        last_wake_coefficients = coefficients[
                            wake_id - 1, length_id, ar_id, chord_id, :, :
                        ]
                        max_wake_pc = np.max(
                            100
                            * np.abs(
                                (these_coefficients - last_wake_coefficients)
                                / last_wake_coefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from wake type: "
                            + str(round(max_wake_pc, 2))
                            + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from wake type: "
                            + str(max_wake_pc)
                        )

                    # If this isn't the first wake length, calculate the wake state APE.
                    if length_id > 0:
                        last_length_coefficients = coefficients[
                            wake_id, length_id - 1, ar_id, chord_id, :, :
                        ]
                        max_length_pc = np.max(
                            100
                            * np.abs(
                                (these_coefficients - last_length_coefficients)
                                / last_length_coefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from wake length: "
                            + str(round(max_length_pc, 2))
                            + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from wake length: "
                            + str(max_length_pc)
                        )

                    # If this isn't the first panel aspect ratio, calculate the panel
                    # aspect ratio APE.
                    if ar_id > 0:
                        last_ar_coefficients = coefficients[
                            wake_id, length_id, ar_id - 1, chord_id, :, :
                        ]
                        max_ar_pc = np.max(
                            100
                            * np.abs(
                                (these_coefficients - last_ar_coefficients)
                                / last_ar_coefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from panel aspect "
                            "ratio: " + str(round(max_ar_pc, 2)) + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from panel aspect "
                            "ratio: " + str(max_ar_pc)
                        )

                    # If this isn't the first number of chordwise panels, calculate
                    # the number of chordwise panels APE.
                    if chord_id > 0:
                        last_chord_coefficients = coefficients[
                            wake_id, length_id, ar_id, chord_id - 1, :, :
                        ]
                        max_chord_pc = np.max(
                            100
                            * np.abs(
                                (these_coefficients - last_chord_coefficients)
                                / last_chord_coefficients
                            )
                        )

                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from chordwise panels: "
                            + str(round(max_chord_pc, 2))
                            + "%"
                        )
                    else:
                        convergence_logger.info(
                            "\t\t\t\t\tMaximum coefficient change from chordwise panels: "
                            + str(max_chord_pc)
                        )

                    # Consider the panel aspect ratio value to be saturated if it is
                    # equal to 1. This is because a panel aspect ratio of 1 is
                    # considered the maximum degree of fineness. Consider the wake
                    # state to be saturated if it False (which corresponds to a free
                    # wake), as this is considered to be the most accurate wake state.
                    wake_saturated = not wake
                    ar_saturated = panel_aspect_ratio == 1

                    # Check if the user only specified one value for any of the four
                    # convergence parameters.
                    single_wake = len(wake_list) == 1
                    single_length = len(wake_lengths_list) == 1
                    single_ar = len(panel_aspect_ratios_list) == 1
                    single_chord = len(num_chordwise_panels_list) == 1

                    # Check if the iteration calculated that it is converged with
                    # respect to any of the four convergence parameters.
                    wake_converged = max_wake_pc < convergence_criteria
                    length_converged = max_length_pc < convergence_criteria
                    ar_converged = max_ar_pc < convergence_criteria
                    chord_converged = max_chord_pc < convergence_criteria

                    # Consider each convergence parameter to have passed it is
                    # converged, single, or saturated.
                    wake_passed = wake_converged or single_wake or wake_saturated
                    length_passed = length_converged or single_length
                    ar_passed = ar_converged or single_ar or ar_saturated
                    chord_passed = chord_converged or single_chord

                    # If all four convergence parameters have passed, then the solver
                    # has found a converged or semi-converged value and will return
                    # the converged parameters.
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
                            # We've tested more than one panel aspect ratio.
                            if ar_converged:
                                # There is no big difference between this panel aspect
                                # ratio and the last (coarser) panel aspect ratio.
                                # Therefore, the last (coarser) panel aspect ratio is
                                # converged.
                                converged_ar_id = ar_id - 1
                            else:
                                # There is a big difference between this panel aspect
                                # ratio and the last (coarser) panel aspect ratio.
                                # However, the panel aspect ratio is one, so it's
                                # saturated. Therefore, this panel aspect ratio is
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
                        converged_iter_time = iter_times[
                            converged_wake_id,
                            converged_length_id,
                            converged_ar_id,
                            converged_chord_id,
                        ]

                        if single_wake or single_length or single_ar or single_chord:
                            convergence_logger.info(
                                "The analysis found a semi-converged mesh:"
                            )
                            if single_wake:
                                convergence_logger.warning(
                                    "Wake type convergence not checked."
                                )
                            if single_length:
                                convergence_logger.warning(
                                    "Wake length convergence not checked."
                                )
                            if single_ar:
                                convergence_logger.warning(
                                    "Panel aspect ratio convergence not checked."
                                )
                            if single_chord:
                                convergence_logger.warning(
                                    "Chordwise panels convergence not checked."
                                )
                        else:
                            convergence_logger.info(
                                "The analysis found a converged mesh:"
                            )

                        if converged_wake:
                            convergence_logger.info("\tWake type: prescribed")
                        else:
                            convergence_logger.info("\tWake type: free")

                        if is_static:
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
                            "\tChordwise panels: " + str(converged_chordwise_panels)
                        )
                        convergence_logger.info(
                            "\tIteration time: "
                            + str(round(converged_iter_time, 3))
                            + " s"
                        )

                        return [
                            converged_wake,
                            converged_wake_length,
                            converged_aspect_ratio,
                            converged_chordwise_panels,
                        ]

    # If all iterations have been checked and none of them resulted in all
    # convergence parameters passing, then indicate that no converged solution was
    # found and return values of None for the converged parameters.
    convergence_logger.info("The analysis did not find a converged mesh.")
    return [None, None, None, None]
