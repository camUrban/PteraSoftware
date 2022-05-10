# ToDo: Document this module.
"""This module contains functions to analyze the convergence of steady and unsteady
solvers.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    analyze_steady_convergence:

    analyze_unsteady_convergence:"""
import logging
import math
import time

import numpy as np

from . import geometry
from . import problems
from . import movement

from . import unsteady_ring_vortex_lattice_method
from . import steady_horseshoe_vortex_lattice_method
from . import output
from . import functions

convergence_logger = logging.getLogger("convergence")


# ToDo: Document this function.
def analyze_steady_convergence(
    base_problem,
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(3, 12),
    convergence_criteria=1.0,
    logging_level="Debug",
):
    logging_level_value = functions.convert_logging_level_name_to_value(logging_level)
    convergence_logger.setLevel(logging_level_value)

    logging.basicConfig()

    convergence_logger.info("Beginning convergence analysis.")

    base_operating_point = base_problem.operating_point

    base_airplanes = base_problem.airplanes

    panel_aspect_ratios_list = list(
        range(panel_aspect_ratio_bounds[0], panel_aspect_ratio_bounds[1] - 1, -1)
    )

    num_chordwise_panels_list = list(
        range(num_chordwise_panels_bounds[0], num_chordwise_panels_bounds[1] + 1)
    )

    iter_times = np.zeros(
        (len(panel_aspect_ratios_list), len(num_chordwise_panels_list))
    )
    force_coefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(base_airplanes),
        )
    )
    moment_coefficients = np.zeros(
        (
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(base_airplanes),
        )
    )

    iteration = 0
    num_iterations = len(panel_aspect_ratios_list) * len(num_chordwise_panels_list)

    for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):

        ar_msg = "Panel aspect ratio: " + str(panel_aspect_ratio)
        convergence_logger.debug(msg=ar_msg)

        for chord_id, num_chordwise_panels in enumerate(num_chordwise_panels_list):

            chordwise_msg = "\tNumber of chordwise panels: " + str(num_chordwise_panels)
            convergence_logger.debug(msg=chordwise_msg)

            iteration += 1
            iteration_msg = (
                "\t\tIteration Number: " + str(iteration) + "/" + str(num_iterations)
            )
            convergence_logger.debug(msg=iteration_msg)

            these_airplanes = []
            for base_airplane in base_airplanes:

                base_wings = base_airplane.wings
                these_wings = []
                for base_wing in base_wings:

                    base_wing_cross_sections = base_wing.wing_cross_sections
                    these_wing_cross_sections = []
                    for (
                        base_wing_cross_section_id,
                        base_wing_cross_section,
                    ) in enumerate(base_wing_cross_sections):

                        if base_wing_cross_section_id < (
                            len(base_wing_cross_sections) - 1
                        ):
                            next_base_wing_cross_section = base_wing_cross_sections[
                                base_wing_cross_section_id + 1
                            ]
                            section_length = (
                                next_base_wing_cross_section.y_le
                                - base_wing_cross_section.y_le
                            )
                            root_chord = base_wing_cross_section.chord
                            tip_chord = next_base_wing_cross_section.chord
                            section_area = section_length * (root_chord + tip_chord) / 2
                            section_standard_mean_chord = section_area / section_length

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
                            geometry.WingCrossSection(
                                # These values are copied from the base wing cross section.
                                x_le=base_wing_cross_section.x_le,
                                y_le=base_wing_cross_section.y_le,
                                z_le=base_wing_cross_section.z_le,
                                chord=base_wing_cross_section.chord,
                                twist=base_wing_cross_section.twist,
                                control_surface_type=base_wing_cross_section.control_surface_type,
                                control_surface_hinge_point=base_wing_cross_section.control_surface_hinge_point,
                                control_surface_deflection=base_wing_cross_section.control_surface_deflection,
                                spanwise_spacing=base_wing_cross_section.spanwise_spacing,
                                # These values change.
                                num_spanwise_panels=this_num_spanwise_panels,
                                airfoil=geometry.Airfoil(
                                    name=base_wing_cross_section.airfoil.name,
                                    coordinates=base_wing_cross_section.airfoil.coordinates,
                                    repanel=base_wing_cross_section.airfoil.repanel,
                                    n_points_per_side=base_wing_cross_section.airfoil.n_points_per_side,
                                ),
                            )
                        )

                    these_wings.append(
                        geometry.Wing(
                            # These values are copied from this base wing.
                            name=base_wing.name,
                            x_le=base_wing.x_le,
                            y_le=base_wing.y_le,
                            z_le=base_wing.z_le,
                            symmetric=base_wing.symmetric,
                            chordwise_spacing=base_wing.chordwise_spacing,
                            # These values change.
                            num_chordwise_panels=num_chordwise_panels,
                            wing_cross_sections=these_wing_cross_sections,
                        )
                    )

                these_airplanes.append(
                    geometry.Airplane(
                        # These values are copied from the base airplane.
                        name=base_airplane.name,
                        x_ref=base_airplane.x_ref,
                        y_ref=base_airplane.y_ref,
                        z_ref=base_airplane.z_ref,
                        weight=base_airplane.weight,
                        # These are kept as None so that they are recalculated with this
                        # airplane's mesh.
                        s_ref=None,
                        c_ref=None,
                        b_ref=None,
                        # This value changes.
                        wings=these_wings,
                    )
                )

            this_problem = problems.SteadyProblem(
                airplanes=these_airplanes, operating_point=base_operating_point
            )

            # ToDo: Have this function be capable of running either steady solver (
            #  not just the horseshoe vortex solver).
            this_solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                steady_problem=this_problem,
            )

            del this_problem

            iter_start = time.time()
            this_solver.run(logging_level="Critical")
            iter_stop = time.time()

            this_iter_time = iter_stop - iter_start

            these_force_coefficients = np.zeros(len(these_airplanes))
            these_moment_coefficients = np.zeros(len(these_airplanes))
            for airplane_id, airplane in enumerate(these_airplanes):
                these_force_coefficients[airplane_id] = np.linalg.norm(
                    airplane.total_near_field_force_coefficients_wind_axes
                )
                these_moment_coefficients[airplane_id] = np.linalg.norm(
                    airplane.total_near_field_moment_coefficients_wind_axes
                )

            force_coefficients[ar_id, chord_id, :] = these_force_coefficients
            moment_coefficients[ar_id, chord_id, :] = these_moment_coefficients

            iter_times[ar_id, chord_id] = this_iter_time

            time_msg = "\t\tIteration Time: " + str(round(this_iter_time, 3)) + " s"
            convergence_logger.debug(msg=time_msg)

            max_ar_pc = np.inf
            max_chord_pc = np.inf

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

                max_ar_pc_msg = (
                    "\t\tMaximum coefficient change from the panel aspect ratio: "
                    + str(round(max_ar_pc, 2))
                    + "%"
                )
                convergence_logger.debug(msg=max_ar_pc_msg)
            else:
                max_ar_pc_msg = (
                    "\t\tMaximum coefficient change from the panel aspect ratio: "
                    + str(max_ar_pc)
                )
                convergence_logger.debug(msg=max_ar_pc_msg)

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

                max_chord_pc_msg = (
                    "\t\tMaximum coefficient change from the number of chordwise panels: "
                    + str(round(max_chord_pc, 2))
                    + "%"
                )
                convergence_logger.debug(msg=max_chord_pc_msg)
            else:
                max_chord_pc_msg = (
                    "\t\tMaximum coefficient change from the number of chordwise panels: "
                    + str(max_chord_pc)
                )
                convergence_logger.debug(msg=max_chord_pc_msg)

            single_ar = len(panel_aspect_ratios_list) == 1
            single_chord = len(num_chordwise_panels_list) == 1

            ar_converged = max_ar_pc < convergence_criteria
            chord_converged = max_chord_pc < convergence_criteria

            ar_passed = ar_converged or single_ar
            chord_passed = chord_converged or single_chord

            if ar_passed and chord_passed:
                if single_ar:
                    converged_ar_id = ar_id
                else:
                    converged_ar_id = ar_id - 1
                if single_chord:
                    converged_chord_id = chord_id
                else:
                    converged_chord_id = chord_id - 1

                converged_chordwise_panels = num_chordwise_panels_list[
                    converged_chord_id
                ]
                converged_aspect_ratio = panel_aspect_ratios_list[converged_ar_id]
                converged_iter_time = iter_times[converged_ar_id, converged_chord_id]

                convergence_logger.info("The analysis found a converged mesh:")
                convergence_logger.info(
                    "\tConverged panel aspect ratio: " + str(converged_aspect_ratio)
                )
                convergence_logger.info(
                    "\tConverged number of chordwise panels: "
                    + str(converged_chordwise_panels)
                )
                convergence_logger.info(
                    "\tConverged iteration time: "
                    + str(round(converged_iter_time, 3))
                    + " s"
                )

                # ToDo: BUG: This is displaying the "super-converged" solver.
                if convergence_logger.level == logging.DEBUG:
                    output.draw(
                        solver=this_solver,
                        scalar_type="lift",
                    )

                return [
                    converged_chordwise_panels,
                    converged_aspect_ratio,
                ]

    convergence_logger.info("The analysis did not find a converged mesh.")
    return [None, None]


# ToDo: Add the ability for this function to deal with static geometry unsteady
#  problems.
# ToDo: Document this function.
def analyze_unsteady_convergence(
    ref_airplane_movements,
    ref_operating_point_movement,
    prescribed_wake=True,
    free_wake=False,
    num_cycles_bounds=(1, 1),
    panel_aspect_ratio_bounds=(2, 1),
    num_chordwise_panels_bounds=(3, 5),
    convergence_criteria=1.0,
    logging_level="Debug",
):
    logging_level_value = functions.convert_logging_level_name_to_value(logging_level)
    convergence_logger.setLevel(logging_level_value)

    logging.basicConfig()

    convergence_logger.info("Beginning convergence analysis.")

    wake_list = []
    if prescribed_wake:
        wake_list.append(True)
    if free_wake:
        wake_list.append(False)
    if not (free_wake or prescribed_wake):
        convergence_logger.critical(
            "Your solver must have at least one type of wake equal to True."
        )

    num_cycles_list = list(range(num_cycles_bounds[0], num_cycles_bounds[1] + 1))

    panel_aspect_ratios_list = list(
        range(panel_aspect_ratio_bounds[0], panel_aspect_ratio_bounds[1] - 1, -1)
    )

    num_chordwise_panels_list = list(
        range(num_chordwise_panels_bounds[0], num_chordwise_panels_bounds[1] + 1)
    )

    iter_times = np.zeros(
        (
            len(wake_list),
            len(num_cycles_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
        )
    )

    force_coefficients = np.zeros(
        (
            len(wake_list),
            len(num_cycles_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplane_movements),
        )
    )
    moment_coefficients = np.zeros(
        (
            len(wake_list),
            len(num_cycles_list),
            len(panel_aspect_ratios_list),
            len(num_chordwise_panels_list),
            len(ref_airplane_movements),
        )
    )

    iteration = 0
    num_iterations = (
        len(wake_list)
        * len(num_cycles_list)
        * len(panel_aspect_ratios_list)
        * len(num_chordwise_panels_list)
    )

    for wake_id, wake in enumerate(wake_list):
        wake_msg = "Wake type: "
        if wake:
            wake_msg += "prescribed"
        else:
            wake_msg += "free"
        convergence_logger.debug(msg=wake_msg)

        for cycle_id, num_cycles in enumerate(num_cycles_list):
            cycle_msg = "\tNumber of cycles: " + str(num_cycles)
            convergence_logger.debug(msg=cycle_msg)

            for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):

                ar_msg = "\t\tPanel aspect ratio: " + str(panel_aspect_ratio)
                convergence_logger.debug(msg=ar_msg)

                for chord_id, num_chordwise_panels in enumerate(
                    num_chordwise_panels_list
                ):

                    chordwise_msg = "\t\t\tNumber of chordwise panels: " + str(
                        num_chordwise_panels
                    )
                    convergence_logger.debug(msg=chordwise_msg)

                    iteration += 1
                    iteration_msg = (
                        "\t\t\t\tIteration Number: "
                        + str(iteration)
                        + "/"
                        + str(num_iterations)
                    )
                    convergence_logger.debug(msg=iteration_msg)

                    these_base_airplanes = []
                    these_airplane_movements = []
                    for ref_airplane_movement in ref_airplane_movements:
                        # 1: Reference this movement's base.
                        ref_base_airplane = ref_airplane_movement.base_airplane

                        # 2: Reference this movement's list of sub-movements.
                        ref_wing_movements = ref_airplane_movement.wing_movements

                        # 3: Create an empty list for the sub-movement base objects.
                        these_base_wings = []

                        # 4: Create an empty list for the sub-movement copies.
                        these_wing_movements = []

                        # 5: Iterate over the sub-movements.
                        for ref_wing_movement in ref_wing_movements:
                            # 1: Reference this movement's base.
                            ref_base_wing = ref_wing_movement.base_wing

                            # 2: Reference this movement's list of sub-movements.
                            ref_wing_cross_section_movements = (
                                ref_wing_movement.wing_cross_section_movements
                            )

                            # 3: Create an empty list for the sub-movement base objects.
                            these_base_wing_cross_sections = []

                            # 4: Create an empty list for the sub-movement copies.
                            these_wing_cross_section_movements = []

                            # 5: Iterate over the sub-movements.
                            for (
                                ref_wing_cross_section_movement_id,
                                ref_wing_cross_section_movement,
                            ) in enumerate(ref_wing_cross_section_movements):
                                # 1: Reference this movement's base.
                                ref_base_wing_cross_section = (
                                    ref_wing_cross_section_movement.base_wing_cross_section
                                )

                                if ref_wing_cross_section_movement_id < (
                                    len(ref_wing_cross_section_movements) - 1
                                ):
                                    next_ref_base_wing_cross_section = (
                                        ref_wing_cross_section_movements[
                                            ref_wing_cross_section_movement_id + 1
                                        ].base_wing_cross_section
                                    )
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

                                # 2: Reference this movement's list of sub-movements.
                                # N/A

                                # 3: Create an empty list for the sub-movement base
                                # objects.
                                # N/A

                                # 4: Create an empty list for the sub-movement copies.
                                # N/A

                                # 5: Iterate over the sub-movements.
                                # N/A

                                # 6: Create a copy of the base.
                                this_base_wing_cross_section = geometry.WingCrossSection(
                                    # These values are copied from the reference
                                    # base wing cross section.
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
                                    airfoil=geometry.Airfoil(
                                        name=ref_base_wing_cross_section.airfoil.name,
                                        coordinates=ref_base_wing_cross_section.airfoil.coordinates,
                                        repanel=ref_base_wing_cross_section.airfoil.repanel,
                                        n_points_per_side=ref_base_wing_cross_section.airfoil.n_points_per_side,
                                    ),
                                )

                                # 7. Create a copy of the new movement.
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

                                # 9. Append the new movement to the list of new
                                # movements.
                                these_wing_cross_section_movements.append(
                                    this_wing_cross_section_movement
                                )

                            # 6: Create a copy of the base.
                            this_base_wing = geometry.Wing(
                                # These values are copied from this reference base wing.
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

                            # 7. Create a copy of the new movement.
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

                            # 9. Append the new movement to the list of new movements.
                            these_wing_movements.append(this_wing_movement)

                        # 6: Create a copy of the base.
                        this_base_airplane = geometry.Airplane(
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
                            b_ref=None,
                            # This value changes.
                            wings=these_base_wings,
                        )

                        # 7. Create a copy of the new movement.
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

                        # 9. Append the new movement to the list of new movements.
                        these_airplane_movements.append(this_airplane_movement)

                    this_movement = movement.Movement(
                        airplane_movements=these_airplane_movements,
                        operating_point_movement=ref_operating_point_movement,
                        num_cycles=num_cycles,
                    )

                    this_problem = problems.UnsteadyProblem(
                        movement=this_movement,
                        only_final_results=True,
                    )

                    this_solver = unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        unsteady_problem=this_problem
                    )

                    iter_start = time.time()
                    this_solver.run(
                        logging_level="Warning",
                        prescribed_wake=wake,
                        calculate_streamlines=False,
                    )
                    iter_stop = time.time()

                    this_iter_time = iter_stop - iter_start

                    these_force_coefficients = np.zeros(len(these_airplane_movements))
                    these_moment_coefficients = np.zeros(len(these_airplane_movements))
                    for airplane_id, airplane in enumerate(these_airplane_movements):
                        these_force_coefficients[airplane_id] = np.linalg.norm(
                            this_problem.final_total_near_field_force_coefficients_wind_axes[
                                airplane_id
                            ]
                        )
                        these_moment_coefficients[airplane_id] = np.linalg.norm(
                            this_problem.final_total_near_field_moment_coefficients_wind_axes[
                                airplane_id
                            ]
                        )

                    force_coefficients[
                        wake_id, cycle_id, ar_id, chord_id, :
                    ] = these_force_coefficients
                    moment_coefficients[
                        wake_id, cycle_id, ar_id, chord_id, :
                    ] = these_moment_coefficients

                    iter_times[wake_id, cycle_id, ar_id, chord_id] = this_iter_time

                    time_msg = (
                        "\t\t\t\tIteration Time: "
                        + str(round(this_iter_time, 3))
                        + " s"
                    )
                    convergence_logger.debug(msg=time_msg)

                    max_wake_pc = np.inf
                    max_cycle_pc = np.inf
                    max_ar_pc = np.inf
                    max_chord_pc = np.inf

                    if wake_id > 0:
                        last_wake_force_coefficients = force_coefficients[
                            wake_id - 1, cycle_id, ar_id, chord_id, :
                        ]
                        last_wake_moment_coefficients = moment_coefficients[
                            wake_id - 1, cycle_id, ar_id, chord_id, :
                        ]
                        max_wake_force_pc = max(
                            100
                            * np.abs(
                                (
                                    these_force_coefficients
                                    - last_wake_force_coefficients
                                )
                                / last_wake_force_coefficients
                            )
                        )
                        max_wake_moment_pc = max(
                            100
                            * np.abs(
                                (
                                    these_moment_coefficients
                                    - last_wake_moment_coefficients
                                )
                                / last_wake_moment_coefficients
                            )
                        )
                        max_wake_pc = max(max_wake_force_pc, max_wake_moment_pc)

                        max_wake_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the wake type: "
                            + str(round(max_wake_pc, 2))
                            + "%"
                        )
                        convergence_logger.debug(msg=max_wake_pc_msg)
                    else:
                        max_wake_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the wake type: "
                            + str(max_wake_pc)
                        )
                        convergence_logger.debug(msg=max_wake_pc_msg)

                    if cycle_id > 0:
                        last_cycle_force_coefficients = force_coefficients[
                            wake_id, cycle_id - 1, ar_id, chord_id, :
                        ]
                        last_cycle_moment_coefficients = moment_coefficients[
                            wake_id, cycle_id - 1, ar_id, chord_id, :
                        ]
                        max_cycle_force_pc = max(
                            100
                            * np.abs(
                                (
                                    these_force_coefficients
                                    - last_cycle_force_coefficients
                                )
                                / last_cycle_force_coefficients
                            )
                        )
                        max_cycle_moment_pc = max(
                            100
                            * np.abs(
                                (
                                    these_moment_coefficients
                                    - last_cycle_moment_coefficients
                                )
                                / last_cycle_moment_coefficients
                            )
                        )
                        max_cycle_pc = max(max_cycle_force_pc, max_cycle_moment_pc)

                        max_cycle_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the number of cycles: "
                            + str(round(max_cycle_pc, 2))
                            + "%"
                        )
                        convergence_logger.debug(msg=max_cycle_pc_msg)
                    else:
                        max_cycle_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the number of cycles: "
                            + str(max_cycle_pc)
                        )
                        convergence_logger.debug(msg=max_cycle_pc_msg)

                    if ar_id > 0:
                        last_ar_force_coefficients = force_coefficients[
                            wake_id, cycle_id, ar_id - 1, chord_id, :
                        ]
                        last_ar_moment_coefficients = moment_coefficients[
                            wake_id, cycle_id, ar_id - 1, chord_id, :
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
                                (
                                    these_moment_coefficients
                                    - last_ar_moment_coefficients
                                )
                                / last_ar_moment_coefficients
                            )
                        )
                        max_ar_pc = max(max_ar_force_pc, max_ar_moment_pc)

                        max_ar_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the panel aspect ratio: "
                            + str(round(max_ar_pc, 2))
                            + "%"
                        )
                        convergence_logger.debug(msg=max_ar_pc_msg)
                    else:
                        max_ar_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the panel aspect ratio: "
                            + str(max_ar_pc)
                        )
                        convergence_logger.debug(msg=max_ar_pc_msg)

                    if chord_id > 0:
                        last_chord_force_coefficients = force_coefficients[
                            wake_id, cycle_id, ar_id, chord_id - 1, :
                        ]
                        last_chord_moment_coefficients = moment_coefficients[
                            wake_id, cycle_id, ar_id, chord_id - 1, :
                        ]
                        max_chord_force_pc = max(
                            100
                            * np.abs(
                                (
                                    these_force_coefficients
                                    - last_chord_force_coefficients
                                )
                                / last_chord_force_coefficients
                            )
                        )
                        max_chord_moment_pc = max(
                            100
                            * np.abs(
                                (
                                    these_moment_coefficients
                                    - last_chord_moment_coefficients
                                )
                                / last_chord_moment_coefficients
                            )
                        )
                        max_chord_pc = max(max_chord_force_pc, max_chord_moment_pc)

                        max_chord_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the number of chordwise panels: "
                            + str(round(max_chord_pc, 2))
                            + "%"
                        )
                        convergence_logger.debug(msg=max_chord_pc_msg)
                    else:
                        max_chord_pc_msg = (
                            "\t\t\t\tMaximum coefficient change from the number of chordwise panels: "
                            + str(max_chord_pc)
                        )
                        convergence_logger.debug(msg=max_chord_pc_msg)

                    single_wake = len(wake_list) == 1
                    single_cycle = len(num_cycles_list) == 1
                    single_ar = len(panel_aspect_ratios_list) == 1
                    single_chord = len(num_chordwise_panels_list) == 1

                    wake_converged = max_wake_pc < convergence_criteria
                    cycle_converged = max_cycle_pc < convergence_criteria
                    ar_converged = max_ar_pc < convergence_criteria
                    chord_converged = max_chord_pc < convergence_criteria

                    wake_saturated = not wake
                    ar_saturated = panel_aspect_ratio == 1

                    wake_passed = wake_converged or single_wake or wake_saturated
                    cycle_passed = cycle_converged or single_cycle
                    ar_passed = ar_converged or single_ar or ar_saturated
                    chord_passed = chord_converged or single_chord

                    if wake_passed and cycle_passed and ar_passed and chord_passed:
                        if single_wake or (not wake_converged):
                            converged_wake_id = wake_id
                        else:
                            converged_wake_id = wake_id - 1
                        if single_cycle:
                            converged_cycle_id = cycle_id
                        else:
                            converged_cycle_id = cycle_id - 1
                        if single_ar or (not ar_converged):
                            converged_ar_id = ar_id
                        else:
                            converged_ar_id = ar_id - 1
                        if single_chord:
                            converged_chord_id = chord_id
                        else:
                            converged_chord_id = chord_id - 1

                        converged_wake = wake_list[converged_wake_id]
                        converged_num_cycles = num_cycles_list[converged_cycle_id]
                        converged_chordwise_panels = num_chordwise_panels_list[
                            converged_chord_id
                        ]
                        converged_aspect_ratio = panel_aspect_ratios_list[
                            converged_ar_id
                        ]
                        converged_iter_time = iter_times[
                            converged_wake_id,
                            converged_cycle_id,
                            converged_ar_id,
                            converged_chord_id,
                        ]

                        convergence_logger.info("The analysis found a converged mesh:")
                        convergence_logger.info(
                            "\tConverged wake type: " + str(converged_wake)
                        )
                        convergence_logger.info(
                            "\tConverged number of cycles: " + str(converged_num_cycles)
                        )
                        convergence_logger.info(
                            "\tConverged panel aspect ratio: "
                            + str(converged_aspect_ratio)
                        )
                        convergence_logger.info(
                            "\tConverged number of chordwise panels: "
                            + str(converged_chordwise_panels)
                        )
                        convergence_logger.info(
                            "\tConverged iteration time: "
                            + str(round(converged_iter_time, 3))
                            + " s"
                        )

                        # ToDo: BUG: This is displaying the "super-converged" solver.
                        if convergence_logger.level == logging.DEBUG:
                            output.draw(
                                solver=this_solver,
                                scalar_type="lift",
                                show_wake_vortices=True,
                            )

                        return [
                            converged_wake,
                            converged_num_cycles,
                            converged_chordwise_panels,
                            converged_aspect_ratio,
                        ]

    convergence_logger.info("The analysis did not find a converged mesh.")
    return [None, None, None, None]
