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
import time

import numpy as np

from . import geometry
from . import problems
from . import steady_horseshoe_vortex_lattice_method
from . import output
from . import functions


# ToDo: Document this function.
def analyze_steady_convergence(
    base_problem,
    panel_aspect_ratio_bounds=(4, 1),
    num_chordwise_panels_bounds=(5, 14),
    convergence_criteria=0.5,
    logging_level="Debug",
):
    convergence_logger = logging.getLogger("convergence")

    logging_level_value = functions.convert_logging_level_name_to_value(logging_level)
    convergence_logger.setLevel(logging_level_value)

    logging.basicConfig()

    convergence_logger.info("Beginning convergence analysis.")

    base_operating_point = base_problem.operating_point

    base_airplane = base_problem.airplanes[0]
    base_wing = base_airplane.wings[0]
    base_wing_cross_sections = base_wing.wing_cross_sections

    if len(base_problem.airplanes) != 1:
        err_msg = "The problem for convergence analyses must have only one airplane."
        convergence_logger.error(msg=err_msg)

    if len(base_airplane.wings) != 1:
        err_msg = "The airplane for convergence analyses must have only one wing."
        convergence_logger.error(msg=err_msg)

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
        (len(panel_aspect_ratios_list), len(num_chordwise_panels_list))
    )
    moment_coefficients = np.zeros(
        (len(panel_aspect_ratios_list), len(num_chordwise_panels_list))
    )

    iteration = 0
    num_iterations = len(panel_aspect_ratios_list) * len(num_chordwise_panels_list)

    for ar_id, panel_aspect_ratio in enumerate(panel_aspect_ratios_list):

        ar_msg = "Panel aspect ratio: " + str(panel_aspect_ratio)
        convergence_logger.debug(msg=ar_msg)

        for chord_id, num_chordwise_panels in enumerate(num_chordwise_panels_list):

            num_spanwise_panels = round(
                (base_airplane.b_ref * num_chordwise_panels)
                / (base_airplane.c_ref * panel_aspect_ratio)
            )

            chordwise_msg = (
                "\tNumber of chordwise panels: "
                + str(num_chordwise_panels)
                + " (Number of spanwise panels: "
                + str(num_spanwise_panels)
                + ")"
            )
            convergence_logger.debug(msg=chordwise_msg)

            iteration += 1
            iteration_msg = (
                "\t\tIteration Number: " + str(iteration) + "/" + str(num_iterations)
            )
            convergence_logger.debug(msg=iteration_msg)

            these_wing_cross_sections = []

            for base_wing_cross_section in base_wing_cross_sections:
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
                        num_spanwise_panels=num_spanwise_panels,
                        airfoil=geometry.Airfoil(
                            name=base_wing_cross_section.airfoil.name,
                            coordinates=base_wing_cross_section.airfoil.coordinates,
                            repanel=base_wing_cross_section.airfoil.repanel,
                            n_points_per_side=base_wing_cross_section.airfoil.n_points_per_side,
                        ),
                    )
                )

            this_airplane = geometry.Airplane(
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
                wings=[
                    geometry.Wing(
                        # These values are copied from the base wing.
                        name=base_wing.name,
                        x_le=base_wing.x_le,
                        y_le=base_wing.y_le,
                        z_le=base_wing.z_le,
                        symmetric=base_wing.symmetric,
                        chordwise_spacing=base_wing.chordwise_spacing,
                        # These values change.
                        num_chordwise_panels=num_chordwise_panels,
                        wing_cross_sections=these_wing_cross_sections,
                    ),
                ],
            )

            this_problem = problems.SteadyProblem(
                airplanes=[this_airplane], operating_point=base_operating_point
            )

            this_solver = steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                steady_problem=this_problem,
            )

            del this_problem

            iter_start = time.time()
            this_solver.run(logging_level="Critical")
            iter_stop = time.time()

            this_iter_time = iter_stop - iter_start
            this_force_coefficient = np.linalg.norm(
                this_airplane.total_near_field_force_coefficients_wind_axes
            )
            this_moment_coefficient = np.linalg.norm(
                this_airplane.total_near_field_moment_coefficients_wind_axes
            )

            iter_times[ar_id, chord_id] = this_iter_time
            force_coefficients[ar_id, chord_id] = this_force_coefficient
            moment_coefficients[ar_id, chord_id] = this_moment_coefficient

            time_msg = "\t\tIteration Time: " + str(round(this_iter_time, 3)) + " s"
            convergence_logger.debug(msg=time_msg)

            max_ar_pc = np.inf
            max_chord_pc = np.inf

            if ar_id > 0:
                last_ar_force_coefficient = force_coefficients[ar_id - 1, chord_id]
                last_ar_moment_coefficient = moment_coefficients[ar_id - 1, chord_id]
                ar_force_pc = 100 * np.abs(
                    (this_force_coefficient - last_ar_force_coefficient)
                    / last_ar_force_coefficient
                )
                ar_moment_pc = 100 * np.abs(
                    (this_moment_coefficient - last_ar_moment_coefficient)
                    / last_ar_moment_coefficient
                )
                max_ar_pc = max(ar_force_pc, ar_moment_pc)

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
                last_chord_force_coefficient = force_coefficients[ar_id, chord_id - 1]
                last_chord_moment_coefficient = moment_coefficients[ar_id, chord_id - 1]
                chord_force_pc = 100 * np.abs(
                    (this_force_coefficient - last_chord_force_coefficient)
                    / last_chord_force_coefficient
                )
                chord_moment_pc = 100 * np.abs(
                    (this_moment_coefficient - last_chord_moment_coefficient)
                    / last_chord_moment_coefficient
                )
                max_chord_pc = max(chord_force_pc, chord_moment_pc)

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
                converged_chordwise_panels = num_chordwise_panels_list[chord_id - 1]
                converged_aspect_ratio = panel_aspect_ratios_list[ar_id - 1]
                converged_spanwise_panels = round(
                    (base_airplane.b_ref * converged_chordwise_panels)
                    / (base_airplane.c_ref * converged_aspect_ratio)
                )
                converged_iter_time = iter_times[ar_id - 1, chord_id - 1]

                convergence_logger.info("The analysis found a converged mesh:")
                convergence_logger.info(
                    "\tConverged panel aspect ratio: " + str(converged_aspect_ratio)
                )
                convergence_logger.info(
                    "\tConverged number of chordwise panels: "
                    + str(converged_chordwise_panels)
                    + " (Converged number of spanwise panels: "
                    + str(converged_spanwise_panels)
                    + ")"
                )
                convergence_logger.info(
                    "\tConverged iteration time: "
                    + str(round(converged_iter_time, 3))
                    + " s"
                )

                if convergence_logger.level == logging.DEBUG:
                    output.draw(
                        solver=this_solver,
                        scalar_type="lift",
                    )
                return [
                    converged_chordwise_panels,
                    converged_aspect_ratio,
                    converged_spanwise_panels,
                ]

    convergence_logger.info("The analysis did not find a converged mesh.")
    return [None, None, None]
