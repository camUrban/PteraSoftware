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

    base_airplanes = base_problem.airplanes
    # base_airplane = base_problem.airplanes[0]
    # base_wings = base_airplane.wings

    if len(base_problem.airplanes) != 1:
        err_msg = "The problem for convergence analyses must have only one airplane."
        convergence_logger.critical(msg=err_msg)

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

                this_num_spanwise_panels = round(
                    (base_airplane.b_ref * num_chordwise_panels)
                    / (base_airplane.c_ref * panel_aspect_ratio)
                )

                base_wings = base_airplane.wings
                these_wings = []
                for base_wing in base_wings:

                    base_wing_cross_sections = base_wing.wing_cross_sections
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
                converged_chordwise_panels = num_chordwise_panels_list[chord_id - 1]
                converged_aspect_ratio = panel_aspect_ratios_list[ar_id - 1]
                converged_iter_time = iter_times[ar_id - 1, chord_id - 1]

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

                if convergence_logger.level == logging.DEBUG:
                    output.draw(
                        solver=this_solver,
                        scalar_type="lift",
                    )

                converged_spanwise_panels = []
                for base_airplane in base_airplanes:
                    converged_spanwise_panels.append(
                        round(
                            (base_airplane.b_ref * converged_chordwise_panels)
                            / (base_airplane.c_ref * converged_aspect_ratio)
                        )
                    )

                return [
                    converged_chordwise_panels,
                    converged_aspect_ratio,
                    converged_spanwise_panels,
                ]

    convergence_logger.info("The analysis did not find a converged mesh.")
    return [None, None, None]
