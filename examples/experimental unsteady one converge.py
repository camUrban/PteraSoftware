"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import time
import gc

import numpy as np

import src

start_time = time.time()

alpha = 5
x_spacing = 0.5
y_spacing = 0.5

min_num_flaps = 2
max_num_flaps = 4
min_num_chord = 10
max_num_chord = 15

wake_state_list = [True, False]
num_flaps_list = [i for i in range(min_num_flaps, max_num_flaps + 1)]
num_chord_list = [i for i in range(min_num_chord, max_num_chord + 1)]

del min_num_flaps
del max_num_flaps
del min_num_chord
del max_num_chord

all_mean_forces = np.zeros(
    (
        len(wake_state_list),
        len(num_flaps_list),
        len(num_chord_list),
        3,
    )
)
all_mean_moments = np.zeros(
    (
        len(wake_state_list),
        len(num_flaps_list),
        len(num_chord_list),
        3,
    )
)

iteration = 1
num_iterations = len(wake_state_list) * len(num_flaps_list) * len(num_chord_list)

this_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
    velocity=1.0,
    alpha=0.0,
    nu=15.06e-6,
)
this_operating_point_movement = src.movement.OperatingPointMovement(
    base_operating_point=this_operating_point,
)
del this_operating_point

root_to_mid_span = 0.2275
root_chord = 0.1094
mid_to_tip_span = 0.350 - 0.2275
tip_chord = 0.0219

root_to_mid_chord = root_chord
mid_to_tip_chord = (root_chord + tip_chord) / 2

for wake_state_id, wake_state in enumerate(wake_state_list):
    print("Prescribed Wake: ", wake_state)
    for num_flaps_id, num_flaps in enumerate(num_flaps_list):
        print("\tNumber of flaps: ", num_flaps)
        for num_chord_id, num_chord in enumerate(num_chord_list):
            print("\t\tNumber of chordwise panels: ", num_chord)
            print("\t\t\t\t", iteration, "/", num_iterations, sep="")

            root_to_mid_panel_chord = root_to_mid_chord / num_chord
            mid_to_tip_panel_chord = mid_to_tip_chord / num_chord

            root_to_mid_num_span = round(
                root_to_mid_span / (5 * root_to_mid_panel_chord)
            )
            mid_to_tip_num_span = round(mid_to_tip_span / (5 * mid_to_tip_panel_chord))

            lead_airplane = src.geometry.Airplane(
                name="Lead Airplane",
                wings=[
                    src.geometry.Wing(
                        name="Main Wing",
                        symmetric=True,
                        chordwise_spacing="uniform",
                        num_chordwise_panels=num_chord,
                        wing_cross_sections=[
                            src.geometry.WingCrossSection(
                                twist=alpha,
                                chord=0.1094,
                                airfoil=src.geometry.Airfoil(name="naca0012"),
                                num_spanwise_panels=root_to_mid_num_span,
                                spanwise_spacing="uniform",
                            ),
                            src.geometry.WingCrossSection(
                                twist=alpha,
                                y_le=0.2275,
                                chord=0.1094,
                                airfoil=src.geometry.Airfoil(name="naca0012"),
                                num_spanwise_panels=mid_to_tip_num_span,
                                spanwise_spacing="uniform",
                            ),
                            src.geometry.WingCrossSection(
                                twist=alpha,
                                y_le=0.350,
                                chord=0.0219,
                                airfoil=src.geometry.Airfoil(name="naca0012"),
                            ),
                        ],
                    ),
                ],
            )

            lead_airplane_movement = src.movement.AirplaneMovement(
                base_airplane=lead_airplane,
                wing_movements=[
                    src.movement.WingMovement(
                        base_wing=lead_airplane.wings[0],
                        wing_cross_sections_movements=[
                            src.movement.WingCrossSectionMovement(
                                base_wing_cross_section=lead_airplane.wings[
                                    0
                                ].wing_cross_sections[0],
                            ),
                            src.movement.WingCrossSectionMovement(
                                base_wing_cross_section=lead_airplane.wings[
                                    0
                                ].wing_cross_sections[1],
                                sweeping_amplitude=15.0,
                                sweeping_period=1,
                                sweeping_spacing="sine",
                            ),
                            src.movement.WingCrossSectionMovement(
                                base_wing_cross_section=lead_airplane.wings[
                                    0
                                ].wing_cross_sections[2],
                                sweeping_amplitude=15.0,
                                sweeping_period=1,
                                sweeping_spacing="sine",
                            ),
                        ],
                    )
                ],
            )

            del lead_airplane

            this_movement = src.movement.Movement(
                airplane_movements=[
                    lead_airplane_movement,
                ],
                operating_point_movement=this_operating_point_movement,
                num_steps=None,
                num_cycles=num_flaps,
                delta_time=None,
            )

            del lead_airplane_movement

            this_problem = src.problems.UnsteadyProblem(
                movement=this_movement,
                only_final_results=True,
            )

            del this_movement

            this_solver = src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                unsteady_problem=this_problem,
            )

            del this_problem

            iter_start = time.time()

            this_solver.run(
                logging_level="Critical",
                prescribed_wake=wake_state,
                calculate_streamlines=False,
            )

            iter_stop = time.time()
            iter_time = round((iter_stop - iter_start), 2)
            print("\t\t\t\t", iter_time, " s", sep="")

            first_results_step = this_solver.first_results_step
            num_steps = this_solver.num_steps
            num_results_steps = num_steps - first_results_step

            total_forces = np.zeros((3, num_results_steps))
            total_moments = np.zeros((3, num_results_steps))

            mean_forces = np.zeros(3)
            mean_moments = np.zeros(3)

            results_step = 0
            for step in range(first_results_step, num_steps):
                airplane = this_solver.steady_problems[step].airplanes[0]
                total_forces[
                    :, results_step
                ] = airplane.total_near_field_force_wind_axes
                total_moments[
                    :, results_step
                ] = airplane.total_near_field_moment_wind_axes
                results_step += 1

            del results_step
            del step
            del airplane

            this_induced_drag = np.mean(total_forces[0, :])
            this_side_force = np.mean(total_forces[1, :])
            this_lift = np.mean(total_forces[2, :])
            this_rolling_moment = np.mean(total_moments[0, :])
            this_pitching_moment = np.mean(total_moments[1, :])
            this_yawing_moment = np.mean(total_moments[2, :])

            mean_forces[0] = this_induced_drag
            mean_forces[1] = this_side_force
            mean_forces[2] = this_lift
            mean_moments[0] = this_rolling_moment
            mean_moments[1] = this_pitching_moment
            mean_moments[2] = this_yawing_moment

            del this_induced_drag
            del this_side_force
            del this_lift
            del this_rolling_moment
            del this_pitching_moment
            del this_yawing_moment

            del total_forces
            del total_moments

            all_mean_forces[wake_state_id, num_flaps_id, num_chord_id] = mean_forces
            all_mean_moments[wake_state_id, num_flaps_id, num_chord_id] = mean_moments

            max_wake_percent_error = np.inf
            max_flap_percent_error = np.inf
            max_chord_percent_error = np.inf

            if wake_state_id > 0:
                last_wake_forces = all_mean_forces[
                    wake_state_id - 1, num_flaps_id, num_chord_id
                ]
                last_wake_moments = all_mean_moments[
                    wake_state_id - 1, num_flaps_id, num_chord_id
                ]
                wake_forces_percent_errors = 100 * np.abs(
                    (mean_forces - last_wake_forces) / np.max(np.abs(last_wake_forces))
                )
                wake_moments_percent_errors = 100 * np.abs(
                    (mean_moments - last_wake_moments)
                    / np.max(np.abs(last_wake_moments))
                )
                max_wake_percent_error = max(
                    np.max(wake_forces_percent_errors),
                    np.max(wake_moments_percent_errors),
                )
                del last_wake_forces
                del last_wake_moments
                del wake_forces_percent_errors
                del wake_moments_percent_errors
                print(
                    "\t\t\t\tMax Wake Percent Error: ",
                    round(max_wake_percent_error),
                    "%",
                    sep="",
                )

            if num_flaps_id > 0:
                last_flap_forces = all_mean_forces[
                    wake_state_id, num_flaps_id - 1, num_chord_id
                ]
                last_flap_moments = all_mean_moments[
                    wake_state_id, num_flaps_id - 1, num_chord_id
                ]
                flap_forces_percent_errors = 100 * np.abs(
                    (mean_forces - last_flap_forces) / np.max(np.abs(last_flap_forces))
                )
                flap_moments_percent_errors = 100 * np.abs(
                    (mean_moments - last_flap_moments)
                    / np.max(np.abs(last_flap_moments))
                )
                max_flap_percent_error = max(
                    np.max(flap_forces_percent_errors),
                    np.max(flap_moments_percent_errors),
                )
                del last_flap_forces
                del last_flap_moments
                del flap_forces_percent_errors
                del flap_moments_percent_errors
                print(
                    "\t\t\t\tMax Flap Percent Error: ",
                    round(max_flap_percent_error),
                    "%",
                    sep="",
                )

            if num_chord_id > 0:
                last_chord_forces = all_mean_forces[
                    wake_state_id, num_flaps_id, num_chord_id - 1
                ]
                last_chord_moments = all_mean_moments[
                    wake_state_id, num_flaps_id, num_chord_id - 1
                ]
                chord_forces_percent_errors = 100 * np.abs(
                    (mean_forces - last_chord_forces)
                    / np.max(np.abs(last_chord_forces))
                )
                chord_moments_percent_errors = 100 * np.abs(
                    (mean_moments - last_chord_moments)
                    / np.max(np.abs(last_chord_moments))
                )
                max_chord_percent_error = max(
                    np.max(chord_forces_percent_errors),
                    np.max(chord_moments_percent_errors),
                )
                del last_chord_forces
                del last_chord_moments
                del chord_forces_percent_errors
                del chord_moments_percent_errors
                print(
                    "\t\t\t\tMax Chord Percent Error: ",
                    round(max_chord_percent_error),
                    "%",
                    sep="",
                )

            max_percent_error = max(
                max_wake_percent_error,
                max_flap_percent_error,
                max_chord_percent_error,
            )
            if max_percent_error == np.inf:
                print(
                    "\t\t\t\tMax Percent Error: ",
                    max_percent_error,
                    sep="",
                )
            else:
                print(
                    "\t\t\t\tMax Percent Error: ",
                    round(max_percent_error),
                    "%",
                    sep="",
                )

            del num_results_steps
            del num_steps
            del first_results_step
            del iter_start
            del iter_stop
            del iter_time
            del mean_forces
            del mean_moments
            del max_wake_percent_error
            del max_flap_percent_error
            del max_chord_percent_error
            del max_percent_error
            del this_solver
            gc.collect()
            iteration += 1

stop_time = time.time()
elapsed_time = round(stop_time - start_time)
print(elapsed_time, " s", sep="")
