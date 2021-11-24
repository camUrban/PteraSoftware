"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import time

import numpy as np

import src

start_time = time.time()

convergence = 1.0

num_airplanes = 3

aspect_ratio = 5.0
speed = 1.0
alpha = 0.0
x_spacing = 0.5
y_spacing = 0.5
root_to_mid_span = 0.2275
root_chord = 0.1094
mid_to_tip_span = 0.350 - 0.2275
tip_chord = 0.0219
flapping_amplitude = 15.0

period = x_spacing / speed
root_to_mid_chord = root_chord
mid_to_tip_chord = (root_chord + tip_chord) / 2

min_num_flaps = 2
max_num_flaps = 4
min_num_chord = 6
max_num_chord = 9

wake_state_list = [True, False]
num_flaps_list = [i for i in range(min_num_flaps, max_num_flaps + 1)]
num_chord_list = [i for i in range(min_num_chord, max_num_chord + 1)]

rms_lifts = np.zeros(
    (len(wake_state_list), len(num_flaps_list), len(num_chord_list), num_airplanes)
)
rms_drags = np.zeros(
    (len(wake_state_list), len(num_flaps_list), len(num_chord_list), num_airplanes)
)
iter_times = np.zeros((len(wake_state_list), len(num_flaps_list), len(num_chord_list)))

iteration = 0
num_iterations = len(wake_state_list) * len(num_flaps_list) * len(num_chord_list)

max_rmspe = None
max_non_wake_rmspe = None
wake_state = None
num_flaps = None
num_chord = None
iter_time = None
wake_state_id = None
num_flaps_id = None
num_chord_id = None

this_operating_point = src.operating_point.OperatingPoint(
    velocity=speed,
    alpha=0.0,
)
this_operating_point_movement = src.movement.OperatingPointMovement(
    base_operating_point=this_operating_point,
)
del this_operating_point

for wake_state_id, wake_state in enumerate(wake_state_list):
    print("Prescribed Wake:", wake_state)
    for num_flaps_id, num_flaps in enumerate(num_flaps_list):
        print("\tNumber of flaps:", num_flaps)
        for num_chord_id, num_chord in enumerate(num_chord_list):
            print("\t\tNumber of chordwise panels:", num_chord)

            iteration += 1
            print("\t\t\t\t", iteration, "/", num_iterations, sep="")

            root_to_mid_panel_chord = root_to_mid_chord / num_chord
            mid_to_tip_panel_chord = mid_to_tip_chord / num_chord

            root_to_mid_num_span = round(
                root_to_mid_span / (aspect_ratio * root_to_mid_panel_chord)
            )
            mid_to_tip_num_span = round(
                mid_to_tip_span / (aspect_ratio * mid_to_tip_panel_chord)
            )

            these_airplane_movements = []
            row = None
            position = None
            offset_sign = None
            for airplane_id in range(num_airplanes):
                if airplane_id == 0:
                    row = 1
                    position = ""
                    offset_sign = 1
                elif airplane_id % 2 != 0:
                    row += 1
                    position = "Right "
                    offset_sign = 1
                else:
                    position = "Left "
                    offset_sign = -1

                this_name = "Airplane " + position + str(row)

                offset = row - 1

                this_airplane = src.geometry.Airplane(
                    name=this_name,
                    x_ref=offset * x_spacing,
                    y_ref=offset_sign * offset * y_spacing,
                    wings=[
                        src.geometry.Wing(
                            name="Main Wing",
                            symmetric=True,
                            chordwise_spacing="uniform",
                            x_le=offset * x_spacing,
                            y_le=offset_sign * offset * y_spacing,
                            num_chordwise_panels=num_chord,
                            wing_cross_sections=[
                                src.geometry.WingCrossSection(
                                    twist=alpha,
                                    chord=root_chord,
                                    airfoil=src.geometry.Airfoil(name="naca0012"),
                                    num_spanwise_panels=root_to_mid_num_span,
                                    spanwise_spacing="uniform",
                                ),
                                src.geometry.WingCrossSection(
                                    twist=alpha,
                                    y_le=root_to_mid_span,
                                    chord=root_chord,
                                    airfoil=src.geometry.Airfoil(name="naca0012"),
                                    num_spanwise_panels=mid_to_tip_num_span,
                                    spanwise_spacing="uniform",
                                ),
                                src.geometry.WingCrossSection(
                                    twist=alpha,
                                    y_le=root_to_mid_span + mid_to_tip_span,
                                    chord=tip_chord,
                                    airfoil=src.geometry.Airfoil(name="naca0012"),
                                ),
                            ],
                        ),
                    ],
                )

                this_airplane_movement = src.movement.AirplaneMovement(
                    base_airplane=this_airplane,
                    wing_movements=[
                        src.movement.WingMovement(
                            base_wing=this_airplane.wings[0],
                            wing_cross_sections_movements=[
                                src.movement.WingCrossSectionMovement(
                                    base_wing_cross_section=this_airplane.wings[
                                        0
                                    ].wing_cross_sections[0],
                                ),
                                src.movement.WingCrossSectionMovement(
                                    base_wing_cross_section=this_airplane.wings[
                                        0
                                    ].wing_cross_sections[1],
                                    sweeping_amplitude=flapping_amplitude,
                                    sweeping_period=period,
                                    sweeping_spacing="sine",
                                ),
                                src.movement.WingCrossSectionMovement(
                                    base_wing_cross_section=this_airplane.wings[
                                        0
                                    ].wing_cross_sections[2],
                                    sweeping_amplitude=flapping_amplitude,
                                    sweeping_period=period,
                                    sweeping_spacing="sine",
                                ),
                            ],
                        )
                    ],
                )

                these_airplane_movements.append(this_airplane_movement)

                del this_airplane
                del this_airplane_movement

            this_movement = src.movement.Movement(
                airplane_movements=these_airplane_movements,
                operating_point_movement=this_operating_point_movement,
                num_steps=None,
                num_cycles=num_flaps,
                delta_time=None,
            )

            del these_airplane_movements

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
            iter_times[wake_state_id, num_flaps_id, num_chord_id] = iter_time

            first_results_step = this_solver.first_results_step
            num_steps = this_solver.num_steps
            num_results_steps = num_steps - first_results_step

            total_forces = np.zeros((num_airplanes, 3, num_results_steps))
            total_moments = np.zeros((num_airplanes, 3, num_results_steps))

            results_step = 0
            for step in range(first_results_step, num_steps):
                airplanes = this_solver.steady_problems[step].airplanes
                for airplane_id, airplane in enumerate(airplanes):
                    total_forces[
                        airplane_id, :, results_step
                    ] = airplane.total_near_field_force_wind_axes
                    total_moments[
                        airplane_id, :, results_step
                    ] = airplane.total_near_field_moment_wind_axes
                results_step += 1

            rms_drag = np.mean(total_forces[:, 0, :] ** 2) ** 0.5
            rms_lift = np.mean(total_forces[:, 2, :] ** 2) ** 0.5

            rms_drags[wake_state_id, num_flaps_id, num_chord_id] = rms_drag
            rms_lifts[wake_state_id, num_flaps_id, num_chord_id] = rms_lift

            max_wake_rmspe = np.inf
            max_flap_rmspe = np.inf
            max_chord_rmspe = np.inf

            if wake_state_id > 0:
                last_wake_rms_lifts = rms_lifts[
                    wake_state_id - 1, num_flaps_id, num_chord_id
                ]
                last_wake_rms_drags = rms_drags[
                    wake_state_id - 1, num_flaps_id, num_chord_id
                ]
                wake_lift_rmspes = 100 * np.abs(
                    (rms_lift - last_wake_rms_lifts) / rms_lift
                )
                wake_drag_rmspes = 100 * np.abs(
                    (rms_drag - last_wake_rms_drags) / rms_drag
                )
                max_wake_lift_rmspe = np.max(wake_lift_rmspes)
                max_wake_drag_rmspe = np.max(wake_drag_rmspes)
                max_wake_rmspe = max(max_wake_lift_rmspe, max_wake_drag_rmspe)

                print(
                    "\t\t\t\tMax Wake RMSPE: ",
                    round(max_wake_rmspe, 2),
                    "%",
                    sep="",
                )
            else:
                print("\t\t\t\tMax Wake RMSPE:", max_wake_rmspe)

            if num_flaps_id > 0:
                last_flap_rms_lifts = rms_lifts[
                    wake_state_id, num_flaps_id - 1, num_chord_id
                ]
                last_flap_rms_drags = rms_drags[
                    wake_state_id, num_flaps_id - 1, num_chord_id
                ]
                flap_lift_rmspes = 100 * np.abs(
                    (rms_lift - last_flap_rms_lifts) / rms_lift
                )
                flap_drag_rmspes = 100 * np.abs(
                    (rms_drag - last_flap_rms_drags) / rms_drag
                )
                max_flap_lift_rmspe = np.max(flap_lift_rmspes)
                max_flap_drag_rmspe = np.max(flap_drag_rmspes)
                max_flap_rmspe = max(max_flap_lift_rmspe, max_flap_drag_rmspe)

                print(
                    "\t\t\t\tMax Flap RMSPE: ",
                    round(max_flap_rmspe, 2),
                    "%",
                    sep="",
                )
            else:
                print("\t\t\t\tMax Flap RMSPE:", max_flap_rmspe)

            if num_chord_id > 0:
                last_chord_rms_lifts = rms_lifts[
                    wake_state_id, num_flaps_id, num_chord_id - 1
                ]
                last_chord_rms_drags = rms_drags[
                    wake_state_id, num_flaps_id, num_chord_id - 1
                ]
                chord_lift_rmspes = 100 * np.abs(
                    (rms_lift - last_chord_rms_lifts) / rms_lift
                )
                chord_drag_rmspes = 100 * np.abs(
                    (rms_drag - last_chord_rms_drags) / rms_drag
                )
                max_chord_lift_rmspe = np.max(chord_lift_rmspes)
                max_chord_drag_rmspe = np.max(chord_drag_rmspes)
                max_chord_rmspe = max(max_chord_lift_rmspe, max_chord_drag_rmspe)

                print(
                    "\t\t\t\tMax Chord RMSPE: ",
                    round(max_chord_rmspe, 2),
                    "%",
                    sep="",
                )
            else:
                print("\t\t\t\tMax Chord RMSPE:", max_chord_rmspe)

            max_rmspe = max(
                max_wake_rmspe,
                max_flap_rmspe,
                max_chord_rmspe,
            )
            max_non_wake_rmspe = max(
                max_flap_rmspe,
                max_chord_rmspe,
            )

            if max_rmspe == np.inf:
                print(
                    "\t\t\t\t\tMax RMSPE: ",
                    max_rmspe,
                    sep="",
                )
            else:
                print(
                    "\t\t\t\t\tMax RMSPE: ",
                    round(max_rmspe, 2),
                    "%",
                    sep="",
                )

            if max_non_wake_rmspe == np.inf:
                print(
                    "\t\t\t\t\tMax Non-Wake RMSPE: ",
                    max_non_wake_rmspe,
                    sep="",
                )
            else:
                print(
                    "\t\t\t\t\tMax Non-Wake RMSPE: ",
                    round(max_non_wake_rmspe, 2),
                    "%",
                    sep="",
                )

            converged = False

            if max_rmspe < convergence:
                converged = True
            elif max_non_wake_rmspe < convergence and wake_state is False:
                converged = True

            if converged:
                break
        else:
            continue
        break
    else:
        continue
    break

if max_rmspe < convergence or (
    max_non_wake_rmspe < convergence and wake_state is False
):
    if max_rmspe < convergence:
        converged_wake_state = True
        this_iter_time = iter_times[
            wake_state_id - 1, num_flaps_id - 1, num_chord_id - 1
        ]
    else:
        converged_wake_state = False
        this_iter_time = iter_times[wake_state_id, num_flaps_id - 1, num_chord_id - 1]

    converged_num_flaps = num_flaps - 1
    converged_num_chord = num_chord - 1

    print("\nThe simulation found a converged solution:")
    print("\tPrescribed wake:\t", converged_wake_state, sep="")
    print("\tNumber of flaps:\t", converged_num_flaps, sep="")
    print("\tNumber of chordwise panels:\t", converged_num_chord, sep="")
    print("\tSimulation time:\t", this_iter_time, " s", sep="")
else:
    print("\nThe simulation could not find a converged solution.")

stop_time = time.time()
elapsed_time = round(stop_time - start_time)
print("\nTotal Time: ", elapsed_time, " s", sep="")
