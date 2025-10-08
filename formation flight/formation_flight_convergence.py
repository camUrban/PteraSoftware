# REFACTOR: I haven't yet started refactoring this module.
# DOCUMENT: Document this script.
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import pterasoftware as ps

start_time = time.time()

# Known Converged Values (0.5% convergence, 0 degrees angle of attack):
#   1 Airplane (flaps, chordwise panels, prescribed wake):
#       2, 7, False
#   3 Airplanes (flaps, chordwise panels, prescribed wake):
#       3, 7, False
#   5 Airplanes (flaps, chordwise panels, prescribed wake):
#       3, 7, False
convergence = 0.5
num_airplanes = 5
min_num_flaps = 1
max_num_flaps = 4
min_num_chord = 5
max_num_chord = 13
wake_state_list = [True, False]

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

converged = None
prescribed_wake = None
num_flaps = None
num_chord = None
iter_time = None
wake_state_id = None
num_flaps_id = None
num_chord_id = None
single_wake = None
single_flap = None
single_chord = None
wake_saturated = None
this_solver = None

this_operating_point = ps.operating_point.OperatingPoint(vCg__E=speed, alpha=0.0)
this_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=this_operating_point
    )
)
del this_operating_point

for wake_state_id, prescribed_wake in enumerate(wake_state_list):
    print("Prescribed Wake:", prescribed_wake)
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

                this_airplane = ps.geometry.airplane.Airplane(
                    wings=[
                        ps.geometry.wing.Wing(
                            wing_cross_sections=[
                                ps.geometry.wing_cross_section.WingCrossSection(
                                    twist=alpha,
                                    chord=root_chord,
                                    airfoil=ps.geometry.airfoil.Airfoil(
                                        name="naca0012"
                                    ),
                                    num_spanwise_panels=root_to_mid_num_span,
                                    spanwise_spacing="cosine",
                                ),
                                ps.geometry.wing_cross_section.WingCrossSection(
                                    twist=alpha,
                                    y_le=root_to_mid_span,
                                    chord=root_chord,
                                    airfoil=ps.geometry.airfoil.Airfoil(
                                        name="naca0012"
                                    ),
                                    num_spanwise_panels=mid_to_tip_num_span,
                                    spanwise_spacing="cosine",
                                ),
                                ps.geometry.wing_cross_section.WingCrossSection(
                                    twist=alpha,
                                    y_le=root_to_mid_span + mid_to_tip_span,
                                    chord=tip_chord,
                                    airfoil=ps.geometry.airfoil.Airfoil(
                                        name="naca0012"
                                    ),
                                ),
                            ],
                            name="Main Wing",
                            symmetric=True,
                            num_chordwise_panels=num_chord,
                            chordwise_spacing="uniform",
                        ),
                    ],
                    name=this_name,
                )

                this_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
                    base_airplane=this_airplane,
                    wing_movements=[
                        ps.movements.wing_movement.WingMovement(
                            base_wing=this_airplane.wings[0],
                            wing_cross_section_movements=[
                                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                                    base_wing_cross_section=this_airplane.wings[
                                        0
                                    ].wing_cross_sections[0]
                                ),
                                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                                    base_wing_cross_section=this_airplane.wings[
                                        0
                                    ].wing_cross_sections[1]
                                ),
                                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                                    base_wing_cross_section=this_airplane.wings[
                                        0
                                    ].wing_cross_sections[2]
                                ),
                            ],
                        )
                    ],
                )

                these_airplane_movements.append(this_airplane_movement)

                del this_airplane
                del this_airplane_movement

            this_movement = ps.movements.movement.Movement(
                airplane_movements=these_airplane_movements,
                operating_point_movement=this_operating_point_movement,
                num_steps=None,
                num_cycles=num_flaps,
                delta_time=None,
            )

            del these_airplane_movements

            this_problem = ps.problems.UnsteadyProblem(
                movement=this_movement,
                only_final_results=True,
            )

            del this_movement

            this_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                unsteady_problem=this_problem,
            )

            del this_problem

            iter_start = time.time()

            this_solver.run(
                logging_level="Critical",
                prescribed_wake=prescribed_wake,
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
                    total_forces[airplane_id, :, results_step] = (
                        airplane.total_near_field_force_wind_axes
                    )
                    total_moments[airplane_id, :, results_step] = (
                        airplane.total_near_field_moment_wind_axes
                    )
                results_step += 1

            these_s_drags = total_forces[:, 0, :] ** 2
            these_s_lifts = total_forces[:, 2, :] ** 2

            these_ms_drags = np.mean(these_s_drags, axis=-1)
            these_ms_lifts = np.mean(these_s_lifts, axis=-1)

            these_rms_drags = these_ms_drags**0.5
            these_rms_lifts = these_ms_lifts**0.5

            rms_drags[wake_state_id, num_flaps_id, num_chord_id, :] = these_rms_drags
            rms_lifts[wake_state_id, num_flaps_id, num_chord_id, :] = these_rms_lifts

            max_wake_rmspc = np.inf
            max_flap_rmspc = np.inf
            max_chord_rmspc = np.inf

            if wake_state_id > 0:
                last_wake_rms_lifts = rms_lifts[
                    wake_state_id - 1, num_flaps_id, num_chord_id, :
                ]
                last_wake_rms_drags = rms_drags[
                    wake_state_id - 1, num_flaps_id, num_chord_id, :
                ]
                wake_lift_rmspcs = 100 * np.abs(
                    (these_rms_lifts - last_wake_rms_lifts) / last_wake_rms_lifts
                )
                wake_drag_rmspcs = 100 * np.abs(
                    (these_rms_drags - last_wake_rms_drags) / last_wake_rms_drags
                )
                max_wake_lift_rmspc = np.max(wake_lift_rmspcs)
                max_wake_drag_rmspc = np.max(wake_drag_rmspcs)
                max_wake_rmspc = max(max_wake_lift_rmspc, max_wake_drag_rmspc)

                print(
                    "\t\t\t\tMax Wake RMSPC: ",
                    round(max_wake_rmspc, 2),
                    "%",
                    sep="",
                )
            else:
                print("\t\t\t\tMax Wake RMSPC:", max_wake_rmspc)

            if num_flaps_id > 0:
                last_flap_rms_lifts = rms_lifts[
                    wake_state_id, num_flaps_id - 1, num_chord_id, :
                ]
                last_flap_rms_drags = rms_drags[
                    wake_state_id, num_flaps_id - 1, num_chord_id, :
                ]
                flap_lift_rmspcs = 100 * np.abs(
                    (these_rms_lifts - last_flap_rms_lifts) / last_flap_rms_lifts
                )
                flap_drag_rmspcs = 100 * np.abs(
                    (these_rms_drags - last_flap_rms_drags) / last_flap_rms_drags
                )
                max_flap_lift_rmspc = np.max(flap_lift_rmspcs)
                max_flap_drag_rmspc = np.max(flap_drag_rmspcs)
                max_flap_rmspc = max(max_flap_lift_rmspc, max_flap_drag_rmspc)

                print(
                    "\t\t\t\tMax Flap RMSPC: ",
                    round(max_flap_rmspc, 2),
                    "%",
                    sep="",
                )
            else:
                print("\t\t\t\tMax Flap RMSPC:", max_flap_rmspc)

            if num_chord_id > 0:
                last_chord_rms_lifts = rms_lifts[
                    wake_state_id, num_flaps_id, num_chord_id - 1, :
                ]
                last_chord_rms_drags = rms_drags[
                    wake_state_id, num_flaps_id, num_chord_id - 1, :
                ]
                chord_lift_rmspcs = 100 * np.abs(
                    (these_rms_lifts - last_chord_rms_lifts) / last_chord_rms_lifts
                )
                chord_drag_rmspcs = 100 * np.abs(
                    (these_rms_drags - last_chord_rms_drags) / last_chord_rms_drags
                )
                max_chord_lift_rmspc = np.max(chord_lift_rmspcs)
                max_chord_drag_rmspc = np.max(chord_drag_rmspcs)
                max_chord_rmspc = max(max_chord_lift_rmspc, max_chord_drag_rmspc)

                print(
                    "\t\t\t\tMax Chord RMSPC: ",
                    round(max_chord_rmspc, 2),
                    "%",
                    sep="",
                )
            else:
                print("\t\t\t\tMax Chord RMSPC:", max_chord_rmspc)

            single_wake = len(wake_state_list) == 1
            single_flap = len(num_flaps_list) == 1
            single_chord = len(num_chord_list) == 1

            wake_converged = max_wake_rmspc < convergence
            flap_converged = max_flap_rmspc < convergence
            chord_converged = max_chord_rmspc < convergence

            wake_passed = wake_converged or single_wake
            flap_passed = flap_converged or single_flap
            chord_passed = chord_converged or single_chord

            wake_saturated = prescribed_wake is False

            converged = False
            if wake_passed and flap_passed and chord_passed:
                converged = True
            elif wake_saturated and flap_passed and chord_passed:
                converged = True

            if converged:
                break
        else:
            continue
        break
    else:
        continue
    break

plot_wake_state_id = len(wake_state_list) - 1
plot_num_flaps_id = len(num_flaps_list) - 1
plot_num_chord_id = len(num_chord_list) - 1

if converged:
    print("\nThe simulation found a converged solution:")

    if single_wake and wake_saturated:
        converged_wake_state_id = wake_state_list.index(False)
    elif single_wake:
        print("\tWarning: Wake state convergence not checked.")
        converged_wake_state_id = wake_state_id
    elif wake_saturated:
        converged_wake_state_id = wake_state_list.index(False)
    else:
        converged_wake_state_id = wake_state_id - 1

    if single_flap:
        print("\tWarning: Flap number convergence not checked.")
        converged_num_flaps_id = num_flaps_id
    else:
        converged_num_flaps_id = num_flaps_id - 1

    if single_chord:
        print("\tWarning: Chordwise panel number convergence not checked.")
        converged_num_chord_id = num_chord_id
    else:
        converged_num_chord_id = num_chord_id - 1

    converged_wake_state = wake_state_list[converged_wake_state_id]
    converged_num_flaps = num_flaps_list[converged_num_flaps_id]
    converged_num_chord = num_chord_list[converged_num_chord_id]
    this_iter_time = iter_times[
        converged_wake_state_id,
        converged_num_flaps_id,
        converged_num_chord_id,
    ]

    plot_wake_state_id = converged_wake_state_id
    plot_num_flaps_id = converged_num_flaps_id
    plot_num_chord_id = converged_num_chord_id

    print("\tPrescribed wake:\t", converged_wake_state, sep="")
    print("\tNumber of flaps:\t", converged_num_flaps, sep="")
    print("\tNumber of chordwise panels:\t", converged_num_chord, sep="")
    print("\tSimulation time:\t", this_iter_time, " s", sep="")
else:
    print("\nThe simulation could not find a converged solution.")

if not single_chord:
    lift_figure, lift_axes = plt.subplots()
    drag_figure, drag_axes = plt.subplots()

    row = None
    for airplane_id in range(num_airplanes):
        if airplane_id == 0:
            row = 1
        elif airplane_id % 2 == 0:
            row += 1
        else:
            continue

        these_rms_lifts = rms_lifts[
            plot_wake_state_id,
            plot_num_flaps_id,
            : (plot_num_chord_id + 2),
            airplane_id,
        ]
        these_rms_drags = rms_drags[
            plot_wake_state_id,
            plot_num_flaps_id,
            : (plot_num_chord_id + 2),
            airplane_id,
        ]

        lift_axes.plot(
            num_chord_list[: (plot_num_chord_id + 2)],
            these_rms_lifts,
            label="Row " + str(row),
            marker="o",
            linestyle="--",
        )
        drag_axes.plot(
            num_chord_list[: (plot_num_chord_id + 2)],
            these_rms_drags,
            label="Row " + str(row),
            marker="o",
            linestyle="--",
        )

    lift_axes.set_xlabel("Number of Chordwise Panels")
    drag_axes.set_xlabel("Number of Chordwise Panels")

    lift_axes.set_ylabel("RMS Lift (N)")
    drag_axes.set_ylabel("RMS Drag (N)")

    lift_axes.set_title("Number of Chordwise Panels\nLift Convergence")
    drag_axes.set_title("Number of Chordwise Panels\nDrag Convergence")

    lift_axes.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
    drag_axes.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))

    lift_axes.set_ylim(bottom=0, top=0.0625)
    drag_axes.set_ylim(bottom=0, top=0.0250)

    lift_axes.legend(loc="lower left")
    drag_axes.legend(loc="lower left")

    lift_figure.show()
    drag_figure.show()

if not single_flap:
    lift_figure, lift_axes = plt.subplots()
    drag_figure, drag_axes = plt.subplots()

    row = None
    for airplane_id in range(num_airplanes):
        if airplane_id == 0:
            row = 1
        elif airplane_id % 2 == 0:
            row += 1
        else:
            continue

        these_rms_lifts = rms_lifts[
            plot_wake_state_id,
            : (plot_num_flaps_id + 2),
            plot_num_chord_id,
            airplane_id,
        ]
        these_rms_drags = rms_drags[
            plot_wake_state_id,
            : (plot_num_flaps_id + 2),
            plot_num_chord_id,
            airplane_id,
        ]

        lift_axes.plot(
            num_flaps_list[: (plot_num_flaps_id + 2)],
            these_rms_lifts,
            label="Row " + str(row),
            marker="o",
            linestyle="--",
        )
        drag_axes.plot(
            num_flaps_list[: (plot_num_flaps_id + 2)],
            these_rms_drags,
            label="Row " + str(row),
            marker="o",
            linestyle="--",
        )

    lift_axes.set_xlabel("Number of Flap Cycles")
    drag_axes.set_xlabel("Number of Flap Cycles")

    lift_axes.set_ylabel("RMS Lift (N)")
    drag_axes.set_ylabel("RMS Drag (N)")

    lift_axes.set_title("Number of Flap Cycles\nLift Convergence")
    drag_axes.set_title("Number of Flap Cycles\nDrag Convergence")

    lift_axes.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
    drag_axes.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))

    lift_axes.set_ylim(bottom=0, top=0.0625)
    drag_axes.set_ylim(bottom=0, top=0.0250)

    lift_axes.legend(loc="lower left")
    drag_axes.legend(loc="lower left")

    lift_figure.show()
    drag_figure.show()

np.save("rms_lifts", rms_lifts)
np.save("rms_drags", rms_drags)
np.save("iter_times", iter_times)

stop_time = time.time()
elapsed_time = round(stop_time - start_time)
print("\nTotal Time: ", elapsed_time, " s", sep="")
