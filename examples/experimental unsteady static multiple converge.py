"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import time

import numpy as np
from tqdm import tqdm

import src

start_time = time.time()

alpha = 5
x_spacing = 0.5
y_spacing = 0.5

min_num_flaps = 2
max_num_flaps = 4
min_num_chord = 4
max_num_chord = 6
min_num_span = 8
max_num_span = 10

wake_state_list = [True, False]
num_flaps_list = [i for i in range(min_num_flaps, max_num_flaps + 1)]
num_chord_list = [i for i in range(min_num_chord, max_num_chord + 1)]
num_span_list = [i for i in range(min_num_span, max_num_span + 1)]

all_mean_forces = np.zeros(
    (
        len(wake_state_list),
        len(num_flaps_list),
        len(num_chord_list),
        len(num_span_list),
        3,
        3,
    )
)
all_mean_moments = np.zeros(
    (
        len(wake_state_list),
        len(num_flaps_list),
        len(num_chord_list),
        len(num_span_list),
        3,
        3,
    )
)

iteration = 0

times = []
for wake_state_val in range(1, len(wake_state_list) + 1):
    for num_flaps_val in range(1, len(num_flaps_list) + 1):
        for num_chord_val in range(1, len(num_chord_list) + 1):
            for num_span_val in range(1, len(num_span_list) + 1):
                this_time = (
                    wake_state_val ** 1.25
                    * num_flaps_val ** 1.5
                    * num_chord_val ** 2
                    * num_span_val ** 1
                )
                times.append(this_time)
total_time = sum(times)

with tqdm(
    total=total_time,
    unit="",
    unit_scale=True,
    ncols=100,
    desc="Simulating",
    bar_format="{desc}:{percentage:3.0f}% |{bar}| Elapsed: {elapsed}, Remaining: {remaining}",
) as bar:
    for wake_state_id, wake_state in enumerate(wake_state_list):
        # print("Prescribed Wake: ", wake_state)
        for num_flaps_id, num_flaps in enumerate(num_flaps_list):
            # print("\tNumber of flaps: ", num_flaps)
            for num_chord_id, num_chord in enumerate(num_chord_list):
                # print("\t\tNumber of chordwise panels: ", num_chord)
                for num_span_id, num_span in enumerate(num_span_list):
                    # print("\t\t\tNumber of spanwise panels: ", num_span)
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
                                        num_spanwise_panels=num_span,
                                    ),
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        y_le=0.2275,
                                        chord=0.1094,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        y_le=0.350,
                                        chord=0.0219,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                ],
                            ),
                        ],
                    )
                    right_airplane = src.geometry.Airplane(
                        name="Right Airplane",
                        x_ref=x_spacing,
                        y_ref=y_spacing,
                        wings=[
                            src.geometry.Wing(
                                name="Main Wing",
                                symmetric=True,
                                chordwise_spacing="uniform",
                                x_le=x_spacing,
                                y_le=y_spacing,
                                num_chordwise_panels=num_chord,
                                wing_cross_sections=[
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        chord=0.1094,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        y_le=0.2275,
                                        chord=0.1094,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        y_le=0.350,
                                        chord=0.0219,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                ],
                            ),
                        ],
                    )
                    left_airplane = src.geometry.Airplane(
                        name="Left Airplane",
                        x_ref=x_spacing,
                        y_ref=-y_spacing,
                        wings=[
                            src.geometry.Wing(
                                name="Main Wing",
                                symmetric=True,
                                chordwise_spacing="uniform",
                                x_le=x_spacing,
                                y_le=-y_spacing,
                                num_chordwise_panels=num_chord,
                                wing_cross_sections=[
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        chord=0.1094,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        y_le=0.2275,
                                        chord=0.1094,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                    src.geometry.WingCrossSection(
                                        twist=alpha,
                                        y_le=0.350,
                                        chord=0.0219,
                                        airfoil=src.geometry.Airfoil(name="naca0012"),
                                        num_spanwise_panels=num_span,
                                    ),
                                ],
                            ),
                        ],
                    )
                    this_operating_point = src.operating_point.OperatingPoint(
                        density=1.225,
                        velocity=1.0,
                        alpha=0.0,
                        nu=15.06e-6,
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
                                ],
                            )
                        ],
                    )
                    right_airplane_movement = src.movement.AirplaneMovement(
                        base_airplane=right_airplane,
                        wing_movements=[
                            src.movement.WingMovement(
                                base_wing=right_airplane.wings[0],
                                wing_cross_sections_movements=[
                                    src.movement.WingCrossSectionMovement(
                                        base_wing_cross_section=right_airplane.wings[
                                            0
                                        ].wing_cross_sections[0],
                                    ),
                                    src.movement.WingCrossSectionMovement(
                                        base_wing_cross_section=right_airplane.wings[
                                            0
                                        ].wing_cross_sections[1],
                                        sweeping_amplitude=15.0,
                                        sweeping_period=1,
                                        sweeping_spacing="sine",
                                    ),
                                ],
                            )
                        ],
                    )
                    left_airplane_movement = src.movement.AirplaneMovement(
                        base_airplane=left_airplane,
                        wing_movements=[
                            src.movement.WingMovement(
                                base_wing=left_airplane.wings[0],
                                wing_cross_sections_movements=[
                                    src.movement.WingCrossSectionMovement(
                                        base_wing_cross_section=left_airplane.wings[
                                            0
                                        ].wing_cross_sections[0],
                                    ),
                                    src.movement.WingCrossSectionMovement(
                                        base_wing_cross_section=left_airplane.wings[
                                            0
                                        ].wing_cross_sections[1],
                                        sweeping_amplitude=15.0,
                                        sweeping_period=1,
                                        sweeping_spacing="sine",
                                    ),
                                ],
                            )
                        ],
                    )
                    this_operating_point_movement = src.movement.OperatingPointMovement(
                        base_operating_point=this_operating_point,
                    )
                    this_movement = src.movement.Movement(
                        airplane_movements=[
                            lead_airplane_movement,
                            right_airplane_movement,
                            left_airplane_movement,
                        ],
                        operating_point_movement=this_operating_point_movement,
                        num_steps=None,
                        num_cycles=num_flaps,
                        delta_time=None,
                    )

                    del lead_airplane_movement
                    del right_airplane_movement
                    del left_airplane_movement
                    del this_operating_point_movement

                    this_problem = src.problems.UnsteadyProblem(
                        movement=this_movement,
                        only_final_results=True,
                    )
                    this_solver = src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        unsteady_problem=this_problem,
                    )

                    del this_problem

                    this_solver.run(
                        logging_level="Critical",
                        prescribed_wake=wake_state,
                        calculate_streamlines=False,
                    )

                    first_results_step = this_solver.first_results_step
                    num_steps = this_solver.num_steps
                    num_airplanes = this_solver.num_airplanes
                    num_results_steps = num_steps - first_results_step

                    total_forces = np.zeros((num_airplanes, 3, num_results_steps))
                    total_moments = np.zeros((num_airplanes, 3, num_results_steps))

                    mean_forces = np.zeros((num_airplanes, 3))
                    mean_moments = np.zeros((num_airplanes, 3))

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

                    for airplane_id in range(num_airplanes):
                        this_induced_drag = np.mean(total_forces[airplane_id, 0, :])
                        this_side_force = np.mean(total_forces[airplane_id, 1, :])
                        this_lift = np.mean(total_forces[airplane_id, 2, :])
                        this_rolling_moment = np.mean(total_moments[airplane_id, 0, :])
                        this_pitching_moment = np.mean(total_moments[airplane_id, 1, :])
                        this_yawing_moment = np.mean(total_moments[airplane_id, 2, :])

                        mean_forces[airplane_id, 0] = this_induced_drag
                        mean_forces[airplane_id, 1] = this_side_force
                        mean_forces[airplane_id, 2] = this_lift
                        mean_moments[airplane_id, 0] = this_rolling_moment
                        mean_moments[airplane_id, 1] = this_pitching_moment
                        mean_moments[airplane_id, 2] = this_yawing_moment

                    all_mean_forces[
                        wake_state_id, num_flaps_id, num_chord_id, num_span_id
                    ] = mean_forces
                    all_mean_moments[
                        wake_state_id, num_flaps_id, num_chord_id, num_span_id
                    ] = mean_moments

                    bar.update(times[iteration])
                    iteration += 1

for wake_state_id, wake_state in enumerate(wake_state_list):
    print("Prescribed Wake: ", wake_state)
    for num_flaps_id, num_flaps in enumerate(num_flaps_list):
        print("\tNumber of flaps: ", num_flaps)
        for num_chord_id, num_chord in enumerate(num_chord_list):
            print("\t\tNumber of chordwise panels: ", num_chord)
            for num_span_id, num_span in enumerate(num_span_list):
                print("\t\t\tNumber of spanwise panels: ", num_span)
                these_forces = all_mean_forces[
                    wake_state_id, num_flaps_id, num_chord_id, num_span_id
                ]
                these_moments = all_mean_moments[
                    wake_state_id, num_flaps_id, num_chord_id, num_span_id
                ]

                max_wake_percent_error = np.inf
                max_flap_percent_error = np.inf
                max_chord_percent_error = np.inf
                max_span_percent_error = np.inf

                if wake_state_id > 0:
                    last_wake_forces = all_mean_forces[
                        wake_state_id - 1, num_flaps_id, num_chord_id, num_span_id
                    ]
                    last_wake_moments = all_mean_moments[
                        wake_state_id - 1, num_flaps_id, num_chord_id, num_span_id
                    ]
                    wake_forces_percent_errors = 100 * np.abs(
                        (these_forces - last_wake_forces) / np.mean(last_wake_forces)
                    )
                    wake_moments_percent_errors = 100 * np.abs(
                        (these_moments - last_wake_moments) / np.mean(last_wake_moments)
                    )
                    max_wake_percent_error = max(
                        np.max(wake_forces_percent_errors),
                        np.max(wake_moments_percent_errors),
                    )
                    print(
                        "\t\t\t\tMax Wake Percent Error: ",
                        round(max_wake_percent_error),
                        "%",
                        sep="",
                    )

                if num_flaps_id > 0:
                    last_flap_forces = all_mean_forces[
                        wake_state_id, num_flaps_id - 1, num_chord_id, num_span_id
                    ]
                    last_flap_moments = all_mean_moments[
                        wake_state_id, num_flaps_id - 1, num_chord_id, num_span_id
                    ]
                    flap_forces_percent_errors = 100 * np.abs(
                        (these_forces - last_flap_forces) / np.mean(last_flap_forces)
                    )
                    flap_moments_percent_errors = 100 * np.abs(
                        (these_moments - last_flap_moments) / np.mean(last_flap_moments)
                    )
                    max_flap_percent_error = max(
                        np.max(flap_forces_percent_errors),
                        np.max(flap_moments_percent_errors),
                    )
                    print(
                        "\t\t\t\tMax Flap Percent Error: ",
                        round(max_flap_percent_error),
                        "%",
                        sep="",
                    )

                if num_chord_id > 0:
                    last_chord_forces = all_mean_forces[
                        wake_state_id, num_flaps_id, num_chord_id - 1, num_span_id
                    ]
                    last_chord_moments = all_mean_moments[
                        wake_state_id, num_flaps_id, num_chord_id - 1, num_span_id
                    ]
                    chord_forces_percent_errors = 100 * np.abs(
                        (these_forces - last_chord_forces) / np.mean(last_chord_forces)
                    )
                    chord_moments_percent_errors = 100 * np.abs(
                        (these_moments - last_chord_moments)
                        / np.mean(last_chord_moments)
                    )
                    max_chord_percent_error = max(
                        np.max(chord_forces_percent_errors),
                        np.max(chord_moments_percent_errors),
                    )
                    print(
                        "\t\t\t\tMax Chord Percent Error: ",
                        round(max_chord_percent_error),
                        "%",
                        sep="",
                    )

                if num_span_id > 0:
                    last_span_forces = all_mean_forces[
                        wake_state_id, num_flaps_id, num_chord_id, num_span_id - 1
                    ]
                    last_span_moments = all_mean_moments[
                        wake_state_id, num_flaps_id, num_chord_id, num_span_id - 1
                    ]
                    span_forces_percent_errors = 100 * np.abs(
                        (these_forces - last_span_forces) / np.mean(last_span_forces)
                    )
                    span_moments_percent_errors = 100 * np.abs(
                        (these_moments - last_span_moments) / np.mean(last_span_moments)
                    )
                    max_span_percent_error = max(
                        np.max(span_forces_percent_errors),
                        np.max(span_moments_percent_errors),
                    )
                    print(
                        "\t\t\t\tMax Span Percent Error: ",
                        round(max_span_percent_error),
                        "%",
                        sep="",
                    )

                max_percent_error = max(
                    max_wake_percent_error,
                    max_flap_percent_error,
                    max_chord_percent_error,
                    max_span_percent_error,
                )
                if max_percent_error is not np.inf:
                    print(
                        "\t\t\t\t\tMax Percent Error: ",
                        round(max_percent_error),
                        "%",
                        sep="",
                    )

stop_time = time.time()
elapsed_time = round(stop_time - start_time)
print("Elapsed Time: ", elapsed_time, " s", sep="")
