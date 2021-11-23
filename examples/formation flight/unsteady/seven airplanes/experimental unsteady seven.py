"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import src

wake_state = True
num_flaps = 2
num_chord = 4

speed = 1.0
alpha = 5
x_spacing = 0.5
y_spacing = 0.5
period = x_spacing / speed

root_to_mid_span = 0.2275
root_chord = 0.1094
mid_to_tip_span = 0.350 - 0.2275
tip_chord = 0.0219

root_to_mid_chord = root_chord
mid_to_tip_chord = (root_chord + tip_chord) / 2

root_to_mid_panel_chord = root_to_mid_chord / num_chord
mid_to_tip_panel_chord = mid_to_tip_chord / num_chord

root_to_mid_num_span = round(root_to_mid_span / (5 * root_to_mid_panel_chord))
mid_to_tip_num_span = round(mid_to_tip_span / (5 * mid_to_tip_panel_chord))

airplane_1 = src.geometry.Airplane(
    name="Airplane 1",
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

airplane_right_2 = src.geometry.Airplane(
    name="Airplane Right 2",
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

airplane_left_2 = src.geometry.Airplane(
    name="Airplane Left 2",
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

airplane_right_3 = src.geometry.Airplane(
    name="Airplane Right 3",
    x_ref=2 * x_spacing,
    y_ref=2 * y_spacing,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            x_le=2 * x_spacing,
            y_le=2 * y_spacing,
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

airplane_left_3 = src.geometry.Airplane(
    name="Airplane Left 3",
    x_ref=2 * x_spacing,
    y_ref=-2 * y_spacing,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            x_le=2 * x_spacing,
            y_le=-2 * y_spacing,
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

airplane_right_4 = src.geometry.Airplane(
    name="Airplane Right 4",
    x_ref=3 * x_spacing,
    y_ref=3 * y_spacing,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            x_le=3 * x_spacing,
            y_le=3 * y_spacing,
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

airplane_left_4 = src.geometry.Airplane(
    name="Airplane Left 4",
    x_ref=3 * x_spacing,
    y_ref=-3 * y_spacing,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            x_le=3 * x_spacing,
            y_le=-3 * y_spacing,
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

this_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
    velocity=speed,
    alpha=0.0,
    nu=15.06e-6,
)

airplane_1_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_1,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_1.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_1.wings[0].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_1.wings[0].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_1.wings[0].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

airplane_right_2_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_right_2,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_right_2.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_2.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_2.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_2.wings[
                        0
                    ].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

airplane_left_2_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_left_2,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_left_2.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_2.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_2.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_2.wings[
                        0
                    ].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

airplane_right_3_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_right_3,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_right_3.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_3.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_3.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_3.wings[
                        0
                    ].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

airplane_left_3_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_left_3,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_left_3.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_3.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_3.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_3.wings[
                        0
                    ].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

airplane_right_4_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_right_4,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_right_4.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_4.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_4.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_right_4.wings[
                        0
                    ].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
            ],
        )
    ],
)

airplane_left_4_movement = src.movement.AirplaneMovement(
    base_airplane=airplane_left_4,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=airplane_left_4.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_4.wings[
                        0
                    ].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_4.wings[
                        0
                    ].wing_cross_sections[1],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
                    sweeping_spacing="sine",
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=airplane_left_4.wings[
                        0
                    ].wing_cross_sections[2],
                    sweeping_amplitude=15.0,
                    sweeping_period=period,
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
        airplane_1_movement,
        airplane_right_2_movement,
        airplane_left_2_movement,
        airplane_right_3_movement,
        airplane_left_3_movement,
        airplane_right_4_movement,
        airplane_left_4_movement,
    ],
    operating_point_movement=this_operating_point_movement,
    num_cycles=num_flaps,
)

del airplane_1_movement
del airplane_right_2_movement
del airplane_left_2_movement
del airplane_right_3_movement
del airplane_left_3_movement
del airplane_right_4_movement
del airplane_left_4_movement
del this_operating_point_movement

this_problem = src.problems.UnsteadyProblem(
    movement=this_movement,
    only_final_results=False,
)

this_solver = (
    src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=this_problem,
    )
)

del this_problem

this_solver.run(
    prescribed_wake=wake_state,
    calculate_streamlines=False,
)

src.output.draw(
    solver=this_solver,
    show_delta_pressures=True,
    show_wake_vortices=True,
)

# src.output.plot_results_versus_time(unsteady_solver=this_solver, testing=False)

src.output.print_unsteady_results(unsteady_solver=this_solver)

src.output.animate(
    unsteady_solver=this_solver,
    show_delta_pressures=True,
    show_wake_vortices=True,
    keep_file=True,
)
