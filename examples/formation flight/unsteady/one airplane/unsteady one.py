"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import src

wake_state = True
num_cycles = 2
num_chord = 11

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
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                    num_spanwise_panels=root_to_mid_num_span,
                    spanwise_spacing="uniform",
                ),
                src.geometry.WingCrossSection(
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                    num_spanwise_panels=mid_to_tip_num_span,
                    spanwise_spacing="uniform",
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0875,
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
    alpha=5.0,
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

this_operating_point_movement = src.movement.OperatingPointMovement(
    base_operating_point=this_operating_point,
)

this_movement = src.movement.Movement(
    airplane_movements=[airplane_1_movement],
    operating_point_movement=this_operating_point_movement,
    num_cycles=num_cycles,
)

del airplane_1_movement
del this_operating_point_movement

this_problem = src.problems.UnsteadyProblem(
    movement=this_movement, only_final_results=False
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

# src.output.animate(
#     unsteady_solver=this_solver,
#     show_delta_pressures=True,
#     show_wake_vortices=True,
#     keep_file=False,
# )
