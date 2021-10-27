"""This is script is an example of running Ptera Software's steady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """
import src

alpha = 15
x_spacing = 0.5
y_spacing = 0.5

lead_airplane = src.geometry.Airplane(
    name="Lead Airplane",
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    twist=alpha,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
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
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    twist=alpha,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
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
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    twist=alpha,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    twist=alpha,
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
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
                    base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[
                        1
                    ],
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
                    base_wing_cross_section=right_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=right_airplane.wings[0].wing_cross_sections[
                        1
                    ],
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
                    base_wing_cross_section=left_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                ),
                src.movement.WingCrossSectionMovement(
                    base_wing_cross_section=left_airplane.wings[0].wing_cross_sections[
                        1
                    ],
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
    num_steps=100,
)

del lead_airplane_movement
del right_airplane_movement
del left_airplane_movement
del this_operating_point_movement

this_problem = src.problems.UnsteadyProblem(
    movement=this_movement,
    only_final_results=True,
)

this_solver = (
    src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=this_problem,
    )
)

del this_problem

this_solver.run(
    prescribed_wake=False,
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
#     keep_file=True,
# )
