"""This is script is an example of running Ptera Software's unsteady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """

import src

example_airplane = src.geometry.Airplane(
    name="Example Airplane",
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            x_le=0.0,
            y_le=0.0,
            z_le=0.0,
            symmetric=True,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    num_spanwise_panels=8,
                    spanwise_spacing="cosine",
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.2275,
                    z_le=0.0,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0875,
                    y_le=0.350,
                    z_le=0.0,
                    chord=0.0219,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
    ],
)

main_wing_root_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
)

main_wing_mid_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
    sweeping_amplitude=30.0,
    sweeping_period=1.0,
)

main_wing_tip_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[2],
    sweeping_amplitude=30.0,
    sweeping_period=1.0,
)

main_wing_movement = src.movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_sections_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_mid_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
)

del main_wing_root_wing_cross_section_movement
del main_wing_mid_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement

airplane_movement = src.movement.AirplaneMovement(  # Define the base airplane object.
    base_airplane=example_airplane,  # Add the list of wing movement objects.
    wing_movements=[main_wing_movement],
)

del main_wing_movement

example_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
    beta=0.0,
    velocity=1.0,
    alpha=1.0,
    nu=15.06e-6,
)

operating_point_movement = src.movement.OperatingPointMovement(
    base_operating_point=example_operating_point,
)

movement = src.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    num_steps=80,
    delta_time=0.025,
)

del airplane_movement
del operating_point_movement

example_problem = src.problems.UnsteadyProblem(
    movement=movement,
)

example_solver = (
    src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)

del example_problem

example_solver.run(
    logging_level="Warning",
    prescribed_wake=True,
)

src.output.animate(
    unsteady_solver=example_solver,
    show_delta_pressures=True,
    show_wake_vortices=True,
    keep_file=False,
)

src.output.plot_results_versus_time(
    unsteady_solver=example_solver,
    testing=False,
)
