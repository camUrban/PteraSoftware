"""This is script is an example of how to run Ptera Software's unsteady ring vortex
lattice method solver on a custom airplane with variable geometry."""

import src

steady = False

example_airplane = src.geometry.Airplane(
    name="Example Airplane",
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
    wings=[
        src.geometry.Wing(
            name="Lead Wing",
            x_le=0.0,
            y_le=0.0,
            z_le=0.0,
            symmetric=False,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    x_le=0.75,
                    y_le=-6.0,
                    z_le=1.0,
                    chord=1.5,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    num_spanwise_panels=8,
                    spanwise_spacing="cosine",
                    chord=1.75,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.75,
                    y_le=6.0,
                    z_le=1.0,
                    chord=1.5,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
        src.geometry.Wing(
            name="Right Wing",
            x_le=15.0,
            y_le=10.0,
            z_le=0.0,
            symmetric=False,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    x_le=0.75,
                    y_le=-6.0,
                    z_le=1.0,
                    chord=1.5,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    num_spanwise_panels=8,
                    spanwise_spacing="cosine",
                    chord=1.75,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.75,
                    y_le=6.0,
                    z_le=1.0,
                    chord=1.5,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
        src.geometry.Wing(
            name="Left Wing",
            x_le=15.0,
            y_le=-10.0,
            z_le=0.0,
            symmetric=False,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    x_le=0.75,
                    y_le=-6.0,
                    z_le=1.0,
                    chord=1.5,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    num_spanwise_panels=8,
                    spanwise_spacing="cosine",
                    chord=1.75,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    x_le=0.75,
                    y_le=6.0,
                    z_le=1.0,
                    chord=1.5,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
            ],
        ),
    ],
)

example_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
    velocity=10.0,
    alpha=5.0,
    nu=15.06e-6,
)

if steady:
    steady_problem = src.problems.SteadyProblem(
        airplane=example_airplane, operating_point=example_operating_point
    )

    steady_solver = src.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=steady_problem,
    )

    del steady_problem

    steady_solver.run(
        logging_level="Warning",
    )

    src.output.draw(
        solver=steady_solver,
        show_delta_pressures=True,
    )
else:
    lead_wing_left_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
        sweeping_amplitude=30.0,
        sweeping_period=1.0,
        sweeping_spacing="sine",
    )
    lead_wing_root_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
    )
    lead_wing_right_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[2],
        sweeping_amplitude=30.0,
        sweeping_period=1.0,
        sweeping_spacing="sine",
    )

    right_wing_left_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[0],
        sweeping_amplitude=30.0,
        sweeping_period=1.0,
        sweeping_spacing="sine",
    )
    right_wing_root_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[1],
    )
    right_wing_right_wing_cross_section_movement = (
        src.movement.WingCrossSectionMovement(
            base_wing_cross_section=example_airplane.wings[1].wing_cross_sections[2],
            sweeping_amplitude=30.0,
            sweeping_period=1.0,
            sweeping_spacing="sine",
        )
    )

    left_wing_left_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[2].wing_cross_sections[0],
        sweeping_amplitude=30.0,
        sweeping_period=1.0,
        sweeping_spacing="sine",
    )
    left_wing_root_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[2].wing_cross_sections[1],
    )
    left_wing_right_wing_cross_section_movement = src.movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[2].wing_cross_sections[2],
        sweeping_amplitude=30.0,
        sweeping_period=1.0,
        sweeping_spacing="sine",
    )

    lead_wing_movement = src.movement.WingMovement(
        base_wing=example_airplane.wings[0],
        wing_cross_sections_movements=[
            lead_wing_left_wing_cross_section_movement,
            lead_wing_root_wing_cross_section_movement,
            lead_wing_right_wing_cross_section_movement,
        ],
    )

    right_wing_movement = src.movement.WingMovement(
        base_wing=example_airplane.wings[1],
        wing_cross_sections_movements=[
            right_wing_left_wing_cross_section_movement,
            right_wing_root_wing_cross_section_movement,
            right_wing_right_wing_cross_section_movement,
        ],
    )

    left_wing_movement = src.movement.WingMovement(
        base_wing=example_airplane.wings[2],
        wing_cross_sections_movements=[
            left_wing_left_wing_cross_section_movement,
            left_wing_root_wing_cross_section_movement,
            left_wing_right_wing_cross_section_movement,
        ],
    )

    del lead_wing_left_wing_cross_section_movement
    del lead_wing_root_wing_cross_section_movement
    del lead_wing_right_wing_cross_section_movement
    del right_wing_left_wing_cross_section_movement
    del right_wing_root_wing_cross_section_movement
    del right_wing_right_wing_cross_section_movement
    del left_wing_left_wing_cross_section_movement
    del left_wing_root_wing_cross_section_movement
    del left_wing_right_wing_cross_section_movement

    airplane_movement = src.movement.AirplaneMovement(
        base_airplane=example_airplane,
        wing_movements=[lead_wing_movement, right_wing_movement, left_wing_movement],
    )

    del lead_wing_movement
    del right_wing_movement
    del left_wing_movement

    operating_point_movement = src.movement.OperatingPointMovement(
        base_operating_point=example_operating_point,
    )

    movement = src.movement.Movement(
        airplane_movement=airplane_movement,
        operating_point_movement=operating_point_movement,
        num_steps=None,
        delta_time=None,
    )

    del airplane_movement
    del operating_point_movement

    unsteady_problem = src.problems.UnsteadyProblem(
        movement=movement,
    )

    unsteady_solver = (
        src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_problem=unsteady_problem,
        )
    )

    del unsteady_problem

    unsteady_solver.run(
        logging_level="Warning",
        prescribed_wake=True,
    )

    src.output.animate(
        unsteady_solver=unsteady_solver,
        show_delta_pressures=True,
        show_wake_vortices=True,
        keep_file=False,
    )

    src.output.plot_results_versus_time(
        unsteady_solver=unsteady_solver,
        testing=False,
    )
