"""This is script is an example of how to run Ptera Software's unsteady ring vortex
lattice method solver on a custom airplane with variable geometry."""

import src

steady = True

lead_airplane = src.geometry.Airplane(
    name="Lead Airplane",
    x_ref=0.0,
    y_ref=0.0,
    z_ref=0.0,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
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

right_airplane = src.geometry.Airplane(
    name="Right Airplane",
    x_ref=-15.0,
    y_ref=10.0,
    z_ref=0.0,
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=6,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
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
    steady_solver = src.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem=src.problems.SteadyProblem(
            airplanes=[lead_airplane, right_airplane],
            operating_point=example_operating_point,
        ),
    )

    steady_solver.run()

    src.output.draw(
        solver=steady_solver,
        show_delta_pressures=True,
    )
else:
    lead_airplane_main_wing_root_wing_cross_section_movement = (
        src.movement.WingCrossSectionMovement(
            base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[0]
        )
    )
    lead_airplane_main_wing_tip_wing_cross_section_movement = (
        src.movement.WingCrossSectionMovement(
            base_wing_cross_section=lead_airplane.wings[0].wing_cross_sections[1],
            sweeping_amplitude=30.0,
            sweeping_period=1.0,
            sweeping_spacing="sine",
        )
    )

    right_airplane_main_wing_root_wing_cross_section_movement = (
        src.movement.WingCrossSectionMovement(
            base_wing_cross_section=right_airplane.wings[0].wing_cross_sections[0]
        )
    )
    right_airplane_main_wing_tip_wing_cross_section_movement = (
        src.movement.WingCrossSectionMovement(
            base_wing_cross_section=right_airplane.wings[0].wing_cross_sections[1],
            sweeping_amplitude=30.0,
            sweeping_period=1.0,
            sweeping_spacing="sine",
        )
    )

    lead_airplane_main_wing_movement = src.movement.WingMovement(
        base_wing=lead_airplane.wings[0],
        wing_cross_sections_movements=[
            lead_airplane_main_wing_root_wing_cross_section_movement,
            lead_airplane_main_wing_tip_wing_cross_section_movement,
        ],
    )
    right_airplane_main_wing_movement = src.movement.WingMovement(
        base_wing=right_airplane.wings[0],
        wing_cross_sections_movements=[
            right_airplane_main_wing_root_wing_cross_section_movement,
            right_airplane_main_wing_tip_wing_cross_section_movement,
        ],
    )

    del lead_airplane_main_wing_root_wing_cross_section_movement
    del lead_airplane_main_wing_tip_wing_cross_section_movement
    del right_airplane_main_wing_root_wing_cross_section_movement
    del right_airplane_main_wing_tip_wing_cross_section_movement

    lead_airplane_movement = src.movement.AirplaneMovement(
        base_airplane=lead_airplane,
        wing_movements=[lead_airplane_main_wing_movement],
    )
    right_airplane_movement = src.movement.AirplaneMovement(
        base_airplane=lead_airplane,
        wing_movements=[lead_airplane_main_wing_movement],
    )

    del lead_airplane_main_wing_movement
    del right_airplane_main_wing_movement

    operating_point_movement = src.movement.OperatingPointMovement(
        base_operating_point=example_operating_point,
    )

    movement = src.movement.Movement(
        airplane_movements=[lead_airplane_movement, right_airplane_movement],
        operating_point_movement=operating_point_movement,
        num_steps=None,
        delta_time=None,
    )

    del lead_airplane_movement
    del right_airplane_movement

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
