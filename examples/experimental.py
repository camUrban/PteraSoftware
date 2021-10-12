"""This is script is an example of running Ptera Software's unsteady ring vortex
lattice method solver on an airplane with geometry similar to the NMT experimental
setup. """

import src

example_airplane = src.geometry.Airplane(
    name="Experimental Airplane",
    wings=[
        src.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                src.geometry.WingCrossSection(
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
                ),
                src.geometry.WingCrossSection(
                    y_le=0.2275,
                    chord=0.1094,
                    airfoil=src.geometry.Airfoil(name="naca0012"),
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

airplane_movement = src.movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[
        src.movement.WingMovement(
            base_wing=example_airplane.wings[0],
            wing_cross_sections_movements=[
                src.movement.WingCrossSectionMovement(
                    example_airplane.wings[0].wing_cross_sections[0],
                ),
                src.movement.WingCrossSectionMovement(
                    example_airplane.wings[0].wing_cross_sections[1],
                    sweeping_amplitude=30.0,
                    sweeping_period=1.0,
                ),
                src.movement.WingCrossSectionMovement(
                    example_airplane.wings[0].wing_cross_sections[2],
                    sweeping_amplitude=30.0,
                    sweeping_period=1.0,
                ),
            ],
        )
    ],
)

example_operating_point = src.operating_point.OperatingPoint(
    density=1.225,
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
)

del airplane_movement
del operating_point_movement

example_problem = src.problems.UnsteadyProblem(movement=movement)

example_solver = (
    src.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)

del example_problem

example_solver.run()

src.output.animate(
    unsteady_solver=example_solver,
    show_delta_pressures=True,
    show_wake_vortices=True,
    keep_file=False,
)

src.output.plot_results_versus_time(unsteady_solver=example_solver)
