"""This script is a single-run version of simulation in unsteady_benchmark_timed.py.
It is useful for profiling the unsteady solver, which cannot be done with
unsteady_benchmark_timed.py. This script doesn't have any expected output images in
the docs directory. Do not commit any changes to this file.
"""

import pterasoftware as ps

flapping_frequency = 1
num_chordwise_panels = 5
num_spanwise_panels = 20

example_airplane = ps.geometry.airplane.Airplane(
    name="Example Airplane",
    wings=[
        ps.geometry.wing.Wing(
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    chord=1.75,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0000",
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    x_le=0.625,
                    y_le=5.0,
                    chord=0.5,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0000",
                    ),
                ),
            ],
        ),
    ],
)

upper_wing_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
)

upper_wing_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
    sweeping_amplitude=15.0,
    sweeping_period=1 / flapping_frequency,
    sweeping_spacing="sine",
    pitching_amplitude=5.0,
    pitching_period=1 / flapping_frequency,
    pitching_spacing="sine",
    heaving_amplitude=5.0,
    heaving_period=1 / flapping_frequency,
    heaving_spacing="sine",
)

upper_wing_movement = ps.movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_sections_movements=[
        upper_wing_root_wing_cross_section_movement,
        upper_wing_tip_wing_cross_section_movement,
    ],
)

del upper_wing_root_wing_cross_section_movement
del upper_wing_tip_wing_cross_section_movement

airplane_movement = ps.movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[upper_wing_movement],
)

del upper_wing_movement

example_operating_point = ps.operating_point.OperatingPoint(
    density=1.225,
    beta=0.0,
    velocity=10.0,
    alpha=0.0,
)

operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=example_operating_point,
)

movement = ps.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
)

del airplane_movement
del operating_point_movement

example_problem = ps.problems.UnsteadyProblem(
    movement=movement, only_final_results=True
)

example_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)

del example_problem

example_solver.run(
    prescribed_wake=True,
    calculate_streamlines=False,
)
