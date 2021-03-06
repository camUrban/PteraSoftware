import timeit

import numpy as np

n_repeat = 3
n_execute = 1

print("\tBenchmarking unsteady solver...")

setup = """
import math

import pterasoftware as ps

flapping_frequency = 1
num_chordwise_panels = 5
num_spanwise_panels = 20

example_airplane = ps.geometry.Airplane(
    name="Example Airplane",
    wings=[
        ps.geometry.Wing(
            name="Main Wing",
            symmetric=True,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    chord=1.75,
                    airfoil=ps.geometry.Airfoil(name="naca0000",),
                ),
                ps.geometry.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    x_le=0.625,
                    y_le=5.0,
                    chord=0.5,
                    airfoil=ps.geometry.Airfoil(name="naca0000",),
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
    base_airplane=example_airplane, wing_movements=[upper_wing_movement],
)

del upper_wing_movement

example_operating_point = ps.operating_point.OperatingPoint(
    density=1.225, beta=0.0, velocity=10.0, alpha=0.0,
)

operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=example_operating_point,
)

movement = ps.movement.Movement(
    airplane_movement=airplane_movement,
    operating_point_movement=operating_point_movement,
)

del airplane_movement
del operating_point_movement

example_problem = ps.problems.UnsteadyProblem(movement=movement, only_final_results=True)

example_solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
    unsteady_problem=example_problem,
)

del example_problem
"""
statement = """
example_solver.run(
    verbose=True, prescribed_wake=True,
)
"""

times = timeit.repeat(repeat=n_repeat, stmt=statement, setup=setup, number=n_execute)
best_time = min(times) / n_execute
best_time_pretty = np.format_float_scientific(best_time, 2)
print("\t\tAverage Time per Loop: " + best_time_pretty + " s")
