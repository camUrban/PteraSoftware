"""This script is used to benchmark the speed of the unsteady solver with a typical
use case. Avoid committing changes to this file."""

import timeit

import numpy as np

n_repeat = 3
n_execute = 1

print("\tBenchmarking unsteady solver...")

setup = """
import pterasoftware as ps

num_chordwise_panels = 5
num_spanwise_panels = 6

flapping_frequency = 1

benchmark_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=num_spanwise_panels,
                    spanwise_spacing="uniform",
                    chord=1.75,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    control_surface_symmetry_type="symmetric",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    spanwise_spacing=None,
                    Lp_Wcsp_Lpp=(0.625, 5.0, 0.0),
                    chord=0.5,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            name="Main Wing",
            Ler_Gs_Cgs=(0.0, 0.25, 0.0),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
        ),
    ],
    name="Unsteady Benchmark Airplane",
)

main_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=benchmark_airplane.wings[0].wing_cross_sections[0]
    )
)
main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=benchmark_airplane.wings[0].wing_cross_sections[1],
    )
)
reflected_main_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=benchmark_airplane.wings[1].wing_cross_sections[0]
    )
)
reflected_main_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=benchmark_airplane.wings[1].wing_cross_sections[1],
    )
)

main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=benchmark_airplane.wings[0],
    wing_cross_section_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    ampAngles_Gs_to_Wn_ixyz=(15.0, 5.0, 5.0),
    periodAngles_Gs_to_Wn_ixyz=(
        1 / flapping_frequency,
        1 / flapping_frequency,
        1 / flapping_frequency,
    ),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
)
reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=benchmark_airplane.wings[1],
    wing_cross_section_movements=[
        reflected_main_wing_root_wing_cross_section_movement,
        reflected_main_wing_tip_wing_cross_section_movement,
    ],
    ampAngles_Gs_to_Wn_ixyz=(15.0, 5.0, 5.0),
    periodAngles_Gs_to_Wn_ixyz=(
        1 / flapping_frequency,
        1 / flapping_frequency,
        1 / flapping_frequency,
    ),
    spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
)

del main_wing_root_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement
del reflected_main_wing_root_wing_cross_section_movement
del reflected_main_wing_tip_wing_cross_section_movement

benchmark_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=benchmark_airplane,
    wing_movements=[main_wing_movement, reflected_main_wing_movement],
)

del main_wing_movement
del reflected_main_wing_movement

benchmark_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0
)

benchmark_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=benchmark_operating_point
    )
)

benchmark_movement = ps.movements.movement.Movement(
    airplane_movements=[benchmark_airplane_movement],
    operating_point_movement=benchmark_operating_point_movement,
    num_cycles=3,
)

del benchmark_airplane_movement
del benchmark_operating_point_movement

benchmark_problem = ps.problems.UnsteadyProblem(
    movement=benchmark_movement, only_final_results=True
)

benchmark_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=benchmark_problem,
    )
)

del benchmark_problem
"""
statement = """
benchmark_solver.run(
    prescribed_wake=True,
    calculate_streamlines=False,
)
"""

times = timeit.repeat(repeat=n_repeat, stmt=statement, setup=setup, number=n_execute)
best_time = min(times) / n_execute
best_time_pretty = np.format_float_scientific(best_time, 2)
print("\t\tAverage Time per Loop: " + best_time_pretty + " s")
