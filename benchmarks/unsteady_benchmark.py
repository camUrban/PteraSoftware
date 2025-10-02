"""This script is a single-run version of simulation in unsteady_benchmark_timed.py.
It is useful for profiling the unsteady solver, which cannot be done with
unsteady_benchmark_timed.py. This script doesn't have any expected output images in
the docs directory. Do not commit any changes to this file."""

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
            symmetry_normal_Wn=(0.0, 1.0, 0.0),
            symmetry_point_Wn_Ler=(0.0, 0.0, 0.0),
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing="uniform",
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
        ),
    ],
)

upper_wing_root_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
    )
)

upper_wing_tip_wing_cross_section_movement = (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
        ampAngles_Wcsp_to_Wcs_izyx=(15.0, 5.0, 5.0),
        periodAngles_Wcsp_to_Wcs_izyx=(
            1 / flapping_frequency,
            1 / flapping_frequency,
            1 / flapping_frequency,
        ),
        spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
    )
)

upper_wing_movement = ps.movements.wing_movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_section_movements=[
        upper_wing_root_wing_cross_section_movement,
        upper_wing_tip_wing_cross_section_movement,
    ],
)

del upper_wing_root_wing_cross_section_movement
del upper_wing_tip_wing_cross_section_movement

airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[upper_wing_movement],
)

del upper_wing_movement

example_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=10.0, alpha=0.0, beta=0.0
)

operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
    base_operating_point=example_operating_point
)

movement = ps.movements.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    num_cycles=3,
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
