"""This is script is an example of how to run Ptera Software's unsteady ring vortex
lattice method solver on a custom airplane with variable geometry. """

import pterasoftware as ps
import numpy as np

example_airplane = ps.geometry.Airplane(
    wings=[
        ps.geometry.Wing(
            name="Caudal Fin",
            symmetric=True,
            num_chordwise_panels=3,
            chordwise_spacing="uniform",
            symmetry_unit_normal_vector=np.array([0, 1, 0]),
            unit_chordwise_vector=np.array([1, 0, 0]),
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=0.0,
                    z_le=0.0,
                    twist=0.0,
                    unit_normal_vector=np.array([0, 1, 0]),
                    num_spanwise_panels=3,
                    spanwise_spacing="cosine",
                    chord=1.0,
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
                ps.geometry.WingCrossSection(
                    x_le=0.0,
                    y_le=5.0,
                    z_le=0.0,
                    chord=1.0,
                    twist=0.0,
                    unit_normal_vector=np.array([0, 1, 0]),
                    airfoil=ps.geometry.Airfoil(
                        name="naca0012",
                    ),
                ),
            ],
        ),
    ],
)

main_wing_root_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[0],
    sweeping_amplitude=0.0,
    sweeping_period=0.0,
    sweeping_spacing="sine",
    pitching_amplitude=0.0,
    pitching_period=0.0,
    pitching_spacing="sine",
    heaving_amplitude=0.0,
    heaving_period=0.0,
    heaving_spacing="sine",
)

main_wing_tip_wing_cross_section_movement = ps.movement.WingCrossSectionMovement(
    base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[1],
    sweeping_amplitude=0.0,
    sweeping_period=0.0,
    sweeping_spacing="sine",
    pitching_amplitude=0.0,
    pitching_period=0.0,
    pitching_spacing="sine",
    heaving_amplitude=0.0,
    heaving_period=0.0,
    heaving_spacing="sine",
)

main_wing_movement = ps.movement.WingMovement(
    base_wing=example_airplane.wings[0],
    wing_cross_sections_movements=[
        main_wing_root_wing_cross_section_movement,
        main_wing_tip_wing_cross_section_movement,
    ],
    x_le_amplitude=0.0,
    x_le_period=0.0,
    x_le_spacing="sine",
    y_le_amplitude=0.0,
    y_le_period=0.0,
    y_le_spacing="sine",
    z_le_amplitude=0.0,
    z_le_period=0.0,
    z_le_spacing="sine",
)

del main_wing_root_wing_cross_section_movement
del main_wing_tip_wing_cross_section_movement

airplane_movement = ps.movement.AirplaneMovement(
    base_airplane=example_airplane,
    wing_movements=[main_wing_movement],
    x_ref_amplitude=0.0,
    x_ref_period=0.0,
    x_ref_spacing="sine",
    y_ref_amplitude=0.0,
    y_ref_period=0.0,
    y_ref_spacing="sine",
    z_ref_amplitude=0.0,
    z_ref_period=0.0,
    z_ref_spacing="sine",
)

del main_wing_movement

example_operating_point = ps.operating_point.OperatingPoint(
    density=1.225,
    beta=0.0,
    velocity=10.0,
    alpha=0.0,
    nu=15.06e-6,
)

operating_point_movement = ps.movement.OperatingPointMovement(
    base_operating_point=example_operating_point,
    velocity_amplitude=0.0,
    velocity_period=0.0,
    velocity_spacing="sine",
)

movement = ps.movement.Movement(
    airplane_movements=[airplane_movement],
    operating_point_movement=operating_point_movement,
    num_steps=None,
    delta_time=None,
)

del airplane_movement
del operating_point_movement

example_problem = ps.problems.UnsteadyProblem(
    movement=movement,
)

example_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=example_problem,
    )
)

del example_problem

example_solver.run(
    logging_level="Warning",
    prescribed_wake=True,
)

# ps.output.animate(
#     unsteady_solver=example_solver,
#     scalar_type="lift",
#     show_wake_vortices=True,
#     save=False,
# )

ps.output.print_unsteady_results(unsteady_solver=example_solver)

ps.output.draw(
    solver=example_solver,
    scalar_type="lift",
    show_wake_vortices=True,
    save=False,
)

# ps.output.plot_results_versus_time(
#     unsteady_solver=example_solver,
#     show=True,
#     save=False,
# )
