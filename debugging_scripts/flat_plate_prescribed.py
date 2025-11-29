"""This is a debugging script for running a prescribed unsteady simulation on a flat
plate with one Panel."""

import pterasoftware as ps

trim_alpha = 5.0
trim_beta = 0.0

flat_plate_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=1,
                    control_surface_symmetry_type=None,
                    spanwise_spacing="uniform",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type=None,
                    spanwise_spacing=None,
                ),
            ],
            chordwise_spacing="uniform",
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=1,
        ),
    ],
    weight=100,
)

flat_plate_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=flat_plate_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=flat_plate_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=flat_plate_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=flat_plate_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

del flat_plate_airplane

flat_plate_operating_point = ps.operating_point.OperatingPoint(
    vCg__E=10.0,
    alpha=trim_alpha,
    beta=trim_beta,
    externalFX_W=0.0,
)

flat_plate_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=flat_plate_operating_point
    )
)

del flat_plate_operating_point

flat_plate_movement = ps.movements.movement.Movement(
    airplane_movements=[flat_plate_airplane_movement],
    operating_point_movement=flat_plate_operating_point_movement,
    delta_time=0.1,
    num_steps=10,
)

del flat_plate_airplane_movement
del flat_plate_operating_point_movement

flat_plate_problem = ps.problems.UnsteadyProblem(
    movement=flat_plate_movement,
    only_final_results=False,
)

del flat_plate_movement

flat_plate_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=flat_plate_problem,
    )
)

del flat_plate_problem

flat_plate_solver.run(
    logging_level="Info",
    prescribed_wake=True,
    calculate_streamlines=False,
)

ps.output.draw(
    solver=flat_plate_solver,
    scalar_type="lift",
    show_streamlines=False,
    show_wake_vortices=True,
    save=False,
    testing=False,
)

ps.output.plot_results_versus_time(
    unsteady_solver=flat_plate_solver, show=True, save=False
)

# delta_time = 0.011074
# converged_num_steps = 91
