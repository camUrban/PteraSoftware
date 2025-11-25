"""This script is for running a prescribed unsteady simulation on the simple glider
after finding its converged and trimmed parameters."""

import pterasoftware as ps

converged_prescribed_wake = True
converged_num_chords = 13
converged_num_chordwise_panels = 7
wing_1_converged_num_spanwise_panels = 12
wing_2_converged_num_spanwise_panels = 2
wing_3_converged_num_spanwise_panels = 5

trim_vCg__E = 12.9
trim_alpha = 3.3
trim_beta = 0.0
trim_externalFX_W = 6.5

simple_glider_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=wing_1_converged_num_spanwise_panels,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca2412",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            chordwise_spacing="uniform",
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=wing_2_converged_num_spanwise_panels,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 0.5),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=wing_3_converged_num_spanwise_panels,
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                    ),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 2.0, 0.0),
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 1.0),
            angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
    ],
    weight=420,
)

simple_glider_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=simple_glider_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[2],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

del simple_glider_airplane

simple_glider_operating_point = ps.operating_point.OperatingPoint(
    vCg__E=trim_vCg__E,
    alpha=trim_alpha,
    beta=trim_beta,
    externalFX_W=trim_externalFX_W,
)

simple_glider_operating_point_movement = (
    ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=simple_glider_operating_point
    )
)

del simple_glider_operating_point

simple_glider_movement = ps.movements.movement.Movement(
    airplane_movements=[simple_glider_airplane_movement],
    operating_point_movement=simple_glider_operating_point_movement,
    num_chords=converged_num_chords,
)

del simple_glider_airplane_movement
del simple_glider_operating_point_movement

simple_glider_problem = ps.problems.UnsteadyProblem(
    movement=simple_glider_movement,
    only_final_results=False,
)

print(f"Delta time: {simple_glider_problem.delta_time:.6f} s")
print(f"Converged number of time steps: {simple_glider_problem.num_steps}")

del simple_glider_movement

simple_glider_solver = (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_problem=simple_glider_problem,
    )
)

del simple_glider_problem

simple_glider_solver.run(
    logging_level="Warning",
    prescribed_wake=converged_prescribed_wake,
    calculate_streamlines=False,
)

ps.output.draw(
    solver=simple_glider_solver,
    scalar_type="lift",
    show_streamlines=False,
    show_wake_vortices=True,
    save=False,
    testing=False,
)

ps.output.plot_results_versus_time(
    unsteady_solver=simple_glider_solver, show=True, save=False
)

# delta_time = 0.011074
# converged_num_steps = 91
