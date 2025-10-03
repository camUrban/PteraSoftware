import pterasoftware as ps

solver_type = "unsteady_ring"

flat_plate_airplane = ps.geometry.airplane.Airplane(
    name="Flat Plate Airplane",
    Cgi_E_I=(0.0, 0.0, 0.0),
    angles_E_to_B_izyx=(0.0, 0.0, 0.0),
    weight=0.0,
    s_ref=None,
    c_ref=None,
    b_ref=None,
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=5,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing="uniform",
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    num_spanwise_panels=None,
                    chord=1.0,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    control_surface_hinge_point=0.75,
                    control_surface_deflection=0.0,
                    spanwise_spacing=None,
                    airfoil=ps.geometry.airfoil.Airfoil(
                        name="naca0012",
                        outline_A_lp=None,
                        resample=True,
                        n_points_per_side=400,
                    ),
                ),
            ],
            name="Right Wing",
            prelimLer_G_Cg=(-0.5, 1.0, 0.0),
            angles_G_to_prelimWn_izyx=(0, 0, 0),
            symmetric=True,
            mirror_only=False,
            symmetry_normal_G=(0, 1, 0),
            symmetry_point_G_Cg=(0, -1.0, 0),
            num_chordwise_panels=5,
            chordwise_spacing="uniform",
        ),
    ],
)

flat_plate_operating_point = ps.operating_point.OperatingPoint(
    rho=1.225, vCg__E=1.0, alpha=15.0, beta=0.0, nu=15.06e-6
)


def get_solver(this_solver_type, this_airplane, this_operating_point):
    match this_solver_type:
        case "steady_horseshoe":
            this_problem = ps.problems.SteadyProblem(
                airplanes=[this_airplane],
                operating_point=this_operating_point,
            )

            del this_airplane
            del this_operating_point

            return ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                steady_problem=this_problem,
            )
        case "steady_ring":
            this_problem = ps.problems.SteadyProblem(
                airplanes=[this_airplane],
                operating_point=this_operating_point,
            )

            del this_airplane
            del this_operating_point

            return ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
                steady_problem=this_problem,
            )
        case "unsteady_ring":
            this_wing_root_wing_cross_section_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=this_airplane.wings[0].wing_cross_sections[
                        0
                    ],
                    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    ampAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    periodAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
                    phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                )
            )

            this_wing_tip_wing_cross_section_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=this_airplane.wings[0].wing_cross_sections[
                        1
                    ],
                    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    ampAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    periodAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
                    phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                )
            )

            this_wing_movement = ps.movements.wing_movement.WingMovement(
                base_wing=this_airplane.wings[0],
                wing_cross_section_movements=[
                    this_wing_root_wing_cross_section_movement,
                    this_wing_tip_wing_cross_section_movement,
                ],
                ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
                phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
                ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
                phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            )

            reflected_wing_root_wing_cross_section_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=this_airplane.wings[1].wing_cross_sections[
                        0
                    ],
                    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    ampAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    periodAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
                    phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                )
            )

            reflected_wing_tip_wing_cross_section_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=this_airplane.wings[1].wing_cross_sections[
                        1
                    ],
                    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                    phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    ampAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    periodAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                    spacingAngles_Wcsp_to_Wcs_izyx=("sine", "sine", "sine"),
                    phaseAngles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
                )
            )

            reflected_wing_movement = ps.movements.wing_movement.WingMovement(
                base_wing=this_airplane.wings[1],
                wing_cross_section_movements=[
                    reflected_wing_root_wing_cross_section_movement,
                    reflected_wing_tip_wing_cross_section_movement,
                ],
                ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
                phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
                ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
                phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            )

            del this_wing_root_wing_cross_section_movement
            del this_wing_tip_wing_cross_section_movement

            this_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=this_airplane,
                wing_movements=[this_wing_movement, reflected_wing_movement],
                ampCgi_E_I=(0.0, 0.0, 0.0),
                periodCgi_E_I=(0.0, 0.0, 0.0),
                spacingCgi_E_I=("sine", "sine", "sine"),
                phaseCgi_E_I=(0.0, 0.0, 0.0),
                ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
                periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
                spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
                phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            )

            del this_airplane
            del this_wing_movement

            this_operating_point_movement = (
                ps.movements.operating_point_movement.OperatingPointMovement(
                    base_operating_point=this_operating_point,
                    ampVCg__E=0.0,
                    periodVCg__E=0.0,
                    spacingVCg__E="sine",
                    phaseVCg__E=0.0,
                )
            )

            del this_operating_point

            this_movement = ps.movements.movement.Movement(
                airplane_movements=[this_airplane_movement],
                operating_point_movement=this_operating_point_movement,
                num_steps=None,
                delta_time=None,
                num_chords=5,
            )

            del this_airplane_movement
            del this_operating_point_movement

            this_problem = ps.problems.UnsteadyProblem(
                movement=this_movement,
                only_final_results=False,
            )

            del this_movement

            return ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                unsteady_problem=this_problem,
            )
        case _:
            raise ValueError(f"Unknown solver type: {this_solver_type}")


flat_plate_solver = get_solver(
    this_solver_type=solver_type,
    this_airplane=flat_plate_airplane,
    this_operating_point=flat_plate_operating_point,
)

flat_plate_solver.run(
    logging_level="Warning",
)

if solver_type in ("steady_horseshoe", "steady_ring"):
    ps.output.draw(
        solver=flat_plate_solver,
        scalar_type="lift",
        show_streamlines=True,
        show_wake_vortices=False,
        save=False,
    )

    ps.output.print_results()
else:
    ps.output.draw(
        solver=flat_plate_solver,
        scalar_type="lift",
        show_streamlines=False,
        show_wake_vortices=True,
        save=False,
    )

    ps.output.animate(
        unsteady_solver=flat_plate_solver,
        scalar_type="lift",
        show_wake_vortices=True,
        save=False,
        testing=False,
    )
