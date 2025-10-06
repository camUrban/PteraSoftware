"""This module creates Airplanes to be used as fixtures."""

import pterasoftware as ps


def make_steady_validation_airplane():
    """This function creates an Airplane to be used as a fixture for testing steady
    solvers.

    The parameters of this Airplane were found to be converged based on the following
    call to analyze_steady_convergence:
    converged_parameters = ps.convergence.analyze_steady_convergence(
        ref_problem=steady_validation_problem,
        solver_type="steady horseshoe vortex lattice method",
        panel_aspect_ratio_bounds=(4, 1),
        num_chordwise_panels_bounds=(3, 20),
        convergence_criteria=0.1,
    ).

    :return steady_validation_airplane: Airplane
        This is the Airplane fixture.
    """
    steady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=20,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(1.0, 5.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                prelimLer_G_Cg=(0.0, 0.0, 0.0),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=14,
                chordwise_spacing="cosine",
            )
        ],
        name="Steady Validation Airplane",
        Cg_E_CgP1=(0.0, 0.0, 0.0),
        angles_E_to_B_izyx=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return steady_validation_airplane


def make_multiple_wing_steady_validation_airplane():
    """This function creates an Airplane with multiple Wings to be used as a fixture
    for testing steady solvers.

    The parameters of this Airplane were found to be converged based on the following
    call to analyze_steady_convergence:
    converged_parameters = ps.convergence.analyze_steady_convergence(
        ref_problem=steady_validation_problem,
        solver_type="steady horseshoe vortex lattice method",
        panel_aspect_ratio_bounds=(4, 1),
        num_chordwise_panels_bounds=(3, 20),
        convergence_criteria=0.1,
    ).

    :return multiple_wing_steady_validation_airplane: Airplane
        This is the Airplane fixture.
    """
    multiple_wing_steady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca23012",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=69,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="uniform",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca23012",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(1.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                prelimLer_G_Cg=(0.0, 0.0, 0.0),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=12,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=16,
                        chord=1.00,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="uniform",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(1.0, 1.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Horizontal Stabilizer",
                prelimLer_G_Cg=(5.0, 0.0, 0.0),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=12,
                chordwise_spacing="uniform",
            ),
        ],
        name="Multiple Wing Steady Validation Airplane",
        Cg_E_CgP1=(0.0, 0.0, 0.0),
        angles_E_to_B_izyx=(0.0, 0.0, 0.0),
        weight=1 * 9.81,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return multiple_wing_steady_validation_airplane


def make_symmetric_unsteady_validation_airplane():
    """This function creates a symmetric Airplane to be used as a fixture for testing
    unsteady solvers.

    The parameters of this Airplane were found to be converged based on the following
    call to analyze_unsteady_convergence:
    converged_parameters = ps.convergence.analyze_unsteady_convergence(
        ref_problem=unsteady_validation_problem,
        prescribed_wake=True,
        free_wake=True,
        num_chords_bounds=(3, 9),
        panel_aspect_ratio_bounds=(4, 1),
        num_chordwise_panels_bounds=(4, 11),
        coefficient_mask=[True, False, True, False, True, False],
        convergence_criteria=1.0,
    ).

    :return symmetric_unsteady_validation_airplane: Airplane
        This is the Airplane fixture.
    """
    symmetric_unsteady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=18,
                        chord=2.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=2.0,
                        Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                prelimLer_G_Cg=(0.0, 0.0, 0.0),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=7,
                chordwise_spacing="uniform",
            ),
        ],
        name="Symmetric Unsteady Validation Airplane",
        Cg_E_CgP1=(0.0, 0.0, 0.0),
        angles_E_to_B_izyx=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return symmetric_unsteady_validation_airplane


# TODO: Check that this test case has converged characteristics.
def make_symmetric_multiple_wing_unsteady_validation_airplane():
    """This function creates a multi-wing, symmetric Airplane to be used as a fixture
    for testing unsteady solvers.

    :return symmetric_multiple_wing_steady_validation_airplane: Airplane
        This is the Airplane fixture.
    """
    symmetric_multiple_wing_steady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=8,
                        chord=1.5,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.5, 5.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                prelimLer_G_Cg=(0.0, 0.0, 0.0),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=8,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(0.25, 1.5, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Horizontal Stabilizer",
                prelimLer_G_Cg=(6.25, 0.0, 1.75),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=8,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type=None,
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(0.25, 0.0, 1.5),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type=None,
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Vertical Stabilizer",
                prelimLer_G_Cg=(6.25, 0.0, 0.125),
                angles_G_to_prelimWn_ixyz=(0.0, 0.0, 0.0),
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=None,
                symmetryPoint_G_Cg=None,
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
        ],
        name="Symmetric Multiple Wing Unsteady Validation Airplane",
        Cg_E_CgP1=(0.0, 0.0, 0.0),
        angles_E_to_B_izyx=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return symmetric_multiple_wing_steady_validation_airplane
