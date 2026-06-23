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
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=14,
                chordwise_spacing="cosine",
            )
        ],
        name="Steady Validation Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
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
                        Lp_Wcsp_Lpp=(1.0, 5.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
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
                Ler_Gs_Cgs=(5.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=12,
                chordwise_spacing="uniform",
            ),
        ],
        name="Multiple Wing Steady Validation Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
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
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=7,
                chordwise_spacing="uniform",
            ),
        ],
        name="Symmetric Unsteady Validation Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return symmetric_unsteady_validation_airplane


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
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
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
                Ler_Gs_Cgs=(6.25, 0.0, 1.75),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
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
                        Lp_Wcsp_Lpp=(0.25, 1.5, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type=None,
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Vertical Stabilizer",
                Ler_Gs_Cgs=(6.25, 0.0, 0.125),
                angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=None,
                symmetryPoint_G_Cg=None,
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
        ],
        name="Symmetric Multiple Wing Unsteady Validation Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return symmetric_multiple_wing_steady_validation_airplane


def make_simple_glider_airplane():
    """This function creates the simple glider Airplane used for free flight testing.

    The simple glider is a three-wing aircraft (cambered main wing, symmetric
    horizontal stabilizer with negative incidence, and a single vertical stabilizer)
    whose planform geometry, center of gravity, and inertia were tuned for static pitch
    and yaw stability and verified in XFLR5. The negative horizontal stabilizer
    incidence relative to the main wing provides the restoring pitch moment that makes
    the glider statically stable, so that it flies a bounded, damped free flight
    trajectory rather than diverging. The center of gravity is left at the geometry
    origin (the default), which is the reference point for the matching inertia matrix
    in make_simple_glider_free_flight_problem.

    The mesh densities here are coarser than the converged values from the original
    convergence study. The static stability that this fixture relies on is a property of
    the continuous planform, center of gravity, and inertia, not of the discretization,
    so the coarser mesh keeps the free flight integration test affordable while
    preserving the stable behavior.

    :return simple_glider_airplane: Airplane
        This is the simple glider Airplane fixture.
    """
    simple_glider_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                        num_spanwise_panels=10,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                        num_spanwise_panels=None,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=4,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=6,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=None,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        spanwise_spacing=None,
                    ),
                ],
                name="Horizontal Stabilizer",
                Ler_Gs_Cgs=(5.0, 0.0, 0.5),
                angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=4,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=6,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=None,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 2.0, 0.0),
                        spanwise_spacing=None,
                    ),
                ],
                name="Vertical Stabilizer",
                Ler_Gs_Cgs=(5.0, 0.0, 1.0),
                angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=None,
                symmetryPoint_G_Cg=None,
                num_chordwise_panels=4,
                chordwise_spacing="uniform",
            ),
        ],
        name="Simple Glider",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=420.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return simple_glider_airplane


def make_flapping_free_flight_airplane():
    """This function creates the flapping-wing Airplane used for free flight testing.

    This is the same airframe as the flapping-wing free flight example: a cambered main
    wing (which flaps symmetrically, defined in the matching movement fixture) plus a
    symmetric V-tail. Unlike the simple glider, this airframe is not tuned for static
    stability or trim; it exists to exercise the strongly coupled free flight solver
    under the large, oscillatory loads that flapping produces.

    The main wing's root is offset from the symmetry plane (its Ler_Gs_Cgs has a nonzero
    y component), so it has type 5 symmetry and is split into two separate, mirrored
    Wings; the matching movement fixture flaps both halves with the same amplitude to
    keep the flapping symmetric. The V-tail's root lies on the symmetry plane, so it
    stays a single mirrored Wing.

    The mesh densities here are coarser than the example's, and the airfoils use the
    default resampling rather than the example's dense override, which keeps the strongly
    coupled flapping integration test affordable. The coupling behavior under test is a
    property of the continuous geometry and mass properties, not of the discretization.

    :return flapping_free_flight_airplane: Airplane
        This is the flapping-wing Airplane fixture.
    """
    flapping_free_flight_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                        num_spanwise_panels=4,
                        chord=1.75,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                        num_spanwise_panels=None,
                        chord=1.5,
                        Lp_Wcsp_Lpp=(0.75, 6.0, 1.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.5, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=4,
                        chord=1.5,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        spanwise_spacing="uniform",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        num_spanwise_panels=None,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.5, 2.0, 1.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        spanwise_spacing=None,
                    ),
                ],
                name="V-Tail",
                Ler_Gs_Cgs=(5.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
            ),
        ],
        name="Flapping Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=420.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return flapping_free_flight_airplane


def make_surface_effect_airplane():
    """This function creates a simple single-wing Airplane for surface effect testing.

    The Airplane uses a NACA 0010 symmetric airfoil with zero twist and zero dihedral.

    :return surface_effect_airplane: Airplane
        This is the Airplane fixture.
    """
    surface_effect_airplane = ps.geometry.airplane.Airplane(
        wings=[
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
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
            )
        ],
        name="Surface Effect Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return surface_effect_airplane
