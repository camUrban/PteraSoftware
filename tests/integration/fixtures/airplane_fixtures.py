"""This module creates airplane objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_steady_validation_airplane: This function creates an airplane object to be
    used as a fixture for testing steady solvers.

    make_multiple_wing_steady_validation_airplane: This function creates a multi-wing
    airplane object to be used as a fixture for testing steady solvers.

    make_symmetric_unsteady_validation_airplane: This function creates a symmetric
    airplane object to be used as a fixture for testing unsteady solvers.

    make_symmetric_multiple_wing_unsteady_validation_airplane: This function creates
    a multi-wing, symmetric airplane object to be used as a fixture for testing
    unsteady solvers.
"""

import numpy as np

import pterasoftware as ps


def make_steady_validation_airplane():
    """This function creates an airplane object to be used as a fixture for testing
    steady solvers.

    The parameters of this airplane were found to be converged based on the following
    call to analyze_steady_convergence:
    converged_parameters = ps.convergence.analyze_steady_convergence(
        ref_problem=steady_validation_problem,
        solver_type="steady horseshoe vortex lattice method",
        panel_aspect_ratio_bounds=(4, 1),
        num_chordwise_panels_bounds=(3, 20),
        convergence_criteria=0.1,
    ).

    :return steady_validation_airplane: Airplane
        This is the airplane fixture.
    """
    # Create and return the airplane object.
    steady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=1.0,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=20,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=1.0,
                        y_le=5.0,
                        z_le=0.0,
                        chord=0.75,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=5.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Main Wing",
                x_le=0.0,
                y_le=0.0,
                z_le=0.0,
                unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                symmetric=True,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=14,
                chordwise_spacing="cosine",
            )
        ],
        name="Steady Validation Airplane",
        x_ref=0.0,
        y_ref=0.0,
        z_ref=0.0,
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return steady_validation_airplane


def make_multiple_wing_steady_validation_airplane():
    """This function creates a multi-wing airplane object to be used as a fixture
    for testing steady solvers.

    The parameters of this airplane were found to be converged based on the following
    call to analyze_steady_convergence:
    converged_parameters = ps.convergence.analyze_steady_convergence(
        ref_problem=steady_validation_problem,
        solver_type="steady horseshoe vortex lattice method",
        panel_aspect_ratio_bounds=(4, 1),
        num_chordwise_panels_bounds=(3, 20),
        convergence_criteria=0.1,
    ).

    :return multiple_wing_steady_validation_airplane: Airplane
        This is the airplane fixture.
    """
    # Create and return the airplane object.
    multiple_wing_steady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca23012",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=1.0,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=69,
                        spanwise_spacing="uniform",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca23012",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=1.0,
                        y_le=5.0,
                        z_le=0.0,
                        chord=0.75,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=69,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Main Wing",
                x_le=0.0,
                y_le=0.0,
                z_le=0.0,
                unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                symmetric=True,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=12,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=1.00,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=-5.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=16,
                        spanwise_spacing="uniform",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=1.0,
                        y_le=1.0,
                        z_le=0.0,
                        chord=0.75,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=-5.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=16,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Horizontal Stabilizer",
                x_le=5.0,
                y_le=0.0,
                z_le=0.0,
                unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                symmetric=True,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=12,
                chordwise_spacing="uniform",
            ),
        ],
        name="Multiple Wing Steady Validation Airplane",
        x_ref=0.0,
        y_ref=0.0,
        z_ref=0.0,
        weight=1 * 9.81,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return multiple_wing_steady_validation_airplane


def make_symmetric_unsteady_validation_airplane():
    """This function creates a symmetric airplane object to be used as a fixture for
    testing unsteady solvers.

    The parameters of this airplane were found to be converged based on the following
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
        This is the airplane fixture.
    """
    # Create and return the airplane object.
    symmetric_unsteady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=2.0,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=18,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=5.0,
                        z_le=0.0,
                        chord=2.0,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=18,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Main Wing",
                x_le=0.0,
                y_le=0.0,
                z_le=0.0,
                unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                symmetric=True,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=7,
                chordwise_spacing="uniform",
            ),
        ],
        name="Symmetric Unsteady Validation Airplane",
        x_ref=0.0,
        y_ref=0.0,
        z_ref=0.0,
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return symmetric_unsteady_validation_airplane


# TODO: Check that this test case has converged characteristics.
def make_symmetric_multiple_wing_unsteady_validation_airplane():
    """This function creates a multi-wing, symmetric airplane object to be used as a
    fixture for testing unsteady solvers.

    :return symmetric_multiple_wing_steady_validation_airplane: Airplane
        This is the airplane fixture.
    """
    # Create and return the airplane object.
    symmetric_multiple_wing_steady_validation_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=1.5,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.5,
                        y_le=5.0,
                        z_le=0.0,
                        chord=1.0,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Main Wing",
                x_le=0.0,
                y_le=0.0,
                z_le=0.0,
                unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                symmetric=True,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=1.0,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=-5.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.25,
                        y_le=1.5,
                        z_le=0.0,
                        chord=0.75,
                        unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                        twist=-5.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Horizontal Stabilizer",
                x_le=6.25,
                y_le=0.0,
                z_le=1.75,
                unit_normal_vector=np.array([0.0, 1.0, 0.0]),
                symmetric=True,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.0,
                        y_le=0.0,
                        z_le=0.0,
                        chord=1.0,
                        unit_normal_vector=np.array([0.0, 0.0, 1.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            coordinates=None,
                            repanel=True,
                            n_points_per_side=50,
                        ),
                        x_le=0.25,
                        y_le=0.0,
                        z_le=1.5,
                        chord=0.75,
                        unit_normal_vector=np.array([0.0, 0.0, 1.0]),
                        twist=0.0,
                        control_surface_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        num_spanwise_panels=8,
                        spanwise_spacing="cosine",
                    ),
                ],
                name="Vertical Stabilizer",
                x_le=6.25,
                y_le=0.0,
                z_le=0.125,
                unit_normal_vector=np.array([0.0, 0.0, 1.0]),
                symmetric=False,
                unit_chordwise_vector=np.array([1.0, 0.0, 0.0]),
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
            ),
        ],
        name="Symmetric Multiple Wing Unsteady Validation Airplane",
        x_ref=0.0,
        y_ref=0.0,
        z_ref=0.0,
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )
    return symmetric_multiple_wing_steady_validation_airplane
