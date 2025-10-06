"""This module contains functions to create geometry objects for use in tests."""

import numpy as np
import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _panel


def make_test_airfoil_fixture():
    """This method makes a fixture that is an Airfoil for testing purposes.

    :return test_airfoil_fixture: Airfoil
        This is the Airfoil configured for testing.
    """
    test_airfoil_fixture = ps.geometry.airfoil.Airfoil(name="naca2412")

    return test_airfoil_fixture


def make_basic_wing_cross_section_fixture(airfoil=None):
    """This method makes a fixture that is a WingCrossSection with typical
    parameters for general testing.

    :param airfoil: Airfoil, optional
        This is the Airfoil to use for the WingCrossSection. If None, a new
        test airfoil fixture will be created. The default is None.

    :return basic_wing_cross_section_fixture: WingCrossSection
        This is the WingCrossSection configured for general testing.
    """
    # Use provided airfoil or create a new one.
    if airfoil is None:
        test_airfoil_fixture = make_test_airfoil_fixture()
    else:
        test_airfoil_fixture = airfoil

    # Create the basic WingCrossSection.
    basic_wing_cross_section_fixture = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=8,
        chord=1.5,
        Lp_Wcsp_Lpp=[0.2, 0.5, 0.1],
        angles_Wcsp_to_Wcs_ixyz=[5.0, -2.0, 3.0],
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=5.0,
        spanwise_spacing="cosine",
    )

    # Return the WingCrossSection fixture.
    return basic_wing_cross_section_fixture


def make_root_wing_cross_section_fixture():
    """This method makes a fixture that is a root WingCrossSection (with
    zero vectors as required by constraints).

    :return root_wing_cross_section_fixture: WingCrossSection
        This is the root WingCrossSection.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the root WingCrossSection.
    root_wing_cross_section_fixture = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=10,
        chord=2.0,
        Lp_Wcsp_Lpp=[0.0, 0.0, 0.0],
        angles_Wcsp_to_Wcs_ixyz=[0.0, 0.0, 0.0],
    )

    # Return the root WingCrossSection fixture.
    return root_wing_cross_section_fixture


def make_tip_wing_cross_section_fixture():
    """This method makes a fixture that is a tip WingCrossSection (with
    None values for spanwise parameters).

    :return tip_wing_cross_section_fixture: WingCrossSection
        This is the tip WingCrossSection.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the tip WingCrossSection.
    tip_wing_cross_section_fixture = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=None,
        chord=0.8,
        Lp_Wcsp_Lpp=[0.5, 2.0, 0.2],
        angles_Wcsp_to_Wcs_ixyz=[10.0, -5.0, 8.0],
        spanwise_spacing=None,
    )

    # Return the tip WingCrossSection fixture.
    return tip_wing_cross_section_fixture


def make_minimal_wing_cross_section_fixture():
    """This method makes a fixture that is a WingCrossSection with minimal
    valid parameters.

    :return minimal_wing_cross_section_fixture: WingCrossSection
        This is the minimal WingCrossSection.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the minimal WingCrossSection.
    minimal_wing_cross_section_fixture = (
        ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil_fixture,
            num_spanwise_panels=1,  # Minimum valid value
        )
    )

    # Return the minimal WingCrossSection fixture.
    return minimal_wing_cross_section_fixture


def make_asymmetric_control_surface_wing_cross_section_fixture():
    """This method makes a fixture that is a WingCrossSection with asymmetric
    control surface configuration.

    :return asymmetric_wing_cross_section_fixture: WingCrossSection
        This is the WingCrossSection with asymmetric control surface.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the asymmetric control surface WingCrossSection.
    asymmetric_wing_cross_section_fixture = (
        ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil_fixture,
            num_spanwise_panels=None,
            chord=1.2,
            Lp_Wcsp_Lpp=[0.1, 1.0, 0.05],
            angles_Wcsp_to_Wcs_ixyz=[2.0, 0.0, -1.0],
            control_surface_symmetry_type="asymmetric",
            control_surface_hinge_point=0.8,
            control_surface_deflection=-5.0,
            spanwise_spacing=None,
        )
    )

    # Return the asymmetric WingCrossSection fixture.
    return asymmetric_wing_cross_section_fixture


def make_origin_wing_fixture():
    """This method makes a fixture that is a Wing positioned at the origin,
    suitable for movement testing.

    :return origin_wing_fixture: Wing
        This is the Wing positioned at the origin for movement testing.
    """
    # Create WingCrossSections for the wing.
    root_wcs = make_root_wing_cross_section_fixture()
    tip_wcs = make_tip_wing_cross_section_fixture()

    # Create wing at origin (for movement testing).
    origin_wing_fixture = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Origin Wing",
        prelimLer_G_Cg=[0.0, 0.0, 0.0],
        angles_G_to_prelimWn_ixyz=[0.0, 0.0, 0.0],
        symmetric=False,
        mirror_only=False,
        symmetryNormal_G=None,
        symmetryPoint_G_Cg=None,
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )

    return origin_wing_fixture


def make_type_1_wing_fixture():
    """This method makes a fixture that is a Wing with type 1 symmetry
    (symmetric=False, mirror_only=False).

    :return type_1_wing_fixture: Wing
        This is the Wing configured for type 1 symmetry testing.
    """
    # Create WingCrossSections for the wing
    root_wcs = make_root_wing_cross_section_fixture()
    tip_wcs = make_tip_wing_cross_section_fixture()

    # Create type 1 wing (no symmetry, no mirroring)
    type_1_wing_fixture = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Type 1 Test Wing",
        prelimLer_G_Cg=[1.0, 0.0, 0.5],
        angles_G_to_prelimWn_ixyz=[0.0, 5.0, 0.0],
        symmetric=False,
        mirror_only=False,
        symmetryNormal_G=None,
        symmetryPoint_G_Cg=None,
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )

    return type_1_wing_fixture


def make_type_2_wing_fixture():
    """This method makes a fixture that is a Wing with type 2 symmetry
    (symmetric=False, mirror_only=True, coincident_symmetry_plane=True).

    :return type_2_wing_fixture: Wing
        This is the Wing configured for type 2 symmetry testing.
    """
    # Create WingCrossSections for the wing
    root_wcs = make_root_wing_cross_section_fixture()
    tip_wcs = make_tip_wing_cross_section_fixture()

    # Create type 2 wing (mirror_only=True, coincident xz-plane symmetry)
    type_2_wing_fixture = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Type 2 Test Wing",
        prelimLer_G_Cg=[1.0, 0.0, 0.5],
        angles_G_to_prelimWn_ixyz=[0.0, 0.0, 0.0],
        symmetric=False,
        mirror_only=True,
        symmetryNormal_G=[0.0, 1.0, 0.0],
        symmetryPoint_G_Cg=[1.0, 0.0, 0.5],
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )

    return type_2_wing_fixture


def make_type_3_wing_fixture():
    """This method makes a fixture that is a Wing with type 3 symmetry
    (symmetric=False, mirror_only=True, coincident_symmetry_plane=False).

    :return type_3_wing_fixture: Wing
        This is the Wing configured for type 3 symmetry testing.
    """
    # Create WingCrossSections for the wing
    root_wcs = make_root_wing_cross_section_fixture()
    tip_wcs = make_tip_wing_cross_section_fixture()

    # Create type 3 wing (mirror_only=True, non-coincident symmetry plane)
    type_3_wing_fixture = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Type 3 Test Wing",
        prelimLer_G_Cg=[1.0, 0.0, 0.5],
        angles_G_to_prelimWn_ixyz=[0.0, 5.0, 0.0],
        symmetric=False,
        mirror_only=True,
        symmetryNormal_G=[0.0, 0.707, 0.707],
        symmetryPoint_G_Cg=[0.5, 0.0, 0.0],
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )

    return type_3_wing_fixture


def make_type_4_wing_fixture():
    """This method makes a fixture that is a Wing with type 4 symmetry
    (symmetric=True, coincident_symmetry_plane=True).

    :return type_4_wing_fixture: Wing
        This is the Wing configured for type 4 symmetry testing.
    """
    # Create WingCrossSections for the wing with symmetric control surfaces
    root_wcs = make_root_wing_cross_section_fixture()
    root_wcs.control_surface_symmetry_type = "symmetric"

    tip_wcs = make_tip_wing_cross_section_fixture()
    # Set chord and control surface for tip (since it was None)
    tip_wcs.chord = 1.0
    tip_wcs.control_surface_symmetry_type = "symmetric"
    tip_wcs.control_surface_hinge_point = 0.75
    tip_wcs.control_surface_deflection = 2.0

    # Create type 4 wing (symmetric=True, coincident xz-plane symmetry)
    type_4_wing_fixture = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Type 4 Test Wing",
        prelimLer_G_Cg=[1.0, 0.0, 0.5],
        angles_G_to_prelimWn_ixyz=[0.0, 0.0, 0.0],
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=[0.0, 1.0, 0.0],
        symmetryPoint_G_Cg=[1.0, 0.0, 0.5],
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )

    return type_4_wing_fixture


def make_type_5_wing_fixture():
    """This method makes a fixture that is a Wing with type 5 symmetry
    (symmetric=True, coincident_symmetry_plane=False).

    Note: Type 5 wings are automatically processed by Airplane.process_wing_symmetry()
    into type 1 and type 3 wings, so this fixture represents the initial state
    before processing.

    :return type_5_wing_fixture: Wing
        This is the Wing configured for type 5 symmetry testing.
    """
    # Create WingCrossSections for the wing with symmetric control surfaces
    root_wcs = make_root_wing_cross_section_fixture()
    root_wcs.control_surface_symmetry_type = "symmetric"

    tip_wcs = make_tip_wing_cross_section_fixture()
    # Set chord and control surface for tip (since it was None)
    tip_wcs.chord = 1.0
    tip_wcs.control_surface_symmetry_type = "symmetric"
    tip_wcs.control_surface_hinge_point = 0.75
    tip_wcs.control_surface_deflection = 2.0

    # Create type 5 wing (symmetric=True, non-coincident symmetry plane)
    type_5_wing_fixture = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Type 5 Test Wing",
        prelimLer_G_Cg=[1.0, 0.0, 0.5],
        angles_G_to_prelimWn_ixyz=[0.0, 5.0, 0.0],
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=[0.0, 0.707, 0.707],
        symmetryPoint_G_Cg=[0.5, 0.0, 0.0],
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )

    return type_5_wing_fixture


def make_basic_airplane_fixture():
    """This method makes a fixture that is an Airplane with basic configuration for
    general testing.

    :return basic_airplane_fixture: Airplane
        This is the Airplane configured for general testing.
    """
    # Create a basic Wing for the Airplane
    wing = make_type_1_wing_fixture()

    # Create the basic Airplane
    basic_airplane_fixture = ps.geometry.airplane.Airplane(
        wings=[wing],
        name="Basic Test Airplane",
        Cg_E_CgP1=[1.0, 0.5, -0.2],
        angles_E_to_B_izyx=[10.0, -5.0, 15.0],
        weight=1000.0,
    )

    return basic_airplane_fixture


def make_first_airplane_fixture():
    """This method makes a fixture that is an Airplane suitable for use
    as the first Airplane in a simulation (with Cg_E_CgP1 set to zeros).

    :return first_airplane_fixture: Airplane
        This is the Airplane configured as the first in a simulation.
    """
    # Create a Wing for the Airplane
    wing = make_type_4_wing_fixture()

    # Create the first Airplane
    first_airplane_fixture = ps.geometry.airplane.Airplane(
        wings=[wing],
        name="First Test Airplane",
        Cg_E_CgP1=[0.0, 0.0, 0.0],
        angles_E_to_B_izyx=[0.0, 0.0, 0.0],
        weight=1500.0,
    )

    return first_airplane_fixture


def make_multi_wing_airplane_fixture():
    """This method makes a fixture that is an Airplane with multiple Wings for
    testing multi-wing configurations.

    :return multi_wing_airplane_fixture: Airplane
        This is the Airplane with multiple Wings.
    """
    # Create multiple Wings for the Airplane
    main_wing = make_type_4_wing_fixture()
    main_wing.name = "Main Wing"

    tail_wing = make_type_2_wing_fixture()
    tail_wing.name = "Tail Wing"
    tail_wing.prelimLer_G_Cg = [0.0, 0.0, 5.0]  # Move tail wing back

    # Create the multi-wing Airplane
    multi_wing_airplane_fixture = ps.geometry.airplane.Airplane(
        wings=[main_wing, tail_wing],
        name="Multi-Wing Test Airplane",
        Cg_E_CgP1=[2.0, -1.0, 0.5],
        angles_E_to_B_izyx=[-5.0, 2.0, -10.0],
        weight=2000.0,
        s_ref=20.0,
        c_ref=1.5,
        b_ref=8.0,
    )

    return multi_wing_airplane_fixture


def make_type_5_wing_airplane_fixture():
    """This method makes a fixture that is an Airplane with a type 5 Wing to test
    Wing symmetry processing (type 5 symmetry gets split into two Wings).

    :return type_5_wing_airplane_fixture: Airplane
        This is the Airplane with a type 5 Wing that will be processed.
    """
    # Create a type 5 Wing for the Airplane
    wing = make_type_5_wing_fixture()

    # Create the Airplane (this will process the type 5 Wing)
    type_5_wing_airplane_fixture = ps.geometry.airplane.Airplane(
        wings=[wing],
        name="Type 5 Wing Test Airplane",
        Cg_E_CgP1=[0.0, 0.0, 0.0],
        angles_E_to_B_izyx=[0.0, 0.0, 0.0],
        weight=1200.0,
    )

    return type_5_wing_airplane_fixture


def make_custom_reference_airplane_fixture():
    """This method makes a fixture that is an Airplane with custom reference
    dimensions set explicitly.

    :return custom_reference_airplane_fixture: Airplane
        This is the Airplane with custom reference dimensions.
    """
    # Create a Wing for the Airplane
    wing = make_type_1_wing_fixture()

    # Create the Airplane with custom reference dimensions
    custom_reference_airplane_fixture = ps.geometry.airplane.Airplane(
        wings=[wing],
        name="Custom Reference Test Airplane",
        Cg_E_CgP1=[0.5, 1.0, -0.8],
        angles_E_to_B_izyx=[20.0, -10.0, 5.0],
        weight=800.0,
        s_ref=15.0,
        c_ref=2.0,
        b_ref=10.0,
    )

    return custom_reference_airplane_fixture


def make_naca0012_airfoil_fixture():
    """This method makes a fixture that is a symmetric NACA 0012 Airfoil for
    testing purposes.

    :return naca0012_airfoil_fixture: Airfoil
        This is the symmetric NACA 0012 Airfoil configured for testing.
    """
    naca0012_airfoil_fixture = ps.geometry.airfoil.Airfoil(name="naca0012")

    return naca0012_airfoil_fixture


def make_naca2412_airfoil_fixture():
    """This method makes a fixture that is a cambered NACA 2412 Airfoil for
    testing purposes.

    :return naca2412_airfoil_fixture: Airfoil
        This is the cambered NACA 2412 Airfoil configured for testing.
    """
    naca2412_airfoil_fixture = ps.geometry.airfoil.Airfoil(name="naca2412")

    return naca2412_airfoil_fixture


def make_custom_outline_airfoil_fixture():
    """This method makes a fixture that is an Airfoil with custom outline
    coordinates for testing purposes.

    :return custom_outline_airfoil_fixture: Airfoil
        This is the Airfoil with custom outline coordinates.
    """
    custom_outline = np.array(
        [
            [1.00, 0.00],  # Outline upper trailing point
            [0.75, 0.05],
            [0.50, 0.03],
            [0.25, 0.01],
            [0.00, 0.00],  # Outline leading point
            [0.25, -0.01],
            [0.50, -0.03],
            [0.75, -0.05],
            [1.00, -0.00],  # Outline lower trailing point
        ]
    )

    custom_outline_airfoil_fixture = ps.geometry.airfoil.Airfoil(
        name="Custom Test Airfoil",
        outline_A_lp=custom_outline,
        resample=False,
        n_points_per_side=400,
    )

    return custom_outline_airfoil_fixture


def make_resampled_airfoil_fixture():
    """This method makes a fixture that is an Airfoil with resampling enabled
    for testing purposes.

    :return resampled_airfoil_fixture: Airfoil
        This is the Airfoil with resampling enabled.
    """
    resampled_airfoil_fixture = ps.geometry.airfoil.Airfoil(
        name="naca3210",
        resample=True,
        n_points_per_side=100,
    )

    return resampled_airfoil_fixture


def make_non_resampled_airfoil_fixture():
    """This method makes a fixture that is an Airfoil with resampling disabled
    for testing purposes.

    :return non_resampled_airfoil_fixture: Airfoil
        This is the Airfoil with resampling disabled.
    """
    simple_outline = np.array(
        [
            [1.00, 0.00],  # Outline upper trailing point
            [0.50, 0.03],
            [0.00, 0.00],  # Outline leading point
            [0.50, -0.03],
            [1.00, -0.00],  # Outline lower trailing point
        ]
    )

    non_resampled_airfoil_fixture = ps.geometry.airfoil.Airfoil(
        name="Non-Resampled Test Airfoil",
        outline_A_lp=simple_outline,
        resample=False,
        n_points_per_side=10,
    )

    return non_resampled_airfoil_fixture


def make_basic_panel_fixture():
    """This method makes a fixture that is a Panel for testing purposes.

    :return basic_panel_fixture: Panel
        This is the Panel configured for testing.
    """
    # Define panel corner points (in geometry axes, relative to the CG) for a
    # simple rectangular panel with 1.0 m chord and 0.5 m span.
    basic_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    return basic_panel_fixture
