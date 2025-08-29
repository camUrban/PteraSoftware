"""This module contains functions to create geometry objects for use in tests.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_test_airfoil_fixture: This method makes a fixture that is an Airfoil
    object for testing purposes.

    make_basic_wing_cross_section_fixture: This method makes a fixture that is a
    WingCrossSection object with typical parameters for general testing.

    make_root_wing_cross_section_fixture: This method makes a fixture that is a
    root WingCrossSection object (with zero vectors as required by constraints).

    make_tip_wing_cross_section_fixture: This method makes a fixture that is a
    tip WingCrossSection object (with None values for spanwise parameters).

    make_minimal_wing_cross_section_fixture: This method makes a fixture that is
    a WingCrossSection object with minimal valid parameters.

    make_asymmetric_control_surface_wing_cross_section_fixture: This method makes
    a fixture that is a WingCrossSection object with asymmetric control surface
    configuration.
"""

# ToDo: Uncomment if we end up using numpy again in this module.
# import numpy as np

import pterasoftware as ps


def make_test_airfoil_fixture():
    """This method makes a fixture that is an Airfoil object for testing purposes.

    :return test_airfoil_fixture: Airfoil
        This is the Airfoil object configured for testing.
    """
    test_airfoil_fixture = ps.geometry.Airfoil(name="naca2412")

    return test_airfoil_fixture


def make_basic_wing_cross_section_fixture():
    """This method makes a fixture that is a WingCrossSection object with typical
    parameters for general testing.

    :return basic_wing_cross_section_fixture: WingCrossSection
        This is the WingCrossSection object configured for general testing.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the basic WingCrossSection object.
    basic_wing_cross_section_fixture = ps.geometry.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=8,
        chord=1.5,
        Lp_Wcsp_Lpp=[0.2, 0.5, 0.1],
        angles_Wcsp_to_Wcs_i321=[5.0, -2.0, 3.0],
        control_surface_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=5.0,
        spanwise_spacing="cosine",
    )

    # Delete the constructing fixture.
    del test_airfoil_fixture

    # Return the WingCrossSection fixture.
    return basic_wing_cross_section_fixture


def make_root_wing_cross_section_fixture():
    """This method makes a fixture that is a root WingCrossSection object (with
    zero vectors as required by constraints).

    :return root_wing_cross_section_fixture: WingCrossSection
        This is the root WingCrossSection object.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the root WingCrossSection object.
    root_wing_cross_section_fixture = ps.geometry.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=10,
        chord=2.0,
        Lp_Wcsp_Lpp=[0.0, 0.0, 0.0],
        angles_Wcsp_to_Wcs_i321=[0.0, 0.0, 0.0],
    )

    # Delete the constructing fixture.
    del test_airfoil_fixture

    # Return the root WingCrossSection fixture.
    return root_wing_cross_section_fixture


def make_tip_wing_cross_section_fixture():
    """This method makes a fixture that is a tip WingCrossSection object (with
    None values for spanwise parameters).

    :return tip_wing_cross_section_fixture: WingCrossSection
        This is the tip WingCrossSection object.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the tip WingCrossSection object.
    tip_wing_cross_section_fixture = ps.geometry.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=None,
        chord=0.8,
        Lp_Wcsp_Lpp=[0.5, 2.0, 0.2],
        angles_Wcsp_to_Wcs_i321=[10.0, -5.0, 8.0],
        spanwise_spacing=None,
    )

    # Delete the constructing fixture.
    del test_airfoil_fixture

    # Return the tip WingCrossSection fixture.
    return tip_wing_cross_section_fixture


def make_minimal_wing_cross_section_fixture():
    """This method makes a fixture that is a WingCrossSection object with minimal
    valid parameters.

    :return minimal_wing_cross_section_fixture: WingCrossSection
        This is the minimal WingCrossSection object.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the minimal WingCrossSection object.
    minimal_wing_cross_section_fixture = ps.geometry.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=1,  # Minimum valid value
    )

    # Delete the constructing fixture.
    del test_airfoil_fixture

    # Return the minimal WingCrossSection fixture.
    return minimal_wing_cross_section_fixture


def make_asymmetric_control_surface_wing_cross_section_fixture():
    """This method makes a fixture that is a WingCrossSection object with asymmetric
    control surface configuration.

    :return asymmetric_wing_cross_section_fixture: WingCrossSection
        This is the WingCrossSection object with asymmetric control surface.
    """
    # Initialize the constructing fixture.
    test_airfoil_fixture = make_test_airfoil_fixture()

    # Create the asymmetric control surface WingCrossSection object.
    asymmetric_wing_cross_section_fixture = ps.geometry.WingCrossSection(
        airfoil=test_airfoil_fixture,
        num_spanwise_panels=6,
        chord=1.2,
        Lp_Wcsp_Lpp=[0.1, 1.0, 0.05],
        angles_Wcsp_to_Wcs_i321=[2.0, 0.0, -1.0],
        control_surface_type="asymmetric",
        control_surface_hinge_point=0.8,
        control_surface_deflection=-10.0,
        spanwise_spacing="uniform",
    )

    # Delete the constructing fixture.
    del test_airfoil_fixture

    # Return the asymmetric WingCrossSection fixture.
    return asymmetric_wing_cross_section_fixture
