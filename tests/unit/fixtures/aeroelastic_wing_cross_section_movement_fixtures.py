"""This module contains functions to create AeroelasticWingCrossSectionMovements for
use in tests."""

import pterasoftware as ps

from . import geometry_fixtures


def make_static_aeroelastic_wing_cross_section_movement_fixture():
    """This method makes a fixture that is an AeroelasticWingCrossSectionMovement with
    all parameters zero (no prescribed movement).

    :return static_aeroelastic_wing_cross_section_movement_fixture:
        AeroelasticWingCrossSectionMovement
        This is the AeroelasticWingCrossSectionMovement with no prescribed movement.
        Its base WingCrossSection is a tip WingCrossSection with non zero base angles,
        so that adding structural deformation produces a predictable result.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the static AeroelasticWingCrossSectionMovement.
    static_aeroelastic_wing_cross_section_movement_fixture = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the AeroelasticWingCrossSectionMovement fixture.
    return static_aeroelastic_wing_cross_section_movement_fixture


def make_basic_aeroelastic_wing_cross_section_movement_fixture():
    """This method makes a fixture that is an AeroelasticWingCrossSectionMovement with
    general-purpose moderate values.

    :return basic_aeroelastic_wing_cross_section_movement_fixture:
        AeroelasticWingCrossSectionMovement
        This is the AeroelasticWingCrossSectionMovement with general-purpose values.
    """
    # Initialize the constructing fixture.
    # Use the tip fixture to ensure Lp values stay non negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the basic AeroelasticWingCrossSectionMovement.
    basic_aeroelastic_wing_cross_section_movement_fixture = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.4, 0.3, 0.15),
        periodLp_Wcsp_Lpp=(2.0, 2.0, 2.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(15.0, 10.0, 5.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(2.0, 2.0, 2.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the AeroelasticWingCrossSectionMovement fixture.
    return basic_aeroelastic_wing_cross_section_movement_fixture
