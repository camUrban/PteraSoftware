"""This module contains functions to create WingCrossSectionMovements for use in
tests."""

import pterasoftware as ps

from . import geometry_fixtures


def make_static_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with all
    parameters zero (no movement).

    :return static_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with no movement.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the static WingCrossSectionMovement.
    static_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
    )

    # Return the WingCrossSectionMovement fixture.
    return static_wing_cross_section_movement_fixture


def make_static_tip_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with all
    parameters zero (no movement), using a tip WingCrossSection as the base.

    :return static_tip_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with no movement for a tip cross
        section.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the static tip WingCrossSectionMovement.
    static_tip_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
    )

    # Return the WingCrossSectionMovement fixture.
    return static_tip_wing_cross_section_movement_fixture


def make_basic_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a WingCrossSectionMovement with
    general-purpose moderate values.

    :return basic_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with general-purpose values.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the basic WingCrossSectionMovement.
    basic_wing_cross_section_movement_fixture = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
    )

    # Return the WingCrossSectionMovement fixture.
    return basic_wing_cross_section_movement_fixture
