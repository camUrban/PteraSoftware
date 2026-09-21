"""This module contains functions to create WingCrossSectionMovements for use in
tests."""

import pterasoftware as ps

from . import geometry_fixtures


def make_static_wing_cross_section_movement_fixture(
    base_wing_cross_section: (
        ps.geometry.wing_cross_section.WingCrossSection | None
    ) = None,
) -> ps.movements.wing_cross_section_movement.WingCrossSectionMovement:
    """This method makes a fixture that is a WingCrossSectionMovement with all
    parameters zero (no movement).

    :param base_wing_cross_section: WingCrossSection, optional This is the base
        WingCrossSection to build the movement around. If None, a new root
        WingCrossSection fixture will be created. The default is None.
    :return static_wing_cross_section_movement_fixture: WingCrossSectionMovement This is
        the WingCrossSectionMovement with no movement.
    """
    # Use the provided WingCrossSection or create a new one.
    if base_wing_cross_section is None:
        base_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )

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


def make_static_tip_wing_cross_section_movement_fixture(
    base_wing_cross_section: (
        ps.geometry.wing_cross_section.WingCrossSection | None
    ) = None,
) -> ps.movements.wing_cross_section_movement.WingCrossSectionMovement:
    """This method makes a fixture that is a WingCrossSectionMovement with all
    parameters zero (no movement), using a tip WingCrossSection as the base.

    :param base_wing_cross_section: WingCrossSection, optional This is the base
        WingCrossSection to build the movement around. If None, a new tip
        WingCrossSection fixture will be created. The default is None.
    :return static_tip_wing_cross_section_movement_fixture: WingCrossSectionMovement
        This is the WingCrossSectionMovement with no movement for a tip cross section.
    """
    # Use the provided WingCrossSection or create a new one.
    if base_wing_cross_section is None:
        base_wing_cross_section = (
            geometry_fixtures.make_tip_wing_cross_section_fixture()
        )

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


def make_basic_wing_cross_section_movement_fixture(
    base_wing_cross_section: (
        ps.geometry.wing_cross_section.WingCrossSection | None
    ) = None,
) -> ps.movements.wing_cross_section_movement.WingCrossSectionMovement:
    """This method makes a fixture that is a WingCrossSectionMovement with general-
    purpose moderate values.

    :param base_wing_cross_section: WingCrossSection, optional This is the base
        WingCrossSection to build the movement around. If None, a new tip
        WingCrossSection fixture will be created. The default is None.
    :return basic_wing_cross_section_movement_fixture: WingCrossSectionMovement This is
        the WingCrossSectionMovement with general-purpose values.
    """
    # Use the provided WingCrossSection or create a new one. The tip fixture ensures Lp
    # values stay non-negative during oscillation.
    if base_wing_cross_section is None:
        base_wing_cross_section = (
            geometry_fixtures.make_tip_wing_cross_section_fixture()
        )

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
