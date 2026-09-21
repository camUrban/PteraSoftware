"""This module contains functions to create WingMovements for use in tests."""

import pterasoftware as ps

from . import geometry_fixtures, wing_cross_section_movement_fixtures


def make_static_wing_movement_fixture(
    base_wing: ps.geometry.wing.Wing | None = None,
) -> ps.movements.wing_movement.WingMovement:
    """This method makes a fixture that is a WingMovement with all parameters zero (no
    movement).

    :param base_wing: Wing, optional This is the base Wing to build the movement around.
        If None, a new origin Wing fixture will be created. The default is None.
    :return static_wing_movement_fixture: WingMovement This is the WingMovement with no
        movement.
    """
    # Use the provided Wing or create a new one.
    if base_wing is None:
        base_wing = geometry_fixtures.make_origin_wing_fixture()
    wing_cross_section_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(
            base_wing.wing_cross_sections[0]
        ),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(
            base_wing.wing_cross_sections[1]
        ),
    ]

    # Create the static WingMovement.
    static_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wing_cross_section_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return static_wing_movement_fixture


def make_basic_wing_movement_fixture(
    base_wing: ps.geometry.wing.Wing | None = None,
) -> ps.movements.wing_movement.WingMovement:
    """This method makes a fixture that is a WingMovement with general-purpose moderate
    values.

    :param base_wing: Wing, optional This is the base Wing to build the movement around.
        If None, a new origin Wing fixture will be created. The default is None.
    :return basic_wing_movement_fixture: WingMovement This is the WingMovement with
        general-purpose values.
    """
    # Use the provided Wing or create a new one.
    if base_wing is None:
        base_wing = geometry_fixtures.make_origin_wing_fixture()
    wing_cross_section_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(
            base_wing.wing_cross_sections[0]
        ),
        wing_cross_section_movement_fixtures.make_basic_wing_cross_section_movement_fixture(
            base_wing.wing_cross_sections[1]
        ),
    ]

    # Create the basic WingMovement.
    basic_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wing_cross_section_movements,
        ampLer_Gs_Cgs=(0.1, 0.05, 0.08),
        periodLer_Gs_Cgs=(2.0, 2.0, 2.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(5.0, 3.0, 2.0),
        periodAngles_Gs_to_Wn_ixyz=(2.0, 2.0, 2.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return basic_wing_movement_fixture


def make_periodic_geometry_wing_movement_fixture(
    base_wing: ps.geometry.wing.Wing | None = None,
) -> ps.movements.wing_movement.WingMovement:
    """This method makes a fixture that is a WingMovement with periodic geometry motion,
    suitable for testing the variable geometry optimization. The fixture uses a 0.1s
    period which aligns well with common delta_time values like 0.01s (10 steps per
    period) and 0.02s (5 steps per period).

    :param base_wing: Wing, optional This is the base Wing to build the movement around.
        If None, a new origin Wing fixture will be created. The default is None.
    :return periodic_geometry_wing_movement_fixture: WingMovement This is the
        WingMovement with periodic geometry motion.
    """
    # Use the provided Wing or create a new one.
    if base_wing is None:
        base_wing = geometry_fixtures.make_origin_wing_fixture()
    wing_cross_section_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(
            base_wing.wing_cross_sections[0]
        ),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(
            base_wing.wing_cross_sections[1]
        ),
    ]

    # Create the periodic-geometry WingMovement.
    periodic_geometry_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wing_cross_section_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(5.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.1, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return periodic_geometry_wing_movement_fixture
