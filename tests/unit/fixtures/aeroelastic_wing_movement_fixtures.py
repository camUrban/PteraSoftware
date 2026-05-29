"""This module contains functions to create AeroelasticWingMovements for use in
tests."""

import pterasoftware as ps

from . import geometry_fixtures


def make_static_aeroelastic_wing_movement_fixture():
    """This method makes a fixture that is an AeroelasticWingMovement with all
    parameters zero (no prescribed movement).

    :return static_aeroelastic_wing_movement_fixture: AeroelasticWingMovement
        This is the AeroelasticWingMovement with no prescribed movement. Its
        AeroelasticWingCrossSectionMovements are built from the base Wing's own
        WingCrossSections so that the generated Wings stay valid.
    """
    # Initialize the constructing fixtures. Build a static
    # AeroelasticWingCrossSectionMovement for each of the base Wing's
    # WingCrossSections, which keeps the generated Wings valid (in particular, the
    # first WingCrossSection keeps its required zero Lp_Wcsp_Lpp).
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section
        )
        for wing_cross_section in base_wing.wing_cross_sections
    ]

    # Create the static AeroelasticWingMovement.
    static_aeroelastic_wing_movement_fixture = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
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
    )

    # Return the AeroelasticWingMovement fixture.
    return static_aeroelastic_wing_movement_fixture


def make_basic_aeroelastic_wing_movement_fixture():
    """This method makes a fixture that is an AeroelasticWingMovement with
    general-purpose moderate values.

    :return basic_aeroelastic_wing_movement_fixture: AeroelasticWingMovement
        This is the AeroelasticWingMovement with general-purpose values. The root
        WingCrossSection's movement is static, while the tip WingCrossSection's
        movement oscillates its angles, so that structural deformation can be tested on
        top of a prescribed oscillation.
    """
    # Initialize the constructing fixtures. The root WingCrossSection's movement stays
    # static (keeping its required zero Lp_Wcsp_Lpp), while the tip's movement
    # oscillates its angles.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    root_wing_cross_section, tip_wing_cross_section = base_wing.wing_cross_sections
    wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=root_wing_cross_section,
        ),
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=tip_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_ixyz=(15.0, 10.0, 5.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(2.0, 2.0, 2.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        ),
    ]

    # Create the basic AeroelasticWingMovement.
    basic_aeroelastic_wing_movement_fixture = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
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
    )

    # Return the AeroelasticWingMovement fixture.
    return basic_aeroelastic_wing_movement_fixture


def make_rotation_offset_aeroelastic_wing_movement_fixture():
    """This method makes a fixture that is an AeroelasticWingMovement identical to the
    basic fixture but with a non zero rotationPointOffset_Gs_Ler.

    :return rotation_offset_aeroelastic_wing_movement_fixture: AeroelasticWingMovement
        This is the AeroelasticWingMovement with a non zero rotationPointOffset_Gs_Ler,
        for testing that angular motion about an offset rotation point shifts the Wing's
        Ler_Gs_Cgs. All other parameters match the basic fixture so the offset effect
        can be isolated.
    """
    # Initialize the constructing fixtures, matching the basic fixture's
    # WingCrossSectionMovements.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    root_wing_cross_section, tip_wing_cross_section = base_wing.wing_cross_sections
    wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=root_wing_cross_section,
        ),
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=tip_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_ixyz=(15.0, 10.0, 5.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(2.0, 2.0, 2.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        ),
    ]

    # Create the AeroelasticWingMovement with a non zero rotation point offset.
    rotation_offset_aeroelastic_wing_movement_fixture = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
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
            rotationPointOffset_Gs_Ler=(0.1, 0.2, -0.1),
        )
    )

    # Return the AeroelasticWingMovement fixture.
    return rotation_offset_aeroelastic_wing_movement_fixture
