"""This module contains functions to create AeroelasticAirplaneMovements for use in
tests."""

import pterasoftware as ps

from . import geometry_fixtures


def make_static_aeroelastic_airplane_movement_fixture():
    """This method makes a fixture that is an AeroelasticAirplaneMovement with all
    parameters zero (no prescribed movement).

    :return static_aeroelastic_airplane_movement_fixture: AeroelasticAirplaneMovement
        This is the AeroelasticAirplaneMovement with no prescribed movement. Its single
        AeroelasticWingMovement child's AeroelasticWingCrossSectionMovements are built
        from the base Airplane's own WingCrossSections so that the generated Airplanes
        stay valid.
    """
    # Initialize the constructing fixtures. Build the base Airplane first, then build
    # the movements from the Airplane's own Wing and WingCrossSections so that they
    # match after the Airplane's symmetry processing.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    base_airplane = ps.geometry.airplane.Airplane(
        wings=[base_wing],
        name="Origin Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )

    # Build a static AeroelasticWingCrossSectionMovement for each of the base Wing's
    # WingCrossSections, which keeps the generated Wings valid (in particular, the
    # first WingCrossSection keeps its required zero Lp_Wcsp_Lpp).
    wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section
        )
        for wing_cross_section in base_airplane.wings[0].wing_cross_sections
    ]

    # Create the static AeroelasticWingMovement child.
    wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
        base_wing=base_airplane.wings[0],
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

    # Create the static AeroelasticAirplaneMovement.
    static_aeroelastic_airplane_movement_fixture = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AeroelasticAirplaneMovement fixture.
    return static_aeroelastic_airplane_movement_fixture


def make_basic_aeroelastic_airplane_movement_fixture():
    """This method makes a fixture that is an AeroelasticAirplaneMovement with
    general-purpose moderate values.

    :return basic_aeroelastic_airplane_movement_fixture: AeroelasticAirplaneMovement
        This is the AeroelasticAirplaneMovement with general-purpose values. Its single
        AeroelasticWingMovement child oscillates its Wing's position and angles, the
        root WingCrossSection's movement stays static, and the tip WingCrossSection's
        movement oscillates its angles, so that structural deformation can be tested on
        top of a prescribed oscillation. Because the base Airplane is the first Airplane
        in a simulation, its Cg_GP1_CgP1 motion stays zero.
    """
    # Initialize the constructing fixtures. Build the base Airplane first, then build
    # the movements from the Airplane's own Wing and WingCrossSections so that they
    # match after the Airplane's symmetry processing.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    base_airplane = ps.geometry.airplane.Airplane(
        wings=[base_wing],
        name="Origin Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )

    # The root WingCrossSection's movement stays static (keeping its required zero
    # Lp_Wcsp_Lpp), while the tip's movement oscillates its angles.
    root_wing_cross_section, tip_wing_cross_section = base_airplane.wings[
        0
    ].wing_cross_sections
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

    # Create the basic AeroelasticWingMovement child.
    wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
        base_wing=base_airplane.wings[0],
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

    # Create the basic AeroelasticAirplaneMovement.
    basic_aeroelastic_airplane_movement_fixture = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AeroelasticAirplaneMovement fixture.
    return basic_aeroelastic_airplane_movement_fixture


def make_mixed_wing_aeroelastic_airplane_movement_fixture():
    """This method makes a fixture that is an AeroelasticAirplaneMovement whose first
    Wing is backed by an AeroelasticWingMovement and whose second Wing is backed by a
    standard WingMovement.

    This exercises both branches of generate_airplane_at_time_step at once: the
    AeroelasticWingMovement child receives its per Wing deformation, while the standard
    WingMovement child is advanced without deformation.

    :return mixed_wing_aeroelastic_airplane_movement_fixture: AeroelasticAirplaneMovement
        This is the AeroelasticAirplaneMovement with one AeroelasticWingMovement child
        and one standard WingMovement child. Both Wings are static so that any change in
        a generated Wing's angles can be attributed to applied deformation.
    """
    # Initialize the constructing fixtures. Build the base Airplane first, then build
    # the movements from the Airplane's own Wings and WingCrossSections so that they
    # match after the Airplane's symmetry processing.
    aeroelastic_wing = geometry_fixtures.make_origin_wing_fixture()
    standard_wing = geometry_fixtures.make_origin_wing_fixture()
    base_airplane = ps.geometry.airplane.Airplane(
        wings=[aeroelastic_wing, standard_wing],
        name="Origin Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )

    # Build the AeroelasticWingMovement child for the first Wing.
    aeroelastic_wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section
        )
        for wing_cross_section in base_airplane.wings[0].wing_cross_sections
    ]
    aeroelastic_wing_movement = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=base_airplane.wings[0],
            wing_cross_section_movements=aeroelastic_wing_cross_section_movements,
        )
    )

    # Build the standard WingMovement child for the second Wing.
    standard_wing_cross_section_movements = [
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section
        )
        for wing_cross_section in base_airplane.wings[1].wing_cross_sections
    ]
    standard_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=base_airplane.wings[1],
        wing_cross_section_movements=standard_wing_cross_section_movements,
    )

    # Create the mixed AeroelasticAirplaneMovement.
    mixed_wing_aeroelastic_airplane_movement_fixture = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[aeroelastic_wing_movement, standard_wing_movement],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AeroelasticAirplaneMovement fixture.
    return mixed_wing_aeroelastic_airplane_movement_fixture
