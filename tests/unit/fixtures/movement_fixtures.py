"""This module contains functions to create Movements for use in tests."""

import pterasoftware as ps

from . import airplane_movement_fixtures, operating_point_fixtures


def make_basic_aeroelastic_movement_fixture():
    """This method makes a fixture that is an AeroelasticMovement for testing.

    :return basic_aeroelastic_movement_fixture: AeroelasticMovement
        This is the AeroelasticMovement configured for general testing.
    """
    # Create a shared airfoil for both wing cross sections.
    airfoil = ps.geometry.airfoil.Airfoil(name="naca2412")

    # Create the root WCS. The first WCS of a Wing must have Lp_Wcsp_Lpp=(0,0,0).
    root_wcs = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=airfoil,
        num_spanwise_panels=1,
        chord=1.0,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spanwise_spacing="uniform",
    )

    # Create the tip WCS.
    tip_wcs = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=airfoil,
        num_spanwise_panels=None,
        chord=0.5,
        Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spanwise_spacing=None,
    )

    # Create the wing.
    wing = ps.geometry.wing.Wing(
        wing_cross_sections=[root_wcs, tip_wcs],
        name="Test Wing",
        Ler_Gs_Cgs=(0.0, 0.0, 0.0),
        angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        symmetric=False,
        mirror_only=False,
        single_step_wing=False,
        num_chordwise_panels=2,
        chordwise_spacing="uniform",
    )

    # Create the airplane.
    airplane = ps.geometry.airplane.Airplane(
        wings=[wing],
        name="Test Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )

    # Create WCS movements using the airplane's WCS objects.
    root_wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
        base_wing_cross_section=airplane.wings[0].wing_cross_sections[0],
    )
    tip_wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
        base_wing_cross_section=airplane.wings[0].wing_cross_sections[1],
    )

    # Create a wing movement with sinusoidal flapping (non-zero period required
    # by generate_inertial_torque_function when spacing is "sine").
    wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=[root_wcs_movement, tip_wcs_movement],
        ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Create the airplane movement.
    airplane_movement = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=airplane,
            wing_movements=[wing_movement],
        )
    )

    # Create the operating point movement.
    op_point_movement = ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
    )

    # Create the AeroelasticMovement.
    basic_aeroelastic_movement_fixture = (
        ps.movements.aeroelastic_movement.AeroelasticMovement(
            airplane_movements=[airplane_movement],
            operating_point_movement=op_point_movement,
            delta_time=0.1,
            num_steps=3,
        )
    )

    return basic_aeroelastic_movement_fixture


def make_static_movement_fixture():
    """This method makes a fixture that is a Movement with all static components.

    :return static_movement_fixture: Movement
        This is the Movement with no motion.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_static_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the static Movement.
    static_movement_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_chords=3,
    )

    # Return the Movement fixture.
    return static_movement_fixture


def make_basic_movement_fixture():
    """This method makes a fixture that is a Movement with general-purpose values.

    :return basic_movement_fixture: Movement
        This is the Movement with general-purpose values for testing.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the basic Movement.
    basic_movement_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_cycles=1,
    )

    # Return the Movement fixture.
    return basic_movement_fixture


def make_static_movement_with_explicit_num_steps_fixture():
    """This method makes a fixture that is a Movement with static motion and
    explicitly set num_steps.

    :return static_movement_with_explicit_num_steps_fixture: Movement
        This is the Movement with static motion and explicit num_steps.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_static_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with explicit num_steps.
    static_movement_with_explicit_num_steps_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_steps=5,
    )

    # Return the Movement fixture.
    return static_movement_with_explicit_num_steps_fixture


def make_non_static_movement_with_explicit_num_steps_fixture():
    """This method makes a fixture that is a Movement with non static motion and
    explicitly set num_steps.

    :return non_static_movement_with_explicit_num_steps_fixture: Movement
        This is the Movement with non static motion and explicit num_steps.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with explicit num_steps.
    non_static_movement_with_explicit_num_steps_fixture = (
        ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=10,
        )
    )

    # Return the Movement fixture.
    return non_static_movement_with_explicit_num_steps_fixture


def make_movement_with_custom_delta_time_fixture():
    """This method makes a fixture that is a Movement with custom delta_time.

    :return movement_with_custom_delta_time_fixture: Movement
        This is the Movement with custom delta_time parameter.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with custom delta_time.
    movement_with_custom_delta_time_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        delta_time=0.05,
        num_cycles=1,
    )

    # Return the Movement fixture.
    return movement_with_custom_delta_time_fixture


def make_movement_with_multiple_airplanes_fixture():
    """This method makes a fixture that is a Movement with multiple AirplaneMovements.

    :return movement_with_multiple_airplanes_fixture: Movement
        This is the Movement with multiple AirplaneMovements.
    """
    # Initialize the constructing fixtures.
    airplane_movements = [
        airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
    ]
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the Movement with multiple AirplaneMovements.
    movement_with_multiple_airplanes_fixture = ps.movements.movement.Movement(
        airplane_movements=airplane_movements,
        operating_point_movement=operating_point_movement,
        num_cycles=1,
    )

    # Return the Movement fixture.
    return movement_with_multiple_airplanes_fixture
