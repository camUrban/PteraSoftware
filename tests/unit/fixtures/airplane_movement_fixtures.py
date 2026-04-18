"""This module contains functions to create AirplaneMovements for use in tests."""

import pterasoftware as ps

from . import geometry_fixtures, wing_movement_fixtures


def make_static_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with all parameters
    zero (no movement).

    :return static_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with no movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the static AirplaneMovement.
    static_airplane_movement_fixture = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the AirplaneMovement fixture.
    return static_airplane_movement_fixture


def make_basic_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with general-purpose
    moderate values.

    :return basic_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with general-purpose values.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_basic_wing_movement_fixture()]

    # Create the basic AirplaneMovement.
    basic_airplane_movement_fixture = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the AirplaneMovement fixture.
    return basic_airplane_movement_fixture


def make_periodic_geometry_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with periodic geometry
    motion suitable for testing the variable geometry optimization.

    The fixture uses a 0.1s period which aligns well with common delta_time values
    like 0.01s (10 steps per period) and 0.02s (5 steps per period).

    :return periodic_geometry_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with periodic geometry motion.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        wing_movement_fixtures.make_periodic_geometry_wing_movement_fixture()
    ]

    # Create the periodic geometry AirplaneMovement.
    periodic_geometry_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return periodic_geometry_airplane_movement_fixture
