"""This module contains functions to create AirplaneMovements for use in tests."""

import numpy as np
import pterasoftware as ps

from . import geometry_fixtures
from . import wing_movement_fixtures


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


def make_sine_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with sine spacing
    for Cg_GP1_CgP1.

    :return sine_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with sine spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with sine spacing for Cg_GP1_CgP1.
    sine_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.1, 0.0, 0.0),
            periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return sine_spacing_Cg_airplane_movement_fixture


def make_uniform_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with uniform spacing
    for Cg_GP1_CgP1.

    :return uniform_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with uniform spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with uniform spacing for Cg_GP1_CgP1.
    uniform_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.1, 0.0, 0.0),
            periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("uniform", "uniform", "uniform"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return uniform_spacing_Cg_airplane_movement_fixture


def make_mixed_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with mixed spacing
    for Cg_GP1_CgP1.

    :return mixed_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with mixed spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with mixed spacing for Cg_GP1_CgP1.
    mixed_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.1, 0.08, 0.06),
            periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
            spacingCg_GP1_CgP1=("sine", "uniform", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return mixed_spacing_Cg_airplane_movement_fixture


def make_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement where Cg_GP1_CgP1 moves.

    :return Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with Cg_GP1_CgP1 movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the moving Cg AirplaneMovement.
    Cg_airplane_movement_fixture = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.08, 0.06, 0.05),
        periodCg_GP1_CgP1=(1.5, 1.5, 1.5),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the AirplaneMovement fixture.
    return Cg_airplane_movement_fixture


def make_phase_offset_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with non-zero phase
    offset for Cg_GP1_CgP1.

    :return phase_offset_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with phase offset for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the phase-offset AirplaneMovement.
    phase_offset_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.08, 0.06, 0.05),
            periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(90.0, -45.0, 60.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return phase_offset_Cg_airplane_movement_fixture


def make_multiple_periods_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with different
    periods for different dimensions.

    :return multiple_periods_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with different periods.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        wing_movement_fixtures.make_multiple_periods_wing_movement_fixture()
    ]

    # Create the multiple-periods AirplaneMovement.
    multiple_periods_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.06, 0.05, 0.04),
            periodCg_GP1_CgP1=(1.0, 2.0, 3.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return multiple_periods_airplane_movement_fixture


def make_custom_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with a custom
    spacing function for Cg_GP1_CgP1.

    :return custom_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with custom spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Define a custom harmonic spacing function.
    def custom_harmonic(x):
        """Custom harmonic spacing function: normalized combination of harmonics.

        This function satisfies all requirements: starts at 0, returns to 0 at
        2*pi, has zero mean, has amplitude of 1, and is periodic.

        :param x: (N,) ndarray of floats
            The input angles in radians.

        :return: (N,) ndarray of floats
            The output values.
        """
        return (3.0 / (2.0 * np.sqrt(2.0))) * (
            np.sin(x) + (1.0 / 3.0) * np.sin(3.0 * x)
        )

    # Create the custom-spacing AirplaneMovement.
    custom_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.08, 0.0, 0.0),
            periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=(custom_harmonic, "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return custom_spacing_Cg_airplane_movement_fixture


def make_mixed_custom_and_standard_spacing_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with mixed custom
    and standard spacing functions.

    :return mixed_custom_and_standard_spacing_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with mixed custom and standard spacing.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        wing_movement_fixtures.make_mixed_custom_and_standard_spacing_wing_movement_fixture()
    ]

    # Define a custom harmonic spacing function.
    def custom_harmonic(x):
        """Custom harmonic spacing function: normalized combination of harmonics.

        :param x: (N,) ndarray of floats
            The input angles in radians.

        :return: (N,) ndarray of floats
            The output values.
        """
        return (3.0 / (2.0 * np.sqrt(2.0))) * (
            np.sin(x) + (1.0 / 3.0) * np.sin(3.0 * x)
        )

    # Create the mixed-spacing AirplaneMovement.
    mixed_custom_and_standard_spacing_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.06, 0.05, 0.04),
            periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
            spacingCg_GP1_CgP1=(custom_harmonic, "uniform", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return mixed_custom_and_standard_spacing_airplane_movement_fixture
