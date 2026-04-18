"""This module contains functions to create CoreAirplaneMovements for use in tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware._core import CoreAirplaneMovement

from . import core_wing_movement_fixtures, geometry_fixtures


def make_static_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with all parameters
    zero (no movement).

    :return static_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with no movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

    # Create the static CoreAirplaneMovement.
    static_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return static_core_airplane_movement_fixture


def make_basic_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with general-purpose
    moderate values.

    :return basic_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with general-purpose values.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_basic_core_wing_movement_fixture()
    ]

    # Create the basic CoreAirplaneMovement.
    basic_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return basic_core_airplane_movement_fixture


def make_sine_spacing_Cg_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with sine spacing
    for Cg_GP1_CgP1.

    :return sine_spacing_Cg_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with sine spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

    # Create the CoreAirplaneMovement with sine spacing for Cg_GP1_CgP1.
    sine_spacing_Cg_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.1, 0.0, 0.0),
        periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return sine_spacing_Cg_core_airplane_movement_fixture


def make_uniform_spacing_Cg_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with uniform spacing
    for Cg_GP1_CgP1.

    :return uniform_spacing_Cg_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with uniform spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

    # Create the CoreAirplaneMovement with uniform spacing for Cg_GP1_CgP1.
    uniform_spacing_Cg_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.1, 0.0, 0.0),
        periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("uniform", "uniform", "uniform"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return uniform_spacing_Cg_core_airplane_movement_fixture


def make_mixed_spacing_Cg_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with mixed spacing
    for Cg_GP1_CgP1.

    :return mixed_spacing_Cg_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with mixed spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

    # Create the CoreAirplaneMovement with mixed spacing for Cg_GP1_CgP1.
    mixed_spacing_Cg_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.1, 0.08, 0.06),
        periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
        spacingCg_GP1_CgP1=("sine", "uniform", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return mixed_spacing_Cg_core_airplane_movement_fixture


def make_Cg_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement where Cg_GP1_CgP1 moves.

    :return Cg_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with Cg_GP1_CgP1 movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

    # Create the moving Cg CoreAirplaneMovement.
    Cg_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.08, 0.06, 0.05),
        periodCg_GP1_CgP1=(1.5, 1.5, 1.5),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return Cg_core_airplane_movement_fixture


def make_phase_offset_Cg_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with non-zero phase
    offset for Cg_GP1_CgP1.

    :return phase_offset_Cg_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with phase offset for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

    # Create the phase-offset CoreAirplaneMovement.
    phase_offset_Cg_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.08, 0.06, 0.05),
        periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(90.0, -45.0, 60.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return phase_offset_Cg_core_airplane_movement_fixture


def make_multiple_periods_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with different
    periods for different dimensions.

    :return multiple_periods_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with different periods.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_multiple_periods_core_wing_movement_fixture()
    ]

    # Create the multiple-periods CoreAirplaneMovement.
    multiple_periods_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.06, 0.05, 0.04),
        periodCg_GP1_CgP1=(1.0, 2.0, 3.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return multiple_periods_core_airplane_movement_fixture


def make_custom_spacing_Cg_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with a custom
    spacing function for Cg_GP1_CgP1.

    :return custom_spacing_Cg_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with custom spacing for Cg_GP1_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
    ]

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

    # Create the custom spacing CoreAirplaneMovement.
    custom_spacing_Cg_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.08, 0.0, 0.0),
        periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=(custom_harmonic, "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return custom_spacing_Cg_core_airplane_movement_fixture


def make_mixed_custom_and_standard_spacing_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with mixed custom
    and standard spacing functions.

    :return mixed_custom_and_standard_spacing_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with mixed custom and standard spacing.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_mixed_custom_and_standard_spacing_core_wing_movement_fixture()
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

    # Create the mixed-spacing CoreAirplaneMovement.
    mixed_custom_and_standard_spacing_core_airplane_movement_fixture = (
        CoreAirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.06, 0.05, 0.04),
            periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
            spacingCg_GP1_CgP1=(custom_harmonic, "uniform", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreAirplaneMovement fixture.
    return mixed_custom_and_standard_spacing_core_airplane_movement_fixture


def make_periodic_geometry_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with periodic geometry
    motion suitable for testing the variable geometry optimization.

    The fixture uses a 0.1s period which aligns well with common delta_time values
    like 0.01s (10 steps per period) and 0.02s (5 steps per period).

    :return periodic_geometry_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with periodic geometry motion.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_periodic_geometry_core_wing_movement_fixture()
    ]

    # Create the periodic geometry CoreAirplaneMovement.
    periodic_geometry_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return periodic_geometry_core_airplane_movement_fixture


def make_angles_only_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement where only Wing angles
    move (no position movement).

    This is useful for testing geometry matching code that compares Wing angles
    separately from Wing positions.

    :return angles_only_core_airplane_movement_fixture: CoreAirplaneMovement
        This is the CoreAirplaneMovement with only Wing angle movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_angles_only_core_wing_movement_fixture()
    ]

    # Create the angles-only CoreAirplaneMovement.
    angles_only_core_airplane_movement_fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Return the CoreAirplaneMovement fixture.
    return angles_only_core_airplane_movement_fixture


def make_2_chordwise_panels_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with a Wing that has
    2 chordwise panels.

    This is useful for testing panel shape comparison in geometry matching.

    :return: CoreAirplaneMovement with a Wing that has 2 chordwise panels.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_2_chordwise_panels_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_2_chordwise_panels_core_wing_movement_fixture()
    ]

    # Create the CoreAirplaneMovement.
    fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    return fixture


def make_3_chordwise_panels_core_airplane_movement_fixture():
    """This method makes a fixture that is an CoreAirplaneMovement with a Wing that has
    3 chordwise panels.

    This is useful for testing panel shape comparison in geometry matching.

    :return: CoreAirplaneMovement with a Wing that has 3 chordwise panels.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_3_chordwise_panels_airplane_fixture()
    wing_movements = [
        core_wing_movement_fixtures.make_3_chordwise_panels_core_wing_movement_fixture()
    ]

    # Create the CoreAirplaneMovement.
    fixture = CoreAirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    return fixture
