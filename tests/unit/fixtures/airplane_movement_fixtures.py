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
        ampCg_E_CgP1=(0.0, 0.0, 0.0),
        periodCg_E_CgP1=(0.0, 0.0, 0.0),
        spacingCg_E_CgP1=("sine", "sine", "sine"),
        phaseCg_E_CgP1=(0.0, 0.0, 0.0),
        ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
        phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
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
        ampCg_E_CgP1=(0.05, 0.03, 0.04),
        periodCg_E_CgP1=(2.0, 2.0, 2.0),
        spacingCg_E_CgP1=("sine", "sine", "sine"),
        phaseCg_E_CgP1=(0.0, 0.0, 0.0),
        ampAngles_E_to_B_izyx=(3.0, 2.0, 1.5),
        periodAngles_E_to_B_izyx=(2.0, 2.0, 2.0),
        spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
        phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
    )

    # Return the AirplaneMovement fixture.
    return basic_airplane_movement_fixture


def make_sine_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with sine spacing
    for Cg_E_CgP1.

    :return sine_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with sine spacing for Cg_E_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with sine spacing for Cg_E_CgP1.
    sine_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.1, 0.0, 0.0),
            periodCg_E_CgP1=(1.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return sine_spacing_Cg_airplane_movement_fixture


def make_uniform_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with uniform spacing
    for Cg_E_CgP1.

    :return uniform_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with uniform spacing for Cg_E_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with uniform spacing for Cg_E_CgP1.
    uniform_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.1, 0.0, 0.0),
            periodCg_E_CgP1=(1.0, 0.0, 0.0),
            spacingCg_E_CgP1=("uniform", "uniform", "uniform"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return uniform_spacing_Cg_airplane_movement_fixture


def make_mixed_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with mixed spacing
    for Cg_E_CgP1.

    :return mixed_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with mixed spacing for Cg_E_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with mixed spacing for Cg_E_CgP1.
    mixed_spacing_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.1, 0.08, 0.06),
            periodCg_E_CgP1=(1.0, 1.0, 1.0),
            spacingCg_E_CgP1=("sine", "uniform", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return mixed_spacing_Cg_airplane_movement_fixture


def make_sine_spacing_angles_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with sine spacing
    for angles_E_to_B_izyx.

    :return sine_spacing_angles_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with sine spacing for angles_E_to_B_izyx.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with sine spacing for angles_E_to_B_izyx.
    sine_spacing_angles_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.0, 0.0, 0.0),
            periodCg_E_CgP1=(0.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(10.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(1.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return sine_spacing_angles_airplane_movement_fixture


def make_uniform_spacing_angles_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with uniform spacing
    for angles_E_to_B_izyx.

    :return uniform_spacing_angles_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with uniform spacing for angles_E_to_B_izyx.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with uniform spacing for angles_E_to_B_izyx.
    uniform_spacing_angles_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.0, 0.0, 0.0),
            periodCg_E_CgP1=(0.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(10.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(1.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("uniform", "uniform", "uniform"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return uniform_spacing_angles_airplane_movement_fixture


def make_mixed_spacing_angles_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with mixed spacing
    for angles_E_to_B_izyx.

    :return mixed_spacing_angles_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with mixed spacing for angles_E_to_B_izyx.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the AirplaneMovement with mixed spacing for angles_E_to_B_izyx.
    mixed_spacing_angles_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.0, 0.0, 0.0),
            periodCg_E_CgP1=(0.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(10.0, 15.0, 8.0),
            periodAngles_E_to_B_izyx=(1.0, 1.0, 1.0),
            spacingAngles_E_to_B_izyx=("sine", "uniform", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return mixed_spacing_angles_airplane_movement_fixture


def make_Cg_only_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement where only Cg_E_CgP1
    moves.

    :return Cg_only_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with only Cg_E_CgP1 movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the Cg-only AirplaneMovement.
    Cg_only_airplane_movement_fixture = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=wing_movements,
        ampCg_E_CgP1=(0.08, 0.06, 0.05),
        periodCg_E_CgP1=(1.5, 1.5, 1.5),
        spacingCg_E_CgP1=("sine", "sine", "sine"),
        phaseCg_E_CgP1=(0.0, 0.0, 0.0),
        ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
        phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
    )

    # Return the AirplaneMovement fixture.
    return Cg_only_airplane_movement_fixture


def make_angles_only_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement where only
    angles_E_to_B_izyx moves.

    :return angles_only_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with only angles_E_to_B_izyx movement.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the angles-only AirplaneMovement.
    angles_only_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.0, 0.0, 0.0),
            periodCg_E_CgP1=(0.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(12.0, 8.0, 6.0),
            periodAngles_E_to_B_izyx=(1.5, 1.5, 1.5),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return angles_only_airplane_movement_fixture


def make_phase_offset_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with non-zero phase
    offset for Cg_E_CgP1.

    :return phase_offset_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with phase offset for Cg_E_CgP1.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the phase-offset AirplaneMovement.
    phase_offset_Cg_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.08, 0.06, 0.05),
            periodCg_E_CgP1=(1.0, 1.0, 1.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(90.0, -45.0, 60.0),
            ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return phase_offset_Cg_airplane_movement_fixture


def make_phase_offset_angles_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with non-zero phase
    offset for angles_E_to_B_izyx.

    :return phase_offset_angles_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with phase offset for angles_E_to_B_izyx.
    """
    # Initialize the constructing fixtures.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

    # Create the phase-offset AirplaneMovement.
    phase_offset_angles_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.0, 0.0, 0.0),
            periodCg_E_CgP1=(0.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(10.0, 12.0, 8.0),
            periodAngles_E_to_B_izyx=(1.0, 1.0, 1.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(45.0, 90.0, -30.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return phase_offset_angles_airplane_movement_fixture


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
            ampCg_E_CgP1=(0.06, 0.05, 0.04),
            periodCg_E_CgP1=(1.0, 2.0, 3.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(8.0, 10.0, 12.0),
            periodAngles_E_to_B_izyx=(0.5, 1.5, 2.5),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return multiple_periods_airplane_movement_fixture


def make_custom_spacing_Cg_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with a custom
    spacing function for Cg_E_CgP1.

    :return custom_spacing_Cg_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with custom spacing for Cg_E_CgP1.
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
            ampCg_E_CgP1=(0.08, 0.0, 0.0),
            periodCg_E_CgP1=(1.0, 0.0, 0.0),
            spacingCg_E_CgP1=(custom_harmonic, "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=("sine", "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return custom_spacing_Cg_airplane_movement_fixture


def make_custom_spacing_angles_airplane_movement_fixture():
    """This method makes a fixture that is an AirplaneMovement with a custom
    spacing function for angles_E_to_B_izyx.

    :return custom_spacing_angles_airplane_movement_fixture: AirplaneMovement
        This is the AirplaneMovement with custom spacing for angles_E_to_B_izyx.
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
    custom_spacing_angles_airplane_movement_fixture = (
        ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_E_CgP1=(0.0, 0.0, 0.0),
            periodCg_E_CgP1=(0.0, 0.0, 0.0),
            spacingCg_E_CgP1=("sine", "sine", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(10.0, 0.0, 0.0),
            periodAngles_E_to_B_izyx=(1.0, 0.0, 0.0),
            spacingAngles_E_to_B_izyx=(custom_harmonic, "sine", "sine"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return custom_spacing_angles_airplane_movement_fixture


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
            ampCg_E_CgP1=(0.06, 0.05, 0.04),
            periodCg_E_CgP1=(1.0, 1.0, 1.0),
            spacingCg_E_CgP1=(custom_harmonic, "uniform", "sine"),
            phaseCg_E_CgP1=(0.0, 0.0, 0.0),
            ampAngles_E_to_B_izyx=(8.0, 10.0, 6.0),
            periodAngles_E_to_B_izyx=(1.0, 1.0, 1.0),
            spacingAngles_E_to_B_izyx=("sine", custom_harmonic, "uniform"),
            phaseAngles_E_to_B_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the AirplaneMovement fixture.
    return mixed_custom_and_standard_spacing_airplane_movement_fixture
