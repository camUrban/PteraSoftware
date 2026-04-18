"""This module contains functions to create CoreWingMovements for use in tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware._core import CoreWingMovement

from . import core_wing_cross_section_movement_fixtures, geometry_fixtures


def make_static_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with all parameters
    zero (no movement).

    :return static_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with no movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the static CoreWingMovement.
    static_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return static_core_wing_movement_fixture


def make_basic_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with general-purpose
    moderate values.

    :return basic_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with general-purpose values.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_basic_core_wing_cross_section_movement_fixture(),
    ]

    # Create the basic CoreWingMovement.
    basic_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.1, 0.05, 0.08),
        periodLer_Gs_Cgs=(2.0, 2.0, 2.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(5.0, 3.0, 2.0),
        periodAngles_Gs_to_Wn_ixyz=(2.0, 2.0, 2.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return basic_core_wing_movement_fixture


def make_sine_spacing_Ler_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with sine spacing for
    Ler_Gs_Cgs.

    :return sine_spacing_Ler_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with sine spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with sine spacing for Ler_Gs_Cgs.
    sine_spacing_Ler_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.2, 0.0, 0.0),
        periodLer_Gs_Cgs=(1.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return sine_spacing_Ler_core_wing_movement_fixture


def make_uniform_spacing_Ler_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with uniform spacing for
    Ler_Gs_Cgs.

    :return uniform_spacing_Ler_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with uniform spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with uniform spacing for Ler_Gs_Cgs.
    uniform_spacing_Ler_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.2, 0.0, 0.0),
        periodLer_Gs_Cgs=(1.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("uniform", "uniform", "uniform"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return uniform_spacing_Ler_core_wing_movement_fixture


def make_mixed_spacing_Ler_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with mixed spacing for
    Ler_Gs_Cgs.

    :return mixed_spacing_Ler_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with mixed spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with mixed spacing for Ler_Gs_Cgs.
    mixed_spacing_Ler_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.2, 0.15, 0.1),
        periodLer_Gs_Cgs=(1.0, 1.0, 1.0),
        spacingLer_Gs_Cgs=("sine", "uniform", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return mixed_spacing_Ler_core_wing_movement_fixture


def make_sine_spacing_angles_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with sine spacing for
    angles_Gs_to_Wn_ixyz.

    :return sine_spacing_angles_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with sine spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with sine spacing for angles_Gs_to_Wn_ixyz.
    sine_spacing_angles_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return sine_spacing_angles_core_wing_movement_fixture


def make_uniform_spacing_angles_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with uniform spacing for
    angles_Gs_to_Wn_ixyz.

    :return uniform_spacing_angles_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with uniform spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with uniform spacing for angles_Gs_to_Wn_ixyz.
    uniform_spacing_angles_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("uniform", "uniform", "uniform"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return uniform_spacing_angles_core_wing_movement_fixture


def make_mixed_spacing_angles_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with mixed spacing for
    angles_Gs_to_Wn_ixyz.

    :return mixed_spacing_angles_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with mixed spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with mixed spacing for angles_Gs_to_Wn_ixyz.
    mixed_spacing_angles_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(10.0, 15.0, 8.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 1.0, 1.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "uniform", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return mixed_spacing_angles_core_wing_movement_fixture


def make_Ler_only_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement where only Ler_Gs_Cgs
    moves.

    :return Ler_only_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with only Ler_Gs_Cgs movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the Ler-only CoreWingMovement.
    Ler_only_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.15, 0.1, 0.08),
        periodLer_Gs_Cgs=(1.5, 1.5, 1.5),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return Ler_only_core_wing_movement_fixture


def make_angles_only_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement where only
    angles_Gs_to_Wn_ixyz moves.

    :return angles_only_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with only angles_Gs_to_Wn_ixyz movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the angles-only CoreWingMovement.
    angles_only_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(12.0, 8.0, 5.0),
        periodAngles_Gs_to_Wn_ixyz=(1.5, 1.5, 1.5),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return angles_only_core_wing_movement_fixture


def make_phase_offset_Ler_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with non-zero phase offset
    for Ler_Gs_Cgs.

    :return phase_offset_Ler_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with phase offset for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the phase-offset CoreWingMovement.
    phase_offset_Ler_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.1, 0.08, 0.06),
        periodLer_Gs_Cgs=(1.0, 1.0, 1.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(90.0, -45.0, 60.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return phase_offset_Ler_core_wing_movement_fixture


def make_phase_offset_angles_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with non-zero phase offset
    for angles_Gs_to_Wn_ixyz.

    :return phase_offset_angles_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with phase offset for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the phase-offset CoreWingMovement.
    phase_offset_angles_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(10.0, 12.0, 8.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 1.0, 1.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(45.0, 90.0, -30.0),
    )

    # Return the CoreWingMovement fixture.
    return phase_offset_angles_core_wing_movement_fixture


def make_multiple_periods_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with different periods
    for different dimensions.

    :return multiple_periods_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with different periods.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_multiple_periods_core_wing_cross_section_movement_fixture(),
    ]

    # Create the multiple-periods CoreWingMovement.
    multiple_periods_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.1, 0.08, 0.06),
        periodLer_Gs_Cgs=(1.0, 2.0, 3.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(8.0, 10.0, 12.0),
        periodAngles_Gs_to_Wn_ixyz=(0.5, 1.5, 2.5),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return multiple_periods_core_wing_movement_fixture


def make_custom_spacing_Ler_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with a custom spacing
    function for Ler_Gs_Cgs.

    :return custom_spacing_Ler_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with custom spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
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

    # Create the custom-spacing CoreWingMovement.
    custom_spacing_Ler_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.15, 0.0, 0.0),
        periodLer_Gs_Cgs=(1.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=(custom_harmonic, "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return custom_spacing_Ler_core_wing_movement_fixture


def make_custom_spacing_angles_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with a custom spacing
    function for angles_Gs_to_Wn_ixyz.

    :return custom_spacing_angles_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with custom spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
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

    # Create the custom-spacing CoreWingMovement.
    custom_spacing_angles_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=(custom_harmonic, "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return custom_spacing_angles_core_wing_movement_fixture


def make_mixed_custom_and_standard_spacing_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with mixed custom and
    standard spacing functions.

    :return mixed_custom_and_standard_spacing_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with mixed custom and standard spacing.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_mixed_custom_and_standard_spacing_core_wing_cross_section_movement_fixture(),
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

    # Create the mixed-spacing CoreWingMovement.
    mixed_custom_and_standard_spacing_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.1, 0.08, 0.06),
        periodLer_Gs_Cgs=(1.0, 1.0, 1.0),
        spacingLer_Gs_Cgs=(custom_harmonic, "uniform", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(8.0, 10.0, 6.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 1.0, 1.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", custom_harmonic, "uniform"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return mixed_custom_and_standard_spacing_core_wing_movement_fixture


def make_rotation_point_offset_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with a non zero rotation
    point offset.

    :return rotation_point_offset_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with a non zero rotationPointOffset_Gs_Ler.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the CoreWingMovement with rotation point offset.
    # The offset is in y direction (0.0, 0.5, 0.0), and we rotate about the x axis.
    # This causes the wing to trace an arc in the yz plane as it rotates.
    rotation_point_offset_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        rotationPointOffset_Gs_Ler=(0.0, 0.5, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return rotation_point_offset_core_wing_movement_fixture


def make_periodic_geometry_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement with periodic geometry motion
    suitable for testing the variable geometry optimization.

    The fixture uses a 0.1s period which aligns well with common delta_time values
    like 0.01s (10 steps per period) and 0.02s (5 steps per period).

    :return periodic_geometry_core_wing_movement_fixture: CoreWingMovement
        This is the CoreWingMovement with periodic geometry motion.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    # Create the periodic geometry CoreWingMovement.
    # Use a 0.1s period for angular motion (plunging wing motion).
    periodic_geometry_core_wing_movement_fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(5.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.1, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingMovement fixture.
    return periodic_geometry_core_wing_movement_fixture


def make_2_chordwise_panels_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement for a Wing with
    2 chordwise panels.

    :return: CoreWingMovement for a Wing with 2 chordwise panels.
    """
    base_wing = geometry_fixtures.make_wing_with_2_chordwise_panels()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    return fixture


def make_3_chordwise_panels_core_wing_movement_fixture():
    """This method makes a fixture that is a CoreWingMovement for a Wing with
    3 chordwise panels.

    :return: CoreWingMovement for a Wing with 3 chordwise panels.
    """
    base_wing = geometry_fixtures.make_wing_with_3_chordwise_panels()
    wcs_movements = [
        core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
        core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
    ]

    fixture = CoreWingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    return fixture
