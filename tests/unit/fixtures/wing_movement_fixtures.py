"""This module contains functions to create WingMovements for use in tests."""

import numpy as np

import pterasoftware as ps

from . import geometry_fixtures, wing_cross_section_movement_fixtures


def make_static_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with all parameters
    zero (no movement).

    :return static_wing_movement_fixture: WingMovement
        This is the WingMovement with no movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the static WingMovement.
    static_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return static_wing_movement_fixture


def make_basic_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with general-purpose
    moderate values.

    :return basic_wing_movement_fixture: WingMovement
        This is the WingMovement with general-purpose values.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_basic_wing_cross_section_movement_fixture(),
    ]

    # Create the basic WingMovement.
    basic_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return basic_wing_movement_fixture


def make_sine_spacing_Ler_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with sine spacing for
    Ler_Gs_Cgs.

    :return sine_spacing_Ler_wing_movement_fixture: WingMovement
        This is the WingMovement with sine spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with sine spacing for Ler_Gs_Cgs.
    sine_spacing_Ler_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return sine_spacing_Ler_wing_movement_fixture


def make_uniform_spacing_Ler_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with uniform spacing for
    Ler_Gs_Cgs.

    :return uniform_spacing_Ler_wing_movement_fixture: WingMovement
        This is the WingMovement with uniform spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with uniform spacing for Ler_Gs_Cgs.
    uniform_spacing_Ler_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return uniform_spacing_Ler_wing_movement_fixture


def make_mixed_spacing_Ler_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with mixed spacing for
    Ler_Gs_Cgs.

    :return mixed_spacing_Ler_wing_movement_fixture: WingMovement
        This is the WingMovement with mixed spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with mixed spacing for Ler_Gs_Cgs.
    mixed_spacing_Ler_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return mixed_spacing_Ler_wing_movement_fixture


def make_sine_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with sine spacing for
    angles_Gs_to_Wn_ixyz.

    :return sine_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with sine spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with sine spacing for angles_Gs_to_Wn_ixyz.
    sine_spacing_angles_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return sine_spacing_angles_wing_movement_fixture


def make_uniform_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with uniform spacing for
    angles_Gs_to_Wn_ixyz.

    :return uniform_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with uniform spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with uniform spacing for angles_Gs_to_Wn_ixyz.
    uniform_spacing_angles_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
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
    )

    # Return the WingMovement fixture.
    return uniform_spacing_angles_wing_movement_fixture


def make_mixed_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with mixed spacing for
    angles_Gs_to_Wn_ixyz.

    :return mixed_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with mixed spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with mixed spacing for angles_Gs_to_Wn_ixyz.
    mixed_spacing_angles_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
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
    )

    # Return the WingMovement fixture.
    return mixed_spacing_angles_wing_movement_fixture


def make_Ler_only_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement where only Ler_Gs_Cgs
    moves.

    :return Ler_only_wing_movement_fixture: WingMovement
        This is the WingMovement with only Ler_Gs_Cgs movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the Ler-only WingMovement.
    Ler_only_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return Ler_only_wing_movement_fixture


def make_angles_only_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement where only
    angles_Gs_to_Wn_ixyz moves.

    :return angles_only_wing_movement_fixture: WingMovement
        This is the WingMovement with only angles_Gs_to_Wn_ixyz movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the angles-only WingMovement.
    angles_only_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return angles_only_wing_movement_fixture


def make_phase_offset_Ler_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with non-zero phase offset
    for Ler_Gs_Cgs.

    :return phase_offset_Ler_wing_movement_fixture: WingMovement
        This is the WingMovement with phase offset for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the phase-offset WingMovement.
    phase_offset_Ler_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return phase_offset_Ler_wing_movement_fixture


def make_phase_offset_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with non-zero phase offset
    for angles_Gs_to_Wn_ixyz.

    :return phase_offset_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with phase offset for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the phase-offset WingMovement.
    phase_offset_angles_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return phase_offset_angles_wing_movement_fixture


def make_multiple_periods_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with different periods
    for different dimensions.

    :return multiple_periods_wing_movement_fixture: WingMovement
        This is the WingMovement with different periods.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_multiple_periods_wing_cross_section_movement_fixture(),
    ]

    # Create the multiple-periods WingMovement.
    multiple_periods_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return multiple_periods_wing_movement_fixture


def make_custom_spacing_Ler_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with a custom spacing
    function for Ler_Gs_Cgs.

    :return custom_spacing_Ler_wing_movement_fixture: WingMovement
        This is the WingMovement with custom spacing for Ler_Gs_Cgs.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
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

    # Create the custom-spacing WingMovement.
    custom_spacing_Ler_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
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

    # Return the WingMovement fixture.
    return custom_spacing_Ler_wing_movement_fixture


def make_custom_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with a custom spacing
    function for angles_Gs_to_Wn_ixyz.

    :return custom_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with custom spacing for angles_Gs_to_Wn_ixyz.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
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

    # Create the custom-spacing WingMovement.
    custom_spacing_angles_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
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
    )

    # Return the WingMovement fixture.
    return custom_spacing_angles_wing_movement_fixture


def make_mixed_custom_and_standard_spacing_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with mixed custom and
    standard spacing functions.

    :return mixed_custom_and_standard_spacing_wing_movement_fixture: WingMovement
        This is the WingMovement with mixed custom and standard spacing.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_mixed_custom_and_standard_spacing_wing_cross_section_movement_fixture(),
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

    # Create the mixed-spacing WingMovement.
    mixed_custom_and_standard_spacing_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
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
    )

    # Return the WingMovement fixture.
    return mixed_custom_and_standard_spacing_wing_movement_fixture
