"""This module contains functions to create WingMovements for use in tests."""

import numpy as np
import pterasoftware as ps

from . import geometry_fixtures
from . import wing_cross_section_movement_fixtures


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
        ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
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
        ampPrelimLer_G_Cg=(0.1, 0.05, 0.08),
        periodPrelimLer_G_Cg=(2.0, 2.0, 2.0),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(5.0, 3.0, 2.0),
        periodAngles_G_to_prelimWn_izyx=(2.0, 2.0, 2.0),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return basic_wing_movement_fixture


def make_sine_spacing_prelimLer_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with sine spacing for
    prelimLer_G_Cg.

    :return sine_spacing_prelimLer_wing_movement_fixture: WingMovement
        This is the WingMovement with sine spacing for prelimLer_G_Cg.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with sine spacing for prelimLer_G_Cg.
    sine_spacing_prelimLer_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.2, 0.0, 0.0),
            periodPrelimLer_G_Cg=(1.0, 0.0, 0.0),
            spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return sine_spacing_prelimLer_wing_movement_fixture


def make_uniform_spacing_prelimLer_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with uniform spacing for
    prelimLer_G_Cg.

    :return uniform_spacing_prelimLer_wing_movement_fixture: WingMovement
        This is the WingMovement with uniform spacing for prelimLer_G_Cg.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with uniform spacing for prelimLer_G_Cg.
    uniform_spacing_prelimLer_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.2, 0.0, 0.0),
            periodPrelimLer_G_Cg=(1.0, 0.0, 0.0),
            spacingPrelimLer_G_Cg=("uniform", "uniform", "uniform"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return uniform_spacing_prelimLer_wing_movement_fixture


def make_mixed_spacing_prelimLer_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with mixed spacing for
    prelimLer_G_Cg.

    :return mixed_spacing_prelimLer_wing_movement_fixture: WingMovement
        This is the WingMovement with mixed spacing for prelimLer_G_Cg.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with mixed spacing for prelimLer_G_Cg.
    mixed_spacing_prelimLer_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.2, 0.15, 0.1),
            periodPrelimLer_G_Cg=(1.0, 1.0, 1.0),
            spacingPrelimLer_G_Cg=("sine", "uniform", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return mixed_spacing_prelimLer_wing_movement_fixture


def make_sine_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with sine spacing for
    angles_G_to_prelimWn_izyx.

    :return sine_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with sine spacing for angles_G_to_prelimWn_izyx.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with sine spacing for angles_G_to_prelimWn_izyx.
    sine_spacing_angles_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(10.0, 0.0, 0.0),
        periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return sine_spacing_angles_wing_movement_fixture


def make_uniform_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with uniform spacing for
    angles_G_to_prelimWn_izyx.

    :return uniform_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with uniform spacing for angles_G_to_prelimWn_izyx.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with uniform spacing for angles_G_to_prelimWn_izyx.
    uniform_spacing_angles_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
            periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
            spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(10.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=("uniform", "uniform", "uniform"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return uniform_spacing_angles_wing_movement_fixture


def make_mixed_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with mixed spacing for
    angles_G_to_prelimWn_izyx.

    :return mixed_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with mixed spacing for angles_G_to_prelimWn_izyx.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the WingMovement with mixed spacing for angles_G_to_prelimWn_izyx.
    mixed_spacing_angles_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
            periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
            spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(10.0, 15.0, 8.0),
            periodAngles_G_to_prelimWn_izyx=(1.0, 1.0, 1.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", "uniform", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return mixed_spacing_angles_wing_movement_fixture


def make_prelimLer_only_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement where only prelimLer_G_Cg
    moves.

    :return prelimLer_only_wing_movement_fixture: WingMovement
        This is the WingMovement with only prelimLer_G_Cg movement.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the prelimLer-only WingMovement.
    prelimLer_only_wing_movement_fixture = ps.movements.wing_movement.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wcs_movements,
        ampPrelimLer_G_Cg=(0.15, 0.1, 0.08),
        periodPrelimLer_G_Cg=(1.5, 1.5, 1.5),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return prelimLer_only_wing_movement_fixture


def make_angles_only_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement where only
    angles_G_to_prelimWn_izyx moves.

    :return angles_only_wing_movement_fixture: WingMovement
        This is the WingMovement with only angles_G_to_prelimWn_izyx movement.
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
        ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(12.0, 8.0, 5.0),
        periodAngles_G_to_prelimWn_izyx=(1.5, 1.5, 1.5),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return angles_only_wing_movement_fixture


def make_phase_offset_prelimLer_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with non-zero phase offset
    for prelimLer_G_Cg.

    :return phase_offset_prelimLer_wing_movement_fixture: WingMovement
        This is the WingMovement with phase offset for prelimLer_G_Cg.
    """
    # Initialize the constructing fixtures.
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wcs_movements = [
        wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
        wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
    ]

    # Create the phase-offset WingMovement.
    phase_offset_prelimLer_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.1, 0.08, 0.06),
            periodPrelimLer_G_Cg=(1.0, 1.0, 1.0),
            spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
            phasePrelimLer_G_Cg=(90.0, -45.0, 60.0),
            ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return phase_offset_prelimLer_wing_movement_fixture


def make_phase_offset_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with non-zero phase offset
    for angles_G_to_prelimWn_izyx.

    :return phase_offset_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with phase offset for angles_G_to_prelimWn_izyx.
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
        ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(10.0, 12.0, 8.0),
        periodAngles_G_to_prelimWn_izyx=(1.0, 1.0, 1.0),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(45.0, 90.0, -30.0),
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
        ampPrelimLer_G_Cg=(0.1, 0.08, 0.06),
        periodPrelimLer_G_Cg=(1.0, 2.0, 3.0),
        spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
        phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
        ampAngles_G_to_prelimWn_izyx=(8.0, 10.0, 12.0),
        periodAngles_G_to_prelimWn_izyx=(0.5, 1.5, 2.5),
        spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
        phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
    )

    # Return the WingMovement fixture.
    return multiple_periods_wing_movement_fixture


def make_custom_spacing_prelimLer_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with a custom spacing
    function for prelimLer_G_Cg.

    :return custom_spacing_prelimLer_wing_movement_fixture: WingMovement
        This is the WingMovement with custom spacing for prelimLer_G_Cg.
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
    custom_spacing_prelimLer_wing_movement_fixture = (
        ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.15, 0.0, 0.0),
            periodPrelimLer_G_Cg=(1.0, 0.0, 0.0),
            spacingPrelimLer_G_Cg=(custom_harmonic, "sine", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", "sine", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return custom_spacing_prelimLer_wing_movement_fixture


def make_custom_spacing_angles_wing_movement_fixture():
    """This method makes a fixture that is a WingMovement with a custom spacing
    function for angles_G_to_prelimWn_izyx.

    :return custom_spacing_angles_wing_movement_fixture: WingMovement
        This is the WingMovement with custom spacing for angles_G_to_prelimWn_izyx.
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
            ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
            periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
            spacingPrelimLer_G_Cg=("sine", "sine", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(10.0, 0.0, 0.0),
            periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
            spacingAngles_G_to_prelimWn_izyx=(custom_harmonic, "sine", "sine"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
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
            ampPrelimLer_G_Cg=(0.1, 0.08, 0.06),
            periodPrelimLer_G_Cg=(1.0, 1.0, 1.0),
            spacingPrelimLer_G_Cg=(custom_harmonic, "uniform", "sine"),
            phasePrelimLer_G_Cg=(0.0, 0.0, 0.0),
            ampAngles_G_to_prelimWn_izyx=(8.0, 10.0, 6.0),
            periodAngles_G_to_prelimWn_izyx=(1.0, 1.0, 1.0),
            spacingAngles_G_to_prelimWn_izyx=("sine", custom_harmonic, "uniform"),
            phaseAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
        )
    )

    # Return the WingMovement fixture.
    return mixed_custom_and_standard_spacing_wing_movement_fixture
