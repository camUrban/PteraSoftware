"""This module contains functions to create CoreWingCrossSectionMovements for use in
tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware._core import CoreWingCrossSectionMovement

from . import geometry_fixtures


def make_sine_spacing_Lp_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with sine
    spacing for Lp_Wcsp_Lpp.

    :return sine_spacing_Lp_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with sine spacing for Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the CoreWingCrossSectionMovement with sine spacing.
    sine_spacing_Lp_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return sine_spacing_Lp_core_wing_cross_section_movement_fixture


def make_uniform_spacing_Lp_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with uniform
    spacing for Lp_Wcsp_Lpp.

    :return uniform_spacing_Lp_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with uniform spacing for
        Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the CoreWingCrossSectionMovement with uniform spacing.
    uniform_spacing_Lp_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("uniform", "uniform", "uniform"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return uniform_spacing_Lp_core_wing_cross_section_movement_fixture


def make_mixed_spacing_Lp_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with mixed
    spacing for Lp_Wcsp_Lpp.

    :return mixed_spacing_Lp_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with mixed spacing for Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    # Use tip fixture which has Lp_Wcsp_Lpp[1] = 2.0, allowing for amplitude of 1.5.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the CoreWingCrossSectionMovement with mixed spacing.
    mixed_spacing_Lp_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(1.0, 1.5, 0.5),
            periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
            spacingLp_Wcsp_Lpp=("sine", "uniform", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return mixed_spacing_Lp_core_wing_cross_section_movement_fixture


def make_sine_spacing_angles_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with sine
    spacing for angles_Wcsp_to_Wcs_ixyz.

    :return sine_spacing_angles_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with sine spacing for
        angles_Wcsp_to_Wcs_ixyz.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the CoreWingCrossSectionMovement with sine spacing for angles.
    sine_spacing_angles_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return sine_spacing_angles_core_wing_cross_section_movement_fixture


def make_uniform_spacing_angles_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with uniform
    spacing for angles_Wcsp_to_Wcs_ixyz.

    :return uniform_spacing_angles_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with uniform spacing for
        angles_Wcsp_to_Wcs_ixyz.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the CoreWingCrossSectionMovement with uniform spacing for angles.
    uniform_spacing_angles_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("uniform", "uniform", "uniform"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return uniform_spacing_angles_core_wing_cross_section_movement_fixture


def make_mixed_spacing_angles_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with mixed
    spacing for angles_Wcsp_to_Wcs_ixyz.

    :return mixed_spacing_angles_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with mixed spacing for
        angles_Wcsp_to_Wcs_ixyz.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the CoreWingCrossSectionMovement with mixed spacing for angles.
    mixed_spacing_angles_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 20.0, 5.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 1.0, 1.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "uniform", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return mixed_spacing_angles_core_wing_cross_section_movement_fixture


def make_static_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with all
    parameters zero (no movement).

    :return static_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with no movement.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the static CoreWingCrossSectionMovement.
    static_core_wing_cross_section_movement_fixture = CoreWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return static_core_wing_cross_section_movement_fixture


def make_static_tip_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with all
    parameters zero (no movement), using a tip WingCrossSection as the base.

    :return static_tip_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with no movement for a tip cross
        section.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the static tip CoreWingCrossSectionMovement.
    static_tip_core_wing_cross_section_movement_fixture = CoreWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return static_tip_core_wing_cross_section_movement_fixture


def make_basic_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with
    general-purpose moderate values.

    :return basic_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with general-purpose values.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the basic CoreWingCrossSectionMovement.
    basic_core_wing_cross_section_movement_fixture = CoreWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.4, 0.3, 0.15),
        periodLp_Wcsp_Lpp=(2.0, 2.0, 2.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(15.0, 10.0, 5.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(2.0, 2.0, 2.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return basic_core_wing_cross_section_movement_fixture


def make_Lp_only_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement where only
    Lp_Wcsp_Lpp moves.

    :return Lp_only_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with only Lp_Wcsp_Lpp movement.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the Lp-only CoreWingCrossSectionMovement.
    Lp_only_core_wing_cross_section_movement_fixture = CoreWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.4, 0.5, 0.15),
        periodLp_Wcsp_Lpp=(1.5, 1.5, 1.5),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return Lp_only_core_wing_cross_section_movement_fixture


def make_angles_only_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement where only
    angles_Wcsp_to_Wcs_ixyz moves.

    :return angles_only_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with only angles_Wcsp_to_Wcs_ixyz
        movement.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the angles-only CoreWingCrossSectionMovement.
    angles_only_core_wing_cross_section_movement_fixture = CoreWingCrossSectionMovement(
        base_wing_cross_section=base_wing_cross_section,
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(20.0, 15.0, 10.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(1.5, 1.5, 1.5),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return angles_only_core_wing_cross_section_movement_fixture


def make_phase_offset_Lp_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with
    non-zero phase offset for Lp_Wcsp_Lpp.

    :return phase_offset_Lp_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with phase offset for Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the phase-offset CoreWingCrossSectionMovement.
    phase_offset_Lp_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(0.4, 0.3, 0.15),
            periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(90.0, -90.0, 45.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return phase_offset_Lp_core_wing_cross_section_movement_fixture


def make_phase_offset_angles_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with
    non-zero phase offset for angles_Wcsp_to_Wcs_ixyz.

    :return phase_offset_angles_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with phase offset for
        angles_Wcsp_to_Wcs_ixyz.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

    # Create the phase-offset CoreWingCrossSectionMovement.
    phase_offset_angles_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 15.0, 20.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 1.0, 1.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(45.0, 90.0, -45.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return phase_offset_angles_core_wing_cross_section_movement_fixture


def make_multiple_periods_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with
    different periods for different dimensions.

    :return multiple_periods_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with different periods.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

    # Create the multiple-periods CoreWingCrossSectionMovement.
    multiple_periods_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(0.4, 0.4, 0.15),
            periodLp_Wcsp_Lpp=(1.0, 2.0, 3.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 15.0, 20.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.5, 1.5, 2.5),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return multiple_periods_core_wing_cross_section_movement_fixture


def make_custom_spacing_Lp_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with a
    custom spacing function for Lp_Wcsp_Lpp.

    :return custom_spacing_Lp_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with custom spacing for Lp_Wcsp_Lpp.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

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

    # Create the custom-spacing CoreWingCrossSectionMovement.
    custom_spacing_Lp_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(0.4, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=(custom_harmonic, "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return custom_spacing_Lp_core_wing_cross_section_movement_fixture


def make_custom_spacing_angles_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with a
    custom spacing function for angles_Wcsp_to_Wcs_ixyz.

    :return custom_spacing_angles_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with custom spacing for
        angles_Wcsp_to_Wcs_ixyz.
    """
    # Initialize the constructing fixture.
    base_wing_cross_section = geometry_fixtures.make_root_wing_cross_section_fixture()

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

    # Create the custom-spacing CoreWingCrossSectionMovement.
    custom_spacing_angles_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=(custom_harmonic, "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return custom_spacing_angles_core_wing_cross_section_movement_fixture


def make_mixed_custom_and_standard_spacing_core_wing_cross_section_movement_fixture():
    """This method makes a fixture that is a CoreWingCrossSectionMovement with mixed
    custom and standard spacing functions.

    :return mixed_custom_and_standard_spacing_core_wing_cross_section_movement_fixture: CoreWingCrossSectionMovement
        This is the CoreWingCrossSectionMovement with mixed custom and standard
        spacing.
    """
    # Initialize the constructing fixture.
    # Use tip fixture to ensure Lp values stay non-negative during oscillation.
    base_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()

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

    # Create the mixed-spacing CoreWingCrossSectionMovement.
    mixed_custom_and_standard_spacing_core_wing_cross_section_movement_fixture = (
        CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=(0.4, 0.3, 0.15),
            periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
            spacingLp_Wcsp_Lpp=(custom_harmonic, "uniform", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ampAngles_Wcsp_to_Wcs_ixyz=(15.0, 10.0, 5.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 1.0, 1.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", custom_harmonic, "uniform"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )
    )

    # Return the CoreWingCrossSectionMovement fixture.
    return mixed_custom_and_standard_spacing_core_wing_cross_section_movement_fixture
