"""This module contains functions to create HorseshoeVortices for use in tests."""

import numpy as np
import pterasoftware as ps


def make_basic_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with basic
    configuration for testing purposes.

    :return basic_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex configured for basic testing. The finite leg is 1
        meter long in the y-direction (in geometry axes), and the quasi-infinite legs
        extend 20 meters in the positive x-direction (in geometry axes) with unit
        strength.
    """
    basic_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )

    return basic_horseshoe_vortex_fixture


def make_short_legs_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with short
    quasi-infinite legs for testing purposes.

    :return short_legs_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with relatively short quasi-infinite legs
        (5 meters) for testing edge cases.
    """
    short_legs_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=5.0,
        strength=1.0,
    )

    return short_legs_horseshoe_vortex_fixture


def make_long_legs_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with very long
    quasi-infinite legs for testing purposes.

    :return long_legs_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with very long quasi-infinite legs (100 meters)
        to better approximate infinite legs.
    """
    long_legs_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=100.0,
        strength=1.0,
    )

    return long_legs_horseshoe_vortex_fixture


def make_tilted_legs_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with tilted
    quasi-infinite legs for testing purposes.

    :return tilted_legs_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with quasi-infinite legs tilted at an angle
        and extending in a direction with x, y, and z components (in geometry axes).
    """
    tilted_legs_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([0.707, 0.0, 0.707], dtype=float),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )

    return tilted_legs_horseshoe_vortex_fixture


def make_wide_finite_leg_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with a wide finite
    leg for testing purposes.

    :return wide_finite_leg_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with a 5-meter finite leg span.
    """
    wide_finite_leg_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 2.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -2.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )

    return wide_finite_leg_horseshoe_vortex_fixture


def make_narrow_finite_leg_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with a narrow finite
    leg for testing purposes.

    :return narrow_finite_leg_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with a 0.1-meter finite leg span.
    """
    narrow_finite_leg_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.05, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.05, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )

    return narrow_finite_leg_horseshoe_vortex_fixture


def make_zero_strength_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with zero strength
    for testing edge cases.

    :return zero_strength_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with zero strength.
    """
    zero_strength_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=0.0,
    )

    return zero_strength_horseshoe_vortex_fixture


def make_negative_strength_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with negative strength
    for testing purposes.

    :return negative_strength_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with negative strength of -1.0.
    """
    negative_strength_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=-1.0,
    )

    return negative_strength_horseshoe_vortex_fixture


def make_high_strength_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex with high strength
    for testing purposes.

    :return high_strength_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex with high strength of 100.0.
    """
    high_strength_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=100.0,
    )

    return high_strength_horseshoe_vortex_fixture


def make_offset_horseshoe_vortex_fixture():
    """This method makes a fixture that is a HorseshoeVortex offset from the
    origin for testing purposes.

    :return offset_horseshoe_vortex_fixture: HorseshoeVortex
        This is the HorseshoeVortex offset by [10.0, 5.0, 3.0] (in geometry axes,
        relative to the CG) with unit strength.
    """
    offset = np.array([10.0, 5.0, 3.0], dtype=float)
    offset_horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        Frhvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float) + offset,
        Flhvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float) + offset,
        leftLegVector_G=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )

    return offset_horseshoe_vortex_fixture
