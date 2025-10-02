"""This module contains functions to create RingVortices for use in tests."""

import numpy as np
import pterasoftware as ps


def make_basic_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex with basic rectangular
    geometry for testing purposes.

    :return basic_ring_vortex_fixture: RingVortex
        This is the RingVortex configured for basic testing. The vortex forms a unit
        square in the xz-plane (in geometry axes) with unit strength.
    """
    basic_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
        Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
        strength=1.0,
    )

    return basic_ring_vortex_fixture


def make_unit_square_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex forming a unit square
    in the xy-plane for testing purposes.

    :return unit_square_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a 1x1 square in the xy-plane centered at the
        origin (in geometry axes, relative to the CG) with unit strength.
    """
    unit_square_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.5, 0.5, 0.0], dtype=float),
        Flrvp_G_Cg=np.array([0.5, -0.5, 0.0], dtype=float),
        Blrvp_G_Cg=np.array([-0.5, -0.5, 0.0], dtype=float),
        Brrvp_G_Cg=np.array([-0.5, 0.5, 0.0], dtype=float),
        strength=1.0,
    )

    return unit_square_ring_vortex_fixture


def make_rectangular_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex forming a rectangle
    for testing purposes.

    :return rectangular_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a 2x1 rectangle in the xz-plane (in geometry
        axes) with strength of 2.5.
    """
    rectangular_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        Blrvp_G_Cg=np.array([2.0, -0.5, 0.0], dtype=float),
        Brrvp_G_Cg=np.array([2.0, 0.5, 0.0], dtype=float),
        strength=2.5,
    )

    return rectangular_ring_vortex_fixture


def make_tilted_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex tilted out of the
    principal planes for testing purposes.

    :return tilted_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a tilted quadrilateral with strength of 1.5.
    """
    tilted_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.0, 0.5, 0.2], dtype=float),
        Flrvp_G_Cg=np.array([0.0, -0.5, 0.1], dtype=float),
        Blrvp_G_Cg=np.array([1.0, -0.5, -0.1], dtype=float),
        Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
        strength=1.5,
    )

    return tilted_ring_vortex_fixture


def make_zero_strength_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex with zero strength for
    testing edge cases.

    :return zero_strength_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a unit square with zero strength.
    """
    zero_strength_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
        Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
        strength=0.0,
    )

    return zero_strength_ring_vortex_fixture


def make_negative_strength_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex with negative strength
    for testing purposes.

    :return negative_strength_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a unit square with negative strength of
        -1.0.
    """
    negative_strength_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
        Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
        Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
        strength=-1.0,
    )

    return negative_strength_ring_vortex_fixture


def make_offset_ring_vortex_fixture():
    """This method makes a fixture that is a RingVortex offset from the origin
    for testing purposes.

    :return offset_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a unit square offset by [5.0, 3.0, 2.0]
        (in geometry axes, relative to the CG) with unit strength.
    """
    offset = np.array([5.0, 3.0, 2.0], dtype=float)
    offset_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float) + offset,
        Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float) + offset,
        Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float) + offset,
        Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float) + offset,
        strength=1.0,
    )

    return offset_ring_vortex_fixture


def make_small_ring_vortex_fixture():
    """This method makes a fixture that is a very small RingVortex for testing
    purposes.

    :return small_ring_vortex_fixture: RingVortex
        This is the RingVortex forming a 0.01x0.01 square in the xy-plane (in
        geometry axes) with unit strength.
    """
    small_ring_vortex_fixture = ps.aerodynamics.RingVortex(
        Frrvp_G_Cg=np.array([0.005, 0.005, 0.0], dtype=float),
        Flrvp_G_Cg=np.array([0.005, -0.005, 0.0], dtype=float),
        Blrvp_G_Cg=np.array([-0.005, -0.005, 0.0], dtype=float),
        Brrvp_G_Cg=np.array([-0.005, 0.005, 0.0], dtype=float),
        strength=1.0,
    )

    return small_ring_vortex_fixture
