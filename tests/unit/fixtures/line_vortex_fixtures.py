"""This module contains functions to create LineVortices for use in tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware._vortices import _line_vortex


def make_basic_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex with basic configuration for
    testing purposes.

    :return basic_line_vortex_fixture: LineVortex
        This is the LineVortex configured for basic testing. The vortex extends from
        [0.0, 0.0, 0.0] to [1.0, 0.0, 0.0] (in the first Airplane's geometry axes,
        relative to the first Airplane's CG) with unit strength.
    """
    basic_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
        strength=1.0,
    )

    return basic_line_vortex_fixture


def make_diagonal_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex oriented diagonally for testing
    purposes.

    :return diagonal_line_vortex_fixture: LineVortex
        This is the LineVortex extending diagonally from [0.0, 0.0, 0.0] to [1.0, 1.0,
        1.0] (in the first Airplane's geometry axes, relative to the first Airplane's
        CG) with unit strength.
    """
    diagonal_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1.0, 1.0, 1.0], dtype=float),
        strength=1.0,
    )

    return diagonal_line_vortex_fixture


def make_y_aligned_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex aligned with the y axis for
    testing purposes.

    :return y_aligned_line_vortex_fixture: LineVortex
        This is the LineVortex extending along the y axis from [0.0, -0.5, 0.0] to
        [0.0, 0.5, 0.0] (in the first Airplane's geometry axes, relative to the first
        Airplane's CG) with unit strength.
    """
    y_aligned_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, -0.5, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([0.0, 0.5, 0.0], dtype=float),
        strength=1.0,
    )

    return y_aligned_line_vortex_fixture


def make_z_aligned_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex aligned with the z axis for
    testing purposes.

    :return z_aligned_line_vortex_fixture: LineVortex
        This is the LineVortex extending along the z axis from [0.0, 0.0, 0.0] to
        [0.0, 0.0, 2.0] (in the first Airplane's geometry axes, relative to the first
        Airplane's CG) with strength of 3.0.
    """
    z_aligned_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([0.0, 0.0, 2.0], dtype=float),
        strength=3.0,
    )

    return z_aligned_line_vortex_fixture


def make_zero_strength_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex with zero strength for testing
    edge cases.

    :return zero_strength_line_vortex_fixture: LineVortex
        This is the LineVortex with zero strength.
    """
    zero_strength_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
        strength=0.0,
    )

    return zero_strength_line_vortex_fixture


def make_negative_strength_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex with negative strength for
    testing purposes.

    :return negative_strength_line_vortex_fixture: LineVortex
        This is the LineVortex with negative strength of -2.5.
    """
    negative_strength_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
        strength=-2.5,
    )

    return negative_strength_line_vortex_fixture


def make_large_strength_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex with very large strength for
    testing purposes.

    :return large_strength_line_vortex_fixture: LineVortex
        This is the LineVortex with large strength of 1e6.
    """
    large_strength_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
        strength=1e6,
    )

    return large_strength_line_vortex_fixture


def make_offset_line_vortex_fixture():
    """This method makes a fixture that is a LineVortex offset from the origin for
    testing purposes.

    :return offset_line_vortex_fixture: LineVortex
        This is the LineVortex offset by [5.0, 3.0, 2.0] (in the first Airplane's
        geometry axes, relative to the first Airplane's CG) with unit strength.
    """
    offset = np.array([5.0, 3.0, 2.0], dtype=float)
    offset_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float) + offset,
        Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float) + offset,
        strength=1.0,
    )

    return offset_line_vortex_fixture


def make_small_line_vortex_fixture():
    """This method makes a fixture that is a very small LineVortex for testing
    purposes.

    :return small_line_vortex_fixture: LineVortex
        This is the LineVortex with length of 0.001 meters.
    """
    small_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([0.001, 0.0, 0.0], dtype=float),
        strength=1.0,
    )

    return small_line_vortex_fixture


def make_long_line_vortex_fixture():
    """This method makes a fixture that is a very long LineVortex for testing
    purposes.

    :return long_line_vortex_fixture: LineVortex
        This is the LineVortex with length of 1000 meters.
    """
    long_line_vortex_fixture = _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1000.0, 0.0, 0.0], dtype=float),
        strength=1.0,
    )

    return long_line_vortex_fixture
