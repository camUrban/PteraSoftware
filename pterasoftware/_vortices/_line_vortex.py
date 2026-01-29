"""Contains the LineVortex class."""

from __future__ import annotations

from typing import cast

import numpy as np


class LineVortex:
    """A class used to contain line vortices.

    **Contains the following methods:**

    vector_GP1: The LineVortex's vector from start to end point (in the first Airplane's
    geometry axes).

    Clvp_GP1_CgP1: The position of the LineVortex's center point (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    **Notes:**

    Computed geometric properties (vector and center point) are lazily evaluated and
    cached.
    """

    def __init__(
        self, Slvp_GP1_CgP1: np.ndarray, Elvp_GP1_CgP1: np.ndarray, strength: float
    ) -> None:
        """The initialization method.

        :param Slvp_GP1_CgP1: A (3,) ndarray of 3 floats representing the position of
            the LineVortex's start point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The units are in meters.
        :param Elvp_GP1_CgP1: A (3,) ndarray of 3 floats representing the position of
            the LineVortex's end point (in the first Airplane's geometry axes, relative
            to the first Airplane's CG). The units are in meters.
        :param strength: The strength of the LineVortex. Its units are in meters squared
            per second.
        :return: None
        """
        # Initialize the immutable attributes. Set those that are numpy arrays to be
        # read only.
        self._Slvp_GP1_CgP1 = Slvp_GP1_CgP1
        self._Slvp_GP1_CgP1.flags.writeable = False
        self._Elvp_GP1_CgP1 = Elvp_GP1_CgP1
        self._Elvp_GP1_CgP1.flags.writeable = False

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._vector_GP1: np.ndarray | None = None
        self._Clvp_GP1_CgP1: np.ndarray | None = None

        # Initialize a mutable attribute to hold the LineVortex's strength.
        self.strength = strength

    # --- Immutable: read only properties ---
    @property
    def Slvp_GP1_CgP1(self) -> np.ndarray:
        return self._Slvp_GP1_CgP1

    @property
    def Elvp_GP1_CgP1(self) -> np.ndarray:
        return self._Elvp_GP1_CgP1

    # --- Immutable derived: manual lazy caching ---
    @property
    def vector_GP1(self) -> np.ndarray:
        """The LineVortex's vector from start to end point (in the first Airplane's
        geometry axes).

        :return: A (3,) ndarray of floats representing the vector from the LineVortex's
            start point to its end point (in the first Airplane's geometry axes). The
            units are in meters.
        """
        if self._vector_GP1 is None:
            self._vector_GP1 = cast(
                np.ndarray, self._Elvp_GP1_CgP1 - self._Slvp_GP1_CgP1
            )
            self._vector_GP1.flags.writeable = False
        return self._vector_GP1

    @property
    def Clvp_GP1_CgP1(self) -> np.ndarray:
        """The position of the LineVortex's center point (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the LineVortex's
            center point (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters.
        """
        if self._Clvp_GP1_CgP1 is None:
            self._Clvp_GP1_CgP1 = self._Slvp_GP1_CgP1 + 0.5 * self.vector_GP1
            self._Clvp_GP1_CgP1.flags.writeable = False
        return self._Clvp_GP1_CgP1
