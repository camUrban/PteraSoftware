"""Contains the LineVortex class."""

from __future__ import annotations

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
    cached. Setting either endpoint position invalidates all cached values, ensuring
    consistency while avoiding redundant computation.
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
        # Declare type annotations and initialize the private cache variables.
        self._Slvp_GP1_CgP1: np.ndarray
        self._Elvp_GP1_CgP1: np.ndarray
        self._vector_GP1: np.ndarray | None = None
        self._Clvp_GP1_CgP1: np.ndarray | None = None

        # Initialize the attributes.
        self.Slvp_GP1_CgP1 = Slvp_GP1_CgP1
        self.Elvp_GP1_CgP1 = Elvp_GP1_CgP1
        self.strength = strength

    @property
    def Slvp_GP1_CgP1(self) -> np.ndarray:
        return self._Slvp_GP1_CgP1

    @Slvp_GP1_CgP1.setter
    def Slvp_GP1_CgP1(self, newSlvp_GP1_CgP1: np.ndarray) -> None:
        self._vector_GP1 = None
        self._Clvp_GP1_CgP1 = None

        self._Slvp_GP1_CgP1 = newSlvp_GP1_CgP1

    @property
    def Elvp_GP1_CgP1(self) -> np.ndarray:
        return self._Elvp_GP1_CgP1

    @Elvp_GP1_CgP1.setter
    def Elvp_GP1_CgP1(self, newElvp_GP1_CgP1: np.ndarray) -> None:
        self._vector_GP1 = None
        self._Clvp_GP1_CgP1 = None

        self._Elvp_GP1_CgP1 = newElvp_GP1_CgP1

    @property
    def vector_GP1(self) -> np.ndarray:
        """The LineVortex's vector from start to end point (in the first Airplane's
        geometry axes).

        :return: A (3,) ndarray of floats representing the vector from the LineVortex's
            start point to its end point (in the first Airplane's geometry axes). The
            units are in meters.
        """
        if self._vector_GP1 is None:
            self._vector_GP1 = self.Elvp_GP1_CgP1 - self.Slvp_GP1_CgP1
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
            self._Clvp_GP1_CgP1 = self.Slvp_GP1_CgP1 + 0.5 * self.vector_GP1
        return self._Clvp_GP1_CgP1
