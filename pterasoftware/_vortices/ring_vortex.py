"""Contains the RingVortex class."""

from __future__ import annotations

import numpy as np

from pterasoftware import _functions

from . import _line_vortex


class RingVortex:
    """A class used to contain ring vortices.

    **Contains the following methods:**

    Crvp_GP1_CgP1: The position of the RingVortex's centroid (in the first Airplane's
    geometry axes, relative to the first Airplane's CG).

    front_leg: The LineVortex representing this RingVortex's front leg.

    left_leg: The LineVortex representing this RingVortex's left leg.

    back_leg: The LineVortex representing this RingVortex's back leg.

    right_leg: The LineVortex representing this RingVortex's right leg.

    area: An estimate of this RingVortex's area.

    **Notes:**

    Computed geometric properties (LineVortex legs, centroid, and area) are lazily
    evaluated and cached. Changing the strength propagates the new strength to the four
    LineVortex legs (if they've be created and cached).
    """

    __slots__ = (
        "_Flrvp_GP1_CgP1",
        "_Frrvp_GP1_CgP1",
        "_Blrvp_GP1_CgP1",
        "_Brrvp_GP1_CgP1",
        "_Crvp_GP1_CgP1",
        "_front_leg",
        "_left_leg",
        "_back_leg",
        "_right_leg",
        "_area",
        "_strength",
        "age",
    )

    def __init__(
        self,
        Frrvp_GP1_CgP1: np.ndarray,
        Flrvp_GP1_CgP1: np.ndarray,
        Blrvp_GP1_CgP1: np.ndarray,
        Brrvp_GP1_CgP1: np.ndarray,
        strength: float,
    ) -> None:
        """The initialization method.

        :param Frrvp_GP1_CgP1: A (3,) ndarray of floats representing the position of the
            RingVortex's front right point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The front right point is defined as
            the end point of the RingVortex's right leg and the start point of its front
            leg. The units are in meters.
        :param Flrvp_GP1_CgP1: A (3,) ndarray of floats representing the position of the
            RingVortex's front left point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The front left point is defined as the
            end point of the RingVortex's front leg and the start point of its left leg.
            The units are in meters.
        :param Blrvp_GP1_CgP1: A (3,) ndarray of floats representing the position of the
            RingVortex's back left point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The back left point is defined as the
            end point of the RingVortex's left leg and the start point of its back leg.
            The units are in meters.
        :param Brrvp_GP1_CgP1: A (3,) ndarray of floats representing the position of the
            RingVortex's back right point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The back right point is defined as the
            end point of the RingVortex's back leg and the start point of its right leg.
            The units are in meters.
        :param strength: The strength of the RingVortex. Its units are in meters squared
            per second.
        :return: None
        """
        # Initialize the immutable attributes. Set those that are numpy arrays to be
        # read only.
        self._Flrvp_GP1_CgP1 = Flrvp_GP1_CgP1
        self._Flrvp_GP1_CgP1.flags.writeable = False
        self._Frrvp_GP1_CgP1 = Frrvp_GP1_CgP1
        self._Frrvp_GP1_CgP1.flags.writeable = False
        self._Blrvp_GP1_CgP1 = Blrvp_GP1_CgP1
        self._Blrvp_GP1_CgP1.flags.writeable = False
        self._Brrvp_GP1_CgP1 = Brrvp_GP1_CgP1
        self._Brrvp_GP1_CgP1.flags.writeable = False

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._Crvp_GP1_CgP1: np.ndarray | None = None
        self._front_leg: _line_vortex.LineVortex | None = None
        self._left_leg: _line_vortex.LineVortex | None = None
        self._back_leg: _line_vortex.LineVortex | None = None
        self._right_leg: _line_vortex.LineVortex | None = None
        self._area: float | None = None

        # Initialize the mutable attributes to hold the HorseshoeVortex's strength and
        # the age of the RingVortex in seconds (in simulation time).
        self._strength = strength
        self.age = 0.0

    # --- Immutable: read only properties ---
    @property
    def Frrvp_GP1_CgP1(self) -> np.ndarray:
        return self._Frrvp_GP1_CgP1

    @property
    def Flrvp_GP1_CgP1(self) -> np.ndarray:
        return self._Flrvp_GP1_CgP1

    @property
    def Blrvp_GP1_CgP1(self) -> np.ndarray:
        return self._Blrvp_GP1_CgP1

    @property
    def Brrvp_GP1_CgP1(self) -> np.ndarray:
        return self._Brrvp_GP1_CgP1

    # --- Immutable derived: manual lazy caching ---
    @property
    def Crvp_GP1_CgP1(self) -> np.ndarray:
        """The position of the RingVortex's centroid (in the first Airplane's geometry
        axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the RingVortex's
            centroid (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters.
        """
        if self._Crvp_GP1_CgP1 is None:
            self._Crvp_GP1_CgP1 = _functions.numba_centroid_of_quadrilateral(
                self._Flrvp_GP1_CgP1,
                self._Frrvp_GP1_CgP1,
                self._Blrvp_GP1_CgP1,
                self._Brrvp_GP1_CgP1,
            )
            self._Crvp_GP1_CgP1.flags.writeable = False
        return self._Crvp_GP1_CgP1

    @property
    def front_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this RingVortex's front leg.

        :return: A LineVortex representing this RingVortex's front leg. The front leg
            goes from the front right point to the front left point.
        """
        if self._front_leg is None:
            self._front_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self._Frrvp_GP1_CgP1,
                Elvp_GP1_CgP1=self._Flrvp_GP1_CgP1,
                strength=self._strength,
            )
        return self._front_leg

    @property
    def left_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this RingVortex's left leg.

        :return: A LineVortex representing this RingVortex's left leg. The left leg goes
            from the front left point to the back left point.
        """
        if self._left_leg is None:
            self._left_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self._Flrvp_GP1_CgP1,
                Elvp_GP1_CgP1=self._Blrvp_GP1_CgP1,
                strength=self._strength,
            )
        return self._left_leg

    @property
    def back_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this RingVortex's back leg.

        :return: A LineVortex representing this RingVortex's back leg. The back leg goes
            from the back left point to the back right point.
        """
        if self._back_leg is None:
            self._back_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self._Blrvp_GP1_CgP1,
                Elvp_GP1_CgP1=self._Brrvp_GP1_CgP1,
                strength=self._strength,
            )
        return self._back_leg

    @property
    def right_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this RingVortex's right leg.

        :return: A LineVortex representing this RingVortex's right leg. The right leg
            goes from the back right point to the front right point.
        """
        if self._right_leg is None:
            self._right_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self._Brrvp_GP1_CgP1,
                Elvp_GP1_CgP1=self._Frrvp_GP1_CgP1,
                strength=self._strength,
            )
        return self._right_leg

    @property
    def area(self) -> float:
        """An estimate of this RingVortex's area.

        This is only an estimate because the surface defined by four line segments in
        3-space is a hyperboloid, and there doesn't seem to be a closed-form equation
        for the surface area of a hyperboloid between four points. Instead, we estimate
        the area using the cross product of RingVortex's diagonal vectors, which should
        be relatively accurate if the RingVortex can be approximated as a planar, convex
        quadrilateral.

        :return: An estimate of the RingVortex's area. The units are square meters.
        """
        if self._area is None:
            firstDiagonal_GP1 = self._Frrvp_GP1_CgP1 - self._Blrvp_GP1_CgP1
            secondDiagonal_GP1 = self._Flrvp_GP1_CgP1 - self._Brrvp_GP1_CgP1

            self._area = float(
                np.linalg.norm(np.cross(firstDiagonal_GP1, secondDiagonal_GP1)) / 2.0
            )
        return self._area

    # --- Mutable: property with propagation to children ---
    @property
    def strength(self) -> float:
        return self._strength

    @strength.setter
    def strength(self, new_strength: float) -> None:
        self._strength = new_strength

        if self._front_leg is not None:
            self._front_leg.strength = self._strength
        if self._left_leg is not None:
            self._left_leg.strength = self._strength
        if self._back_leg is not None:
            self._back_leg.strength = self._strength
        if self._right_leg is not None:
            self._right_leg.strength = self._strength
