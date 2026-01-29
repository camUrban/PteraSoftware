"""Contains the HorseshoeVortex class."""

from __future__ import annotations

import numpy as np

from . import _line_vortex


class HorseshoeVortex:
    """A class used to contain horseshoe vortices.

    **Contains the following methods:**

    Brhvp_GP1_CgP1: The position of the HorseshoeVortex's back right point (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    Blhvp_GP1_CgP1: The position of the HorseshoeVortex's back left point (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    right_leg: The LineVortex representing this HorseshoeVortex's right leg.

    finite_leg: The LineVortex representing this HorseshoeVortex's finite leg.

    left_leg: The LineVortex representing this HorseshoeVortex's left leg.

    **Notes:**

    Computed geometric properties (back right point, back left point, and LineVortex
    legs) are lazily evaluated and cached. Setting front point positions, the left leg
    vector, or the left/right leg lengths invalidates dependent cached values, ensuring
    consistency while avoiding redundant computation.
    """

    def __init__(
        self,
        Frhvp_GP1_CgP1: np.ndarray,
        Flhvp_GP1_CgP1: np.ndarray,
        leftLegVector_GP1: np.ndarray,
        left_right_leg_lengths: float,
        strength: float,
    ) -> None:
        """The initialization method.

        :param Frhvp_GP1_CgP1: A (3,) ndarray of floats representing the position of the
            HorseshoeVortex's front right point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The front right point is defined as
            the start point of the HorseshoeVortex's front leg, which is also its one
            finite leg. The units are in meters.
        :param Flhvp_GP1_CgP1: A (3,) ndarray of floats representing the position of the
            HorseshoeVortex's front left point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The front left point is defined as the
            end point of the HorseshoeVortex's front leg, which is also its one finite
            leg. The units are in meters.
        :param leftLegVector_GP1: A (3,) ndarray of floats representing the direction
            vector of the HorseshoeVortex's left leg (in the first Airplane's geometry
            axes). The left leg starts from the front left point and ends at the back
            left point. It is one of the HorseshoeVortex's two quasi infinite legs, the
            other being the right leg. The right leg's vector (in the first Airplane's
            geometry axes) is defined as -1.0 times this vector. It will be normalized
            to a unit vector during initialization. The units are in meters.
        :param left_right_leg_lengths: The length of the HorseshoeVortex's left and
            right quasi infinite legs. I recommend setting it to at least 20 times the
            length of the finite leg. The units are in meters.
        :param strength: The strength of the HorseshoeVortex. Its units are in meters
            squared per second.
        :return: None
        """
        # Declare type annotations and initialize the private cache variables.
        self._Frhvp_GP1_CgP1: np.ndarray
        self._Flhvp_GP1_CgP1: np.ndarray
        self._leftLegVector_GP1: np.ndarray
        self._left_right_leg_lengths: float
        self._strength: float
        self._Brhvp_GP1_CgP1: np.ndarray | None = None
        self._Blhvp_GP1_CgP1: np.ndarray | None = None
        self._right_leg: _line_vortex.LineVortex | None = None
        self._finite_leg: _line_vortex.LineVortex | None = None
        self._left_leg: _line_vortex.LineVortex | None = None

        # Initialize the attributes.
        self.Frhvp_GP1_CgP1 = Frhvp_GP1_CgP1
        self.Flhvp_GP1_CgP1 = Flhvp_GP1_CgP1
        self.leftLegVector_GP1 = leftLegVector_GP1
        self.left_right_leg_lengths = left_right_leg_lengths
        self.strength = strength

    @property
    def Frhvp_GP1_CgP1(self) -> np.ndarray:
        return self._Frhvp_GP1_CgP1

    @Frhvp_GP1_CgP1.setter
    def Frhvp_GP1_CgP1(self, newFrhvp_GP1_CgP1: np.ndarray) -> None:
        self._Brhvp_GP1_CgP1 = None

        self._Frhvp_GP1_CgP1 = newFrhvp_GP1_CgP1

        if self._right_leg is not None:
            self._right_leg.Elvp_GP1_CgP1 = self._Frhvp_GP1_CgP1
        if self._finite_leg is not None:
            self._finite_leg.Slvp_GP1_CgP1 = self._Frhvp_GP1_CgP1

    @property
    def Flhvp_GP1_CgP1(self) -> np.ndarray:
        return self._Flhvp_GP1_CgP1

    @Flhvp_GP1_CgP1.setter
    def Flhvp_GP1_CgP1(self, newFlhvp_GP1_CgP1: np.ndarray) -> None:
        self._Blhvp_GP1_CgP1 = None

        self._Flhvp_GP1_CgP1 = newFlhvp_GP1_CgP1

        if self._finite_leg is not None:
            self._finite_leg.Elvp_GP1_CgP1 = self._Flhvp_GP1_CgP1
        if self._left_leg is not None:
            self._left_leg.Slvp_GP1_CgP1 = self._Flhvp_GP1_CgP1

    @property
    def leftLegVector_GP1(self) -> np.ndarray:
        return self._leftLegVector_GP1

    @leftLegVector_GP1.setter
    def leftLegVector_GP1(self, newLeftLegVector_GP1: np.ndarray) -> None:
        self._Brhvp_GP1_CgP1 = None
        self._Blhvp_GP1_CgP1 = None

        self._leftLegVector_GP1 = newLeftLegVector_GP1 / np.linalg.norm(
            newLeftLegVector_GP1
        )

        if self._right_leg is not None:
            self._right_leg.Slvp_GP1_CgP1 = self.Brhvp_GP1_CgP1
        if self._left_leg is not None:
            self._left_leg.Elvp_GP1_CgP1 = self.Blhvp_GP1_CgP1

    @property
    def left_right_leg_lengths(self) -> float:
        return self._left_right_leg_lengths

    @left_right_leg_lengths.setter
    def left_right_leg_lengths(self, new_left_right_leg_lengths: float) -> None:
        self._Brhvp_GP1_CgP1 = None
        self._Blhvp_GP1_CgP1 = None

        self._left_right_leg_lengths = new_left_right_leg_lengths

        if self._right_leg is not None:
            self._right_leg.Slvp_GP1_CgP1 = self.Brhvp_GP1_CgP1
        if self._left_leg is not None:
            self._left_leg.Elvp_GP1_CgP1 = self.Blhvp_GP1_CgP1

    @property
    def strength(self) -> float:
        return self._strength

    @strength.setter
    def strength(self, new_strength: float) -> None:
        self._strength = new_strength

        if self._right_leg is not None:
            self._right_leg.strength = self._strength
        if self._finite_leg is not None:
            self._finite_leg.strength = self._strength
        if self._left_leg is not None:
            self._left_leg.strength = self._strength

    @property
    def Brhvp_GP1_CgP1(self) -> np.ndarray:
        """The position of the HorseshoeVortex's back right point (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the
            HorseshoeVortex's back right point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The units are in meters.
        """
        if self._Brhvp_GP1_CgP1 is None:
            self._Brhvp_GP1_CgP1 = (
                self.Frhvp_GP1_CgP1
                + self.leftLegVector_GP1 * self.left_right_leg_lengths
            )
        return self._Brhvp_GP1_CgP1

    @property
    def Blhvp_GP1_CgP1(self) -> np.ndarray:
        """The position of the HorseshoeVortex's back left point (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the
            HorseshoeVortex's back left point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). The units are in meters.
        """
        if self._Blhvp_GP1_CgP1 is None:
            self._Blhvp_GP1_CgP1 = (
                self.Flhvp_GP1_CgP1
                + self.leftLegVector_GP1 * self.left_right_leg_lengths
            )
        return self._Blhvp_GP1_CgP1

    @property
    def right_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this HorseshoeVortex's right leg.

        :return: A LineVortex representing this HorseshoeVortex's right leg. The right
            leg goes from the back right point to the front right point.
        """
        if self._right_leg is None:
            self._right_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self.Brhvp_GP1_CgP1,
                Elvp_GP1_CgP1=self.Frhvp_GP1_CgP1,
                strength=self.strength,
            )
        return self._right_leg

    @property
    def finite_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this HorseshoeVortex's finite leg.

        :return: A LineVortex representing this HorseshoeVortex's finite leg. The finite
            leg goes from the front right point to the front left point.
        """
        if self._finite_leg is None:
            self._finite_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self.Frhvp_GP1_CgP1,
                Elvp_GP1_CgP1=self.Flhvp_GP1_CgP1,
                strength=self.strength,
            )
        return self._finite_leg

    @property
    def left_leg(self) -> _line_vortex.LineVortex:
        """The LineVortex representing this HorseshoeVortex's left leg.

        :return: A LineVortex representing this HorseshoeVortex's left leg. The left leg
            goes from the front left point to the back left point.
        """
        if self._left_leg is None:
            self._left_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=self.Flhvp_GP1_CgP1,
                Elvp_GP1_CgP1=self.Blhvp_GP1_CgP1,
                strength=self.strength,
            )
        return self._left_leg
