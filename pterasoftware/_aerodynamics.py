"""Contains the vortex classes and functions for calculating induced velocities."""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from . import _functions

# Squire's parameter relates to the size of the vortex cores and the rate at which they
# grow. The value of this parameter is slightly controversial. It dramatically affects
# the stability of the result. I'm using this value, as cited for use in flapping-wing
# vehicles in "Role of Filament Strain in the Free-Vortex Modeling of Rotor Wakes"
# (Ananthan and Leishman, 2004). It is unitless.
_squire = 10**-4

# Lamb's constant relates to the size of the vortex cores and the rate at which they
# grow. The value of this parameter is well agreed upon, and published in "Extended
# Unsteady Vortex-Lattice Method for Insect Flapping Wings" (Nguyen et al., 2016). It is
# unitless.
_lamb = 1.25643

# The local machine error is used to fix removable discontinuities in the induced
# velocity functions.
_eps = np.finfo(float).eps


class RingVortex:
    """A class used to contain ring vortices.

    **Contains the following methods:**

    area: An estimate of this RingVortex's area.

    update_strength: Updates the strength of this RingVortex and of its four LineVortex
    legs.

    update_position: Updates the position of the RingVortex and of its attributes.
    """

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
        self.Flrvp_GP1_CgP1 = Flrvp_GP1_CgP1
        self.Frrvp_GP1_CgP1 = Frrvp_GP1_CgP1
        self.Blrvp_GP1_CgP1 = Blrvp_GP1_CgP1
        self.Brrvp_GP1_CgP1 = Brrvp_GP1_CgP1
        self.strength = strength

        # Initialize the LineVortices that make up the RingVortex.
        self.front_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Frrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Flrvp_GP1_CgP1,
            strength=self.strength,
        )
        self.left_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Flrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Blrvp_GP1_CgP1,
            strength=self.strength,
        )
        self.back_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Blrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Brrvp_GP1_CgP1,
            strength=self.strength,
        )
        self.right_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Brrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Frrvp_GP1_CgP1,
            strength=self.strength,
        )

        # Initialize a variable to hold the centroid of the RingVortex.
        self.Crvp_GP1_CgP1 = _functions.numba_centroid_of_quadrilateral(
            self.Flrvp_GP1_CgP1,
            self.Frrvp_GP1_CgP1,
            self.Blrvp_GP1_CgP1,
            self.Brrvp_GP1_CgP1,
        )

        # Initialize a variable to hold the age of the RingVortex in seconds (in
        # simulation time).
        self.age: float = 0.0

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
        firstDiagonal_GP1 = self.Frrvp_GP1_CgP1 - self.Blrvp_GP1_CgP1
        secondDiagonal_GP1 = self.Flrvp_GP1_CgP1 - self.Brrvp_GP1_CgP1

        return float(
            np.linalg.norm(np.cross(firstDiagonal_GP1, secondDiagonal_GP1)) / 2.0
        )

    def update_strength(self, strength: float) -> None:
        """Updates the strength of this RingVortex and of its four LineVortex legs.

        :param strength: The new strength of the RingVortex. Its units are in meters
            squared per second.
        :return: None
        """
        self.strength = strength
        self.right_leg.strength = strength
        self.front_leg.strength = strength
        self.left_leg.strength = strength
        self.back_leg.strength = strength

    def update_position(
        self,
        Frrvp_GP1_CgP1: np.ndarray,
        Flrvp_GP1_CgP1: np.ndarray,
        Blrvp_GP1_CgP1: np.ndarray,
        Brrvp_GP1_CgP1: np.ndarray,
    ) -> None:
        """Updates the position of the RingVortex and of its attributes.

        :param Frrvp_GP1_CgP1: A (3,) ndarray of floats representing the new position of
            the RingVortex's front right point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG).
        :param Flrvp_GP1_CgP1: A (3,) ndarray of floats representing the new position of
            the RingVortex's front left point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG).
        :param Blrvp_GP1_CgP1: A (3,) ndarray of floats representing the new position of
            the RingVortex's back left point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG).
        :param Brrvp_GP1_CgP1: A (3,) ndarray of floats representing the new position of
            the RingVortex's back right point (in the first Airplane's geometry axes,
            relative to the first Airplane's CG).
        :return: None
        """
        self.Flrvp_GP1_CgP1 = Flrvp_GP1_CgP1
        self.Frrvp_GP1_CgP1 = Frrvp_GP1_CgP1
        self.Blrvp_GP1_CgP1 = Blrvp_GP1_CgP1
        self.Brrvp_GP1_CgP1 = Brrvp_GP1_CgP1

        # TODO: Is it really faster to make new LineVortices instead of writing an
        #  update_position method inside of LineVortex and calling that from this
        #  method? Initialize the line vortices that make up the ring vortex. Or
        #  maybe I should make decorate the LineVortices as properties of RingVortex?
        self.front_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Frrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Flrvp_GP1_CgP1,
            strength=self.strength,
        )
        self.left_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Flrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Blrvp_GP1_CgP1,
            strength=self.strength,
        )
        self.back_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Blrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Brrvp_GP1_CgP1,
            strength=self.strength,
        )
        self.right_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Brrvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Frrvp_GP1_CgP1,
            strength=self.strength,
        )

        # Update the variable with the centroid of the RingVortex (in the first
        # Airplane's geometry axes, relative to the first Airplane's CG).
        self.Crvp_GP1_CgP1 = _functions.numba_centroid_of_quadrilateral(
            self.Flrvp_GP1_CgP1,
            self.Frrvp_GP1_CgP1,
            self.Blrvp_GP1_CgP1,
            self.Brrvp_GP1_CgP1,
        )


class HorseshoeVortex:
    """A class used to contain horseshoe vortices.

    **Contains the following methods:**

    update_strength: Updates the strength of this HorseshoeVortex and of its three
    LineVortex legs.
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
        self.Frhvp_GP1_CgP1 = Frhvp_GP1_CgP1
        self.Flhvp_GP1_CgP1 = Flhvp_GP1_CgP1
        self.leftLegVector_GP1 = leftLegVector_GP1 / np.linalg.norm(leftLegVector_GP1)
        self.left_right_leg_lengths = left_right_leg_lengths
        self.strength = strength

        # TODO: Consider making these properties.
        # Save attributes for the back right and back left horseshoe vortex points (in
        # the first Airplane's geometry axes, relative to the first Airplane's CG).
        self.Brhvp_GP1_CgP1 = (
            self.Frhvp_GP1_CgP1 + self.leftLegVector_GP1 * self.left_right_leg_lengths
        )
        self.Blhvp_GP1_CgP1 = (
            self.Flhvp_GP1_CgP1 + self.leftLegVector_GP1 * self.left_right_leg_lengths
        )

        # TODO: Consider making these properties.
        # Initialize LineVortices to represent the HorseshoeVortex's legs.
        self.right_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Brhvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Frhvp_GP1_CgP1,
            strength=self.strength,
        )
        self.finite_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Frhvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Flhvp_GP1_CgP1,
            strength=self.strength,
        )
        self.left_leg = _LineVortex(
            Slvp_GP1_CgP1=self.Flhvp_GP1_CgP1,
            Elvp_GP1_CgP1=self.Blhvp_GP1_CgP1,
            strength=self.strength,
        )

    # TODO: If we make the LineVortices strengths properties, then we can get rid of
    #  this method.
    def update_strength(self, strength: float) -> None:
        """Updates the strength of this HorseshoeVortex and of its three LineVortex
        legs.

        :param strength: The new strength of the HorseshoeVortex. Its units are in
            meters squared per second.
        :return: None
        """
        self.strength = strength
        self.right_leg.strength = strength
        self.finite_leg.strength = strength
        self.left_leg.strength = strength


class _LineVortex:
    """A class used to contain line vortices.

    **Contains the following methods:**

    None
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
        self.Slvp_GP1_CgP1 = Slvp_GP1_CgP1
        self.Elvp_GP1_CgP1 = Elvp_GP1_CgP1
        self.strength = strength

        # TODO: Consider making these properties.
        # Initialize variables to hold the vector from the LineVortex's start to end
        # point (in the first Airplane's geometry axes), and its center point (in the
        # first Airplane's geometry axes, relative to the first Airplane's CG).
        self.vector_GP1 = self.Elvp_GP1_CgP1 - self.Slvp_GP1_CgP1
        self.Clvp_GP1_CgP1 = self.Slvp_GP1_CgP1 + 0.5 * self.vector_GP1


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_ring_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrrvp_GP1_CgP1: np.ndarray,
    stackFrrvp_GP1_CgP1: np.ndarray,
    stackFlrvp_GP1_CgP1: np.ndarray,
    stackBlrvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of RingVortices and
    finds the cumulative induced velocity at every point.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M RingVortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' front right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' front left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackBlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' back left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        RingVortices. The units are in meters squared per second.
    :param ages: For bound RingVortices, this must be None. For RingVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M RingVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrrvp_GP1_CgP1,
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
        stackBrrvp_GP1_CgP1,
    ]

    stackVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], 3))

    # Get the velocity induced by each leg of each RingVortex (in the first Airplane's
    # geometry axes, observed from the Earth frame).
    for i in range(4):
        stackVInd_GP1__E += _collapsed_velocities_from_line_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
            stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return stackVInd_GP1__E


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_ring_vortices_chordwise_segments(
    stackP_GP1_CgP1: np.ndarray,
    stackBrrvp_GP1_CgP1: np.ndarray,
    stackFrrvp_GP1_CgP1: np.ndarray,
    stackFlrvp_GP1_CgP1: np.ndarray,
    stackBlrvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of RingVortices and
    finds the cumulative induced velocity at every point due to the RingVortices' left
    and right LineVortex legs.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M RingVortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' front right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' front left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackBlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' back left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        RingVortices. The units are in meters squared per second.
    :param ages: For bound RingVortices, this must be None. For RingVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M RingVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame) due to the M RingVortices' left and right LineVortex legs. The units are
        in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
    ]

    stackVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], 3))

    # Get the velocity induced by the left and right legs of each RingVortex (in the
    # first Airplane's geometry axes, observed from the Earth frame).
    for i in range(2):
        stackVInd_GP1__E += _collapsed_velocities_from_line_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
            stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return stackVInd_GP1__E


@njit(cache=True, fastmath=False)
def expanded_velocities_from_ring_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrrvp_GP1_CgP1: np.ndarray,
    stackFrrvp_GP1_CgP1: np.ndarray,
    stackFlrvp_GP1_CgP1: np.ndarray,
    stackBlrvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of RingVortices and
    finds the induced velocity at every point due to each RingVortex.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M RingVortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' front right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' front left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackBlrvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M RingVortices' back left vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        RingVortices. The units are in meters squared per second.
    :param ages: For bound RingVortices, this must be None. For RingVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M RingVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M RingVortices. The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrrvp_GP1_CgP1,
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrrvp_GP1_CgP1,
        stackFlrvp_GP1_CgP1,
        stackBlrvp_GP1_CgP1,
        stackBrrvp_GP1_CgP1,
    ]

    gridVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], strengths.shape[0], 3))

    # Get the velocity induced by each leg of each RingVortex (in the first Airplane's
    # geometry axes, observed from the Earth frame).
    for i in range(4):
        gridVInd_GP1__E += _expanded_velocities_from_line_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
            stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return gridVInd_GP1__E


# TODO: Remove the ability to specify HorseshoeVortex ages they are never used in
#  unsteady simulations.
@njit(cache=True, fastmath=False)
def collapsed_velocities_from_horseshoe_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrhvp_GP1_CgP1: np.ndarray,
    stackFrhvp_GP1_CgP1: np.ndarray,
    stackFlhvp_GP1_CgP1: np.ndarray,
    stackBlhvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of HorseshoeVortices and
    finds the cumulative induced velocity at every point.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M HorseshoeVortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M HorseshoeVortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M HorseshoeVortices' front left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackBlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M HorseshoeVortices' back left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        HorseshoeVortices. The units are in meters squared per second.
    :param ages: For bound HorseshoeVortices, this must be None. For HorseshoeVortices
        that have been shed into the wake, it must be a (M,) ndarray of floats
        representing the ages of the M HorseshoeVortices in seconds. The default is
        None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrhvp_GP1_CgP1,
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
        stackBlhvp_GP1_CgP1,
    ]

    stackVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], 3))

    # Get the velocity induced by each leg of each HorseshoeVortex (in the first
    # Airplane's geometry axes, observed from the Earth frame).
    for i in range(3):
        stackVInd_GP1__E += _collapsed_velocities_from_line_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
            stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return stackVInd_GP1__E


# TODO: Remove the ability to specify HorseshoeVortex ages they are never used in
#  unsteady simulations.
@njit(cache=True, fastmath=False)
def expanded_velocities_from_horseshoe_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackBrhvp_GP1_CgP1: np.ndarray,
    stackFrhvp_GP1_CgP1: np.ndarray,
    stackFlhvp_GP1_CgP1: np.ndarray,
    stackBlhvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of HorseshoeVortices and
    finds the induced velocity at every point due to each HorseshoeVortex.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackBrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        M HorseshoeVortices' back right vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackFrhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M HorseshoeVortices' front right vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackFlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M HorseshoeVortices' front left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param stackBlhvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M HorseshoeVortices' back left vertices (in the first Airplane's geometry
        axes, relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of M
        HorseshoeVortices. The units are in meters squared per second.
    :param ages: For bound HorseshoeVortices, this must be None. For HorseshoeVortices
        that have been shed into the wake, it must be a (M,) ndarray of floats
        representing the ages of the M HorseshoeVortices in seconds. The default is
        None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M HorseshoeVortices. The units are in meters per second.
    """
    listStackSlvp_GP1_CgP1 = [
        stackBrhvp_GP1_CgP1,
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
    ]
    listStackElvp_GP1_CgP1 = [
        stackFrhvp_GP1_CgP1,
        stackFlhvp_GP1_CgP1,
        stackBlhvp_GP1_CgP1,
    ]

    gridVInd_GP1__E = np.zeros((stackP_GP1_CgP1.shape[0], strengths.shape[0], 3))

    # Get the velocity induced by each leg of each HorseshoeVortex (in the first
    # Airplane's geometry axes, observed from the Earth frame).
    for i in range(3):
        gridVInd_GP1__E += _expanded_velocities_from_line_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackSlvp_GP1_CgP1=listStackSlvp_GP1_CgP1[i],
            stackElvp_GP1_CgP1=listStackElvp_GP1_CgP1[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return gridVInd_GP1__E


@njit(cache=True, fastmath=False)
def _collapsed_velocities_from_line_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackSlvp_GP1_CgP1: np.ndarray,
    stackElvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of LineVortices and
    finds the cumulative induced velocity at every point.

    This function uses a modified version of the Bio-Savart law to create a smooth
    induced velocity decay based on a LineVortex's core radius. The radius is determined
    based on the LineVortex's age and the kinematic viscosity. If the age of the
    LineVortex is 0.0 seconds, the radius is set to 0.0 meters. The age of a LineVortex
    in only relevant for vortices that have been shed into the wake.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    **Citation:**

    Equation adapted from: "Extended Unsteady Vortex-Lattice Method for Insect Flapping
    Wings"

    Authors: Anh Tuan Nguyen, Joong-Kwan Kim, Jong-Seob Han, and Jae-Hung Han

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackSlvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M LineVortices' starting vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackElvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M LineVortices' ending vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        LineVortices. The units are in meters squared per second.
    :param ages: For bound LineVortices, this must be None. For LineVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M LineVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    num_vortices = stackSlvp_GP1_CgP1.shape[0]
    num_points = stackP_GP1_CgP1.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in
    # the first Airplane's geometry axes, observed from the Earth frame).
    stackVInd_GP1__E = np.zeros((num_points, 3))

    # If the user didn't specify any ages, set the age of each LineVortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
        Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]
        strength = strengths[vortex_id]
        age = ages[vortex_id]

        # Calculate the radius of the LineVortex's core. If the age is 0.0 seconds,
        # this will evaluate to be 0.0 meters.
        r_c = 2 * math.sqrt(_lamb * (nu + _squire * abs(strength)) * age)

        # The r0_GP1 vector goes from the LineVortex's start point to its end point (in
        # the first Airplane's geometry axes).
        r0X_GP1 = Elvp_GP1_CgP1[0] - Slvp_GP1_CgP1[0]
        r0Y_GP1 = Elvp_GP1_CgP1[1] - Slvp_GP1_CgP1[1]
        r0Z_GP1 = Elvp_GP1_CgP1[2] - Slvp_GP1_CgP1[2]

        # Find r0_GP1's length.
        r0 = math.sqrt(r0X_GP1**2 + r0Y_GP1**2 + r0Z_GP1**2)

        c_1 = strength / (4 * math.pi)
        c_2 = r0**2 * r_c**2

        for point_id in range(num_points):
            P_GP1_CgP1 = stackP_GP1_CgP1[point_id]

            # The r1_GP1 vector goes from P_GP1_CgP1 to the LineVortex's start point (in
            # the first Airplane's geometry axes).
            r1X_GP1 = Slvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r1Y_GP1 = Slvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r1Z_GP1 = Slvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r2_GP1 vector goes from P_GP1_CgP1 to the LineVortex's end point (in
            # the first Airplane's geometry axes).
            r2X_GP1 = Elvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r2Y_GP1 = Elvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r2Z_GP1 = Elvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r3_GP1 vector is the cross product of r1_GP1 and r2_GP1 (in the first
            # Airplane's geometry axes).
            r3X_GP1 = r1Y_GP1 * r2Z_GP1 - r1Z_GP1 * r2Y_GP1
            r3Y_GP1 = r1Z_GP1 * r2X_GP1 - r1X_GP1 * r2Z_GP1
            r3Z_GP1 = r1X_GP1 * r2Y_GP1 - r1Y_GP1 * r2X_GP1

            # Find the lengths of r1_GP1, r2_GP1, and r3_GP1.
            r1 = math.sqrt(r1X_GP1**2 + r1Y_GP1**2 + r1Z_GP1**2)
            r2 = math.sqrt(r2X_GP1**2 + r2Y_GP1**2 + r2Z_GP1**2)
            r3 = math.sqrt(r3X_GP1**2 + r3Y_GP1**2 + r3Z_GP1**2)

            c_3 = r1X_GP1 * r2X_GP1 + r1Y_GP1 * r2Y_GP1 + r1Z_GP1 * r2Z_GP1

            # If part of the LineVortex is so close to P_GP1_CgP1 that they are touching
            # (within machine epsilon), there is a removable discontinuity. In this
            # case, continue to the next point because there is no velocity induced by
            # the current LineVortex at this point.
            if r1 < _eps or r2 < _eps or r3**2 < _eps:
                continue
            else:
                c_4 = c_1 * (r1 + r2) * (r1 * r2 - c_3) / (r1 * r2 * (r3**2 + c_2))
                stackVInd_GP1__E[point_id, 0] += c_4 * r3X_GP1
                stackVInd_GP1__E[point_id, 1] += c_4 * r3Y_GP1
                stackVInd_GP1__E[point_id, 2] += c_4 * r3Z_GP1
    return stackVInd_GP1__E


@njit(cache=True, fastmath=False)
def _expanded_velocities_from_line_vortices(
    stackP_GP1_CgP1: np.ndarray,
    stackSlvp_GP1_CgP1: np.ndarray,
    stackElvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of LineVortices and
    finds the induced velocity at every point due to each LineVortex.

    This function uses a modified version of the Bio-Savart law to create a smooth
    induced velocity decay based on a LineVortex's core radius. The radius is determined
    based on a LineVortex's age and the kinematic viscosity. If the age of the
    LineVortex is 0.0 seconds, the radius is set to 0.0 meters. The age of a vortex in
    only relevant for LineVortices that have been shed into the wake.

    This function's performance has been highly optimized for unsteady simulations via
    Numba. While using Numba dramatically increases unsteady simulation performance, it
    does cause a performance drop for the less intense steady simulations.

    **Citation:**

    Equation adapted from: "Extended Unsteady Vortex-Lattice Method for Insect Flapping
    Wings"

    Authors: Anh Tuan Nguyen, Joong-Kwan Kim, Jong-Seob Han, and Jae-Hung Han

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackSlvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of M
        LineVortices' starting vertices (in the first Airplane's geometry axes, relative
        to the first Airplane's CG). The units are in meters.
    :param stackElvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M LineVortices' ending vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        LineVortices. The units are in meters squared per second.
    :param ages: For bound LineVortices, this must be None. For LineVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M LineVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M LineVortices. The units are in meters per second.
    """
    num_vortices = stackSlvp_GP1_CgP1.shape[0]
    num_points = stackP_GP1_CgP1.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in the
    # first Airplane's geometry axes, observed from the Earth frame).
    gridVInd_GP1__E = np.zeros((num_points, num_vortices, 3))

    # If the user didn't specify any ages, set the age of each LineVortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
        Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]
        strength = strengths[vortex_id]
        age = ages[vortex_id]

        # Calculate the radius of the LineVortex's core. If the age is 0.0 seconds,
        # this will evaluate to be 0.0 meters.
        r_c = 2 * math.sqrt(_lamb * (nu + _squire * abs(strength)) * age)

        # The r0_GP1 vector goes from the LineVortex's start point to its end point (in
        # the first Airplane's geometry axes).
        r0X_GP1 = Elvp_GP1_CgP1[0] - Slvp_GP1_CgP1[0]
        r0Y_GP1 = Elvp_GP1_CgP1[1] - Slvp_GP1_CgP1[1]
        r0Z_GP1 = Elvp_GP1_CgP1[2] - Slvp_GP1_CgP1[2]

        # Find r0_GP1's length.
        r0 = math.sqrt(r0X_GP1**2 + r0Y_GP1**2 + r0Z_GP1**2)

        c_1 = strength / (4 * math.pi)
        c_2 = r0**2 * r_c**2

        for point_id in range(num_points):
            P_GP1_CgP1 = stackP_GP1_CgP1[point_id]

            # The r1_GP1 vector goes from P_GP1_CgP1 to the LineVortex's start point (in
            # the first Airplane's geometry axes).
            r1X_GP1 = Slvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r1Y_GP1 = Slvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r1Z_GP1 = Slvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r2_GP1 vector goes from P_GP1_CgP1 to the LineVortex's end point (in
            # the first Airplane's geometry axes).
            r2X_GP1 = Elvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r2Y_GP1 = Elvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r2Z_GP1 = Elvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r3_GP1 vector is the cross product of r1_GP1 and r2_GP1 (in the first
            # Airplane's geometry axes).
            r3X_GP1 = r1Y_GP1 * r2Z_GP1 - r1Z_GP1 * r2Y_GP1
            r3Y_GP1 = r1Z_GP1 * r2X_GP1 - r1X_GP1 * r2Z_GP1
            r3Z_GP1 = r1X_GP1 * r2Y_GP1 - r1Y_GP1 * r2X_GP1

            # Find the lengths of r1_GP1, r2_GP1, and r3_GP1.
            r1 = math.sqrt(r1X_GP1**2 + r1Y_GP1**2 + r1Z_GP1**2)
            r2 = math.sqrt(r2X_GP1**2 + r2Y_GP1**2 + r2Z_GP1**2)
            r3 = math.sqrt(r3X_GP1**2 + r3Y_GP1**2 + r3Z_GP1**2)

            c_3 = r1X_GP1 * r2X_GP1 + r1Y_GP1 * r2Y_GP1 + r1Z_GP1 * r2Z_GP1

            # If part of the LineVortex is so close to P_GP1_CgP1 that they are touching
            # (within machine epsilon), there is a removable discontinuity. In this
            # case, continue to the next point because there is no velocity induced by
            # the current LineVortex at this point.
            if r1 < _eps or r2 < _eps or r3**2 < _eps:
                continue
            else:
                c_4 = c_1 * (r1 + r2) * (r1 * r2 - c_3) / (r1 * r2 * (r3**2 + c_2))
                gridVInd_GP1__E[point_id, vortex_id, 0] = c_4 * r3X_GP1
                gridVInd_GP1__E[point_id, vortex_id, 1] = c_4 * r3Y_GP1
                gridVInd_GP1__E[point_id, vortex_id, 2] = c_4 * r3Z_GP1
    return gridVInd_GP1__E
