# TODO: Consider making this module private (renaming it with a _ prefix).
"""This module contains vortex class definitions, and functions for calculating
induced velocities.

This module contains the following classes:
    RingVortex: This class is used to contain ring vortices.

    HorseshoeVortex: This class is used to contain horseshoe vortices.

This module contains the following exceptions:
    None

This module contains the following functions:
    collapsed_velocities_from_ring_vortices: This function takes in a group of
    points, and the attributes of a group of RingVortices. At every point, it finds
    the cumulative induced velocity due to all the RingVortices.

    collapsed_velocities_from_ring_vortices_chordwise_segments: This function takes
    in a group of points, and the attributes of a group of RingVortices. At every
    point, it finds the cumulative induced velocity due to all the RingVortices'
    chordwise segments.

    expanded_velocities_from_ring_vortices: This function takes in a group of points,
    and the attributes of a group of RingVortices. At every point, it finds the
    induced velocity due to each RingVortex.

    collapsed_velocities_from_horseshoe_vortices: This function takes in a group of
    points, and the attributes of a group of HorseshoeVortices. At every point,
    it finds the cumulative induced velocity due to all the HorseshoeVortices.

    expanded_velocities_from_horseshoe_vortices: This function takes in a group of
    points, and the attributes of a group of HorseshoeVortices. At every point,
    it finds the induced velocity due to each HorseshoeVortex.
"""

import math

import numpy as np
from numba import njit

from . import functions
from . import parameter_validation

# Set the value of Squire's parameter that will be used by the induced velocity
# functions. Squire's parameter relates to the size of the vortex cores and the rate
# at which they grow. The value of this parameter is slightly controversial. It
# dramatically affects the stability of the result. I'm using this value, as cited
# for use in flapping-wing vehicles in "Role of Filament Strain in the Free-Vortex
# Modeling of Rotor Wakes" (Ananthan and Leishman, 2004). It is unitless.
_squire = 10**-4

# Set the value of Lamb's constant that will be used by the induced velocity
# functions. Lamb's constant relates to the size of the vortex cores and the rate at
# which they grow. The value of this parameter is well agreed upon, and published in
# "Extended Unsteady Vortex-Lattice Method for Insect Flapping Wings" (Nguyen et al.,
# 2016). It is unitless.
_lamb = 1.25643

# Set the value of the local machine error. This will be used to fix removable
# discontinuities in the induced velocity functions.
_eps = np.finfo(float).eps


class RingVortex:
    """This class is used to contain ring vortices.

    This class contains the following public methods:
        update_strength: This method updates the strength of this RingVortex,
        and the strength of its four legs' LineVortices.

        update_position: This method updates the position of the RingVortex, and the
        positions of all its attributes.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        Frrvp_G_Cg,
        Flrvp_G_Cg,
        Blrvp_G_Cg,
        Brrvp_G_Cg,
        strength,
    ):
        """This is the initialization method.

        :param Frrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's front-right point (in geometry
            axes, relative to the CG). The front-right point is defined as the end
            point of the RingVortex's right leg and the start point of its front leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param Flrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's front-left point (in geometry
            axes, relative to the CG). The front-left point is defined as the end
            point of the RingVortex's front leg and the start point of its left leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param Blrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's back-left point (in geometry axes,
            relative to the CG). The back-left point is defined as the end point of
            the RingVortex's left leg and the start point of its back leg. Can be a
            list, tuple, or numpy array of numbers (int or float). Values are
            converted to floats internally. The units are in meters.

        :param Brrvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the RingVortex's back-right point (in geometry
            axes, relative to the CG). The back-right point is defined as the end
            point of the RingVortex's back leg and the start point of its right leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param strength: number

            This is the strength of the RingVortex It must be a positive number and
            will be converted internally to a float. Its units are in meters squared
            per second.
        """
        self.Flrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Flrvp_G_Cg, "Flrvp_G_Cg"
        )
        self.Frrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Frrvp_G_Cg, "Frrvp_G_Cg"
        )
        self.Blrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Blrvp_G_Cg, "Blrvp_G_Cg"
        )
        self.Brrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Brrvp_G_Cg, "Brrvp_G_Cg"
        )
        if strength is not None:
            strength = parameter_validation.number_return_float(strength, "strength")
        self.strength = strength

        # Initialize the LineVortices that make up the RingVortex.
        self.front_leg = _LineVortex(
            Slvp_G_Cg=self.Frrvp_G_Cg,
            Elvp_G_Cg=self.Flrvp_G_Cg,
            strength=self.strength,
        )
        self.left_leg = _LineVortex(
            Slvp_G_Cg=self.Flrvp_G_Cg,
            Elvp_G_Cg=self.Blrvp_G_Cg,
            strength=self.strength,
        )
        self.back_leg = _LineVortex(
            Slvp_G_Cg=self.Blrvp_G_Cg,
            Elvp_G_Cg=self.Brrvp_G_Cg,
            strength=self.strength,
        )
        self.right_leg = _LineVortex(
            Slvp_G_Cg=self.Brrvp_G_Cg,
            Elvp_G_Cg=self.Frrvp_G_Cg,
            strength=self.strength,
        )

        # TODO: Consider making this a property.
        # Initialize a variable to hold the centroid of the RingVortex.
        self.Crvp_G_Cg = functions.numba_centroid_of_quadrilateral(
            self.Flrvp_G_Cg,
            self.Frrvp_G_Cg,
            self.Blrvp_G_Cg,
            self.Brrvp_G_Cg,
        )

        # Initialize a variable to hold the age of the RingVortex in seconds (in
        # simulation time).
        self.age = 0

    def update_strength(self, strength):
        """This method updates the strength of this RingVortex, and the strength
        of its four legs' LineVortices.

        :param strength: float

            This is the strength of this vortex, and of its four
            legs' LineVortices. Its units are meters squared per second.

        :return: None
        """
        self.strength = strength
        self.right_leg.strength = strength
        self.front_leg.strength = strength
        self.left_leg.strength = strength
        self.back_leg.strength = strength

    def update_position(self, Frrvp_G_Cg, Flrvp_G_Cg, Blrvp_G_Cg, Brrvp_G_Cg):
        """This method updates the position of the RingVortex, and the positions of all
        its attributes.

        :param Frrvp_G_Cg: array-like of 3 numbers

            New position [x, y, z] of the RingVortex's front-right point (in geometry
            axes, relative to the CG). The front-right point is defined as the end
            point of the RingVortex's right leg and the start point of its front leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param Flrvp_G_Cg: array-like of 3 numbers

            New position [x, y, z] of the RingVortex's front-left point (in geometry
            axes, relative to the CG). The front-left point is defined as the end
            point of the RingVortex's front leg and the start point of its left leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :param Blrvp_G_Cg: array-like of 3 numbers

            New position [x, y, z] of the RingVortex's back-left point (in geometry axes,
            relative to the CG). The back-left point is defined as the end point of
            the RingVortex's left leg and the start point of its back leg. Can be a
            list, tuple, or numpy array of numbers (int or float). Values are
            converted to floats internally. The units are in meters.

        :param Brrvp_G_Cg: array-like of 3 numbers

            New position [x, y, z] of the RingVortex's back-right point (in geometry
            axes, relative to the CG). The back-right point is defined as the end
            point of the RingVortex's back leg and the start point of its right leg.
            Can be a list, tuple, or numpy array of numbers (int or float). Values
            are converted to floats internally. The units are in meters.

        :return: None
        """
        self.Flrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Flrvp_G_Cg, "Flrvp_G_Cg"
        )
        self.Frrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Frrvp_G_Cg, "Frrvp_G_Cg"
        )
        self.Blrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Blrvp_G_Cg, "Blrvp_G_Cg"
        )
        self.Brrvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Brrvp_G_Cg, "Brrvp_G_Cg"
        )

        # TODO: Is it really faster to make new LineVortices instead of writing an
        #  update_position method inside of LineVortex and calling that from this
        #  method? Initialize the line vortices that make up the ring vortex. Or
        #  maybe I should make decorate the LineVortices as properties of RingVortex?
        self.front_leg = _LineVortex(
            Slvp_G_Cg=self.Frrvp_G_Cg,
            Elvp_G_Cg=self.Flrvp_G_Cg,
            strength=self.strength,
        )
        self.left_leg = _LineVortex(
            Slvp_G_Cg=self.Flrvp_G_Cg,
            Elvp_G_Cg=self.Blrvp_G_Cg,
            strength=self.strength,
        )
        self.back_leg = _LineVortex(
            Slvp_G_Cg=self.Blrvp_G_Cg,
            Elvp_G_Cg=self.Brrvp_G_Cg,
            strength=self.strength,
        )
        self.right_leg = _LineVortex(
            Slvp_G_Cg=self.Brrvp_G_Cg,
            Elvp_G_Cg=self.Frrvp_G_Cg,
            strength=self.strength,
        )

        # Update the variable with the centroid of the RingVortex.
        self.Crvp_G_Cg = functions.numba_centroid_of_quadrilateral(
            self.Flrvp_G_Cg,
            self.Frrvp_G_Cg,
            self.Blrvp_G_Cg,
            self.Brrvp_G_Cg,
        )


class HorseshoeVortex:
    """This class is used to contain horseshoe vortices.

    :param Frhvp_G_Cg: array-like of 3 numbers

        Position [x, y, z] of the HorseshoeVortex's front-right point (in geometry
        axes, relative to the CG). The front-right point is defined as the start
        point of the HorseshoeVortex's front leg, which is also its one finite leg.
        Can be a list, tuple, or numpy array of numbers (int or float). Values are
        converted to floats internally. The units are in meters.

    :param Flhvp_G_Cg: array-like of 3 numbers

        Position [x, y, z] of the HorseshoeVortex's front-left point (in geometry
        axes, relative to the CG). The front-left point is defined as the end point
        of the HorseshoeVortex's front leg, which is also its one finite leg. Can be
        a list, tuple, or numpy array of numbers (int or float). Values are converted
        to floats internally. The units are in meters.

    :param leftLegVector_G: array-like of 3 numbers

        Direction vector of the HorseshoeVortex's left leg (in geometry axes). The
        left leg starts from the front-left point and ends at the back-left point. It
        is one of the HorseshoeVortex's two quasi-infinite legs, the other being the
        right leg. It can be a list, tuple, or numpy array of numbers (int or float).
        Values are converted to floats internally. If this isn't already a unit
        vector, it will be converted to one during initialization. The right leg's
        vector (in geometry axes) is defined as -1 times this vector. The units are
        in meters.

    :param left_right_leg_lengths: number

        This is the length of the HorseshoeVortex's left and right quasi-infinite
        legs. It must be a positive number and will be converted internally to a
        float. I recommend setting it to at least 20 times the length of the finite
        leg. The units are in meters.

    :param strength: number

        This is the strength of the HorseshoeVortex. It must be a number and will be
        converted internally to a float. Its units are in meters squared per second.
    """

    def __init__(
        self,
        Frhvp_G_Cg,
        Flhvp_G_Cg,
        leftLegVector_G,
        left_right_leg_lengths,
        strength,
    ):
        """This is the constructor method."""
        self.Frhvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Frhvp_G_Cg, "Frhvp_G_Cg"
        )
        self.Flhvp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Flhvp_G_Cg, "Flhvp_G_Cg"
        )
        self.leftLegVector_G = (
            parameter_validation.threeD_number_vectorLike_return_float_unit_vector(
                leftLegVector_G, "leftLegVector_G"
            )
        )
        self.left_right_leg_lengths = parameter_validation.positive_number_return_float(
            left_right_leg_lengths, "left_right_leg_lengths"
        )
        if strength is not None:
            strength = parameter_validation.number_return_float(strength, "strength")
        self.strength = strength

        # TODO: Consider making these properties.
        # Save attributes for the back right and back left horseshoe vortex points (in
        # geometry axes, relative to the CG).
        self.Brhvp_G_Cg = (
            self.Frhvp_G_Cg + self.leftLegVector_G * self.left_right_leg_lengths
        )
        self.Blhvp_G_Cg = (
            self.Flhvp_G_Cg + self.leftLegVector_G * self.left_right_leg_lengths
        )

        # TODO: Consider making these properties.
        # Initialize LineVortices to represent the HorseshoeVortex's legs.
        self.right_leg = _LineVortex(
            Slvp_G_Cg=self.Brhvp_G_Cg,
            Elvp_G_Cg=self.Frhvp_G_Cg,
            strength=self.strength,
        )
        self.finite_leg = _LineVortex(
            Slvp_G_Cg=self.Frhvp_G_Cg,
            Elvp_G_Cg=self.Flhvp_G_Cg,
            strength=self.strength,
        )
        self.left_leg = _LineVortex(
            Slvp_G_Cg=self.Flhvp_G_Cg,
            Elvp_G_Cg=self.Blhvp_G_Cg,
            strength=self.strength,
        )

    # TODO: If we make the LineVortices strengths, then we can get rid of this method.
    def update_strength(self, strength):
        """This method updates the strength of this HorseshoeVortex object, and the
        strength of its leg LineVortices.

        :param strength: float

            This is the strength of this vortex, and of its LineVortex legs. Its
            units are meters squared per second.
        """
        self.strength = strength
        self.right_leg.strength = strength
        self.finite_leg.strength = strength
        self.left_leg.strength = strength


class _LineVortex:
    """This class is used to contain line vortices.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, Slvp_G_Cg, Elvp_G_Cg, strength):
        """This is the initialization method.

        Note: This is a private class, so it doesn't validate its parameters.
        Validate them before passing them in.

        :param Slvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the LineVortex's start point (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. The units are in
            meters.

        :param Elvp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the LineVortex's end point (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. The units are in
            meters.

        :param strength: number

            This is the strength of the LineVortex. It must be a number and will be
            converted internally to a float. Its units are in meters squared per second.
        """
        self.Slvp_G_Cg = Slvp_G_Cg
        self.Elvp_G_Cg = Elvp_G_Cg
        self.strength = strength

        # TODO: Consider making these properties.
        # Initialize variables to hold the vector from the LineVortex's start to end
        # point (in geometry axes), and its center point (in geometry axes, relative
        # to the CG).
        self.vector_G = self.Elvp_G_Cg - self.Slvp_G_Cg
        self.Clvp_G_Cg = self.Slvp_G_Cg + 0.5 * self.vector_G


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_ring_vortices(
    stackP_G_Cg,
    stackBrrvp_G_Cg,
    stackFrrvp_G_Cg,
    stackFlrvp_G_Cg,
    stackBlrvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    RingVortices. At every point, it finds the cumulative induced velocity due to all
    the RingVortices.

    Note: See the docstring for collapsed_velocities_from_line_vortices for
    information on usage.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    # TODO: Update docstring parameters' descriptions and types.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackBrrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's back right vertex's position in meters.
    :param stackFrrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front right vertex's position in meters.
    :param stackFlrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param stackBlrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of ring
        vortices. Each holds the strength of that ring vortex in meters squared per
        second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that ring vortex in seconds. This
        is only relevant for vortices that have been shed into the wake. The default
        value is None. If the age of a specific vortex is 0.0 seconds, then the
        vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 2D array of floats
        This is an array of shape (N x 3), and it holds the cumulative induced
        velocity at each of the N points due to all the ring vortices. The units are
        meters per second.
    """
    listStackSlvp_G_Cg = [
        stackBrrvp_G_Cg,
        stackFrrvp_G_Cg,
        stackFlrvp_G_Cg,
        stackBlrvp_G_Cg,
    ]
    listStackElvp_G_Cg = [
        stackFrrvp_G_Cg,
        stackFlrvp_G_Cg,
        stackBlrvp_G_Cg,
        stackBrrvp_G_Cg,
    ]

    stackVInd_G__E = np.zeros((stackP_G_Cg.shape[0], 3))

    # Get the velocity induced by each leg of the RingVortex.
    for i in range(4):
        stackVInd_G__E += _collapsed_velocities_from_line_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackSlvp_G_Cg=listStackSlvp_G_Cg[i],
            stackElvp_G_Cg=listStackElvp_G_Cg[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return stackVInd_G__E


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_ring_vortices_chordwise_segments(
    stackP_G_Cg,
    stackBrrvp_G_Cg,
    stackFrrvp_G_Cg,
    stackFlrvp_G_Cg,
    stackBlrvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    RingVortices. At every point, it finds the cumulative induced velocity due to all
    the RingVortices' chordwise segments.

    Note: See the docstring for collapsed_velocities_from_line_vortices for
    information on usage.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    # TODO: Update docstring parameters' descriptions and types.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackBrrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's back right vertex's position in meters.
    :param stackFrrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front right vertex's position in meters.
    :param stackFlrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param stackBlrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of ring
        vortices. Each holds the strength of that ring vortex in meters squared per
        second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that ring vortex in seconds. This
        is only relevant for vortices that have been shed into the wake. The default
        value is None. If the age of a specific vortex is 0.0 seconds, then the
        vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 2D array of floats
        This is an array of shape (N x 3), and it holds the cumulative induced
        velocity at each of the N points due to all the ring vortices' spanwise
        segments. The units are meters per second.
    """
    listStackSlvp_G_Cg = [
        stackBrrvp_G_Cg,
        stackFlrvp_G_Cg,
    ]
    listStackElvp_G_Cg = [
        stackFrrvp_G_Cg,
        stackBlrvp_G_Cg,
    ]

    stackVInd_G__E = np.zeros((stackP_G_Cg.shape[0], 3))

    # Get the velocity induced by each leg of the RingVortex.
    for i in range(2):
        stackVInd_G__E += _collapsed_velocities_from_line_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackSlvp_G_Cg=listStackSlvp_G_Cg[i],
            stackElvp_G_Cg=listStackElvp_G_Cg[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return stackVInd_G__E


@njit(cache=True, fastmath=False)
def expanded_velocities_from_ring_vortices(
    stackP_G_Cg,
    stackBrrvp_G_Cg,
    stackFrrvp_G_Cg,
    stackFlrvp_G_Cg,
    stackBlrvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    RingVortices. At every point, it finds the induced velocity due to each RingVortex.

    Note: See the docstring for collapsed_velocities_from_line_vortices for
    information on usage.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    # TODO: Update docstring parameters' descriptions and types.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackBrrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's back right vertex's position in meters.
    :param stackFrrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front right vertex's position in meters.
    :param stackFlrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param stackBlrvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of ring
        vortices. Each holds the strength of that ring vortex in meters squared per
        second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that ring vortex in seconds. This
        is only relevant for vortices that have been shed into the wake. The default
        value is None. If the age of a specific vortex is 0.0 seconds, then the
        vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 3D array of floats
        This is an array of shape (N x M x 3), where each row/column pair identifies
        the velocity induced at one point by one of the ring vortices. The units are
        meters per second.
    """
    listStackSlvp_G_Cg = [
        stackBrrvp_G_Cg,
        stackFrrvp_G_Cg,
        stackFlrvp_G_Cg,
        stackBlrvp_G_Cg,
    ]
    listStackElvp_G_Cg = [
        stackFrrvp_G_Cg,
        stackFlrvp_G_Cg,
        stackBlrvp_G_Cg,
        stackBrrvp_G_Cg,
    ]

    gridVInd_G__E = np.zeros((stackP_G_Cg.shape[0], strengths.shape[0], 3))

    # Get the velocity induced by each leg of the RingVortex.
    for i in range(4):
        gridVInd_G__E += _expanded_velocities_from_line_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackSlvp_G_Cg=listStackSlvp_G_Cg[i],
            stackElvp_G_Cg=listStackElvp_G_Cg[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return gridVInd_G__E


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_horseshoe_vortices(
    stackP_G_Cg,
    stackBrhvp_G_Cg,
    stackFrhvp_G_Cg,
    stackFlhvp_G_Cg,
    stackBlhvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    HorseshoeVortices. At every point, it finds the cumulative induced velocity due to
    all the HorseshoeVortices.

    Note: See the docstring for collapsed_velocities_from_line_vortices for
    information on usage.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    # TODO: Update docstring parameters' descriptions and types.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackBrhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's back right vertex's position in meters.
    :param stackFrhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front right vertex's position in meters.
    :param stackFlhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front left vertex's position in meters.
    :param stackBlhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front left vertex's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of horseshoe
        vortices. Each holds the strength of that horseshoe vortex in meters squared
        per second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that horseshoe vortex in seconds.
        This is only relevant for vortices that have been shed into the wake. The
        default value is None. If the age of a specific vortex is 0.0 seconds,
        then the vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 2D array of floats
        This is an array of shape (N x 3), and it holds the cumulative induced
        velocity at each of the N points due to all the horseshoe vortices. The units
        are meters per second.
    """
    listStackSlvp_G_Cg = [
        stackBrhvp_G_Cg,
        stackFrhvp_G_Cg,
        stackFlhvp_G_Cg,
    ]
    listStackElvp_G_Cg = [
        stackFrhvp_G_Cg,
        stackFlhvp_G_Cg,
        stackBlhvp_G_Cg,
    ]

    stackVInd_G__E = np.zeros((stackP_G_Cg.shape[0], 3))

    # Get the velocity induced by each leg of the RingVortex.
    for i in range(3):
        stackVInd_G__E += _collapsed_velocities_from_line_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackSlvp_G_Cg=listStackSlvp_G_Cg[i],
            stackElvp_G_Cg=listStackElvp_G_Cg[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return stackVInd_G__E


@njit(cache=True, fastmath=False)
def expanded_velocities_from_horseshoe_vortices(
    stackP_G_Cg,
    stackBrhvp_G_Cg,
    stackFrhvp_G_Cg,
    stackFlhvp_G_Cg,
    stackBlhvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    HorseshoeVortices. At every point, it finds the induced velocity due to each
    HorseshoeVortex.

    Note: See the docstring for collapsed_velocities_from_line_vortices for
    information on usage.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    # TODO: Update docstring parameters' descriptions and types.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackBrhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's back right vertex's position in meters.
    :param stackFrhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front right vertex's position in meters.
    :param stackFlhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front left vertex's position in meters.
    :param stackBlhvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front left vertex's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of horseshoe
        vortices. Each holds the strength of that horseshoe vortex in meters squared
        per second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that horseshoe vortex in seconds.
        This is only relevant for vortices that have been shed into the wake. The
        default value is None. If the age of a specific vortex is 0.0 seconds,
        then the vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 2D array of floats
        This is an array of shape (N x M x 3), where each row/column pair identifies
        the velocity induced at one point by one of the horseshoe vortices. The units
        are meters per second.
    """
    listStackSlvp_G_Cg = [
        stackBrhvp_G_Cg,
        stackFrhvp_G_Cg,
        stackFlhvp_G_Cg,
    ]
    listStackElvp_G_Cg = [
        stackFrhvp_G_Cg,
        stackFlhvp_G_Cg,
        stackBlhvp_G_Cg,
    ]

    gridVInd_G__E = np.zeros((stackP_G_Cg.shape[0], strengths.shape[0], 3))

    # Get the velocity induced by each leg of the RingVortex.
    for i in range(3):
        gridVInd_G__E += _expanded_velocities_from_line_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackSlvp_G_Cg=listStackSlvp_G_Cg[i],
            stackElvp_G_Cg=listStackElvp_G_Cg[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return gridVInd_G__E


@njit(cache=True, fastmath=False)
def _collapsed_velocities_from_line_vortices(
    stackP_G_Cg,
    stackSlvp_G_Cg,
    stackElvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    LineVortices. At every point, it finds the cumulative induced velocity due to all
    the LineVortices.


    Example Usage (single point case):
    ------------------------------------------------------------------------------------

    Returns vIndP_G__E: the velocity induced at P_G_Cg, a point on the body (in
    geometry axes, relative to the CG), due solely to the supplied LineVortex.
    vIndP_G__E is in geometry axes, and observed from the Earth frame. The velocity
    returned is as observed from the Earth frame because the Biot-Savart law
    calculates the velocity induced by the LineVortex, as observed from an inertial
    frame where the induced velocity approaches zero as distance from the LineVortex
    approaches infinity. We choose this frame to be the inertial Earth frame.

    To obtain the apparent velocity "felt" by the point P_G_Cg:

    1.) Calculate flow velocity at P (in geometry axes, observed from the Earth
    frame), by summing the velocity induced by all LineVortices:

        vFlowP_G__E = ∑[vIndP_G__E(P_G_Cg, Slvp_G_Cg, Elvp_G_Cg, strength)]

    2.) Find the material-point velocity at P (in geometry axes, observed from the
    Earth frame):

        vMatP_G__E = vCg_G__E + omegaCg_G__E x P_G_Cg + vDeformP_G__E

    3.) Combine for apparent ("felt") velocity at P (in geometry axes, observed from
    the Earth frame):

        vFeltP_G__E = vFlowP_G__E - vMatP_G__E

    Conventions:

    - All vectors are 3-vectors expressed in geometry axes.

    - "x" denotes the 3D cross product.

    Variables (inputs, geometry, and constants)

    - G: The Airplane's geometry axis system (geometry axes), with basis directions
    +x = aft, +y = right, and +z = up.

    - Cg: The Airplane's CG, which serves as the reference point used for relative
    positions.

    - E: Inertial Earth reference frame where the observed induced velocity from the
    LineVortex approaches zero as distance from the LineVortex approaches infinity.

    - Slvp_G_Cg: start point of the specified LineVortex segment (in geometry axes,
    relative to the CG) (m).

    - Elvp_G_Cg: end point of the specified LineVortex segment (in geometry axes,
    relative to the CG) (m).

    - P_G_Cg: evaluation point on the body (in geometry axes, relative to the CG) (
    m). It is also where vFlowP_G__E, vMatP_G__E, vFeltP_G__E are evaluated.

    - strength: LineVortex strength (circulation). It is a frame-independent scalar
    with units m^2/s.

    - LineVortex: a finite straight vortex segment defined by Slvp_G_Cg, Elvp_G_Cg,
    and strength (circulation).

    Velocities (all in geometry axes, observed from the Earth frame)

    - vIndP_G__E: induced velocity at P_G_Cg due only to the specified LineVortex (m/s).

    - vFlowP_G__E: total flow velocity at P_G_Cg (sum of all LineVortices'
    vIndP_G__E) (m/s). We assume that the external fluid is quiescent (no wind or
    current), so this is only composed of the induced velocities.

    - vCg_G__E: rigid translational velocity of the body at the CG (m/s).

    - omegaCg_G__E: rigid angular velocity of the body (about the CG) (rad/s).

    - vDeformP_G__E: additional local velocity of point P due to deformation/flapping
    relative to the body’s rigid motion (m/s).

    - vMatP_G__E: inertial velocity of the material point P (rigid + deform) =
    vCg_G__E + omegaCg_G__E x P_G_Cg + vDeformP_G__E (m/s).

    - vFeltP_G__E: apparent velocity "felt" by the material point P, equal to
    vFlowP_G__E - vMatP_G__E (m/s).

    ------------------------------------------------------------------------------------

    Citation: The equations in this function are from "Extended Unsteady
    Vortex-Lattice Method for Insect Flapping Wings" (Nguyen et al., 2016)

    Note: This function uses a modified version of the Bio-Savart law to create a
    smooth induced velocity decay based on a vortex's core radius. The radius is
    determined based on a vortex's age and kinematic viscosity. If the age of the
    vortex is 0.0 seconds, the radius is set to 0.0 meters. The age of a vortex in
    only relevant for vortices that have been shed into the wake.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    # TODO: Update docstring parameters' types and descriptions.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackSlvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's Slvp_G_Cg's position in meters.
    :param stackElvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's Elvp_G_Cg's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the strength of that line vortex in meters
        squared per second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that line vortex in seconds. This
        is only relevant for vortices that have been shed into the wake. The default
        value is None. If the age of a specific vortex is 0.0 seconds, then the
        vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 2D array of floats
        This is an array of shape (N x 3), and it holds the cumulative induced
        velocity at each of the N points due to all the line vortices. The units are
        meters per second.
    """
    num_vortices = stackSlvp_G_Cg.shape[0]
    num_points = stackP_G_Cg.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in
    # geometry axes, observed from the Earth frame).
    stackVInd_G__E = np.zeros((num_points, 3))

    # If the user didn't specify any ages, set the age of each LineVortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_G_Cg = stackSlvp_G_Cg[vortex_id]
        Elvp_G_Cg = stackElvp_G_Cg[vortex_id]
        strength = strengths[vortex_id]
        age = ages[vortex_id]

        # Calculate the radius of the LineVortex's core. If the age is 0.0 seconds,
        # this will evaluate to be 0.0 meters.
        r_c = 2 * math.sqrt(_lamb * (nu + _squire * abs(strength)) * age)

        # The r0_G vector goes from the LineVortex's start point to its end point (in
        # geometry axes).
        r0X_G = Elvp_G_Cg[0] - Slvp_G_Cg[0]
        r0Y_G = Elvp_G_Cg[1] - Slvp_G_Cg[1]
        r0Z_G = Elvp_G_Cg[2] - Slvp_G_Cg[2]

        # Find r0_G's length.
        r0 = math.sqrt(r0X_G**2 + r0Y_G**2 + r0Z_G**2)

        c_1 = strength / (4 * math.pi)
        c_2 = r0**2 * r_c**2

        for point_id in range(num_points):
            P_G_Cg = stackP_G_Cg[point_id]

            # The r1_G vector goes from P_G_Cg to the LineVortex's start point (in
            # geometry axes).
            r1X_G = Slvp_G_Cg[0] - P_G_Cg[0]
            r1Y_G = Slvp_G_Cg[1] - P_G_Cg[1]
            r1Z_G = Slvp_G_Cg[2] - P_G_Cg[2]

            # The r2_G vector goes from P_G_Cg to the LineVortex's end point (in
            # geometry axes).
            r2X_G = Elvp_G_Cg[0] - P_G_Cg[0]
            r2Y_G = Elvp_G_Cg[1] - P_G_Cg[1]
            r2Z_G = Elvp_G_Cg[2] - P_G_Cg[2]

            # The r3_G vector is the cross product of r1_G and r2_G (in geometry axes).
            r3X_G = r1Y_G * r2Z_G - r1Z_G * r2Y_G
            r3Y_G = r1Z_G * r2X_G - r1X_G * r2Z_G
            r3Z_G = r1X_G * r2Y_G - r1Y_G * r2X_G

            # Find the lengths of r1_G, r2_G, and r3_G.
            r1 = math.sqrt(r1X_G**2 + r1Y_G**2 + r1Z_G**2)
            r2 = math.sqrt(r2X_G**2 + r2Y_G**2 + r2Z_G**2)
            r3 = math.sqrt(r3X_G**2 + r3Y_G**2 + r3Z_G**2)

            c_3 = r1X_G * r2X_G + r1Y_G * r2Y_G + r1Z_G * r2Z_G

            # If part of the LineVortex is so close to P_G_Cg that they are touching
            # (within machine epsilon), there is a removable discontinuity. In this
            # case, continue to the next point because there is no velocity induced
            # by the current LineVortex at this point.
            if r1 < _eps or r2 < _eps or r3**2 < _eps:
                continue
            else:
                c_4 = c_1 * (r1 + r2) * (r1 * r2 - c_3) / (r1 * r2 * (r3**2 + c_2))
                stackVInd_G__E[point_id, 0] += c_4 * r3X_G
                stackVInd_G__E[point_id, 1] += c_4 * r3Y_G
                stackVInd_G__E[point_id, 2] += c_4 * r3Z_G
    return stackVInd_G__E


@njit(cache=True, fastmath=False)
def _expanded_velocities_from_line_vortices(
    stackP_G_Cg,
    stackSlvp_G_Cg,
    stackElvp_G_Cg,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    LineVortices. At every point, it finds the induced velocity due to each LineVortex.

    Note: See the docstring for collapsed_velocities_from_line_vortices for
    information on usage.

    Citation: The equations in this function are from "Extended Unsteady
    Vortex-Lattice Method for Insect Flapping Wings" (Nguyen et al., 2016)

    Note: This function uses a modified version of the Bio-Savart law to create a
    smooth induced velocity decay based on a vortex's core radius. The radius is
    determined based on a vortex's age and kinematic viscosity. If the age of the
    vortex is 0.0 seconds, the radius is set to 0.0 meters. The age of a vortex in
    only relevant for vortices that have been shed into the wake.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

     # TODO: Update docstring parameters' types and descriptions.
    :param stackP_G_Cg: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param stackSlvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's Slvp_G_Cg's position in meters.
    :param stackElvp_G_Cg: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's Elvp_G_Cg's position in meters.
    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the strength of that line vortex in meters
        squared per second.
    :param ages: 1D array of floats, optional
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the age of that line vortex in seconds. This
        is only relevant for vortices that have been shed into the wake. The default
        value is None. If the age of a specific vortex is 0.0 seconds, then the
        vortex core radius is set to 0.0 meters.
    :param nu: float, optional
        This variable is a float that represents the kinematic viscosity of the fluid
        in meters squared per second. The default value is 0.0 meters squared per
        second.
    :return velocities: 3D array of floats
        This is an array of shape (N x M x 3), where each row/column pair identifies
        the velocity induced at one point by one of the line vortices. The units are
        meters per second.
    """
    num_vortices = stackSlvp_G_Cg.shape[0]
    num_points = stackP_G_Cg.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in
    # geometry axes, observed from the Earth frame).
    gridVInd_G__E = np.zeros((num_points, num_vortices, 3))

    # If the user didn't specify any ages, set the age of each LineVortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_G_Cg = stackSlvp_G_Cg[vortex_id]
        Elvp_G_Cg = stackElvp_G_Cg[vortex_id]
        strength = strengths[vortex_id]
        age = ages[vortex_id]

        # Calculate the radius of the LineVortex's core. If the age is 0.0 seconds,
        # this will evaluate to be 0.0 meters.
        r_c = 2 * math.sqrt(_lamb * (nu + _squire * abs(strength)) * age)

        # The r0_G vector goes from the LineVortex's start point to its end point (in
        # geometry axes).
        r0X_G = Elvp_G_Cg[0] - Slvp_G_Cg[0]
        r0Y_G = Elvp_G_Cg[1] - Slvp_G_Cg[1]
        r0Z_G = Elvp_G_Cg[2] - Slvp_G_Cg[2]

        # Find r0_G's length.
        r0 = math.sqrt(r0X_G**2 + r0Y_G**2 + r0Z_G**2)

        c_1 = strength / (4 * math.pi)
        c_2 = r0**2 * r_c**2

        for point_id in range(num_points):
            P_G_Cg = stackP_G_Cg[point_id]

            # The r1_G vector goes from P_G_Cg to the LineVortex's start point (in
            # geometry axes).
            r1X_G = Slvp_G_Cg[0] - P_G_Cg[0]
            r1Y_G = Slvp_G_Cg[1] - P_G_Cg[1]
            r1Z_G = Slvp_G_Cg[2] - P_G_Cg[2]

            # The r2_G vector goes from P_G_Cg to the LineVortex's end point (in
            # geometry axes).
            r2X_G = Elvp_G_Cg[0] - P_G_Cg[0]
            r2Y_G = Elvp_G_Cg[1] - P_G_Cg[1]
            r2Z_G = Elvp_G_Cg[2] - P_G_Cg[2]

            # The r3_G vector is the cross product of r1_G and r2_G (in geometry axes).
            r3X_G = r1Y_G * r2Z_G - r1Z_G * r2Y_G
            r3Y_G = r1Z_G * r2X_G - r1X_G * r2Z_G
            r3Z_G = r1X_G * r2Y_G - r1Y_G * r2X_G

            # Find the lengths of r1_G, r2_G, and r3_G.
            r1 = math.sqrt(r1X_G**2 + r1Y_G**2 + r1Z_G**2)
            r2 = math.sqrt(r2X_G**2 + r2Y_G**2 + r2Z_G**2)
            r3 = math.sqrt(r3X_G**2 + r3Y_G**2 + r3Z_G**2)

            c_3 = r1X_G * r2X_G + r1Y_G * r2Y_G + r1Z_G * r2Z_G

            # If part of the LineVortex is so close to P_G_Cg that they are touching
            # (within machine epsilon), there is a removable discontinuity. In this
            # case, continue to the next point because there is no velocity induced
            # by the current LineVortex at this point.
            if r1 < _eps or r2 < _eps or r3**2 < _eps:
                continue
            else:
                c_4 = c_1 * (r1 + r2) * (r1 * r2 - c_3) / (r1 * r2 * (r3**2 + c_2))
                gridVInd_G__E[point_id, vortex_id, 0] = c_4 * r3X_G
                gridVInd_G__E[point_id, vortex_id, 1] = c_4 * r3Y_G
                gridVInd_G__E[point_id, vortex_id, 2] = c_4 * r3Z_G
    return gridVInd_G__E
