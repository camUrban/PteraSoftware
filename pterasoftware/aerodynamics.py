"""This module contains vortex class definitions, and useful aerodynamic functions.

This module contains the following classes:
    LineVortex: This class is used to contain line vortices.

    HorseshoeVortex: This class is used to contain horseshoe vortices.

    RingVortex: This class is used to contain ring vortices.

This module contains the following exceptions:
    None

This module contains the following functions:
    collapsed_velocities_from_horseshoe_vortices: This function takes in a group of
    points, and the attributes of a group of horseshoe vortices. At every point,
    it finds the cumulative induced velocity due to all the horseshoe vortices.

    expanded_velocities_from_horseshoe_vortices: This function takes in a group of
    points, and the attributes of a group of horseshoe vortices. At every point,
    it finds the induced velocity due to each horseshoe vortex.

    collapsed_velocities_from_ring_vortices: This function takes in a group of
    points, and the attributes of a group of ring vortices. At every point, it finds
    the cumulative induced velocity due to all the ring vortices.

    collapsed_velocities_from_ring_vortices_chordwise_segments: This function takes
    in a group of points, and the attributes of a group of ring vortices. At every
    point, it finds the cumulative induced velocity due to all the ring vortices'
    chordwise segments.

    expanded_velocities_from_ring_vortices: This function takes in a group of points,
    and the attributes of a group of ring vortices. At every point, it finds the
    induced velocity due to each ring vortex.

    collapsed_velocities_from_line_vortices: This function takes in a group of
    points, and the attributes of a group of line vortices. At every point, it finds
    the cumulative induced velocity due to all the line vortices.

    expanded_velocities_from_line_vortices: This function takes in a group of points,
    and the attributes of a group of line vortices. At every point, it finds the
    induced velocity due to each line vortex.
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
squire = 10**-4

# Set the value of Lamb's constant that will be used by the induced velocity
# functions. Lamb's constant relates to the size of the vortex cores and the rate at
# which they grow. The value of this parameter is well agreed upon, and published in
# "Extended Unsteady Vortex-Lattice Method for Insect Flapping Wings" (Nguyen et al.,
# 2016). It is unitless.
lamb = 1.25643

# Set the value of the local machine error. This will be used to fix removable
# discontinuities in the induced velocity functions.
eps = np.finfo(float).eps


class LineVortex:
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

        :param Slvp_G_Cg: array-like of 3 numbers Position [x, y, z] of the LineVortex's
            start point (in geometry axes, relative to the CG). Can be a list, tuple, or
            numpy array of numbers (int or float). Values are converted to floats
            internally. The units are in meters.
        :param Elvp_G_Cg: array-like of 3 numbers Position [x, y, z] of the LineVortex's
            end point (in geometry axes, relative to the CG). Can be a list, tuple, or
            numpy array of numbers (int or float). Values are converted to floats
            internally. The units are in meters.
        :param strength: number This is the strength of the LineVortex. It can be a
            positive int or float and will be converted to a float internally. The units
            are in meters squared per second.
        """
        self.Slvp_G_Cg = parameter_validation.positive_number_return_float(
            Slvp_G_Cg, "Slvp_G_Cg"
        )
        self.Elvp_G_Cg = parameter_validation.positive_number_return_float(
            Elvp_G_Cg, "Elvp_G_Cg"
        )
        self.strength = parameter_validation.positive_number_return_float(
            strength, "strength"
        )

        # Initialize variables to hold the vector from the LineVortex's start to end
        # point (in geometry axes), and its center point (in geometry axes, relative
        # to the CG).
        self.vector_G = self.Elvp_G_Cg - self.Slvp_G_Cg
        self.Clvp_G_Cg = self.Slvp_G_Cg + 0.5 * self.vector_G


# NOTE: I'm in the process of refactoring this class.
class HorseshoeVortex:
    """This class is used to contain horseshoe vortices.

    :param Frhvp_G_Cg: Position [x, y, z] of the HorseshoeVortex's front-right point (in
        geometry axes, relative to the CG). The front-right point is defined as the
        start point of the HorseshoeVortex's front leg, which is also its one finite
        leg. Can be a list, tuple, or numpy array of numbers (int or float). Values are
        converted to floats internally. The units are in meters.
    :type Frhvp_G_Cg: array-like of 3 numbers
    :param Flhvp_G_Cg: Position [x, y, z] of the HorseshoeVortex's front-left point (in
        geometry axes, relative to the CG). The front-left point is defined as the end
        point of the HorseshoeVortex's front leg, which is also its one finite leg. Can
        be a list, tuple, or numpy array of numbers (int or float). Values are converted
        to floats internally. The units are in meters.
    :type Flhvp_G_Cg: array-like of 3 numbers
    :param leftLegVector_G: Direction vector of the HorseshoeVortex's left leg (in
        geometry axes). The left leg starts from the front-left point and ends at the
        back-left point. It is one of the HorseshoeVortex's two quasi-infinite legs, the
        other being the right-leg. It can be a list, tuple, or numpy array of numbers
        (int or float). Values are converted to floats internally. If this isn't already
        a unit vector, it will be converted to one during initialization. The right-
        leg's vector (in geometry axes) is defined as -1 times this vector. The units
        are in meters.
    :type leftLegVector_G: array-like of 3 numbers
    :param left_right_leg_lengths: This is the length of the HorseshoeVortex's left and
        right quasi-infinite legs. It must be a positive number and will be converted
        internally to a float. I recommend setting it to at least 20 times the length of
        the finite leg. The units are in meters.
    :type left_right_leg_lengths: number
    :param strength: This is the strength of the HorseshoeVortex It must be a positive
        number and will be converted internally to a float. Its units are in meters
        squared per second.
    :type strength: number :raises [ErrorType]: [ErrorDescription]
    :return: [ ReturnDescription]
    :rtype: [ReturnType]
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
        self.strength = parameter_validation.positive_number_return_float(
            strength, "strength"
        )

        self.right_leg_origin = (
            self.Frhvp_G_Cg + leftLegVector_G * left_right_leg_lengths
        )
        self.left_leg_termination = (
            self.Flhvp_G_Cg + leftLegVector_G * left_right_leg_lengths
        )

        # Initialize LineVortices to represent the HorseshoeVortex's legs.
        self.right_leg = LineVortex(
            Slvp_G_Cg=self.right_leg_origin,
            Elvp_G_Cg=self.Frhvp_G_Cg,
            strength=self.strength,
        )
        self.finite_leg = LineVortex(
            Slvp_G_Cg=self.Frhvp_G_Cg,
            Elvp_G_Cg=self.Flhvp_G_Cg,
            strength=self.strength,
        )
        self.left_leg = LineVortex(
            Slvp_G_Cg=self.Flhvp_G_Cg,
            Elvp_G_Cg=self.left_leg_termination,
            strength=self.strength,
        )

    def update_strength(self, strength):
        """This method updates the strength of this HorseshoeVortex object, and the
        strength of its leg LineVortices.

        :param strength: This is the strength of this vortex, and of its line vortex
            legs. Its units are meters squared per second.
        :type strength: float
        """
        self.strength = strength
        self.right_leg.strength = strength
        self.finite_leg.strength = strength
        self.left_leg.strength = strength


class RingVortex:
    """This class is used to contain ring vortices.

    This class contains the following public methods:
        update_strength: This method updates the strength of this ring vortex object,
        and the strength of its four legs' line vortex objects.

        update_position: This method updates the position of the ring vortex, and the
        positions of all its attributes.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        front_left_vertex,
        front_right_vertex,
        back_left_vertex,
        back_right_vertex,
        strength,
    ):
        """This is the initialization method.

        :param front_left_vertex: 1D array This is a vector_G containing the x, y, and z
            coordinates of the vortex's front left point. It's a (,3) array with units
            of meters.
        :param front_right_vertex: 1D array This is a vector_G containing the x, y, and
            z coordinates of the vortex's front right point. It's a (,3) array with
            units of meters.
        :param back_left_vertex: 1D array This is a vector_G containing the x, y, and z
            coordinates of the vortex's back left point. It's a (,3) array with units of
            meters.
        :param back_right_vertex: 1D array This is a vector_G containing the x, y, and z
            coordinates of the vortex's back right point. It's a (,3) array with units
            of meters.
        :param strength: float This is the strength of the vortex in meters squared per
            second.
        """
        self.front_left_vertex = front_left_vertex
        self.front_right_vertex = front_right_vertex
        self.back_left_vertex = back_left_vertex
        self.back_right_vertex = back_right_vertex
        self.strength = strength

        # Initialize the line vortices that make up the ring vortex.
        self.front_leg = LineVortex(
            origin=self.front_right_vertex,
            termination=self.front_left_vertex,
            strength=self.strength,
        )
        self.left_leg = LineVortex(
            origin=self.front_left_vertex,
            termination=self.back_left_vertex,
            strength=self.strength,
        )
        self.back_leg = LineVortex(
            origin=self.back_left_vertex,
            termination=self.back_right_vertex,
            strength=self.strength,
        )
        self.right_leg = LineVortex(
            origin=self.back_right_vertex,
            termination=self.front_right_vertex,
            strength=self.strength,
        )

        # Initialize a variable to hold the centroid of the ring vortex.
        self.center = functions.numba_centroid_of_quadrilateral(
            self.front_left_vertex,
            self.front_right_vertex,
            self.back_left_vertex,
            self.back_right_vertex,
        )

        # Initialize a variable to hold the age of the ring vortex in seconds.
        self.age = 0

    def update_strength(self, strength):
        """This method updates the strength of this ring vortex object, and the strength
        of its four legs' line vortex objects.

        :param strength: float This is the strength of this vortex, and of its four
            legs' line vortices. Its units are meters squared per second.
        :return: None
        """
        self.strength = strength
        self.right_leg.strength = strength
        self.front_leg.strength = strength
        self.left_leg.strength = strength
        self.back_leg.strength = strength

    def update_position(
        self, front_left_vertex, front_right_vertex, back_left_vertex, back_right_vertex
    ):
        """This method updates the position of the ring vortex, and the positions of all
        its attributes.

        :param front_left_vertex: 1D array This is the new position of the x, y, and z
            coordinates of the front left vertex. It is a (,3) array with units of
            meters.
        :param front_right_vertex: 1D array This is the new position of the x, y, and z
            coordinates of the front right vertex. It is a (,3) array with units of
            meters.
        :param back_left_vertex: 1D array This is the new position of the x, y, and z
            coordinates of the back left vertex. It is a (,3) array with units of
            meters.
        :param back_right_vertex: 1D array This is the new position of the x, y, and z
            coordinates of the back right vertex. It is a (,3) array with units of
            meters.
        :return: None
        """
        self.front_left_vertex = front_left_vertex
        self.front_right_vertex = front_right_vertex
        self.back_left_vertex = back_left_vertex
        self.back_right_vertex = back_right_vertex

        # Initialize the line vortices that make up the ring vortex.
        self.front_leg = LineVortex(
            origin=self.front_right_vertex,
            termination=self.front_left_vertex,
            strength=self.strength,
        )
        self.left_leg = LineVortex(
            origin=self.front_left_vertex,
            termination=self.back_left_vertex,
            strength=self.strength,
        )
        self.back_leg = LineVortex(
            origin=self.back_left_vertex,
            termination=self.back_right_vertex,
            strength=self.strength,
        )
        self.right_leg = LineVortex(
            origin=self.back_right_vertex,
            termination=self.front_right_vertex,
            strength=self.strength,
        )

        # Initialize a variable to hold the centroid of the ring vortex.
        self.center = functions.numba_centroid_of_quadrilateral(
            self.front_left_vertex,
            self.front_right_vertex,
            self.back_left_vertex,
            self.back_right_vertex,
        )


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_horseshoe_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    horseshoe vortices. At every point, it finds the cumulative induced velocity due to
    all the horseshoe vortices.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param back_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D array of floats
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
    origins_list = [
        back_right_vortex_vertices,
        front_right_vortex_vertices,
        front_left_vortex_vertices,
    ]
    terminations_list = [
        front_right_vortex_vertices,
        front_left_vortex_vertices,
        back_left_vortex_vertices,
    ]
    induced_velocities = np.zeros((points.shape[0], 3))

    # Get the velocity induced by each leg of the ring vortex.
    for i in range(3):
        induced_velocities += collapsed_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=False)
def expanded_velocities_from_horseshoe_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    horseshoe vortices. At every point, it finds the induced velocity due to each
    horseshoe vortex.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param back_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of
        horseshoe vortices. Each row contains the x, y, and z float coordinates of
        that horseshoe vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D array of floats
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
    origins_list = [
        back_right_vortex_vertices,
        front_right_vortex_vertices,
        front_left_vortex_vertices,
    ]
    terminations_list = [
        front_right_vortex_vertices,
        front_left_vortex_vertices,
        back_left_vortex_vertices,
    ]
    induced_velocities = np.zeros((points.shape[0], strengths.shape[0], 3))

    # Get the velocity induced by each leg of the ring vortex.
    for i in range(3):
        induced_velocities += expanded_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_ring_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of ring
    vortices. At every point, it finds the cumulative induced velocity due to all the
    ring vortices.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param back_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D array of floats
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
    origins_list = [
        back_right_vortex_vertices,
        front_right_vortex_vertices,
        front_left_vortex_vertices,
        back_left_vortex_vertices,
    ]
    terminations_list = [
        front_right_vortex_vertices,
        front_left_vortex_vertices,
        back_left_vortex_vertices,
        back_right_vortex_vertices,
    ]
    induced_velocities = np.zeros((points.shape[0], 3))

    # Get the velocity induced by each leg of the ring vortex.
    for i in range(4):
        induced_velocities += collapsed_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_ring_vortices_chordwise_segments(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of ring
    vortices. At every point, it finds the cumulative induced velocity due to all the
    ring vortices' chordwise segments.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param back_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D array of floats
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
    origins_list = [
        back_right_vortex_vertices,
        front_left_vortex_vertices,
    ]
    terminations_list = [
        front_right_vortex_vertices,
        back_left_vortex_vertices,
    ]
    induced_velocities = np.zeros((points.shape[0], 3))

    # Get the velocity induced by each leg of the ring vortex.
    for i in range(2):
        induced_velocities += collapsed_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=False)
def expanded_velocities_from_ring_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of ring
    vortices. At every point, it finds the induced velocity due to each ring vortex.

    Note: This function's performance has been highly optimized for unsteady
    simulations via Numba. While using Numba dramatically increases unsteady
    simulation performance, it does cause a performance drop for the less intense
    steady simulations.

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param back_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of ring
        vortices. Each row contains the x, y, and z float coordinates of that ring
        vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D array of floats
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
    origins_list = [
        back_right_vortex_vertices,
        front_right_vortex_vertices,
        front_left_vortex_vertices,
        back_left_vortex_vertices,
    ]
    terminations_list = [
        front_right_vortex_vertices,
        front_left_vortex_vertices,
        back_left_vortex_vertices,
        back_right_vortex_vertices,
    ]
    induced_velocities = np.zeros((points.shape[0], strengths.shape[0], 3))

    # Get the velocity induced by each leg of the ring vortex.
    for i in range(4):
        induced_velocities += expanded_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=False)
def collapsed_velocities_from_line_vortices(
    points,
    origins,
    terminations,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of line
    vortices. At every point, it finds the cumulative induced velocity due to all the
    line vortices.

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

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param origins: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's Slvp_G_Cg's position in meters.
    :param terminations: 2D array of floats
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
    num_vortices = origins.shape[0]
    num_points = points.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities.
    velocities = np.zeros((num_points, 3))

    # If the user didn't specify any ages, set the age of each vortex to 0.0 seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        origin = origins[vortex_id]
        termination = terminations[vortex_id]
        strength = strengths[vortex_id]
        age = ages[vortex_id]

        # Calculate the radius of the vortex's core. If the age is 0.0 seconds,
        # this will evaluate to be 0.0 meters.
        r_c = 2 * math.sqrt(lamb * (nu + squire * abs(strength)) * age)

        # The r_0 vector_G goes from the line vortex's Slvp_G_Cg to its Elvp_G_Cg.
        r_0_x = termination[0] - origin[0]
        r_0_y = termination[1] - origin[1]
        r_0_z = termination[2] - origin[2]

        # Find the r_0 vector_G's length.
        r_0 = math.sqrt(r_0_x**2 + r_0_y**2 + r_0_z**2)

        c_1 = strength / (4 * math.pi)
        c_2 = r_0**2 * r_c**2

        for point_id in range(num_points):
            point = points[point_id]

            # The r_1 vector_G goes from the point to the line vortex's Slvp_G_Cg.
            r_1_x = origin[0] - point[0]
            r_1_y = origin[1] - point[1]
            r_1_z = origin[2] - point[2]

            # The r_2 vector_G goes from the point to the line vortex's Elvp_G_Cg.
            r_2_x = termination[0] - point[0]
            r_2_y = termination[1] - point[1]
            r_2_z = termination[2] - point[2]

            # The r_3 vector_G is the cross product of the r_1 and r_2 vectors.
            r_3_x = r_1_y * r_2_z - r_1_z * r_2_y
            r_3_y = r_1_z * r_2_x - r_1_x * r_2_z
            r_3_z = r_1_x * r_2_y - r_1_y * r_2_x

            # Find the r_1, r_2, and r_3 vectors' lengths.
            r_1 = math.sqrt(r_1_x**2 + r_1_y**2 + r_1_z**2)
            r_2 = math.sqrt(r_2_x**2 + r_2_y**2 + r_2_z**2)
            r_3 = math.sqrt(r_3_x**2 + r_3_y**2 + r_3_z**2)

            c_3 = r_1_x * r_2_x + r_1_y * r_2_y + r_1_z * r_2_z

            # If part of the vortex is so close to the point that they are touching (
            # within machine epsilon), there is a removable discontinuity. In this
            # case, continue to the next point because there is no velocity induced
            # by the current vortex at this point.
            if r_1 < eps or r_2 < eps or r_3**2 < eps:
                continue
            else:
                c_4 = (
                    c_1 * (r_1 + r_2) * (r_1 * r_2 - c_3) / (r_1 * r_2 * (r_3**2 + c_2))
                )
                velocities[point_id, 0] += c_4 * r_3_x
                velocities[point_id, 1] += c_4 * r_3_y
                velocities[point_id, 2] += c_4 * r_3_z

    return velocities


@njit(cache=True, fastmath=False)
def expanded_velocities_from_line_vortices(
    points,
    origins,
    terminations,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of line
    vortices. At every point, it finds the induced velocity due to each line vortex.

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

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of points.
        Each row contains the x, y, and z float coordinates of that point's position
        in meters.
    :param origins: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's Slvp_G_Cg's position in meters.
    :param terminations: 2D array of floats
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
    num_vortices = origins.shape[0]
    num_points = points.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities.
    velocities = np.zeros((num_points, num_vortices, 3))

    # If the user didn't specify any ages, set the age of each vortex to 0.0 seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        origin = origins[vortex_id]
        termination = terminations[vortex_id]
        strength = strengths[vortex_id]
        age = ages[vortex_id]

        # Calculate the radius of the vortex's core. If the age is 0.0 seconds,
        # this will evaluate to be 0.0 meters.
        r_c = 2 * math.sqrt(lamb * (nu + squire * abs(strength)) * age)

        # The r_0 vector_G goes from the line vortex's Slvp_G_Cg to its Elvp_G_Cg.
        r_0_x = termination[0] - origin[0]
        r_0_y = termination[1] - origin[1]
        r_0_z = termination[2] - origin[2]

        # Find the r_0 vector_G's length.
        r_0 = math.sqrt(r_0_x**2 + r_0_y**2 + r_0_z**2)

        c_1 = strength / (4 * math.pi)
        c_2 = r_0**2 * r_c**2

        for point_id in range(num_points):
            point = points[point_id]

            # The r_1 vector_G goes from the point to the line vortex's Slvp_G_Cg.
            r_1_x = origin[0] - point[0]
            r_1_y = origin[1] - point[1]
            r_1_z = origin[2] - point[2]

            # The r_2 vector_G goes from the point to the line vortex's Elvp_G_Cg.
            r_2_x = termination[0] - point[0]
            r_2_y = termination[1] - point[1]
            r_2_z = termination[2] - point[2]

            # The r_3 vector_G is the cross product of the r_1 and r_2 vectors.
            r_3_x = r_1_y * r_2_z - r_1_z * r_2_y
            r_3_y = r_1_z * r_2_x - r_1_x * r_2_z
            r_3_z = r_1_x * r_2_y - r_1_y * r_2_x

            # Find the r_1, r_2, and r_3 vectors' lengths.
            r_1 = math.sqrt(r_1_x**2 + r_1_y**2 + r_1_z**2)
            r_2 = math.sqrt(r_2_x**2 + r_2_y**2 + r_2_z**2)
            r_3 = math.sqrt(r_3_x**2 + r_3_y**2 + r_3_z**2)

            c_3 = r_1_x * r_2_x + r_1_y * r_2_y + r_1_z * r_2_z

            # If part of the vortex is so close to the point that they are touching (
            # within machine epsilon), there is a removable discontinuity. In this
            # case, set the velocity components to their true values, which are 0.0
            # meters per second.
            if r_1 < eps or r_2 < eps or r_3**2 < eps:
                continue
            else:
                c_4 = (
                    c_1 * (r_1 + r_2) * (r_1 * r_2 - c_3) / (r_1 * r_2 * (r_3**2 + c_2))
                )
                velocities[point_id, vortex_id, 0] = c_4 * r_3_x
                velocities[point_id, vortex_id, 1] = c_4 * r_3_y
                velocities[point_id, vortex_id, 2] = c_4 * r_3_z

    return velocities
