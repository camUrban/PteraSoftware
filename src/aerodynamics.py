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
    it finds the cumulative induced velocity due to all of the horseshoe vortices.

    expanded_velocities_from_horseshoe_vortices: This function takes in a group of
    points, and the attributes of a group of horseshoe vortices. At every point,
    it finds the induced velocity due to each horseshoe vortex.

    collapsed_velocities_from_ring_vortices: This function takes in a group of
    points, and the attributes of a group of ring vortices. At every point, it finds
    the cumulative induced velocity due to all of the ring vortices.

    expanded_velocities_from_ring_vortices: This function takes in a group of points,
    and the attributes of a group of ring vortices. At every point, it finds the
    induced velocity due to each ring vortex.

    collapsed_velocities_from_line_vortices: This function takes in a group of
    points, and the attributes of a group of line vortices. At every point, it finds
    the cumulative induced velocity due to all of the line vortices.

    expanded_velocities_from_line_vortices: This function takes in a group of points,
    and the attributes of a group of line vortices. At every point, it finds the
    induced velocity due to each line vortex.
"""
import math

import numpy as np
from numba import njit, prange

from . import functions


# Set the value of Squire's parameter that will be used by the induced velocity
# functions. Squire's parameter relates to the size of the vortex cores and the rate
# at which they grow. The value of this parameter is slightly controversial. It
# dramatically affect the stability of the result. I'm using this value, as cited for
# use in flapping-wing vehicles in "Role of Filament Strain in the Free-Vortex
# Modeling of Rotor Wakes" (Ananthan and Leishman, 2004). It is unitless.
squire = 10 ** -4

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

    def __init__(self, origin, termination, strength):
        """This is the initialization method.

        :param origin: 1D array
            This is a vector containing the x, y, and z coordinates of the origin of
            the line vortex. It's a (3,) array. Its units are meters.
        :param termination: 1D array
            This is a vector containing the x, y, and z coordinates of the
            termination of the line vortex. It's a (3,) array. Its units are meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        """
        self.origin = origin
        self.termination = termination
        self.strength = strength

        # Initialize variables to hold the vector from the vortex's origin to
        # termination, and the point halfway between the origin and termination.
        self.vector = self.termination - self.origin
        self.center = self.origin + 0.5 * self.vector


class HorseshoeVortex:
    """This class is used to contain horseshoe vortices.

    This class contains the following public methods:
        update_strength: This method updates the strength of this horseshoe vortex
        object, and the strength of its legs line vortex objects.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        finite_leg_origin,
        finite_leg_termination,
        strength,
        infinite_leg_direction,
        infinite_leg_length,
    ):
        """This is the initialization method.

        :param finite_leg_origin: 1D array
            This is a vector containing the x, y, and z coordinates of the origin of
            the vortex's finite leg. It's a (,3) array. It's units are meters.
        :param finite_leg_termination: 1D array
            This is a vector containing the x, y, and z coordinates of the
            termination of the vortex's finite leg. It's a (,3) array. It's units are
            meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        :param infinite_leg_direction: 1D array
            This is a unit vector containing the direction that the infinite legs
            extend towards. It's a (,3) array. It's units are meters. It's default
            value is the unit vector in the positive x direction.
        :param infinite_leg_length: float
            This is the length back to extend the quasi-infinite legs of the
            horseshoe vortex. It's units are meters.
        """
        self.finite_leg_origin = finite_leg_origin
        self.finite_leg_termination = finite_leg_termination
        self.strength = strength
        self.infinite_leg_direction = infinite_leg_direction
        self.infinite_leg_length = infinite_leg_length
        self.right_leg_origin = (
            self.finite_leg_origin + infinite_leg_direction * infinite_leg_length
        )
        self.left_leg_termination = (
            self.finite_leg_termination + infinite_leg_direction * infinite_leg_length
        )

        # Initialize a line vortex to represent the horseshoe's finite leg.
        self.right_leg = LineVortex(
            origin=self.right_leg_origin,
            termination=self.finite_leg_origin,
            strength=self.strength,
        )
        self.finite_leg = LineVortex(
            origin=self.finite_leg_origin,
            termination=self.finite_leg_termination,
            strength=self.strength,
        )
        self.left_leg = LineVortex(
            origin=self.finite_leg_termination,
            termination=self.left_leg_termination,
            strength=self.strength,
        )

    def update_strength(self, strength):
        """This method updates the strength of this horseshoe vortex object, and the
        strength of its legs line vortex objects.

        :param strength: float
            This is the strength of this vortex, and of its line vortex legs. Its
            units are meters squared per second.
        :return: None
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

        :param front_left_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            front left point. It's a (,3) array with units of meters.
        :param front_right_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            front right point. It's a (,3) array with units of meters.
        :param back_left_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            back left point. It's a (,3) array with units of meters.
        :param back_right_vertex: 1D array
            This is a vector containing the x, y, and z coordinates of the vortex's
            back right point. It's a (,3) array with units of meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
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
        """This method updates the strength of this ring vortex object, and the
        strength of its four legs' line vortex objects.

        :param strength: float
            This is the strength of this vortex, and of its four legs' line vortices.
            Its units are meters squared per second.
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
        """This method updates the position of the ring vortex, and the positions of
        all its attributes.

        :param front_left_vertex: 1D array
            This is the new position of the x, y, and z coordinates of the front left
            vertex. It is a (,3) array with units of meters.
        :param front_right_vertex: 1D array
            This is the new position of the x, y, and z coordinates of the front
            right vertex. It is a (,3) array with units of meters.
        :param back_left_vertex: 1D array
            This is the new position of the x, y, and z coordinates of the back left
            vertex. It is a (,3) array with units of meters.
        :param back_right_vertex: 1D array
            This is the new position of the x, y, and z coordinates of the back right
            vertex. It is a (,3) array with units of meters.
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


@njit(cache=True, fastmath=True, parallel=True)
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
    horseshoe vortices. At every point, it finds the cumulative induced velocity due
    to all of the horseshoe vortices.

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
        velocity at each of the N points due to all of the horseshoe vortices. The
        units are meters per second.
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
    for i in prange(3):
        induced_velocities += collapsed_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=True, parallel=True)
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
    for i in prange(3):
        induced_velocities += expanded_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=True, parallel=True)
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
    """This function takes in a group of points, and the attributes of a group of
    ring vortices. At every point, it finds the cumulative induced velocity due to
    all of the ring vortices.

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
        velocity at each of the N points due to all of the ring vortices. The units
        are meters per second.
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
    for i in prange(4):
        induced_velocities += collapsed_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=True, parallel=True)
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
    """This function takes in a group of points, and the attributes of a group of
    ring vortices. At every point, it finds the induced velocity due to each ring
    vortex.

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
    for i in prange(4):
        induced_velocities += expanded_velocities_from_line_vortices(
            points=points,
            origins=origins_list[i],
            terminations=terminations_list[i],
            strengths=strengths,
            ages=ages,
            nu=nu,
        )
    return induced_velocities


@njit(cache=True, fastmath=True)
def collapsed_velocities_from_line_vortices(
    points,
    origins,
    terminations,
    strengths,
    ages=None,
    nu=0.0,
):
    """This function takes in a group of points, and the attributes of a group of
    line vortices. At every point, it finds the cumulative induced velocity due to
    all of the line vortices.

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
        vortex's origin's position in meters.
    :param terminations: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's termination's position in meters.
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
        velocity at each of the N points due to all of the line vortices. The units
        are meters per second.
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

        # The r_0 vector goes from the line vortex's origin to its termination.
        r_0_x = termination[0] - origin[0]
        r_0_y = termination[1] - origin[1]
        r_0_z = termination[2] - origin[2]

        # Find the r_0 vector's length.
        r_0 = math.sqrt(r_0_x ** 2 + r_0_y ** 2 + r_0_z ** 2)

        c_1 = strength / (4 * math.pi)
        c_2 = r_0 ** 2 * r_c ** 2

        for point_id in range(num_points):
            point = points[point_id]

            # The r_1 vector goes from the point to the line vortex's origin.
            r_1_x = origin[0] - point[0]
            r_1_y = origin[1] - point[1]
            r_1_z = origin[2] - point[2]

            # The r_2 vector goes from the point to the line vortex's termination.
            r_2_x = termination[0] - point[0]
            r_2_y = termination[1] - point[1]
            r_2_z = termination[2] - point[2]

            # The r_3 vector is the cross product of the r_1 and r_2 vectors.
            r_3_x = r_1_y * r_2_z - r_1_z * r_2_y
            r_3_y = r_1_z * r_2_x - r_1_x * r_2_z
            r_3_z = r_1_x * r_2_y - r_1_y * r_2_x

            # Find the r_1, r_2, and r_3 vectors' lengths.
            r_1 = math.sqrt(r_1_x ** 2 + r_1_y ** 2 + r_1_z ** 2)
            r_2 = math.sqrt(r_2_x ** 2 + r_2_y ** 2 + r_2_z ** 2)
            r_3 = math.sqrt(r_3_x ** 2 + r_3_y ** 2 + r_3_z ** 2)

            c_3 = r_1_x * r_2_x + r_1_y * r_2_y + r_1_z * r_2_z

            # If part of the vortex is so close to the point that they are touching (
            # within machine epsilon), there is a removable discontinuity. In this
            # case, continue to the next point because there is no velocity induced
            # by the current vortex at this point.
            if r_1 < eps or r_2 < eps or r_3 ** 2 < eps:
                continue
            else:
                c_4 = (
                    c_1
                    * (r_1 + r_2)
                    * (r_1 * r_2 - c_3)
                    / (r_1 * r_2 * (r_3 ** 2 + c_2))
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
    """This function takes in a group of points, and the attributes of a group of
    line vortices. At every point, it finds the induced velocity due to each line
    vortex.

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
        vortex's origin's position in meters.
    :param terminations: 2D array of floats
        This variable is an array of shape (M x 3), where M is the number of line
        vortices. Each row contains the x, y, and z float coordinates of that line
        vortex's termination's position in meters.
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
    velocities = np.empty((num_points, num_vortices, 3))

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

        # The r_0 vector goes from the line vortex's origin to its termination.
        r_0_x = termination[0] - origin[0]
        r_0_y = termination[1] - origin[1]
        r_0_z = termination[2] - origin[2]

        # Find the r_0 vector's length.
        r_0 = math.sqrt(r_0_x ** 2 + r_0_y ** 2 + r_0_z ** 2)

        c_1 = strength / (4 * math.pi)
        c_2 = r_0 ** 2 * r_c ** 2

        for point_id in range(num_points):
            point = points[point_id]

            # The r_1 vector goes from the point to the line vortex's origin.
            r_1_x = origin[0] - point[0]
            r_1_y = origin[1] - point[1]
            r_1_z = origin[2] - point[2]

            # The r_2 vector goes from the point to the line vortex's termination.
            r_2_x = termination[0] - point[0]
            r_2_y = termination[1] - point[1]
            r_2_z = termination[2] - point[2]

            # The r_3 vector is the cross product of the r_1 and r_2 vectors.
            r_3_x = r_1_y * r_2_z - r_1_z * r_2_y
            r_3_y = r_1_z * r_2_x - r_1_x * r_2_z
            r_3_z = r_1_x * r_2_y - r_1_y * r_2_x

            # Find the r_1, r_2, and r_3 vectors' lengths.
            r_1 = math.sqrt(r_1_x ** 2 + r_1_y ** 2 + r_1_z ** 2)
            r_2 = math.sqrt(r_2_x ** 2 + r_2_y ** 2 + r_2_z ** 2)
            r_3 = math.sqrt(r_3_x ** 2 + r_3_y ** 2 + r_3_z ** 2)

            c_3 = r_1_x * r_2_x + r_1_y * r_2_y + r_1_z * r_2_z

            # If part of the vortex is so close to the point that they are touching (
            # within machine epsilon), there is a removable discontinuity. In this
            # case, set the velocity components to their true values, which are 0.0
            # meters per second.
            if r_1 < eps or r_2 < eps or r_3 ** 2 < eps:
                velocities[point_id, vortex_id, 0] = 0.0
                velocities[point_id, vortex_id, 1] = 0.0
                velocities[point_id, vortex_id, 2] = 0.0
            else:
                c_4 = (
                    c_1
                    * (r_1 + r_2)
                    * (r_1 * r_2 - c_3)
                    / (r_1 * r_2 * (r_3 ** 2 + c_2))
                )
                velocities[point_id, vortex_id, 0] = c_4 * r_3_x
                velocities[point_id, vortex_id, 1] = c_4 * r_3_y
                velocities[point_id, vortex_id, 2] = c_4 * r_3_z

    return velocities
