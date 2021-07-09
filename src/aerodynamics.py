"""This module contains vortex class definitions, and useful aerodynamic functions.

This module contains the following classes:
    LineVortex: This class is used to contain line vortices.

    HorseshoeVortex: This class is used to contain horseshoe vortices.

    RingVortex: This class is used to contain ring vortices.

This module contains the following exceptions:
    None

This module contains the following functions:
    calculate_velocity_induced_by_horseshoe_vortices: This function takes in a group of
    points, and the attributes of a group of horseshoe vortices. At every point,
    it finds the induced velocity due to every horseshoe vortex, which are
    characterized by groups of back right vertices, front right vertices, front left
    vertices, back left vertices, and strengths.

    calculate_velocity_induced_by_ring_vortices: This function takes in a group of
    points, and the attributes of a group of ring vortices. At every point, it finds
    the induced velocity due to every ring vortex, which are characterized by groups
    of back right vertices, front right vertices, front left vertices, back left
    vertices, and strengths.

    calculate_velocity_induced_by_line_vortices: This function takes in a group of
    points, and the attributes of a group of line vortices. At every point, it finds
    the induced velocity due to every line vortex, which are characterized by groups
    of origins, terminations, and strengths.

    numba_subtract: This function takes in two arrays, one containing N vectors of 3
    components, and one containing M vectors of 3 components. The function then
    calculates and returns a matrix containing the vectors that go from each of the N
    vectors to each of the M vectors.

    numba_explicit_norm: This function takes in an array of shape (M x N x 3),
    which represents M x N vectors, each having 3 components. The function calculates
    and returns the lengths of each vector.
    
    numba_2d_explicit_cross: This function takes in two arrays, which each contain N 
    x M vectors of 3 components. The function then calculates and returns the cross 
    product of the two vectors at each position.

    numba_vector_absolute_magnitude: This function takes in an array of floats which
    is of size (M x N x 3). This represents N x M vectors, each with 3 components.
    The function then calculates and returns the vector absolute magnitude, which is
    the sum of the square of each component.

    numba_discard_singularities: This function takes in an 3d array and replaces any
    inf or nan values with float 0s.

    numba_compute_k_constants: This is a helper function which computes the k
    constants used by the calculate_velocity_induced_by_line_vortices function.

    numba_compute_dot_constants: This function computes the constants used by the
    calculate_velocity_induced_by_line_vortices function that are related to the dot
    product.

    numba_collapse: This function is a helper for the
    calculate_velocity_induced_by_line_vortices function. It uses Numba to speed up
    the summed effects from each of the M line vortices on each of the N points.
"""
import numpy as np
from numba import njit, prange

from . import functions


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


def calculate_velocity_induced_by_horseshoe_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    collapse=True,
):
    """This function takes in a group of points, and the attributes of a group of
    horseshoe vortices. At every point, it finds the induced velocity due to every
    horseshoe vortex, which are characterized by groups of back right vertices,
    front right vertices, front left vertices, back left vertices, and strengths.

    Note: This method uses vectorization, and therefore is much faster for batch
    operations than using the vortex objects' class methods for calculating induced
    velocity.

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
    :param collapse: bool, optional
        This variable determines whether or not the user would like the output to be
        of shape (N x M x 3) or of shape (N x 3). If true, than the effect from every
        horseshoe vortex on a given point will be summed, so the result will be of
        shape (N x 3), where each row identifies the summed effects on a point. If
        false, than the effect from every horseshoe vortex will remain distinct,
        and the shape will be (N x M x 3), where each row/column pair identifies the
        effect on a point by one of the horseshoe vortices.
    :return induced_velocities: either a 2D array of floats or a 3D array of floats
        If collapse is true, the output is the summed effects from every horseshoe
        vortex on a given point. The result will be of shape (N x 3), where each row
        identifies the effects on a point. If false, than the effect from every
        horseshoe vortex will remain distinct, and the shape will be (N x M x 3),
        where each row/column pair identifies the effect on one point by one of the
        horseshoe vortices. Either way, the results units are meters per second.
    """

    # Get the velocity induced by each leg of the horseshoe vortex.
    right_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=back_right_vortex_vertices,
        terminations=front_right_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    finite_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_right_vortex_vertices,
        terminations=front_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    left_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_left_vortex_vertices,
        terminations=back_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )

    # Calculate the total induced velocity by summing the velocities induced by each
    # leg.
    induced_velocities = (
        right_leg_velocities + finite_leg_velocities + left_leg_velocities
    )

    # Return the induced velocity.
    return induced_velocities


def calculate_velocity_induced_by_ring_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    collapse=True,
):
    """This function takes in a group of points, and the attributes of a group of
    ring vortices. At every point, it finds the induced velocity due to every ring
    vortex, which are characterized by groups of back right vertices, front right
    vertices, front left vertices, back left vertices, and strengths.

    Note: This method uses vectorization, and therefore is much faster for batch
    operations than using the vortex objects' class methods for calculating induced
    velocity.

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
    :param collapse: bool, optional
        This variable determines whether or not the user would like the output to be
        of shape (N x M x 3) or of shape (N x 3). If true, than the effect from every
        ring vortex on a given point will be summed, so the result will be of shape (
        N x 3), where each row identifies the summed effects on a point. If false,
        than the effect from every ring vortex will remain distinct, and the shape
        will be (N x M x 3), where each row/column pair identifies the effect on a
        point by one of the ring vortices.
    :return induced_velocities: either a 2D array of floats or a 3D array of floats
        be of shape (N x 3), where each row identifies the effects on a point. If
        false, than the effect from every ring vortex will remain distinct, and the
        shape will be (N x M x 3), where each row/column pair identifies the effect
        on one point by one of the ring vortices. Either way, the results units are
        meters per second.
    """

    # Get the velocity induced by each leg of the ring vortex.
    right_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=back_right_vortex_vertices,
        terminations=front_right_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    front_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_right_vortex_vertices,
        terminations=front_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    left_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_left_vortex_vertices,
        terminations=back_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    back_leg_velocities = calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=back_left_vortex_vertices,
        terminations=back_right_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )

    # Calculate the total induced velocity by summing the velocities induced by each
    # leg.
    induced_velocities = (
        right_leg_velocities
        + front_leg_velocities
        + left_leg_velocities
        + back_leg_velocities
    )

    # Return the induced velocity.
    return induced_velocities


def calculate_velocity_induced_by_line_vortices(
    points, origins, terminations, strengths, collapse=True
):
    """This function takes in a group of points, and the attributes of a group of
    line vortices. At every point, it finds the induced velocity due to every line
    vortex, which are characterized by groups of origins, terminations, and strengths.

    Note: This function uses methodology described on pp. 251-255 of the second
    edition of "Low-Speed Aerodynamics" by Joseph Katz and Allen Plotkin.

    Citation: Some of the following code was adapted by Jérôme Richard as a response
    to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param points: 2D array of floats
        This variable is an array of shape (N x 3), where N is the number of
        points. Each row contains the x, y, and z float coordinates of that point's
        position in meters.
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
    :param collapse: bool, optional
        This variable determines whether or not the user would like the output to be
        of shape (N x M x 3) or of shape (N x 3). If true, than the effect from every
        line vortex on a given point will be summed, so the result will be of shape
        (N x 3), where each row identifies the summed effects on a point. If false,
        than the effect from every line vortex will remain distinct, and the shape will
        be (N x M x 3), where each row/column pair identifies the effect on a point by
        one of the line vortices.
    :return induced_velocities: either a 2D array of floats or a 3D array of floats
        If collapse is true, the output is the summed effects from every line vortex
        on a given point. The result will be of shape (N x 3), where each row
        identifies the effects on a point. If false, than the effect from every line
        vortex will remain distinct, and the shape will be (N x M x 3), where each
        row/column pair identifies the effect on one point by one of the line
        vortices. Either way, the results units are meters per second.
    """
    r_1 = numba_subtract(points, origins)
    r_2 = numba_subtract(points, terminations)

    r_3 = numba_2d_explicit_cross(r_1, r_2)

    r_3_abs_mag = numba_vector_absolute_magnitude(r_3)

    r_1_len = numba_explicit_norm(r_1)
    r_2_len = numba_explicit_norm(r_2)

    radius = 3.0e-16
    r_1_len[r_1_len < radius] = 0
    r_2_len[r_2_len < radius] = 0
    r_3_abs_mag[r_3_abs_mag < radius] = 0

    dot_const_1, dot_const_2 = numba_compute_dot_constants(r_1, r_2)

    with np.errstate(divide="ignore", invalid="ignore"):
        k = numba_compute_k_constants(
            strengths,
            r_3_abs_mag,
            dot_const_1,
            r_1_len,
            dot_const_2,
            r_2_len,
        )
        k = np.expand_dims(k, axis=2)
        induced_velocities = k * r_3

    numba_discard_singularities(induced_velocities)

    if collapse:
        induced_velocities = numba_collapse(induced_velocities)

    return induced_velocities


@njit(parallel=True, cache=True)
def numba_subtract(vectors_1, vectors_2):
    """This function takes in two arrays, one containing N vectors of 3 components,
    and one containing M vectors of 3 components. The function then calculates and
    returns a matrix containing the vectors that go from each of the N vectors to
    each of the M vectors.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param vectors_1: array of floats of size (N x 3)
        This is the first array of N vectors.
    :param vectors_2: array of floats of size (M x 3)
        This is the second array of M vectors.
    :return crosses: array of floats of size (N x M x 3)
        This is array of vectors going between each of the N vectors to each of the M
        vectors.
    """
    diffs = np.empty((vectors_1.shape[0], vectors_2.shape[0], 3))
    for i in prange(diffs.shape[0]):
        for j in range(diffs.shape[1]):
            for k in range(3):
                diffs[i, j, k] = vectors_1[i, k] - vectors_2[j, k]
    return diffs


@njit(parallel=True, cache=True)
def numba_explicit_norm(vectors):
    """This function takes in an array of shape (M x N x 3), which represents M x N
    vectors, each having 3 components. The function calculates and returns the
    lengths of each vector.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param vectors: array of floats of shape (M x N x 3)
        This is the array of vectors.
    :return norms: array of floats of shape (M x N)
        This is the array of vector lengths.
    """
    norms = np.empty((vectors.shape[0], vectors.shape[1]))
    for i in prange(norms.shape[0]):
        for j in range(norms.shape[1]):
            norms[i, j] = np.sqrt(
                vectors[i, j, 0] ** 2 + vectors[i, j, 1] ** 2 + vectors[i, j, 2] ** 2
            )
    return norms


@njit(parallel=True, cache=True)
def numba_2d_explicit_cross(vectors_1, vectors_2):
    """This function takes in two arrays, which each contain N x M vectors of 3
    components. The function then calculates and returns the cross product of the two
    vectors at each position.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param vectors_1: array of floats of size (N x M x 3)
        This is the first array of N x M vectors.
    :param vectors_2: array of floats of size (N x M x 3)
        This is the second array of N x M vectors.
    :return crosses: array of floats of size (N x M x 3)
        This is the cross product of the two inputted vectors at each of the N x M
        positions.
    """
    crosses = np.empty(vectors_1.shape)
    for i in prange(crosses.shape[0]):
        for j in range(crosses.shape[1]):
            crosses[i, j, 0] = (
                vectors_1[i, j, 1] * vectors_2[i, j, 2]
                - vectors_1[i, j, 2] * vectors_2[i, j, 1]
            )
            crosses[i, j, 1] = (
                vectors_1[i, j, 2] * vectors_2[i, j, 0]
                - vectors_1[i, j, 0] * vectors_2[i, j, 2]
            )
            crosses[i, j, 2] = (
                vectors_1[i, j, 0] * vectors_2[i, j, 1]
                - vectors_1[i, j, 1] * vectors_2[i, j, 0]
            )
    return crosses


@njit(parallel=True, cache=True)
def numba_vector_absolute_magnitude(crosses):
    """This function takes in an array of floats which is of size (M x N x 3). This
    represents N x M vectors, each with 3 components. The function then calculates
    and returns the vector absolute magnitude, which is the sum of the square of each
    component.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param crosses: 3d array of floats
        This is the (M x N x 3) array, which holds the M x N vectors.
    :return: 2d array of floats
        This is the (M x N) array, which hold each vector's vector absolute magnitude.
    """
    return crosses[:, :, 0] ** 2 + crosses[:, :, 1] ** 2 + crosses[:, :, 2] ** 2


@njit(parallel=True, cache=True)
def numba_discard_singularities(induced_velocities):
    """This function takes in an 3d array and replaces any inf or nan values with
    float 0s.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param induced_velocities: 3d array of floats
        This is the array of values to sanitize.
    :return: None
    """
    for i in prange(induced_velocities.shape[0]):
        for j in range(induced_velocities.shape[1]):
            for k in range(3):
                if np.isinf(induced_velocities[i, j, k]) or np.isnan(
                    induced_velocities[i, j, k]
                ):
                    induced_velocities[i, j, k] = 0.0


@njit(parallel=True, cache=True)
def numba_compute_k_constants(
    strengths,
    r_3_abs_mag,
    dot_const_1,
    r_1_len,
    dot_const_2,
    r_2_len,
):
    """This is a helper function which computes the k constants used by the
    calculate_velocity_induced_by_line_vortices function.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param strengths: 1D array of floats
        This variable is an array of shape (, M), where M is the number of line
        vortices. Each position contains the strength of that line vortex in meters
        squared per second.
    :param r_3_abs_mag: 2D array of floats
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :param dot_const_1: 2D array of floats
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :param r_1_len: 2D array of floats
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :param dot_const_2: 2D array of floats
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :param r_2_len: 2D array of floats
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :return: 2D array of floats
        This is the array containing the k constants which will be used by the
        calculate_velocity_induced_by_line_vortices function.
    """
    return (
        strengths
        / (4 * np.pi * r_3_abs_mag)
        * (dot_const_1 / r_1_len - dot_const_2 / r_2_len)
    )


@njit(parallel=True, cache=True)
def numba_compute_dot_constants(r_1, r_2):
    """This function computes the constants used by the
    calculate_velocity_induced_by_line_vortices function that are related to the dot
    product.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param r_1: array of floats of shape (N x M x 3)
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :param r_2: array of floats of shape (N x M x 3)
        This variable is an array of shape (N x M), where N is the number of points
        and M is the number of line vortices. The value at each position is
        calculated by the calculate_velocity_induced_by_line_vortices function and is
        related to the relative positions of the points and vortices.
    :return: two arrays of floats, each of shape (N x M)
        These two arrays contain constants used by the
        calculate_velocity_induced_by_line_vortices function that are related to the
        dot product.
    """
    n, m = r_1.shape[0], r_1.shape[1]
    dot_const_1 = np.empty((n, m))
    dot_const_2 = np.empty((n, m))
    for i in prange(n):
        for j in range(m):
            dot_const_1[i, j] = 0.0
            dot_const_2[i, j] = 0.0
            for k in range(3):
                r_0 = r_1[i, j, k] - r_2[i, j, k]
                dot_const_1[i, j] += r_0 * r_1[i, j, k]
                dot_const_2[i, j] += r_0 * r_2[i, j, k]
    return dot_const_1, dot_const_2


@njit(parallel=True, cache=True)
def numba_collapse(induced_velocities):
    """This function is a helper for the calculate_velocity_induced_by_line_vortices
    function. It uses Numba to speed up the summed effects from each of the M line
    vortices on each of the N points.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param induced_velocities: 3D array of floats
        This matrix contains the velocity induced by every line vortex at every
        point. The units are meters per second.
    :return collapses: 2D array of floats
        This matrix is the total induced velocity at each point due to all of the
        line vortices. The units are meters per second.
    """
    n, m = induced_velocities.shape[0], induced_velocities.shape[1]
    collapses = np.empty((n, 3))
    for i in prange(n):
        collapses[i, 0] = np.sum(induced_velocities[i, :, 0])
        collapses[i, 1] = np.sum(induced_velocities[i, :, 1])
        collapses[i, 2] = np.sum(induced_velocities[i, :, 2])
    return collapses
