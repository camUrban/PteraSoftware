"""This module contains vortex class definitions, and useful aerodynamic functions.

This module contains the following classes:
    LineVortex: This class is used to contain line vortices.
    HorseshoeVortex: This class is used to contain horseshoe vortices.
    RingVortex: This class is used to contain ring vortices.

This module contains the following exceptions:
    None

This module contains the following functions:
    calculate_velocity_induced_by_line_vortices: This function takes in a group of points, origins, terminations and
                                                 strengths. At every point, it finds the induced velocity due to every
                                                 line vortex, which are characterized by the groups of origins,
                                                 terminations, and strengths.
    calculate_velocity_induced_by_horseshoe_vortices: This function takes in a group of points, and the attributes of a
                                                      group of horseshoe vortices. At every point, it finds the induced
                                                      velocity due to every horseshoe vortex, which are characterized by
                                                      groups of back right vertices, front right vertices, front left
                                                      vertices, back left vertices, and strengths.
    calculate_velocity_induced_by_ring_vortices: This function takes in a group of points, and the attributes of a group
                                                 of ring vortices. At every point, it finds the induced velocity due to
                                                 every ring vortex, which are characterized by groups of back right
                                                 vertices, front right vertices, front left vertices, back left
                                                 vertices, and strengths.
"""

import numpy as np

import pterasoftware as ps


class LineVortex:
    """This class is used to contain line vortices.

    This class contains the following public methods:
        calculate_normalized_induced_velocity: This method calculates the velocity induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_induced_velocity: This method calculates the velocity induced at a point by this vortex with its given
                                    vortex strength.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, origin, termination, strength):
        """This is the initialization method.

        :param origin: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the origin of the line vortex. It's a (3,)
            ndarray. Its units are meters.
        :param termination: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the termination of the line vortex. It's a (3,)
            ndarray. Its units are meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        """

        self.origin = origin
        self.termination = termination
        self.strength = strength

        # Initialize variables to hold the vector from the vortex's origin to termination, and the point halfway between
        # the origin and termination.
        self.vector = self.termination - self.origin
        self.center = self.origin + 0.5 * self.vector

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with a unit vortex strength.

        :param point: 1D ndarray
            This parameter is the point where the induced velocity is to be calculated. It's a (3,) ndarray. Its units
            are meters.
        :return normalized_induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity. It's a (3,) ndarray. Its
            units are meters per second.
        """

        normalized_induced_velocity = self.calculate_induced_velocity(
            point=point, overriding_strength=1
        )

        # Return the normalized induced velocity.
        return normalized_induced_velocity

    def calculate_induced_velocity(self, point, overriding_strength=None):
        """This method calculates the velocity induced at a point by this vortex with its given vortex strength.

        This method uses methodology described on pp. 251-255 of the second edition of "Low-Speed Aerodynamics" by
        Joseph Katz and Allen Plotkin.

        :param point: 1D ndarray
            This parameter is the point where the induced velocity is to be calculated. It's a (3,) ndarray. Its units
            are meters.
        :param overriding_strength: float
            This is the magnitude of the vorticity. It's sign is given by the using the right hand rule on the vector
            from the line vortex's origin to termination. Its units are meters squared per second. It's default value is
            None. If None, then this method will use the vortex's assigned strength. Otherwise, it will use this
            strength.
        :return induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity. It's a (3,) ndarray. Its
            units are meters per second.
        """

        # If the overriding strength is not None, then calculate the induced velocity using the overriding value.
        # Otherwise, use the object's strength attribute's value.
        if overriding_strength is None:
            strength = self.strength
        else:
            strength = overriding_strength

        # Define the vectors from the vortex to the point.
        r_1 = point - self.origin
        r_2 = point - self.termination

        # Define the vector from the vortex origin to the vortex termination.
        r_0 = r_1 - r_2

        # Calculate the vector cross product.
        r_1_cross_r_2 = np.cross(r_1, r_2)

        # Calculate the cross product's absolute magnitude.
        r_1_cross_r_2_absolute_magnitude = (
            r_1_cross_r_2[0] ** 2 + r_1_cross_r_2[1] ** 2 + r_1_cross_r_2[2] ** 2
        )

        # Calculate the vector lengths.
        r_1_length = np.linalg.norm(r_1)
        r_2_length = np.linalg.norm(r_2)

        # Check for singularities.
        line_vortex_radius = 3.0e-16
        if (
            r_1_length < line_vortex_radius
            or r_2_length < line_vortex_radius
            or r_1_cross_r_2_absolute_magnitude < line_vortex_radius
        ):
            # If there is a singularity, the induced velocity is zero.
            return np.array([0, 0, 0])

        # Calculate the vector dot products.
        r_0_dot_r_1 = np.dot(r_0, r_1)
        r_0_dot_r_2 = np.dot(r_0, r_2)

        # Calculate the k coefficient.
        k = (
            strength
            / (4 * np.pi * r_1_cross_r_2_absolute_magnitude)
            * (r_0_dot_r_1 / r_1_length - r_0_dot_r_2 / r_2_length)
        )

        # Calculate the induced velocity components, and combine them into the induced velocity ndarray.
        u = k * r_1_cross_r_2[0]
        v = k * r_1_cross_r_2[1]
        w = k * r_1_cross_r_2[2]
        return np.array([u, v, w])


class HorseshoeVortex:
    """This class is used to contain horseshoe vortices.

    This class contains the following public methods:
        calculate_normalized_induced_velocity: This method calculates the velocity induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_induced_velocity: This method calculates the velocity induced at a point by this vortex with its given
                                    vortex strength.
        update_strength: This method updates the strength of this horseshoe vortex object, and the strength of its legs
                         line vortex objects.

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

        :param finite_leg_origin: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the origin of the vortex's finite leg. It's a
            (,3) ndarray. It's units are meters.
        :param finite_leg_termination: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the termination of the vortex's finite leg. It's
            a (,3) ndarray. It's units are meters.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        :param infinite_leg_direction: 1D ndarray
            This is a unit vector containing the direction that the infinite legs extend towards. It's a (,3) ndarray.
            It's units are meters. It's default value is the unit vector in the positive x direction.
        :param infinite_leg_length: float
            This is the length back to extend the quasi-infinite legs of the horseshoe vortex. It's units are meters.
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

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with a unit vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at. It's a
            (,3) ndarray with units of meters.
        :return normalized_induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity. It's a (,3) ndarray with
            units of meters per second.
        """

        normalized_velocity_induced_by_right_leg = self.right_leg.calculate_normalized_induced_velocity(
            point=point
        )
        normalized_velocity_induced_by_finite_leg = self.finite_leg.calculate_normalized_induced_velocity(
            point=point
        )
        normalized_velocity_induced_by_left_leg = self.left_leg.calculate_normalized_induced_velocity(
            point=point
        )

        # Sum the velocities induced by each leg to get the velocity induced by the entire horseshoe.
        normalized_induced_velocity = (
            normalized_velocity_induced_by_right_leg
            + normalized_velocity_induced_by_finite_leg
            + normalized_velocity_induced_by_left_leg
        )

        # Return the velocity induced by the vortex.
        return normalized_induced_velocity

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with its given vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at. It's a
            (,3) ndarray with units of meters.
        :return induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity. It's a (,3) ndarray with
            units of meters per second.
        """

        velocity_induced_by_right_leg = self.right_leg.calculate_induced_velocity(
            point=point
        )
        velocity_induced_by_finite_leg = self.finite_leg.calculate_induced_velocity(
            point=point
        )
        velocity_induced_by_left_leg = self.left_leg.calculate_induced_velocity(
            point=point
        )

        # Sum the velocities induced by each leg to get the velocity induced by the entire horseshoe.
        induced_velocity = (
            velocity_induced_by_right_leg
            + velocity_induced_by_finite_leg
            + velocity_induced_by_left_leg
        )

        # Return the velocity induced by the vortex.
        return induced_velocity

    def update_strength(self, strength):
        """This method updates the strength of this horseshoe vortex object, and the strength of its legs line vortex
        objects.

        :param strength: float
            This is the strength of this vortex, and of its line vortex legs. Its units are meters squared per second.
        :return: None
        """

        self.strength = strength
        self.right_leg.strength = strength
        self.finite_leg.strength = strength
        self.left_leg.strength = strength


class RingVortex:
    """This class is used to contain ring vortices.

    This class contains the following public methods:
        calculate_normalized_induced_velocity: This method calculates the velocity induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_induced_velocity: This method calculates the velocity induced at a point by this vortex with its given
                                    vortex strength.
        update_strength: This method updates the strength of this ring vortex object, and the strength of its
                         four legs' line vortex objects.
        update_position: This method updates the position of the ring vortex, and the positions of all its attributes.

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

        :param front_left_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's front left point. It's a (,3)
            ndarray with units of meters.
        :param front_right_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's front right point. It's a (,3)
            ndarray with units of meters.
        :param back_left_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's back left point. It's a (,3)
            ndarray with units of meters.
        :param back_right_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's back right point. It's a (,3)
            ndarray with units of meters.
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
        self.center = ps.geometry.centroid_of_quadrilateral(
            self.front_left_vertex,
            self.front_right_vertex,
            self.back_left_vertex,
            self.back_right_vertex,
        )

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with a unit vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at. It's a
            (,3) ndarray with units of meters.
        :return normalized_induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity. It's a (,3) ndarray with
            units of meters per second.
        """

        normalized_velocity_induced_by_front_leg = self.front_leg.calculate_normalized_induced_velocity(
            point=point
        )
        normalized_velocity_induced_by_left_leg = self.left_leg.calculate_normalized_induced_velocity(
            point=point
        )
        normalized_velocity_induced_by_back_leg = self.back_leg.calculate_normalized_induced_velocity(
            point=point
        )
        normalized_velocity_induced_by_right_leg = self.right_leg.calculate_normalized_induced_velocity(
            point=point
        )

        # Sum the velocities induced by each leg to get the velocity induced by the entire ring vortex.
        normalized_induced_velocity = (
            normalized_velocity_induced_by_front_leg
            + normalized_velocity_induced_by_left_leg
            + normalized_velocity_induced_by_back_leg
            + normalized_velocity_induced_by_right_leg
        )

        return normalized_induced_velocity

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with its given vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at. It's a
            (,3) ndarray with units of meters.
        :return induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity. It's a (,3) ndarray with
            units of meters per second.
        """

        velocity_induced_by_front_leg = self.front_leg.calculate_induced_velocity(
            point=point
        )
        velocity_induced_by_left_leg = self.left_leg.calculate_induced_velocity(
            point=point
        )
        velocity_induced_by_back_leg = self.back_leg.calculate_induced_velocity(
            point=point
        )
        velocity_induced_by_right_leg = self.right_leg.calculate_induced_velocity(
            point=point
        )

        # Sum the velocities induced by each leg to get the velocity induced by the entire ring vortex.
        induced_velocity = (
            velocity_induced_by_front_leg
            + velocity_induced_by_left_leg
            + velocity_induced_by_back_leg
            + velocity_induced_by_right_leg
        )

        return induced_velocity

    def update_strength(self, strength):
        """This method updates the strength of this ring vortex object, and the strength of its four legs' line vortex
        objects.

        :param strength: float
            This is the strength of this vortex, and of its four legs' line vortices. Its units are meters squared per
            second.
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
        """This method updates the position of the ring vortex, and the positions of all its attributes.

        :param front_left_vertex: 1D ndarray
            This is the new position of the x, y, and z coordinates of the front left vertex. It is a (,3) ndarray with
            units of meters.
        :param front_right_vertex: 1D ndarray
            This is the new position of the x, y, and z coordinates of the front right vertex. It is a (,3) ndarray with
            units of meters.
        :param back_left_vertex: 1D ndarray
            This is the new position of the x, y, and z coordinates of the back left vertex. It is a (,3) ndarray with
            units of meters.
        :param back_right_vertex: 1D ndarray
            This is the new position of the x, y, and z coordinates of the back right vertex. It is a (,3) ndarray with
            units of meters.
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
        self.center = ps.geometry.centroid_of_quadrilateral(
            self.front_left_vertex,
            self.front_right_vertex,
            self.back_left_vertex,
            self.back_right_vertex,
        )


def calculate_velocity_induced_by_line_vortices(
    points, origins, terminations, strengths, collapse=True
):
    """ This function takes in a group of points, origins, terminations and strengths. At every point, it finds the
    induced velocity due to every line vortex, which are characterized by the groups of origins, terminations, and
    strengths.

        This method uses vectorization, and therefore is much faster for batch operations than using the vortex objects'
        class methods for calculating induced velocity.

        This function uses methodology described on pp. 251-255 of the second edition of "Low-Speed Aerodynamics" by
        Joseph Katz and Allen Plotkin.

    :param points: 2D ndarray of floats
        This variable is an ndarray of shape (N x 3), where N is the number of points. Each row contains the x, y, and z
        float coordinates of that point's position in meters.
    :param origins: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of line vortices. Each row contains the x,
        y, and z float coordinates of that line vortex's origin's position in meters.
    :param terminations: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of line vortices. Each row contains the x,
        y, and z float coordinates of that line vortex's termination's position in meters.
    :param strengths: 1D ndarray of floats
        This variable is an ndarray of shape (, M), where M is the number of line vortices. Each position contains the
        strength of that line vortex in meters squared per second.
    :param collapse: bool, optional
        This variable determines whether or not the user would like the output to be of shape (N x M x 3) or of shape
        (N x 3). If true, than the effect from every line vortex on a given point will be summed, so the result will be
        of shape (N x 3), where each row identifies the summed effects on a point. If false, than the effect from every
        line vortex will remain distinct, and the shape will be (N x M x 3), where each row/column pair identifies the
        effect on a point by one of the line vortices.
    :return induced_velocities: either a 2D ndarray of floats or a 3D ndarray of floats
        If collapse is true, the output is the summed effects from every line vortex on a given point. The result will
        be of shape (N x 3), where each row identifies the effects on a point. If false, than the effect from every line
        vortex will remain distinct, and the shape will be (N x M x 3), where each row/column pair identifies the effect
        on one point by one of the line vortices. Either way, the results units are meters per second.
    """

    # Validate the input. These steps significantly reduces the complexity of debugging.

    # Check that all the inputs are ndarrays.
    assert isinstance(points, np.ndarray)
    assert isinstance(origins, np.ndarray)
    assert isinstance(terminations, np.ndarray)
    assert isinstance(strengths, np.ndarray)

    # Get the shapes of the input ndarrays.
    points_shape = points.shape
    origins_shape = origins.shape
    terminations_shape = terminations.shape
    strengths_shape = strengths.shape

    # Get the number of dimensions of the input ndarrays.
    points_num_dims = points.ndim
    origins_num_dims = origins.ndim
    terminations_num_dims = terminations.ndim
    strengths_num_dims = strengths.ndim

    # Validate that each input has the correct dimensionality.
    assert points_num_dims == 2
    assert origins_num_dims == 2
    assert terminations_num_dims == 2
    assert strengths_num_dims == 1

    # Find the number of line vortices, according to the origins, terminations, and strengths inputs.
    num_origins = origins_shape[0]
    num_terminations = terminations_shape[0]
    num_strengths = strengths_shape[0]

    # Check that the origins, terminations, and strengths inputs agree on the number of line vortices.
    assert num_origins == num_terminations == num_strengths

    # Get the number of coordinates for each row in the points, origins, and terminations inputs.
    points_num_coordinates = points_shape[-1]
    origin_num_coordinates = origins_shape[-1]
    termination_num_coordinates = terminations_shape[-1]

    # Check that each row of the points, origins, and terminations inputs all have three coordinates (x, y, and z).
    assert points_num_coordinates == 3
    assert origin_num_coordinates == 3
    assert termination_num_coordinates == 3

    # Check that every position in every input holds a float.
    for coordinate in np.ravel(points):
        assert isinstance(coordinate, float)
    for coordinate in np.ravel(origins):
        assert isinstance(coordinate, float)
    for coordinate in np.ravel(terminations):
        assert isinstance(coordinate, float)
    for strength in np.ravel(strengths):
        assert isinstance(strength, float)

    # We have now verified that:
    #   1.  points is a 2D ndarray of floats with shape Nx3. Each row stores one of N points. The columns store each
    #       point's x, y, and z coordinates.
    #   2.  origins is a 2D ndarray of floats with shape Mx3. Each row stores one of M origins. The columns store each
    #       origin's x, y, and z coordinates.
    #   3.  terminations is a 2D ndarray of floats with shape Mx3. Each row stores one of M terminations. The columns
    #       store each termination's x, y, and z coordinates
    #   4.  strengths is a 1D ndarray of floats with shape M. Each positions stores one of M strengths. The columns
    #       store each strength's value.

    # Get the number of vortices, and the number of points.
    num_vortices = num_strengths
    num_points = points_shape[0]

    # Expand the dimensionality of the points input. It is now of shape (N x 1 x 3). This will allow numpy to
    # broadcast the upcoming subtractions.
    points = np.expand_dims(points, axis=1)

    # Define the vectors from the vortex to the points. r_1 and r_2 now both are of shape (N x M x 3). Each row/column
    # pair holds the vector associated with each point/vortex pair.
    r_1 = points - origins
    r_2 = points - terminations

    # Define the vector from the vortex origins to the vortex terminations. This is of shape (N x M x 3).
    r_0 = r_1 - r_2

    # Calculate the vector cross product. This is of shape (N x M x 3).
    r_1_cross_r_2 = np.cross(r_1, r_2, axis=-1)

    # Calculate the cross product's absolute magnitude. This is of shape (N x M).
    r_1_cross_r_2_absolute_magnitude = (
        r_1_cross_r_2[:, :, 0] ** 2
        + r_1_cross_r_2[:, :, 1] ** 2
        + r_1_cross_r_2[:, :, 2] ** 2
    )

    # Calculate the vector lengths. These are of shape (N x M).
    r_1_length = np.linalg.norm(r_1, axis=-1)
    r_2_length = np.linalg.norm(r_2, axis=-1)

    # Define the radius of the line vortices. This is used to get rid of any singularities.
    radius = 3.0e-16

    # Set the lengths and the absolute magnitudes to zero, at the places where the lengths and absolute magnitudes are
    # less than the vortex radius. This insures that the calculation for the constant k will produce np.inf or np.nan
    # values at the locations where there are singularities.
    r_1_length[r_1_length < radius] = 0
    r_2_length[r_2_length < radius] = 0
    r_1_cross_r_2_absolute_magnitude[r_1_cross_r_2_absolute_magnitude < radius] = 0

    # Calculate the vector dot products. This uses numpy's einsum function for speed.
    r_0_dot_r_1 = np.einsum("ijk,ijk->ij", r_0, r_1)
    r_0_dot_r_2 = np.einsum("ijk,ijk->ij", r_0, r_2)

    # Calculate k and then the induced velocity, ignoring any divide-by-zero or nan errors. k is of shape (N x M)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (
            strengths
            / (4 * np.pi * r_1_cross_r_2_absolute_magnitude)
            * (r_0_dot_r_1 / r_1_length - r_0_dot_r_2 / r_2_length)
        )

        # Set the shape of k to be (N x M x 1) to support numpy broadcasting in the subsequent multiplication.
        k = np.expand_dims(k, axis=2)

        # Multiple k by the cross products of r_1 and r_2 to get the non-collapsed matrix of induced velocities. This is
        # of shape (M x N x 3).
        induced_velocities = k * r_1_cross_r_2

    # Get the shape of the induced velocity matrix.
    induced_velocities_shape = induced_velocities.shape

    # Check that the calculations produced the expected shape.
    assert induced_velocities_shape == (num_points, num_vortices, 3)

    # Set the values of the induced velocity to zero where there are singularities.
    induced_velocities[np.isinf(induced_velocities)] = 0
    induced_velocities[np.isnan(induced_velocities)] = 0

    if collapse:
        induced_velocities = np.sum(induced_velocities, axis=1)
        induced_velocities_shape = induced_velocities.shape
        assert induced_velocities_shape == (num_points, 3)

    return induced_velocities


def calculate_velocity_induced_by_horseshoe_vortices(
    points,
    back_right_vortex_vertices,
    front_right_vortex_vertices,
    front_left_vortex_vertices,
    back_left_vortex_vertices,
    strengths,
    collapse=True,
):
    """ This function takes in a group of points, and the attributes of a group of horseshoe vortices. At every point,
    it finds the induced velocity due to every horseshoe vortex, which are characterized by groups of back right
    vertices, front right vertices, front left vertices, back left vertices, and strengths.

        This method uses vectorization, and therefore is much faster for batch operations than using the vortex objects'
        class methods for calculating induced velocity.

    :param points: 2D ndarray of floats
        This variable is an ndarray of shape (N x 3), where N is the number of points. Each row contains the x, y, and z
        float coordinates of that point's position in meters.
    :param back_right_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of horseshoe vortices. Each row contains the
        x, y, and z float coordinates of that horseshoe vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of horseshoe vortices. Each row contains the
        x, y, and z float coordinates of that horseshoe vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of horseshoe vortices. Each row contains the
        x, y, and z float coordinates of that horseshoe vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of horseshoe vortices. Each row contains the
        x, y, and z float coordinates of that horseshoe vortex's front left vertex's position in meters.
    :param strengths: 1D ndarray of floats
        This variable is an ndarray of shape (, M), where M is the number of horseshoe vortices. Each holds the strength
        of that horseshoe vortex in meters squared per second.
    :param collapse: bool, optional
        This variable determines whether or not the user would like the output to be of shape (N x M x 3) or of shape
        (N x 3). If true, than the effect from every horseshoe vortex on a given point will be summed, so the result
        will be of shape (N x 3), where each row identifies the summed effects on a point. If false, than the effect
        from every horseshoe vortex will remain distinct, and the shape will be (N x M x 3), where each row/column pair
        identifies the effect on a point by one of the horseshoe vortices.
    :return induced_velocities: either a 2D ndarray of floats or a 3D ndarray of floats
        If collapse is true, the output is the summed effects from every horseshoe vortex on a given point. The result
        will be of shape (N x 3), where each row identifies the effects on a point. If false, than the effect from every
        horseshoe vortex will remain distinct, and the shape will be (N x M x 3), where each row/column pair identifies
        the effect on one point by one of the horseshoe vortices. Either way, the results units are meters per second.
    """

    # Get the velocity induced by each leg of the horseshoe vortex.
    right_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=back_right_vortex_vertices,
        terminations=front_right_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    finite_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_right_vortex_vertices,
        terminations=front_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    left_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_left_vortex_vertices,
        terminations=back_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )

    # Calculate the total induced velocity by summing the velocities induced by each leg.
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
    """ This function takes in a group of points, and the attributes of a group of ring vortices. At every point, it
    finds the induced velocity due to every ring vortex, which are characterized by groups of back right vertices, front
    right vertices, front left vertices, back left vertices, and strengths.

        This method uses vectorization, and therefore is much faster for batch operations than using the vortex objects'
        class methods for calculating induced velocity.

    :param points: 2D ndarray of floats
        This variable is an ndarray of shape (N x 3), where N is the number of points. Each row contains the x, y, and z
        float coordinates of that point's position in meters.
    :param back_right_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of ring vortices. Each row contains the x,
        y, and z float coordinates of that ring vortex's back right vertex's position in meters.
    :param front_right_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of ring vortices. Each row contains the x,
        y, and z float coordinates of that ring vortex's front right vertex's position in meters.
    :param front_left_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of ring vortices. Each row contains the x,
        y, and z float coordinates of that ring vortex's front left vertex's position in meters.
    :param back_left_vortex_vertices: 2D ndarray of floats
        This variable is an ndarray of shape (M x 3), where M is the number of ring vortices. Each row contains the x,
        y, and z float coordinates of that ring vortex's front left vertex's position in meters.
    :param strengths: 1D ndarray of floats
        This variable is an ndarray of shape (, M), where M is the number of ring vortices. Each holds the strength of
        that ring vortex in meters squared per second.
    :param collapse: bool, optional
        This variable determines whether or not the user would like the output to be of shape (N x M x 3) or of shape
        (N x 3). If true, than the effect from every ring vortex on a given point will be summed, so the result will be
        of shape (N x 3), where each row identifies the summed effects on a point. If false, than the effect from every
        ring vortex will remain distinct, and the shape will be (N x M x 3), where each row/column pair identifies the
        effect on a point by one of the ring vortices.
    :return induced_velocities: either a 2D ndarray of floats or a 3D ndarray of floats
        If collapse is true, the output is the summed effects from every ring vortex on a given point. The result will
        be of shape (N x 3), where each row identifies the effects on a point. If false, than the effect from every ring
        vortex will remain distinct, and the shape will be (N x M x 3), where each row/column pair identifies the effect
        on one point by one of the ring vortices. Either way, the results units are meters per second.
    """

    # Get the velocity induced by each leg of the ring vortex.
    right_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=back_right_vortex_vertices,
        terminations=front_right_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    front_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_right_vortex_vertices,
        terminations=front_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    left_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=front_left_vortex_vertices,
        terminations=back_left_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )
    back_leg_velocities = ps.aerodynamics.calculate_velocity_induced_by_line_vortices(
        points=points,
        origins=back_left_vortex_vertices,
        terminations=back_right_vortex_vertices,
        strengths=strengths,
        collapse=collapse,
    )

    # Calculate the total induced velocity by summing the velocities induced by each leg.
    induced_velocities = (
        right_leg_velocities
        + front_leg_velocities
        + left_leg_velocities
        + back_leg_velocities
    )

    # Return the induced velocity.
    return induced_velocities
