"""This module contains vortex class definitions.

This module contains the following classes:
    LineVortex: This class is used to contain line vortices.
    HorseshoeVortex: This class is used to contain horseshoe vortices.
    RingVortex: This class is used to contain ring vortices.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


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

        This function uses methodology described on pp. 251-255 of the second edition of "Low-Speed Aerodynamics" by
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
        self.center = asmvp.geometry.centroid_of_quadrilateral(
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
        self.center = asmvp.geometry.centroid_of_quadrilateral(
            self.front_left_vertex,
            self.front_right_vertex,
            self.back_left_vertex,
            self.back_right_vertex,
        )
