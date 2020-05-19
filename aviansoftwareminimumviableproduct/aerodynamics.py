"""This module contains useful aerodynamics functions, and the vortex class definitions.

This module contains the following classes:
    LineVortex: This class is used to contain line vortices.
    HorseshoeVortex: This class is used to contain horseshoe vortices.
    RingVortex: This class is used to contain ring vortices.

This module contains the following exceptions:
    None

This module contains the following functions:
    induced_velocity_from_line_vortex: Calculates the velocity induced at a specified point by a line vortex of a
                                       specified strength, and at a specified location.
    induced_velocity_from_horseshoe_vortex: Calculates the velocity induced at a specified point by a horseshoe vortex
                                            of a specified strength, and at a specified location.
    induced_velocity_from_ring_vortex: Calculates the velocity induced at a specified point by a ring vortex of a
                                       specified strength, and at a specified location.
"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


def calculate_velocity_induced_from_line_vortex(point, vortex_origin, vortex_termination, vortex_strength):
    """This function calculates the velocity induced at a specified point by a line vortex of a specified strength,
    and at a specified location.

    This function uses methodology described on pp. 251-255 of the second edition of "Low-Speed Aerodynamics" by Joseph
    Katz and Allen Plotkin.

    :param point: This parameter is the point where the induced velocity is to be calculated. It's a (3,) ndarray.
                  Its units are meters.
    :param vortex_origin: This parameter is the point where the line vortex begins. It's a (3,) ndarray. Its units
                          are meters.
    :param vortex_termination: This parameter is the point where the line vortex ends. It's a (3,) ndarray. Its
                               units are meters.
    :param vortex_strength: This is the magnitude of the vorticity. It's sign is given by the using the right hand rule
                            on the vector from the line vortex's origin to termination. Its units are meters squared
                            per second.
    :return velocity_induced_by_vortex: This is the velocity induced at the point by the line vortex. It's a (3,) numpy
                                        array. Its units are meters per second.
    """

    # Define the vectors from the vortex to the point.
    r_1 = point - vortex_origin
    r_2 = point - vortex_termination

    # Define the vector from the vortex origin to the vortex termination.
    r_0 = r_1 - r_2

    # Calculate the vector cross product.
    r_1_cross_r_2 = np.cross(r_1, r_2)

    # Calculate the cross product's absolute value.
    r_1_cross_r_2_absolute_value = r_1_cross_r_2[0] ** 2 + r_1_cross_r_2[1] ** 2 + r_1_cross_r_2[2] ** 2

    # Calculate the vector lengths.
    r_1_length = np.linalg.norm(r_1)
    r_2_length = np.linalg.norm(r_2)

    # Check for singularities.
    line_vortex_radius = 3.0e-16
    if (r_1_length < line_vortex_radius
            or r_2_length < line_vortex_radius
            or r_1_cross_r_2_absolute_value < line_vortex_radius):
        # If there is a singularity, the induced velocity is zero.
        velocity_induced_by_vortex = np.array([0, 0, 0])
    else:
        # Calculate the vector dot products.
        r_0_dot_r_1 = np.dot(r_0, r_1)
        r_0_dot_r_2 = np.dot(r_0, r_2)

        # Calculate the k coefficient.
        k = (vortex_strength / (4 * np.pi * r_1_cross_r_2_absolute_value)
             * (r_0_dot_r_1 / r_1_length - r_0_dot_r_2 / r_2_length))

        # Calculate the induced velocity components, and combine them into the induced velocity ndarray.
        u = k * r_1_cross_r_2[0]
        v = k * r_1_cross_r_2[1]
        w = k * r_1_cross_r_2[2]
        velocity_induced_by_vortex = np.array([u, v, w])

    # Return the induced velocity.
    return velocity_induced_by_vortex


def calculate_velocity_induced_from_horseshoe_vortex(point, finite_leg_origin, finite_leg_termination, vortex_strength):
    """This function calculates the velocity induced at a specified point by a horseshoe vortex of a specified strength,
    and at a specified location.

    This function uses methodology described on pp. 131-132 of "Flight Vehicle Aerodynamics" by Mark Drela.

    Notes:
        1. This function assumes that the horseshoe vortex's infinite legs, which are rays originating from the finite
        leg's origin and termination points, go straight back in the positive x direction.

    :param point: This parameter is the point where the induced velocity is to be calculated. It's a (3,) ndarray.
                  Its units are meters.
    :param finite_leg_origin: This parameter is the point where the horseshoe vortex's finite leg begins. It's a (3,)
                              ndarray. Its units are meters.
    :param finite_leg_termination: This parameter is the point where the horseshoe vortex's finite leg ends. It's a (3,)
                                   ndarray. Its units are meters.
    :param vortex_strength: This is the magnitude of the vorticity. It's sign is given by the using the right hand rule
                            on the vector from the horseshoe vortex's finite leg's origin to termination. Its units are
                            meters squared per second.
    :return velocity_induced_by_vortex: This is the velocity induced at the point by the horseshoe vortex. It's a (3,)
                                        ndarray. Its units are meters per second.
    """
    velocity_induced_by_right_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=finite_leg_origin + np.array([20, 0, 0]),
        vortex_termination=finite_leg_origin,
        vortex_strength=vortex_strength
    )
    velocity_induced_by_front_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=finite_leg_origin,
        vortex_termination=finite_leg_termination,
        vortex_strength=vortex_strength
    )
    velocity_induced_by_left_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=finite_leg_termination,
        vortex_termination=finite_leg_termination + np.array([20, 0, 0]),
        vortex_strength=vortex_strength
    )

    # Sum the velocities induced by each leg to get the velocity induced by the entire horseshoe.
    velocity_induced_by_vortex = (
        velocity_induced_by_right_leg
        + velocity_induced_by_front_leg
        + velocity_induced_by_left_leg
    )

    # Return the velocity induced by the vortex.
    return velocity_induced_by_vortex

    # # Define vectors from the vortex's finite leg to the point.
    # a = point - finite_leg_origin
    # b = point - finite_leg_termination
    #
    # # Calculate the vector lengths.
    # a_length = np.linalg.norm(a)
    # b_length = np.linalg.norm(b)
    #
    # # Define the x direction unit vector.
    # x_unit = np.array([1, 0, 0])
    #
    # # Calculate the vector dot and cross products.
    # a_cross_b = np.cross(a, b)
    # a_dot_b = np.dot(a, b)
    # a_cross_x_unit = np.cross(a, x_unit)
    # b_cross_x_unit = np.cross(b, x_unit)
    # a_dot_x_unit = np.dot(a, x_unit)
    # b_dot_x_unit = np.dot(b, x_unit)
    #
    # # Check for singularities on each of the three legs
    # horseshoe_radius = 3.0e-16
    # if ((abs(a_length * b_length + a_dot_b) < horseshoe_radius)
    #         or (abs(a_length) < horseshoe_radius)
    #         or (abs(b_length) < horseshoe_radius)):
    #     term1 = 0
    # else:
    #     term1 = (a_cross_b / (a_length * b_length + a_dot_b)) * (1 / a_length + 1 / b_length)
    #
    # if ((abs(a_length - a_dot_x_unit) < horseshoe_radius)
    #         or (abs(a_length) < horseshoe_radius)):
    #     term2 = 0
    # else:
    #     term2 = (a_cross_x_unit / (a_length - a_dot_x_unit)) * (1 / a_length)
    #
    # if ((abs(b_length - b_dot_x_unit) < horseshoe_radius)
    #         or (abs(b_length) < horseshoe_radius)):
    #     term3 = 0
    # else:
    #     term3 = (b_cross_x_unit / (b_length - b_dot_x_unit)) * (1 / b_length)
    #
    # velocity_induced_by_vortex = (vortex_strength / (4 * np.pi) * (term1 + term2 - term3))
    #
    # # Return the velocity induced by the vortex.
    # return velocity_induced_by_vortex


def calculate_velocity_induced_from_ring_vortex(point, front_left_vortex_vertex, front_right_vortex_vertex,
                                                back_left_vortex_vertex, back_right_vortex_vertex, vortex_strength):
    """This function calculates the velocity induced at a specified point by a ring vortex of a specified strength,
    and at a specified location.

    This function uses methodology described on pp. 255-256 of the second edition of "Low-Speed
    Aerodynamics" by Joseph Katz and Allen Plotkin.

    :param point: This parameter is the point where the induced velocity is to be calculated. It's a (3,) ndarray.
                  Its units are meters.
    :param front_left_vortex_vertex: This parameter is the ring vortex's front left point. It's a (3,) ndarray. Its
                                     units are meters.
    :param front_right_vortex_vertex: This parameter is the ring vortex's front right point. It's a (3,) ndarray.
                                      Its units are meters.
    :param back_left_vortex_vertex: This parameter is the ring vortex's back left point. It's a (3,) ndarray. Its
                                    units are meters.
    :param back_right_vortex_vertex: This parameter is the ring vortex's back right point. It's a (3,) ndarray. Its
                                     units are meters.
    :param vortex_strength: This is the magnitude of the vorticity. It's sign is given by the using the right hand rule
                            clockwise around the ring vortex when viewed from above. Its units are meters squared per
                            second.
    :return velocity_induced_by_vortex: This is the velocity induced at the point by the line vortex. It's a (3,) numpy
                                        array. Its units are meters per second.
    """

    # Calculate the velocity induced by each leg of the ring vortex. These function calls will return zero if there is a
    # singularity detected on that leg.
    velocity_induced_by_front_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=front_right_vortex_vertex,
        vortex_termination=front_left_vortex_vertex,
        vortex_strength=vortex_strength
    )
    velocity_induced_by_left_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=front_left_vortex_vertex,
        vortex_termination=back_left_vortex_vertex,
        vortex_strength=vortex_strength
    )
    velocity_induced_by_back_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=back_left_vortex_vertex,
        vortex_termination=back_right_vortex_vertex,
        vortex_strength=vortex_strength
    )
    velocity_induced_by_right_leg = calculate_velocity_induced_from_line_vortex(
        point=point,
        vortex_origin=back_right_vortex_vertex,
        vortex_termination=front_right_vortex_vertex,
        vortex_strength=vortex_strength
    )

    # Sum the velocities induced by each leg to get the velocity induced by the entire ring.
    velocity_induced_by_vortex = (velocity_induced_by_front_leg
                                  + velocity_induced_by_right_leg
                                  + velocity_induced_by_back_leg
                                  + velocity_induced_by_left_leg)

    # Return the velocity induced by the vortex.
    return velocity_induced_by_vortex


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
            This is a vector containing the x, y, and z coordinates of the origin of the line vortex.
        :param termination: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the termination of the line vortex.
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

        # Initialize variables to hold the near field force and moment on the vortex. These will be set to None until a
        # solution is found.
        self.near_field_force = None
        self.near_field_moment = None

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with a unit vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return normalized_induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        normalized_induced_velocity = calculate_velocity_induced_from_line_vortex(point=point,
                                                                                  vortex_origin=self.origin,
                                                                                  vortex_termination=self.termination,
                                                                                  vortex_strength=1)
        return normalized_induced_velocity

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with its given vortex strength.

        :param point:
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        induced_velocity = calculate_velocity_induced_from_line_vortex(point=point, vortex_origin=self.origin,
                                                                       vortex_termination=self.termination,
                                                                       vortex_strength=self.strength)
        return induced_velocity


class HorseshoeVortex:
    """This class is used to contain horseshoe vortices.

    This class contains the following public methods:
        calculate_normalized_induced_velocity: This method calculates the velocity induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_induced_velocity: This method calculates the velocity induced at a point by this vortex with its given
                                    vortex strength.
        calculate_normalized_induced_downwash: This method calculates the downwash induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_induced_downwash: This method calculates the downwash induced at a point by this vortex with its given
                                    vortex strength.
        update_strength: This method updates the strength of this horseshoe vortex object, and the strength of its
                         finite leg's line vortex object.
        update_force_and_moment: This method updates the force and moment on this horseshoe vortex based on the force
                                 on its finite leg's line vortex object.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, finite_leg_origin, finite_leg_termination, strength):
        """This is the initialization method.

        :param finite_leg_origin: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the origin of the vortex's finite leg.
        :param finite_leg_termination: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the termination of the vortex's finite leg.
        :param strength: float
            This is the strength of the vortex in meters squared per second.
        """

        self.finite_leg_origin = finite_leg_origin
        self.finite_leg_termination = finite_leg_termination
        self.strength = strength

        # Initialize a line vortex to represent the horseshoe's finite leg.
        self.finite_leg = LineVortex(origin=self.finite_leg_origin, termination=self.finite_leg_termination,
                                     strength=self.strength)

        # Initialize variables to hold the near field force and moment on the vortex. These will be set to None until a
        # solution is found.
        self.near_field_force = None
        self.near_field_moment = None

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with a unit vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return normalized_induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        normalized_induced_velocity = calculate_velocity_induced_from_horseshoe_vortex(
            point=point, finite_leg_origin=self.finite_leg_origin, finite_leg_termination=self.finite_leg_termination,
            vortex_strength=1)
        return normalized_induced_velocity

    def calculate_normalized_induced_downwash(self, point):
        """This method calculates the downwash induced at a point by this vortex with a unit vortex strength.

        The induced downwash of a horseshoe vortex is the velocity induced by the infinite legs. This is found by
        subtracting the velocity induced by just the finite leg from the total horseshoe vortex's induced velocity.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced downwash at.
        :return normalized_induced_downwash: 1D ndarray
            This is a vector containing the x, y, and z components of the induced downwash.
        """

        # Find this vortex's total normalized induced velocity, and the normalized induced velocity from its finite leg.
        normalized_induced_velocity = self.calculate_normalized_induced_velocity(point=point)
        normalized_induced_velocity_from_finite_leg = self.finite_leg.calculate_normalized_induced_velocity(point=point)

        # Find the normalized downwash.
        normalized_induced_downwash = normalized_induced_velocity - normalized_induced_velocity_from_finite_leg

        return normalized_induced_downwash

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with its given vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        induced_velocity = calculate_velocity_induced_from_horseshoe_vortex(
            point=point, finite_leg_origin=self.finite_leg_origin, finite_leg_termination=self.finite_leg_termination,
            vortex_strength=self.strength)
        return induced_velocity

    def calculate_induced_downwash(self, point):
        """This method calculates the downwash induced at a point by this vortex with its given vortex strength.

        The induced downwash of a horseshoe vortex is the velocity induced by the infinite legs. This is found by
        subtracting the velocity induced by just the finite leg from the total horseshoe vortex's induced velocity.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced downwash at.
        :return induced_downwash: 1D ndarray
            This is a vector containing the x, y, and z components of the induced downwash.
        """

        # Find this vortex's total induced velocity, and the induced velocity from its finite leg.
        induced_velocity = self.calculate_induced_velocity(point=point)
        induced_velocity_from_finite_leg = self.finite_leg.calculate_induced_velocity(point=point)

        # Find the downwash.
        induced_downwash = induced_velocity - induced_velocity_from_finite_leg

        return induced_downwash

    def update_strength(self, strength):
        """This method updates the strength of this horseshoe vortex object, and the strength of its finite leg's line
        vortex object.

        :param strength: float
            This is the strength of this vortex, and of its finite leg's line vortex. Its units are meters squared per
            second.
        :return: None
        """

        self.strength = strength
        self.finite_leg.strength = strength

    def update_force_and_moment(self):
        """This method updates the force and moment on this horseshoe vortex based on the force on its finite leg's line
        vortex object.

        :return: None
        """

        self.near_field_force = self.finite_leg.near_field_force
        self.near_field_moment = self.finite_leg.near_field_moment


class RingVortex:
    """This class is used to contain ring vortices.

    This class contains the following public methods:
        calculate_normalized_induced_velocity: This method calculates the velocity induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_normalized_induced_downwash: This method calculates the downwash induced at a point by this vortex
                                               with a unit vortex strength.
        calculate_induced_velocity: This method calculates the velocity induced at a point by this vortex with its given
                                    vortex strength.
        calculate_induced_downwash: This method calculates the downwash induced at a point by this vortex with its given
                                    vortex strength.
        update_strength: This method updates the strength of this ring vortex object, and the strength of its
                         four legs' line vortex objects.
        update_force_and_moment: This method updates the force and moment on this ring vortex based on the force
                                 on its four legs' line vortex objects.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, front_left_vertex, front_right_vertex, back_left_vertex, back_right_vertex, strength):
        """This is the initialization method.

        :param front_left_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's front left point.
        :param front_right_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's front right point.
        :param back_left_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's back left point.
        :param back_right_vertex: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the vortex's back right point.
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
            strength=self.strength
        )
        self.left_leg = LineVortex(
            origin=self.front_left_vertex,
            termination=self.back_left_vertex,
            strength=self.strength
        )
        self.back_leg = LineVortex(
            origin=self.back_left_vertex,
            termination=self.back_right_vertex,
            strength=self.strength
        )
        self.right_leg = LineVortex(
            origin=self.back_right_vertex,
            termination=self.front_right_vertex,
            strength=self.strength
        )

        # Initialize a variable to hold the centroid of the ring vortex.
        self.center = asmvp.geometry.centroid_of_quadrilateral(front_left_vertex, front_right_vertex, back_left_vertex,
                                                               back_right_vertex)

        # Initialize variables to hold the near field force and moment on the vortex. These will be set to None until a
        # solution is found.
        self.near_field_force = None
        self.near_field_moment = None

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with a unit vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return normalized_induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        normalized_induced_velocity = calculate_velocity_induced_from_ring_vortex(
            point=point,
            front_left_vortex_vertex=self.front_left_vertex,
            front_right_vortex_vertex=self.front_right_vertex,
            back_left_vortex_vertex=self.back_left_vertex,
            back_right_vortex_vertex=self.back_right_vertex,
            vortex_strength=1
        )
        return normalized_induced_velocity

    def calculate_normalized_induced_downwash(self, point):
        """This method calculates the downwash induced at a point by this vortex with a unit vortex strength.

        The induced downwash of a ring vortex is the velocity induced by it's left and right legs.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced downwash at.
        :return normalized_induced_downwash: 1D ndarray
            This is a vector containing the x, y, and z components of the induced downwash.
        """

        # Find the normalized velocity induced by the ring's left and right legs.
        normalized_velocity_induced_by_left_leg = self.left_leg.calculate_normalized_induced_velocity(point)
        normalized_velocity_induced_by_right_leg = self.right_leg.calculate_normalized_induced_velocity(point)

        # Find the normalized induced downwash.
        normalized_induced_downwash = normalized_velocity_induced_by_left_leg + normalized_velocity_induced_by_right_leg

        return normalized_induced_downwash

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this vortex with its given vortex strength.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return induced_velocity: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        induced_velocity = calculate_velocity_induced_from_ring_vortex(
            point=point,
            front_left_vortex_vertex=self.front_left_vertex,
            front_right_vortex_vertex=self.front_right_vertex,
            back_left_vortex_vertex=self.back_left_vertex,
            back_right_vortex_vertex=self.back_right_vertex,
            vortex_strength=self.strength
        )
        return induced_velocity

    def calculate_induced_downwash(self, point):
        """This method calculates the downwash induced at a point by this vortex with its given vortex strength.

        The induced downwash of a ring vortex is the velocity induced by it's left and right legs.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced downwash at.
        :return induced_downwash: 1D ndarray
            This is a vector containing the x, y, and z components of the induced downwash.
        """

        # Find the normalized velocity induced by the ring's left and right legs.
        velocity_induced_by_left_leg = self.left_leg.calculate_induced_velocity(point)
        velocity_induced_by_right_leg = self.right_leg.calculate_induced_velocity(point)

        # Find the normalized induced downwash.
        induced_downwash = velocity_induced_by_left_leg + velocity_induced_by_right_leg

        return induced_downwash

    def update_strength(self, strength):
        """This method updates the strength of this ring vortex object, and the strength of its four legs' line vortex
        objects.

        :param strength: float
            This is the strength of this vortex, and of its four legs' line vortices. Its units are meters squared per
            second.
        :return: None
        """

        self.strength = strength
        self.front_leg.strength = strength
        self.left_leg.strength = strength
        self.back_leg.strength = strength
        self.right_leg.strength = strength

    def update_force_and_moment(self):
        """This method updates the force and moment on this ring vortex based on the force on its four legs' line
        objects.

        :return: None
        """

        self.near_field_force = (self.front_leg.near_field_force + self.left_leg.near_field_force
                                 + self.back_leg.near_field_force + self.right_leg.near_field_force)
        self.near_field_moment = np.cross(self.near_field_force, self.center)
