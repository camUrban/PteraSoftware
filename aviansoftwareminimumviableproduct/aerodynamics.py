"""This module contains useful aerodynamics functions.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    induced_velocity_from_line_vortex: Calculates the velocity induced at a specified point by a line vortex of a
                                       specified strength, and at a specified location.
    induced_velocity_from_horseshoe_vortex: Calculates the velocity induced at a specified point by a line vortex of a
                                            specified strength, and at a specified location.
    induced_velocity_from_ring_vortex: Calculates the velocity induced at a specified point by a line vortex of a
                                       specified strength, and at a specified location.
"""

import numpy as np


def induced_velocity_from_line_vortex(point, vortex_origin, vortex_termination, vortex_strength):
    """This function calculates the velocity induced at a specified point by a line vortex of a specified strength,
    and at a specified location. It uses methodology described on pp. 251-255 of the second edition of "Low-Speed
    Aerodynamics" by Joseph Katz and Allen Plotkin.

    :param point: This parameter is the point where the induced velocity is to be calculated. It's a (3,) numpy array.
                  It's units are meters.
    :param vortex_origin: This parameter is the point where the line vortex begins. It's a (3,) numpy array. It's units
                          are meters.
    :param vortex_termination: This parameter is the point where the line vortex ends. It's a (3,) numpy array. It's
                               units are meters.
    :param vortex_strength: This is the magnitude of the vorticity. It's sign is given by the using the right hand rule
                            on the vector from the line vortex's origin to termination. It's units are meters squared
                            per second.
    :return velocity_induced_by_vortex: This is the velocity induced at the point by the line vortex. It's a (3,) numpy
                                        array. It's units are meters per second.
    """

    # Define vectors from the vortex to the point.
    r_1 = point - vortex_origin
    r_2 = point - vortex_termination

    # Define vector from vortex origin to vortex termination.
    r_0 = r_1 - r_2

    # Calculate vector cross product.
    r_1_cross_r_2 = np.cross(r_1, r_2)

    # Calculate cross product absolute value.
    r_1_cross_r_2_absolute_value = r_1_cross_r_2[0] ** 2 + r_1_cross_r_2[1] ** 2 + r_1_cross_r_2[2] ** 2

    # Calculate vector lengths.
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

        # Calculate the induced velocity components, and combine them into the induced velocity numpy array.
        u = k * r_1_cross_r_2[0]
        v = k * r_1_cross_r_2[1]
        w = k * r_1_cross_r_2[2]
        velocity_induced_by_vortex = np.array([u, v, w])

    # Return the induced velocity.
    return velocity_induced_by_vortex


def induced_velocity_from_horseshoe_vortex(point, finite_leg_origin, finite_leg_termination, vortex_strength):
    """This function calculates the velocity induced at a specified point by a horseshoe vortex of a specified strength,
    and at a specified location. It uses methodology described on pp. 131-132 of "Flight Vehicle Aerodynamics" by Mark
    Drela.

    Notes:
        1. This function assumes that the horseshoe vortex's infinite legs, which are rays originating from the finite
        leg's origin and termination points, go straight back in the positive x direction.

    :param point: This parameter is the point where the induced velocity is to be calculated. It's a (3,) numpy array.
                  It's units are meters.
    :param finite_leg_origin: This parameter is the point where the horseshoe vortex's finite leg begins. It's a (3,)
                              numpy array. It's units are meters.
    :param finite_leg_termination: This parameter is the point where the horseshoe vortex's finite leg ends. It's a (3,)
                                   numpy array. It's units are meters.
    :param vortex_strength: This is the magnitude of the vorticity. It's sign is given by the using the right hand rule
                            on the vector from the horseshoe vortex's finite leg's origin to termination. It's units are
                            meters squared per second.
    :return velocity_induced_by_vortex: This is the velocity induced at the point by the horseshoe vortex. It's a (3,)
                                        numpy array. It's units are meters per second.
    """

    # Define vectors from the vortex's finite leg to the point.
    a = point - finite_leg_origin
    b = point - finite_leg_termination

    # Calculate vector lengths.
    a_length = np.linalg.norm(a)
    b_length = np.linalg.norm(b)

    # Define the x direction unit vector.
    x_unit = np.array([1, 0, 0])

    # Find a point on an infinite line that encompasses the first infinite leg, which is a ray. This assumes the leg
    # goes straight back in the positive x direction.
    point_on_first_infinite_leg = finite_leg_origin + x_unit
    # Find the distance between the point and this line. Uses: Weisstein, Eric W. "Point-Line
    # Distance--3-Dimensional." From MathWorld--A Wolfram Web Resource.
    # https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    distance_along_first_infinite_leg = (np.linalg.norm(np.cross((point - finite_leg_origin),
                                                                 (point - point_on_first_infinite_leg)))
                                         / np.linalg.norm(point_on_first_infinite_leg - finite_leg_origin))

    # Find a point on an infinite line that encompasses the second infinite leg, which is a ray. This assumes the leg
    # goes straight back in the positive x direction.
    point_on_second_infinite_leg = finite_leg_termination + x_unit
    # Find the distance between the point and this line. Uses: Weisstein, Eric W. "Point-Line
    # Distance--3-Dimensional." From MathWorld--A Wolfram Web Resource.
    # https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    distance_along_second_infinite_leg = (np.linalg.norm(np.cross((point - finite_leg_termination),
                                                                  (point - point_on_second_infinite_leg)))
                                          / np.linalg.norm(point_on_second_infinite_leg
                                                           - finite_leg_termination))

    # Check for singularities.
    horseshoe_vortex_radius = 3.0e-16
    if (a_length < horseshoe_vortex_radius or b_length < horseshoe_vortex_radius
            or distance_along_first_infinite_leg < horseshoe_vortex_radius
            or distance_along_second_infinite_leg < horseshoe_vortex_radius):
        # If there is a singularity, the induced velocity is zero.
        velocity_induced_by_vortex = np.array([0, 0, 0])
    else:
        # Calculate the vector dot and cross products.
        a_cross_b = np.cross(a, b)
        a_dot_b = np.dot(a, b)
        a_cross_x_unit = np.cross(a, x_unit)
        b_cross_x_unit = np.cross(b, x_unit)
        a_dot_x_unit = np.dot(a, x_unit)
        b_dot_x_unit = np.dot(b, x_unit)

        # Calculate the velocity induced by the vortex.
        velocity_induced_by_vortex = (vortex_strength / (4 * np.pi)
                                      * ((a_cross_b / (a_length * b_length + a_dot_b))
                                         * (1 / a_length + 1 / b_length)
                                         + (a_cross_x_unit / (a_length - a_dot_x_unit))
                                         * (1 / a_length)
                                         - (b_cross_x_unit / (b_length - b_dot_x_unit))
                                         * (1 / b_length)))

    # Return the velocity induced by the vortex.
    return velocity_induced_by_vortex


def induced_velocity_from_ring_vortex(point, front_left_vortex_point, front_right_vortex_point, back_left_vortex_point,
                                      back_right_vortex_point, vortex_strength):
    """This function calculates the velocity induced at a specified point by a ring vortex of a specified strength,
    and at a specified location. It uses methodology described on pp. 255-256 of the second edition of "Low-Speed
    Aerodynamics" by Joseph Katz and Allen Plotkin.

    :param point: This parameter is the point where the induced velocity is to be calculated. It's a (3,) numpy array.
                  It's units are meters.
    :param front_left_vortex_point: This parameter is the ring vortex's front left point. It's a (3,) numpy array. It's
                                    units are meters.
    :param front_right_vortex_point: This parameter is the ring vortex's front right point. It's a (3,) numpy array.
                                     It's units are meters.
    :param back_left_vortex_point: This parameter is the ring vortex's back left point. It's a (3,) numpy array. It's
                                   units are meters.
    :param back_right_vortex_point: This parameter is the ring vortex's back right point. It's a (3,) numpy array. It's
                                    units are meters.
    :param vortex_strength: This is the magnitude of the vorticity. It's sign is given by the using the right hand rule
                            clockwise around the ring vortex when viewed from above. It's units are meters squared per
                            second.
    :return velocity_induced_by_vortex: This is the velocity induced at the point by the line vortex. It's a (3,) numpy
                                        array. It's units are meters per second.
    """

    # Calculate the velocity induced by each leg of the ring vortex. These function calls will return zero if there is a
    # singularity detected on that leg.
    velocity_induced_by_upper_leg = induced_velocity_from_line_vortex(point=point,
                                                                      vortex_origin=front_right_vortex_point,
                                                                      vortex_termination=front_left_vortex_point,
                                                                      vortex_strength=vortex_strength)
    velocity_induced_by_left_leg = induced_velocity_from_line_vortex(point=point,
                                                                     vortex_origin=front_left_vortex_point,
                                                                     vortex_termination=back_left_vortex_point,
                                                                     vortex_strength=vortex_strength)
    velocity_induced_by_lower_leg = induced_velocity_from_line_vortex(point=point,
                                                                      vortex_origin=back_left_vortex_point,
                                                                      vortex_termination=back_right_vortex_point,
                                                                      vortex_strength=vortex_strength)
    velocity_induced_by_right_leg = induced_velocity_from_line_vortex(point=point,
                                                                      vortex_origin=back_right_vortex_point,
                                                                      vortex_termination=front_right_vortex_point,
                                                                      vortex_strength=vortex_strength)

    # Sum the velocities induced by each leg to get the velocity induced by the entire ring.
    velocity_induced_by_vortex = (velocity_induced_by_upper_leg
                                  + velocity_induced_by_left_leg
                                  + velocity_induced_by_lower_leg
                                  + velocity_induced_by_right_leg)

    # Return the velocity induced by the vortex.
    return velocity_induced_by_vortex
