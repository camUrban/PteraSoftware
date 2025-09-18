"""This module contains functions shared by other modules in the pterasoftware package.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    cosspace: This function is used to create an array containing a specified number
    of values between a specified minimum and maximum value that are spaced via a
    cosine function.

    angle_axis_rotation_matrix: This function is used to find the rotation matrix for
    a given axis and angle.

    numba_centroid_of_quadrilateral: This function is used to find the centroid of a
    quadrilateral. It has been optimized for JIT compilation using Numba.

    calculate_streamlines: This function calculates the location of the streamlines
    coming off the back of the wings.

    convert_logging_level_name_to_value: This function takes in a string that
    represents the logging level and returns the integer that can be used to set the
    logger to this level.

    process_steady_solver_forces: This function uses the forces and moments a solver
    has found on its panels to find the forces, moments, and associated coefficients
    on each airplane in the solver.

    process_unsteady_solver_forces: This function uses the forces and moments a
    solver has found on its panels to find the forces, moments, and associated
    coefficients on each airplane in the solver.

    update_ring_vortex_solvers_panel_attributes: This function populates a ring
    vortex solver's attributes with the attributes of a given panel.

    calculate_steady_freestream_wing_influences: This function finds the vector of
    freestream-wing influence coefficients associated with this problem.

    numba_1d_explicit_cross: This function takes in two arrays, each of which contain
    N vectors of 3 components. The function then calculates and returns the cross
    product of the two vectors at each position.

    reflect_point_across_plane: This function finds the coordinates of the reflection
    of a point across a plane in 3D space.

    interp_between_points: This function finds the MxN points between M pairs of
    points in 3D space given an array of N normalized spacings.
"""

import logging

import numpy as np
from numba import njit


def cosspace(minimum, maximum, n_points=50, endpoint=True):
    """This function is used to create an array containing a specified number of
    values between a specified minimum and maximum value that are spaced via a cosine
    function.

    Citation:
        Adapted from:         geometry.cosspace in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    :param minimum: float
        This is the minimum value of the range of numbers you would like spaced.
    :param maximum: float
        This is the maximum value of the range of numbers you would like spaced.
    :param n_points: int, optional
        This is the number of points to space. The default is 50.
    :param endpoint: bool, optional
        This determines if the maximum value will be included in the output. The
        default is True.
    :return cosine_spaced_points: 1D array
        This is a 1D array of the points, ranging from the minimum to the maximum
        value (inclusive), spaced via a cosine function.
    """

    # Find the mean and the amplitude of the cosine function.
    mean = (maximum + minimum) / 2
    amp = (maximum - minimum) / 2

    # Space the points by applying cosine to the linspace function. Then return the
    # points.
    cosine_spaced_points = mean + amp * np.cos(
        np.linspace(np.pi, 0, n_points, endpoint=endpoint)
    )
    return cosine_spaced_points


def angle_axis_rotation_matrix(angle, axis, axis_already_normalized=False):
    """This function is used to find the rotation matrix for a given axis and angle.

    Citation:
        Adapted from:         geometry.angle_axis_rotation_matrix in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    For more information on the math behind this method, see:
    https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    :param angle: float or 1D array of 3 angles
        This is the angle. Provide it in radians.
    :param axis: 1D array of 3 elements.
        This is the axis.
    :param axis_already_normalized: bool, optional
        Set this as True if the axis given is already normalized, which will skip
        internal normalization and speed up the method. If not, set as False. The
        default is False.
    :return rotation matrix: array
        This is the rotation matrix. If the given angle is a scalar, this will be a 3
        x 3 array. If the given angle is a vector, this will be a 3 x 3 x N array.
    """

    # Normalize the axis is it is not already normalized.
    if not axis_already_normalized:
        axis = axis / np.linalg.norm(axis)

    # Define constants to use when calculating the rotation matrix.
    sin_theta = np.sin(angle)
    cos_theta = np.cos(angle)
    u_x = axis[0]
    u_y = axis[1]
    u_z = axis[2]

    # Calculate the rotation matrix.
    rotation_matrix = np.array(
        [
            [
                cos_theta + u_x**2 * (1 - cos_theta),
                u_x * u_y * (1 - cos_theta) - u_z * sin_theta,
                u_x * u_z * (1 - cos_theta) + u_y * sin_theta,
            ],
            [
                u_y * u_x * (1 - cos_theta) + u_z * sin_theta,
                cos_theta + u_y**2 * (1 - cos_theta),
                u_y * u_z * (1 - cos_theta) - u_x * sin_theta,
            ],
            [
                u_z * u_x * (1 - cos_theta) - u_y * sin_theta,
                u_z * u_y * (1 - cos_theta) + u_x * sin_theta,
                cos_theta + u_z**2 * (1 - cos_theta),
            ],
        ]
    )

    # Return the rotation matrix.
    return rotation_matrix


@njit(cache=True, fastmath=False)
def numba_centroid_of_quadrilateral(
    front_left_vertex, front_right_vertex, back_left_vertex, back_right_vertex
):
    """This function is used to find the centroid of a quadrilateral. It has been
    optimized for JIT compilation using Numba.

    :param front_left_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the front left
        vertex of the quadrilateral.
    :param front_right_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the front right
        vertex of the quadrilateral.
    :param back_left_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the back left
        vertex of the quadrilateral.
    :param back_right_vertex: 1D array of floats
        This is an array containing the x, y, and z components of the back right
        vertex of the quadrilateral.
    :return: 1D array of floats
        This is an array containing the x, y, and z components of the centroid of the
        quadrilateral.
    """

    x_average = (
        front_left_vertex[0]
        + front_right_vertex[0]
        + back_left_vertex[0]
        + back_right_vertex[0]
    ) / 4
    y_average = (
        front_left_vertex[1]
        + front_right_vertex[1]
        + back_left_vertex[1]
        + back_right_vertex[1]
    ) / 4
    z_average = (
        front_left_vertex[2]
        + front_right_vertex[2]
        + back_left_vertex[2]
        + back_right_vertex[2]
    ) / 4

    return np.array([x_average, y_average, z_average])


def calculate_streamlines(solver, num_steps=25, delta_time=0.02):
    """This function calculates the location of the streamlines coming off the back
    of the wings.

    This method is vectorized to increase performance.

    :param solver: Solver
        This is the solver object for which to calculate the streamlines.
    :param num_steps: int, optional
        This is the integer number of points along each streamline (not including the
        initial points). It can be increased for higher fidelity visuals. The default
        value is 25.
    :param delta_time: float, optional
        This is the time in seconds between each time step It can be
        decreased for higher fidelity visuals or to make the streamlines shorter.
        Its default value is 0.02 seconds.
    :return: None
    """
    # Initialize an array to hold this solver's matrix of streamline points.
    solver.streamline_points = np.expand_dims(solver.seed_points, axis=0)

    # Iterate through the streamline steps.
    for step in range(num_steps):
        # Get the last row of streamline points.
        last_row_streamline_points = solver.streamline_points[-1, :, :]

        # Add the freestream velocity to the induced velocity to get the total
        # velocity at each of the last row of streamline points.
        total_velocities = solver.calculate_solution_velocity(
            points=last_row_streamline_points
        )

        # Interpolate the positions on a new row of streamline points.
        new_row_streamline_points = (
            last_row_streamline_points + total_velocities * delta_time
        )

        # Stack the new row of streamline points to the bottom of the matrix of
        # streamline points.
        solver.streamline_points = np.vstack(
            (
                solver.streamline_points,
                np.expand_dims(new_row_streamline_points, axis=0),
            )
        )


def convert_logging_level_name_to_value(name):
    """This function takes in a string that represents the logging level and returns
    the integer that can be used to set the logger to this level.

    :param name: str
        This is the string representation of the logging level. The options are
        "Debug", "Info", "Warning", "Error", and "Critical".
    :return: int
        This is the integer value that can used to set the appropriate logging level.
    """
    logging_levels = {
        "Debug": logging.DEBUG,
        "Info": logging.INFO,
        "Warning": logging.WARNING,
        "Error": logging.ERROR,
        "Critical": logging.CRITICAL,
    }
    try:
        return logging_levels[name]
    except KeyError:
        raise Exception("The name of the logging level provided is not a valid option.")


def process_steady_solver_forces(
    steady_solver, near_field_forces_geometry_axes, near_field_moments_geometry_axes
):
    """This function uses the forces and moments a solver has found on its panels to
    find the forces, moments, and associated coefficients on each airplane in the
    solver.

    :param steady_solver: SteadySolver
        This is the solver whose forces will be processed.
    :param near_field_forces_geometry_axes: Nx3 array of floats
        This is an array of the near field forces in geometry axes on each of the
        solver's panels. The array is size Nx3, where N is the number of panels. The
        units are Newtons.
    :param near_field_moments_geometry_axes: Nx3 array of floats
        This is an array of the near field moments in geometry axes on each of the
        solver's panels. The array is size Nx3, where N is the number of panels. The
        units are Newton-meters.
    :return:
    """
    # Find this operating point's dynamic pressure. The units are Pascals.
    dynamic_pressure = steady_solver.operating_point.calculate_dynamic_pressure()

    # Find the rotation matrix that will be used to convert the geometry frame values
    # to wind frame values.
    rotation_matrix = np.transpose(
        steady_solver.operating_point.calculate_rotation_matrix_wind_to_geometry()
    )

    # Iterate through this solver's panels.
    for panel_num, panel in enumerate(steady_solver.panels):
        # Get this panel's near field forces and moments in geometry axes and wind axes.
        this_force_geometry_axes = near_field_forces_geometry_axes[panel_num, :]
        this_moment_geometry_axes = near_field_moments_geometry_axes[panel_num, :]
        this_force_wind_axes = rotation_matrix @ this_force_geometry_axes
        this_moment_wind_axes = rotation_matrix @ this_moment_geometry_axes

        # Update the force and moment on this panel.
        panel.forces_G = this_force_geometry_axes
        panel.moments_G_Cg = this_moment_geometry_axes
        panel.forces_W = this_force_wind_axes
        panel.moments_W_Cg = this_moment_wind_axes

        # Update the force coefficients this panel.
        panel.update_coefficients(dynamic_pressure)

    # Initialize arrays to hold each airplane's total force and moment in geometry
    # axes.
    total_near_field_forces_geometry_axes = np.zeros((steady_solver.num_airplanes, 3))
    total_near_field_moments_geometry_axes = np.zeros((steady_solver.num_airplanes, 3))

    # Iterate through each airplane and find the total force and moment experienced
    # by each by summing up the contribution's from its panels.
    for airplane_num, airplane in enumerate(steady_solver.airplanes):
        for wing in airplane.wings:
            for panel in np.ravel(wing.panels):
                total_near_field_forces_geometry_axes[
                    airplane_num, :
                ] += panel.near_field_force_geometry_axes
                total_near_field_moments_geometry_axes[
                    airplane_num, :
                ] += panel.near_field_moment_geometry_axes

    # For each airplane, find the total force and moment it experiences in wind axes
    # from the rotation matrix and the total force and moment it experiences in
    # geometry axes.
    for airplane_num, airplane in enumerate(steady_solver.airplanes):
        airplane.total_near_field_force_wind_axes = (
            rotation_matrix @ total_near_field_forces_geometry_axes[airplane_num]
        )
        airplane.total_near_field_moment_wind_axes = (
            rotation_matrix @ total_near_field_moments_geometry_axes[airplane_num]
        )

    # Iterate through the airplanes and calculate each one's coefficients.
    for airplane in steady_solver.airplanes:
        # Calculate this airplane's force coefficients.
        induced_drag_coefficient = (
            -airplane.total_near_field_force_wind_axes[0]
            / dynamic_pressure
            / airplane.s_ref
        )
        side_force_coefficient = (
            airplane.total_near_field_force_wind_axes[1]
            / dynamic_pressure
            / airplane.s_ref
        )
        lift_coefficient = (
            -airplane.total_near_field_force_wind_axes[2]
            / dynamic_pressure
            / airplane.s_ref
        )

        # Calculate this airplane's moment coefficients.
        rolling_moment_coefficient = (
            airplane.total_near_field_moment_wind_axes[0]
            / dynamic_pressure
            / airplane.s_ref
            / airplane.b_ref
        )
        pitching_moment_coefficient = (
            airplane.total_near_field_moment_wind_axes[1]
            / dynamic_pressure
            / airplane.s_ref
            / airplane.c_ref
        )
        yawing_moment_coefficient = (
            airplane.total_near_field_moment_wind_axes[2]
            / dynamic_pressure
            / airplane.s_ref
            / airplane.b_ref
        )

        # Populate this airplane's force and moment coefficient attributes.
        airplane.total_near_field_force_coefficients_wind_axes = np.array(
            [
                induced_drag_coefficient,
                side_force_coefficient,
                lift_coefficient,
            ]
        )
        airplane.total_near_field_moment_coefficients_wind_axes = np.array(
            [
                rolling_moment_coefficient,
                pitching_moment_coefficient,
                yawing_moment_coefficient,
            ]
        )


def process_unsteady_solver_forces(
    unsteady_solver, near_field_forces_geometry_axes, near_field_moments_geometry_axes
):
    """This function uses the forces and moments a solver has found on its panels to
    find the forces, moments, and associated coefficients on each airplane in the
    solver.

    :param unsteady_solver: UnsteadySolver
        This is the solver whose forces will be processed.
    :param near_field_forces_geometry_axes: Nx3 array of floats
        This is an array of the near field forces in geometry axes on each of the
        solver's panels. The array is size Nx3, where N is the number of panels. The
        units are Newtons.
    :param near_field_moments_geometry_axes: Nx3 array of floats
        This is an array of the near field moments in geometry axes on each of the
        solver's panels. The array is size Nx3, where N is the number of panels. The
        units are Newton-meters.
    :return:
    """
    operating_point = unsteady_solver.current_operating_point

    # Find this operating point's dynamic pressure. The units are Pascals.
    dynamic_pressure = operating_point.calculate_dynamic_pressure()

    # Find the rotation matrix that will be used to convert the geometry frame values
    # to wind frame values.
    rotation_matrix = np.transpose(
        operating_point.calculate_rotation_matrix_wind_to_geometry()
    )

    # Iterate through this solver's panels.
    for panel_num, panel in enumerate(unsteady_solver.panels):
        # Get this panel's near field forces and moments in geometry axes and wind
        # axes.
        this_force_geometry_axes = near_field_forces_geometry_axes[panel_num, :]
        this_moment_geometry_axes = near_field_moments_geometry_axes[panel_num, :]
        this_force_wind_axes = rotation_matrix @ this_force_geometry_axes
        this_moment_wind_axes = rotation_matrix @ this_moment_geometry_axes

        # Update the force and moment on this panel.
        panel.forces_G = this_force_geometry_axes
        panel.moments_G_Cg = this_moment_geometry_axes
        panel.forces_W = this_force_wind_axes
        panel.moments_W_Cg = this_moment_wind_axes

        # Update the force coefficients on this panel.
        panel.update_coefficients(dynamic_pressure=dynamic_pressure)

    # Initialize arrays for each airplane's total force and moment in geometry axes.
    num_airplanes = unsteady_solver.num_airplanes
    total_near_field_forces_geometry_axes = np.zeros((num_airplanes, 3))
    total_near_field_moments_geometry_axes = np.zeros((num_airplanes, 3))

    # Iterate through each airplane and find the total force and moment experienced
    # by each by summing up the contribution's from its panels.
    for airplane_num, airplane in enumerate(unsteady_solver.current_airplanes):
        for wing in airplane.wings:
            for panel in np.ravel(wing.panels):
                total_near_field_forces_geometry_axes[
                    airplane_num, :
                ] += panel.near_field_force_geometry_axes
                total_near_field_moments_geometry_axes[
                    airplane_num, :
                ] += panel.near_field_moment_geometry_axes

    # For each airplane, find the total force and moment it experiences in wind axes
    # from the rotation matrix and the total force and moment it experiences in
    # geometry axes.
    for airplane_num, airplane in enumerate(unsteady_solver.current_airplanes):
        airplane.total_near_field_force_wind_axes = (
            rotation_matrix @ total_near_field_forces_geometry_axes[airplane_num]
        )
        airplane.total_near_field_moment_wind_axes = (
            rotation_matrix @ total_near_field_moments_geometry_axes[airplane_num]
        )

    # Iterate through the airplanes and calculate each one's coefficients.
    for airplane in unsteady_solver.current_airplanes:
        # Calculate this airplane's force coefficients.
        induced_drag_coefficient = (
            -airplane.total_near_field_force_wind_axes[0]
            / dynamic_pressure
            / airplane.s_ref
        )
        side_force_coefficient = (
            airplane.total_near_field_force_wind_axes[1]
            / dynamic_pressure
            / airplane.s_ref
        )
        lift_coefficient = (
            -airplane.total_near_field_force_wind_axes[2]
            / dynamic_pressure
            / airplane.s_ref
        )

        # Calculate this airplane's moment coefficients.
        rolling_moment_coefficient = (
            airplane.total_near_field_moment_wind_axes[0]
            / dynamic_pressure
            / airplane.s_ref
            / airplane.b_ref
        )
        pitching_moment_coefficient = (
            airplane.total_near_field_moment_wind_axes[1]
            / dynamic_pressure
            / airplane.s_ref
            / airplane.c_ref
        )
        yawing_moment_coefficient = (
            airplane.total_near_field_moment_wind_axes[2]
            / dynamic_pressure
            / airplane.s_ref
            / airplane.b_ref
        )

        # Populate this airplane's force and moment coefficient attributes.
        airplane.total_near_field_force_coefficients_wind_axes = np.array(
            [
                induced_drag_coefficient,
                side_force_coefficient,
                lift_coefficient,
            ]
        )
        airplane.total_near_field_moment_coefficients_wind_axes = np.array(
            [
                rolling_moment_coefficient,
                pitching_moment_coefficient,
                yawing_moment_coefficient,
            ]
        )


def update_ring_vortex_solvers_panel_attributes(
    solver, global_panel_position, panel, airplane
):
    """This function populates a ring vortex solver's attributes with the attributes
    of a given panel.

    :param solver: SteadySolver or UnsteadySolver
        This is the solver object whose attributes are to be updated. It should be a
        solver that uses ring vortices.
    :param global_panel_position: int
        This is the position of the panel with respect to the global array of all
        panels.
    :param panel: Panel
        This is the panel object whose attributes will be used to update the solver's
        attributes.
    :param airplane: Airplane
        This is the Airplane object to which the Panel object belongs.
    :return:
    """

    # Update the solver's list of attributes with this panel's attributes.
    solver.panels[global_panel_position] = panel
    solver.panel_normal_directions[global_panel_position, :] = panel.unitNormal_G
    solver.panel_areas[global_panel_position] = panel.area
    solver.panel_collocation_points[global_panel_position, :] = panel.Cpp_G_Cg
    solver.panel_back_right_vortex_vertices[global_panel_position, :] = (
        panel.ring_vortex.rightLeg_G.origin
    )
    solver.panel_front_right_vortex_vertices[global_panel_position, :] = (
        panel.ring_vortex.rightLeg_G.termination
    )
    solver.panel_front_left_vortex_vertices[global_panel_position, :] = (
        panel.ring_vortex.leftLeg_G.origin
    )
    solver.panel_back_left_vortex_vertices[global_panel_position, :] = (
        panel.ring_vortex.leftLeg_G.termination
    )
    solver.panel_right_vortex_centers[global_panel_position, :] = (
        panel.ring_vortex.rightLeg_G.Crvp_G_Cg
    )
    solver.panel_right_vortex_vectors[global_panel_position, :] = (
        panel.ring_vortex.rightLeg_G.vector
    )
    solver.panel_front_vortex_centers[global_panel_position, :] = (
        panel.ring_vortex.frontLeg_G.Crvp_G_Cg
    )
    solver.panel_front_vortex_vectors[global_panel_position, :] = (
        panel.ring_vortex.frontLeg_G.vector
    )
    solver.panel_left_vortex_centers[global_panel_position, :] = (
        panel.ring_vortex.leftLeg_G.Crvp_G_Cg
    )
    solver.panel_left_vortex_vectors[global_panel_position, :] = (
        panel.ring_vortex.leftLeg_G.vector
    )
    solver.panel_back_vortex_centers[global_panel_position, :] = (
        panel.ring_vortex.backLeg_G.Crvp_G_Cg
    )
    solver.panel_back_vortex_vectors[global_panel_position, :] = (
        panel.ring_vortex.backLeg_G.vector
    )
    solver.panel_is_trailing_edge[global_panel_position] = panel.is_trailing_edge
    solver.panel_is_leading_edge[global_panel_position] = panel.is_leading_edge
    solver.panel_is_right_edge[global_panel_position] = panel.is_right_edge
    solver.panel_is_left_edge[global_panel_position] = panel.is_left_edge
    solver.panel_moment_references[global_panel_position, :] = airplane.xyz_ref

    # Check if this panel is on the trailing edge. If it is, calculate its
    # streamline seed point and add it to the solver's # array of seed points.
    if panel.is_trailing_edge:
        solver.seed_points = np.vstack(
            (
                solver.seed_points,
                panel.Blpp_G_Cg + 0.5 * (panel.Brpp_G_Cg - panel.Blpp_G_Cg),
            )
        )


def calculate_steady_freestream_wing_influences(steady_solver):
    """This function finds the vector of freestream-wing influence coefficients
    associated with this problem.

    :return: None
    """
    # Take the batch dot product of the freestream velocity with each panel's
    # normal direction. This is now the problem's 1D array of freestream-wing
    # influence coefficients.
    steady_solver.freestream_wing_influences = np.einsum(
        "ij,j->i",
        steady_solver.panel_normal_directions,
        steady_solver.freestream_velocity,
    )


@njit(cache=True, fastmath=False)
def numba_1d_explicit_cross(vectors_1, vectors_2):
    """This function takes in two arrays, each of which contain N vectors of 3
    components. The function then calculates and returns the cross product of the two
    vectors at each position.

    Note: This function has been optimized for JIT compilation and parallel
    computation using Numba.

    Citation: Some or all of the following code was written by Jérôme Richard as a
    response to a question on Stack Overflow. The original response is here:
    https://stackoverflow.com/a/66757029/13240504.

    :param vectors_1: array of floats of size (N x 3)
        This is the first array of N vectors.
    :param vectors_2: array of floats of size (N x 3)
        This is the second array of N vectors.
    :return crosses: array of floats of size (N x 3)
        This is the cross product of the two inputted vectors at each of the N
        positions.
    """
    crosses = np.zeros(vectors_1.shape)
    for i in range(crosses.shape[0]):
        crosses[i, 0] = (
            vectors_1[i, 1] * vectors_2[i, 2] - vectors_1[i, 2] * vectors_2[i, 1]
        )
        crosses[i, 1] = (
            vectors_1[i, 2] * vectors_2[i, 0] - vectors_1[i, 0] * vectors_2[i, 2]
        )
        crosses[i, 2] = (
            vectors_1[i, 0] * vectors_2[i, 1] - vectors_1[i, 1] * vectors_2[i, 0]
        )
    return crosses


def reflect_point_across_plane(point, plane_unit_normal, plane_point):
    """This function finds the coordinates of the reflection of a point across a
    plane in 3D space.

    As the plane doesn't necessarily include the origin, this is an affine
    transformation. The transformation matrix and details are discussed in the
    following source: https://en.wikipedia.org/wiki/Transformation_matrix#Reflection_2

    :param point: (3,) array of floats
        This is a (3,) array of the coordinates of the point to reflect. The units
        are meters.
    :param plane_unit_normal: (3,) array of floats
        This is a (3,) array of the components of the plane's unit normal vector. It
        must have unit magnitude. The units are meters.
    :param plane_point: (3,) array of floats
        This is a (3,) array of the coordinates of a point on the plane. The units
        are meters.
    :return: (3,) array of floats
        This is a (3,) array of the components of the reflected point. The units are
        meters.
    """
    [x, y, z] = point
    [a, b, c] = plane_unit_normal

    # Find the dot product between the point on the plane and unit normal.
    d = np.dot(-plane_point, plane_unit_normal)

    # To make the transformation matrix readable, define new variables for the
    # products of the vector components and the dot product.
    ab = a * b
    ac = a * c
    ad = a * d
    bc = b * c
    bd = b * d
    cd = c * d
    a2 = a**2
    b2 = b**2
    c2 = c**2

    # Define the affine transformation matrix. See the citation in the docstring for
    # details.
    transformation_matrix = np.array(
        [
            [1 - 2 * a2, -2 * ab, -2 * ac, -2 * ad],
            [-2 * ab, 1 - 2 * b2, -2 * bc, -2 * bd],
            [-2 * ac, -2 * bc, 1 - 2 * c2, -2 * cd],
            [0, 0, 0, 1],
        ]
    )

    expanded_coordinates = np.array([[x], [y], [z], [1]])

    transformed_expanded_coordinates = transformation_matrix @ expanded_coordinates

    # Extract the reflected coordinates in 3D space.
    [x_prime, y_prime, z_prime] = transformed_expanded_coordinates[:-1, 0]

    return np.array([x_prime, y_prime, z_prime])


@njit(cache=True, fastmath=False)
def interp_between_points(start_points, end_points, norm_spacings):
    """This function finds the MxN points between M pairs of points in 3D space given
    an array of N normalized spacings.

    :param start_points: (M, 3) array of floats
        This is the (M, 3) array containing the coordinates of the M starting points.
        The units are meters.
    :param end_points: (M, 3) array of floats
        This is the (M, 3) array containing the coordinates of the M ending points.
        The units are meters.
    :param norm_spacings: (N,) array of floats
        This is the (N,) array of the N spacings between the starting points and
        ending points. The values are unitless and must be normalized from 0 to 1.
    :return points: (M, N, 3) array of floats
        This is the (M, N, 3) array of the coordinates of the MxN interpolated
        points. The units are meters.
    """
    m = start_points.shape[0]
    n = norm_spacings.size

    points = np.zeros((m, n, 3))

    for i in range(m):
        start_point = start_points[i, :]
        end_point = end_points[i, :]

        vector = end_point - start_point

        for j in range(n):
            norm_spacing = norm_spacings[j]

            spacing = norm_spacing * vector

            points[i, j, :] = start_point + spacing

    return points
