# ToDo: Properly document this module
"""This module contains functions used by other modules in the src package.

"""
import logging

import numpy as np
from numba import njit, prange


def cosspace(
    minimum,
    maximum,
    n_points=50,
    endpoint=True,
):
    """This function is used to create a array containing a specified number of
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
        This sets whether or not the maximum value will be included in the output.
        The default is True.
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


def reflect_over_xz_plane(input_vector):
    """This function is used to flip a the y coordinate of a coordinate vector.

    Citation:
        Adapted from:         geometry.reflect_over_xz_plane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    :param input_vector: array
        This can either be a 1D array of three items, a M x 3 2D array, or a M x
        N x 3 3D array. N and
        represent arbitrary numbers of rows or columns.
    :return output vector: array
        This is a array with each vertex's y variable flipped.
    """

    # Initialize the output vector.
    output_vector = input_vector

    # Find the shape of the input vector.
    shape = np.shape(output_vector)

    # Characterize the input vector.
    if len(shape) == 1 and shape[0] == 3:
        # The input vector is a 1D array of 3 items. Flip the vertex's y variable.
        output_vector = output_vector * np.array([1, -1, 1])
    elif len(shape) == 2 and shape[1] == 3:
        # The input vector is a 2D array of shape M x 3. Where M is some arbitrary
        # number of rows. Flip each
        # vertex's y variable.
        output_vector = output_vector * np.array([1, -1, 1])
    elif len(shape) == 3 and shape[2] == 3:  # 3D MxNx3 vector
        # The input vector is a 3D array of shape M x N x 3. Where M is some
        # arbitrary number of rows, and N is
        # some arbitrary number of columns. Flip each vertex's y variable.
        output_vector = output_vector * np.array([1, -1, 1])
    else:
        # The input vector is an unacceptable shape. Throw an error.
        raise Exception("Invalid input for reflect_over_xz_plane.")

    # Return the output vector.
    return output_vector


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
                cos_theta + u_x ** 2 * (1 - cos_theta),
                u_x * u_y * (1 - cos_theta) - u_z * sin_theta,
                u_x * u_z * (1 - cos_theta) + u_y * sin_theta,
            ],
            [
                u_y * u_x * (1 - cos_theta) + u_z * sin_theta,
                cos_theta + u_y ** 2 * (1 - cos_theta),
                u_y * u_z * (1 - cos_theta) - u_x * sin_theta,
            ],
            [
                u_z * u_x * (1 - cos_theta) - u_y * sin_theta,
                u_z * u_y * (1 - cos_theta) + u_x * sin_theta,
                cos_theta + u_z ** 2 * (1 - cos_theta),
            ],
        ]
    )

    # Return the rotation matrix.
    return rotation_matrix


@njit(cache=True)
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


def calculate_streamlines(solver, num_steps=10, delta_time=0.1):
    """Calculates the location of the streamlines coming off the back of the wings.

    This method is vectorized to increase performance.

    :param solver: Solver
        This is the solver object for which to calculate the streamlines.
    :param num_steps: int, optional
        This is the integer number of points along each streamline (not including the
        initial points). It can be increased for higher fidelity visuals. The default
        value is 10.
    :param delta_time: float, optional
        This is the time in seconds between each time current_step It can be
        decreased for higher fidelity visuals or to make the streamlines shorter.
        It's default value is 0.1 seconds.
    :return: None
    """
    # Initialize a array to hold this solver's matrix of streamline points.
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


# ToDo: Update this function's documentation.
def convert_logging_level_name_to_value(name):
    """

    :param name:
    :return:
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


# ToDo: Update this function's documentation.
def process_steady_solver_forces(
    steady_solver, near_field_forces_geometry_axes, near_field_moments_geometry_axes
):
    """

    :param steady_solver:
    :param near_field_forces_geometry_axes:
    :param near_field_moments_geometry_axes:
    :return:
    """
    # Initialize a variable to hold the global panel position.
    global_panel_position = 0

    # Iterate through this solver's panels.
    for panel in steady_solver.panels:
        # Update the force and moment on this panel.
        panel.near_field_force_geometry_axes = near_field_forces_geometry_axes[
            global_panel_position, :
        ]
        panel.near_field_moment_geometry_axes = near_field_moments_geometry_axes[
            global_panel_position, :
        ]

        # Update the pressure on this panel.
        panel.update_pressure()

        # Increment the global panel position.
        global_panel_position += 1

    # Sum up the near field forces and moments on every panel to find the total
    # force and moment on the geometry.
    total_near_field_force_geometry_axes = np.sum(
        near_field_forces_geometry_axes, axis=0
    )
    total_near_field_moment_geometry_axes = np.sum(
        near_field_moments_geometry_axes, axis=0
    )

    # Find the total near field force in wind axes from the rotation matrix and
    # the total near field force in geometry axes.
    steady_solver.airplane.total_near_field_force_wind_axes = (
        np.transpose(
            steady_solver.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes()
        )
        @ total_near_field_force_geometry_axes
    )

    # Find the total near field moment in wind axes from the rotation matrix and
    # the total near field moment in geometry axes.
    steady_solver.airplane.total_near_field_moment_wind_axes = (
        np.transpose(
            steady_solver.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes()
        )
        @ total_near_field_moment_geometry_axes
    )

    # Calculate the current_airplane's induced drag coefficient
    induced_drag_coefficient = (
        -steady_solver.airplane.total_near_field_force_wind_axes[0]
        / steady_solver.operating_point.calculate_dynamic_pressure()
        / steady_solver.airplane.s_ref
    )

    # Calculate the current_airplane's side force coefficient.
    side_force_coefficient = (
        steady_solver.airplane.total_near_field_force_wind_axes[1]
        / steady_solver.operating_point.calculate_dynamic_pressure()
        / steady_solver.airplane.s_ref
    )

    # Calculate the current_airplane's lift coefficient.
    lift_coefficient = (
        -steady_solver.airplane.total_near_field_force_wind_axes[2]
        / steady_solver.operating_point.calculate_dynamic_pressure()
        / steady_solver.airplane.s_ref
    )

    # Calculate the current_airplane's rolling moment coefficient.
    rolling_moment_coefficient = (
        steady_solver.airplane.total_near_field_moment_wind_axes[0]
        / steady_solver.operating_point.calculate_dynamic_pressure()
        / steady_solver.airplane.s_ref
        / steady_solver.airplane.b_ref
    )

    # Calculate the current_airplane's pitching moment coefficient.
    pitching_moment_coefficient = (
        steady_solver.airplane.total_near_field_moment_wind_axes[1]
        / steady_solver.operating_point.calculate_dynamic_pressure()
        / steady_solver.airplane.s_ref
        / steady_solver.airplane.c_ref
    )

    # Calculate the current_airplane's yawing moment coefficient.
    yawing_moment_coefficient = (
        steady_solver.airplane.total_near_field_moment_wind_axes[2]
        / steady_solver.operating_point.calculate_dynamic_pressure()
        / steady_solver.airplane.s_ref
        / steady_solver.airplane.b_ref
    )

    steady_solver.airplane.total_near_field_force_coefficients_wind_axes = np.array(
        [induced_drag_coefficient, side_force_coefficient, lift_coefficient]
    )
    steady_solver.airplane.total_near_field_moment_coefficients_wind_axes = np.array(
        [
            rolling_moment_coefficient,
            pitching_moment_coefficient,
            yawing_moment_coefficient,
        ]
    )


# ToDo: Update this function's documentation.
def update_ring_vortex_solvers_panel_attributes(solver, global_panel_position, panel):
    """

    :param solver:
    :param global_panel_position:
    :param panel:
    :return:
    """
    # Update the solver's list of attributes with this panel's attributes.
    solver.panels[global_panel_position] = panel
    solver.panel_normal_directions[global_panel_position, :] = panel.normal_direction
    solver.panel_areas[global_panel_position] = panel.area
    solver.panel_centers[global_panel_position] = panel.center
    solver.panel_collocation_points[global_panel_position, :] = panel.collocation_point
    solver.panel_back_right_vortex_vertices[
        global_panel_position, :
    ] = panel.ring_vortex.right_leg.origin
    solver.panel_front_right_vortex_vertices[
        global_panel_position, :
    ] = panel.ring_vortex.right_leg.termination
    solver.panel_front_left_vortex_vertices[
        global_panel_position, :
    ] = panel.ring_vortex.left_leg.origin
    solver.panel_back_left_vortex_vertices[
        global_panel_position, :
    ] = panel.ring_vortex.left_leg.termination
    solver.panel_right_vortex_centers[
        global_panel_position, :
    ] = panel.ring_vortex.right_leg.center
    solver.panel_right_vortex_vectors[
        global_panel_position, :
    ] = panel.ring_vortex.right_leg.vector
    solver.panel_front_vortex_centers[
        global_panel_position, :
    ] = panel.ring_vortex.front_leg.center
    solver.panel_front_vortex_vectors[
        global_panel_position, :
    ] = panel.ring_vortex.front_leg.vector
    solver.panel_left_vortex_centers[
        global_panel_position, :
    ] = panel.ring_vortex.left_leg.center
    solver.panel_left_vortex_vectors[
        global_panel_position, :
    ] = panel.ring_vortex.left_leg.vector
    solver.panel_back_vortex_centers[
        global_panel_position, :
    ] = panel.ring_vortex.back_leg.center
    solver.panel_back_vortex_vectors[
        global_panel_position, :
    ] = panel.ring_vortex.back_leg.vector
    solver.panel_is_trailing_edge[global_panel_position] = panel.is_trailing_edge
    solver.panel_is_leading_edge[global_panel_position] = panel.is_leading_edge
    solver.panel_is_right_edge[global_panel_position] = panel.is_right_edge
    solver.panel_is_left_edge[global_panel_position] = panel.is_left_edge

    # Check if this panel is on the trailing edge.
    if panel.is_trailing_edge:
        # If it is, calculate it's streamline seed point and add it to
        # the solver's array of seed points.
        solver.seed_points = np.vstack(
            (
                solver.seed_points,
                panel.back_left_vertex
                + 0.5 * (panel.back_right_vertex - panel.back_left_vertex),
            )
        )


def calculate_steady_freestream_wing_influences(steady_solver):
    """This method finds the vector of freestream-wing influence coefficients
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


@njit(parallel=True, cache=True, fastmath=True)
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
    for i in prange(crosses.shape[0]):
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
