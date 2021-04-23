# ToDo: Properly document this module
import numpy as np
from numba import njit


def cosspace(
    minimum=0.0,
    maximum=1.0,
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

    :param minimum: float, optional
        This is the minimum value of the range of numbers you would like spaced. The
        default is 0.0.
    :param maximum: float, optional
        This is the maximum value of the range of numbers you would like spaced. The
        default is 1.0.
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
