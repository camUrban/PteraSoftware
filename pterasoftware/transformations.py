"""This module contains functions used for geometric transformations.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    generate_homog: This function converts a 3D vector to homogeneous coordinates
    for use with 4x4 transformation matrices.

    generate_R: This function generates either a passive or active rotation matrix
    using an angle vector, and parameters specifying if the matrix should be passive
    or active, intrinsic or extrinsic, and the order by which to apply the rotations
    about each axis.

    generate_T_rot: This function converts a rotation matrix to a rotational
    transformation matrix for use with homogeneous coordinates.

    generate_T_trans: This function generates either a passive or active
    translational transformation matrix using a vector and a parameter specifying if
    the transformation should be passive or active.

    generate_T_reflect: This function generates either a passive or active
    reflectional transformation matrix about a plane which is defined by a point (in
    "A" axes, relative to point "a") and a normal unit vector (in "A" axes)."""

import numpy as np

from . import parameter_validation


def generate_homog(vector_A, has_point):
    """This function converts a 3D vector to homogeneous coordinates for use with
    4x4 transformation matrices.

    Homogeneous coordinates extend 3D vectors to 4D by adding a fourth component. For
    position vectors (points), the fourth component is 1. For free vectors (such as
    velocity or force vectors), the fourth component is 0. This allows 4x4
    transformation matrices to handle both translations and rotations in a unified
    framework.

    :param vector_A: array-like of 3 numbers

        This is the 3D vector to convert to homogeneous coordinates. The units depend
        on the type of vector (meters for position, meters/second for velocity,
        etc.). Can be a tuple, list, or numpy array of numbers (int or float). Values
        are converted to floats internally.

    :param has_point: bool

        If True, treats the vector as a position vector (point) and sets the fourth
        component to 1. If False, treats the vector as a free vector and sets the
        fourth component to 0.

    :return: (4,) ndarray of floats

        This is the vector in homogeneous coordinates. For position vectors,
        the fourth component is 1. For free vectors, the fourth component is 0. The
        units are the same as the input vector.
    """
    vector_A = parameter_validation.validate_3d_vector_float(vector_A, "vector_A")
    has_point = parameter_validation.validate_boolean(has_point, "has_point")

    vectorHomog_A = np.zeros(4, dtype=float)

    vectorHomog_A[:3] = vector_A

    if has_point:
        vectorHomog_A[-1] = 1.0

    return vectorHomog_A


def generate_R(angles, passive, intrinsic, order):
    """This function generates either a passive or active rotation matrix using an
    angle vector, and parameters specifying if the matrix should be passive or
    active, intrinsic or extrinsic, and the order by which to apply the rotations
    about each axis.

    Notes: Angles are in degrees with signs defined using the right-hand rule.

    Passive Use-Case: Let `r_A` be a vector in "A" axes, but we want to find `r_B`,
    which is the same vector, but expressed in "B" axes. The orientation of "B" axes
    relative to "A" axes is defined by the angle vector `angles` (with rotations in
    order and type defined by the variables `order` and `intrinsic`). Then:

    ```
     R_pas_A_to_B = generate_R(angles, True, intrinsic, order)
    r_B =  R_pas_A_to_B @ r_A
    ```

    Active Use-Case: Let `r_A` be a vector in "A" axes, but we want to find
    `r_A_prime`, which is `r_A` rotated by the specified sequence: about the fixed
    "A" axes if `intrinsic=False`, or about the current, newly-rotated axes if
    `intrinsic=True`, with angles given by `angles` and the sequence defined by
    `order`. Then:

    ```
    rotate_R_act = generate_R(angles, False, intrinsic, order)
    r_A_prime = rotate_R_act @ r_A
    ```

    :param angles: array-like of 3 numbers

        This is the angle vector, in degrees, with signs defined using the right-hand
        rule. Can be a tuple, list, or numpy array of numbers (int or float). Values
        are converted to floats internally. For `passive=True`, it describes the
        orientation of "B" axes with respect to "A" axes. For `passive=False`,
        it prescribes the angles by which to rotate a vector in "A" axes. In both
        cases, the rotations' type is specified by the intrinsic parameter. Angles
        are always listed as [about x-axis, about y-axis, about z-axis], but are
        applied in the sequence given by `order` (e.g., order="zxy" applies angles[
        2], angles[0], angles[1]).

    :param passive: bool

        If True, returns a matrix that changes coordinates from "A" to "B" axes (`r_B
        = R @ r_A`). If False, returns a matrix that rotates vectors in "A" axes (
        `r_A_prime = R @ r_A`).

    :param intrinsic: bool

        If True, each subsequent rotation is applied to the current, newly-rotated
        axes. If False, all rotations are performed about the original, non-rotated
        "A" axes.

    :param order: string of length 3

        This is the string of three characters that represents the rotation order.
        Each character can be 'x', 'y', or 'z'. Only Tait-Bryan angles are accepted
        so all accepted characters must be distinct.

    :return: (3, 3) ndarray of floats

        This is the rotation matrix.
    """
    angles = parameter_validation.validate_3d_vector_float(angles, "angles")
    passive = parameter_validation.validate_boolean(passive, "passive")
    intrinsic = parameter_validation.validate_boolean(intrinsic, "intrinsic")
    order = parameter_validation.validate_rotation_order(order, "order")

    angleX_rad, angleY_rad, angleZ_rad = np.radians(angles)

    R_act_x = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(angleX_rad), -np.sin(angleX_rad)],
            [0.0, np.sin(angleX_rad), np.cos(angleX_rad)],
        ]
    )
    R_act_y = np.array(
        [
            [np.cos(angleY_rad), 0.0, np.sin(angleY_rad)],
            [0.0, 1.0, 0.0],
            [-np.sin(angleY_rad), 0.0, np.cos(angleY_rad)],
        ]
    )
    R_act_z = np.array(
        [
            [np.cos(angleZ_rad), -np.sin(angleZ_rad), 0],
            [np.sin(angleZ_rad), np.cos(angleZ_rad), 0],
            [0.0, 0.0, 1.0],
        ]
    )

    R_act_components = [R_act_x, R_act_y, R_act_z]

    order_nums = [{"x": 1, "y": 2, "z": 3}[order_char] for order_char in order]

    if intrinsic:
        order_nums.reverse()

    order_ids = [order_num - 1 for order_num in order_nums]

    R_act = np.eye(3, dtype=float)
    for order_id in order_ids:
        R_act = R_act_components[order_id] @ R_act

    if passive:
        return R_act.T
    return R_act


def generate_T_rot(R):
    """This function converts a rotation matrix to a rotation transformation matrix
    using homogeneous coordinates.

    :param R: 3x3 array-like of numbers

        This is the rotation matrix. It can be any 3x3 array-like object and contain
        either floats or ints. Internally, it will be converted to a (3, 3) ndarray
        of floats. The rotation matrix can be active or passive, intrinsic or
        extrinsic, and perform rotations in any order. The generated transformation
        matrix will retain these same attributes.

    :return: (4, 4) ndarray of floats

        This is the transformation matrix.
    """
    R = parameter_validation.validate_3_by_3_matrix_float(R, "R")
    T_rot = np.eye(4, dtype=float)

    T_rot[:3, :3] = R
    return T_rot


def generate_T_trans(translations, passive):
    """This function generates either a passive or active translational
    transformation matrix using a vector and a parameter specifying if the
    transformation should be passive or active.

    Passive Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `c_A_b`,
    which describes the location of "c", relative to the "b" point. The position of
    "b" is defined by `translations` (in "A" axes, relative to the point a). Then:

    ```
    T_pas_A_a_to_A_b = generate_T_trans(translations, True)
    cHomog_A_b = T_pas_A_a_to_A_b @ generate_homog(c_A_a, True)
    c_A_b = cHomog_A_b[:3]
    ```

    Active Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `cPrime_A_a`,
    which is the position of "cPrime", which is point "c" offset by `translations` (
    in "A" axes). Then:

    ```
    translate_T_act = generate_T_trans(translations, False)
    cPrimeHomog_A_a = translate_T_act @ generate_homog(c_A_a, True)
    cPrime_A_a = cPrimeHomog_A_a[:3]
    ```

    :param translations: array-like of 3 numbers

        Can be a tuple, list, or numpy array of numbers (int or float). Values are
        converted to floats internally. For `passive=True`, this is the position of
        the "b" point (in "A" axes, relative to the "a" point). For `passive=False`,
        this is the position (in "A" axes) of the offset point "cPrime" relative to
        the original "c" point.

    :param passive: bool

        If True, returns a matrix that changes reference point of a vector in
        homogeneous coordinates (`rHomog_A_b = T_trans @ rHomog_A_a`). If False,
        returns the position vector of point offset from an original position,
        still relative to the same point (`cPrimeHomog_A_a = T_trans @ cHomog_A_a`).

    :return: (4, 4) ndarray of floats

        This is the transformation matrix.
    """
    p = parameter_validation.validate_3d_vector_float(translations, "translations")
    passive = parameter_validation.validate_boolean(passive, "passive")
    T_trans = np.eye(4, dtype=float)

    T_trans[:3, 3] = -p if passive else p
    return T_trans


def generate_T_reflect(plane_point_A_a, plane_normal_A, passive):
    """This function generates either a passive or active reflectional transformation
    matrix about a plane which is defined by a point (in "A" axes, relative to point
    "a") and a normal vector (in "A" axes).

    Note: This function generates identical matrices for both passive and active
    cases, which is correct. However, it retains the `passive` flag for API
    consistency and also as a reminder to think about what the final matrix represents.

    Passive Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `c_B_b`,
    which describes the location of "c" in "B" axes, relative to the "b" point. The
    orientation of "B" is "A" reflected across the plane defined by plane_point_A_a
    and plane_normal_A. The "b" point is located at the "a" point's position,
    reflected across the same plane. Then:

    ```
    T_pas_A_a_to_B_b = generate_T_reflect(plane_point_A_a, plane_normal_A, True)
    cHomog_B_b = T_pas_A_a_to_B_b @ generate_homog(c_A_a, True)
    c_B_b = cHomog_B_b[:3]
    ```

    Active Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `cPrime_A_a`,
    which is the position of "cPrime", point "c" reflected across the plane defined
    by plane_point_A_a and plane_normal_A. Then:

    ```
    reflect_T_act = generate_T_reflect(plane_point_A_a, plane_normal_A, False)
    cPrimeHomog_A_a = reflect_T_act @ generate_homog(c_A_a, True)
    cPrime_A_a = cPrimeHomog_A_a[:3]
    ```

    :param plane_point_A_a: array-like of 3 numbers

        This is a point on the reflection plane (in "A" axes, relative to the "a"
        point). Can be a tuple, list, or numpy array of numbers (int or float).
        Values are converted to floats internally. It is in units of meters.

    :param plane_normal_A: array-like of 3 numbers

        This is the normal vector of the reflection plane (in "A" axes). It is
        normalized internally. Can be a tuple, list, or numpy array of numbers (int
        or float). Values are converted to floats internally. It must have a non-zero
        length.

    :param passive: bool

        If True, returns a matrix that changes reference point and axes of a vector
        in homogeneous coordinates to a reference point and axes reflected about the
        specified plane (`cHomog_B_b = T_reflect @ cHomog_A_a`). If False, it returns
        a matrix that reflects a vector (in its original axes) about a specified
        plane (`cPrimeHomog_A_a = T_reflect @ cHomog_A_a`).

    :return: (4, 4) ndarray of floats

        This is the transformation matrix.
    """
    p = parameter_validation.validate_3d_vector_float(
        plane_point_A_a, "plane_point_A_a"
    )
    n_hat = parameter_validation.validate_3d_unit_vector_norm_float(
        plane_normal_A, "plane_normal_A"
    )
    passive = parameter_validation.validate_boolean(passive, "passive")

    T_reflect = np.eye(4, dtype=float)

    S = np.eye(3, dtype=float) - 2 * np.outer(n_hat, n_hat)
    d = 2 * (np.dot(p, n_hat)) * n_hat

    T_reflect[:3, :3] = S
    T_reflect[:3, 3] = d

    return T_reflect
