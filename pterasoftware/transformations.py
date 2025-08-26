"""This module contains functions use for geometric transformations.

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
"""

import numpy as np


def generate_homog(vector_A, has_point):
    """This function converts a 3D vector to homogeneous coordinates for use with
    4x4 transformation matrices.

    Homogeneous coordinates extend 3D vectors to 4D by adding a fourth component.
    For position vectors (points), the fourth component is 1. For direction
    vectors (such as velocity or force vectors), the fourth component is 0. This
    allows 4x4 transformation matrices to handle both translations and rotations
    in a unified framework.

    :param vector_A: (3,) ndarray of floats
        This is the 3D vector to convert to homogeneous coordinates. The units
        depend on the type of vector (meters for position, meters/second for
        velocity, etc.).
    :param has_point: bool
        If True, treats the vector as a position vector (point) and sets the
        fourth component to 1. If False, treats the vector as a direction vector
        and sets the fourth component to 0.
    :return: (4,) ndarray of floats
        This is the vector in homogeneous coordinates. For position vectors, the
        fourth component is 1. For direction vectors, the fourth component is 0.
        The units are the same as the input vector.
    """
    vector_A_homo = np.zeros(4, dtype=float)
    vector_A_homo[:3] = vector_A

    if has_point:
        vector_A_homo[-1] = 1

    return vector_A_homo


def generate_R(angles, passive, intrinsic, order):
    """This function generates either a passive or active rotation matrix using an
    angle vector, and parameters specifying if the matrix should be passive or
    active, intrinsic or extrinsic, and the order by which to apply the rotations
    about each axis.

    Notes: This function doesn't validate inputs or handle gimbal lock. Angles are in
    degrees with signs defined using the right-hand rule.

    Passive Use-Case: Let `r_A` be a vector in "A" axes, but we want to find `r_B`,
    which is the same vector, but expressed in "B" axes. The orientation of "B" axes
    relative to "A" axes is defined by the angle vector `angles` (with rotations in
    order and type defined by the variables `order` and `intrinsic`). Then:
    `R_pas_A_to_B = generate_R(angles, True, intrinsic, order)` and `r_B =
    R_pas_A_to_B @ r_A`.

    Active Use-Case: Let `r_A` be a vector in "A" axes, but we want to find
    `r_A_prime`, which is `r_A` rotated by the specified sequence: about the fixed
    "A" axes if `intrinsic=False`, or about the current, newly-rotated axes if
    `intrinsic=True`, with angles given by `angles` and the sequence defined by
    `order`. Then: `R_act = generate_R(angles, False, intrinsic, order)` and
    `r_A_prime = R_act @ r_A`.

    :param angles: (3,) ndarray of floats
        This is the angle vector in with signs defined using the right-hand rule. For
        `passive=True`, it describes the orientation of "B" axes with respect to "A"
        axes. For `passive=False`, it prescribes the angles by which to rotate a
        vector in "A" axes. In both cases, the rotations' type is specified by the
        intrinsic parameter. Angles are always listed as [about axis 1, about axis 2,
        about axis 3], but are applied in the sequence given by `order` (e.g.,
        order="312" applies angles[2], angles[0], angles[1]).
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
        Each character can be '1', '2', or '3'. Only Tait-Bryan angles are accepted
        so all accepted characters must be distinct.
    :return: (3, 3) ndarray of floats
        This is the rotation matrix.
    """
    angles = np.array(angles, dtype=float)

    angle_1_rad, angle_2_rad, angle_3_rad = np.radians(angles)

    R_act_1 = np.array(
        [
            [1, 0, 0],
            [0, np.cos(angle_1_rad), -np.sin(angle_1_rad)],
            [0, np.sin(angle_1_rad), np.cos(angle_1_rad)],
        ]
    )
    R_act_2 = np.array(
        [
            [np.cos(angle_2_rad), 0, np.sin(angle_2_rad)],
            [0, 1, 0],
            [-np.sin(angle_2_rad), 0, np.cos(angle_2_rad)],
        ]
    )
    R_act_3 = np.array(
        [
            [np.cos(angle_3_rad), -np.sin(angle_3_rad), 0],
            [np.sin(angle_3_rad), np.cos(angle_3_rad), 0],
            [0, 0, 1],
        ]
    )

    R_act_components = [R_act_1, R_act_2, R_act_3]

    order_nums = [ord(order_num) - ord("0") for order_num in order]

    if intrinsic:
        order_nums.reverse()

    order_ids = [order_num - 1 for order_num in order_nums]

    R_act = np.eye(3)
    for order_id in order_ids:
        R_act = R_act_components[order_id] @ R_act

    if passive:
        return R_act.T
    return R_act


def generate_T_rot(R):
    """This function converts a rotation matrix to a rotational transformation matrix
    for use with homogeneous coordinates.

    Notes: This function doesn't validate inputs.

    :param (3, 3) ndarray of floats
        This is the rotation matrix. It can be active or passive, intrinsic or
        extrinsic, and perform rotations in any order. The generated transformation
        matrix will retain these same attributes.
    :return: (4, 4) ndarray of floats
        This is the transformation matrix.
    """
    T_rot = np.eye(4)
    T_rot[:3, :3] = R
    return T_rot


def generate_T_trans(translations, passive):
    """This function generates either a passive or active translational
    transformation matrix using a vector and a parameter specifying if the
    transformation should be passive or active.

    Passive Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c", (in "A" axes, relative to the "a" point). We want to find `c_A_b`,
    which describes the location of "c", relative to the "b" point. The position of
    "b" is defined by 'translations' (in "A" axes, relative to the point a). Then:
    `T_trans_pas_a_to_b = generate_T_trans(translate, True)`, `c_A_b_homo =
    T_trans_pas_a_to_b @ generate_homog(c_A_b, True)`, and `c_A_b = c_A_b_homo[:3]`.

    Active Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `cPrime_A_a`,
    which is the position of "cPrime", which is point "c" offset by `translations` (
    in "A" axes). Then: `T_trans_act = generate_T_trans( translations, False)`,
    `cPrime_A_a_homo = T_trans_act @ generate_homog(c_A_a, True)`, and `cPrime_A_a =
    cPrime_A_a_homo[:3]`.

    :param translations: (3,) ndarray of floats
        For `passive=True`, this is the position of the "b" point (in "A" axes,
        relative to the "a" point). For `passive=False`, this is the position (in "A"
        axes) of the offset point "cPrime" relative to the original "c" point.
    :param passive: bool
        If True, returns a matrix that changes reference point of a vector in
        homogeneous coordinates (`r_A_b_homo = T_trans @ r_A_a_homo`). If False,
        returns the position vector of point offset from an original position,
        still relative to the same point (`cPrime_A_a_homo = T_trans @ c_A_a_homo`).
    :return: (4, 4) ndarray of floats
        This is the transformation matrix.
    """
    T_trans = np.eye(4)
    T_trans[:3, 3] = -translations if passive else translations

    return T_trans
