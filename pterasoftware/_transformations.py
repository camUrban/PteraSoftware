"""Contains shared functions for geometric transformations."""

from __future__ import annotations

import numpy as np

# Threshold below which R_to_quat_wxyz's per component branch checks trigger the
# catastrophic-cancellation-resistant alternative formulas for extracting the
# corresponding quaternion component. Sarabandi and Thomas (2019) Section 4 Table 1
# identifies eta = 0 as the experimentally optimal value across worst-case error,
# average error, and standard deviation.
_QUAT_BRANCH_THRESHOLD = 0.0

# Threshold above which |sin(angleY)| triggers the gimbal-lock fallback decomposition
# in R_to_angles_izyx, equivalent to angleY within about 8.1e-5 degrees of the +/- 90
# degree pole. Picked to sit four orders of magnitude above float64 epsilon (1e-16): the
# standard atan2 formula keeps cos(angleY) as a shared factor that cancels in the ratio,
# so it works correctly until cos(angleY) approaches machine epsilon (cos at this
# threshold is sqrt(2e-12) ~ 1.4e-6). The worst-case round-trip angle error just below
# the threshold is roughly 1 nanodegree.
_GIMBAL_LOCK_THRESHOLD = 1.0 - 1.0e-12


def _generate_homogs(vectors_A: np.ndarray, is_position: bool) -> np.ndarray:
    """Converts 3D vector(s) to homogeneous coordinates for use with (4,4)
    transformation matrices.

    Homogeneous coordinates extend 3D vectors to 4D by adding a fourth component. For
    position vectors (is_position=True), which represent a point, the fourth component
    is 1.0 so the transform's translation applies. Otherwise (is_position=False, for
    example a velocity, force, or moment), the fourth component is 0.0 so only the
    rotation applies. This allows (4,4) transformation matrices to handle both
    translations and rotations in a unified framework.

    This function handles both single vectors and arrays of vectors efficiently.

    :param vectors_A: A (...,3) ndarray of floats representing the vector(s) to convert
        to homogeneous coordinates.
    :param is_position: True if the vector is a position (a point), False otherwise (for
        example a velocity, force, or moment).
    :return: A (...,4) ndarray of floats (with same leading dimensions as the input)
        representing the vector(s) in homogeneous coordinates.
    """
    # Create a homogeneous ndarray with one extra dimension.
    vectorsHomog_A = np.zeros(vectors_A.shape[:-1] + (4,), dtype=float)

    # Copy the vectors' three components to the homogeneous ndarray.
    vectorsHomog_A[..., :3] = vectors_A

    # Set the homogeneous coordinate.
    if is_position:
        vectorsHomog_A[..., -1] = 1.0

    return vectorsHomog_A


def generate_rot_T(
    angles: np.ndarray,
    passive: bool,
    intrinsic: bool,
    order: str,
) -> np.ndarray:
    """Generates a rotational transformation matrix.

    **Passive Use-Case:**

    Let ``r_A`` be a non-position vector in "A" axes, but we want to find ``r_B``, which
    is the same vector, but expressed in "B" axes. The orientation of "B" axes relative
    to "A" axes is defined by the angle vector ``angles`` (with rotations in order and
    type defined by the variables ``order`` and ``intrinsic``). Then:

    | ``T_pas_A_to_B=generate_rot_T(angles,True,intrinsic,order)``

    | ``r_B=apply_T_to_vectors(T_pas_A_to_B,r_A,is_position=False)``

    **Active Use-Case:**

    Let ``r_A`` be a non-position vector in "A" axes, but we want to find ``rPrime_A``,
    which is ``r_A`` rotated by the specified sequence: about the fixed "A" axes if
    ``intrinsic=False``, or about the current, newly-rotated axes if ``intrinsic=True``,
    with angles given by ``angles`` and the sequence defined by ``order``. Then:

    | ``rot_T_act=generate_rot_T(angles,False,intrinsic,order)``

    | ``rPrime_A=apply_T_to_vectors(rot_T_act,r_A,is_position=False)``

    :param angles: A (3,) ndarray of floats representing the rotation angles, with signs
        defined using the right-hand rule. For `passive=True`, it describes the
        orientation of "B" axes with respect to "A" axes. For `passive=False`, it
        prescribes the angles by which to rotate a vector in "A" axes. In both cases,
        the rotations' type is specified by the intrinsic parameter. Angles are always
        listed as [about x axis, about y axis, about z axis], but are applied in the
        sequence given by `order` (e.g., order="zxy" applies angles[2], angles[0],
        angles[1]). The units are in degrees.
    :param passive: Set this to True to return a matrix that changes coordinates from
        "A" to "B" axes (``r_B=R@r_A``). Set this to False to return a matrix that
        rotates vectors in "A" axes (``rPrime_A=R@r_A``).
    :param intrinsic: Set this to True to return a transformation matrix where each
        subsequent rotation is applied to the current, newly-rotated axes. Set this to
        False to return a transformation matrix where rotations are performed about the
        original, non rotated "A" axes.
    :param order: A str of three chars that represents the rotation order. Each char can
        be 'x', 'y', or 'z'. Only Tait-Bryan angles are accepted so all accepted chars
        must be distinct.
    :return: The transformation matrix as a (4,4) ndarray of floats.
    """
    angleXRad, angleYRad, angleZRad = np.deg2rad(angles)

    x_R_act = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(angleXRad), -np.sin(angleXRad)],
            [0.0, np.sin(angleXRad), np.cos(angleXRad)],
        ]
    )
    y_R_act = np.array(
        [
            [np.cos(angleYRad), 0.0, np.sin(angleYRad)],
            [0.0, 1.0, 0.0],
            [-np.sin(angleYRad), 0.0, np.cos(angleYRad)],
        ]
    )
    z_R_act = np.array(
        [
            [np.cos(angleZRad), -np.sin(angleZRad), 0.0],
            [np.sin(angleZRad), np.cos(angleZRad), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    R_act_components = [x_R_act, y_R_act, z_R_act]

    order_nums = [{"x": 1, "y": 2, "z": 3}[order_char] for order_char in order]

    if intrinsic:
        order_nums.reverse()

    order_ids = [order_num - 1 for order_num in order_nums]

    R_act = np.eye(3, dtype=float)
    for order_id in order_ids:
        R_act = R_act_components[order_id] @ R_act

    R = R_act
    if passive:
        R = R.T

    T = np.eye(4, dtype=float)

    T[:3, :3] = R
    return T


def compute_offset_rotation_adjustment(
    rotation_matrix: np.ndarray,
    offset: np.ndarray,
) -> np.ndarray:
    """Compute the positional adjustment required to rotate about an offset point.

    When a rotation is applied about a point other than the origin, the effective
    position must be adjusted. This function computes that adjustment vector.

    The adjustment is given by: (I - R) @ offset

    where R is the active rotation matrix and offset is the rotation point.

    :param rotation_matrix: A (3,3) ndarray representing the active rotation matrix.
    :param offset: A (3,) ndarray representing the rotation point offset.
    :return: A (3,) ndarray representing the position adjustment.
    """
    return np.asarray((np.eye(3, dtype=float) - rotation_matrix) @ offset, dtype=float)


def generate_2D_rot_R(
    angle: float,
    passive: bool,
) -> np.ndarray:
    """Generates a 2D rotational matrix.

    **Passive Use-Case:**

    Let ``r_A`` be a 2D non-position vector in "A" axes, but we want to find ``r_B``,
    which is the same vector, but expressed in "B" axes. The orientation of "B" axes
    relative to "A" axes is defined by the angle ``angle``. Then:

    | ``R_pas_A_to_B=generate_2D_rot_R(angle,True)``

    | ``r_B=R_pas_A_to_B@r_A``

    **Active Use-Case:**

    Let ``r_A`` be a 2D non-position vector in "A" axes, but we want to find
    ``rPrime_A``, which is ``r_A`` rotated by ``angle``. Then:

    | ``rot_R_act=generate_2D_rot_R(angle,False)``

    | ``rPrime_A=rot_R_act@r_A``

    :param angle: A float representing the rotation angle, with signs defined using the
        right-hand rule. For ``passive=True``, it describes the orientation of "B" axes
        with respect to "A" axes. For ``passive=False``, it prescribes the angle by
        which to rotate a vector in "A" axes. The units are in degrees.
    :param passive: Set this to True to return a matrix that changes coordinates from
        "A" to "B" axes (``r_B=R@r_A``). Set this to False to return a matrix that
        rotates vectors in "A" axes (``rPrime_A=R@r_A``).
    :return: The rotation matrix as a (2,2) ndarray of floats.
    """
    angleRad = np.deg2rad(angle)

    R_act = np.array(
        [
            [np.cos(angleRad), -np.sin(angleRad)],
            [np.sin(angleRad), np.cos(angleRad)],
        ],
        dtype=float,
    )

    if passive:
        return R_act.T

    return R_act


def generate_trans_T(
    translations: np.ndarray,
    passive: bool,
) -> np.ndarray:
    """Generates a translational transformation matrix.

    **Passive Use-Case:**

    Let ``c_A_a`` be a vector which describes the location of point "c" (in "A" axes,
    relative to the "a" point). We want to find ``c_A_b``, which describes the location
    of "c", relative to the "b" point. The position of "b" is defined by
    ``translations`` (in "A" axes, relative to the point a). Then:

    | ``T_pas_A_a_to_A_b=generate_trans_T(translations,True)``

    | ``c_A_b=apply_T_to_vectors(T_pas_A_a_to_A_b,c_A_a,is_position=True)``

    **Active Use-Case:**

    Let ``c_A_a`` be a vector which describes the location of point "c" (in "A" axes,
    relative to the "a" point). We want to find ``cPrime_A_a``, which is the position of
    "cPrime", which is point "c" offset by `translations` (in "A" axes). Then:

    | ``translate_T_act=generate_trans_T(translations,False)``

    | ``cPrime_A_a=apply_T_to_vectors(translate_T_act,c_A_a,is_position=True)``

    :param translations: A (3,) ndarray of floats representing the translations. For
        ``passive=True``, this is the position of the "b" point (in "A" axes, relative
        to the "a" point). For ``passive=False``, this is the position (in "A" axes) of
        the offset point "cPrime" relative to the original "c" point. The units are in
        meters.
    :param passive: Set this to True to return a matrix that changes the reference point
        of a vector in homogeneous coordinates (``rHomog_A_b=T_trans@rHomog_A_a``). Set
        this to False to return a matrix that finds the new position vector of point
        after translating it from its original position
        (``cPrimeHomog_A_a=T_trans@cHomog_A_a``).
    :return: The transformation matrix as a (4,4) ndarray of floats.
    """
    p = translations
    T_trans = np.eye(4, dtype=float)
    T_trans[:3, 3] = -p if passive else p
    return T_trans


# noinspection PyUnusedLocal
def generate_reflect_T(
    plane_point_A_a: np.ndarray,
    plane_normal_A: np.ndarray,
    passive: bool,
) -> np.ndarray:
    """Generates a reflectional transformation matrix about a plane defined by a point
    (in "A" axes, relative to point "a") and a normal vector (in "A" axes).

    **Passive Use-Case:**

    Let ``c_A_a`` be a vector which describes the location of point "c" (in "A" axes,
    relative to the "a" point). We want to find ``c_B_b``, which describes the location
    of "c" in "B" axes, relative to the "b" point. The orientation of "B" is "A"
    reflected across the plane defined by ``plane_point_A_a`` and ``plane_normal_A``.
    The "b" point is located at the "a" point's position, reflected across the same
    plane. Then:

    | ``T_pas_A_a_to_B_b=generate_reflect_T(plane_point_A_a,plane_normal_A,True)``

    | ``c_B_b=apply_T_to_vectors(T_pas_A_a_to_B_b,c_A_a,is_position=True)``

    **Active Use-Case:**

    Let ``c_A_a`` be a vector which describes the location of point "c" (in "A" axes,
    relative to the "a" point). We want to find ``cPrime_A_a``, which is the position of
    "cPrime", point "c" reflected across the plane defined by ``plane_point_A_a`` and
    ``plane_normal_A``. Then:

    | ``reflect_T_act=generate_reflect_T(plane_point_A_a,plane_normal_A,False)``

    | ``c_A_a=apply_T_to_vectors(reflect_T_act,c_A_a,is_position=True)``

    **Notes:**

    This function generates identical matrices for both passive and active cases, which
    is correct. However, it retains the `passive` flag for API consistency and as a
    reminder to consider what the final matrix represents.

    **Warning:**

    A reflection is an improper transform. Pseudovectors (axial vectors, such as
    moments, angular velocities, vorticities, or surface normals) do not reflect the way
    position, velocity, or force vectors do: under a reflection a pseudovector picks up
    an extra sign change. Applying this matrix to a pseudovector with
    apply_T_to_vectors() therefore returns the negative of the correctly reflected
    vector.

    :param plane_point_A_a: A (3,) ndarray of floats representing a point on the
        reflection plane (in "A" axes, relative to the "a" point). The units are in
        meters.
    :param plane_normal_A: A (3,) ndarray of floats representing a vector (in "A" axes)
        normal to the reflection plane. It must have a non zero magnitude, and will be
        normalized to a unit vector.
    :param passive: Set this to True to return a matrix that changes reference point and
        axes of a vector in homogeneous coordinates to a reference point and axes
        reflected about the specified plane (``cHomog_B_b=T_reflect@cHomog_A_a``). Set
        this to False to return a matrix that reflects a vector (in its original axes,
        relative to its original reference point) about a specified plane
        (``cPrimeHomog_A_a=T_reflect@cHomog_A_a``).
    :return: The transformation matrix as a (4,4) ndarray of floats.
    """
    p = plane_point_A_a
    norm = np.linalg.norm(plane_normal_A)
    if norm == 0:
        raise ValueError("plane_normal_A must have a non zero length.")
    n_hat = plane_normal_A / norm
    T_reflect = np.eye(4, dtype=float)

    S = np.eye(3, dtype=float) - 2 * np.outer(n_hat, n_hat)
    d = 2 * (np.dot(p, n_hat)) * n_hat

    T_reflect[:3, :3] = S
    T_reflect[:3, 3] = d

    return T_reflect


def _left_compose_T(valid_T_chain: list[np.ndarray]) -> np.ndarray:
    """Left-compose a list of homogeneous transformations.

    :param valid_T_chain: A list of ndarrays of floats, each with shape (4,4),
        representing the series of transformations.
    :return: A single (4,4) ndarray of floats representing the composed transformation.
        For example, if ``valid_T_chain=[T_1,T_2,...,T_n]``, this function will return
        ``T_n@...@T_2@T_1``.
    """
    if len(valid_T_chain) == 1:
        return valid_T_chain[0]
    composed_T = np.eye(4, dtype=float)
    for valid_T in valid_T_chain:
        composed_T = valid_T @ composed_T
    return composed_T


def compose_T_pas(
    *T_pas_chain: np.ndarray,
) -> np.ndarray:
    """Compose a chain of passive homogeneous transformations.

    **Use-Case:**

    | ``T_pas_A_a_to_C_c=compose_T_pas(T_pas_A_a_to_B_b,T_pas_B_b_to_C_c)``

    :param T_pas_chain: One or more (4,4) ndarrays of floats representing the passive
        homogeneous transforms along the path from the original axes and reference point
        to the final axes and reference point.
    :return: The composed transformation matrix as a (4,4) ndarray of floats.
    """
    if not T_pas_chain:
        raise ValueError("At least one transform must be provided.")

    return _left_compose_T(list(T_pas_chain))


def compose_T_act(
    *T_act_chain: np.ndarray,
) -> np.ndarray:
    """Compose a chain of active homogeneous transformations.

    **Use-Case:**

    | ``composed_T_act=compose_T_act(reflect_T_act,rot_T_act,trans_T_act)``

    **Notes:**

    This function left-composes the supplied active transforms: given
    ``compose_T_act(T1,T2,...,Tn)`` it returns ``Tn@...@T2@T1``. Interpreting these as
    active transformations, this implies that they occur in the order in which they are
    passed.

    Active translations created with ``generate_trans_T(...,passive=False)`` interpret
    the components in the same axes the vector is expressed in (e.g., geometry axes).
    Therefore: ```T_act=compose_T_act(rot_T_act,trans_T_act)``` applies a rotation first
    and then a *world-fixed* translation. This can seem counter-intuitive.

    If you instead want a *body-fixed* translation (e.g., "+10 along x' after the
    rotation"), either pre-rotate the components before building the translation:

    | ``R=rot_T_act[:3,:3]``

    | ``tPrime_A=R@t_A``

    | ``transBodyFixed_T_act=generate_trans_T(tPrime_A,False)``

    | ``T_act=compose_T_act(rot_T_act,transBodyFixed_T_act)``

    or, pass the translation before the rotation:

    | ``T_act=compose_T_act(trans_T_act,rot_T_act)``

    :param T_act_chain: One or more (4,4) ndarrays of floats representing the active
        homogeneous transforms applied in order.
    :return: The composed transformation matrix as a (4,4) ndarray of floats.
    """
    if not T_act_chain:
        raise ValueError("At least one transform must be provided.")

    return _left_compose_T(list(T_act_chain))


def _invert_T_rigid(valid_T: np.ndarray) -> np.ndarray:
    """Invert a rigid homogeneous transform.

    **Notes:**

    A valid rigid homogeneous transform can be broken down into two components:

    | ``R=valid_T[:3,:3]``

    | ``t=valid_T[:3,3]``

    This function uses these components to return the inverse of ``valid_T``:

    | ``[[R.T,-R.T@t];[0,1]]``

    :param valid_T: A (4,4) ndarray of floats representing a valid, rigid homogeneous
        transform.
    :return: A (4,4) ndarray of floats representing the inverse of the input transform.
    """
    valid_R = valid_T[:3, :3]
    valid_t = valid_T[:3, 3]

    valid_transpose_R = valid_R.T

    valid_T_inv = np.eye(4, dtype=float)

    valid_T_inv[:3, :3] = valid_transpose_R
    valid_T_inv[:3, 3] = -valid_transpose_R @ valid_t

    return valid_T_inv


def invert_T_pas(T_pas: np.ndarray) -> np.ndarray:
    """Inverts a passive homogeneous transform.

    A passive transform maps components of the same physical quantity between an initial
    axis system and reference point and a target axis system and reference point. For
    example, if ``T_pas_A_a_to_B_b`` maps components from "A" axes (relative to point
    "a") to "B" axes (relative to point "b"), then:

    | ``T_pas_B_b_to_A_a=invert_T_pas(T_pas_A_a_to_B_b)``

    | ``rHomog_A_a=T_pas_B_b_to_A_a@rHomog_B_b``

    **Notes:**

    For position vectors (``is_position=True``), the translation component matters.
    Otherwise (``is_position=False``, for example a velocity, force, or moment),
    translation has no effect because the homogeneous last coordinate is 0.0.

    :param T_pas: A (4,4) ndarray of floats representing a passive homogeneous transform
        mapping from source axes and reference point to target axes and reference point.
    :return: A (4,4) ndarray of floats representing the passive transform that maps back
        from the target axes and reference point to the original axes and reference
        point.
    """
    return _invert_T_rigid(T_pas)


def invert_T_act(T_act: np.ndarray) -> np.ndarray:
    """Inverts an active homogeneous transform.

    An active transform re-orients and optionally translates a quantity within the same
    axis system. For example, if ``T_act`` transforms the non-position vector ``q_A``
    (in "A" axes) to the non-position vector ``qPrime_A`` (in "A" axes), then:

    | ``q_A=apply_T_to_vectors(invert_T_act(T_act),qPrime_A,is_position=False)``

    **Notes:**

    For position vectors (``is_position=True``), both orientation and translation are
    undone. Otherwise (``is_position=False``, for example a velocity, force, or moment),
    only the orientation is undone; translation has no effect because the homogeneous
    last coordinate is 0.0.

    :param T_act: A (4,4) ndarray of floats representing an active homogeneous transform
        that operated within the current axis system.
    :return: A (4,4) ndarray of floats representing the active transform that exactly
        undoes T_act.
    """
    return _invert_T_rigid(T_act)


def convert_T_pas_to_T_act(
    T_pas: np.ndarray,
) -> np.ndarray:
    """Converts a passive transformation matrix to an active transformation matrix.

    A passive transform describes how to re-express the same physical quantity in
    different axes or relative to a different reference point. An active transform
    describes how to change the physical quantity itself within the same axes. This
    function converts between these interpretations by inverting the matrix.

    :param T_pas: A (4,4) ndarray of floats representing a passive transformation
        matrix.
    :return: A (4,4) ndarray of floats representing the converted active transformation
        matrix.
    """
    return np.linalg.inv(T_pas)


def convert_T_act_to_T_pas(
    T_act: np.ndarray,
) -> np.ndarray:
    """Converts an active transformation matrix to a passive transformation matrix.

    An active transform describes how to change the physical quantity itself within the
    same axes. A passive transform describes how to re-express the same physical
    quantity in different axes or relative to a different reference point. This function
    converts between these interpretations by inverting the matrix.

    :param T_act: A (4,4) ndarray of floats representing an active transformation
        matrix.
    :return: A (4,4) ndarray of floats representing the converted passive transformation
        matrix.
    """
    return np.linalg.inv(T_act)


def apply_T_to_vectors(
    T: np.ndarray,
    vectors_A: np.ndarray,
    is_position: bool,
) -> np.ndarray:
    """Applies a homogeneous transform to 3-element vector(s) and returns 3-element
    vector(s).

    For passive T, this function maps components from source axes and reference point to
    target axes and reference point. For active T, this function re-orients or
    translates the vector(s) within the same axes.

    This function handles both single vectors and arrays of vectors efficiently using
    einsum operations.

    **Warning:**

    When T is an improper transform (such as a reflection), this correctly reflects
    position, velocity, and force vectors but returns the negative of the correct result
    for pseudovectors (axial vectors such as moments, angular velocities, vorticities,
    or surface normals).

    :param T: A (4,4) ndarray of floats representing a homogeneous transform (active or
        passive).
    :param vectors_A: A (...,3) ndarray of floats representing the vector(s) to
        transform. Can be a single (3,) vector or a (...,3) array of vectors.
    :param is_position: True if the vector is a position (a point), so the transform's
        translation applies. False otherwise (for example a velocity, force, or moment),
        so only the rotation applies.
    :return: A ndarray of floats with same shape as ``vectors_A`` representing the
        transformed vector(s).
    """
    vectorsHomog_A = _generate_homogs(vectors_A, is_position)
    return np.asarray(
        np.einsum("ij,...j->...i", T, vectorsHomog_A)[..., :3], dtype=float
    )


def R_to_quat_wxyz(R: np.ndarray) -> np.ndarray:
    """Converts a rotation matrix to a unit quaternion.

    **Citation:**

    Equation adapted from: "Accurate Computation of Quaternions from Rotation Matrices"

    Authors: Soheil Sarabandi and Federico Thomas

    Date retrieved: 11/25/2025

    :param R: A (3,3) ndarray of floats representing a rotation matrix.
    :return: A (4,) ndarray of floats representing the unit quaternion.
    """
    r_11, r_12, r_13 = R[0]
    r_21, r_22, r_23 = R[1]
    r_31, r_32, r_33 = R[2]

    q_1: float
    check_1 = r_11 + r_22 + r_33
    if check_1 > _QUAT_BRANCH_THRESHOLD:
        q_1 = 0.5 * np.sqrt(1 + check_1)
    else:
        num_1 = (r_32 - r_23) ** 2 + (r_13 - r_31) ** 2 + (r_21 - r_12) ** 2
        den_1 = 3 - check_1
        q_1 = 0.5 * np.sqrt(num_1 / den_1)

    q_2_abs: float
    check_2 = r_11 - r_22 - r_33
    if check_2 > _QUAT_BRANCH_THRESHOLD:
        q_2_abs = 0.5 * np.sqrt(1 + check_2)
    else:
        num_2 = (r_32 - r_23) ** 2 + (r_12 + r_21) ** 2 + (r_31 + r_13) ** 2
        den_2 = 3 - check_2
        q_2_abs = 0.5 * np.sqrt(num_2 / den_2)
    q_2_sign = 1.0 if (r_32 - r_23) >= 0.0 else -1.0
    q_2 = q_2_sign * q_2_abs

    q_3_abs: float
    check_3 = -r_11 + r_22 - r_33
    if check_3 > _QUAT_BRANCH_THRESHOLD:
        q_3_abs = 0.5 * np.sqrt(1 + check_3)
    else:
        num_3 = (r_13 - r_31) ** 2 + (r_12 + r_21) ** 2 + (r_23 + r_32) ** 2
        den_3 = 3 - check_3
        q_3_abs = 0.5 * np.sqrt(num_3 / den_3)
    q_3_sign = 1.0 if (r_13 - r_31) >= 0.0 else -1.0
    q_3 = q_3_sign * q_3_abs

    q_4_abs: float
    check_4 = -r_11 - r_22 + r_33
    if check_4 > _QUAT_BRANCH_THRESHOLD:
        q_4_abs = 0.5 * np.sqrt(1 + check_4)
    else:
        num_4 = (r_21 - r_12) ** 2 + (r_31 + r_13) ** 2 + (r_32 + r_23) ** 2
        den_4 = 3 - check_4
        q_4_abs = 0.5 * np.sqrt(num_4 / den_4)
    q_4_sign = 1.0 if (r_21 - r_12) >= 0.0 else -1.0
    q_4 = q_4_sign * q_4_abs

    q = np.asarray((q_1, q_2, q_3, q_4), dtype=float)
    q_mag = float(np.linalg.norm(q))

    return q / q_mag


def R_to_angles_izyx(R: np.ndarray) -> np.ndarray:
    """Converts a rotation matrix to intrinsic z-y'-x" Euler angles in degrees.

    The returned ``[angleX, angleY, angleZ]`` always satisfies the matrix round-trip
    property: ``generate_rot_T(angles=[angleX, angleY, angleZ], passive=True,
    intrinsic=True, order="zyx")`` produces a (4,4) homogeneous transform whose
    rotation block equals R.

    The angle round-trip property (recovering the same angles that built R) only
    holds when ``|sin(angleY)| <= 1 - 1e-12``, equivalently when ``|angleY|`` is more
    than about 8e-5 degrees from the +/- 90 degree gimbal-lock pole. Inside the pole
    band, the standard atan2 decomposition becomes ill conditioned because
    cos(angleY) approaches machine epsilon, so the helper falls back to an
    alternative decomposition with ``angleX = 0`` and ``angleZ`` absorbing the
    indeterminate rotation. That fallback is one of infinitely many valid
    decompositions at the pole, so the recovered angles will differ from the
    originals there even though the matrix is reconstructed correctly.

    :param R: A (3,3) ndarray of floats representing a rotation matrix.
    :return: A (3,) ndarray of floats representing [angleX, angleY, angleZ] in
        degrees.
    """
    sin_angleY = float(np.clip(-R[0, 2], -1.0, 1.0))
    angleY = float(np.rad2deg(np.arcsin(sin_angleY)))
    if abs(sin_angleY) > _GIMBAL_LOCK_THRESHOLD:
        angleX = 0.0
        angleZ = float(np.rad2deg(np.arctan2(-R[1, 0], R[1, 1])))
    else:
        angleX = float(np.rad2deg(np.arctan2(R[1, 2], R[2, 2])))
        angleZ = float(np.rad2deg(np.arctan2(R[0, 1], R[0, 0])))
    return np.asarray([angleX, angleY, angleZ], dtype=float)


def alpha_and_beta_from_vInf_BP1(
    vInf_BP1__E: np.ndarray,
    vCg__E: float,
) -> tuple[float, float]:
    """Extracts the angle of attack and angle of sideslip from the freestream velocity
    in the first Airplane's body axes.

    Returns ``(nan, nan)`` when ``vCg__E`` is exactly zero. The freestream has no
    preferred direction at zero speed, so alpha and beta are physically undefined;
    NaN is the IEEE representation of an undefined real value. Callers that need to
    handle the zero-speed case must check for NaN in the result. The broader free
    flight simulation pipeline is already mathematically degenerate at zero speed
    (load coefficients are forces / qInf, and qInf = 0.5 * rho * vCg__E**2 = 0), so
    this function does not attempt to paper over the degeneracy with a substitute
    value.

    :param vInf_BP1__E: A (3,) ndarray of floats representing the freestream velocity
        (in the first Airplane's body axes, observed from the Earth frame) in meters per
        second.
    :param vCg__E: A float representing the speed of the first Airplane's CG (observed
        from the Earth frame) in meters per second.
    :return: A tuple (alpha, beta) where alpha is the angle of attack in degrees and
        beta is the angle of sideslip in degrees. Both are NaN if vCg__E is zero.
    """
    if vCg__E == 0.0:
        return float("nan"), float("nan")

    vInfX_BP1__E, vInfY_BP1__E, vInfZ_BP1__E = vInf_BP1__E

    # Invert the wind axes construction defined in docs/AXES_POINTS_AND_FRAMES.md and
    # implemented by OperatingPoint. In that convention the CG velocity in body axes (the
    # negated freestream, vCg_BP1__E = -vInf_BP1__E) has components
    #   x: vCg__E * cos(alpha) * cos(beta)
    #   y: vCg__E * cos(alpha) * sin(beta)
    #   z: vCg__E * sin(alpha)
    # so alpha follows from the body z component (arcsin) and beta from the body x and y
    # components (arctan2). Extracting them in this order, rather than the more common
    # textbook order that swaps which angle uses arcsin, keeps this function the exact
    # inverse of OperatingPoint. Deriving alpha and beta here and storing them back on an
    # OperatingPoint then reproduces the original freestream.
    sin_alpha = float(np.clip(-vInfZ_BP1__E / vCg__E, -1.0, 1.0))
    alpha = float(np.rad2deg(np.arcsin(sin_alpha)))
    beta = float(np.rad2deg(np.arctan2(-vInfY_BP1__E, -vInfX_BP1__E)))
    return alpha, beta
