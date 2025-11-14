"""This module contains functions used for geometric transformations."""

import numpy as np

from . import _parameter_validation


def _generate_homogs(vectors_A: np.ndarray, has_point: bool) -> np.ndarray:
    """This function converts 3D vector(s) to homogeneous coordinates for use with
    (4,4) transformation matrices.

    Homogeneous coordinates extend 3D vectors to 4D by adding a fourth component. For
    vectors relative to a reference point (has_point=True), such as position vectors,
    the fourth component is 1.0. For free vectors (has_point=False), such as velocity
    or force vectors, the fourth component is 0.0. This allows (4,4) transformation
    matrices to handle both translations and rotations in a unified framework.

    This function handles both single vectors and arrays of vectors efficiently.

    This is a private function, so it assumes all parameters have been pre-validated.

    :param vectors_A: ndarray of shape (..., 3)
        Array of 3D vectors or a single 3D vector
    :param has_point: bool
        True if vectors are relative to a reference point, False for free vectors
    :return: ndarray of shape (..., 4)
        Array of homogeneous vectors with same leading dimensions as input
    """
    # Create homogeneous array with one extra dimension
    vectorsHomog_A = np.zeros(vectors_A.shape[:-1] + (4,), dtype=float)

    # Copy the 3D components
    vectorsHomog_A[..., :3] = vectors_A

    # Set the homogeneous coordinate
    if has_point:
        vectorsHomog_A[..., -1] = 1.0

    return vectorsHomog_A


def generate_rot_T(angles, passive, intrinsic, order) -> np.ndarray:
    """This function generates a transformation matrix representing a passive or
    active rotation using an angle vector, and parameters specifying if the matrix
    should be passive or active, intrinsic or extrinsic, and the order by which to
    apply the rotations about each axis.

    Notes: Angles are in degrees with signs defined using the right-hand rule.

    Passive Use-Case: Let `r_A` be a free vector in "A" axes, but we want to find
    `r_B`, which is the same vector, but expressed in "B" axes. The orientation of
    "B" axes relative to "A" axes is defined by the angle vector `angles` (with
    rotations in order and type defined by the variables `order` and `intrinsic`). Then:

    ```
    T_pas_A_to_B = generate_rot_T(angles, True, intrinsic, order)
    r_B =  apply_T_to_vectors(T_pas_A_to_B, r_A, has_point=False)
    ```

    Active Use-Case: Let `r_A` be a free vector in "A" axes, but we want to find
    `rPrime_A`, which is `r_A` rotated by the specified sequence: about the fixed "A"
    axes if `intrinsic=False`, or about the current, newly-rotated axes if
    `intrinsic=True`, with angles given by `angles` and the sequence defined by
    `order`. Then:

    ```
    rot_T_act = generate_rot_T(angles, False, intrinsic, order)
    rPrime_A = apply_T_to_vectors(rot_T_act, r_A, has_point=False)
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

    :param passive: boolLike

        If True, returns a matrix that changes coordinates from "A" to "B" axes (`r_B
        = R @ r_A`). If False, returns a matrix that rotates vectors in "A" axes (
        `rPrime_A = R @ r_A`). It can be a boolean or a NumPy boolean and will be
        converted internally to a boolean.

    :param intrinsic: boolLike

        If True, each subsequent rotation is applied to the current, newly-rotated
        axes. If False, all rotations are performed about the original, non-rotated
        "A" axes. It can be a boolean or a NumPy boolean and will be converted
        internally to a boolean.

    :param order: string of length 3

        This is the string of three characters that represents the rotation order.
        Each character can be 'x', 'y', or 'z'. Only Tait-Bryan angles are accepted
        so all accepted characters must be distinct.

    :return: (4, 4) ndarray of floats

        This is the transformation matrix.
    """
    angles = _parameter_validation.threeD_number_vectorLike_return_float(
        angles, "angles"
    )
    passive = _parameter_validation.boolLike_return_bool(passive, "passive")
    intrinsic = _parameter_validation.boolLike_return_bool(intrinsic, "intrinsic")
    order = _parameter_validation.rotation_order_return_string(order, "order")

    angleX_rad, angleY_rad, angleZ_rad = np.radians(angles)

    x_R_act = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(angleX_rad), -np.sin(angleX_rad)],
            [0.0, np.sin(angleX_rad), np.cos(angleX_rad)],
        ]
    )
    y_R_act = np.array(
        [
            [np.cos(angleY_rad), 0.0, np.sin(angleY_rad)],
            [0.0, 1.0, 0.0],
            [-np.sin(angleY_rad), 0.0, np.cos(angleY_rad)],
        ]
    )
    z_R_act = np.array(
        [
            [np.cos(angleZ_rad), -np.sin(angleZ_rad), 0.0],
            [np.sin(angleZ_rad), np.cos(angleZ_rad), 0.0],
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


def generate_trans_T(translations, passive) -> np.ndarray:
    """This function generates either a passive or active translational
    transformation matrix using a vector and a parameter specifying if the
    transformation should be passive or active.

    Passive Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `c_A_b`,
    which describes the location of "c", relative to the "b" point. The position of
    "b" is defined by `translations` (in "A" axes, relative to the point a). Then:

    ```
    T_pas_A_a_to_A_b = generate_trans_T(translations, True)
    c_A_b = apply_T_to_vectors(T_pas_A_a_to_A_b, c_A_a, has_point=True)
    ```

    Active Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `cPrime_A_a`,
    which is the position of "cPrime", which is point "c" offset by `translations` (
    in "A" axes). Then:

    ```
    translate_T_act = generate_trans_T(translations, False)
    cPrime_A_a = apply_T_to_vectors(translate_T_act, c_A_a, has_point=True)
    ```

    :param translations: array-like of 3 numbers

        Can be a tuple, list, or numpy array of numbers (int or float). Values are
        converted to floats internally. For `passive=True`, this is the position of
        the "b" point (in "A" axes, relative to the "a" point). For `passive=False`,
        this is the position (in "A" axes) of the offset point "cPrime" relative to
        the original "c" point.

    :param passive: boolLike

        If True, returns a matrix that changes reference point of a vector in
        homogeneous coordinates (`rHomog_A_b = T_trans @ rHomog_A_a`). If False,
        returns the position vector of point offset from an original position,
        still relative to the same point (`cPrimeHomog_A_a = T_trans @ cHomog_A_a`).
        It can be a boolean or a NumPy boolean and will be converted internally to a
        boolean.

    :return: (4, 4) ndarray of floats

        This is the transformation matrix.
    """
    p = _parameter_validation.threeD_number_vectorLike_return_float(
        translations, "translations"
    )
    passive = _parameter_validation.boolLike_return_bool(passive, "passive")
    T_trans = np.eye(4, dtype=float)

    T_trans[:3, 3] = -p if passive else p
    return T_trans


def generate_reflect_T(plane_point_A_a, plane_normal_A, passive) -> np.ndarray:
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
    T_pas_A_a_to_B_b = generate_reflect_T(plane_point_A_a, plane_normal_A, True)
    c_B_b = apply_T_to_vectors(T_pas_A_a_to_B_b, c_A_a, has_point=True)
    ```

    Active Use-Case: Let `c_A_a` be a vector which describes the location of point
    "c" (in "A" axes, relative to the "a" point). We want to find `cPrime_A_a`,
    which is the position of "cPrime", point "c" reflected across the plane defined
    by plane_point_A_a and plane_normal_A. Then:

    ```
    reflect_T_act = generate_reflect_T(plane_point_A_a, plane_normal_A, False)
    c_A_a = apply_T_to_vectors(reflect_T_act, c_A_a, has_point=True)
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

    :param passive: boolLike

        If True, returns a matrix that changes reference point and axes of a vector
        in homogeneous coordinates to a reference point and axes reflected about the
        specified plane (`cHomog_B_b = T_reflect @ cHomog_A_a`). If False, it returns
        a matrix that reflects a vector (in its original axes) about a specified
        plane (`cPrimeHomog_A_a = T_reflect @ cHomog_A_a`). It can be a boolean or a
        NumPy boolean and will be converted internally to a boolean.

    :return: (4, 4) ndarray of floats

        This is the transformation matrix.
    """
    p = _parameter_validation.threeD_number_vectorLike_return_float(
        plane_point_A_a, "plane_point_A_a"
    )
    n_hat = _parameter_validation.threeD_number_vectorLike_return_float_unit_vector(
        plane_normal_A, "plane_normal_A"
    )
    # noinspection PyUnusedLocal
    passive = _parameter_validation.boolLike_return_bool(passive, "passive")

    T_reflect = np.eye(4, dtype=float)

    S = np.eye(3, dtype=float) - 2 * np.outer(n_hat, n_hat)
    d = 2 * (np.dot(p, n_hat)) * n_hat

    T_reflect[:3, :3] = S
    T_reflect[:3, 3] = d

    return T_reflect


def _left_compose_T(valid_T_chain) -> np.ndarray:
    """Left-compose a list of homogeneous transformations

    Given: `valid_T_chain=[T_1, T_2, ..., T_n]`
    Internally returns: `T_n @ ... @ T_2 @ T_1`

    This is a private function, so it assumes all parameters have been pre-validated.

    :param valid_T_chain: (4,4) ndarray of floats
    :return: (4,4) ndarray of floats
    """
    if len(valid_T_chain) == 1:
        return valid_T_chain[0]
    composed_T = np.eye(4, dtype=float)
    for valid_T in valid_T_chain:
        composed_T = valid_T @ composed_T
    return composed_T


def compose_T_pas(*T_pas_chain) -> np.ndarray:
    """Compose a chain of passive homogeneous transforms.

    Pass transforms in the "path order" from the original axes/point to
    the final axes/point:

    T_pas_A_a_to_C_c = compose_T_pas(T_pas_A_a_to_B_b, T_pas_B_b_to_C_c)

    Internally returns: T_pas_B_b_to_C_c @ T_pas_A_a_to_B_b

    :param T_pas_chain: sequence of (4,4) array-like objects filled with numbers
        The passive homogeneous transforms along the path.

    :return: (4,4) ndarray of floats
        The composed passive homogeneous transform.
    """
    if not T_pas_chain:
        raise ValueError("At least one transform must be provided.")

    valid_T_pas_chain = []
    for T_pas_id, T_pas in enumerate(T_pas_chain):
        valid_T_pas_chain.append(
            _parameter_validation.fourByFour_number_arrayLike_return_float(
                T_pas, f"T_pas_chain[{T_pas_id}]"
            )
        )
    return _left_compose_T(valid_T_pas_chain)


def compose_T_act(*T_act_chain) -> np.ndarray:
    """Compose a chain of active homogeneous transforms.

    Pass transforms in the *chronological order* you want to apply them to the vector:

    composed_T_act = compose_T_act(reflect_T_act, rot_T_act, trans_T_act)

    Internally returns: translate_T_act @ rot_T_act @ reflect_T_act

    Note on body-fixed vs world-fixed interpretations:
        This function left-composes the supplied active transforms: given
        compose_T_act( T1, T2, ..., Tn) it returns Tn @ ... @ T2 @ T1. Pass
        transforms in the chronological order you want them applied.

        Active translations created with generate_trans_T(..., passive=False)
        interpret the components in the same axes the vector is expressed in (e.g.,
        G axes). Therefore: ``` T_act = compose_T_act(rot_T_act, trans_T_act) ```
        applies a rotation first and then a *world-fixed* translation.

        If you instead want a *body-fixed* translation (e.g., “+10 along x′ after the
        rotation”), either:
            (a) pre-rotate the components before building the translation:
            ```
            R = rot_T_act[:3, :3]
            tPrime_A = R @ t_A
            transBodyFixed_T_act = generate_trans_T(tPrime_A, False)
            T_act = compose_T_act(rot_T_act, transBodyFixed_T_act)
            ```
            or
            (b) pass the translation before the rotation:
            ```
            T_act = compose_T_act(trans_T_act, rot_T_act)
            ```

        Both yield the same homogeneous matrix: [[R, R @ t_body], [0, 0, 0, 1]].

    :param T_act_chain: sequence of (4,4) array-like objects filled with numbers

        The active homogeneous transforms in chronological order.

    :return: (4,4) ndarray of floats

        The composed active homogeneous transform.
    """
    if not T_act_chain:
        raise ValueError("At least one transform must be provided.")

    valid_T_act_chain = []
    for T_act_id, T_act in enumerate(T_act_chain):
        valid_T_act_chain.append(
            _parameter_validation.fourByFour_number_arrayLike_return_float(
                T_act, f"T_act_chain[{T_act_id}]"
            )
        )
    return _left_compose_T(valid_T_act_chain)


def _invert_T_rigid(valid_T: np.ndarray) -> np.ndarray:
    """Invert a rigid homogeneous transform.

    Given:
    ```
    R = valid_T[:3, :3]
    t = valid_T[:3, 3]
    ```
    Internally returns:
    [ R.T   -R.T @ t ]
    [  0       1     ]

    This is a private function, so it assumes all parameters have been pre-validated.

    :param valid_T: (4,4) ndarray of floats
    :return: (4,4) ndarray of floats
    """
    valid_R = valid_T[:3, :3]
    valid_t = valid_T[:3, 3]

    valid_transpose_R = valid_R.T

    valid_T_inv = np.eye(4, dtype=float)

    valid_T_inv[:3, :3] = valid_transpose_R
    valid_T_inv[:3, 3] = -valid_transpose_R @ valid_t

    return valid_T_inv


def invert_T_pas(T_pas) -> np.ndarray:
    """Invert a passive homogeneous transform.

    A passive transform maps components of the same physical quantity between an
    initial axis system and reference point and a target axis system and reference
    point. For example, if T_pas_A_a_to_B_b maps components from A axes relative to
    point a to B axes relative to point b, then:

    ```
    T_pas_B_b_to_A_a = invert_T_pas(T_pas_A_a_to_B_b)
    rHomog_A_a = T_pas_B_b_to_A_a @ rHomog_B_b
    ```

    Notes:
    - For vectors relative to a reference point (has_point=True), such as position
    vectors, the translation component matters.
    - For free vectors (has_point=False), translation has no effect because the
    homogeneous last coordinate is 0.0.

    Example 1:
    # Given T_pas_A_a_to_B_b, get the inverse mapper:
    T_pas_B_b_to_A_a = invert_T_pas(T_pas_A_a_to_B_b)

    # Apply to the position vector r (in B axes, relative to point b):
    r_A_a = apply_T_to_vectors(T_pas_B_b_to_A_a, r_B_b, has_point=True)

    Example 2:
    # Apply to the free vector v (in B axes) (translation will be ignored):
    v_A = apply_T_to_vectors(T_pas_B_b_to_A_a, v_B)

    :param T_pas: (4,4) array-like of numbers

        Passive homogeneous transform mapping from source axes+reference point to
        target axes + reference point.  Can be any (4,4) array-like object of size of
        numbers (int or float). Values are converted to floats internally.

    :return: (4,4) ndarray of floats

        The passive transform that maps back from the target axes + reference point
        to the original axes + reference point.
    """
    valid_T_pas = _parameter_validation.fourByFour_number_arrayLike_return_float(
        T_pas, "T_pas"
    )
    return _invert_T_rigid(valid_T_pas)


def invert_T_act(T_act) -> np.ndarray:
    """Invert an active homogeneous transform.

    An active transform re-orients and optionally translates a quantity within the
    same axis system. For example, if T_act transforms the free vector q (in A axes)
    to the free vector qPrime (in A axes), then:

    ```
    q_A = apply_T_to_vectors(invert_T_act(T_act), qPrime_A, has_point=False)
    ```

    Notes:
    - For vectors relative to a reference point (has_point=True), such as position
    vectors, both orientation and translation are undone.
    - For free vectors (has_point=False), only the orientation is undone; translation
    has no effect because the homogeneous last coordinate is 0.

    Example 1:
    # Undo an active sequence applied to a position:
    rPrime_A = apply_T_to_vectors(T_act, r_A, has_point=True)
    undo_T_act = invert_T_act(T_act)

    recovered_r_A = apply_T_to_vectors(undo_T_act, rPrime_A, has_point=True)

    Example 2:
    # For a free vector, translation has no effect:
    fPrime_A = apply_T_to_vectors(T_act, f_A, has_point=False)
    undo_T_act = invert_T_act(T_act)

    recovered_f_A = apply_T_to_vectors(undo_T_act, fPrime_A, has_point=False)

    :param T_act: (4,4) array-like of numbers

        Active homogeneous transform that operated within the current axis system.
        Can be any (4,4) array-like object of size of numbers (int or float). Values
        are converted to floats internally.

    :return: (4,4) ndarray of floats

        The active transform that exactly undoes T_act.
    """
    valid_T_act = _parameter_validation.fourByFour_number_arrayLike_return_float(
        T_act, "T_act"
    )
    return _invert_T_rigid(valid_T_act)


def convert_T_pas_to_T_act(T_pas) -> np.ndarray:
    """Convert a passive transformation matrix to an active transformation matrix.

    :param T_pas: (4,4) array-like of numbers

        The passive transformation matrix.

    :return: (4,4) ndarray of floats

        The converted active transformation matrix.
    """
    valid_T_pas = _parameter_validation.fourByFour_number_arrayLike_return_float(
        T_pas, "T_pas"
    )
    return np.linalg.inv(valid_T_pas)


def convert_T_act_to_T_pas(T_act) -> np.ndarray:
    """Convert an active transformation matrix to a passive transformation matrix.

    :param T_act: (4,4) array-like of numbers

        The active transformation matrix.

    :return: (4,4) ndarray of floats

        The converted passive transformation matrix.
    """
    valid_T_act = _parameter_validation.fourByFour_number_arrayLike_return_float(
        T_act, "T_act"
    )
    return np.linalg.inv(valid_T_act)


def apply_T_to_vectors(T, vectors_A, has_point) -> np.ndarray:
    """Apply a homogeneous transform to 3-element vector(s) and return 3-element
    vector(s).

    For passive T: interpret as mapping components from source axes/point to target
    axes/point.

    For active T: interpret as re-orienting/translating the vector(s) within the same
    axes.

    This function handles both single vectors and arrays of vectors efficiently using
    einsum operations.

    :param T: (4,4) array-like of numbers

        Homogeneous transform (active or passive). Can be any (4,4) array-like object
        of size of numbers (int or float). Values are converted to floats internally.

    :param vectors_A: array-like of shape (..., 3)

        The 3-vector(s) to transform. Can be a single vector (3,) or array of vectors
        (..., 3). Can be a tuple, list, or numpy array of numbers (int or float).
        Values are converted to floats internally.

    :param has_point: boolLike

        True for vectors relative to a reference point such as position vectors.
        False for free vectors. It can be a boolean or a NumPy boolean and will be
        converted internally to a boolean.

    :return: ndarray of floats with same shape as input vectors_A

        The transformed vector(s).
    """
    T = _parameter_validation.fourByFour_number_arrayLike_return_float(T, "T")
    vectors_A = (
        _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
            vectors_A, "vectors_A"
        )
    )
    has_point = _parameter_validation.boolLike_return_bool(has_point, "has_point")

    vectorsHomog_A = _generate_homogs(vectors_A, has_point)

    return np.einsum("ij,...j->...i", T, vectorsHomog_A)[..., :3]
