"""This module contains functions to create fixtures for parameter validation tests."""

import numpy as np


# Valid string fixtures.
def make_valid_str_fixture():
    """Makes a fixture that is a valid string for str_return_str.

    :return: A valid string.
    """
    return "test_string"


def make_empty_str_fixture():
    """Makes a fixture that is an empty string for str_return_str.

    :return: An empty string (still valid as a str type).
    """
    return ""


# Valid bool fixtures.
def make_valid_bool_fixture():
    """Makes a fixture that is a valid bool for boolLike_return_bool.

    :return: A valid bool.
    """
    return True


def make_valid_numpy_bool_fixture():
    """Makes a fixture that is a valid numpy bool for boolLike_return_bool.

    :return: A valid numpy bool.
    """
    return np.bool_(False)


# Valid int in range fixtures.
def make_valid_int_fixture():
    """Makes a fixture that is a valid int for int_in_range_return_int.

    :return: A valid int.
    """
    return 5


def make_valid_int_at_lower_bound_inclusive_fixture():
    """Makes a fixture that is a valid int at the lower bound (inclusive).

    :return: A valid int at lower bound.
    """
    return 0


def make_valid_int_at_upper_bound_inclusive_fixture():
    """Makes a fixture that is a valid int at the upper bound (inclusive).

    :return: A valid int at upper bound.
    """
    return 10


# Valid number in range fixtures.
def make_valid_float_fixture():
    """Makes a fixture that is a valid float for number_in_range_return_float.

    :return: A valid float.
    """
    return 3.14


def make_valid_int_as_number_fixture():
    """Makes a fixture that is a valid int for number_in_range_return_float.

    :return: A valid int (acceptable as a number).
    """
    return 7


# Valid array-like fixtures.
def make_valid_1d_array_fixture():
    """Makes a fixture that is a valid 1D array for array-like validation functions.

    :return: A (5,) ndarray of floats.
    """
    return np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)


def make_valid_list_fixture():
    """Makes a fixture that is a valid list for array-like validation functions.

    :return: A list of floats.
    """
    return [1.0, 2.0, 3.0]


def make_valid_tuple_fixture():
    """Makes a fixture that is a valid tuple for array-like validation functions.

    :return: A tuple of floats.
    """
    return (1.0, 2.0, 3.0)


def make_valid_nested_list_fixture():
    """Makes a fixture that is a valid nested list for 2D array-like validation.

    :return: A nested list of floats.
    """
    return [[1.0, 2.0], [3.0, 4.0]]


# Valid 3D vector fixtures.
def make_valid_3d_vector_fixture():
    """Makes a fixture that is a valid 3D vector for threeD_number_vectorLike_return_float.

    :return: A (3,) ndarray of floats.
    """
    return np.array([1.0, 2.0, 3.0], dtype=float)


def make_valid_3d_list_fixture():
    """Makes a fixture that is a valid 3D list for threeD_number_vectorLike_return_float.

    :return: A list with 3 elements.
    """
    return [1.0, 2.0, 3.0]


def make_valid_3d_tuple_fixture():
    """Makes a fixture that is a valid 3D tuple for threeD_number_vectorLike_return_float.

    :return: A tuple with 3 elements.
    """
    return (1.0, 2.0, 3.0)


def make_valid_unit_vector_fixture():
    """Makes a fixture that is a valid unit vector (already normalized).

    :return: A (3,) ndarray of floats with magnitude 1.0.
    """
    return np.array([1.0, 0.0, 0.0], dtype=float)


def make_valid_non_unit_vector_fixture():
    """Makes a fixture that is a valid vector that needs normalization.

    :return: A (3,) ndarray of floats with magnitude != 1.0.
    """
    return np.array([3.0, 4.0, 0.0], dtype=float)


# Valid 2D vector fixtures.
def make_valid_2d_vector_fixture():
    """Makes a fixture that is a valid 2D vector for arrayLike_of_twoD_number_vectorLikes_return_float.

    :return: A (2,) ndarray of floats.
    """
    return np.array([1.0, 2.0], dtype=float)


def make_valid_2d_vectors_array_fixture():
    """Makes a fixture that is a valid array of 2D vectors.

    :return: A (3, 2) ndarray of floats.
    """
    return np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=float)


# Valid 3D vectors array fixtures.
def make_valid_3d_vectors_array_fixture():
    """Makes a fixture that is a valid array of 3D vectors.

    :return: A (4, 3) ndarray of floats.
    """
    return np.array(
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]],
        dtype=float,
    )


# Valid spacing vector fixtures.
def make_valid_spacing_vector_all_strings_fixture():
    """Makes a fixture that is a valid spacing vector with all string elements.

    :return: A tuple with 3 string elements.
    """
    return ("sine", "uniform", "sine")


def make_valid_spacing_vector_all_callables_fixture():
    """Makes a fixture that is a valid spacing vector with all callable elements.

    :return: A tuple with 3 callable elements.
    """

    def custom_spacing(x):
        return np.sin(x)

    return (custom_spacing, custom_spacing, custom_spacing)


def make_valid_spacing_vector_mixed_fixture():
    """Makes a fixture that is a valid spacing vector with mixed string and callable elements.

    :return: A tuple with mixed string and callable elements.
    """

    def custom_spacing(x):
        return np.sin(x)

    return ("sine", custom_spacing, "uniform")


# Valid 4x4 matrix fixtures.
def make_valid_4x4_matrix_fixture():
    """Makes a fixture that is a valid 4x4 matrix for fourByFour_number_arrayLike_return_float.

    :return: A (4, 4) ndarray of floats.
    """
    return np.eye(4, dtype=float)


def make_valid_4x4_list_fixture():
    """Makes a fixture that is a valid 4x4 nested list for fourByFour_number_arrayLike_return_float.

    :return: A nested list with shape (4, 4).
    """
    return [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]


# Valid N-D vector fixtures.
def make_valid_nd_vector_fixture():
    """Makes a fixture that is a valid N-D vector for nD_number_vectorLike_return_float.

    :return: A (7,) ndarray of floats.
    """
    return np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], dtype=float)


# Valid non-empty list fixtures.
def make_valid_non_empty_list_fixture():
    """Makes a fixture that is a valid non-empty list for non_empty_list_return_list.

    :return: A non-empty list.
    """
    return [1, 2, 3]


def make_valid_single_element_list_fixture():
    """Makes a fixture that is a valid single-element list.

    :return: A list with one element.
    """
    return ["single"]


# Valid rotation order fixtures.
def make_valid_rotation_order_xyz_fixture():
    """Makes a fixture that is a valid rotation order string "xyz".

    :return: The string "xyz".
    """
    return "xyz"


def make_valid_rotation_order_zyx_fixture():
    """Makes a fixture that is a valid rotation order string "zyx".

    :return: The string "zyx".
    """
    return "zyx"


# Invalid fixtures (for testing error conditions).
def make_invalid_non_str_fixture():
    """Makes a fixture that is not a string for testing str_return_str errors.

    :return: An int (invalid type for str validation).
    """
    return 123


def make_invalid_non_bool_fixture():
    """Makes a fixture that is not a bool for testing boolLike_return_bool errors.

    :return: A string (invalid type for bool validation).
    """
    return "true"


def make_invalid_float_for_int_fixture():
    """Makes a fixture that is a float for testing int_in_range_return_int errors.

    :return: A float (invalid type for int validation).
    """
    return 5.5


def make_invalid_nan_fixture():
    """Makes a fixture that is NaN for testing finite value validation errors.

    :return: np.nan.
    """
    return np.nan


def make_invalid_inf_fixture():
    """Makes a fixture that is infinity for testing finite value validation errors.

    :return: np.inf.
    """
    return np.inf


def make_invalid_neg_inf_fixture():
    """Makes a fixture that is negative infinity for testing finite value validation errors.

    :return: -np.inf.
    """
    return -np.inf


def make_invalid_2d_vector_as_3d_fixture():
    """Makes a fixture that is a 2D vector for testing 3D vector validation errors.

    :return: A (2,) ndarray (invalid shape for 3D validation).
    """
    return np.array([1.0, 2.0], dtype=float)


def make_invalid_4d_vector_as_3d_fixture():
    """Makes a fixture that is a 4D vector for testing 3D vector validation errors.

    :return: A (4,) ndarray (invalid shape for 3D validation).
    """
    return np.array([1.0, 2.0, 3.0, 4.0], dtype=float)


def make_invalid_zero_vector_fixture():
    """Makes a fixture that is a zero vector for testing unit vector validation errors.

    :return: A (3,) ndarray with all zeros.
    """
    return np.array([0.0, 0.0, 0.0], dtype=float)


def make_invalid_3x3_matrix_fixture():
    """Makes a fixture that is a 3x3 matrix for testing 4x4 matrix validation errors.

    :return: A (3, 3) ndarray (invalid shape for 4x4 validation).
    """
    return np.eye(3, dtype=float)


def make_invalid_empty_list_fixture():
    """Makes a fixture that is an empty list for testing non-empty list validation errors.

    :return: An empty list.
    """
    return []


def make_invalid_rotation_order_fixture():
    """Makes a fixture that is an invalid rotation order string.

    :return: An invalid rotation order string.
    """
    return "abc"


def make_invalid_spacing_string_fixture():
    """Makes a fixture that is an invalid spacing string.

    :return: A tuple with an invalid spacing string.
    """
    return ("sine", "invalid", "uniform")


def make_invalid_spacing_type_fixture():
    """Makes a fixture that is an invalid spacing type (not string or callable).

    :return: A tuple with an invalid type (int).
    """
    return ("sine", 123, "uniform")


def make_scalar_fixture():
    """Makes a fixture that is a scalar for testing array-like validation errors.

    :return: A scalar value.
    """
    return 5.0


def make_invalid_array_with_nan_fixture():
    """Makes a fixture that is an array containing NaN.

    :return: A (3,) ndarray containing NaN.
    """
    return np.array([1.0, np.nan, 3.0], dtype=float)


def make_invalid_array_with_inf_fixture():
    """Makes a fixture that is an array containing infinity.

    :return: A (3,) ndarray containing infinity.
    """
    return np.array([1.0, np.inf, 3.0], dtype=float)
