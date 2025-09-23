# TODO: Consider making this module private (renaming it with a _ prefix).
"""This module contains common parameter validation functions."""

import numpy as np


def string_return_string(string, name):
    """Validates that a value is a string and returns it as a string. name must also
    be a string."""
    if not isinstance(name, str):
        raise TypeError("name must be a string.")

    if not isinstance(string, str):
        raise TypeError(f"{name} must be a string.")
    return string


def boolLike_return_bool(value, name):
    """Validates that a value is a boolean or NumPy boolean and returns it as a
    boolean. name must be a string."""
    name = string_return_string(name, "name")

    if not isinstance(value, (bool, np.bool)):
        raise TypeError(
            f"{name} is a {type(value)} equal to {value} but it must be a boolean."
        )

    return bool(value)


def int_return_int(value, name):
    """Validates that a value is an int and returns it as an int. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test.
    """
    name = string_return_string(name, "name")

    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")

    return value


# TODO: Consider getting rid of this function as its functionality can be replicated
#  with validate_scalar_int_in_range.
def non_negative_int_return_int(value, name):
    """Validates that a value is an int with a value greater than or equal to zero
    and returns it as an int. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")

    if value < 0:
        raise ValueError(f"{name} must be greater than or equal to zero.")

    return int(value)


# TODO: Consider getting rid of this function as its functionality can be replicated
#  with validate_scalar_int_in_range.
def positive_int_return_int(value, name):
    """Validates that a value is an int with a value greater than zero and returns it
    as an int. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    if not isinstance(value, int):
        raise TypeError(f"{name} must be an int.")

    if value <= 0:
        raise ValueError(f"{name} must be positive.")

    return value


def int_in_range_return_int(
    value, name, min_val, min_inclusive, max_val, max_inclusive
):
    """Validates that a value is an int and in a custom range and returns it as an
    int. If min_val or max_val is None, then value's magnitude isn't checked relative
    to that parameter. If min_val or max_val is None, the corresponding '*_inclusive'
    parameter must also be None. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")
    if min_val is not None:
        min_val = number_return_float(min_val, "min_val")
        min_inclusive = boolLike_return_bool(min_inclusive, "min_inclusive")
    else:
        if min_inclusive is not None:
            raise ValueError("min_inclusive must be None if min_val is None")
    if max_val is not None:
        max_val = number_return_float(max_val, "max_val")
        max_inclusive = boolLike_return_bool(max_inclusive, "max_inclusive")
    else:
        if max_inclusive is not None:
            raise ValueError("max_inclusive must be None if max_val is None")

    if min_val is not None and max_val is not None:
        if min_val >= max_val:
            raise ValueError("min_val must be less than max_val")

    if not isinstance(value, int):
        raise TypeError(f"{name} must be an int.")

    if min_val is not None:
        if min_inclusive:
            if not value >= min_val:
                raise ValueError(f"{name} must be greater than or equal to {min_val}.")
        else:
            if not value > min_val:
                raise ValueError(f"{name} must be greater than {min_val}.")

    if max_val is not None:
        if max_inclusive:
            if not value <= max_val:
                raise ValueError(f"{name} must be less than or equal to {max_val}.")
        else:
            if not value < max_val:
                raise ValueError(f"{name} must be less than {min_val}.")

    return int(value)


def number_return_float(value, name):
    """Validates a value is a number and returns it as a float. name must be a
    string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    return float(value)


# TODO: Consider getting rid of this function as its functionality can be replicated
#  with number_in_range_return_float.
def non_negative_number_return_float(value, name):
    """Validates a value is a number and is greater than or equal to zero and returns
    it as a float. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if value < 0:
        raise ValueError(f"{name} must be greater than or equal to zero.")

    return float(value)


# TODO: Consider getting rid of this function as its functionality can be replicated
#  with number_in_range_return_float.
def positive_number_return_float(value, name):
    """Validates a value is a number and is greater than zero and returns it as a
    float. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if value <= 0:
        raise ValueError(f"{name} must be positive.")

    return float(value)


def number_in_range_return_float(
    value, name, min_val, min_inclusive, max_val, max_inclusive
):
    """Validates a value is a number and is in a custom range and returns it as a
    float. If min_val or max_val is None, then value's magnitude isn't checked
    relative to that parameter. If min_val or max_val is None, the corresponding
    '*_inclusive' parameter must also be None. If not None, these parameters must be
    booleans. If neither min_val nor max_val are None, then min_val must be less than
    max_val. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")
    if min_val is not None:
        min_val = number_return_float(min_val, "min_val")
        min_inclusive = boolLike_return_bool(min_inclusive, "min_inclusive")
    else:
        if min_inclusive is not None:
            raise ValueError("min_inclusive must be None if min_val is None")
    if max_val is not None:
        max_val = number_return_float(max_val, "max_val")
        max_inclusive = boolLike_return_bool(max_inclusive, "max_inclusive")
    else:
        if max_inclusive is not None:
            raise ValueError("max_inclusive must be None if max_val is None")

    if min_val is not None and max_val is not None:
        if min_val >= max_val:
            raise ValueError("min_val must be less than max_val")

    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if min_val is not None:
        if min_inclusive:
            if not value >= min_val:
                raise ValueError(f"{name} must be greater than or equal to {min_val}.")
        else:
            if not value > min_val:
                raise ValueError(f"{name} must be greater than {min_val}.")

    if max_val is not None:
        if max_inclusive:
            if not value <= max_val:
                raise ValueError(f"{name} must be less than or equal to {max_val}.")
        else:
            if not value < max_val:
                raise ValueError(f"{name} must be less than {min_val}.")

    return float(value)


def arrayLike_of_twoD_number_vectorLikes_return_float(vectors, name):
    """Validates a value is an array-like object of 2D number vector-like objects (
    array-like objects with shape (2,)). It then returns it as a (...,2) numpy array
    of floats. name must be a string.

    Accepts both single vectors of shape (2,) and arrays of vectors with shape
    (..., 2) where the last dimension must be 2.

    Note: np.nan, np.inf, and -np.inf won't pass this test.
    """
    name = string_return_string(name, "name")

    try:
        vectors = np.asarray(vectors, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vectors.ndim == 0:
        raise ValueError(f"{name} cannot be a scalar.")

    if vectors.shape[-1] != 2:
        raise ValueError(f"{name} must have 2 elements in the last dimension.")

    if not np.isfinite(vectors).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return vectors


def threeD_number_vectorLike_return_float(vector, name):
    """Validates a value is a 3D vector-like object (array-like object with shape (3,
    )). It then returns it as a (3,) numpy array of floats. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        vector = np.asarray(vector, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    if not np.isfinite(vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return vector


def arrayLike_of_threeD_number_vectorLikes_return_float(vectors, name):
    """Validates a value is an array-like object of 3D number vector-like objects (
    array-like objects with shape (3,)). It then returns it as a (...,3) numpy array
    of floats. name must be a string.

    Accepts both single vectors of shape (2,) and arrays of vectors with shape
    (..., 2) where the last dimension must be 2.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        vectors = np.asarray(vectors, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vectors.ndim == 0:
        raise ValueError(f"{name} cannot be a scalar.")

    if vectors.shape[-1] != 3:
        raise ValueError(f"{name} must have 3 elements in the last dimension.")

    if not np.isfinite(vectors).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return vectors


def threeD_number_vectorLike_return_float_unit_vector(vector, name):
    """Validates a value is a 3D vector-like object (array-like object with shape (3,
    )). It then returns it as a (3,) numpy array of floats, normalized to have a
    magnitude of 1.0. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        vector = np.asarray(vector, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    if not np.isfinite(vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError(f"{name} must have a non-zero length.")
    elif not np.isclose(norm, 1.0):
        return vector / norm
    return vector


def fourD_homog_number_vectorLike_return_float(vector, name):
    """Validates a value is a 4D homogeneous vector-like object (array-like object
    with shape (4,) with a final value equal to 0.0 or 1.0). It then returns it as a
    (4,) numpy array of floats. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        vector = np.asarray(vector, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vector.shape != (4,):
        raise ValueError(f"{name} must be a 4-element vector.")

    if not np.isfinite(vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    last_val = vector[-1]
    if last_val != 0 and last_val != 1:
        raise ValueError(f"{name}'s last element must be 0.0 or 1.0.")

    return vector


def nD_number_vectorLike_return_float(vector, name):
    """Validates a value is an ND vector-like object (array-like object with shape (N,
    )). It then returns it as an (N,) numpy array of floats. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        vector = np.asarray(vector, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vector.ndim != 1:
        raise ValueError(f"{name} must be an N-element vector.")

    if not np.isfinite(vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return vector


def threeByThree_number_arrayLike_return_float(matrix, name):
    """Validates a value is a (3,3) array-like object. It then returns it as a (3,
    3) numpy array of floats. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        matrix = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix.")

    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return matrix


def fourByFour_number_arrayLike_return_float(matrix, name):
    """Validates a value is a (4,4) array-like object. It then returns it as a (4,
    4) numpy array of floats. name must be a string.

    Note: np.nan, np.inf, and -np.inf won't pass this test."""
    name = string_return_string(name, "name")

    try:
        matrix = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if matrix.shape != (4, 4):
        raise ValueError(f"{name} must be a 4x4 matrix.")

    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return matrix


def list_return_list(list_parameter, name):
    """Validates a list and returns it. name must be a string."""
    name = string_return_string(name, "name")

    if not isinstance(list_parameter, list):
        raise TypeError(f"{name} must be a list.")
    return list_parameter


def non_empty_list_return_list(list_parameter, name):
    """Validates a non-empty list and returns it. name must be a string."""
    name = string_return_string(name, "name")

    if not isinstance(list_parameter, list):
        raise TypeError(f"{name} must be a list.")

    if len(list_parameter) < 1:
        raise ValueError(f"{name} must have at least one element.")

    return list_parameter


def rotation_order_return_string(order, name):
    """Validates string representing a Tait-Bryan rotation sequence, and returns it.
    name must be a string."""
    name = string_return_string(name, "name")

    if not isinstance(order, str):
        raise TypeError(f"{name} must be a string.")

    if len(order) != 3:
        raise ValueError(f"{name} must have 3 characters.")

    valid_orders = ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]
    if order not in valid_orders:
        raise ValueError(f"{name} must be one of {valid_orders}.")

    return order
