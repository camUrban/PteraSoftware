"""This module contains common parameter validation functions."""

import numpy as np


def validate_boolean(value, name):
    """Validates a boolean parameter."""
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a boolean.")

    return value


def validate_scalar_int(value, name):
    """Validates a scalar is an int and returns it."""
    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    return value


def validate_non_negative_scalar_int(value, name):
    """Validates a scalar is an int with a value greater than or equal to zero and
    returns it."""

    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if value < 0:
        raise ValueError(f"{name} must be greater than or equal to zero.")

    return value


def validate_positive_scalar_int(value, name):
    """Validates a scalar is an int with a value greater than zero and returns it."""
    if not isinstance(value, int):
        raise TypeError(f"{name} must be an int.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if value <= 0:
        raise ValueError(f"{name} must be positive.")

    return value


def validate_scalar_float(value, name):
    """Validates a scalar value and returns it as a float."""
    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    return float(value)


def validate_non_negative_scalar_float(value, name):
    """Validates a scalar value greater than or equal to zero and returns it as a
    float."""
    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if value < 0:
        raise ValueError(f"{name} must be greater than or equal to zero.")

    return float(value)


def validate_positive_scalar_float(value, name):
    """Validates a scalar value greater than zero and returns it as a float."""
    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if value <= 0:
        raise ValueError(f"{name} must be positive.")

    return float(value)


def validate_scalar_in_range_float(
    value, name, min_val, min_inclusive, max_val, max_inclusive
):
    """Validates a scalar value in a custom range and returns it as a float."""
    if not isinstance(value, (int, float, np.number)):
        raise TypeError(f"{name} must be numeric.")

    if not np.isfinite(value).all():
        raise ValueError(f"{name} can't be nan, inf, or -inf.")

    if min_inclusive:
        if not value >= min_val:
            raise ValueError(f"{name} must be greater than or equal to {min_val}.")
    else:
        if not value > min_val:
            raise ValueError(f"{name} must be greater than {min_val}.")

    if max_inclusive:
        if not value <= max_val:
            raise ValueError(f"{name} must be less than or equal to {max_val}.")
    else:
        if not value < max_val:
            raise ValueError(f"{name} must be less than {min_val}.")

    return float(value)


def validate_3d_vector_float(vector, name):
    """Validates a 3D vector and returns it as a vector of floats."""
    try:
        vector = np.asarray(vector, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    if not np.isfinite(vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return vector


def validate_3d_vectors_array_float(vectors, name):
    """Validates an array of 3D vectors and returns it as an array of floats.
    
    Accepts both single vectors of shape (3,) and arrays of vectors with shape 
    (..., 3) where the last dimension must be 3.
    """
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


def validate_3d_unit_vector_norm_float(vector, name):
    """Validates a 3D unit vector and returns it as a normalized vector of floats."""
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


def validate_4d_homog_vector_float(vector, name):
    """Validates a 4D vector in homogeneous coordinates and returns it as a vector of
    floats."""
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


def validate_3_by_3_matrix_float(matrix, name):
    """Validates a 3x3 matrix and returns it as a matrix of floats."""
    try:
        matrix = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix.")

    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return matrix


def validate_4_by_4_matrix_float(matrix, name):
    """Validates a 4x4 matrix and returns it as a matrix of floats."""
    try:
        matrix = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and numeric.")

    if matrix.shape != (4, 4):
        raise ValueError(f"{name} must be a 4x4 matrix.")

    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return matrix


def validate_string(string, name):
    """Validates a string and returns it."""
    if not isinstance(string, str):
        raise TypeError(f"{name} must be a string.")
    return string


def validate_list(list_parameter, name):
    """Validates a list and returns it."""
    if not isinstance(list_parameter, list):
        raise TypeError(f"{name} must be a list.")
    return list_parameter


def validate_non_empty_list(list_parameter, name):
    """Validates a non-empty list and returns it."""
    if not isinstance(list_parameter, list):
        raise TypeError(f"{name} must be a list.")

    if len(list_parameter) < 1:
        raise ValueError(f"{name} must have at least one element.")

    return list_parameter


def validate_rotation_order(order, name):
    """Validates string representing a Tait-Bryan rotation sequence, and returns it."""
    if not isinstance(order, str):
        raise TypeError(f"{name} must be a string.")

    if len(order) != 3:
        raise ValueError(f"{name} must have 3 characters.")

    valid_orders = ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]
    if order not in valid_orders:
        raise ValueError(f"{name} must be one of {valid_orders}.")

    return order
