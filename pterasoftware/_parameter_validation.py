"""Contains shared parameter validation functions."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, cast

import numpy as np


def str_return_str(value: Any, name: str) -> str:
    """Validates that a value is a str and returns it as a str.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value.
    """
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a str.")
    return value


def boolLike_return_bool(value: Any, name: str) -> bool:
    """Validates that a value is a bool or a numpy bool and returns it as a bool.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value.
    """
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    else:
        raise TypeError(
            f"{name} is a {type(value)} equal to {value} but it must be a bool or a "
            f"numpy bool."
        )


def int_in_range_return_int(
    value: Any,
    name: str,
    min_val: int | float | None = None,
    min_inclusive: bool | None = None,
    max_val: int | float | None = None,
    max_inclusive: bool | None = None,
) -> int:
    """Validates that a value is an int, that it lies in a custom range, and returns it
    as an int.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :param min_val: The minimum of the accepted range. If there is no minimum, leave as
        None. The default is None.
    :param min_inclusive: Determines whether the minimum of the range is a valid value.
        This must be None if min_val is None and a bool if min_val isn't None. The
        default is None.
    :param max_val: The maximum of the accepted range. If there is no maximum, leave as
        None. If both min_val and max_val are not None, then min_val must be less than
        max_val. The default is None.
    :param max_inclusive: Determines whether the maximum of the range is a valid value.
        This must be None if max_val is None and a bool if max_val isn't None. The
        default is None.
    :return: The validated value.
    """
    if min_val is not None:
        if min_inclusive is None:
            raise ValueError("min_inclusive can't be None if min_value isn't None.")
    else:
        if min_inclusive is not None:
            raise ValueError("min_inclusive must be None if min_val is None")

    if max_val is not None:
        if max_inclusive is None:
            raise ValueError("max_inclusive can't be None if max_value isn't None.")
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
                raise ValueError(f"{name} must be less than {max_val}.")

    return int(value)


def number_in_range_return_float(
    value: Any,
    name: str,
    min_val: int | float | None = None,
    min_inclusive: bool | None = None,
    max_val: int | float | None = None,
    max_inclusive: bool | None = None,
) -> float:
    """Validates that a value is a number (an int or a float), that it lies in a custom
    range, and returns it as a float.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :param min_val: The minimum of the accepted range. If there is no minimum, leave as
        None. The default is None.
    :param min_inclusive: Determines whether the minimum of the range is a valid value.
        This must be None if min_val is None and a bool if min_val isn't None. The
        default is None.
    :param max_val: The maximum of the accepted range. If there is no maximum, leave as
        None. If both min_val and max_val are not None, then min_val must be less than
        max_val. The default is None.
    :param max_inclusive: Determines whether the maximum of the range is a valid value.
        This must be None if max_val is None and a bool if max_val isn't None. The
        default is None.
    :return: The validated value.
    """
    if min_val is not None:
        if min_inclusive is None:
            raise ValueError("min_inclusive can't be None if min_value isn't None.")
    else:
        if min_inclusive is not None:
            raise ValueError("min_inclusive must be None if min_val is None")

    if max_val is not None:
        if max_inclusive is None:
            raise ValueError("max_inclusive can't be None if max_value isn't None.")
    else:
        if max_inclusive is not None:
            raise ValueError("max_inclusive must be None if max_val is None")

    if min_val is not None and max_val is not None:
        if min_val >= max_val:
            raise ValueError("min_val must be less than max_val")

    if not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be an int or a float.")

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
                raise ValueError(f"{name} must be less than {max_val}.")

    return float(value)


def arrayLike_of_numbers_in_range_return_float(
    value: Any,
    name: str,
    min_val: int | float | None = None,
    min_inclusive: bool | None = None,
    max_val: int | float | None = None,
    max_inclusive: bool | None = None,
) -> np.ndarray:
    """Validates a value is a number falling in a custom range, or an array-like object
    of numbers with every element falling in a custom range, and returns the input as a
    ndarray of floats.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :param min_val: The minimum of the accepted range. If there is no minimum, leave as
        None. The default is None.
    :param min_inclusive: Determines whether the minimum of the range is a valid value.
        This must be None if min_val is None and a bool if min_val isn't None. The
        default is None.
    :param max_val: The maximum of the accepted range. If there is no maximum, leave as
        None. If both min_val and max_val are not None, then min_val must be less than
        max_val. The default is None.
    :param max_inclusive: Determines whether the maximum of the range is a valid value.
        This must be None if max_val is None and a bool if max_val isn't None. The
        default is None.
    :return: The validated value. If the input value is a ndarray, this will be a
        ndarray of floats with same shape. If it is a number, this will be ndarray of
        floats with shape ().
    """
    if min_val is not None:
        if min_inclusive is None:
            raise ValueError("min_inclusive can't be None if min_value isn't None.")
    else:
        if min_inclusive is not None:
            raise ValueError("min_inclusive must be None if min_val is None")

    if max_val is not None:
        if max_inclusive is None:
            raise ValueError("max_inclusive can't be None if max_value isn't None.")
    else:
        if max_inclusive is not None:
            raise ValueError("max_inclusive must be None if max_val is None")

    if min_val is not None and max_val is not None:
        if min_val >= max_val:
            raise ValueError("min_val must be less than max_val")

    try:
        validated_value = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if not np.isfinite(validated_value).all():
        raise ValueError(f"{name} (or all its elements) can't be nan, inf, or -inf.")

    if min_val is not None:
        if min_inclusive:
            if not np.all(validated_value >= min_val):
                raise ValueError(
                    f"{name} (or all its elements) must be greater than or equal to {min_val}."
                )
        else:
            if not np.all(validated_value > min_val):
                raise ValueError(
                    f"{name} (or all its elements) must be greater than {min_val}."
                )

    if max_val is not None:
        if max_inclusive:
            if not np.all(validated_value <= max_val):
                raise ValueError(
                    f"{name} (or all its elements) must be less than or equal to {max_val}."
                )
        else:
            if not np.all(validated_value < max_val):
                raise ValueError(
                    f"{name} (or all its elements) must be less than {max_val}."
                )

    return validated_value


def arrayLike_of_twoD_number_vectorLikes_return_float(
    value: Any,
    name: str,
) -> np.ndarray:
    """Validates a value is an array-like object of 2D number vector-like objects
    (array-like objects with shape (2,)). It then returns it as a (...,2) ndarray of
    floats.

    Accepts both single vector-like objects of shape (2,) and array-like objects with
    shape (...,2). np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value as a ndarray of floats with the same shape as the input
        value.
    """
    try:
        validated_vectors = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if validated_vectors.ndim == 0:
        raise ValueError(f"{name} cannot be a scalar.")

    if validated_vectors.shape[-1] != 2:
        raise ValueError(f"{name} must have 2 elements in the last dimension.")

    if not np.isfinite(validated_vectors).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return validated_vectors


def threeD_number_vectorLike_return_float(value: Any, name: str) -> np.ndarray:
    """Validates a value is a 3D vector-like object (array-like object with shape (3,)),
    and returns it as a (3,) ndarray of floats.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value as a ndarray of floats with the same shape as the input
        value.
    """
    try:
        validated_vector = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if validated_vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    if not np.isfinite(validated_vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return validated_vector


def arrayLike_of_threeD_number_vectorLikes_return_float(
    value: Any, name: str
) -> np.ndarray:
    """Validates a value is an array-like object of 3D vector-like objects of numbers
    (array-like objects with shape (3,)) of numbers, and returns it as a (...,3) ndarray
    of floats.

    Accepts both single vector-like objects of shape (3,) and array-like objects with
    shape (...,3). np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value as a ndarray of floats with the same shape as the input
        value.
    """
    try:
        validated_vectors = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if validated_vectors.ndim == 0:
        raise ValueError(f"{name} cannot be a scalar.")

    if validated_vectors.shape[-1] != 3:
        raise ValueError(f"{name} must have 3 elements in the last dimension.")

    if not np.isfinite(validated_vectors).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return validated_vectors


def threeD_number_vectorLike_return_float_unit_vector(
    value: Any, name: str
) -> np.ndarray:
    """Validates a value is a 3D vector-like object (array-like object with shape (3,))
    of numbers, and returns it as a (3,) ndarray of floats, normalized to have a
    magnitude of 1.0.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value as a ndarray of floats with the same shape as the input
        value.
    """
    try:
        validated_vector = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if validated_vector.shape != (3,):
        raise ValueError(f"{name} must be a 3-element vector.")

    if not np.isfinite(validated_vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    norm = np.linalg.norm(validated_vector)
    if norm == 0:
        raise ValueError(f"{name} must have a non zero length.")
    elif not np.isclose(norm, 1.0):
        return cast(np.ndarray, validated_vector / norm)
    return validated_vector


def threeD_spacing_vectorLike_return_tuple(value: Any, name: str) -> tuple[
    str | Callable[[np.ndarray], np.ndarray],
    str | Callable[[np.ndarray], np.ndarray],
    str | Callable[[np.ndarray], np.ndarray],
]:
    """Validates a value is a 3D vector-like object (array-like object with shape (3,))
    of spacing specifications, and then returns it as a tuple of 3 spacing
    specifications.

    Each element can either be a str ("sine" or "uniform") or a callable (custom spacing
    function).

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value.
    """
    # Convert to list if ndarray, or validate it's a list/tuple.
    if isinstance(value, np.ndarray):
        if value.ndim != 1 or value.shape[0] != 3:
            raise ValueError(f"{name} must be a 3-element vector.")
        value = value.tolist()
    elif isinstance(value, (list, tuple)):
        if len(value) != 3:
            raise ValueError(f"{name} must be a 3-element vector.")
    else:
        raise TypeError(f"{name} must be array-like (tuple, list, or ndarray).")

    # Check each element is either a valid str or a callable.
    validated_list = []
    for i, elem in enumerate(value):
        if isinstance(elem, str):
            if elem not in ["sine", "uniform"]:
                raise ValueError(
                    f"Element {i} of {name} must be 'sine', 'uniform', or a callable, "
                    f"got str '{elem}'."
                )
            validated_list.append(elem)
        elif callable(elem):
            validated_list.append(elem)
        else:
            raise TypeError(
                f"Element {i} of {name} must be a str ('sine' or 'uniform') or a "
                f"callable, got {type(elem).__name__}."
            )

    validated_value = tuple(validated_list)
    return cast(
        tuple[
            str | Callable[[np.ndarray], np.ndarray],
            str | Callable[[np.ndarray], np.ndarray],
            str | Callable[[np.ndarray], np.ndarray],
        ],
        validated_value,
    )


def nD_number_vectorLike_return_float(value: Any, name: str) -> np.ndarray:
    """Validates a value is an ND vector-like object (array-like object with shape (N,))
    of numbers, and returns it as an (N,) ndarray of floats.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value as a ndarray of floats with the same shape as the input
        value.
    """
    try:
        validated_vector = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if validated_vector.ndim != 1:
        raise ValueError(f"{name} must be an N-element vector.")

    if not np.isfinite(validated_vector).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return validated_vector


def fourByFour_number_arrayLike_return_float(value: Any, name: str) -> np.ndarray:
    """Validates a value is a (4,4) array-like object. It then returns it as a (4,4)
    ndarray of floats.

    np.nan, np.inf, and -np.inf aren't valid values.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value as a ndarray of floats with the same shape as the input
        value.
    """
    try:
        validated_matrix = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be array-like and contain ints or floats.")

    if validated_matrix.shape != (4, 4):
        raise ValueError(f"{name} must be a 4x4 matrix.")

    if not np.isfinite(validated_matrix).all():
        raise ValueError(f"{name} can't contain any nan, inf, or -inf elements.")

    return validated_matrix


def non_empty_list_return_list(value: Any, name: str) -> list[Any]:
    """Validates a value is a non empty list and returns it.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value.
    """
    if not isinstance(value, list):
        raise TypeError(f"{name} must be a list.")

    if len(value) < 1:
        raise ValueError(f"{name} must have at least one element.")

    return value


def rotation_order_return_str(value: Any, name: str) -> str:
    """Validates a value is a str representing a Tait-Bryan rotation sequence, and
    returns it as a str.

    :param value: The value to validate.
    :param name: The name of the value.
    :return: The validated value.
    """
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a str.")

    if len(value) != 3:
        raise ValueError(f"{name} must have 3 characters.")

    valid_orders = ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]
    if value not in valid_orders:
        raise ValueError(f"{name} must be one of {valid_orders}.")

    return value
