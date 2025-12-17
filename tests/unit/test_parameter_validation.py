"""This module contains a class to test parameter validation functions."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _parameter_validation as pv
from tests.unit.fixtures import parameter_validation_fixtures as pvf


class TestStrReturnStr(unittest.TestCase):
    """A class with functions to test str_return_str."""

    def test_valid_str(self):
        """Test str_return_str with a valid string."""
        valid_str = pvf.make_valid_str_fixture()
        result = pv.str_return_str(valid_str, "test_param")
        self.assertEqual(result, valid_str)
        self.assertIsInstance(result, str)

    def test_valid_empty_str(self):
        """Test str_return_str with an empty string (still valid type)."""
        empty_str = pvf.make_empty_str_fixture()
        result = pv.str_return_str(empty_str, "test_param")
        self.assertEqual(result, "")
        self.assertIsInstance(result, str)

    def test_invalid_int(self):
        """Test str_return_str raises TypeError with an int."""
        invalid_value = pvf.make_invalid_non_str_fixture()
        with self.assertRaises(TypeError) as context:
            pv.str_return_str(invalid_value, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("str", str(context.exception))

    def test_invalid_none(self):
        """Test str_return_str raises TypeError with None."""
        with self.assertRaises(TypeError):
            pv.str_return_str(None, "test_param")

    def test_invalid_list(self):
        """Test str_return_str raises TypeError with a list."""
        with self.assertRaises(TypeError):
            pv.str_return_str(["a", "b"], "test_param")


class TestBoolLikeReturnBool(unittest.TestCase):
    """A class with functions to test boolLike_return_bool."""

    def test_valid_bool_true(self):
        """Test boolLike_return_bool with True."""
        result = pv.boolLike_return_bool(True, "test_param")
        self.assertTrue(result)
        self.assertIsInstance(result, bool)

    def test_valid_bool_false(self):
        """Test boolLike_return_bool with False."""
        result = pv.boolLike_return_bool(False, "test_param")
        self.assertFalse(result)
        self.assertIsInstance(result, bool)

    def test_valid_numpy_bool_true(self):
        """Test boolLike_return_bool with numpy bool True."""
        np_bool = np.bool_(True)
        result = pv.boolLike_return_bool(np_bool, "test_param")
        self.assertTrue(result)
        self.assertIsInstance(result, bool)

    def test_valid_numpy_bool_false(self):
        """Test boolLike_return_bool with numpy bool False."""
        np_bool = pvf.make_valid_numpy_bool_fixture()
        result = pv.boolLike_return_bool(np_bool, "test_param")
        self.assertFalse(result)
        self.assertIsInstance(result, bool)

    def test_invalid_string(self):
        """Test boolLike_return_bool raises TypeError with a string."""
        invalid_value = pvf.make_invalid_non_bool_fixture()
        with self.assertRaises(TypeError) as context:
            pv.boolLike_return_bool(invalid_value, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("bool", str(context.exception))

    def test_invalid_int(self):
        """Test boolLike_return_bool raises TypeError with an int."""
        with self.assertRaises(TypeError):
            pv.boolLike_return_bool(1, "test_param")

    def test_invalid_none(self):
        """Test boolLike_return_bool raises TypeError with None."""
        with self.assertRaises(TypeError):
            pv.boolLike_return_bool(None, "test_param")


class TestIntInRangeReturnInt(unittest.TestCase):
    """A class with functions to test int_in_range_return_int."""

    def test_valid_int_no_bounds(self):
        """Test int_in_range_return_int with a valid int and no bounds."""
        valid_int = pvf.make_valid_int_fixture()
        result = pv.int_in_range_return_int(valid_int, "test_param")
        self.assertEqual(result, valid_int)
        self.assertIsInstance(result, int)

    def test_valid_int_with_min_inclusive(self):
        """Test int_in_range_return_int with min bound inclusive."""
        result = pv.int_in_range_return_int(
            5, "test_param", min_val=0, min_inclusive=True
        )
        self.assertEqual(result, 5)

    def test_valid_int_at_min_inclusive(self):
        """Test int_in_range_return_int at min bound (inclusive)."""
        result = pv.int_in_range_return_int(
            0, "test_param", min_val=0, min_inclusive=True
        )
        self.assertEqual(result, 0)

    def test_valid_int_with_min_exclusive(self):
        """Test int_in_range_return_int with min bound exclusive."""
        result = pv.int_in_range_return_int(
            1, "test_param", min_val=0, min_inclusive=False
        )
        self.assertEqual(result, 1)

    def test_valid_int_with_max_inclusive(self):
        """Test int_in_range_return_int with max bound inclusive."""
        result = pv.int_in_range_return_int(
            5, "test_param", max_val=10, max_inclusive=True
        )
        self.assertEqual(result, 5)

    def test_valid_int_at_max_inclusive(self):
        """Test int_in_range_return_int at max bound (inclusive)."""
        result = pv.int_in_range_return_int(
            10, "test_param", max_val=10, max_inclusive=True
        )
        self.assertEqual(result, 10)

    def test_valid_int_with_max_exclusive(self):
        """Test int_in_range_return_int with max bound exclusive."""
        result = pv.int_in_range_return_int(
            9, "test_param", max_val=10, max_inclusive=False
        )
        self.assertEqual(result, 9)

    def test_valid_int_with_both_bounds_inclusive(self):
        """Test int_in_range_return_int with both bounds inclusive."""
        result = pv.int_in_range_return_int(
            5,
            "test_param",
            min_val=0,
            min_inclusive=True,
            max_val=10,
            max_inclusive=True,
        )
        self.assertEqual(result, 5)

    def test_valid_int_with_both_bounds_exclusive(self):
        """Test int_in_range_return_int with both bounds exclusive."""
        result = pv.int_in_range_return_int(
            5,
            "test_param",
            min_val=0,
            min_inclusive=False,
            max_val=10,
            max_inclusive=False,
        )
        self.assertEqual(result, 5)

    def test_invalid_float(self):
        """Test int_in_range_return_int raises TypeError with a float."""
        invalid_value = pvf.make_invalid_float_for_int_fixture()
        with self.assertRaises(TypeError) as context:
            pv.int_in_range_return_int(invalid_value, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("int", str(context.exception))

    def test_invalid_string(self):
        """Test int_in_range_return_int raises TypeError with a string."""
        with self.assertRaises(TypeError):
            pv.int_in_range_return_int("5", "test_param")

    def test_invalid_below_min_inclusive(self):
        """Test int_in_range_return_int raises ValueError below min (inclusive)."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(-1, "test_param", min_val=0, min_inclusive=True)
        self.assertIn("test_param", str(context.exception))
        self.assertIn("greater than or equal to", str(context.exception))

    def test_invalid_at_min_exclusive(self):
        """Test int_in_range_return_int raises ValueError at min (exclusive)."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(0, "test_param", min_val=0, min_inclusive=False)
        self.assertIn("test_param", str(context.exception))
        self.assertIn("greater than", str(context.exception))

    def test_invalid_above_max_inclusive(self):
        """Test int_in_range_return_int raises ValueError above max (inclusive)."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(11, "test_param", max_val=10, max_inclusive=True)
        self.assertIn("test_param", str(context.exception))
        self.assertIn("less than or equal to", str(context.exception))

    def test_invalid_at_max_exclusive(self):
        """Test int_in_range_return_int raises ValueError at max (exclusive)."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(
                10, "test_param", max_val=10, max_inclusive=False
            )
        self.assertIn("test_param", str(context.exception))
        self.assertIn("less than", str(context.exception))

    def test_invalid_min_inclusive_without_min_val(self):
        """Test int_in_range_return_int raises ValueError with min_inclusive but no min_val."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(5, "test_param", min_inclusive=True)
        self.assertIn(
            "min_inclusive must be None if min_val is None", str(context.exception)
        )

    def test_invalid_min_val_without_min_inclusive(self):
        """Test int_in_range_return_int raises ValueError with min_val but no min_inclusive."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(5, "test_param", min_val=0)
        self.assertIn(
            "min_inclusive can't be None if min_value isn't None",
            str(context.exception),
        )

    def test_invalid_max_inclusive_without_max_val(self):
        """Test int_in_range_return_int raises ValueError with max_inclusive but no max_val."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(5, "test_param", max_inclusive=True)
        self.assertIn(
            "max_inclusive must be None if max_val is None", str(context.exception)
        )

    def test_invalid_max_val_without_max_inclusive(self):
        """Test int_in_range_return_int raises ValueError with max_val but no max_inclusive."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(5, "test_param", max_val=10)
        self.assertIn(
            "max_inclusive can't be None if max_value isn't None",
            str(context.exception),
        )

    def test_invalid_min_val_greater_than_max_val(self):
        """Test int_in_range_return_int raises ValueError when min_val >= max_val."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(
                5,
                "test_param",
                min_val=10,
                min_inclusive=True,
                max_val=5,
                max_inclusive=True,
            )
        self.assertIn("min_val must be less than max_val", str(context.exception))

    def test_invalid_min_val_equals_max_val(self):
        """Test int_in_range_return_int raises ValueError when min_val == max_val."""
        with self.assertRaises(ValueError) as context:
            pv.int_in_range_return_int(
                5,
                "test_param",
                min_val=5,
                min_inclusive=True,
                max_val=5,
                max_inclusive=True,
            )
        self.assertIn("min_val must be less than max_val", str(context.exception))


class TestNumberInRangeReturnFloat(unittest.TestCase):
    """A class with functions to test number_in_range_return_float."""

    def test_valid_float_no_bounds(self):
        """Test number_in_range_return_float with a valid float and no bounds."""
        valid_float = pvf.make_valid_float_fixture()
        result = pv.number_in_range_return_float(valid_float, "test_param")
        self.assertEqual(result, valid_float)
        self.assertIsInstance(result, float)

    def test_valid_int_no_bounds(self):
        """Test number_in_range_return_float with a valid int (converted to float)."""
        valid_int = pvf.make_valid_int_as_number_fixture()
        result = pv.number_in_range_return_float(valid_int, "test_param")
        self.assertEqual(result, float(valid_int))
        self.assertIsInstance(result, float)

    def test_valid_float_with_min_inclusive(self):
        """Test number_in_range_return_float with min bound inclusive."""
        result = pv.number_in_range_return_float(
            5.5, "test_param", min_val=0.0, min_inclusive=True
        )
        self.assertEqual(result, 5.5)

    def test_valid_float_at_min_inclusive(self):
        """Test number_in_range_return_float at min bound (inclusive)."""
        result = pv.number_in_range_return_float(
            0.0, "test_param", min_val=0.0, min_inclusive=True
        )
        self.assertEqual(result, 0.0)

    def test_valid_float_with_max_inclusive(self):
        """Test number_in_range_return_float with max bound inclusive."""
        result = pv.number_in_range_return_float(
            5.5, "test_param", max_val=10.0, max_inclusive=True
        )
        self.assertEqual(result, 5.5)

    def test_valid_float_at_max_inclusive(self):
        """Test number_in_range_return_float at max bound (inclusive)."""
        result = pv.number_in_range_return_float(
            10.0, "test_param", max_val=10.0, max_inclusive=True
        )
        self.assertEqual(result, 10.0)

    def test_valid_float_with_both_bounds(self):
        """Test number_in_range_return_float with both bounds."""
        result = pv.number_in_range_return_float(
            5.5,
            "test_param",
            min_val=0.0,
            min_inclusive=True,
            max_val=10.0,
            max_inclusive=True,
        )
        self.assertEqual(result, 5.5)

    def test_invalid_string(self):
        """Test number_in_range_return_float raises TypeError with a string."""
        with self.assertRaises(TypeError) as context:
            pv.number_in_range_return_float("5.5", "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("int or a float", str(context.exception))

    def test_invalid_nan(self):
        """Test number_in_range_return_float raises ValueError with NaN."""
        invalid_value = pvf.make_invalid_nan_fixture()
        with self.assertRaises(ValueError) as context:
            pv.number_in_range_return_float(invalid_value, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("nan", str(context.exception).lower())

    def test_invalid_inf(self):
        """Test number_in_range_return_float raises ValueError with infinity."""
        invalid_value = pvf.make_invalid_inf_fixture()
        with self.assertRaises(ValueError) as context:
            pv.number_in_range_return_float(invalid_value, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("inf", str(context.exception).lower())

    def test_invalid_neg_inf(self):
        """Test number_in_range_return_float raises ValueError with negative infinity."""
        invalid_value = pvf.make_invalid_neg_inf_fixture()
        with self.assertRaises(ValueError) as context:
            pv.number_in_range_return_float(invalid_value, "test_param")
        self.assertIn("test_param", str(context.exception))

    def test_invalid_below_min(self):
        """Test number_in_range_return_float raises ValueError below min."""
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(
                -0.1, "test_param", min_val=0.0, min_inclusive=True
            )

    def test_invalid_above_max(self):
        """Test number_in_range_return_float raises ValueError above max."""
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(
                10.1, "test_param", max_val=10.0, max_inclusive=True
            )

    def test_invalid_parameter_combinations(self):
        """Test number_in_range_return_float raises ValueError with invalid parameter combinations."""
        # min_inclusive without min_val
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(5.0, "test_param", min_inclusive=True)

        # min_val without min_inclusive
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(5.0, "test_param", min_val=0.0)

        # max_inclusive without max_val
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(5.0, "test_param", max_inclusive=True)

        # max_val without max_inclusive
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(5.0, "test_param", max_val=10.0)

        # min_val >= max_val
        with self.assertRaises(ValueError):
            pv.number_in_range_return_float(
                5.0,
                "test_param",
                min_val=10.0,
                min_inclusive=True,
                max_val=5.0,
                max_inclusive=True,
            )


class TestArrayLikeOfNumbersInRangeReturnFloat(unittest.TestCase):
    """A class with functions to test arrayLike_of_numbers_in_range_return_float."""

    def test_valid_1d_array(self):
        """Test arrayLike_of_numbers_in_range_return_float with a valid 1D array."""
        valid_array = pvf.make_valid_1d_array_fixture()
        result = pv.arrayLike_of_numbers_in_range_return_float(
            valid_array, "test_param"
        )
        npt.assert_array_equal(result, valid_array)
        self.assertEqual(result.dtype, float)

    def test_valid_list(self):
        """Test arrayLike_of_numbers_in_range_return_float with a valid list."""
        valid_list = pvf.make_valid_list_fixture()
        result = pv.arrayLike_of_numbers_in_range_return_float(valid_list, "test_param")
        npt.assert_array_equal(result, np.array(valid_list, dtype=float))
        self.assertEqual(result.dtype, float)

    def test_valid_tuple(self):
        """Test arrayLike_of_numbers_in_range_return_float with a valid tuple."""
        valid_tuple = pvf.make_valid_tuple_fixture()
        result = pv.arrayLike_of_numbers_in_range_return_float(
            valid_tuple, "test_param"
        )
        npt.assert_array_equal(result, np.array(valid_tuple, dtype=float))
        self.assertEqual(result.dtype, float)

    def test_valid_nested_list(self):
        """Test arrayLike_of_numbers_in_range_return_float with a valid nested list."""
        valid_nested = pvf.make_valid_nested_list_fixture()
        result = pv.arrayLike_of_numbers_in_range_return_float(
            valid_nested, "test_param"
        )
        npt.assert_array_equal(result, np.array(valid_nested, dtype=float))

    def test_valid_scalar(self):
        """Test arrayLike_of_numbers_in_range_return_float with a scalar."""
        result = pv.arrayLike_of_numbers_in_range_return_float(5.0, "test_param")
        self.assertEqual(result, 5.0)
        self.assertEqual(result.shape, ())

    def test_valid_with_range(self):
        """Test arrayLike_of_numbers_in_range_return_float with valid range."""
        result = pv.arrayLike_of_numbers_in_range_return_float(
            [1.0, 2.0, 3.0],
            "test_param",
            min_val=0.0,
            min_inclusive=True,
            max_val=10.0,
            max_inclusive=True,
        )
        npt.assert_array_equal(result, np.array([1.0, 2.0, 3.0], dtype=float))

    def test_invalid_string_array(self):
        """Test arrayLike_of_numbers_in_range_return_float raises TypeError with string array."""
        with self.assertRaises(TypeError) as context:
            pv.arrayLike_of_numbers_in_range_return_float(["a", "b", "c"], "test_param")
        self.assertIn("test_param", str(context.exception))

    def test_invalid_array_with_nan(self):
        """Test arrayLike_of_numbers_in_range_return_float raises ValueError with NaN."""
        invalid_array = pvf.make_invalid_array_with_nan_fixture()
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_numbers_in_range_return_float(invalid_array, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("nan", str(context.exception).lower())

    def test_invalid_array_with_inf(self):
        """Test arrayLike_of_numbers_in_range_return_float raises ValueError with infinity."""
        invalid_array = pvf.make_invalid_array_with_inf_fixture()
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_numbers_in_range_return_float(invalid_array, "test_param")
        self.assertIn("test_param", str(context.exception))

    def test_invalid_below_min(self):
        """Test arrayLike_of_numbers_in_range_return_float raises ValueError below min."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_numbers_in_range_return_float(
                [-1.0, 2.0, 3.0],
                "test_param",
                min_val=0.0,
                min_inclusive=True,
            )
        self.assertIn("test_param", str(context.exception))

    def test_invalid_above_max(self):
        """Test arrayLike_of_numbers_in_range_return_float raises ValueError above max."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_numbers_in_range_return_float(
                [1.0, 2.0, 11.0],
                "test_param",
                max_val=10.0,
                max_inclusive=True,
            )
        self.assertIn("test_param", str(context.exception))


class TestArrayLikeOfTwoDNumberVectorLikesReturnFloat(unittest.TestCase):
    """A class with functions to test arrayLike_of_twoD_number_vectorLikes_return_float."""

    def test_valid_single_2d_vector(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float with a single 2D vector."""
        valid_vector = pvf.make_valid_2d_vector_fixture()
        result = pv.arrayLike_of_twoD_number_vectorLikes_return_float(
            valid_vector, "test_param"
        )
        npt.assert_array_equal(result, valid_vector)
        self.assertEqual(result.dtype, float)

    def test_valid_array_of_2d_vectors(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float with array of 2D vectors."""
        valid_array = pvf.make_valid_2d_vectors_array_fixture()
        result = pv.arrayLike_of_twoD_number_vectorLikes_return_float(
            valid_array, "test_param"
        )
        npt.assert_array_equal(result, valid_array)
        self.assertEqual(result.shape, (3, 2))

    def test_valid_list_of_2d_vectors(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float with list of 2D vectors."""
        result = pv.arrayLike_of_twoD_number_vectorLikes_return_float(
            [[1.0, 2.0], [3.0, 4.0]], "test_param"
        )
        npt.assert_array_equal(result, np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float))

    def test_invalid_scalar(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float raises ValueError with scalar."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_twoD_number_vectorLikes_return_float(5.0, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("scalar", str(context.exception).lower())

    def test_invalid_3d_vector(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float raises ValueError with 3D vector."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_twoD_number_vectorLikes_return_float(
                [1.0, 2.0, 3.0], "test_param"
            )
        self.assertIn("test_param", str(context.exception))
        self.assertIn("2 elements", str(context.exception))

    def test_invalid_with_nan(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float raises ValueError with NaN."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_twoD_number_vectorLikes_return_float(
                [1.0, np.nan], "test_param"
            )
        self.assertIn("test_param", str(context.exception))
        self.assertIn("nan", str(context.exception).lower())

    def test_invalid_string(self):
        """Test arrayLike_of_twoD_number_vectorLikes_return_float raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.arrayLike_of_twoD_number_vectorLikes_return_float(
                "invalid", "test_param"
            )


class TestThreeDNumberVectorLikeReturnFloat(unittest.TestCase):
    """A class with functions to test threeD_number_vectorLike_return_float."""

    def test_valid_3d_array(self):
        """Test threeD_number_vectorLike_return_float with a valid 3D array."""
        valid_vector = pvf.make_valid_3d_vector_fixture()
        result = pv.threeD_number_vectorLike_return_float(valid_vector, "test_param")
        npt.assert_array_equal(result, valid_vector)
        self.assertEqual(result.dtype, float)
        self.assertEqual(result.shape, (3,))

    def test_valid_3d_list(self):
        """Test threeD_number_vectorLike_return_float with a valid 3D list."""
        valid_list = pvf.make_valid_3d_list_fixture()
        result = pv.threeD_number_vectorLike_return_float(valid_list, "test_param")
        npt.assert_array_equal(result, np.array(valid_list, dtype=float))
        self.assertEqual(result.shape, (3,))

    def test_valid_3d_tuple(self):
        """Test threeD_number_vectorLike_return_float with a valid 3D tuple."""
        valid_tuple = pvf.make_valid_3d_tuple_fixture()
        result = pv.threeD_number_vectorLike_return_float(valid_tuple, "test_param")
        npt.assert_array_equal(result, np.array(valid_tuple, dtype=float))
        self.assertEqual(result.shape, (3,))

    def test_valid_with_ints(self):
        """Test threeD_number_vectorLike_return_float with integers (converted to float)."""
        result = pv.threeD_number_vectorLike_return_float([1, 2, 3], "test_param")
        npt.assert_array_equal(result, np.array([1.0, 2.0, 3.0], dtype=float))
        self.assertEqual(result.dtype, float)

    def test_invalid_2d_vector(self):
        """Test threeD_number_vectorLike_return_float raises ValueError with 2D vector."""
        invalid_vector = pvf.make_invalid_2d_vector_as_3d_fixture()
        with self.assertRaises(ValueError) as context:
            pv.threeD_number_vectorLike_return_float(invalid_vector, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("3-element", str(context.exception))

    def test_invalid_4d_vector(self):
        """Test threeD_number_vectorLike_return_float raises ValueError with 4D vector."""
        invalid_vector = pvf.make_invalid_4d_vector_as_3d_fixture()
        with self.assertRaises(ValueError) as context:
            pv.threeD_number_vectorLike_return_float(invalid_vector, "test_param")
        self.assertIn("test_param", str(context.exception))

    def test_invalid_with_nan(self):
        """Test threeD_number_vectorLike_return_float raises ValueError with NaN."""
        with self.assertRaises(ValueError) as context:
            pv.threeD_number_vectorLike_return_float([1.0, np.nan, 3.0], "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("nan", str(context.exception).lower())

    def test_invalid_with_inf(self):
        """Test threeD_number_vectorLike_return_float raises ValueError with infinity."""
        with self.assertRaises(ValueError):
            pv.threeD_number_vectorLike_return_float([1.0, np.inf, 3.0], "test_param")

    def test_invalid_string(self):
        """Test threeD_number_vectorLike_return_float raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.threeD_number_vectorLike_return_float("invalid", "test_param")


class TestArrayLikeOfThreeDNumberVectorLikesReturnFloat(unittest.TestCase):
    """A class with functions to test arrayLike_of_threeD_number_vectorLikes_return_float."""

    def test_valid_single_3d_vector(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float with single 3D vector."""
        valid_vector = pvf.make_valid_3d_vector_fixture()
        result = pv.arrayLike_of_threeD_number_vectorLikes_return_float(
            valid_vector, "test_param"
        )
        npt.assert_array_equal(result, valid_vector)
        self.assertEqual(result.shape, (3,))

    def test_valid_array_of_3d_vectors(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float with array of 3D vectors."""
        valid_array = pvf.make_valid_3d_vectors_array_fixture()
        result = pv.arrayLike_of_threeD_number_vectorLikes_return_float(
            valid_array, "test_param"
        )
        npt.assert_array_equal(result, valid_array)
        self.assertEqual(result.shape, (4, 3))

    def test_valid_list_of_3d_vectors(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float with list of 3D vectors."""
        result = pv.arrayLike_of_threeD_number_vectorLikes_return_float(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], "test_param"
        )
        expected = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=float)
        npt.assert_array_equal(result, expected)

    def test_invalid_scalar(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float raises ValueError with scalar."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_threeD_number_vectorLikes_return_float(5.0, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("scalar", str(context.exception).lower())

    def test_invalid_2d_vectors(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float raises ValueError with 2D vectors."""
        with self.assertRaises(ValueError) as context:
            pv.arrayLike_of_threeD_number_vectorLikes_return_float(
                [[1.0, 2.0], [3.0, 4.0]], "test_param"
            )
        self.assertIn("test_param", str(context.exception))
        self.assertIn("3 elements", str(context.exception))

    def test_invalid_with_nan(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float raises ValueError with NaN."""
        with self.assertRaises(ValueError):
            pv.arrayLike_of_threeD_number_vectorLikes_return_float(
                [[1.0, np.nan, 3.0]], "test_param"
            )

    def test_invalid_string(self):
        """Test arrayLike_of_threeD_number_vectorLikes_return_float raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.arrayLike_of_threeD_number_vectorLikes_return_float(
                "invalid", "test_param"
            )


class TestThreeDNumberVectorLikeReturnFloatUnitVector(unittest.TestCase):
    """A class with functions to test threeD_number_vectorLike_return_float_unit_vector."""

    def test_valid_unit_vector(self):
        """Test threeD_number_vectorLike_return_float_unit_vector with already unit vector."""
        valid_vector = pvf.make_valid_unit_vector_fixture()
        result = pv.threeD_number_vectorLike_return_float_unit_vector(
            valid_vector, "test_param"
        )
        npt.assert_allclose(result, valid_vector, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(np.linalg.norm(result), 1.0, rtol=1e-10, atol=1e-14)

    def test_valid_non_unit_vector_normalized(self):
        """Test threeD_number_vectorLike_return_float_unit_vector normalizes non-unit vector."""
        valid_vector = pvf.make_valid_non_unit_vector_fixture()
        result = pv.threeD_number_vectorLike_return_float_unit_vector(
            valid_vector, "test_param"
        )
        expected = valid_vector / np.linalg.norm(valid_vector)
        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(np.linalg.norm(result), 1.0, rtol=1e-10, atol=1e-14)

    def test_valid_list_normalized(self):
        """Test threeD_number_vectorLike_return_float_unit_vector normalizes list input."""
        result = pv.threeD_number_vectorLike_return_float_unit_vector(
            [0.0, 3.0, 4.0], "test_param"
        )
        expected = np.array([0.0, 0.6, 0.8], dtype=float)
        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_invalid_zero_vector(self):
        """Test threeD_number_vectorLike_return_float_unit_vector raises ValueError with zero vector."""
        zero_vector = pvf.make_invalid_zero_vector_fixture()
        with self.assertRaises(ValueError) as context:
            pv.threeD_number_vectorLike_return_float_unit_vector(
                zero_vector, "test_param"
            )
        self.assertIn("test_param", str(context.exception))
        self.assertIn("non zero length", str(context.exception))

    def test_invalid_2d_vector(self):
        """Test threeD_number_vectorLike_return_float_unit_vector raises ValueError with 2D vector."""
        invalid_vector = pvf.make_invalid_2d_vector_as_3d_fixture()
        with self.assertRaises(ValueError):
            pv.threeD_number_vectorLike_return_float_unit_vector(
                invalid_vector, "test_param"
            )

    def test_invalid_with_nan(self):
        """Test threeD_number_vectorLike_return_float_unit_vector raises ValueError with NaN."""
        with self.assertRaises(ValueError):
            pv.threeD_number_vectorLike_return_float_unit_vector(
                [1.0, np.nan, 3.0], "test_param"
            )

    def test_invalid_string(self):
        """Test threeD_number_vectorLike_return_float_unit_vector raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.threeD_number_vectorLike_return_float_unit_vector(
                "invalid", "test_param"
            )


class TestThreeDSpacingVectorLikeReturnTuple(unittest.TestCase):
    """A class with functions to test threeD_spacing_vectorLike_return_tuple."""

    def test_valid_all_strings(self):
        """Test threeD_spacing_vectorLike_return_tuple with all string elements."""
        valid_spacing = pvf.make_valid_spacing_vector_all_strings_fixture()
        result = pv.threeD_spacing_vectorLike_return_tuple(valid_spacing, "test_param")
        self.assertEqual(result, valid_spacing)
        self.assertIsInstance(result, tuple)

    def test_valid_all_callables(self):
        """Test threeD_spacing_vectorLike_return_tuple with all callable elements."""
        valid_spacing = pvf.make_valid_spacing_vector_all_callables_fixture()
        result = pv.threeD_spacing_vectorLike_return_tuple(valid_spacing, "test_param")
        self.assertEqual(len(result), 3)
        for elem in result:
            self.assertTrue(callable(elem))

    def test_valid_mixed(self):
        """Test threeD_spacing_vectorLike_return_tuple with mixed string and callable."""
        valid_spacing = pvf.make_valid_spacing_vector_mixed_fixture()
        result = pv.threeD_spacing_vectorLike_return_tuple(valid_spacing, "test_param")
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], "sine")
        self.assertTrue(callable(result[1]))
        self.assertEqual(result[2], "uniform")

    def test_valid_list_input(self):
        """Test threeD_spacing_vectorLike_return_tuple with list input."""
        result = pv.threeD_spacing_vectorLike_return_tuple(
            ["sine", "uniform", "sine"], "test_param"
        )
        self.assertEqual(result, ("sine", "uniform", "sine"))

    def test_valid_ndarray_input(self):
        """Test threeD_spacing_vectorLike_return_tuple with ndarray input."""
        result = pv.threeD_spacing_vectorLike_return_tuple(
            np.array(["sine", "uniform", "sine"]), "test_param"
        )
        self.assertEqual(result, ("sine", "uniform", "sine"))

    def test_invalid_string_value(self):
        """Test threeD_spacing_vectorLike_return_tuple raises ValueError with invalid string."""
        invalid_spacing = pvf.make_invalid_spacing_string_fixture()
        with self.assertRaises(ValueError) as context:
            pv.threeD_spacing_vectorLike_return_tuple(invalid_spacing, "test_param")
        self.assertIn("Element 1", str(context.exception))
        self.assertIn("sine", str(context.exception))
        self.assertIn("uniform", str(context.exception))

    def test_invalid_type(self):
        """Test threeD_spacing_vectorLike_return_tuple raises TypeError with invalid type."""
        invalid_spacing = pvf.make_invalid_spacing_type_fixture()
        with self.assertRaises(TypeError) as context:
            pv.threeD_spacing_vectorLike_return_tuple(invalid_spacing, "test_param")
        self.assertIn("Element 1", str(context.exception))

    def test_invalid_wrong_length(self):
        """Test threeD_spacing_vectorLike_return_tuple raises ValueError with wrong length."""
        with self.assertRaises(ValueError) as context:
            pv.threeD_spacing_vectorLike_return_tuple(("sine", "uniform"), "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("3-element", str(context.exception))

    def test_invalid_string_input(self):
        """Test threeD_spacing_vectorLike_return_tuple raises TypeError with string input."""
        with self.assertRaises(TypeError) as context:
            pv.threeD_spacing_vectorLike_return_tuple("sine", "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("array-like", str(context.exception).lower())


class TestNDNumberVectorLikeReturnFloat(unittest.TestCase):
    """A class with functions to test nD_number_vectorLike_return_float."""

    def test_valid_nd_array(self):
        """Test nD_number_vectorLike_return_float with valid N-D array."""
        valid_vector = pvf.make_valid_nd_vector_fixture()
        result = pv.nD_number_vectorLike_return_float(valid_vector, "test_param")
        npt.assert_array_equal(result, valid_vector)
        self.assertEqual(result.dtype, float)

    def test_valid_list(self):
        """Test nD_number_vectorLike_return_float with valid list."""
        result = pv.nD_number_vectorLike_return_float([1, 2, 3, 4, 5], "test_param")
        npt.assert_array_equal(result, np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=float))

    def test_valid_single_element(self):
        """Test nD_number_vectorLike_return_float with single element."""
        result = pv.nD_number_vectorLike_return_float([5.0], "test_param")
        npt.assert_array_equal(result, np.array([5.0], dtype=float))
        self.assertEqual(result.shape, (1,))

    def test_invalid_2d_array(self):
        """Test nD_number_vectorLike_return_float raises ValueError with 2D array."""
        with self.assertRaises(ValueError) as context:
            pv.nD_number_vectorLike_return_float([[1.0, 2.0], [3.0, 4.0]], "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("N-element", str(context.exception))

    def test_invalid_scalar(self):
        """Test nD_number_vectorLike_return_float raises ValueError with scalar."""
        with self.assertRaises(ValueError):
            pv.nD_number_vectorLike_return_float(5.0, "test_param")

    def test_invalid_with_nan(self):
        """Test nD_number_vectorLike_return_float raises ValueError with NaN."""
        with self.assertRaises(ValueError):
            pv.nD_number_vectorLike_return_float([1.0, np.nan, 3.0], "test_param")

    def test_invalid_string(self):
        """Test nD_number_vectorLike_return_float raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.nD_number_vectorLike_return_float("invalid", "test_param")


class TestFourByFourNumberArrayLikeReturnFloat(unittest.TestCase):
    """A class with functions to test fourByFour_number_arrayLike_return_float."""

    def test_valid_4x4_array(self):
        """Test fourByFour_number_arrayLike_return_float with valid 4x4 array."""
        valid_matrix = pvf.make_valid_4x4_matrix_fixture()
        result = pv.fourByFour_number_arrayLike_return_float(valid_matrix, "test_param")
        npt.assert_array_equal(result, valid_matrix)
        self.assertEqual(result.dtype, float)
        self.assertEqual(result.shape, (4, 4))

    def test_valid_4x4_list(self):
        """Test fourByFour_number_arrayLike_return_float with valid 4x4 nested list."""
        valid_list = pvf.make_valid_4x4_list_fixture()
        result = pv.fourByFour_number_arrayLike_return_float(valid_list, "test_param")
        npt.assert_array_equal(result, np.array(valid_list, dtype=float))
        self.assertEqual(result.shape, (4, 4))

    def test_valid_with_ints(self):
        """Test fourByFour_number_arrayLike_return_float converts ints to floats."""
        int_matrix = np.eye(4, dtype=int)
        result = pv.fourByFour_number_arrayLike_return_float(int_matrix, "test_param")
        npt.assert_array_equal(result, np.eye(4, dtype=float))
        self.assertEqual(result.dtype, float)

    def test_invalid_3x3_matrix(self):
        """Test fourByFour_number_arrayLike_return_float raises ValueError with 3x3 matrix."""
        invalid_matrix = pvf.make_invalid_3x3_matrix_fixture()
        with self.assertRaises(ValueError) as context:
            pv.fourByFour_number_arrayLike_return_float(invalid_matrix, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("4x4", str(context.exception))

    def test_invalid_4x3_matrix(self):
        """Test fourByFour_number_arrayLike_return_float raises ValueError with 4x3 matrix."""
        with self.assertRaises(ValueError):
            pv.fourByFour_number_arrayLike_return_float(
                np.zeros((4, 3), dtype=float), "test_param"
            )

    def test_invalid_with_nan(self):
        """Test fourByFour_number_arrayLike_return_float raises ValueError with NaN."""
        matrix_with_nan = np.eye(4, dtype=float)
        matrix_with_nan[0, 0] = np.nan
        with self.assertRaises(ValueError):
            pv.fourByFour_number_arrayLike_return_float(matrix_with_nan, "test_param")

    def test_invalid_string(self):
        """Test fourByFour_number_arrayLike_return_float raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.fourByFour_number_arrayLike_return_float("invalid", "test_param")


class TestNonEmptyListReturnList(unittest.TestCase):
    """A class with functions to test non_empty_list_return_list."""

    def test_valid_non_empty_list(self):
        """Test non_empty_list_return_list with a valid non-empty list."""
        valid_list = pvf.make_valid_non_empty_list_fixture()
        result = pv.non_empty_list_return_list(valid_list, "test_param")
        self.assertEqual(result, valid_list)
        self.assertIsInstance(result, list)

    def test_valid_single_element_list(self):
        """Test non_empty_list_return_list with a single-element list."""
        valid_list = pvf.make_valid_single_element_list_fixture()
        result = pv.non_empty_list_return_list(valid_list, "test_param")
        self.assertEqual(result, valid_list)

    def test_valid_mixed_types_list(self):
        """Test non_empty_list_return_list with a list of mixed types."""
        mixed_list = [1, "two", 3.0, None]
        result = pv.non_empty_list_return_list(mixed_list, "test_param")
        self.assertEqual(result, mixed_list)

    def test_invalid_empty_list(self):
        """Test non_empty_list_return_list raises ValueError with empty list."""
        empty_list = pvf.make_invalid_empty_list_fixture()
        with self.assertRaises(ValueError) as context:
            pv.non_empty_list_return_list(empty_list, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("at least one element", str(context.exception))

    def test_invalid_tuple(self):
        """Test non_empty_list_return_list raises TypeError with tuple."""
        with self.assertRaises(TypeError) as context:
            pv.non_empty_list_return_list((1, 2, 3), "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("list", str(context.exception))

    def test_invalid_string(self):
        """Test non_empty_list_return_list raises TypeError with string."""
        with self.assertRaises(TypeError):
            pv.non_empty_list_return_list("string", "test_param")

    def test_invalid_ndarray(self):
        """Test non_empty_list_return_list raises TypeError with ndarray."""
        with self.assertRaises(TypeError):
            pv.non_empty_list_return_list(np.array([1, 2, 3]), "test_param")


class TestRotationOrderReturnStr(unittest.TestCase):
    """A class with functions to test rotation_order_return_str."""

    def test_valid_xyz(self):
        """Test rotation_order_return_str with valid 'xyz' order."""
        valid_order = pvf.make_valid_rotation_order_xyz_fixture()
        result = pv.rotation_order_return_str(valid_order, "test_param")
        self.assertEqual(result, "xyz")
        self.assertIsInstance(result, str)

    def test_valid_zyx(self):
        """Test rotation_order_return_str with valid 'zyx' order."""
        valid_order = pvf.make_valid_rotation_order_zyx_fixture()
        result = pv.rotation_order_return_str(valid_order, "test_param")
        self.assertEqual(result, "zyx")

    def test_valid_all_orders(self):
        """Test rotation_order_return_str with all valid orders."""
        valid_orders = ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]
        for order in valid_orders:
            with self.subTest(order=order):
                result = pv.rotation_order_return_str(order, "test_param")
                self.assertEqual(result, order)

    def test_invalid_order(self):
        """Test rotation_order_return_str raises ValueError with invalid order."""
        invalid_order = pvf.make_invalid_rotation_order_fixture()
        with self.assertRaises(ValueError) as context:
            pv.rotation_order_return_str(invalid_order, "test_param")
        self.assertIn("test_param", str(context.exception))

    def test_invalid_wrong_length(self):
        """Test rotation_order_return_str raises ValueError with wrong length."""
        with self.assertRaises(ValueError) as context:
            pv.rotation_order_return_str("xy", "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("3 characters", str(context.exception))

    def test_invalid_uppercase(self):
        """Test rotation_order_return_str raises ValueError with uppercase."""
        with self.assertRaises(ValueError):
            pv.rotation_order_return_str("XYZ", "test_param")

    def test_invalid_int(self):
        """Test rotation_order_return_str raises TypeError with int."""
        with self.assertRaises(TypeError) as context:
            pv.rotation_order_return_str(123, "test_param")
        self.assertIn("test_param", str(context.exception))
        self.assertIn("str", str(context.exception))

    def test_invalid_none(self):
        """Test rotation_order_return_str raises TypeError with None."""
        with self.assertRaises(TypeError):
            pv.rotation_order_return_str(None, "test_param")


if __name__ == "__main__":
    unittest.main()
