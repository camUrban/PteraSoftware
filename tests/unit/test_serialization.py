"""This module contains classes to test functions in the serialization module."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware._oscillation import (
    oscillating_lin_at_time,
    oscillating_sin_at_time,
)

# noinspection PyProtectedMember
from pterasoftware._panel import Panel

# noinspection PyProtectedMember
from pterasoftware._serialization import (
    _FORMAT_VERSION,
    _deserialize_value,
    _ndarray_from_dict,
    _ndarray_to_dict,
    _object_from_dict,
    _object_to_dict,
    _serialize_value,
    load,
    save,
)
from pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method import (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver,
)

# noinspection PyProtectedMember
from pterasoftware.geometry.airfoil import Airfoil
from pterasoftware.geometry.airplane import Airplane
from pterasoftware.geometry.wing import Wing
from pterasoftware.geometry.wing_cross_section import WingCrossSection
from pterasoftware.movements.aeroelastic_airplane_movement import (
    AeroelasticAirplaneMovement,
)
from pterasoftware.movements.aeroelastic_movement import AeroelasticMovement
from pterasoftware.movements.aeroelastic_operating_point_movement import (
    AeroelasticOperatingPointMovement,
)
from pterasoftware.movements.aeroelastic_wing_cross_section_movement import (
    AeroelasticWingCrossSectionMovement,
)
from pterasoftware.movements.aeroelastic_wing_movement import AeroelasticWingMovement
from pterasoftware.movements.airplane_movement import AirplaneMovement
from pterasoftware.movements.movement import Movement
from pterasoftware.movements.operating_point_movement import OperatingPointMovement
from pterasoftware.movements.wing_cross_section_movement import (
    WingCrossSectionMovement,
)
from pterasoftware.movements.wing_movement import WingMovement
from pterasoftware.operating_point import OperatingPoint
from pterasoftware.problems import (
    AeroelasticUnsteadyProblem,
    SteadyProblem,
    UnsteadyProblem,
)
from pterasoftware.steady_horseshoe_vortex_lattice_method import (
    SteadyHorseshoeVortexLatticeMethodSolver,
)
from pterasoftware.steady_ring_vortex_lattice_method import (
    SteadyRingVortexLatticeMethodSolver,
)
from pterasoftware.unsteady_ring_vortex_lattice_method import (
    UnsteadyRingVortexLatticeMethodSolver,
)
from tests.unit.fixtures import (
    problem_fixtures,
    serialization_fixtures,
    solver_fixtures,
)


class TestNdarrayRoundTrip(unittest.TestCase):
    """This class contains methods for testing _ndarray_to_dict and
    _ndarray_from_dict round trips.
    """

    def test_float64_1d(self):
        """Tests round trip for a 1D float64 array.

        :return: None
        """
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.dtype, np.float64)

    def test_float64_2d(self):
        """Tests round trip for a 2D float64 array.

        :return: None
        """
        arr = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.shape, (3, 2))

    def test_float64_3d(self):
        """Tests round trip for a 3D float64 array.

        :return: None
        """
        arr = np.arange(24.0, dtype=np.float64).reshape((2, 4, 3))
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.shape, (2, 4, 3))

    def test_int64_1d(self):
        """Tests round trip for a 1D int64 array.

        :return: None
        """
        arr = np.array([10, 20, 30], dtype=np.int64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.dtype, np.int64)

    def test_bool_1d(self):
        """Tests round trip for a 1D bool array.

        :return: None
        """
        arr = np.array([True, False, True, False], dtype=np.bool_)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.dtype, np.bool_)

    def test_empty_float64(self):
        """Tests round trip for an empty float64 array.

        :return: None
        """
        arr = np.array([], dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.shape, (0,))
        self.assertEqual(result.dtype, np.float64)

    def test_empty_2d(self):
        """Tests round trip for an empty 2D array with a non trivial second dimension.

        :return: None
        """
        arr = np.empty((0, 3), dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertEqual(result.shape, (0, 3))
        self.assertEqual(result.dtype, np.float64)

    def test_writeable_preserved(self):
        """Tests that a writable array remains writable after round trip.

        :return: None
        """
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        self.assertTrue(arr.flags.writeable)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertTrue(result.flags.writeable)

    def test_read_only_preserved(self):
        """Tests that a read only array remains read only after round trip.

        :return: None
        """
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        arr.flags.writeable = False
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertFalse(result.flags.writeable)

    def test_missing_writeable_defaults_to_writeable(self):
        """Tests that a missing writeable field defaults to a writable array.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        del serialized_dict["writeable"]
        result = _ndarray_from_dict(serialized_dict)
        self.assertTrue(result.flags.writeable)

    def test_dtype_object_with_none(self):
        """Tests round trip for a dtype=object array containing None values.

        :return: None
        """
        arr = np.empty(3, dtype=object)
        arr[0] = None
        arr[1] = None
        arr[2] = None
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertEqual(result.dtype, object)
        self.assertEqual(result.shape, (3,))
        for i in range(3):
            self.assertIsNone(result[i])

    def test_dtype_object_2d_shape(self):
        """Tests that a 2D dtype=object array preserves its shape after round trip.

        :return: None
        """
        arr = np.empty((2, 3), dtype=object)
        for i in range(2):
            for j in range(3):
                arr[i, j] = None
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertEqual(result.shape, (2, 3))
        self.assertEqual(result.dtype, object)

    def test_dtype_object_read_only_preserved(self):
        """Tests that a read only dtype=object array remains read only after round trip.

        :return: None
        """
        arr = np.empty(2, dtype=object)
        arr[0] = None
        arr[1] = None
        arr.flags.writeable = False
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertFalse(result.flags.writeable)

    def test_dtype_object_writeable_preserved(self):
        """Tests that a writable dtype=object array remains writable after round trip.

        :return: None
        """
        arr = np.empty(2, dtype=object)
        arr[0] = None
        arr[1] = None
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertTrue(result.flags.writeable)


class TestNdarrayToDict(unittest.TestCase):
    """This class contains methods for testing _ndarray_to_dict output structure."""

    def test_numeric_dict_keys(self):
        """Tests that a numeric array produces a dict with the expected keys.

        :return: None
        """
        arr = np.array([1.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        self.assertEqual(serialized_dict["_type"], "ndarray")
        self.assertEqual(serialized_dict["dtype"], "float64")
        self.assertEqual(serialized_dict["shape"], [1])
        self.assertIn("data", serialized_dict)
        self.assertIn("writeable", serialized_dict)
        self.assertNotIn("items", serialized_dict)

    def test_object_dict_keys(self):
        """Tests that a dtype=object array produces a dict with the expected keys.

        :return: None
        """
        arr = np.empty(1, dtype=object)
        arr[0] = None
        serialized_dict = _ndarray_to_dict(arr)
        self.assertEqual(serialized_dict["_type"], "ndarray")
        self.assertEqual(serialized_dict["dtype"], "object")
        self.assertEqual(serialized_dict["shape"], [1])
        self.assertIn("items", serialized_dict)
        self.assertIn("writeable", serialized_dict)
        self.assertNotIn("data", serialized_dict)

    def test_base64_data_is_string(self):
        """Tests that the base64 encoded data is a string.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        self.assertIsInstance(serialized_dict["data"], str)


class TestSerializeValue(unittest.TestCase):
    """This class contains methods for testing _serialize_value."""

    def test_none(self):
        """Tests that None serializes to None.

        :return: None
        """
        self.assertIsNone(_serialize_value(None))

    def test_bool_true(self):
        """Tests that True serializes to True.

        :return: None
        """
        result = _serialize_value(True)
        self.assertIs(result, True)

    def test_bool_false(self):
        """Tests that False serializes to False.

        :return: None
        """
        result = _serialize_value(False)
        self.assertIs(result, False)

    def test_np_bool(self):
        """Tests that a numpy bool serializes to a Python bool.

        :return: None
        """
        result = _serialize_value(np.bool_(True))
        self.assertIs(result, True)
        self.assertIsInstance(result, bool)

    def test_bool_not_wrapped_as_int(self):
        """Tests that bool values are not wrapped as int dicts.

        :return: None
        """
        result = _serialize_value(True)
        self.assertNotIsInstance(result, dict)

    def test_int(self):
        """Tests that a Python int serializes to a wrapped dict.

        :return: None
        """
        result = _serialize_value(42)
        self.assertEqual(result, {"_type": "int", "value": 42})

    def test_np_int64(self):
        """Tests that a numpy int64 serializes to a wrapped dict with a Python int.

        :return: None
        """
        result = _serialize_value(np.int64(7))
        self.assertEqual(result, {"_type": "int", "value": 7})
        assert isinstance(result, dict)
        self.assertIsInstance(result["value"], int)

    def test_float(self):
        """Tests that a Python float serializes to a wrapped dict.

        :return: None
        """
        result = _serialize_value(3.14)
        self.assertEqual(result, {"_type": "float", "value": 3.14})

    def test_np_float64(self):
        """Tests that a numpy float64 serializes to a wrapped dict with a Python float.

        :return: None
        """
        result = _serialize_value(np.float64(2.718))
        self.assertEqual(result, {"_type": "float", "value": 2.718})
        assert isinstance(result, dict)
        self.assertIsInstance(result["value"], float)

    def test_float_inf_raises(self):
        """Tests that serializing inf raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(float("inf"))

    def test_float_negative_inf_raises(self):
        """Tests that serializing negative inf raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(float("-inf"))

    def test_float_nan_raises(self):
        """Tests that serializing NaN raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(float("nan"))

    def test_str(self):
        """Tests that a string serializes to itself.

        :return: None
        """
        result = _serialize_value("hello")
        self.assertEqual(result, "hello")

    def test_ndarray(self):
        """Tests that a numpy array delegates to _ndarray_to_dict.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        result = _serialize_value(arr)
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "ndarray")

    def test_tuple(self):
        """Tests that a tuple serializes to a dict with items.

        :return: None
        """
        result = _serialize_value((1, 2.0, "three"))
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "tuple")
        self.assertEqual(len(result["items"]), 3)

    def test_tuple_empty(self):
        """Tests that an empty tuple serializes correctly.

        :return: None
        """
        result = _serialize_value(())
        self.assertEqual(result, {"_type": "tuple", "items": []})

    def test_list(self):
        """Tests that a list serializes to a dict with items.

        :return: None
        """
        result = _serialize_value([1, 2.0, "three"])
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "list")
        self.assertEqual(len(result["items"]), 3)

    def test_list_empty(self):
        """Tests that an empty list serializes correctly.

        :return: None
        """
        result = _serialize_value([])
        self.assertEqual(result, {"_type": "list", "items": []})

    def test_nested_tuple(self):
        """Tests that a nested tuple serializes recursively.

        :return: None
        """
        result = _serialize_value((1, (2, 3)))
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "tuple")
        inner = result["items"][1]
        assert isinstance(inner, dict)
        self.assertEqual(inner["_type"], "tuple")
        self.assertEqual(len(inner["items"]), 2)

    def test_callable_sine(self):
        """Tests that the oscillating_sin_at_time function serializes by name.

        :return: None
        """
        result = _serialize_value(oscillating_sin_at_time)
        self.assertEqual(result, {"_type": "callable", "name": "sine"})

    def test_callable_uniform(self):
        """Tests that the oscillating_lin_at_time function serializes by name.

        :return: None
        """
        result = _serialize_value(oscillating_lin_at_time)
        self.assertEqual(result, {"_type": "callable", "name": "uniform"})

    def test_callable_custom_raises(self):
        """Tests that a custom callable raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(lambda x: x)

    def test_unsupported_type_raises(self):
        """Tests that an unsupported type raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _serialize_value(set())


class TestDeserializeValue(unittest.TestCase):
    """This class contains methods for testing _deserialize_value."""

    def test_none(self):
        """Tests that None deserializes to None.

        :return: None
        """
        self.assertIsNone(_deserialize_value(None))

    def test_bool_true(self):
        """Tests that True deserializes to True.

        :return: None
        """
        self.assertIs(_deserialize_value(True), True)

    def test_bool_false(self):
        """Tests that False deserializes to False.

        :return: None
        """
        self.assertIs(_deserialize_value(False), False)

    def test_str(self):
        """Tests that a string deserializes to itself.

        :return: None
        """
        self.assertEqual(_deserialize_value("hello"), "hello")

    def test_int(self):
        """Tests that a wrapped int dict deserializes to an int.

        :return: None
        """
        result = _deserialize_value({"_type": "int", "value": 42})
        self.assertEqual(result, 42)
        self.assertIsInstance(result, int)

    def test_float(self):
        """Tests that a wrapped float dict deserializes to a float.

        :return: None
        """
        result = _deserialize_value({"_type": "float", "value": 3.14})
        self.assertEqual(result, 3.14)
        self.assertIsInstance(result, float)

    def test_ndarray(self):
        """Tests that an ndarray dict deserializes to a numpy array.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        result = _deserialize_value(serialized_dict)
        npt.assert_array_equal(result, arr)

    def test_tuple(self):
        """Tests that a tuple dict deserializes to a tuple.

        :return: None
        """
        data = {
            "_type": "tuple",
            "items": [
                {"_type": "int", "value": 1},
                {"_type": "float", "value": 2.0},
                "three",
            ],
        }
        result = _deserialize_value(data)
        self.assertEqual(result, (1, 2.0, "three"))
        self.assertIsInstance(result, tuple)

    def test_list(self):
        """Tests that a list dict deserializes to a list.

        :return: None
        """
        data = {
            "_type": "list",
            "items": [
                {"_type": "int", "value": 1},
                {"_type": "float", "value": 2.0},
            ],
        }
        result = _deserialize_value(data)
        self.assertEqual(result, [1, 2.0])
        self.assertIsInstance(result, list)

    def test_callable_sine(self):
        """Tests that a callable dict with name "sine" deserializes to
        oscillating_sin_at_time.

        :return: None
        """
        result = _deserialize_value({"_type": "callable", "name": "sine"})
        self.assertIs(result, oscillating_sin_at_time)

    def test_callable_uniform(self):
        """Tests that a callable dict with name "uniform" deserializes to
        oscillating_lin_at_time.

        :return: None
        """
        result = _deserialize_value({"_type": "callable", "name": "uniform"})
        self.assertIs(result, oscillating_lin_at_time)

    def test_callable_unknown_name_raises(self):
        """Tests that an unknown callable name raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value({"_type": "callable", "name": "unknown"})

    def test_bare_int_raises(self):
        """Tests that a bare JSON int raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value(42)

    def test_bare_float_raises(self):
        """Tests that a bare JSON float raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value(3.14)

    def test_dict_without_type_raises(self):
        """Tests that a dict without a _type key raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value({"key": "value"})

    def test_unknown_type_tag_raises(self):
        """Tests that an unknown _type tag raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _deserialize_value({"_type": "unknown"})


class TestValueRoundTrip(unittest.TestCase):
    """This class contains methods for testing _serialize_value and _deserialize_value
    round trips.
    """

    def test_none(self):
        """Tests round trip for None.

        :return: None
        """
        self.assertIsNone(_deserialize_value(_serialize_value(None)))

    def test_bool(self):
        """Tests round trip for bool values.

        :return: None
        """
        for value in [True, False]:
            with self.subTest(value=value):
                self.assertIs(_deserialize_value(_serialize_value(value)), value)

    def test_int(self):
        """Tests round trip for int values.

        :return: None
        """
        for value in [0, 1, -1, 42, -999]:
            with self.subTest(value=value):
                result = _deserialize_value(_serialize_value(value))
                self.assertEqual(result, value)
                self.assertIsInstance(result, int)

    def test_np_int64(self):
        """Tests round trip for numpy int64 values.

        :return: None
        """
        result = _deserialize_value(_serialize_value(np.int64(7)))
        self.assertEqual(result, 7)
        self.assertIsInstance(result, int)

    def test_float(self):
        """Tests round trip for float values.

        :return: None
        """
        for value in [0.0, 1.5, -3.14, 1e-10, 1e10]:
            with self.subTest(value=value):
                result = _deserialize_value(_serialize_value(value))
                self.assertEqual(result, value)
                self.assertIsInstance(result, float)

    def test_np_float64(self):
        """Tests round trip for numpy float64 values.

        :return: None
        """
        result = _deserialize_value(_serialize_value(np.float64(2.718)))
        self.assertEqual(result, 2.718)
        self.assertIsInstance(result, float)

    def test_str(self):
        """Tests round trip for string values.

        :return: None
        """
        for value in ["", "hello", "cosine", "uniform"]:
            with self.subTest(value=value):
                self.assertEqual(_deserialize_value(_serialize_value(value)), value)

    def test_tuple(self):
        """Tests round trip for a tuple with mixed types.

        :return: None
        """
        value = (1, 2.0, "three", None, True)
        result = _deserialize_value(_serialize_value(value))
        self.assertEqual(result, value)
        self.assertIsInstance(result, tuple)

    def test_list(self):
        """Tests round trip for a list with mixed types.

        :return: None
        """
        value = [1, 2.0, "three", None, False]
        result = _deserialize_value(_serialize_value(value))
        self.assertEqual(result, value)
        self.assertIsInstance(result, list)

    def test_nested_containers(self):
        """Tests round trip for nested tuples and lists.

        :return: None
        """
        value = ([1, 2], (3.0, "four"), [None, True])
        result = _deserialize_value(_serialize_value(value))
        assert isinstance(result, tuple)
        self.assertEqual(result[0], [1, 2])
        self.assertIsInstance(result[0], list)
        self.assertEqual(result[1], (3.0, "four"))
        self.assertIsInstance(result[1], tuple)

    def test_callable_sine(self):
        """Tests round trip for the oscillating_sin_at_time function.

        :return: None
        """
        self.assertIs(
            _deserialize_value(_serialize_value(oscillating_sin_at_time)),
            oscillating_sin_at_time,
        )

    def test_callable_uniform(self):
        """Tests round trip for the oscillating_lin_at_time function.

        :return: None
        """
        self.assertIs(
            _deserialize_value(_serialize_value(oscillating_lin_at_time)),
            oscillating_lin_at_time,
        )


class TestObjectToDict(unittest.TestCase):
    """This class contains methods for testing _object_to_dict."""

    def test_unregistered_class_raises(self):
        """Tests that an unregistered class raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _object_to_dict("not a Ptera Software object")


class TestObjectFromDict(unittest.TestCase):
    """This class contains methods for testing _object_from_dict."""

    def test_unknown_class_raises(self):
        """Tests that an unknown class name raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _object_from_dict({"_type": "UnknownClass"})


class TestSaveLoad(unittest.TestCase):
    """This class contains methods for testing save and load."""

    def test_file_contains_format_version(self):
        """Tests that the saved file contains the format version.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.json"
            save(path, operating_point)
            with open(path) as f:
                data = json.load(f)
        self.assertEqual(data["_format_version"], _FORMAT_VERSION)

    def test_file_contains_provenance(self):
        """Tests that the saved file contains provenance metadata.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.json"
            save(path, operating_point)
            with open(path) as f:
                data = json.load(f)
        self.assertIn("_saved_at", data)
        self.assertIn("_pterasoftware_version", data)
        self.assertIn("_commit", data)
        self.assertIn("_dirty", data)

    def test_format_version_mismatch_raises(self):
        """Tests that loading a file with a mismatched format version raises.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.json"
            save(path, operating_point)
            with open(path) as f:
                data = json.load(f)
            data["_format_version"] = 9999
            with open(path, "w") as f:
                json.dump(data, f)
            with self.assertRaises(ValueError):
                load(path)

    def test_save_accepts_string_path(self):
        """Tests that save accepts a string path in addition to a Path object.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / "test.json")
            save(path, operating_point)
            result = load(path)
        assert isinstance(result, OperatingPoint)

    def test_save_invalid_extension_raises(self):
        """Tests that save raises a ValueError for an unsupported file extension.

        :return: None
        """
        operating_point = OperatingPoint()
        with self.assertRaises(ValueError):
            save("test.txt", operating_point)

    def test_load_invalid_extension_raises(self):
        """Tests that load raises a ValueError for an unsupported file extension.

        :return: None
        """
        with self.assertRaises(ValueError):
            load("test.dat")

    def test_gzip_bomb_protection(self):
        """Tests that the max_size parameter on load controls the size limit.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.json.gz"
            save(path, operating_point)
            with self.assertRaises(ValueError):
                load(path, max_size=10)


class TestAirfoilRoundTrip(unittest.TestCase):
    """This class contains methods for testing Airfoil serialization round trips."""

    def test_round_trip(self):
        """Tests that an Airfoil survives a full round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA0012")
        result = _deserialize_value(_serialize_value(airfoil))
        assert isinstance(result, Airfoil)
        self.assertEqual(result.name, "NACA0012")
        npt.assert_array_equal(result.outline_A_lp, airfoil.outline_A_lp)
        self.assertEqual(result.resample, airfoil.resample)
        self.assertEqual(result.n_points_per_side, airfoil.n_points_per_side)

    def test_mcl_round_trip(self):
        """Tests that the mean camber line array survives round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA2412")
        result = _deserialize_value(_serialize_value(airfoil))
        assert isinstance(result, Airfoil)
        assert result.mcl_A_lp is not None
        assert airfoil.mcl_A_lp is not None
        npt.assert_array_equal(result.mcl_A_lp, airfoil.mcl_A_lp)

    def test_writeable_flags(self):
        """Tests that deserialized Airfoil arrays preserve their writeable flags.

        :return: None
        """
        airfoil = Airfoil(name="NACA0012")
        result = _deserialize_value(_serialize_value(airfoil))
        assert isinstance(result, Airfoil)
        self.assertFalse(result.outline_A_lp.flags.writeable)

    def test_save_load_round_trip(self):
        """Tests that an Airfoil survives a save/load round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA0012")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "airfoil.json"
            save(path, airfoil)
            result = load(path)
        assert isinstance(result, Airfoil)
        self.assertEqual(result.name, "NACA0012")
        npt.assert_array_equal(result.outline_A_lp, airfoil.outline_A_lp)


class TestOperatingPointRoundTrip(unittest.TestCase):
    """This class contains methods for testing OperatingPoint serialization round
    trips.
    """

    def test_round_trip(self):
        """Tests that an OperatingPoint survives a full round trip.

        :return: None
        """
        operating_point = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0, beta=0.0)
        result = _deserialize_value(_serialize_value(operating_point))
        assert isinstance(result, OperatingPoint)
        self.assertEqual(result.rho, 1.225)
        self.assertEqual(result.vCg__E, 10.0)
        self.assertEqual(result.alpha, 5.0)
        self.assertEqual(result.beta, 0.0)

    def test_with_surface_effect(self):
        """Tests that an OperatingPoint with surface effect parameters survives round
        trip.

        :return: None
        """
        operating_point = OperatingPoint(
            surfaceNormal_E=(0.0, 0.0, 1.0),
            surfacePoint_E_Eo=(0.0, 0.0, -1.0),
        )
        result = _deserialize_value(_serialize_value(operating_point))
        assert isinstance(result, OperatingPoint)
        assert result.surfaceNormal_E is not None
        assert operating_point.surfaceNormal_E is not None
        npt.assert_array_equal(result.surfaceNormal_E, operating_point.surfaceNormal_E)
        assert result.surfacePoint_E_Eo is not None
        assert operating_point.surfacePoint_E_Eo is not None
        npt.assert_array_equal(
            result.surfacePoint_E_Eo, operating_point.surfacePoint_E_Eo
        )

    def test_none_surface_params(self):
        """Tests that None surface parameters remain None after round trip.

        :return: None
        """
        operating_point = OperatingPoint()
        result = _deserialize_value(_serialize_value(operating_point))
        assert isinstance(result, OperatingPoint)
        self.assertIsNone(result.surfaceNormal_E)
        self.assertIsNone(result.surfacePoint_E_Eo)

    def test_none_caches_round_trip(self):
        """Tests that uncomputed caches remain None after round trip.

        :return: None
        """
        operating_point = OperatingPoint()
        result = _deserialize_value(_serialize_value(operating_point))
        assert isinstance(result, OperatingPoint)
        self.assertIsNone(object.__getattribute__(result, "_qInf__E"))
        self.assertIsNone(object.__getattribute__(result, "_T_pas_GP1_CgP1_to_W_CgP1"))
        self.assertIsNone(object.__getattribute__(result, "_vInf_GP1__E"))

    def test_save_load_round_trip(self):
        """Tests that an OperatingPoint survives a save/load round trip.

        :return: None
        """
        operating_point = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "operating_point.json"
            save(path, operating_point)
            result = load(path)
        assert isinstance(result, OperatingPoint)
        self.assertEqual(result.rho, 1.225)
        self.assertEqual(result.alpha, 5.0)


class TestWingCrossSectionRoundTrip(unittest.TestCase):
    """This class contains methods for testing WingCrossSection serialization round
    trips.
    """

    def test_round_trip(self):
        """Tests that a WingCrossSection survives a full round trip.

        :return: None
        """
        wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=8,
            chord=1.0,
        )
        result = _deserialize_value(_serialize_value(wing_cross_section))
        assert isinstance(result, WingCrossSection)
        self.assertEqual(result.airfoil.name, "NACA0012")
        self.assertEqual(result.num_spanwise_panels, 8)
        self.assertEqual(result.chord, 1.0)

    def test_tip_wing_cross_section(self):
        """Tests that a tip WingCrossSection (num_spanwise_panels=None) survives round
        trip.

        :return: None
        """
        wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=None,
        )
        result = _deserialize_value(_serialize_value(wing_cross_section))
        assert isinstance(result, WingCrossSection)
        self.assertIsNone(result.num_spanwise_panels)

    def test_nested_airfoil(self):
        """Tests that the nested Airfoil object survives round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA2412")
        wing_cross_section = WingCrossSection(
            airfoil=airfoil,
            num_spanwise_panels=8,
        )
        result = _deserialize_value(_serialize_value(wing_cross_section))
        assert isinstance(result, WingCrossSection)
        npt.assert_array_equal(result.airfoil.outline_A_lp, airfoil.outline_A_lp)

    def test_save_load_round_trip(self):
        """Tests that a WingCrossSection survives a save/load round trip.

        :return: None
        """
        wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=8,
            chord=1.0,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing_cross_section.json"
            save(path, wing_cross_section)
            result = load(path)
        assert isinstance(result, WingCrossSection)
        self.assertEqual(result.chord, 1.0)
        self.assertEqual(result.num_spanwise_panels, 8)


class TestPanelRoundTrip(unittest.TestCase):
    """This class contains methods for testing Panel serialization round trips."""

    def test_round_trip(self):
        """Tests that a Panel survives a full round trip.

        :return: None
        """
        panel = serialization_fixtures.make_basic_panel_fixture()
        result = _deserialize_value(_serialize_value(panel))
        assert isinstance(result, Panel)
        npt.assert_array_equal(result.Frpp_G_Cg, np.array([1.0, 0.0, 0.0]))
        npt.assert_array_equal(result.Flpp_G_Cg, np.array([1.0, 1.0, 0.0]))
        self.assertTrue(result.is_leading_edge)
        self.assertFalse(result.is_trailing_edge)

    def test_none_mutable_attrs(self):
        """Tests that None mutable attributes remain None after round trip.

        :return: None
        """
        panel = serialization_fixtures.make_basic_panel_fixture()
        result = _deserialize_value(_serialize_value(panel))
        assert isinstance(result, Panel)
        self.assertIsNone(result.forces_GP1)

    def test_writeable_flags(self):
        """Tests that deserialized Panel arrays preserve their writeable flags.

        :return: None
        """
        panel = serialization_fixtures.make_basic_panel_fixture()
        result = _deserialize_value(_serialize_value(panel))
        assert isinstance(result, Panel)
        self.assertFalse(result.Frpp_G_Cg.flags.writeable)


class TestWingRoundTrip(unittest.TestCase):
    """This class contains methods for testing Wing serialization round trips."""

    def test_round_trip(self):
        """Tests that a meshed Wing survives a full round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        wing = airplane.wings[0]
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        self.assertEqual(result.name, wing.name)
        self.assertEqual(result.num_chordwise_panels, wing.num_chordwise_panels)
        self.assertEqual(result.chordwise_spacing, wing.chordwise_spacing)

    def test_panels_dtype_object_round_trip(self):
        """Tests that the dtype=object _panels array survives round trip with Panel
        objects.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        wing = airplane.wings[0]
        assert wing.panels is not None
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        assert result.panels is not None
        self.assertEqual(result.panels.shape, wing.panels.shape)
        self.assertEqual(result.panels.dtype, object)
        for i in range(result.panels.shape[0]):
            for j in range(result.panels.shape[1]):
                assert isinstance(result.panels[i, j], Panel)
                npt.assert_array_equal(
                    result.panels[i, j].Frpp_G_Cg,
                    wing.panels[i, j].Frpp_G_Cg,
                )

    def test_wing_cross_sections_tuple_round_trip(self):
        """Tests that the _wing_cross_sections tuple survives round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        wing = airplane.wings[0]
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        self.assertEqual(len(result.wing_cross_sections), len(wing.wing_cross_sections))
        for orig, loaded in zip(wing.wing_cross_sections, result.wing_cross_sections):
            assert isinstance(loaded, WingCrossSection)
            self.assertEqual(loaded.chord, orig.chord)

    def test_save_load_round_trip(self):
        """Tests that a meshed Wing survives a save/load round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        wing = airplane.wings[0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing.json"
            save(path, wing)
            result = load(path)
        assert isinstance(result, Wing)
        self.assertEqual(result.name, wing.name)
        self.assertEqual(result.num_chordwise_panels, wing.num_chordwise_panels)


class TestAirplaneRoundTrip(unittest.TestCase):
    """This class contains methods for testing Airplane serialization round trips."""

    def test_round_trip(self):
        """Tests that a meshed Airplane survives a full round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        result = _deserialize_value(_serialize_value(airplane))
        assert isinstance(result, Airplane)
        self.assertEqual(result.name, airplane.name)
        self.assertEqual(len(result.wings), len(airplane.wings))
        self.assertEqual(result.s_ref, airplane.s_ref)
        self.assertEqual(result.c_ref, airplane.c_ref)
        self.assertEqual(result.b_ref, airplane.b_ref)

    def test_nested_wing_panels(self):
        """Tests that an Airplane's nested Wing Panel arrays survive round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        result = _deserialize_value(_serialize_value(airplane))
        assert isinstance(result, Airplane)
        for orig_wing, loaded_wing in zip(airplane.wings, result.wings):
            assert isinstance(loaded_wing, Wing)
            assert loaded_wing.panels is not None
            assert orig_wing.panels is not None
            self.assertEqual(loaded_wing.panels.shape, orig_wing.panels.shape)

    def test_save_load_round_trip(self):
        """Tests that a meshed Airplane survives a save/load round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "airplane.json"
            save(path, airplane)
            result = load(path)
        assert isinstance(result, Airplane)
        self.assertEqual(result.name, airplane.name)
        self.assertEqual(len(result.wings), len(airplane.wings))


class TestSteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing SteadyProblem serialization round
    trips.
    """

    def test_round_trip(self):
        """Tests that a SteadyProblem survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, SteadyProblem)
        self.assertEqual(len(result.airplanes), 1)
        assert isinstance(result.airplanes[0], Airplane)
        assert isinstance(result.operating_point, OperatingPoint)
        self.assertEqual(result.operating_point.rho, 1.225)
        self.assertEqual(result.operating_point.alpha, 5.0)

    def test_panels_have_gp1_coordinates(self):
        """Tests that deserialized Panels retain their GP1 coordinates set by the
        SteadyProblem.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        orig_panel = problem.airplanes[0].wings[0].panels[0, 0]
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, SteadyProblem)
        loaded_panel = result.airplanes[0].wings[0].panels[0, 0]
        assert orig_panel.Frpp_GP1_CgP1 is not None
        assert loaded_panel.Frpp_GP1_CgP1 is not None
        npt.assert_array_equal(loaded_panel.Frpp_GP1_CgP1, orig_panel.Frpp_GP1_CgP1)

    def test_save_load_round_trip(self):
        """Tests that a SteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "problem.json"
            save(path, problem)
            result = load(path)
        assert isinstance(result, SteadyProblem)
        self.assertEqual(len(result.airplanes), 1)
        assert isinstance(result.operating_point, OperatingPoint)

    def test_formation_flight_round_trip(self):
        """Tests that a SteadyProblem with two Airplanes survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_formation_steady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, SteadyProblem)
        self.assertEqual(len(result.airplanes), 2)
        for i in range(2):
            assert isinstance(result.airplanes[i], Airplane)
            assert result.airplanes[i].wings[0].panels is not None

    def test_formation_flight_gp1_coordinates(self):
        """Tests that both Airplanes' Panels have correct GP1 coordinates after round
        trip in a formation flight SteadyProblem.

        :return: None
        """
        problem = serialization_fixtures.make_formation_steady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, SteadyProblem)
        for i in range(2):
            orig_panel = problem.airplanes[i].wings[0].panels[0, 0]
            loaded_panel = result.airplanes[i].wings[0].panels[0, 0]
            assert orig_panel.Frpp_GP1_CgP1 is not None
            assert loaded_panel.Frpp_GP1_CgP1 is not None
            npt.assert_array_equal(loaded_panel.Frpp_GP1_CgP1, orig_panel.Frpp_GP1_CgP1)


class TestSteadyHorseshoeSolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing SteadyHorseshoeVortexLatticeMethodSolver
    serialization round trips.
    """

    def test_solved_round_trip(self):
        """Tests that a solved steady horseshoe solver survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        self.assertTrue(result.ran)
        assert result.airplanes[0].forces_W is not None
        assert solver.airplanes[0].forces_W is not None
        npt.assert_array_equal(
            result.airplanes[0].forces_W, solver.airplanes[0].forces_W
        )

    def test_shared_reference_identity(self):
        """Tests that the solver's shared references point into the SteadyProblem
        graph after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        self.assertIs(result.airplanes, result._steady_problem.airplanes)
        self.assertIs(result.operating_point, result._steady_problem.operating_point)

    def test_pre_run_round_trip(self):
        """Tests that a pre run solver survives round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_save_load_round_trip(self):
        """Tests that a solved steady horseshoe solver survives a save/load round
        trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.json"
            save(path, solver)
            result = load(path)
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestSteadyRingSolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing SteadyRingVortexLatticeMethodSolver
    serialization round trips.
    """

    def test_solved_round_trip(self):
        """Tests that a solved steady ring solver survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)
        assert result.airplanes[0].forces_W is not None
        assert solver.airplanes[0].forces_W is not None
        npt.assert_array_equal(
            result.airplanes[0].forces_W, solver.airplanes[0].forces_W
        )

    def test_shared_reference_identity(self):
        """Tests that the solver's shared references point into the SteadyProblem
        graph after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        self.assertIs(result.airplanes, result._steady_problem.airplanes)
        self.assertIs(result.operating_point, result._steady_problem.operating_point)

    def test_pre_run_round_trip(self):
        """Tests that a pre run solver survives round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_save_load_round_trip(self):
        """Tests that a solved steady ring solver survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.json"
            save(path, solver)
            result = load(path)
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestMovementClassesRoundTrip(unittest.TestCase):
    """This class contains methods for testing movement class serialization round
    trips.
    """

    def test_operating_point_movement(self):
        """Tests that an OperatingPointMovement survives a full round trip.

        :return: None
        """
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
        )
        result = _deserialize_value(_serialize_value(operating_point_movement))
        assert isinstance(result, OperatingPointMovement)
        self.assertEqual(result.base_operating_point.vCg__E, 10.0)

    def test_wing_cross_section_movement(self):
        """Tests that a WingCrossSectionMovement survives a full round trip.

        :return: None
        """
        wing_cross_section_movement = WingCrossSectionMovement(
            base_wing_cross_section=WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=4,
                chord=1.0,
            ),
        )
        result = _deserialize_value(_serialize_value(wing_cross_section_movement))
        assert isinstance(result, WingCrossSectionMovement)
        self.assertEqual(result.base_wing_cross_section.chord, 1.0)

    def test_wing_movement(self):
        """Tests that a WingMovement survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        wing_movement = problem.movement.airplane_movements[0].wing_movements[0]
        result = _deserialize_value(_serialize_value(wing_movement))
        assert isinstance(result, WingMovement)
        self.assertEqual(result.base_wing.name, wing_movement.base_wing.name)
        self.assertEqual(
            len(result.wing_cross_section_movements),
            len(wing_movement.wing_cross_section_movements),
        )

    def test_airplane_movement(self):
        """Tests that an AirplaneMovement survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        airplane_movement = problem.movement.airplane_movements[0]
        result = _deserialize_value(_serialize_value(airplane_movement))
        assert isinstance(result, AirplaneMovement)
        self.assertEqual(
            result.base_airplane.name, airplane_movement.base_airplane.name
        )

    def test_bare_movement(self):
        """Tests that a bare Movement (not inside UnsteadyProblem) serializes all
        slots including _airplanes and _operating_points.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        movement = problem.movement
        result = _deserialize_value(_serialize_value(movement))
        assert isinstance(result, Movement)
        self.assertEqual(len(result.airplanes), len(movement.airplanes))
        self.assertEqual(len(result.operating_points), len(movement.operating_points))

    def test_save_load_operating_point_movement(self):
        """Tests that an OperatingPointMovement survives a save/load round trip.

        :return: None
        """
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "operating_point_movement.json"
            save(path, operating_point_movement)
            result = load(path)
        assert isinstance(result, OperatingPointMovement)
        self.assertEqual(result.base_operating_point.vCg__E, 10.0)

    def test_save_load_wing_cross_section_movement(self):
        """Tests that a WingCrossSectionMovement survives a save/load round trip.

        :return: None
        """
        wing_cross_section_movement = WingCrossSectionMovement(
            base_wing_cross_section=WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=4,
                chord=1.0,
            ),
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing_cross_section_movement.json"
            save(path, wing_cross_section_movement)
            result = load(path)
        assert isinstance(result, WingCrossSectionMovement)
        self.assertEqual(result.base_wing_cross_section.chord, 1.0)

    def test_save_load_wing_movement(self):
        """Tests that a WingMovement survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        wing_movement = problem.movement.airplane_movements[0].wing_movements[0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing_movement.json"
            save(path, wing_movement)
            result = load(path)
        assert isinstance(result, WingMovement)
        self.assertEqual(result.base_wing.name, wing_movement.base_wing.name)

    def test_save_load_airplane_movement(self):
        """Tests that an AirplaneMovement survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        airplane_movement = problem.movement.airplane_movements[0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "airplane_movement.json"
            save(path, airplane_movement)
            result = load(path)
        assert isinstance(result, AirplaneMovement)
        self.assertEqual(
            result.base_airplane.name, airplane_movement.base_airplane.name
        )

    def test_save_load_movement(self):
        """Tests that a bare Movement survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        movement = problem.movement
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "movement.json"
            save(path, movement)
            result = load(path)
        assert isinstance(result, Movement)
        self.assertEqual(len(result.airplanes), len(movement.airplanes))


class TestUnsteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing UnsteadyProblem serialization round
    trips.
    """

    def test_round_trip(self):
        """Tests that an UnsteadyProblem survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, UnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)
        self.assertEqual(len(result.steady_problems), len(problem.steady_problems))

    def test_movement_airplanes_identity(self):
        """Tests that Movement._airplanes and SteadyProblem.airplanes point to the
        same objects after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, UnsteadyProblem)
        for step in range(result.num_steps):
            for airplane_movement_index in range(
                len(result.movement.airplane_movements)
            ):
                self.assertIs(
                    result.movement.airplanes[airplane_movement_index][step],
                    result.steady_problems[step].airplanes[airplane_movement_index],
                )

    def test_movement_operating_points_identity(self):
        """Tests that Movement._operating_points and SteadyProblem.operating_point
        point to the same objects after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, UnsteadyProblem)
        for step in range(result.num_steps):
            self.assertIs(
                result.movement.operating_points[step],
                result.steady_problems[step].operating_point,
            )

    def test_save_load_round_trip(self):
        """Tests that an UnsteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "unsteady.json"
            save(path, problem)
            result = load(path)
        assert isinstance(result, UnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)

    def test_formation_flight_round_trip(self):
        """Tests that an UnsteadyProblem with two Airplanes survives a full round
        trip with correct DAG identity.

        :return: None
        """
        problem = serialization_fixtures.make_formation_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, UnsteadyProblem)
        self.assertEqual(len(result.movement.airplane_movements), 2)

        # Verify DAG identity for both airplanes across all time steps.
        for step in range(result.num_steps):
            for airplane_movement_index in range(2):
                self.assertIs(
                    result.movement.airplanes[airplane_movement_index][step],
                    result.steady_problems[step].airplanes[airplane_movement_index],
                )
            self.assertIs(
                result.movement.operating_points[step],
                result.steady_problems[step].operating_point,
            )


class TestUnsteadySolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing UnsteadyRingVortexLatticeMethodSolver
    serialization round trips.
    """

    def test_solved_round_trip(self):
        """Tests that a solved unsteady solver survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)
        self.assertEqual(result.num_steps, solver.num_steps)

    def test_shared_reference_identity(self):
        """Tests that the solver's shared references point into the UnsteadyProblem
        graph after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        self.assertIs(result.steady_problems, result.unsteady_problem.steady_problems)

    def test_movement_dag_identity(self):
        """Tests that the Movement <-> SteadyProblem DAG identity is preserved
        through the unsteady solver round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        unsteady_problem = result.unsteady_problem
        for step in range(unsteady_problem.num_steps):
            for airplane_movement_index in range(
                len(unsteady_problem.movement.airplane_movements)
            ):
                self.assertIs(
                    unsteady_problem.movement.airplanes[airplane_movement_index][step],
                    unsteady_problem.steady_problems[step].airplanes[
                        airplane_movement_index
                    ],
                )

    def test_pre_run_round_trip(self):
        """Tests that a pre run unsteady solver survives round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_save_load_round_trip(self):
        """Tests that a solved unsteady solver survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.json"
            save(path, solver)
            result = load(path)
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestAeroelasticMovementClassesRoundTrip(unittest.TestCase):
    """This class contains methods for testing aeroelastic movement class
    serialization round trips.
    """

    def setUp(self):
        """Build a shared AeroelasticUnsteadyProblem to source the movement graph.

        :return: None
        """
        self.problem = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )

    def test_aeroelastic_movement(self):
        """Tests that an AeroelasticMovement survives a full round trip.

        :return: None
        """
        movement = self.problem.movement
        result = _deserialize_value(_serialize_value(movement))
        assert isinstance(result, AeroelasticMovement)
        self.assertEqual(len(result.operating_points), len(movement.operating_points))

    def test_aeroelastic_airplane_movement(self):
        """Tests that an AeroelasticAirplaneMovement survives a full round trip.

        :return: None
        """
        airplane_movement = self.problem.movement.airplane_movements[0]
        result = _deserialize_value(_serialize_value(airplane_movement))
        assert isinstance(result, AeroelasticAirplaneMovement)
        self.assertEqual(
            result.base_airplane.name, airplane_movement.base_airplane.name
        )

    def test_aeroelastic_wing_movement(self):
        """Tests that an AeroelasticWingMovement survives a full round trip, including
        its all-None second-derivative companion.

        :return: None
        """
        wing_movement = self.problem.movement.airplane_movements[0].wing_movements[0]
        result = _deserialize_value(_serialize_value(wing_movement))
        assert isinstance(result, AeroelasticWingMovement)
        self.assertEqual(result.base_wing.name, wing_movement.base_wing.name)
        self.assertEqual(
            result.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
        )

    def test_aeroelastic_wing_cross_section_movement(self):
        """Tests that an AeroelasticWingCrossSectionMovement survives a full round
        trip.

        :return: None
        """
        wing_cross_section_movement = (
            self.problem.movement.airplane_movements[0]
            .wing_movements[0]
            .wing_cross_section_movements[0]
        )
        result = _deserialize_value(_serialize_value(wing_cross_section_movement))
        assert isinstance(result, AeroelasticWingCrossSectionMovement)

    def test_aeroelastic_operating_point_movement(self):
        """Tests that an AeroelasticOperatingPointMovement survives a full round trip.

        :return: None
        """
        operating_point_movement = self.problem.movement.operating_point_movement
        result = _deserialize_value(_serialize_value(operating_point_movement))
        assert isinstance(result, AeroelasticOperatingPointMovement)

    def test_custom_callable_spacing_is_not_serializable(self):
        """Tests that an AeroelasticWingMovement with a custom callable angular spacing
        raises on serialization, matching the standard movement-class behavior.

        :return: None
        """
        wing_movement = self.problem.movement.airplane_movements[0].wing_movements[0]
        custom_wing_movement = AeroelasticWingMovement(
            base_wing=wing_movement.base_wing,
            wing_cross_section_movements=list(
                wing_movement.wing_cross_section_movements
            ),
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=(lambda time: 0.0, "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[lambda time: 0.0, None, None],
        )
        with self.assertRaises(ValueError):
            _serialize_value(custom_wing_movement)


class TestAeroelasticUnsteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing AeroelasticUnsteadyProblem
    serialization round trips.
    """

    def test_round_trip(self):
        """Tests that an AeroelasticUnsteadyProblem survives a full round trip.

        :return: None
        """
        problem = problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, AeroelasticUnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)
        self.assertEqual(len(result.steady_problems), len(problem.steady_problems))
        self.assertEqual(result.wing_density, problem.wing_density)
        self.assertEqual(result.spring_constant, problem.spring_constant)
        self.assertEqual(result.damping_constant, problem.damping_constant)
        self.assertEqual(result.step_discards, problem.step_discards)

    def test_save_load_round_trip(self):
        """Tests that an AeroelasticUnsteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aeroelastic.json"
            save(path, problem)
            result = load(path)
        assert isinstance(result, AeroelasticUnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)


class TestAeroelasticUnsteadySolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing
    AeroelasticUnsteadyRingVortexLatticeMethodSolver serialization round trips.
    """

    def test_pre_run_round_trip(self):
        """Tests that a pre run aeroelastic solver survives round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_solved_round_trip(self):
        """Tests that a solved aeroelastic solver survives a full round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)
        self.assertEqual(result.num_steps, solver.num_steps)

    def test_shared_reference_identity(self):
        """Tests that the solver's reconstructed steady problems are the same objects
        as those reachable through its AeroelasticUnsteadyProblem after round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        for reconstructed, problem_side in zip(
            result.steady_problems, result.unsteady_problem.steady_problems
        ):
            self.assertIs(reconstructed, problem_side)

    def test_per_wing_state_round_trip(self):
        """Tests that the AeroelasticUnsteadyProblem's per-wing deformation state
        survives the solver round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        original = solver.unsteady_problem
        reconstructed = result.unsteady_problem
        self.assertEqual(
            len(reconstructed.net_deformation_per_wing),
            len(original.net_deformation_per_wing),
        )
        npt.assert_array_equal(
            reconstructed.net_deformation_per_wing[0],
            original.net_deformation_per_wing[0],
        )

    def test_save_load_round_trip(self):
        """Tests that a solved aeroelastic solver survives a save/load round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aeroelastic_solver.json"
            save(path, solver)
            result = load(path)
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)
