"""This module contains classes to test functions in the serialization module."""

import functools
import hashlib
import inspect
import json
import tempfile
import textwrap
import unittest
import zipfile
from pathlib import Path
from typing import Any, cast
from unittest import mock

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware._mujoco_model import MuJoCoModel

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
    UnboundCallable,
    _deserialize_value,
    _ndarray_from_dict,
    _ndarray_to_dict,
    _object_from_dict,
    _object_to_dict,
    _serialize_value,
    hash_object,
    load,
    save,
)
from pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method import (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver,
)
from pterasoftware.free_flight_unsteady_ring_vortex_lattice_method import (
    FreeFlightUnsteadyRingVortexLatticeMethodSolver,
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
from pterasoftware.movements.aeroelastic_wing_cross_section_movement import (
    AeroelasticWingCrossSectionMovement,
)
from pterasoftware.movements.aeroelastic_wing_movement import AeroelasticWingMovement
from pterasoftware.movements.airplane_movement import AirplaneMovement
from pterasoftware.movements.free_flight_movement import FreeFlightMovement
from pterasoftware.movements.free_flight_operating_point_movement import (
    FreeFlightOperatingPointMovement,
)
from pterasoftware.movements.movement import Movement
from pterasoftware.movements.operating_point_movement import OperatingPointMovement
from pterasoftware.movements.wing_cross_section_movement import (
    WingCrossSectionMovement,
)
from pterasoftware.movements.wing_movement import WingMovement
from pterasoftware.operating_point import OperatingPoint
from pterasoftware.problems import (
    AeroelasticUnsteadyProblem,
    FreeFlightUnsteadyProblem,
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
    mujoco_model_fixtures,
    problem_fixtures,
    serialization_fixtures,
    solver_fixtures,
)


def rewrite_archive(path: Path, replacements: dict[str, bytes | None]) -> None:
    """Rewrites a saved archive in place with some of its members replaced, dropped, or
    added.

    :param path: The path of the archive to rewrite.
    :param replacements: A dict mapping member names to their new contents, or to None
        to drop the member. Members that are not named keep their contents, the
        archive's member order is preserved, and a named member that the archive does
        not hold is appended after the existing members.
    :return: None
    """
    with zipfile.ZipFile(path) as archive:
        members = [(info.filename, archive.read(info)) for info in archive.infolist()]
    existing_names = {name for name, _ in members}
    for name, content in replacements.items():
        if name not in existing_names and content is not None:
            members.append((name, content))
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in members:
            if name in replacements:
                replacement = replacements[name]
                if replacement is None:
                    continue
                content = replacement
            archive.writestr(name, content)


def read_member(path: Path, name: str) -> dict[str, Any]:
    """Reads one JSON member of a saved archive.

    :param path: The path of the archive to read.
    :param name: The name of the member to read.
    :return: The dict parsed from the member.
    """
    with zipfile.ZipFile(path) as archive:
        member = json.loads(archive.read(name))
    assert isinstance(member, dict)
    return member


def read_header(path: Path) -> dict[str, Any]:
    """Reads the header member of a saved archive.

    :param path: The path of the archive to read.
    :return: The dict parsed from the archive's header member.
    """
    return read_member(path, "header.json")


def count_records(data: object, type_name: str) -> int:
    """Counts the full records of one class in parsed member data, ignoring refs.

    :param data: The parsed JSON of one archive member.
    :param type_name: The class name to count records of.
    :return: The number of dicts tagged with the class name that carry an "_id" key.
    """
    count = 0
    pending: list[object] = [data]
    while pending:
        node = pending.pop()
        if isinstance(node, dict):
            if node.get("_type") == type_name and "_id" in node:
                count += 1
            pending.extend(node.values())
        elif isinstance(node, list):
            pending.extend(node)
    return count


class TestNdarrayRoundTrip(unittest.TestCase):
    """This class contains methods for testing _ndarray_to_dict and _ndarray_from_dict
    round trips."""

    def test_float64_1d(self) -> None:
        """Tests round trip for a 1D float64 array.

        :return: None
        """
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.dtype, np.float64)

    def test_float64_2d(self) -> None:
        """Tests round trip for a 2D float64 array.

        :return: None
        """
        arr = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.shape, (3, 2))

    def test_float64_3d(self) -> None:
        """Tests round trip for a 3D float64 array.

        :return: None
        """
        arr = np.arange(24.0, dtype=np.float64).reshape((2, 4, 3))
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.shape, (2, 4, 3))

    def test_int64_1d(self) -> None:
        """Tests round trip for a 1D int64 array.

        :return: None
        """
        arr = np.array([10, 20, 30], dtype=np.int64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.dtype, np.int64)

    def test_bool_1d(self) -> None:
        """Tests round trip for a 1D bool array.

        :return: None
        """
        arr = np.array([True, False, True, False], dtype=np.bool_)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.dtype, np.bool_)

    def test_empty_float64(self) -> None:
        """Tests round trip for an empty float64 array.

        :return: None
        """
        arr = np.array([], dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertEqual(result.shape, (0,))
        self.assertEqual(result.dtype, np.float64)

    def test_empty_2d(self) -> None:
        """Tests round trip for an empty 2D array with a non trivial second dimension.

        :return: None
        """
        arr = np.empty((0, 3), dtype=np.float64)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertEqual(result.shape, (0, 3))
        self.assertEqual(result.dtype, np.float64)

    def test_writeable_preserved(self) -> None:
        """Tests that a writable array remains writable after round trip.

        :return: None
        """
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        self.assertTrue(arr.flags.writeable)
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertTrue(result.flags.writeable)

    def test_read_only_preserved(self) -> None:
        """Tests that a read only array remains read only after round trip.

        :return: None
        """
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        arr.flags.writeable = False
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        npt.assert_array_equal(result, arr)
        self.assertFalse(result.flags.writeable)

    def test_missing_writeable_defaults_to_writeable(self) -> None:
        """Tests that a missing writeable field defaults to a writable array.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        del serialized_dict["writeable"]
        result = _ndarray_from_dict(serialized_dict)
        self.assertTrue(result.flags.writeable)

    def test_dtype_object_with_none(self) -> None:
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

    def test_dtype_object_2d_shape(self) -> None:
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

    def test_dtype_object_read_only_preserved(self) -> None:
        """Tests that a read only dtype=object array remains read only after round trip.

        :return: None
        """
        arr = np.empty(2, dtype=object)
        arr[0] = None
        arr[1] = None
        arr.flags.writeable = False
        result = _ndarray_from_dict(_ndarray_to_dict(arr))
        self.assertFalse(result.flags.writeable)

    def test_dtype_object_writeable_preserved(self) -> None:
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

    def test_numeric_dict_keys(self) -> None:
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

    def test_object_dict_keys(self) -> None:
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

    def test_base64_data_is_string(self) -> None:
        """Tests that the base64 encoded data is a string.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        self.assertIsInstance(serialized_dict["data"], str)


class TestSerializeValue(unittest.TestCase):
    """This class contains methods for testing _serialize_value."""

    def test_none(self) -> None:
        """Tests that None serializes to None.

        :return: None
        """
        self.assertIsNone(_serialize_value(None))

    def test_bool_true(self) -> None:
        """Tests that True serializes to True.

        :return: None
        """
        result = _serialize_value(True)
        self.assertIs(result, True)

    def test_bool_false(self) -> None:
        """Tests that False serializes to False.

        :return: None
        """
        result = _serialize_value(False)
        self.assertIs(result, False)

    def test_np_bool(self) -> None:
        """Tests that a numpy bool serializes to a Python bool.

        :return: None
        """
        result = _serialize_value(np.bool_(True))
        self.assertIs(result, True)
        self.assertIsInstance(result, bool)

    def test_bool_not_wrapped_as_int(self) -> None:
        """Tests that bool values are not wrapped as int dicts.

        :return: None
        """
        result = _serialize_value(True)
        self.assertNotIsInstance(result, dict)

    def test_int(self) -> None:
        """Tests that a Python int serializes to a wrapped dict.

        :return: None
        """
        result = _serialize_value(42)
        self.assertEqual(result, {"_type": "int", "value": 42})

    def test_np_int64(self) -> None:
        """Tests that a numpy int64 serializes to a wrapped dict with a Python int.

        :return: None
        """
        result = _serialize_value(np.int64(7))
        self.assertEqual(result, {"_type": "int", "value": 7})
        assert isinstance(result, dict)
        self.assertIsInstance(result["value"], int)

    def test_float(self) -> None:
        """Tests that a Python float serializes to a wrapped dict.

        :return: None
        """
        result = _serialize_value(3.14)
        self.assertEqual(result, {"_type": "float", "value": 3.14})

    def test_np_float64(self) -> None:
        """Tests that a numpy float64 serializes to a wrapped dict with a Python float.

        :return: None
        """
        result = _serialize_value(np.float64(2.718))
        self.assertEqual(result, {"_type": "float", "value": 2.718})
        assert isinstance(result, dict)
        self.assertIsInstance(result["value"], float)

    def test_float_inf_raises(self) -> None:
        """Tests that serializing inf raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(float("inf"))

    def test_float_negative_inf_raises(self) -> None:
        """Tests that serializing negative inf raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(float("-inf"))

    def test_float_nan_raises(self) -> None:
        """Tests that serializing NaN raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _serialize_value(float("nan"))

    def test_str(self) -> None:
        """Tests that a string serializes to itself.

        :return: None
        """
        result = _serialize_value("hello")
        self.assertEqual(result, "hello")

    def test_bytes(self) -> None:
        """Tests that bytes serialize to a dict with a base64 string.

        :return: None
        """
        result = _serialize_value(b"\x00\x01binary")
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "bytes")
        self.assertIsInstance(result["data"], str)

    def test_ndarray(self) -> None:
        """Tests that a numpy array delegates to _ndarray_to_dict.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        result = _serialize_value(arr)
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "ndarray")

    def test_tuple(self) -> None:
        """Tests that a tuple serializes to a dict with items.

        :return: None
        """
        result = _serialize_value((1, 2.0, "three"))
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "tuple")
        self.assertEqual(len(result["items"]), 3)

    def test_tuple_empty(self) -> None:
        """Tests that an empty tuple serializes correctly.

        :return: None
        """
        result = _serialize_value(())
        self.assertEqual(result, {"_type": "tuple", "items": []})

    def test_list(self) -> None:
        """Tests that a list serializes to a dict with items.

        :return: None
        """
        result = _serialize_value([1, 2.0, "three"])
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "list")
        self.assertEqual(len(result["items"]), 3)

    def test_list_empty(self) -> None:
        """Tests that an empty list serializes correctly.

        :return: None
        """
        result = _serialize_value([])
        self.assertEqual(result, {"_type": "list", "items": []})

    def test_dict(self) -> None:
        """Tests that a dict serializes to a dict with items.

        :return: None
        """
        result = _serialize_value({"a": 1, "b": 2.0, "c": "three"})
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "dict")
        self.assertEqual(len(result["items"]), 3)

    def test_dict_empty(self) -> None:
        """Tests that an empty dict serializes correctly.

        :return: None
        """
        result = _serialize_value({})
        self.assertEqual(result, {"_type": "dict", "items": {}})

    def test_dict_non_str_key_raises(self) -> None:
        """Tests that a dict with a non-str key raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _serialize_value({1: "one"})

    def test_nested_tuple(self) -> None:
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

    def test_callable_sine(self) -> None:
        """Tests that the oscillating_sin_at_time function serializes by name.

        :return: None
        """
        result = _serialize_value(oscillating_sin_at_time)
        self.assertEqual(result, {"_type": "callable", "name": "sine"})

    def test_callable_uniform(self) -> None:
        """Tests that the oscillating_lin_at_time function serializes by name.

        :return: None
        """
        result = _serialize_value(oscillating_lin_at_time)
        self.assertEqual(result, {"_type": "callable", "name": "uniform"})

    def test_custom_callable_serializes_as_marker(self) -> None:
        """Tests that a custom callable serializes as an inert marker holding its
        qualified name, its dedented source text, and the SHA-256 hash of that text.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        source = textwrap.dedent(inspect.getsource(custom_spacing))
        result = _serialize_value(custom_spacing)
        self.assertEqual(
            result,
            {
                "_type": "custom_callable",
                "qualname": (
                    "tests.unit.fixtures.serialization_fixtures."
                    "make_custom_spacing_fixture.<locals>.custom_spacing"
                ),
                "source": source,
                "source_hash": hashlib.sha256(source.encode("utf-8")).hexdigest(),
            },
        )

    def test_custom_callable_lambda_serializes_as_marker(self) -> None:
        """Tests that a lambda serializes as a marker whose source is the enclosing
        statement.

        :return: None
        """
        result = _serialize_value(lambda time: 0.0)
        assert isinstance(result, dict)
        self.assertEqual(result["_type"], "custom_callable")
        self.assertTrue(result["qualname"].endswith(".<locals>.<lambda>"))
        self.assertIn("lambda time: 0.0", result["source"])
        self.assertEqual(
            result["source_hash"],
            hashlib.sha256(result["source"].encode("utf-8")).hexdigest(),
        )

    def test_custom_callable_without_source_serializes_as_marker(self) -> None:
        """Tests that a callable whose source cannot be retrieved serializes as a marker
        with null source and source hash rather than raising.

        :return: None
        """
        result = _serialize_value(len)
        self.assertEqual(
            result,
            {
                "_type": "custom_callable",
                "qualname": "builtins.len",
                "source": None,
                "source_hash": None,
            },
        )

    def test_custom_callable_partial_serializes_as_marker(self) -> None:
        """Tests that a functools.partial, which has no __qualname__ of its own,
        serializes as a marker named after its type.

        :return: None
        """
        result = _serialize_value(functools.partial(max, 0.0))
        self.assertEqual(
            result,
            {
                "_type": "custom_callable",
                "qualname": "functools.partial",
                "source": None,
                "source_hash": None,
            },
        )

    def test_unbound_callable_serializes_as_its_marker(self) -> None:
        """Tests that an UnboundCallable serializes back to the marker it was loaded
        from, without re-inspecting anything.

        :return: None
        """
        unbound = serialization_fixtures.make_unbound_callable_fixture()
        result = _serialize_value(unbound)
        self.assertEqual(
            result,
            {
                "_type": "custom_callable",
                "qualname": unbound.qualname,
                "source": unbound.source,
                "source_hash": unbound.source_hash,
            },
        )

    def test_unsupported_type_raises(self) -> None:
        """Tests that an unsupported type raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _serialize_value(set())

    def test_shared_object_serializes_as_ref(self) -> None:
        """Tests that the second encounter of one object serializes as a ref dict.

        :return: None
        """
        operating_point = OperatingPoint()
        result = _serialize_value((operating_point, operating_point))
        assert isinstance(result, dict)
        first, second = result["items"]
        self.assertEqual(first["_type"], "OperatingPoint")
        self.assertEqual(second, {"_type": "ref", "id": first["_id"]})


class TestDeserializeValue(unittest.TestCase):
    """This class contains methods for testing _deserialize_value."""

    def test_none(self) -> None:
        """Tests that None deserializes to None.

        :return: None
        """
        self.assertIsNone(_deserialize_value(None))

    def test_bool_true(self) -> None:
        """Tests that True deserializes to True.

        :return: None
        """
        self.assertIs(_deserialize_value(True), True)

    def test_bool_false(self) -> None:
        """Tests that False deserializes to False.

        :return: None
        """
        self.assertIs(_deserialize_value(False), False)

    def test_str(self) -> None:
        """Tests that a string deserializes to itself.

        :return: None
        """
        self.assertEqual(_deserialize_value("hello"), "hello")

    def test_bytes(self) -> None:
        """Tests that a bytes dict deserializes to bytes.

        :return: None
        """
        result = _deserialize_value({"_type": "bytes", "data": "AAFiaW5hcnk="})
        self.assertEqual(result, b"\x00\x01binary")
        self.assertIsInstance(result, bytes)

    def test_int(self) -> None:
        """Tests that a wrapped int dict deserializes to an int.

        :return: None
        """
        result = _deserialize_value({"_type": "int", "value": 42})
        self.assertEqual(result, 42)
        self.assertIsInstance(result, int)

    def test_float(self) -> None:
        """Tests that a wrapped float dict deserializes to a float.

        :return: None
        """
        result = _deserialize_value({"_type": "float", "value": 3.14})
        self.assertEqual(result, 3.14)
        self.assertIsInstance(result, float)

    def test_ndarray(self) -> None:
        """Tests that an ndarray dict deserializes to a numpy array.

        :return: None
        """
        arr = np.array([1.0, 2.0], dtype=np.float64)
        serialized_dict = _ndarray_to_dict(arr)
        result = _deserialize_value(serialized_dict)
        npt.assert_array_equal(result, arr)

    def test_tuple(self) -> None:
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

    def test_list(self) -> None:
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

    def test_dict(self) -> None:
        """Tests that a dict dict deserializes to a dict.

        :return: None
        """
        data = {
            "_type": "dict",
            "items": {
                "a": {"_type": "int", "value": 1},
                "b": {"_type": "float", "value": 2.0},
            },
        }
        result = _deserialize_value(data)
        self.assertEqual(result, {"a": 1, "b": 2.0})
        self.assertIsInstance(result, dict)

    def test_callable_sine(self) -> None:
        """Tests that a callable dict with name "sine" deserializes to
        oscillating_sin_at_time.

        :return: None
        """
        result = _deserialize_value({"_type": "callable", "name": "sine"})
        self.assertIs(result, oscillating_sin_at_time)

    def test_callable_uniform(self) -> None:
        """Tests that a callable dict with name "uniform" deserializes to
        oscillating_lin_at_time.

        :return: None
        """
        result = _deserialize_value({"_type": "callable", "name": "uniform"})
        self.assertIs(result, oscillating_lin_at_time)

    def test_callable_unknown_name_raises(self) -> None:
        """Tests that an unknown callable name raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value({"_type": "callable", "name": "unknown"})

    def test_custom_callable_deserializes_to_unbound_callable(self) -> None:
        """Tests that a custom callable marker deserializes to an UnboundCallable
        carrying the marker's fields.

        :return: None
        """
        source = "def my_spacing(time: float) -> float:\n    return 0.0\n"
        source_hash = hashlib.sha256(source.encode("utf-8")).hexdigest()
        result = _deserialize_value(
            {
                "_type": "custom_callable",
                "qualname": "my_module.my_spacing",
                "source": source,
                "source_hash": source_hash,
            }
        )
        assert isinstance(result, UnboundCallable)
        self.assertEqual(result.qualname, "my_module.my_spacing")
        self.assertEqual(result.source, source)
        self.assertEqual(result.source_hash, source_hash)

    def test_custom_callable_rebinds_from_callables(self) -> None:
        """Tests that a custom callable marker whose qualified name is in callables
        deserializes to the supplied function itself.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        marker = _serialize_value(custom_spacing)
        assert isinstance(marker, dict)
        result = _deserialize_value(
            marker, callables={marker["qualname"]: custom_spacing}
        )
        self.assertIs(result, custom_spacing)

    def test_custom_callable_not_in_callables_stays_unbound(self) -> None:
        """Tests that a custom callable marker whose qualified name is not in callables
        still deserializes to an UnboundCallable, and that the value level deserializer
        does not reject the unmatched key, since only the top level load can know
        whether a key went unused.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        marker = _serialize_value(custom_spacing)
        result = _deserialize_value(
            marker, callables={"some_module.other_function": custom_spacing}
        )
        self.assertIsInstance(result, UnboundCallable)

    def test_custom_callable_without_source_deserializes(self) -> None:
        """Tests that a custom callable marker with null source and source hash
        deserializes to an UnboundCallable with None for both.

        :return: None
        """
        result = _deserialize_value(
            {
                "_type": "custom_callable",
                "qualname": "builtins.len",
                "source": None,
                "source_hash": None,
            }
        )
        assert isinstance(result, UnboundCallable)
        self.assertEqual(result.qualname, "builtins.len")
        self.assertIsNone(result.source)
        self.assertIsNone(result.source_hash)

    def test_bare_int_raises(self) -> None:
        """Tests that a bare JSON int raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value(42)

    def test_bare_float_raises(self) -> None:
        """Tests that a bare JSON float raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value(3.14)

    def test_dict_without_type_raises(self) -> None:
        """Tests that a dict without a _type key raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value({"key": "value"})

    def test_unknown_type_tag_raises(self) -> None:
        """Tests that an unknown _type tag raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _deserialize_value({"_type": "unknown"})

    def test_chunked_placeholder_without_chunk_values_raises(self) -> None:
        """Tests that a chunked slot placeholder raises a ValueError when no chunk
        values are supplied, since only load holds the step members it splices from.

        :return: None
        """
        placeholder = {
            "_type": "chunked",
            "key": "list_num_wake_vortices",
            "container": "list",
            "length": 0,
        }
        with self.assertRaises(ValueError):
            _deserialize_value(placeholder)

    def test_unknown_ref_id_raises(self) -> None:
        """Tests that a ref to an id absent from the reference table raises a
        ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            _deserialize_value({"_type": "ref", "id": 0})


class TestValueRoundTrip(unittest.TestCase):
    """This class contains methods for testing _serialize_value and _deserialize_value
    round trips."""

    def test_none(self) -> None:
        """Tests round trip for None.

        :return: None
        """
        self.assertIsNone(_deserialize_value(_serialize_value(None)))

    def test_bool(self) -> None:
        """Tests round trip for bool values.

        :return: None
        """
        for value in [True, False]:
            with self.subTest(value=value):
                self.assertIs(_deserialize_value(_serialize_value(value)), value)

    def test_int(self) -> None:
        """Tests round trip for int values.

        :return: None
        """
        for value in [0, 1, -1, 42, -999]:
            with self.subTest(value=value):
                result = _deserialize_value(_serialize_value(value))
                self.assertEqual(result, value)
                self.assertIsInstance(result, int)

    def test_np_int64(self) -> None:
        """Tests round trip for numpy int64 values.

        :return: None
        """
        result = _deserialize_value(_serialize_value(np.int64(7)))
        self.assertEqual(result, 7)
        self.assertIsInstance(result, int)

    def test_float(self) -> None:
        """Tests round trip for float values.

        :return: None
        """
        for value in [0.0, 1.5, -3.14, 1e-10, 1e10]:
            with self.subTest(value=value):
                result = _deserialize_value(_serialize_value(value))
                self.assertEqual(result, value)
                self.assertIsInstance(result, float)

    def test_np_float64(self) -> None:
        """Tests round trip for numpy float64 values.

        :return: None
        """
        result = _deserialize_value(_serialize_value(np.float64(2.718)))
        self.assertEqual(result, 2.718)
        self.assertIsInstance(result, float)

    def test_str(self) -> None:
        """Tests round trip for string values.

        :return: None
        """
        for value in ["", "hello", "cosine", "uniform"]:
            with self.subTest(value=value):
                self.assertEqual(_deserialize_value(_serialize_value(value)), value)

    def test_bytes(self) -> None:
        """Tests round trip for bytes values.

        :return: None
        """
        for value in [b"", b"hello", b"\x00\x01\xff"]:
            with self.subTest(value=value):
                result = _deserialize_value(_serialize_value(value))
                self.assertEqual(result, value)
                self.assertIsInstance(result, bytes)

    def test_tuple(self) -> None:
        """Tests round trip for a tuple with mixed types.

        :return: None
        """
        value = (1, 2.0, "three", None, True)
        result = _deserialize_value(_serialize_value(value))
        self.assertEqual(result, value)
        self.assertIsInstance(result, tuple)

    def test_list(self) -> None:
        """Tests round trip for a list with mixed types.

        :return: None
        """
        value = [1, 2.0, "three", None, False]
        result = _deserialize_value(_serialize_value(value))
        self.assertEqual(result, value)
        self.assertIsInstance(result, list)

    def test_dict(self) -> None:
        """Tests round trip for a dict with mixed value types.

        :return: None
        """
        value = {"a": 1, "b": 2.0, "c": "three", "d": None, "e": b"bytes"}
        result = _deserialize_value(_serialize_value(value))
        self.assertEqual(result, value)
        self.assertIsInstance(result, dict)

    def test_nested_containers(self) -> None:
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

    def test_callable_sine(self) -> None:
        """Tests round trip for the oscillating_sin_at_time function.

        :return: None
        """
        self.assertIs(
            _deserialize_value(_serialize_value(oscillating_sin_at_time)),
            oscillating_sin_at_time,
        )

    def test_callable_uniform(self) -> None:
        """Tests round trip for the oscillating_lin_at_time function.

        :return: None
        """
        self.assertIs(
            _deserialize_value(_serialize_value(oscillating_lin_at_time)),
            oscillating_lin_at_time,
        )

    def test_custom_callable(self) -> None:
        """Tests that a custom callable round trips to an UnboundCallable whose own
        serialization is the identical marker, so a loaded object can be saved again.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        marker = _serialize_value(custom_spacing)
        result = _deserialize_value(marker)
        assert isinstance(result, UnboundCallable)
        self.assertEqual(_serialize_value(result), marker)


class TestUnboundCallable(unittest.TestCase):
    """This class contains methods for testing UnboundCallable."""

    def test_is_callable(self) -> None:
        """Tests that an UnboundCallable passes the callable checks the movement and
        problem classes dispatch on.

        :return: None
        """
        self.assertTrue(
            callable(serialization_fixtures.make_unbound_callable_fixture())
        )

    def test_properties(self) -> None:
        """Tests that the constructor arguments are exposed as properties.

        :return: None
        """
        unbound = serialization_fixtures.make_unbound_callable_fixture()
        self.assertEqual(unbound.qualname, "my_module.my_spacing")
        self.assertEqual(
            unbound.source, "def my_spacing(time: float) -> float:\n    return 0.0\n"
        )
        assert unbound.source is not None
        self.assertEqual(
            unbound.source_hash,
            hashlib.sha256(unbound.source.encode("utf-8")).hexdigest(),
        )

    def test_properties_are_read_only(self) -> None:
        """Tests that the properties cannot be reassigned.

        :return: None
        """
        unbound = serialization_fixtures.make_unbound_callable_fixture()
        for name in ("qualname", "source", "source_hash"):
            with self.subTest(name=name):
                with self.assertRaises(AttributeError):
                    setattr(unbound, name, "changed")

    def test_qualname_validation(self) -> None:
        """Tests that a non-str qualname is rejected.

        :return: None
        """
        bad_qualname: Any = 42
        with self.assertRaises(TypeError):
            UnboundCallable(qualname=bad_qualname, source=None, source_hash=None)

    def test_source_validation(self) -> None:
        """Tests that a source that is neither a str nor None is rejected.

        :return: None
        """
        bad_source: Any = 42
        with self.assertRaises(TypeError):
            UnboundCallable(qualname="my_module.f", source=bad_source, source_hash=None)

    def test_source_hash_validation(self) -> None:
        """Tests that a source hash that is neither a str nor None is rejected.

        :return: None
        """
        bad_source_hash: Any = 42
        with self.assertRaises(TypeError):
            UnboundCallable(
                qualname="my_module.f", source=None, source_hash=bad_source_hash
            )

    def test_call_raises_naming_function_and_source(self) -> None:
        """Tests that invoking an UnboundCallable raises a RuntimeError that names the
        original function and includes its stored source text.

        :return: None
        """
        unbound = serialization_fixtures.make_unbound_callable_fixture()
        assert unbound.source is not None
        with self.assertRaises(RuntimeError) as context:
            unbound(0.5)
        message = str(context.exception)
        self.assertIn(unbound.qualname, message)
        self.assertIn(unbound.source, message)

    def test_call_raises_without_source(self) -> None:
        """Tests that invoking an UnboundCallable with no stored source still raises a
        RuntimeError that names the original function.

        :return: None
        """
        unbound = serialization_fixtures.make_sourceless_unbound_callable_fixture()
        with self.assertRaises(RuntimeError) as context:
            unbound(0.5)
        self.assertIn(unbound.qualname, str(context.exception))

    def test_repr_names_function(self) -> None:
        """Tests that the repr names the original function.

        :return: None
        """
        unbound = serialization_fixtures.make_unbound_callable_fixture()
        self.assertIn(unbound.qualname, repr(unbound))


class TestObjectToDict(unittest.TestCase):
    """This class contains methods for testing _object_to_dict."""

    def test_unregistered_class_raises(self) -> None:
        """Tests that an unregistered class raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _object_to_dict("not a Ptera Software object")


class TestObjectFromDict(unittest.TestCase):
    """This class contains methods for testing _object_from_dict."""

    def test_unknown_class_raises(self) -> None:
        """Tests that an unknown class name raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            _object_from_dict({"_type": "UnknownClass"})


class TestHashObject(unittest.TestCase):
    """This class contains methods for testing hash_object."""

    def test_returns_hex_string(self) -> None:
        """Tests that hash_object returns a 64 character lowercase hex string.

        :return: None
        """
        result = hash_object(OperatingPoint())
        self.assertIsInstance(result, str)
        self.assertEqual(len(result), 64)
        self.assertTrue(all(character in "0123456789abcdef" for character in result))

    def test_deterministic_for_same_object(self) -> None:
        """Tests that hashing the same object twice returns the same digest.

        :return: None
        """
        operating_point = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
        self.assertEqual(hash_object(operating_point), hash_object(operating_point))

    def test_equal_for_distinct_but_equal_objects(self) -> None:
        """Tests that two distinct objects with equal content hash identically.

        :return: None
        """
        first = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
        second = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
        self.assertIsNot(first, second)
        self.assertEqual(hash_object(first), hash_object(second))

    def test_differs_for_different_content(self) -> None:
        """Tests that objects differing in one attribute hash differently.

        :return: None
        """
        first = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
        second = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=6.0)
        self.assertNotEqual(hash_object(first), hash_object(second))

    def test_differs_across_classes(self) -> None:
        """Tests that objects of different classes hash differently.

        :return: None
        """
        self.assertNotEqual(
            hash_object(OperatingPoint()),
            hash_object(Airfoil(name="NACA0012")),
        )

    def test_differs_for_different_aliasing(self) -> None:
        """Tests that graphs with equal content but different sharing hash differently.

        :return: None
        """
        shared_airfoil = Airfoil(name="naca2412")
        shared_wing = Wing(
            wing_cross_sections=[
                WingCrossSection(
                    airfoil=shared_airfoil,
                    num_spanwise_panels=3,
                    chord=1.0,
                    spanwise_spacing="uniform",
                ),
                WingCrossSection(
                    airfoil=shared_airfoil,
                    num_spanwise_panels=None,
                    chord=0.5,
                    Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
                ),
            ],
        )
        distinct_wing = Wing(
            wing_cross_sections=[
                WingCrossSection(
                    airfoil=Airfoil(name="naca2412"),
                    num_spanwise_panels=3,
                    chord=1.0,
                    spanwise_spacing="uniform",
                ),
                WingCrossSection(
                    airfoil=Airfoil(name="naca2412"),
                    num_spanwise_panels=None,
                    chord=0.5,
                    Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
                ),
            ],
        )
        self.assertNotEqual(hash_object(shared_wing), hash_object(distinct_wing))

    def test_stable_across_save_load_round_trip(self) -> None:
        """Tests that a save and load round trip preserves an object's digest.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        original_hash = hash_object(problem)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "problem.psz"
            save(path, problem)
            loaded = load(path)
        assert isinstance(loaded, SteadyProblem)
        self.assertEqual(hash_object(loaded), original_hash)

    def test_folds_in_format_version(self) -> None:
        """Tests that changing the format version changes the digest.

        :return: None
        """
        operating_point = OperatingPoint()
        base_hash = hash_object(operating_point)
        with mock.patch(
            "pterasoftware._serialization._FORMAT_VERSION", _FORMAT_VERSION + 1
        ):
            bumped_hash = hash_object(operating_point)
        self.assertNotEqual(base_hash, bumped_hash)

    def test_equal_for_custom_callables_with_same_source(self) -> None:
        """Tests that two distinct function objects with identical source text hash
        identically.

        :return: None
        """
        first = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=serialization_fixtures.make_custom_spacing_fixture(),
        )
        second = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=serialization_fixtures.make_custom_spacing_fixture(),
        )
        self.assertIsNot(first.spacingVCg__E, second.spacingVCg__E)
        self.assertEqual(hash_object(first), hash_object(second))

    def test_differs_for_custom_callables_with_different_source(self) -> None:
        """Tests that custom callables with different source text hash differently, so
        editing a custom spacing function invalidates cached results.

        :return: None
        """
        first = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=serialization_fixtures.make_custom_spacing_fixture(),
        )
        second = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=serialization_fixtures.make_other_custom_spacing_fixture(),
        )
        self.assertNotEqual(hash_object(first), hash_object(second))

    def test_stable_across_save_load_round_trip_with_custom_callable(self) -> None:
        """Tests that a save and load round trip preserves the digest of an object
        holding a custom callable, even though the callable loads as an UnboundCallable.

        :return: None
        """
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=serialization_fixtures.make_custom_spacing_fixture(),
        )
        original_hash = hash_object(operating_point_movement)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "operating_point_movement.psz"
            save(path, operating_point_movement)
            loaded = load(path)
        assert isinstance(loaded, OperatingPointMovement)
        self.assertEqual(hash_object(loaded), original_hash)

    def test_unregistered_object_raises(self) -> None:
        """Tests that hashing an unregistered object raises a TypeError.

        :return: None
        """
        with self.assertRaises(TypeError):
            hash_object("not a Ptera Software object")

    def test_stable_across_save_load_round_trip_for_unsteady_solver(self) -> None:
        """Tests that a save and load round trip preserves a solved unsteady solver's
        digest, even though its per step data is split across step members.

        :return: None
        """
        solver = UnsteadyRingVortexLatticeMethodSolver(
            serialization_fixtures.make_unsteady_problem_fixture()
        )
        solver.run()
        original_hash = hash_object(solver)
        self.assertEqual(hash_object(solver), original_hash)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.psz"
            save(path, solver)
            loaded = load(path)
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)
        self.assertEqual(hash_object(loaded), original_hash)

    def test_differs_for_different_step_content(self) -> None:
        """Tests that changing a value inside one time step's chunked data changes the
        digest.

        :return: None
        """
        solver = UnsteadyRingVortexLatticeMethodSolver(
            serialization_fixtures.make_unsteady_problem_fixture()
        )
        solver.run()
        original_hash = hash_object(solver)
        solver.listStackBrwrvp_GP1_CgP1[-1][0, 0] += 1.0
        self.assertNotEqual(hash_object(solver), original_hash)


class TestSaveLoad(unittest.TestCase):
    """This class contains methods for testing save and load."""

    def test_file_contains_format_version(self) -> None:
        """Tests that the saved archive's header contains the format version.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
        self.assertEqual(header["_format_version"], _FORMAT_VERSION)

    def test_file_contains_provenance(self) -> None:
        """Tests that the saved archive's header contains provenance metadata.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
        self.assertIn("_saved_at", header)
        self.assertIn("_pterasoftware_version", header)
        self.assertIn("_commit", header)
        self.assertIn("_dirty", header)

    def test_archive_layout(self) -> None:
        """Tests that a saved archive holds a header member followed by a root member,
        and that the header names the saved class and lists the root member as the only
        data member.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            with zipfile.ZipFile(path) as archive:
                names = archive.namelist()
                root = json.loads(archive.read("root.json"))
            header = read_header(path)
        self.assertEqual(names, ["header.json", "root.json"])
        self.assertEqual(header["_type"], "OperatingPoint")
        self.assertEqual(header["num_chunks"], 0)
        self.assertEqual(header["members"], ["root.json"])
        self.assertEqual(root["_type"], "OperatingPoint")

    def test_format_version_mismatch_raises(self) -> None:
        """Tests that loading a file with a mismatched format version raises, reporting
        both format versions and naming no package version.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
            header["_format_version"] = 9999
            header["_pterasoftware_version"] = "4.0.0"
            rewrite_archive(path, {"header.json": json.dumps(header).encode()})
            with self.assertRaises(ValueError) as context:
                load(path)
        message = str(context.exception)
        self.assertIn("9999", message)
        self.assertIn(str(_FORMAT_VERSION), message)
        self.assertNotIn("4.0.0", message)

    def test_non_public_header_type_raises(self) -> None:
        """Tests that loading an archive whose header names a class that is not a public
        saveable class raises a TypeError before any data is read.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
            header["_type"] = "Panel"
            rewrite_archive(path, {"header.json": json.dumps(header).encode()})
            with self.assertRaises(TypeError):
                load(path)

    def test_root_type_mismatch_raises(self) -> None:
        """Tests that loading an archive whose root member holds a different class than
        the header names raises a ValueError.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
            header["_type"] = "Airfoil"
            rewrite_archive(path, {"header.json": json.dumps(header).encode()})
            with self.assertRaises(ValueError):
                load(path)

    def test_missing_manifest_member_raises(self) -> None:
        """Tests that loading an archive that lacks a member its header lists raises a
        ValueError.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            rewrite_archive(path, {"root.json": None})
            with self.assertRaises(ValueError):
                load(path)

    def test_manifest_not_ending_in_root_raises(self) -> None:
        """Tests that loading an archive whose manifest does not end with the root
        member raises a ValueError.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
            header["members"] = ["root.json", "extra.json"]
            rewrite_archive(path, {"header.json": json.dumps(header).encode()})
            with self.assertRaises(ValueError):
                load(path)

    def test_empty_manifest_raises(self) -> None:
        """Tests that loading an archive whose manifest is empty raises a ValueError.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
            header["members"] = []
            rewrite_archive(path, {"header.json": json.dumps(header).encode()})
            with self.assertRaises(ValueError):
                load(path)

    def test_duplicate_manifest_entries_raise(self) -> None:
        """Tests that loading an archive whose manifest lists a member twice raises a
        ValueError.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            header = read_header(path)
            header["members"] = ["root.json", "root.json"]
            rewrite_archive(path, {"header.json": json.dumps(header).encode()})
            with self.assertRaises(ValueError):
                load(path)

    def test_unlisted_members_are_ignored(self) -> None:
        """Tests that an archive member the manifest does not list does not affect a
        load.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            rewrite_archive(path, {"extra.json": b"{not valid json"})
            result = load(path)
        assert isinstance(result, OperatingPoint)

    def test_non_zip_file_raises(self) -> None:
        """Tests that loading a .psz file that is not a zip archive raises a ValueError.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            path.write_bytes(b"not a zip archive")
            with self.assertRaises(ValueError):
                load(path)

    def test_save_accepts_string_path(self) -> None:
        """Tests that save accepts a string path in addition to a Path object.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / "test.psz")
            save(path, operating_point)
            result = load(path)
        assert isinstance(result, OperatingPoint)

    def test_save_invalid_extension_raises(self) -> None:
        """Tests that save raises a ValueError for an unsupported file extension,
        including the retired JSON extensions.

        :return: None
        """
        operating_point = OperatingPoint()
        for name in ["test.txt", "test.json", "test.json.gz"]:
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    save(name, operating_point)

    def test_load_invalid_extension_raises(self) -> None:
        """Tests that load raises a ValueError for an unsupported file extension,
        including the retired JSON extensions.

        :return: None
        """
        for name in ["test.dat", "test.json", "test.json.gz"]:
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    load(name)

    def test_save_directory_raises(self) -> None:
        """Tests that save raises a ValueError when the path is an existing directory.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            path.mkdir()
            with self.assertRaises(ValueError):
                save(path, operating_point)

    def test_load_directory_raises(self) -> None:
        """Tests that load raises a ValueError when the path is an existing directory.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            path.mkdir()
            with self.assertRaises(ValueError):
                load(path)

    def test_max_size_caps_decompressed_size(self) -> None:
        """Tests that the max_size parameter on load caps the decompressed size.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            with self.assertRaises(ValueError):
                load(path, max_size=10)

    def test_max_size_is_cumulative_across_members(self) -> None:
        """Tests that max_size caps the total decompressed size of every member read, so
        members that each fit under the cap still raise when their total exceeds it, and
        a cap equal to their total loads.

        :return: None
        """
        operating_point = OperatingPoint()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.psz"
            save(path, operating_point)
            with zipfile.ZipFile(path) as archive:
                header_size = len(archive.read("header.json"))
                root_size = len(archive.read("root.json"))
            self.assertLess(header_size, root_size)
            with self.assertRaises(ValueError):
                load(path, max_size=header_size + root_size - 1)
            result = load(path, max_size=header_size + root_size)
        assert isinstance(result, OperatingPoint)


class TestChunkedSaveLoad(unittest.TestCase):
    """This class contains methods for testing how save splits an unsteady solver's per
    step data across the archive's step members, and how load splices it back."""

    solver_slot_keys = [
        "_listStackBrbrvp_GP1_CgP1",
        "_listStackFrbrvp_GP1_CgP1",
        "_listStackFlbrvp_GP1_CgP1",
        "_listStackBlbrvp_GP1_CgP1",
        "list_num_wake_vortices",
        "_list_wake_vortex_strengths",
        "listStackBrwrvp_GP1_CgP1",
        "listStackFrwrvp_GP1_CgP1",
        "listStackFlwrvp_GP1_CgP1",
        "listStackBlwrvp_GP1_CgP1",
    ]
    problem_key = "unsteady_problem/_steady_problems"

    def setUp(self) -> None:
        """Solve an unsteady solver and save it to a temporary file.

        :return: None
        """
        self.solver = UnsteadyRingVortexLatticeMethodSolver(
            serialization_fixtures.make_unsteady_problem_fixture()
        )
        self.solver.run()
        self.num_steps = self.solver.num_steps
        self.step_names = [f"steps/{step:08d}.json" for step in range(self.num_steps)]
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "solver.psz"
        save(self.path, self.solver)

    def test_unsteady_solver_archive_layout(self) -> None:
        """Tests that a solved unsteady solver's archive holds a header, one step member
        per time step, and a root member, in that order, and that the header describes
        that layout.

        :return: None
        """
        with zipfile.ZipFile(self.path) as archive:
            names = archive.namelist()
        header = read_header(self.path)
        self.assertEqual(names, ["header.json", *self.step_names, "root.json"])
        self.assertEqual(header["_type"], "UnsteadyRingVortexLatticeMethodSolver")
        self.assertEqual(header["num_chunks"], self.num_steps)
        self.assertEqual(header["members"], [*self.step_names, "root.json"])

    def test_steady_solver_archive_layout(self) -> None:
        """Tests that a steady solver's archive has no step members and a manifest
        holding only the root member.

        :return: None
        """
        solver = SteadyRingVortexLatticeMethodSolver(
            serialization_fixtures.make_steady_problem_fixture()
        )
        solver.run()
        path = Path(self.tmp.name) / "steady_solver.psz"
        save(path, solver)
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
        header = read_header(path)
        self.assertEqual(names, ["header.json", "root.json"])
        self.assertEqual(header["num_chunks"], 0)
        self.assertEqual(header["members"], ["root.json"])

    def test_step_member_keys(self) -> None:
        """Tests that every step member holds exactly one entry per chunked sequence,
        keyed by the slot name for the solver's own sequences and by the qualified key
        for its UnsteadyProblem's SteadyProblems.

        :return: None
        """
        for name in self.step_names:
            with self.subTest(member=name):
                member = read_member(self.path, name)
                self.assertEqual(
                    list(member), [*self.solver_slot_keys, self.problem_key]
                )
                self.assertEqual(member[self.problem_key]["_type"], "SteadyProblem")

    def test_root_placeholders(self) -> None:
        """Tests that the root member replaces each chunked slot with a placeholder that
        names its key, its container type, and its length.

        :return: None
        """
        root = read_member(self.path, "root.json")
        for key in self.solver_slot_keys:
            with self.subTest(slot=key):
                self.assertEqual(
                    root[key],
                    {
                        "_type": "chunked",
                        "key": key,
                        "container": "list",
                        "length": self.num_steps,
                    },
                )
        self.assertEqual(
            root["unsteady_problem"]["_steady_problems"],
            {
                "_type": "chunked",
                "key": self.problem_key,
                "container": "tuple",
                "length": self.num_steps,
            },
        )

    def test_per_step_records_written_once(self) -> None:
        """Tests that every distinct Airplane is written as a full record exactly once
        across all members, so the root member's other references to the per step
        Airplanes are refs into the step members.

        :return: None
        """
        distinct_airplanes = {
            id(airplane)
            for steady_problem in self.solver.steady_problems
            for airplane in steady_problem.airplanes
        }
        for (
            airplane_movement
        ) in self.solver.unsteady_problem.movement.airplane_movements:
            distinct_airplanes.add(id(airplane_movement.base_airplane))
        total = sum(
            count_records(read_member(self.path, name), "Airplane")
            for name in [*self.step_names, "root.json"]
        )
        self.assertEqual(total, len(distinct_airplanes))
        for name in self.step_names:
            with self.subTest(member=name):
                self.assertEqual(
                    count_records(read_member(self.path, name), "Airplane"), 1
                )

    def test_containers_restored(self) -> None:
        """Tests that the UnsteadyProblem's SteadyProblems load as a tuple and the
        solver's per step sequences load as lists.

        :return: None
        """
        loaded = load(self.path)
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)
        self.assertIsInstance(loaded.unsteady_problem.steady_problems, tuple)
        for key in self.solver_slot_keys:
            with self.subTest(slot=key):
                self.assertIsInstance(getattr(loaded, key), list)
                self.assertEqual(len(getattr(loaded, key)), self.num_steps)

    def test_step_values_round_trip(self) -> None:
        """Tests that the per step wake data loads with the same values it was saved
        with.

        :return: None
        """
        loaded = load(self.path)
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)
        self.assertEqual(
            loaded.list_num_wake_vortices, self.solver.list_num_wake_vortices
        )
        for step in range(self.num_steps):
            npt.assert_array_equal(
                loaded.listStackBrwrvp_GP1_CgP1[step],
                self.solver.listStackBrwrvp_GP1_CgP1[step],
            )

    def test_solver_aliases_hold_after_load(self) -> None:
        """Tests that the solver's current Airplanes and OperatingPoint are the same
        objects as the current step's SteadyProblem's after load, even though they are
        written in different members.

        :return: None
        """
        loaded = load(self.path)
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)
        current_problem = loaded.steady_problems[loaded._current_step]
        self.assertEqual(len(loaded.current_airplanes), len(current_problem.airplanes))
        for loaded_airplane, problem_airplane in zip(
            loaded.current_airplanes, current_problem.airplanes
        ):
            self.assertIs(loaded_airplane, problem_airplane)
        self.assertIs(loaded.current_operating_point, current_problem.operating_point)
        self.assertIs(loaded.steady_problems, loaded.unsteady_problem.steady_problems)

    def test_movement_refs_resolve_into_step_members(self) -> None:
        """Tests that the Movement's per step Airplanes and OperatingPoints are the same
        objects as the SteadyProblems' after load.

        :return: None
        """
        loaded = load(self.path)
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)
        unsteady_problem = loaded.unsteady_problem
        assert isinstance(unsteady_problem, UnsteadyProblem)
        movement = unsteady_problem.movement
        for step in range(self.num_steps):
            self.assertIs(
                movement.airplanes[0][step], loaded.steady_problems[step].airplanes[0]
            )
            self.assertIs(
                movement.operating_points[step],
                loaded.steady_problems[step].operating_point,
            )

    def test_ragged_sequences_round_trip(self) -> None:
        """Tests that a coupled problem whose SteadyProblems list is shorter than the
        solver's per step lists is saved with the shorter sequence absent from the later
        step members and loads with both lengths preserved.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        num_problems = len(solver.unsteady_problem.steady_problems)
        self.assertLess(num_problems, solver.num_steps)
        path = Path(self.tmp.name) / "ragged.psz"
        save(path, solver)
        self.assertEqual(read_header(path)["num_chunks"], solver.num_steps)
        for step in range(solver.num_steps):
            with self.subTest(step=step):
                member = read_member(path, f"steps/{step:08d}.json")
                self.assertEqual(self.problem_key in member, step < num_problems)
        loaded = load(path)
        assert isinstance(loaded, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertEqual(len(loaded.unsteady_problem.steady_problems), num_problems)
        self.assertEqual(len(loaded.list_num_wake_vortices), solver.num_steps)

    def test_shared_airfoil_written_once_across_members(self) -> None:
        """Tests that an Airfoil shared by every step's WingCrossSections and by the
        base geometry is written as a full record in the first step member only, and
        loads as one object shared across every step and the base geometry.

        :return: None
        """
        solver = UnsteadyRingVortexLatticeMethodSolver(
            serialization_fixtures.make_variable_unsteady_problem_fixture()
        )
        solver.run()
        path = Path(self.tmp.name) / "variable.psz"
        save(path, solver)
        step_names = [f"steps/{step:08d}.json" for step in range(solver.num_steps)]
        for name in [*step_names, "root.json"]:
            with self.subTest(member=name):
                self.assertEqual(
                    count_records(read_member(path, name), "Airfoil"),
                    2 if name == step_names[0] else 0,
                )
        loaded = load(path)
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)
        unsteady_problem = loaded.unsteady_problem
        assert isinstance(unsteady_problem, UnsteadyProblem)
        base_wing = unsteady_problem.movement.airplane_movements[0].base_airplane.wings[
            0
        ]
        for step in range(solver.num_steps):
            wing = loaded.steady_problems[step].airplanes[0].wings[0]
            for section, base_section in zip(
                wing.wing_cross_sections, base_wing.wing_cross_sections
            ):
                self.assertIs(section.airfoil, base_section.airfoil)

    def test_reordered_manifest_raises(self) -> None:
        """Tests that a manifest listing a step member ahead of the member holding a
        record it refs raises a ValueError instead of loading.

        :return: None
        """
        solver = UnsteadyRingVortexLatticeMethodSolver(
            serialization_fixtures.make_variable_unsteady_problem_fixture()
        )
        path = Path(self.tmp.name) / "variable.psz"
        save(path, solver)
        header = read_header(path)
        members = header["members"]
        members[0], members[1] = members[1], members[0]
        rewrite_archive(path, {"header.json": json.dumps(header).encode()})
        with self.assertRaises(ValueError):
            load(path)

    def test_max_size_is_cumulative_across_step_members(self) -> None:
        """Tests that a max_size every member fits under on its own still raises when
        the members' total exceeds it, and that a max_size equal to the total loads.

        :return: None
        """
        with zipfile.ZipFile(self.path) as archive:
            sizes = [len(archive.read(name)) for name in archive.namelist()]
        with self.assertRaises(ValueError):
            load(self.path, max_size=max(sizes))
        with self.assertRaises(ValueError):
            load(self.path, max_size=sum(sizes) - 1)
        loaded = load(self.path, max_size=sum(sizes))
        assert isinstance(loaded, UnsteadyRingVortexLatticeMethodSolver)

    def test_corrupted_placeholder_length_raises(self) -> None:
        """Tests that a placeholder whose length disagrees with the number of elements
        the step members hold raises a ValueError.

        :return: None
        """
        root = read_member(self.path, "root.json")
        root["list_num_wake_vortices"]["length"] += 1
        rewrite_archive(self.path, {"root.json": json.dumps(root).encode()})
        with self.assertRaises(ValueError):
            load(self.path)

    def test_unknown_step_key_raises(self) -> None:
        """Tests that a step member key that no chunked slot claims raises a ValueError
        naming the key.

        :return: None
        """
        member = read_member(self.path, self.step_names[0])
        member["bogus"] = None
        rewrite_archive(self.path, {self.step_names[0]: json.dumps(member).encode()})
        with self.assertRaises(ValueError) as context:
            load(self.path)
        self.assertIn("bogus", str(context.exception))

    def test_missing_step_member_raises(self) -> None:
        """Tests that an archive lacking a step member its manifest lists raises a
        ValueError.

        :return: None
        """
        rewrite_archive(self.path, {self.step_names[-1]: None})
        with self.assertRaises(ValueError):
            load(self.path)


class TestLoadCallables(unittest.TestCase):
    """This class contains methods for testing load's callables argument, which rebinds
    custom callables to the functions a saved file recorded by name."""

    def setUp(self) -> None:
        """Save an OperatingPointMovement with a custom spacing function to a temporary
        file and record the function's qualified name.

        :return: None
        """
        self.custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        self.qualname = (
            self.custom_spacing.__module__ + "." + self.custom_spacing.__qualname__
        )
        self.operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=self.custom_spacing,
        )
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "operating_point_movement.psz"
        save(self.path, self.operating_point_movement)

    def test_rebinds_custom_callable(self) -> None:
        """Tests that load substitutes the supplied function for its marker, so the
        loaded movement regenerates the same motion as the original.

        :return: None
        """
        result = load(self.path, callables={self.qualname: self.custom_spacing})
        assert isinstance(result, OperatingPointMovement)
        self.assertIs(result.spacingVCg__E, self.custom_spacing)
        self.assertEqual(
            result.generate_operating_point_at_time_step(1, 0.1).vCg__E,
            self.operating_point_movement.generate_operating_point_at_time_step(
                1, 0.1
            ).vCg__E,
        )

    def test_default_leaves_callable_unbound(self) -> None:
        """Tests that load without callables still restores the marker as an
        UnboundCallable.

        :return: None
        """
        result = load(self.path)
        assert isinstance(result, OperatingPointMovement)
        self.assertIsInstance(result.spacingVCg__E, UnboundCallable)

    def test_matching_source_does_not_warn(self) -> None:
        """Tests that rebinding a function whose source hash matches the recorded one
        logs no warning.

        :return: None
        """
        with self.assertNoLogs("pterasoftware._serialization", level="WARNING"):
            load(self.path, callables={self.qualname: self.custom_spacing})

    def test_different_source_warns_and_still_rebinds(self) -> None:
        """Tests that rebinding a function whose source hash differs from the recorded
        one logs a warning naming the function, and that the rebind still happens.

        :return: None
        """
        other_spacing = serialization_fixtures.make_other_custom_spacing_fixture()
        with self.assertLogs(
            "pterasoftware._serialization", level="WARNING"
        ) as context:
            result = load(self.path, callables={self.qualname: other_spacing})
        self.assertTrue(any(self.qualname in line for line in context.output))
        assert isinstance(result, OperatingPointMovement)
        self.assertIs(result.spacingVCg__E, other_spacing)

    def test_missing_recorded_hash_does_not_warn(self) -> None:
        """Tests that rebinding a marker whose source could not be retrieved at save
        time skips the hash comparison and logs no warning.

        :return: None
        """
        # A builtin has no retrievable source. Its signature does not match a spacing
        # function's, which is fine because nothing here generates motion.
        sourceless_spacing: Any = len
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=sourceless_spacing,
        )
        path = Path(self.tmp.name) / "sourceless.psz"
        save(path, operating_point_movement)
        with self.assertNoLogs("pterasoftware._serialization", level="WARNING"):
            result = load(path, callables={"builtins.len": self.custom_spacing})
        assert isinstance(result, OperatingPointMovement)
        self.assertIs(result.spacingVCg__E, self.custom_spacing)

    def test_unhashable_rebound_function_warns(self) -> None:
        """Tests that rebinding a function whose own source cannot be retrieved, to a
        marker that did record a hash, logs a warning naming the function, since the
        recorded reference cannot be checked against it.

        :return: None
        """
        # A partial has no retrievable source. Its signature does not match a spacing
        # function's, which is fine because nothing here generates motion.
        sourceless_spacing: Any = functools.partial(max, 0.0)
        with self.assertLogs(
            "pterasoftware._serialization", level="WARNING"
        ) as context:
            result = load(self.path, callables={self.qualname: sourceless_spacing})
        self.assertTrue(any(self.qualname in line for line in context.output))
        assert isinstance(result, OperatingPointMovement)
        self.assertIs(result.spacingVCg__E, sourceless_spacing)

    def test_unused_key_raises(self) -> None:
        """Tests that a callables key matching no marker in the file raises a ValueError
        naming the key, so a typo cannot silently leave a placeholder in place.

        :return: None
        """
        with self.assertRaises(ValueError) as context:
            load(
                self.path,
                callables={
                    self.qualname: self.custom_spacing,
                    "some_module.misspelled": self.custom_spacing,
                },
            )
        self.assertIn("some_module.misspelled", str(context.exception))

    def test_non_mapping_callables_raises(self) -> None:
        """Tests that a callables argument that is not a mapping raises a TypeError.

        :return: None
        """
        bad_callables: Any = [self.custom_spacing]
        with self.assertRaises(TypeError):
            load(self.path, callables=bad_callables)

    def test_non_str_key_raises(self) -> None:
        """Tests that a callables key that is not a str raises a TypeError.

        :return: None
        """
        bad_callables: Any = {42: self.custom_spacing}
        with self.assertRaises(TypeError):
            load(self.path, callables=bad_callables)

    def test_non_callable_value_raises(self) -> None:
        """Tests that a callables value that is not callable raises a TypeError.

        :return: None
        """
        bad_callables: Any = {self.qualname: "not callable"}
        with self.assertRaises(TypeError):
            load(self.path, callables=bad_callables)


class TestAirfoilRoundTrip(unittest.TestCase):
    """This class contains methods for testing Airfoil serialization round trips."""

    def test_round_trip(self) -> None:
        """Tests that an Airfoil survives a full round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA0012")
        result = _deserialize_value(_serialize_value(airfoil))
        assert isinstance(result, Airfoil)
        self.assertEqual(result.name, "NACA0012")
        npt.assert_array_equal(result.outline_A_Lp, airfoil.outline_A_Lp)
        self.assertEqual(result.resample, airfoil.resample)
        self.assertEqual(result.n_points_per_side, airfoil.n_points_per_side)

    def test_mcl_round_trip(self) -> None:
        """Tests that the mean camber line array survives round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA2412")
        result = _deserialize_value(_serialize_value(airfoil))
        assert isinstance(result, Airfoil)
        assert result.mcl_A_Lp is not None
        assert airfoil.mcl_A_Lp is not None
        npt.assert_array_equal(result.mcl_A_Lp, airfoil.mcl_A_Lp)

    def test_writeable_flags(self) -> None:
        """Tests that deserialized Airfoil arrays preserve their writeable flags.

        :return: None
        """
        airfoil = Airfoil(name="NACA0012")
        result = _deserialize_value(_serialize_value(airfoil))
        assert isinstance(result, Airfoil)
        self.assertFalse(result.outline_A_Lp.flags.writeable)

    def test_save_load_round_trip(self) -> None:
        """Tests that an Airfoil survives a save/load round trip.

        :return: None
        """
        airfoil = Airfoil(name="NACA0012")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "airfoil.psz"
            save(path, airfoil)
            result = load(path)
        assert isinstance(result, Airfoil)
        self.assertEqual(result.name, "NACA0012")
        npt.assert_array_equal(result.outline_A_Lp, airfoil.outline_A_Lp)


class TestOperatingPointRoundTrip(unittest.TestCase):
    """This class contains methods for testing OperatingPoint serialization round
    trips."""

    def test_round_trip(self) -> None:
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

    def test_with_surface_effect(self) -> None:
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

    def test_none_surface_params(self) -> None:
        """Tests that None surface parameters remain None after round trip.

        :return: None
        """
        operating_point = OperatingPoint()
        result = _deserialize_value(_serialize_value(operating_point))
        assert isinstance(result, OperatingPoint)
        self.assertIsNone(result.surfaceNormal_E)
        self.assertIsNone(result.surfacePoint_E_Eo)

    def test_none_caches_round_trip(self) -> None:
        """Tests that uncomputed caches remain None after round trip.

        :return: None
        """
        operating_point = OperatingPoint()
        result = _deserialize_value(_serialize_value(operating_point))
        assert isinstance(result, OperatingPoint)
        self.assertIsNone(object.__getattribute__(result, "_qInf__E"))
        self.assertIsNone(object.__getattribute__(result, "_T_pas_GP1_CgP1_to_W_CgP1"))
        self.assertIsNone(object.__getattribute__(result, "_vInf_GP1__E"))

    def test_save_load_round_trip(self) -> None:
        """Tests that an OperatingPoint survives a save/load round trip.

        :return: None
        """
        operating_point = OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "operating_point.psz"
            save(path, operating_point)
            result = load(path)
        assert isinstance(result, OperatingPoint)
        self.assertEqual(result.rho, 1.225)
        self.assertEqual(result.alpha, 5.0)


class TestWingCrossSectionRoundTrip(unittest.TestCase):
    """This class contains methods for testing WingCrossSection serialization round
    trips."""

    def test_round_trip(self) -> None:
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

    def test_tip_wing_cross_section(self) -> None:
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

    def test_nested_airfoil(self) -> None:
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
        npt.assert_array_equal(result.airfoil.outline_A_Lp, airfoil.outline_A_Lp)

    def test_save_load_round_trip(self) -> None:
        """Tests that a WingCrossSection survives a save/load round trip.

        :return: None
        """
        wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=8,
            chord=1.0,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing_cross_section.psz"
            save(path, wing_cross_section)
            result = load(path)
        assert isinstance(result, WingCrossSection)
        self.assertEqual(result.chord, 1.0)
        self.assertEqual(result.num_spanwise_panels, 8)


class TestPanelRoundTrip(unittest.TestCase):
    """This class contains methods for testing Panel serialization round trips."""

    def test_round_trip(self) -> None:
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

    def test_none_mutable_attrs(self) -> None:
        """Tests that None mutable attributes remain None after round trip.

        :return: None
        """
        panel = serialization_fixtures.make_basic_panel_fixture()
        result = _deserialize_value(_serialize_value(panel))
        assert isinstance(result, Panel)
        self.assertIsNone(result.forces_GP1)

    def test_writeable_flags(self) -> None:
        """Tests that deserialized Panel arrays preserve their writeable flags.

        :return: None
        """
        panel = serialization_fixtures.make_basic_panel_fixture()
        result = _deserialize_value(_serialize_value(panel))
        assert isinstance(result, Panel)
        self.assertFalse(result.Frpp_G_Cg.flags.writeable)


class TestWingRoundTrip(unittest.TestCase):
    """This class contains methods for testing Wing serialization round trips."""

    def test_round_trip(self) -> None:
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
        self.assertEqual(result.spanwise_mesh, wing.spanwise_mesh)

    def test_exploded_wing_spanwise_mesh_round_trip(self) -> None:
        """Tests that an exploded Wing's spanwise mesh marker survives round trip.

        :return: None
        """
        root_wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="naca2412"),
            num_spanwise_panels=3,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing="uniform",
        )
        tip_wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="naca2412"),
            num_spanwise_panels=None,
            chord=0.5,
            Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing=None,
        )
        wing = Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
            explode_into_strips=True,
        )
        self.assertEqual(wing.spanwise_mesh, "exploded")
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        self.assertEqual(result.spanwise_mesh, "exploded")

    def test_edge_defined_wing_round_trip(self) -> None:
        """Tests that a from_edge_points Wing's marker, edge curves, and tip trim
        fraction survive a round trip.

        :return: None
        """
        ys = np.linspace(0.0, 1.0, 11)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.5 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        wing = Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=Airfoil(name="naca0012"),
            tip_trim_fraction=0.1,
        )
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        self.assertEqual(result.spanwise_mesh, "edge_defined")
        self.assertEqual(result.tip_trim_fraction, 0.1)
        assert result.leadingEdgePoints_Wn_Ler is not None
        assert result.trailingEdgePoints_Wn_Ler is not None
        npt.assert_array_equal(result.leadingEdgePoints_Wn_Ler, leading)
        npt.assert_array_equal(result.trailingEdgePoints_Wn_Ler, trailing)
        self.assertFalse(result.leadingEdgePoints_Wn_Ler.flags.writeable)
        self.assertFalse(result.trailingEdgePoints_Wn_Ler.flags.writeable)

    def test_panels_dtype_object_round_trip(self) -> None:
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

    def test_wing_cross_sections_tuple_round_trip(self) -> None:
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

    def test_save_load_round_trip(self) -> None:
        """Tests that a meshed Wing survives a save/load round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        wing = airplane.wings[0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing.psz"
            save(path, wing)
            result = load(path)
        assert isinstance(result, Wing)
        self.assertEqual(result.name, wing.name)
        self.assertEqual(result.num_chordwise_panels, wing.num_chordwise_panels)

    def test_shared_airfoil_identity_round_trip(self) -> None:
        """Tests that one Airfoil shared by both WingCrossSections is still one shared
        Airfoil after a round trip.

        :return: None
        """
        airfoil = Airfoil(name="naca2412")
        root_wing_cross_section = WingCrossSection(
            airfoil=airfoil,
            num_spanwise_panels=3,
            chord=1.0,
            spanwise_spacing="uniform",
        )
        tip_wing_cross_section = WingCrossSection(
            airfoil=airfoil,
            num_spanwise_panels=None,
            chord=0.5,
            Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
        )
        wing = Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
        )
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        self.assertIs(
            result.wing_cross_sections[0].airfoil,
            result.wing_cross_sections[1].airfoil,
        )

    def test_distinct_airfoils_identity_round_trip(self) -> None:
        """Tests that equal but distinct Airfoils are still distinct after a round trip.

        :return: None
        """
        root_wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="naca2412"),
            num_spanwise_panels=3,
            chord=1.0,
            spanwise_spacing="uniform",
        )
        tip_wing_cross_section = WingCrossSection(
            airfoil=Airfoil(name="naca2412"),
            num_spanwise_panels=None,
            chord=0.5,
            Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
        )
        wing = Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
        )
        result = _deserialize_value(_serialize_value(wing))
        assert isinstance(result, Wing)
        self.assertIsNot(
            result.wing_cross_sections[0].airfoil,
            result.wing_cross_sections[1].airfoil,
        )


class TestAirplaneRoundTrip(unittest.TestCase):
    """This class contains methods for testing Airplane serialization round trips."""

    def test_round_trip(self) -> None:
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

    def test_nested_wing_panels(self) -> None:
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

    def test_save_load_round_trip(self) -> None:
        """Tests that a meshed Airplane survives a save/load round trip.

        :return: None
        """
        airplane = serialization_fixtures.make_meshed_airplane_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "airplane.psz"
            save(path, airplane)
            result = load(path)
        assert isinstance(result, Airplane)
        self.assertEqual(result.name, airplane.name)
        self.assertEqual(len(result.wings), len(airplane.wings))


class TestSteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing SteadyProblem serialization round
    trips."""

    def test_round_trip(self) -> None:
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

    def test_panels_have_gp1_coordinates(self) -> None:
        """Tests that deserialized Panels retain their GP1 coordinates set by the
        SteadyProblem.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        orig_panels = problem.airplanes[0].wings[0].panels
        assert orig_panels is not None
        orig_panel = orig_panels[0, 0]
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, SteadyProblem)
        loaded_panels = result.airplanes[0].wings[0].panels
        assert loaded_panels is not None
        loaded_panel = loaded_panels[0, 0]
        assert orig_panel.Frpp_GP1_CgP1 is not None
        assert loaded_panel.Frpp_GP1_CgP1 is not None
        npt.assert_array_equal(loaded_panel.Frpp_GP1_CgP1, orig_panel.Frpp_GP1_CgP1)

    def test_save_load_round_trip(self) -> None:
        """Tests that a SteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "problem.psz"
            save(path, problem)
            result = load(path)
        assert isinstance(result, SteadyProblem)
        self.assertEqual(len(result.airplanes), 1)
        assert isinstance(result.operating_point, OperatingPoint)

    def test_formation_flight_round_trip(self) -> None:
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

    def test_formation_flight_gp1_coordinates(self) -> None:
        """Tests that both Airplanes' Panels have correct GP1 coordinates after round
        trip in a formation flight SteadyProblem.

        :return: None
        """
        problem = serialization_fixtures.make_formation_steady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, SteadyProblem)
        for i in range(2):
            orig_panels = problem.airplanes[i].wings[0].panels
            loaded_panels = result.airplanes[i].wings[0].panels
            assert orig_panels is not None
            assert loaded_panels is not None
            orig_panel = orig_panels[0, 0]
            loaded_panel = loaded_panels[0, 0]
            assert orig_panel.Frpp_GP1_CgP1 is not None
            assert loaded_panel.Frpp_GP1_CgP1 is not None
            npt.assert_array_equal(loaded_panel.Frpp_GP1_CgP1, orig_panel.Frpp_GP1_CgP1)


class TestSteadyHorseshoeSolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing SteadyHorseshoeVortexLatticeMethodSolver
    serialization round trips."""

    def test_solved_round_trip(self) -> None:
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

    def test_shared_reference_identity(self) -> None:
        """Tests that the solver's shared references point into the SteadyProblem graph
        after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        for solver_airplane, problem_airplane in zip(
            result.airplanes, result._steady_problem.airplanes, strict=True
        ):
            self.assertIs(solver_airplane, problem_airplane)
        self.assertIs(result.operating_point, result._steady_problem.operating_point)

    def test_pre_run_round_trip(self) -> None:
        """Tests that a pre run solver survives round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_save_load_round_trip(self) -> None:
        """Tests that a solved steady horseshoe solver survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.psz"
            save(path, solver)
            result = load(path)
        assert isinstance(result, SteadyHorseshoeVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestSteadyRingSolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing SteadyRingVortexLatticeMethodSolver
    serialization round trips."""

    def test_solved_round_trip(self) -> None:
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

    def test_shared_reference_identity(self) -> None:
        """Tests that the solver's shared references point into the SteadyProblem graph
        after round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        for solver_airplane, problem_airplane in zip(
            result.airplanes, result._steady_problem.airplanes, strict=True
        ):
            self.assertIs(solver_airplane, problem_airplane)
        self.assertIs(result.operating_point, result._steady_problem.operating_point)

    def test_pre_run_round_trip(self) -> None:
        """Tests that a pre run solver survives round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_save_load_round_trip(self) -> None:
        """Tests that a solved steady ring solver survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_steady_problem_fixture()
        solver = SteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.psz"
            save(path, solver)
            result = load(path)
        assert isinstance(result, SteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestMovementClassesRoundTrip(unittest.TestCase):
    """This class contains methods for testing movement class serialization round
    trips."""

    def test_operating_point_movement(self) -> None:
        """Tests that an OperatingPointMovement survives a full round trip.

        :return: None
        """
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
        )
        result = _deserialize_value(_serialize_value(operating_point_movement))
        assert isinstance(result, OperatingPointMovement)
        self.assertEqual(result.base_operating_point.vCg__E, 10.0)

    def test_wing_cross_section_movement(self) -> None:
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

    def test_wing_movement(self) -> None:
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

    def test_airplane_movement(self) -> None:
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

    def test_bare_movement(self) -> None:
        """Tests that a bare Movement (not inside UnsteadyProblem) serializes all slots
        including _airplanes and _operating_points.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        movement = problem.movement
        result = _deserialize_value(_serialize_value(movement))
        assert isinstance(result, Movement)
        self.assertEqual(len(result.airplanes), len(movement.airplanes))
        self.assertEqual(len(result.operating_points), len(movement.operating_points))

    def test_save_load_operating_point_movement(self) -> None:
        """Tests that an OperatingPointMovement survives a save/load round trip.

        :return: None
        """
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "operating_point_movement.psz"
            save(path, operating_point_movement)
            result = load(path)
        assert isinstance(result, OperatingPointMovement)
        self.assertEqual(result.base_operating_point.vCg__E, 10.0)

    def test_custom_spacing_operating_point_movement(self) -> None:
        """Tests that an OperatingPointMovement with a custom spacing function round
        trips with the function replaced by an UnboundCallable.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=custom_spacing,
        )
        result = _deserialize_value(_serialize_value(operating_point_movement))
        assert isinstance(result, OperatingPointMovement)
        assert isinstance(result.spacingVCg__E, UnboundCallable)
        self.assertEqual(
            result.spacingVCg__E.qualname,
            custom_spacing.__module__ + "." + custom_spacing.__qualname__,
        )
        self.assertEqual(result.ampVCg__E, 1.0)
        self.assertEqual(result.periodVCg__E, 1.0)

    def test_custom_spacing_operating_point_movement_raises_when_regenerated(
        self,
    ) -> None:
        """Tests that a round tripped OperatingPointMovement with a custom spacing
        function raises only when it is asked to generate motion, which invokes the
        UnboundCallable. The custom spacing validator wraps the UnboundCallable's
        RuntimeError in a ValueError, so the ValueError's message names the function and
        its chained context is the RuntimeError.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=custom_spacing,
        )
        result = _deserialize_value(_serialize_value(operating_point_movement))
        assert isinstance(result, OperatingPointMovement)
        with self.assertRaises(ValueError) as context:
            result.generate_operating_point_at_time_step(1, 0.1)
        self.assertIn(custom_spacing.__qualname__, str(context.exception))
        self.assertIsInstance(context.exception.__context__, RuntimeError)

    def test_save_load_custom_spacing_operating_point_movement(self) -> None:
        """Tests that an OperatingPointMovement with a custom spacing function survives
        a save/load round trip with the function replaced by an UnboundCallable.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        operating_point_movement = OperatingPointMovement(
            base_operating_point=OperatingPoint(),
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            spacingVCg__E=custom_spacing,
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "operating_point_movement.psz"
            save(path, operating_point_movement)
            result = load(path)
        assert isinstance(result, OperatingPointMovement)
        assert isinstance(result.spacingVCg__E, UnboundCallable)
        self.assertEqual(
            result.spacingVCg__E.qualname,
            custom_spacing.__module__ + "." + custom_spacing.__qualname__,
        )

    def test_save_load_wing_cross_section_movement(self) -> None:
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
            path = Path(tmp) / "wing_cross_section_movement.psz"
            save(path, wing_cross_section_movement)
            result = load(path)
        assert isinstance(result, WingCrossSectionMovement)
        self.assertEqual(result.base_wing_cross_section.chord, 1.0)

    def test_save_load_wing_movement(self) -> None:
        """Tests that a WingMovement survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        wing_movement = problem.movement.airplane_movements[0].wing_movements[0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wing_movement.psz"
            save(path, wing_movement)
            result = load(path)
        assert isinstance(result, WingMovement)
        self.assertEqual(result.base_wing.name, wing_movement.base_wing.name)

    def test_save_load_airplane_movement(self) -> None:
        """Tests that an AirplaneMovement survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        airplane_movement = problem.movement.airplane_movements[0]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "airplane_movement.psz"
            save(path, airplane_movement)
            result = load(path)
        assert isinstance(result, AirplaneMovement)
        self.assertEqual(
            result.base_airplane.name, airplane_movement.base_airplane.name
        )

    def test_save_load_movement(self) -> None:
        """Tests that a bare Movement survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        movement = problem.movement
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "movement.psz"
            save(path, movement)
            result = load(path)
        assert isinstance(result, Movement)
        self.assertEqual(len(result.airplanes), len(movement.airplanes))


class TestUnsteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing UnsteadyProblem serialization round
    trips."""

    def test_round_trip(self) -> None:
        """Tests that an UnsteadyProblem survives a full round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, UnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)
        self.assertEqual(len(result.steady_problems), len(problem.steady_problems))

    def test_movement_airplanes_identity(self) -> None:
        """Tests that Movement._airplanes and SteadyProblem.airplanes point to the same
        objects after round trip.

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

    def test_movement_operating_points_identity(self) -> None:
        """Tests that Movement._operating_points and SteadyProblem.operating_point point
        to the same objects after round trip.

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

    def test_save_load_round_trip(self) -> None:
        """Tests that an UnsteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "unsteady.psz"
            save(path, problem)
            result = load(path)
        assert isinstance(result, UnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)

    def test_formation_flight_round_trip(self) -> None:
        """Tests that an UnsteadyProblem with two Airplanes survives a full round trip
        with correct DAG identity.

        :return: None
        """
        problem = serialization_fixtures.make_formation_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, UnsteadyProblem)
        self.assertEqual(len(result.movement.airplane_movements), 2)

        # Verify DAG identity for both Airplanes across all time steps.
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
    serialization round trips."""

    def test_solved_round_trip(self) -> None:
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

    def test_shared_reference_identity(self) -> None:
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

    def test_movement_dag_identity(self) -> None:
        """Tests that the Movement <-> SteadyProblem DAG identity is preserved through
        the unsteady solver round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        unsteady_problem = result.unsteady_problem
        assert isinstance(unsteady_problem, UnsteadyProblem)
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

    def test_pre_run_round_trip(self) -> None:
        """Tests that a pre run unsteady solver survives round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_save_load_round_trip(self) -> None:
        """Tests that a solved unsteady solver survives a save/load round trip.

        :return: None
        """
        problem = serialization_fixtures.make_unsteady_problem_fixture()
        solver = UnsteadyRingVortexLatticeMethodSolver(problem)
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "solver.psz"
            save(path, solver)
            result = load(path)
        assert isinstance(result, UnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestAeroelasticMovementClassesRoundTrip(unittest.TestCase):
    """This class contains methods for testing aeroelastic movement class serialization
    round trips."""

    def setUp(self) -> None:
        """Build a shared AeroelasticUnsteadyProblem to source the movement graph.

        :return: None
        """
        self.problem = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )

    def test_aeroelastic_movement(self) -> None:
        """Tests that an AeroelasticMovement survives a full round trip.

        :return: None
        """
        movement = self.problem.movement
        assert isinstance(movement, AeroelasticMovement)
        result = _deserialize_value(_serialize_value(movement))
        assert isinstance(result, AeroelasticMovement)
        self.assertEqual(len(result.operating_points), len(movement.operating_points))

    def test_aeroelastic_airplane_movement(self) -> None:
        """Tests that an AeroelasticAirplaneMovement survives a full round trip.

        :return: None
        """
        airplane_movement = self.problem.movement.airplane_movements[0]
        result = _deserialize_value(_serialize_value(airplane_movement))
        assert isinstance(result, AeroelasticAirplaneMovement)
        self.assertEqual(
            result.base_airplane.name, airplane_movement.base_airplane.name
        )

    def test_aeroelastic_wing_movement(self) -> None:
        """Tests that an AeroelasticWingMovement survives a full round trip, including
        its all-None second-derivative companion.

        :return: None
        """
        wing_movement = self.problem.movement.airplane_movements[0].wing_movements[0]
        assert isinstance(wing_movement, AeroelasticWingMovement)
        result = _deserialize_value(_serialize_value(wing_movement))
        assert isinstance(result, AeroelasticWingMovement)
        self.assertEqual(result.base_wing.name, wing_movement.base_wing.name)
        self.assertEqual(
            result.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
        )

    def test_aeroelastic_wing_cross_section_movement(self) -> None:
        """Tests that an AeroelasticWingCrossSectionMovement survives a full round trip.

        :return: None
        """
        wing_cross_section_movement = (
            self.problem.movement.airplane_movements[0]
            .wing_movements[0]
            .wing_cross_section_movements[0]
        )
        result = _deserialize_value(_serialize_value(wing_cross_section_movement))
        assert isinstance(result, AeroelasticWingCrossSectionMovement)

    def test_custom_callable_spacing_round_trip(self) -> None:
        """Tests that an AeroelasticWingMovement with a custom callable angular spacing
        and its second-derivative companion round trips with both callables replaced by
        UnboundCallables and the named spacings untouched.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        custom_second_derivative = (
            serialization_fixtures.make_custom_second_derivative_fixture()
        )
        wing_movement = self.problem.movement.airplane_movements[0].wing_movements[0]
        custom_wing_movement = AeroelasticWingMovement(
            base_wing=wing_movement.base_wing,
            wing_cross_section_movements=cast(
                list[AeroelasticWingCrossSectionMovement],
                list(wing_movement.wing_cross_section_movements),
            ),
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=(custom_spacing, "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[
                custom_second_derivative,
                None,
                None,
            ],
        )
        result = _deserialize_value(_serialize_value(custom_wing_movement))
        assert isinstance(result, AeroelasticWingMovement)

        spacing = result.spacingAngles_Gs_to_Wn_ixyz[0]
        assert isinstance(spacing, UnboundCallable)
        self.assertEqual(
            spacing.qualname,
            custom_spacing.__module__ + "." + custom_spacing.__qualname__,
        )
        self.assertEqual(result.spacingAngles_Gs_to_Wn_ixyz[1:], ("sine", "sine"))

        second_derivative = result.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz[0]
        assert isinstance(second_derivative, UnboundCallable)
        self.assertEqual(
            second_derivative.qualname,
            custom_second_derivative.__module__
            + "."
            + custom_second_derivative.__qualname__,
        )
        self.assertEqual(
            result.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz[1:], (None, None)
        )

    def test_partial_rebind_leaves_other_callable_unbound(self) -> None:
        """Tests that rebinding only the custom spacing of an AeroelasticWingMovement
        restores that function and leaves its second-derivative companion as an
        UnboundCallable.

        :return: None
        """
        custom_spacing = serialization_fixtures.make_custom_spacing_fixture()
        custom_second_derivative = (
            serialization_fixtures.make_custom_second_derivative_fixture()
        )
        wing_movement = self.problem.movement.airplane_movements[0].wing_movements[0]
        custom_wing_movement = AeroelasticWingMovement(
            base_wing=wing_movement.base_wing,
            wing_cross_section_movements=cast(
                list[AeroelasticWingCrossSectionMovement],
                list(wing_movement.wing_cross_section_movements),
            ),
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=(custom_spacing, "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[
                custom_second_derivative,
                None,
                None,
            ],
        )
        qualname = custom_spacing.__module__ + "." + custom_spacing.__qualname__
        result = _deserialize_value(
            _serialize_value(custom_wing_movement),
            callables={qualname: custom_spacing},
        )
        assert isinstance(result, AeroelasticWingMovement)
        self.assertIs(result.spacingAngles_Gs_to_Wn_ixyz[0], custom_spacing)
        self.assertIsInstance(
            result.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz[0], UnboundCallable
        )


class TestAeroelasticUnsteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing AeroelasticUnsteadyProblem serialization
    round trips."""

    def test_round_trip(self) -> None:
        """Tests that an AeroelasticUnsteadyProblem survives a full round trip.

        :return: None
        """
        problem = problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, AeroelasticUnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)
        self.assertEqual(len(result.steady_problems), len(problem.steady_problems))
        self.assertEqual(result.wing_density, problem.wing_density)
        self.assertEqual(result.spring_constant_rad, problem.spring_constant_rad)
        self.assertEqual(result.damping_constant_rad, problem.damping_constant_rad)
        self.assertEqual(result.step_discards, problem.step_discards)

    def test_save_load_round_trip(self) -> None:
        """Tests that an AeroelasticUnsteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aeroelastic.psz"
            save(path, problem)
            result = load(path)
        assert isinstance(result, AeroelasticUnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)


class TestAeroelasticUnsteadySolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing
    AeroelasticUnsteadyRingVortexLatticeMethodSolver serialization round trips."""

    def test_pre_run_round_trip(self) -> None:
        """Tests that a pre run aeroelastic solver survives round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)

    def test_solved_round_trip(self) -> None:
        """Tests that a solved aeroelastic solver survives a full round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)
        self.assertEqual(result.num_steps, solver.num_steps)

    def test_shared_reference_identity(self) -> None:
        """Tests that the solver's reconstructed steady problems are the same objects as
        those reachable through its AeroelasticUnsteadyProblem after round trip.

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

    def test_per_wing_state_round_trip(self) -> None:
        """Tests that the AeroelasticUnsteadyProblem's per-wing deformation state (both
        the angle and angle derivative time series) survives the solver round trip.

        The spring-damper ODEs are re-seeded from both time series, so both must survive
        for a reconstructed problem's structural state to be usable.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        original = solver.unsteady_problem
        reconstructed = result.unsteady_problem
        assert isinstance(original, AeroelasticUnsteadyProblem)
        assert isinstance(reconstructed, AeroelasticUnsteadyProblem)
        self.assertEqual(
            len(reconstructed.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz),
            len(original.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz),
        )
        self.assertEqual(
            len(reconstructed.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]),
            len(original.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0]),
        )
        npt.assert_array_equal(
            reconstructed.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0][-1],
            original.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0][-1],
        )
        self.assertEqual(
            len(reconstructed._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz),
            len(original._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz),
        )
        self.assertEqual(
            len(reconstructed._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[0]),
            len(original._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[0]),
        )
        npt.assert_array_equal(
            reconstructed._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[0][-1],
            original._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[0][-1],
        )

    def test_save_load_round_trip(self) -> None:
        """Tests that a solved aeroelastic solver survives a save/load round trip.

        :return: None
        """
        solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()
        solver.run()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aeroelastic_solver.psz"
            save(path, solver)
            result = load(path)
        assert isinstance(result, AeroelasticUnsteadyRingVortexLatticeMethodSolver)
        self.assertTrue(result.ran)


class TestMuJoCoModelRoundTrip(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel serialization round trips."""

    def setUp(self) -> None:
        """Build a shared FreeFlightUnsteadyProblem to source a constructed MuJoCoModel.

        :return: None
        """
        self.problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )

    def test_round_trip(self) -> None:
        """Tests that a MuJoCoModel survives a round trip, preserving its serializable
        slots.

        :return: None
        """
        model = self.problem._mujoco_model
        result = _deserialize_value(_serialize_value(model))
        assert isinstance(result, MuJoCoModel)
        self.assertEqual(result.xml_str, model.xml_str)
        self.assertEqual(result.body_id, model.body_id)
        npt.assert_array_equal(result.initial_qpos, model.initial_qpos)
        npt.assert_array_equal(result.initial_qvel, model.initial_qvel)

    def test_mujoco_assets_model_round_trip(self) -> None:
        """Tests that a MuJoCoModel built with mujoco_assets survives a round trip, with
        the rebuilt engine resolving the XML string's mesh reference from the restored
        assets dict.

        :return: None
        """
        model = mujoco_model_fixtures.make_render_geometry_mujoco_model_fixture()
        result = _deserialize_value(_serialize_value(model))
        assert isinstance(result, MuJoCoModel)
        self.assertEqual(result.xml_str, model.xml_str)
        self.assertEqual(result._mujoco_assets, model._mujoco_assets)
        # Stepping is another verification that the deserialized MuJoCoModel is usable.
        result.step()

    def test_rebuilt_engine_is_functional(self) -> None:
        """Tests that a round-tripped MuJoCoModel can be queried and stepped, confirming
        its native model and data objects were rebuilt from the XML string.

        :return: None
        """
        result = _deserialize_value(_serialize_value(self.problem._mujoco_model))
        assert isinstance(result, MuJoCoModel)
        state = result.get_state()
        self.assertEqual(
            set(state.keys()),
            {
                "position_E_Eo",
                "R_pas_E_to_BP1",
                "velocity_E__E",
                "omegas_BP1__E",
                "time",
            },
        )
        # Stepping exercises the rebuilt native model and data objects.
        result.step()


class TestFreeFlightMovementClassesRoundTrip(unittest.TestCase):
    """This class contains methods for testing free flight movement class serialization
    round trips."""

    def setUp(self) -> None:
        """Build a shared FreeFlightUnsteadyProblem to source the movement graph.

        :return: None
        """
        self.problem = (
            problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        )

    def test_free_flight_movement(self) -> None:
        """Tests that a FreeFlightMovement survives a full round trip, including its
        pregenerated airplanes.

        :return: None
        """
        movement = self.problem.movement
        assert isinstance(movement, FreeFlightMovement)
        result = _deserialize_value(_serialize_value(movement))
        assert isinstance(result, FreeFlightMovement)
        self.assertEqual(result.num_steps, movement.num_steps)
        self.assertEqual(len(result.airplanes[0]), len(movement.airplanes[0]))

    def test_free_flight_operating_point_movement(self) -> None:
        """Tests that a FreeFlightOperatingPointMovement survives a full round trip.

        :return: None
        """
        operating_point_movement = self.problem.movement.operating_point_movement
        assert isinstance(operating_point_movement, FreeFlightOperatingPointMovement)
        result = _deserialize_value(_serialize_value(operating_point_movement))
        assert isinstance(result, FreeFlightOperatingPointMovement)
        self.assertEqual(
            len(result.operating_points),
            len(operating_point_movement.operating_points),
        )


class TestFreeFlightUnsteadyProblemRoundTrip(unittest.TestCase):
    """This class contains methods for testing FreeFlightUnsteadyProblem serialization
    round trips."""

    def test_round_trip(self) -> None:
        """Tests that a FreeFlightUnsteadyProblem survives a full round trip, with its
        MuJoCoModel rebuilt and functional.

        :return: None
        """
        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, FreeFlightUnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)
        self.assertEqual(len(result.steady_problems), len(problem.steady_problems))
        npt.assert_array_equal(result.I_BP1_CgP1, problem.I_BP1_CgP1)
        self.assertEqual(result.mass, problem.mass)
        self.assertEqual(result.k_max, problem.k_max)
        self.assertIsNone(result.external_loads_fn)
        self.assertIsInstance(result._mujoco_model, MuJoCoModel)
        # The rebuilt MuJoCoModel is functional.
        result._mujoco_model.step()

    def test_non_default_k_max_round_trip(self) -> None:
        """Tests that a non-default k_max survives a round trip.

        :return: None
        """
        fixture = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        problem = FreeFlightUnsteadyProblem(
            movement=fixture._free_flight_movement,
            mass=fixture.mass,
            I_BP1_CgP1=fixture.I_BP1_CgP1,
            k_max=5,
        )
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, FreeFlightUnsteadyProblem)
        self.assertEqual(result.k_max, 5)

    def test_custom_external_loads_fn_round_trip(self) -> None:
        """Tests that a FreeFlightUnsteadyProblem with a custom external_loads_fn round
        trips with the function replaced by an UnboundCallable and its MuJoCoModel
        rebuilt and functional.

        :return: None
        """

        # noinspection PyUnusedLocal
        def external_loads_fn(
            operating_point: OperatingPoint, airplane: Airplane
        ) -> tuple[np.ndarray, np.ndarray]:
            return np.zeros(3, dtype=float), np.zeros(3, dtype=float)

        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture(
            external_loads_fn=external_loads_fn
        )
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, FreeFlightUnsteadyProblem)
        assert isinstance(result.external_loads_fn, UnboundCallable)
        self.assertEqual(
            result.external_loads_fn.qualname,
            external_loads_fn.__module__ + "." + external_loads_fn.__qualname__,
        )
        # The rebuilt MuJoCoModel is functional.
        result._mujoco_model.step()

    def test_custom_external_loads_fn_rebinds(self) -> None:
        """Tests that a FreeFlightUnsteadyProblem's custom external_loads_fn is restored
        as the supplied function when its qualified name is in callables.

        :return: None
        """

        # noinspection PyUnusedLocal
        def external_loads_fn(
            operating_point: OperatingPoint, airplane: Airplane
        ) -> tuple[np.ndarray, np.ndarray]:
            return np.zeros(3, dtype=float), np.zeros(3, dtype=float)

        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture(
            external_loads_fn=external_loads_fn
        )
        qualname = external_loads_fn.__module__ + "." + external_loads_fn.__qualname__
        result = _deserialize_value(
            _serialize_value(problem), callables={qualname: external_loads_fn}
        )
        assert isinstance(result, FreeFlightUnsteadyProblem)
        self.assertIs(result.external_loads_fn, external_loads_fn)

    def test_mujoco_assets_problem_round_trip(self) -> None:
        """Tests that a FreeFlightUnsteadyProblem whose MuJoCoModel was built with
        mujoco_assets survives a full round trip, preserving the assets dict.

        :return: None
        """
        mujoco_assets = {
            "tetrahedron.stl": (
                mujoco_model_fixtures.make_tetrahedron_stl_bytes_fixture()
            )
        }
        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture(
            mujoco_assets=mujoco_assets
        )
        result = _deserialize_value(_serialize_value(problem))
        assert isinstance(result, FreeFlightUnsteadyProblem)
        self.assertEqual(result._mujoco_model._mujoco_assets, mujoco_assets)
        # The rebuilt MuJoCoModel is functional.
        result._mujoco_model.step()

    def test_save_load_round_trip(self) -> None:
        """Tests that a FreeFlightUnsteadyProblem survives a save/load round trip.

        :return: None
        """
        problem = problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "free_flight.psz"
            save(path, problem)
            result = load(path)
        assert isinstance(result, FreeFlightUnsteadyProblem)
        self.assertEqual(result.num_steps, problem.num_steps)


class TestFreeFlightUnsteadySolverRoundTrip(unittest.TestCase):
    """This class contains methods for testing
    FreeFlightUnsteadyRingVortexLatticeMethodSolver serialization round trips."""

    def test_pre_run_round_trip(self) -> None:
        """Tests that a pre run free flight solver survives a round trip.

        :return: None
        """
        solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, FreeFlightUnsteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)
        self.assertEqual(result.num_steps, solver.num_steps)

    def test_shared_reference_identity(self) -> None:
        """Tests that the solver's steady problems are the same objects as those
        reachable through its FreeFlightUnsteadyProblem after a round trip.

        :return: None
        """
        solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()
        result = _deserialize_value(_serialize_value(solver))
        assert isinstance(result, FreeFlightUnsteadyRingVortexLatticeMethodSolver)
        for reconstructed, problem_side in zip(
            result.steady_problems, result.unsteady_problem.steady_problems
        ):
            self.assertIs(reconstructed, problem_side)

    def test_save_load_round_trip(self) -> None:
        """Tests that a free flight solver survives a save/load round trip.

        :return: None
        """
        solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "free_flight_solver.psz"
            save(path, solver)
            result = load(path)
        assert isinstance(result, FreeFlightUnsteadyRingVortexLatticeMethodSolver)
        self.assertFalse(result.ran)
