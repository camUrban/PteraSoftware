"""Contains the functions for serializing and deserializing Ptera Software objects."""

from __future__ import annotations

import base64
import gzip
import json
import subprocess
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np

from . import _logging
from ._oscillation import oscillating_lin_at_time, oscillating_sin_at_time

# This module is inherently coupled to the internals of every class in the package
# (it reads __slots__, knows class structure, and imports all classes into its
# registry), so importing from a sibling private module is acceptable here.
# noinspection PyProtectedMember
from ._panel import Panel

# noinspection PyProtectedMember
from .aeroelastic_unsteady_ring_vortex_lattice_method import (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver,
)
from .geometry.airfoil import Airfoil
from .geometry.airplane import Airplane
from .geometry.wing import Wing
from .geometry.wing_cross_section import WingCrossSection
from .movements.aeroelastic_airplane_movement import AeroelasticAirplaneMovement
from .movements.aeroelastic_movement import AeroelasticMovement
from .movements.aeroelastic_operating_point_movement import (
    AeroelasticOperatingPointMovement,
)
from .movements.aeroelastic_wing_cross_section_movement import (
    AeroelasticWingCrossSectionMovement,
)
from .movements.aeroelastic_wing_movement import AeroelasticWingMovement
from .movements.airplane_movement import AirplaneMovement
from .movements.movement import Movement
from .movements.operating_point_movement import OperatingPointMovement
from .movements.wing_cross_section_movement import WingCrossSectionMovement
from .movements.wing_movement import WingMovement
from .operating_point import OperatingPoint
from .problems import AeroelasticUnsteadyProblem, SteadyProblem, UnsteadyProblem
from .steady_horseshoe_vortex_lattice_method import (
    SteadyHorseshoeVortexLatticeMethodSolver,
)
from .steady_ring_vortex_lattice_method import SteadyRingVortexLatticeMethodSolver
from .unsteady_ring_vortex_lattice_method import (
    UnsteadyRingVortexLatticeMethodSolver,
)

_logger = _logging.get_logger("_serialization")

# Maps serializable callable names to their function objects and vice versa.
_CALLABLE_NAME_TO_FUNC = {
    "sine": oscillating_sin_at_time,
    "uniform": oscillating_lin_at_time,
}
_CALLABLE_FUNC_TO_NAME = {func: name for name, func in _CALLABLE_NAME_TO_FUNC.items()}

# Increments only when the serialization structure changes (slots added/removed/
# renamed, class registry changed, encoding strategy changed).
_FORMAT_VERSION = 4


def _all_slots(cls: type) -> list[str]:
    """Collects all __slots__ from a class and its parents via the MRO.

    Walks the method resolution order so that inherited slots (e.g., those on
    CoreMovement) are included alongside the class's own slots.

    :param cls: The class to inspect.
    :return: A list of slot names in MRO order (parent slots first).
    """
    slots: list[str] = []
    for klass in reversed(cls.__mro__):
        for slot in getattr(klass, "__slots__", ()):
            if slot not in slots:
                slots.append(slot)
    return slots


# Default maximum decompressed size in bytes when reading gzip files. Prevents gzip
# bombs from exhausting memory. Users can override this via the max_size parameter on
# load().
_DEFAULT_MAX_DECOMPRESSED_SIZE = 4_000_000_000  # 4 GB

# Maps class names to their types for deserialization dispatch.
_CLASS_REGISTRY: dict[str, type] = {
    "Airfoil": Airfoil,
    "OperatingPoint": OperatingPoint,
    "WingCrossSection": WingCrossSection,
    "Panel": Panel,
    "Wing": Wing,
    "Airplane": Airplane,
    "SteadyProblem": SteadyProblem,
    "SteadyHorseshoeVortexLatticeMethodSolver": SteadyHorseshoeVortexLatticeMethodSolver,
    "SteadyRingVortexLatticeMethodSolver": SteadyRingVortexLatticeMethodSolver,
    "Movement": Movement,
    "AirplaneMovement": AirplaneMovement,
    "WingMovement": WingMovement,
    "WingCrossSectionMovement": WingCrossSectionMovement,
    "OperatingPointMovement": OperatingPointMovement,
    "UnsteadyProblem": UnsteadyProblem,
    "UnsteadyRingVortexLatticeMethodSolver": UnsteadyRingVortexLatticeMethodSolver,
    "AeroelasticMovement": AeroelasticMovement,
    "AeroelasticAirplaneMovement": AeroelasticAirplaneMovement,
    "AeroelasticWingMovement": AeroelasticWingMovement,
    "AeroelasticWingCrossSectionMovement": AeroelasticWingCrossSectionMovement,
    "AeroelasticOperatingPointMovement": AeroelasticOperatingPointMovement,
    "AeroelasticUnsteadyProblem": AeroelasticUnsteadyProblem,
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver": AeroelasticUnsteadyRingVortexLatticeMethodSolver,
}

# Classes that can be saved and loaded as top level objects via save() and load().
# Internal classes (e.g., Panel) are excluded because they are not part of the public
# API and their structure may change without a format version bump. They are still
# serializable as nested objects within public classes.
_PUBLIC_SAVEABLE_CLASSES: frozenset[str] = frozenset(
    {
        "Airfoil",
        "OperatingPoint",
        "WingCrossSection",
        "Wing",
        "Airplane",
        "SteadyProblem",
        "SteadyHorseshoeVortexLatticeMethodSolver",
        "SteadyRingVortexLatticeMethodSolver",
        "Movement",
        "AirplaneMovement",
        "WingMovement",
        "WingCrossSectionMovement",
        "OperatingPointMovement",
        "UnsteadyProblem",
        "UnsteadyRingVortexLatticeMethodSolver",
        "AeroelasticMovement",
        "AeroelasticAirplaneMovement",
        "AeroelasticWingMovement",
        "AeroelasticWingCrossSectionMovement",
        "AeroelasticOperatingPointMovement",
        "AeroelasticUnsteadyProblem",
        "AeroelasticUnsteadyRingVortexLatticeMethodSolver",
    }
)

# Slots on steady solvers that are aliases into the SteadyProblem graph.
_STEADY_SOLVER_SKIP_SLOTS: frozenset[str] = frozenset(
    {"airplanes", "operating_point", "reynolds_numbers", "vInf_GP1__E", "panels"}
)

# Slots on the unsteady solver that are aliases into the UnsteadyProblem graph.
_UNSTEADY_SOLVER_SKIP_SLOTS: frozenset[str] = frozenset(
    {"steady_problems", "current_airplanes", "current_operating_point", "panels"}
)

# Slots on Movement that are redundant when serialized inside an UnsteadyProblem.
# The canonical data lives in the SteadyProblem objects; these are reconstructed
# on deserialization to preserve object identity.
_MOVEMENT_SKIP_SLOTS: frozenset[str] = frozenset({"_airplanes", "_operating_points"})


def save(path: str | Path, obj: object) -> None:
    """Saves a Ptera Software object to a JSON file.

    If the path ends with ".json.gz", the output is gzip compressed automatically. Using
    ".json.gz" is highly recommended for all but the smallest unmeshed geometry objects.
    Gzip compression typically reduces file sizes by 10x or more, and plain ".json"
    files use indented formatting that further increases their size.

    :param path: The file path to save to. Should end with ".json" or ".json.gz".
    :param obj: The Ptera Software object to save. Must be a public Ptera Software class
        (e.g., Airplane, SteadyProblem, or a solver). Internal classes such as Panel
        cannot be saved directly.
    :return: None
    """
    path = Path(path)
    if not path.name.endswith(".json") and not path.name.endswith(".json.gz"):
        raise ValueError(
            f"Path must end with '.json' or '.json.gz', got '{path.name}'."
        )

    class_name = type(obj).__name__
    if class_name not in _PUBLIC_SAVEABLE_CLASSES:
        raise TypeError(
            f"{class_name} is not a public saveable class. Only public Ptera Software "
            f"classes can be saved via save()."
        )

    _logger.info("Saving %s to %s.", class_name, path)

    data = _object_to_dict(obj)

    # Add the metadata header fields to the serialized dict.
    provenance = _get_provenance()
    header = {"_format_version": _FORMAT_VERSION, **provenance}
    data = {**header, **data}

    if path.name.endswith(".json.gz"):
        # Use compact format for gzip since readability does not matter and
        # whitespace would increase the pre-compression size.
        json_bytes = json.dumps(data).encode("utf-8")
        with gzip.open(path, "wb") as f:
            f.write(json_bytes)
    else:
        # Use indented format for plain JSON for human readability.
        json_bytes = json.dumps(data, indent=2).encode("utf-8")
        with open(path, "wb") as f:
            f.write(json_bytes)

    file_size = path.stat().st_size
    _logger.info("Saved %s to %s (%d bytes).", type(obj).__name__, path, file_size)


def load(path: str | Path, max_size: int | None = None) -> object:
    """Loads a Ptera Software object from a JSON file.

    If the path ends with ".json.gz", the input is gzip decompressed automatically.

    :param path: The file path to load from.
    :param max_size: The maximum decompressed size in bytes for gzip files. If None, the
        default of 4 GB is used. Set this to a larger value if loading very large
        simulation results. Only applies to ".json.gz" files.
    :return: The deserialized Ptera Software object.
    """
    path = Path(path)
    if not path.name.endswith(".json") and not path.name.endswith(".json.gz"):
        raise ValueError(
            f"Path must end with '.json' or '.json.gz', got '{path.name}'."
        )
    _logger.info("Loading from %s.", path)

    if max_size is None:
        max_size = _DEFAULT_MAX_DECOMPRESSED_SIZE

    if path.name.endswith(".json.gz"):
        with gzip.open(path, "rb") as f:
            raw = f.read(max_size + 1)
            if len(raw) > max_size:
                raise ValueError(
                    f"Decompressed file exceeds the maximum allowed size of "
                    f"{max_size} bytes."
                )
    else:
        with open(path, "rb") as f:
            raw = f.read()

    data = json.loads(raw)

    # Check format version compatibility.
    file_version = data.get("_format_version")
    if file_version != _FORMAT_VERSION:
        raise ValueError(
            f"Format version mismatch: file has version {file_version}, but the "
            f"current code expects version {_FORMAT_VERSION}."
        )

    # Validate that the top level type is a public saveable class.
    top_level_type = data.get("_type")
    if top_level_type not in _PUBLIC_SAVEABLE_CLASSES:
        raise TypeError(
            f"'{top_level_type}' is not a public saveable class. Only files containing "
            f"public Ptera Software classes can be loaded via load()."
        )

    # Log provenance warnings.
    _log_load_warnings(data)

    obj = _object_from_dict(data)
    _logger.info("Loaded %s from %s.", type(obj).__name__, path)
    return obj


def _get_provenance() -> dict[str, str | bool | None]:
    """Returns a dict of provenance metadata for the serialized file.

    The provenance fields are informational only and are never checked at load time. The
    git derived fields (_commit and _dirty) are best effort and are set to None if the
    code is running outside a git repository.

    :return: A dict with provenance metadata.
    """
    try:
        pkg_version = version("PteraSoftware")
    except PackageNotFoundError:  # pragma: no cover
        pkg_version = None

    commit = None
    dirty = None
    try:
        commit = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
            )
            .decode("ascii")
            .strip()
        )
        status = (
            subprocess.check_output(
                ["git", "status", "--porcelain"], stderr=subprocess.DEVNULL
            )
            .decode("ascii")
            .strip()
        )
        dirty = len(status) > 0
    except (FileNotFoundError, subprocess.CalledProcessError):  # pragma: no cover
        _logger.warning("Git is not available. Provenance fields will be null.")

    return {
        "_pterasoftware_version": pkg_version,
        "_commit": commit,
        "_dirty": dirty,
        "_saved_at": datetime.now(timezone.utc).isoformat(),
    }


def _log_load_warnings(data: dict[str, Any]) -> None:
    """Logs warnings about provenance metadata during deserialization.

    :param data: The top level dict loaded from the JSON file.
    :return: None
    """
    if data.get("_dirty"):  # pragma: no branch
        _logger.warning(
            "The file was saved with uncommitted changes (_dirty=true). The _commit "
            "hash may not fully represent the code state at save time."
        )

    file_commit = data.get("_commit")
    if file_commit is not None:
        try:
            current_commit = (
                subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
                )
                .decode("ascii")
                .strip()
            )
            if file_commit != current_commit:  # pragma: no cover
                _logger.warning(
                    "The file was saved at commit %s, but the current HEAD is %s.",
                    file_commit[:12],
                    current_commit[:12],
                )
            current_status = (
                subprocess.check_output(
                    ["git", "status", "--porcelain"], stderr=subprocess.DEVNULL
                )
                .decode("ascii")
                .strip()
            )
            if len(current_status) > 0:  # pragma: no branch
                _logger.warning("The current working tree has uncommitted changes.")
        except (FileNotFoundError, subprocess.CalledProcessError):  # pragma: no cover
            pass


def _object_to_dict(
    obj: object, *, _skip_slots: frozenset[str] = frozenset()
) -> dict[str, Any]:
    """Generically serializes a Ptera Software object to a dict.

    Iterates over the class's __slots__ to discover all attributes, including cache
    slots. JSON keys are the slot names themselves (e.g., "_rho", "_name", "forces_W").
    Slots listed in _skip_slots are serialized as null (used by the shared reference
    optimization for UnsteadyProblem and solver classes).

    :param obj: The Ptera Software object to serialize.
    :param _skip_slots: A frozenset of slot names to serialize as null.
    :return: A dict with "_type" set to the class name and one key per slot.
    """
    cls = type(obj)
    class_name = cls.__name__
    if class_name not in _CLASS_REGISTRY:
        raise TypeError(f"_object_to_dict does not handle {class_name}.")

    _logger.debug("Serializing %s.", class_name)

    # Solver skip slots always apply because the aliases are always redundant
    # with solver._steady_problem or solver.unsteady_problem.
    if isinstance(
        obj,
        (SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver),
    ):
        _skip_slots = _skip_slots | _STEADY_SOLVER_SKIP_SLOTS
    elif isinstance(obj, UnsteadyRingVortexLatticeMethodSolver):
        _skip_slots = _skip_slots | _UNSTEADY_SOLVER_SKIP_SLOTS

    result: dict[str, Any] = {"_type": class_name}
    is_unsteady_problem = isinstance(obj, UnsteadyProblem)
    for slot_name in _all_slots(cls):
        if slot_name in _skip_slots:
            result[slot_name] = None
        elif is_unsteady_problem and slot_name == "_movement":
            # Pass _MOVEMENT_SKIP_SLOTS to the Movement child so that
            # _airplanes and _operating_points are serialized as null.
            result[slot_name] = _object_to_dict(
                getattr(obj, slot_name), _skip_slots=_MOVEMENT_SKIP_SLOTS
            )
        else:
            result[slot_name] = _serialize_value(getattr(obj, slot_name))
    return result


def _object_from_dict(data: dict[str, Any]) -> object:
    """Generically deserializes a Ptera Software object from a dict.

    Uses the "_type" tag to look up the class in the registry. Creates an empty instance
    via object.__new__(cls), bypassing __init__ entirely. Then restores every serialized
    attribute (including caches) via object.__setattr__(). After restoring all slots,
    calls _reconstruct_shared_references() to rebuild any slots that were skipped during
    serialization.

    :param data: The dict produced by _object_to_dict.
    :return: The reconstructed Ptera Software object.
    """
    type_tag = data["_type"]
    cls = _CLASS_REGISTRY.get(type_tag)
    if cls is None:
        raise TypeError(f"Unknown class in _object_from_dict: '{type_tag}'.")

    _logger.debug("Deserializing %s.", type_tag)

    obj: object = object.__new__(cls)
    for slot_name in _all_slots(cls):
        object.__setattr__(obj, slot_name, _deserialize_value(data[slot_name]))
    _reconstruct_shared_references(obj)
    return obj


def _reconstruct_shared_references(obj: object) -> None:
    """Rebuilds shared references that were skipped during serialization.

    Called after _object_from_dict restores all serialized slots. Handles
    UnsteadyProblem (Movement <-> SteadyProblem aliases), steady solvers (solver ->
    SteadyProblem aliases), and the unsteady solver (solver -> UnsteadyProblem aliases).

    :param obj: The deserialized Ptera Software object.
    :return: None
    """
    if isinstance(
        obj,
        (SteadyHorseshoeVortexLatticeMethodSolver, SteadyRingVortexLatticeMethodSolver),
    ):
        # This module is inherently coupled to class internals (it reads __slots__
        # and writes private backing stores for all classes), so accessing a private
        # attribute directly is acceptable here.
        # noinspection PyProtectedMember
        problem = obj._steady_problem
        object.__setattr__(obj, "airplanes", problem.airplanes)
        object.__setattr__(obj, "operating_point", problem.operating_point)
        object.__setattr__(obj, "reynolds_numbers", problem.reynolds_numbers)
        object.__setattr__(obj, "vInf_GP1__E", problem.operating_point.vInf_GP1__E)

        if obj.ran:
            # Reconstruct the flattened panels array (same order as
            # _collapse_geometry).
            panels_list: list[Panel] = []
            for airplane in problem.airplanes:
                for wing in airplane.wings:
                    assert wing.panels is not None
                    panels_list.extend(np.ravel(wing.panels))
            object.__setattr__(obj, "panels", np.array(panels_list, dtype=object))
        else:
            # Pre run: panels is an uninitialized object array sized to
            # num_panels.
            object.__setattr__(obj, "panels", np.empty(obj.num_panels, dtype=object))

    if isinstance(obj, UnsteadyProblem):
        movement = obj.movement
        steady_problems = obj.steady_problems
        num_steps = len(steady_problems)
        num_airplane_movements = len(movement.airplane_movements)

        # Reconstruct Movement._airplanes as a tuple[tuple[Airplane, ...], ...].
        # Outer index = airplane movement, inner index = time step.
        airplanes = tuple(
            tuple(
                steady_problems[step].airplanes[airplane_movement_index]
                for step in range(num_steps)
            )
            for airplane_movement_index in range(num_airplane_movements)
        )
        object.__setattr__(movement, "_airplanes", airplanes)

        # Reconstruct Movement._operating_points: tuple[OperatingPoint, ...].
        # Indexed by time step.
        operating_points = tuple(
            steady_problems[step].operating_point for step in range(num_steps)
        )
        object.__setattr__(movement, "_operating_points", operating_points)

    if isinstance(obj, UnsteadyRingVortexLatticeMethodSolver):
        unsteady_problem = obj.unsteady_problem
        object.__setattr__(obj, "steady_problems", unsteady_problem.steady_problems)

        # This module is inherently coupled to class internals (it reads __slots__
        # and writes private backing stores for all classes), so accessing a private
        # attribute directly is acceptable here.
        # noinspection PyProtectedMember
        current_step = obj._current_step
        current_steady_problem = unsteady_problem.steady_problems[current_step]
        object.__setattr__(
            obj, "current_operating_point", current_steady_problem.operating_point
        )

        if obj.ran:
            object.__setattr__(
                obj, "current_airplanes", current_steady_problem.airplanes
            )

            # Reconstruct flattened panels from current step's airplanes.
            current_panels_list: list[Panel] = []
            for airplane in current_steady_problem.airplanes:
                for wing in airplane.wings:
                    assert wing.panels is not None
                    current_panels_list.extend(np.ravel(wing.panels))
            object.__setattr__(
                obj, "panels", np.array(current_panels_list, dtype=object)
            )
        else:
            # Pre run: current_airplanes is an empty tuple and panels is an
            # empty object array.
            object.__setattr__(obj, "current_airplanes", ())
            object.__setattr__(obj, "panels", np.empty(0, dtype=object))


def _ndarray_to_dict(arr: np.ndarray) -> dict[str, Any]:
    """Converts a NumPy ndarray to a JSON serializable dict.

    For numeric and bool dtypes, the array data is encoded as a base64 string with dtype
    and shape metadata. For dtype=object arrays (e.g., Wing._panels,
    Wing.wake_ring_vortices), each element is serialized individually via
    _serialize_value. The writeable flag is recorded so that deserialization can restore
    the original mutability.

    :param arr: The ndarray to serialize.
    :return: A dict representing the serialized ndarray.
    """
    if arr.dtype == object:
        return {
            "_type": "ndarray",
            "dtype": "object",
            "shape": list(arr.shape),
            "items": [_serialize_value(item) for item in arr.ravel()],
            "writeable": bool(arr.flags.writeable),
        }

    return {
        "_type": "ndarray",
        "dtype": str(arr.dtype),
        "shape": list(arr.shape),
        "data": base64.b64encode(arr.tobytes()).decode("ascii"),
        "writeable": bool(arr.flags.writeable),
    }


def _ndarray_from_dict(array_dict: dict[str, Any]) -> np.ndarray:
    """Reconstructs a NumPy ndarray from a dict produced by _ndarray_to_dict.

    Dispatches on dtype: base64 decode for numeric and bool dtypes, element by element
    deserialization for dtype=object (reshaping to the original shape). After
    reconstruction, restores the writeable flag from the dict's "writeable" field. If
    the field is absent, the array defaults to writeable.

    :param array_dict: The dict produced by _ndarray_to_dict.
    :return: The reconstructed ndarray.
    """
    shape = array_dict["shape"]
    writeable = array_dict.get("writeable", True)

    if array_dict["dtype"] == "object":
        items = [_deserialize_value(item) for item in array_dict["items"]]
        arr = np.empty(len(items), dtype=object)
        for i, item in enumerate(items):
            arr[i] = item
        arr = arr.reshape(shape)
    else:
        raw = base64.b64decode(array_dict["data"])
        arr = (
            np.frombuffer(raw, dtype=np.dtype(array_dict["dtype"]))
            .reshape(shape)
            .copy()
        )

    if not writeable:
        arr.flags.writeable = False

    return arr


def _serialize_value(value: object) -> object:
    """Serializes a single value based on its runtime type.

    Dispatch order: None, bool (before int since bool is a subclass of int), int, float,
    str, ndarray, tuple, list, callable. int and float are wrapped in {"_type": ...,
    "value": ...} dicts rather than serialized as bare JSON numbers to eliminate the
    int/float ambiguity that arises because JSON has a single number type. None, bool,
    and str remain bare JSON values because they map to unambiguous JSON types (null,
    boolean, string).

    :param value: The value to serialize.
    :return: The JSON serializable representation of the value.
    """
    if value is None:
        return None

    # Check bool before int because bool is a subclass of int.
    if isinstance(value, (bool, np.bool_)):
        return bool(value)

    if isinstance(value, (int, np.integer)):
        return {"_type": "int", "value": int(value)}

    if isinstance(value, (float, np.floating)):
        float_value = float(value)
        if not np.isfinite(float_value):
            raise ValueError(f"Cannot serialize non finite float: {float_value}.")
        return {"_type": "float", "value": float_value}

    if isinstance(value, str):
        return value

    if isinstance(value, np.ndarray):
        return _ndarray_to_dict(value)

    if isinstance(value, tuple):
        return {
            "_type": "tuple",
            "items": [_serialize_value(item) for item in value],
        }

    if isinstance(value, list):
        return {
            "_type": "list",
            "items": [_serialize_value(item) for item in value],
        }

    if type(value).__name__ in _CLASS_REGISTRY:
        return _object_to_dict(value)

    if callable(value):
        name = _CALLABLE_FUNC_TO_NAME.get(value)
        if name is None:
            raise ValueError(
                "Only 'sine' and 'uniform' spacing functions are serializable. Custom "
                "callables cannot be serialized."
            )
        return {"_type": "callable", "name": name}

    raise TypeError(f"_serialize_value does not handle {type(value).__name__}.")


def _deserialize_value(data: object) -> object:
    """Deserializes a single value from its JSON representation.

    The format is self describing via _type tags. None, bool, and str are bare JSON
    values. All other types are wrapped in dicts with a _type key. Bare JSON numbers
    (int or float without a _type wrapper) raise a ValueError because all numeric values
    are wrapped during serialization.

    :param data: The serialized data.
    :return: The deserialized value.
    """
    if data is None:
        return None

    if isinstance(data, bool):
        return data

    if isinstance(data, str):
        return data

    if isinstance(data, dict):
        type_tag = data.get("_type")
        if type_tag is None:
            raise ValueError(
                "Dict without '_type' key encountered during deserialization."
            )
        if type_tag == "int":
            return int(data["value"])
        if type_tag == "float":
            return float(data["value"])
        if type_tag == "ndarray":
            return _ndarray_from_dict(data)
        if type_tag == "tuple":
            return tuple(_deserialize_value(item) for item in data["items"])
        if type_tag == "list":
            return [_deserialize_value(item) for item in data["items"]]
        if type_tag == "callable":
            name = data["name"]
            func = _CALLABLE_NAME_TO_FUNC.get(name)
            if func is None:
                raise ValueError(f"Unknown callable name: '{name}'.")
            return func
        if type_tag in _CLASS_REGISTRY:
            return _object_from_dict(data)
        raise TypeError(f"Unknown _type tag: '{type_tag}'.")

    if isinstance(data, (int, float)):
        raise ValueError(
            f"Bare JSON number {data} encountered during deserialization. All numeric "
            "values should be wrapped in _type dicts."
        )

    raise TypeError(  # pragma: no cover
        f"_deserialize_value does not handle {type(data).__name__}."
    )
