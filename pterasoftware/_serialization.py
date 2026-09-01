"""Contains the functions for serializing and deserializing Ptera Software objects."""

from __future__ import annotations

import base64
import hashlib
import inspect
import json
import subprocess
import textwrap
import zipfile
from collections.abc import Callable, Iterator, Mapping
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np

from . import _logging, _parameter_validation

# This module is inherently coupled to the internals of every class in the package (it
# reads __slots__, knows class structure, and imports all classes into its registry), so
# importing from a sibling private module is acceptable here.
# noinspection PyProtectedMember
from ._mujoco_model import MuJoCoModel
from ._oscillation import oscillating_lin_at_time, oscillating_sin_at_time

# noinspection PyProtectedMember
from ._panel import Panel

# noinspection PyProtectedMember
from .aeroelastic_unsteady_ring_vortex_lattice_method import (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver,
)

# noinspection PyProtectedMember
from .free_flight_unsteady_ring_vortex_lattice_method import (
    FreeFlightUnsteadyRingVortexLatticeMethodSolver,
)
from .geometry.airfoil import Airfoil
from .geometry.airplane import Airplane
from .geometry.wing import Wing
from .geometry.wing_cross_section import WingCrossSection
from .movements.aeroelastic_airplane_movement import AeroelasticAirplaneMovement
from .movements.aeroelastic_movement import AeroelasticMovement
from .movements.aeroelastic_wing_cross_section_movement import (
    AeroelasticWingCrossSectionMovement,
)
from .movements.aeroelastic_wing_movement import AeroelasticWingMovement
from .movements.airplane_movement import AirplaneMovement
from .movements.free_flight_movement import FreeFlightMovement
from .movements.free_flight_operating_point_movement import (
    FreeFlightOperatingPointMovement,
)
from .movements.movement import Movement
from .movements.operating_point_movement import OperatingPointMovement
from .movements.wing_cross_section_movement import WingCrossSectionMovement
from .movements.wing_movement import WingMovement
from .operating_point import OperatingPoint

# noinspection PyProtectedMember
from .problems import (
    AeroelasticUnsteadyProblem,
    FreeFlightUnsteadyProblem,
    SteadyProblem,
    UnsteadyProblem,
    _CoupledUnsteadyProblem,
)
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
_FORMAT_VERSION = 28


class UnboundCallable:
    """A placeholder standing in for a custom callable that could not be rebuilt when a
    saved object was loaded.

    **Contains the following methods:**

    None

    **Notes:**

    Saved files store a custom callable (a custom spacing function, an
    AeroelasticWingMovement's second-derivative function, or a
    FreeFlightUnsteadyProblem's external_loads_fn) as an inert marker holding the
    callable's qualified name, its source text when that could be retrieved, and a
    SHA-256 hash of that text. Nothing in the file is ever executed, so loading rebuilds
    the marker as an UnboundCallable rather than as the original function.

    An UnboundCallable passes callable checks, so a loaded object keeps the same
    structure as the one that was saved, and saving it again writes the same marker back
    out. It cannot be invoked: calling it raises a RuntimeError that names the original
    function and prints its recorded source. Loading a solved simulation never invokes
    these callables, so a file holding them loads and can be visualized without any of
    them. A user who needs to re-solve or regenerate motion passes the original
    functions to load's callables argument, keyed by the qualified names recorded in the
    file, and load restores them in place of the placeholders.
    """

    __slots__ = ("_qualname", "_source", "_source_hash")

    def __init__(
        self, qualname: str, source: str | None, source_hash: str | None
    ) -> None:
        """The initialization method.

        :param qualname: The original callable's qualified name, including its module.
        :param source: The original callable's dedented source text, or None if it could
            not be retrieved when the object was saved.
        :param source_hash: The lowercase hexadecimal SHA-256 digest of source encoded
            as UTF-8, or None if source is None.
        :return: None
        """
        if not isinstance(qualname, str):
            raise TypeError(f"qualname must be a str, got {type(qualname).__name__}.")
        self._qualname = qualname

        if source is not None and not isinstance(source, str):
            raise TypeError(
                f"source must be a str or None, got {type(source).__name__}."
            )
        self._source = source

        if source_hash is not None and not isinstance(source_hash, str):
            raise TypeError(
                f"source_hash must be a str or None, got {type(source_hash).__name__}."
            )
        self._source_hash = source_hash

    # --- Immutable: read only properties ---
    @property
    def qualname(self) -> str:
        return self._qualname

    @property
    def source(self) -> str | None:
        return self._source

    @property
    def source_hash(self) -> str | None:
        return self._source_hash

    def __call__(self, *args: object, **kwargs: object) -> object:
        """Raises a RuntimeError naming the original callable and its source.

        :param args: Ignored positional arguments, accepted so the placeholder can stand
            wherever the original callable stood.
        :param kwargs: Ignored keyword arguments, accepted for the same reason.
        :return: Never returns.
        """
        if self._source is None:
            source_note = (
                "Its source text could not be retrieved when the object was saved."
            )
        else:
            source_note = "Its recorded source is:\n\n" + self._source
        raise RuntimeError(
            f"The custom callable {self._qualname} was replaced by a placeholder when "
            "the object was loaded, because saved files store only a callable's name "
            "and source text and never execute anything they contain. To use it, load "
            "the file again and pass the function as load(path, callables="
            f"{{{self._qualname!r}: function}}). {source_note}"
        )

    def __repr__(self) -> str:
        """Returns a string naming the original callable.

        :return: The repr string.
        """
        return f"UnboundCallable(qualname={self._qualname!r})"


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


# This is the default maximum decompressed size in bytes when reading archives. The cap
# is cumulative across every member read during one load(). Prevents zip bombs from
# exhausting memory. Users can override this via the max_size parameter on load().
_DEFAULT_MAX_DECOMPRESSED_SIZE = 4_000_000_000  # 4 GB

# This is the maximum read size in bytes for decompressing archive members during
# load(). Each read is sized to this value or to the remaining allowance under the
# maximum allowed size, whichever is smaller. Reading in bounded chunks lets the
# cumulative decompressed size be checked as each chunk arrives, so a member that
# decompresses far beyond the cap is abandoned after at most one extra chunk instead of
# being decompressed in full.
_MEMBER_READ_CHUNK_SIZE = 16_777_216  # 16 MiB

# These are the names of the archive members that every saved file holds. The header is
# written first and the root object last. Between them, an object with per step
# sequences gets one step member per time step, named by the zero padded step index.
_HEADER_MEMBER_NAME = "header.json"
_ROOT_MEMBER_NAME = "root.json"
_STEP_MEMBER_NAME_FORMAT = "steps/{:08d}.json"

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
    "AeroelasticUnsteadyProblem": AeroelasticUnsteadyProblem,
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver": AeroelasticUnsteadyRingVortexLatticeMethodSolver,
    "FreeFlightMovement": FreeFlightMovement,
    "FreeFlightOperatingPointMovement": FreeFlightOperatingPointMovement,
    "FreeFlightUnsteadyProblem": FreeFlightUnsteadyProblem,
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver": FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    "MuJoCoModel": MuJoCoModel,
}

# These are the classes that can be saved and loaded as top level objects via save() and
# load(). Internal classes (e.g., Panel) are excluded because they are not part of the
# public API and their structure may change without a format version bump. They are
# still serializable as nested objects within public classes.
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
        "AeroelasticUnsteadyProblem",
        "AeroelasticUnsteadyRingVortexLatticeMethodSolver",
        "FreeFlightMovement",
        "FreeFlightOperatingPointMovement",
        "FreeFlightUnsteadyProblem",
        "FreeFlightUnsteadyRingVortexLatticeMethodSolver",
    }
)

# These are the slots on MuJoCoModel that are serialized as null. _model and _data wrap
# native MuJoCo state and are rebuilt from the serialized XML string and the serialized
# _mujoco_assets dict on deserialization via MuJoCoModel._rebuild_engine.
_MUJOCO_MODEL_SKIP_SLOTS: frozenset[str] = frozenset({"_model", "_data"})

# These are the per step sequences that save() splits across the step members instead of
# writing inline in the root member, keyed by the class that owns them. Owners are
# matched by isinstance, so subclasses inherit their parent's entries. The unsteady
# solver's ten lists hold one ndarray or int per time step, and an unsteady problem's
# _steady_problems holds one SteadyProblem per time step, which is where the per step
# Airplanes and their Panels live. UnsteadyProblem and _CoupledUnsteadyProblem are
# siblings that each declare their own _steady_problems slot, so both need an entry.
# Everything else stays inline in the root member on purpose: an UnsteadyProblem's per
# Airplane final load lists and the aeroelastic deformation lists are small, and the
# coupled movements' per step Airplane and OperatingPoint lists serialize as refs into
# the step members because the root member is written last.
_CHUNKED_SLOTS: tuple[tuple[type, tuple[str, ...]], ...] = (
    (
        UnsteadyRingVortexLatticeMethodSolver,
        (
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
        ),
    ),
    (UnsteadyProblem, ("_steady_problems",)),
    (_CoupledUnsteadyProblem, ("_steady_problems",)),
)


def save(path: str | Path, obj: object) -> None:
    """Saves a Ptera Software object to a .psz file.

    A .psz file is a zip archive of compressed JSON members. The first member is a
    header holding the serialization format version, provenance metadata, the saved
    object's class name, and a manifest of the members that follow. The last member
    holds the object itself. For an unsteady solver or an unsteady problem, the per time
    step data (each step's SteadyProblem and the solver's per step bound and wake vortex
    arrays) is split out of the last member into one member per time step, so that
    saving never holds more than one time step's worth of serialized data in memory at
    once on top of the object itself.

    Shared references are preserved. An object referenced from many places (e.g., an
    Airfoil shared by every per step WingCrossSection) is written once, every other
    reference is written as a pointer to it, and load() restores the sharing.

    Custom callables (custom spacing functions, an AeroelasticWingMovement's second-
    derivative functions, and a FreeFlightUnsteadyProblem's external_loads_fn) cannot be
    written to JSON as code. Each is stored as an inert marker holding its qualified
    name, its source text when that can be retrieved, and a hash of that text. The built
    in "sine" and "uniform" spacings are stored by name and restored exactly. See load()
    for how the markers are restored.

    The archive records an internal serialization format version. load() accepts a file
    only when that format version matches the running code's exactly, and there is no
    migration of files written under a different format version. The format version is
    independent of the package version and bumps whenever the serialized structure
    changes. To read a saved file later, run a build of Ptera Software whose format
    version matches the file's.

    :param path: The file path to save to. Should end with ".psz".
    :param obj: The Ptera Software object to save. Must be a public Ptera Software class
        (e.g., Airplane, SteadyProblem, or a solver). Internal classes such as Panel
        cannot be saved directly.
    :return: None
    """
    path = _parameter_validation.pathLike_return_path(path, "path", (".psz",))

    class_name = type(obj).__name__
    if class_name not in _PUBLIC_SAVEABLE_CLASSES:
        raise TypeError(
            f"{class_name} is not a public saveable class. Only public Ptera Software "
            f"classes can be saved via save()."
        )

    _logger.info(_logging.indent() + "Saving %s to %s", class_name, path)

    num_chunks = max(
        (
            len(getattr(owner, slot_name))
            for _, owner, slot_name in _collect_chunked_sequences(obj)
        ),
        default=0,
    )
    header = {
        "_format_version": _FORMAT_VERSION,
        **_get_provenance(),
        "_type": class_name,
        "num_chunks": num_chunks,
        "members": [
            *(_STEP_MEMBER_NAME_FORMAT.format(step) for step in range(num_chunks)),
            _ROOT_MEMBER_NAME,
        ],
    }

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(_HEADER_MEMBER_NAME, json.dumps(header))
        for member_name, payload in _emit_members(obj):
            # Members are dumped without sort_keys and without indent. load() resolves
            # refs in the order the keys appear in the member's text, and sorting the
            # keys could move a ref ahead of the record it points at. Only hash_object
            # sorts keys, because it never resolves refs. Indentation would only inflate
            # the pre-compression size.
            archive.writestr(member_name, json.dumps(payload))
            # Drop the payload before the emitter builds the next one so that only one
            # member's dict tree is alive at a time.
            del payload

    file_size = path.stat().st_size
    _logger.info(
        _logging.indent() + "Saved %s to %s (%d bytes)",
        type(obj).__name__,
        path,
        file_size,
    )


def load(
    path: str | Path,
    max_size: int | None = None,
    callables: Mapping[str, Callable[..., object]] | None = None,
) -> object:
    """Loads a Ptera Software object from a .psz file.

    The archive's header is read first, and its format version and class name are
    checked before any of the object's data is decompressed. The members are then read
    one at a time in the order the header's manifest lists them, and each is
    decompressed in bounded chunks so that the cumulative decompressed size can be
    checked against max_size as it grows. Each step member's data is deserialized and
    the member's text dropped before the next member is read, so loading never holds
    more than one member's worth of text and parsed JSON in memory at once on top of the
    object being rebuilt.

    Nothing in the file is ever executed. A custom callable stored by save() as an inert
    marker is restored as an UnboundCallable, a placeholder that passes callable checks
    but raises a RuntimeError naming the original function and printing its recorded
    source if it is ever invoked. Loading a solved simulation never invokes these
    callables, so such a file loads and can be visualized as usual. To re-solve or
    regenerate motion, pass the original functions in callables, keyed by the qualified
    names the markers recorded (the names the RuntimeError and the file report), and
    they are restored in place of the placeholders. Every key must match a marker in the
    file, so a misspelled key raises rather than silently leaving a placeholder behind.
    When a marker recorded a source hash and the supplied function's own source does not
    hash to the same value (or cannot be retrieved), a warning is logged naming the
    function, because the file was saved with a different definition.

    The file records an internal serialization format version. A file is accepted only
    when that format version matches the running code's exactly, and there is no
    migration of files written under a different format version. A mismatch raises a
    ValueError reporting the file's format version and the running code's. To read the
    file, run a build of Ptera Software whose format version matches the file's.

    :param path: The file path to load from. Should end with ".psz".
    :param max_size: The maximum cumulative decompressed size in bytes of every member
        read from the archive. If None, the default of 4 GB is used. Set this to a
        larger value if loading very large simulation results.
    :param callables: A mapping from the qualified names of custom callables recorded in
        the file to the functions to restore in their place, or None to restore every
        custom callable as an UnboundCallable. The default is None.
    :return: The deserialized Ptera Software object.
    """
    path = Path(path)
    if not path.name.lower().endswith(".psz"):
        raise ValueError(f"Path must end with '.psz', got '{path.name}'.")
    if path.is_dir():
        raise ValueError(f"Path must be a file path, got directory '{path}'.")
    if callables is not None:
        if not isinstance(callables, Mapping):
            raise TypeError(
                f"callables must be a mapping or None, got {type(callables).__name__}."
            )
        for key, func in callables.items():
            if not isinstance(key, str):
                raise TypeError(
                    f"callables keys must be str, got {type(key).__name__}."
                )
            if not callable(func):
                raise TypeError(
                    f"callables values must be callable, but the value for '{key}' is "
                    f"{type(func).__name__}."
                )
    _logger.info(_logging.indent() + "Loading from %s", path)

    if max_size is None:
        max_size = _DEFAULT_MAX_DECOMPRESSED_SIZE

    # The reference table, the chunked sequence elements, and the custom callable hashes
    # accumulate across every member, because a step member's records are referenced
    # from later members and the root member's placeholders splice the elements back in.
    table: dict[int, object] = {}
    chunk_values: dict[str, list[object]] = {}
    recorded_hashes: dict[str, set[str | None]] = {}

    try:
        with zipfile.ZipFile(path) as archive:
            header_bytes = _read_member_bytes(archive, _HEADER_MEMBER_NAME, max_size, 0)
            consumed = len(header_bytes)
            header = json.loads(header_bytes)
            del header_bytes
            if not isinstance(header, dict):
                raise ValueError(
                    f"The file's {_HEADER_MEMBER_NAME} member must hold a JSON "
                    f"object, got {type(header).__name__}."
                )

            # Check format version compatibility. A file loads only when its stored
            # format version matches the running code's, and there is no migration path,
            # so a mismatch is unrecoverable under the running code. The error reports
            # both format versions only. It names no package version to install. The
            # gate keys on the format integer, not the package version, and the stored
            # _pterasoftware_version is provenance that does not reliably identify a
            # build at the file's format version.
            file_version = header.get("_format_version")
            if file_version != _FORMAT_VERSION:
                raise ValueError(
                    f"Format version mismatch: file has format version {file_version}, "
                    f"but the current code uses format version {_FORMAT_VERSION}. A "
                    f"file loads only under a build of Ptera Software whose format "
                    f"version matches the file's."
                )

            # Validate that the top level type is a public saveable class. The header
            # carries the class name, so this check happens before any of the object's
            # data is decompressed.
            top_level_type = header.get("_type")
            if top_level_type not in _PUBLIC_SAVEABLE_CLASSES:
                raise TypeError(
                    f"'{top_level_type}' is not a public saveable class. Only files "
                    f"containing public Ptera Software classes can be loaded via "
                    f"load()."
                )

            # Log provenance warnings.
            _log_load_warnings(header)

            # Validate the manifest. Archive entries that the manifest does not list are
            # ignored, so an archive can carry extra members without affecting a load.
            member_names = header.get("members")
            if (
                not isinstance(member_names, list)
                or not member_names
                or not all(isinstance(name, str) for name in member_names)
            ):
                raise ValueError(
                    "The file's header must list its members as a non empty list of "
                    "member names."
                )
            if len(set(member_names)) != len(member_names):
                raise ValueError("The file's header lists a member more than once.")
            if member_names[-1] != _ROOT_MEMBER_NAME:
                raise ValueError(
                    f"The file's header must list {_ROOT_MEMBER_NAME} as its last "
                    f"member, got '{member_names[-1]}'."
                )

            # Stream the step members in manifest order. Each member maps the key of
            # every chunked sequence with an element at that step to that element, and
            # JSON preserves the key order, so refs within a member resolve in the order
            # they were written.
            for member_name in member_names[:-1]:
                member_bytes = _read_member_bytes(
                    archive, member_name, max_size, consumed
                )
                consumed += len(member_bytes)
                member_data = json.loads(member_bytes)
                del member_bytes
                if not isinstance(member_data, dict):
                    raise ValueError(
                        f"The file's {member_name} member must hold a JSON object, got "
                        f"{type(member_data).__name__}."
                    )
                if callables is not None:
                    _collect_custom_callable_hashes(member_data, recorded_hashes)
                for key, value in member_data.items():
                    chunk_values.setdefault(key, []).append(
                        _deserialize_value(
                            value,
                            table=table,
                            callables=callables,
                            chunk_values=chunk_values,
                        )
                    )
                del member_data

            root_bytes = _read_member_bytes(
                archive, _ROOT_MEMBER_NAME, max_size, consumed
            )
            root_data = json.loads(root_bytes)
            del root_bytes
    except zipfile.BadZipFile as error:
        raise ValueError(f"'{path}' is not a valid .psz file: {error}") from error

    if not isinstance(root_data, dict) or root_data.get("_type") != top_level_type:
        raise ValueError(
            f"The file's {_ROOT_MEMBER_NAME} member must hold a '{top_level_type}' "
            f"record to match its header."
        )

    if callables is not None:
        _collect_custom_callable_hashes(root_data, recorded_hashes)
        _check_callables_against_markers(recorded_hashes, callables)

    obj = _object_from_dict(
        root_data, table=table, callables=callables, chunk_values=chunk_values
    )

    # Splicing a chunked slot removes its key, so any key left over came from a step
    # member and was claimed by no slot in the root member.
    if chunk_values:
        raise ValueError(
            "The file's step members hold keys that no chunked slot claimed: "
            + ", ".join(f"'{key}'" for key in sorted(chunk_values))
            + "."
        )

    _logger.info(_logging.indent() + "Loaded %s from %s", type(obj).__name__, path)
    return obj


def _read_member_bytes(
    archive: zipfile.ZipFile, member_name: str, max_size: int, consumed: int
) -> bytearray:
    """Reads and decompresses one archive member under a cumulative size cap.

    The member is decompressed in bounded chunks, and the read stops as soon as the
    bytes read so far, together with the bytes consumed by earlier members, exceed
    max_size. The sizes recorded in the archive's directory are never trusted, because a
    crafted archive can misstate them, so the bounded reads are the only enforcement.

    :param archive: The open archive to read from.
    :param member_name: The name of the member to read.
    :param max_size: The maximum cumulative decompressed size in bytes across every
        member read during this load.
    :param consumed: The decompressed size in bytes of the members read so far.
    :return: The member's decompressed bytes.
    """
    try:
        member = archive.open(member_name)
    except KeyError:
        raise ValueError(f"The file is missing the member '{member_name}'.") from None

    decompressed = bytearray()
    with member:
        while True:
            # Each read is also capped at the remaining allowance so that the total
            # decompressed data never exceeds the cap by more than one byte, even when
            # the allowance is smaller than the chunk size.
            remaining = max_size + 1 - consumed - len(decompressed)
            chunk = member.read(min(_MEMBER_READ_CHUNK_SIZE, remaining))
            if not chunk:
                break
            decompressed += chunk
            if consumed + len(decompressed) > max_size:
                raise ValueError(
                    f"The file's decompressed size exceeds the maximum allowed size "
                    f"of {max_size} bytes."
                )
    return decompressed


def _collect_chunked_sequences(obj: object) -> list[tuple[str, object, str]]:
    """Collects the per step sequences that are split across step members when a Ptera
    Software object is saved.

    The owners are the object itself and, when the object is an unsteady solver, its
    UnsteadyProblem. Each owner contributes the slots that _CHUNKED_SLOTS declares for
    its class. The order is fixed (the object's own slots first, then its
    UnsteadyProblem's) because it decides the order the sequences' elements are
    serialized in, and with it the memo ids that refs carry and the digest hash_object
    computes.

    :param obj: The Ptera Software object being saved.
    :return: A list of (key, owner, slot name) tuples, one per sequence. The key is the
        slot name when the owner is the object itself, and the owner's path from the
        object joined to the slot name with a slash otherwise. Step members use the key
        to label each sequence's element, and the root member's placeholder for the slot
        carries the same key.
    """
    owners: list[tuple[str, object]] = [("", obj)]
    if isinstance(obj, UnsteadyRingVortexLatticeMethodSolver):
        owners.append(("unsteady_problem", obj.unsteady_problem))

    sequences: list[tuple[str, object, str]] = []
    for owner_path, owner in owners:
        for owner_type, slot_names in _CHUNKED_SLOTS:
            if isinstance(owner, owner_type):
                for slot_name in slot_names:
                    key = slot_name if owner_path == "" else f"{owner_path}/{slot_name}"
                    sequences.append((key, owner, slot_name))
    return sequences


def _emit_members(obj: object) -> Iterator[tuple[str, dict[str, Any]]]:
    """Yields the data members of an archive holding a Ptera Software object, in write
    order.

    The step members come first, one per time step, each mapping the key of every
    chunked sequence that has an element at that step to that element's serialized form.
    A sequence shorter than the longest one simply has no key in the later step members.
    The root member comes last and holds the object's own record, with each chunked slot
    replaced by a placeholder naming its key, its container type, and its length.

    One identity memo threads across every member, so a ref always points at a record in
    the same member or an earlier one. Writing the root member last is what lets the
    object's other references to per step data (a solver's current Airplanes, or a
    Movement's per step Airplane lists) serialize as cheap refs into the step members.
    save() writes the members in this order, and hash_object() folds them into its
    digest in this order, so the two agree on the memo ids that the refs carry.

    :param obj: The Ptera Software object to serialize.
    :return: An iterator of (member name, member payload) tuples ending with the root
        member, which holds the object's own record.
    """
    sequences = _collect_chunked_sequences(obj)
    named_sequences = [
        (key, owner, slot_name, getattr(owner, slot_name))
        for key, owner, slot_name in sequences
    ]
    num_chunks = max((len(sequence) for *_, sequence in named_sequences), default=0)

    memo: dict[int, tuple[int, object]] = {}
    for step in range(num_chunks):
        yield _STEP_MEMBER_NAME_FORMAT.format(step), {
            key: _serialize_value(sequence[step], memo=memo)
            for key, _, _, sequence in named_sequences
            if step < len(sequence)
        }

    # Nothing reachable from a per step element refers back to a chunked slot's owner,
    # so no owner can have been serialized inside a step member. If one had been, its
    # chunked slots would have been written inline there, and the root member's
    # placeholders for them would splice nothing.
    chunk_placeholders: dict[int, dict[str, dict[str, Any]]] = {}
    for key, owner, slot_name, sequence in named_sequences:
        if id(owner) in memo:
            raise RuntimeError(
                f"The owner of the chunked slot '{key}' was reached from a step member "
                "before the root member was written."
            )
        if isinstance(sequence, tuple):
            container = "tuple"
        elif isinstance(sequence, list):
            container = "list"
        else:
            raise TypeError(
                f"The chunked slot '{key}' must hold a tuple or a list, got "
                f"{type(sequence).__name__}."
            )
        chunk_placeholders.setdefault(id(owner), {})[slot_name] = {
            "_type": "chunked",
            "key": key,
            "container": container,
            "length": len(sequence),
        }
    yield _ROOT_MEMBER_NAME, _object_to_dict(
        obj, memo=memo, chunk_placeholders=chunk_placeholders
    )


def _collect_custom_callable_hashes(
    data: object, recorded_hashes: dict[str, set[str | None]]
) -> None:
    """Collects the recorded source hashes of every custom callable marker in one parsed
    member, keyed by qualified name.

    Walks the parsed JSON without constructing anything. The hashes accumulate across
    members so that _check_callables_against_markers can see every marker in the file
    even though each member's parsed data is dropped after it is deserialized.

    :param data: The parsed JSON of one archive member.
    :param recorded_hashes: The dict mapping each marker's qualified name to the set of
        source hashes recorded for it, which this call adds to.
    :return: None
    """
    pending: list[object] = [data]
    while pending:
        node = pending.pop()
        if isinstance(node, dict):
            if node.get("_type") == "custom_callable":
                recorded_hashes.setdefault(node["qualname"], set()).add(
                    node["source_hash"]
                )
            else:
                pending.extend(node.values())
        elif isinstance(node, list):
            pending.extend(node)


def _check_callables_against_markers(
    recorded_hashes: dict[str, set[str | None]],
    callables: Mapping[str, Callable[..., object]],
) -> None:
    """Checks a callables mapping against the custom callable markers a file recorded.

    A key that matches no marker raises, so a typo fails fast instead of leaving a
    placeholder in place. A key whose function does not hash to a recorded hash
    (including a function whose source cannot be retrieved) logs one warning naming the
    function. A marker that recorded no hash gives nothing to check against, so it never
    warns.

    :param recorded_hashes: The dict mapping each marker's qualified name to the set of
        source hashes recorded for it across every member of the file, as built by
        _collect_custom_callable_hashes.
    :param callables: The validated mapping passed to load().
    :return: None
    """
    unused = sorted(key for key in callables if key not in recorded_hashes)
    if unused:
        raise ValueError(
            "callables has keys that match no custom callable in the file: "
            + ", ".join(f"'{key}'" for key in unused)
            + ". The file's custom callables are: "
            + ", ".join(f"'{key}'" for key in sorted(recorded_hashes))
            + "."
        )

    for qualname, func in callables.items():
        rebound_hash = _custom_callable_to_dict(func)["source_hash"]
        recorded = recorded_hashes[qualname]
        if any(
            this_hash is not None and this_hash != rebound_hash
            for this_hash in recorded
        ):
            _logger.warning(
                _logging.indent()
                + "The function supplied for %s does not match the source the file "
                "recorded for it, so the file was saved with a different definition",
                qualname,
            )


def hash_object(obj: object) -> str:
    """Returns a hex digest identifying a Ptera Software object's content.

    The digest is the SHA-256 hash of the object's JSON serialization, computed over
    every serialized slot (including cache slots) and folded together with the current
    serialization format version. The serialization is folded in member by member, in
    the same order save() writes the members, with each member's name delimiting its
    payload. Provenance metadata is never part of the digest, because the header member
    is not a data member. Two object graphs with identical content and identical
    internal sharing structure produce the same digest regardless of instance identity,
    and any difference in serialized content produces a different digest. Because shared
    references serialize as refs, two graphs whose contents are equal but whose aliasing
    differs (e.g., one Airfoil shared by every WingCrossSection versus an equal but
    distinct Airfoil per WingCrossSection) produce different digests. The digest
    survives a save and load round trip, but it is only stable within a single
    serialization format version, because a format version change is folded into the
    hash and shifts every digest.

    :param obj: The Ptera Software object to hash. Must be an instance of a class
        registered for serialization; an unregistered type raises a TypeError.
    :return: A 64 character lowercase hexadecimal string holding the SHA-256 digest of
        the object's content.
    """
    digest = hashlib.sha256()
    digest.update(json.dumps({"_format_version": _FORMAT_VERSION}).encode("utf-8"))
    for member_name, payload in _emit_members(obj):
        # Sorting the keys here is safe because hashing never resolves refs, unlike
        # save(), which must dump each member in insertion order so that every ref
        # follows the record it points at.
        digest.update(member_name.encode("utf-8"))
        digest.update(json.dumps(payload, sort_keys=True).encode("utf-8"))
        del payload
    return digest.hexdigest()


def _get_provenance() -> dict[str, str | bool | None]:
    """Returns a dict of provenance metadata for the serialized file.

    The provenance fields are informational only and are never checked at load time. The
    git derived fields (_commit and _dirty) are best effort and are set to None unless
    the package directory sits inside a Ptera Software development checkout (its parent
    is the repository toplevel), which covers both installs outside any git repository
    and installs whose venv lives inside an unrelated repository.

    :return: A dict with provenance metadata.
    """
    try:
        pkg_version = version("PteraSoftware")
    except PackageNotFoundError:  # pragma: no cover
        pkg_version = None

    commit = None
    dirty = None
    package_dir = Path(__file__).resolve().parent
    try:
        toplevel = (
            subprocess.check_output(
                ["git", "rev-parse", "--show-toplevel"],
                cwd=package_dir,
                stderr=subprocess.DEVNULL,
            )
            .decode("ascii")
            .strip()
        )
        if Path(toplevel).resolve() != package_dir.parent.resolve():
            return {
                "_pterasoftware_version": pkg_version,
                "_commit": None,
                "_dirty": None,
                "_saved_at": datetime.now(timezone.utc).isoformat(),
            }
        commit = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=package_dir,
                stderr=subprocess.DEVNULL,
            )
            .decode("ascii")
            .strip()
        )
        status = (
            subprocess.check_output(
                ["git", "status", "--porcelain"],
                cwd=package_dir,
                stderr=subprocess.DEVNULL,
            )
            .decode("ascii")
            .strip()
        )
        dirty = len(status) > 0
    except (
        FileNotFoundError,
        subprocess.CalledProcessError,
        UnicodeDecodeError,
    ):  # pragma: no cover
        _logger.debug(
            _logging.indent()
            + "The package git state could not be read, so the provenance fields "
            "will be null"
        )

    return {
        "_pterasoftware_version": pkg_version,
        "_commit": commit,
        "_dirty": dirty,
        "_saved_at": datetime.now(timezone.utc).isoformat(),
    }


def _log_load_warnings(data: dict[str, Any]) -> None:
    """Logs warnings about provenance metadata during deserialization.

    :param data: The dict loaded from the archive's header member.
    :return: None
    """
    if data.get("_dirty"):  # pragma: no branch
        _logger.warning(
            _logging.indent()
            + "The file was saved with uncommitted changes, so the hash may not "
            "fully represent the code state"
        )

    file_commit = data.get("_commit")
    if file_commit is not None:
        try:
            package_dir = Path(__file__).resolve().parent
            toplevel = (
                subprocess.check_output(
                    ["git", "rev-parse", "--show-toplevel"],
                    cwd=package_dir,
                    stderr=subprocess.DEVNULL,
                )
                .decode("ascii")
                .strip()
            )
            if Path(toplevel).resolve() != package_dir.parent.resolve():
                return
            current_commit = (
                subprocess.check_output(
                    ["git", "rev-parse", "HEAD"],
                    cwd=package_dir,
                    stderr=subprocess.DEVNULL,
                )
                .decode("ascii")
                .strip()
            )
            if file_commit != current_commit:  # pragma: no cover
                _logger.warning(
                    _logging.indent()
                    + "The file was saved at commit %s, but the current HEAD is %s",
                    file_commit[:12],
                    current_commit[:12],
                )
            current_status = (
                subprocess.check_output(
                    ["git", "status", "--porcelain"],
                    cwd=package_dir,
                    stderr=subprocess.DEVNULL,
                )
                .decode("ascii")
                .strip()
            )
            if len(current_status) > 0:  # pragma: no branch
                _logger.warning(
                    _logging.indent()
                    + "The current working tree has uncommitted changes"
                )
        except (
            FileNotFoundError,
            subprocess.CalledProcessError,
            UnicodeDecodeError,
        ):  # pragma: no cover
            pass


def _object_to_dict(
    obj: object,
    *,
    memo: dict[int, tuple[int, object]] | None = None,
    chunk_placeholders: dict[int, dict[str, dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    """Generically serializes a Ptera Software object to a dict.

    Iterates over the class's __slots__ to discover all attributes, including cache
    slots. JSON keys are the slot names themselves (e.g., "_rho", "_name", "forces_W").
    The MuJoCoModel's native engine objects are serialized as null and rebuilt on
    deserialization. A slot that chunk_placeholders names for this object is written as
    its placeholder instead of being serialized, because its elements were already
    written to the step members.

    Shared references are preserved via an identity memo. The first encounter of each
    object serializes it fully and stamps the dict with an "_id" key holding a small
    integer. Every later encounter of the same object (by identity) serializes as a
    {"_type": "ref", "id": ...} dict pointing at that integer, so an object that is
    referenced from many places (e.g., an Airfoil shared by every per step
    WingCrossSection) is stored only once and its sharing survives the round trip.

    :param obj: The Ptera Software object to serialize.
    :param memo: The identity memo mapping an object's id() to its "_id" integer and the
        object itself (held to keep its id() from being recycled during the walk). If
        None, a fresh memo is created, making this a top level call.
    :param chunk_placeholders: A dict mapping the id() of each object that owns chunked
        sequences to a dict mapping each chunked slot's name to the placeholder dict to
        write in its place. It is threaded through every recursion so that an owner
        nested anywhere in the graph is found. If None, no slot is replaced, which is
        the case for every in memory conversion outside of save() and hash_object().
    :return: A dict with "_type" set to the class name, "_id" set to the memo integer,
        and one key per slot, or a {"_type": "ref", "id": ...} dict if the object was
        already serialized under this memo.
    """
    cls = type(obj)
    class_name = cls.__name__
    if class_name not in _CLASS_REGISTRY:
        raise TypeError(f"_object_to_dict does not handle {class_name}.")

    if memo is None:
        memo = {}

    memoized = memo.get(id(obj))
    if memoized is not None:
        return {"_type": "ref", "id": memoized[0]}

    _logger.debug(_logging.indent() + "Serializing %s", class_name)

    # The MuJoCoModel's native model and data objects cannot be serialized, so they are
    # rebuilt from the XML string and the assets dict on deserialization.
    skip_slots: frozenset[str] = frozenset()
    if isinstance(obj, MuJoCoModel):
        skip_slots = _MUJOCO_MODEL_SKIP_SLOTS

    placeholders: dict[str, dict[str, Any]] = {}
    if chunk_placeholders is not None:
        placeholders = chunk_placeholders.get(id(obj), {})

    # Register the object in the memo before serializing its slots so that any reference
    # back to it from within its own subtree serializes as a ref.
    ref_id = len(memo)
    memo[id(obj)] = (ref_id, obj)

    result: dict[str, Any] = {"_type": class_name, "_id": ref_id}
    for slot_name in _all_slots(cls):
        if slot_name in placeholders:
            result[slot_name] = placeholders[slot_name]
        elif slot_name in skip_slots:
            result[slot_name] = None
        else:
            result[slot_name] = _serialize_value(
                getattr(obj, slot_name),
                memo=memo,
                chunk_placeholders=chunk_placeholders,
            )
    return result


def _object_from_dict(
    data: dict[str, Any],
    *,
    table: dict[int, object] | None = None,
    callables: Mapping[str, Callable[..., object]] | None = None,
    chunk_values: dict[str, list[object]] | None = None,
) -> object:
    """Generically deserializes a Ptera Software object from a dict.

    Uses the "_type" tag to look up the class in the registry. Creates an empty instance
    via object.__new__(cls), bypassing __init__ entirely. Then restores every serialized
    attribute (including caches) via object.__setattr__(). For a MuJoCoModel, rebuilds
    the native engine objects that were serialized as null.

    The instance is registered in the table under the dict's "_id" integer before its
    slots are restored, so {"_type": "ref", "id": ...} dicts elsewhere in the data
    resolve to the same instance and shared references survive the round trip.

    :param data: The dict produced by _object_to_dict.
    :param table: The reference table mapping "_id" integers to their reconstructed
        instances. If None, a fresh table is created, making this a top level call.
    :param callables: The mapping from qualified names to functions threaded through to
        _deserialize_value so custom callable markers can be rebound. If None, every
        marker deserializes to an UnboundCallable.
    :param chunk_values: The dict mapping each chunked sequence's key to the elements
        already deserialized from the step members, threaded through to
        _deserialize_value so a chunked slot's placeholder can be spliced. If None, a
        placeholder raises, which is the case for every in memory conversion outside of
        load().
    :return: The reconstructed Ptera Software object.
    """
    type_tag = data["_type"]
    cls = _CLASS_REGISTRY.get(type_tag)
    if cls is None:
        raise TypeError(f"Unknown class in _object_from_dict: '{type_tag}'.")

    if table is None:
        table = {}

    _logger.debug(_logging.indent() + "Deserializing %s", type_tag)

    obj: object = object.__new__(cls)
    table[data["_id"]] = obj
    for slot_name in _all_slots(cls):
        object.__setattr__(
            obj,
            slot_name,
            _deserialize_value(
                data[slot_name],
                table=table,
                callables=callables,
                chunk_values=chunk_values,
            ),
        )

    if isinstance(obj, MuJoCoModel):
        # The native model and data objects were skipped during serialization. Rebuild
        # them from the deserialized XML string. This module is inherently coupled to
        # class internals, so calling a private method directly is acceptable here.
        # noinspection PyProtectedMember
        obj._rebuild_engine()
    return obj


def _ndarray_to_dict(
    arr: np.ndarray,
    *,
    memo: dict[int, tuple[int, object]] | None = None,
    chunk_placeholders: dict[int, dict[str, dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    """Converts a NumPy ndarray to a JSON serializable dict.

    For numeric and bool dtypes, the array data is encoded as a base64 string with dtype
    and shape metadata. For dtype=object arrays (e.g., Wing._panels,
    Wing.wake_ring_vortices), each element is serialized individually via
    _serialize_value. The writeable flag is recorded so that deserialization can restore
    the original mutability.

    :param arr: The ndarray to serialize.
    :param memo: The identity memo threaded through to _serialize_value for dtype=object
        elements. If None, a fresh memo is created.
    :param chunk_placeholders: The chunked slot placeholders threaded through to
        _serialize_value for dtype=object elements. If None, no slot is replaced.
    :return: A dict representing the serialized ndarray.
    """
    if arr.dtype == object:
        if memo is None:
            memo = {}
        return {
            "_type": "ndarray",
            "dtype": "object",
            "shape": list(arr.shape),
            "items": [
                _serialize_value(item, memo=memo, chunk_placeholders=chunk_placeholders)
                for item in arr.ravel()
            ],
            "writeable": bool(arr.flags.writeable),
        }

    return {
        "_type": "ndarray",
        "dtype": str(arr.dtype),
        "shape": list(arr.shape),
        "data": base64.b64encode(arr.tobytes()).decode("ascii"),
        "writeable": bool(arr.flags.writeable),
    }


def _ndarray_from_dict(
    array_dict: dict[str, Any],
    *,
    table: dict[int, object] | None = None,
    callables: Mapping[str, Callable[..., object]] | None = None,
    chunk_values: dict[str, list[object]] | None = None,
) -> np.ndarray:
    """Reconstructs a NumPy ndarray from a dict produced by _ndarray_to_dict.

    Dispatches on dtype: base64 decode for numeric and bool dtypes, element by element
    deserialization for dtype=object (reshaping to the original shape). After
    reconstruction, restores the writeable flag from the dict's "writeable" field. If
    the field is absent, the array defaults to writeable.

    :param array_dict: The dict produced by _ndarray_to_dict.
    :param table: The reference table threaded through to _deserialize_value for
        dtype=object elements. If None, a fresh table is created.
    :param callables: The mapping from qualified names to functions threaded through to
        _deserialize_value for dtype=object elements. If None, every custom callable
        marker deserializes to an UnboundCallable.
    :param chunk_values: The chunked sequence elements threaded through to
        _deserialize_value for dtype=object elements. If None, a chunked slot
        placeholder raises.
    :return: The reconstructed ndarray.
    """
    shape = array_dict["shape"]
    writeable = array_dict.get("writeable", True)

    if array_dict["dtype"] == "object":
        if table is None:
            table = {}
        items = [
            _deserialize_value(
                item, table=table, callables=callables, chunk_values=chunk_values
            )
            for item in array_dict["items"]
        ]
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


def _serialize_value(
    value: object,
    *,
    memo: dict[int, tuple[int, object]] | None = None,
    chunk_placeholders: dict[int, dict[str, dict[str, Any]]] | None = None,
) -> object:
    """Serializes a single value based on its runtime type.

    Dispatch order: None, bool (before int since bool is a subclass of int), int, float,
    str, bytes, ndarray, tuple, list, dict, registered class instance, UnboundCallable,
    callable. The built in "sine" and "uniform" spacing functions serialize by name, and
    every other callable serializes as an inert custom callable marker (see
    _custom_callable_to_dict). int and float are wrapped in {"_type": ..., "value": ...}
    dicts rather than serialized as bare JSON numbers to eliminate the int/float
    ambiguity that arises because JSON has a single number type. None, bool, and str
    remain bare JSON values because they map to unambiguous JSON types (null, boolean,
    string). bytes are encoded as base64 strings, the same encoding _ndarray_to_dict
    uses for array buffers. dicts must have str keys because JSON object keys are
    strings, and a lossy key coercion would break the round trip.

    :param value: The value to serialize.
    :param memo: The identity memo threaded through to _object_to_dict for registered
        class instances, so shared references serialize once. If None, a fresh memo is
        created, making this a top level call.
    :param chunk_placeholders: The chunked slot placeholders threaded through to
        _object_to_dict for registered class instances and through every container
        recursion. If None, no slot is replaced.
    :return: The JSON serializable representation of the value.
    """
    if memo is None:
        memo = {}

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

    if isinstance(value, bytes):
        return {
            "_type": "bytes",
            "data": base64.b64encode(value).decode("ascii"),
        }

    if isinstance(value, np.ndarray):
        return _ndarray_to_dict(value, memo=memo, chunk_placeholders=chunk_placeholders)

    if isinstance(value, tuple):
        return {
            "_type": "tuple",
            "items": [
                _serialize_value(item, memo=memo, chunk_placeholders=chunk_placeholders)
                for item in value
            ],
        }

    if isinstance(value, list):
        return {
            "_type": "list",
            "items": [
                _serialize_value(item, memo=memo, chunk_placeholders=chunk_placeholders)
                for item in value
            ],
        }

    if isinstance(value, dict):
        for key in value:
            if not isinstance(key, str):
                raise TypeError(
                    "_serialize_value only handles dicts with str keys, not "
                    f"{type(key).__name__} keys."
                )
        return {
            "_type": "dict",
            "items": {
                key: _serialize_value(
                    item, memo=memo, chunk_placeholders=chunk_placeholders
                )
                for key, item in value.items()
            },
        }

    if type(value).__name__ in _CLASS_REGISTRY:
        return _object_to_dict(value, memo=memo, chunk_placeholders=chunk_placeholders)

    # An UnboundCallable re-emits the marker it was loaded from, so a loaded object
    # saves and hashes identically to the original without re-inspecting anything.
    if isinstance(value, UnboundCallable):
        return {
            "_type": "custom_callable",
            "qualname": value.qualname,
            "source": value.source,
            "source_hash": value.source_hash,
        }

    if callable(value):
        name = _CALLABLE_FUNC_TO_NAME.get(value)
        if name is None:
            return _custom_callable_to_dict(value)
        return {"_type": "callable", "name": name}

    raise TypeError(f"_serialize_value does not handle {type(value).__name__}.")


def _custom_callable_to_dict(value: Callable[..., object]) -> dict[str, Any]:
    """Serializes a custom callable as an inert marker.

    The marker records the callable's qualified name, its dedented source text when
    inspect.getsource can retrieve it, and the SHA-256 hash of that text. A callable
    without a __qualname__ of its own (a functools.partial or an instance of a class
    defining __call__) is named after its type. Source retrieval fails for builtins,
    partials, callable instances, and functions defined in a REPL or by exec, in which
    case the source and its hash are both null. The marker is deserialized as an
    UnboundCallable, and nothing stored in it is ever executed.

    :param value: The custom callable to serialize.
    :return: The marker dict.
    """
    qualname = getattr(value, "__qualname__", None)
    if isinstance(qualname, str):
        module = getattr(value, "__module__", None)
    else:
        qualname = type(value).__qualname__
        module = type(value).__module__
    if isinstance(module, str):
        qualname = f"{module}.{qualname}"

    try:
        source: str | None = textwrap.dedent(inspect.getsource(value))
    except (OSError, TypeError):
        source = None

    source_hash = None
    if source is not None:
        source_hash = hashlib.sha256(source.encode("utf-8")).hexdigest()

    return {
        "_type": "custom_callable",
        "qualname": qualname,
        "source": source,
        "source_hash": source_hash,
    }


def _deserialize_value(
    data: object,
    *,
    table: dict[int, object] | None = None,
    callables: Mapping[str, Callable[..., object]] | None = None,
    chunk_values: dict[str, list[object]] | None = None,
) -> object:
    """Deserializes a single value from its JSON representation.

    The format is self describing via _type tags. None, bool, and str are bare JSON
    values. All other types are wrapped in dicts with a _type key. Bare JSON numbers
    (int or float without a _type wrapper) raise a ValueError because all numeric values
    are wrapped during serialization. A "chunked" placeholder, which stands in for a per
    step sequence that save() split across the step members, is spliced back together
    from chunk_values.

    :param data: The serialized data.
    :param table: The reference table mapping "_id" integers to their reconstructed
        instances, used to resolve {"_type": "ref", "id": ...} dicts. If None, a fresh
        table is created, making this a top level call.
    :param callables: A mapping from qualified names to functions. A custom callable
        marker whose qualified name is a key deserializes to that key's function, and
        every other marker deserializes to an UnboundCallable. Keys that match no marker
        are ignored here, because only load() sees the whole file and can tell that a
        key went unused. If None, every marker deserializes to an UnboundCallable.
    :param chunk_values: A dict mapping each chunked sequence's key to the elements
        already deserialized from the step members, in step order. Splicing a
        placeholder removes its key, so load() can tell afterwards whether any step
        member held a key that no slot claimed. If None, a placeholder raises, because
        only load() holds the step members it splices from.
    :return: The deserialized value.
    """
    if table is None:
        table = {}

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
        if type_tag == "bytes":
            return base64.b64decode(data["data"])
        if type_tag == "ndarray":
            return _ndarray_from_dict(
                data, table=table, callables=callables, chunk_values=chunk_values
            )
        if type_tag == "ref":
            ref_id = data["id"]
            if ref_id not in table:
                raise ValueError(
                    f"Reference to unknown shared object id {ref_id} encountered "
                    "during deserialization."
                )
            return table[ref_id]
        if type_tag == "tuple":
            return tuple(
                _deserialize_value(
                    item, table=table, callables=callables, chunk_values=chunk_values
                )
                for item in data["items"]
            )
        if type_tag == "list":
            return [
                _deserialize_value(
                    item, table=table, callables=callables, chunk_values=chunk_values
                )
                for item in data["items"]
            ]
        if type_tag == "dict":
            return {
                key: _deserialize_value(
                    item, table=table, callables=callables, chunk_values=chunk_values
                )
                for key, item in data["items"].items()
            }
        if type_tag == "chunked":
            if chunk_values is None:
                raise ValueError(
                    "A chunked slot placeholder was encountered outside of load(), "
                    "which is the only reader that holds the step members it splices "
                    "from."
                )
            key = data["key"]
            items = chunk_values.pop(key, [])
            if len(items) != data["length"]:
                raise ValueError(
                    f"The chunked slot '{key}' should hold {data['length']} elements, "
                    f"but the file's step members hold {len(items)}."
                )
            container = data["container"]
            if container == "tuple":
                return tuple(items)
            if container == "list":
                return items
            raise ValueError(f"Unknown chunked container: '{container}'.")
        if type_tag == "callable":
            name = data["name"]
            func = _CALLABLE_NAME_TO_FUNC.get(name)
            if func is None:
                raise ValueError(f"Unknown callable name: '{name}'.")
            return func
        if type_tag == "custom_callable":
            qualname = data["qualname"]
            if callables is not None and qualname in callables:
                return callables[qualname]
            return UnboundCallable(
                qualname=qualname,
                source=data["source"],
                source_hash=data["source_hash"],
            )
        if type_tag in _CLASS_REGISTRY:
            return _object_from_dict(
                data, table=table, callables=callables, chunk_values=chunk_values
            )
        raise TypeError(f"Unknown _type tag: '{type_tag}'.")

    if isinstance(data, (int, float)):
        raise ValueError(
            f"Bare JSON number {data} encountered during deserialization. All numeric "
            "values should be wrapped in _type dicts."
        )

    raise TypeError(  # pragma: no cover
        f"_deserialize_value does not handle {type(data).__name__}."
    )
