"""Contains functions for caching convergence solves and mesh-resolution memos."""

from __future__ import annotations

import json
import os
import time
from collections.abc import Callable
from pathlib import Path

import numpy as np

from . import _logging

_logger = _logging.get_logger("_convergence_cache")

# The schema version of the JSON solve cache. It is bumped only when the on-disk
# structure of the cache file changes, for example how the load coefficients are stored
# or the addition of the memo section. A cache file whose version does not match is
# ignored and rebuilt from scratch. This is independent of the serialization format
# version that hash_object folds into each key, which guards against changes to the
# reference problem's serialized structure.
_SOLVE_CACHE_VERSION = 3


def solve_cache_key(ref_problem_hash: str, *components: object) -> str:
    """Builds a stable solve-cache key from a reference problem's content hash and the
    parameters that determine a single solve.

    The content hash identifies the geometry, operating point, and motion being
    converged, and the remaining components identify one mesh in the sweep, for example
    the solver type, the Panel aspect ratio, and the number of chordwise Panels. The
    components are joined into one string, so their order is significant and every
    parameter that changes the solve must be included.

    :param ref_problem_hash: The content hash of the reference problem, as returned by
        the serialization module's hash_object.
    :param components: The remaining parameters that determine the solve, in a fixed
        order. Each is converted to a string.
    :return: A string that uniquely identifies the solve within the cache.
    """
    return "|".join([ref_problem_hash, *(str(component) for component in components)])


def load_solve_cache(
    cache_path: Path | None,
) -> dict[str, tuple[np.ndarray, float]]:
    """Loads a JSON solve cache from disk into a dictionary of load coefficients and
    solve times.

    The cache is a pure optimization, so any file that cannot be read as a current-
    schema solve cache yields an empty cache rather than an error and the sweep rebuilds
    it by solving each mesh. This covers a missing file, an unreadable file (for example
    a directory or a permission error), a file that is not valid JSON, a file whose
    schema version does not match the current one, and a current-version file whose
    entries section is missing or malformed. An unreadable or invalid file is logged as
    a warning so a silently degraded cache is still observable.

    :param cache_path: The path to the JSON solve cache file, or None to skip caching
        and return an empty cache.
    :return: A dict mapping each cache key to a tuple of a (num_airplanes, 6) ndarray of
        floats holding that mesh's per-Airplane load coefficients (the three force
        coefficients (in wind axes) followed by the three moment coefficients (in wind
        axes, relative to the first Airplane's CG)) and a float holding the time, in
        seconds, its solve took when it ran.
    """
    if cache_path is None or not cache_path.exists():
        return {}

    try:
        with open(cache_path) as cache_file:
            data = json.load(cache_file)
    except (OSError, json.JSONDecodeError) as error:
        _logger.warning(
            "Ignoring unreadable solve cache at %s (%s); rebuilding it.",
            cache_path,
            error,
        )
        return {}

    if data.get("_cache_version") != _SOLVE_CACHE_VERSION:
        return {}

    entries = data.get("entries", {})
    if not isinstance(entries, dict):
        _logger.warning(
            "Ignoring solve cache at %s with a malformed entries section; rebuilding it.",
            cache_path,
        )
        return {}

    return {
        key: (
            np.array(entry["coefficients"], dtype=float),
            float(entry["solve_time"]),
        )
        for key, entry in entries.items()
    }


def write_cache(
    cache_path: Path,
    solve_cache: dict[str, tuple[np.ndarray, float]],
    memo_cache: dict[str, float],
) -> None:
    """Writes the whole cache (solved coefficients, solve times, and memos) to disk as
    JSON, replacing any existing file atomically.

    Both sections are written together as one file image, so neither the solved
    coefficients nor the memos can clobber the other. The file is written to a temporary
    file in the same directory and then moved into place with os.replace, so a run
    interrupted mid-write leaves the existing cache intact rather than a half-written
    file. There is no file locking, so concurrent convergence runs sharing one cache
    file are not supported.

    :param cache_path: The path to the JSON cache file.
    :param solve_cache: A dict mapping each solve-cache key to a tuple of a
        (num_airplanes, 6) ndarray of floats holding that mesh's per-Airplane load
        coefficients (the three force coefficients (in wind axes) followed by the three
        moment coefficients (in wind axes, relative to the first Airplane's CG)) and a
        float holding the time, in seconds, its solve took when it ran.
    :param memo_cache: A dict mapping each absolute-valued memo key to its value, as
        returned by memos_to_disk.
    :return: None
    """
    data = {
        "_cache_version": _SOLVE_CACHE_VERSION,
        "entries": {
            key: {
                "coefficients": coefficients.tolist(),
                "solve_time": solve_time,
            }
            for key, (coefficients, solve_time) in solve_cache.items()
        },
        "memos": memo_cache,
    }

    temp_path = cache_path.parent / (cache_path.name + ".tmp")
    with open(temp_path, "w") as cache_file:
        json.dump(data, cache_file)
    os.replace(temp_path, cache_path)


def cached_solve(
    cache: dict[str, tuple[np.ndarray, float]],
    key: str,
    solve: Callable[[], np.ndarray],
    persist: Callable[[], None] | None = None,
) -> tuple[np.ndarray, float]:
    """Returns a mesh's load coefficients and solve time from the cache, solving,
    timing, and storing them on a miss.

    On a hit, the stored coefficients and solve time are returned and the expensive
    solve is skipped. On a miss, solve is called and timed, its result is added to the
    in-memory cache, and, when a persist callback is given, it is invoked to write the
    whole cache through to disk so an interrupted study keeps every iteration solved so
    far. The returned solve time is the time the solve took when it actually ran, which
    on a hit may have been during an earlier run or on another machine. The returned
    array is the one held in the cache and must be treated as read-only by the caller.

    :param cache: The in-memory cache mapping each key to a tuple of a (num_airplanes,
        6) ndarray of floats and a float solve time in seconds, as returned by
        load_solve_cache. It is updated in place on a miss.
    :param key: The cache key identifying this solve, as returned by solve_cache_key.
    :param solve: A callable that builds and solves the mesh and returns its
        (num_airplanes, 6) ndarray of floats of per-Airplane load coefficients (the
        three force coefficients (in wind axes) followed by the three moment
        coefficients (in wind axes, relative to the first Airplane's CG)). It is called
        only on a cache miss.
    :param persist: An optional callback that writes the whole cache through to disk. It
        is called only on a cache miss, and only when caching to a file. When None, the
        cache is kept in memory only. The default is None.
    :return: A tuple of a (num_airplanes, 6) ndarray of floats holding the mesh's per-
        Airplane load coefficients (the three force coefficients (in wind axes) followed
        by the three moment coefficients (in wind axes, relative to the first Airplane's
        CG)) and a float holding the time, in seconds, its solve took when it ran.
    """
    if key in cache:
        return cache[key]

    solve_start = time.time()
    coefficients = solve()
    solve_time = time.time() - solve_start
    cache[key] = (coefficients, solve_time)

    if persist is not None:
        persist()

    return cache[key]


def memos_to_disk(
    ref_problem_hash: str,
    panel_aspect_ratios_list: list[int],
    num_chordwise_panels_list: list[int],
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int],
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int],
    delta_time_cache: dict[tuple[int, int], float],
) -> dict[str, float]:
    """Translates the in-loop memo caches into a single dictionary keyed on absolute
    mesh values for storage on disk.

    The in-loop caches are keyed on bounds-relative sweep indices (the Panel aspect
    ratio index and the number of chordwise Panels index), which are meaningful only
    within the run that produced them. This converts each key's first two components to
    the absolute Panel aspect ratio and number of chordwise Panels, so a later run with
    different bounds can recover the memo by its absolute mesh rather than by an index
    that would point elsewhere. The remaining key components (the Airplane, Wing, and
    WingCrossSection indices) are structural indices into the reference geometry and are
    stable across runs, so they pass through unchanged. The absolute keys reuse the
    solve-cache key scheme, anchored to the reference problem's content hash, so one
    cache file can hold memos for multiple reference problems.

    :param ref_problem_hash: The content hash of the reference problem, as returned by
        the serialization module's hash_object. It anchors each memo key to this
        geometry, operating point, and motion.
    :param panel_aspect_ratios_list: The list of Panel aspect ratios being tested, used
        to convert each Panel aspect ratio index to its absolute value.
    :param num_chordwise_panels_list: The list of numbers of chordwise Panels being
        tested, used to convert each number of chordwise Panels index to its absolute
        value.
    :param num_spanwise_panels_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index, WingCrossSection
        index) tuple to a number of spanwise Panels, for trapezoidal Wings.
    :param num_wing_cross_sections_cache: The cache mapping a (Panel aspect ratio index,
        number of chordwise Panels index, Airplane index, Wing index) tuple to a number
        of WingCrossSections, for edge-defined Wings.
    :param delta_time_cache: The cache mapping a (Panel aspect ratio index, number of
        chordwise Panels index) tuple to an optimized delta_time.
    :return: A dict mapping each absolute-valued memo key to its value, ready to be
        written to the cache file's memo section.
    """
    memos: dict[str, float] = {}

    for (
        ar_id,
        chord_id,
        airplane_id,
        wing_id,
        section_id,
    ), num_spanwise_panels in num_spanwise_panels_cache.items():
        memos[
            solve_cache_key(
                ref_problem_hash,
                "spanwise",
                panel_aspect_ratios_list[ar_id],
                num_chordwise_panels_list[chord_id],
                airplane_id,
                wing_id,
                section_id,
            )
        ] = num_spanwise_panels

    for (
        ar_id,
        chord_id,
        airplane_id,
        wing_id,
    ), num_wing_cross_sections in num_wing_cross_sections_cache.items():
        memos[
            solve_cache_key(
                ref_problem_hash,
                "cross_sections",
                panel_aspect_ratios_list[ar_id],
                num_chordwise_panels_list[chord_id],
                airplane_id,
                wing_id,
            )
        ] = num_wing_cross_sections

    for (ar_id, chord_id), delta_time in delta_time_cache.items():
        memos[
            solve_cache_key(
                ref_problem_hash,
                "delta_time",
                panel_aspect_ratios_list[ar_id],
                num_chordwise_panels_list[chord_id],
            )
        ] = delta_time

    return memos


def memos_from_disk(
    memos: dict[str, float],
    ref_problem_hash: str,
    panel_aspect_ratios_list: list[int],
    num_chordwise_panels_list: list[int],
) -> tuple[
    dict[tuple[int, int, int, int, int], int],
    dict[tuple[int, int, int, int], int],
    dict[tuple[int, int], float],
]:
    """Translates a disk memo dictionary back into the three in-loop memo caches keyed
    on this run's bounds-relative sweep indices.

    This inverts memos_to_disk against this run's parameter lists. Each memo whose
    absolute Panel aspect ratio and number of chordwise Panels both fall within this
    run's bounds is mapped onto that run's Panel aspect ratio index and number of
    chordwise Panels index, so the caches pre-populate exactly as if this run had
    already resolved that mesh. Memos written for a different reference problem, or
    whose absolute mesh lies outside this run's bounds, are dropped. Spanwise and
    WingCrossSection counts are returned as ints and delta_times as floats, regardless
    of how JSON loaded them.

    :param memos: The disk memo dictionary mapping each absolute-valued memo key to its
        value, as returned by load_memo_cache.
    :param ref_problem_hash: The content hash of the reference problem, as returned by
        the serialization module's hash_object. Only memos anchored to this hash are
        kept.
    :param panel_aspect_ratios_list: The list of Panel aspect ratios being tested, used
        to convert each absolute Panel aspect ratio to its index and to reject out-of-
        bounds memos.
    :param num_chordwise_panels_list: The list of numbers of chordwise Panels being
        tested, used to convert each absolute number of chordwise Panels to its index
        and to reject out-of-bounds memos.
    :return: A tuple of the three in-loop memo caches, in order the number of spanwise
        Panels cache, the number of WingCrossSections cache, and the delta_time cache,
        each keyed on this run's sweep indices.
    """
    num_spanwise_panels_cache: dict[tuple[int, int, int, int, int], int] = {}
    num_wing_cross_sections_cache: dict[tuple[int, int, int, int], int] = {}
    delta_time_cache: dict[tuple[int, int], float] = {}

    for key, value in memos.items():
        parts = key.split("|")

        # Skip memos written for a different reference problem sharing this cache file.
        if parts[0] != ref_problem_hash:
            continue

        kind = parts[1]
        panel_aspect_ratio = int(parts[2])
        num_chordwise_panels = int(parts[3])

        # Skip memos whose absolute mesh falls outside this run's bounds, where the
        # sweep index would be meaningless.
        if (
            panel_aspect_ratio not in panel_aspect_ratios_list
            or num_chordwise_panels not in num_chordwise_panels_list
        ):
            continue

        ar_id = panel_aspect_ratios_list.index(panel_aspect_ratio)
        chord_id = num_chordwise_panels_list.index(num_chordwise_panels)

        if kind == "spanwise":
            num_spanwise_panels_cache[
                (ar_id, chord_id, int(parts[4]), int(parts[5]), int(parts[6]))
            ] = int(value)
        elif kind == "cross_sections":
            num_wing_cross_sections_cache[
                (ar_id, chord_id, int(parts[4]), int(parts[5]))
            ] = int(value)
        elif kind == "delta_time":
            delta_time_cache[(ar_id, chord_id)] = float(value)

    return (
        num_spanwise_panels_cache,
        num_wing_cross_sections_cache,
        delta_time_cache,
    )


def load_memo_cache(cache_path: Path | None) -> dict[str, float]:
    """Loads the memo section of a JSON cache file into a dictionary keyed on absolute
    mesh values.

    The memo cache is a pure optimization, so any file that cannot be read as a current-
    schema cache yields an empty dictionary rather than an error and the sweep resolves
    each mesh as usual. This covers a None path, a missing file, an unreadable file (for
    example a directory or a permission error), a file that is not valid JSON, a file
    whose schema version does not match the current one, and a current-version file
    whose memo section is missing or malformed. An unreadable or invalid file is logged
    as a warning so a silently degraded cache is still observable.

    :param cache_path: The path to the JSON cache file, or None to skip caching and
        return an empty dictionary.
    :return: A dict mapping each absolute-valued memo key to its value, as stored in the
        cache file's memo section.
    """
    if cache_path is None or not cache_path.exists():
        return {}

    try:
        with open(cache_path) as cache_file:
            data = json.load(cache_file)
    except (OSError, json.JSONDecodeError) as error:
        _logger.warning(
            "Ignoring unreadable memo cache at %s (%s); resolving each mesh.",
            cache_path,
            error,
        )
        return {}

    if data.get("_cache_version") != _SOLVE_CACHE_VERSION:
        return {}

    memos: dict[str, float] = data.get("memos", {})
    if not isinstance(memos, dict):
        _logger.warning(
            "Ignoring cache at %s with a malformed memo section; resolving each mesh.",
            cache_path,
        )
        return {}

    return memos
