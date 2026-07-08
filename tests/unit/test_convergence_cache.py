"""This module contains classes to test the convergence cache functions."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _convergence_cache


class TestSolveCacheKey(unittest.TestCase):
    """This class contains methods for testing _convergence_cache.solve_cache_key."""

    def test_returns_str(self) -> None:
        """Test that the key is returned as a string."""
        key = _convergence_cache.solve_cache_key("abc123", "steady ring", 4, 3)
        self.assertIsInstance(key, str)

    def test_deterministic(self) -> None:
        """Test that identical inputs produce identical keys."""
        first = _convergence_cache.solve_cache_key("abc123", "steady ring", 4, 3)
        second = _convergence_cache.solve_cache_key("abc123", "steady ring", 4, 3)
        self.assertEqual(first, second)

    def test_different_component_differs(self) -> None:
        """Test that changing one component changes the key."""
        first = _convergence_cache.solve_cache_key("abc123", "steady ring", 4, 3)
        second = _convergence_cache.solve_cache_key("abc123", "steady ring", 4, 4)
        self.assertNotEqual(first, second)

    def test_different_hash_differs(self) -> None:
        """Test that changing the reference problem hash changes the key."""
        first = _convergence_cache.solve_cache_key("abc123", "steady ring", 4, 3)
        second = _convergence_cache.solve_cache_key("def456", "steady ring", 4, 3)
        self.assertNotEqual(first, second)

    def test_component_order_matters(self) -> None:
        """Test that two components in a different order produce different keys."""
        first = _convergence_cache.solve_cache_key("abc123", 4, 3)
        second = _convergence_cache.solve_cache_key("abc123", 3, 4)
        self.assertNotEqual(first, second)


class TestSolveCache(unittest.TestCase):
    """This class contains methods for testing _convergence_cache.cached_solve,
    _convergence_cache.load_solve_cache, and _convergence_cache.write_cache.
    """

    def test_miss_calls_solve_and_returns_result(self) -> None:
        """Test that a cache miss calls solve and returns its coefficients and a
        non-negative solve time.
        """
        cache: dict[str, tuple[np.ndarray, float]] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        result_coefficients, result_solve_time = _convergence_cache.cached_solve(
            cache, "key", lambda: coefficients
        )
        npt.assert_array_equal(result_coefficients, coefficients)
        self.assertIsInstance(result_solve_time, float)
        self.assertGreaterEqual(result_solve_time, 0.0)

    def test_miss_stores_result_in_cache(self) -> None:
        """Test that a cache miss adds the coefficients and solve time to the in-memory
        cache.
        """
        cache: dict[str, tuple[np.ndarray, float]] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        result = _convergence_cache.cached_solve(cache, "key", lambda: coefficients)
        self.assertIn("key", cache)
        npt.assert_array_equal(cache["key"][0], coefficients)
        self.assertEqual(cache["key"][1], result[1])

    def test_hit_returns_cached_without_solving(self) -> None:
        """Test that a cache hit returns the stored coefficients and solve time and
        skips solve.
        """
        cached = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        cache: dict[str, tuple[np.ndarray, float]] = {"key": (cached, 42.5)}
        solve_calls = 0

        def solve() -> np.ndarray:
            nonlocal solve_calls
            solve_calls += 1
            return np.zeros((1, 6), dtype=float)

        result_coefficients, result_solve_time = _convergence_cache.cached_solve(
            cache, "key", solve
        )
        npt.assert_array_equal(result_coefficients, cached)
        self.assertEqual(result_solve_time, 42.5)
        self.assertEqual(solve_calls, 0)

    def test_no_persist_callback_does_not_write(self) -> None:
        """Test that a miss with no persist callback writes nothing to disk."""
        cache: dict[str, tuple[np.ndarray, float]] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            _convergence_cache.cached_solve(cache, "key", lambda: coefficients)
            self.assertFalse(path.exists())

    def test_hit_does_not_persist(self) -> None:
        """Test that a cache hit does not invoke the persist callback."""
        cached = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        cache: dict[str, tuple[np.ndarray, float]] = {"key": (cached, 42.5)}
        persist_calls = 0

        def persist() -> None:
            nonlocal persist_calls
            persist_calls += 1

        _convergence_cache.cached_solve(
            cache, "key", lambda: np.zeros((1, 6), dtype=float), persist
        )
        self.assertEqual(persist_calls, 0)

    def test_write_through_persists_to_disk(self) -> None:
        """Test that a miss invokes the persist callback so the entry reloads."""
        cache: dict[str, tuple[np.ndarray, float]] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            _convergence_cache.cached_solve(
                cache,
                "key",
                lambda: coefficients,
                lambda: _convergence_cache.write_cache(path, cache, {}),
            )
            reloaded = _convergence_cache.load_solve_cache(path)
        self.assertIn("key", reloaded)
        npt.assert_array_equal(reloaded["key"][0], coefficients)
        self.assertEqual(reloaded["key"][1], cache["key"][1])

    def test_load_missing_file_returns_empty(self) -> None:
        """Test that loading a missing cache file returns an empty cache."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "missing.json"
            self.assertEqual(_convergence_cache.load_solve_cache(path), {})

    def test_load_none_path_returns_empty(self) -> None:
        """Test that loading a None cache path returns an empty cache."""
        self.assertEqual(_convergence_cache.load_solve_cache(None), {})

    def test_load_version_mismatch_returns_empty(self) -> None:
        """Test that a cache file with a mismatched schema version is ignored."""
        cache = {"key": (np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float), 42.5)}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            _convergence_cache.write_cache(path, cache, {})
            with open(path) as cache_file:
                data = json.load(cache_file)
            data["_cache_version"] = _convergence_cache._SOLVE_CACHE_VERSION + 1
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(_convergence_cache.load_solve_cache(path), {})

    def test_round_trip_preserves_coefficients_and_solve_times(self) -> None:
        """Test that writing then loading a cache preserves its coefficients and solve
        times.
        """
        cache = {
            "first": (np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float), 1.25),
            "second": (
                np.array(
                    [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]],
                    dtype=float,
                ),
                600.0,
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            _convergence_cache.write_cache(path, cache, {})
            reloaded = _convergence_cache.load_solve_cache(path)
        self.assertEqual(set(reloaded), set(cache))
        for key in cache:
            npt.assert_array_equal(reloaded[key][0], cache[key][0])
            self.assertEqual(reloaded[key][1], cache[key][1])

    def test_load_invalid_json_returns_empty(self) -> None:
        """Test that a cache file that is not valid JSON is ignored with a warning."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            with open(path, "w") as cache_file:
                cache_file.write("{not valid json")
            with self.assertLogs("pterasoftware._convergence_cache", level="WARNING"):
                self.assertEqual(_convergence_cache.load_solve_cache(path), {})

    def test_load_missing_entries_section_returns_empty(self) -> None:
        """Test that a current-version file with no entries section loads as empty."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                "memos": {},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(_convergence_cache.load_solve_cache(path), {})

    def test_load_malformed_entries_section_returns_empty(self) -> None:
        """Test that a file whose entries section is not a dict is ignored with a
        warning.
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                "entries": [1, 2, 3],
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            with self.assertLogs("pterasoftware._convergence_cache", level="WARNING"):
                self.assertEqual(_convergence_cache.load_solve_cache(path), {})

    def test_load_malformed_entry_returns_empty(self) -> None:
        """Test that a current-version file whose entries section holds a malformed
        entry is ignored with a warning.
        """
        # The variants cover an entry that is not a dict, an entry missing its
        # coefficients, an entry missing its solve time, and an entry with non-numeric
        # coefficients.
        malformed_entries: tuple[dict[str, object], ...] = (
            {"key": [1.0, 2.0, 3.0]},
            {"key": {"solve_time": 42.5}},
            {"key": {"coefficients": [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]}},
            {"key": {"coefficients": [["junk"] * 6], "solve_time": 42.5}},
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            for entries in malformed_entries:
                with self.subTest(entries=entries):
                    data = {
                        "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                        "entries": entries,
                    }
                    with open(path, "w") as cache_file:
                        json.dump(data, cache_file)
                    with self.assertLogs(
                        "pterasoftware._convergence_cache", level="WARNING"
                    ):
                        self.assertEqual(_convergence_cache.load_solve_cache(path), {})


class TestMemoTranslation(unittest.TestCase):
    """This class contains methods for testing _convergence_cache.memos_to_disk and
    _convergence_cache.memos_from_disk, which translate the in-loop memo caches between their
    in-memory bounds-relative index keys and the absolute-valued keys stored on disk.
    """

    def test_to_disk_keys_on_absolute_values(self) -> None:
        """Test that translating to disk keys each memo on its absolute Panel aspect
        ratio and number of chordwise Panels rather than its sweep index.

        :return: None
        """
        ar_list = [4, 3, 2, 1]
        chord_list = [3, 4]

        # The sweep index (1, 0) resolves to Panel aspect ratio 3 and 3 chordwise Panels.
        num_spanwise_panels_cache = {(1, 0, 0, 0, 0): 5}
        num_wing_cross_sections_cache = {(1, 0, 0, 0): 7}
        delta_time_cache = {(1, 0): 0.0125}

        memos = _convergence_cache.memos_to_disk(
            "abc123",
            ar_list,
            chord_list,
            num_spanwise_panels_cache,
            num_wing_cross_sections_cache,
            delta_time_cache,
        )

        self.assertEqual(
            memos[
                _convergence_cache.solve_cache_key("abc123", "spanwise", 3, 3, 0, 0, 0)
            ],
            5,
        )
        self.assertEqual(
            memos[
                _convergence_cache.solve_cache_key(
                    "abc123", "cross_sections", 3, 3, 0, 0
                )
            ],
            7,
        )
        self.assertEqual(
            memos[_convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3)],
            0.0125,
        )

    def test_from_disk_maps_absolute_to_new_indices(self) -> None:
        """Test that translating from disk maps each in-bounds memo onto the sweep
        indices of a run whose bounds differ from the run that wrote it.

        :return: None
        """
        memos = {
            _convergence_cache.solve_cache_key("abc123", "spanwise", 3, 3, 0, 0, 0): 5,
            _convergence_cache.solve_cache_key(
                "abc123", "cross_sections", 3, 3, 0, 0
            ): 7,
            _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3): 0.0125,
        }

        # A later run over Panel aspect ratios (3, 2) and a single chordwise value of 3.
        ar_list = [3, 2]
        chord_list = [3]

        spanwise, cross_sections, delta_time = _convergence_cache.memos_from_disk(
            memos, "abc123", ar_list, chord_list
        )

        # Panel aspect ratio 3 is index 0 and 3 chordwise Panels is index 0 in this run.
        self.assertEqual(spanwise, {(0, 0, 0, 0, 0): 5})
        self.assertEqual(cross_sections, {(0, 0, 0, 0): 7})
        self.assertEqual(delta_time, {(0, 0): 0.0125})

    def test_from_disk_drops_out_of_bounds(self) -> None:
        """Test that translating from disk drops any memo whose absolute Panel aspect
        ratio or number of chordwise Panels lies outside the new run's bounds.

        :return: None
        """
        memos = {
            _convergence_cache.solve_cache_key("abc123", "delta_time", 4, 3): 0.02,
            _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3): 0.0125,
        }

        ar_list = [3, 2]
        chord_list = [3]

        _, _, delta_time = _convergence_cache.memos_from_disk(
            memos, "abc123", ar_list, chord_list
        )

        # Panel aspect ratio 4 is outside the (3, 2) bounds, so only 3 survives.
        self.assertEqual(delta_time, {(0, 0): 0.0125})

    def test_from_disk_drops_other_reference_problems(self) -> None:
        """Test that translating from disk ignores memos written for a different
        reference problem sharing the same cache file.

        :return: None
        """
        memos = {
            _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3): 0.0125,
            _convergence_cache.solve_cache_key("def456", "delta_time", 3, 3): 0.5,
        }

        ar_list = [3]
        chord_list = [3]

        _, _, delta_time = _convergence_cache.memos_from_disk(
            memos, "abc123", ar_list, chord_list
        )

        self.assertEqual(delta_time, {(0, 0): 0.0125})

    def test_from_disk_preserves_value_types(self) -> None:
        """Test that spanwise and cross-section counts come back as ints and delta_time
        as a float, even when JSON loading widened an integral delta_time.

        :return: None
        """
        memos: dict[str, float] = {
            _convergence_cache.solve_cache_key("abc123", "spanwise", 3, 3, 0, 0, 0): 5,
            _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3): 1,
        }

        ar_list = [3]
        chord_list = [3]

        spanwise, _, delta_time = _convergence_cache.memos_from_disk(
            memos, "abc123", ar_list, chord_list
        )

        self.assertIsInstance(spanwise[(0, 0, 0, 0, 0)], int)
        self.assertIsInstance(delta_time[(0, 0)], float)

    def test_round_trip_across_different_bounds(self) -> None:
        """Test that memos written by one run are recovered at the correct indices by a
        later run whose bounds overlap only partly.

        :return: None
        """
        write_ar_list = [4, 3, 2, 1]
        write_chord_list = [3, 4]

        # Index (0, 0) is Panel aspect ratio 4 with 3 chordwise Panels, and index (1, 1)
        # is Panel aspect ratio 3 with 4 chordwise Panels.
        num_spanwise_panels_cache = {
            (0, 0, 0, 0, 0): 4,
            (1, 1, 0, 0, 0): 6,
        }
        # Index (1, 0) is Panel aspect ratio 3 with 3 chordwise Panels.
        num_wing_cross_sections_cache = {(1, 0, 0, 0): 7}
        # Index (1, 1) is Panel aspect ratio 3 with 4 chordwise Panels.
        delta_time_cache = {(1, 1): 0.0125}

        memos = _convergence_cache.memos_to_disk(
            "abc123",
            write_ar_list,
            write_chord_list,
            num_spanwise_panels_cache,
            num_wing_cross_sections_cache,
            delta_time_cache,
        )

        # A later run over Panel aspect ratios (3, 2) and a single chordwise value of 4.
        read_ar_list = [3, 2]
        read_chord_list = [4]

        spanwise, cross_sections, delta_time = _convergence_cache.memos_from_disk(
            memos, "abc123", read_ar_list, read_chord_list
        )

        # Only the Panel aspect ratio 3, 4-chordwise entries fall inside the new bounds,
        # where Panel aspect ratio 3 is index 0 and 4 chordwise Panels is index 0.
        self.assertEqual(spanwise, {(0, 0, 0, 0, 0): 6})
        # The cross-section memo sat at 3 chordwise Panels, which is out of bounds now.
        self.assertEqual(cross_sections, {})
        self.assertEqual(delta_time, {(0, 0): 0.0125})


class TestMemoCacheDisk(unittest.TestCase):
    """This class contains methods for testing _convergence_cache.load_memo_cache, which reads
    the memo section of a JSON cache file.
    """

    def test_load_missing_file_returns_empty(self) -> None:
        """Test that loading memos from a missing cache file returns an empty dict.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "missing.json"
            self.assertEqual(_convergence_cache.load_memo_cache(path), {})

    def test_load_none_path_returns_empty(self) -> None:
        """Test that loading memos from a None cache path returns an empty dict.

        :return: None
        """
        self.assertEqual(_convergence_cache.load_memo_cache(None), {})

    def test_load_version_mismatch_returns_empty(self) -> None:
        """Test that memos in a cache file with a mismatched schema version are ignored.

        :return: None
        """
        key = _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION + 1,
                "entries": {},
                "memos": {key: 0.0125},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(_convergence_cache.load_memo_cache(path), {})

    def test_load_reads_memos_section(self) -> None:
        """Test that loading returns the stored memo section keyed by absolute value.

        :return: None
        """
        key = _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                "entries": {},
                "memos": {key: 0.0125},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(_convergence_cache.load_memo_cache(path), {key: 0.0125})

    def test_load_missing_memos_section_returns_empty(self) -> None:
        """Test that a current-version cache file with no memo section loads as an empty
        memo dict rather than raising.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                "entries": {},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(_convergence_cache.load_memo_cache(path), {})

    def test_load_invalid_json_returns_empty(self) -> None:
        """Test that a memo cache file that is not valid JSON is ignored with a warning.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            with open(path, "w") as cache_file:
                cache_file.write("{not valid json")
            with self.assertLogs("pterasoftware._convergence_cache", level="WARNING"):
                self.assertEqual(_convergence_cache.load_memo_cache(path), {})

    def test_load_malformed_memos_section_returns_empty(self) -> None:
        """Test that a current-version file whose memo section is not a dict is ignored
        with a warning.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                "entries": {},
                "memos": [1, 2, 3],
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            with self.assertLogs("pterasoftware._convergence_cache", level="WARNING"):
                self.assertEqual(_convergence_cache.load_memo_cache(path), {})

    def test_load_non_numeric_memo_value_returns_empty(self) -> None:
        """Test that a current-version file whose memo section holds a non-numeric
        value is ignored with a warning.

        :return: None
        """
        key = _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            for value in ("junk", None, [0.0125]):
                with self.subTest(value=value):
                    data = {
                        "_cache_version": _convergence_cache._SOLVE_CACHE_VERSION,
                        "entries": {},
                        "memos": {key: value},
                    }
                    with open(path, "w") as cache_file:
                        json.dump(data, cache_file)
                    with self.assertLogs(
                        "pterasoftware._convergence_cache", level="WARNING"
                    ):
                        self.assertEqual(_convergence_cache.load_memo_cache(path), {})

    def test_write_cache_round_trips_both_sections(self) -> None:
        """Test that the unified writer stores both the solve entries and the memos so
        each section reloads independently without clobbering the other.

        :return: None
        """
        solve_cache = {
            "solve": (np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float), 42.5)
        }
        memo_key = _convergence_cache.solve_cache_key("abc123", "delta_time", 3, 3)
        memo_cache = {memo_key: 0.0125}

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            _convergence_cache.write_cache(path, solve_cache, memo_cache)
            reloaded_solve = _convergence_cache.load_solve_cache(path)
            reloaded_memos = _convergence_cache.load_memo_cache(path)

        self.assertEqual(set(reloaded_solve), {"solve"})
        npt.assert_array_equal(reloaded_solve["solve"][0], solve_cache["solve"][0])
        self.assertEqual(reloaded_solve["solve"][1], solve_cache["solve"][1])
        self.assertEqual(reloaded_memos, memo_cache)
