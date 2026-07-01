"""This module contains classes to test the convergence analysis functions."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import numpy.testing as npt

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import convergence
from tests.unit.fixtures import geometry_fixtures


class TestConvergedParameterId(unittest.TestCase):
    """This class contains methods for testing
    convergence._converged_parameter_id.
    """

    def test_single_returns_this_id(self) -> None:
        """Test that a single tested value returns this iteration's own index."""
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=0, single=True, converged=False
            ),
            0,
        )

    def test_single_takes_precedence_over_converged(self) -> None:
        """Test that a single tested value returns this iteration's own index even
        when the converged flag is set.
        """
        self.assertEqual(
            convergence._converged_parameter_id(this_id=3, single=True, converged=True),
            3,
        )

    def test_converged_returns_coarser_id(self) -> None:
        """Test that a converged iteration returns the incrementally coarser
        index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=4, single=False, converged=True
            ),
            3,
        )

    def test_saturated_returns_this_id(self) -> None:
        """Test that an iteration that passed without converging returns this
        iteration's own index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=5, single=False, converged=False
            ),
            5,
        )

    def test_converged_id_depends_on_this_id(self) -> None:
        """Test that the coarser index tracks this iteration's index rather than a
        fixed value.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=7, single=False, converged=True
            ),
            6,
        )

    def test_returns_int(self) -> None:
        """Test that the converged index is returned as an int."""
        result = convergence._converged_parameter_id(
            this_id=2, single=False, converged=True
        )
        self.assertIsInstance(result, int)


class TestGetWingSectionNumSpanwisePanels(unittest.TestCase):
    """This class contains methods for testing
    convergence._get_wing_section_num_spanwise_panels, the search that picks a
    non-edge-defined Wing section's number of spanwise Panels.
    """

    def setUp(self) -> None:
        """Set up a root and tip WingCrossSection defining one wing section.

        :return: None
        """
        self.root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        self.tip_wing_cross_section = (
            geometry_fixtures.make_tip_wing_cross_section_fixture()
        )
        self.num_chordwise_panels = 4
        self.chordwise_spacing = "uniform"

    def _average_panel_aspect_ratio(self, num_spanwise_panels) -> float:
        """Meshes the wing section at a number of spanwise Panels and returns its average
        Panel aspect ratio.
        """
        return convergence._get_wing_section_average_panel_aspect_ratio(
            self.num_chordwise_panels,
            self.chordwise_spacing,
            self.root_wing_cross_section,
            self.tip_wing_cross_section,
            num_spanwise_panels,
        )

    def _num_spanwise_panels_for(self, target) -> int:
        """Searches for the number of spanwise Panels that hits a target average Panel
        aspect ratio, starting from the smallest valid count.
        """
        return convergence._get_wing_section_num_spanwise_panels(
            desired_average_panel_aspect_ratio=target,
            num_chordwise_panels=self.num_chordwise_panels,
            chordwise_spacing=self.chordwise_spacing,
            ref_root_wing_cross_section=self.root_wing_cross_section,
            ref_tip_wing_cross_section=self.tip_wing_cross_section,
            start_val=1,
        )

    def test_returns_int(self) -> None:
        """Test that the number of spanwise Panels is returned as an int."""
        self.assertIsInstance(self._num_spanwise_panels_for(4), int)

    def test_result_is_closest_to_target(self) -> None:
        """Test that the chosen number of spanwise Panels gives an average Panel aspect
        ratio at least as close to the target as either neighboring count.
        """
        target = 4
        result = self._num_spanwise_panels_for(target)

        chosen_difference = abs(self._average_panel_aspect_ratio(result) - target)
        upper_difference = abs(self._average_panel_aspect_ratio(result + 1) - target)
        self.assertLessEqual(chosen_difference, upper_difference)

        # The smallest valid number of spanwise Panels is one, so the lower neighbor only
        # exists above that.
        if result > 1:
            lower_difference = abs(
                self._average_panel_aspect_ratio(result - 1) - target
            )
            self.assertLessEqual(chosen_difference, lower_difference)

    def test_smaller_target_needs_at_least_as_many_panels(self) -> None:
        """Test that a smaller target Panel aspect ratio, a finer mesh, needs at least as
        many spanwise Panels as a larger one.
        """
        self.assertGreaterEqual(
            self._num_spanwise_panels_for(2), self._num_spanwise_panels_for(4)
        )


class TestGetNumWingCrossSectionsForPanelAr(unittest.TestCase):
    """This class contains methods for testing
    convergence._get_num_wing_cross_sections_for_panel_ar, the search that picks an
    edge-defined Wing's number of WingCrossSections.
    """

    @staticmethod
    def _tapered_edge_points() -> tuple[np.ndarray, np.ndarray]:
        """Builds straight, tapered leading and trailing edge curves.

        The leading edge sweeps back from the origin to (0.5, 2.0, 0.0) and the trailing
        edge is unswept at x = 1.0, so the chord tapers from a unit root chord to a 0.5
        tip chord over a two meter half span.
        """
        ys = np.linspace(0.0, 2.0, 20)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.25 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        return leading, trailing

    def _make_edge_wing(self, symmetric=False) -> ps.geometry.wing.Wing:
        """Builds an edge-defined Wing from the tapered edge curves."""
        leading, trailing = self._tapered_edge_points()
        return ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            name="Edge Wing",
            symmetric=symmetric,
            symmetryNormal_G=(0.0, 1.0, 0.0) if symmetric else None,
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0) if symmetric else None,
            num_chordwise_panels=4,
        )

    def setUp(self) -> None:
        """Set up a reference edge-defined Wing.

        :return: None
        """
        self.ref_wing = self._make_edge_wing()
        self.num_chordwise_panels = 4

    def _average_panel_aspect_ratio(self, num_wing_cross_sections) -> float:
        """Rebuilds and meshes the edge-defined Wing at a number of WingCrossSections and
        returns its average Panel aspect ratio.
        """
        refined_wing = convergence._build_edge_defined_wing(
            self.ref_wing, self.num_chordwise_panels, num_wing_cross_sections
        )
        airplane = ps.geometry.airplane.Airplane(wings=[refined_wing])
        average_panel_aspect_ratio = airplane.wings[0].average_panel_aspect_ratio
        assert average_panel_aspect_ratio is not None
        return average_panel_aspect_ratio

    def _num_wing_cross_sections_for(self, target, ref_wing=None) -> int:
        """Searches for the number of WingCrossSections that hits a target average Panel
        aspect ratio, starting from the smallest valid count.
        """
        return convergence._get_num_wing_cross_sections_for_panel_ar(
            desired_average_panel_aspect_ratio=target,
            num_chordwise_panels=self.num_chordwise_panels,
            ref_wing=self.ref_wing if ref_wing is None else ref_wing,
            start_val=2,
        )

    def test_returns_int(self) -> None:
        """Test that the number of WingCrossSections is returned as an int."""
        self.assertIsInstance(self._num_wing_cross_sections_for(4), int)

    def test_result_is_closest_to_target(self) -> None:
        """Test that the chosen number of WingCrossSections gives an average Panel aspect
        ratio at least as close to the target as either neighboring count.
        """
        target = 4
        result = self._num_wing_cross_sections_for(target)

        chosen_difference = abs(self._average_panel_aspect_ratio(result) - target)
        upper_difference = abs(self._average_panel_aspect_ratio(result + 1) - target)
        self.assertLessEqual(chosen_difference, upper_difference)

        # An edge-defined Wing needs at least two WingCrossSections, so the lower neighbor
        # only exists above that.
        if result > 2:
            lower_difference = abs(
                self._average_panel_aspect_ratio(result - 1) - target
            )
            self.assertLessEqual(chosen_difference, lower_difference)

    def test_smaller_target_needs_at_least_as_many_wing_cross_sections(self) -> None:
        """Test that a smaller target Panel aspect ratio, a finer mesh, needs at least as
        many WingCrossSections as a larger one.
        """
        self.assertGreaterEqual(
            self._num_wing_cross_sections_for(2), self._num_wing_cross_sections_for(4)
        )

    def test_symmetric_matches_asymmetric(self) -> None:
        """Test that the count is measured on the half span, so a symmetric Wing and an
        asymmetric Wing built from the same half-span curves need the same count.
        """
        asymmetric_result = self._num_wing_cross_sections_for(
            4, ref_wing=self._make_edge_wing(symmetric=False)
        )
        symmetric_result = self._num_wing_cross_sections_for(
            4, ref_wing=self._make_edge_wing(symmetric=True)
        )
        self.assertEqual(asymmetric_result, symmetric_result)


class TestValidateCoefficientMask(unittest.TestCase):
    """This class contains methods for testing
    convergence._validate_coefficient_mask.
    """

    def test_none_returns_all_true(self) -> None:
        """Test that None returns a (6,) mask of all True."""
        result = convergence._validate_coefficient_mask(None)
        self.assertTrue(np.array_equal(result, np.ones(6, dtype=bool)))

    def test_tuple_is_returned_as_bool_array(self) -> None:
        """Test that a valid tuple is returned as an equivalent (6,) bool ndarray."""
        result = convergence._validate_coefficient_mask(
            (True, False, True, False, True, False)
        )
        self.assertTrue(
            np.array_equal(result, np.array([1, 0, 1, 0, 1, 0], dtype=bool))
        )

    def test_non_tuple_raises_type_error(self) -> None:
        """Test that a non-tuple, non-None mask raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_coefficient_mask([True] * 6)  # type: ignore[arg-type]

    def test_wrong_length_raises_value_error(self) -> None:
        """Test that a mask without exactly six elements raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_coefficient_mask(
                (True, True, True)  # type: ignore[arg-type]
            )

    def test_non_bool_element_raises_type_error(self) -> None:
        """Test that a mask with a non-bool element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_coefficient_mask(
                (True, 1, True, True, True, True)  # type: ignore[arg-type]
            )

    def test_all_false_raises_value_error(self) -> None:
        """Test that a mask with no True element raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_coefficient_mask((False,) * 6)


class TestCheckCoefficientConvergence(unittest.TestCase):
    """This class contains methods for testing
    convergence._check_coefficient_convergence.
    """

    def setUp(self) -> None:
        """Set up a full mask and the tolerances used across the tests.

        :return: None
        """
        self.mask = np.ones(6, dtype=bool)
        self.rtol = 0.05
        self.atol = 0.001

    def test_identical_coefficients_converge(self) -> None:
        """Test that identical coefficients converge with a perfect metric."""
        these = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        converged, metric, _ = convergence._check_coefficient_convergence(
            these, these.copy(), self.rtol, self.atol, self.mask
        )
        self.assertTrue(converged)
        self.assertEqual(metric, 100.0)

    def test_large_relative_change_does_not_converge(self) -> None:
        """Test that a coefficient changing by more than the relative tolerance does not
        converge and is reported as the limiting coefficient.
        """
        coarser = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        these = coarser.copy()
        these[0, 2] = 2.5
        converged, _, limiting_id = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertFalse(converged)
        self.assertEqual(limiting_id, 2)

    def test_masked_out_coefficient_is_ignored(self) -> None:
        """Test that masking out the only offending coefficient makes the check
        converge.
        """
        coarser = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        these = coarser.copy()
        these[0, 2] = 2.5
        mask = np.array([True, True, False, True, True, True], dtype=bool)
        converged, _, _ = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, mask
        )
        self.assertTrue(converged)

    def test_absolute_tolerance_floors_near_zero(self) -> None:
        """Test that a coefficient near zero converges via the absolute tolerance floor
        even though its relative change is large.
        """
        these = np.zeros((1, 6), dtype=float)
        coarser = np.zeros((1, 6), dtype=float)
        coarser[0, 0] = 0.5 * self.atol
        converged, _, _ = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertTrue(converged)

    def test_all_airplanes_must_converge(self) -> None:
        """Test that the check fails when any one Airplane has an unconverged
        coefficient.
        """
        coarser = np.array(
            [
                [1.0, 0.0, 2.0, 0.1, 0.0, 0.05],
                [1.0, 0.0, 2.0, 0.1, 0.0, 0.05],
            ]
        )
        these = coarser.copy()
        these[1, 0] = 2.0
        converged, _, limiting_id = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertFalse(converged)
        self.assertEqual(limiting_id, 0)

    def test_returns_bool_float_int(self) -> None:
        """Test that the result is a bool, a float, and an int."""
        these = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        converged, metric, limiting_id = convergence._check_coefficient_convergence(
            these, these.copy(), self.rtol, self.atol, self.mask
        )
        self.assertIsInstance(converged, bool)
        self.assertIsInstance(metric, float)
        self.assertIsInstance(limiting_id, int)


class TestSolveCacheKey(unittest.TestCase):
    """This class contains methods for testing convergence._solve_cache_key."""

    def test_returns_str(self) -> None:
        """Test that the key is returned as a string."""
        key = convergence._solve_cache_key("abc123", "steady ring", 4, 3)
        self.assertIsInstance(key, str)

    def test_deterministic(self) -> None:
        """Test that identical inputs produce identical keys."""
        first = convergence._solve_cache_key("abc123", "steady ring", 4, 3)
        second = convergence._solve_cache_key("abc123", "steady ring", 4, 3)
        self.assertEqual(first, second)

    def test_different_component_differs(self) -> None:
        """Test that changing one component changes the key."""
        first = convergence._solve_cache_key("abc123", "steady ring", 4, 3)
        second = convergence._solve_cache_key("abc123", "steady ring", 4, 4)
        self.assertNotEqual(first, second)

    def test_different_hash_differs(self) -> None:
        """Test that changing the reference problem hash changes the key."""
        first = convergence._solve_cache_key("abc123", "steady ring", 4, 3)
        second = convergence._solve_cache_key("def456", "steady ring", 4, 3)
        self.assertNotEqual(first, second)

    def test_component_order_matters(self) -> None:
        """Test that two components in a different order produce different keys."""
        first = convergence._solve_cache_key("abc123", 4, 3)
        second = convergence._solve_cache_key("abc123", 3, 4)
        self.assertNotEqual(first, second)


class TestSolveCache(unittest.TestCase):
    """This class contains methods for testing convergence._cached_solve,
    convergence._load_solve_cache, and convergence._write_cache.
    """

    def test_miss_calls_solve_and_returns_result(self) -> None:
        """Test that a cache miss calls solve and returns its coefficients."""
        cache: dict[str, np.ndarray] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        result = convergence._cached_solve(cache, "key", lambda: coefficients)
        npt.assert_array_equal(result, coefficients)

    def test_miss_stores_result_in_cache(self) -> None:
        """Test that a cache miss adds the coefficients to the in-memory cache."""
        cache: dict[str, np.ndarray] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        convergence._cached_solve(cache, "key", lambda: coefficients)
        self.assertIn("key", cache)
        npt.assert_array_equal(cache["key"], coefficients)

    def test_hit_returns_cached_without_solving(self) -> None:
        """Test that a cache hit returns the stored coefficients and skips solve."""
        cached = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        cache: dict[str, np.ndarray] = {"key": cached}
        solve_calls = 0

        def solve() -> np.ndarray:
            nonlocal solve_calls
            solve_calls += 1
            return np.zeros((1, 6), dtype=float)

        result = convergence._cached_solve(cache, "key", solve)
        npt.assert_array_equal(result, cached)
        self.assertEqual(solve_calls, 0)

    def test_no_persist_callback_does_not_write(self) -> None:
        """Test that a miss with no persist callback writes nothing to disk."""
        cache: dict[str, np.ndarray] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            convergence._cached_solve(cache, "key", lambda: coefficients)
            self.assertFalse(path.exists())

    def test_hit_does_not_persist(self) -> None:
        """Test that a cache hit does not invoke the persist callback."""
        cached = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        cache: dict[str, np.ndarray] = {"key": cached}
        persist_calls = 0

        def persist() -> None:
            nonlocal persist_calls
            persist_calls += 1

        convergence._cached_solve(
            cache, "key", lambda: np.zeros((1, 6), dtype=float), persist
        )
        self.assertEqual(persist_calls, 0)

    def test_write_through_persists_to_disk(self) -> None:
        """Test that a miss invokes the persist callback so the entry reloads."""
        cache: dict[str, np.ndarray] = {}
        coefficients = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            convergence._cached_solve(
                cache,
                "key",
                lambda: coefficients,
                lambda: convergence._write_cache(path, cache, {}),
            )
            reloaded = convergence._load_solve_cache(path)
        self.assertIn("key", reloaded)
        npt.assert_array_equal(reloaded["key"], coefficients)

    def test_load_missing_file_returns_empty(self) -> None:
        """Test that loading a missing cache file returns an empty cache."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "missing.json"
            self.assertEqual(convergence._load_solve_cache(path), {})

    def test_load_none_path_returns_empty(self) -> None:
        """Test that loading a None cache path returns an empty cache."""
        self.assertEqual(convergence._load_solve_cache(None), {})

    def test_load_version_mismatch_returns_empty(self) -> None:
        """Test that a cache file with a mismatched schema version is ignored."""
        cache = {"key": np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            convergence._write_cache(path, cache, {})
            with open(path) as cache_file:
                data = json.load(cache_file)
            data["_cache_version"] = convergence._SOLVE_CACHE_VERSION + 1
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(convergence._load_solve_cache(path), {})

    def test_round_trip_preserves_coefficients(self) -> None:
        """Test that writing then loading a cache preserves its coefficients."""
        cache = {
            "first": np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float),
            "second": np.array(
                [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]],
                dtype=float,
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            convergence._write_cache(path, cache, {})
            reloaded = convergence._load_solve_cache(path)
        self.assertEqual(set(reloaded), set(cache))
        for key in cache:
            npt.assert_array_equal(reloaded[key], cache[key])


class TestMemoTranslation(unittest.TestCase):
    """This class contains methods for testing convergence._memos_to_disk and
    convergence._memos_from_disk, which translate the in-loop memo caches between their
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

        memos = convergence._memos_to_disk(
            "abc123",
            ar_list,
            chord_list,
            num_spanwise_panels_cache,
            num_wing_cross_sections_cache,
            delta_time_cache,
        )

        self.assertEqual(
            memos[convergence._solve_cache_key("abc123", "spanwise", 3, 3, 0, 0, 0)], 5
        )
        self.assertEqual(
            memos[convergence._solve_cache_key("abc123", "cross_sections", 3, 3, 0, 0)],
            7,
        )
        self.assertEqual(
            memos[convergence._solve_cache_key("abc123", "delta_time", 3, 3)], 0.0125
        )

    def test_from_disk_maps_absolute_to_new_indices(self) -> None:
        """Test that translating from disk maps each in-bounds memo onto the sweep
        indices of a run whose bounds differ from the run that wrote it.

        :return: None
        """
        memos = {
            convergence._solve_cache_key("abc123", "spanwise", 3, 3, 0, 0, 0): 5,
            convergence._solve_cache_key("abc123", "cross_sections", 3, 3, 0, 0): 7,
            convergence._solve_cache_key("abc123", "delta_time", 3, 3): 0.0125,
        }

        # A later run over Panel aspect ratios (3, 2) and a single chordwise value of 3.
        ar_list = [3, 2]
        chord_list = [3]

        spanwise, cross_sections, delta_time = convergence._memos_from_disk(
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
            convergence._solve_cache_key("abc123", "delta_time", 4, 3): 0.02,
            convergence._solve_cache_key("abc123", "delta_time", 3, 3): 0.0125,
        }

        ar_list = [3, 2]
        chord_list = [3]

        _, _, delta_time = convergence._memos_from_disk(
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
            convergence._solve_cache_key("abc123", "delta_time", 3, 3): 0.0125,
            convergence._solve_cache_key("def456", "delta_time", 3, 3): 0.5,
        }

        ar_list = [3]
        chord_list = [3]

        _, _, delta_time = convergence._memos_from_disk(
            memos, "abc123", ar_list, chord_list
        )

        self.assertEqual(delta_time, {(0, 0): 0.0125})

    def test_from_disk_preserves_value_types(self) -> None:
        """Test that spanwise and cross-section counts come back as ints and delta_time
        as a float, even when JSON loading widened an integral delta_time.

        :return: None
        """
        memos: dict[str, float] = {
            convergence._solve_cache_key("abc123", "spanwise", 3, 3, 0, 0, 0): 5,
            convergence._solve_cache_key("abc123", "delta_time", 3, 3): 1,
        }

        ar_list = [3]
        chord_list = [3]

        spanwise, _, delta_time = convergence._memos_from_disk(
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

        memos = convergence._memos_to_disk(
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

        spanwise, cross_sections, delta_time = convergence._memos_from_disk(
            memos, "abc123", read_ar_list, read_chord_list
        )

        # Only the Panel aspect ratio 3, 4-chordwise entries fall inside the new bounds,
        # where Panel aspect ratio 3 is index 0 and 4 chordwise Panels is index 0.
        self.assertEqual(spanwise, {(0, 0, 0, 0, 0): 6})
        # The cross-section memo sat at 3 chordwise Panels, which is out of bounds now.
        self.assertEqual(cross_sections, {})
        self.assertEqual(delta_time, {(0, 0): 0.0125})


class TestMemoCacheDisk(unittest.TestCase):
    """This class contains methods for testing convergence._load_memo_cache, which reads
    the memo section of a JSON cache file.
    """

    def test_load_missing_file_returns_empty(self) -> None:
        """Test that loading memos from a missing cache file returns an empty dict.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "missing.json"
            self.assertEqual(convergence._load_memo_cache(path), {})

    def test_load_none_path_returns_empty(self) -> None:
        """Test that loading memos from a None cache path returns an empty dict.

        :return: None
        """
        self.assertEqual(convergence._load_memo_cache(None), {})

    def test_load_version_mismatch_returns_empty(self) -> None:
        """Test that memos in a cache file with a mismatched schema version are ignored.

        :return: None
        """
        key = convergence._solve_cache_key("abc123", "delta_time", 3, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": convergence._SOLVE_CACHE_VERSION + 1,
                "entries": {},
                "memos": {key: 0.0125},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(convergence._load_memo_cache(path), {})

    def test_load_reads_memos_section(self) -> None:
        """Test that loading returns the stored memo section keyed by absolute value.

        :return: None
        """
        key = convergence._solve_cache_key("abc123", "delta_time", 3, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": convergence._SOLVE_CACHE_VERSION,
                "entries": {},
                "memos": {key: 0.0125},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(convergence._load_memo_cache(path), {key: 0.0125})

    def test_load_missing_memos_section_returns_empty(self) -> None:
        """Test that a current-version cache file with no memo section loads as an empty
        memo dict rather than raising.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            data = {
                "_cache_version": convergence._SOLVE_CACHE_VERSION,
                "entries": {},
            }
            with open(path, "w") as cache_file:
                json.dump(data, cache_file)
            self.assertEqual(convergence._load_memo_cache(path), {})

    def test_write_cache_round_trips_both_sections(self) -> None:
        """Test that the unified writer stores both the solve entries and the memos so
        each section reloads independently without clobbering the other.

        :return: None
        """
        solve_cache = {"solve": np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]], dtype=float)}
        memo_key = convergence._solve_cache_key("abc123", "delta_time", 3, 3)
        memo_cache = {memo_key: 0.0125}

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cache.json"
            convergence._write_cache(path, solve_cache, memo_cache)
            reloaded_solve = convergence._load_solve_cache(path)
            reloaded_memos = convergence._load_memo_cache(path)

        self.assertEqual(set(reloaded_solve), {"solve"})
        npt.assert_array_equal(reloaded_solve["solve"], solve_cache["solve"])
        self.assertEqual(reloaded_memos, memo_cache)
