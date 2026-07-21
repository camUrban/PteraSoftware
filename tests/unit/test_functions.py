"""This module contains classes to test shared utility functions."""

import threading
import unittest

import numpy as np
import numpy.testing as npt
import threadpoolctl

# noinspection PyProtectedMember
from pterasoftware import _functions


class TestCosspace(unittest.TestCase):
    """Tests for the cosspace function."""

    def test_includes_endpoints(self) -> None:
        """Cosspace should include minimum and maximum by default."""
        points = _functions.cosspace(0.0, 10.0, n_points=5)

        self.assertAlmostEqual(points[0], 0.0)
        self.assertAlmostEqual(points[-1], 10.0)

    def test_excludes_endpoint_when_requested(self) -> None:
        """Cosspace should exclude maximum when endpoint=False."""
        points = _functions.cosspace(
            0.0,
            10.0,
            n_points=5,
            endpoint=False,
        )

        self.assertAlmostEqual(points[0], 0.0)
        self.assertLess(points[-1], 10.0)

    def test_returns_correct_number_of_points(self) -> None:
        """Cosspace should return the requested number of points."""
        points = _functions.cosspace(0.0, 1.0, n_points=17)

        self.assertEqual(len(points), 17)

    def test_spacing_is_symmetric(self) -> None:
        """Cosspace output should be symmetric about the midpoint."""
        points = _functions.cosspace(-1.0, 1.0, n_points=9)

        npt.assert_allclose(points, -points[::-1], atol=1e-14)

    def test_single_point(self) -> None:
        """Single-point cosspace should return the minimum value."""
        points = _functions.cosspace(2.0, 8.0, n_points=1)

        npt.assert_allclose(points, np.array([2.0]), atol=1e-14)


class TestNumbaCentroidOfQuadrilateral(unittest.TestCase):
    """Tests for the numba_centroid_of_quadrilateral function."""

    def test_unit_square_in_xy_plane(self) -> None:
        """Centroid of a unit square should be at its geometric center."""
        front_left = np.array([0.0, 1.0, 0.0])
        front_right = np.array([1.0, 1.0, 0.0])
        back_left = np.array([0.0, 0.0, 0.0])
        back_right = np.array([1.0, 0.0, 0.0])

        centroid = _functions.numba_centroid_of_quadrilateral(
            front_left,
            front_right,
            back_left,
            back_right,
        )

        expected = np.array([0.5, 0.5, 0.0])

        npt.assert_allclose(centroid, expected, atol=1e-14)

    def test_translated_quadrilateral(self) -> None:
        """Centroid should translate consistently with the input points."""
        front_left = np.array([2.0, 3.0, 4.0])
        front_right = np.array([4.0, 3.0, 4.0])
        back_left = np.array([2.0, 1.0, 4.0])
        back_right = np.array([4.0, 1.0, 4.0])

        centroid = _functions.numba_centroid_of_quadrilateral(
            front_left,
            front_right,
            back_left,
            back_right,
        )

        expected = np.array([3.0, 2.0, 4.0])

        npt.assert_allclose(centroid, expected, atol=1e-14)

    def test_nonplanar_quadrilateral(self) -> None:
        """Centroid should equal the average of all vertex coordinates."""
        front_left = np.array([0.0, 0.0, 0.0])
        front_right = np.array([2.0, 0.0, 1.0])
        back_left = np.array([0.0, 2.0, 2.0])
        back_right = np.array([2.0, 2.0, 3.0])

        centroid = _functions.numba_centroid_of_quadrilateral(
            front_left,
            front_right,
            back_left,
            back_right,
        )

        expected = np.array([1.0, 1.0, 1.5])

        npt.assert_allclose(centroid, expected, atol=1e-14)

    def test_identical_points(self) -> None:
        """Centroid of identical points should equal that point."""
        point = np.array([1.5, -2.0, 3.25])

        centroid = _functions.numba_centroid_of_quadrilateral(
            point,
            point,
            point,
            point,
        )

        npt.assert_allclose(centroid, point, atol=1e-14)


class TestNumba1dExplicitCross(unittest.TestCase):
    """Tests for the numba_1d_explicit_cross function."""

    def test_known_cross_products(self) -> None:
        """Cross products should match analytically known results."""
        vectors1 = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )

        vectors2 = np.array(
            [
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

        cross_products = _functions.numba_1d_explicit_cross(
            vectors1,
            vectors2,
        )

        expected = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ]
        )

        npt.assert_allclose(cross_products, expected, atol=1e-14)

    def test_parallel_vectors(self) -> None:
        """Parallel vectors should produce zero cross products."""
        vectors1 = np.array(
            [
                [1.0, 2.0, 3.0],
                [2.0, 4.0, 6.0],
            ]
        )

        vectors2 = np.array(
            [
                [2.0, 4.0, 6.0],
                [1.0, 2.0, 3.0],
            ]
        )

        cross_products = _functions.numba_1d_explicit_cross(
            vectors1,
            vectors2,
        )

        expected = np.zeros((2, 3))

        npt.assert_allclose(cross_products, expected, atol=1e-14)

    def test_matches_numpy_cross(self) -> None:
        """Cross products should match numpy.cross."""
        vectors1 = np.array(
            [
                [1.0, 2.0, 3.0],
                [-1.0, 4.0, 0.5],
                [2.5, -3.0, 1.0],
            ]
        )

        vectors2 = np.array(
            [
                [4.0, 5.0, 6.0],
                [2.0, -1.0, 3.0],
                [0.0, 2.0, -4.0],
            ]
        )

        cross_products = _functions.numba_1d_explicit_cross(
            vectors1,
            vectors2,
        )

        expected = np.cross(vectors1, vectors2)

        npt.assert_allclose(cross_products, expected, atol=1e-14)

    def test_empty_input(self) -> None:
        """Empty vector stacks should return an empty array."""
        vectors1 = np.empty((0, 3))
        vectors2 = np.empty((0, 3))

        cross_products = _functions.numba_1d_explicit_cross(
            vectors1,
            vectors2,
        )

        expected = np.empty((0, 3))

        npt.assert_allclose(cross_products, expected, atol=1e-14)


class TestInterpBetweenPoints(unittest.TestCase):
    """Tests for the interp_between_points function."""

    def test_single_pair_linear_interpolation(self) -> None:
        """Interpolated points should lie linearly between a point pair."""
        start_points = np.array([[0.0, 0.0, 0.0]])
        end_points = np.array([[10.0, 0.0, 0.0]])

        norm_spacings = np.array([0.0, 0.5, 1.0])

        interpolated_points = _functions.interp_between_points(
            start_points,
            end_points,
            norm_spacings,
        )

        expected = np.array(
            [
                [
                    [0.0, 0.0, 0.0],
                    [5.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                ]
            ]
        )

        npt.assert_allclose(interpolated_points, expected, atol=1e-14)

    def test_multiple_point_pairs(self) -> None:
        """Interpolation should work independently for multiple point pairs."""
        start_points = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ]
        )

        end_points = np.array(
            [
                [2.0, 0.0, 0.0],
                [3.0, 3.0, 3.0],
            ]
        )

        norm_spacings = np.array([0.0, 0.5, 1.0])

        interpolated_points = _functions.interp_between_points(
            start_points,
            end_points,
            norm_spacings,
        )

        expected = np.array(
            [
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                ],
                [
                    [1.0, 1.0, 1.0],
                    [2.0, 2.0, 2.0],
                    [3.0, 3.0, 3.0],
                ],
            ]
        )

        npt.assert_allclose(interpolated_points, expected, atol=1e-14)

    def test_zero_spacing_returns_start_points(self) -> None:
        """Zero normalized spacing should return the start points."""
        start_points = np.array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
            ]
        )

        end_points = np.array(
            [
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
            ]
        )

        norm_spacings = np.array([0.0])

        interpolated_points = _functions.interp_between_points(
            start_points,
            end_points,
            norm_spacings,
        )

        expected = start_points[:, np.newaxis, :]

        npt.assert_allclose(interpolated_points, expected, atol=1e-14)

    def test_unit_spacing_returns_end_points(self) -> None:
        """Unit normalized spacing should return the end points."""
        start_points = np.array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
            ]
        )

        end_points = np.array(
            [
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
            ]
        )

        norm_spacings = np.array([1.0])

        interpolated_points = _functions.interp_between_points(
            start_points,
            end_points,
            norm_spacings,
        )

        expected = end_points[:, np.newaxis, :]

        npt.assert_allclose(interpolated_points, expected, atol=1e-14)


class TestFormatDuration(unittest.TestCase):
    """Tests for the format_duration function."""

    def test_sub_minute_duration_shows_seconds_only(self) -> None:
        """A duration under one minute should format as unpadded seconds alone."""
        self.assertEqual(_functions.format_duration(2.345), "2.35 s")

    def test_zero_duration(self) -> None:
        """A zero duration should format as zero seconds."""
        self.assertEqual(_functions.format_duration(0.0), "0.00 s")

    def test_sub_hour_duration_drops_hours_and_comma(self) -> None:
        """A duration under one hour should drop the hours and the comma."""
        self.assertEqual(_functions.format_duration(1234.5), "20 min and 34.5 s")

    def test_exactly_one_minute(self) -> None:
        """A duration of exactly one minute should roll over to the minutes form."""
        self.assertEqual(_functions.format_duration(60.0), "1 min and 0.00 s")

    def test_multi_hour_duration_shows_all_portions(self) -> None:
        """A duration over one hour should show hours, minutes, and seconds."""
        self.assertEqual(
            _functions.format_duration(45296.789), "12 hr, 34 min, and 56.8 s"
        )

    def test_exactly_one_hour(self) -> None:
        """A duration of exactly one hour should roll over to the hours form."""
        self.assertEqual(_functions.format_duration(3600.0), "1 hr, 0 min, and 0.00 s")

    def test_one_thousand_hours_shows_hours_only(self) -> None:
        """A duration of at least 1000 hours should format as hours alone."""
        self.assertEqual(_functions.format_duration(3.6e6), "1.00E+03 hr")

    def test_negative_sub_minute_duration(self) -> None:
        """A negative duration should carry its sign on the seconds form."""
        self.assertEqual(_functions.format_duration(-2.345), "-2.35 s")

    def test_negative_compound_duration_signs_largest_unit(self) -> None:
        """A negative duration should carry its sign on the largest unit."""
        self.assertEqual(
            _functions.format_duration(-45296.789), "-12 hr, 34 min, and 56.8 s"
        )

    def test_seconds_rounding_carries_into_minutes(self) -> None:
        """A seconds remainder that rounds to 60 should carry into the minutes."""
        self.assertEqual(_functions.format_duration(59.999), "1 min and 0.00 s")

    def test_seconds_rounding_carries_into_hours(self) -> None:
        """A carried seconds remainder should cascade into the hours."""
        self.assertEqual(
            _functions.format_duration(3599.999), "1 hr, 0 min, and 0.00 s"
        )

    def test_seconds_rounding_carries_into_hours_only_form(self) -> None:
        """A carried seconds remainder should cascade into the hours-only form."""
        self.assertEqual(_functions.format_duration(3.6e6 - 0.001), "1.00E+03 hr")

    def test_tiny_seconds_remainder_uses_scientific_notation(self) -> None:
        """A tiny seconds remainder should keep three significant figures."""
        self.assertEqual(_functions.format_duration(60.0000001), "1 min and 1.00E-07 s")

    def test_left_pad_gives_every_form_a_constant_width(self) -> None:
        """Left-padded results should share a constant width across all forms."""
        durations = [0.0, 2.345, 1234.5, 45296.789, 3.6e6, -45296.789]
        for duration in durations:
            unpadded = _functions.format_duration(duration)
            padded = _functions.format_duration(duration, left_pad=True)
            self.assertEqual(padded, unpadded.rjust(_functions._DURATION_PAD_WIDTH))
            self.assertEqual(len(padded), _functions._DURATION_PAD_WIDTH)


class TestSolveLoopThreadLimits(unittest.TestCase):
    """Tests for the solve_loop_thread_limits function."""

    @staticmethod
    def _blas_num_threads() -> dict[str, int]:
        """Returns each controllable BLAS library's current thread count, keyed by its
        file path."""
        return {
            library["filepath"]: library["num_threads"]
            for library in threadpoolctl.threadpool_info()
            if library["user_api"] == "blas"
        }

    def test_below_threshold_limits_blas_to_one_thread(self) -> None:
        """Inside the below-threshold context, every controllable BLAS library should be
        limited to 1 thread."""
        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            for library in threadpoolctl.threadpool_info():
                if library["user_api"] == "blas":
                    self.assertEqual(library["num_threads"], 1)

    def test_below_threshold_leaves_non_blas_pools_alone(self) -> None:
        """The below-threshold context should not limit non BLAS thread pools, such as
        the one behind Numba's threading layer."""
        outside_num_threads = {
            library["filepath"]: library["num_threads"]
            for library in threadpoolctl.threadpool_info()
            if library["user_api"] != "blas"
        }

        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            inside_num_threads = {
                library["filepath"]: library["num_threads"]
                for library in threadpoolctl.threadpool_info()
                if library["user_api"] != "blas"
            }

        self.assertEqual(inside_num_threads, outside_num_threads)

    def test_below_threshold_restores_blas_after_exiting(self) -> None:
        """The below-threshold context should restore each BLAS library's original
        thread count when it exits."""
        outside_num_threads = self._blas_num_threads()

        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            pass

        self.assertEqual(self._blas_num_threads(), outside_num_threads)

    def test_at_threshold_leaves_blas_at_full_width(self) -> None:
        """A Panel count at the threshold should leave every BLAS library's thread count
        untouched."""
        outside_num_threads = self._blas_num_threads()

        with _functions.solve_loop_thread_limits(_functions._SOLVE_THREAD_THRESHOLD):
            inside_num_threads = self._blas_num_threads()

        self.assertEqual(inside_num_threads, outside_num_threads)

    def test_above_threshold_leaves_blas_at_full_width(self) -> None:
        """A Panel count above the threshold should leave every BLAS library's thread
        count untouched."""
        outside_num_threads = self._blas_num_threads()

        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD + 1
        ):
            inside_num_threads = self._blas_num_threads()

        self.assertEqual(inside_num_threads, outside_num_threads)

    def test_concurrent_solve_loop_in_another_thread_raises(self) -> None:
        """A second solve loop entered while one is already running should raise, rather
        than corrupt the process-wide BLAS thread count that both would save and
        restore."""
        first_entered = threading.Event()
        second_may_exit = threading.Event()
        second_error = []

        def run_second() -> None:
            first_entered.wait(5)
            try:
                with _functions.solve_loop_thread_limits(
                    _functions._SOLVE_THREAD_THRESHOLD - 1
                ):
                    pass
            except RuntimeError as error:
                second_error.append(error)
            finally:
                second_may_exit.set()

        second_thread = threading.Thread(target=run_second, name="SecondSolveLoop")
        second_thread.start()

        # Hold the first solve loop open until the second one has tried to start its
        # own, so the two genuinely overlap.
        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            first_entered.set()
            second_may_exit.wait(5)

        second_thread.join(5)

        self.assertEqual(len(second_error), 1)
        self.assertIn("already in progress", str(second_error[0]))
        self.assertIn("separate processes", str(second_error[0]))

    def test_concurrent_solve_loop_leaves_blas_uncorrupted(self) -> None:
        """After a rejected concurrent solve loop, the first solve loop should still
        restore each BLAS library's original thread count."""
        outside_num_threads = self._blas_num_threads()

        first_entered = threading.Event()
        second_may_exit = threading.Event()

        def run_second() -> None:
            first_entered.wait(5)
            try:
                with _functions.solve_loop_thread_limits(
                    _functions._SOLVE_THREAD_THRESHOLD - 1
                ):
                    pass
            except RuntimeError:
                pass
            finally:
                second_may_exit.set()

        second_thread = threading.Thread(target=run_second, name="SecondSolveLoop")
        second_thread.start()

        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            first_entered.set()
            second_may_exit.wait(5)

        second_thread.join(5)

        self.assertEqual(self._blas_num_threads(), outside_num_threads)

    def test_sequential_solve_loops_in_different_threads_are_allowed(self) -> None:
        """The guard should reject only overlapping solve loops, so a solve loop that
        starts after another has finished should run even from a different thread."""
        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            pass

        errors = []

        def run_after() -> None:
            try:
                with _functions.solve_loop_thread_limits(
                    _functions._SOLVE_THREAD_THRESHOLD - 1
                ):
                    pass
            except RuntimeError as error:
                errors.append(error)

        later_thread = threading.Thread(target=run_after, name="LaterSolveLoop")
        later_thread.start()
        later_thread.join(5)

        self.assertEqual(errors, [])

    def test_solve_loop_releases_the_guard_when_its_body_raises(self) -> None:
        """A solver run that raises should still release the guard, so the next run is
        not rejected forever."""
        with self.assertRaises(ZeroDivisionError):
            with _functions.solve_loop_thread_limits(
                _functions._SOLVE_THREAD_THRESHOLD - 1
            ):
                raise ZeroDivisionError

        # The guard is free again, so this must not raise.
        with _functions.solve_loop_thread_limits(
            _functions._SOLVE_THREAD_THRESHOLD - 1
        ):
            pass
