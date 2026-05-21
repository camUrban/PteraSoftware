"""This module contains classes to test shared utility functions."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _functions


class TestCosspace(unittest.TestCase):
    """Tests for the cosspace function."""

    def test_includes_endpoints(self):
        """Cosspace should include minimum and maximum by default."""
        points = _functions.cosspace(0.0, 10.0, n_points=5)

        self.assertAlmostEqual(points[0], 0.0)
        self.assertAlmostEqual(points[-1], 10.0)

    def test_excludes_endpoint_when_requested(self):
        """Cosspace should exclude maximum when endpoint=False."""
        points = _functions.cosspace(
            0.0,
            10.0,
            n_points=5,
            endpoint=False,
        )

        self.assertAlmostEqual(points[0], 0.0)
        self.assertLess(points[-1], 10.0)

    def test_returns_correct_number_of_points(self):
        """Cosspace should return the requested number of points."""
        points = _functions.cosspace(0.0, 1.0, n_points=17)

        self.assertEqual(len(points), 17)

    def test_spacing_is_symmetric(self):
        """Cosspace output should be symmetric about the midpoint."""
        points = _functions.cosspace(-1.0, 1.0, n_points=9)

        npt.assert_allclose(points, -points[::-1], atol=1e-14)

    def test_single_point(self):
        """Single-point cosspace should return the minimum value."""
        points = _functions.cosspace(2.0, 8.0, n_points=1)

        npt.assert_allclose(points, np.array([2.0]), atol=1e-14)


class TestNumbaCentroidOfQuadrilateral(unittest.TestCase):
    """Tests for the numba_centroid_of_quadrilateral function."""

    def test_unit_square_in_xy_plane(self):
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

    def test_translated_quadrilateral(self):
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

    def test_nonplanar_quadrilateral(self):
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

    def test_identical_points(self):
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

    def test_known_cross_products(self):
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

    def test_parallel_vectors(self):
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

    def test_matches_numpy_cross(self):
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

    def test_empty_input(self):
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

    def test_single_pair_linear_interpolation(self):
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

    def test_multiple_point_pairs(self):
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

    def test_zero_spacing_returns_start_points(self):
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

    def test_unit_spacing_returns_end_points(self):
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


if __name__ == "__main__":
    unittest.main()
