"""This module contains a class to test LineVortices."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware._vortices import _line_vortex
from tests.unit.fixtures import line_vortex_fixtures


class TestLineVortex(unittest.TestCase):
    """This is a class with functions to test LineVortices."""

    def setUp(self):
        """Set up test fixtures for LineVortex tests."""
        # Create fixtures using line_vortex_fixtures module.
        self.basic_line_vortex = line_vortex_fixtures.make_basic_line_vortex_fixture()
        self.diagonal_line_vortex = (
            line_vortex_fixtures.make_diagonal_line_vortex_fixture()
        )
        self.y_aligned_line_vortex = (
            line_vortex_fixtures.make_y_aligned_line_vortex_fixture()
        )
        self.z_aligned_line_vortex = (
            line_vortex_fixtures.make_z_aligned_line_vortex_fixture()
        )
        self.zero_strength_line_vortex = (
            line_vortex_fixtures.make_zero_strength_line_vortex_fixture()
        )
        self.negative_strength_line_vortex = (
            line_vortex_fixtures.make_negative_strength_line_vortex_fixture()
        )
        self.large_strength_line_vortex = (
            line_vortex_fixtures.make_large_strength_line_vortex_fixture()
        )
        self.offset_line_vortex = line_vortex_fixtures.make_offset_line_vortex_fixture()
        self.small_line_vortex = line_vortex_fixtures.make_small_line_vortex_fixture()
        self.long_line_vortex = line_vortex_fixtures.make_long_line_vortex_fixture()

    def test_initialization_valid_parameters(self):
        """Test LineVortex initialization with valid parameters."""
        # Test that basic LineVortex initializes correctly.
        self.assertIsInstance(self.basic_line_vortex, _line_vortex.LineVortex)
        npt.assert_array_equal(
            self.basic_line_vortex.Slvp_GP1_CgP1,
            np.array([0.0, 0.0, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.basic_line_vortex.Elvp_GP1_CgP1,
            np.array([1.0, 0.0, 0.0], dtype=float),
        )
        self.assertEqual(self.basic_line_vortex.strength, 1.0)

    def test_initialization_with_different_strengths(self):
        """Test LineVortex initialization with various strength values."""
        # Test with zero strength.
        self.assertEqual(self.zero_strength_line_vortex.strength, 0.0)

        # Test with negative strength.
        self.assertEqual(self.negative_strength_line_vortex.strength, -2.5)

        # Test with large strength.
        self.assertEqual(self.large_strength_line_vortex.strength, 1e6)

    def test_initialization_with_different_orientations(self):
        """Test LineVortex initialization with various orientations."""
        # Test y aligned.
        npt.assert_array_equal(
            self.y_aligned_line_vortex.Slvp_GP1_CgP1,
            np.array([0.0, -0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.y_aligned_line_vortex.Elvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.0], dtype=float),
        )

        # Test z aligned.
        npt.assert_array_equal(
            self.z_aligned_line_vortex.Slvp_GP1_CgP1,
            np.array([0.0, 0.0, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.z_aligned_line_vortex.Elvp_GP1_CgP1,
            np.array([0.0, 0.0, 2.0], dtype=float),
        )

        # Test diagonal.
        npt.assert_array_equal(
            self.diagonal_line_vortex.Slvp_GP1_CgP1,
            np.array([0.0, 0.0, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.diagonal_line_vortex.Elvp_GP1_CgP1,
            np.array([1.0, 1.0, 1.0], dtype=float),
        )

    def test_vector_property_basic(self):
        """Test that vector_GP1 property returns correct value for basic vortex."""
        # For basic LineVortex, vector should be [1.0, 0.0, 0.0].
        expected_vector = np.array([1.0, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.basic_line_vortex.vector_GP1,
            expected_vector,
            decimal=10,
        )

    def test_vector_property_diagonal(self):
        """Test that vector_GP1 property returns correct value for diagonal vortex."""
        # For diagonal LineVortex, vector should be [1.0, 1.0, 1.0].
        expected_vector = np.array([1.0, 1.0, 1.0], dtype=float)
        npt.assert_array_almost_equal(
            self.diagonal_line_vortex.vector_GP1,
            expected_vector,
            decimal=10,
        )

    def test_vector_property_y_aligned(self):
        """Test that vector_GP1 property returns correct value for y aligned vortex."""
        # For y aligned LineVortex, vector should be [0.0, 1.0, 0.0].
        expected_vector = np.array([0.0, 1.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.y_aligned_line_vortex.vector_GP1,
            expected_vector,
            decimal=10,
        )

    def test_vector_property_z_aligned(self):
        """Test that vector_GP1 property returns correct value for z aligned vortex."""
        # For z aligned LineVortex, vector should be [0.0, 0.0, 2.0].
        expected_vector = np.array([0.0, 0.0, 2.0], dtype=float)
        npt.assert_array_almost_equal(
            self.z_aligned_line_vortex.vector_GP1,
            expected_vector,
            decimal=10,
        )

    def test_center_point_property_basic(self):
        """Test that Clvp_GP1_CgP1 property returns correct value for basic vortex."""
        # For basic LineVortex from [0, 0, 0] to [1, 0, 0], center should be [0.5, 0,
        # 0].
        expected_center = np.array([0.5, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.basic_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_center_point_property_diagonal(self):
        """Test that Clvp_GP1_CgP1 property returns correct value for diagonal
        vortex."""
        # For diagonal LineVortex from [0, 0, 0] to [1, 1, 1], center should be
        # [0.5, 0.5, 0.5].
        expected_center = np.array([0.5, 0.5, 0.5], dtype=float)
        npt.assert_array_almost_equal(
            self.diagonal_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_center_point_property_y_aligned(self):
        """Test that Clvp_GP1_CgP1 property returns correct value for y aligned
        vortex."""
        # For y aligned LineVortex from [0, -0.5, 0] to [0, 0.5, 0], center should be
        # [0, 0, 0].
        expected_center = np.array([0.0, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.y_aligned_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_center_point_property_z_aligned(self):
        """Test that Clvp_GP1_CgP1 property returns correct value for z aligned
        vortex."""
        # For z aligned LineVortex from [0, 0, 0] to [0, 0, 2], center should be
        # [0, 0, 1].
        expected_center = np.array([0.0, 0.0, 1.0], dtype=float)
        npt.assert_array_almost_equal(
            self.z_aligned_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_center_point_property_offset(self):
        """Test that Clvp_GP1_CgP1 property returns correct value for offset vortex."""
        # For offset LineVortex, center should be [5.5, 3.0, 2.0].
        expected_center = np.array([5.5, 3.0, 2.0], dtype=float)
        npt.assert_array_almost_equal(
            self.offset_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_small_line_vortex(self):
        """Test LineVortex with very small dimensions."""
        # Should handle very small LineVortices without issues.
        expected_vector = np.array([0.001, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.small_line_vortex.vector_GP1,
            expected_vector,
            decimal=10,
        )

        expected_center = np.array([0.0005, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.small_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_long_line_vortex(self):
        """Test LineVortex with very large dimensions."""
        # Should handle very long LineVortices without issues.
        expected_vector = np.array([1000.0, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.long_line_vortex.vector_GP1,
            expected_vector,
            decimal=10,
        )

        expected_center = np.array([500.0, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.long_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

    def test_offset_vortex(self):
        """Test LineVortex positioned far from origin."""
        # Verify offset was applied correctly.
        expected_offset = np.array([5.0, 3.0, 2.0], dtype=float)
        npt.assert_array_almost_equal(
            self.offset_line_vortex.Slvp_GP1_CgP1,
            np.array([0.0, 0.0, 0.0], dtype=float) + expected_offset,
            decimal=10,
        )
        npt.assert_array_almost_equal(
            self.offset_line_vortex.Elvp_GP1_CgP1,
            np.array([1.0, 0.0, 0.0], dtype=float) + expected_offset,
            decimal=10,
        )

    def test_strength_can_be_modified(self):
        """Test that strength can be modified after initialization."""
        # Create a fresh fixture.
        line_vortex = line_vortex_fixtures.make_basic_line_vortex_fixture()

        # Modify strength.
        line_vortex.strength = 5.0
        self.assertEqual(line_vortex.strength, 5.0)

        # Modify to negative.
        line_vortex.strength = -3.0
        self.assertEqual(line_vortex.strength, -3.0)

        # Modify to zero.
        line_vortex.strength = 0.0
        self.assertEqual(line_vortex.strength, 0.0)

    def test_vector_length(self):
        """Test that vector length matches distance between endpoints."""
        # For basic LineVortex, length should be 1.0.
        vector_length = np.linalg.norm(self.basic_line_vortex.vector_GP1)
        npt.assert_almost_equal(vector_length, 1.0, decimal=10)

        # For diagonal LineVortex, length should be sqrt(3).
        vector_length = np.linalg.norm(self.diagonal_line_vortex.vector_GP1)
        npt.assert_almost_equal(vector_length, np.sqrt(3), decimal=10)

        # For z aligned LineVortex, length should be 2.0.
        vector_length = np.linalg.norm(self.z_aligned_line_vortex.vector_GP1)
        npt.assert_almost_equal(vector_length, 2.0, decimal=10)

    def test_center_is_midpoint_of_endpoints(self):
        """Test that center point is exactly the midpoint of start and end points."""
        # For all fixtures, center should be (start + end) / 2.
        expected_center = (
            self.basic_line_vortex.Slvp_GP1_CgP1 + self.basic_line_vortex.Elvp_GP1_CgP1
        ) / 2.0
        npt.assert_array_almost_equal(
            self.basic_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

        expected_center = (
            self.diagonal_line_vortex.Slvp_GP1_CgP1
            + self.diagonal_line_vortex.Elvp_GP1_CgP1
        ) / 2.0
        npt.assert_array_almost_equal(
            self.diagonal_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )

        expected_center = (
            self.offset_line_vortex.Slvp_GP1_CgP1
            + self.offset_line_vortex.Elvp_GP1_CgP1
        ) / 2.0
        npt.assert_array_almost_equal(
            self.offset_line_vortex.Clvp_GP1_CgP1,
            expected_center,
            decimal=10,
        )


if __name__ == "__main__":
    unittest.main()
