"""This module contains a class to test RingVortices."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _aerodynamics
from tests.unit.fixtures import ring_vortex_fixtures


class TestRingVortex(unittest.TestCase):
    """This is a class with functions to test RingVortices."""

    def setUp(self):
        """Set up test fixtures for RingVortex tests."""
        # Create fixtures using ring_vortex_fixtures module.
        self.basic_ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()
        self.unit_square_ring_vortex = (
            ring_vortex_fixtures.make_unit_square_ring_vortex_fixture()
        )
        self.rectangular_ring_vortex = (
            ring_vortex_fixtures.make_rectangular_ring_vortex_fixture()
        )
        self.tilted_ring_vortex = ring_vortex_fixtures.make_tilted_ring_vortex_fixture()
        self.zero_strength_ring_vortex = (
            ring_vortex_fixtures.make_zero_strength_ring_vortex_fixture()
        )
        self.negative_strength_ring_vortex = (
            ring_vortex_fixtures.make_negative_strength_ring_vortex_fixture()
        )
        self.offset_ring_vortex = ring_vortex_fixtures.make_offset_ring_vortex_fixture()
        self.small_ring_vortex = ring_vortex_fixtures.make_small_ring_vortex_fixture()
        self.large_strength_ring_vortex = (
            ring_vortex_fixtures.make_large_strength_ring_vortex_fixture()
        )
        self.large_coordinate_ring_vortex = (
            ring_vortex_fixtures.make_large_coordinate_ring_vortex_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test RingVortex initialization with valid parameters."""
        # Test that basic RingVortex initializes correctly.
        self.assertIsInstance(self.basic_ring_vortex, _aerodynamics.RingVortex)
        npt.assert_array_equal(
            self.basic_ring_vortex.Frrvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.Flrvp_GP1_CgP1,
            np.array([0.0, -0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.Blrvp_GP1_CgP1,
            np.array([1.0, -0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.Brrvp_GP1_CgP1,
            np.array([1.0, 0.5, 0.0], dtype=float),
        )
        self.assertEqual(self.basic_ring_vortex.strength, 1.0)
        self.assertEqual(self.basic_ring_vortex.age, 0)

    def test_initialization_with_different_strengths(self):
        """Test RingVortex initialization with various strength values."""
        # Test with zero strength.
        self.assertEqual(self.zero_strength_ring_vortex.strength, 0.0)

        # Test with negative strength.
        self.assertEqual(self.negative_strength_ring_vortex.strength, -1.0)

        # Test with positive strength.
        self.assertEqual(self.rectangular_ring_vortex.strength, 2.5)

    def test_initialization_with_different_geometries(self):
        """Test RingVortex initialization with various geometries."""
        # Test unit square geometry.
        npt.assert_array_equal(
            self.unit_square_ring_vortex.Frrvp_GP1_CgP1,
            np.array([0.5, 0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.unit_square_ring_vortex.Flrvp_GP1_CgP1,
            np.array([0.5, -0.5, 0.0], dtype=float),
        )

        # Test rectangular geometry.
        npt.assert_array_equal(
            self.rectangular_ring_vortex.Brrvp_GP1_CgP1,
            np.array([2.0, 0.5, 0.0], dtype=float),
        )

        # Test tilted geometry.
        npt.assert_array_equal(
            self.tilted_ring_vortex.Frrvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.2], dtype=float),
        )

    def test_line_vortex_legs_created(self):
        """Test that RingVortex creates four LineVortex legs correctly."""
        # Check that all four legs are created.
        self.assertIsNotNone(self.basic_ring_vortex.front_leg)
        self.assertIsNotNone(self.basic_ring_vortex.left_leg)
        self.assertIsNotNone(self.basic_ring_vortex.back_leg)
        self.assertIsNotNone(self.basic_ring_vortex.right_leg)

        # Check that leg strengths match the RingVortex strength.
        self.assertEqual(
            self.basic_ring_vortex.front_leg.strength, self.basic_ring_vortex.strength
        )
        self.assertEqual(
            self.basic_ring_vortex.left_leg.strength, self.basic_ring_vortex.strength
        )
        self.assertEqual(
            self.basic_ring_vortex.back_leg.strength, self.basic_ring_vortex.strength
        )
        self.assertEqual(
            self.basic_ring_vortex.right_leg.strength, self.basic_ring_vortex.strength
        )

    def test_line_vortex_legs_positions(self):
        """Test that LineVortex leg positions match RingVortex corners."""
        # Test front leg (from front right to front left).
        npt.assert_array_equal(
            self.basic_ring_vortex.front_leg.Slvp_GP1_CgP1,
            self.basic_ring_vortex.Frrvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.front_leg.Elvp_GP1_CgP1,
            self.basic_ring_vortex.Flrvp_GP1_CgP1,
        )

        # Test left leg (from front left to back left).
        npt.assert_array_equal(
            self.basic_ring_vortex.left_leg.Slvp_GP1_CgP1,
            self.basic_ring_vortex.Flrvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.left_leg.Elvp_GP1_CgP1,
            self.basic_ring_vortex.Blrvp_GP1_CgP1,
        )

        # Test back leg (from back left to back right).
        npt.assert_array_equal(
            self.basic_ring_vortex.back_leg.Slvp_GP1_CgP1,
            self.basic_ring_vortex.Blrvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.back_leg.Elvp_GP1_CgP1,
            self.basic_ring_vortex.Brrvp_GP1_CgP1,
        )

        # Test right leg (from back right to front right).
        npt.assert_array_equal(
            self.basic_ring_vortex.right_leg.Slvp_GP1_CgP1,
            self.basic_ring_vortex.Brrvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.right_leg.Elvp_GP1_CgP1,
            self.basic_ring_vortex.Frrvp_GP1_CgP1,
        )

    def test_centroid_calculation(self):
        """Test that RingVortex centroid is calculated correctly."""
        # For a unit square centered at origin, centroid should be at origin.
        expectedCentroid_GP1_CgP1 = np.array([0.0, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.unit_square_ring_vortex.Crvp_GP1_CgP1,
            expectedCentroid_GP1_CgP1,
            decimal=10,
        )

        # For the basic RingVortex fixture, the centroid should be at [0.5, 0.0,
        # 0.0] (in the first Airplane's geometry axes, relative to the first
        # Airplane's CG).
        expectedCentroid_GP1_CgP1 = np.array([0.5, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.basic_ring_vortex.Crvp_GP1_CgP1, expectedCentroid_GP1_CgP1, decimal=10
        )

    def test_initial_age_is_zero(self):
        """Test that RingVortex age is initialized to zero."""
        self.assertEqual(self.basic_ring_vortex.age, 0)
        self.assertEqual(self.unit_square_ring_vortex.age, 0)
        self.assertEqual(self.tilted_ring_vortex.age, 0)

    def test_very_small_vortex(self):
        """Test RingVortex with very small dimensions."""
        # Should handle very small RingVortices without issues.
        self.assertIsNotNone(self.small_ring_vortex.Crvp_GP1_CgP1)
        self.assertIsNotNone(self.small_ring_vortex.front_leg)

    def test_offset_vortex(self):
        """Test RingVortex positioned far from origin."""
        # Verify offset was applied correctly.
        expected_offset = np.array([5.0, 3.0, 2.0], dtype=float)
        npt.assert_array_almost_equal(
            self.offset_ring_vortex.Frrvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.0], dtype=float) + expected_offset,
            decimal=10,
        )

    def test_area_property_unit_square(self):
        """Test that area property returns correct value for unit square."""
        # For unit square RingVortex, area should be 1.0 square meter.
        expected_area = 1.0
        npt.assert_almost_equal(
            self.unit_square_ring_vortex.area, expected_area, decimal=10
        )

    def test_area_property_rectangular(self):
        """Test that area property returns correct value for rectangular vortex."""
        # For rectangular RingVortex (2x1), area should be 2.0 square meters.
        expected_area = 2.0
        npt.assert_almost_equal(
            self.rectangular_ring_vortex.area, expected_area, decimal=10
        )

    def test_area_property_basic(self):
        """Test that area property returns correct value for basic vortex."""
        # For basic RingVortex (1x1), area should be 1.0 square meter.
        expected_area = 1.0
        npt.assert_almost_equal(self.basic_ring_vortex.area, expected_area, decimal=10)

    def test_area_property_small(self):
        """Test that area property returns correct value for small vortex."""
        # For small RingVortex (0.01x0.01), area should be 0.0001 square meters.
        expected_area = 0.0001
        npt.assert_almost_equal(self.small_ring_vortex.area, expected_area, decimal=10)

    def test_area_property_positive(self):
        """Test that area property is always positive."""
        # All areas should be positive.
        self.assertGreater(self.basic_ring_vortex.area, 0)
        self.assertGreater(self.unit_square_ring_vortex.area, 0)
        self.assertGreater(self.rectangular_ring_vortex.area, 0)
        self.assertGreater(self.tilted_ring_vortex.area, 0)
        self.assertGreater(self.small_ring_vortex.area, 0)
        self.assertGreater(self.offset_ring_vortex.area, 0)

    def test_line_vortex_legs_are_line_vortex_type(self):
        """Test that LineVortex legs are of the correct type."""
        # Verify all legs are LineVortex instances.
        self.assertIsInstance(
            self.basic_ring_vortex.front_leg, _aerodynamics._LineVortex
        )
        self.assertIsInstance(
            self.basic_ring_vortex.left_leg, _aerodynamics._LineVortex
        )
        self.assertIsInstance(
            self.basic_ring_vortex.back_leg, _aerodynamics._LineVortex
        )
        self.assertIsInstance(
            self.basic_ring_vortex.right_leg, _aerodynamics._LineVortex
        )

    def test_large_strength_vortex(self):
        """Test RingVortex with very large strength."""
        # Verify the large strength RingVortex initializes correctly.
        self.assertEqual(self.large_strength_ring_vortex.strength, 1e6)
        self.assertEqual(self.large_strength_ring_vortex.front_leg.strength, 1e6)
        self.assertEqual(self.large_strength_ring_vortex.left_leg.strength, 1e6)
        self.assertEqual(self.large_strength_ring_vortex.back_leg.strength, 1e6)
        self.assertEqual(self.large_strength_ring_vortex.right_leg.strength, 1e6)

    def test_large_coordinate_vortex(self):
        """Test RingVortex with very large coordinates."""
        # Verify the large coordinate RingVortex initializes correctly.
        expected_offset = np.array([1e6, 1e6, 1e6], dtype=float)
        npt.assert_array_almost_equal(
            self.large_coordinate_ring_vortex.Frrvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.0], dtype=float) + expected_offset,
            decimal=4,
        )

        # Verify centroid is calculated correctly for large coordinates.
        expected_centroid = np.array([0.5, 0.0, 0.0], dtype=float) + expected_offset
        npt.assert_array_almost_equal(
            self.large_coordinate_ring_vortex.Crvp_GP1_CgP1,
            expected_centroid,
            decimal=4,
        )


if __name__ == "__main__":
    unittest.main()
