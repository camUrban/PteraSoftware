"""This module contains classes to test RingVortices."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _vortices

# noinspection PyProtectedMember
from pterasoftware._vortices import _line_vortex
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
        self.assertIsInstance(self.basic_ring_vortex, _vortices.ring_vortex.RingVortex)
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
        self.assertIsInstance(self.basic_ring_vortex.front_leg, _line_vortex.LineVortex)
        self.assertIsInstance(self.basic_ring_vortex.left_leg, _line_vortex.LineVortex)
        self.assertIsInstance(self.basic_ring_vortex.back_leg, _line_vortex.LineVortex)
        self.assertIsInstance(self.basic_ring_vortex.right_leg, _line_vortex.LineVortex)

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

    def test_strength_propagation_to_legs(self):
        """Test that changing RingVortex strength propagates to all legs."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access all legs to create them.
        _ = ring_vortex.front_leg
        _ = ring_vortex.left_leg
        _ = ring_vortex.back_leg
        _ = ring_vortex.right_leg

        # Change strength.
        new_strength = 5.0
        ring_vortex.strength = new_strength

        # Verify all legs have the new strength.
        self.assertEqual(ring_vortex.front_leg.strength, new_strength)
        self.assertEqual(ring_vortex.left_leg.strength, new_strength)
        self.assertEqual(ring_vortex.back_leg.strength, new_strength)
        self.assertEqual(ring_vortex.right_leg.strength, new_strength)

    def test_strength_propagation_to_uncreated_legs(self):
        """Test that changing strength before accessing legs still works."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Change strength before accessing legs.
        new_strength = 7.5
        ring_vortex.strength = new_strength

        # Access legs now and verify they have the new strength.
        self.assertEqual(ring_vortex.front_leg.strength, new_strength)
        self.assertEqual(ring_vortex.left_leg.strength, new_strength)
        self.assertEqual(ring_vortex.back_leg.strength, new_strength)
        self.assertEqual(ring_vortex.right_leg.strength, new_strength)

    def test_strength_zero_to_nonzero(self):
        """Test changing strength from zero to nonzero."""
        # Create a zero strength fixture.
        ring_vortex = ring_vortex_fixtures.make_zero_strength_ring_vortex_fixture()

        # Access legs.
        _ = ring_vortex.front_leg

        # Change to nonzero.
        ring_vortex.strength = 10.0

        # Verify.
        self.assertEqual(ring_vortex.strength, 10.0)
        self.assertEqual(ring_vortex.front_leg.strength, 10.0)


class TestRingVortexImmutability(unittest.TestCase):
    """Tests for RingVortex attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.basic_ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

    def test_immutable_Frrvp_GP1_CgP1_property(self):
        """Test that Frrvp_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.Frrvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Frrvp_GP1_CgP1_array_read_only(self):
        """Test that Frrvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_ring_vortex.Frrvp_GP1_CgP1[0] = 999.0

    def test_immutable_Flrvp_GP1_CgP1_property(self):
        """Test that Flrvp_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.Flrvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Flrvp_GP1_CgP1_array_read_only(self):
        """Test that Flrvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_ring_vortex.Flrvp_GP1_CgP1[0] = 999.0

    def test_immutable_Blrvp_GP1_CgP1_property(self):
        """Test that Blrvp_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.Blrvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Blrvp_GP1_CgP1_array_read_only(self):
        """Test that Blrvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_ring_vortex.Blrvp_GP1_CgP1[0] = 999.0

    def test_immutable_Brrvp_GP1_CgP1_property(self):
        """Test that Brrvp_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.Brrvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Brrvp_GP1_CgP1_array_read_only(self):
        """Test that Brrvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_ring_vortex.Brrvp_GP1_CgP1[0] = 999.0

    def test_immutable_Crvp_GP1_CgP1_property(self):
        """Test that Crvp_GP1_CgP1 property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_ring_vortex.Crvp_GP1_CgP1
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.Crvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Crvp_GP1_CgP1_array_read_only(self):
        """Test that Crvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_ring_vortex.Crvp_GP1_CgP1[0] = 999.0

    def test_immutable_front_leg_property(self):
        """Test that front_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_ring_vortex.front_leg
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.front_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_immutable_left_leg_property(self):
        """Test that left_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_ring_vortex.left_leg
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.left_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_immutable_back_leg_property(self):
        """Test that back_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_ring_vortex.back_leg
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.back_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_immutable_right_leg_property(self):
        """Test that right_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_ring_vortex.right_leg
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.right_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_immutable_area_property(self):
        """Test that area property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_ring_vortex.area
        with self.assertRaises(AttributeError):
            self.basic_ring_vortex.area = 5.0

    def test_mutable_strength_property(self):
        """Test that strength property is mutable."""
        # Strength should be modifiable.
        self.basic_ring_vortex.strength = 5.0
        self.assertEqual(self.basic_ring_vortex.strength, 5.0)

        self.basic_ring_vortex.strength = -3.0
        self.assertEqual(self.basic_ring_vortex.strength, -3.0)

    def test_mutable_age_attribute(self):
        """Test that age attribute is mutable."""
        # Age should be modifiable.
        self.basic_ring_vortex.age = 5.0
        self.assertEqual(self.basic_ring_vortex.age, 5.0)

        self.basic_ring_vortex.age = 10.0
        self.assertEqual(self.basic_ring_vortex.age, 10.0)


class TestRingVortexLazyCaching(unittest.TestCase):
    """Tests for RingVortex lazy caching behavior."""

    def test_centroid_lazy_evaluation(self):
        """Test that centroid is lazily evaluated."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Verify that _Crvp_GP1_CgP1 is None before access.
        self.assertIsNone(ring_vortex._Crvp_GP1_CgP1)

        # Access the property.
        _ = ring_vortex.Crvp_GP1_CgP1

        # Verify that the cache is now populated.
        self.assertIsNotNone(ring_vortex._Crvp_GP1_CgP1)

    def test_centroid_caching(self):
        """Test that centroid returns the same cached object on repeated access."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access centroid twice.
        centroid_first = ring_vortex.Crvp_GP1_CgP1
        centroid_second = ring_vortex.Crvp_GP1_CgP1

        # Verify they are the same object (not just equal).
        self.assertIs(centroid_first, centroid_second)

    def test_front_leg_lazy_evaluation(self):
        """Test that front_leg is lazily evaluated."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Verify that _front_leg is None before access.
        self.assertIsNone(ring_vortex._front_leg)

        # Access the property.
        _ = ring_vortex.front_leg

        # Verify that the cache is now populated.
        self.assertIsNotNone(ring_vortex._front_leg)

    def test_front_leg_caching(self):
        """Test that front_leg returns the same cached object on repeated access."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access front_leg twice.
        front_leg_first = ring_vortex.front_leg
        front_leg_second = ring_vortex.front_leg

        # Verify they are the same object.
        self.assertIs(front_leg_first, front_leg_second)

    def test_left_leg_lazy_evaluation(self):
        """Test that left_leg is lazily evaluated."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Verify that _left_leg is None before access.
        self.assertIsNone(ring_vortex._left_leg)

        # Access the property.
        _ = ring_vortex.left_leg

        # Verify that the cache is now populated.
        self.assertIsNotNone(ring_vortex._left_leg)

    def test_left_leg_caching(self):
        """Test that left_leg returns the same cached object on repeated access."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access left_leg twice.
        left_leg_first = ring_vortex.left_leg
        left_leg_second = ring_vortex.left_leg

        # Verify they are the same object.
        self.assertIs(left_leg_first, left_leg_second)

    def test_back_leg_lazy_evaluation(self):
        """Test that back_leg is lazily evaluated."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Verify that _back_leg is None before access.
        self.assertIsNone(ring_vortex._back_leg)

        # Access the property.
        _ = ring_vortex.back_leg

        # Verify that the cache is now populated.
        self.assertIsNotNone(ring_vortex._back_leg)

    def test_back_leg_caching(self):
        """Test that back_leg returns the same cached object on repeated access."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access back_leg twice.
        back_leg_first = ring_vortex.back_leg
        back_leg_second = ring_vortex.back_leg

        # Verify they are the same object.
        self.assertIs(back_leg_first, back_leg_second)

    def test_right_leg_lazy_evaluation(self):
        """Test that right_leg is lazily evaluated."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Verify that _right_leg is None before access.
        self.assertIsNone(ring_vortex._right_leg)

        # Access the property.
        _ = ring_vortex.right_leg

        # Verify that the cache is now populated.
        self.assertIsNotNone(ring_vortex._right_leg)

    def test_right_leg_caching(self):
        """Test that right_leg returns the same cached object on repeated access."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access right_leg twice.
        right_leg_first = ring_vortex.right_leg
        right_leg_second = ring_vortex.right_leg

        # Verify they are the same object.
        self.assertIs(right_leg_first, right_leg_second)

    def test_area_lazy_evaluation(self):
        """Test that area is lazily evaluated."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Verify that _area is None before access.
        self.assertIsNone(ring_vortex._area)

        # Access the property.
        _ = ring_vortex.area

        # Verify that the cache is now populated.
        self.assertIsNotNone(ring_vortex._area)

    def test_area_caching(self):
        """Test that area returns the same cached value on repeated access."""
        # Create a fresh fixture.
        ring_vortex = ring_vortex_fixtures.make_basic_ring_vortex_fixture()

        # Access area twice.
        area_first = ring_vortex.area
        area_second = ring_vortex.area

        # Verify they are equal (floats are immutable so no identity check).
        self.assertEqual(area_first, area_second)


if __name__ == "__main__":
    unittest.main()
