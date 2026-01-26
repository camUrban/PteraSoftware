"""This module contains a class to test HorseshoeVortices."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _aerodynamics
from tests.unit.fixtures import horseshoe_vortex_fixtures


class TestHorseshoeVortex(unittest.TestCase):
    """This is a class with functions to test HorseshoeVortices."""

    def setUp(self):
        """Set up test fixtures for HorseshoeVortex tests."""
        # Create fixtures using horseshoe_vortex_fixtures module.
        self.basic_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )
        self.short_legs_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_short_legs_horseshoe_vortex_fixture()
        )
        self.long_legs_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_long_legs_horseshoe_vortex_fixture()
        )
        self.tilted_legs_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_tilted_legs_horseshoe_vortex_fixture()
        )
        self.wide_finite_leg_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_wide_finite_leg_horseshoe_vortex_fixture()
        )
        self.narrow_finite_leg_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_narrow_finite_leg_horseshoe_vortex_fixture()
        )
        self.zero_strength_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_zero_strength_horseshoe_vortex_fixture()
        )
        self.negative_strength_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_negative_strength_horseshoe_vortex_fixture()
        )
        self.high_strength_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_high_strength_horseshoe_vortex_fixture()
        )
        self.offset_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_offset_horseshoe_vortex_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test HorseshoeVortex initialization with valid parameters."""
        # Test that basic HorseshoeVortex initializes correctly.
        self.assertIsInstance(
            self.basic_horseshoe_vortex, _aerodynamics.HorseshoeVortex
        )
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.Frhvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.Flhvp_GP1_CgP1,
            np.array([0.0, -0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.leftLegVector_GP1,
            np.array([1.0, 0.0, 0.0], dtype=float),
        )
        self.assertEqual(self.basic_horseshoe_vortex.left_right_leg_lengths, 20.0)
        self.assertEqual(self.basic_horseshoe_vortex.strength, 1.0)

    def test_initialization_with_different_strengths(self):
        """Test HorseshoeVortex initialization with various strength values."""
        # Test with zero strength.
        self.assertEqual(self.zero_strength_horseshoe_vortex.strength, 0.0)

        # Test with negative strength.
        self.assertEqual(self.negative_strength_horseshoe_vortex.strength, -1.0)

        # Test with high positive strength.
        self.assertEqual(self.high_strength_horseshoe_vortex.strength, 100.0)

    def test_initialization_with_different_leg_lengths(self):
        """Test HorseshoeVortex initialization with various leg lengths."""
        # Test with short legs.
        self.assertEqual(self.short_legs_horseshoe_vortex.left_right_leg_lengths, 5.0)

        # Test with long legs.
        self.assertEqual(self.long_legs_horseshoe_vortex.left_right_leg_lengths, 100.0)

        # Test with basic legs.
        self.assertEqual(self.basic_horseshoe_vortex.left_right_leg_lengths, 20.0)

    def test_initialization_with_different_finite_leg_widths(self):
        """Test HorseshoeVortex initialization with various finite leg widths."""
        # Test wide finite leg.
        npt.assert_array_equal(
            self.wide_finite_leg_horseshoe_vortex.Frhvp_GP1_CgP1,
            np.array([0.0, 2.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.wide_finite_leg_horseshoe_vortex.Flhvp_GP1_CgP1,
            np.array([0.0, -2.5, 0.0], dtype=float),
        )

        # Test narrow finite leg.
        npt.assert_array_equal(
            self.narrow_finite_leg_horseshoe_vortex.Frhvp_GP1_CgP1,
            np.array([0.0, 0.05, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.narrow_finite_leg_horseshoe_vortex.Flhvp_GP1_CgP1,
            np.array([0.0, -0.05, 0.0], dtype=float),
        )

    def test_left_leg_vector_normalization(self):
        """Test that leftLegVector_GP1 is normalized to unit length."""
        # The tilted legs fixture uses [0.707, 0.0, 0.707] which should be normalized.
        left_leg_length = np.linalg.norm(
            self.tilted_legs_horseshoe_vortex.leftLegVector_GP1
        )
        npt.assert_almost_equal(left_leg_length, 1.0, decimal=10)

        # Basic fixture uses [1.0, 0.0, 0.0] which is already unit length.
        left_leg_length = np.linalg.norm(self.basic_horseshoe_vortex.leftLegVector_GP1)
        npt.assert_almost_equal(left_leg_length, 1.0, decimal=10)

    def test_back_points_calculation(self):
        """Test that back points (Brhvp_GP1_CgP1 and Blhvp_GP1_CgP1) are calculated
        correctly."""
        # For basic HorseshoeVortex, back points should be 20 meters in the
        # positive x-direction from front points.
        expected_Brhvp_GP1_CgP1 = np.array([20.0, 0.5, 0.0], dtype=float)
        expected_Blhvp_GP1_CgP1 = np.array([20.0, -0.5, 0.0], dtype=float)

        npt.assert_array_almost_equal(
            self.basic_horseshoe_vortex.Brhvp_GP1_CgP1,
            expected_Brhvp_GP1_CgP1,
            decimal=10,
        )
        npt.assert_array_almost_equal(
            self.basic_horseshoe_vortex.Blhvp_GP1_CgP1,
            expected_Blhvp_GP1_CgP1,
            decimal=10,
        )

    def test_back_points_with_different_leg_lengths(self):
        """Test back point calculation with different leg lengths."""
        # For short legs (5 meters).
        expected_Brhvp_GP1_CgP1 = np.array([5.0, 0.5, 0.0], dtype=float)
        expected_Blhvp_GP1_CgP1 = np.array([5.0, -0.5, 0.0], dtype=float)

        npt.assert_array_almost_equal(
            self.short_legs_horseshoe_vortex.Brhvp_GP1_CgP1,
            expected_Brhvp_GP1_CgP1,
            decimal=10,
        )
        npt.assert_array_almost_equal(
            self.short_legs_horseshoe_vortex.Blhvp_GP1_CgP1,
            expected_Blhvp_GP1_CgP1,
            decimal=10,
        )

        # For long legs (100 meters).
        expected_Brhvp_GP1_CgP1 = np.array([100.0, 0.5, 0.0], dtype=float)
        expected_Blhvp_GP1_CgP1 = np.array([100.0, -0.5, 0.0], dtype=float)

        npt.assert_array_almost_equal(
            self.long_legs_horseshoe_vortex.Brhvp_GP1_CgP1,
            expected_Brhvp_GP1_CgP1,
            decimal=10,
        )
        npt.assert_array_almost_equal(
            self.long_legs_horseshoe_vortex.Blhvp_GP1_CgP1,
            expected_Blhvp_GP1_CgP1,
            decimal=10,
        )

    def test_line_vortex_legs_created(self):
        """Test that HorseshoeVortex creates three LineVortex legs correctly."""
        # Check that all three legs are created.
        self.assertIsNotNone(self.basic_horseshoe_vortex.right_leg)
        self.assertIsNotNone(self.basic_horseshoe_vortex.finite_leg)
        self.assertIsNotNone(self.basic_horseshoe_vortex.left_leg)

        # Check that leg strengths match the HorseshoeVortex strength.
        self.assertEqual(
            self.basic_horseshoe_vortex.right_leg.strength,
            self.basic_horseshoe_vortex.strength,
        )
        self.assertEqual(
            self.basic_horseshoe_vortex.finite_leg.strength,
            self.basic_horseshoe_vortex.strength,
        )
        self.assertEqual(
            self.basic_horseshoe_vortex.left_leg.strength,
            self.basic_horseshoe_vortex.strength,
        )

    def test_line_vortex_legs_positions(self):
        """Test that LineVortex leg positions match HorseshoeVortex points."""
        # Test right leg (from back right to front right).
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.right_leg.Slvp_GP1_CgP1,
            self.basic_horseshoe_vortex.Brhvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.right_leg.Elvp_GP1_CgP1,
            self.basic_horseshoe_vortex.Frhvp_GP1_CgP1,
        )

        # Test finite leg (from front right to front left).
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.finite_leg.Slvp_GP1_CgP1,
            self.basic_horseshoe_vortex.Frhvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.finite_leg.Elvp_GP1_CgP1,
            self.basic_horseshoe_vortex.Flhvp_GP1_CgP1,
        )

        # Test left leg (from front left to back left).
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.left_leg.Slvp_GP1_CgP1,
            self.basic_horseshoe_vortex.Flhvp_GP1_CgP1,
        )
        npt.assert_array_equal(
            self.basic_horseshoe_vortex.left_leg.Elvp_GP1_CgP1,
            self.basic_horseshoe_vortex.Blhvp_GP1_CgP1,
        )

    def test_tilted_leg_vector(self):
        """Test HorseshoeVortex with tilted leg vector."""
        # The tilted legs fixture should have a normalized vector.
        npt.assert_almost_equal(
            np.linalg.norm(self.tilted_legs_horseshoe_vortex.leftLegVector_GP1),
            1.0,
            decimal=10,
        )

        # Back points should be displaced along the tilted direction.
        expected_displacement = (
            self.tilted_legs_horseshoe_vortex.leftLegVector_GP1
            * self.tilted_legs_horseshoe_vortex.left_right_leg_lengths
        )
        expected_Blhvp_GP1_CgP1 = (
            self.tilted_legs_horseshoe_vortex.Flhvp_GP1_CgP1 + expected_displacement
        )

        npt.assert_array_almost_equal(
            self.tilted_legs_horseshoe_vortex.Blhvp_GP1_CgP1,
            expected_Blhvp_GP1_CgP1,
            decimal=10,
        )

    def test_offset_horseshoe_vortex(self):
        """Test HorseshoeVortex positioned far from origin."""
        # Verify offset was applied correctly.
        expected_offset = np.array([10.0, 5.0, 3.0], dtype=float)
        npt.assert_array_almost_equal(
            self.offset_horseshoe_vortex.Frhvp_GP1_CgP1,
            np.array([0.0, 0.5, 0.0], dtype=float) + expected_offset,
            decimal=10,
        )
        npt.assert_array_almost_equal(
            self.offset_horseshoe_vortex.Flhvp_GP1_CgP1,
            np.array([0.0, -0.5, 0.0], dtype=float) + expected_offset,
            decimal=10,
        )

    def test_cache_invalidation_back_points_on_front_right_change(self):
        """Test that back right point cache is invalidated when front right changes."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access back right point to populate cache.
        _ = horseshoe_vortex.Brhvp_GP1_CgP1

        # Change front right point.
        new_front_right = np.array([0.5, 0.5, 0.0], dtype=float)
        horseshoe_vortex.Frhvp_GP1_CgP1 = new_front_right

        # Verify that back right point is recalculated.
        expected_back_right = new_front_right + (
            horseshoe_vortex.leftLegVector_GP1 * horseshoe_vortex.left_right_leg_lengths
        )
        npt.assert_array_almost_equal(
            horseshoe_vortex.Brhvp_GP1_CgP1,
            expected_back_right,
            decimal=10,
        )

    def test_cache_invalidation_back_points_on_front_left_change(self):
        """Test that back left point cache is invalidated when front left changes."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access back left point to populate cache.
        _ = horseshoe_vortex.Blhvp_GP1_CgP1

        # Change front left point.
        new_front_left = np.array([0.5, -0.5, 0.0], dtype=float)
        horseshoe_vortex.Flhvp_GP1_CgP1 = new_front_left

        # Verify that back left point is recalculated.
        expected_back_left = new_front_left + (
            horseshoe_vortex.leftLegVector_GP1 * horseshoe_vortex.left_right_leg_lengths
        )
        npt.assert_array_almost_equal(
            horseshoe_vortex.Blhvp_GP1_CgP1,
            expected_back_left,
            decimal=10,
        )

    def test_cache_invalidation_back_points_on_leg_vector_change(self):
        """Test that back points are invalidated when left leg vector changes."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access back points to populate cache.
        initial_back_right = horseshoe_vortex.Brhvp_GP1_CgP1.copy()
        initial_back_left = horseshoe_vortex.Blhvp_GP1_CgP1.copy()

        # Change leg vector.
        new_leg_vector = np.array([0.707, 0.0, 0.707], dtype=float)
        horseshoe_vortex.leftLegVector_GP1 = new_leg_vector

        # Verify that back points are recalculated (different from initial).
        self.assertFalse(
            np.allclose(horseshoe_vortex.Brhvp_GP1_CgP1, initial_back_right)
        )
        self.assertFalse(
            np.allclose(horseshoe_vortex.Blhvp_GP1_CgP1, initial_back_left)
        )

    def test_cache_invalidation_back_points_on_leg_length_change(self):
        """Test that back points are invalidated when leg lengths change."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access back points to populate cache.
        initial_back_right = horseshoe_vortex.Brhvp_GP1_CgP1.copy()

        # Change leg length.
        new_leg_length = 50.0
        horseshoe_vortex.left_right_leg_lengths = new_leg_length

        # Verify that back points are recalculated.
        expected_back_right = horseshoe_vortex.Frhvp_GP1_CgP1 + (
            horseshoe_vortex.leftLegVector_GP1 * new_leg_length
        )
        npt.assert_array_almost_equal(
            horseshoe_vortex.Brhvp_GP1_CgP1,
            expected_back_right,
            decimal=10,
        )
        self.assertFalse(
            np.allclose(horseshoe_vortex.Brhvp_GP1_CgP1, initial_back_right)
        )

    def test_leg_position_updates_on_front_right_change(self):
        """Test that LineVortex leg positions update when front right changes."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access legs to create them.
        _ = horseshoe_vortex.right_leg
        _ = horseshoe_vortex.finite_leg

        # Change front right point.
        new_front_right = np.array([0.25, 0.75, 0.1], dtype=float)
        horseshoe_vortex.Frhvp_GP1_CgP1 = new_front_right

        # Verify that right leg end point is updated.
        npt.assert_array_equal(
            horseshoe_vortex.right_leg.Elvp_GP1_CgP1,
            new_front_right,
        )

        # Verify that finite leg start point is updated.
        npt.assert_array_equal(
            horseshoe_vortex.finite_leg.Slvp_GP1_CgP1,
            new_front_right,
        )

    def test_leg_position_updates_on_front_left_change(self):
        """Test that LineVortex leg positions update when front left changes."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access legs to create them.
        _ = horseshoe_vortex.finite_leg
        _ = horseshoe_vortex.left_leg

        # Change front left point.
        new_front_left = np.array([0.25, -0.75, 0.1], dtype=float)
        horseshoe_vortex.Flhvp_GP1_CgP1 = new_front_left

        # Verify that finite leg end point is updated.
        npt.assert_array_equal(
            horseshoe_vortex.finite_leg.Elvp_GP1_CgP1,
            new_front_left,
        )

        # Verify that left leg start point is updated.
        npt.assert_array_equal(
            horseshoe_vortex.left_leg.Slvp_GP1_CgP1,
            new_front_left,
        )

    def test_leg_position_updates_on_leg_vector_change(self):
        """Test that LineVortex leg positions update when leg vector changes."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access legs to create them.
        _ = horseshoe_vortex.right_leg
        _ = horseshoe_vortex.left_leg

        # Change leg vector.
        new_leg_vector = np.array([0.0, 0.0, 1.0], dtype=float)
        horseshoe_vortex.leftLegVector_GP1 = new_leg_vector

        # Verify that right leg start point is updated to new back right.
        npt.assert_array_almost_equal(
            horseshoe_vortex.right_leg.Slvp_GP1_CgP1,
            horseshoe_vortex.Brhvp_GP1_CgP1,
            decimal=10,
        )

        # Verify that left leg end point is updated to new back left.
        npt.assert_array_almost_equal(
            horseshoe_vortex.left_leg.Elvp_GP1_CgP1,
            horseshoe_vortex.Blhvp_GP1_CgP1,
            decimal=10,
        )

    def test_leg_position_updates_on_leg_length_change(self):
        """Test that LineVortex leg positions update when leg lengths change."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access legs to create them.
        _ = horseshoe_vortex.right_leg
        _ = horseshoe_vortex.left_leg

        # Change leg length.
        new_leg_length = 50.0
        horseshoe_vortex.left_right_leg_lengths = new_leg_length

        # Verify that right leg start point is updated to new back right.
        npt.assert_array_almost_equal(
            horseshoe_vortex.right_leg.Slvp_GP1_CgP1,
            horseshoe_vortex.Brhvp_GP1_CgP1,
            decimal=10,
        )

        # Verify that left leg end point is updated to new back left.
        npt.assert_array_almost_equal(
            horseshoe_vortex.left_leg.Elvp_GP1_CgP1,
            horseshoe_vortex.Blhvp_GP1_CgP1,
            decimal=10,
        )

    def test_strength_propagation_to_legs(self):
        """Test that changing HorseshoeVortex strength propagates to all legs."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access all legs to create them.
        _ = horseshoe_vortex.right_leg
        _ = horseshoe_vortex.finite_leg
        _ = horseshoe_vortex.left_leg

        # Change strength.
        new_strength = 5.0
        horseshoe_vortex.strength = new_strength

        # Verify all legs have the new strength.
        self.assertEqual(horseshoe_vortex.right_leg.strength, new_strength)
        self.assertEqual(horseshoe_vortex.finite_leg.strength, new_strength)
        self.assertEqual(horseshoe_vortex.left_leg.strength, new_strength)

    def test_strength_propagation_to_uncreated_legs(self):
        """Test that changing strength before accessing legs still works."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Change strength before accessing legs.
        new_strength = 7.5
        horseshoe_vortex.strength = new_strength

        # Access legs now and verify they have the new strength.
        self.assertEqual(horseshoe_vortex.right_leg.strength, new_strength)
        self.assertEqual(horseshoe_vortex.finite_leg.strength, new_strength)
        self.assertEqual(horseshoe_vortex.left_leg.strength, new_strength)

    def test_multiple_front_point_updates(self):
        """Test that multiple sequential front point updates work correctly."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Update both front points.
        new_front_right = np.array([1.0, 1.0, 0.0], dtype=float)
        new_front_left = np.array([1.0, -1.0, 0.0], dtype=float)
        horseshoe_vortex.Frhvp_GP1_CgP1 = new_front_right
        horseshoe_vortex.Flhvp_GP1_CgP1 = new_front_left

        # Verify back points are recalculated.
        expected_back_right = new_front_right + (
            horseshoe_vortex.leftLegVector_GP1 * horseshoe_vortex.left_right_leg_lengths
        )
        expected_back_left = new_front_left + (
            horseshoe_vortex.leftLegVector_GP1 * horseshoe_vortex.left_right_leg_lengths
        )

        npt.assert_array_almost_equal(
            horseshoe_vortex.Brhvp_GP1_CgP1,
            expected_back_right,
            decimal=10,
        )
        npt.assert_array_almost_equal(
            horseshoe_vortex.Blhvp_GP1_CgP1,
            expected_back_left,
            decimal=10,
        )

    def test_strength_zero_to_nonzero(self):
        """Test changing strength from zero to nonzero."""
        # Create a zero strength fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_zero_strength_horseshoe_vortex_fixture()
        )

        # Access legs.
        _ = horseshoe_vortex.right_leg

        # Change to nonzero.
        horseshoe_vortex.strength = 10.0

        # Verify.
        self.assertEqual(horseshoe_vortex.strength, 10.0)
        self.assertEqual(horseshoe_vortex.right_leg.strength, 10.0)

    def test_line_vortex_legs_are_line_vortex_type(self):
        """Test that LineVortex legs are of the correct type."""
        # Verify all legs are LineVortex instances.
        self.assertIsInstance(
            self.basic_horseshoe_vortex.right_leg, _aerodynamics._LineVortex
        )
        self.assertIsInstance(
            self.basic_horseshoe_vortex.finite_leg, _aerodynamics._LineVortex
        )
        self.assertIsInstance(
            self.basic_horseshoe_vortex.left_leg, _aerodynamics._LineVortex
        )


if __name__ == "__main__":
    unittest.main()
