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

    def test_update_strength_method(self):
        """Test the update_strength method."""
        # Get initial strength.
        initial_strength = self.basic_horseshoe_vortex.strength

        # Update to new strength.
        new_strength = 5.0
        self.basic_horseshoe_vortex.update_strength(new_strength)

        # Verify HorseshoeVortex strength was updated.
        self.assertEqual(self.basic_horseshoe_vortex.strength, new_strength)
        self.assertNotEqual(self.basic_horseshoe_vortex.strength, initial_strength)

        # Verify all leg strengths were updated.
        self.assertEqual(self.basic_horseshoe_vortex.right_leg.strength, new_strength)
        self.assertEqual(self.basic_horseshoe_vortex.finite_leg.strength, new_strength)
        self.assertEqual(self.basic_horseshoe_vortex.left_leg.strength, new_strength)

    def test_update_strength_to_zero(self):
        """Test updating strength to zero."""
        self.basic_horseshoe_vortex.update_strength(0.0)
        self.assertEqual(self.basic_horseshoe_vortex.strength, 0.0)
        self.assertEqual(self.basic_horseshoe_vortex.right_leg.strength, 0.0)
        self.assertEqual(self.basic_horseshoe_vortex.finite_leg.strength, 0.0)
        self.assertEqual(self.basic_horseshoe_vortex.left_leg.strength, 0.0)

    def test_update_strength_to_negative(self):
        """Test updating strength to negative value."""
        self.basic_horseshoe_vortex.update_strength(-3.5)
        self.assertEqual(self.basic_horseshoe_vortex.strength, -3.5)
        self.assertEqual(self.basic_horseshoe_vortex.right_leg.strength, -3.5)
        self.assertEqual(self.basic_horseshoe_vortex.finite_leg.strength, -3.5)
        self.assertEqual(self.basic_horseshoe_vortex.left_leg.strength, -3.5)

    def test_none_strength_handling(self):
        """Test that None strength is handled correctly."""
        horseshoe_vortex = _aerodynamics.HorseshoeVortex(
            Frhvp_GP1_CgP1=np.array([0.0, 0.5, 0.0], dtype=float),
            Flhvp_GP1_CgP1=np.array([0.0, -0.5, 0.0], dtype=float),
            leftLegVector_GP1=np.array([1.0, 0.0, 0.0], dtype=float),
            left_right_leg_lengths=20.0,
            strength=None,
        )

        # Verify strength is None.
        self.assertIsNone(horseshoe_vortex.strength)

        # Verify leg strengths are also None.
        self.assertIsNone(horseshoe_vortex.right_leg.strength)
        self.assertIsNone(horseshoe_vortex.finite_leg.strength)
        self.assertIsNone(horseshoe_vortex.left_leg.strength)

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


if __name__ == "__main__":
    unittest.main()
