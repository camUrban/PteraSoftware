"""This module contains classes to test HorseshoeVortices."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _vortices

# noinspection PyProtectedMember
from pterasoftware._vortices import _line_vortex
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
            self.basic_horseshoe_vortex, _vortices.horseshoe_vortex.HorseshoeVortex
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

    def test_tilted_leg_vector_back_point_displacement(self):
        """Test that back points are displaced along a tilted leg vector direction."""
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

    def test_line_vortex_legs_are_line_vortex_type(self):
        """Test that LineVortex legs are of the correct type."""
        # Verify all legs are LineVortex instances.
        self.assertIsInstance(
            self.basic_horseshoe_vortex.right_leg, _line_vortex.LineVortex
        )
        self.assertIsInstance(
            self.basic_horseshoe_vortex.finite_leg, _line_vortex.LineVortex
        )
        self.assertIsInstance(
            self.basic_horseshoe_vortex.left_leg, _line_vortex.LineVortex
        )


class TestHorseshoeVortexImmutability(unittest.TestCase):
    """Tests for HorseshoeVortex attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.basic_horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

    def test_immutable_Frhvp_GP1_CgP1_property(self):
        """Test that Frhvp_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.Frhvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Frhvp_GP1_CgP1_array_read_only(self):
        """Test that Frhvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_horseshoe_vortex.Frhvp_GP1_CgP1[0] = 999.0

    def test_immutable_Flhvp_GP1_CgP1_property(self):
        """Test that Flhvp_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.Flhvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Flhvp_GP1_CgP1_array_read_only(self):
        """Test that Flhvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_horseshoe_vortex.Flhvp_GP1_CgP1[0] = 999.0

    def test_immutable_leftLegVector_GP1_property(self):
        """Test that leftLegVector_GP1 property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.leftLegVector_GP1 = np.array([0.0, 1.0, 0.0])

    def test_immutable_leftLegVector_GP1_array_read_only(self):
        """Test that leftLegVector_GP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_horseshoe_vortex.leftLegVector_GP1[0] = 999.0

    def test_immutable_left_right_leg_lengths_property(self):
        """Test that left_right_leg_lengths property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.left_right_leg_lengths = 50.0

    def test_immutable_Brhvp_GP1_CgP1_property(self):
        """Test that Brhvp_GP1_CgP1 property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_horseshoe_vortex.Brhvp_GP1_CgP1
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.Brhvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Brhvp_GP1_CgP1_array_read_only(self):
        """Test that Brhvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_horseshoe_vortex.Brhvp_GP1_CgP1[0] = 999.0

    def test_immutable_Blhvp_GP1_CgP1_property(self):
        """Test that Blhvp_GP1_CgP1 property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_horseshoe_vortex.Blhvp_GP1_CgP1
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.Blhvp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])

    def test_immutable_Blhvp_GP1_CgP1_array_read_only(self):
        """Test that Blhvp_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_horseshoe_vortex.Blhvp_GP1_CgP1[0] = 999.0

    def test_immutable_right_leg_property(self):
        """Test that right_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_horseshoe_vortex.right_leg
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.right_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_immutable_finite_leg_property(self):
        """Test that finite_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_horseshoe_vortex.finite_leg
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.finite_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_immutable_left_leg_property(self):
        """Test that left_leg property is read only."""
        # Access to trigger lazy computation.
        _ = self.basic_horseshoe_vortex.left_leg
        with self.assertRaises(AttributeError):
            self.basic_horseshoe_vortex.left_leg = _line_vortex.LineVortex(
                Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
                Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
                strength=1.0,
            )

    def test_mutable_strength_property(self):
        """Test that strength property is mutable."""
        # Strength should be modifiable.
        self.basic_horseshoe_vortex.strength = 5.0
        self.assertEqual(self.basic_horseshoe_vortex.strength, 5.0)

        self.basic_horseshoe_vortex.strength = -3.0
        self.assertEqual(self.basic_horseshoe_vortex.strength, -3.0)


class TestHorseshoeVortexLazyCaching(unittest.TestCase):
    """Tests for HorseshoeVortex lazy caching behavior."""

    def test_back_points_lazy_evaluation(self):
        """Test that back points are lazily evaluated."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Verify that _Brhvp_GP1_CgP1 is None before access.
        self.assertIsNone(horseshoe_vortex._Brhvp_GP1_CgP1)
        self.assertIsNone(horseshoe_vortex._Blhvp_GP1_CgP1)

        # Access the properties.
        _ = horseshoe_vortex.Brhvp_GP1_CgP1
        _ = horseshoe_vortex.Blhvp_GP1_CgP1

        # Verify that the caches are now populated.
        self.assertIsNotNone(horseshoe_vortex._Brhvp_GP1_CgP1)
        self.assertIsNotNone(horseshoe_vortex._Blhvp_GP1_CgP1)

    def test_back_points_caching(self):
        """Test that back points return the same cached object on repeated access."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access back right point twice.
        back_right_first = horseshoe_vortex.Brhvp_GP1_CgP1
        back_right_second = horseshoe_vortex.Brhvp_GP1_CgP1

        # Verify they are the same object (not just equal).
        self.assertIs(back_right_first, back_right_second)

        # Access back left point twice.
        back_left_first = horseshoe_vortex.Blhvp_GP1_CgP1
        back_left_second = horseshoe_vortex.Blhvp_GP1_CgP1

        # Verify they are the same object.
        self.assertIs(back_left_first, back_left_second)

    def test_line_vortex_legs_lazy_evaluation(self):
        """Test that LineVortex legs are lazily evaluated."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Verify that leg caches are None before access.
        self.assertIsNone(horseshoe_vortex._right_leg)
        self.assertIsNone(horseshoe_vortex._finite_leg)
        self.assertIsNone(horseshoe_vortex._left_leg)

        # Access the properties.
        _ = horseshoe_vortex.right_leg
        _ = horseshoe_vortex.finite_leg
        _ = horseshoe_vortex.left_leg

        # Verify that the caches are now populated.
        self.assertIsNotNone(horseshoe_vortex._right_leg)
        self.assertIsNotNone(horseshoe_vortex._finite_leg)
        self.assertIsNotNone(horseshoe_vortex._left_leg)

    def test_line_vortex_legs_caching(self):
        """Test that LineVortex legs return the same cached object on repeated access."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Access right leg twice.
        right_leg_first = horseshoe_vortex.right_leg
        right_leg_second = horseshoe_vortex.right_leg

        # Verify they are the same object.
        self.assertIs(right_leg_first, right_leg_second)

        # Access finite leg twice.
        finite_leg_first = horseshoe_vortex.finite_leg
        finite_leg_second = horseshoe_vortex.finite_leg

        # Verify they are the same object.
        self.assertIs(finite_leg_first, finite_leg_second)

        # Access left leg twice.
        left_leg_first = horseshoe_vortex.left_leg
        left_leg_second = horseshoe_vortex.left_leg

        # Verify they are the same object.
        self.assertIs(left_leg_first, left_leg_second)

    def test_back_points_computed_before_legs_when_accessing_legs(self):
        """Test that accessing legs triggers computation of back points."""
        # Create a fresh fixture.
        horseshoe_vortex = (
            horseshoe_vortex_fixtures.make_basic_horseshoe_vortex_fixture()
        )

        # Verify caches start empty.
        self.assertIsNone(horseshoe_vortex._Brhvp_GP1_CgP1)
        self.assertIsNone(horseshoe_vortex._Blhvp_GP1_CgP1)

        # Access right leg (which needs back right point).
        _ = horseshoe_vortex.right_leg

        # Verify back right point is now computed.
        self.assertIsNotNone(horseshoe_vortex._Brhvp_GP1_CgP1)

        # Back left point should still be None since it wasn't needed.
        self.assertIsNone(horseshoe_vortex._Blhvp_GP1_CgP1)

        # Access left leg (which needs back left point).
        _ = horseshoe_vortex.left_leg

        # Verify back left point is now computed.
        self.assertIsNotNone(horseshoe_vortex._Blhvp_GP1_CgP1)


if __name__ == "__main__":
    unittest.main()
