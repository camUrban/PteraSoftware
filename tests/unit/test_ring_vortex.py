"""This module contains a class to test RingVortices."""

import unittest
import numpy as np
import numpy.testing as npt

import pterasoftware as ps
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

    def test_initialization_valid_parameters(self):
        """Test RingVortex initialization with valid parameters."""
        # Test that basic RingVortex initializes correctly.
        self.assertIsInstance(self.basic_ring_vortex, ps.aerodynamics.RingVortex)
        npt.assert_array_equal(
            self.basic_ring_vortex.Frrvp_G_Cg, np.array([0.0, 0.5, 0.0], dtype=float)
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.Flrvp_G_Cg, np.array([0.0, -0.5, 0.0], dtype=float)
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.Blrvp_G_Cg, np.array([1.0, -0.5, 0.0], dtype=float)
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.Brrvp_G_Cg, np.array([1.0, 0.5, 0.0], dtype=float)
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
            self.unit_square_ring_vortex.Frrvp_G_Cg,
            np.array([0.5, 0.5, 0.0], dtype=float),
        )
        npt.assert_array_equal(
            self.unit_square_ring_vortex.Flrvp_G_Cg,
            np.array([0.5, -0.5, 0.0], dtype=float),
        )

        # Test rectangular geometry.
        npt.assert_array_equal(
            self.rectangular_ring_vortex.Brrvp_G_Cg,
            np.array([2.0, 0.5, 0.0], dtype=float),
        )

        # Test tilted geometry.
        npt.assert_array_equal(
            self.tilted_ring_vortex.Frrvp_G_Cg,
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
            self.basic_ring_vortex.front_leg.Slvp_G_Cg,
            self.basic_ring_vortex.Frrvp_G_Cg,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.front_leg.Elvp_G_Cg,
            self.basic_ring_vortex.Flrvp_G_Cg,
        )

        # Test left leg (from front left to back left).
        npt.assert_array_equal(
            self.basic_ring_vortex.left_leg.Slvp_G_Cg,
            self.basic_ring_vortex.Flrvp_G_Cg,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.left_leg.Elvp_G_Cg,
            self.basic_ring_vortex.Blrvp_G_Cg,
        )

        # Test back leg (from back left to back right).
        npt.assert_array_equal(
            self.basic_ring_vortex.back_leg.Slvp_G_Cg,
            self.basic_ring_vortex.Blrvp_G_Cg,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.back_leg.Elvp_G_Cg,
            self.basic_ring_vortex.Brrvp_G_Cg,
        )

        # Test right leg (from back right to front right).
        npt.assert_array_equal(
            self.basic_ring_vortex.right_leg.Slvp_G_Cg,
            self.basic_ring_vortex.Brrvp_G_Cg,
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.right_leg.Elvp_G_Cg,
            self.basic_ring_vortex.Frrvp_G_Cg,
        )

    def test_centroid_calculation(self):
        """Test that RingVortex centroid is calculated correctly."""
        # For a unit square centered at origin, centroid should be at origin.
        expected_centroid = np.array([0.0, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.unit_square_ring_vortex.Crvp_G_Cg, expected_centroid, decimal=10
        )

        # For the basic RingVortex fixture, the centroid should be at [0.5, 0.0,
        # 0.0] (in geometry axes, relative to the CG).
        expected_centroid = np.array([0.5, 0.0, 0.0], dtype=float)
        npt.assert_array_almost_equal(
            self.basic_ring_vortex.Crvp_G_Cg, expected_centroid, decimal=10
        )

    def test_initial_age_is_zero(self):
        """Test that RingVortex age is initialized to zero."""
        self.assertEqual(self.basic_ring_vortex.age, 0)
        self.assertEqual(self.unit_square_ring_vortex.age, 0)
        self.assertEqual(self.tilted_ring_vortex.age, 0)

    def test_update_strength_method(self):
        """Test the update_strength method."""
        # Get initial strength.
        initial_strength = self.basic_ring_vortex.strength

        # Update to new strength.
        new_strength = 5.0
        self.basic_ring_vortex.update_strength(new_strength)

        # Verify RingVortex strength was updated.
        self.assertEqual(self.basic_ring_vortex.strength, new_strength)
        self.assertNotEqual(self.basic_ring_vortex.strength, initial_strength)

        # Verify all leg strengths were updated.
        self.assertEqual(self.basic_ring_vortex.front_leg.strength, new_strength)
        self.assertEqual(self.basic_ring_vortex.left_leg.strength, new_strength)
        self.assertEqual(self.basic_ring_vortex.back_leg.strength, new_strength)
        self.assertEqual(self.basic_ring_vortex.right_leg.strength, new_strength)

    def test_update_strength_to_zero(self):
        """Test updating strength to zero."""
        self.basic_ring_vortex.update_strength(0.0)
        self.assertEqual(self.basic_ring_vortex.strength, 0.0)
        self.assertEqual(self.basic_ring_vortex.front_leg.strength, 0.0)
        self.assertEqual(self.basic_ring_vortex.left_leg.strength, 0.0)
        self.assertEqual(self.basic_ring_vortex.back_leg.strength, 0.0)
        self.assertEqual(self.basic_ring_vortex.right_leg.strength, 0.0)

    def test_update_strength_to_negative(self):
        """Test updating strength to negative value."""
        self.basic_ring_vortex.update_strength(-3.5)
        self.assertEqual(self.basic_ring_vortex.strength, -3.5)
        self.assertEqual(self.basic_ring_vortex.front_leg.strength, -3.5)
        self.assertEqual(self.basic_ring_vortex.left_leg.strength, -3.5)
        self.assertEqual(self.basic_ring_vortex.back_leg.strength, -3.5)
        self.assertEqual(self.basic_ring_vortex.right_leg.strength, -3.5)

    def test_update_position_method(self):
        """Test the update_position method."""
        # Define new positions.
        new_Frrvp_G_Cg = np.array([1.0, 1.0, 0.5], dtype=float)
        new_Flrvp_G_Cg = np.array([1.0, -1.0, 0.5], dtype=float)
        new_Blrvp_G_Cg = np.array([2.0, -1.0, 0.5], dtype=float)
        new_Brrvp_G_Cg = np.array([2.0, 1.0, 0.5], dtype=float)

        # Update position.
        self.basic_ring_vortex.update_position(
            new_Frrvp_G_Cg, new_Flrvp_G_Cg, new_Blrvp_G_Cg, new_Brrvp_G_Cg
        )

        # Verify corner positions were updated.
        npt.assert_array_equal(self.basic_ring_vortex.Frrvp_G_Cg, new_Frrvp_G_Cg)
        npt.assert_array_equal(self.basic_ring_vortex.Flrvp_G_Cg, new_Flrvp_G_Cg)
        npt.assert_array_equal(self.basic_ring_vortex.Blrvp_G_Cg, new_Blrvp_G_Cg)
        npt.assert_array_equal(self.basic_ring_vortex.Brrvp_G_Cg, new_Brrvp_G_Cg)

    def test_update_position_updates_line_vortices(self):
        """Test that update_position updates LineVortex legs."""
        # Define new positions.
        new_Frrvp_G_Cg = np.array([1.0, 1.0, 0.5], dtype=float)
        new_Flrvp_G_Cg = np.array([1.0, -1.0, 0.5], dtype=float)
        new_Blrvp_G_Cg = np.array([2.0, -1.0, 0.5], dtype=float)
        new_Brrvp_G_Cg = np.array([2.0, 1.0, 0.5], dtype=float)

        # Update position.
        self.basic_ring_vortex.update_position(
            new_Frrvp_G_Cg, new_Flrvp_G_Cg, new_Blrvp_G_Cg, new_Brrvp_G_Cg
        )

        # Verify front leg was updated.
        npt.assert_array_equal(
            self.basic_ring_vortex.front_leg.Slvp_G_Cg, new_Frrvp_G_Cg
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.front_leg.Elvp_G_Cg, new_Flrvp_G_Cg
        )

        # Verify left leg was updated.
        npt.assert_array_equal(
            self.basic_ring_vortex.left_leg.Slvp_G_Cg, new_Flrvp_G_Cg
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.left_leg.Elvp_G_Cg, new_Blrvp_G_Cg
        )

        # Verify back leg was updated.
        npt.assert_array_equal(
            self.basic_ring_vortex.back_leg.Slvp_G_Cg, new_Blrvp_G_Cg
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.back_leg.Elvp_G_Cg, new_Brrvp_G_Cg
        )

        # Verify right leg was updated.
        npt.assert_array_equal(
            self.basic_ring_vortex.right_leg.Slvp_G_Cg, new_Brrvp_G_Cg
        )
        npt.assert_array_equal(
            self.basic_ring_vortex.right_leg.Elvp_G_Cg, new_Frrvp_G_Cg
        )

    def test_update_position_updates_centroid(self):
        """Test that update_position updates the centroid."""
        # Define new positions that form a square centered at [1.5, 0.0, 0.5] (in
        # geometry axes, relative to the CG).
        new_Frrvp_G_Cg = np.array([1.0, 1.0, 0.5], dtype=float)
        new_Flrvp_G_Cg = np.array([1.0, -1.0, 0.5], dtype=float)
        new_Blrvp_G_Cg = np.array([2.0, -1.0, 0.5], dtype=float)
        new_Brrvp_G_Cg = np.array([2.0, 1.0, 0.5], dtype=float)

        # Update position.
        self.basic_ring_vortex.update_position(
            new_Frrvp_G_Cg, new_Flrvp_G_Cg, new_Blrvp_G_Cg, new_Brrvp_G_Cg
        )

        # Expected centroid.
        expected_centroid = np.array([1.5, 0.0, 0.5], dtype=float)
        npt.assert_array_almost_equal(
            self.basic_ring_vortex.Crvp_G_Cg, expected_centroid, decimal=10
        )

    def test_update_position_preserves_strength(self):
        """Test that update_position preserves strength."""
        original_strength = self.basic_ring_vortex.strength

        # Define new positions.
        new_Frrvp_G_Cg = np.array([1.0, 1.0, 0.5], dtype=float)
        new_Flrvp_G_Cg = np.array([1.0, -1.0, 0.5], dtype=float)
        new_Blrvp_G_Cg = np.array([2.0, -1.0, 0.5], dtype=float)
        new_Brrvp_G_Cg = np.array([2.0, 1.0, 0.5], dtype=float)

        # Update position.
        self.basic_ring_vortex.update_position(
            new_Frrvp_G_Cg, new_Flrvp_G_Cg, new_Blrvp_G_Cg, new_Brrvp_G_Cg
        )

        # Verify strength was preserved.
        self.assertEqual(self.basic_ring_vortex.strength, original_strength)
        self.assertEqual(self.basic_ring_vortex.front_leg.strength, original_strength)
        self.assertEqual(self.basic_ring_vortex.left_leg.strength, original_strength)
        self.assertEqual(self.basic_ring_vortex.back_leg.strength, original_strength)
        self.assertEqual(self.basic_ring_vortex.right_leg.strength, original_strength)

    def test_parameter_validation_invalid_types(self):
        """Test parameter validation with invalid types."""
        # Test with string instead of array.
        with self.assertRaises((TypeError, ValueError)):
            ps.aerodynamics.RingVortex(
                Frrvp_G_Cg="not_an_array",
                Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
                Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
                Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
                strength=1.0,
            )

        # Test with wrong-sized ndarray.
        with self.assertRaises((TypeError, ValueError)):
            ps.aerodynamics.RingVortex(
                Frrvp_G_Cg=np.array([0.0, 0.5], dtype=float),  # Only 2 elements
                Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
                Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
                Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
                strength=1.0,
            )

        # Test with string strength.
        with self.assertRaises((TypeError, ValueError)):
            ps.aerodynamics.RingVortex(
                Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
                Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
                Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
                Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
                strength="not_a_number",
            )

    def test_none_strength_handling(self):
        """Test that None strength is handled correctly."""
        ring_vortex = ps.aerodynamics.RingVortex(
            Frrvp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
            Flrvp_G_Cg=np.array([0.0, -0.5, 0.0], dtype=float),
            Blrvp_G_Cg=np.array([1.0, -0.5, 0.0], dtype=float),
            Brrvp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
            strength=None,
        )

        # Verify strength is None.
        self.assertIsNone(ring_vortex.strength)

        # Verify leg strengths are also None.
        self.assertIsNone(ring_vortex.front_leg.strength)
        self.assertIsNone(ring_vortex.left_leg.strength)
        self.assertIsNone(ring_vortex.back_leg.strength)
        self.assertIsNone(ring_vortex.right_leg.strength)

        # Clean up.
        del ring_vortex

    def test_very_small_vortex(self):
        """Test RingVortex with very small dimensions."""
        # Should handle very small RingVortices without issues.
        self.assertIsNotNone(self.small_ring_vortex.Crvp_G_Cg)
        self.assertIsNotNone(self.small_ring_vortex.front_leg)

    def test_offset_vortex(self):
        """Test RingVortex positioned far from origin."""
        # Verify offset was applied correctly.
        expected_offset = np.array([5.0, 3.0, 2.0], dtype=float)
        npt.assert_array_almost_equal(
            self.offset_ring_vortex.Frrvp_G_Cg,
            np.array([0.0, 0.5, 0.0], dtype=float) + expected_offset,
            decimal=10,
        )


if __name__ == "__main__":
    unittest.main()
