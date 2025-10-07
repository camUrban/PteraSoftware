"""This module contains a class to test Airfoils."""

import unittest
import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestAirfoil(unittest.TestCase):
    """This class contains unit tests for the Airfoil class."""

    def setUp(self):
        """Set up test fixtures for Airfoil tests."""
        # Create fixtures for different Airfoil types
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()
        self.naca2412_airfoil = geometry_fixtures.make_naca2412_airfoil_fixture()
        self.custom_outline_airfoil = (
            geometry_fixtures.make_custom_outline_airfoil_fixture()
        )
        self.resampled_airfoil = geometry_fixtures.make_resampled_airfoil_fixture()
        self.non_resampled_airfoil = (
            geometry_fixtures.make_non_resampled_airfoil_fixture()
        )
        self.named_airfoil = geometry_fixtures.make_named_airfoil_fixture()

    def test_initialization_naca_airfoils(self):
        """Test Airfoil initialization with NACA airfoil names."""
        # Test NACA 0012 initialization
        self.assertEqual(self.naca0012_airfoil.name, "naca0012")
        self.assertTrue(self.naca0012_airfoil.resample)
        self.assertEqual(self.naca0012_airfoil.n_points_per_side, 400)
        self.assertIsNotNone(self.naca0012_airfoil.outline_A_lp)

        # Test NACA 2412 initialization
        self.assertEqual(self.naca2412_airfoil.name, "naca2412")
        self.assertTrue(self.naca2412_airfoil.resample)
        self.assertEqual(self.naca2412_airfoil.n_points_per_side, 400)
        self.assertIsNotNone(self.naca2412_airfoil.outline_A_lp)

    def test_initialization_custom_outline(self):
        """Test Airfoil initialization with custom outline."""
        self.assertEqual(self.custom_outline_airfoil.name, "Custom Test Airfoil")
        self.assertFalse(self.custom_outline_airfoil.resample)
        self.assertIsInstance(self.custom_outline_airfoil.n_points_per_side, int)
        self.assertIsNotNone(self.custom_outline_airfoil.outline_A_lp)

    def test_initialization_resampling_parameters(self):
        """Test Airfoil initialization with different resampling parameters."""
        # Test resampled airfoil
        self.assertTrue(self.resampled_airfoil.resample)
        self.assertEqual(self.resampled_airfoil.n_points_per_side, 100)

        # Test non-resampled airfoil
        self.assertFalse(self.non_resampled_airfoil.resample)
        self.assertIsInstance(self.non_resampled_airfoil.n_points_per_side, int)

    def test_outline_shape_and_bounds(self):
        """Test that airfoil outline has correct shape and coordinate bounds."""
        for airfoil in [
            self.naca0012_airfoil,
            self.naca2412_airfoil,
            self.custom_outline_airfoil,
            self.resampled_airfoil,
            self.non_resampled_airfoil,
        ]:
            # Check outline shape
            self.assertEqual(len(airfoil.outline_A_lp.shape), 2)
            self.assertEqual(airfoil.outline_A_lp.shape[1], 2)

            # Check x-coordinates are roughly within the [0, 1] range
            x_coords = airfoil.outline_A_lp[:, 0]
            self.assertTrue(np.all(x_coords >= 0.0 - 0.01))
            self.assertTrue(np.all(x_coords <= 1.0 + 0.01))

            # Check that outline contains multiple points
            self.assertGreater(airfoil.outline_A_lp.shape[0], 3)

    def test_airfoil_closure(self):
        """Test that airfoil outlines form approximately closed loops."""
        for airfoil in [
            self.naca0012_airfoil,
            self.naca2412_airfoil,
            self.resampled_airfoil,
        ]:
            outline = airfoil.outline_A_lp
            # Check that first and last points are close (closed loop)
            npt.assert_allclose(outline[0], outline[-1], rtol=2e-2, atol=2e-2)

    def test_leading_edge_location(self):
        """Test that leading edge is properly located at approximately x=0."""
        for airfoil in [
            self.naca0012_airfoil,
            self.naca2412_airfoil,
            self.resampled_airfoil,
        ]:
            outline = airfoil.outline_A_lp
            x_coords = outline[:, 0]
            # Check that minimum x-coordinate is approximately 0
            self.assertAlmostEqual(np.min(x_coords), 0.0, places=2)

    def test_trailing_edge_location(self):
        """Test that trailing edge is properly located at approximately x=1."""
        for airfoil in [
            self.naca0012_airfoil,
            self.naca2412_airfoil,
            self.resampled_airfoil,
        ]:
            outline = airfoil.outline_A_lp
            x_coords = outline[:, 0]
            # Check that maximum x-coordinate is approximately 1
            self.assertAlmostEqual(np.max(x_coords), 1.0, places=2)

    def test_symmetric_airfoil_properties(self):
        """Test properties specific to symmetric airfoils (NACA 0012)."""
        outline = self.naca0012_airfoil.outline_A_lp
        y_coords = outline[:, 1]

        # For a symmetric airfoil, expect roughly equal positive and negative y values
        positive_y = y_coords[y_coords > 0]
        negative_y = y_coords[y_coords < 0]
        self.assertGreater(len(positive_y), 0)
        self.assertGreater(len(negative_y), 0)

    def test_mcl_A_lp_attribute(self):
        """Test that mcl_A_lp attribute is properly set."""
        for airfoil in [
            self.naca0012_airfoil,
            self.naca2412_airfoil,
            self.resampled_airfoil,
        ]:
            self.assertIsNotNone(airfoil.mcl_A_lp)
            self.assertEqual(len(airfoil.mcl_A_lp.shape), 2)
            self.assertEqual(airfoil.mcl_A_lp.shape[1], 2)

            # MCL x-values should span approximately [0, 1]
            mclX_A_lp = airfoil.mcl_A_lp[:, 0]
            self.assertAlmostEqual(np.min(mclX_A_lp), 0.0, places=2)
            self.assertAlmostEqual(np.max(mclX_A_lp), 1.0, places=2)

    def test_get_resampled_mcl_method(self):
        """Test the get_resampled_mcl method."""
        num_points = 50
        mcl_fractions = np.linspace(0, 1, num_points)
        mcl_A_lp = self.naca2412_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)
        self.assertEqual(mcl_A_lp.shape[0], num_points)

        # Check x-values are still approximately in the range [0, 1]
        mclX_A_lp = mcl_A_lp[:, 0]
        self.assertAlmostEqual(np.min(mclX_A_lp), 0.0, places=2)
        self.assertAlmostEqual(np.max(mclX_A_lp), 1.0, places=2)

    def test_add_control_surface_method(self):
        """Test the add_control_surface method."""
        # Test adding control surface at 75% chord
        hinge_point = 0.75
        deflection = 5.0

        modified_airfoil = self.naca0012_airfoil.add_control_surface(
            deflection=deflection,
            hinge_point=hinge_point,
        )

        # Check that we get a new Airfoil
        self.assertIsInstance(modified_airfoil, ps.geometry.airfoil.Airfoil)
        self.assertIsNot(modified_airfoil, self.naca0012_airfoil)

        # Check that outline is modified
        self.assertIsNotNone(modified_airfoil.outline_A_lp)

        # Control surface should change the trailing edge shape
        original_outline = self.naca0012_airfoil.outline_A_lp
        modified_outline = modified_airfoil.outline_A_lp
        self.assertFalse(np.allclose(original_outline, modified_outline))

    def test_parameter_validation_invalid_inputs(self):
        """Test that invalid parameters raise appropriate errors."""
        # Test outlines with x-values significantly less than 0.0.
        invalid_outline = np.array([[-0.1, 0.0], [0.5, 0.1], [1.0, 0.0], [0.5, -0.1]])
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(
                name="Invalid", outline_A_lp=invalid_outline, resample=False
            )

        # Test outlines with x-values significantly greater than 1.0.
        invalid_outline2 = np.array([[0.0, 0.0], [0.5, 0.1], [1.1, 0.0], [0.5, -0.1]])
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(
                name="Invalid2", outline_A_lp=invalid_outline2, resample=False
            )

    def test_naca_airfoil_thickness(self):
        """Test that the generated NACA0012 and NACA2412 Airfoils have approximately
        the correct maximum thickness."""
        for airfoil in [self.naca0012_airfoil, self.naca2412_airfoil]:
            outline_A_lp = airfoil.outline_A_lp
            outlineY_A_lp = outline_A_lp[:, 1]

            # Find upper and lower surfaces roughly
            # For properly ordered airfoils, we expect some variation in y
            max_thickness = np.max(outlineY_A_lp) - np.min(outlineY_A_lp)
            self.assertAlmostEqual(max_thickness, 0.12, places=2)

    def test_initialization_named_airfoil(self):
        """Test Airfoil initialization with a named airfoil from the _airfoils data
        directory."""
        # Test that the named airfoil loads correctly from the _airfoils directory
        self.assertEqual(self.named_airfoil.name, "a18")
        self.assertTrue(self.named_airfoil.resample)
        self.assertEqual(self.named_airfoil.n_points_per_side, 400)
        self.assertIsNotNone(self.named_airfoil.outline_A_lp)

        # Test that the outline has the correct shape and bounds
        self.assertEqual(len(self.named_airfoil.outline_A_lp.shape), 2)
        self.assertEqual(self.named_airfoil.outline_A_lp.shape[1], 2)

        # Check x-coordinates are roughly within the [0, 1] range
        x_coords = self.named_airfoil.outline_A_lp[:, 0]
        self.assertTrue(np.all(x_coords >= 0.0 - 0.01))
        self.assertTrue(np.all(x_coords <= 1.0 + 0.01))

        # Check that outline contains multiple points
        self.assertGreater(self.named_airfoil.outline_A_lp.shape[0], 3)

        # Test that the MCL is populated
        self.assertIsNotNone(self.named_airfoil.mcl_A_lp)
        self.assertEqual(len(self.named_airfoil.mcl_A_lp.shape), 2)
        self.assertEqual(self.named_airfoil.mcl_A_lp.shape[1], 2)


if __name__ == "__main__":
    unittest.main()
