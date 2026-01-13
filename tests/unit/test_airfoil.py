"""This module contains a class to test Airfoils."""

import importlib.resources
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

    def test_add_control_surface_zero_deflection_returns_self(self):
        """Verify that zero deflection returns the same Airfoil instance (optimization)."""
        result = self.naca0012_airfoil.add_control_surface(
            deflection=0.0,
            hinge_point=0.75,
        )
        self.assertIs(result, self.naca0012_airfoil)

    def test_parameter_validation_invalid_inputs(self):
        """Test that invalid parameters raise appropriate errors."""
        # Test outlines with insufficient points in upper portion. These outlines have
        # only 4 points total, and the point ordering places only 1 point in the upper
        # portion (the leading point itself), which fails the requirement for at least
        # 3 unique points in each portion.
        invalid_outline = np.array([[-0.1, 0.0], [0.5, 0.1], [1.0, 0.0], [0.5, -0.1]])
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(
                name="Invalid", outline_A_lp=invalid_outline, resample=False
            )

        invalid_outline2 = np.array([[0.0, 0.0], [0.5, 0.1], [1.1, 0.0], [0.5, -0.1]])
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(
                name="Invalid2", outline_A_lp=invalid_outline2, resample=False
            )

    def test_naca_4_series_validation(self):
        """Test validation of NACA 4 series airfoil parameters."""
        # Test that NACA 4 series airfoils with thickness above 30% raise a ValueError.
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA0031")

        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA2499")

        # Test that NACA0000 (zero thickness) raises a ValueError.
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA0000")

        # Test that NACA 4 series at exactly 30% thickness is valid.
        airfoil_30_percent = ps.geometry.airfoil.Airfoil(name="NACA0030")
        self.assertIsNotNone(airfoil_30_percent.outline_A_lp)

        # Test that airfoils with inconsistent camber parameters (first two digits must
        # both be zero or both be non zero) raise a ValueError.
        # Case 1: Camber without camber location (e.g., NACA1012, NACA9012).
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA1012")

        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA9012")

        # Case 2: Camber location without camber (e.g., NACA0112, NACA0912).
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA0112")

        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA0912")

        # Test that symmetric airfoils (both first digits are 0) are valid.
        airfoil_symmetric = ps.geometry.airfoil.Airfoil(name="NACA0012")
        self.assertIsNotNone(airfoil_symmetric.outline_A_lp)

        # Test that cambered airfoils (both first digits are non zero) are valid.
        airfoil_cambered = ps.geometry.airfoil.Airfoil(name="NACA2412")
        self.assertIsNotNone(airfoil_cambered.outline_A_lp)

        # Test that airfoils with position of maximum camber too close to the leading
        # edge relative to camber and thickness raise a ValueError. The constraint is:
        # camber_loc >= max_camber + thickness / 2.
        # NACA9130: camber_loc=0.1, max_camber=0.09, thickness=0.30
        # 0.1 < 0.09 + 0.15 = 0.24, so this should fail.
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA9130")

        # NACA5115: camber_loc=0.1, max_camber=0.05, thickness=0.15
        # 0.1 < 0.05 + 0.075 = 0.125, so this should fail.
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="NACA5115")

        # NACA4212: camber_loc=0.2, max_camber=0.04, thickness=0.12
        # 0.2 >= 0.04 + 0.06 = 0.10, so this should pass.
        airfoil_valid_camber_pos = ps.geometry.airfoil.Airfoil(name="NACA4212")
        self.assertIsNotNone(airfoil_valid_camber_pos.outline_A_lp)

        # NACA2110: camber_loc=0.1, max_camber=0.02, thickness=0.10
        # 0.1 >= 0.02 + 0.05 = 0.07, so this should pass.
        airfoil_boundary_camber_pos = ps.geometry.airfoil.Airfoil(name="NACA2110")
        self.assertIsNotNone(airfoil_boundary_camber_pos.outline_A_lp)

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

    def test_excessive_rotation_rejection(self):
        """Test that outlines with excessive rotation (>15 degrees) are rejected."""
        # Create a valid airfoil outline
        valid_outline = np.array(
            [
                [1.00, 0.00],
                [0.75, 0.06],
                [0.50, 0.08],
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, -0.08],
                [0.75, -0.06],
                [1.00, 0.00],
            ]
        )

        # Rotate the outline by 20 degrees (exceeds 15 degree limit)
        angle_rad = np.radians(20.0)
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        rotated_outline = (rotation_matrix @ valid_outline.T).T

        with self.assertRaises(ValueError) as context:
            ps.geometry.airfoil.Airfoil(
                name="Excessively Rotated",
                outline_A_lp=rotated_outline,
                resample=False,
            )
        self.assertIn("excessive rotation", str(context.exception).lower())

    def test_minor_rotation_correction(self):
        """Test that minor rotation offsets (<15 degrees) are auto-corrected."""
        # Create a valid airfoil outline
        valid_outline = np.array(
            [
                [1.00, 0.00],
                [0.75, 0.06],
                [0.50, 0.08],
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, -0.08],
                [0.75, -0.06],
                [1.00, 0.00],
            ]
        )

        # Rotate the outline by 10 degrees (within 15 degree limit)
        angle_rad = np.radians(10.0)
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
        rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        rotated_outline = (rotation_matrix @ valid_outline.T).T

        # Should succeed and normalize the rotation
        airfoil = ps.geometry.airfoil.Airfoil(
            name="Minor Rotation",
            outline_A_lp=rotated_outline,
            resample=False,
        )

        # After normalization, leading point should be at ~[0, 0]
        lp_index = int(np.argmin(airfoil.outline_A_lp[:, 0]))
        lp = airfoil.outline_A_lp[lp_index, :]
        npt.assert_allclose(lp, [0.0, 0.0], atol=1e-6)

        # Trailing edge should be at ~x=1.0
        self.assertAlmostEqual(airfoil.outline_A_lp[0, 0], 1.0, places=5)
        self.assertAlmostEqual(airfoil.outline_A_lp[-1, 0], 1.0, places=5)

    def test_auto_normalization_translation_and_scale(self):
        """Test that outlines are auto-normalized for position and scale."""
        # Create an outline that is translated and scaled (not at origin, chord != 1)
        # This is a simple diamond-like airfoil shape scaled by 2x and translated
        base_outline = np.array(
            [
                [1.00, 0.00],
                [0.75, 0.06],
                [0.50, 0.08],
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, -0.08],
                [0.75, -0.06],
                [1.00, 0.00],
            ]
        )

        # Scale by 2x and translate by (5, 3)
        scaled_translated = base_outline * 2.0 + np.array([5.0, 3.0])

        airfoil = ps.geometry.airfoil.Airfoil(
            name="Scaled and Translated",
            outline_A_lp=scaled_translated,
            resample=False,
        )

        # After normalization, leading point should be at ~[0, 0]
        lp_index = int(np.argmin(airfoil.outline_A_lp[:, 0]))
        lp = airfoil.outline_A_lp[lp_index, :]
        npt.assert_allclose(lp, [0.0, 0.0], atol=1e-6)

        # Trailing edge should be at ~x=1.0 (unit chord)
        self.assertAlmostEqual(airfoil.outline_A_lp[0, 0], 1.0, places=5)
        self.assertAlmostEqual(airfoil.outline_A_lp[-1, 0], 1.0, places=5)

        # x values should be in [0, 1]
        self.assertTrue(np.all(airfoil.outline_A_lp[:, 0] >= 0.0 - 1e-9))
        self.assertTrue(np.all(airfoil.outline_A_lp[:, 0] <= 1.0 + 1e-9))

    def test_self_intersection_rejection(self):
        """Test that self-intersecting outlines are rejected."""
        # Create an outline where upper surface crosses below lower surface
        # This is a "twisted" airfoil that would self-intersect
        self_intersecting_outline = np.array(
            [
                [1.00, 0.00],
                [0.75, 0.06],
                [0.50, -0.05],  # Upper surface goes below y=0
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, 0.05],  # Lower surface goes above y=0
                [0.75, -0.06],
                [1.00, 0.00],
            ]
        )

        with self.assertRaises(ValueError) as context:
            ps.geometry.airfoil.Airfoil(
                name="Self Intersecting",
                outline_A_lp=self_intersecting_outline,
                resample=False,
            )
        # The error should mention self-intersection or upper/lower y values
        error_msg = str(context.exception).lower()
        self.assertTrue(
            "upper" in error_msg or "lower" in error_msg or "intersect" in error_msg
        )

    def test_open_trailing_edge_support(self):
        """Test that open (blunt) trailing edges are supported."""
        # Create an outline with an open trailing edge (upper TE y != lower TE y)
        open_te_outline = np.array(
            [
                [1.00, 0.01],  # Upper TE slightly above centerline
                [0.75, 0.06],
                [0.50, 0.08],
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, -0.08],
                [0.75, -0.06],
                [1.00, -0.01],  # Lower TE slightly below centerline
            ]
        )

        # Should succeed - open trailing edges are allowed
        airfoil = ps.geometry.airfoil.Airfoil(
            name="Open TE",
            outline_A_lp=open_te_outline,
            resample=False,
        )

        # Verify airfoil was created successfully
        self.assertIsNotNone(airfoil.outline_A_lp)
        self.assertIsNotNone(airfoil.mcl_A_lp)

        # Upper TE y should be >= lower TE y
        upper_te_y = airfoil.outline_A_lp[0, 1]
        lower_te_y = airfoil.outline_A_lp[-1, 1]
        self.assertGreaterEqual(upper_te_y, lower_te_y)

    def test_upper_x_non_increasing_validation(self):
        """Test that upper outline must have non-increasing x values."""
        # Create an outline where upper portion has increasing x (invalid)
        invalid_upper_outline = np.array(
            [
                [1.00, 0.00],
                [0.50, 0.06],  # x decreases (OK)
                [0.60, 0.08],  # x increases (invalid for upper portion)
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, -0.08],
                [0.75, -0.06],
                [1.00, 0.00],
            ]
        )

        with self.assertRaises(ValueError) as context:
            ps.geometry.airfoil.Airfoil(
                name="Invalid Upper",
                outline_A_lp=invalid_upper_outline,
                resample=False,
            )
        self.assertIn("upper", str(context.exception).lower())

    def test_lower_x_non_decreasing_validation(self):
        """Test that lower outline must have non-decreasing x values."""
        # Create an outline where lower portion has decreasing x (invalid)
        invalid_lower_outline = np.array(
            [
                [1.00, 0.00],
                [0.75, 0.06],
                [0.50, 0.08],
                [0.25, 0.06],
                [0.00, 0.00],
                [0.25, -0.06],
                [0.50, -0.08],
                [0.40, -0.06],  # x decreases (invalid for lower portion)
                [1.00, 0.00],
            ]
        )

        with self.assertRaises(ValueError) as context:
            ps.geometry.airfoil.Airfoil(
                name="Invalid Lower",
                outline_A_lp=invalid_lower_outline,
                resample=False,
            )
        self.assertIn("lower", str(context.exception).lower())

    def test_minimum_points_validation(self):
        """Test that outline must have at least 5 points."""
        # Create an outline with only 4 points
        too_few_points = np.array(
            [
                [1.00, 0.00],
                [0.50, 0.05],
                [0.00, 0.00],
                [0.50, -0.05],
            ]
        )

        with self.assertRaises(ValueError) as context:
            ps.geometry.airfoil.Airfoil(
                name="Too Few Points",
                outline_A_lp=too_few_points,
                resample=False,
            )
        # Should fail for either "at least five points" or "at least three unique points"
        error_msg = str(context.exception).lower()
        self.assertTrue(
            "five" in error_msg or "three" in error_msg or "unique" in error_msg
        )

    def test_all_valid_naca4_airfoils_load(self):
        """Test that all valid NACA 4 series airfoils load without errors.

        NACA 4 series constraints:
        1. Cannot be "0000" (zero thickness)
        2. Thickness (last two digits) must be <= 30%
        3. First two digits must either both be zero (symmetric) or both be
           non zero (cambered)
        4. For cambered airfoils: camber_loc >= max_camber + thickness/2
        """
        failed_airfoils = []

        for first_digit in range(10):  # 0-9: max camber (%)
            for second_digit in range(10):  # 0-9: camber location (x10%)
                for thickness in range(1, 31):  # 01-30: thickness (%)
                    # Constraint: First two digits must both be zero or both be non-zero
                    max_camber = first_digit * 0.01
                    camber_loc = second_digit * 0.1
                    if (max_camber > 0) != (camber_loc > 0):
                        continue

                    # Constraint: For cambered airfoils, camber_loc >= max_camber + thickness/2
                    if max_camber > 0:
                        thickness_fraction = thickness * 0.01
                        if camber_loc < max_camber + thickness_fraction / 2:
                            continue

                    # Build the NACA name and try to load
                    name = f"NACA{first_digit}{second_digit}{thickness:02d}"
                    try:
                        airfoil = ps.geometry.airfoil.Airfoil(name=name)
                        # Basic sanity checks
                        self.assertIsNotNone(airfoil.outline_A_lp)
                        self.assertIsNotNone(airfoil.mcl_A_lp)
                    except Exception as e:
                        failed_airfoils.append((name, str(e)))

        # Report all failures at once for easier debugging
        if failed_airfoils:
            failure_msg = "\n".join(
                [f"  {name}: {error}" for name, error in failed_airfoils]
            )
            self.fail(
                f"{len(failed_airfoils)} NACA 4 series airfoils failed:\n{failure_msg}"
            )

    def test_all_database_airfoils_load(self):
        """Test that all airfoils in the database load without errors."""
        # Get all airfoil names from the database
        airfoils_dir = importlib.resources.files("pterasoftware.geometry").joinpath(
            "_airfoils"
        )
        airfoil_names = []
        for item in airfoils_dir.iterdir():
            if item.name.endswith(".dat"):
                airfoil_names.append(item.name[:-4])  # Remove .dat extension

        self.assertGreater(len(airfoil_names), 0, "No airfoils found in database")

        failed_airfoils = []
        for name in sorted(airfoil_names):
            try:
                airfoil = ps.geometry.airfoil.Airfoil(name=name)
                # Basic sanity checks
                self.assertIsNotNone(airfoil.outline_A_lp)
                self.assertIsNotNone(airfoil.mcl_A_lp)
            except Exception as e:
                failed_airfoils.append((name, str(e)))

        # Report all failures at once for easier debugging
        if failed_airfoils:
            failure_msg = "\n".join(
                [f"  {name}: {error}" for name, error in failed_airfoils]
            )
            self.fail(
                f"{len(failed_airfoils)} database airfoils failed:\n{failure_msg}"
            )


if __name__ == "__main__":
    unittest.main()
