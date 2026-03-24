"""This module contains classes to test Airfoils."""

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


class TestAirfoilImmutability(unittest.TestCase):
    """Tests for Airfoil attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()
        self.naca2412_airfoil = geometry_fixtures.make_naca2412_airfoil_fixture()

    def test_immutable_name_property(self):
        """Test that name property is read only."""
        with self.assertRaises(AttributeError):
            self.naca0012_airfoil.name = "new_name"

    def test_immutable_resample_property(self):
        """Test that resample property is read only."""
        with self.assertRaises(AttributeError):
            self.naca0012_airfoil.resample = False

    def test_immutable_n_points_per_side_property(self):
        """Test that n_points_per_side property is read only."""
        with self.assertRaises(AttributeError):
            self.naca0012_airfoil.n_points_per_side = 100

    def test_immutable_outline_A_lp_property(self):
        """Test that outline_A_lp property is read only."""
        new_outline = np.zeros((10, 2))
        with self.assertRaises(AttributeError):
            self.naca0012_airfoil.outline_A_lp = new_outline

    def test_immutable_mcl_A_lp_property(self):
        """Test that mcl_A_lp property is read only."""
        new_mcl = np.zeros((10, 2))
        with self.assertRaises(AttributeError):
            self.naca0012_airfoil.mcl_A_lp = new_mcl

    def test_outline_A_lp_array_read_only(self):
        """Test that outline_A_lp array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.outline_A_lp[0, 0] = 999.0

    def test_mcl_A_lp_array_read_only(self):
        """Test that mcl_A_lp array cannot be modified in place."""
        self.assertIsNotNone(self.naca0012_airfoil.mcl_A_lp)
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.mcl_A_lp[0, 0] = 999.0


class TestAirfoilDeepCopy(unittest.TestCase):
    """Tests for Airfoil.__deepcopy__ method."""

    def setUp(self):
        """Set up test fixtures for deepcopy tests."""
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()
        self.naca2412_airfoil = geometry_fixtures.make_naca2412_airfoil_fixture()
        self.custom_outline_airfoil = (
            geometry_fixtures.make_custom_outline_airfoil_fixture()
        )
        self.resampled_airfoil = geometry_fixtures.make_resampled_airfoil_fixture()
        self.non_resampled_airfoil = (
            geometry_fixtures.make_non_resampled_airfoil_fixture()
        )

    def test_deepcopy_creates_new_instance(self):
        """Test that deepcopy creates a new Airfoil instance."""
        import copy

        original = self.naca0012_airfoil
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.geometry.airfoil.Airfoil)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_all_attributes(self):
        """Test that deepcopy preserves all attribute values."""
        import copy

        original = self.naca0012_airfoil
        copied = copy.deepcopy(original)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.resample, original.resample)
        self.assertEqual(copied.n_points_per_side, original.n_points_per_side)
        npt.assert_array_equal(copied.outline_A_lp, original.outline_A_lp)
        npt.assert_array_equal(copied.mcl_A_lp, original.mcl_A_lp)

    def test_deepcopy_creates_independent_outline_array(self):
        """Test that deepcopy creates an independent copy of outline_A_lp."""
        import copy

        original = self.naca0012_airfoil
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.outline_A_lp, original.outline_A_lp)
        npt.assert_array_equal(copied.outline_A_lp, original.outline_A_lp)

    def test_deepcopy_creates_independent_mcl_array(self):
        """Test that deepcopy creates an independent copy of mcl_A_lp."""
        import copy

        original = self.naca0012_airfoil
        copied = copy.deepcopy(original)

        self.assertIsNotNone(copied.mcl_A_lp)
        self.assertIsNotNone(original.mcl_A_lp)
        self.assertIsNot(copied.mcl_A_lp, original.mcl_A_lp)
        npt.assert_array_equal(copied.mcl_A_lp, original.mcl_A_lp)

    def test_deepcopy_arrays_remain_read_only(self):
        """Test that deepcopied arrays are still read only."""
        import copy

        original = self.naca0012_airfoil
        copied = copy.deepcopy(original)

        with self.assertRaises(ValueError):
            copied.outline_A_lp[0, 0] = 999.0

        with self.assertRaises(ValueError):
            copied.mcl_A_lp[0, 0] = 999.0

    def test_deepcopy_non_resampled_airfoil(self):
        """Test that deepcopy works correctly for non-resampled Airfoils."""
        import copy

        original = self.non_resampled_airfoil
        copied = copy.deepcopy(original)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.resample, original.resample)
        self.assertFalse(copied.resample)
        npt.assert_array_equal(copied.outline_A_lp, original.outline_A_lp)

    def test_deepcopy_custom_outline_airfoil(self):
        """Test that deepcopy works correctly for custom outline Airfoils."""
        import copy

        original = self.custom_outline_airfoil
        copied = copy.deepcopy(original)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.resample, original.resample)
        self.assertEqual(copied.n_points_per_side, original.n_points_per_side)
        npt.assert_array_equal(copied.outline_A_lp, original.outline_A_lp)
        npt.assert_array_equal(copied.mcl_A_lp, original.mcl_A_lp)

    def test_deepcopy_handles_memo_correctly(self):
        """Test that deepcopy uses memo dictionary correctly to prevent duplicates."""
        import copy

        original = self.naca0012_airfoil
        memo = {}
        copied1 = copy.deepcopy(original, memo)
        copied2 = copy.deepcopy(original, memo)

        # Both copies should reference the same object in memo
        self.assertIs(copied1, copied2)


class TestAirfoilGetPlottableData(unittest.TestCase):
    """Tests for Airfoil.get_plottable_data method."""

    def setUp(self):
        """Set up test fixtures for get_plottable_data tests."""
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()
        self.naca2412_airfoil = geometry_fixtures.make_naca2412_airfoil_fixture()

    def test_get_plottable_data_returns_list_when_show_false(self):
        """Test that get_plottable_data returns a list when show is False."""
        result = self.naca0012_airfoil.get_plottable_data(show=False)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_returns_outline_and_mcl(self):
        """Test that get_plottable_data returns outline and MCL data."""
        result = self.naca0012_airfoil.get_plottable_data(show=False)

        outline_data = result[0]
        mcl_data = result[1]

        self.assertIsInstance(outline_data, np.ndarray)
        self.assertIsInstance(mcl_data, np.ndarray)

        # Verify shape is (N, 2)
        self.assertEqual(len(outline_data.shape), 2)
        self.assertEqual(outline_data.shape[1], 2)
        self.assertEqual(len(mcl_data.shape), 2)
        self.assertEqual(mcl_data.shape[1], 2)

    def test_get_plottable_data_outline_matches_attribute(self):
        """Test that returned outline matches outline_A_lp attribute."""
        result = self.naca0012_airfoil.get_plottable_data(show=False)

        npt.assert_array_equal(result[0], self.naca0012_airfoil.outline_A_lp)

    def test_get_plottable_data_mcl_matches_attribute(self):
        """Test that returned MCL matches mcl_A_lp attribute."""
        result = self.naca0012_airfoil.get_plottable_data(show=False)

        npt.assert_array_equal(result[1], self.naca0012_airfoil.mcl_A_lp)

    def test_get_plottable_data_default_show_is_false(self):
        """Test that the default value for show parameter is False."""
        # Call without any parameters
        result = self.naca0012_airfoil.get_plottable_data()

        # Should return a list (not None) when show defaults to False
        self.assertIsInstance(result, list)

    def test_get_plottable_data_accepts_numpy_bool(self):
        """Test that get_plottable_data accepts numpy bool for show parameter."""
        result = self.naca0012_airfoil.get_plottable_data(show=np.bool_(False))

        self.assertIsInstance(result, list)

    def test_get_plottable_data_cambered_airfoil(self):
        """Test get_plottable_data with a cambered airfoil."""
        result = self.naca2412_airfoil.get_plottable_data(show=False)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

        _ = result[0]
        mcl_data = result[1]

        # For cambered airfoil, MCL should not be at y=0 everywhere
        mcl_y_values = mcl_data[:, 1]
        # MCL for cambered airfoil has non-zero y values (camber)
        self.assertGreater(np.max(np.abs(mcl_y_values)), 0.0)


class TestAirfoilAddControlSurface(unittest.TestCase):
    """Additional tests for Airfoil.add_control_surface method."""

    def setUp(self):
        """Set up test fixtures for add_control_surface tests."""
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()
        self.naca2412_airfoil = geometry_fixtures.make_naca2412_airfoil_fixture()

    def test_add_control_surface_positive_deflection(self):
        """Test add_control_surface with positive (downward) deflection."""
        modified = self.naca0012_airfoil.add_control_surface(
            deflection=5.0, hinge_point=0.75
        )

        self.assertIsInstance(modified, ps.geometry.airfoil.Airfoil)
        self.assertIn("flapped", modified.name)

    def test_add_control_surface_negative_deflection(self):
        """Test add_control_surface with negative (upward) deflection."""
        modified = self.naca0012_airfoil.add_control_surface(
            deflection=-5.0, hinge_point=0.75
        )

        self.assertIsInstance(modified, ps.geometry.airfoil.Airfoil)
        self.assertIn("flapped", modified.name)

    def test_add_control_surface_various_hinge_points(self):
        """Test add_control_surface with various hinge point locations."""
        hinge_points = [0.5, 0.6, 0.7, 0.8, 0.9]
        for hinge_point in hinge_points:
            with self.subTest(hinge_point=hinge_point):
                modified = self.naca0012_airfoil.add_control_surface(
                    deflection=3.0, hinge_point=hinge_point
                )
                self.assertIsInstance(modified, ps.geometry.airfoil.Airfoil)

    def test_add_control_surface_invalid_deflection(self):
        """Test that invalid deflection values raise errors."""
        # Deflection too large
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.add_control_surface(deflection=6.0, hinge_point=0.75)

        with self.assertRaises(ValueError):
            self.naca0012_airfoil.add_control_surface(deflection=-6.0, hinge_point=0.75)

    def test_add_control_surface_invalid_hinge_point(self):
        """Test that invalid hinge point values raise errors."""
        # Hinge point at edge (exclusive boundaries)
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.add_control_surface(deflection=3.0, hinge_point=0.0)

        with self.assertRaises(ValueError):
            self.naca0012_airfoil.add_control_surface(deflection=3.0, hinge_point=1.0)

        # Hinge point outside valid range
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.add_control_surface(deflection=3.0, hinge_point=-0.1)

        with self.assertRaises(ValueError):
            self.naca0012_airfoil.add_control_surface(deflection=3.0, hinge_point=1.1)

    def test_add_control_surface_accepts_int_parameters(self):
        """Test that add_control_surface accepts integer parameters."""
        modified = self.naca0012_airfoil.add_control_surface(
            deflection=3, hinge_point=0.75  # int for deflection
        )
        self.assertIsInstance(modified, ps.geometry.airfoil.Airfoil)

    def test_add_control_surface_preserves_name_base(self):
        """Test that modified Airfoil name includes original name."""
        original_name = self.naca0012_airfoil.name
        modified = self.naca0012_airfoil.add_control_surface(
            deflection=3.0, hinge_point=0.75
        )
        self.assertIn(original_name, modified.name)

    def test_add_control_surface_cambered_airfoil(self):
        """Test add_control_surface on a cambered airfoil."""
        modified = self.naca2412_airfoil.add_control_surface(
            deflection=4.0, hinge_point=0.8
        )
        self.assertIsInstance(modified, ps.geometry.airfoil.Airfoil)
        self.assertIsNotNone(modified.outline_A_lp)
        self.assertIsNotNone(modified.mcl_A_lp)


class TestAirfoilGetResampledMcl(unittest.TestCase):
    """Additional tests for Airfoil.get_resampled_mcl method."""

    def setUp(self):
        """Set up test fixtures for get_resampled_mcl tests."""
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()
        self.naca2412_airfoil = geometry_fixtures.make_naca2412_airfoil_fixture()

    def test_get_resampled_mcl_returns_correct_shape(self):
        """Test that get_resampled_mcl returns correct shape."""
        mcl_fractions = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
        result = self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

        self.assertEqual(result.shape, (5, 2))

    def test_get_resampled_mcl_endpoints(self):
        """Test that get_resampled_mcl endpoints are at x=0 and x=1."""
        mcl_fractions = np.array([0.0, 0.5, 1.0])
        result = self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

        # First point should be at x ~ 0.0
        self.assertAlmostEqual(result[0, 0], 0.0, places=4)

        # Last point should be at x ~ 1.0
        self.assertAlmostEqual(result[-1, 0], 1.0, places=4)

    def test_get_resampled_mcl_monotonic_x(self):
        """Test that resampled MCL has monotonically increasing x values."""
        mcl_fractions = np.linspace(0, 1, 20)
        result = self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

        x_values = result[:, 0]
        x_diffs = np.diff(x_values)

        # All differences should be positive (monotonically increasing)
        self.assertTrue(np.all(x_diffs >= 0))

    def test_get_resampled_mcl_invalid_first_value(self):
        """Test that invalid first value in mcl_fractions raises error."""
        mcl_fractions = np.array([0.1, 0.5, 1.0])  # First value is not 0.0
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

    def test_get_resampled_mcl_invalid_last_value(self):
        """Test that invalid last value in mcl_fractions raises error."""
        mcl_fractions = np.array([0.0, 0.5, 0.9])  # Last value is not 1.0
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

    def test_get_resampled_mcl_values_outside_range(self):
        """Test that values outside [0, 1] range raise error."""
        mcl_fractions_low = np.array([-0.1, 0.5, 1.0])
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions_low)

        mcl_fractions_high = np.array([0.0, 0.5, 1.1])
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions_high)

    def test_get_resampled_mcl_non_ascending_values(self):
        """Test that non-ascending values raise error."""
        mcl_fractions = np.array([0.0, 0.7, 0.5, 1.0])  # 0.7 > 0.5
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

    def test_get_resampled_mcl_duplicate_values(self):
        """Test that duplicate values raise error."""
        mcl_fractions = np.array([0.0, 0.5, 0.5, 1.0])  # Duplicate 0.5
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

    def test_get_resampled_mcl_too_few_values(self):
        """Test that too few values raise error."""
        mcl_fractions = np.array([0.5])  # Only one value
        with self.assertRaises(ValueError):
            self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

    def test_get_resampled_mcl_accepts_list(self):
        """Test that get_resampled_mcl accepts list input."""
        mcl_fractions = [0.0, 0.25, 0.5, 0.75, 1.0]
        result = self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

        self.assertEqual(result.shape, (5, 2))

    def test_get_resampled_mcl_accepts_tuple(self):
        """Test that get_resampled_mcl accepts tuple input."""
        mcl_fractions = (0.0, 0.5, 1.0)
        result = self.naca0012_airfoil.get_resampled_mcl(mcl_fractions=mcl_fractions)

        self.assertEqual(result.shape, (3, 2))


class TestAirfoilEdgeCases(unittest.TestCase):
    """Tests for edge cases and special scenarios in Airfoil."""

    def setUp(self):
        """Set up test fixtures for edge case tests."""
        self.naca0012_airfoil = geometry_fixtures.make_naca0012_airfoil_fixture()

    def test_minimum_n_points_per_side(self):
        """Test Airfoil with minimum valid n_points_per_side (3)."""
        min_points_airfoil = (
            geometry_fixtures.make_minimum_n_points_per_side_airfoil_fixture()
        )

        self.assertEqual(min_points_airfoil.n_points_per_side, 3)
        self.assertIsNotNone(min_points_airfoil.outline_A_lp)
        self.assertIsNotNone(min_points_airfoil.mcl_A_lp)

    def test_invalid_n_points_per_side_too_low(self):
        """Test that n_points_per_side below 3 raises error."""
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="naca0012", n_points_per_side=2)

        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="naca0012", n_points_per_side=1)

    def test_thick_naca_airfoil(self):
        """Test a thick NACA airfoil (30% thickness)."""
        thick_airfoil = geometry_fixtures.make_thick_naca_airfoil_fixture()

        self.assertIsNotNone(thick_airfoil.outline_A_lp)
        self.assertIsNotNone(thick_airfoil.mcl_A_lp)

        # Check thickness is approximately 30%
        outline = thick_airfoil.outline_A_lp
        y_coords = outline[:, 1]
        max_thickness = np.max(y_coords) - np.min(y_coords)
        self.assertAlmostEqual(max_thickness, 0.30, places=2)

    def test_blunt_trailing_edge_airfoil(self):
        """Test an Airfoil with blunt trailing edge."""
        blunt_te_airfoil = geometry_fixtures.make_blunt_trailing_edge_airfoil_fixture()

        self.assertIsNotNone(blunt_te_airfoil.outline_A_lp)
        self.assertIsNotNone(blunt_te_airfoil.mcl_A_lp)

        # Upper TE y should be greater than lower TE y
        upper_te_y = blunt_te_airfoil.outline_A_lp[0, 1]
        lower_te_y = blunt_te_airfoil.outline_A_lp[-1, 1]
        self.assertGreater(upper_te_y, lower_te_y)

    def test_case_insensitive_naca_name(self):
        """Test that NACA airfoil names are case insensitive."""
        mixed_case_airfoil = (
            geometry_fixtures.make_case_insensitive_naca_airfoil_fixture()
        )

        self.assertIsNotNone(mixed_case_airfoil.outline_A_lp)
        self.assertIsNotNone(mixed_case_airfoil.mcl_A_lp)

    def test_whitespace_in_name(self):
        """Test that leading and trailing whitespace in name is handled."""
        airfoil = ps.geometry.airfoil.Airfoil(name="  naca0012  ")

        self.assertIsNotNone(airfoil.outline_A_lp)
        self.assertIsNotNone(airfoil.mcl_A_lp)

    def test_invalid_airfoil_name(self):
        """Test that invalid airfoil name raises error."""
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(name="not_a_valid_airfoil_name")

    def test_invalid_name_type(self):
        """Test that invalid name type raises error."""
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airfoil.Airfoil(name=12345)

    def test_invalid_resample_type(self):
        """Test that invalid resample type raises error."""
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airfoil.Airfoil(name="naca0012", resample="yes")

    def test_invalid_n_points_per_side_type(self):
        """Test that invalid n_points_per_side type raises error."""
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airfoil.Airfoil(name="naca0012", n_points_per_side=100.5)

    def test_outline_with_different_array_types(self):
        """Test that outline_A_lp accepts different array like types."""
        # Test with list of lists
        outline_list = [
            [1.00, 0.00],
            [0.75, 0.05],
            [0.50, 0.06],
            [0.25, 0.04],
            [0.00, 0.00],
            [0.25, -0.04],
            [0.50, -0.06],
            [0.75, -0.05],
            [1.00, 0.00],
        ]
        airfoil_list = ps.geometry.airfoil.Airfoil(
            name="list_outline", outline_A_lp=outline_list, resample=False
        )
        self.assertIsNotNone(airfoil_list.outline_A_lp)

        # Test with tuple of tuples
        outline_tuple = tuple(tuple(row) for row in outline_list)
        airfoil_tuple = ps.geometry.airfoil.Airfoil(
            name="tuple_outline", outline_A_lp=outline_tuple, resample=False
        )
        self.assertIsNotNone(airfoil_tuple.outline_A_lp)

    def test_outline_with_int_coordinates(self):
        """Test that outline_A_lp accepts integer coordinates."""
        outline_int = np.array(
            [
                [1, 0],
                [0, 0],  # Minimal points
                [1, 0],
            ],
            dtype=int,
        )
        # This will fail validation due to too few points, but type should be accepted
        with self.assertRaises(ValueError):
            ps.geometry.airfoil.Airfoil(
                name="int_outline", outline_A_lp=outline_int, resample=False
            )


if __name__ == "__main__":
    unittest.main()
