"""This module contains classes to test Panels."""

import copy
import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _panel
from tests.unit.fixtures import panel_fixtures


class TestPanel(unittest.TestCase):
    """This class contains unit tests for the Panel class."""

    def setUp(self):
        """Set up test fixtures for Panel tests."""
        self.basic_panel = panel_fixtures.make_basic_panel_fixture()

    def test_initialization_valid_parameters(self):
        """Test Panel initialization with valid parameters."""
        panel = self.basic_panel

        # Test that Panel initializes correctly
        self.assertIsInstance(panel, _panel.Panel)

        # Test that corner points are correctly stored
        npt.assert_array_equal(panel.Frpp_G_Cg, np.array([0.0, 0.5, 0.0]))
        npt.assert_array_equal(panel.Flpp_G_Cg, np.array([0.0, 0.0, 0.0]))
        npt.assert_array_equal(panel.Blpp_G_Cg, np.array([1.0, 0.0, 0.0]))
        npt.assert_array_equal(panel.Brpp_G_Cg, np.array([1.0, 0.5, 0.0]))

        # Test that edge bools are correctly stored
        self.assertFalse(panel.is_leading_edge)
        self.assertFalse(panel.is_trailing_edge)

    def test_initial_attribute_values(self):
        """Test that optional attributes start as None."""
        panel = self.basic_panel

        # Test that position attributes are None initially
        self.assertIsNone(panel.is_right_edge)
        self.assertIsNone(panel.is_left_edge)
        self.assertIsNone(panel.local_chordwise_position)
        self.assertIsNone(panel.local_spanwise_position)

        # Test that force and moment attributes are None initially
        self.assertIsNone(panel.forces_GP1)
        self.assertIsNone(panel.moments_GP1_CgP1)
        self.assertIsNone(panel.forces_W)
        self.assertIsNone(panel.moments_W_CgP1)

    def test_rightLeg_G_property(self):
        """Test right leg vector calculation."""
        panel = self.basic_panel

        rightLeg_G = panel.rightLeg_G

        # Right leg should go from back-right to front-right
        expected_rightLeg_G = np.array([-1.0, 0.0, 0.0])
        npt.assert_array_almost_equal(rightLeg_G, expected_rightLeg_G)

    def test_frontLeg_G_property(self):
        """Test front leg vector calculation."""
        panel = self.basic_panel

        frontLeg_G = panel.frontLeg_G

        # Front leg should go from front-right to front-left
        expected_frontLeg_G = np.array([0.0, -0.5, 0.0])
        npt.assert_array_almost_equal(frontLeg_G, expected_frontLeg_G)

    def test_leftLeg_G_property(self):
        """Test left leg vector calculation."""
        panel = self.basic_panel

        leftLeg_G = panel.leftLeg_G

        # Left leg should go from front-left to back-left
        expected_leftLeg_G = np.array([1.0, 0.0, 0.0])
        npt.assert_array_almost_equal(leftLeg_G, expected_leftLeg_G)

    def test_backLeg_G_property(self):
        """Test back leg vector calculation."""
        panel = self.basic_panel

        backLeg_G = panel.backLeg_G

        # Back leg should go from back-left to back-right
        expected_backLeg_G = np.array([0.0, 0.5, 0.0])
        npt.assert_array_almost_equal(backLeg_G, expected_backLeg_G)

    def test_Frbvp_G_Cg_property(self):
        """Test front-right bound vortex point at 75% chord."""
        panel = self.basic_panel

        Frbvp_G_Cg = panel.Frbvp_G_Cg

        # Should be at back-right plus 75% of right leg (towards front)
        # Brpp_G_Cg = [1.0, 0.5, 0.0]
        # rightLeg_G = [-1.0, 0.0, 0.0]
        # Expected: [1.0, 0.5, 0.0] + 0.75 * [-1.0, 0.0, 0.0] = [0.25, 0.5, 0.0]
        expected_Frbvp_G_Cg = np.array([0.25, 0.5, 0.0])
        npt.assert_array_almost_equal(Frbvp_G_Cg, expected_Frbvp_G_Cg)

    def test_Flbvp_G_Cg_property(self):
        """Test front-left bound vortex point at 25% chord."""
        panel = self.basic_panel

        Flbvp_G_Cg = panel.Flbvp_G_Cg

        # Should be at front-left plus 25% of left leg (towards back)
        # Flpp_G_Cg = [0.0, 0.0, 0.0]
        # leftLeg_G = [1.0, 0.0, 0.0]
        # Expected: [0.0, 0.0, 0.0] + 0.25 * [1.0, 0.0, 0.0] = [0.25, 0.0, 0.0]
        expected_Flbvp_G_Cg = np.array([0.25, 0.0, 0.0])
        npt.assert_array_almost_equal(Flbvp_G_Cg, expected_Flbvp_G_Cg)

    def test_Cpp_G_Cg_property(self):
        """Test collocation point at 75% chord midspan."""
        panel = self.basic_panel

        Cpp_G_Cg = panel.Cpp_G_Cg

        # Should be at 75% chord, midspan
        # Expected: [0.75, 0.25, 0.0]
        expected_Cpp_G_Cg = np.array([0.75, 0.25, 0.0])
        npt.assert_array_almost_equal(Cpp_G_Cg, expected_Cpp_G_Cg)

    def test_area_property(self):
        """Test area calculation for rectangular panel."""
        panel = self.basic_panel

        area = panel.area

        # For a rectangular panel with chord 1.0 m and span 0.5 m
        expected_area = 0.5  # square meters
        self.assertAlmostEqual(area, expected_area, places=10)

    def test_unitNormal_G_property(self):
        """Test unit normal vector calculation."""
        panel = self.basic_panel

        unitNormal_G = panel.unitNormal_G

        # For a flat panel in the xy-plane, normal should point in +z direction
        expected_unitNormal_G = np.array([0.0, 0.0, 1.0])
        npt.assert_array_almost_equal(unitNormal_G, expected_unitNormal_G)

        # Verify it's a unit vector
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)

    def test_calculate_projected_area_aligned(self):
        """Test projected area when normal is aligned with panel normal."""
        panel = self.basic_panel

        # Project onto xy-plane (panel is in xy-plane)
        normal_G = np.array([0.0, 0.0, 1.0])
        projected_area = panel.calculate_projected_area(normal_G)

        # Should equal full area
        self.assertAlmostEqual(projected_area, panel.area, places=10)

    def test_calculate_projected_area_perpendicular(self):
        """Test projected area when normal is perpendicular to panel."""
        panel = self.basic_panel

        # Project onto xz-plane (perpendicular to panel in xy-plane)
        normal_G = np.array([0.0, 1.0, 0.0])
        projected_area = panel.calculate_projected_area(normal_G)

        # Should be approximately zero
        self.assertAlmostEqual(projected_area, 0.0, places=10)

    def test_calculate_projected_area_45_degrees(self):
        """Test projected area at a 45 degree angle."""
        panel = self.basic_panel

        # Project at 45 degrees
        normal_G = np.array([0.0, 1.0, 1.0])  # Will be normalized internally
        projected_area = panel.calculate_projected_area(normal_G)

        # Should be area * cos(45 deg) = area / sqrt(2)
        expected_projected_area = panel.area / np.sqrt(2)
        self.assertAlmostEqual(projected_area, expected_projected_area, places=10)

    def test_calculate_projected_area_normalizes_input(self):
        """Test that calculate_projected_area normalizes non-unit normal vectors."""
        panel = self.basic_panel

        # Use a non-unit vector that points in z-direction
        normal_G = np.array([0.0, 0.0, 5.0])
        projected_area = panel.calculate_projected_area(normal_G)

        # Should still equal full area (since direction is same as panel normal)
        self.assertAlmostEqual(projected_area, panel.area, places=10)

    def test_nearly_planar_panel(self):
        """Test with a nearly planar panel."""
        # Create a panel with very slight twist
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 1.0, 0.001]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([2.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([2.0, 1.0, 0.002]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should still calculate area and normal without issues
        area = panel.area
        self.assertGreater(area, 0.0)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)

    def test_twisted_panel(self):
        """Test with a non-planar (twisted) panel."""
        # Create a twisted panel
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 1.0, 0.5]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([2.0, 0.0, -0.5]),
            Brpp_G_Cg=np.array([2.0, 1.0, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should still calculate area and normal without issues
        area = panel.area
        self.assertGreater(area, 0.0)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)

    def test_small_panel(self):
        """Test with very small panel dimensions."""
        # Create a very small panel (0.01 m x 0.01 m)
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.00, 0.01, 0.0]),
            Flpp_G_Cg=np.array([0.00, 0.00, 0.0]),
            Blpp_G_Cg=np.array([0.01, 0.00, 0.0]),
            Brpp_G_Cg=np.array([0.01, 0.01, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should still calculate properties correctly
        area = panel.area
        expected_area = 0.0001
        self.assertAlmostEqual(area, expected_area, places=10)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)

    def test_aspect_ratio_property(self):
        """Test aspect ratio calculation for rectangular panel."""
        panel = self.basic_panel

        aspect_ratio = panel.aspect_ratio

        # For basic panel: span = 0.5 m, chord = 1.0 m
        # Aspect ratio = span / chord = 0.5
        expected_aspect_ratio = 0.5
        self.assertAlmostEqual(aspect_ratio, expected_aspect_ratio, places=10)

    def test_aspect_ratio_square_panel(self):
        """Test aspect ratio for a square panel."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 1.0, 0.0]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1.0, 1.0, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        aspect_ratio = panel.aspect_ratio

        # For a square panel, aspect ratio should be 1.0
        expected_aspect_ratio = 1.0
        self.assertAlmostEqual(aspect_ratio, expected_aspect_ratio, places=10)

    def test_aspect_ratio_high_aspect_panel(self):
        """Test aspect ratio for a high aspect ratio panel."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 4.0, 0.0]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1.0, 4.0, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        aspect_ratio = panel.aspect_ratio

        # For a panel with span = 4.0 m, chord = 1.0 m
        # Aspect ratio = 4.0
        expected_aspect_ratio = 4.0
        self.assertAlmostEqual(aspect_ratio, expected_aspect_ratio, places=10)


class TestPanelGP1Properties(unittest.TestCase):
    """Tests for Panel properties in global (GP1) coordinates."""

    def setUp(self):
        """Set up test fixtures for GP1 property tests."""
        self.basic_panel = panel_fixtures.make_panel_with_gp1_positions_fixture()
        self.offset = np.array([10.0, 20.0, 5.0])

    def test_rightLeg_GP1_property(self):
        """Test right leg vector calculation in GP1 coordinates."""
        panel = self.basic_panel

        rightLeg_GP1 = panel.rightLeg_GP1

        # Right leg should be same as local since it's a vector (translation invariant)
        expected_rightLeg_GP1 = np.array([-1.0, 0.0, 0.0])
        npt.assert_array_almost_equal(rightLeg_GP1, expected_rightLeg_GP1)

    def test_rightLeg_GP1_returns_none_when_not_set(self):
        """Test that rightLeg_GP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        # GP1 positions are None by default
        self.assertIsNone(panel.rightLeg_GP1)

    def test_frontLeg_GP1_property(self):
        """Test front leg vector calculation in GP1 coordinates."""
        panel = self.basic_panel

        frontLeg_GP1 = panel.frontLeg_GP1

        # Front leg should be same as local since it's a vector
        expected_frontLeg_GP1 = np.array([0.0, -0.5, 0.0])
        npt.assert_array_almost_equal(frontLeg_GP1, expected_frontLeg_GP1)

    def test_frontLeg_GP1_returns_none_when_not_set(self):
        """Test that frontLeg_GP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.frontLeg_GP1)

    def test_leftLeg_GP1_property(self):
        """Test left leg vector calculation in GP1 coordinates."""
        panel = self.basic_panel

        leftLeg_GP1 = panel.leftLeg_GP1

        # Left leg should be same as local since it's a vector
        expected_leftLeg_GP1 = np.array([1.0, 0.0, 0.0])
        npt.assert_array_almost_equal(leftLeg_GP1, expected_leftLeg_GP1)

    def test_leftLeg_GP1_returns_none_when_not_set(self):
        """Test that leftLeg_GP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.leftLeg_GP1)

    def test_backLeg_GP1_property(self):
        """Test back leg vector calculation in GP1 coordinates."""
        panel = self.basic_panel

        backLeg_GP1 = panel.backLeg_GP1

        # Back leg should be same as local since it's a vector
        expected_backLeg_GP1 = np.array([0.0, 0.5, 0.0])
        npt.assert_array_almost_equal(backLeg_GP1, expected_backLeg_GP1)

    def test_backLeg_GP1_returns_none_when_not_set(self):
        """Test that backLeg_GP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.backLeg_GP1)

    def test_Frbvp_GP1_CgP1_property(self):
        """Test front right bound vortex point in GP1 coordinates."""
        panel = self.basic_panel

        Frbvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1

        # Should be local value plus offset
        expected_Frbvp_GP1_CgP1 = np.array([0.25, 0.5, 0.0]) + self.offset
        npt.assert_array_almost_equal(Frbvp_GP1_CgP1, expected_Frbvp_GP1_CgP1)

    def test_Frbvp_GP1_CgP1_returns_none_when_not_set(self):
        """Test that Frbvp_GP1_CgP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.Frbvp_GP1_CgP1)

    def test_Flbvp_GP1_CgP1_property(self):
        """Test front left bound vortex point in GP1 coordinates."""
        panel = self.basic_panel

        Flbvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1

        # Should be local value plus offset
        expected_Flbvp_GP1_CgP1 = np.array([0.25, 0.0, 0.0]) + self.offset
        npt.assert_array_almost_equal(Flbvp_GP1_CgP1, expected_Flbvp_GP1_CgP1)

    def test_Flbvp_GP1_CgP1_returns_none_when_not_set(self):
        """Test that Flbvp_GP1_CgP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.Flbvp_GP1_CgP1)

    def test_Cpp_GP1_CgP1_property(self):
        """Test collocation point in GP1 coordinates."""
        panel = self.basic_panel

        Cpp_GP1_CgP1 = panel.Cpp_GP1_CgP1

        # Should be local value plus offset
        expected_Cpp_GP1_CgP1 = np.array([0.75, 0.25, 0.0]) + self.offset
        npt.assert_array_almost_equal(Cpp_GP1_CgP1, expected_Cpp_GP1_CgP1)

    def test_Cpp_GP1_CgP1_returns_none_when_not_set(self):
        """Test that Cpp_GP1_CgP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.Cpp_GP1_CgP1)

    def test_unitNormal_GP1_property(self):
        """Test unit normal vector in GP1 coordinates."""
        panel = self.basic_panel

        unitNormal_GP1 = panel.unitNormal_GP1

        # Should be same as local since it's a unit vector (rotation invariant for
        # pure translation)
        expected_unitNormal_GP1 = np.array([0.0, 0.0, 1.0])
        npt.assert_array_almost_equal(unitNormal_GP1, expected_unitNormal_GP1)

        # Verify it's a unit vector
        self.assertAlmostEqual(np.linalg.norm(unitNormal_GP1), 1.0, places=10)

    def test_unitNormal_GP1_returns_none_when_not_set(self):
        """Test that unitNormal_GP1 returns None when GP1 positions not set."""
        panel = panel_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.unitNormal_GP1)

    def test_GP1_properties_with_rotated_panel(self):
        """Test GP1 properties when global position includes rotation effect."""
        panel = panel_fixtures.make_basic_panel_fixture()

        # Simulate a -90 degree rotation about z-axis in global coordinates
        # Rotation matrix for -90 degrees about z: [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
        # Original: Frpp = [0, 0.5, 0] -> [0.5, 0, 0]
        # Original: Flpp = [0, 0, 0] -> [0, 0, 0]
        # Original: Blpp = [1, 0, 0] -> [0, -1, 0]
        # Original: Brpp = [1, 0.5, 0] -> [0.5, -1, 0]
        panel.Frpp_GP1_CgP1 = np.array([0.5, 0.0, 0.0])
        panel.Flpp_GP1_CgP1 = np.array([0.0, 0.0, 0.0])
        panel.Blpp_GP1_CgP1 = np.array([0.0, -1.0, 0.0])
        panel.Brpp_GP1_CgP1 = np.array([0.5, -1.0, 0.0])

        # The leg vectors should now be rotated
        # Original rightLeg_G = Frpp - Brpp = [0, 0.5, 0] - [1, 0.5, 0] = [-1, 0, 0]
        # After -90 degree rotation: [0, 1, 0]
        rightLeg_GP1 = panel.rightLeg_GP1
        expected_rightLeg_GP1 = np.array([0.0, 1.0, 0.0])
        npt.assert_array_almost_equal(rightLeg_GP1, expected_rightLeg_GP1)

        # Unit normal should still point in +z direction for flat panel (rotation
        # around z axis preserves the z component of vectors in the xy plane)
        unitNormal_GP1 = panel.unitNormal_GP1
        expected_unitNormal_GP1 = np.array([0.0, 0.0, 1.0])
        npt.assert_array_almost_equal(unitNormal_GP1, expected_unitNormal_GP1)


class TestPanelDegenerateCases(unittest.TestCase):
    """Tests for Panel behavior with degenerate or edge case geometries."""

    def test_very_large_panel(self):
        """Test with very large panel dimensions."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 1000.0, 0.0]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1000.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1000.0, 1000.0, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should calculate properties correctly
        area = panel.area
        expected_area = 1000000.0
        self.assertAlmostEqual(area, expected_area, places=5)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)
        npt.assert_array_almost_equal(unitNormal_G, np.array([0.0, 0.0, 1.0]))

    def test_panel_with_negative_coordinates(self):
        """Test with panel in negative coordinate space."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([-1.0, -0.5, -1.0]),
            Flpp_G_Cg=np.array([-1.0, -1.0, -1.0]),
            Blpp_G_Cg=np.array([-2.0, -1.0, -1.0]),
            Brpp_G_Cg=np.array([-2.0, -0.5, -1.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should calculate properties correctly
        area = panel.area
        expected_area = 0.5
        self.assertAlmostEqual(area, expected_area, places=10)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)

        # Collocation point should be in negative space
        Cpp_G_Cg = panel.Cpp_G_Cg
        self.assertTrue(all(coord < 0 for coord in Cpp_G_Cg))

    def test_panel_spanning_origin(self):
        """Test with panel that spans across the origin."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([-0.5, 0.5, 0.0]),
            Flpp_G_Cg=np.array([-0.5, -0.5, 0.0]),
            Blpp_G_Cg=np.array([0.5, -0.5, 0.0]),
            Brpp_G_Cg=np.array([0.5, 0.5, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should calculate properties correctly
        area = panel.area
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=10)

        # Collocation point should be near origin
        Cpp_G_Cg = panel.Cpp_G_Cg
        npt.assert_array_almost_equal(Cpp_G_Cg, np.array([0.25, 0.0, 0.0]))

    def test_nearly_zero_area_panel(self):
        """Test with a panel that has very small (nearly zero) area."""
        # Create a very thin panel
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 1e-10, 0.0]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1.0, 1e-10, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Area should be very small but computable
        area = panel.area
        self.assertGreater(area, 0.0)
        self.assertLess(area, 1e-8)

    def test_collinear_points_panel(self):
        """Test with a degenerate panel where points are collinear (zero area)."""
        # Create a panel with collinear points (all on x-axis)
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Area should be zero
        area = panel.area
        self.assertAlmostEqual(area, 0.0, places=10)

    def test_coincident_points_panel(self):
        """Test with a degenerate panel where all points are coincident."""
        point = np.array([1.0, 2.0, 3.0])
        panel = _panel.Panel(
            Frpp_G_Cg=point.copy(),
            Flpp_G_Cg=point.copy(),
            Blpp_G_Cg=point.copy(),
            Brpp_G_Cg=point.copy(),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Area should be zero
        area = panel.area
        self.assertAlmostEqual(area, 0.0, places=10)

        # Leg vectors should be zero
        npt.assert_array_almost_equal(panel.rightLeg_G, np.array([0.0, 0.0, 0.0]))
        npt.assert_array_almost_equal(panel.frontLeg_G, np.array([0.0, 0.0, 0.0]))
        npt.assert_array_almost_equal(panel.leftLeg_G, np.array([0.0, 0.0, 0.0]))
        npt.assert_array_almost_equal(panel.backLeg_G, np.array([0.0, 0.0, 0.0]))

    def test_panel_in_yz_plane(self):
        """Test with panel in the yz plane (perpendicular to x-axis)."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 0.0, 0.5]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([0.0, 1.0, 0.0]),
            Brpp_G_Cg=np.array([0.0, 1.0, 0.5]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        area = panel.area
        expected_area = 0.5
        self.assertAlmostEqual(area, expected_area, places=10)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)
        # Normal should point in x direction (or negative x)
        self.assertAlmostEqual(abs(unitNormal_G[0]), 1.0, places=10)
        self.assertAlmostEqual(unitNormal_G[1], 0.0, places=10)
        self.assertAlmostEqual(unitNormal_G[2], 0.0, places=10)

    def test_panel_in_xz_plane(self):
        """Test with panel in the xz plane (perpendicular to y-axis)."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 0.0, 0.5]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1.0, 0.0, 0.5]),
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        area = panel.area
        expected_area = 0.5
        self.assertAlmostEqual(area, expected_area, places=10)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)
        # Normal should point in y direction (or negative y)
        self.assertAlmostEqual(unitNormal_G[0], 0.0, places=10)
        self.assertAlmostEqual(abs(unitNormal_G[1]), 1.0, places=10)
        self.assertAlmostEqual(unitNormal_G[2], 0.0, places=10)

    def test_edge_flags_true(self):
        """Test panel with is_leading_edge and is_trailing_edge set to True."""
        panel = _panel.Panel(
            Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
            Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
            Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
            Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
            is_leading_edge=True,
            is_trailing_edge=True,
        )

        self.assertTrue(panel.is_leading_edge)
        self.assertTrue(panel.is_trailing_edge)

        # Geometric properties should still work
        area = panel.area
        expected_area = 0.5
        self.assertAlmostEqual(area, expected_area, places=10)


class TestPanelImmutability(unittest.TestCase):
    """Tests for Panel attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.basic_panel = panel_fixtures.make_basic_panel_fixture()

    def test_Frpp_G_Cg_is_read_only_property(self):
        """Test that Frpp_G_Cg property cannot be reassigned."""
        with self.assertRaises(AttributeError):
            self.basic_panel.Frpp_G_Cg = np.array([99.0, 99.0, 99.0])

    def test_Frpp_G_Cg_array_is_read_only(self):
        """Test that Frpp_G_Cg array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_panel.Frpp_G_Cg[0] = 999.0

    def test_Flpp_G_Cg_is_read_only_property(self):
        """Test that Flpp_G_Cg property cannot be reassigned."""
        with self.assertRaises(AttributeError):
            self.basic_panel.Flpp_G_Cg = np.array([99.0, 99.0, 99.0])

    def test_Flpp_G_Cg_array_is_read_only(self):
        """Test that Flpp_G_Cg array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_panel.Flpp_G_Cg[0] = 999.0

    def test_Blpp_G_Cg_is_read_only_property(self):
        """Test that Blpp_G_Cg property cannot be reassigned."""
        with self.assertRaises(AttributeError):
            self.basic_panel.Blpp_G_Cg = np.array([99.0, 99.0, 99.0])

    def test_Blpp_G_Cg_array_is_read_only(self):
        """Test that Blpp_G_Cg array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_panel.Blpp_G_Cg[0] = 999.0

    def test_Brpp_G_Cg_is_read_only_property(self):
        """Test that Brpp_G_Cg property cannot be reassigned."""
        with self.assertRaises(AttributeError):
            self.basic_panel.Brpp_G_Cg = np.array([99.0, 99.0, 99.0])

    def test_Brpp_G_Cg_array_is_read_only(self):
        """Test that Brpp_G_Cg array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_panel.Brpp_G_Cg[0] = 999.0

    def test_is_leading_edge_is_read_only_property(self):
        """Test that is_leading_edge property cannot be reassigned."""
        with self.assertRaises(AttributeError):
            self.basic_panel.is_leading_edge = True

    def test_is_trailing_edge_is_read_only_property(self):
        """Test that is_trailing_edge property cannot be reassigned."""
        with self.assertRaises(AttributeError):
            self.basic_panel.is_trailing_edge = True


class TestPanelSetOnceProperties(unittest.TestCase):
    """Tests for Panel set once properties."""

    def setUp(self):
        """Set up test fixtures for set once property tests."""
        self.basic_panel = panel_fixtures.make_basic_panel_fixture()

    def test_Frpp_GP1_CgP1_can_be_set_once(self):
        """Test that Frpp_GP1_CgP1 can be set once."""
        self.assertIsNone(self.basic_panel.Frpp_GP1_CgP1)
        self.basic_panel.Frpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        npt.assert_array_equal(
            self.basic_panel.Frpp_GP1_CgP1, np.array([1.0, 2.0, 3.0])
        )

    def test_Frpp_GP1_CgP1_cannot_be_set_twice(self):
        """Test that Frpp_GP1_CgP1 cannot be set twice."""
        self.basic_panel.Frpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(AttributeError):
            self.basic_panel.Frpp_GP1_CgP1 = np.array([4.0, 5.0, 6.0])

    def test_Frpp_GP1_CgP1_array_is_read_only_after_set(self):
        """Test that Frpp_GP1_CgP1 array is read only after being set."""
        self.basic_panel.Frpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            self.basic_panel.Frpp_GP1_CgP1[0] = 999.0

    def test_Flpp_GP1_CgP1_can_be_set_once(self):
        """Test that Flpp_GP1_CgP1 can be set once."""
        self.assertIsNone(self.basic_panel.Flpp_GP1_CgP1)
        self.basic_panel.Flpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        npt.assert_array_equal(
            self.basic_panel.Flpp_GP1_CgP1, np.array([1.0, 2.0, 3.0])
        )

    def test_Flpp_GP1_CgP1_cannot_be_set_twice(self):
        """Test that Flpp_GP1_CgP1 cannot be set twice."""
        self.basic_panel.Flpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(AttributeError):
            self.basic_panel.Flpp_GP1_CgP1 = np.array([4.0, 5.0, 6.0])

    def test_Blpp_GP1_CgP1_can_be_set_once(self):
        """Test that Blpp_GP1_CgP1 can be set once."""
        self.assertIsNone(self.basic_panel.Blpp_GP1_CgP1)
        self.basic_panel.Blpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        npt.assert_array_equal(
            self.basic_panel.Blpp_GP1_CgP1, np.array([1.0, 2.0, 3.0])
        )

    def test_Blpp_GP1_CgP1_cannot_be_set_twice(self):
        """Test that Blpp_GP1_CgP1 cannot be set twice."""
        self.basic_panel.Blpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(AttributeError):
            self.basic_panel.Blpp_GP1_CgP1 = np.array([4.0, 5.0, 6.0])

    def test_Brpp_GP1_CgP1_can_be_set_once(self):
        """Test that Brpp_GP1_CgP1 can be set once."""
        self.assertIsNone(self.basic_panel.Brpp_GP1_CgP1)
        self.basic_panel.Brpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        npt.assert_array_equal(
            self.basic_panel.Brpp_GP1_CgP1, np.array([1.0, 2.0, 3.0])
        )

    def test_Brpp_GP1_CgP1_cannot_be_set_twice(self):
        """Test that Brpp_GP1_CgP1 cannot be set twice."""
        self.basic_panel.Brpp_GP1_CgP1 = np.array([1.0, 2.0, 3.0])
        with self.assertRaises(AttributeError):
            self.basic_panel.Brpp_GP1_CgP1 = np.array([4.0, 5.0, 6.0])

    def test_is_right_edge_can_be_set_once(self):
        """Test that is_right_edge can be set once."""
        self.assertIsNone(self.basic_panel.is_right_edge)
        self.basic_panel.is_right_edge = True
        self.assertTrue(self.basic_panel.is_right_edge)

    def test_is_right_edge_cannot_be_set_twice(self):
        """Test that is_right_edge cannot be set twice."""
        self.basic_panel.is_right_edge = True
        with self.assertRaises(AttributeError):
            self.basic_panel.is_right_edge = False

    def test_is_left_edge_can_be_set_once(self):
        """Test that is_left_edge can be set once."""
        self.assertIsNone(self.basic_panel.is_left_edge)
        self.basic_panel.is_left_edge = False
        self.assertFalse(self.basic_panel.is_left_edge)

    def test_is_left_edge_cannot_be_set_twice(self):
        """Test that is_left_edge cannot be set twice."""
        self.basic_panel.is_left_edge = False
        with self.assertRaises(AttributeError):
            self.basic_panel.is_left_edge = True

    def test_local_chordwise_position_can_be_set_once(self):
        """Test that local_chordwise_position can be set once."""
        self.assertIsNone(self.basic_panel.local_chordwise_position)
        self.basic_panel.local_chordwise_position = 5
        self.assertEqual(self.basic_panel.local_chordwise_position, 5)

    def test_local_chordwise_position_cannot_be_set_twice(self):
        """Test that local_chordwise_position cannot be set twice."""
        self.basic_panel.local_chordwise_position = 5
        with self.assertRaises(AttributeError):
            self.basic_panel.local_chordwise_position = 10

    def test_local_spanwise_position_can_be_set_once(self):
        """Test that local_spanwise_position can be set once."""
        self.assertIsNone(self.basic_panel.local_spanwise_position)
        self.basic_panel.local_spanwise_position = 3
        self.assertEqual(self.basic_panel.local_spanwise_position, 3)

    def test_local_spanwise_position_cannot_be_set_twice(self):
        """Test that local_spanwise_position cannot be set twice."""
        self.basic_panel.local_spanwise_position = 3
        with self.assertRaises(AttributeError):
            self.basic_panel.local_spanwise_position = 7


class TestPanelDeepCopy(unittest.TestCase):
    """Tests for Panel.__deepcopy__ method."""

    def test_deepcopy_creates_new_instance(self):
        """Test that deepcopy creates a new Panel instance."""
        original = panel_fixtures.make_basic_panel_fixture()
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, _panel.Panel)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_corner_positions(self):
        """Test that deepcopy preserves corner position values."""
        original = panel_fixtures.make_basic_panel_fixture()
        copied = copy.deepcopy(original)

        npt.assert_array_equal(copied.Frpp_G_Cg, original.Frpp_G_Cg)
        npt.assert_array_equal(copied.Flpp_G_Cg, original.Flpp_G_Cg)
        npt.assert_array_equal(copied.Blpp_G_Cg, original.Blpp_G_Cg)
        npt.assert_array_equal(copied.Brpp_G_Cg, original.Brpp_G_Cg)

    def test_deepcopy_creates_independent_corner_arrays(self):
        """Test that deepcopy creates independent copies of corner position arrays."""
        original = panel_fixtures.make_basic_panel_fixture()
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.Frpp_G_Cg, original.Frpp_G_Cg)
        self.assertIsNot(copied.Flpp_G_Cg, original.Flpp_G_Cg)
        self.assertIsNot(copied.Blpp_G_Cg, original.Blpp_G_Cg)
        self.assertIsNot(copied.Brpp_G_Cg, original.Brpp_G_Cg)

    def test_deepcopy_corner_arrays_are_read_only(self):
        """Test that deepcopied corner position arrays are read only."""
        original = panel_fixtures.make_basic_panel_fixture()
        copied = copy.deepcopy(original)

        with self.assertRaises(ValueError):
            copied.Frpp_G_Cg[0] = 999.0
        with self.assertRaises(ValueError):
            copied.Flpp_G_Cg[0] = 999.0
        with self.assertRaises(ValueError):
            copied.Blpp_G_Cg[0] = 999.0
        with self.assertRaises(ValueError):
            copied.Brpp_G_Cg[0] = 999.0

    def test_deepcopy_preserves_edge_flags(self):
        """Test that deepcopy preserves is_leading_edge and is_trailing_edge."""
        original = panel_fixtures.make_leading_and_trailing_edge_panel_fixture()
        copied = copy.deepcopy(original)

        self.assertEqual(copied.is_leading_edge, original.is_leading_edge)
        self.assertEqual(copied.is_trailing_edge, original.is_trailing_edge)

    def test_deepcopy_preserves_mesh_position_attributes(self):
        """Test that deepcopy preserves mesh position attributes."""
        original = panel_fixtures.make_panel_with_set_once_attributes_fixture()
        copied = copy.deepcopy(original)

        self.assertEqual(copied.is_right_edge, original.is_right_edge)
        self.assertEqual(copied.is_left_edge, original.is_left_edge)
        self.assertEqual(
            copied.local_chordwise_position, original.local_chordwise_position
        )
        self.assertEqual(
            copied.local_spanwise_position, original.local_spanwise_position
        )

    def test_deepcopy_preserves_cached_local_properties(self):
        """Test that deepcopy preserves cached local geometric properties."""
        original = panel_fixtures.make_panel_with_cached_properties_fixture()
        copied = copy.deepcopy(original)

        npt.assert_array_equal(copied.rightLeg_G, original.rightLeg_G)
        npt.assert_array_equal(copied.frontLeg_G, original.frontLeg_G)
        npt.assert_array_equal(copied.leftLeg_G, original.leftLeg_G)
        npt.assert_array_equal(copied.backLeg_G, original.backLeg_G)
        npt.assert_array_equal(copied.Frbvp_G_Cg, original.Frbvp_G_Cg)
        npt.assert_array_equal(copied.Flbvp_G_Cg, original.Flbvp_G_Cg)
        npt.assert_array_equal(copied.Cpp_G_Cg, original.Cpp_G_Cg)
        npt.assert_array_equal(copied.unitNormal_G, original.unitNormal_G)
        self.assertEqual(copied.area, original.area)
        self.assertEqual(copied.aspect_ratio, original.aspect_ratio)

    def test_deepcopy_resets_gp1_positions_to_none(self):
        """Test that deepcopy resets GP1 positions to None."""
        original = panel_fixtures.make_panel_with_gp1_positions_fixture()

        # Verify original has GP1 positions set.
        self.assertIsNotNone(original.Frpp_GP1_CgP1)
        self.assertIsNotNone(original.Flpp_GP1_CgP1)
        self.assertIsNotNone(original.Blpp_GP1_CgP1)
        self.assertIsNotNone(original.Brpp_GP1_CgP1)

        copied = copy.deepcopy(original)

        # Verify copy has GP1 positions reset to None.
        self.assertIsNone(copied.Frpp_GP1_CgP1)
        self.assertIsNone(copied.Flpp_GP1_CgP1)
        self.assertIsNone(copied.Blpp_GP1_CgP1)
        self.assertIsNone(copied.Brpp_GP1_CgP1)

    def test_deepcopy_resets_gp1_derived_properties_to_none(self):
        """Test that deepcopy resets GP1 derived properties to None."""
        original = panel_fixtures.make_panel_with_gp1_positions_fixture()

        # Access GP1 derived properties to cache them.
        _ = original.rightLeg_GP1
        _ = original.frontLeg_GP1
        _ = original.leftLeg_GP1
        _ = original.backLeg_GP1
        _ = original.Frbvp_GP1_CgP1
        _ = original.Flbvp_GP1_CgP1
        _ = original.Cpp_GP1_CgP1
        _ = original.unitNormal_GP1

        copied = copy.deepcopy(original)

        # Verify copy has GP1 derived properties reset to None.
        self.assertIsNone(copied.rightLeg_GP1)
        self.assertIsNone(copied.frontLeg_GP1)
        self.assertIsNone(copied.leftLeg_GP1)
        self.assertIsNone(copied.backLeg_GP1)
        self.assertIsNone(copied.Frbvp_GP1_CgP1)
        self.assertIsNone(copied.Flbvp_GP1_CgP1)
        self.assertIsNone(copied.Cpp_GP1_CgP1)
        self.assertIsNone(copied.unitNormal_GP1)

    def test_deepcopy_resets_loads_to_none(self):
        """Test that deepcopy resets force and moment attributes to None."""
        original = panel_fixtures.make_basic_panel_fixture()

        # Set mock values for loads.
        original.forces_GP1 = np.array([1.0, 2.0, 3.0])
        original.moments_GP1_CgP1 = np.array([4.0, 5.0, 6.0])
        original.forces_W = np.array([7.0, 8.0, 9.0])
        original.moments_W_CgP1 = np.array([10.0, 11.0, 12.0])

        copied = copy.deepcopy(original)

        self.assertIsNone(copied.forces_GP1)
        self.assertIsNone(copied.moments_GP1_CgP1)
        self.assertIsNone(copied.forces_W)
        self.assertIsNone(copied.moments_W_CgP1)

    def test_deepcopy_allows_setting_gp1_positions_on_copy(self):
        """Test that deepcopy allows setting GP1 positions on the copy."""
        original = panel_fixtures.make_panel_with_gp1_positions_fixture()
        copied = copy.deepcopy(original)

        # Should be able to set GP1 positions on copy since they were reset.
        new_position = np.array([100.0, 200.0, 300.0])
        copied.Frpp_GP1_CgP1 = new_position.copy()
        npt.assert_array_equal(copied.Frpp_GP1_CgP1, new_position)

    def test_deepcopy_with_uncached_properties(self):
        """Test deepcopy when no cached properties have been accessed."""
        original = panel_fixtures.make_basic_panel_fixture()
        copied = copy.deepcopy(original)

        # Access properties on copy and verify they are correct.
        self.assertIsNotNone(copied.rightLeg_G)
        self.assertIsNotNone(copied.area)
        npt.assert_array_equal(copied.rightLeg_G, original.rightLeg_G)
        self.assertEqual(copied.area, original.area)


if __name__ == "__main__":
    unittest.main()
