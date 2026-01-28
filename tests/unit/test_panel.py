"""This module contains classes to test Panels."""

import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _panel
from tests.unit.fixtures import geometry_fixtures


class TestPanel(unittest.TestCase):
    """This class contains unit tests for the Panel class."""

    def setUp(self):
        """Set up test fixtures for Panel tests."""
        self.basic_panel = geometry_fixtures.make_basic_panel_fixture()

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

        # Test that vortex attributes are None initially
        self.assertIsNone(panel.ring_vortex)
        self.assertIsNone(panel.horseshoe_vortex)

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

        # Should be area * cos(45°) = area / sqrt(2)
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
        self.basic_panel = geometry_fixtures.make_basic_panel_fixture()

        # Set up global positions (simulating what a solver would do)
        # Use an offset to distinguish from local coordinates
        self.offset = np.array([10.0, 20.0, 5.0])
        self.basic_panel.Frpp_GP1_CgP1 = self.basic_panel.Frpp_G_Cg + self.offset
        self.basic_panel.Flpp_GP1_CgP1 = self.basic_panel.Flpp_G_Cg + self.offset
        self.basic_panel.Blpp_GP1_CgP1 = self.basic_panel.Blpp_G_Cg + self.offset
        self.basic_panel.Brpp_GP1_CgP1 = self.basic_panel.Brpp_G_Cg + self.offset

    def test_rightLeg_GP1_property(self):
        """Test right leg vector calculation in GP1 coordinates."""
        panel = self.basic_panel

        rightLeg_GP1 = panel.rightLeg_GP1

        # Right leg should be same as local since it's a vector (translation invariant)
        expected_rightLeg_GP1 = np.array([-1.0, 0.0, 0.0])
        npt.assert_array_almost_equal(rightLeg_GP1, expected_rightLeg_GP1)

    def test_rightLeg_GP1_returns_none_when_not_set(self):
        """Test that rightLeg_GP1 returns None when GP1 positions not set."""
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

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
        panel = geometry_fixtures.make_basic_panel_fixture()

        self.assertIsNone(panel.unitNormal_GP1)

    def test_GP1_properties_with_rotated_panel(self):
        """Test GP1 properties when global position includes rotation effect."""
        panel = geometry_fixtures.make_basic_panel_fixture()

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


if __name__ == "__main__":
    unittest.main()
