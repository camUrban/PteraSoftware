"""This module contains a class to test Panels.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    TestPanel: This is a class with functions to test Panels.
"""

import unittest
import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestPanel(unittest.TestCase):
    """This class contains unit tests for the Panel class."""

    def setUp(self):
        """Set up test fixtures for Panel tests."""
        self.basic_panel = geometry_fixtures.make_basic_panel_fixture()

    def tearDown(self):
        """Clean up test fixtures."""
        del self.basic_panel

    def test_initialization_valid_parameters(self):
        """Test Panel initialization with valid parameters."""
        panel = self.basic_panel

        # Test that Panel initializes correctly
        self.assertIsInstance(panel, ps.geometry.panel.Panel)

        # Test that corner points are correctly stored
        npt.assert_array_equal(panel.Frpp_G_Cg, np.array([0.0, 0.5, 0.0]))
        npt.assert_array_equal(panel.Flpp_G_Cg, np.array([0.0, 0.0, 0.0]))
        npt.assert_array_equal(panel.Blpp_G_Cg, np.array([1.0, 0.0, 0.0]))
        npt.assert_array_equal(panel.Brpp_G_Cg, np.array([1.0, 0.5, 0.0]))

        # Test that edge booleans are correctly stored
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
        self.assertIsNone(panel.forces_G)
        self.assertIsNone(panel.moments_G_Cg)
        self.assertIsNone(panel.forces_W)
        self.assertIsNone(panel.moments_W_Cg)

        # Test that coefficient attributes are None initially
        self.assertIsNone(panel.induced_drag_coefficient)
        self.assertIsNone(panel.side_force_coefficient)
        self.assertIsNone(panel.lift_coefficient)

    def test_parameter_validation_corner_points(self):
        """Test parameter validation for corner point inputs."""
        # Test invalid Frpp_G_Cg type
        with self.assertRaises(TypeError):
            ps.geometry.panel.Panel(
                Frpp_G_Cg="invalid",
                Flpp_G_Cg=[0.0, 0.0, 0.0],
                Blpp_G_Cg=[1.0, 0.0, 0.0],
                Brpp_G_Cg=[1.0, 0.5, 0.0],
                is_leading_edge=False,
                is_trailing_edge=False,
            )

        # Test invalid Frpp_G_Cg size
        with self.assertRaises((ValueError, TypeError)):
            ps.geometry.panel.Panel(
                Frpp_G_Cg=[0.0, 0.5],  # Only 2 elements
                Flpp_G_Cg=[0.0, 0.0, 0.0],
                Blpp_G_Cg=[1.0, 0.0, 0.0],
                Brpp_G_Cg=[1.0, 0.5, 0.0],
                is_leading_edge=False,
                is_trailing_edge=False,
            )

    def test_parameter_validation_edge_booleans(self):
        """Test parameter validation for edge boolean inputs."""
        # Test invalid is_leading_edge type
        with self.assertRaises(TypeError):
            ps.geometry.panel.Panel(
                Frpp_G_Cg=[0.0, 0.5, 0.0],
                Flpp_G_Cg=[0.0, 0.0, 0.0],
                Blpp_G_Cg=[1.0, 0.0, 0.0],
                Brpp_G_Cg=[1.0, 0.5, 0.0],
                is_leading_edge="invalid",
                is_trailing_edge=False,
            )

        # Test invalid is_trailing_edge type
        with self.assertRaises(TypeError):
            ps.geometry.panel.Panel(
                Frpp_G_Cg=[0.0, 0.5, 0.0],
                Flpp_G_Cg=[0.0, 0.0, 0.0],
                Blpp_G_Cg=[1.0, 0.0, 0.0],
                Brpp_G_Cg=[1.0, 0.5, 0.0],
                is_leading_edge=False,
                is_trailing_edge=123,
            )

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

    def test_unitSpanwise_G_property(self):
        """Test unit spanwise vector calculation."""
        panel = self.basic_panel

        unitSpanwise_G = panel.unitSpanwise_G

        # For this panel, spanwise should point in +y direction
        expected_unitSpanwise_G = np.array([0.0, 1.0, 0.0])
        npt.assert_array_almost_equal(unitSpanwise_G, expected_unitSpanwise_G)

        # Verify it's a unit vector
        self.assertAlmostEqual(np.linalg.norm(unitSpanwise_G), 1.0, places=10)

    def test_unitChordwise_G_property(self):
        """Test unit chordwise vector calculation."""
        panel = self.basic_panel

        unitChordwise_G = panel.unitChordwise_G

        # For this panel, chordwise should point in +x direction
        expected_unitChordwise_G = np.array([1.0, 0.0, 0.0])
        npt.assert_array_almost_equal(unitChordwise_G, expected_unitChordwise_G)

        # Verify it's a unit vector
        self.assertAlmostEqual(np.linalg.norm(unitChordwise_G), 1.0, places=10)

    def test_average_span_property(self):
        """Test average span calculation."""
        panel = self.basic_panel

        average_span = panel.average_span

        # For a rectangular panel with span 0.5 m on both front and back
        expected_average_span = 0.5
        self.assertAlmostEqual(average_span, expected_average_span, places=10)

    def test_average_chord_property(self):
        """Test average chord calculation."""
        panel = self.basic_panel

        average_chord = panel.average_chord

        # For a rectangular panel with chord 1.0 m on both left and right
        expected_average_chord = 1.0
        self.assertAlmostEqual(average_chord, expected_average_chord, places=10)

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
        """Test projected area at 45 degree angle."""
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

    def test_calculate_projected_area_validation(self):
        """Test validation of normal vector parameter."""
        panel = self.basic_panel

        # Test invalid normal_G type
        with self.assertRaises(TypeError):
            panel.calculate_projected_area("invalid")

        # Test invalid normal_G size
        with self.assertRaises((ValueError, TypeError)):
            panel.calculate_projected_area([0.0, 1.0])  # Only 2 elements

    def test_update_coefficients(self):
        """Test force coefficient updates."""
        panel = self.basic_panel

        # Set forces in wind axes for testing
        panel.forces_W = np.array([-10.0, 2.0, -50.0])

        # Calculate coefficients with a test dynamic pressure
        freestream_q = 100.0  # Pascals
        panel.update_coefficients(freestream_q)

        # Test induced drag coefficient: -forces_W[0] / area / q
        expected_induced_drag_coeff = 10.0 / panel.area / freestream_q
        self.assertAlmostEqual(
            panel.induced_drag_coefficient, expected_induced_drag_coeff, places=10
        )

        # Test side force coefficient: forces_W[1] / area / q
        expected_side_force_coeff = 2.0 / panel.area / freestream_q
        self.assertAlmostEqual(
            panel.side_force_coefficient, expected_side_force_coeff, places=10
        )

        # Test lift coefficient: -forces_W[2] / area / q
        expected_lift_coeff = 50.0 / panel.area / freestream_q
        self.assertAlmostEqual(
            panel.lift_coefficient, expected_lift_coeff, places=10
        )

    def test_update_coefficients_validation(self):
        """Test validation of freestream_q parameter."""
        panel = self.basic_panel
        panel.forces_W = np.array([1.0, 0.0, 1.0])

        # Test invalid freestream_q type
        with self.assertRaises(TypeError):
            panel.update_coefficients("invalid")

        # Test negative freestream_q
        with self.assertRaises((ValueError, TypeError)):
            panel.update_coefficients(-100.0)

        # Test zero freestream_q
        with self.assertRaises((ValueError, TypeError, ZeroDivisionError)):
            panel.update_coefficients(0.0)

    def test_nearly_planar_panel(self):
        """Test with a nearly planar panel."""
        # Create a panel with very slight twist
        panel = ps.geometry.panel.Panel(
            Frpp_G_Cg=[0.0, 1.0, 0.001],
            Flpp_G_Cg=[0.0, 0.0, 0.0],
            Blpp_G_Cg=[2.0, 0.0, 0.0],
            Brpp_G_Cg=[2.0, 1.0, 0.002],
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
        panel = ps.geometry.panel.Panel(
            Frpp_G_Cg=[0.0, 1.0, 0.5],
            Flpp_G_Cg=[0.0, 0.0, 0.0],
            Blpp_G_Cg=[2.0, 0.0, -0.5],
            Brpp_G_Cg=[2.0, 1.0, 0.0],
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
        panel = ps.geometry.panel.Panel(
            Frpp_G_Cg=[0.00, 0.01, 0.0],
            Flpp_G_Cg=[0.00, 0.00, 0.0],
            Blpp_G_Cg=[0.01, 0.00, 0.0],
            Brpp_G_Cg=[0.01, 0.01, 0.0],
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Should still calculate properties correctly
        area = panel.area
        expected_area = 0.0001
        self.assertAlmostEqual(area, expected_area, places=10)

        unitNormal_G = panel.unitNormal_G
        self.assertAlmostEqual(np.linalg.norm(unitNormal_G), 1.0, places=10)

    def test_tapered_panel(self):
        """Test with a tapered (trapezoidal) panel."""
        # Create a trapezoidal panel with different front and back spans
        panel = ps.geometry.panel.Panel(
            Frpp_G_Cg=[0.0, 1.0, 0.0],  # Front span: 1.0 m
            Flpp_G_Cg=[0.0, 0.0, 0.0],
            Blpp_G_Cg=[2.0, 0.0, 0.0],
            Brpp_G_Cg=[2.0, 0.5, 0.0],  # Back span: 0.5 m
            is_leading_edge=False,
            is_trailing_edge=False,
        )

        # Average span should be (1.0 + 0.5) / 2 = 0.75
        average_span = panel.average_span
        expected_average_span = 0.75
        self.assertAlmostEqual(average_span, expected_average_span, places=10)

        # Area should be approximately 1.5 square meters (trapezoid area)
        area = panel.area
        self.assertGreater(area, 0.0)


if __name__ == "__main__":
    unittest.main()