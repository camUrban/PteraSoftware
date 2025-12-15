"""This module contains a class to test WingCrossSections."""

import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestWingCrossSection(unittest.TestCase):
    """This is a class with functions to test WingCrossSections."""

    def setUp(self):
        """Set up test fixtures for WingCrossSection tests."""
        # Create fixtures using geometry_fixtures module
        self.test_airfoil = geometry_fixtures.make_test_airfoil_fixture()
        self.basic_wing_cross_section = (
            geometry_fixtures.make_basic_wing_cross_section_fixture(self.test_airfoil)
        )
        self.root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        self.tip_wing_cross_section = (
            geometry_fixtures.make_tip_wing_cross_section_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test WingCrossSection initialization with valid parameters."""
        # Test that basic WingCrossSection initializes correctly
        self.assertIsInstance(
            self.basic_wing_cross_section,
            ps.geometry.wing_cross_section.WingCrossSection,
        )
        self.assertEqual(self.basic_wing_cross_section.airfoil, self.test_airfoil)
        self.assertEqual(self.basic_wing_cross_section.num_spanwise_panels, 8)
        self.assertEqual(self.basic_wing_cross_section.chord, 1.5)
        np.testing.assert_array_equal(
            self.basic_wing_cross_section.Lp_Wcsp_Lpp, np.array([0.2, 0.5, 0.1])
        )
        np.testing.assert_array_equal(
            self.basic_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            np.array([5.0, -2.0, 3.0]),
        )
        self.assertEqual(
            self.basic_wing_cross_section.control_surface_symmetry_type, "symmetric"
        )
        self.assertEqual(
            self.basic_wing_cross_section.control_surface_hinge_point, 0.75
        )
        self.assertEqual(self.basic_wing_cross_section.control_surface_deflection, 5.0)
        self.assertEqual(self.basic_wing_cross_section.spanwise_spacing, "cosine")
        self.assertFalse(self.basic_wing_cross_section.validated)

    def test_initialization_default_parameters(self):
        """Test WingCrossSection initialization with default parameters."""
        wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=self.test_airfoil,
            num_spanwise_panels=5,
        )

        self.assertEqual(wing_cross_section.chord, 1.0)
        np.testing.assert_array_equal(
            wing_cross_section.Lp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0])
        )
        np.testing.assert_array_equal(
            wing_cross_section.angles_Wcsp_to_Wcs_ixyz, np.array([0.0, 0.0, 0.0])
        )
        self.assertEqual(wing_cross_section.control_surface_symmetry_type, None)
        self.assertEqual(wing_cross_section.control_surface_hinge_point, 0.75)
        self.assertEqual(wing_cross_section.control_surface_deflection, 0.0)
        self.assertIsNone(wing_cross_section.spanwise_spacing)

    def test_airfoil_parameter_validation(self):
        """Test that airfoil parameter is properly validated."""
        # Test with invalid airfoil type
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing_cross_section.WingCrossSection(
                airfoil="not_an_airfoil",
                num_spanwise_panels=8,
            )

    def test_num_spanwise_panels_validation(self):
        """Test num_spanwise_panels parameter validation."""
        # Test with valid positive integer
        wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=self.test_airfoil,
            num_spanwise_panels=15,
        )
        self.assertEqual(wing_cross_section.num_spanwise_panels, 15)

        # Test with None (valid for tip cross sections)
        wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=self.test_airfoil,
            num_spanwise_panels=None,
        )
        self.assertIsNone(wing_cross_section.num_spanwise_panels)

        # Test with invalid values
        invalid_values = [0, -5, 2.5, "eight"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=invalid_value,
                    )

    def test_chord_validation(self):
        """Test chord parameter validation."""
        # Test with valid positive values
        valid_chords = [0.1, 1.0, 2.5, 10.0]
        for chord in valid_chords:
            with self.subTest(chord=chord):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8,
                    chord=chord,
                )
                self.assertEqual(wing_cross_section.chord, chord)

        # Test with invalid values
        invalid_chords = [0.0, -1.0, -0.5, "one"]
        for invalid_chord in invalid_chords:
            with self.subTest(invalid_chord=invalid_chord):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        chord=invalid_chord,
                    )

    def test_Lp_Wcsp_Lpp_validation(self):
        """Test Lp_Wcsp_Lpp parameter validation."""
        # Test with valid array-like inputs
        valid_vectors = [
            np.array([0.0, 0.0, 0.0]),  # numpy array of floats
            [1.0, 2.0, 0.5],  # list of floats
            [1, 2, 0],  # list of ints
            (-0.5, 1.0, -0.2),  # tuple of floats
            np.array([1, 2, 0]),  # numpy array of ints
        ]
        for vector in valid_vectors:
            with self.subTest(vector=vector):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8,
                    Lp_Wcsp_Lpp=vector,
                )
                np.testing.assert_array_equal(wing_cross_section.Lp_Wcsp_Lpp, vector)

        # Test that second component must be non-negative
        with self.assertRaises(ValueError):
            ps.geometry.wing_cross_section.WingCrossSection(
                airfoil=self.test_airfoil,
                num_spanwise_panels=8,
                Lp_Wcsp_Lpp=np.array([1.0, -0.5, 0.0]),  # Negative y-component
            )

        # Test with invalid shapes/types
        invalid_vectors = [
            np.array([1.0, 2.0]),  # Wrong size
            np.array([1.0, 2.0, 3.0, 4.0]),  # Wrong size
            "not_a_vector",  # String
        ]
        for invalid_vector in invalid_vectors:
            with self.subTest(invalid_vector=invalid_vector):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        Lp_Wcsp_Lpp=invalid_vector,
                    )

    def test_angles_Wcsp_to_Wcs_ixyz_validation(self):
        """Test angles_Wcsp_to_Wcs_ixyz parameter validation."""
        # Test with valid array-like angles (within [-90, 90] range)
        valid_angles = [
            np.array([0.0, 0.0, 0.0]),  # numpy array of floats
            [45.0, -30.0, 60.0],  # list of floats
            [45, -30, 60],  # list of ints
            (89.9, -89.9, 0.0),  # tuple of floats
            np.array([30, -15, 45]),  # numpy array of ints
            [90.0, -90.0, 0.0],  # boundary values (list)
            np.array([90.0, 0.0, -90.0]),  # boundary values (array)
        ]
        for angles in valid_angles:
            with self.subTest(angles=angles):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8,
                    angles_Wcsp_to_Wcs_ixyz=angles,
                )
                np.testing.assert_array_equal(
                    wing_cross_section.angles_Wcsp_to_Wcs_ixyz, angles
                )

        # Test with angles outside valid range (using various array-like formats)
        invalid_angles = [
            [90.1, 0.0, 0.0],  # Greater than 90 (list)
            np.array([0.0, -90.1, 0.0]),  # Less than -90 (array)
            [95, 0, 0],  # Greater than 90 (list of ints)
            (0.0, -100.0, 0.0),  # Less than -90 (tuple)
        ]
        for invalid_angle in invalid_angles:
            with self.subTest(invalid_angle=invalid_angle):
                with self.assertRaises(ValueError):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        angles_Wcsp_to_Wcs_ixyz=invalid_angle,
                    )

    def test_control_surface_symmetry_type_validation(self):
        """Test control_surface_symmetry_type parameter validation."""
        # Test with valid types
        valid_types = ["symmetric", "asymmetric"]
        for control_type in valid_types:
            with self.subTest(control_type=control_type):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8,
                    control_surface_symmetry_type=control_type,
                )
                self.assertEqual(
                    wing_cross_section.control_surface_symmetry_type, control_type
                )

        # Test with invalid types
        invalid_types = ["invalid", "flap", "", 123]
        for invalid_type in invalid_types:
            with self.subTest(invalid_type=invalid_type):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        control_surface_symmetry_type=invalid_type,
                    )

    def test_control_surface_hinge_point_validation(self):
        """Test control_surface_hinge_point parameter validation."""
        # Test with valid values (in range 0 < x < 1)
        valid_hinge_points = [0.1, 0.5, 0.75, 0.9, 0.999]
        for hinge_point in valid_hinge_points:
            with self.subTest(hinge_point=hinge_point):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8,
                    control_surface_hinge_point=hinge_point,
                )
                self.assertEqual(
                    wing_cross_section.control_surface_hinge_point, hinge_point
                )

        # Test with invalid values (outside range or edge values)
        invalid_hinge_points = [0.0, 1.0, -0.1, 1.1, "point"]
        for invalid_hinge_point in invalid_hinge_points:
            with self.subTest(invalid_hinge_point=invalid_hinge_point):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        control_surface_hinge_point=invalid_hinge_point,
                    )

    def test_control_surface_deflection_validation(self):
        """Test control_surface_deflection parameter validation."""
        # Test with valid values (in range -5 <= x <= 5)
        valid_deflections = [0.0, 1.0, -5.0, 5.0, -4.1, 4]
        for deflection in valid_deflections:
            with self.subTest(deflection=deflection):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8,
                    control_surface_deflection=deflection,
                )
                self.assertEqual(
                    wing_cross_section.control_surface_deflection, deflection
                )

        # Test with invalid values (outside range or edge values)
        invalid_deflections = [-90.0, 90.0, -100.0, 120.0, "deflection"]
        for invalid_deflection in invalid_deflections:
            with self.subTest(invalid_deflection=invalid_deflection):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        control_surface_deflection=invalid_deflection,
                    )

    def test_spanwise_spacing_validation(self):
        """Test spanwise_spacing parameter validation."""
        # Test with valid values
        valid_spacings = ["cosine", "uniform", None]
        for spacing in valid_spacings:
            with self.subTest(spacing=spacing):
                wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=self.test_airfoil,
                    num_spanwise_panels=8 if spacing is not None else None,
                    spanwise_spacing=spacing,
                )
                self.assertEqual(wing_cross_section.spanwise_spacing, spacing)

        # Test with invalid values
        invalid_spacings = ["linear", "exponential", "", 123]
        for invalid_spacing in invalid_spacings:
            with self.subTest(invalid_spacing=invalid_spacing):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=self.test_airfoil,
                        num_spanwise_panels=8,
                        spanwise_spacing=invalid_spacing,
                    )

    def test_validate_root_constraints(self):
        """Test validate_root_constraints method."""
        # Test that root WingCrossSection passes validation
        try:
            self.root_wing_cross_section.validate_root_constraints()
        except Exception as e:
            self.fail(f"Root WingCrossSection validation failed unexpectedly: {e}")

        # Test that non-root WingCrossSections fail validation
        with self.assertRaises(ValueError):
            self.basic_wing_cross_section.validate_root_constraints()

        with self.assertRaises(ValueError):
            self.tip_wing_cross_section.validate_root_constraints()

        # Test that root WingCrossSection with num_spanwise_panels=None fails
        invalid_root_wcs = (
            geometry_fixtures.make_invalid_root_wing_cross_section_fixture()
        )
        with self.assertRaises(ValueError):
            invalid_root_wcs.validate_root_constraints()

    def test_validate_tip_constraints(self):
        """Test validate_tip_constraints method."""
        # Test that tip WingCrossSection passes validation
        try:
            self.tip_wing_cross_section.validate_tip_constraints()
        except Exception as e:
            self.fail(f"Tip WingCrossSection validation failed unexpectedly: {e}")

        # Test that non-tip WingCrossSections fail validation
        with self.assertRaises(ValueError):
            self.basic_wing_cross_section.validate_tip_constraints()

        with self.assertRaises(ValueError):
            self.root_wing_cross_section.validate_tip_constraints()

    def test_validate_mid_constraints(self):
        """Test validate_mid_constraints method."""
        # Test that valid middle WingCrossSection passes validation
        middle_wcs = geometry_fixtures.make_middle_wing_cross_section_fixture()
        try:
            middle_wcs.validate_mid_constraints()
        except Exception as e:
            self.fail(f"Middle WingCrossSection validation failed unexpectedly: {e}")

        # Test that basic WingCrossSection (has num_spanwise_panels) passes validation
        try:
            self.basic_wing_cross_section.validate_mid_constraints()
        except Exception as e:
            self.fail(
                f"Basic WingCrossSection validation as middle failed unexpectedly: {e}"
            )

        # Test that middle WingCrossSection with num_spanwise_panels=None fails
        invalid_middle_wcs = (
            geometry_fixtures.make_invalid_middle_wing_cross_section_fixture()
        )
        with self.assertRaises(ValueError):
            invalid_middle_wcs.validate_mid_constraints()

        # Test that tip WingCrossSection (has num_spanwise_panels=None) fails
        with self.assertRaises(ValueError):
            self.tip_wing_cross_section.validate_mid_constraints()

    def test_transformation_matrices_not_validated(self):
        """Test that transformation matrices return None when not validated."""
        # Test with unvalidated WingCrossSection
        self.assertIsNone(self.basic_wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp)
        self.assertIsNone(self.basic_wing_cross_section.T_pas_Wcs_Lp_to_Wcsp_Lpp)

    def test_transformation_matrices_validated(self):
        """Test that transformation matrices work correctly when validated."""
        # Manually set validated flag to True
        self.basic_wing_cross_section.validated = True

        # Test that matrices are returned
        T_forward = self.basic_wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp
        T_inverse = self.basic_wing_cross_section.T_pas_Wcs_Lp_to_Wcsp_Lpp

        self.assertIsInstance(T_forward, np.ndarray)
        self.assertIsInstance(T_inverse, np.ndarray)
        self.assertEqual(T_forward.shape, (4, 4))
        self.assertEqual(T_inverse.shape, (4, 4))

        # Test that they are inverses of each other
        identity = T_forward @ T_inverse
        expected_identity = np.eye(4)
        np.testing.assert_allclose(identity, expected_identity, rtol=1e-10, atol=1e-14)

        # Test with root WingCrossSection (should be identity matrices)
        self.root_wing_cross_section.validated = True
        T_root_forward = self.root_wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp
        T_root_inverse = self.root_wing_cross_section.T_pas_Wcs_Lp_to_Wcsp_Lpp

        np.testing.assert_allclose(T_root_forward, np.eye(4), rtol=1e-10, atol=1e-14)
        np.testing.assert_allclose(T_root_inverse, np.eye(4), rtol=1e-10, atol=1e-14)

    def test_comprehensive_initialization_edge_cases(self):
        """Test edge cases in initialization that combine multiple parameters."""
        # Test with minimal valid parameters
        minimal_wing_cross_section = (
            geometry_fixtures.make_minimal_wing_cross_section_fixture()
        )
        self.assertEqual(minimal_wing_cross_section.num_spanwise_panels, 1)

        # Test with maximum reasonable values
        max_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=self.test_airfoil,
            num_spanwise_panels=100,
            chord=50.0,
            Lp_Wcsp_Lpp=np.array([10.0, 20.0, 5.0]),
            angles_Wcsp_to_Wcs_ixyz=np.array([90.0, -90.0, 90.0]),
            control_surface_hinge_point=0.001,  # Very small but valid
            control_surface_deflection=5.0,  # Maximum valid deflection
        )
        self.assertEqual(max_wing_cross_section.num_spanwise_panels, 100)
        self.assertEqual(max_wing_cross_section.chord, 50.0)


if __name__ == "__main__":
    unittest.main()
