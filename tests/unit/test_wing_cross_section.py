"""This module contains classes to test WingCrossSections."""

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


class TestWingCrossSectionImmutability(unittest.TestCase):
    """Tests for WingCrossSection attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.test_airfoil = geometry_fixtures.make_test_airfoil_fixture()
        self.basic_wing_cross_section = (
            geometry_fixtures.make_basic_wing_cross_section_fixture(self.test_airfoil)
        )

    def test_immutable_airfoil_property(self):
        """Test that airfoil property is read only."""
        new_airfoil = geometry_fixtures.make_test_airfoil_fixture()
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.airfoil = new_airfoil

    def test_immutable_num_spanwise_panels_property(self):
        """Test that num_spanwise_panels property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.num_spanwise_panels = 20

    def test_immutable_chord_property(self):
        """Test that chord property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.chord = 2.5

    def test_immutable_Lp_Wcsp_Lpp_property(self):
        """Test that Lp_Wcsp_Lpp property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.Lp_Wcsp_Lpp = np.array([1.0, 2.0, 3.0])

    def test_immutable_Lp_Wcsp_Lpp_array_read_only(self):
        """Test that Lp_Wcsp_Lpp array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_wing_cross_section.Lp_Wcsp_Lpp[0] = 999.0

    def test_immutable_angles_Wcsp_to_Wcs_ixyz_property(self):
        """Test that angles_Wcsp_to_Wcs_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.angles_Wcsp_to_Wcs_ixyz = np.array(
                [1.0, 2.0, 3.0]
            )

    def test_immutable_angles_Wcsp_to_Wcs_ixyz_array_read_only(self):
        """Test that angles_Wcsp_to_Wcs_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_wing_cross_section.angles_Wcsp_to_Wcs_ixyz[0] = 999.0

    def test_immutable_control_surface_hinge_point_property(self):
        """Test that control_surface_hinge_point property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.control_surface_hinge_point = 0.5

    def test_immutable_control_surface_deflection_property(self):
        """Test that control_surface_deflection property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.control_surface_deflection = 2.0

    def test_immutable_spanwise_spacing_property(self):
        """Test that spanwise_spacing property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.spanwise_spacing = "uniform"

    def test_mutable_control_surface_symmetry_type(self):
        """Test that control_surface_symmetry_type remains mutable."""
        # This attribute must remain mutable for type 5 symmetry processing
        self.basic_wing_cross_section.control_surface_symmetry_type = "asymmetric"
        self.assertEqual(
            self.basic_wing_cross_section.control_surface_symmetry_type, "asymmetric"
        )

        self.basic_wing_cross_section.control_surface_symmetry_type = None
        self.assertIsNone(self.basic_wing_cross_section.control_surface_symmetry_type)

    def test_set_once_validated_property(self):
        """Test that validated can only be set once."""
        # First set should succeed
        self.assertFalse(self.basic_wing_cross_section.validated)
        self.basic_wing_cross_section.validated = True
        self.assertTrue(self.basic_wing_cross_section.validated)

        # Second set should raise AttributeError
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.validated = True

    def test_set_once_symmetry_type_property(self):
        """Test that symmetry_type can only be set once."""
        # First set should succeed
        self.assertIsNone(self.basic_wing_cross_section.symmetry_type)
        self.basic_wing_cross_section.symmetry_type = 4
        self.assertEqual(self.basic_wing_cross_section.symmetry_type, 4)

        # Second set should raise AttributeError
        with self.assertRaises(AttributeError):
            self.basic_wing_cross_section.symmetry_type = 3

    def test_validated_false_to_false_succeeds(self):
        """Test that setting validated from False to False succeeds.

        The set once logic only blocks when validated is already True. Setting
        from False to False is allowed. This test documents this behavior.
        """
        self.assertFalse(self.basic_wing_cross_section.validated)
        self.basic_wing_cross_section.validated = False
        self.assertFalse(self.basic_wing_cross_section.validated)

    def test_T_pas_Wcsp_Lpp_to_Wcs_Lp_array_read_only(self):
        """Test that T_pas_Wcsp_Lpp_to_Wcs_Lp array cannot be modified in place."""
        self.basic_wing_cross_section.validated = True
        T = self.basic_wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp
        self.assertIsNotNone(T)
        with self.assertRaises(ValueError):
            T[0, 0] = 999.0

    def test_T_pas_Wcs_Lp_to_Wcsp_Lpp_array_read_only(self):
        """Test that T_pas_Wcs_Lp_to_Wcsp_Lpp array cannot be modified in place."""
        self.basic_wing_cross_section.validated = True
        T = self.basic_wing_cross_section.T_pas_Wcs_Lp_to_Wcsp_Lpp
        self.assertIsNotNone(T)
        with self.assertRaises(ValueError):
            T[0, 0] = 999.0

    def test_T_pas_Wcsp_Lpp_to_Wcs_Lp_caching_returns_same_object(self):
        """Test that repeated access to T_pas_Wcsp_Lpp_to_Wcs_Lp returns the same
        cached object.
        """
        self.basic_wing_cross_section.validated = True
        T1 = self.basic_wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp
        T2 = self.basic_wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp
        self.assertIs(T1, T2)

    def test_T_pas_Wcs_Lp_to_Wcsp_Lpp_caching_returns_same_object(self):
        """Test that repeated access to T_pas_Wcs_Lp_to_Wcsp_Lpp returns the same
        cached object.
        """
        self.basic_wing_cross_section.validated = True
        T1 = self.basic_wing_cross_section.T_pas_Wcs_Lp_to_Wcsp_Lpp
        T2 = self.basic_wing_cross_section.T_pas_Wcs_Lp_to_Wcsp_Lpp
        self.assertIs(T1, T2)


class TestWingCrossSectionDeepCopy(unittest.TestCase):
    """Tests for WingCrossSection.__deepcopy__ method."""

    def setUp(self):
        """Set up test fixtures for deepcopy tests."""
        self.test_airfoil = geometry_fixtures.make_test_airfoil_fixture()
        self.basic_wing_cross_section = (
            geometry_fixtures.make_basic_wing_cross_section_fixture(self.test_airfoil)
        )

    def test_deepcopy_creates_new_instance(self):
        """Test that deepcopy creates a new WingCrossSection instance."""
        import copy

        original = self.basic_wing_cross_section
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.geometry.wing_cross_section.WingCrossSection)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_all_attributes(self):
        """Test that deepcopy preserves all attribute values."""
        import copy

        original = self.basic_wing_cross_section
        original.validated = True
        original.symmetry_type = 4

        copied = copy.deepcopy(original)

        self.assertEqual(copied.num_spanwise_panels, original.num_spanwise_panels)
        self.assertEqual(copied.chord, original.chord)
        np.testing.assert_array_equal(copied.Lp_Wcsp_Lpp, original.Lp_Wcsp_Lpp)
        np.testing.assert_array_equal(
            copied.angles_Wcsp_to_Wcs_ixyz, original.angles_Wcsp_to_Wcs_ixyz
        )
        self.assertEqual(
            copied.control_surface_symmetry_type, original.control_surface_symmetry_type
        )
        self.assertEqual(
            copied.control_surface_hinge_point, original.control_surface_hinge_point
        )
        self.assertEqual(
            copied.control_surface_deflection, original.control_surface_deflection
        )
        self.assertEqual(copied.spanwise_spacing, original.spanwise_spacing)
        self.assertEqual(copied.validated, original.validated)
        self.assertEqual(copied.symmetry_type, original.symmetry_type)

    def test_deepcopy_creates_independent_arrays(self):
        """Test that deepcopy creates independent copies of numpy arrays."""
        import copy

        original = self.basic_wing_cross_section
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.Lp_Wcsp_Lpp, original.Lp_Wcsp_Lpp)
        self.assertIsNot(
            copied.angles_Wcsp_to_Wcs_ixyz, original.angles_Wcsp_to_Wcs_ixyz
        )

    def test_deepcopy_independence_modifying_copy(self):
        """Test that immutable attributes cannot be modified on the copy."""
        import copy

        original = self.basic_wing_cross_section
        copied = copy.deepcopy(original)

        # Verify that immutable properties cannot be set
        with self.assertRaises(AttributeError):
            copied.chord = 999.0

        # Verify that numpy arrays are read only
        with self.assertRaises(ValueError):
            copied.Lp_Wcsp_Lpp[0] = 100.0

    def test_deepcopy_independence_modifying_original(self):
        """Test that immutable attributes cannot be modified on the original."""
        import copy

        original = self.basic_wing_cross_section
        copied = copy.deepcopy(original)

        # Verify that immutable properties cannot be set on the original
        with self.assertRaises(AttributeError):
            original.chord = 999.0

        # Verify that numpy arrays are read only on the original
        with self.assertRaises(ValueError):
            original.Lp_Wcsp_Lpp[0] = 100.0

        # Verify that the copy still has the original values
        self.assertEqual(copied.chord, original.chord)
        np.testing.assert_array_equal(copied.Lp_Wcsp_Lpp, original.Lp_Wcsp_Lpp)

    def test_deepcopy_with_none_values(self):
        """Test that deepcopy handles None values correctly."""
        import copy

        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        copied = copy.deepcopy(tip_wcs)

        self.assertIsNone(copied.num_spanwise_panels)
        self.assertIsNone(copied.spanwise_spacing)
        self.assertIsNone(copied.symmetry_type)

    def test_deepcopy_preserves_transformation_matrix_behavior(self):
        """Test that copied WingCrossSection has correct transformation matrix behavior."""
        import copy

        original = self.basic_wing_cross_section
        original.validated = True

        copied = copy.deepcopy(original)

        original_T = original.T_pas_Wcsp_Lpp_to_Wcs_Lp
        copied_T = copied.T_pas_Wcsp_Lpp_to_Wcs_Lp

        self.assertIsNotNone(original_T)
        self.assertIsNotNone(copied_T)
        np.testing.assert_array_almost_equal(copied_T, original_T)

    def test_deepcopy_unvalidated_cross_section(self):
        """Test that deepcopy handles unvalidated WingCrossSections correctly."""
        import copy

        original = self.basic_wing_cross_section
        self.assertFalse(original.validated)

        copied = copy.deepcopy(original)

        self.assertFalse(copied.validated)
        self.assertIsNone(copied.T_pas_Wcsp_Lpp_to_Wcs_Lp)

    def test_deepcopy_airfoil_independence(self):
        """Test that deepcopied Airfoil is a separate instance with immutable attributes."""
        import copy

        original = self.basic_wing_cross_section
        copied = copy.deepcopy(original)

        # Verify that the Airfoil is a separate instance
        self.assertIsNot(copied.airfoil, original.airfoil)

        # Verify that the Airfoil's name is immutable (no setter)
        with self.assertRaises(AttributeError):
            copied.airfoil.name = "modified_name"

        # Verify that the Airfoil arrays are independent and read only
        self.assertIsNot(copied.airfoil.outline_A_lp, original.airfoil.outline_A_lp)
        with self.assertRaises(ValueError):
            copied.airfoil.outline_A_lp[0, 0] = 999.0

    def test_deepcopy_cached_transformation_matrices_read_only(self):
        """Test that deepcopied cached transformation matrices are read only."""
        import copy

        original = self.basic_wing_cross_section
        original.validated = True

        # Trigger caching on original.
        _ = original.T_pas_Wcsp_Lpp_to_Wcs_Lp
        _ = original.T_pas_Wcs_Lp_to_Wcsp_Lpp

        copied = copy.deepcopy(original)

        # Verify that the copied matrices are read only.
        with self.assertRaises(ValueError):
            copied.T_pas_Wcsp_Lpp_to_Wcs_Lp[0, 0] = 999.0

        with self.assertRaises(ValueError):
            copied.T_pas_Wcs_Lp_to_Wcsp_Lpp[0, 0] = 999.0

    def test_deepcopy_cached_transformation_matrices_are_independent(self):
        """Test that deepcopied cached transformation matrices are independent copies."""
        import copy

        original = self.basic_wing_cross_section
        original.validated = True

        # Trigger caching on original.
        original_T_forward = original.T_pas_Wcsp_Lpp_to_Wcs_Lp
        original_T_inverse = original.T_pas_Wcs_Lp_to_Wcsp_Lpp

        copied = copy.deepcopy(original)

        # Verify that the copied matrices are different objects.
        self.assertIsNot(copied.T_pas_Wcsp_Lpp_to_Wcs_Lp, original_T_forward)
        self.assertIsNot(copied.T_pas_Wcs_Lp_to_Wcsp_Lpp, original_T_inverse)

        # Verify that the values are equal.
        np.testing.assert_array_equal(
            copied.T_pas_Wcsp_Lpp_to_Wcs_Lp, original_T_forward
        )
        np.testing.assert_array_equal(
            copied.T_pas_Wcs_Lp_to_Wcsp_Lpp, original_T_inverse
        )


class TestWingCrossSectionGetPlottableData(unittest.TestCase):
    """Tests for WingCrossSection.get_plottable_data method."""

    def setUp(self):
        """Set up test fixtures for get_plottable_data tests."""
        self.test_airfoil = geometry_fixtures.make_test_airfoil_fixture()
        self.basic_wing_cross_section = (
            geometry_fixtures.make_basic_wing_cross_section_fixture(self.test_airfoil)
        )
        self.root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )

    def test_get_plottable_data_returns_none_when_not_validated(self):
        """Test that get_plottable_data returns None when not validated."""
        # Set symmetry_type but not validated.
        self.basic_wing_cross_section.symmetry_type = 1
        self.assertFalse(self.basic_wing_cross_section.validated)

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        self.assertIsNone(result)

    def test_get_plottable_data_returns_none_when_symmetry_type_not_set(self):
        """Test that get_plottable_data returns None when symmetry_type not set."""
        # Set validated but not symmetry_type.
        self.basic_wing_cross_section.validated = True
        self.assertIsNone(self.basic_wing_cross_section.symmetry_type)

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        self.assertIsNone(result)

    def test_get_plottable_data_returns_none_when_neither_set(self):
        """Test that get_plottable_data returns None when neither validated nor
        symmetry_type is set.
        """
        self.assertFalse(self.basic_wing_cross_section.validated)
        self.assertIsNone(self.basic_wing_cross_section.symmetry_type)

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        self.assertIsNone(result)

    def test_get_plottable_data_returns_list_when_valid(self):
        """Test that get_plottable_data returns a list when validated and
        symmetry_type is set.
        """
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_returns_ndarrays(self):
        """Test that get_plottable_data returns ndarrays for outline and MCL."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        self.assertIsInstance(result[0], np.ndarray)
        self.assertIsInstance(result[1], np.ndarray)

    def test_get_plottable_data_returns_3d_points(self):
        """Test that get_plottable_data returns arrays with 3 columns (x, y, z)."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        # Both outline and MCL should have 3 columns (x, y, z).
        self.assertEqual(result[0].shape[1], 3)
        self.assertEqual(result[1].shape[1], 3)

    def test_get_plottable_data_y_components_are_zero(self):
        """Test that get_plottable_data returns points with zero y components.

        The points are in wing cross section axes relative to the leading point,
        so the y components should all be zero (the cross section is in the xz
        plane).
        """
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        # All y components should be zero.
        np.testing.assert_array_equal(result[0][:, 1], 0.0)
        np.testing.assert_array_equal(result[1][:, 1], 0.0)

    def test_get_plottable_data_scaled_by_chord(self):
        """Test that get_plottable_data returns points scaled by chord."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1
        chord = self.basic_wing_cross_section.chord

        result = self.basic_wing_cross_section.get_plottable_data(show=False)
        outline = result[0]

        # The x range should be approximately [0, chord].
        x_min = np.min(outline[:, 0])
        x_max = np.max(outline[:, 0])

        self.assertAlmostEqual(x_min, 0.0, places=5)
        self.assertAlmostEqual(x_max, chord, places=5)

    def test_get_plottable_data_default_show_is_false(self):
        """Test that get_plottable_data default for show is False."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        # Call without show parameter, should return data (not None).
        result = self.basic_wing_cross_section.get_plottable_data()

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_accepts_numpy_bool(self):
        """Test that get_plottable_data accepts numpy bool for show parameter."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=np.bool_(False))

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_with_symmetry_type_1(self):
        """Test get_plottable_data with symmetry type 1 (no symmetry)."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_with_root_cross_section(self):
        """Test get_plottable_data with a root WingCrossSection (identity transform)."""
        self.root_wing_cross_section.validated = True
        self.root_wing_cross_section.symmetry_type = 1

        result = self.root_wing_cross_section.get_plottable_data(show=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result), 2)

        # For root cross section, chord is 2.0.
        chord = self.root_wing_cross_section.chord
        outline = result[0]
        x_max = np.max(outline[:, 0])
        self.assertAlmostEqual(x_max, chord, places=5)

    def test_get_plottable_data_with_unit_chord(self):
        """Test get_plottable_data with a WingCrossSection with chord=1.0."""
        # Create a WingCrossSection with chord=1.0.
        wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=self.test_airfoil,
            num_spanwise_panels=8,
            chord=1.0,
        )
        wing_cross_section.validated = True
        wing_cross_section.symmetry_type = 1

        result = wing_cross_section.get_plottable_data(show=False)

        self.assertIsNotNone(result)
        outline = result[0]

        # With chord=1.0, the x range should be [0, 1].
        x_min = np.min(outline[:, 0])
        x_max = np.max(outline[:, 0])

        self.assertAlmostEqual(x_min, 0.0, places=5)
        self.assertAlmostEqual(x_max, 1.0, places=5)

    def test_get_plottable_data_mcl_within_outline_bounds(self):
        """Test that MCL points are within the outline x bounds."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        result = self.basic_wing_cross_section.get_plottable_data(show=False)
        outline = result[0]
        mcl = result[1]

        outline_x_min = np.min(outline[:, 0])
        outline_x_max = np.max(outline[:, 0])
        mcl_x_min = np.min(mcl[:, 0])
        mcl_x_max = np.max(mcl[:, 0])

        # MCL should be within outline x bounds.
        self.assertGreaterEqual(mcl_x_min, outline_x_min - 1e-10)
        self.assertLessEqual(mcl_x_max, outline_x_max + 1e-10)

    def test_get_plottable_data_invalid_show_type_raises(self):
        """Test that get_plottable_data raises error for invalid show type."""
        self.basic_wing_cross_section.validated = True
        self.basic_wing_cross_section.symmetry_type = 1

        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            self.basic_wing_cross_section.get_plottable_data(show="invalid")


if __name__ == "__main__":
    unittest.main()
