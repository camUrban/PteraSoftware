"""This module contains classes to test Wings."""

import copy
import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestWing(unittest.TestCase):
    """This class contains unit tests for the Wing class."""

    def setUp(self):
        """Set up test fixtures for Wing tests."""
        # Create fixtures for all Wing types
        self.type_1_wing = geometry_fixtures.make_type_1_wing_fixture()
        self.type_2_wing = geometry_fixtures.make_type_2_wing_fixture()
        self.type_3_wing = geometry_fixtures.make_type_3_wing_fixture()
        self.type_4_wing = geometry_fixtures.make_type_4_wing_fixture()
        self.type_5_wing = geometry_fixtures.make_type_5_wing_fixture()

        # Create additional test fixtures
        self.root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

    def test_initialization_valid_parameters(self):
        """Test Wing initialization with valid parameters for all types."""
        # Test that all Wing types initialize correctly
        wings_to_test = [
            (self.type_1_wing, "type 1"),
            (self.type_2_wing, "type 2"),
            (self.type_3_wing, "type 3"),
            (self.type_4_wing, "type 4"),
            (self.type_5_wing, "type 5"),
        ]

        for wing, wing_type in wings_to_test:
            with self.subTest(wing_type=wing_type):
                self.assertIsInstance(wing, ps.geometry.wing.Wing)
                self.assertIsInstance(wing.wing_cross_sections, tuple)
                self.assertEqual(len(wing.wing_cross_sections), 2)
                self.assertIsInstance(wing.name, str)
                self.assertEqual(len(wing.Ler_Gs_Cgs), 3)
                self.assertEqual(len(wing.angles_Gs_to_Wn_ixyz), 3)
                self.assertIsInstance(wing.symmetric, bool)
                self.assertIsInstance(wing.mirror_only, bool)
                self.assertEqual(wing.num_chordwise_panels, 8)
                self.assertEqual(wing.chordwise_spacing, "cosine")

    def test_wing_cross_sections_validation(self):
        """Test that wing_cross_sections parameter validation works correctly."""
        # Test empty list raises error
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(wing_cross_sections=[])

        # Test non-list raises error
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(wing_cross_sections="not a list")

        # Test single WingCrossSection raises error (need at least 2)
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(wing_cross_sections=[self.root_wcs])

        # Test non-WingCrossSection objects raise error
        with self.assertRaises(TypeError):
            ps.geometry.wing.Wing(wing_cross_sections=[self.root_wcs, "not a wcs"])

    def test_symmetry_parameter_validation(self):
        """Test symmetry parameter validation logic."""
        # Test that symmetric and mirror_only cannot both be True
        # Create fresh fixtures since WingCrossSections can only be validated once
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                symmetric=True,
                mirror_only=True,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )

        # Test that symmetry parameters must be None when no symmetry
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=[0.0, 1.0, 0.0],
            )

        # Test that symmetry parameters must be provided when symmetric=True
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=None,
            )

    def test_wing_type_1_properties(self):
        """Test type 1 wing (symmetric=False, mirror_only=False) properties."""
        wing = self.type_1_wing

        # Test basic properties
        self.assertFalse(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        self.assertIsNone(wing.symmetryNormal_G)
        self.assertIsNone(wing.symmetryPoint_G_Cg)

        # Test that symmetry_type is None before meshing
        self.assertIsNone(wing.symmetry_type)

    def test_wing_type_2_properties(self):
        """Test type 2 wing (mirror_only=True, coincident symmetry plane) properties."""
        wing = self.type_2_wing

        # Test basic properties
        self.assertFalse(wing.symmetric)
        self.assertTrue(wing.mirror_only)
        npt.assert_array_equal(wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([1.0, 0.0, 0.5]))

    def test_wing_type_3_properties(self):
        """Test type 3 wing (mirror_only=True, non-coincident symmetry plane)
        properties."""
        wing = self.type_3_wing

        # Test basic properties
        self.assertFalse(wing.symmetric)
        self.assertTrue(wing.mirror_only)
        npt.assert_array_equal(wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([0.0, -0.5, 0.0]))

    def test_wing_type_4_properties(self):
        """Test type 4 wing (symmetric=True, coincident symmetry plane) properties."""
        wing = self.type_4_wing

        # Test basic properties
        self.assertTrue(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        npt.assert_array_equal(wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([1.0, 0.0, 0.5]))

        # Test that WingCrossSections have control surface symmetry types
        for wcs in wing.wing_cross_sections:
            self.assertEqual(wcs.control_surface_symmetry_type, "symmetric")

    def test_wing_type_5_properties(self):
        """Test type 5 wing (symmetric=True, non-coincident symmetry plane)
        properties."""
        wing = self.type_5_wing

        # Test basic properties (before Airplane processing)
        self.assertTrue(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        npt.assert_array_equal(
            wing.symmetryNormal_G, np.array([0.0, np.sqrt(2) / 2, np.sqrt(2) / 2])
        )
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([0.5, 0.0, 0.0]))

    def test_generate_mesh_symmetry_type_1(self):
        """Test generate_mesh method with type 1 symmetry."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)
        self.assertEqual(wing.symmetry_type, 1)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_symmetry_type_2(self):
        """Test generate_mesh method with type 2 symmetry."""
        wing = geometry_fixtures.make_type_2_wing_fixture()
        wing.generate_mesh(2)
        self.assertEqual(wing.symmetry_type, 2)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_symmetry_type_3(self):
        """Test generate_mesh method with type 3 symmetry."""
        wing = geometry_fixtures.make_type_3_wing_fixture()
        wing.generate_mesh(3)
        self.assertEqual(wing.symmetry_type, 3)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_symmetry_type_4(self):
        """Test generate_mesh method with type 4 symmetry."""
        wing = geometry_fixtures.make_type_4_wing_fixture()
        wing.generate_mesh(4)
        self.assertEqual(wing.symmetry_type, 4)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_invalid_symmetry_type(self):
        """Test generate_mesh method with invalid symmetry types."""
        wing = self.type_1_wing

        invalid_types = [0, 5, -1, 10, "invalid", None]
        for invalid_type in invalid_types:
            with self.subTest(invalid_type=invalid_type):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    wing.generate_mesh(invalid_type)

    def test_transformation_matrices_before_meshing(self):
        """Test that transformation matrices are None before meshing."""
        wing = self.type_1_wing

        # All transformation matrices should be None before meshing
        self.assertIsNone(wing.T_pas_G_Cg_to_Wn_Ler)
        self.assertIsNone(wing.T_pas_Wn_Ler_to_G_Cg)

    def test_transformation_matrices_after_meshing(self):
        """Test that transformation matrices are available after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Transformation matrices should be available after meshing
        self.assertIsNotNone(wing.T_pas_G_Cg_to_Wn_Ler)
        self.assertIsNotNone(wing.T_pas_Wn_Ler_to_G_Cg)

        # Should be 4x4 matrices
        self.assertEqual(wing.T_pas_G_Cg_to_Wn_Ler.shape, (4, 4))
        self.assertEqual(wing.T_pas_Wn_Ler_to_G_Cg.shape, (4, 4))

        # Should be inverses of each other
        identity = wing.T_pas_G_Cg_to_Wn_Ler @ wing.T_pas_Wn_Ler_to_G_Cg
        npt.assert_allclose(identity, np.eye(4), atol=1e-14)

    def test_wing_axes_vectors_after_meshing(self):
        """Test wing axes vectors after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Wing axes basis vectors should be available
        self.assertIsNotNone(wing.WnX_G)
        self.assertIsNotNone(wing.WnY_G)
        self.assertIsNotNone(wing.WnZ_G)

        # Should be 3-element vectors
        self.assertEqual(len(wing.WnX_G), 3)
        self.assertEqual(len(wing.WnY_G), 3)
        self.assertEqual(len(wing.WnZ_G), 3)

        # Should form an orthonormal basis
        npt.assert_allclose(np.linalg.norm(wing.WnX_G), 1.0, atol=1e-14)
        npt.assert_allclose(np.linalg.norm(wing.WnY_G), 1.0, atol=1e-14)
        npt.assert_allclose(np.linalg.norm(wing.WnZ_G), 1.0, atol=1e-14)

        # Should be orthogonal
        npt.assert_allclose(np.dot(wing.WnX_G, wing.WnY_G), 0.0, atol=1e-14)
        npt.assert_allclose(wing.WnY_G @ wing.WnZ_G, 0.0, atol=1e-14)
        npt.assert_allclose(wing.WnZ_G @ wing.WnX_G, 0.0, atol=1e-14)

    def test_wing_cross_section_transformations_after_meshing(self):
        """Test WingCrossSection transformation matrices after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Should have transformation lists for each WingCrossSection
        num_wcs = len(wing.wing_cross_sections)
        self.assertEqual(len(wing.children_T_pas_Wn_Ler_to_Wcs_Lp), num_wcs)
        self.assertEqual(len(wing.children_T_pas_Wcs_Lp_to_Wn_Ler), num_wcs)
        self.assertEqual(len(wing.children_T_pas_G_Cg_to_Wcs_Lp), num_wcs)
        self.assertEqual(len(wing.children_T_pas_Wcs_Lp_to_G_Cg), num_wcs)

        # Each transformation should be a 4x4 matrix
        for i in range(num_wcs):
            self.assertEqual(wing.children_T_pas_Wn_Ler_to_Wcs_Lp[i].shape, (4, 4))
            self.assertEqual(wing.children_T_pas_Wcs_Lp_to_Wn_Ler[i].shape, (4, 4))
            self.assertEqual(wing.children_T_pas_G_Cg_to_Wcs_Lp[i].shape, (4, 4))
            self.assertEqual(wing.children_T_pas_Wcs_Lp_to_G_Cg[i].shape, (4, 4))

    def test_geometric_properties_after_meshing(self):
        """Test geometric property calculations after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Test that geometric properties are available and positive
        self.assertGreater(wing.projected_area, 0.0)
        self.assertGreater(wing.wetted_area, 0.0)
        self.assertGreater(wing.span, 0.0)
        self.assertGreater(wing.standard_mean_chord, 0.0)
        self.assertGreater(wing.mean_aerodynamic_chord, 0.0)

        # Test that wetted area is greater than projected area (both sides)
        self.assertGreaterEqual(wing.wetted_area, wing.projected_area)

    def test_geometric_properties_before_meshing_return_none(self):
        """Test that geometric properties return None before meshing."""
        wing = self.type_1_wing  # Not meshed

        properties_to_test = [
            "projected_area",
            "wetted_area",
            "span",
            "standard_mean_chord",
            "mean_aerodynamic_chord",
        ]

        for prop in properties_to_test:
            with self.subTest(property=prop):
                result = getattr(wing, prop)
                self.assertIsNone(result)

    def test_wing_with_different_chordwise_spacing(self):
        """Test Wing creation with different chordwise spacing options."""
        spacing_options = ["uniform", "cosine"]

        for spacing in spacing_options:
            with self.subTest(spacing=spacing):
                # Create fresh fixtures for each iteration since WingCrossSections
                # can only be validated once
                root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
                tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[root_wcs, tip_wcs],
                    chordwise_spacing=spacing,
                )
                self.assertEqual(wing.chordwise_spacing, spacing)

    def test_wing_with_different_chordwise_panels(self):
        """Test Wing creation with different numbers of chordwise panels."""
        panel_counts = [1, 4, 8, 16, 32]

        for count in panel_counts:
            with self.subTest(count=count):
                # Create fresh fixtures for each iteration since WingCrossSections
                # can only be validated once
                root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
                tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[root_wcs, tip_wcs],
                    num_chordwise_panels=count,
                )
                self.assertEqual(wing.num_chordwise_panels, count)

    def test_wing_parameter_validation(self):
        """Test parameter validation for Wing initialization."""
        # Test invalid Ler position
        # Create fresh fixtures since WingCrossSections can only be validated once
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs], Ler_Gs_Cgs="invalid"
            )

        # Test invalid angles
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                angles_Gs_to_Wn_ixyz="invalid",
            )

        # Test invalid num_chordwise_panels
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                num_chordwise_panels=0,
            )

        # Test invalid chordwise_spacing
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                chordwise_spacing="invalid_spacing",
            )

        # Test invalid symmetric
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs], symmetric="invalid"
            )

        # Test invalid mirror_only
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs], mirror_only="invalid"
            )

        # Test invalid explode_into_strips
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs], explode_into_strips="invalid"
            )

    def test_wing_name_validation(self):
        """Test Wing name parameter validation."""
        # Test valid string name
        # Create fresh fixtures since WingCrossSections can only be validated once
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wcs, tip_wcs], name="Test Wing Name"
        )
        self.assertEqual(wing.name, "Test Wing Name")

        # Test invalid name type
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(wing_cross_sections=[root_wcs, tip_wcs], name=123)

    def test_symmetry_normal_normalization(self):
        """Test that symmetry normal vectors are properly normalized."""
        # Create fresh fixtures since WingCrossSections can only be validated once
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        # Create Wing with non-unit normal vector
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wcs, tip_wcs],
            symmetric=False,
            mirror_only=True,
            symmetryNormal_G=[0.0, 5.0, 0.0],
            symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
        )

        # Should be normalized to unit vector
        npt.assert_allclose(np.linalg.norm(wing.symmetryNormal_G), 1.0, atol=1e-14)
        npt.assert_allclose(
            wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]), atol=1e-14
        )

    def test_three_section_wing_validation(self):
        """Test Wing with 3 WingCrossSections validates correctly."""
        # Test that valid 3-WingCrossSection Wing initializes correctly
        wing = geometry_fixtures.make_three_section_wing_fixture()
        self.assertIsInstance(wing, ps.geometry.wing.Wing)
        self.assertEqual(len(wing.wing_cross_sections), 3)

        # Verify all WingCrossSections are validated
        for wcs in wing.wing_cross_sections:
            self.assertTrue(wcs.validated)

    def test_four_section_wing_validation(self):
        """Test Wing with 4 WingCrossSections validates correctly."""
        # Test that valid 4-WingCrossSection Wing initializes correctly
        wing = geometry_fixtures.make_four_section_wing_fixture()
        self.assertIsInstance(wing, ps.geometry.wing.Wing)
        self.assertEqual(len(wing.wing_cross_sections), 4)

        # Verify all WingCrossSections are validated
        for wcs in wing.wing_cross_sections:
            self.assertTrue(wcs.validated)

    def test_invalid_middle_wing_cross_section_raises_error(self):
        """Test that Wing with invalid middle WingCrossSection raises ValueError."""
        # Test that Wing with middle WingCrossSection having num_spanwise_panels=None fails
        with self.assertRaises(ValueError):
            geometry_fixtures.make_invalid_three_section_wing_fixture()

    def test_invalid_root_wing_cross_section_raises_error(self):
        """Test that Wing with invalid root WingCrossSection raises ValueError."""
        # Test that Wing with root WingCrossSection having num_spanwise_panels=None fails
        with self.assertRaises(ValueError):
            geometry_fixtures.make_invalid_root_wing_fixture()

    def test_symmetry_point_none_when_symmetric_raises_value_error(self):
        """Test that symmetryPoint_G_Cg=None with symmetric=True raises ValueError."""
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=None,
            )

    def test_symmetry_point_not_none_when_no_symmetry_raises_value_error(self):
        """Test that symmetryPoint_G_Cg not None with no symmetry raises ValueError."""
        root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wcs, tip_wcs],
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=None,
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )

    def test_span_simple_rectangular_wing(self):
        """Test span calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (from y=0 to y=2.0)
        expected_span = 2.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_simple_tapered_wing(self):
        """Test span calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 3.0 meters (from y=0 to y=3.0)
        expected_span = 3.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_symmetric_continuous_rectangular_wing(self):
        """Test span calculation for symmetric continuous rectangular Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # Expected span: 5.0 meters (2 * 2.5, due to symmetry)
        expected_span = 5.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_three_section_tapered_wing(self):
        """Test span calculation for three section tapered Wing."""
        wing = geometry_fixtures.make_three_section_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 4.0 meters (from y=0 to y=4.0)
        expected_span = 4.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_projected_area_simple_rectangular_wing(self):
        """Test projected area calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Expected projected area: 2.0 square meters (1.0 m chord * 2.0 m span)
        expected_area = 2.0
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_projected_area_simple_tapered_wing(self):
        """Test projected area calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected projected area: 4.5 square meters
        # Trapezoid: (chord_root + chord_tip) / 2 * span = (2.0 + 1.0) / 2 * 3.0 = 4.5
        expected_area = 4.5
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_projected_area_symmetric_continuous_rectangular_wing(self):
        """Test projected area calculation for symmetric continuous rectangular
        Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # Expected projected area: 7.5 square meters
        # Rectangle: chord * span = 1.5 * 5.0 = 7.5
        expected_area = 7.5
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_projected_area_three_section_tapered_wing(self):
        """Test projected area calculation for three section tapered Wing."""
        wing = geometry_fixtures.make_three_section_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected projected area: 8.0 square meters
        # Section 1 (root to middle): (3.0 + 2.0) / 2 * 2.0 = 5.0
        # Section 2 (middle to tip): (2.0 + 1.0) / 2 * 2.0 = 3.0
        # Total: 5.0 + 3.0 = 8.0
        expected_area = 8.0
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_wetted_area_greater_than_projected_area(self):
        """Test that wetted area is greater than or equal to projected area for all
        Wings."""
        wings = [
            geometry_fixtures.make_simple_rectangular_wing_fixture(),
            geometry_fixtures.make_simple_tapered_wing_fixture(),
            geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture(),
            geometry_fixtures.make_three_section_tapered_wing_fixture(),
        ]

        symmetry_types = [1, 1, 4, 1]

        for wing, symmetry_type in zip(wings, symmetry_types):
            with self.subTest(wing=wing.name):
                wing.generate_mesh(symmetry_type)

                projected_area = wing.projected_area
                wetted_area = wing.wetted_area

                self.assertIsNotNone(projected_area)
                self.assertIsNotNone(wetted_area)
                self.assertGreaterEqual(wetted_area, projected_area)

    def test_standard_mean_chord_simple_rectangular_wing(self):
        """Test standard mean chord calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Standard mean chord = projected_area / span = 2.0 / 2.0 = 1.0
        expected_smc = 1.0
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_simple_tapered_wing(self):
        """Test standard mean chord calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Standard mean chord = projected_area / span = 4.5 / 3.0 = 1.5
        expected_smc = 1.5
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_symmetric_continuous_rectangular_wing(self):
        """Test standard mean chord calculation for symmetric continuous rectangular
        Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # Standard mean chord = projected_area / span = 7.5 / 5.0 = 1.5
        expected_smc = 1.5
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_three_section_tapered_wing(self):
        """Test standard mean chord calculation for three section tapered Wing."""
        wing = geometry_fixtures.make_three_section_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Standard mean chord = projected_area / span = 8.0 / 4.0 = 2.0
        expected_smc = 2.0
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_mean_aerodynamic_chord_simple_rectangular_wing(self):
        """Test mean aerodynamic chord calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # For a rectangular wing (constant chord), MAC = chord = 1.0
        expected_mac = 1.0
        actual_mac = wing.mean_aerodynamic_chord

        self.assertIsNotNone(actual_mac)
        npt.assert_allclose(actual_mac, expected_mac, rtol=1e-10, atol=1e-14)

    def test_mean_aerodynamic_chord_simple_tapered_wing(self):
        """Test mean aerodynamic chord calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # For a linearly tapered wing:
        # MAC = (2/3) * (c_root + c_tip - c_root * c_tip / (c_root + c_tip))
        # With c_root=2.0, c_tip=1.0:
        # MAC = (2/3) * (2.0 + 1.0 - 2.0 * 1.0 / (2.0 + 1.0))
        # MAC = (2/3) * (3.0 - 2.0/3.0) = (2/3) * (7.0/3.0) = 14.0/9.0 = 1.555...
        c_root = 2.0
        c_tip = 1.0
        expected_mac = (2.0 / 3.0) * (
            c_root + c_tip - c_root * c_tip / (c_root + c_tip)
        )

        actual_mac = wing.mean_aerodynamic_chord

        self.assertIsNotNone(actual_mac)
        npt.assert_allclose(actual_mac, expected_mac, rtol=1e-10, atol=1e-14)

    def test_mean_aerodynamic_chord_symmetric_continuous_rectangular_wing(self):
        """Test mean aerodynamic chord calculation for symmetric continuous
        rectangular Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # For a rectangular wing (constant chord), MAC = chord = 1.5
        expected_mac = 1.5
        actual_mac = wing.mean_aerodynamic_chord

        self.assertIsNotNone(actual_mac)
        npt.assert_allclose(actual_mac, expected_mac, rtol=1e-10, atol=1e-14)

    def test_geometric_properties_consistency(self):
        """Test that geometric properties are internally consistent across different
        Wings."""
        wings_data = [
            (geometry_fixtures.make_simple_rectangular_wing_fixture(), 1),
            (geometry_fixtures.make_simple_tapered_wing_fixture(), 1),
            (
                geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture(),
                4,
            ),
            (geometry_fixtures.make_three_section_tapered_wing_fixture(), 1),
        ]

        for wing, symmetry_type in wings_data:
            with self.subTest(wing=wing.name):
                wing.generate_mesh(symmetry_type)

                # Get properties
                span = wing.span
                projected_area = wing.projected_area
                standard_mean_chord = wing.standard_mean_chord

                # Verify consistency: projected_area = span * standard_mean_chord
                self.assertIsNotNone(span)
                self.assertIsNotNone(projected_area)
                self.assertIsNotNone(standard_mean_chord)

                calculated_area = span * standard_mean_chord
                npt.assert_allclose(
                    projected_area, calculated_area, rtol=1e-10, atol=1e-14
                )

    def test_span_rotated_wing_x_axis(self):
        """Test span calculation invariance for Wing rotated about x axis."""
        # Create a Wing rotated 45 degrees about x axis
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([45.0, 0.0, 0.0])
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotation about x axis does not affect y extent)
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_rotated_wing_y_axis(self):
        """Test span calculation invariance for Wing rotated about y axis."""
        # Create a Wing rotated 30 degrees about y axis
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([0.0, 30.0, 0.0])
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotation about y axis does not affect y extent)
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_rotated_wing_z_axis(self):
        """Test span calculation invariance for Wing rotated about z axis."""
        # Create a Wing rotated 60 degrees about z axis
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([0.0, 0.0, 60.0])
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotation about z axis does not affect y extent)
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_rotated_wing_combined_rotations(self):
        """Test span calculation invariance for Wing with combined rotations."""
        # Create a Wing with combined rotations
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture(
            [15.0, 25.0, 35.0]
        )
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotations do not affect y extent in wing axes)
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_wing_with_rotated_cross_sections(self):
        """Test span calculation for Wing with rotated WingCrossSections."""
        wing = geometry_fixtures.make_wing_with_rotated_cross_sections_fixture()
        wing.generate_mesh(1)

        # Expected span: 5.0 meters (rotation about y axis does not affect y position)
        expected_span = 5.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_swept_wing(self):
        """Test span calculation for swept Wing."""
        wing = geometry_fixtures.make_swept_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 3.0 meters (sweep in x direction does not affect y extent)
        expected_span = 3.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_dihedral_wing(self):
        """Test span calculation for Wing with dihedral."""
        wing = geometry_fixtures.make_dihedral_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 3.0 meters (dihedral in z direction does not affect y extent in wing axes)
        expected_span = 3.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_rotated_wing(self):
        """Test standard mean chord calculation for rotated Wing."""
        # Create a Wing rotated 45 degrees about x axis
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([45.0, 0.0, 0.0])
        wing.generate_mesh(1)

        # Expected standard mean chord: 1.0 (projected_area / span = 2.0 / 2.0)
        expected_smc = 1.0

        actual_smc = wing.standard_mean_chord
        self.assertIsNotNone(actual_smc)

        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_swept_wing(self):
        """Test standard mean chord calculation for swept Wing."""
        wing = geometry_fixtures.make_swept_wing_fixture()
        wing.generate_mesh(1)

        # Expected standard mean chord: projected_area / span = 4.5 / 3.0 = 1.5
        expected_smc = 1.5

        actual_smc = wing.standard_mean_chord
        self.assertIsNotNone(actual_smc)

        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_average_panel_aspect_ratio_returns_none_before_meshing(self):
        """Test that average_panel_aspect_ratio returns None before meshing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()

        self.assertIsNone(wing.average_panel_aspect_ratio)

    def test_average_panel_aspect_ratio_returns_positive_after_meshing(self):
        """Test that average_panel_aspect_ratio returns a positive value after meshing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        average_aspect_ratio = wing.average_panel_aspect_ratio
        self.assertIsNotNone(average_aspect_ratio)
        self.assertGreater(average_aspect_ratio, 0.0)

    def test_average_panel_aspect_ratio_simple_rectangular_wing(self):
        """Test average_panel_aspect_ratio calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        average_aspect_ratio = wing.average_panel_aspect_ratio
        self.assertIsNotNone(average_aspect_ratio)

        # For a rectangular wing with uniform spacing:
        # Panel chord = wing_chord / num_chordwise_panels = 1.0 / 4 = 0.25
        # Panel span = wing_span / num_spanwise_panels = 2.0 / 8 = 0.25
        # Expected aspect ratio ~ 1.0 for roughly square panels.
        # This is an approximate check since actual Panel aspect ratios depend on
        # the meshing algorithm's implementation details.
        self.assertGreater(average_aspect_ratio, 0.0)
        self.assertLess(average_aspect_ratio, 100.0)

    def test_average_panel_aspect_ratio_type_4_symmetric_wing(self):
        """Test average_panel_aspect_ratio for type 4 symmetric Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        average_aspect_ratio = wing.average_panel_aspect_ratio
        self.assertIsNotNone(average_aspect_ratio)
        self.assertGreater(average_aspect_ratio, 0.0)

    def test_average_panel_aspect_ratio_caching(self):
        """Test that average_panel_aspect_ratio is cached after first access."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Access twice.
        first_access = wing.average_panel_aspect_ratio
        second_access = wing.average_panel_aspect_ratio

        # Values should be identical.
        self.assertEqual(first_access, second_access)

    def test_angles_Gs_to_Wn_ixyz_validation_boundary_values(self):
        """Test angles_Gs_to_Wn_ixyz validation with boundary values."""
        # Test with boundary values (should be valid).
        valid_angles_sets = [
            [90.0, 0.0, 0.0],
            [0.0, 90.0, 0.0],
            [0.0, 0.0, 90.0],
            [-90.0, 0.0, 0.0],
            [0.0, -90.0, 0.0],
            [0.0, 0.0, -90.0],
            [90.0, 90.0, 90.0],
            [-90.0, -90.0, -90.0],
        ]

        for angles in valid_angles_sets:
            with self.subTest(angles=angles):
                root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
                tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[root_wcs, tip_wcs],
                    angles_Gs_to_Wn_ixyz=angles,
                )
                npt.assert_array_equal(wing.angles_Gs_to_Wn_ixyz, np.array(angles))

    def test_angles_Gs_to_Wn_ixyz_validation_outside_range(self):
        """Test angles_Gs_to_Wn_ixyz validation with values outside valid range."""
        # Test with values outside valid range (should raise ValueError).
        invalid_angles_sets = [
            [90.1, 0.0, 0.0],
            [0.0, 90.1, 0.0],
            [0.0, 0.0, 90.1],
            [-90.1, 0.0, 0.0],
            [0.0, -90.1, 0.0],
            [0.0, 0.0, -90.1],
            [100.0, 0.0, 0.0],
            [0.0, -100.0, 0.0],
        ]

        for angles in invalid_angles_sets:
            with self.subTest(angles=angles):
                root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
                tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
                with self.assertRaises(ValueError):
                    ps.geometry.wing.Wing(
                        wing_cross_sections=[root_wcs, tip_wcs],
                        angles_Gs_to_Wn_ixyz=angles,
                    )

    def test_Ler_Gs_Cgs_accepts_various_input_types(self):
        """Test that Ler_Gs_Cgs accepts various array-like input types."""
        input_formats = [
            np.array([1.0, 2.0, 3.0]),  # ndarray
            [1.0, 2.0, 3.0],  # list
            (1.0, 2.0, 3.0),  # tuple
            [1, 2, 3],  # list of ints
            (1, 2, 3),  # tuple of ints
        ]

        for input_val in input_formats:
            with self.subTest(input_format=type(input_val).__name__):
                root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
                tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[root_wcs, tip_wcs],
                    Ler_Gs_Cgs=input_val,
                )
                expected = np.array([1.0, 2.0, 3.0])
                npt.assert_array_equal(wing.Ler_Gs_Cgs, expected)

    def test_angles_Gs_to_Wn_ixyz_accepts_various_input_types(self):
        """Test that angles_Gs_to_Wn_ixyz accepts various array-like input types."""
        input_formats = [
            np.array([10.0, 20.0, 30.0]),  # ndarray
            [10.0, 20.0, 30.0],  # list
            (10.0, 20.0, 30.0),  # tuple
            [10, 20, 30],  # list of ints
            (10, 20, 30),  # tuple of ints
        ]

        for input_val in input_formats:
            with self.subTest(input_format=type(input_val).__name__):
                root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
                tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[root_wcs, tip_wcs],
                    angles_Gs_to_Wn_ixyz=input_val,
                )
                expected = np.array([10.0, 20.0, 30.0])
                npt.assert_array_equal(wing.angles_Gs_to_Wn_ixyz, expected)


class TestWingDeepCopy(unittest.TestCase):
    """Tests for Wing.__deepcopy__ method."""

    def setUp(self):
        """Set up test fixtures for deepcopy tests."""
        self.type_1_wing = geometry_fixtures.make_type_1_wing_fixture()
        self.type_4_wing = geometry_fixtures.make_type_4_wing_fixture()
        self.root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

    def test_deepcopy_creates_new_instance(self):
        """Test that deepcopy creates a new Wing instance."""
        import copy

        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.geometry.wing.Wing)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_wing_parameters(self):
        """Test that deepcopy preserves Wing parameters."""
        import copy

        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.symmetric, original.symmetric)
        self.assertEqual(copied.mirror_only, original.mirror_only)
        self.assertEqual(copied.num_chordwise_panels, original.num_chordwise_panels)
        self.assertEqual(copied.chordwise_spacing, original.chordwise_spacing)
        npt.assert_array_equal(copied.Ler_Gs_Cgs, original.Ler_Gs_Cgs)
        npt.assert_array_equal(
            copied.angles_Gs_to_Wn_ixyz, original.angles_Gs_to_Wn_ixyz
        )

    def test_deepcopy_creates_independent_arrays(self):
        """Test that deepcopy creates independent copies of numpy arrays."""
        import copy

        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.Ler_Gs_Cgs, original.Ler_Gs_Cgs)
        self.assertIsNot(copied.angles_Gs_to_Wn_ixyz, original.angles_Gs_to_Wn_ixyz)

    def test_deepcopy_creates_independent_wing_cross_sections(self):
        """Test that deepcopy creates independent WingCrossSection copies."""
        import copy

        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertEqual(
            len(copied.wing_cross_sections), len(original.wing_cross_sections)
        )
        for orig_wcs, copied_wcs in zip(
            original.wing_cross_sections, copied.wing_cross_sections
        ):
            self.assertIsNot(orig_wcs, copied_wcs)
            self.assertEqual(copied_wcs.chord, orig_wcs.chord)

    def test_deepcopy_preserves_symmetry_attributes(self):
        """Test that deepcopy preserves symmetry attributes correctly."""
        import copy

        original = self.type_4_wing
        copied = copy.deepcopy(original)

        self.assertEqual(copied.symmetric, original.symmetric)
        self.assertEqual(copied.mirror_only, original.mirror_only)
        npt.assert_array_equal(copied.symmetryNormal_G, original.symmetryNormal_G)
        npt.assert_array_equal(copied.symmetryPoint_G_Cg, original.symmetryPoint_G_Cg)
        self.assertIsNot(copied.symmetryNormal_G, original.symmetryNormal_G)
        self.assertIsNot(copied.symmetryPoint_G_Cg, original.symmetryPoint_G_Cg)

    def test_deepcopy_preserves_none_symmetry_attributes(self):
        """Test that deepcopy handles None symmetry attributes correctly."""
        import copy

        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertIsNone(copied.symmetryNormal_G)
        self.assertIsNone(copied.symmetryPoint_G_Cg)

    def test_deepcopy_unmeshed_wing(self):
        """Test that deepcopy handles unmeshed Wings correctly."""
        import copy

        original = self.type_1_wing
        self.assertIsNone(original.panels)

        copied = copy.deepcopy(original)

        self.assertIsNone(copied.symmetry_type)
        self.assertIsNone(copied.num_spanwise_panels)
        self.assertIsNone(copied.num_panels)
        self.assertIsNone(copied.panels)
        self.assertIsNone(copied.gridWrvp_GP1_CgP1)

    def test_deepcopy_meshed_wing_preserves_mesh_metadata(self):
        """Test that deepcopy preserves mesh metadata for meshed Wings."""
        import copy

        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertEqual(copied.symmetry_type, original.symmetry_type)
        self.assertEqual(copied.num_spanwise_panels, original.num_spanwise_panels)
        self.assertEqual(copied.num_panels, original.num_panels)

    def test_deepcopy_meshed_wing_preserves_panels(self):
        """Test that deepcopy preserves Panels for meshed Wings."""
        import copy

        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertIsNotNone(copied.panels)
        self.assertEqual(copied.panels.shape, original.panels.shape)

        for i in range(original.panels.shape[0]):
            for j in range(original.panels.shape[1]):
                orig_panel = original.panels[i, j]
                copied_panel = copied.panels[i, j]
                self.assertIsNot(orig_panel, copied_panel)
                npt.assert_array_equal(copied_panel.Frpp_G_Cg, orig_panel.Frpp_G_Cg)

    def test_deepcopy_resets_wake_state(self):
        """Test that deepcopy resets wake state to an empty array."""
        import copy

        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertIsNotNone(copied.gridWrvp_GP1_CgP1)
        self.assertEqual(copied.gridWrvp_GP1_CgP1.shape[0], 0)
        self.assertEqual(
            copied.gridWrvp_GP1_CgP1.shape[1], original.num_spanwise_panels + 1
        )

    def test_deepcopy_independence_modifying_copy_mutable_attrs(self):
        """Test that modifying mutable attributes on the copy does not affect the
        original."""
        import copy

        original = self.type_4_wing
        original.generate_mesh(4)
        original_symmetric = original.symmetric
        original_symmetryNormal = original.symmetryNormal_G.copy()

        copied = copy.deepcopy(original)

        # Modify mutable attributes on the copy.
        copied.symmetric = not original_symmetric
        copied.symmetryNormal_G[0] = 999.0

        # Verify the original is unchanged.
        self.assertEqual(original.symmetric, original_symmetric)
        npt.assert_array_equal(original.symmetryNormal_G, original_symmetryNormal)

    def test_deepcopy_independence_modifying_original_mutable_attrs(self):
        """Test that modifying mutable attributes on the original does not affect the
        copy."""
        import copy

        original = self.type_4_wing
        original.generate_mesh(4)

        copied = copy.deepcopy(original)
        copied_symmetric = copied.symmetric
        copied_symmetryNormal = copied.symmetryNormal_G.copy()

        # Modify mutable attributes on the original.
        original.symmetric = not copied_symmetric
        original.symmetryNormal_G[0] = 999.0

        # Verify the copy is unchanged.
        self.assertEqual(copied.symmetric, copied_symmetric)
        npt.assert_array_equal(copied.symmetryNormal_G, copied_symmetryNormal)

    def test_immutable_attributes_are_read_only(self):
        """Test that immutable Wing attributes cannot be modified."""
        wing = self.type_1_wing
        wing.generate_mesh(1)

        # Test that setting immutable attributes raises AttributeError.
        with self.assertRaises(AttributeError):
            wing.name = "New Name"

        with self.assertRaises(AttributeError):
            wing.num_chordwise_panels = 16

        with self.assertRaises(AttributeError):
            wing.chordwise_spacing = "uniform"

        # Test that numpy arrays are read-only (raises ValueError on in-place mutation).
        with self.assertRaises(ValueError):
            wing.Ler_Gs_Cgs[0] = 999.0

        with self.assertRaises(ValueError):
            wing.angles_Gs_to_Wn_ixyz[0] = 45.0

    def test_set_once_attributes_cannot_be_reassigned(self):
        """Test that set-once attributes raise AttributeError on second assignment."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Test that set-once attributes raise AttributeError when set again.
        with self.assertRaises(AttributeError):
            wing.symmetry_type = 2

        with self.assertRaises(AttributeError):
            wing.num_spanwise_panels = 100

        with self.assertRaises(AttributeError):
            wing.num_panels = 1000

        with self.assertRaises(AttributeError):
            wing.panels = np.empty((1, 1), dtype=object)

    def test_wing_cross_sections_tuple_immutability(self):
        """Test that wing_cross_sections tuple cannot be mutated."""
        wing = self.type_1_wing

        # Verify it's a tuple (tuples don't have append method).
        self.assertIsInstance(wing.wing_cross_sections, tuple)

        # Attempting to call append should raise AttributeError.
        with self.assertRaises(AttributeError):
            # noinspection PyUnresolvedReferences
            wing.wing_cross_sections.append(self.root_wcs)

    def test_deepcopy_preserves_geometric_properties(self):
        """Test that deepcopy preserves geometric property calculations."""
        import copy

        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertAlmostEqual(copied.span, original.span, places=10)
        self.assertAlmostEqual(
            copied.projected_area, original.projected_area, places=10
        )
        self.assertAlmostEqual(copied.wetted_area, original.wetted_area, places=10)

    def test_deepcopy_type_4_wing(self):
        """Test that deepcopy works correctly for type 4 symmetric Wings."""
        import copy

        original = self.type_4_wing
        original.generate_mesh(4)

        copied = copy.deepcopy(original)

        self.assertEqual(copied.symmetry_type, 4)
        self.assertEqual(copied.symmetric, True)
        self.assertIsNotNone(copied.panels)
        self.assertAlmostEqual(copied.span, original.span, places=10)

    def test_deepcopy_copied_wing_is_functional(self):
        """Test that copied Wings are fully functional."""
        import copy

        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        span = copied.span
        projected_area = copied.projected_area
        standard_mean_chord = copied.standard_mean_chord
        mean_aerodynamic_chord = copied.mean_aerodynamic_chord

        self.assertIsNotNone(span)
        self.assertIsNotNone(projected_area)
        self.assertIsNotNone(standard_mean_chord)
        self.assertIsNotNone(mean_aerodynamic_chord)
        self.assertGreater(span, 0.0)
        self.assertGreater(projected_area, 0.0)

    def test_deepcopy_meshed_wing_with_populated_axis_caches(self):
        """Test deepcopy of a meshed Wing whose WnX_G, WnY_G, and WnZ_G caches
        have been populated by accessing those properties.

        Accessing the axis-vector properties forces the lazy-computed arrays to
        be cached as non-None, which exercises the not-None copy branches inside
        __deepcopy__ that are skipped when only calling generate_mesh.
        """
        import copy

        import numpy.testing as npt

        original = self.type_1_wing
        original.generate_mesh(1)

        # Access cached properties to populate internal caches before deepcopy.
        _wn_x = original.WnX_G
        _wn_y = original.WnY_G
        _wn_z = original.WnZ_G
        _children_to_wcs = original.children_T_pas_Wn_Ler_to_Wcs_Lp
        _children_from_wcs = original.children_T_pas_Wcs_Lp_to_Wn_Ler

        copied = copy.deepcopy(original)

        # Verify the copied wing has matching axis vectors.
        npt.assert_array_equal(copied.WnX_G, original.WnX_G)
        npt.assert_array_equal(copied.WnY_G, original.WnY_G)
        npt.assert_array_equal(copied.WnZ_G, original.WnZ_G)

        # Verify the children transformation matrices are copied correctly.
        for i in range(len(original.wing_cross_sections)):
            npt.assert_array_equal(
                copied.children_T_pas_Wn_Ler_to_Wcs_Lp[i],
                original.children_T_pas_Wn_Ler_to_Wcs_Lp[i],
            )
            npt.assert_array_equal(
                copied.children_T_pas_Wcs_Lp_to_Wn_Ler[i],
                original.children_T_pas_Wcs_Lp_to_Wn_Ler[i],
            )


class TestWingGetPlottableData(unittest.TestCase):
    """Tests for Wing.get_plottable_data method."""

    def setUp(self):
        """Set up test fixtures for get_plottable_data tests."""

    def test_get_plottable_data_returns_none_when_symmetry_type_not_set(self):
        """Test that get_plottable_data returns None when symmetry_type not set."""
        wing = geometry_fixtures.make_type_1_wing_fixture()

        # Symmetry type not set (Wing not meshed).
        self.assertIsNone(wing.symmetry_type)

        result = wing.get_plottable_data(show=False)

        self.assertIsNone(result)

    def test_get_plottable_data_returns_list_when_meshed(self):
        """Test that get_plottable_data returns a list when meshed."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_returns_list_of_lists(self):
        """Test that get_plottable_data returns two lists of ndarrays."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)

        # First element is list of Airfoil outlines.
        self.assertIsInstance(result[0], list)
        # Second element is list of Airfoil mean camber lines.
        self.assertIsInstance(result[1], list)

    def test_get_plottable_data_returns_ndarrays_for_each_cross_section(self):
        """Test that get_plottable_data returns one ndarray per WingCrossSection."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)

        num_wcs = len(wing.wing_cross_sections)

        # Should have one outline array per WingCrossSection.
        self.assertEqual(len(result[0]), num_wcs)
        # Should have one MCL array per WingCrossSection.
        self.assertEqual(len(result[1]), num_wcs)

        for outline in result[0]:
            self.assertIsInstance(outline, np.ndarray)
            self.assertEqual(outline.shape[1], 3)  # 3D points

        for mcl in result[1]:
            self.assertIsInstance(mcl, np.ndarray)
            self.assertEqual(mcl.shape[1], 3)  # 3D points

    def test_get_plottable_data_three_section_wing(self):
        """Test get_plottable_data for Wing with 3 WingCrossSections."""
        wing = geometry_fixtures.make_three_section_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)

        # Should have 3 outlines and 3 MCLs.
        self.assertEqual(len(result[0]), 3)
        self.assertEqual(len(result[1]), 3)

    def test_get_plottable_data_type_4_symmetric_wing(self):
        """Test get_plottable_data for type 4 symmetric Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        result = wing.get_plottable_data(show=False)

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_default_show_is_false(self):
        """Test that get_plottable_data default for show is False."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Call without show parameter, should return data (not None).
        result = wing.get_plottable_data()

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_accepts_numpy_bool(self):
        """Test that get_plottable_data accepts numpy bool for show parameter."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=np.bool_(False))

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_invalid_show_type_raises(self):
        """Test that get_plottable_data raises error for invalid show type."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            wing.get_plottable_data(show="invalid")

    def test_get_plottable_data_with_panels_meshed(self):
        """Test get_plottable_data with meshed Panels."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Verify Panels are meshed.
        self.assertIsNotNone(wing.panels)

        result = wing.get_plottable_data(show=False)

        self.assertIsNotNone(result)


class TestWingTransformationMatrixCaching(unittest.TestCase):
    """Tests for Wing transformation matrix caching behavior."""

    def test_T_pas_G_Cg_to_Wn_Ler_returns_same_object_on_repeated_access(self):
        """Test that T_pas_G_Cg_to_Wn_Ler returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        T1 = wing.T_pas_G_Cg_to_Wn_Ler
        T2 = wing.T_pas_G_Cg_to_Wn_Ler

        self.assertIs(T1, T2)

    def test_T_pas_Wn_Ler_to_G_Cg_returns_same_object_on_repeated_access(self):
        """Test that T_pas_Wn_Ler_to_G_Cg returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        T1 = wing.T_pas_Wn_Ler_to_G_Cg
        T2 = wing.T_pas_Wn_Ler_to_G_Cg

        self.assertIs(T1, T2)

    def test_WnX_G_returns_same_object_on_repeated_access(self):
        """Test that WnX_G returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        v1 = wing.WnX_G
        v2 = wing.WnX_G

        self.assertIs(v1, v2)

    def test_WnY_G_returns_same_object_on_repeated_access(self):
        """Test that WnY_G returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        v1 = wing.WnY_G
        v2 = wing.WnY_G

        self.assertIs(v1, v2)

    def test_WnZ_G_returns_same_object_on_repeated_access(self):
        """Test that WnZ_G returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        v1 = wing.WnZ_G
        v2 = wing.WnZ_G

        self.assertIs(v1, v2)

    def test_children_T_pas_Wn_Ler_to_Wcs_Lp_returns_same_object_on_repeated_access(
        self,
    ):
        """Test that children_T_pas_Wn_Ler_to_Wcs_Lp returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_Wn_Ler_to_Wcs_Lp
        list2 = wing.children_T_pas_Wn_Ler_to_Wcs_Lp

        self.assertIs(list1, list2)

    def test_children_T_pas_Wcs_Lp_to_Wn_Ler_returns_same_object_on_repeated_access(
        self,
    ):
        """Test that children_T_pas_Wcs_Lp_to_Wn_Ler returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_Wcs_Lp_to_Wn_Ler
        list2 = wing.children_T_pas_Wcs_Lp_to_Wn_Ler

        self.assertIs(list1, list2)

    def test_children_T_pas_G_Cg_to_Wcs_Lp_returns_same_object_on_repeated_access(self):
        """Test that children_T_pas_G_Cg_to_Wcs_Lp returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_G_Cg_to_Wcs_Lp
        list2 = wing.children_T_pas_G_Cg_to_Wcs_Lp

        self.assertIs(list1, list2)

    def test_children_T_pas_Wcs_Lp_to_G_Cg_returns_same_object_on_repeated_access(self):
        """Test that children_T_pas_Wcs_Lp_to_G_Cg returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_Wcs_Lp_to_G_Cg
        list2 = wing.children_T_pas_Wcs_Lp_to_G_Cg

        self.assertIs(list1, list2)

    def test_transformation_matrices_are_read_only(self):
        """Test that transformation matrices are read only."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        with self.assertRaises(ValueError):
            wing.T_pas_G_Cg_to_Wn_Ler[0, 0] = 999.0

        with self.assertRaises(ValueError):
            wing.T_pas_Wn_Ler_to_G_Cg[0, 0] = 999.0

    def test_basis_vectors_are_read_only(self):
        """Test that basis vectors are read only."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        with self.assertRaises(ValueError):
            wing.WnX_G[0] = 999.0

        with self.assertRaises(ValueError):
            wing.WnY_G[0] = 999.0

        with self.assertRaises(ValueError):
            wing.WnZ_G[0] = 999.0

    def test_children_transformation_matrices_are_read_only(self):
        """Test that children transformation matrices are read only."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        for T in wing.children_T_pas_Wn_Ler_to_Wcs_Lp:
            with self.assertRaises(ValueError):
                T[0, 0] = 999.0

        for T in wing.children_T_pas_Wcs_Lp_to_Wn_Ler:
            with self.assertRaises(ValueError):
                T[0, 0] = 999.0

        for T in wing.children_T_pas_G_Cg_to_Wcs_Lp:
            with self.assertRaises(ValueError):
                T[0, 0] = 999.0

        for T in wing.children_T_pas_Wcs_Lp_to_G_Cg:
            with self.assertRaises(ValueError):
                T[0, 0] = 999.0


class TestSingleStepWingMethods(unittest.TestCase):
    """This class contains unit tests for Wing.explode_wing,
    Wing.interpolate_between_wing_cross_sections, and the explode_into_strips
    parameter."""

    def _make_wcs_3span(self):
        """Create a root WCS with num_spanwise_panels=3."""
        return ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
            num_spanwise_panels=3,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing="uniform",
        )

    def _make_tip_wcs(self):
        """Create a tip WCS with num_spanwise_panels=None."""
        return ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
            num_spanwise_panels=None,
            chord=0.5,
            Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing=None,
        )

    def _make_plain_wing(self, explode_into_strips=False):
        """Create a minimal 2-WCS wing."""
        return ps.geometry.wing.Wing(
            wing_cross_sections=[self._make_wcs_3span(), self._make_tip_wcs()],
            name="Test Wing",
            Ler_Gs_Cgs=(0.0, 0.0, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            explode_into_strips=explode_into_strips,
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
        )

    def test_explode_into_strips_false_wcs_count_unchanged(self):
        """Test that explode_into_strips=False keeps the original two WCS."""
        wing = self._make_plain_wing(explode_into_strips=False)
        self.assertEqual(len(wing.wing_cross_sections), 2)

    def test_explode_into_strips_true_correct_wcs_count(self):
        """Test that explode_into_strips=True with root num_spanwise=3 produces 4 WCS
        (root copy plus 3 interpolated including the tip)."""
        wing = self._make_plain_wing(explode_into_strips=True)
        self.assertEqual(len(wing.wing_cross_sections), 4)

    def test_explode_into_strips_true_non_tip_have_num_spanwise_one(self):
        """Test that all non-tip WCS have num_spanwise_panels=1 after explode."""
        wing = self._make_plain_wing(explode_into_strips=True)
        for wcs in wing.wing_cross_sections[:-1]:
            with self.subTest(wcs=wcs):
                self.assertEqual(wcs.num_spanwise_panels, 1)

    def test_explode_into_strips_true_tip_has_none_spanwise(self):
        """Test that the last WCS (tip) has num_spanwise_panels=None after explode."""
        wing = self._make_plain_wing(explode_into_strips=True)
        self.assertIsNone(wing.wing_cross_sections[-1].num_spanwise_panels)

    def test_interpolate_returns_correct_count_with_first_wcs(self):
        """Test that interpolate_between_wing_cross_sections with first_wcs=True
        returns N+1 WCS (root copy plus N interpolated)."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.interpolate_between_wing_cross_sections(
            self._make_wcs_3span(), self._make_tip_wcs(), first_wcs=True
        )
        # N=3: root copy + 3 interpolated = 4 total
        self.assertEqual(len(result), 4)

    def test_interpolate_returns_correct_count_without_first_wcs(self):
        """Test that interpolate_between_wing_cross_sections with first_wcs=False
        returns N WCS (no root copy)."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.interpolate_between_wing_cross_sections(
            self._make_wcs_3span(), self._make_tip_wcs(), first_wcs=False
        )
        # N=3: no root copy => 3 WCS total
        self.assertEqual(len(result), 3)

    def test_interpolate_root_copy_chord(self):
        """Test that the first WCS in the result has the root chord when
        first_wcs=True."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.interpolate_between_wing_cross_sections(
            self._make_wcs_3span(), self._make_tip_wcs(), first_wcs=True
        )
        self.assertAlmostEqual(result[0].chord, 1.0)

    def test_interpolate_tip_chord(self):
        """Test that the last WCS in the result has the tip chord."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.interpolate_between_wing_cross_sections(
            self._make_wcs_3span(), self._make_tip_wcs(), first_wcs=True
        )
        self.assertAlmostEqual(result[-1].chord, 0.5)

    def test_interpolate_intermediate_chord_linearly_interpolated(self):
        """Test that intermediate WCS chords are linearly interpolated between
        root (1.0) and tip (0.5)."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.interpolate_between_wing_cross_sections(
            self._make_wcs_3span(), self._make_tip_wcs(), first_wcs=True
        )
        # N=3: t=1/3, 2/3, 3/3
        # chord at t=1/3: (1 - 1/3)*1.0 + (1/3)*0.5 = 5/6
        expected_first_interp = (1.0 - 1.0 / 3.0) * 1.0 + (1.0 / 3.0) * 0.5
        self.assertAlmostEqual(result[1].chord, expected_first_interp, places=10)

    def test_interpolate_lp_y_divided_by_n(self):
        """Test that the Lp_Wcsp_Lpp y-component of each interpolated WCS is
        tip_Lp_y / N."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.interpolate_between_wing_cross_sections(
            self._make_wcs_3span(), self._make_tip_wcs(), first_wcs=False
        )
        # tip has Lp_y = 0.5; N=3 => each section Lp_y = 0.5/3
        expected_lp_y = 0.5 / 3.0
        for wcs in result:
            with self.subTest(wcs=wcs):
                self.assertAlmostEqual(
                    float(wcs.Lp_Wcsp_Lpp[1]), expected_lp_y, places=10
                )

    def test_explode_wing_with_two_wcs_returns_correct_count(self):
        """Test that explode_wing with a 2-WCS input (root: num=3, tip) returns 4
        WCS."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.explode_wing([self._make_wcs_3span(), self._make_tip_wcs()])
        self.assertEqual(len(result), 4)

    def test_explode_wing_all_non_tip_have_num_spanwise_one(self):
        """Test that explode_wing produces WCS where every non-tip entry has
        num_spanwise_panels=1."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.explode_wing([self._make_wcs_3span(), self._make_tip_wcs()])
        for wcs in result[:-1]:
            with self.subTest(wcs=wcs):
                self.assertEqual(wcs.num_spanwise_panels, 1)

    def test_explode_wing_last_wcs_is_tip(self):
        """Test that explode_wing produces a final WCS with num_spanwise_panels=None."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing.explode_wing([self._make_wcs_3span(), self._make_tip_wcs()])
        self.assertIsNone(result[-1].num_spanwise_panels)

    def test_explode_wing_rejects_non_uniform_spanwise_spacing(self):
        """Test that explode_wing raises ValueError when a non tip WCS uses cosine
        spanwise spacing, since the explosion assumes uniformly distributed
        intermediates."""
        wing = self._make_plain_wing(explode_into_strips=False)
        cosine_root = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
            num_spanwise_panels=3,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing="cosine",
        )
        with self.assertRaises(ValueError):
            wing.explode_wing([cosine_root, self._make_tip_wcs()])


if __name__ == "__main__":
    unittest.main()
