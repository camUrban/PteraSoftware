"""This module contains a class to test Wings."""

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
        self.test_airfoil = geometry_fixtures.make_test_airfoil_fixture()
        self.root_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        self.tip_wcs = geometry_fixtures.make_tip_wing_cross_section_fixture()

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
                self.assertIsInstance(wing.wing_cross_sections, list)
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
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                symmetric=True,
                mirror_only=True,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )

        # Test that symmetry parameters must be None when no symmetry
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=[0.0, 1.0, 0.0],
            )

        # Test that symmetry parameters must be provided when symmetric=True
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
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
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[self.root_wcs, self.tip_wcs],
                    chordwise_spacing=spacing,
                )
                self.assertEqual(wing.chordwise_spacing, spacing)

    def test_wing_with_different_chordwise_panels(self):
        """Test Wing creation with different numbers of chordwise panels."""
        panel_counts = [1, 4, 8, 16, 32]

        for count in panel_counts:
            with self.subTest(count=count):
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[self.root_wcs, self.tip_wcs],
                    num_chordwise_panels=count,
                )
                self.assertEqual(wing.num_chordwise_panels, count)

    def test_wing_parameter_validation(self):
        """Test parameter validation for Wing initialization."""
        # Test invalid Ler position
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs], Ler_Gs_Cgs="invalid"
            )

        # Test invalid angles
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                angles_Gs_to_Wn_ixyz="invalid",
            )

        # Test invalid num_chordwise_panels
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                num_chordwise_panels=0,
            )

        # Test invalid chordwise_spacing
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                chordwise_spacing="invalid_spacing",
            )

    def test_wing_name_validation(self):
        """Test Wing name parameter validation."""
        # Test valid string name
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[self.root_wcs, self.tip_wcs], name="Test Wing Name"
        )
        self.assertEqual(wing.name, "Test Wing Name")

        # Test invalid name type
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs], name=123
            )

    def test_symmetry_normal_normalization(self):
        """Test that symmetry normal vectors are properly normalized."""
        # Create Wing with non-unit normal vector
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[self.root_wcs, self.tip_wcs],
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

    def test_properties_none_before_meshing(self):
        """Test that span, standard_mean_chord, and mean_aerodynamic_chord return
        None before meshing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()

        # Properties should return None before meshing
        self.assertIsNone(wing.span)
        self.assertIsNone(wing.projected_area)
        self.assertIsNone(wing.wetted_area)
        self.assertIsNone(wing.standard_mean_chord)
        self.assertIsNone(wing.mean_aerodynamic_chord)

    def test_properties_available_after_meshing(self):
        """Test that span, standard_mean_chord, and mean_aerodynamic_chord are
        available after meshing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Properties should be available and positive after meshing
        self.assertIsNotNone(wing.span)
        self.assertIsNotNone(wing.projected_area)
        self.assertIsNotNone(wing.wetted_area)
        self.assertIsNotNone(wing.standard_mean_chord)
        self.assertIsNotNone(wing.mean_aerodynamic_chord)

        self.assertGreater(wing.span, 0.0)
        self.assertGreater(wing.projected_area, 0.0)
        self.assertGreater(wing.wetted_area, 0.0)
        self.assertGreater(wing.standard_mean_chord, 0.0)
        self.assertGreater(wing.mean_aerodynamic_chord, 0.0)

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


if __name__ == "__main__":
    unittest.main()
