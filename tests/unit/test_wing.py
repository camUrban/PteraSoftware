"""This module contains a class to test Wings.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    TestWings: This is a class with functions to test Wings.
"""

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

    def tearDown(self):
        """Clean up test fixtures."""
        del self.type_1_wing
        del self.type_2_wing
        del self.type_3_wing
        del self.type_4_wing
        del self.type_5_wing
        del self.test_airfoil
        del self.root_wcs
        del self.tip_wcs

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
                self.assertEqual(len(wing.prelimLer_G_Cg), 3)
                self.assertEqual(len(wing.angles_G_to_prelimWn), 3)
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
                symmetry_normal_Wn=[0.0, 1.0, 0.0],
                symmetry_point_Wn_Ler=[0.0, 0.0, 0.0],
            )

        # Test that symmetry parameters must be None when no symmetry
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                symmetric=False,
                mirror_only=False,
                symmetry_normal_Wn=[0.0, 1.0, 0.0],  # Should be None
            )

        # Test that symmetry parameters must be provided when symmetric=True
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                symmetric=True,
                mirror_only=False,
                symmetry_normal_Wn=None,  # Should be provided
            )

    def test_wing_type_1_properties(self):
        """Test type 1 wing (symmetric=False, mirror_only=False) properties."""
        wing = self.type_1_wing

        # Test basic properties
        self.assertFalse(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        self.assertIsNone(wing.symmetry_normal_Wn)
        self.assertIsNone(wing.symmetry_point_Wn_Ler)

        # Test that symmetry_type is None before meshing
        self.assertIsNone(wing.symmetry_type)

    def test_wing_type_2_properties(self):
        """Test type 2 wing (mirror_only=True, coincident symmetry plane) properties."""
        wing = self.type_2_wing

        # Test basic properties
        self.assertFalse(wing.symmetric)
        self.assertTrue(wing.mirror_only)
        npt.assert_array_equal(wing.symmetry_normal_Wn, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetry_point_Wn_Ler, np.array([0.0, 0.0, 0.0]))

    def test_wing_type_3_properties(self):
        """Test type 3 wing (mirror_only=True, non-coincident symmetry plane) properties."""
        wing = self.type_3_wing

        # Test basic properties
        self.assertFalse(wing.symmetric)
        self.assertTrue(wing.mirror_only)
        # Check that symmetry plane is non-coincident
        self.assertFalse(
            np.allclose(wing.symmetry_normal_Wn, np.array([0.0, 1.0, 0.0]))
        )
        self.assertFalse(
            np.allclose(wing.symmetry_point_Wn_Ler, np.array([0.0, 0.0, 0.0]))
        )

    def test_wing_type_4_properties(self):
        """Test type 4 wing (symmetric=True, coincident symmetry plane) properties."""
        wing = self.type_4_wing

        # Test basic properties
        self.assertTrue(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        npt.assert_array_equal(wing.symmetry_normal_Wn, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetry_point_Wn_Ler, np.array([0.0, 0.0, 0.0]))

        # Test that WingCrossSections have control surface symmetry types
        for wcs in wing.wing_cross_sections:
            self.assertEqual(wcs.control_surface_symmetry_type, "symmetric")

    def test_wing_type_5_properties(self):
        """Test type 5 wing (symmetric=True, non-coincident symmetry plane) properties."""
        wing = self.type_5_wing

        # Test basic properties (before Airplane processing)
        self.assertTrue(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        # Check that symmetry plane is non-coincident
        self.assertFalse(
            np.allclose(wing.symmetry_normal_Wn, np.array([0.0, 1.0, 0.0]))
        )
        self.assertFalse(
            np.allclose(wing.symmetry_point_Wn_Ler, np.array([0.0, 0.0, 0.0]))
        )

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
        # Test invalid preliminary position
        with self.assertRaises(TypeError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                prelimLer_G_Cg="invalid",
            )

        # Test invalid angles
        with self.assertRaises(TypeError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                angles_G_to_prelimWn="invalid",
            )

        # Test invalid num_chordwise_panels
        with self.assertRaises((ValueError, TypeError)):
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                num_chordwise_panels=0,  # Must be positive
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
            ps.geometry.wing.Wing(
                wing_cross_sections=[self.root_wcs, self.tip_wcs],
                name=123,  # Should be string
            )

    def test_symmetry_normal_normalization(self):
        """Test that symmetry normal vectors are properly normalized."""
        # Create Wing with non-unit normal vector
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[self.root_wcs, self.tip_wcs],
            symmetric=False,
            mirror_only=True,
            symmetry_normal_Wn=[0.0, 5.0, 0.0],  # Non-unit vector
            symmetry_point_Wn_Ler=[0.0, 0.0, 0.0],
        )

        # Should be normalized to unit vector
        npt.assert_allclose(np.linalg.norm(wing.symmetry_normal_Wn), 1.0, atol=1e-14)
        npt.assert_allclose(
            wing.symmetry_normal_Wn, np.array([0.0, 1.0, 0.0]), atol=1e-14
        )

    # TODO: Finalize Wing's get_plottable_data testing.
    # def test_wing_get_plottable_data(self):
    #     """Test that the get_plottable_data method works correctly."""
    #     airplane = geometry_fixtures.make_basic_airplane_fixture()
    #     airplane = geometry_fixtures.make_type_5_wing_airplane_fixture()
    #     this_wing = airplane.wings[1]
    #     this_wing.get_plottable_data(show=True)
    #     # airplane.draw()


if __name__ == "__main__":
    unittest.main()
