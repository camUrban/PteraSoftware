"""This module contains classes to test Wings."""

import copy
import unittest
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestWing(unittest.TestCase):
    """This class contains unit tests for the Wing class."""

    def setUp(self) -> None:
        """Set up test fixtures for Wing tests."""
        # Create fixtures for all Wing types.
        self.type_1_wing = geometry_fixtures.make_type_1_wing_fixture()
        self.type_2_wing = geometry_fixtures.make_type_2_wing_fixture()
        self.type_3_wing = geometry_fixtures.make_type_3_wing_fixture()
        self.type_4_wing = geometry_fixtures.make_type_4_wing_fixture()
        self.type_5_wing = geometry_fixtures.make_type_5_wing_fixture()

        # Create additional test fixtures.
        self.root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )

    def test_initialization_valid_parameters(self) -> None:
        """Test Wing initialization with valid parameters for all types."""
        # Test that all Wing types initialize correctly.
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

    def test_wing_cross_sections_validation(self) -> None:
        """Test that wing_cross_sections parameter validation works correctly."""
        # Test empty list raises error.
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(wing_cross_sections=[])

        # Test non-list raises error.
        bad_wing_cross_sections: Any = "not a list"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(wing_cross_sections=bad_wing_cross_sections)

        # Test single WingCrossSection raises error (need at least 2).
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(wing_cross_sections=[self.root_wing_cross_section])

        # Test non-WingCrossSection objects raise error.
        bad_wing_cross_section: Any = "not a wing_cross_section"
        with self.assertRaises(TypeError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    self.root_wing_cross_section,
                    bad_wing_cross_section,
                ]
            )

    def test_symmetry_parameter_validation(self) -> None:
        """Test symmetry parameter validation logic."""
        # Test that symmetric and mirror_only cannot both be True. Create fresh fixtures
        # since WingCrossSections can only be validated once.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                symmetric=True,
                mirror_only=True,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )

        # Test that symmetry parameters must be None when no symmetry.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=[0.0, 1.0, 0.0],
            )

        # Test that symmetry parameters must be provided when symmetric=True.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=None,
            )

    def test_wing_type_1_properties(self) -> None:
        """Test type 1 wing (symmetric=False, mirror_only=False) properties."""
        wing = self.type_1_wing

        # Test basic properties.
        self.assertFalse(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        self.assertIsNone(wing.symmetryNormal_G)
        self.assertIsNone(wing.symmetryPoint_G_Cg)

        # Test that symmetry_type is None before meshing.
        self.assertIsNone(wing.symmetry_type)

    def test_wing_type_2_properties(self) -> None:
        """Test type 2 wing (mirror_only=True, coincident symmetry plane) properties."""
        wing = self.type_2_wing

        # Test basic properties.
        self.assertFalse(wing.symmetric)
        self.assertTrue(wing.mirror_only)
        npt.assert_array_equal(wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([1.0, 0.0, 0.5]))

    def test_wing_type_3_properties(self) -> None:
        """Test type 3 wing (mirror_only=True, non-coincident symmetry plane)
        properties."""
        wing = self.type_3_wing

        # Test basic properties.
        self.assertFalse(wing.symmetric)
        self.assertTrue(wing.mirror_only)
        npt.assert_array_equal(wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([0.0, -0.5, 0.0]))

    def test_wing_type_4_properties(self) -> None:
        """Test type 4 wing (symmetric=True, coincident symmetry plane) properties."""
        wing = self.type_4_wing

        # Test basic properties.
        self.assertTrue(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        npt.assert_array_equal(wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]))
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([1.0, 0.0, 0.5]))

        # Test that WingCrossSections have control surface symmetry types.
        for wing_cross_section in wing.wing_cross_sections:
            self.assertEqual(
                wing_cross_section.control_surface_symmetry_type, "symmetric"
            )

    def test_wing_type_5_properties(self) -> None:
        """Test type 5 wing (symmetric=True, non-coincident symmetry plane)
        properties."""
        wing = self.type_5_wing

        # Test basic properties (before Airplane processing).
        self.assertTrue(wing.symmetric)
        self.assertFalse(wing.mirror_only)
        npt.assert_array_equal(
            wing.symmetryNormal_G, np.array([0.0, np.sqrt(2) / 2, np.sqrt(2) / 2])
        )
        npt.assert_array_equal(wing.symmetryPoint_G_Cg, np.array([0.5, 0.0, 0.0]))

    def test_symmetric_wing_panel_normals_mirror_across_symmetry_plane(self) -> None:
        """Test that a symmetric Wing's panel normals mirror correctly across its xz
        symmetry plane.

        Each panel at +y has a partner at -y whose unit normal is the source normal with
        only its y component negated. This pins the inner-outer corner swap in the
        mirrored-wing meshing, which keeps the recomputed normals correctly oriented;
        reflecting the corner points without that swap would negate the normal's x and z
        components instead.
        """
        wing = geometry_fixtures.make_symmetric_dihedral_wing_fixture()
        wing.generate_mesh(4)
        assert wing.panels is not None
        panels = wing.panels.flatten()

        centroids = np.array(
            [
                (panel.Frpp_G_Cg + panel.Flpp_G_Cg + panel.Brpp_G_Cg + panel.Blpp_G_Cg)
                / 4.0
                for panel in panels
            ]
        )
        normals = np.array([panel.unitNormal_G for panel in panels])

        # The dihedral should give the normals a non-trivial spanwise component, so the
        # y-flip below is a meaningful check rather than an identity.
        self.assertGreater(np.max(np.abs(normals[:, 1])), 0.1)

        # Reflection across the xz (y = 0) symmetry plane negates only the y component.
        y_flip = np.array([1.0, -1.0, 1.0])
        for centroid, normal in zip(centroids, normals):
            distances = np.linalg.norm(centroids - centroid * y_flip, axis=1)
            partner = int(np.argmin(distances))

            # A true mirror partner exists across the y = 0 plane.
            self.assertLess(distances[partner], 1e-9)

            # The mirror partner's normal is this normal with its y component negated.
            npt.assert_allclose(normals[partner], normal * y_flip, atol=1e-10)

    def test_generate_mesh_symmetry_type_1(self) -> None:
        """Test generate_mesh method with type 1 symmetry."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)
        self.assertEqual(wing.symmetry_type, 1)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        assert wing.panels is not None
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_symmetry_type_2(self) -> None:
        """Test generate_mesh method with type 2 symmetry."""
        wing = geometry_fixtures.make_type_2_wing_fixture()
        wing.generate_mesh(2)
        self.assertEqual(wing.symmetry_type, 2)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        assert wing.panels is not None
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_symmetry_type_3(self) -> None:
        """Test generate_mesh method with type 3 symmetry."""
        wing = geometry_fixtures.make_type_3_wing_fixture()
        wing.generate_mesh(3)
        self.assertEqual(wing.symmetry_type, 3)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        assert wing.panels is not None
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_symmetry_type_4(self) -> None:
        """Test generate_mesh method with type 4 symmetry."""
        wing = geometry_fixtures.make_type_4_wing_fixture()
        wing.generate_mesh(4)
        self.assertEqual(wing.symmetry_type, 4)
        self.assertTrue(hasattr(wing, "panels"))
        self.assertIsInstance(wing.panels, np.ndarray)
        assert wing.panels is not None
        self.assertEqual(wing.panels.ndim, 2)

    def test_generate_mesh_invalid_symmetry_type(self) -> None:
        """Test generate_mesh method with invalid symmetry types."""
        wing = self.type_1_wing

        invalid_types: list[Any] = [0, 5, -1, 10, "invalid", None]
        for invalid_type in invalid_types:
            with self.subTest(invalid_type=invalid_type):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    wing.generate_mesh(invalid_type)

    def test_transformation_matrices_before_meshing(self) -> None:
        """Test that transformation matrices are None before meshing."""
        wing = self.type_1_wing

        # All transformation matrices should be None before meshing.
        self.assertIsNone(wing.T_pas_G_Cg_to_Wn_Ler)
        self.assertIsNone(wing.T_pas_Wn_Ler_to_G_Cg)

    def test_transformation_matrices_after_meshing(self) -> None:
        """Test that transformation matrices are available after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Transformation matrices should be available after meshing.
        self.assertIsNotNone(wing.T_pas_G_Cg_to_Wn_Ler)
        self.assertIsNotNone(wing.T_pas_Wn_Ler_to_G_Cg)
        assert wing.T_pas_G_Cg_to_Wn_Ler is not None
        assert wing.T_pas_Wn_Ler_to_G_Cg is not None

        # They should be 4 x 4 matrices.
        self.assertEqual(wing.T_pas_G_Cg_to_Wn_Ler.shape, (4, 4))
        self.assertEqual(wing.T_pas_Wn_Ler_to_G_Cg.shape, (4, 4))

        # They should be inverses of each other.
        identity = wing.T_pas_G_Cg_to_Wn_Ler @ wing.T_pas_Wn_Ler_to_G_Cg
        npt.assert_allclose(identity, np.eye(4), atol=1e-14)

    def test_wing_axes_vectors_after_meshing(self) -> None:
        """Test wing axes vectors after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Wing axes basis vectors should be available.
        self.assertIsNotNone(wing.WnX_G)
        self.assertIsNotNone(wing.WnY_G)
        self.assertIsNotNone(wing.WnZ_G)
        assert wing.WnX_G is not None
        assert wing.WnY_G is not None
        assert wing.WnZ_G is not None

        # They should be 3-element vectors.
        self.assertEqual(len(wing.WnX_G), 3)
        self.assertEqual(len(wing.WnY_G), 3)
        self.assertEqual(len(wing.WnZ_G), 3)

        # They should form an orthonormal basis.
        npt.assert_allclose(np.linalg.norm(wing.WnX_G), 1.0, atol=1e-14)
        npt.assert_allclose(np.linalg.norm(wing.WnY_G), 1.0, atol=1e-14)
        npt.assert_allclose(np.linalg.norm(wing.WnZ_G), 1.0, atol=1e-14)

        # They should be orthogonal.
        npt.assert_allclose(np.dot(wing.WnX_G, wing.WnY_G), 0.0, atol=1e-14)
        npt.assert_allclose(wing.WnY_G @ wing.WnZ_G, 0.0, atol=1e-14)
        npt.assert_allclose(wing.WnZ_G @ wing.WnX_G, 0.0, atol=1e-14)

    def test_wing_cross_section_transformations_after_meshing(self) -> None:
        """Test WingCrossSection transformation matrices after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # It should have transformation lists for each WingCrossSection.
        num_wing_cross_sections = len(wing.wing_cross_sections)
        self.assertEqual(
            len(wing.children_T_pas_Wn_Ler_to_Wcs_Lp), num_wing_cross_sections
        )
        self.assertEqual(
            len(wing.children_T_pas_Wcs_Lp_to_Wn_Ler), num_wing_cross_sections
        )
        self.assertEqual(
            len(wing.children_T_pas_G_Cg_to_Wcs_Lp), num_wing_cross_sections
        )
        self.assertEqual(
            len(wing.children_T_pas_Wcs_Lp_to_G_Cg), num_wing_cross_sections
        )

        # Each transformation should be a 4 x 4 matrix.
        for i in range(num_wing_cross_sections):
            self.assertEqual(wing.children_T_pas_Wn_Ler_to_Wcs_Lp[i].shape, (4, 4))
            self.assertEqual(wing.children_T_pas_Wcs_Lp_to_Wn_Ler[i].shape, (4, 4))
            self.assertEqual(wing.children_T_pas_G_Cg_to_Wcs_Lp[i].shape, (4, 4))
            self.assertEqual(wing.children_T_pas_Wcs_Lp_to_G_Cg[i].shape, (4, 4))

    def test_geometric_properties_after_meshing(self) -> None:
        """Test geometric property calculations after meshing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Test that geometric properties are available and positive.
        assert wing.projected_area is not None
        assert wing.wetted_area is not None
        assert wing.span is not None
        assert wing.standard_mean_chord is not None
        assert wing.mean_aerodynamic_chord is not None
        self.assertGreater(wing.projected_area, 0.0)
        self.assertGreater(wing.wetted_area, 0.0)
        self.assertGreater(wing.span, 0.0)
        self.assertGreater(wing.standard_mean_chord, 0.0)
        self.assertGreater(wing.mean_aerodynamic_chord, 0.0)

        # Test that wetted area is greater than projected area (both sides).
        self.assertGreaterEqual(wing.wetted_area, wing.projected_area)

    def test_geometric_properties_before_meshing_return_none(self) -> None:
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

    def test_wing_with_different_chordwise_spacing(self) -> None:
        """Test Wing creation with different chordwise spacing options."""
        spacing_options = ["uniform", "cosine"]

        for spacing in spacing_options:
            with self.subTest(spacing=spacing):
                # Create fresh fixtures for each iteration since WingCrossSections can
                # only be validated once.
                root_wing_cross_section = (
                    geometry_fixtures.make_root_wing_cross_section_fixture()
                )
                tip_wing_cross_section = (
                    geometry_fixtures.make_tip_wing_cross_section_fixture()
                )
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        root_wing_cross_section,
                        tip_wing_cross_section,
                    ],
                    chordwise_spacing=spacing,
                )
                self.assertEqual(wing.chordwise_spacing, spacing)

    def test_wing_with_different_chordwise_panels(self) -> None:
        """Test Wing creation with different numbers of chordwise panels."""
        panel_counts = [1, 4, 8, 16, 32]

        for count in panel_counts:
            with self.subTest(count=count):
                # Create fresh fixtures for each iteration since WingCrossSections can
                # only be validated once.
                root_wing_cross_section = (
                    geometry_fixtures.make_root_wing_cross_section_fixture()
                )
                tip_wing_cross_section = (
                    geometry_fixtures.make_tip_wing_cross_section_fixture()
                )
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        root_wing_cross_section,
                        tip_wing_cross_section,
                    ],
                    num_chordwise_panels=count,
                )
                self.assertEqual(wing.num_chordwise_panels, count)

    def test_wing_parameter_validation(self) -> None:
        """Test parameter validation for Wing initialization."""
        # Test invalid Ler position. Create fresh fixtures since WingCrossSections can
        # only be validated once.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        bad_Ler_Gs_Cgs: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                Ler_Gs_Cgs=bad_Ler_Gs_Cgs,
            )

        # Test invalid angles.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        bad_angles_Gs_to_Wn_ixyz: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                angles_Gs_to_Wn_ixyz=bad_angles_Gs_to_Wn_ixyz,
            )

        # Test invalid num_chordwise_panels.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                num_chordwise_panels=0,
            )

        # Test invalid chordwise_spacing.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                chordwise_spacing="invalid_spacing",
            )

        # Test invalid symmetric.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        bad_symmetric: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                symmetric=bad_symmetric,
            )

        # Test invalid mirror_only.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        bad_mirror_only: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                mirror_only=bad_mirror_only,
            )

        # Test invalid explode_into_strips.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        bad_explode_into_strips: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                explode_into_strips=bad_explode_into_strips,
            )

    def test_wing_name_validation(self) -> None:
        """Test Wing name parameter validation."""
        # Test valid string name. Create fresh fixtures since WingCrossSections can only
        # be validated once.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
            name="Test Wing Name",
        )
        self.assertEqual(wing.name, "Test Wing Name")

        # Test invalid name type.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        bad_name: Any = 123
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                name=bad_name,
            )

    def test_symmetry_normal_normalization(self) -> None:
        """Test that symmetry normal vectors are properly normalized."""
        # Create fresh fixtures since WingCrossSections can only be validated once.
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        # Create Wing with non-unit normal vector.
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
            symmetric=False,
            mirror_only=True,
            symmetryNormal_G=[0.0, 5.0, 0.0],
            symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
        )

        # It should be normalized to unit vector.
        assert wing.symmetryNormal_G is not None
        npt.assert_allclose(np.linalg.norm(wing.symmetryNormal_G), 1.0, atol=1e-14)
        npt.assert_allclose(
            wing.symmetryNormal_G, np.array([0.0, 1.0, 0.0]), atol=1e-14
        )

    def test_three_section_wing_validation(self) -> None:
        """Test Wing with 3 WingCrossSections validates correctly."""
        # Test that valid 3-WingCrossSection Wing initializes correctly.
        wing = geometry_fixtures.make_three_section_wing_fixture()
        self.assertIsInstance(wing, ps.geometry.wing.Wing)
        self.assertEqual(len(wing.wing_cross_sections), 3)

        # Verify all WingCrossSections are validated.
        for wing_cross_section in wing.wing_cross_sections:
            self.assertTrue(wing_cross_section.validated)

    def test_four_section_wing_validation(self) -> None:
        """Test Wing with 4 WingCrossSections validates correctly."""
        # Test that valid 4-WingCrossSection Wing initializes correctly.
        wing = geometry_fixtures.make_four_section_wing_fixture()
        self.assertIsInstance(wing, ps.geometry.wing.Wing)
        self.assertEqual(len(wing.wing_cross_sections), 4)

        # Verify all WingCrossSections are validated.
        for wing_cross_section in wing.wing_cross_sections:
            self.assertTrue(wing_cross_section.validated)

    def test_invalid_middle_wing_cross_section_raises_error(self) -> None:
        """Test that Wing with invalid middle WingCrossSection raises ValueError."""
        # Test that Wing with middle WingCrossSection having num_spanwise_panels=None
        # fails.
        with self.assertRaises(ValueError):
            geometry_fixtures.make_invalid_three_section_wing_fixture()

    def test_invalid_root_wing_cross_section_raises_error(self) -> None:
        """Test that Wing with invalid root WingCrossSection raises ValueError."""
        # Test that Wing with root WingCrossSection having num_spanwise_panels=None
        # fails.
        with self.assertRaises(ValueError):
            geometry_fixtures.make_invalid_root_wing_fixture()

    def test_symmetry_point_none_when_symmetric_raises_value_error(self) -> None:
        """Test that symmetryPoint_G_Cg=None with symmetric=True raises ValueError."""
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=None,
            )

    def test_symmetry_point_not_none_when_no_symmetry_raises_value_error(self) -> None:
        """Test that symmetryPoint_G_Cg not None with no symmetry raises ValueError."""
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        tip_wing_cross_section = geometry_fixtures.make_tip_wing_cross_section_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing(
                wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=None,
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )

    def test_span_simple_rectangular_wing(self) -> None:
        """Test span calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (from y = 0.0 to y = 2.0).
        expected_span = 2.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        assert actual_span is not None
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_simple_tapered_wing(self) -> None:
        """Test span calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 3.0 meters (from y = 0.0 to y = 3.0).
        expected_span = 3.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        assert actual_span is not None
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_symmetric_continuous_rectangular_wing(self) -> None:
        """Test span calculation for symmetric continuous rectangular Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # Expected span: 5.0 meters (2.0 * 2.5, due to symmetry).
        expected_span = 5.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        assert actual_span is not None
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_three_section_tapered_wing(self) -> None:
        """Test span calculation for three section tapered Wing."""
        wing = geometry_fixtures.make_three_section_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 4.0 meters (from y = 0.0 to y = 4.0).
        expected_span = 4.0
        actual_span = wing.span

        self.assertIsNotNone(actual_span)
        assert actual_span is not None
        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_projected_area_simple_rectangular_wing(self) -> None:
        """Test projected area calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Expected projected area: 2.0 square meters (1.0 m chord * 2.0 m span).
        expected_area = 2.0
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        assert actual_area is not None
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_projected_area_simple_tapered_wing(self) -> None:
        """Test projected area calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected projected area: 4.5 square meters. Trapezoid: (chord_root +
        # chord_tip) / 2.0 * span = (2.0 + 1.0) / 2.0 * 3.0 = 4.5.
        expected_area = 4.5
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        assert actual_area is not None
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_projected_area_symmetric_continuous_rectangular_wing(self) -> None:
        """Test projected area calculation for symmetric continuous rectangular
        Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # Expected projected area: 7.5 square meters. Rectangle: chord * span = 1.5 *
        # 5.0 = 7.5.
        expected_area = 7.5
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        assert actual_area is not None
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_projected_area_three_section_tapered_wing(self) -> None:
        """Test projected area calculation for three section tapered Wing."""
        wing = geometry_fixtures.make_three_section_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Expected projected area: 8.0 square meters. Section 1 (root to middle): (3.0 +
        # 2.0) / 2.0 * 2.0 = 5.0. Section 2 (middle to tip): (2.0 + 1.0) / 2.0 * 2.0 =
        # 3.0. Total: 5.0 + 3.0 = 8.0.
        expected_area = 8.0
        actual_area = wing.projected_area

        self.assertIsNotNone(actual_area)
        assert actual_area is not None
        npt.assert_allclose(actual_area, expected_area, rtol=1e-10, atol=1e-14)

    def test_wetted_area_greater_than_projected_area(self) -> None:
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
                assert projected_area is not None
                assert wetted_area is not None
                self.assertGreaterEqual(wetted_area, projected_area)

    def test_standard_mean_chord_simple_rectangular_wing(self) -> None:
        """Test standard mean chord calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Standard mean chord = projected_area / span = 2.0 / 2.0 = 1.0.
        expected_smc = 1.0
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        assert actual_smc is not None
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_simple_tapered_wing(self) -> None:
        """Test standard mean chord calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Standard mean chord = projected_area / span = 4.5 / 3.0 = 1.5.
        expected_smc = 1.5
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        assert actual_smc is not None
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_symmetric_continuous_rectangular_wing(self) -> None:
        """Test standard mean chord calculation for symmetric continuous rectangular
        Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # Standard mean chord = projected_area / span = 7.5 / 5.0 = 1.5.
        expected_smc = 1.5
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        assert actual_smc is not None
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_three_section_tapered_wing(self) -> None:
        """Test standard mean chord calculation for three section tapered Wing."""
        wing = geometry_fixtures.make_three_section_tapered_wing_fixture()
        wing.generate_mesh(1)

        # Standard mean chord = projected_area / span = 8.0 / 4.0 = 2.0.
        expected_smc = 2.0
        actual_smc = wing.standard_mean_chord

        self.assertIsNotNone(actual_smc)
        assert actual_smc is not None
        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_mean_aerodynamic_chord_simple_rectangular_wing(self) -> None:
        """Test mean aerodynamic chord calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # For a rectangular wing (constant chord), MAC = chord = 1.0.
        expected_mac = 1.0
        actual_mac = wing.mean_aerodynamic_chord

        self.assertIsNotNone(actual_mac)
        assert actual_mac is not None
        npt.assert_allclose(actual_mac, expected_mac, rtol=1e-10, atol=1e-14)

    def test_mean_aerodynamic_chord_simple_tapered_wing(self) -> None:
        """Test mean aerodynamic chord calculation for simple tapered Wing."""
        wing = geometry_fixtures.make_simple_tapered_wing_fixture()
        wing.generate_mesh(1)

        # For a linearly tapered wing:
        # MAC = (2.0 / 3.0) * (c_root + c_tip - c_root * c_tip / (c_root + c_tip))
        # With c_root = 2.0, c_tip = 1.0:
        # MAC = (2.0 / 3.0) * (2.0 + 1.0 - 2.0 * 1.0 / (2.0 + 1.0))
        # MAC = (2.0 / 3.0) * (3.0 - 2.0 / 3.0) = (2.0 / 3.0) * (7.0 / 3.0) = 14.0 / 9.0 = 1.555...
        c_root = 2.0
        c_tip = 1.0
        expected_mac = (2.0 / 3.0) * (
            c_root + c_tip - c_root * c_tip / (c_root + c_tip)
        )

        actual_mac = wing.mean_aerodynamic_chord

        self.assertIsNotNone(actual_mac)
        assert actual_mac is not None
        npt.assert_allclose(actual_mac, expected_mac, rtol=1e-10, atol=1e-14)

    def test_mean_aerodynamic_chord_symmetric_continuous_rectangular_wing(self) -> None:
        """Test mean aerodynamic chord calculation for symmetric continuous
        rectangular Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        # For a rectangular wing (constant chord), MAC = chord = 1.5.
        expected_mac = 1.5
        actual_mac = wing.mean_aerodynamic_chord

        self.assertIsNotNone(actual_mac)
        assert actual_mac is not None
        npt.assert_allclose(actual_mac, expected_mac, rtol=1e-10, atol=1e-14)

    def test_geometric_properties_consistency(self) -> None:
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

                # Get properties.
                span = wing.span
                projected_area = wing.projected_area
                standard_mean_chord = wing.standard_mean_chord

                # Verify consistency: projected_area = span * standard_mean_chord.
                self.assertIsNotNone(span)
                self.assertIsNotNone(projected_area)
                self.assertIsNotNone(standard_mean_chord)
                assert span is not None
                assert projected_area is not None
                assert standard_mean_chord is not None

                calculated_area = span * standard_mean_chord
                npt.assert_allclose(
                    projected_area, calculated_area, rtol=1e-10, atol=1e-14
                )

    def test_span_rotated_wing_x_axis(self) -> None:
        """Test span calculation invariance for Wing rotated about x axis."""
        # Create a Wing rotated 45.0 degrees about x axis.
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([45.0, 0.0, 0.0])
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotation about x axis does not affect y extent).
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_rotated_wing_y_axis(self) -> None:
        """Test span calculation invariance for Wing rotated about y axis."""
        # Create a Wing rotated 30.0 degrees about y axis.
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([0.0, 30.0, 0.0])
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotation about y axis does not affect y extent).
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_rotated_wing_z_axis(self) -> None:
        """Test span calculation invariance for Wing rotated about z axis."""
        # Create a Wing rotated 60.0 degrees about z axis.
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([0.0, 0.0, 60.0])
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotation about z axis does not affect y extent).
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_rotated_wing_combined_rotations(self) -> None:
        """Test span calculation invariance for Wing with combined rotations."""
        # Create a Wing with combined rotations.
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture(
            [15.0, 25.0, 35.0]
        )
        wing.generate_mesh(1)

        # Expected span: 2.0 meters (rotations do not affect y extent in wing axes).
        expected_span = 2.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_wing_with_rotated_cross_sections(self) -> None:
        """Test span calculation for Wing with rotated WingCrossSections."""
        wing = geometry_fixtures.make_wing_with_rotated_cross_sections_fixture()
        wing.generate_mesh(1)

        # Expected span: 5.0 meters (rotation about y axis does not affect y position).
        expected_span = 5.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_swept_wing(self) -> None:
        """Test span calculation for swept Wing."""
        wing = geometry_fixtures.make_swept_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 3.0 meters (sweep in x direction does not affect y extent).
        expected_span = 3.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_span_dihedral_wing(self) -> None:
        """Test span calculation for Wing with dihedral."""
        wing = geometry_fixtures.make_dihedral_wing_fixture()
        wing.generate_mesh(1)

        # Expected span: 3.0 meters (dihedral in z direction does not affect y extent in
        # wing axes).
        expected_span = 3.0

        actual_span = wing.span
        self.assertIsNotNone(actual_span)
        assert actual_span is not None

        npt.assert_allclose(actual_span, expected_span, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_rotated_wing(self) -> None:
        """Test standard mean chord calculation for rotated Wing."""
        # Create a Wing rotated 45.0 degrees about x axis.
        wing = geometry_fixtures.make_rotated_rectangular_wing_fixture([45.0, 0.0, 0.0])
        wing.generate_mesh(1)

        # Expected standard mean chord: 1.0 (projected_area / span = 2.0 / 2.0).
        expected_smc = 1.0

        actual_smc = wing.standard_mean_chord
        self.assertIsNotNone(actual_smc)
        assert actual_smc is not None

        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_standard_mean_chord_swept_wing(self) -> None:
        """Test standard mean chord calculation for swept Wing."""
        wing = geometry_fixtures.make_swept_wing_fixture()
        wing.generate_mesh(1)

        # Expected standard mean chord: projected_area / span = 4.5 / 3.0 = 1.5.
        expected_smc = 1.5

        actual_smc = wing.standard_mean_chord
        self.assertIsNotNone(actual_smc)
        assert actual_smc is not None

        npt.assert_allclose(actual_smc, expected_smc, rtol=1e-10, atol=1e-14)

    def test_average_panel_aspect_ratio_returns_none_before_meshing(self) -> None:
        """Test that average_panel_aspect_ratio returns None before meshing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()

        self.assertIsNone(wing.average_panel_aspect_ratio)

    def test_average_panel_aspect_ratio_returns_positive_after_meshing(self) -> None:
        """Test that average_panel_aspect_ratio returns a positive value after meshing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        average_aspect_ratio = wing.average_panel_aspect_ratio
        self.assertIsNotNone(average_aspect_ratio)
        assert average_aspect_ratio is not None
        self.assertGreater(average_aspect_ratio, 0.0)

    def test_average_panel_aspect_ratio_simple_rectangular_wing(self) -> None:
        """Test average_panel_aspect_ratio calculation for simple rectangular Wing."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        average_aspect_ratio = wing.average_panel_aspect_ratio
        self.assertIsNotNone(average_aspect_ratio)
        assert average_aspect_ratio is not None

        # For a rectangular wing with uniform spacing: Panel chord = wing_chord /
        # num_chordwise_panels = 1.0 / 4 = 0.25. Panel span = wing_span /
        # num_spanwise_panels = 2.0 / 8 = 0.25. The expected aspect ratio is
        # approximately 1.0 for roughly square Panels. This is an approximate check
        # since actual Panel aspect ratios depend on the meshing algorithm's
        # implementation details.
        self.assertGreater(average_aspect_ratio, 0.0)
        self.assertLess(average_aspect_ratio, 100.0)

    def test_average_panel_aspect_ratio_type_4_symmetric_wing(self) -> None:
        """Test average_panel_aspect_ratio for type 4 symmetric Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        average_aspect_ratio = wing.average_panel_aspect_ratio
        self.assertIsNotNone(average_aspect_ratio)
        assert average_aspect_ratio is not None
        self.assertGreater(average_aspect_ratio, 0.0)

    def test_average_panel_aspect_ratio_caching(self) -> None:
        """Test that average_panel_aspect_ratio is cached after first access."""
        wing = geometry_fixtures.make_simple_rectangular_wing_fixture()
        wing.generate_mesh(1)

        # Access twice.
        first_access = wing.average_panel_aspect_ratio
        second_access = wing.average_panel_aspect_ratio

        # Values should be identical.
        self.assertEqual(first_access, second_access)

    def test_angles_Gs_to_Wn_ixyz_validation_boundary_values(self) -> None:
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
                root_wing_cross_section = (
                    geometry_fixtures.make_root_wing_cross_section_fixture()
                )
                tip_wing_cross_section = (
                    geometry_fixtures.make_tip_wing_cross_section_fixture()
                )
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        root_wing_cross_section,
                        tip_wing_cross_section,
                    ],
                    angles_Gs_to_Wn_ixyz=angles,
                )
                npt.assert_array_equal(wing.angles_Gs_to_Wn_ixyz, np.array(angles))

    def test_angles_Gs_to_Wn_ixyz_validation_outside_range(self) -> None:
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
                root_wing_cross_section = (
                    geometry_fixtures.make_root_wing_cross_section_fixture()
                )
                tip_wing_cross_section = (
                    geometry_fixtures.make_tip_wing_cross_section_fixture()
                )
                with self.assertRaises(ValueError):
                    ps.geometry.wing.Wing(
                        wing_cross_sections=[
                            root_wing_cross_section,
                            tip_wing_cross_section,
                        ],
                        angles_Gs_to_Wn_ixyz=angles,
                    )

    def test_Ler_Gs_Cgs_accepts_various_input_types(self) -> None:
        """Test that Ler_Gs_Cgs accepts various array-like input types."""
        input_formats: list[np.ndarray | Sequence[float | int]] = [
            np.array([1.0, 2.0, 3.0]),  # ndarray
            [1.0, 2.0, 3.0],  # list
            (1.0, 2.0, 3.0),  # tuple
            [1, 2, 3],  # list of ints
            (1, 2, 3),  # tuple of ints
        ]

        for input_val in input_formats:
            with self.subTest(input_format=type(input_val).__name__):
                root_wing_cross_section = (
                    geometry_fixtures.make_root_wing_cross_section_fixture()
                )
                tip_wing_cross_section = (
                    geometry_fixtures.make_tip_wing_cross_section_fixture()
                )
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        root_wing_cross_section,
                        tip_wing_cross_section,
                    ],
                    Ler_Gs_Cgs=input_val,
                )
                expected = np.array([1.0, 2.0, 3.0])
                npt.assert_array_equal(wing.Ler_Gs_Cgs, expected)

    def test_angles_Gs_to_Wn_ixyz_accepts_various_input_types(self) -> None:
        """Test that angles_Gs_to_Wn_ixyz accepts various array-like input types."""
        input_formats: list[np.ndarray | Sequence[float | int]] = [
            np.array([10.0, 20.0, 30.0]),  # ndarray
            [10.0, 20.0, 30.0],  # list
            (10.0, 20.0, 30.0),  # tuple
            [10, 20, 30],  # list of ints
            (10, 20, 30),  # tuple of ints
        ]

        for input_val in input_formats:
            with self.subTest(input_format=type(input_val).__name__):
                root_wing_cross_section = (
                    geometry_fixtures.make_root_wing_cross_section_fixture()
                )
                tip_wing_cross_section = (
                    geometry_fixtures.make_tip_wing_cross_section_fixture()
                )
                wing = ps.geometry.wing.Wing(
                    wing_cross_sections=[
                        root_wing_cross_section,
                        tip_wing_cross_section,
                    ],
                    angles_Gs_to_Wn_ixyz=input_val,
                )
                expected = np.array([10.0, 20.0, 30.0])
                npt.assert_array_equal(wing.angles_Gs_to_Wn_ixyz, expected)


class TestWingDeepCopy(unittest.TestCase):
    """Tests for Wing.__deepcopy__ method."""

    def setUp(self) -> None:
        """Set up test fixtures for deepcopy tests."""
        self.type_1_wing = geometry_fixtures.make_type_1_wing_fixture()
        self.type_4_wing = geometry_fixtures.make_type_4_wing_fixture()
        self.root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )

    def test_deepcopy_creates_new_instance(self) -> None:
        """Test that deepcopy creates a new Wing instance."""
        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.geometry.wing.Wing)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_wing_parameters(self) -> None:
        """Test that deepcopy preserves Wing parameters."""
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

    def test_deepcopy_creates_independent_arrays(self) -> None:
        """Test that deepcopy creates independent copies of numpy arrays."""
        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.Ler_Gs_Cgs, original.Ler_Gs_Cgs)
        self.assertIsNot(copied.angles_Gs_to_Wn_ixyz, original.angles_Gs_to_Wn_ixyz)

    def test_deepcopy_creates_independent_wing_cross_sections(self) -> None:
        """Test that deepcopy creates independent WingCrossSection copies."""
        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertEqual(
            len(copied.wing_cross_sections), len(original.wing_cross_sections)
        )
        for orig_wing_cross_section, copied_wing_cross_section in zip(
            original.wing_cross_sections, copied.wing_cross_sections
        ):
            self.assertIsNot(orig_wing_cross_section, copied_wing_cross_section)
            self.assertEqual(
                copied_wing_cross_section.chord, orig_wing_cross_section.chord
            )

    def test_deepcopy_preserves_symmetry_attributes(self) -> None:
        """Test that deepcopy preserves symmetry attributes correctly."""
        original = self.type_4_wing
        copied = copy.deepcopy(original)

        self.assertEqual(copied.symmetric, original.symmetric)
        self.assertEqual(copied.mirror_only, original.mirror_only)
        npt.assert_array_equal(copied.symmetryNormal_G, original.symmetryNormal_G)
        npt.assert_array_equal(copied.symmetryPoint_G_Cg, original.symmetryPoint_G_Cg)
        self.assertIsNot(copied.symmetryNormal_G, original.symmetryNormal_G)
        self.assertIsNot(copied.symmetryPoint_G_Cg, original.symmetryPoint_G_Cg)

    def test_deepcopy_preserves_none_symmetry_attributes(self) -> None:
        """Test that deepcopy handles None symmetry attributes correctly."""
        original = self.type_1_wing
        copied = copy.deepcopy(original)

        self.assertIsNone(copied.symmetryNormal_G)
        self.assertIsNone(copied.symmetryPoint_G_Cg)

    def test_deepcopy_unmeshed_wing(self) -> None:
        """Test that deepcopy handles unmeshed Wings correctly."""
        original = self.type_1_wing
        self.assertIsNone(original.panels)

        copied = copy.deepcopy(original)

        self.assertIsNone(copied.symmetry_type)
        self.assertIsNone(copied.num_spanwise_panels)
        self.assertIsNone(copied.num_panels)
        self.assertIsNone(copied.panels)
        self.assertIsNone(copied.gridWrvp_GP1_CgP1)

    def test_deepcopy_meshed_wing_preserves_mesh_metadata(self) -> None:
        """Test that deepcopy preserves mesh metadata for meshed Wings."""
        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertEqual(copied.symmetry_type, original.symmetry_type)
        self.assertEqual(copied.num_spanwise_panels, original.num_spanwise_panels)
        self.assertEqual(copied.num_panels, original.num_panels)

    def test_deepcopy_meshed_wing_preserves_panels(self) -> None:
        """Test that deepcopy preserves Panels for meshed Wings."""
        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertIsNotNone(copied.panels)
        assert copied.panels is not None
        assert original.panels is not None
        self.assertEqual(copied.panels.shape, original.panels.shape)

        for i in range(original.panels.shape[0]):
            for j in range(original.panels.shape[1]):
                orig_panel = original.panels[i, j]
                copied_panel = copied.panels[i, j]
                self.assertIsNot(orig_panel, copied_panel)
                npt.assert_array_equal(copied_panel.Frpp_G_Cg, orig_panel.Frpp_G_Cg)

    def test_deepcopy_resets_wake_state(self) -> None:
        """Test that deepcopy resets wake state to an empty array."""
        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        self.assertIsNotNone(copied.gridWrvp_GP1_CgP1)
        assert copied.gridWrvp_GP1_CgP1 is not None
        assert original.num_spanwise_panels is not None
        self.assertEqual(copied.gridWrvp_GP1_CgP1.shape[0], 0)
        self.assertEqual(
            copied.gridWrvp_GP1_CgP1.shape[1], original.num_spanwise_panels + 1
        )

    def test_deepcopy_independence_modifying_copy_mutable_attrs(self) -> None:
        """Test that modifying mutable attributes on the copy does not affect the
        original."""
        original = self.type_4_wing
        original.generate_mesh(4)
        original_symmetric = original.symmetric
        assert original.symmetryNormal_G is not None
        original_symmetryNormal = original.symmetryNormal_G.copy()

        copied = copy.deepcopy(original)

        # Modify mutable attributes on the copy.
        copied.symmetric = not original_symmetric
        assert copied.symmetryNormal_G is not None
        copied.symmetryNormal_G[0] = 999.0

        # Verify the original is unchanged.
        self.assertEqual(original.symmetric, original_symmetric)
        npt.assert_array_equal(original.symmetryNormal_G, original_symmetryNormal)

    def test_deepcopy_independence_modifying_original_mutable_attrs(self) -> None:
        """Test that modifying mutable attributes on the original does not affect the
        copy."""
        original = self.type_4_wing
        original.generate_mesh(4)

        copied = copy.deepcopy(original)
        copied_symmetric = copied.symmetric
        assert copied.symmetryNormal_G is not None
        copied_symmetryNormal = copied.symmetryNormal_G.copy()

        # Modify mutable attributes on the original.
        original.symmetric = not copied_symmetric
        assert original.symmetryNormal_G is not None
        original.symmetryNormal_G[0] = 999.0

        # Verify the copy is unchanged.
        self.assertEqual(copied.symmetric, copied_symmetric)
        npt.assert_array_equal(copied.symmetryNormal_G, copied_symmetryNormal)

    def test_immutable_attributes_are_read_only(self) -> None:
        """Test that immutable Wing attributes cannot be modified."""
        wing = self.type_1_wing
        wing.generate_mesh(1)

        # Test that setting immutable attributes raises AttributeError.
        with self.assertRaises(AttributeError):
            setattr(wing, "name", "New Name")

        with self.assertRaises(AttributeError):
            setattr(wing, "num_chordwise_panels", 16)

        with self.assertRaises(AttributeError):
            setattr(wing, "chordwise_spacing", "uniform")

        # Test that numpy arrays are read-only (raises ValueError on in-place mutation).
        with self.assertRaises(ValueError):
            wing.Ler_Gs_Cgs[0] = 999.0

        with self.assertRaises(ValueError):
            wing.angles_Gs_to_Wn_ixyz[0] = 45.0

    def test_set_once_attributes_cannot_be_reassigned(self) -> None:
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

    def test_wing_cross_sections_tuple_immutability(self) -> None:
        """Test that wing_cross_sections tuple cannot be mutated."""
        wing = self.type_1_wing

        # Verify it's a tuple (tuples don't have append method).
        self.assertIsInstance(wing.wing_cross_sections, tuple)

        # Attempting to call append should raise AttributeError.
        wing_cross_sections: Any = wing.wing_cross_sections
        with self.assertRaises(AttributeError):
            # noinspection PyUnresolvedReferences
            wing_cross_sections.append(self.root_wing_cross_section)

    def test_deepcopy_preserves_geometric_properties(self) -> None:
        """Test that deepcopy preserves geometric property calculations."""
        original = self.type_1_wing
        original.generate_mesh(1)

        copied = copy.deepcopy(original)

        assert copied.span is not None
        assert original.span is not None
        assert copied.projected_area is not None
        assert original.projected_area is not None
        assert copied.wetted_area is not None
        assert original.wetted_area is not None
        self.assertAlmostEqual(copied.span, original.span, places=10)
        self.assertAlmostEqual(
            copied.projected_area, original.projected_area, places=10
        )
        self.assertAlmostEqual(copied.wetted_area, original.wetted_area, places=10)

    def test_deepcopy_type_4_wing(self) -> None:
        """Test that deepcopy works correctly for type 4 symmetric Wings."""
        original = self.type_4_wing
        original.generate_mesh(4)

        copied = copy.deepcopy(original)

        self.assertEqual(copied.symmetry_type, 4)
        self.assertEqual(copied.symmetric, True)
        self.assertIsNotNone(copied.panels)
        assert copied.span is not None
        assert original.span is not None
        self.assertAlmostEqual(copied.span, original.span, places=10)

    def test_deepcopy_copied_wing_is_functional(self) -> None:
        """Test that copied Wings are fully functional."""
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
        assert span is not None
        assert projected_area is not None
        self.assertGreater(span, 0.0)
        self.assertGreater(projected_area, 0.0)

    def test_deepcopy_meshed_wing_with_populated_axis_caches(self) -> None:
        """Test deepcopy of a meshed Wing whose WnX_G, WnY_G, and WnZ_G caches
        have been populated by accessing those properties.

        Accessing the axis-vector properties forces the lazy-computed arrays to
        be cached as non-None, which exercises the not-None copy branches inside
        __deepcopy__ that are skipped when only calling generate_mesh.
        """
        original = self.type_1_wing
        original.generate_mesh(1)

        # Access cached properties to populate internal caches before deepcopy.
        _WnX_G = original.WnX_G
        _WnY_G = original.WnY_G
        _WnZ_G = original.WnZ_G
        _children_to_wing_cross_section = original.children_T_pas_Wn_Ler_to_Wcs_Lp
        _children_from_wing_cross_section = original.children_T_pas_Wcs_Lp_to_Wn_Ler

        copied = copy.deepcopy(original)

        # Verify the copied Wing has matching axis vectors.
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

    def setUp(self) -> None:
        """Set up test fixtures for get_plottable_data tests."""

    def test_get_plottable_data_returns_none_when_symmetry_type_not_set(self) -> None:
        """Test that get_plottable_data returns None when symmetry_type not set."""
        wing = geometry_fixtures.make_type_1_wing_fixture()

        # Symmetry type not set (Wing not meshed).
        self.assertIsNone(wing.symmetry_type)

        result = wing.get_plottable_data(show=False)

        self.assertIsNone(result)

    def test_get_plottable_data_returns_list_when_meshed(self) -> None:
        """Test that get_plottable_data returns a list when meshed."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)
        assert result is not None

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_returns_list_of_lists(self) -> None:
        """Test that get_plottable_data returns two lists of ndarrays."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)
        assert result is not None

        # First element is list of Airfoil outlines.
        self.assertIsInstance(result[0], list)
        # Second element is list of Airfoil mean camber lines.
        self.assertIsInstance(result[1], list)

    def test_get_plottable_data_returns_ndarrays_for_each_cross_section(self) -> None:
        """Test that get_plottable_data returns one ndarray per WingCrossSection."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)
        assert result is not None

        num_wing_cross_sections = len(wing.wing_cross_sections)

        # It should have one outline array per WingCrossSection.
        self.assertEqual(len(result[0]), num_wing_cross_sections)
        # It should have one MCL array per WingCrossSection.
        self.assertEqual(len(result[1]), num_wing_cross_sections)

        for outline in result[0]:
            self.assertIsInstance(outline, np.ndarray)
            self.assertEqual(outline.shape[1], 3)  # These are 3D points.

        for mcl in result[1]:
            self.assertIsInstance(mcl, np.ndarray)
            self.assertEqual(mcl.shape[1], 3)  # These are 3D points.

    def test_get_plottable_data_three_section_wing(self) -> None:
        """Test get_plottable_data for Wing with 3 WingCrossSections."""
        wing = geometry_fixtures.make_three_section_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=False)
        assert result is not None

        # It should have 3 outlines and 3 MCLs.
        self.assertEqual(len(result[0]), 3)
        self.assertEqual(len(result[1]), 3)

    def test_get_plottable_data_type_4_symmetric_wing(self) -> None:
        """Test get_plottable_data for type 4 symmetric Wing."""
        wing = geometry_fixtures.make_symmetric_continuous_rectangular_wing_fixture()
        wing.generate_mesh(4)

        result = wing.get_plottable_data(show=False)
        assert result is not None

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_default_show_is_false(self) -> None:
        """Test that get_plottable_data default for show is False."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Call without a show parameter. It should return data (not None).
        result = wing.get_plottable_data()

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_accepts_numpy_bool(self) -> None:
        """Test that get_plottable_data accepts numpy bool for show parameter."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        result = wing.get_plottable_data(show=np.bool_(False))

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_invalid_show_type_raises(self) -> None:
        """Test that get_plottable_data raises error for invalid show type."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        bad_show: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            wing.get_plottable_data(show=bad_show)

    def test_get_plottable_data_with_panels_meshed(self) -> None:
        """Test get_plottable_data with meshed Panels."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        # Verify Panels are meshed.
        self.assertIsNotNone(wing.panels)

        result = wing.get_plottable_data(show=False)

        self.assertIsNotNone(result)


class TestWingTransformationMatrixCaching(unittest.TestCase):
    """Tests for Wing transformation matrix caching behavior."""

    def test_T_pas_G_Cg_to_Wn_Ler_returns_same_object_on_repeated_access(self) -> None:
        """Test that T_pas_G_Cg_to_Wn_Ler returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        T1 = wing.T_pas_G_Cg_to_Wn_Ler
        T2 = wing.T_pas_G_Cg_to_Wn_Ler

        self.assertIs(T1, T2)

    def test_T_pas_Wn_Ler_to_G_Cg_returns_same_object_on_repeated_access(self) -> None:
        """Test that T_pas_Wn_Ler_to_G_Cg returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        T1 = wing.T_pas_Wn_Ler_to_G_Cg
        T2 = wing.T_pas_Wn_Ler_to_G_Cg

        self.assertIs(T1, T2)

    def test_WnX_G_returns_same_object_on_repeated_access(self) -> None:
        """Test that WnX_G returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        v1 = wing.WnX_G
        v2 = wing.WnX_G

        self.assertIs(v1, v2)

    def test_WnY_G_returns_same_object_on_repeated_access(self) -> None:
        """Test that WnY_G returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        v1 = wing.WnY_G
        v2 = wing.WnY_G

        self.assertIs(v1, v2)

    def test_WnZ_G_returns_same_object_on_repeated_access(self) -> None:
        """Test that WnZ_G returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        v1 = wing.WnZ_G
        v2 = wing.WnZ_G

        self.assertIs(v1, v2)

    def test_children_T_pas_Wn_Ler_to_Wcs_Lp_returns_same_object_on_repeated_access(
        self,
    ) -> None:
        """Test that children_T_pas_Wn_Ler_to_Wcs_Lp returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_Wn_Ler_to_Wcs_Lp
        list2 = wing.children_T_pas_Wn_Ler_to_Wcs_Lp

        self.assertIs(list1, list2)

    def test_children_T_pas_Wcs_Lp_to_Wn_Ler_returns_same_object_on_repeated_access(
        self,
    ) -> None:
        """Test that children_T_pas_Wcs_Lp_to_Wn_Ler returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_Wcs_Lp_to_Wn_Ler
        list2 = wing.children_T_pas_Wcs_Lp_to_Wn_Ler

        self.assertIs(list1, list2)

    def test_children_T_pas_G_Cg_to_Wcs_Lp_returns_same_object_on_repeated_access(
        self,
    ) -> None:
        """Test that children_T_pas_G_Cg_to_Wcs_Lp returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_G_Cg_to_Wcs_Lp
        list2 = wing.children_T_pas_G_Cg_to_Wcs_Lp

        self.assertIs(list1, list2)

    def test_children_T_pas_Wcs_Lp_to_G_Cg_returns_same_object_on_repeated_access(
        self,
    ) -> None:
        """Test that children_T_pas_Wcs_Lp_to_G_Cg returns the same cached object."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        list1 = wing.children_T_pas_Wcs_Lp_to_G_Cg
        list2 = wing.children_T_pas_Wcs_Lp_to_G_Cg

        self.assertIs(list1, list2)

    def test_transformation_matrices_are_read_only(self) -> None:
        """Test that transformation matrices are read only."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        assert wing.T_pas_G_Cg_to_Wn_Ler is not None
        assert wing.T_pas_Wn_Ler_to_G_Cg is not None

        with self.assertRaises(ValueError):
            wing.T_pas_G_Cg_to_Wn_Ler[0, 0] = 999.0

        with self.assertRaises(ValueError):
            wing.T_pas_Wn_Ler_to_G_Cg[0, 0] = 999.0

    def test_basis_vectors_are_read_only(self) -> None:
        """Test that basis vectors are read only."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        wing.generate_mesh(1)

        assert wing.WnX_G is not None
        assert wing.WnY_G is not None
        assert wing.WnZ_G is not None

        with self.assertRaises(ValueError):
            wing.WnX_G[0] = 999.0

        with self.assertRaises(ValueError):
            wing.WnY_G[0] = 999.0

        with self.assertRaises(ValueError):
            wing.WnZ_G[0] = 999.0

    def test_children_transformation_matrices_are_read_only(self) -> None:
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


class TestExplodeIntoStripsMethods(unittest.TestCase):
    """This class contains unit tests for Wing._explode_wing,
    Wing._interpolate_between_wing_cross_sections, and the explode_into_strips
    parameter."""

    @staticmethod
    def _make_root_wing_cross_section() -> (
        ps.geometry.wing_cross_section.WingCrossSection
    ):
        """Create a root WingCrossSection with num_spanwise_panels=3."""
        return ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
            num_spanwise_panels=3,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing="uniform",
        )

    @staticmethod
    def _make_tip_wing_cross_section() -> (
        ps.geometry.wing_cross_section.WingCrossSection
    ):
        """Create a tip WingCrossSection with num_spanwise_panels=None."""
        return ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
            num_spanwise_panels=None,
            chord=0.5,
            Lp_Wcsp_Lpp=(0.0, 0.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            spanwise_spacing=None,
        )

    def _make_plain_wing(
        self, explode_into_strips: bool = False
    ) -> ps.geometry.wing.Wing:
        """Create a minimal two-WingCrossSection wing."""
        return ps.geometry.wing.Wing(
            wing_cross_sections=[
                self._make_root_wing_cross_section(),
                self._make_tip_wing_cross_section(),
            ],
            name="Test Wing",
            Ler_Gs_Cgs=(0.0, 0.0, 0.0),
            angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            explode_into_strips=explode_into_strips,
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
        )

    def test_explode_into_strips_false_wing_cross_section_count_unchanged(self) -> None:
        """Test that explode_into_strips=False keeps the original two WingCrossSections."""
        wing = self._make_plain_wing(explode_into_strips=False)
        self.assertEqual(len(wing.wing_cross_sections), 2)

    def test_explode_into_strips_true_correct_wing_cross_section_count(self) -> None:
        """Test that explode_into_strips=True with root num_spanwise=3 produces 4 WingCrossSections
        (root copy plus 3 interpolated including the tip)."""
        wing = self._make_plain_wing(explode_into_strips=True)
        self.assertEqual(len(wing.wing_cross_sections), 4)

    def test_explode_into_strips_true_non_tip_have_num_spanwise_one(self) -> None:
        """Test that all non-tip WingCrossSections have num_spanwise_panels=1 after explode."""
        wing = self._make_plain_wing(explode_into_strips=True)
        for wing_cross_section in wing.wing_cross_sections[:-1]:
            with self.subTest(wing_cross_section=wing_cross_section):
                self.assertEqual(wing_cross_section.num_spanwise_panels, 1)

    def test_explode_into_strips_true_tip_has_none_spanwise(self) -> None:
        """Test that the last WingCrossSection (tip) has num_spanwise_panels=None after explode."""
        wing = self._make_plain_wing(explode_into_strips=True)
        self.assertIsNone(wing.wing_cross_sections[-1].num_spanwise_panels)

    def test_interpolate_returns_n_sections(self) -> None:
        """Test that _interpolate_between_wing_cross_sections returns N WingCrossSections (the
        sections downstream of first_wing_cross_section, with no root copy), where N is first_wing_cross_section's spanwise
        panel count."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._interpolate_between_wing_cross_sections(
            self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()
        )
        # N = 3 => 3 interpolated sections
        self.assertEqual(len(result), 3)

    def test_interpolate_last_section_has_tip_chord(self) -> None:
        """Test that the last WingCrossSection in the result has the tip chord."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._interpolate_between_wing_cross_sections(
            self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()
        )
        self.assertAlmostEqual(result[-1].chord, 0.5)

    def test_interpolate_first_section_chord_linearly_interpolated(self) -> None:
        """Test that the first interpolated WingCrossSection chord is linearly interpolated
        between root (1.0) and tip (0.5)."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._interpolate_between_wing_cross_sections(
            self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()
        )
        # N = 3: the first section is at t = 1.0 / 3.0
        # chord at t = 1.0 / 3.0: (1.0 - 1.0 / 3.0) * 1.0 + (1.0 / 3.0) * 0.5 = 5.0 /
        # 6.0
        expected_first_interp = (1.0 - 1.0 / 3.0) * 1.0 + (1.0 / 3.0) * 0.5
        self.assertAlmostEqual(result[0].chord, expected_first_interp, places=10)

    def test_interpolate_Lp_y_divided_by_n(self) -> None:
        """Test that the Lp_Wcsp_Lpp y-component of each interpolated WingCrossSection is
        tip_Lp_y / N."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._interpolate_between_wing_cross_sections(
            self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()
        )
        # tip has Lp_y = 0.5. N = 3 => each section Lp_y = 0.5 / 3.0.
        expected_Lp_y = 0.5 / 3.0
        for wing_cross_section in result:
            with self.subTest(wing_cross_section=wing_cross_section):
                self.assertAlmostEqual(
                    float(wing_cross_section.Lp_Wcsp_Lpp[1]), expected_Lp_y, places=10
                )

    def test_explode_wing_with_two_wing_cross_sections_returns_correct_count(
        self,
    ) -> None:
        """Test that _explode_wing with a two-WingCrossSection input (root: num=3, tip) returns 4
        WingCrossSections."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._explode_wing(
            [self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()]
        )
        self.assertEqual(len(result), 4)

    def test_explode_wing_first_wing_cross_section_is_root(self) -> None:
        """Test that _explode_wing seeds the result with the root WingCrossSection (root chord, a
        single spanwise panel)."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._explode_wing(
            [self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()]
        )
        self.assertAlmostEqual(result[0].chord, 1.0)
        self.assertEqual(result[0].num_spanwise_panels, 1)

    def test_explode_wing_all_non_tip_have_num_spanwise_one(self) -> None:
        """Test that _explode_wing produces WingCrossSections where every non-tip entry has
        num_spanwise_panels=1."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._explode_wing(
            [self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()]
        )
        for wing_cross_section in result[:-1]:
            with self.subTest(wing_cross_section=wing_cross_section):
                self.assertEqual(wing_cross_section.num_spanwise_panels, 1)

    def test_explode_wing_last_wing_cross_section_is_tip(self) -> None:
        """Test that _explode_wing produces a final WingCrossSection with num_spanwise_panels=None."""
        wing = self._make_plain_wing(explode_into_strips=False)
        result = wing._explode_wing(
            [self._make_root_wing_cross_section(), self._make_tip_wing_cross_section()]
        )
        self.assertIsNone(result[-1].num_spanwise_panels)

    def test_explode_wing_rejects_non_uniform_spanwise_spacing(self) -> None:
        """Test that _explode_wing raises ValueError when a non tip WingCrossSection uses cosine
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
            wing._explode_wing([cosine_root, self._make_tip_wing_cross_section()])

    def test_spanwise_mesh_default_is_trapezoidal(self) -> None:
        """Test that a Wing built without explode_into_strips has a trapezoidal spanwise
        mesh marker."""
        wing = self._make_plain_wing(explode_into_strips=False)
        self.assertEqual(wing.spanwise_mesh, "trapezoidal")

    def test_spanwise_mesh_exploded_is_exploded(self) -> None:
        """Test that a Wing built with explode_into_strips has an exploded spanwise mesh
        marker."""
        wing = self._make_plain_wing(explode_into_strips=True)
        self.assertEqual(wing.spanwise_mesh, "exploded")

    def test_spanwise_mesh_is_read_only(self) -> None:
        """Test that the spanwise_mesh marker cannot be reassigned."""
        wing = self._make_plain_wing(explode_into_strips=False)
        with self.assertRaises(AttributeError):
            setattr(wing, "spanwise_mesh", "exploded")

    def test_spanwise_mesh_preserved_by_deepcopy(self) -> None:
        """Test that deep copying a Wing preserves its spanwise mesh marker."""
        wing = self._make_plain_wing(explode_into_strips=True)
        wing_copy = copy.deepcopy(wing)
        self.assertEqual(wing_copy.spanwise_mesh, "exploded")


class TestFromEdgePoints(unittest.TestCase):
    """This class contains unit tests for the Wing.from_edge_points constructor."""

    @staticmethod
    def _straight_edge_points(
        num_input_points: int = 11,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build densely sampled straight leading and trailing edge curves.

        The leading edge is a straight backward sweep from the origin to (0.5, 1.0,
        0.0), so x = 0.5 * y. The trailing edge is a straight, unswept line at x =
        1.0. The chord is therefore 1.0 - 0.5 * y, tapering from a unit root chord to
        a 0.5 tip chord. PCHIP reproduces these linear curves exactly, so the
        resampled WingCrossSections take predictable values.
        """
        ys = np.linspace(0.0, 1.0, num_input_points)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((0.5 * ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        return leading, trailing

    @staticmethod
    def _pointed_tip_edge_points(
        num_input_points: int = 11,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build straight edge curves for a planform that tapers to a point at the tip.

        The leading edge runs from the origin to (1.0, 1.0, 0.0) and the trailing
        edge is unswept at x = 1.0, so the chord is 1.0 - y and falls to zero exactly
        at the tip. Such a planform is invalid without a tip trim.
        """
        ys = np.linspace(0.0, 1.0, num_input_points)
        zeros = np.zeros_like(ys)
        leading = np.column_stack((ys, ys, zeros))
        trailing = np.column_stack((np.ones_like(ys), ys, zeros))
        return leading, trailing

    def _make_edge_wing(
        self, num_wing_cross_sections: int = 5, tip_trim_fraction: float = 0.0
    ) -> ps.geometry.wing.Wing:
        """Build a from_edge_points Wing from the straight tapered edge curves."""
        leading, trailing = self._straight_edge_points()
        return ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=num_wing_cross_sections,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            name="Edge Wing",
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
            tip_trim_fraction=tip_trim_fraction,
        )

    def test_returns_wing(self) -> None:
        """Test that from_edge_points returns a Wing instance."""
        wing = self._make_edge_wing()
        self.assertIsInstance(wing, ps.geometry.wing.Wing)

    def test_wing_cross_section_count(self) -> None:
        """Test that from_edge_points produces num_wing_cross_sections WingCrossSections."""
        wing = self._make_edge_wing(num_wing_cross_sections=5)
        self.assertEqual(len(wing.wing_cross_sections), 5)

    def test_spanwise_mesh_is_edge_defined(self) -> None:
        """Test that a from_edge_points Wing reports an edge_defined spanwise mesh."""
        wing = self._make_edge_wing()
        self.assertEqual(wing.spanwise_mesh, "edge_defined")

    def test_non_tip_wing_cross_sections_have_num_spanwise_one(self) -> None:
        """Test that every non tip WingCrossSection has num_spanwise_panels of one."""
        wing = self._make_edge_wing()
        for wing_cross_section in wing.wing_cross_sections[:-1]:
            with self.subTest(wing_cross_section=wing_cross_section):
                self.assertEqual(wing_cross_section.num_spanwise_panels, 1)

    def test_tip_wing_cross_section_has_none_spanwise(self) -> None:
        """Test that the tip WingCrossSection has num_spanwise_panels of None."""
        wing = self._make_edge_wing()
        self.assertIsNone(wing.wing_cross_sections[-1].num_spanwise_panels)

    def test_chords_match_linear_taper(self) -> None:
        """Test that the resampled chords follow the expected linear taper."""
        wing = self._make_edge_wing(num_wing_cross_sections=5)
        ys = np.linspace(0.0, 1.0, 5)
        expected_chords = 1.0 - 0.5 * ys
        for wing_cross_section, expected_chord in zip(
            wing.wing_cross_sections, expected_chords
        ):
            with self.subTest(expected_chord=expected_chord):
                self.assertAlmostEqual(
                    wing_cross_section.chord, float(expected_chord), places=10
                )

    def test_root_offset_is_zero(self) -> None:
        """Test that the root WingCrossSection has a zero leading point offset."""
        wing = self._make_edge_wing()
        npt.assert_array_equal(
            wing.wing_cross_sections[0].Lp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0])
        )

    def test_leading_point_offsets_match_expected(self) -> None:
        """Test that each non root leading point offset is the spanwise step of the
        straight swept leading edge."""
        wing = self._make_edge_wing(num_wing_cross_sections=5)
        for wing_cross_section in wing.wing_cross_sections[1:]:
            with self.subTest(wing_cross_section=wing_cross_section):
                npt.assert_array_almost_equal(
                    wing_cross_section.Lp_Wcsp_Lpp,
                    np.array([0.125, 0.25, 0.0]),
                    decimal=10,
                )

    def test_all_angle_vectors_are_zero(self) -> None:
        """Test that every WingCrossSection keeps a zero angle vector (untwisted)."""
        wing = self._make_edge_wing()
        for wing_cross_section in wing.wing_cross_sections:
            with self.subTest(wing_cross_section=wing_cross_section):
                npt.assert_array_equal(
                    wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
                    np.array([0.0, 0.0, 0.0]),
                )

    def test_stored_leading_edge_points_match_input(self) -> None:
        """Test that the stored leading edge curve matches the supplied points."""
        leading, _ = self._straight_edge_points()
        wing = self._make_edge_wing()
        npt.assert_array_equal(wing.leadingEdgePoints_Wn_Ler, leading)

    def test_stored_trailing_edge_points_match_input(self) -> None:
        """Test that the stored trailing edge curve matches the supplied points."""
        _, trailing = self._straight_edge_points()
        wing = self._make_edge_wing()
        npt.assert_array_equal(wing.trailingEdgePoints_Wn_Ler, trailing)

    def test_stored_curves_are_read_only(self) -> None:
        """Test that the stored edge curves cannot be mutated in place."""
        wing = self._make_edge_wing()
        assert wing.leadingEdgePoints_Wn_Ler is not None
        assert wing.trailingEdgePoints_Wn_Ler is not None
        self.assertFalse(wing.leadingEdgePoints_Wn_Ler.flags.writeable)
        self.assertFalse(wing.trailingEdgePoints_Wn_Ler.flags.writeable)
        with self.assertRaises(ValueError):
            wing.leadingEdgePoints_Wn_Ler[0, 0] = 1.0
        with self.assertRaises(ValueError):
            wing.trailingEdgePoints_Wn_Ler[0, 0] = 1.0

    def test_edge_properties_are_read_only(self) -> None:
        """Test that the edge curve and tip trim properties cannot be reassigned."""
        wing = self._make_edge_wing()
        with self.assertRaises(AttributeError):
            setattr(wing, "leadingEdgePoints_Wn_Ler", None)
        with self.assertRaises(AttributeError):
            setattr(wing, "trailingEdgePoints_Wn_Ler", None)
        with self.assertRaises(AttributeError):
            setattr(wing, "tip_trim_fraction", 0.5)

    def test_tip_trim_fraction_defaults_to_zero(self) -> None:
        """Test that the stored tip trim fraction defaults to zero."""
        wing = self._make_edge_wing()
        self.assertEqual(wing.tip_trim_fraction, 0.0)

    def test_tip_trim_fraction_stored(self) -> None:
        """Test that a non zero tip trim fraction is stored on the Wing."""
        wing = self._make_edge_wing(tip_trim_fraction=0.25)
        self.assertEqual(wing.tip_trim_fraction, 0.25)

    def test_edge_attributes_none_for_normal_wing(self) -> None:
        """Test that a Wing not built from edge points reports None edge attributes."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        self.assertIsNone(wing.leadingEdgePoints_Wn_Ler)
        self.assertIsNone(wing.trailingEdgePoints_Wn_Ler)
        self.assertIsNone(wing.tip_trim_fraction)

    def test_non_symmetric_wing_cross_sections_have_no_control_surface_type(
        self,
    ) -> None:
        """Test that a non symmetric from_edge_points Wing leaves every control surface
        symmetry type None, as symmetry types 1 through 3 require."""
        wing = self._make_edge_wing()
        for wing_cross_section in wing.wing_cross_sections:
            with self.subTest(wing_cross_section=wing_cross_section):
                self.assertIsNone(wing_cross_section.control_surface_symmetry_type)

    def test_symmetric_wing_cross_sections_have_symmetric_control_surface_type(
        self,
    ) -> None:
        """Test that a symmetric from_edge_points Wing marks every control surface
        symmetry type symmetric, as symmetry types 4 and 5 require."""
        leading, trailing = self._straight_edge_points()
        wing = ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        )
        for wing_cross_section in wing.wing_cross_sections:
            with self.subTest(wing_cross_section=wing_cross_section):
                self.assertEqual(
                    wing_cross_section.control_surface_symmetry_type, "symmetric"
                )

    def test_symmetric_wing_meshes_through_airplane(self) -> None:
        """Test that a symmetric from_edge_points Wing meshes through an Airplane, the
        path that requires a non None control surface symmetry type on every
        WingCrossSection."""
        leading, trailing = self._straight_edge_points()
        wing = ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            symmetric=True,
            symmetryNormal_G=(0.0, 1.0, 0.0),
            symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            num_chordwise_panels=2,
            chordwise_spacing="uniform",
        )
        airplane = ps.geometry.airplane.Airplane(wings=[wing], name="Edge Airplane")
        meshed_wing = airplane.wings[0]
        self.assertEqual(meshed_wing.symmetry_type, 4)
        self.assertIsNotNone(meshed_wing.panels)

    def test_tip_trim_moves_outermost_section_inboard(self) -> None:
        """Test that a tip trim resamples over a shortened span, leaving the outermost
        WingCrossSection with a finite chord inboard of the geometric tip."""
        wing = self._make_edge_wing(num_wing_cross_sections=5, tip_trim_fraction=0.2)
        # The trimmed tip sits at y = 0.8, where the chord is 1.0 - 0.5 * 0.8 = 0.6.
        self.assertAlmostEqual(wing.wing_cross_sections[-1].chord, 0.6, places=10)

    def test_pointed_tip_rejected_without_trim(self) -> None:
        """Test that a planform tapering to a point at the tip is rejected without a
        tip trim, since the tip chord would be zero."""
        leading, trailing = self._pointed_tip_edge_points()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_pointed_tip_accepted_with_trim(self) -> None:
        """Test that a tip trim lets a planform tapering to a point at the tip build
        with a finite tip chord."""
        leading, trailing = self._pointed_tip_edge_points()
        wing = ps.geometry.wing.Wing.from_edge_points(
            leadingEdgePoints_Wn_Ler=leading,
            trailingEdgePoints_Wn_Ler=trailing,
            num_wing_cross_sections=5,
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            tip_trim_fraction=0.2,
        )
        # The trimmed tip sits at y = 0.8, where the chord is 1.0 - 0.8 = 0.2.
        self.assertAlmostEqual(wing.wing_cross_sections[-1].chord, 0.2, places=10)

    def test_deepcopy_preserves_edge_attributes(self) -> None:
        """Test that deep copying a from_edge_points Wing preserves and isolates its
        edge attributes."""
        wing = self._make_edge_wing(tip_trim_fraction=0.1)
        wing_copy = copy.deepcopy(wing)
        self.assertEqual(wing_copy.spanwise_mesh, "edge_defined")
        self.assertEqual(wing_copy.tip_trim_fraction, 0.1)
        npt.assert_array_equal(
            wing_copy.leadingEdgePoints_Wn_Ler, wing.leadingEdgePoints_Wn_Ler
        )
        npt.assert_array_equal(
            wing_copy.trailingEdgePoints_Wn_Ler, wing.trailingEdgePoints_Wn_Ler
        )
        # The copied curves are independent and remain read only.
        self.assertIsNot(
            wing_copy.leadingEdgePoints_Wn_Ler, wing.leadingEdgePoints_Wn_Ler
        )
        assert wing_copy.leadingEdgePoints_Wn_Ler is not None
        assert wing_copy.trailingEdgePoints_Wn_Ler is not None
        self.assertFalse(wing_copy.leadingEdgePoints_Wn_Ler.flags.writeable)
        self.assertFalse(wing_copy.trailingEdgePoints_Wn_Ler.flags.writeable)

    def test_rejects_non_increasing_leading_edge_y(self) -> None:
        """Test that a leading edge with non increasing y components is rejected."""
        leading, trailing = self._straight_edge_points()
        leading[2, 1] = leading[1, 1]
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_non_planar_points(self) -> None:
        """Test that a non zero z component anywhere is rejected."""
        leading, trailing = self._straight_edge_points()
        leading[3, 2] = 0.1
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_leading_edge_not_anchored_at_origin(self) -> None:
        """Test that a leading edge whose first point is not the origin is rejected."""
        leading, trailing = self._straight_edge_points()
        leading[0, 0] = 0.1
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_trailing_edge_root_without_root_chord(self) -> None:
        """Test that a trailing edge whose first point lacks a positive root chord is
        rejected."""
        leading, trailing = self._straight_edge_points()
        trailing[0, 0] = 0.0
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_mismatched_tip_y(self) -> None:
        """Test that curves spanning different maximum y values are rejected."""
        leading, trailing = self._straight_edge_points()
        trailing[-1, 1] = 0.9
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_too_few_points(self) -> None:
        """Test that a curve with fewer than two points is rejected."""
        _, trailing = self._straight_edge_points()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=np.array([[0.0, 0.0, 0.0]]),
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_too_few_wing_cross_sections(self) -> None:
        """Test that fewer than two WingCrossSections is rejected."""
        leading, trailing = self._straight_edge_points()
        with self.assertRaises(ValueError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=1,
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            )

    def test_rejects_non_airfoil(self) -> None:
        """Test that an airfoil argument that is not an Airfoil is rejected."""
        leading, trailing = self._straight_edge_points()
        bad_airfoil: Any = "naca0012"
        with self.assertRaises(TypeError):
            ps.geometry.wing.Wing.from_edge_points(
                leadingEdgePoints_Wn_Ler=leading,
                trailingEdgePoints_Wn_Ler=trailing,
                num_wing_cross_sections=5,
                airfoil=bad_airfoil,
            )

    def test_rejects_tip_trim_fraction_out_of_range(self) -> None:
        """Test that a tip trim fraction outside [0, 1) is rejected."""
        leading, trailing = self._straight_edge_points()
        for bad_fraction in (-0.1, 1.0, 1.5):
            with self.subTest(tip_trim_fraction=bad_fraction):
                with self.assertRaises(ValueError):
                    ps.geometry.wing.Wing.from_edge_points(
                        leadingEdgePoints_Wn_Ler=leading,
                        trailingEdgePoints_Wn_Ler=trailing,
                        num_wing_cross_sections=5,
                        airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                        tip_trim_fraction=bad_fraction,
                    )
