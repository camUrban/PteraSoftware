"""This module contains classes to test Airplanes."""

import copy
import unittest
from collections.abc import Sequence
from typing import Any
from unittest.mock import PropertyMock, patch

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestAirplane(unittest.TestCase):
    """This is a class with functions to test Airplanes."""

    def setUp(self) -> None:
        """Set up test fixtures for Airplane tests."""
        # Create fixtures for different Airplane types
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()
        self.first_airplane = geometry_fixtures.make_first_airplane_fixture()
        self.multi_wing_airplane = geometry_fixtures.make_multi_wing_airplane_fixture()
        self.type_5_wing_airplane = (
            geometry_fixtures.make_type_5_wing_airplane_fixture()
        )
        self.custom_reference_airplane = (
            geometry_fixtures.make_custom_reference_airplane_fixture()
        )

        # Create additional test fixtures
        self.test_wing_type_1 = geometry_fixtures.make_type_1_wing_fixture()

    def test_wings_parameter_validation(self) -> None:
        """Test that wings parameter validation works correctly."""
        # Test empty list raises error
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[])

        # Test non-list raises error
        bad_wings: Any = "not a list"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=bad_wings)

        # Test non-Wing objects raise error
        non_wing_list: Any = ["not a wing"]
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=non_wing_list)

        # Test mixed valid and invalid Wings
        mixed_wing_list: Any = [self.test_wing_type_1, "invalid"]
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=mixed_wing_list)

    def test_name_parameter_validation(self) -> None:
        """Test name parameter validation."""
        # Test a valid string name. Create a fresh fixture since Wings can only be
        # processed once.
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        airplane = ps.geometry.airplane.Airplane(
            wings=[test_wing], name="Valid Test Name"
        )
        self.assertEqual(airplane.name, "Valid Test Name")

        # Test invalid name types
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        bad_name: Any = 123
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], name=bad_name)

        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        none_name: Any = None
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], name=none_name)

    def test_Cg_GP1_CgP1_parameter_validation(self) -> None:
        """Test Cg_GP1_CgP1 parameter validation."""
        # Test valid 3D vectors
        valid_positions: list[np.ndarray | Sequence[float]] = [
            [0.0, 0.0, 0.0],
            [1.0, -2.0, 3.5],
            np.array([0.5, 1.5, -0.8]),
            (2.0, -1.0, 0.0),
        ]

        for position in valid_positions:
            with self.subTest(position=position):
                # Create fresh fixture since Wings can only be processed once
                test_wing = geometry_fixtures.make_type_1_wing_fixture()
                airplane = ps.geometry.airplane.Airplane(
                    wings=[test_wing], Cg_GP1_CgP1=position
                )
                npt.assert_array_equal(airplane.Cg_GP1_CgP1, position)

        # Test invalid positions
        invalid_positions: list[Any] = [
            [1.0, 2.0],  # Wrong size
            [1.0, 2.0, 3.0, 4.0],  # Wrong size
            "not a vector",  # String
            None,  # None
        ]

        for invalid_position in invalid_positions:
            with self.subTest(invalid_position=invalid_position):
                # Create fresh fixture since Wings can only be processed once
                test_wing = geometry_fixtures.make_type_1_wing_fixture()
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.airplane.Airplane(
                        wings=[test_wing], Cg_GP1_CgP1=invalid_position
                    )

    def test_weight_parameter_validation(self) -> None:
        """Test weight parameter validation."""
        # Test valid weights (non-negative)
        valid_weights = [0.0, 0.1, 100.0, 5000.0, 10000.0]

        for weight in valid_weights:
            with self.subTest(weight=weight):
                # Create fresh fixture since Wings can only be processed once
                test_wing = geometry_fixtures.make_type_1_wing_fixture()
                airplane = ps.geometry.airplane.Airplane(
                    wings=[test_wing], weight=weight
                )
                self.assertEqual(airplane.weight, weight)

        # Test invalid weights (negative)
        invalid_weights = [-0.1, -100.0, -1000.0]

        for invalid_weight in invalid_weights:
            with self.subTest(invalid_weight=invalid_weight):
                # Create fresh fixture since Wings can only be processed once
                test_wing = geometry_fixtures.make_type_1_wing_fixture()
                with self.assertRaises(ValueError):
                    ps.geometry.airplane.Airplane(
                        wings=[test_wing], weight=invalid_weight
                    )

        # Test invalid weight types. Create fresh fixture since Wings can only be
        # processed once.
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        bad_weight: Any = "heavy"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], weight=bad_weight)

    def test_reference_dimensions_default_behavior(self) -> None:
        """Test reference dimensions default to first Wing's properties."""
        # Create Airplane with no explicit reference dimensions
        airplane = ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1])

        # Reference dimensions should be populated from the first Wing
        first_wing = airplane.wings[0]
        self.assertEqual(airplane.s_ref, first_wing.projected_area)
        self.assertEqual(airplane.c_ref, first_wing.mean_aerodynamic_chord)
        self.assertEqual(airplane.b_ref, first_wing.span)

    def test_reference_dimensions_explicit_values(self) -> None:
        """Test reference dimensions with explicit values."""
        # Test custom reference Airplane
        self.assertEqual(self.custom_reference_airplane.s_ref, 15.0)
        self.assertEqual(self.custom_reference_airplane.c_ref, 2.0)
        self.assertEqual(self.custom_reference_airplane.b_ref, 10.0)

        # Test validation of reference dimensions. Create fresh fixtures since Wings can
        # only be processed once.
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[test_wing], s_ref=-1.0)

        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[test_wing], c_ref=0.0)

        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        bad_b_ref: Any = "large"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], b_ref=bad_b_ref)

    def test_s_ref_none_with_none_projected_area_raises(self) -> None:
        """Test that s_ref=None raises ValueError when wing's projected_area is None."""
        with patch.object(
            ps.geometry.wing.Wing,
            "projected_area",
            new_callable=PropertyMock,
            return_value=None,
        ):
            test_wing = geometry_fixtures.make_type_1_wing_fixture()
            with self.assertRaises(ValueError):
                ps.geometry.airplane.Airplane(wings=[test_wing])

    def test_c_ref_none_with_none_mean_aerodynamic_chord_raises(self) -> None:
        """Test that c_ref=None raises ValueError when wing's mean_aerodynamic_chord is
        None."""
        with patch.object(
            ps.geometry.wing.Wing,
            "mean_aerodynamic_chord",
            new_callable=PropertyMock,
            return_value=None,
        ):
            test_wing = geometry_fixtures.make_type_1_wing_fixture()
            with self.assertRaises(ValueError):
                ps.geometry.airplane.Airplane(wings=[test_wing], s_ref=2.0)

    def test_b_ref_none_with_none_span_raises(self) -> None:
        """Test that b_ref=None raises ValueError when wing's span is None."""
        with patch.object(
            ps.geometry.wing.Wing,
            "span",
            new_callable=PropertyMock,
            return_value=None,
        ):
            test_wing = geometry_fixtures.make_type_1_wing_fixture()
            with self.assertRaises(ValueError):
                ps.geometry.airplane.Airplane(wings=[test_wing], s_ref=2.0, c_ref=1.0)

    def test_num_panels_calculation(self) -> None:
        """Test that num_panels is calculated correctly from all Wings."""
        # Single-Wing Airplane
        single_wing_panels = self.basic_airplane.wings[0].num_panels
        self.assertEqual(self.basic_airplane.num_panels, single_wing_panels)

        # Multi-Wing Airplane
        wing_panel_counts: list[int] = []
        for wing in self.multi_wing_airplane.wings:
            assert wing.num_panels is not None
            wing_panel_counts.append(wing.num_panels)
        expected_panels = sum(wing_panel_counts)
        self.assertEqual(self.multi_wing_airplane.num_panels, expected_panels)

    def test_force_moment_attributes_initialization(self) -> None:
        """Test that force and moment attributes are initialized to None."""
        airplane = self.basic_airplane

        # All force/moment attributes should be None initially
        self.assertIsNone(airplane.forces_W)
        self.assertIsNone(airplane.forceCoefficients_W)
        self.assertIsNone(airplane.moments_W_CgP1)
        self.assertIsNone(airplane.momentCoefficients_W_CgP1)
        self.assertIsNone(airplane.forces_G)
        self.assertIsNone(airplane.forceCoefficients_G)
        self.assertIsNone(airplane.moments_G_Cg)
        self.assertIsNone(airplane.momentCoefficients_G_Cg)
        self.assertIsNone(airplane.moments_W_Cg)
        self.assertIsNone(airplane.momentCoefficients_W_Cg)

    def test_validate_first_airplane_constraints_valid(self) -> None:
        """Test validate_first_airplane_constraints with valid first Airplane."""
        # First Airplane should pass validation (Cg_GP1_CgP1 is all zeros)
        try:
            self.first_airplane.validate_first_airplane_constraints()
        except Exception as e:
            self.fail(f"First airplane validation failed unexpectedly: {e}")

    def test_validate_first_airplane_constraints_invalid(self) -> None:
        """Test validate_first_airplane_constraints with invalid Airplane."""
        # Basic Airplane should fail validation (Cg_GP1_CgP1 is not all zeros)
        with self.assertRaises(ValueError):
            self.basic_airplane.validate_first_airplane_constraints()

        # Custom reference Airplane should fail validation
        with self.assertRaises(ValueError):
            self.custom_reference_airplane.validate_first_airplane_constraints()

    def test_multi_wing_configuration(self) -> None:
        """Test Airplane with multiple Wings."""
        airplane = self.multi_wing_airplane

        # Should have multiple Wings
        self.assertEqual(len(airplane.wings), 2)
        self.assertIsInstance(airplane.wings[0], ps.geometry.wing.Wing)
        self.assertIsInstance(airplane.wings[1], ps.geometry.wing.Wing)

        # Wings should have different names
        self.assertNotEqual(airplane.wings[0].name, airplane.wings[1].name)

    def test_type_5_wing_processing(self) -> None:
        """Test that type 5 Wings are processed correctly into two Wings."""
        airplane = self.type_5_wing_airplane

        # Type 5 Wing should be split into two Wings during initialization
        self.assertEqual(len(airplane.wings), 2)

        # First Wing should be the original (now type 1)
        first_wing = airplane.wings[0]
        self.assertFalse(first_wing.symmetric)
        self.assertFalse(first_wing.mirror_only)

        # Second Wing should be the reflected Wing (type 3)
        second_wing = airplane.wings[1]
        self.assertFalse(second_wing.symmetric)
        self.assertTrue(second_wing.mirror_only)
        self.assertTrue(second_wing.name.startswith("Reflected"))

    def test_process_wing_symmetry_type_1(self) -> None:
        """Test process_wing_symmetry with type 1 Wing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 1)

    def test_process_wing_symmetry_type_2(self) -> None:
        """Test process_wing_symmetry with type 2 Wing."""
        wing = geometry_fixtures.make_type_2_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 2)

    def test_process_wing_symmetry_type_3(self) -> None:
        """Test process_wing_symmetry with type 3 Wing."""
        wing = geometry_fixtures.make_type_3_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 3)

    def test_process_wing_symmetry_type_4(self) -> None:
        """Test process_wing_symmetry with type 4 Wing."""
        wing = geometry_fixtures.make_type_4_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 4)

    def test_process_wing_symmetry_type_5(self) -> None:
        """Test process_wing_symmetry with type 5 Wing."""
        wing = geometry_fixtures.make_type_5_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with two Wings
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

        # First Wing should be modified to type 1
        first_wing = result[0]
        self.assertEqual(first_wing.symmetry_type, 1)
        self.assertFalse(first_wing.symmetric)
        self.assertFalse(first_wing.mirror_only)
        self.assertIsNone(first_wing.symmetryNormal_G)
        self.assertIsNone(first_wing.symmetryPoint_G_Cg)

        # Second Wing should be reflected type 3
        second_wing = result[1]
        self.assertEqual(second_wing.symmetry_type, 3)
        self.assertFalse(second_wing.symmetric)
        self.assertTrue(second_wing.mirror_only)
        self.assertTrue(second_wing.name.startswith("Reflected"))

    def test_process_wing_symmetry_control_surface_validation_types_1_2_3(self) -> None:
        """Test control surface validation for symmetry types 1, 2, 3."""
        # Type 1: should fail with control surfaces. Create fresh fixtures since
        # WingCrossSections can only be validated once.
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            # This fixture has a control surface.
            geometry_fixtures.make_basic_wing_cross_section_fixture(),
        ]
        with self.assertRaises(ValueError):
            wing_type_1 = ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                symmetric=False,
                mirror_only=False,
            )
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing_type_1)

        # Type 2: should fail with control surfaces. Create fresh fixtures since
        # WingCrossSections can only be validated once.
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            # This fixture has a control surface.
            geometry_fixtures.make_basic_wing_cross_section_fixture(),
        ]
        with self.assertRaises(ValueError):
            wing_type_2 = ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                symmetric=False,
                mirror_only=True,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing_type_2)

    def test_process_wing_symmetry_control_surface_validation_types_4_5(self) -> None:
        """Test control surface validation for symmetry types 4, 5."""
        # Create Wings without control surface configurations
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            # This fixture has no control surface.
            geometry_fixtures.make_minimal_wing_cross_section_fixture(),
        ]

        # Type 4: should fail without control surfaces
        with self.assertRaises(ValueError):
            wing_type_4 = ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=[0.0, 1.0, 0.0],
                symmetryPoint_G_Cg=[0.0, 0.0, 0.0],
            )
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing_type_4)

    def test_process_wing_symmetry_type_5_control_surface_deflections(self) -> None:
        """Test type 5 Wing processing with different control surface deflections."""
        # Create the asymmetric control surface WingCrossSections. Use the root fixture
        # with an asymmetric control surface already configured.
        root_wing_cross_section = (
            geometry_fixtures.make_root_asymmetric_control_surface_wing_cross_section_fixture()
        )
        tip_wing_cross_section = (
            geometry_fixtures.make_asymmetric_control_surface_wing_cross_section_fixture()
        )

        # Create type 5 Wing with asymmetric control surfaces
        wing_cross_sections = [root_wing_cross_section, tip_wing_cross_section]

        wing = ps.geometry.wing.Wing(
            wing_cross_sections=wing_cross_sections,
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=[0.0, 0.707, 0.707],
            symmetryPoint_G_Cg=[0.5, 0.0, 0.0],
        )

        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return two Wings
        self.assertEqual(len(result), 2)

        # Check that asymmetric control surface deflections are negated in the reflected
        # Wing.
        original_wing = result[0]
        reflected_wing = result[1]

        # Find corresponding WingCrossSections with asymmetric control surfaces
        for i, wing_cross_section in enumerate(original_wing.wing_cross_sections):
            reflected_wing_cross_section = reflected_wing.wing_cross_sections[i]
            # Reflected Wing should have None-type control surface symmetry
            self.assertEqual(
                reflected_wing_cross_section.control_surface_symmetry_type, None
            )

    def test_process_wing_symmetry_type_4_asymmetric_root_raises(self) -> None:
        """A type 4 Wing whose root cross section (on the coincident symmetry plane) has
        an asymmetric control surface with a nonzero deflection must raise.

        The original and mirrored halves would deflect that shared cross section in
        opposite directions, tearing the mesh at the centerline seam.
        """
        root_wing_cross_section = (
            geometry_fixtures.make_root_asymmetric_control_surface_wing_cross_section_fixture()
        )
        tip_wing_cross_section = (
            geometry_fixtures.make_tip_wing_cross_section_with_control_surface_fixture()
        )
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
            Ler_Gs_Cgs=[1.0, 0.0, 0.5],
            angles_Gs_to_Wn_ixyz=[0.0, 0.0, 0.0],
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=[0.0, 1.0, 0.0],
            symmetryPoint_G_Cg=[1.0, 0.0, 0.5],
            num_chordwise_panels=8,
            chordwise_spacing="cosine",
        )
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

    def test_process_wing_symmetry_type_4_asymmetric_tip_allowed(self) -> None:
        """A type 4 Wing may carry an asymmetric control surface on an off-plane cross
        section (an outboard aileron), since only the shared root cross section tears
        the mesh.

        The guard must reject the root case without over-restricting this one.
        """
        root_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        root_wing_cross_section.control_surface_symmetry_type = "symmetric"
        tip_wing_cross_section = (
            geometry_fixtures.make_asymmetric_control_surface_wing_cross_section_fixture()
        )
        wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wing_cross_section, tip_wing_cross_section],
            Ler_Gs_Cgs=[1.0, 0.0, 0.5],
            angles_Gs_to_Wn_ixyz=[0.0, 0.0, 0.0],
            symmetric=True,
            mirror_only=False,
            symmetryNormal_G=[0.0, 1.0, 0.0],
            symmetryPoint_G_Cg=[1.0, 0.0, 0.5],
            num_chordwise_panels=8,
            chordwise_spacing="cosine",
        )
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 4)

    def test_airplane_with_various_wing_combinations(self) -> None:
        """Test Airplane with various combinations of Wing types."""
        # Mix of different Wing types
        wings = [
            geometry_fixtures.make_type_1_wing_fixture(),
            geometry_fixtures.make_type_2_wing_fixture(),
            geometry_fixtures.make_type_4_wing_fixture(),
        ]

        airplane = ps.geometry.airplane.Airplane(
            wings=wings, name="Mixed Wing Type Airplane"
        )

        # Should have at least the original number of Wings (type 5 could add more)
        self.assertGreaterEqual(len(airplane.wings), 3)

        # All Wings should be processed and meshed
        for wing in airplane.wings:
            self.assertIsNotNone(wing.symmetry_type)
            self.assertIsNotNone(wing.panels)


class TestAirplaneDeepCopy(unittest.TestCase):
    """Tests for Airplane.__deepcopy__ method."""

    def setUp(self) -> None:
        """Set up test fixtures for deepcopy tests."""
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()
        self.first_airplane = geometry_fixtures.make_first_airplane_fixture()
        self.multi_wing_airplane = geometry_fixtures.make_multi_wing_airplane_fixture()

    def test_deepcopy_creates_new_instance(self) -> None:
        """Test that deepcopy creates a new Airplane instance."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.geometry.airplane.Airplane)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_airplane_parameters(self) -> None:
        """Test that deepcopy preserves Airplane parameters."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.weight, original.weight)
        self.assertEqual(copied.s_ref, original.s_ref)
        self.assertEqual(copied.c_ref, original.c_ref)
        self.assertEqual(copied.b_ref, original.b_ref)
        npt.assert_array_equal(copied.Cg_GP1_CgP1, original.Cg_GP1_CgP1)

    def test_deepcopy_creates_independent_cg_array(self) -> None:
        """Test that deepcopy creates an independent copy of Cg_GP1_CgP1."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.Cg_GP1_CgP1, original.Cg_GP1_CgP1)

    def test_deepcopy_creates_independent_wings(self) -> None:
        """Test that deepcopy creates independent Wing copies."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(len(copied.wings), len(original.wings))
        for orig_wing, copied_wing in zip(original.wings, copied.wings):
            self.assertIsNot(orig_wing, copied_wing)
            self.assertEqual(copied_wing.name, orig_wing.name)

    def test_deepcopy_multi_wing_airplane(self) -> None:
        """Test that deepcopy works correctly for multi wing Airplanes."""
        original = self.multi_wing_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(len(copied.wings), len(original.wings))
        for i, (orig_wing, copied_wing) in enumerate(zip(original.wings, copied.wings)):
            with self.subTest(wing_index=i):
                self.assertIsNot(orig_wing, copied_wing)
                self.assertEqual(copied_wing.symmetry_type, orig_wing.symmetry_type)

    def test_deepcopy_resets_forces_and_moments(self) -> None:
        """Test that deepcopy resets forces and moments to None."""
        original = self.basic_airplane
        original.forces_W = np.array([1.0, 2.0, 3.0])
        original.moments_W_CgP1 = np.array([0.1, 0.2, 0.3])
        original.forceCoefficients_W = np.array([0.01, 0.02, 0.03])
        original.momentCoefficients_W_CgP1 = np.array([0.001, 0.002, 0.003])
        original.forces_G = np.array([1.0, 2.0, 3.0])
        original.forceCoefficients_G = np.array([0.01, 0.02, 0.03])
        original.moments_G_Cg = np.array([0.1, 0.2, 0.3])
        original.momentCoefficients_G_Cg = np.array([0.001, 0.002, 0.003])
        original.moments_W_Cg = np.array([0.1, 0.2, 0.3])
        original.momentCoefficients_W_Cg = np.array([0.001, 0.002, 0.003])

        copied = copy.deepcopy(original)

        self.assertIsNone(copied.forces_W)
        self.assertIsNone(copied.moments_W_CgP1)
        self.assertIsNone(copied.forceCoefficients_W)
        self.assertIsNone(copied.momentCoefficients_W_CgP1)
        self.assertIsNone(copied.forces_G)
        self.assertIsNone(copied.forceCoefficients_G)
        self.assertIsNone(copied.moments_G_Cg)
        self.assertIsNone(copied.momentCoefficients_G_Cg)
        self.assertIsNone(copied.moments_W_Cg)
        self.assertIsNone(copied.momentCoefficients_W_Cg)

    def test_deepcopy_preserves_num_panels(self) -> None:
        """Test that deepcopy preserves the total number of Panels."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(copied.num_panels, original.num_panels)

    def test_deepcopy_preserves_wing_panels(self) -> None:
        """Test that deepcopy preserves Wing Panels."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        for orig_wing, copied_wing in zip(original.wings, copied.wings):
            self.assertIsNotNone(copied_wing.panels)
            assert orig_wing.panels is not None
            assert copied_wing.panels is not None
            self.assertEqual(copied_wing.panels.shape, orig_wing.panels.shape)

            for i in range(orig_wing.panels.shape[0]):
                for j in range(orig_wing.panels.shape[1]):
                    orig_panel = orig_wing.panels[i, j]
                    copied_panel = copied_wing.panels[i, j]
                    self.assertIsNot(orig_panel, copied_panel)
                    npt.assert_array_equal(copied_panel.Frpp_G_Cg, orig_panel.Frpp_G_Cg)

    def test_deepcopy_resets_wing_wake_state(self) -> None:
        """Test that deepcopy resets wake state in all Wings."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        for copied_wing in copied.wings:
            assert copied_wing.gridWrvp_GP1_CgP1 is not None
            self.assertEqual(copied_wing.gridWrvp_GP1_CgP1.shape[0], 0)

    def test_deepcopy_copied_airplane_is_functional(self) -> None:
        """Test that copied Airplanes are fully functional."""
        original = self.basic_airplane
        copied = copy.deepcopy(original)

        num_panels = copied.num_panels
        s_ref = copied.s_ref
        c_ref = copied.c_ref
        b_ref = copied.b_ref

        self.assertGreater(num_panels, 0)
        self.assertGreater(s_ref, 0.0)
        self.assertGreater(c_ref, 0.0)
        self.assertGreater(b_ref, 0.0)

    def test_deepcopy_first_airplane(self) -> None:
        """Test that deepcopy works correctly for first Airplane (Cg at origin)."""
        original = self.first_airplane
        copied = copy.deepcopy(original)

        npt.assert_array_equal(copied.Cg_GP1_CgP1, np.array([0.0, 0.0, 0.0]))
        copied.validate_first_airplane_constraints()


class TestAirplaneImmutability(unittest.TestCase):
    """Tests for Airplane attribute immutability."""

    def setUp(self) -> None:
        """Set up test fixtures for immutability tests."""
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()

    def test_immutable_wings_property(self) -> None:
        """Test that wings property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "wings", ())

    def test_immutable_wings_tuple_cannot_be_modified(self) -> None:
        """Test that wings tuple elements cannot be reassigned."""
        wings: Any = self.basic_airplane.wings
        with self.assertRaises(TypeError):
            wings[0] = None

    def test_immutable_name_property(self) -> None:
        """Test that name property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "name", "New Name")

    def test_immutable_Cg_GP1_CgP1_property(self) -> None:
        """Test that Cg_GP1_CgP1 property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "Cg_GP1_CgP1", np.array([1.0, 2.0, 3.0]))

    def test_immutable_Cg_GP1_CgP1_array_read_only(self) -> None:
        """Test that Cg_GP1_CgP1 array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.basic_airplane.Cg_GP1_CgP1[0] = 999.0

    def test_immutable_weight_property(self) -> None:
        """Test that weight property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "weight", 5000.0)

    def test_immutable_s_ref_property(self) -> None:
        """Test that s_ref property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "s_ref", 100.0)

    def test_immutable_c_ref_property(self) -> None:
        """Test that c_ref property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "c_ref", 5.0)

    def test_immutable_b_ref_property(self) -> None:
        """Test that b_ref property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_airplane, "b_ref", 20.0)

    def test_mutable_forces_W(self) -> None:
        """Test that forces_W attribute remains mutable."""
        self.basic_airplane.forces_W = np.array([1.0, 2.0, 3.0])
        npt.assert_array_equal(self.basic_airplane.forces_W, np.array([1.0, 2.0, 3.0]))

    def test_mutable_forceCoefficients_W(self) -> None:
        """Test that forceCoefficients_W attribute remains mutable."""
        self.basic_airplane.forceCoefficients_W = np.array([0.01, 0.02, 0.03])
        npt.assert_array_equal(
            self.basic_airplane.forceCoefficients_W, np.array([0.01, 0.02, 0.03])
        )

    def test_mutable_moments_W_CgP1(self) -> None:
        """Test that moments_W_CgP1 attribute remains mutable."""
        self.basic_airplane.moments_W_CgP1 = np.array([0.1, 0.2, 0.3])
        npt.assert_array_equal(
            self.basic_airplane.moments_W_CgP1, np.array([0.1, 0.2, 0.3])
        )

    def test_mutable_momentCoefficients_W_CgP1(self) -> None:
        """Test that momentCoefficients_W_CgP1 attribute remains mutable."""
        self.basic_airplane.momentCoefficients_W_CgP1 = np.array([0.001, 0.002, 0.003])
        npt.assert_array_equal(
            self.basic_airplane.momentCoefficients_W_CgP1,
            np.array([0.001, 0.002, 0.003]),
        )

    def test_mutable_forces_G(self) -> None:
        """Test that forces_G attribute remains mutable."""
        self.basic_airplane.forces_G = np.array([1.0, 2.0, 3.0])
        npt.assert_array_equal(self.basic_airplane.forces_G, np.array([1.0, 2.0, 3.0]))

    def test_mutable_forceCoefficients_G(self) -> None:
        """Test that forceCoefficients_G attribute remains mutable."""
        self.basic_airplane.forceCoefficients_G = np.array([0.01, 0.02, 0.03])
        npt.assert_array_equal(
            self.basic_airplane.forceCoefficients_G, np.array([0.01, 0.02, 0.03])
        )

    def test_mutable_moments_G_Cg(self) -> None:
        """Test that moments_G_Cg attribute remains mutable."""
        self.basic_airplane.moments_G_Cg = np.array([0.1, 0.2, 0.3])
        npt.assert_array_equal(
            self.basic_airplane.moments_G_Cg, np.array([0.1, 0.2, 0.3])
        )

    def test_mutable_momentCoefficients_G_Cg(self) -> None:
        """Test that momentCoefficients_G_Cg attribute remains mutable."""
        self.basic_airplane.momentCoefficients_G_Cg = np.array([0.001, 0.002, 0.003])
        npt.assert_array_equal(
            self.basic_airplane.momentCoefficients_G_Cg,
            np.array([0.001, 0.002, 0.003]),
        )

    def test_mutable_moments_W_Cg(self) -> None:
        """Test that moments_W_Cg attribute remains mutable."""
        self.basic_airplane.moments_W_Cg = np.array([0.1, 0.2, 0.3])
        npt.assert_array_equal(
            self.basic_airplane.moments_W_Cg, np.array([0.1, 0.2, 0.3])
        )

    def test_mutable_momentCoefficients_W_Cg(self) -> None:
        """Test that momentCoefficients_W_Cg attribute remains mutable."""
        self.basic_airplane.momentCoefficients_W_Cg = np.array([0.001, 0.002, 0.003])
        npt.assert_array_equal(
            self.basic_airplane.momentCoefficients_W_Cg,
            np.array([0.001, 0.002, 0.003]),
        )


class TestAirplaneDeepCopyWithCgGP1CgP1(unittest.TestCase):
    """Tests for Airplane.deep_copy_with_Cg_GP1_CgP1 method."""

    def setUp(self) -> None:
        """Set up test fixtures for deep_copy_with_Cg_GP1_CgP1 tests."""
        self.first_airplane = geometry_fixtures.make_first_airplane_fixture()
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()
        self.multi_wing_airplane = geometry_fixtures.make_multi_wing_airplane_fixture()

    def test_deep_copy_with_Cg_GP1_CgP1_creates_new_instance(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 creates a new Airplane instance."""
        original = self.first_airplane
        new_position = [5.0, 2.0, -1.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        self.assertIsInstance(copied, ps.geometry.airplane.Airplane)
        self.assertIsNot(original, copied)

    def test_deep_copy_with_Cg_GP1_CgP1_uses_new_position(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 uses the specified position."""
        original = self.first_airplane
        new_position = [5.0, 2.0, -1.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        npt.assert_array_equal(copied.Cg_GP1_CgP1, new_position)
        # Original should be unchanged.
        npt.assert_array_equal(original.Cg_GP1_CgP1, np.array([0.0, 0.0, 0.0]))

    def test_deep_copy_with_Cg_GP1_CgP1_preserves_other_parameters(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 preserves all other parameters."""
        original = self.basic_airplane
        new_position = [10.0, -5.0, 3.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.weight, original.weight)
        self.assertEqual(copied.s_ref, original.s_ref)
        self.assertEqual(copied.c_ref, original.c_ref)
        self.assertEqual(copied.b_ref, original.b_ref)

    def test_deep_copy_with_Cg_GP1_CgP1_creates_independent_wings(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 creates independent Wing copies."""
        original = self.basic_airplane
        new_position = [10.0, -5.0, 3.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        self.assertEqual(len(copied.wings), len(original.wings))
        for orig_wing, copied_wing in zip(original.wings, copied.wings):
            self.assertIsNot(orig_wing, copied_wing)
            self.assertEqual(copied_wing.name, orig_wing.name)

    def test_deep_copy_with_Cg_GP1_CgP1_preserves_num_panels(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 preserves num_panels."""
        original = self.basic_airplane
        new_position = [10.0, -5.0, 3.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        self.assertEqual(copied.num_panels, original.num_panels)

    def test_deep_copy_with_Cg_GP1_CgP1_resets_forces_and_moments(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 resets forces and moments to None."""
        original = self.basic_airplane
        original.forces_W = np.array([1.0, 2.0, 3.0])
        original.moments_W_CgP1 = np.array([0.1, 0.2, 0.3])
        original.forceCoefficients_W = np.array([0.01, 0.02, 0.03])
        original.momentCoefficients_W_CgP1 = np.array([0.001, 0.002, 0.003])
        original.forces_G = np.array([1.0, 2.0, 3.0])
        original.forceCoefficients_G = np.array([0.01, 0.02, 0.03])
        original.moments_G_Cg = np.array([0.1, 0.2, 0.3])
        original.momentCoefficients_G_Cg = np.array([0.001, 0.002, 0.003])
        original.moments_W_Cg = np.array([0.1, 0.2, 0.3])
        original.momentCoefficients_W_Cg = np.array([0.001, 0.002, 0.003])

        new_position = [10.0, -5.0, 3.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        self.assertIsNone(copied.forces_W)
        self.assertIsNone(copied.moments_W_CgP1)
        self.assertIsNone(copied.forceCoefficients_W)
        self.assertIsNone(copied.momentCoefficients_W_CgP1)
        self.assertIsNone(copied.forces_G)
        self.assertIsNone(copied.forceCoefficients_G)
        self.assertIsNone(copied.moments_G_Cg)
        self.assertIsNone(copied.momentCoefficients_G_Cg)
        self.assertIsNone(copied.moments_W_Cg)
        self.assertIsNone(copied.momentCoefficients_W_Cg)

    def test_deep_copy_with_Cg_GP1_CgP1_resets_transformation_cache(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 resets T_pas_G_Cg_to_GP1_CgP1 cache."""
        original = self.basic_airplane
        # Access the transformation matrix to cache it.
        _ = original.T_pas_G_Cg_to_GP1_CgP1

        new_position = [10.0, -5.0, 3.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        # The cached transformation should be None and will be recomputed on access.
        # Access it and verify it uses the new position.
        T = copied.T_pas_G_Cg_to_GP1_CgP1
        # The translation component should match the new position.
        npt.assert_array_almost_equal(T[:3, 3], new_position)

    def test_deep_copy_with_Cg_GP1_CgP1_accepts_various_array_types(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 accepts various array-like inputs."""
        _ = self.first_airplane
        valid_positions: list[np.ndarray | Sequence[float]] = [
            [1.0, 2.0, 3.0],  # list
            (4.0, 5.0, 6.0),  # tuple
            np.array([7.0, 8.0, 9.0]),  # numpy array
            [1, 2, 3],  # list of ints
        ]

        for position in valid_positions:
            with self.subTest(position=position):
                # Create fresh fixture since Wings can only be processed once.
                test_airplane = geometry_fixtures.make_first_airplane_fixture()
                copied = test_airplane.deep_copy_with_Cg_GP1_CgP1(position)
                npt.assert_array_equal(copied.Cg_GP1_CgP1, position)

    def test_deep_copy_with_Cg_GP1_CgP1_validates_position(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 validates the position parameter."""
        original = self.first_airplane

        # Test invalid position shapes.
        invalid_positions: list[Any] = [
            [1.0, 2.0],  # Wrong size.
            [1.0, 2.0, 3.0, 4.0],  # Wrong size.
            "not a vector",  # String.
        ]

        for invalid_position in invalid_positions:
            with self.subTest(invalid_position=invalid_position):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    original.deep_copy_with_Cg_GP1_CgP1(invalid_position)

    def test_deep_copy_with_Cg_GP1_CgP1_multi_wing(self) -> None:
        """Test that deep_copy_with_Cg_GP1_CgP1 works correctly for multi wing
        Airplanes."""
        original = self.multi_wing_airplane
        new_position = [15.0, -10.0, 5.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        self.assertEqual(len(copied.wings), len(original.wings))
        npt.assert_array_equal(copied.Cg_GP1_CgP1, new_position)

        for i, (orig_wing, copied_wing) in enumerate(zip(original.wings, copied.wings)):
            # noinspection PyUnresolvedReferences
            with self.subTest(wing_index=i):
                self.assertIsNot(orig_wing, copied_wing)
                self.assertEqual(copied_wing.symmetry_type, orig_wing.symmetry_type)

    def test_deep_copy_with_Cg_GP1_CgP1_new_position_array_is_read_only(self) -> None:
        """Test that the new position array is read only after copying."""
        original = self.first_airplane
        new_position = [5.0, 2.0, -1.0]
        copied = original.deep_copy_with_Cg_GP1_CgP1(new_position)

        with self.assertRaises(ValueError):
            copied.Cg_GP1_CgP1[0] = 999.0


class TestAirplaneTPasGCgToGP1CgP1(unittest.TestCase):
    """Tests for Airplane.T_pas_G_Cg_to_GP1_CgP1 property."""

    def setUp(self) -> None:
        """Set up test fixtures for T_pas_G_Cg_to_GP1_CgP1 tests."""
        self.first_airplane = geometry_fixtures.make_first_airplane_fixture()
        self.follower_airplane = geometry_fixtures.make_follower_airplane_fixture()
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()

    def test_first_airplane_returns_identity(self) -> None:
        """Test that first Airplane (Cg at origin) returns identity transformation."""
        T = self.first_airplane.T_pas_G_Cg_to_GP1_CgP1
        expected_identity = np.eye(4, dtype=float)
        npt.assert_array_almost_equal(T, expected_identity)

    def test_follower_airplane_returns_translation_matrix(self) -> None:
        """Test that follower Airplane returns correct translation transformation."""
        T = self.follower_airplane.T_pas_G_Cg_to_GP1_CgP1
        position = self.follower_airplane.Cg_GP1_CgP1

        # For a passive translation, the position appears in the last column.
        npt.assert_array_almost_equal(T[:3, 3], position)

        # The rotation part should be identity.
        npt.assert_array_almost_equal(T[:3, :3], np.eye(3, dtype=float))

    def test_transformation_matrix_shape(self) -> None:
        """Test that the transformation matrix has the correct shape."""
        T = self.basic_airplane.T_pas_G_Cg_to_GP1_CgP1

        self.assertEqual(T.shape, (4, 4))

    def test_transformation_matrix_is_ndarray(self) -> None:
        """Test that the transformation matrix is a numpy ndarray."""
        T = self.basic_airplane.T_pas_G_Cg_to_GP1_CgP1

        self.assertIsInstance(T, np.ndarray)

    def test_transformation_matrix_is_read_only(self) -> None:
        """Test that the transformation matrix cannot be modified in place."""
        T = self.basic_airplane.T_pas_G_Cg_to_GP1_CgP1

        with self.assertRaises(ValueError):
            T[0, 0] = 999.0

    def test_transformation_matrix_is_cached(self) -> None:
        """Test that the transformation matrix is cached after first access."""
        T1 = self.basic_airplane.T_pas_G_Cg_to_GP1_CgP1
        T2 = self.basic_airplane.T_pas_G_Cg_to_GP1_CgP1

        self.assertIs(T1, T2)

    def test_transformation_matrix_homogeneous_last_row(self) -> None:
        """Test that the transformation matrix has correct homogeneous last row."""
        T = self.basic_airplane.T_pas_G_Cg_to_GP1_CgP1

        expected_last_row = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T[3, :], expected_last_row)

    def test_transformation_matrix_applies_correct_translation(self) -> None:
        """Test that the transformation matrix applies the correct translation to a
        point."""
        position = np.array([5.0, 2.0, -1.0])
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        airplane = ps.geometry.airplane.Airplane(
            wings=[test_wing], Cg_GP1_CgP1=position
        )

        T = airplane.T_pas_G_Cg_to_GP1_CgP1

        # Transform the origin point in local coordinates.
        local_point_homogeneous = np.array([0.0, 0.0, 0.0, 1.0])
        global_point_homogeneous = T @ local_point_homogeneous

        # The global point should be at the Airplane's CG position.
        npt.assert_array_almost_equal(global_point_homogeneous[:3], position)


class TestAirplaneGetPlottableData(unittest.TestCase):
    """Tests for Airplane.get_plottable_data method."""

    def setUp(self) -> None:
        """Set up test fixtures for get_plottable_data tests."""
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()
        self.multi_wing_airplane = geometry_fixtures.make_multi_wing_airplane_fixture()

    def test_get_plottable_data_returns_list_when_show_is_false(self) -> None:
        """Test that get_plottable_data returns a list when show is False."""
        result = self.basic_airplane.get_plottable_data(show=False)

        self.assertIsInstance(result, list)

    def test_get_plottable_data_returns_two_sub_lists(self) -> None:
        """Test that get_plottable_data returns two sub lists (outlines and MCLs)."""
        result = self.basic_airplane.get_plottable_data(show=False)

        assert result is not None
        self.assertEqual(len(result), 2)

    def test_get_plottable_data_structure_matches_wings_and_cross_sections(
        self,
    ) -> None:
        """Test that the returned data structure matches the number of Wings and
        WingCrossSections."""
        result = self.basic_airplane.get_plottable_data(show=False)
        assert result is not None
        outlines = result[0]
        mcls = result[1]

        # The number of sub lists should match the number of Wings.
        self.assertEqual(len(outlines), len(self.basic_airplane.wings))
        self.assertEqual(len(mcls), len(self.basic_airplane.wings))

        # Each Wing's sub list should have the same number of cross sections.
        for wing_id, wing in enumerate(self.basic_airplane.wings):
            expected_num = len(wing.wing_cross_sections)
            self.assertEqual(len(outlines[wing_id]), expected_num)
            self.assertEqual(len(mcls[wing_id]), expected_num)

    def test_get_plottable_data_returns_ndarrays(self) -> None:
        """Test that get_plottable_data returns ndarrays for each cross section."""
        result = self.basic_airplane.get_plottable_data(show=False)
        assert result is not None
        outlines = result[0]
        mcls = result[1]

        for wing_outlines in outlines:
            for outline in wing_outlines:
                self.assertIsInstance(outline, np.ndarray)

        for wing_mcls in mcls:
            for mcl in wing_mcls:
                self.assertIsInstance(mcl, np.ndarray)

    def test_get_plottable_data_returns_3d_points(self) -> None:
        """Test that get_plottable_data returns arrays with 3 columns (x, y, z)."""
        result = self.basic_airplane.get_plottable_data(show=False)
        assert result is not None
        outlines = result[0]
        mcls = result[1]

        for wing_outlines in outlines:
            for outline in wing_outlines:
                self.assertEqual(outline.shape[1], 3)

        for wing_mcls in mcls:
            for mcl in wing_mcls:
                self.assertEqual(mcl.shape[1], 3)

    def test_get_plottable_data_default_show_is_false(self) -> None:
        """Test that get_plottable_data default for show is False."""
        result = self.basic_airplane.get_plottable_data()

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_accepts_numpy_bool(self) -> None:
        """Test that get_plottable_data accepts numpy bool for show parameter."""
        result = self.basic_airplane.get_plottable_data(show=np.bool_(False))

        self.assertIsNotNone(result)
        self.assertIsInstance(result, list)

    def test_get_plottable_data_multi_wing_airplane(self) -> None:
        """Test get_plottable_data with a multi wing Airplane."""
        result = self.multi_wing_airplane.get_plottable_data(show=False)
        assert result is not None
        outlines = result[0]
        mcls = result[1]

        # Should have data for each Wing.
        self.assertEqual(len(outlines), len(self.multi_wing_airplane.wings))
        self.assertEqual(len(mcls), len(self.multi_wing_airplane.wings))

    def test_get_plottable_data_invalid_show_type_raises(self) -> None:
        """Test that get_plottable_data raises error for invalid show type."""
        bad_show: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            self.basic_airplane.get_plottable_data(show=bad_show)


class TestAirplaneDraw(unittest.TestCase):
    """Tests for Airplane.draw method."""

    def setUp(self) -> None:
        """Set up test fixtures for draw tests."""
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()

    def test_draw_runs_without_error_in_testing_mode(self) -> None:
        """Test that draw runs without error in testing mode."""
        # Use testing=True to avoid blocking on window close.
        try:
            self.basic_airplane.draw(save=False, testing=True)
        except Exception as e:
            self.fail(f"draw() raised {type(e).__name__}: {e}")

    def test_draw_accepts_numpy_bool_for_save(self) -> None:
        """Test that draw accepts numpy bool for save parameter."""
        try:
            self.basic_airplane.draw(save=np.bool_(False), testing=np.bool_(True))
        except Exception as e:
            self.fail(f"draw() raised {type(e).__name__}: {e}")

    def test_draw_invalid_save_type_raises(self) -> None:
        """Test that draw raises error for invalid save type."""
        bad_save: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            self.basic_airplane.draw(save=bad_save, testing=True)

    def test_draw_invalid_testing_type_raises(self) -> None:
        """Test that draw raises error for invalid testing type."""
        bad_testing: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            self.basic_airplane.draw(save=False, testing=bad_testing)
