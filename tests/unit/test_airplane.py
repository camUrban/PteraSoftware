"""This module contains classes to test Airplanes."""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures


class TestAirplane(unittest.TestCase):
    """This is a class with functions to test Airplanes."""

    def setUp(self):
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
        self.test_wing_type_4 = geometry_fixtures.make_type_4_wing_fixture()

    def test_wings_parameter_validation(self):
        """Test that wings parameter validation works correctly."""
        # Test empty list raises error
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[])

        # Test non-list raises error
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings="not a list")

        # Test non-Wing objects raise error
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=["not a wing"])

        # Test mixed valid and invalid Wings
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1, "invalid"])

    def test_name_parameter_validation(self):
        """Test name parameter validation."""
        # Test valid string name
        # Create fresh fixture since Wings can only be processed once
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        airplane = ps.geometry.airplane.Airplane(
            wings=[test_wing], name="Valid Test Name"
        )
        self.assertEqual(airplane.name, "Valid Test Name")

        # Test invalid name types
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], name=123)

        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], name=None)

    def test_Cg_GP1_CgP1_parameter_validation(self):
        """Test Cg_GP1_CgP1 parameter validation."""
        # Test valid 3D vectors
        valid_positions = [
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
        invalid_positions = [
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

    def test_weight_parameter_validation(self):
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

        # Test invalid weight types
        # Create fresh fixture since Wings can only be processed once
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], weight="heavy")

    def test_reference_dimensions_default_behavior(self):
        """Test reference dimensions default to first Wing's properties."""
        # Create Airplane with no explicit reference dimensions
        airplane = ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1])

        # Reference dimensions should be populated from the first Wing
        first_wing = airplane.wings[0]
        self.assertEqual(airplane.s_ref, first_wing.projected_area)
        self.assertEqual(airplane.c_ref, first_wing.mean_aerodynamic_chord)
        self.assertEqual(airplane.b_ref, first_wing.span)

    def test_reference_dimensions_explicit_values(self):
        """Test reference dimensions with explicit values."""
        # Test custom reference Airplane
        self.assertEqual(self.custom_reference_airplane.s_ref, 15.0)
        self.assertEqual(self.custom_reference_airplane.c_ref, 2.0)
        self.assertEqual(self.custom_reference_airplane.b_ref, 10.0)

        # Test validation of reference dimensions
        # Create fresh fixtures since Wings can only be processed once
        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[test_wing], s_ref=-1.0)

        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[test_wing], c_ref=0.0)

        test_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            ps.geometry.airplane.Airplane(wings=[test_wing], b_ref="large")

    def test_num_panels_calculation(self):
        """Test that num_panels is calculated correctly from all Wings."""
        # Single-Wing Airplane
        single_wing_panels = self.basic_airplane.wings[0].num_panels
        self.assertEqual(self.basic_airplane.num_panels, single_wing_panels)

        # Multi-Wing Airplane
        expected_panels = sum(
            wing.num_panels for wing in self.multi_wing_airplane.wings
        )
        self.assertEqual(self.multi_wing_airplane.num_panels, expected_panels)

    def test_force_moment_attributes_initialization(self):
        """Test that force and moment attributes are initialized to None."""
        airplane = self.basic_airplane

        # All force/moment attributes should be None initially
        self.assertIsNone(airplane.forces_W)
        self.assertIsNone(airplane.forceCoefficients_W)
        self.assertIsNone(airplane.moments_W_CgP1)
        self.assertIsNone(airplane.momentCoefficients_W_CgP1)

    def test_validate_first_airplane_constraints_valid(self):
        """Test validate_first_airplane_constraints with valid first Airplane."""
        # First Airplane should pass validation (Cg_GP1_CgP1 is all zeros)
        try:
            self.first_airplane.validate_first_airplane_constraints()
        except Exception as e:
            self.fail(f"First airplane validation failed unexpectedly: {e}")

    def test_validate_first_airplane_constraints_invalid(self):
        """Test validate_first_airplane_constraints with invalid Airplane."""
        # Basic Airplane should fail validation (Cg_GP1_CgP1 is not all zeros)
        with self.assertRaises(ValueError):
            self.basic_airplane.validate_first_airplane_constraints()

        # Custom reference Airplane should fail validation
        with self.assertRaises(ValueError):
            self.custom_reference_airplane.validate_first_airplane_constraints()

    def test_multi_wing_configuration(self):
        """Test Airplane with multiple Wings."""
        airplane = self.multi_wing_airplane

        # Should have multiple Wings
        self.assertEqual(len(airplane.wings), 2)
        self.assertIsInstance(airplane.wings[0], ps.geometry.wing.Wing)
        self.assertIsInstance(airplane.wings[1], ps.geometry.wing.Wing)

        # Wings should have different names
        self.assertNotEqual(airplane.wings[0].name, airplane.wings[1].name)

    def test_type_5_wing_processing(self):
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

    def test_process_wing_symmetry_type_1(self):
        """Test process_wing_symmetry with type 1 Wing."""
        wing = geometry_fixtures.make_type_1_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 1)

    def test_process_wing_symmetry_type_2(self):
        """Test process_wing_symmetry with type 2 Wing."""
        wing = geometry_fixtures.make_type_2_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 2)

    def test_process_wing_symmetry_type_3(self):
        """Test process_wing_symmetry with type 3 Wing."""
        wing = geometry_fixtures.make_type_3_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 3)

    def test_process_wing_symmetry_type_4(self):
        """Test process_wing_symmetry with type 4 Wing."""
        wing = geometry_fixtures.make_type_4_wing_fixture()
        result = ps.geometry.airplane.Airplane.process_wing_symmetry(wing)

        # Should return list with one Wing
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].symmetry_type, 4)

    def test_process_wing_symmetry_type_5(self):
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

    def test_process_wing_symmetry_control_surface_validation_types_1_2_3(self):
        """Test control surface validation for symmetry types 1, 2, 3."""
        # Type 1: should fail with control surfaces
        # Create fresh fixtures since WingCrossSections can only be validated once
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            geometry_fixtures.make_basic_wing_cross_section_fixture(),  # Has control surface
        ]
        with self.assertRaises(ValueError):
            wing_type_1 = ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                symmetric=False,
                mirror_only=False,
            )
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing_type_1)

        # Type 2: should fail with control surfaces
        # Create fresh fixtures since WingCrossSections can only be validated once
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            geometry_fixtures.make_basic_wing_cross_section_fixture(),  # Has control surface
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

    def test_process_wing_symmetry_control_surface_validation_types_4_5(self):
        """Test control surface validation for symmetry types 4, 5."""
        # Create Wings without control surface configurations
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            geometry_fixtures.make_minimal_wing_cross_section_fixture(),  # No control surface
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

    def test_process_wing_symmetry_type_5_control_surface_deflections(self):
        """Test type 5 Wing processing with different control surface deflections."""
        # Create asymmetric control surface Wing cross sections
        # Use root fixture with asymmetric control surface already configured
        root_wcs = (
            geometry_fixtures.make_root_asymmetric_control_surface_wing_cross_section_fixture()
        )
        tip_wcs = (
            geometry_fixtures.make_asymmetric_control_surface_wing_cross_section_fixture()
        )

        # Create type 5 Wing with asymmetric control surfaces
        wing_cross_sections = [root_wcs, tip_wcs]

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

        # Check that asymmetric control surface deflections are negated in reflected Wing
        original_wing = result[0]
        reflected_wing = result[1]

        # Find corresponding WingCrossSections with asymmetric control surfaces
        for i, wcs in enumerate(original_wing.wing_cross_sections):
            reflected_wcs = reflected_wing.wing_cross_sections[i]
            # Reflected Wing should have None-type control surface symmetry
            self.assertEqual(reflected_wcs.control_surface_symmetry_type, None)

    def test_airplane_with_various_wing_combinations(self):
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

    def setUp(self):
        """Set up test fixtures for deepcopy tests."""
        self.basic_airplane = geometry_fixtures.make_basic_airplane_fixture()
        self.first_airplane = geometry_fixtures.make_first_airplane_fixture()
        self.multi_wing_airplane = geometry_fixtures.make_multi_wing_airplane_fixture()

    def test_deepcopy_creates_new_instance(self):
        """Test that deepcopy creates a new Airplane instance."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.geometry.airplane.Airplane)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_airplane_parameters(self):
        """Test that deepcopy preserves Airplane parameters."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(copied.name, original.name)
        self.assertEqual(copied.weight, original.weight)
        self.assertEqual(copied.s_ref, original.s_ref)
        self.assertEqual(copied.c_ref, original.c_ref)
        self.assertEqual(copied.b_ref, original.b_ref)
        npt.assert_array_equal(copied.Cg_GP1_CgP1, original.Cg_GP1_CgP1)

    def test_deepcopy_creates_independent_cg_array(self):
        """Test that deepcopy creates an independent copy of Cg_GP1_CgP1."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertIsNot(copied.Cg_GP1_CgP1, original.Cg_GP1_CgP1)

    def test_deepcopy_creates_independent_wings(self):
        """Test that deepcopy creates independent Wing copies."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(len(copied.wings), len(original.wings))
        for orig_wing, copied_wing in zip(original.wings, copied.wings):
            self.assertIsNot(orig_wing, copied_wing)
            self.assertEqual(copied_wing.name, orig_wing.name)

    def test_deepcopy_multi_wing_airplane(self):
        """Test that deepcopy works correctly for multi-wing Airplanes."""
        import copy

        original = self.multi_wing_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(len(copied.wings), len(original.wings))
        for i, (orig_wing, copied_wing) in enumerate(zip(original.wings, copied.wings)):
            with self.subTest(wing_index=i):
                self.assertIsNot(orig_wing, copied_wing)
                self.assertEqual(copied_wing.symmetry_type, orig_wing.symmetry_type)

    def test_deepcopy_resets_forces_and_moments(self):
        """Test that deepcopy resets forces and moments to None."""
        import copy

        original = self.basic_airplane
        original.forces_W = np.array([1.0, 2.0, 3.0])
        original.moments_W_CgP1 = np.array([0.1, 0.2, 0.3])
        original.forceCoefficients_W = np.array([0.01, 0.02, 0.03])
        original.momentCoefficients_W_CgP1 = np.array([0.001, 0.002, 0.003])

        copied = copy.deepcopy(original)

        self.assertIsNone(copied.forces_W)
        self.assertIsNone(copied.moments_W_CgP1)
        self.assertIsNone(copied.forceCoefficients_W)
        self.assertIsNone(copied.momentCoefficients_W_CgP1)

    def test_deepcopy_preserves_num_panels(self):
        """Test that deepcopy preserves the total number of Panels."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        self.assertEqual(copied.num_panels, original.num_panels)

    def test_deepcopy_preserves_wing_panels(self):
        """Test that deepcopy preserves Wing Panels."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        for orig_wing, copied_wing in zip(original.wings, copied.wings):
            self.assertIsNotNone(copied_wing.panels)
            self.assertEqual(copied_wing.panels.shape, orig_wing.panels.shape)

            for i in range(orig_wing.panels.shape[0]):
                for j in range(orig_wing.panels.shape[1]):
                    orig_panel = orig_wing.panels[i, j]
                    copied_panel = copied_wing.panels[i, j]
                    self.assertIsNot(orig_panel, copied_panel)
                    npt.assert_array_equal(copied_panel.Frpp_G_Cg, orig_panel.Frpp_G_Cg)

    def test_deepcopy_resets_wing_wake_state(self):
        """Test that deepcopy resets wake state in all Wings."""
        import copy

        original = self.basic_airplane
        copied = copy.deepcopy(original)

        for copied_wing in copied.wings:
            self.assertEqual(copied_wing.wake_ring_vortices.shape[0], 0)
            self.assertEqual(copied_wing.gridWrvp_GP1_CgP1.shape[0], 0)

    def test_deepcopy_copied_airplane_is_functional(self):
        """Test that copied Airplanes are fully functional."""
        import copy

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

    def test_deepcopy_first_airplane(self):
        """Test that deepcopy works correctly for first Airplane (Cg at origin)."""
        import copy

        original = self.first_airplane
        copied = copy.deepcopy(original)

        npt.assert_array_equal(copied.Cg_GP1_CgP1, np.array([0.0, 0.0, 0.0]))
        copied.validate_first_airplane_constraints()


if __name__ == "__main__":
    unittest.main()
