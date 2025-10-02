"""This module contains a class to test Airplanes."""

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

    def test_initialization_valid_parameters(self):
        """Test Airplane initialization with valid parameters."""
        # Test that basic Airplane initializes correctly
        self.assertIsInstance(self.basic_airplane, ps.geometry.airplane.Airplane)
        self.assertIsInstance(self.basic_airplane.wings, list)
        self.assertEqual(len(self.basic_airplane.wings), 1)
        self.assertEqual(self.basic_airplane.name, "Basic Test Airplane")
        npt.assert_array_equal(self.basic_airplane.Cgi_E_I, np.array([1.0, 0.5, -0.2]))
        npt.assert_array_equal(
            self.basic_airplane.angles_E_to_B_izyx, np.array([10.0, -5.0, 15.0])
        )
        self.assertEqual(self.basic_airplane.weight, 1000.0)

    def test_wings_parameter_validation(self):
        """Test that wings parameter validation works correctly."""
        # Test empty list raises error
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(wings=[])

        # Test non-list raises error
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings="not a list")

        # Test non-Wing objects raise error
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=["not a wing"])

        # Test mixed valid and invalid Wings
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1, "invalid"])

    def test_name_parameter_validation(self):
        """Test name parameter validation."""
        # Test valid string name
        airplane = ps.geometry.airplane.Airplane(
            wings=[self.test_wing_type_1], name="Valid Test Name"
        )
        self.assertEqual(airplane.name, "Valid Test Name")

        # Test invalid name types
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1], name=123)

        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1], name=None)

    def test_Cgi_E_I_parameter_validation(self):
        """Test Cgi_E_I parameter validation."""
        # Test valid 3D vectors
        valid_positions = [
            [0.0, 0.0, 0.0],
            [1.0, -2.0, 3.5],
            np.array([0.5, 1.5, -0.8]),
            (2.0, -1.0, 0.0),
        ]

        for position in valid_positions:
            with self.subTest(position=position):
                airplane = ps.geometry.airplane.Airplane(
                    wings=[self.test_wing_type_1], Cgi_E_I=position
                )
                npt.assert_array_equal(airplane.Cgi_E_I, position)

        # Test invalid positions
        invalid_positions = [
            [1.0, 2.0],  # Wrong size
            [1.0, 2.0, 3.0, 4.0],  # Wrong size
            "not a vector",  # String
            None,  # None
        ]

        for invalid_position in invalid_positions:
            with self.subTest(invalid_position=invalid_position):
                with self.assertRaises((ValueError, TypeError)):
                    ps.geometry.airplane.Airplane(
                        wings=[self.test_wing_type_1], Cgi_E_I=invalid_position
                    )

    def test_angles_E_to_B_izyx_parameter_validation(self):
        """Test angles_E_to_B_izyx parameter validation."""
        # Test valid angle vectors (within range (-180, 180])
        valid_angles = [
            [0.0, 0.0, 0.0],
            [45.0, -30.0, 90.0],
            [179.9, -179.9, 180.0],
            np.array([15.0, -45.0, 60.0]),
            (-90.0, 45.0, -120.0),
        ]

        for angles in valid_angles:
            with self.subTest(angles=angles):
                airplane = ps.geometry.airplane.Airplane(
                    wings=[self.test_wing_type_1], angles_E_to_B_izyx=angles
                )
                npt.assert_array_equal(airplane.angles_E_to_B_izyx, angles)

        # Test angles outside valid range
        invalid_angles = [
            [180.1, 0.0, 0.0],  # Greater than 180
            [0.0, -180.1, 0.0],  # Less than -180
            [200.0, 0.0, 0.0],  # Much greater than 180
            [0.0, 0.0, -200.0],  # Much less than -180
        ]

        for invalid_angle in invalid_angles:
            with self.subTest(invalid_angle=invalid_angle):
                with self.assertRaises(ValueError):
                    ps.geometry.airplane.Airplane(
                        wings=[self.test_wing_type_1], angles_E_to_B_izyx=invalid_angle
                    )

    def test_weight_parameter_validation(self):
        """Test weight parameter validation."""
        # Test valid weights (non-negative)
        valid_weights = [0.0, 0.1, 100.0, 5000.0, 10000.0]

        for weight in valid_weights:
            with self.subTest(weight=weight):
                airplane = ps.geometry.airplane.Airplane(
                    wings=[self.test_wing_type_1], weight=weight
                )
                self.assertEqual(airplane.weight, weight)

        # Test invalid weights (negative)
        invalid_weights = [-0.1, -100.0, -1000.0]

        for invalid_weight in invalid_weights:
            with self.subTest(invalid_weight=invalid_weight):
                with self.assertRaises(ValueError):
                    ps.geometry.airplane.Airplane(
                        wings=[self.test_wing_type_1], weight=invalid_weight
                    )

        # Test invalid weight types
        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(wings=[self.test_wing_type_1], weight="heavy")

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
        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(
                wings=[self.test_wing_type_1], s_ref=-1.0  # Negative value
            )

        with self.assertRaises(ValueError):
            ps.geometry.airplane.Airplane(
                wings=[self.test_wing_type_1], c_ref=0.0  # Zero value
            )

        with self.assertRaises(TypeError):
            ps.geometry.airplane.Airplane(
                wings=[self.test_wing_type_1], b_ref="large"  # String
            )

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
        self.assertIsNone(airplane.moments_W_Cg)
        self.assertIsNone(airplane.momentCoefficients_W_Cg)

    def test_validate_first_airplane_constraints_valid(self):
        """Test validate_first_airplane_constraints with valid first Airplane."""
        # First Airplane should pass validation (Cgi_E_I is all zeros)
        try:
            self.first_airplane.validate_first_airplane_constraints()
        except Exception as e:
            self.fail(f"First airplane validation failed unexpectedly: {e}")

    def test_validate_first_airplane_constraints_invalid(self):
        """Test validate_first_airplane_constraints with invalid Airplane."""
        # Basic Airplane should fail validation (Cgi_E_I is not all zeros)
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
        self.assertIsNone(first_wing.symmetry_normal_Wn)
        self.assertIsNone(first_wing.symmetry_point_Wn_Ler)

        # Second Wing should be reflected type 3
        second_wing = result[1]
        self.assertEqual(second_wing.symmetry_type, 3)
        self.assertFalse(second_wing.symmetric)
        self.assertTrue(second_wing.mirror_only)
        self.assertTrue(second_wing.name.startswith("Reflected"))

    def test_process_wing_symmetry_control_surface_validation_types_1_2_3(self):
        """Test control surface validation for symmetry types 1, 2, 3."""
        # Create Wings with invalid control surface configurations
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            geometry_fixtures.make_basic_wing_cross_section_fixture(),  # Has control surface
        ]

        # Type 1: should fail with control surfaces
        with self.assertRaises(ValueError):
            wing_type_1 = ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                symmetric=False,
                mirror_only=False,
            )
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing_type_1)

        # Type 2: should fail with control surfaces
        with self.assertRaises(ValueError):
            wing_type_2 = ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                symmetric=False,
                mirror_only=True,
                symmetry_normal_Wn=[0.0, 1.0, 0.0],
                symmetry_point_Wn_Ler=[0.0, 0.0, 0.0],
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
                symmetry_normal_Wn=[0.0, 1.0, 0.0],
                symmetry_point_Wn_Ler=[0.0, 0.0, 0.0],
            )
            ps.geometry.airplane.Airplane.process_wing_symmetry(wing_type_4)

    def test_process_wing_symmetry_type_5_control_surface_deflections(self):
        """Test type 5 Wing processing with different control surface deflections."""
        # Create asymmetric control surface Wing cross section
        asymmetric_wcs = (
            geometry_fixtures.make_asymmetric_control_surface_wing_cross_section_fixture()
        )

        # Create type 5 Wing with asymmetric control surfaces
        wing_cross_sections = [
            geometry_fixtures.make_root_wing_cross_section_fixture(),
            asymmetric_wcs,
        ]
        wing_cross_sections[0].control_surface_symmetry_type = "asymmetric"
        wing_cross_sections[0].control_surface_deflection = 2.5

        wing = ps.geometry.wing.Wing(
            wing_cross_sections=wing_cross_sections,
            symmetric=True,
            mirror_only=False,
            symmetry_normal_Wn=[0.0, 0.707, 0.707],  # Non-coincident plane
            symmetry_point_Wn_Ler=[0.5, 0.0, 0.0],
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

    def test_comprehensive_airplane_properties(self):
        """Test comprehensive property access for different Airplane configurations."""
        airplanes_to_test = [
            (self.basic_airplane, "basic"),
            (self.first_airplane, "first"),
            (self.multi_wing_airplane, "multi-wing"),
            (self.custom_reference_airplane, "custom reference"),
        ]

        for airplane, airplane_type in airplanes_to_test:
            with self.subTest(airplane_type=airplane_type):
                # Test basic properties
                self.assertIsInstance(airplane.wings, list)
                self.assertGreater(len(airplane.wings), 0)
                self.assertIsInstance(airplane.name, str)
                self.assertEqual(len(airplane.Cgi_E_I), 3)
                self.assertEqual(len(airplane.angles_E_to_B_izyx), 3)
                self.assertGreaterEqual(airplane.weight, 0.0)

                # Test reference dimensions
                self.assertGreater(airplane.s_ref, 0.0)
                self.assertGreater(airplane.c_ref, 0.0)
                self.assertGreater(airplane.b_ref, 0.0)

                # Test calculated properties
                self.assertGreaterEqual(airplane.num_panels, 0)

    def test_edge_case_angle_boundaries(self):
        """Test angle validation at boundary values."""
        # Test boundary values that should be valid
        valid_boundary_angles = [
            [180.0, 0.0, 0.0],  # Exactly 180 should be valid
            [0.0, -179.999, 0.0],  # Just above -180 should be valid
            [-179.999, 179.999, 180.0],  # Multiple boundary values
        ]

        for angles in valid_boundary_angles:
            with self.subTest(angles=angles):
                airplane = ps.geometry.airplane.Airplane(
                    wings=[self.test_wing_type_1], angles_E_to_B_izyx=angles
                )
                npt.assert_allclose(airplane.angles_E_to_B_izyx, angles, atol=1e-10)

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

    # TODO: Finalize Airplane's get_plottable_data testing.
    # def test_airplane_get_plottable_data(self):
    #     """Test that the get_plottable_data method works correctly."""
    #     airplane = geometry_fixtures.make_basic_airplane_fixture()
    #     # airplane = geometry_fixtures.make_type_5_wing_airplane_fixture()
    #     # airplane = geometry_fixtures.make_multi_wing_airplane_fixture()
    #     # airplane.get_plottable_data(show=True)
    #     airplane.wings[0].get_plottable_data(show=True)
    #     # airplane.wings[0].wing_cross_sections[0].get_plottable_data(show=True)
    #     # airplane.draw()


if __name__ == "__main__":
    unittest.main()
