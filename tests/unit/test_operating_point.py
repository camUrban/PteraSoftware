"""This module contains a class to test OperatingPoints."""

import unittest
import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import operating_point_fixtures


class TestOperatingPoint(unittest.TestCase):
    """This is a class with functions to test OperatingPoints."""

    def setUp(self):
        """Set up test fixtures for OperatingPoint tests."""
        # Create fixtures for various OperatingPoint configurations.
        self.basic_op = operating_point_fixtures.make_basic_operating_point_fixture()
        self.zero_alpha_beta_op = (
            operating_point_fixtures.make_zero_alpha_beta_operating_point_fixture()
        )
        self.high_alpha_op = (
            operating_point_fixtures.make_high_alpha_operating_point_fixture()
        )
        self.negative_alpha_op = (
            operating_point_fixtures.make_negative_alpha_operating_point_fixture()
        )
        self.nonzero_beta_op = (
            operating_point_fixtures.make_nonzero_beta_operating_point_fixture()
        )
        self.high_speed_op = (
            operating_point_fixtures.make_high_speed_operating_point_fixture()
        )
        self.low_density_op = (
            operating_point_fixtures.make_low_density_operating_point_fixture()
        )
        self.with_external_force_op = (
            operating_point_fixtures.make_with_external_force_operating_point_fixture()
        )
        self.custom_viscosity_op = (
            operating_point_fixtures.make_custom_viscosity_operating_point_fixture()
        )
        self.boundary_alpha_op = (
            operating_point_fixtures.make_boundary_alpha_operating_point_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test OperatingPoint initialization with valid parameters."""
        # Test basic OperatingPoint initialization
        self.assertIsInstance(self.basic_op, ps.operating_point.OperatingPoint)
        self.assertEqual(self.basic_op.rho, 1.225)
        self.assertEqual(self.basic_op.vCg__E, 10.0)
        self.assertEqual(self.basic_op.alpha, 5.0)
        self.assertEqual(self.basic_op.beta, 0.0)
        self.assertEqual(self.basic_op.externalFX_W, 0.0)
        self.assertEqual(self.basic_op.nu, 15.06e-6)

    def test_initialization_with_defaults(self):
        """Test that default values are applied correctly."""
        # Create OperatingPoint with all defaults
        op_default = ps.operating_point.OperatingPoint()

        # Verify default values
        self.assertEqual(op_default.rho, 1.225)
        self.assertEqual(op_default.vCg__E, 10.0)
        self.assertEqual(op_default.alpha, 5.0)
        self.assertEqual(op_default.beta, 0.0)
        self.assertEqual(op_default.externalFX_W, 0.0)
        self.assertEqual(op_default.nu, 15.06e-6)

    def test_rho_parameter_validation(self):
        """Test rho parameter validation."""
        # Test valid positive floats
        valid_rho_values = [0.5, 1.225, 2.0, 10.0]

        for rho in valid_rho_values:
            with self.subTest(rho=rho):
                op = ps.operating_point.OperatingPoint(rho=rho)
                self.assertEqual(op.rho, float(rho))

        # Test invalid values (negative, zero, non-numeric)
        invalid_rho_values = [-1.0, 0.0, -0.5, "invalid", None]

        for invalid_rho in invalid_rho_values:
            with self.subTest(invalid_rho=invalid_rho):
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(rho=invalid_rho)

    def test_vCg__E_parameter_validation(self):
        """Test vCg__E parameter validation."""
        # Test valid positive floats
        valid_vCg_values = [0.01, 1.0, 10.0, 100.0, 300.0]

        for vCg__E in valid_vCg_values:
            with self.subTest(vCg__E=vCg__E):
                op = ps.operating_point.OperatingPoint(vCg__E=vCg__E)
                self.assertEqual(op.vCg__E, float(vCg__E))

        # Test invalid values (negative, zero, non-numeric)
        invalid_vCg_values = [-10.0, 0.0, -0.1, "invalid", None]

        for invalid_vCg in invalid_vCg_values:
            with self.subTest(invalid_vCg=invalid_vCg):
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(vCg__E=invalid_vCg)

    def test_alpha_parameter_validation(self):
        """Test alpha parameter validation."""
        # Test valid range (-180, 180]
        valid_alpha_values = [-179.9, -90.0, 0.0, 45.0, 90.0, 180.0]

        for alpha in valid_alpha_values:
            with self.subTest(alpha=alpha):
                op = ps.operating_point.OperatingPoint(alpha=alpha)
                self.assertEqual(op.alpha, float(alpha))

        # Test invalid values (outside range)
        invalid_alpha_values = [180.1, -180.0, -180.1, 200.0, -200.0]

        for invalid_alpha in invalid_alpha_values:
            with self.subTest(invalid_alpha=invalid_alpha):
                with self.assertRaises(ValueError):
                    ps.operating_point.OperatingPoint(alpha=invalid_alpha)

        # Test non-numeric values
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(alpha="invalid")

    def test_beta_parameter_validation(self):
        """Test beta parameter validation."""
        # Test valid range (-180, 180]
        valid_beta_values = [-179.9, -90.0, 0.0, 10.0, 90.0, 180.0]

        for beta in valid_beta_values:
            with self.subTest(beta=beta):
                op = ps.operating_point.OperatingPoint(beta=beta)
                self.assertEqual(op.beta, float(beta))

        # Test invalid values (outside range)
        invalid_beta_values = [180.1, -180.0, -180.1, 200.0, -200.0]

        for invalid_beta in invalid_beta_values:
            with self.subTest(invalid_beta=invalid_beta):
                with self.assertRaises(ValueError):
                    ps.operating_point.OperatingPoint(beta=invalid_beta)

        # Test non-numeric values
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(beta="invalid")

    def test_externalFX_W_parameter_validation(self):
        """Test externalFX_W parameter validation."""
        # Test valid values (any number including negative, zero, positive)
        valid_external_force_values = [-100.0, -10.0, 0.0, 10.0, 50.0, 1000.0]

        for externalFX_W in valid_external_force_values:
            with self.subTest(externalFX_W=externalFX_W):
                op = ps.operating_point.OperatingPoint(externalFX_W=externalFX_W)
                self.assertEqual(op.externalFX_W, float(externalFX_W))

        # Test invalid types (string, None)
        invalid_external_force_values = ["invalid", None]

        for invalid_external_force in invalid_external_force_values:
            with self.subTest(invalid_external_force=invalid_external_force):
                with self.assertRaises(TypeError):
                    ps.operating_point.OperatingPoint(
                        externalFX_W=invalid_external_force
                    )

    def test_nu_parameter_validation(self):
        """Test nu parameter validation."""
        # Test valid positive floats
        valid_nu_values = [1e-6, 15.06e-6, 20.0e-6, 1e-5]

        for nu in valid_nu_values:
            with self.subTest(nu=nu):
                op = ps.operating_point.OperatingPoint(nu=nu)
                self.assertEqual(op.nu, float(nu))

        # Test invalid values (negative, zero, non-numeric)
        invalid_nu_values = [-1e-6, 0.0, -15.06e-6, "invalid", None]

        for invalid_nu in invalid_nu_values:
            with self.subTest(invalid_nu=invalid_nu):
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(nu=invalid_nu)

    def test_qInf__E_calculation(self):
        """Test qInf__E calculation accuracy."""
        # Verify qInf__E = 0.5 x rho x vCg__E^2
        expected_qInf = 0.5 * self.basic_op.rho * self.basic_op.vCg__E**2
        self.assertAlmostEqual(self.basic_op.qInf__E, expected_qInf, places=10)

        # Test with different configurations
        test_cases = [
            (self.high_speed_op, 0.5 * 1.225 * 100.0**2),
            (self.low_density_op, 0.5 * 0.3 * 10.0**2),
            (self.zero_alpha_beta_op, 0.5 * 1.225 * 10.0**2),
        ]

        for op, expected_qInf in test_cases:
            with self.subTest(op=op):
                self.assertAlmostEqual(op.qInf__E, expected_qInf, places=10)

    def test_qInf__E_scaling_with_velocity(self):
        """Test qInf__E quadratic scaling with velocity."""
        # Create OperatingPoints with different velocities
        op_v10 = ps.operating_point.OperatingPoint(vCg__E=10.0)
        op_v20 = ps.operating_point.OperatingPoint(vCg__E=20.0)

        # Verify quadratic scaling
        ratio = op_v20.qInf__E / op_v10.qInf__E
        self.assertAlmostEqual(ratio, 4.0, places=10)

    def test_qInf__E_scaling_with_density(self):
        """Test qInf__E linear scaling with density."""
        # Create OperatingPoints with different densities
        op_rho1 = ps.operating_point.OperatingPoint(rho=1.0)
        op_rho2 = ps.operating_point.OperatingPoint(rho=2.0)

        # Verify linear scaling
        ratio = op_rho2.qInf__E / op_rho1.qInf__E
        self.assertAlmostEqual(ratio, 2.0, places=10)

    def test_T_pas_G_Cg_to_W_Cg_shape_and_type(self):
        """Test T_pas_G_Cg_to_W_Cg shape and type."""
        T = self.basic_op.T_pas_G_Cg_to_W_Cg

        # Should be 4x4 matrix
        self.assertEqual(T.shape, (4, 4))

        # Should be ndarray of floats
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_W_Cg_to_G_Cg_shape_and_type(self):
        """Test T_pas_W_Cg_to_G_Cg shape and type."""
        T = self.basic_op.T_pas_W_Cg_to_G_Cg

        # Should be 4x4 matrix
        self.assertEqual(T.shape, (4, 4))

        # Should be ndarray of floats
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_transformation_matrices_are_inverses(self):
        """Test that transformation matrices are inverses of each other."""
        # Test for all fixtures
        fixtures = [
            self.basic_op,
            self.zero_alpha_beta_op,
            self.high_alpha_op,
            self.negative_alpha_op,
            self.nonzero_beta_op,
            self.boundary_alpha_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_forward = op.T_pas_G_Cg_to_W_Cg
                T_inverse = op.T_pas_W_Cg_to_G_Cg

                # Verify T_forward @ T_inverse = Identity
                identity = T_forward @ T_inverse
                npt.assert_allclose(identity, np.eye(4), atol=1e-14)

                # Verify T_inverse @ T_forward = Identity
                identity_reverse = T_inverse @ T_forward
                npt.assert_allclose(identity_reverse, np.eye(4), atol=1e-14)

    def test_transformation_zero_alpha_beta(self):
        """Test transformation with zero alpha and beta."""
        op = self.zero_alpha_beta_op

        # With alpha = beta = 0, wind axes align with body axes.
        # The freestream direction in wind axes is [-1, 0, 0].
        # Transform to geometry axes to verify direction.
        vInfHat_G__E = op.vInfHat_G__E

        # Geometry axes: +x = aft, +y = right, +z = up
        # Body axes: +x = forward, +y = right, +z = down
        # With alpha = beta = 0, wind axes align with body axes.
        # So freestream [-1, 0, 0] in wind axes means freestream comes from
        # backward in body axes (forward velocity).
        # In geometry axes, this should be [+1, 0, 0] (forward in geometry is
        # aft direction, opposite of body).
        npt.assert_allclose(vInfHat_G__E[0], 1.0, atol=1e-14)
        npt.assert_allclose(vInfHat_G__E[1], 0.0, atol=1e-14)
        npt.assert_allclose(vInfHat_G__E[2], 0.0, atol=1e-14)

    def test_transformation_alpha_only(self):
        """Test transformation with alpha only (beta = 0)."""
        # Positive alpha means the nose points above the direction of travel,
        # so the relative wind comes from below. In geometry axes (+z = up),
        # the freestream velocity vector should have a positive z-component.
        op_positive_alpha = ps.operating_point.OperatingPoint(alpha=10.0, beta=0.0)
        vInf_G__E_pos = op_positive_alpha.vInfHat_G__E

        # With positive alpha, the freestream should have a positive z-component
        # in geometry axes (wind comes from below).
        self.assertGreater(vInf_G__E_pos[2], 0.0)

        # The y-component should be approximately zero (no sideslip).
        npt.assert_allclose(vInf_G__E_pos[1], 0.0, atol=1e-14)

        # Negative alpha means the nose points below the direction of travel,
        # so the relative wind comes from above.
        op_negative_alpha = ps.operating_point.OperatingPoint(alpha=-10.0, beta=0.0)
        vInf_G__E_neg = op_negative_alpha.vInfHat_G__E

        # With negative alpha, the freestream should have a negative z-component
        # in geometry axes (wind comes from above).
        self.assertLess(vInf_G__E_neg[2], 0.0)

        # The y-component should be approximately zero (no sideslip).
        npt.assert_allclose(vInf_G__E_neg[1], 0.0, atol=1e-14)

    def test_transformation_beta_only(self):
        """Test transformation with beta only (alpha = 0)."""
        # Positive beta means the nose points to the left of the direction of
        # travel, so the airplane moves to the right of where the nose points.
        # The relative wind comes from the right. In geometry axes (+y = right),
        # the freestream velocity vector should have a negative y-component
        # (pointing left, opposite the airplane's rightward motion).
        op_positive_beta = ps.operating_point.OperatingPoint(alpha=0.0, beta=10.0)
        vInf_G__E_pos = op_positive_beta.vInfHat_G__E

        # With positive beta, the freestream should have a negative y-component
        # in geometry axes (wind from right, pointing left).
        self.assertLess(vInf_G__E_pos[1], 0.0)

        # The z-component should be approximately zero (no vertical component).
        npt.assert_allclose(vInf_G__E_pos[2], 0.0, atol=1e-14)

        # Negative beta means the nose points to the right of the direction of
        # travel, so the airplane moves to the left of where the nose points.
        # The relative wind comes from the left.
        op_negative_beta = ps.operating_point.OperatingPoint(alpha=0.0, beta=-10.0)
        vInf_G__E_neg = op_negative_beta.vInfHat_G__E

        # With negative beta, the freestream should have a positive y-component
        # in geometry axes (wind from left, pointing right).
        self.assertGreater(vInf_G__E_neg[1], 0.0)

        # The z-component should be approximately zero (no vertical component).
        npt.assert_allclose(vInf_G__E_neg[2], 0.0, atol=1e-14)

    def test_transformation_both_alpha_beta(self):
        """Test transformation with both alpha and beta non-zero."""
        # Test that transformations remain valid with both angles present.
        op = self.nonzero_beta_op
        T = op.T_pas_G_Cg_to_W_Cg

        # Extract the rotation part (top-left 3x3 submatrix).
        R = T[:3, :3]

        # Verify that the rotation matrix is orthonormal.
        npt.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        npt.assert_allclose(np.linalg.det(R), 1.0, atol=1e-14)

    def test_transformation_boundary_angles(self):
        """Test transformation with boundary angle values."""
        # Test with alpha at boundary
        op_alpha_boundary = ps.operating_point.OperatingPoint(alpha=180.0, beta=0.0)
        T_alpha = op_alpha_boundary.T_pas_G_Cg_to_W_Cg
        R_alpha = T_alpha[:3, :3]

        # Verify no numerical issues
        npt.assert_allclose(R_alpha @ R_alpha.T, np.eye(3), atol=1e-14)
        self.assertFalse(np.any(np.isnan(T_alpha)))
        self.assertFalse(np.any(np.isinf(T_alpha)))

        # Test with beta at boundary
        op_beta_boundary = ps.operating_point.OperatingPoint(alpha=0.0, beta=180.0)
        T_beta = op_beta_boundary.T_pas_G_Cg_to_W_Cg
        R_beta = T_beta[:3, :3]

        # Verify no numerical issues
        npt.assert_allclose(R_beta @ R_beta.T, np.eye(3), atol=1e-14)
        self.assertFalse(np.any(np.isnan(T_beta)))
        self.assertFalse(np.any(np.isinf(T_beta)))

    def test_vInfHat_G__E_is_unit_vector(self):
        """Test that vInfHat_G__E is a unit vector."""
        # Test for all fixtures
        fixtures = [
            self.basic_op,
            self.zero_alpha_beta_op,
            self.high_alpha_op,
            self.negative_alpha_op,
            self.nonzero_beta_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                vInfHat = op.vInfHat_G__E
                magnitude = np.linalg.norm(vInfHat)
                npt.assert_allclose(magnitude, 1.0, atol=1e-14)

    def test_vInfHat_G__E_direction(self):
        """Test vInfHat_G__E direction consistency."""
        # For zero alpha and beta, verify direction aligns as expected.
        # In wind axes, freestream is [-1, 0, 0].
        # Transform to geometry axes and verify it's a unit vector.
        op = self.zero_alpha_beta_op
        vInfHat = op.vInfHat_G__E

        # Should be a 3-element vector
        self.assertEqual(len(vInfHat), 3)

        # Should be unit vector
        npt.assert_allclose(np.linalg.norm(vInfHat), 1.0, atol=1e-14)

    def test_vInf_G__E_magnitude(self):
        """Test that vInf_G__E magnitude equals vCg__E."""
        # Test for all fixtures
        fixtures = [
            self.basic_op,
            self.high_speed_op,
            self.high_alpha_op,
            self.nonzero_beta_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                vInf = op.vInf_G__E
                magnitude = np.linalg.norm(vInf)
                npt.assert_allclose(magnitude, op.vCg__E, atol=1e-14)

    def test_vInf_G__E_equals_vInfHat_times_speed(self):
        """Test that vInf_G__E equals vInfHat_G__E times vCg__E."""
        # Test for all fixtures
        fixtures = [
            self.basic_op,
            self.zero_alpha_beta_op,
            self.high_alpha_op,
            self.negative_alpha_op,
            self.nonzero_beta_op,
            self.high_speed_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                vInf = op.vInf_G__E
                expected_vInf = op.vInfHat_G__E * op.vCg__E
                npt.assert_allclose(vInf, expected_vInf, atol=1e-14)

    def test_vInf_G__E_with_various_alpha_beta(self):
        """Test vInf_G__E direction changes with alpha and beta."""
        # Create OperatingPoints with different alpha/beta combinations
        test_cases = [
            {"alpha": 0.0, "beta": 0.0},
            {"alpha": 10.0, "beta": 0.0},
            {"alpha": 0.0, "beta": 10.0},
            {"alpha": 45.0, "beta": 0.0},
            {"alpha": 0.0, "beta": 45.0},
            {"alpha": 30.0, "beta": 15.0},
        ]

        for params in test_cases:
            with self.subTest(params=params):
                op = ps.operating_point.OperatingPoint(**params)
                vInf = op.vInf_G__E

                # Should be 3-element vector
                self.assertEqual(len(vInf), 3)

                # Should have magnitude equal to vCg__E
                magnitude = np.linalg.norm(vInf)
                npt.assert_allclose(magnitude, op.vCg__E, atol=1e-14)

    def test_multiple_operating_points_independent(self):
        """Test that multiple OperatingPoints are independent."""
        # Create two OperatingPoints with different parameters
        op1 = ps.operating_point.OperatingPoint(alpha=10.0, vCg__E=20.0)
        op2 = ps.operating_point.OperatingPoint(alpha=30.0, vCg__E=50.0)

        # Verify they have different properties
        self.assertNotEqual(op1.alpha, op2.alpha)
        self.assertNotEqual(op1.vCg__E, op2.vCg__E)
        self.assertNotEqual(op1.qInf__E, op2.qInf__E)

        # Verify their transformation matrices are different
        self.assertFalse(np.allclose(op1.T_pas_G_Cg_to_W_Cg, op2.T_pas_G_Cg_to_W_Cg))

        # Verify their velocity vectors are different
        self.assertFalse(np.allclose(op1.vInf_G__E, op2.vInf_G__E))

    def test_comprehensive_operating_point_properties(self):
        """Test all properties on various fixtures."""
        fixtures = [
            (self.basic_op, "basic"),
            (self.zero_alpha_beta_op, "zero_alpha_beta"),
            (self.high_alpha_op, "high_alpha"),
            (self.negative_alpha_op, "negative_alpha"),
            (self.nonzero_beta_op, "nonzero_beta"),
            (self.high_speed_op, "high_speed"),
            (self.low_density_op, "low_density"),
            (self.with_external_force_op, "with_external_force"),
            (self.custom_viscosity_op, "custom_viscosity"),
            (self.boundary_alpha_op, "boundary_alpha"),
        ]

        for op, fixture_name in fixtures:
            with self.subTest(fixture=fixture_name):
                # Test all attributes exist and have correct types
                self.assertIsInstance(op.rho, float)
                self.assertIsInstance(op.vCg__E, float)
                self.assertIsInstance(op.alpha, float)
                self.assertIsInstance(op.beta, float)
                self.assertIsInstance(op.externalFX_W, float)
                self.assertIsInstance(op.nu, float)

                # Test all properties return correct types
                self.assertIsInstance(op.qInf__E, float)
                self.assertIsInstance(op.T_pas_G_Cg_to_W_Cg, np.ndarray)
                self.assertIsInstance(op.T_pas_W_Cg_to_G_Cg, np.ndarray)
                self.assertIsInstance(op.vInfHat_G__E, np.ndarray)
                self.assertIsInstance(op.vInf_G__E, np.ndarray)

                # Test no NaN or Inf values in properties
                self.assertFalse(np.isnan(op.qInf__E))
                self.assertFalse(np.any(np.isnan(op.T_pas_G_Cg_to_W_Cg)))
                self.assertFalse(np.any(np.isnan(op.T_pas_W_Cg_to_G_Cg)))
                self.assertFalse(np.any(np.isnan(op.vInfHat_G__E)))
                self.assertFalse(np.any(np.isnan(op.vInf_G__E)))

    def test_edge_case_angle_boundaries(self):
        """Test angle validation at boundary values."""
        # Test boundary values that should be valid
        valid_boundary_cases = [
            {"alpha": 180.0, "beta": 0.0},
            {"alpha": 0.0, "beta": 180.0},
            {"alpha": -179.999, "beta": 0.0},
            {"alpha": 0.0, "beta": -179.999},
            {"alpha": 180.0, "beta": 180.0},
        ]

        for params in valid_boundary_cases:
            with self.subTest(params=params):
                op = ps.operating_point.OperatingPoint(**params)
                self.assertAlmostEqual(op.alpha, params["alpha"], places=10)
                self.assertAlmostEqual(op.beta, params["beta"], places=10)

    def test_very_low_speed(self):
        """Test with very low but valid speed."""
        op = ps.operating_point.OperatingPoint(vCg__E=0.01)

        # Should still calculate qInf correctly
        expected_qInf = 0.5 * 1.225 * 0.01**2
        self.assertAlmostEqual(op.qInf__E, expected_qInf, places=10)

        # Should still produce valid velocity vectors
        self.assertEqual(len(op.vInf_G__E), 3)
        npt.assert_allclose(np.linalg.norm(op.vInf_G__E), 0.01, atol=1e-14)

    def test_very_high_speed(self):
        """Test with very high speed."""
        op = ps.operating_point.OperatingPoint(vCg__E=300.0)

        # Should still calculate qInf correctly
        expected_qInf = 0.5 * 1.225 * 300.0**2
        self.assertAlmostEqual(op.qInf__E, expected_qInf, places=10)

        # Should still produce valid velocity vectors
        self.assertEqual(len(op.vInf_G__E), 3)
        npt.assert_allclose(np.linalg.norm(op.vInf_G__E), 300.0, atol=1e-12)

    def test_extreme_density_values(self):
        """Test with extreme but valid density values."""
        # Very low density
        op_low = ps.operating_point.OperatingPoint(rho=0.01)
        expected_qInf_low = 0.5 * 0.01 * 10.0**2
        self.assertAlmostEqual(op_low.qInf__E, expected_qInf_low, places=10)

        # Very high density
        op_high = ps.operating_point.OperatingPoint(rho=10.0)
        expected_qInf_high = 0.5 * 10.0 * 10.0**2
        self.assertAlmostEqual(op_high.qInf__E, expected_qInf_high, places=10)

        # Both should produce valid transformations
        self.assertEqual(op_low.T_pas_G_Cg_to_W_Cg.shape, (4, 4))
        self.assertEqual(op_high.T_pas_G_Cg_to_W_Cg.shape, (4, 4))


if __name__ == "__main__":
    unittest.main()
