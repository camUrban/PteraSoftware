"""This module contains a class to test OperatingPoints."""

import copy
import unittest
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import operating_point_fixtures


class TestOperatingPoint(unittest.TestCase):
    """This is a class with functions to test OperatingPoints."""

    def setUp(self) -> None:
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
        self.negative_beta_op = (
            operating_point_fixtures.make_negative_beta_operating_point_fixture()
        )
        self.boundary_beta_op = (
            operating_point_fixtures.make_boundary_beta_operating_point_fixture()
        )
        self.combined_boundary_angles_op = (
            operating_point_fixtures.make_combined_boundary_angles_operating_point_fixture()
        )
        self.very_low_speed_op = (
            operating_point_fixtures.make_very_low_speed_operating_point_fixture()
        )
        self.integer_parameters_op = (
            operating_point_fixtures.make_integer_parameters_operating_point_fixture()
        )
        self.negative_external_force_op = (
            operating_point_fixtures.make_negative_external_force_operating_point_fixture()
        )
        self.near_boundary_alpha_op = (
            operating_point_fixtures.make_near_boundary_alpha_operating_point_fixture()
        )
        self.near_boundary_beta_op = (
            operating_point_fixtures.make_near_boundary_beta_operating_point_fixture()
        )
        self.with_attitude_angles_op = (
            operating_point_fixtures.make_with_attitude_angles_operating_point_fixture()
        )
        self.with_cg_position_op = (
            operating_point_fixtures.make_with_cg_position_operating_point_fixture()
        )
        self.with_ground_surface_op = (
            operating_point_fixtures.make_with_ground_surface_operating_point_fixture()
        )
        self.with_tilted_surface_op = (
            operating_point_fixtures.make_with_tilted_surface_operating_point_fixture()
        )

    def test_initialization_valid_parameters(self) -> None:
        """Test OperatingPoint initialization with valid parameters."""
        # Test basic OperatingPoint initialization
        self.assertIsInstance(self.basic_op, ps.operating_point.OperatingPoint)
        self.assertEqual(self.basic_op.rho, 1.225)
        self.assertEqual(self.basic_op.vCg__E, 10.0)
        self.assertEqual(self.basic_op.alpha, 5.0)
        self.assertEqual(self.basic_op.beta, 0.0)
        self.assertEqual(self.basic_op.externalFX_W, 0.0)
        self.assertEqual(self.basic_op.nu, 15.06e-6)

    def test_initialization_with_defaults(self) -> None:
        """Test that default values are applied correctly."""
        # Create an OperatingPoint with all defaults. Omitting alpha issues a
        # FutureWarning because its implicit default will change in v6.0.0.
        with self.assertWarns(FutureWarning):
            op_default = ps.operating_point.OperatingPoint()

        # Verify default values
        self.assertEqual(op_default.rho, 1.225)
        self.assertEqual(op_default.vCg__E, 10.0)
        self.assertEqual(op_default.alpha, 5.0)
        self.assertEqual(op_default.beta, 0.0)
        self.assertEqual(op_default.externalFX_W, 0.0)
        self.assertEqual(op_default.nu, 15.06e-6)

    def test_rho_parameter_validation(self) -> None:
        """Test rho parameter validation."""
        # Test valid positive floats
        valid_rho_values = [0.5, 1.225, 2.0, 10.0]

        for rho in valid_rho_values:
            with self.subTest(rho=rho):
                op = ps.operating_point.OperatingPoint(alpha=5.0, rho=rho)
                self.assertEqual(op.rho, float(rho))

        # Test invalid values (negative, zero, non-numeric)
        invalid_rho_values: list[Any] = [-1.0, 0.0, -0.5, "invalid", None]

        for invalid_rho in invalid_rho_values:
            with self.subTest(invalid_rho=invalid_rho):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(alpha=5.0, rho=invalid_rho)

    def test_vCg__E_parameter_validation(self) -> None:
        """Test vCg__E parameter validation."""
        # Test valid positive floats
        valid_vCg_values = [0.01, 1.0, 10.0, 100.0, 300.0]

        for vCg__E in valid_vCg_values:
            with self.subTest(vCg__E=vCg__E):
                op = ps.operating_point.OperatingPoint(alpha=5.0, vCg__E=vCg__E)
                self.assertEqual(op.vCg__E, float(vCg__E))

        # Test invalid values (negative, zero, non-numeric)
        invalid_vCg_values: list[Any] = [-10.0, 0.0, -0.1, "invalid", None]

        for invalid_vCg in invalid_vCg_values:
            with self.subTest(invalid_vCg=invalid_vCg):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(alpha=5.0, vCg__E=invalid_vCg)

    def test_alpha_parameter_validation(self) -> None:
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
        bad_alpha: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(alpha=bad_alpha)

    def test_beta_parameter_validation(self) -> None:
        """Test beta parameter validation."""
        # Test valid range (-180, 180]
        valid_beta_values = [-179.9, -90.0, 0.0, 10.0, 90.0, 180.0]

        for beta in valid_beta_values:
            with self.subTest(beta=beta):
                op = ps.operating_point.OperatingPoint(alpha=5.0, beta=beta)
                self.assertEqual(op.beta, float(beta))

        # Test invalid values (outside range)
        invalid_beta_values = [180.1, -180.0, -180.1, 200.0, -200.0]

        for invalid_beta in invalid_beta_values:
            with self.subTest(invalid_beta=invalid_beta):
                with self.assertRaises(ValueError):
                    ps.operating_point.OperatingPoint(alpha=5.0, beta=invalid_beta)

        # Test non-numeric values
        bad_beta: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(alpha=5.0, beta=bad_beta)

    def test_externalFX_W_parameter_validation(self) -> None:
        """Test externalFX_W parameter validation."""
        # Test valid values (any number including negative, zero, positive)
        valid_external_force_values = [-100.0, -10.0, 0.0, 10.0, 50.0, 1000.0]

        for externalFX_W in valid_external_force_values:
            with self.subTest(externalFX_W=externalFX_W):
                op = ps.operating_point.OperatingPoint(
                    alpha=5.0, externalFX_W=externalFX_W
                )
                self.assertEqual(op.externalFX_W, float(externalFX_W))

        # Test invalid types (string, None)
        invalid_external_force_values: list[Any] = ["invalid", None]

        for invalid_external_force in invalid_external_force_values:
            with self.subTest(invalid_external_force=invalid_external_force):
                with self.assertRaises(TypeError):
                    ps.operating_point.OperatingPoint(
                        alpha=5.0, externalFX_W=invalid_external_force
                    )

    def test_nu_parameter_validation(self) -> None:
        """Test nu parameter validation."""
        # Test valid positive floats
        valid_nu_values = [1e-6, 15.06e-6, 20.0e-6, 1e-5]

        for nu in valid_nu_values:
            with self.subTest(nu=nu):
                op = ps.operating_point.OperatingPoint(alpha=5.0, nu=nu)
                self.assertEqual(op.nu, float(nu))

        # Test invalid values (negative, zero, non-numeric)
        invalid_nu_values: list[Any] = [-1e-6, 0.0, -15.06e-6, "invalid", None]

        for invalid_nu in invalid_nu_values:
            with self.subTest(invalid_nu=invalid_nu):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(alpha=5.0, nu=invalid_nu)

    def test_qInf__E_calculation(self) -> None:
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

    def test_qInf__E_scaling_with_velocity(self) -> None:
        """Test qInf__E quadratic scaling with velocity."""
        # Create OperatingPoints with different velocities
        op_v10 = ps.operating_point.OperatingPoint(alpha=5.0, vCg__E=10.0)
        op_v20 = ps.operating_point.OperatingPoint(alpha=5.0, vCg__E=20.0)

        # Verify quadratic scaling
        ratio = op_v20.qInf__E / op_v10.qInf__E
        self.assertAlmostEqual(ratio, 4.0, places=10)

    def test_qInf__E_scaling_with_density(self) -> None:
        """Test qInf__E linear scaling with density."""
        # Create OperatingPoints with different densities
        op_rho1 = ps.operating_point.OperatingPoint(alpha=5.0, rho=1.0)
        op_rho2 = ps.operating_point.OperatingPoint(alpha=5.0, rho=2.0)

        # Verify linear scaling
        ratio = op_rho2.qInf__E / op_rho1.qInf__E
        self.assertAlmostEqual(ratio, 2.0, places=10)

    def test_T_pas_GP1_CgP1_to_W_CgP1_shape_and_type(self) -> None:
        """Test T_pas_GP1_CgP1_to_W_CgP1 shape and type."""
        T = self.basic_op.T_pas_GP1_CgP1_to_W_CgP1

        # Should be 4x4 matrix
        self.assertEqual(T.shape, (4, 4))

        # Should be ndarray of floats
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_W_CgP1_to_GP1_CgP1_shape_and_type(self) -> None:
        """Test T_pas_W_CgP1_to_GP1_CgP1 shape and type."""
        T = self.basic_op.T_pas_W_CgP1_to_GP1_CgP1

        # Should be 4x4 matrix
        self.assertEqual(T.shape, (4, 4))

        # Should be ndarray of floats
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_transformation_matrices_are_inverses(self) -> None:
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
                T_forward = op.T_pas_GP1_CgP1_to_W_CgP1
                T_inverse = op.T_pas_W_CgP1_to_GP1_CgP1

                # Verify T_forward @ T_inverse = Identity
                identity = T_forward @ T_inverse
                npt.assert_allclose(identity, np.eye(4), atol=1e-14)

                # Verify T_inverse @ T_forward = Identity
                identity_reverse = T_inverse @ T_forward
                npt.assert_allclose(identity_reverse, np.eye(4), atol=1e-14)

    def test_transformation_zero_alpha_beta(self) -> None:
        """Test transformation with zero alpha and beta."""
        op = self.zero_alpha_beta_op

        # With alpha = beta = 0, wind axes align with body axes. The freestream
        # direction in wind axes is [-1, 0, 0]. Transform to the first Airplane's
        # geometry axes to verify direction.
        vInfHat_GP1__E = op.vInfHat_GP1__E

        # In geometry axes, +x points aft, +y points right, and +z points up. In body
        # axes, +x points forward, +y points right, and +z points down. With alpha =
        # beta = 0, wind axes align with body axes. So freestream [-1, 0, 0] in wind
        # axes means freestream comes from backward in body axes (forward velocity). In
        # geometry axes, this should be [+1, 0, 0] (forward in geometry axes is the aft
        # direction, opposite of body axes).
        npt.assert_allclose(vInfHat_GP1__E[0], 1.0, atol=1e-14)
        npt.assert_allclose(vInfHat_GP1__E[1], 0.0, atol=1e-14)
        npt.assert_allclose(vInfHat_GP1__E[2], 0.0, atol=1e-14)

    def test_transformation_alpha_only(self) -> None:
        """Test transformation with alpha only (beta = 0)."""
        # Positive alpha means the nose points above the direction of travel, so the
        # relative wind comes from below. In the first Airplane's geometry axes (+z =
        # up), the freestream velocity vector should have a positive z component.
        op_positive_alpha = ps.operating_point.OperatingPoint(alpha=10.0, beta=0.0)
        vInf_GP1__E_pos = op_positive_alpha.vInfHat_GP1__E

        # With positive alpha, the freestream should have a positive z component in the
        # first Airplane's geometry axes (wind comes from below).
        self.assertGreater(vInf_GP1__E_pos[2], 0.0)

        # The y component should be approximately zero (no sideslip).
        npt.assert_allclose(vInf_GP1__E_pos[1], 0.0, atol=1e-14)

        # Negative alpha means the nose points below the direction of travel, so the
        # relative wind comes from above.
        op_negative_alpha = ps.operating_point.OperatingPoint(alpha=-10.0, beta=0.0)
        vInf_GP1__E_neg = op_negative_alpha.vInfHat_GP1__E

        # With negative alpha, the freestream should have a negative z component in the
        # first Airplane's geometry axes (wind comes from above).
        self.assertLess(vInf_GP1__E_neg[2], 0.0)

        # The y component should be approximately zero (no sideslip).
        npt.assert_allclose(vInf_GP1__E_neg[1], 0.0, atol=1e-14)

    def test_transformation_beta_only(self) -> None:
        """Test transformation with beta only (alpha = 0)."""
        # Positive beta means the nose points to the left of the direction of travel, so
        # the airplane moves to the right of where the nose points. The relative wind
        # comes from the right. In the first Airplane's geometry axes (+y = right), the
        # freestream velocity vector should have a negative y component (pointing left,
        # opposite the airplane's rightward motion).
        op_positive_beta = ps.operating_point.OperatingPoint(alpha=0.0, beta=10.0)
        vInf_GP1__E_pos = op_positive_beta.vInfHat_GP1__E

        # With positive beta, the freestream should have a negative y component in the
        # first Airplane's geometry axes (wind from right, pointing left).
        self.assertLess(vInf_GP1__E_pos[1], 0.0)

        # The z component should be approximately zero (no vertical component).
        npt.assert_allclose(vInf_GP1__E_pos[2], 0.0, atol=1e-14)

        # Negative beta means the nose points to the right of the direction of travel,
        # so the airplane moves to the left of where the nose points. The relative wind
        # comes from the left.
        op_negative_beta = ps.operating_point.OperatingPoint(alpha=0.0, beta=-10.0)
        vInf_GP1__E_neg = op_negative_beta.vInfHat_GP1__E

        # With negative beta, the freestream should have a positive y component in the
        # first Airplane's geometry axes (wind from left, pointing right).
        self.assertGreater(vInf_GP1__E_neg[1], 0.0)

        # The z component should be approximately zero (no vertical component).
        npt.assert_allclose(vInf_GP1__E_neg[2], 0.0, atol=1e-14)

    def test_transformation_both_alpha_beta(self) -> None:
        """Test transformation with both alpha and beta non-zero."""
        # Test that transformations remain valid with both angles present.
        op = self.nonzero_beta_op
        T = op.T_pas_GP1_CgP1_to_W_CgP1

        # Extract the rotation part (top-left 3x3 submatrix).
        R = T[:3, :3]

        # Verify that the rotation matrix is orthonormal.
        npt.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        npt.assert_allclose(np.linalg.det(R), 1.0, atol=1e-14)

    def test_transformation_boundary_angles(self) -> None:
        """Test transformation with boundary angle values."""
        # Test with alpha at boundary
        op_alpha_boundary = ps.operating_point.OperatingPoint(alpha=180.0, beta=0.0)
        T_alpha = op_alpha_boundary.T_pas_GP1_CgP1_to_W_CgP1
        R_alpha = T_alpha[:3, :3]

        # Verify no numerical issues
        npt.assert_allclose(R_alpha @ R_alpha.T, np.eye(3), atol=1e-14)
        self.assertFalse(np.any(np.isnan(T_alpha)))
        self.assertFalse(np.any(np.isinf(T_alpha)))

        # Test with beta at boundary
        op_beta_boundary = ps.operating_point.OperatingPoint(alpha=0.0, beta=180.0)
        T_beta = op_beta_boundary.T_pas_GP1_CgP1_to_W_CgP1
        R_beta = T_beta[:3, :3]

        # Verify no numerical issues
        npt.assert_allclose(R_beta @ R_beta.T, np.eye(3), atol=1e-14)
        self.assertFalse(np.any(np.isnan(T_beta)))
        self.assertFalse(np.any(np.isinf(T_beta)))

    def test_vInfHat_GP1__E_is_unit_vector(self) -> None:
        """Test that vInfHat_GP1__E is a unit vector."""
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
                vInfHat_GP1__E = op.vInfHat_GP1__E
                magnitude = np.linalg.norm(vInfHat_GP1__E)
                npt.assert_allclose(magnitude, 1.0, atol=1e-14)

    def test_vInfHat_GP1__E_direction(self) -> None:
        """Test vInfHat_GP1__E direction consistency."""
        # For zero alpha and beta, verify direction aligns as expected. In wind axes,
        # freestream is [-1, 0, 0]. Transform to the first Airplane's geometry axes and
        # verify it's a unit vector.
        op = self.zero_alpha_beta_op
        vInfHat_GP1__E = op.vInfHat_GP1__E

        # Should be a 3-element vector
        self.assertEqual(len(vInfHat_GP1__E), 3)

        # Should be unit vector
        npt.assert_allclose(np.linalg.norm(vInfHat_GP1__E), 1.0, atol=1e-14)

    def test_vInf_GP1__E_magnitude(self) -> None:
        """Test that vInf_GP1__E magnitude equals vCg__E."""
        # Test for all fixtures
        fixtures = [
            self.basic_op,
            self.high_speed_op,
            self.high_alpha_op,
            self.nonzero_beta_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                vInfHat_GP1__E = op.vInf_GP1__E
                magnitude = np.linalg.norm(vInfHat_GP1__E)
                npt.assert_allclose(magnitude, op.vCg__E, atol=1e-14)

    def test_vInf_GP1__E_equals_vInfHat_times_speed(self) -> None:
        """Test that vInf_GP1__E equals vInfHat_GP1__E times vCg__E."""
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
                vInfHat_GP1__E = op.vInf_GP1__E
                expected_vInfHat_GP1__E = op.vInfHat_GP1__E * op.vCg__E
                npt.assert_allclose(vInfHat_GP1__E, expected_vInfHat_GP1__E, atol=1e-14)

    def test_vInf_GP1__E_with_various_alpha_beta(self) -> None:
        """Test vInf_GP1__E direction changes with alpha and beta."""
        # Create OperatingPoints with different alpha/beta combinations
        test_cases: list[dict[str, Any]] = [
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
                vInfHat_GP1__E = op.vInf_GP1__E

                # Should be 3-element vector
                self.assertEqual(len(vInfHat_GP1__E), 3)

                # Should have magnitude equal to vCg__E
                magnitude = np.linalg.norm(vInfHat_GP1__E)
                npt.assert_allclose(magnitude, op.vCg__E, atol=1e-14)

    def test_multiple_operating_points_independent(self) -> None:
        """Test that multiple OperatingPoints are independent."""
        # Create two OperatingPoints with different parameters
        op1 = ps.operating_point.OperatingPoint(alpha=10.0, vCg__E=20.0)
        op2 = ps.operating_point.OperatingPoint(alpha=30.0, vCg__E=50.0)

        # Verify they have different properties
        self.assertNotEqual(op1.alpha, op2.alpha)
        self.assertNotEqual(op1.vCg__E, op2.vCg__E)
        self.assertNotEqual(op1.qInf__E, op2.qInf__E)

        # Verify their transformation matrices are different
        self.assertFalse(
            np.allclose(op1.T_pas_GP1_CgP1_to_W_CgP1, op2.T_pas_GP1_CgP1_to_W_CgP1)
        )

        # Verify their velocity vectors are different
        self.assertFalse(np.allclose(op1.vInf_GP1__E, op2.vInf_GP1__E))

    def test_comprehensive_operating_point_properties(self) -> None:
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
                self.assertIsInstance(op.T_pas_GP1_CgP1_to_W_CgP1, np.ndarray)
                self.assertIsInstance(op.T_pas_W_CgP1_to_GP1_CgP1, np.ndarray)
                self.assertIsInstance(op.vInfHat_GP1__E, np.ndarray)
                self.assertIsInstance(op.vInf_GP1__E, np.ndarray)

                # Test no NaN or Inf values in properties
                self.assertFalse(np.isnan(op.qInf__E))
                self.assertFalse(np.any(np.isnan(op.T_pas_GP1_CgP1_to_W_CgP1)))
                self.assertFalse(np.any(np.isnan(op.T_pas_W_CgP1_to_GP1_CgP1)))
                self.assertFalse(np.any(np.isnan(op.vInfHat_GP1__E)))
                self.assertFalse(np.any(np.isnan(op.vInf_GP1__E)))

    def test_edge_case_angle_boundaries(self) -> None:
        """Test angle validation at boundary values."""
        # Test boundary values that should be valid
        valid_boundary_cases: list[dict[str, Any]] = [
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

    def test_very_high_speed(self) -> None:
        """Test with very high speed."""
        op = ps.operating_point.OperatingPoint(alpha=5.0, vCg__E=300.0)

        # Should still calculate qInf correctly
        expected_qInf = 0.5 * 1.225 * 300.0**2
        self.assertAlmostEqual(op.qInf__E, expected_qInf, places=10)

        # Should still produce valid velocity vectors
        self.assertEqual(len(op.vInf_GP1__E), 3)
        npt.assert_allclose(np.linalg.norm(op.vInf_GP1__E), 300.0, atol=1e-12)

    def test_extreme_density_values(self) -> None:
        """Test with extreme but valid density values."""
        # Very low density
        op_low = ps.operating_point.OperatingPoint(alpha=5.0, rho=0.01)
        expected_qInf_low = 0.5 * 0.01 * 10.0**2
        self.assertAlmostEqual(op_low.qInf__E, expected_qInf_low, places=10)

        # Very high density
        op_high = ps.operating_point.OperatingPoint(alpha=5.0, rho=10.0)
        expected_qInf_high = 0.5 * 10.0 * 10.0**2
        self.assertAlmostEqual(op_high.qInf__E, expected_qInf_high, places=10)

        # Both should produce valid transformations
        self.assertEqual(op_low.T_pas_GP1_CgP1_to_W_CgP1.shape, (4, 4))
        self.assertEqual(op_high.T_pas_GP1_CgP1_to_W_CgP1.shape, (4, 4))

    def test_immutable_attributes_raise_attribute_error(self) -> None:
        """Test that setting read-only properties raises AttributeError."""
        op = self.basic_op

        # Test all immutable scalar attributes
        with self.assertRaises(AttributeError):
            setattr(op, "rho", 2.0)

        with self.assertRaises(AttributeError):
            setattr(op, "vCg__E", 20.0)

        with self.assertRaises(AttributeError):
            setattr(op, "alpha", 10.0)

        with self.assertRaises(AttributeError):
            setattr(op, "beta", 5.0)

        with self.assertRaises(AttributeError):
            setattr(op, "externalFX_W", 100.0)

        with self.assertRaises(AttributeError):
            setattr(op, "nu", 20.0e-6)

    def test_derived_property_caching(self) -> None:
        """Test that derived properties are cached and return same objects."""
        op = self.basic_op

        # Access qInf__E twice, should return same value (float)
        qInf_1 = op.qInf__E
        qInf_2 = op.qInf__E
        self.assertEqual(qInf_1, qInf_2)

        # Access transformation matrices twice, should return same objects
        T_forward_1 = op.T_pas_GP1_CgP1_to_W_CgP1
        T_forward_2 = op.T_pas_GP1_CgP1_to_W_CgP1
        self.assertIs(T_forward_1, T_forward_2)

        T_inverse_1 = op.T_pas_W_CgP1_to_GP1_CgP1
        T_inverse_2 = op.T_pas_W_CgP1_to_GP1_CgP1
        self.assertIs(T_inverse_1, T_inverse_2)

        # Access velocity vectors twice, should return same objects
        vInfHat_1 = op.vInfHat_GP1__E
        vInfHat_2 = op.vInfHat_GP1__E
        self.assertIs(vInfHat_1, vInfHat_2)

        vInf_1 = op.vInf_GP1__E
        vInf_2 = op.vInf_GP1__E
        self.assertIs(vInf_1, vInf_2)

    def test_cached_numpy_arrays_are_read_only(self) -> None:
        """Test that cached numpy arrays cannot be mutated in place."""
        op = self.basic_op

        # Test T_pas_GP1_CgP1_to_W_CgP1 is read-only
        T_forward = op.T_pas_GP1_CgP1_to_W_CgP1
        with self.assertRaises(ValueError):
            T_forward[0, 0] = 999.0

        # Test T_pas_W_CgP1_to_GP1_CgP1 is read-only
        T_inverse = op.T_pas_W_CgP1_to_GP1_CgP1
        with self.assertRaises(ValueError):
            T_inverse[0, 0] = 999.0

        # Test vInfHat_GP1__E is read-only
        vInfHat = op.vInfHat_GP1__E
        with self.assertRaises(ValueError):
            vInfHat[0] = 999.0

        # Test vInf_GP1__E is read-only
        vInf = op.vInf_GP1__E
        with self.assertRaises(ValueError):
            vInf[0] = 999.0

    def test_derived_properties_computed_correctly_after_caching(self) -> None:
        """Test that derived properties return correct values after caching."""
        # Create a fresh OperatingPoint
        op = ps.operating_point.OperatingPoint(
            rho=1.5, vCg__E=25.0, alpha=15.0, beta=5.0
        )

        # Access all derived properties to populate caches
        qInf = op.qInf__E
        T_forward = op.T_pas_GP1_CgP1_to_W_CgP1
        T_inverse = op.T_pas_W_CgP1_to_GP1_CgP1
        vInfHat = op.vInfHat_GP1__E
        vInf = op.vInf_GP1__E

        # Verify qInf__E calculation
        expected_qInf = 0.5 * 1.5 * 25.0**2
        self.assertAlmostEqual(qInf, expected_qInf, places=10)

        # Verify transformation matrices are still inverses
        identity = T_forward @ T_inverse
        npt.assert_allclose(identity, np.eye(4), atol=1e-14)

        # Verify vInfHat is unit vector
        npt.assert_allclose(np.linalg.norm(vInfHat), 1.0, atol=1e-14)

        # Verify vInf magnitude equals vCg__E
        npt.assert_allclose(np.linalg.norm(vInf), 25.0, atol=1e-14)

        # Verify vInf = vInfHat * vCg__E
        npt.assert_allclose(vInf, vInfHat * 25.0, atol=1e-14)

    def test_integer_parameters_converted_to_float(self) -> None:
        """Test that integer parameters are internally converted to floats."""
        op = self.integer_parameters_op

        # Verify all parameters are stored as floats.
        self.assertIsInstance(op.rho, float)
        self.assertIsInstance(op.vCg__E, float)
        self.assertIsInstance(op.alpha, float)
        self.assertIsInstance(op.beta, float)
        self.assertIsInstance(op.externalFX_W, float)
        self.assertIsInstance(op.nu, float)

        # Verify values are correct.
        self.assertEqual(op.rho, 1.0)
        self.assertEqual(op.vCg__E, 10.0)
        self.assertEqual(op.alpha, 5.0)
        self.assertEqual(op.beta, 0.0)
        self.assertEqual(op.externalFX_W, 0.0)
        self.assertEqual(op.nu, 1.0)

    def test_negative_beta_transformation(self) -> None:
        """Test transformation with negative beta."""
        op = self.negative_beta_op

        # Verify beta is stored correctly.
        self.assertEqual(op.beta, -15.0)

        # Verify transformation matrix is still valid.
        T = op.T_pas_GP1_CgP1_to_W_CgP1
        R = T[:3, :3]

        # Rotation matrix should be orthonormal.
        npt.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        npt.assert_allclose(np.linalg.det(R), 1.0, atol=1e-14)

        # Transformation matrices should be inverses.
        T_inverse = op.T_pas_W_CgP1_to_GP1_CgP1
        identity = T @ T_inverse
        npt.assert_allclose(identity, np.eye(4), atol=1e-14)

    def test_boundary_beta_transformation(self) -> None:
        """Test transformation with beta at boundary (180 degrees)."""
        op = self.boundary_beta_op

        # Verify beta is stored correctly.
        self.assertEqual(op.beta, 180.0)

        # Verify transformation matrix is valid.
        T = op.T_pas_GP1_CgP1_to_W_CgP1
        R = T[:3, :3]

        npt.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        self.assertFalse(np.any(np.isnan(T)))
        self.assertFalse(np.any(np.isinf(T)))

    def test_combined_boundary_angles_transformation(self) -> None:
        """Test transformation with both alpha and beta at boundary values."""
        op = self.combined_boundary_angles_op

        # Verify angles are stored correctly.
        self.assertEqual(op.alpha, 180.0)
        self.assertEqual(op.beta, 180.0)

        # Verify transformation matrix is valid.
        T = op.T_pas_GP1_CgP1_to_W_CgP1
        R = T[:3, :3]

        npt.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        self.assertFalse(np.any(np.isnan(T)))
        self.assertFalse(np.any(np.isinf(T)))

        # Transformation matrices should be inverses.
        T_inverse = op.T_pas_W_CgP1_to_GP1_CgP1
        identity = T @ T_inverse
        npt.assert_allclose(identity, np.eye(4), atol=1e-14)

    def test_very_low_speed_fixture(self) -> None:
        """Test the very low speed fixture properties."""
        op = self.very_low_speed_op

        # Verify speed is stored correctly.
        self.assertEqual(op.vCg__E, 0.01)

        # Verify qInf calculation.
        expected_qInf = 0.5 * 1.225 * 0.01**2
        self.assertAlmostEqual(op.qInf__E, expected_qInf, places=14)

        # Verify velocity vector magnitude.
        npt.assert_allclose(np.linalg.norm(op.vInf_GP1__E), 0.01, atol=1e-14)

    def test_negative_external_force_fixture(self) -> None:
        """Test the negative external force fixture properties."""
        op = self.negative_external_force_op

        # Verify external force is stored correctly.
        self.assertEqual(op.externalFX_W, -25.0)
        self.assertIsInstance(op.externalFX_W, float)

    def test_near_boundary_alpha_fixture(self) -> None:
        """Test the near boundary alpha fixture properties."""
        op = self.near_boundary_alpha_op

        # Verify alpha is stored correctly (near -180).
        self.assertEqual(op.alpha, -179.999)

        # Verify transformation is still valid.
        T = op.T_pas_GP1_CgP1_to_W_CgP1
        self.assertFalse(np.any(np.isnan(T)))
        self.assertFalse(np.any(np.isinf(T)))

    def test_near_boundary_beta_fixture(self) -> None:
        """Test the near boundary beta fixture properties."""
        op = self.near_boundary_beta_op

        # Verify beta is stored correctly (near -180).
        self.assertEqual(op.beta, -179.999)

        # Verify transformation is still valid.
        T = op.T_pas_GP1_CgP1_to_W_CgP1
        self.assertFalse(np.any(np.isnan(T)))
        self.assertFalse(np.any(np.isinf(T)))

    def test_all_new_fixtures_comprehensive_properties(self) -> None:
        """Test all properties on the new fixtures."""
        fixtures = [
            (self.negative_beta_op, "negative_beta"),
            (self.boundary_beta_op, "boundary_beta"),
            (self.combined_boundary_angles_op, "combined_boundary_angles"),
            (self.very_low_speed_op, "very_low_speed"),
            (self.integer_parameters_op, "integer_parameters"),
            (self.negative_external_force_op, "negative_external_force"),
            (self.near_boundary_alpha_op, "near_boundary_alpha"),
            (self.near_boundary_beta_op, "near_boundary_beta"),
        ]

        for op, fixture_name in fixtures:
            with self.subTest(fixture=fixture_name):
                # Test all attributes exist and have correct types.
                self.assertIsInstance(op.rho, float)
                self.assertIsInstance(op.vCg__E, float)
                self.assertIsInstance(op.alpha, float)
                self.assertIsInstance(op.beta, float)
                self.assertIsInstance(op.externalFX_W, float)
                self.assertIsInstance(op.nu, float)

                # Test all properties return correct types.
                self.assertIsInstance(op.qInf__E, float)
                self.assertIsInstance(op.T_pas_GP1_CgP1_to_W_CgP1, np.ndarray)
                self.assertIsInstance(op.T_pas_W_CgP1_to_GP1_CgP1, np.ndarray)
                self.assertIsInstance(op.vInfHat_GP1__E, np.ndarray)
                self.assertIsInstance(op.vInf_GP1__E, np.ndarray)

                # Test no NaN or Inf values in properties.
                self.assertFalse(np.isnan(op.qInf__E))
                self.assertFalse(np.any(np.isnan(op.T_pas_GP1_CgP1_to_W_CgP1)))
                self.assertFalse(np.any(np.isnan(op.T_pas_W_CgP1_to_GP1_CgP1)))
                self.assertFalse(np.any(np.isnan(op.vInfHat_GP1__E)))
                self.assertFalse(np.any(np.isnan(op.vInf_GP1__E)))

    # --- Tests for angles_E_to_BP1_izyx ---

    def test_angles_E_to_BP1_izyx_default(self) -> None:
        """Test that an unset angles_E_to_BP1_izyx resolves to the level-flight
        attitude.

        With the default alpha of 5 degrees and no sideslip, the unset attitude resolves
        to a 5 degree pitch (body izyx (0, 5, 0)) so the body flies level along Earth +x
        rather than tilting the flight path relative to Earth.
        """
        op = ps.operating_point.OperatingPoint(alpha=5.0)
        npt.assert_allclose(op.angles_E_to_BP1_izyx, [0.0, 5.0, 0.0], atol=1e-12)

    def test_angles_E_to_BP1_izyx_parameter_validation_valid(self) -> None:
        """Test angles_E_to_BP1_izyx parameter validation with valid values."""
        # Test various valid angles_E_to_BP1_izyx values within (-180, 180].
        valid_angles_values: list[np.ndarray | Sequence[float]] = [
            (0.0, 0.0, 0.0),
            (45.0, 30.0, 15.0),
            (-45.0, -30.0, -15.0),
            (180.0, 180.0, 180.0),
            (-179.999, -179.999, -179.999),
            [90.0, -90.0, 0.0],
            np.array([10.0, 20.0, 30.0]),
        ]

        for angles in valid_angles_values:
            with self.subTest(angles=angles):
                op = ps.operating_point.OperatingPoint(
                    alpha=5.0,
                    angles_E_to_BP1_izyx=angles,
                )
                npt.assert_array_almost_equal(
                    op.angles_E_to_BP1_izyx,
                    np.array(angles, dtype=float),
                )

    def test_angles_E_to_BP1_izyx_parameter_validation_invalid_range(self) -> None:
        """Test angles_E_to_BP1_izyx parameter validation with values outside range."""
        # Test angles outside the valid range (-180, 180].
        invalid_angles_values = [
            (180.1, 0.0, 0.0),
            (0.0, 180.1, 0.0),
            (0.0, 0.0, 180.1),
            (-180.0, 0.0, 0.0),
            (0.0, -180.0, 0.0),
            (0.0, 0.0, -180.0),
            (-180.1, 0.0, 0.0),
            (200.0, 200.0, 200.0),
        ]

        for invalid_angles in invalid_angles_values:
            with self.subTest(invalid_angles=invalid_angles):
                with self.assertRaises(ValueError):
                    ps.operating_point.OperatingPoint(
                        alpha=5.0,
                        angles_E_to_BP1_izyx=invalid_angles,
                    )

    def test_angles_E_to_BP1_izyx_parameter_validation_invalid_type(self) -> None:
        """Test angles_E_to_BP1_izyx parameter validation with invalid types."""
        # Test invalid types for angles_E_to_BP1_izyx. None is excluded because it is
        # the valid default that resolves to the level-flight attitude.
        invalid_angles_values: list[Any] = [
            (0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0),
            "invalid",
            (1.0, "invalid", 0.0),
        ]

        for invalid_angles in invalid_angles_values:
            with self.subTest(invalid_angles=invalid_angles):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(
                        alpha=5.0,
                        angles_E_to_BP1_izyx=invalid_angles,
                    )

    def test_angles_E_to_BP1_izyx_shape_and_type(self) -> None:
        """Test angles_E_to_BP1_izyx shape and type."""
        angles = self.with_attitude_angles_op.angles_E_to_BP1_izyx

        # Should be a 3-element ndarray of floats.
        self.assertEqual(len(angles), 3)
        self.assertIsInstance(angles, np.ndarray)
        self.assertEqual(angles.dtype, float)

    def test_angles_E_to_BP1_izyx_conversion_to_float_array(self) -> None:
        """Test that angles_E_to_BP1_izyx is converted to a float array."""
        # Test with integer values.
        op = ps.operating_point.OperatingPoint(
            alpha=5.0,
            angles_E_to_BP1_izyx=(10, 20, 30),
        )
        self.assertEqual(op.angles_E_to_BP1_izyx.dtype, float)
        npt.assert_array_equal(op.angles_E_to_BP1_izyx, [10.0, 20.0, 30.0])

    def test_angles_E_to_BP1_izyx_immutable(self) -> None:
        """Test that angles_E_to_BP1_izyx is read only."""
        op = self.with_attitude_angles_op

        # The property should not be settable.
        with self.assertRaises(AttributeError):
            setattr(op, "angles_E_to_BP1_izyx", np.array([0.0, 0.0, 0.0]))

        # The underlying array should not be writable.
        with self.assertRaises(ValueError):
            op.angles_E_to_BP1_izyx[0] = 999.0

    # --- Tests for CgP1_E_Eo ---

    def test_CgP1_E_Eo_default(self) -> None:
        """Test that CgP1_E_Eo defaults to (0, 0, 0)."""
        op = ps.operating_point.OperatingPoint(alpha=5.0)
        npt.assert_array_equal(op.CgP1_E_Eo, [0.0, 0.0, 0.0])

    def test_CgP1_E_Eo_parameter_validation_valid(self) -> None:
        """Test CgP1_E_Eo parameter validation with valid values."""
        valid_cg_values: list[np.ndarray | Sequence[float]] = [
            (0.0, 0.0, 0.0),
            (100.0, 200.0, -50.0),
            (-10.0, -20.0, -30.0),
            [1.0, 2.0, 3.0],
            np.array([4.0, 5.0, 6.0]),
        ]

        for cg in valid_cg_values:
            with self.subTest(cg=cg):
                op = ps.operating_point.OperatingPoint(alpha=5.0, CgP1_E_Eo=cg)
                npt.assert_array_almost_equal(op.CgP1_E_Eo, np.array(cg, dtype=float))

    def test_CgP1_E_Eo_parameter_validation_invalid(self) -> None:
        """Test CgP1_E_Eo parameter validation with invalid values."""
        invalid_cg_values: list[Any] = [
            (0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0),
            "invalid",
            None,
            (1.0, "invalid", 0.0),
        ]

        for invalid_cg in invalid_cg_values:
            with self.subTest(invalid_cg=invalid_cg):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.operating_point.OperatingPoint(alpha=5.0, CgP1_E_Eo=invalid_cg)

    def test_CgP1_E_Eo_shape_and_type(self) -> None:
        """Test CgP1_E_Eo shape and type."""
        cg = self.with_cg_position_op.CgP1_E_Eo

        # Should be a 3-element ndarray of floats.
        self.assertEqual(len(cg), 3)
        self.assertIsInstance(cg, np.ndarray)
        self.assertEqual(cg.dtype, float)

    def test_CgP1_E_Eo_conversion_to_float_array(self) -> None:
        """Test that CgP1_E_Eo is converted to a float array."""
        op = ps.operating_point.OperatingPoint(alpha=5.0, CgP1_E_Eo=(10, 20, 30))
        self.assertEqual(op.CgP1_E_Eo.dtype, float)
        npt.assert_array_equal(op.CgP1_E_Eo, [10.0, 20.0, 30.0])

    def test_CgP1_E_Eo_immutable(self) -> None:
        """Test that CgP1_E_Eo is read only."""
        op = self.with_cg_position_op

        # The property should not be settable.
        with self.assertRaises(AttributeError):
            setattr(op, "CgP1_E_Eo", np.array([0.0, 0.0, 0.0]))

        # The underlying array should not be writable.
        with self.assertRaises(ValueError):
            op.CgP1_E_Eo[0] = 999.0

    # --- Tests for surfaceNormal_E and surfacePoint_E_Eo ---

    def test_surface_parameters_default_to_none(self) -> None:
        """Test that surfaceNormal_E and surfacePoint_E_Eo default to None."""
        op = ps.operating_point.OperatingPoint(alpha=5.0)
        self.assertIsNone(op.surfaceNormal_E)
        self.assertIsNone(op.surfacePoint_E_Eo)

    def test_surface_parameters_both_provided(self) -> None:
        """Test that providing both surface parameters succeeds."""
        op = ps.operating_point.OperatingPoint(
            alpha=5.0,
            surfaceNormal_E=(0.0, 0.0, 1.0),
            surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        )
        assert op.surfaceNormal_E is not None
        npt.assert_allclose(op.surfaceNormal_E, [0.0, 0.0, 1.0], atol=1e-14)
        npt.assert_array_equal(op.surfacePoint_E_Eo, [0.0, 0.0, 0.0])

    def test_surface_parameters_mutual_exclusivity(self) -> None:
        """Test that providing only one surface parameter raises ValueError."""
        # surfaceNormal_E without surfacePoint_E_Eo.
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=(0.0, 0.0, 1.0),
                surfacePoint_E_Eo=None,
            )

        # surfacePoint_E_Eo without surfaceNormal_E.
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=None,
                surfacePoint_E_Eo=(0.0, 0.0, 0.0),
            )

    def test_surfaceNormal_E_is_normalized(self) -> None:
        """Test that surfaceNormal_E is normalized to a unit vector."""
        # Provide a non unit normal; it should be normalized internally.
        op = ps.operating_point.OperatingPoint(
            alpha=5.0,
            surfaceNormal_E=(0.0, 0.0, 3.0),
            surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        )
        assert op.surfaceNormal_E is not None
        npt.assert_allclose(np.linalg.norm(op.surfaceNormal_E), 1.0, atol=1e-14)
        npt.assert_allclose(op.surfaceNormal_E, [0.0, 0.0, 1.0], atol=1e-14)

    def test_surfaceNormal_E_validation_invalid(self) -> None:
        """Test surfaceNormal_E validation with invalid values."""
        # Zero vector should be rejected.
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=(0.0, 0.0, 0.0),
                surfacePoint_E_Eo=(0.0, 0.0, 0.0),
            )

        # Wrong length.
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=(0.0, 1.0),
                surfacePoint_E_Eo=(0.0, 0.0, 0.0),
            )

        # Non numeric.
        bad_surfaceNormal_E: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=bad_surfaceNormal_E,
                surfacePoint_E_Eo=(0.0, 0.0, 0.0),
            )

    def test_surfacePoint_E_Eo_validation_invalid(self) -> None:
        """Test surfacePoint_E_Eo validation with invalid values."""
        # Wrong length.
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=(0.0, 0.0, 1.0),
                surfacePoint_E_Eo=(0.0, 0.0),
            )

        # Non numeric.
        bad_surfacePoint_E_Eo: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(
                alpha=5.0,
                surfaceNormal_E=(0.0, 0.0, 1.0),
                surfacePoint_E_Eo=bad_surfacePoint_E_Eo,
            )

    def test_surface_parameters_immutable(self) -> None:
        """Test that surfaceNormal_E and surfacePoint_E_Eo are read only."""
        op = self.with_ground_surface_op

        # The properties should not be settable.
        with self.assertRaises(AttributeError):
            setattr(op, "surfaceNormal_E", np.array([1.0, 0.0, 0.0]))
        with self.assertRaises(AttributeError):
            setattr(op, "surfacePoint_E_Eo", np.array([1.0, 0.0, 0.0]))

        # The underlying arrays should not be writable.
        assert op.surfaceNormal_E is not None
        assert op.surfacePoint_E_Eo is not None
        with self.assertRaises(ValueError):
            op.surfaceNormal_E[0] = 999.0
        with self.assertRaises(ValueError):
            op.surfacePoint_E_Eo[0] = 999.0

    def test_surface_parameters_accept_various_array_likes(self) -> None:
        """Test that surface parameters accept tuples, lists, and ndarrays."""
        input_types: list[
            tuple[np.ndarray | Sequence[float], np.ndarray | Sequence[float]]
        ] = [
            ((0.0, 0.0, 1.0), (1.0, 2.0, 3.0)),
            ([0.0, 0.0, 1.0], [1.0, 2.0, 3.0]),
            (np.array([0.0, 0.0, 1.0]), np.array([1.0, 2.0, 3.0])),
        ]

        for normal, point in input_types:
            with self.subTest(normal_type=type(normal).__name__):
                op = ps.operating_point.OperatingPoint(
                    alpha=5.0,
                    surfaceNormal_E=normal,
                    surfacePoint_E_Eo=point,
                )
                assert op.surfaceNormal_E is not None
                npt.assert_allclose(op.surfaceNormal_E, [0.0, 0.0, 1.0], atol=1e-14)
                npt.assert_array_equal(op.surfacePoint_E_Eo, [1.0, 2.0, 3.0])

    # --- Tests for intermediate transformation properties (GP1 <-> BP1, BP1 <-> W) ---

    def test_T_pas_GP1_CgP1_to_BP1_CgP1_shape_and_type(self) -> None:
        """Test T_pas_GP1_CgP1_to_BP1_CgP1 shape and type."""
        T = self.basic_op.T_pas_GP1_CgP1_to_BP1_CgP1

        self.assertEqual(T.shape, (4, 4))
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_GP1_CgP1_to_BP1_CgP1_is_180_deg_about_y(self) -> None:
        """Test that GP1 to BP1 is a 180 degree rotation about y.

        Geometry axes: +x aft, +y right, +z up. Body axes: +x forward, +y right, +z
        down. So x and z flip sign; y is unchanged.
        """
        T = self.basic_op.T_pas_GP1_CgP1_to_BP1_CgP1
        R = T[:3, :3]

        expected_R = np.diag([-1.0, 1.0, -1.0])
        npt.assert_allclose(R, expected_R, atol=1e-14)

        # Translation part should be zero.
        npt.assert_allclose(T[:3, 3], [0.0, 0.0, 0.0], atol=1e-14)

    def test_T_pas_GP1_CgP1_to_BP1_CgP1_and_inverse_are_inverses(self) -> None:
        """Test that GP1 to BP1 and BP1 to GP1 are inverses."""
        T_forward = self.basic_op.T_pas_GP1_CgP1_to_BP1_CgP1
        T_inverse = self.basic_op.T_pas_BP1_CgP1_to_GP1_CgP1

        npt.assert_allclose(T_forward @ T_inverse, np.eye(4), atol=1e-14)
        npt.assert_allclose(T_inverse @ T_forward, np.eye(4), atol=1e-14)

    def test_T_pas_BP1_CgP1_to_W_CgP1_shape_and_type(self) -> None:
        """Test T_pas_BP1_CgP1_to_W_CgP1 shape and type."""
        T = self.basic_op.T_pas_BP1_CgP1_to_W_CgP1

        self.assertEqual(T.shape, (4, 4))
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_BP1_CgP1_to_W_CgP1_and_inverse_are_inverses(self) -> None:
        """Test that BP1 to W and W to BP1 are inverses."""
        fixtures = [
            self.basic_op,
            self.zero_alpha_beta_op,
            self.high_alpha_op,
            self.nonzero_beta_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_forward = op.T_pas_BP1_CgP1_to_W_CgP1
                T_inverse = op.T_pas_W_CgP1_to_BP1_CgP1

                npt.assert_allclose(T_forward @ T_inverse, np.eye(4), atol=1e-14)

    def test_GP1_to_W_composition_equals_product(self) -> None:
        """Test that GP1 to W equals GP1 to BP1 composed with BP1 to W."""
        fixtures = [
            self.basic_op,
            self.zero_alpha_beta_op,
            self.high_alpha_op,
            self.nonzero_beta_op,
            self.boundary_alpha_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_GP1_to_W = op.T_pas_GP1_CgP1_to_W_CgP1
                T_GP1_to_BP1 = op.T_pas_GP1_CgP1_to_BP1_CgP1
                T_BP1_to_W = op.T_pas_BP1_CgP1_to_W_CgP1

                npt.assert_allclose(T_GP1_to_W, T_BP1_to_W @ T_GP1_to_BP1, atol=1e-14)

    def test_intermediate_transformations_read_only(self) -> None:
        """Test that intermediate transformation matrices are read only."""
        op = self.basic_op

        with self.assertRaises(ValueError):
            op.T_pas_GP1_CgP1_to_BP1_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_BP1_CgP1_to_GP1_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_BP1_CgP1_to_W_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_W_CgP1_to_BP1_CgP1[0, 0] = 999.0

    def test_intermediate_transformations_cached(self) -> None:
        """Test that intermediate transformation properties return the same objects on
        repeated access."""
        op = self.basic_op

        self.assertIs(op.T_pas_GP1_CgP1_to_BP1_CgP1, op.T_pas_GP1_CgP1_to_BP1_CgP1)
        self.assertIs(op.T_pas_BP1_CgP1_to_GP1_CgP1, op.T_pas_BP1_CgP1_to_GP1_CgP1)
        self.assertIs(op.T_pas_BP1_CgP1_to_W_CgP1, op.T_pas_BP1_CgP1_to_W_CgP1)
        self.assertIs(op.T_pas_W_CgP1_to_BP1_CgP1, op.T_pas_W_CgP1_to_BP1_CgP1)

    # --- Tests for Earth transformation properties ---

    def test_T_pas_E_CgP1_to_BP1_CgP1_shape_and_type(self) -> None:
        """Test T_pas_E_CgP1_to_BP1_CgP1 shape and type."""
        T = self.with_attitude_angles_op.T_pas_E_CgP1_to_BP1_CgP1

        self.assertEqual(T.shape, (4, 4))
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_E_CgP1_to_BP1_CgP1_identity_when_zero_angles(self) -> None:
        """Test that E to BP1 is identity when angles_E_to_BP1_izyx is all zeros."""
        op = self.basic_op
        T = op.T_pas_E_CgP1_to_BP1_CgP1

        npt.assert_allclose(T, np.eye(4), atol=1e-14)

    def test_T_pas_E_CgP1_to_BP1_CgP1_and_inverse_are_inverses(self) -> None:
        """Test that E to BP1 and BP1 to E are inverses."""
        fixtures = [
            self.basic_op,
            self.with_attitude_angles_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_forward = op.T_pas_E_CgP1_to_BP1_CgP1
                T_inverse = op.T_pas_BP1_CgP1_to_E_CgP1

                npt.assert_allclose(T_forward @ T_inverse, np.eye(4), atol=1e-14)
                npt.assert_allclose(T_inverse @ T_forward, np.eye(4), atol=1e-14)

    def test_T_pas_E_CgP1_to_BP1_CgP1_orthonormal_rotation(self) -> None:
        """Test that the rotation part of E to BP1 is orthonormal."""
        T = self.with_attitude_angles_op.T_pas_E_CgP1_to_BP1_CgP1
        R = T[:3, :3]

        npt.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        npt.assert_allclose(np.linalg.det(R), 1.0, atol=1e-14)

    def test_T_pas_E_CgP1_to_GP1_CgP1_shape_and_type(self) -> None:
        """Test T_pas_E_CgP1_to_GP1_CgP1 shape and type."""
        T = self.with_attitude_angles_op.T_pas_E_CgP1_to_GP1_CgP1

        self.assertEqual(T.shape, (4, 4))
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_E_CgP1_to_GP1_CgP1_and_inverse_are_inverses(self) -> None:
        """Test that E to GP1 and GP1 to E are inverses."""
        fixtures = [
            self.basic_op,
            self.with_attitude_angles_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_forward = op.T_pas_E_CgP1_to_GP1_CgP1
                T_inverse = op.T_pas_GP1_CgP1_to_E_CgP1

                npt.assert_allclose(T_forward @ T_inverse, np.eye(4), atol=1e-14)
                npt.assert_allclose(T_inverse @ T_forward, np.eye(4), atol=1e-14)

    def test_E_to_GP1_is_180_deg_about_y_when_zero_angles(self) -> None:
        """Test that E to GP1 equals the 180 degree about y rotation when
        angles_E_to_BP1_izyx is all zeros.

        When E to BP1 is identity, E to GP1 is just BP1 to GP1, which is the 180 degree
        rotation about y (flipping x and z).
        """
        op = self.basic_op
        T = op.T_pas_E_CgP1_to_GP1_CgP1
        R = T[:3, :3]

        expected_R = np.diag([-1.0, 1.0, -1.0])
        npt.assert_allclose(R, expected_R, atol=1e-14)

    def test_E_to_GP1_composition_equals_product(self) -> None:
        """Test that E to GP1 equals E to BP1 composed with BP1 to GP1."""
        fixtures = [
            self.basic_op,
            self.with_attitude_angles_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_E_to_GP1 = op.T_pas_E_CgP1_to_GP1_CgP1
                T_E_to_BP1 = op.T_pas_E_CgP1_to_BP1_CgP1
                T_BP1_to_GP1 = op.T_pas_BP1_CgP1_to_GP1_CgP1

                npt.assert_allclose(T_E_to_GP1, T_BP1_to_GP1 @ T_E_to_BP1, atol=1e-14)

    def test_T_pas_W_CgP1_to_E_CgP1_shape_and_type(self) -> None:
        """Test T_pas_W_CgP1_to_E_CgP1 shape and type."""
        T = self.with_attitude_angles_op.T_pas_W_CgP1_to_E_CgP1

        self.assertEqual(T.shape, (4, 4))
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.dtype, float)

    def test_T_pas_W_CgP1_to_E_CgP1_and_inverse_are_inverses(self) -> None:
        """Test that W to E and E to W are inverses."""
        fixtures = [
            self.basic_op,
            self.with_attitude_angles_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_forward = op.T_pas_W_CgP1_to_E_CgP1
                T_inverse = op.T_pas_E_CgP1_to_W_CgP1

                npt.assert_allclose(T_forward @ T_inverse, np.eye(4), atol=1e-14)
                npt.assert_allclose(T_inverse @ T_forward, np.eye(4), atol=1e-14)

    def test_W_to_E_composition_equals_product(self) -> None:
        """Test that W to E equals W to BP1 composed with BP1 to E."""
        fixtures = [
            self.basic_op,
            self.with_attitude_angles_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                T_W_to_E = op.T_pas_W_CgP1_to_E_CgP1
                T_W_to_BP1 = op.T_pas_W_CgP1_to_BP1_CgP1
                T_BP1_to_E = op.T_pas_BP1_CgP1_to_E_CgP1

                npt.assert_allclose(T_W_to_E, T_BP1_to_E @ T_W_to_BP1, atol=1e-14)

    def test_earth_transformations_read_only(self) -> None:
        """Test that Earth transformation matrices are read only."""
        op = self.with_attitude_angles_op

        with self.assertRaises(ValueError):
            op.T_pas_E_CgP1_to_BP1_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_BP1_CgP1_to_E_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_E_CgP1_to_GP1_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_GP1_CgP1_to_E_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_W_CgP1_to_E_CgP1[0, 0] = 999.0
        with self.assertRaises(ValueError):
            op.T_pas_E_CgP1_to_W_CgP1[0, 0] = 999.0

    def test_earth_transformations_cached(self) -> None:
        """Test that Earth transformation properties return the same objects on repeated
        access."""
        op = self.with_attitude_angles_op

        self.assertIs(op.T_pas_E_CgP1_to_BP1_CgP1, op.T_pas_E_CgP1_to_BP1_CgP1)
        self.assertIs(op.T_pas_BP1_CgP1_to_E_CgP1, op.T_pas_BP1_CgP1_to_E_CgP1)
        self.assertIs(op.T_pas_E_CgP1_to_GP1_CgP1, op.T_pas_E_CgP1_to_GP1_CgP1)
        self.assertIs(op.T_pas_GP1_CgP1_to_E_CgP1, op.T_pas_GP1_CgP1_to_E_CgP1)
        self.assertIs(op.T_pas_W_CgP1_to_E_CgP1, op.T_pas_W_CgP1_to_E_CgP1)
        self.assertIs(op.T_pas_E_CgP1_to_W_CgP1, op.T_pas_E_CgP1_to_W_CgP1)

    # --- Tests for derived surface properties ---

    def test_surfaceNormal_GP1_is_none_when_no_surface(self) -> None:
        """Test that surfaceNormal_GP1 returns None when no surface is defined."""
        self.assertIsNone(self.basic_op.surfaceNormal_GP1)

    def test_surfacePoint_GP1_CgP1_is_none_when_no_surface(self) -> None:
        """Test that surfacePoint_GP1_CgP1 returns None when no surface is defined."""
        self.assertIsNone(self.basic_op.surfacePoint_GP1_CgP1)

    def test_surfaceNormal_GP1_is_unit_vector(self) -> None:
        """Test that surfaceNormal_GP1 is a unit vector when a surface is defined."""
        fixtures = [
            self.with_ground_surface_op,
            self.with_tilted_surface_op,
        ]

        for op in fixtures:
            with self.subTest(op=op):
                normal = op.surfaceNormal_GP1
                self.assertIsNotNone(normal)
                assert normal is not None
                npt.assert_allclose(np.linalg.norm(normal), 1.0, atol=1e-14)

    def test_surfaceNormal_GP1_with_zero_angles(self) -> None:
        """Test surfaceNormal_GP1 when angles_E_to_BP1_izyx is all zeros.

        The ground surface fixture has surfaceNormal_E = (0, 0, -1) (down in
        Earth axes) and angles_E_to_BP1_izyx = (0, 0, 0). With zero attitude
        angles, E to GP1 is the 180 degree about y rotation, which flips x and
        z. So the Earth z down vector (0, 0, -1) should map to GP1 z up
        (0, 0, 1).
        """
        normal = self.with_ground_surface_op.surfaceNormal_GP1
        assert normal is not None
        npt.assert_allclose(normal, [0.0, 0.0, 1.0], atol=1e-14)

    def test_surfacePoint_GP1_CgP1_with_zero_angles(self) -> None:
        """Test surfacePoint_GP1_CgP1 when angles_E_to_BP1_izyx is all zeros.

        The ground surface fixture has surfacePoint_E_Eo = (0, 0, 0) and
        CgP1_E_Eo = (0, 0, -10). The surface point relative to CG in Earth axes
        is (0, 0, 0) - (0, 0, -10) = (0, 0, 10). The E to GP1 rotation (180
        degrees about y) flips x and z, giving (0, 0, -10) in GP1 axes. This
        means the ground is 10 meters below the CG in GP1 axes (z up).
        """
        point = self.with_ground_surface_op.surfacePoint_GP1_CgP1
        assert point is not None
        npt.assert_allclose(point, [0.0, 0.0, -10.0], atol=1e-14)

    def test_surfaceNormal_GP1_translation_does_not_affect_normal(self) -> None:
        """Test that changing CgP1_E_Eo does not affect surfaceNormal_GP1.

        The normal is a non-position vector, so it should be independent of the CG
        position.
        """
        op_no_offset = ps.operating_point.OperatingPoint(
            alpha=5.0,
            CgP1_E_Eo=(0.0, 0.0, 0.0),
            surfaceNormal_E=(0.0, 0.0, -1.0),
            surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        )
        op_with_offset = ps.operating_point.OperatingPoint(
            alpha=5.0,
            CgP1_E_Eo=(100.0, 200.0, -50.0),
            surfaceNormal_E=(0.0, 0.0, -1.0),
            surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        )

        assert op_no_offset.surfaceNormal_GP1 is not None
        assert op_with_offset.surfaceNormal_GP1 is not None
        npt.assert_allclose(
            op_no_offset.surfaceNormal_GP1,
            op_with_offset.surfaceNormal_GP1,
            atol=1e-14,
        )

    def test_surfacePoint_GP1_CgP1_depends_on_cg_position(self) -> None:
        """Test that surfacePoint_GP1_CgP1 changes with CgP1_E_Eo.

        With the same surface in Earth axes, changing the CG position should change
        where the surface point is relative to the CG in GP1 axes.
        """
        op_near = ps.operating_point.OperatingPoint(
            alpha=5.0,
            angles_E_to_BP1_izyx=(0.0, 0.0, 0.0),
            CgP1_E_Eo=(0.0, 0.0, -5.0),
            surfaceNormal_E=(0.0, 0.0, -1.0),
            surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        )
        op_far = ps.operating_point.OperatingPoint(
            alpha=5.0,
            angles_E_to_BP1_izyx=(0.0, 0.0, 0.0),
            CgP1_E_Eo=(0.0, 0.0, -20.0),
            surfaceNormal_E=(0.0, 0.0, -1.0),
            surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        )

        # With zero attitude angles, E to GP1 flips z. The surface point relative to CG
        # is (0, 0, -Cg_z) in Earth, which maps to (0, 0, Cg_z) in GP1. But Cg_z is
        # negative (above ground), so the GP1 z component should be negative (below CG
        # in GP1's z up frame).
        point_near = op_near.surfacePoint_GP1_CgP1
        point_far = op_far.surfacePoint_GP1_CgP1
        assert point_near is not None
        assert point_far is not None

        # The airplane further from the ground should have a more negative z component
        # in GP1.
        self.assertLess(point_far[2], point_near[2])

        # Verify specific values. E to GP1 flips x and z. Near: surfacePoint_E_CgP1 =
        # (0,0,0) - (0,0,-5) = (0,0,5). In GP1: (0,0,-5).
        npt.assert_allclose(point_near, [0.0, 0.0, -5.0], atol=1e-14)

        # Far: surfacePoint_E_CgP1 = (0,0,0) - (0,0,-20) = (0,0,20). In GP1: (0,0,-20).
        npt.assert_allclose(point_far, [0.0, 0.0, -20.0], atol=1e-14)

    def test_derived_surface_properties_read_only(self) -> None:
        """Test that derived surface property arrays are read only."""
        op = self.with_ground_surface_op

        assert op.surfaceNormal_GP1 is not None
        assert op.surfacePoint_GP1_CgP1 is not None
        with self.assertRaises(ValueError):
            op.surfaceNormal_GP1[0] = 999.0
        with self.assertRaises(ValueError):
            op.surfacePoint_GP1_CgP1[0] = 999.0

    def test_derived_surface_properties_cached(self) -> None:
        """Test that derived surface properties return the same objects on repeated
        access."""
        op = self.with_ground_surface_op

        self.assertIs(op.surfaceNormal_GP1, op.surfaceNormal_GP1)
        self.assertIs(op.surfacePoint_GP1_CgP1, op.surfacePoint_GP1_CgP1)

    def test_derived_surface_properties_with_tilted_surface(self) -> None:
        """Test derived surface properties with non zero attitude angles.

        The tilted surface fixture has angles_E_to_BP1_izyx = (0, 10, 0),
        which is a 10 degree pitch (about the y axis). The surface normal in
        Earth axes is (0, 0, -1). The combined E to GP1 transformation is
        BP1_to_GP1 (180 degrees about y) composed with E_to_BP1 (10 degrees
        about y). Both rotations are pure y rotations, so they combine to 190
        degrees about y. A pure y rotation by theta maps (0, 0, -1) to
        (sin(theta), 0, -cos(theta)).
        """
        op = self.with_tilted_surface_op
        normal = op.surfaceNormal_GP1
        assert normal is not None

        theta = 190.0
        thetaRad = np.deg2rad(theta)
        expected_normal = np.array([np.sin(thetaRad), 0.0, -np.cos(thetaRad)])
        npt.assert_allclose(normal, expected_normal, atol=1e-14)

        # The normal should still be a unit vector.
        npt.assert_allclose(np.linalg.norm(normal), 1.0, atol=1e-14)

    def test_surfaceReflect_T_act_GP1_CgP1_is_none_when_no_surface(self) -> None:
        """Test that surfaceReflect_T_act_GP1_CgP1 returns None when no surface is
        defined."""
        self.assertIsNone(self.basic_op.surfaceReflect_T_act_GP1_CgP1)

    def test_surfaceReflect_T_act_GP1_CgP1_shape_and_type(self) -> None:
        """Test surfaceReflect_T_act_GP1_CgP1 shape and type."""
        op = self.with_ground_surface_op
        T = op.surfaceReflect_T_act_GP1_CgP1

        self.assertIsNotNone(T)
        assert T is not None
        self.assertIsInstance(T, np.ndarray)
        self.assertEqual(T.shape, (4, 4))
        self.assertEqual(T.dtype, float)

    def test_surfaceReflect_T_act_GP1_CgP1_reflects_point(self) -> None:
        """Test that the reflection matrix correctly reflects a point across the
        image surface.

        The ground surface fixture has a horizontal ground at z = 0 in Earth
        axes with CG at (0, 0, -10). With zero attitude angles, the GP1 to E
        transformation is a 180 degree rotation about y (x and z flip). So
        CgP1 at (0, 0, -10) in Earth becomes (0, 0, 10) in GP1, and the
        ground at z = 0 in Earth becomes z = 0 in GP1 (since the origin maps
        to the origin under 180 degree y rotation). A point at (1, 2, 3) in
        GP1_CgP1 (3 meters above CgP1 in GP1 z, which is 3 meters below CgP1
        in Earth z) should reflect to (1, 2, -23) in GP1_CgP1 when reflected
        about the ground plane.
        """
        from pterasoftware import _transformations

        op = self.with_ground_surface_op
        T = op.surfaceReflect_T_act_GP1_CgP1
        assert T is not None

        # A test point in GP1_CgP1.
        point = np.array([1.0, 2.0, 3.0])

        reflected = _transformations.apply_T_to_vectors(T, point, is_position=True)

        # The surface point in GP1_CgP1 is at z = -10.0 (10.0 meters below CgP1 in the
        # GP1 z direction, corresponding to the ground at z = 0.0 in Earth axes). The
        # surface normal in GP1 axes is (0.0, 0.0, 1.0) (a 180 degree rotation about the
        # y axis flips Earth's (0.0, 0.0, -1.0) to (0.0, 0.0, 1.0)). Reflecting (1.0,
        # 2.0, 3.0) across a plane at z = -10.0 with normal (0.0, 0.0, 1.0) gives (1.0,
        # 2.0, -23.0).
        expected = np.array([1.0, 2.0, -23.0])
        npt.assert_allclose(reflected, expected, atol=1e-12)

    def test_surfaceReflect_T_act_GP1_CgP1_reflects_velocity(self) -> None:
        """Test that the reflection matrix correctly reflects a non-position vector.

        Non-position vector reflection should only negate the component along the
        surface normal, with no translational contribution.
        """
        from pterasoftware import _transformations

        op = self.with_ground_surface_op
        T = op.surfaceReflect_T_act_GP1_CgP1
        assert T is not None

        # The surface normal in GP1 is (0, 0, 1) for this fixture. Reflecting a velocity
        # vector should negate only the z component.
        velocity = np.array([5.0, -3.0, 7.0])
        reflected = _transformations.apply_T_to_vectors(T, velocity, is_position=False)

        expected = np.array([5.0, -3.0, -7.0])
        npt.assert_allclose(reflected, expected, atol=1e-12)

    def test_surfaceReflect_T_act_GP1_CgP1_is_involution(self) -> None:
        """Test that reflecting twice returns the original point."""
        from pterasoftware import _transformations

        op = self.with_ground_surface_op
        T = op.surfaceReflect_T_act_GP1_CgP1
        assert T is not None

        point = np.array([1.0, 2.0, 3.0])
        once = _transformations.apply_T_to_vectors(T, point, is_position=True)
        twice = _transformations.apply_T_to_vectors(T, once, is_position=True)

        npt.assert_allclose(twice, point, atol=1e-12)

    def test_surfaceReflect_T_act_GP1_CgP1_read_only(self) -> None:
        """Test that the reflection matrix is read only."""
        op = self.with_ground_surface_op
        T = op.surfaceReflect_T_act_GP1_CgP1
        assert T is not None

        with self.assertRaises(ValueError):
            T[0, 0] = 999.0

    def test_surfaceReflect_T_act_GP1_CgP1_cached(self) -> None:
        """Test that the reflection matrix returns the same object on repeated
        access."""
        op = self.with_ground_surface_op
        self.assertIs(
            op.surfaceReflect_T_act_GP1_CgP1, op.surfaceReflect_T_act_GP1_CgP1
        )

    def test_surfaceReflect_T_act_GP1_CgP1_with_tilted_surface(self) -> None:
        """Test reflection matrix with non zero attitude angles.

        The tilted surface fixture has angles_E_to_BP1_izyx = (0, 10, 0) and
        CG at (50, 0, -20) in Earth axes. The surface is at z = 0 in Earth
        with normal (0, 0, -1). Reflecting a point on the surface plane
        should return the same point.
        """
        from pterasoftware import _transformations

        op = self.with_tilted_surface_op
        T = op.surfaceReflect_T_act_GP1_CgP1
        assert T is not None

        # The surface point in GP1_CgP1 should lie on the reflection plane. Reflecting
        # it should return the same point.
        surface_point = op.surfacePoint_GP1_CgP1
        assert surface_point is not None
        reflected = _transformations.apply_T_to_vectors(
            T, surface_point, is_position=True
        )
        npt.assert_allclose(reflected, surface_point, atol=1e-12)

    def test_g_E_default(self) -> None:
        """Test that g_E defaults to no gravitational field (the zero vector)."""
        op = ps.operating_point.OperatingPoint(alpha=5.0)
        npt.assert_array_equal(op.g_E, [0.0, 0.0, 0.0])

    def test_omegas_BP1__E_default(self) -> None:
        """Test that omegas_BP1__E defaults to the zero vector."""
        op = ps.operating_point.OperatingPoint(alpha=5.0)
        npt.assert_array_equal(op.omegas_BP1__E, [0.0, 0.0, 0.0])

    def test_g_E_accepts_custom_value(self) -> None:
        """Test that a non default g_E is stored as a ndarray of floats."""
        op = ps.operating_point.OperatingPoint(alpha=5.0, g_E=(1.0, -2.0, 3.5))
        self.assertIsInstance(op.g_E, np.ndarray)
        self.assertEqual(op.g_E.dtype, float)
        npt.assert_array_equal(op.g_E, [1.0, -2.0, 3.5])

    def test_omegas_BP1__E_accepts_custom_value(self) -> None:
        """Test that a non default omegas_BP1__E is stored as a ndarray of floats."""
        op = ps.operating_point.OperatingPoint(
            alpha=5.0, omegas_BP1__E=(0.1, -0.2, 0.3)
        )
        self.assertIsInstance(op.omegas_BP1__E, np.ndarray)
        self.assertEqual(op.omegas_BP1__E.dtype, float)
        npt.assert_array_equal(op.omegas_BP1__E, [0.1, -0.2, 0.3])

    def test_g_E_accepts_zero(self) -> None:
        """Test that an all zero g_E is valid (for zero gravity simulations)."""
        op = ps.operating_point.OperatingPoint(alpha=5.0, g_E=(0.0, 0.0, 0.0))
        npt.assert_array_equal(op.g_E, [0.0, 0.0, 0.0])

    def test_g_E_and_omegas_BP1__E_accept_various_array_likes(self) -> None:
        """Test that both parameters accept tuples, lists, and ndarrays."""
        array_like_pairs: list[
            tuple[np.ndarray | Sequence[float], np.ndarray | Sequence[float]]
        ] = [
            ((1.0, 2.0, 3.0), (0.1, 0.2, 0.3)),
            ([1.0, 2.0, 3.0], [0.1, 0.2, 0.3]),
            (np.array([1.0, 2.0, 3.0]), np.array([0.1, 0.2, 0.3])),
        ]

        for g, omegas in array_like_pairs:
            with self.subTest(input_type=type(g).__name__):
                op = ps.operating_point.OperatingPoint(
                    alpha=5.0, g_E=g, omegas_BP1__E=omegas
                )
                npt.assert_array_equal(op.g_E, [1.0, 2.0, 3.0])
                npt.assert_array_equal(op.omegas_BP1__E, [0.1, 0.2, 0.3])

    def test_g_E_validation_invalid(self) -> None:
        """Test g_E validation with invalid values."""
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(alpha=5.0, g_E=(0.0, 0.0))
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(alpha=5.0, g_E=(0.0, 0.0, float("nan")))
        bad_g_E: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(alpha=5.0, g_E=bad_g_E)

    def test_omegas_BP1__E_validation_invalid(self) -> None:
        """Test omegas_BP1__E validation with invalid values."""
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(alpha=5.0, omegas_BP1__E=(0.0, 0.0))
        with self.assertRaises(ValueError):
            ps.operating_point.OperatingPoint(
                alpha=5.0, omegas_BP1__E=(0.0, 0.0, float("inf"))
            )
        bad_omegas_BP1__E: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.operating_point.OperatingPoint(
                alpha=5.0, omegas_BP1__E=bad_omegas_BP1__E
            )

    def test_g_E_immutable(self) -> None:
        """Test that g_E is read only at both the property and array level."""
        op = self.basic_op
        with self.assertRaises(AttributeError):
            setattr(op, "g_E", np.array([1.0, 0.0, 0.0]))
        with self.assertRaises(ValueError):
            op.g_E[0] = 999.0

    def test_omegas_BP1__E_immutable(self) -> None:
        """Test that omegas_BP1__E is read only at both the property and array level."""
        op = self.basic_op
        with self.assertRaises(AttributeError):
            setattr(op, "omegas_BP1__E", np.array([1.0, 0.0, 0.0]))
        with self.assertRaises(ValueError):
            op.omegas_BP1__E[0] = 999.0

    def test_g_E_converts_integers_to_float(self) -> None:
        """Test that integer inputs for g_E are converted to floats."""
        op = ps.operating_point.OperatingPoint(alpha=5.0, g_E=(1, -2, 3))
        self.assertEqual(op.g_E.dtype, float)

    def test_omegas_BP1__E_converts_integers_to_float(self) -> None:
        """Test that integer inputs for omegas_BP1__E are converted to floats."""
        op = ps.operating_point.OperatingPoint(alpha=5.0, omegas_BP1__E=(1, -2, 3))
        self.assertEqual(op.omegas_BP1__E.dtype, float)


class TestOperatingPointDeepCopy(unittest.TestCase):
    """Tests for OperatingPoint.__deepcopy__ method."""

    def setUp(self) -> None:
        """Set up an OperatingPoint for deepcopy tests."""
        self.operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )

    def test_deepcopy_creates_new_instance(self) -> None:
        """Test that deepcopy creates a new, distinct OperatingPoint."""
        copied = copy.deepcopy(self.operating_point)
        self.assertIsInstance(copied, ps.operating_point.OperatingPoint)
        self.assertIsNot(copied, self.operating_point)

    def test_deepcopy_preserves_scalar_attributes(self) -> None:
        """Test that deepcopy preserves the scalar attribute values."""
        copied = copy.deepcopy(self.operating_point)
        self.assertEqual(copied.rho, self.operating_point.rho)
        self.assertEqual(copied.vCg__E, self.operating_point.vCg__E)
        self.assertEqual(copied.alpha, self.operating_point.alpha)
        self.assertEqual(copied.beta, self.operating_point.beta)
        self.assertEqual(copied.externalFX_W, self.operating_point.externalFX_W)
        self.assertEqual(copied.nu, self.operating_point.nu)

    def test_deepcopy_creates_independent_arrays(self) -> None:
        """Test that deepcopy copies the immutable arrays into new arrays."""
        copied = copy.deepcopy(self.operating_point)
        for name in ("angles_E_to_BP1_izyx", "CgP1_E_Eo", "g_E", "omegas_BP1__E"):
            original_array = getattr(self.operating_point, name)
            copied_array = getattr(copied, name)
            self.assertIsNot(copied_array, original_array)
            npt.assert_array_equal(copied_array, original_array)

    def test_deepcopy_immutable_arrays_remain_read_only(self) -> None:
        """Test that the deepcopied immutable arrays are still read only."""
        copied = copy.deepcopy(self.operating_point)
        for name in ("angles_E_to_BP1_izyx", "CgP1_E_Eo", "g_E", "omegas_BP1__E"):
            with self.assertRaises(ValueError):
                getattr(copied, name)[0] = 999.0

    def test_deepcopy_preserves_populated_caches(self) -> None:
        """Test that caches populated on the original survive the copy, read only."""
        # Populate one scalar cache and two array caches on the original.
        _ = self.operating_point.qInf__E
        _ = self.operating_point.vInf_GP1__E
        _ = self.operating_point.T_pas_GP1_CgP1_to_W_CgP1

        copied = copy.deepcopy(self.operating_point)

        self.assertEqual(
            object.__getattribute__(copied, "_qInf__E"),
            object.__getattribute__(self.operating_point, "_qInf__E"),
        )
        for slot in ("_vInf_GP1__E", "_T_pas_GP1_CgP1_to_W_CgP1"):
            copied_cache = object.__getattribute__(copied, slot)
            original_cache = object.__getattribute__(self.operating_point, slot)
            self.assertIsNotNone(copied_cache)
            npt.assert_array_equal(copied_cache, original_cache)
            self.assertFalse(copied_cache.flags.writeable)

    def test_deepcopy_leaves_uncomputed_caches_none(self) -> None:
        """Test that caches left uncomputed on the original stay None on the copy."""
        copied = copy.deepcopy(self.operating_point)
        for slot in ("_qInf__E", "_vInf_GP1__E", "_T_pas_GP1_CgP1_to_W_CgP1"):
            self.assertIsNone(object.__getattribute__(copied, slot))

    def test_deepcopy_none_surface_params_stay_none(self) -> None:
        """Test that None surface parameters remain None after deepcopy."""
        copied = copy.deepcopy(self.operating_point)
        self.assertIsNone(copied.surfaceNormal_E)
        self.assertIsNone(copied.surfacePoint_E_Eo)

    def test_deepcopy_surface_arrays_remain_read_only(self) -> None:
        """Test that populated surface arrays are copied read only."""
        operating_point = ps.operating_point.OperatingPoint(
            alpha=5.0,
            surfaceNormal_E=(0.0, 0.0, 1.0),
            surfacePoint_E_Eo=(0.0, 0.0, -1.0),
        )
        copied = copy.deepcopy(operating_point)
        self.assertIsNot(copied.surfaceNormal_E, operating_point.surfaceNormal_E)
        assert copied.surfaceNormal_E is not None
        assert copied.surfacePoint_E_Eo is not None
        with self.assertRaises(ValueError):
            copied.surfaceNormal_E[0] = 999.0
        with self.assertRaises(ValueError):
            copied.surfacePoint_E_Eo[0] = 999.0

    def test_deepcopy_handles_memo_correctly(self) -> None:
        """Test that deepcopy uses the memo dict to return the same object twice."""
        memo: dict[int, Any] = {}
        first = copy.deepcopy(self.operating_point, memo)
        second = copy.deepcopy(self.operating_point, memo)
        self.assertIs(first, second)
