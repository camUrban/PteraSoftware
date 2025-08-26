"""This module contains a class to test functions in the functions module.

This module contains the following classes:
    TestGenerateR: This class contains methods for testing the generate_R function.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps


class TestGenerateR(unittest.TestCase):
    """This class contains methods for testing the generate_R function.

    This class contains the following public methods:
        test_identity_rotations: Tests that zero angles produce identity matrices.
        test_passive_vs_active_relationship: Tests the transpose relationship between
        passive and active matrices.
        test_intrinsic_vs_extrinsic_relationship: Tests the relationship between
        intrinsic and extrinsic rotations.
        test_rotation_matrix_properties: Tests that generated matrices are proper
        rotation matrices.
        test_specific_known_active_rotations: Tests specific active rotations with
        known results.
        test_specific_known_passive_rotations: Tests specific passive rotations with
        known results.
        test_composition_properties: Tests composition of multiple rotations.
        test_large_angle_handling: Tests handling of large angles and angle wrapping.
        test_different_orders: Tests all valid Tait-Bryan rotation orders.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_identity_rotations(self):
        """Tests that zero angles produce identity matrices for all configurations.

        :return: None
        """
        zero_angles = np.array([0.0, 0.0, 0.0])

        # Test all combinations of passive/active and intrinsic/extrinsic
        for passive in [True, False]:
            for intrinsic in [True, False]:
                for order in ["123", "132", "213", "231", "312", "321"]:
                    with self.subTest(
                        passive=passive, intrinsic=intrinsic, order=order
                    ):
                        R = ps.functions.generate_R(
                            zero_angles, passive, intrinsic, order
                        )
                        npt.assert_allclose(R, np.eye(3), atol=1e-14)

    def test_passive_vs_active_relationship(self):
        """Tests that passive and active matrices are transposes of each other.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])

        for intrinsic in [True, False]:
            for order in ["123", "132", "213", "231", "312", "321"]:
                with self.subTest(intrinsic=intrinsic, order=order):
                    R_passive = ps.functions.generate_R(angles, True, intrinsic, order)
                    R_active = ps.functions.generate_R(angles, False, intrinsic, order)
                    npt.assert_allclose(R_passive, R_active.T, atol=1e-14)

    def test_intrinsic_vs_extrinsic_relationship(self):
        """Tests the relationship between intrinsic and extrinsic rotations.

        For the same angles, intrinsic rotations with order ABC should equal
        extrinsic rotations with reversed order CBA.

        :return: None
        """
        angles = np.array([20.0, 30.0, 40.0])

        # Test pairs of orders that should be equivalent
        order_pairs = [
            ("123", "321"),
            ("132", "231"),
            ("213", "312"),
            ("231", "132"),
            ("312", "213"),
            ("321", "123"),
        ]

        for passive in [True, False]:
            for intrinsic_order, extrinsic_order in order_pairs:
                with self.subTest(
                    passive=passive,
                    intrinsic_order=intrinsic_order,
                    extrinsic_order=extrinsic_order,
                ):
                    R_intrinsic = ps.functions.generate_R(
                        angles, passive, True, intrinsic_order
                    )
                    R_extrinsic = ps.functions.generate_R(
                        angles, passive, False, extrinsic_order
                    )
                    npt.assert_allclose(R_intrinsic, R_extrinsic, atol=1e-14)

    def test_rotation_matrix_properties(self):
        """Tests that generated matrices satisfy rotation matrix properties.

        Tests: determinant = 1, orthogonality (R.T @ R = I), and proper rotation.

        :return: None
        """
        # Test with various angle combinations
        test_angles = [
            np.array([0.0, 0.0, 0.0]),
            np.array([30.0, 0.0, 0.0]),
            np.array([0.0, 45.0, 0.0]),
            np.array([0.0, 0.0, 60.0]),
            np.array([30.0, 45.0, 60.0]),
            np.array([90.0, 0.0, 90.0]),
            np.array([180.0, 90.0, 45.0]),
            np.array([-30.0, -45.0, -60.0]),
        ]

        for angles in test_angles:
            for passive in [True, False]:
                for intrinsic in [True, False]:
                    for order in ["123", "321"]:  # Test representative orders
                        with self.subTest(
                            angles=angles,
                            passive=passive,
                            intrinsic=intrinsic,
                            order=order,
                        ):
                            R = ps.functions.generate_R(
                                angles, passive, intrinsic, order
                            )

                            # Test determinant = 1 (proper rotation)
                            det = np.linalg.det(R)
                            self.assertAlmostEqual(det, 1.0, places=14)

                            # Test orthogonality: R.T @ R = I
                            identity_test = R.T @ R
                            npt.assert_allclose(identity_test, np.eye(3), atol=1e-14)

                            # Test that R @ R.T = I
                            identity_test2 = R @ R.T
                            npt.assert_allclose(identity_test2, np.eye(3), atol=1e-14)

    def test_specific_known_active_rotations(self):
        """Tests specific active rotations with analytically known results.

        :return: None
        """
        # Test 90-degree rotation about x-axis (order "123", only first angle)
        angles_x90 = np.array([90.0, 0.0, 0.0])
        R_act_x90_expected = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])

        R_act_x90 = ps.functions.generate_R(angles_x90, False, True, "123")
        npt.assert_allclose(R_act_x90, R_act_x90_expected, atol=1e-14)

        # Test 90-degree rotation about y-axis (order "123", only second angle)
        angles_y90 = np.array([0.0, 90.0, 0.0])
        R_act_y90_expected = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])

        R_act_y90 = ps.functions.generate_R(angles_y90, False, True, "123")
        npt.assert_allclose(R_act_y90, R_act_y90_expected, atol=1e-14)

        # Test 90-degree rotation about z-axis (order "123", only third angle)
        angles_z90 = np.array([0.0, 0.0, 90.0])
        R_act_z90_expected = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

        R_act_z90 = ps.functions.generate_R(angles_z90, False, True, "123")
        npt.assert_allclose(R_act_z90, R_act_z90_expected, atol=1e-14)

    def test_specific_known_passive_rotations(self):
        """Tests specific passive rotations with analytically known results.

        :return: None
        """
        # Test 90-degree rotation about x-axis (order "123", only first angle)
        angles_x90 = np.array([90.0, 0.0, 0.0])
        v_A = np.array([0.0, 1.0, 0.0])
        v_B_expected = np.array([0.0, 0.0, -1.0])

        R_pas_x90 = ps.functions.generate_R(angles_x90, True, True, "123")
        v_B = R_pas_x90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

        # Test 90-degree rotation about y-axis (order "123", only second angle)
        angles_y90 = np.array([0.0, 90.0, 0.0])
        v_A = np.array([0.0, 0.0, 1.0])
        v_B_expected = np.array([-1.0, 0.0, 0.0])

        R_pas_y90 = ps.functions.generate_R(angles_y90, True, True, "123")
        v_B = R_pas_y90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

        # Test 90-degree rotation about z-axis (order "123", only third angle)
        angles_z90 = np.array([0.0, 0.0, 90.0])
        v_A = np.array([1.0, 0.0, 0.0])
        v_B_expected = np.array([0.0, -1.0, 0])

        R_pas_z90 = ps.functions.generate_R(angles_z90, True, True, "123")
        v_B = R_pas_z90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

    def test_composition_properties(self):
        """Tests composition properties of rotations.

        Tests that applying rotations in sequence equals the composed matrix.

        :return: None
        """
        # Test vector
        test_vector = np.array([1.0, 2.0, 3.0])

        # Two sets of angles
        angles1 = np.array([30.0, 0.0, 0.0])
        angles2 = np.array([0.0, 45.0, 0.0])

        # Get individual rotation matrices
        R1 = ps.functions.generate_R(angles1, False, True, "123")
        R2 = ps.functions.generate_R(angles2, False, True, "123")

        # Apply rotations sequentially
        v_rotated_sequential = R2 @ (R1 @ test_vector)

        # Apply rotations sequentially without parentheses
        v_rotated_sequential_no_parentheses = R2 @ R1 @ test_vector

        # Apply composed rotation
        R_composed = R2 @ R1
        v_rotated_composed = R_composed @ test_vector

        npt.assert_allclose(
            v_rotated_sequential, v_rotated_sequential_no_parentheses, atol=1e-14
        )
        npt.assert_allclose(v_rotated_sequential, v_rotated_composed, atol=1e-14)

    def test_large_angle_handling(self):
        """Tests handling of large angles beyond ±180 degrees.

        :return: None
        """
        # Large angles that should be equivalent to smaller ones due to periodicity
        large_angles = np.array([450.0, -270.0, 720.0])  # Equivalent to [90, 90, 0]
        equivalent_angles = np.array([90.0, 90.0, 0.0])

        for passive in [True, False]:
            for intrinsic in [True, False]:
                for order in ["123", "321"]:
                    with self.subTest(
                        passive=passive, intrinsic=intrinsic, order=order
                    ):
                        R_large = ps.functions.generate_R(
                            large_angles, passive, intrinsic, order
                        )
                        R_equivalent = ps.functions.generate_R(
                            equivalent_angles, passive, intrinsic, order
                        )

                        # Should produce the same rotation matrix
                        npt.assert_allclose(R_large, R_equivalent, atol=1e-14)

    def test_different_orders(self):
        """Tests all valid Tait-Bryan rotation orders.

        Ensures all six valid orders produce valid rotation matrices.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])
        valid_orders = ["123", "132", "213", "231", "312", "321"]

        for order in valid_orders:
            with self.subTest(order=order):
                for passive in [True, False]:
                    for intrinsic in [True, False]:
                        R = ps.functions.generate_R(angles, passive, intrinsic, order)

                        # Should be a valid rotation matrix
                        det = np.linalg.det(R)
                        self.assertAlmostEqual(det, 1.0, places=14)

                        identity_test = R.T @ R
                        npt.assert_allclose(identity_test, np.eye(3), atol=1e-14)

    def test_edge_case_angles(self):
        """Tests edge case angle values.

        :return: None
        """
        edge_case_angles = [
            np.array([0.0, 0.0, 0.0]),  # All zeros
            np.array([180.0, 0.0, 0.0]),  # 180 degree rotation
            np.array([0.0, 180.0, 0.0]),  # 180 degree rotation
            np.array([0.0, 0.0, 180.0]),  # 180 degree rotation
            np.array([180.0, 180.0, 180.0]),  # All 180 degrees
            np.array([360.0, 0.0, 0.0]),  # Full rotation
            np.array([-180.0, 0.0, 0.0]),  # Negative 180
        ]

        for angles in edge_case_angles:
            with self.subTest(angles=angles):
                for passive in [True, False]:
                    R = ps.functions.generate_R(angles, passive, True, "123")

                    # Should be a valid rotation matrix
                    det = np.linalg.det(R)
                    self.assertAlmostEqual(det, 1.0, places=14)

                    identity_test = R.T @ R
                    npt.assert_allclose(identity_test, np.eye(3), atol=1e-14)


if __name__ == "__main__":
    unittest.main()
