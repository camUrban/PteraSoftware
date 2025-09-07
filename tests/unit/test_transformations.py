"""This module contains classes to test functions in the transformations module.

This module contains the following classes:
    TestGenerateRotT: This class contains methods for testing the generate_rot_T
    function.
    TestGenerateTransT: This class contains methods for testing the generate_trans_T
    function.
    TestGenerateReflectT: This class contains methods for testing the
    generate_reflect_T function.
    TestComposeTPas: This class contains methods for testing the compose_T_pas
    function.
    TestComposeTAct: This class contains methods for testing the compose_T_act
    function.
    TestInvertTPas: This class contains methods for testing the invert_T_pas
    function.
    TestInvertTAct: This class contains methods for testing the invert_T_act
    function.
    TestApplyTToVectors: This class contains methods for testing the apply_T_to_vectors
    function.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps


class TestGenerateRotT(unittest.TestCase):
    """This class contains methods for testing the generate_rot_T function.

    This class contains the following public methods:
        test_identity_transformations: Tests that zero angles produce identity matrices.
        test_transformation_matrix_structure: Tests the structure of transformation matrices.
        test_passive_vs_active_relationship: Tests the transpose relationship between
        passive and active rotation components.
        test_intrinsic_vs_extrinsic_relationship: Tests the relationship between
        intrinsic and extrinsic rotations.
        test_rotation_matrix_properties: Tests that rotation components are proper
        rotation matrices.
        test_specific_known_active_rotations: Tests specific active rotations with
        known results.
        test_specific_known_passive_rotations: Tests specific passive rotations with
        known results.
        test_composition_properties: Tests composition of multiple rotations.
        test_large_angle_handling: Tests handling of large angles and angle wrapping.
        test_different_orders: Tests all valid Tait-Bryan rotation orders.
        test_edge_case_angles: Tests edge case angle values.
        test_homogeneous_coordinate_transformations: Tests transformation of
        homogeneous coordinates.
        test_invalid_rotation_order_rejected: Tests that passing in invalid angle
        orders raises a value error.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_identity_transformations(self):
        """Tests that zero angles produce identity matrices for all configurations.

        :return: None
        """
        zero_angles = np.array([0.0, 0.0, 0.0])

        # Test all combinations of passive/active and intrinsic/extrinsic
        for passive in [True, False]:
            for intrinsic in [True, False]:
                for order in ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]:
                    with self.subTest(
                        passive=passive, intrinsic=intrinsic, order=order
                    ):
                        T = ps.transformations.generate_rot_T(
                            zero_angles, passive, intrinsic, order
                        )
                        npt.assert_allclose(T, np.eye(4), atol=1e-14)

    def test_transformation_matrix_structure(self):
        """Tests that transformation matrices have correct 4x4 structure.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])

        for passive in [True, False]:
            for intrinsic in [True, False]:
                for order in ["xyz", "zyx"]:  # Test representative orders
                    with self.subTest(
                        passive=passive, intrinsic=intrinsic, order=order
                    ):
                        T = ps.transformations.generate_rot_T(
                            angles, passive, intrinsic, order
                        )

                        # Test output shape
                        self.assertEqual(T.shape, (4, 4))

                        # Test that translation part is zero
                        expected_translation = np.array([0.0, 0.0, 0.0])
                        npt.assert_array_equal(T[:3, 3], expected_translation)

                        # Test that bottom row is [0, 0, 0, 1]
                        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
                        npt.assert_array_equal(T[3, :], expected_bottom)

                        # Test determinant is 1
                        det = np.linalg.det(T)
                        self.assertAlmostEqual(det, 1.0, places=14)

    def test_passive_vs_active_relationship(self):
        """Tests that passive and active rotation components are transposes.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])

        for intrinsic in [True, False]:
            for order in ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]:
                with self.subTest(intrinsic=intrinsic, order=order):
                    T_passive = ps.transformations.generate_rot_T(
                        angles, True, intrinsic, order
                    )
                    T_active = ps.transformations.generate_rot_T(
                        angles, False, intrinsic, order
                    )

                    # Extract rotation parts
                    R_passive = T_passive[:3, :3]
                    R_active = T_active[:3, :3]

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
            ("xyz", "zyx"),
            ("xzy", "yzx"),
            ("yxz", "zxy"),
            ("yzx", "xzy"),
            ("zxy", "yxz"),
            ("zyx", "xyz"),
        ]

        for passive in [True, False]:
            for intrinsic_order, extrinsic_order in order_pairs:
                with self.subTest(
                    passive=passive,
                    intrinsic_order=intrinsic_order,
                    extrinsic_order=extrinsic_order,
                ):
                    T_intrinsic = ps.transformations.generate_rot_T(
                        angles, passive, True, intrinsic_order
                    )
                    T_extrinsic = ps.transformations.generate_rot_T(
                        angles, passive, False, extrinsic_order
                    )

                    # Extract rotation parts
                    R_intrinsic = T_intrinsic[:3, :3]
                    R_extrinsic = T_extrinsic[:3, :3]

                    npt.assert_allclose(R_intrinsic, R_extrinsic, atol=1e-14)

    def test_rotation_matrix_properties(self):
        """Tests that rotation components satisfy rotation matrix properties.

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
                    for order in ["xyz", "zyx"]:  # Test representative orders
                        with self.subTest(
                            angles=angles,
                            passive=passive,
                            intrinsic=intrinsic,
                            order=order,
                        ):
                            T = ps.transformations.generate_rot_T(
                                angles, passive, intrinsic, order
                            )
                            R = T[:3, :3]  # Extract rotation part

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
        # Test 90-degree rotation about x-axis (order "xyz", only first angle)
        angles_x90 = np.array([90.0, 0.0, 0.0])
        R_act_x90_expected = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])

        T_act_x90 = ps.transformations.generate_rot_T(angles_x90, False, True, "xyz")
        R_act_x90 = T_act_x90[:3, :3]
        npt.assert_allclose(R_act_x90, R_act_x90_expected, atol=1e-14)

        # Test 90-degree rotation about y-axis (order "xyz", only second angle)
        angles_y90 = np.array([0.0, 90.0, 0.0])
        R_act_y90_expected = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])

        T_act_y90 = ps.transformations.generate_rot_T(angles_y90, False, True, "xyz")
        R_act_y90 = T_act_y90[:3, :3]
        npt.assert_allclose(R_act_y90, R_act_y90_expected, atol=1e-14)

        # Test 90-degree rotation about z-axis (order "xyz", only third angle)
        angles_z90 = np.array([0.0, 0.0, 90.0])
        R_act_z90_expected = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

        T_act_z90 = ps.transformations.generate_rot_T(angles_z90, False, True, "xyz")
        R_act_z90 = T_act_z90[:3, :3]
        npt.assert_allclose(R_act_z90, R_act_z90_expected, atol=1e-14)

    def test_specific_known_passive_rotations(self):
        """Tests specific passive rotations with analytically known results.

        :return: None
        """
        # Test 90-degree rotation about x-axis (order "xyz", only first angle)
        angles_x90 = np.array([90.0, 0.0, 0.0])
        v_A = np.array([0.0, 1.0, 0.0])
        v_B_expected = np.array([0.0, 0.0, -1.0])

        T_pas_x90 = ps.transformations.generate_rot_T(angles_x90, True, True, "xyz")
        R_pas_x90 = T_pas_x90[:3, :3]
        v_B = R_pas_x90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

        # Test 90-degree rotation about y-axis (order "xyz", only second angle)
        angles_y90 = np.array([0.0, 90.0, 0.0])
        v_A = np.array([0.0, 0.0, 1.0])
        v_B_expected = np.array([-1.0, 0.0, 0.0])

        T_pas_y90 = ps.transformations.generate_rot_T(angles_y90, True, True, "xyz")
        R_pas_y90 = T_pas_y90[:3, :3]
        v_B = R_pas_y90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

        # Test 90-degree rotation about z-axis (order "xyz", only third angle)
        angles_z90 = np.array([0.0, 0.0, 90.0])
        v_A = np.array([1.0, 0.0, 0.0])
        v_B_expected = np.array([0.0, -1.0, 0])

        T_pas_z90 = ps.transformations.generate_rot_T(angles_z90, True, True, "xyz")
        R_pas_z90 = T_pas_z90[:3, :3]
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
        T1 = ps.transformations.generate_rot_T(angles1, False, True, "xyz")
        T2 = ps.transformations.generate_rot_T(angles2, False, True, "xyz")
        R1 = T1[:3, :3]
        R2 = T2[:3, :3]

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
                for order in ["xyz", "zyx"]:
                    with self.subTest(
                        passive=passive, intrinsic=intrinsic, order=order
                    ):
                        T_large = ps.transformations.generate_rot_T(
                            large_angles, passive, intrinsic, order
                        )
                        T_equivalent = ps.transformations.generate_rot_T(
                            equivalent_angles, passive, intrinsic, order
                        )

                        # Extract rotation parts
                        R_large = T_large[:3, :3]
                        R_equivalent = T_equivalent[:3, :3]

                        # Should produce the same rotation matrix
                        npt.assert_allclose(R_large, R_equivalent, atol=1e-14)

    def test_different_orders(self):
        """Tests all valid Tait-Bryan rotation orders.

        Ensures all six valid orders produce valid transformation matrices.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])
        valid_orders = ["xyz", "xzy", "yxz", "yzx", "zxy", "zyx"]

        for order in valid_orders:
            with self.subTest(order=order):
                for passive in [True, False]:
                    for intrinsic in [True, False]:
                        T = ps.transformations.generate_rot_T(
                            angles, passive, intrinsic, order
                        )
                        R = T[:3, :3]  # Extract rotation part

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
                    T = ps.transformations.generate_rot_T(angles, passive, True, "xyz")
                    R = T[:3, :3]  # Extract rotation part

                    # Should be a valid rotation matrix
                    det = np.linalg.det(R)
                    self.assertAlmostEqual(det, 1.0, places=14)

                    identity_test = R.T @ R
                    npt.assert_allclose(identity_test, np.eye(3), atol=1e-14)

    def test_homogeneous_coordinate_transformations(self):
        """Tests transformation of homogeneous coordinates.

        :return: None
        """
        # Test rotation of position vectors
        T_rot_act_z90 = ps.transformations.generate_rot_T(
            np.array([0.0, 0.0, 90.0]), False, True, "xyz"
        )

        # Test position vector [1, 0, 0] -> should become [0, 1, 0]
        pos_vec = np.array([1.0, 0.0, 0.0])
        transformed_pos = ps.transformations.apply_T_to_vectors(
            T_rot_act_z90, pos_vec, True
        )

        expected_transformed = np.array([0.0, 1.0, 0.0])
        npt.assert_allclose(transformed_pos, expected_transformed, atol=1e-14)

        # Test direction vector [1, 0, 0] -> should become [0, 1, 0]
        dir_vec = np.array([1.0, 0.0, 0.0])
        transformed_dir = ps.transformations.apply_T_to_vectors(
            T_rot_act_z90, dir_vec, False
        )

        expected_transformed_dir = np.array([0.0, 1.0, 0.0])
        npt.assert_allclose(transformed_dir, expected_transformed_dir, atol=1e-14)

    def test_invalid_rotation_order_rejected(self):
        """Tests that passing in invalid angle orders raises a value error.

        :return: None
        """
        angles = np.array([10.0, 20.0, 30.0])
        for bad in ["xyx", "xxx", "zz", "wxy", "x_y", ""]:
            with self.assertRaises(ValueError):
                ps.transformations.generate_rot_T(angles, True, True, bad)


class TestGenerateTransT(unittest.TestCase):
    """This class contains methods for testing the generate_trans_T function.

    This class contains the following public methods:
        test_passive_transformation_with_homog_coords: Tests passive translation
        transformations on generated homogeneous coordinates.

        test_active_transformation_with_homog_coords: Tests active translation
        transformations on homogeneous coordinates.

        test_zero_translation: Tests with zero translation vector.
        test_transformation_properties: Tests properties of transformation matrices.

        test_direction_vector_transformations: Tests that direction vectors are
        unaffected by translation transformations.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_passive_transformations_with_homog_coords(self):
        """Tests passive translation transformations on generated homogenous
        coordinates.

        :return: None
        """
        # Test passive transformation with position vector
        b_A_a = np.array([1.0, 2.0, 3.0])
        c_A_a = np.array([5.0, 6.0, 7.0])

        T_trans_pas = ps.transformations.generate_trans_T(b_A_a, True)
        c_A_b = ps.transformations.apply_T_to_vectors(T_trans_pas, c_A_a, True)

        # For passive transformation: new_position = old_position - translation
        expected_c_A_b = c_A_a - b_A_a
        npt.assert_array_equal(c_A_b, expected_c_A_b)

    def test_active_transformation_homog_coords(self):
        """Tests active translation transformations on homogeneous coordinates.

        :return: None
        """
        # Test active transformations
        cPrime_A_c = np.array([1.0, 2.0, 3.0])
        c_A_a = np.array([5.0, 6.0, 7.0])

        # Test active transformation with position vector
        T_trans_act = ps.transformations.generate_trans_T(cPrime_A_c, False)
        cPrime_A_a = ps.transformations.apply_T_to_vectors(T_trans_act, c_A_a, True)

        # For active transformation: new_position = old_position + translation
        expected_cPrime_A_a = c_A_a + cPrime_A_c
        npt.assert_array_equal(cPrime_A_a, expected_cPrime_A_a)

    def test_zero_translation(self):
        """Tests with zero translation vector.

        :return: None
        """
        zero_translation = np.array([0.0, 0.0, 0.0])

        # Test passive with zero translation
        T_passive_zero = ps.transformations.generate_trans_T(zero_translation, True)
        npt.assert_array_equal(T_passive_zero, np.eye(4))

        # Test active with zero translation
        T_active_zero = ps.transformations.generate_trans_T(zero_translation, False)
        npt.assert_array_equal(T_active_zero, np.eye(4))

    def test_transformation_properties(self):
        """Tests properties of transformation matrices.

        :return: None
        """
        translations = np.array([1.0, 2.0, 3.0])

        for passive in [True, False]:
            with self.subTest(passive=passive):
                T = ps.transformations.generate_trans_T(translations, passive)

                # Test output shape
                self.assertEqual(T.shape, (4, 4))

                # Test that rotation part is identity
                npt.assert_array_equal(T[:3, :3], np.eye(3))

                # Test that bottom row is [0, 0, 0, 1]
                expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
                npt.assert_array_equal(T[3, :], expected_bottom)

                # Test determinant is 1
                det = np.linalg.det(T)
                self.assertAlmostEqual(det, 1.0, places=14)

    def test_direction_vector_transformations(self):
        """Tests that direction vectors are unaffected by translation transformations.

        :return: None
        """

        # Generate passive and active transformation matrices.
        translation = np.array([1.0, 2.0, 3.0])
        T_trans_pas = ps.transformations.generate_trans_T(translation, True)
        T_trans_act = ps.transformations.generate_trans_T(translation, False)

        # Test that direction vectors are unaffected by translation
        direction = np.array([1.0, 0.0, 0.0])

        passive_transformed_direction = ps.transformations.apply_T_to_vectors(
            T_trans_pas, direction, False
        )
        active_transformed_direction = ps.transformations.apply_T_to_vectors(
            T_trans_act, direction, False
        )

        # Direction vectors should be unchanged by pure translation
        npt.assert_array_equal(passive_transformed_direction, direction)
        npt.assert_array_equal(active_transformed_direction, direction)


class TestGenerateReflectT(unittest.TestCase):
    """This class contains methods for testing the generate_reflect_T function.

    This class contains the following public methods:
        test_reflection_about_origin_planes: Tests reflection about planes through
        origin.
        test_reflection_about_offset_planes: Tests reflection about offset planes.
        test_zero_length_normal_error: Tests error handling for zero-length normals.
        test_transformation_properties: Tests properties of transformation matrices.
        test_symmetric_reflections: Tests that double reflections return original.
        test_normal_vector_normalization: Tests that non-unit normal vectors are
        correctly normalized.
        test_active_passive_equivalence: Tests that passive and active modes produce
        identical matrices.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_reflection_about_origin_planes(self):
        """Tests reflection about planes passing through the origin.

        :return: None
        """
        # Test reflection about xy-plane (z=0, normal=[0,0,1])
        plane_point = np.array([0.0, 0.0, 0.0])
        plane_normal = np.array([0.0, 0.0, 1.0])

        T_reflect_act = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, False
        )

        # Point [1, 2, 3] should reflect to [1, 2, -3]
        test_point = np.array([1.0, 2.0, 3.0])
        reflected_point = ps.transformations.apply_T_to_vectors(
            T_reflect_act, test_point, True
        )

        expected_reflected = np.array([1.0, 2.0, -3.0])
        npt.assert_allclose(reflected_point, expected_reflected, atol=1e-14)

        # Test reflection about xz-plane (y=0, normal=[0,1,0])
        plane_normal_y = np.array([0.0, 1.0, 0.0])
        T_reflect_act_y = ps.transformations.generate_reflect_T(
            plane_point, plane_normal_y, False
        )

        reflected_point_y = ps.transformations.apply_T_to_vectors(
            T_reflect_act_y, test_point, True
        )

        expected_reflected_y = np.array([1.0, -2.0, 3.0])
        npt.assert_allclose(reflected_point_y, expected_reflected_y, atol=1e-14)

    def test_reflection_about_offset_planes(self):
        """Tests reflection about planes not passing through the origin.

        :return: None
        """
        # Test reflection about plane z=2 (normal=[0,0,1], point=[0,0,2])
        plane_point = np.array([0.0, 0.0, 2.0])
        plane_normal = np.array([0.0, 0.0, 1.0])

        T_reflect_act = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, False
        )

        # Point [1, 2, 5] should reflect to [1, 2, -1]
        test_point = np.array([1.0, 2.0, 5.0])
        reflected_point = ps.transformations.apply_T_to_vectors(
            T_reflect_act, test_point, True
        )

        expected_reflected = np.array([1.0, 2.0, -1.0])
        npt.assert_allclose(reflected_point, expected_reflected, atol=1e-14)

        # Point on the plane should remain unchanged
        plane_test_point = np.array([3.0, 4.0, 2.0])
        reflected_plane_point = ps.transformations.apply_T_to_vectors(
            T_reflect_act, plane_test_point, True
        )

        npt.assert_allclose(reflected_plane_point, plane_test_point, atol=1e-14)

    def test_zero_length_normal_error(self):
        """Tests error handling for zero-length normal vectors.

        :return: None
        """
        plane_point = np.array([0.0, 0.0, 0.0])
        zero_normal = np.array([0.0, 0.0, 0.0])

        with self.assertRaises(ValueError) as context:
            ps.transformations.generate_reflect_T(plane_point, zero_normal, False)

        self.assertIn(
            "plane_normal_A must have a non-zero length.", str(context.exception)
        )

    def test_transformation_properties(self):
        """Tests properties of reflection transformation matrices.

        :return: None
        """
        plane_point = np.array([1.0, 1.0, 1.0])
        plane_normal = np.array([1.0, 1.0, 1.0])  # Will be normalized internally

        for passive in [True, False]:
            with self.subTest(passive=passive):
                T = ps.transformations.generate_reflect_T(
                    plane_point, plane_normal, passive
                )

                # Test output shape
                self.assertEqual(T.shape, (4, 4))

                # Test that bottom row is [0, 0, 0, 1]
                expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
                npt.assert_array_equal(T[3, :], expected_bottom)

                # Test determinant is -1 (proper reflection)
                det = np.linalg.det(T)
                self.assertAlmostEqual(det, -1.0, places=12)

                # Reflection part should have determinant -1
                reflection_det = np.linalg.det(T[:3, :3])
                self.assertAlmostEqual(reflection_det, -1.0, places=12)

    def test_symmetric_reflections(self):
        """Tests that double reflections return the original point.

        :return: None
        """
        plane_point = np.array([1.0, 2.0, 3.0])
        plane_normal = np.array([0.0, 1.0, 0.0])

        T_reflect_act = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, False
        )

        test_point = np.array([5.0, 8.0, -2.0])

        # Apply reflection twice
        reflected_once = ps.transformations.apply_T_to_vectors(
            T_reflect_act, test_point, True
        )
        final_point = ps.transformations.apply_T_to_vectors(
            T_reflect_act, reflected_once, True
        )

        npt.assert_allclose(final_point, test_point, atol=1e-14)

    def test_normal_vector_normalization(self):
        """Tests that non-unit normal vectors are correctly normalized.

        :return: None
        """
        plane_point = np.array([0.0, 0.0, 0.0])
        # Use non-unit normal vector - should be normalized internally
        plane_normal = np.array([0.0, 0.0, 5.0])  # Should normalize to [0,0,1]

        T_reflect_act = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, False
        )

        test_point = np.array([1.0, 2.0, 3.0])
        reflected_point = ps.transformations.apply_T_to_vectors(
            T_reflect_act, test_point, True
        )

        # Should give same result as unit normal [0,0,1]
        expected_reflected = np.array([1.0, 2.0, -3.0])
        npt.assert_allclose(reflected_point, expected_reflected, atol=1e-14)

    def test_passive_active_equivalence(self):
        """Tests that passive and active modes produce identical matrices.

        :return: None
        """
        plane_point = np.array([1.0, 2.0, 3.0])
        plane_normal = np.array([1.0, 0.0, 1.0])

        T_reflect_pas = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, True
        )
        T_reflect_act = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, False
        )

        # Matrices should be identical for reflections
        npt.assert_array_equal(T_reflect_pas, T_reflect_act)


class TestComposeTPas(unittest.TestCase):
    """This class contains methods for testing the compose_T_pas function.

    This class contains the following public methods:
        test_basic_passive_composition: Tests basic passive transformation composition.
        test_composition_order: Tests that composition order matters.
        test_identity_composition: Tests composition with identity matrices.
        test_inverse_composition: Tests composition with inverse transformations.
        test_multiple_transformations: Tests composition of multiple transformations.
        test_rotation_translation_composition: Tests composition of rotation and translation.
        test_matrix_properties: Tests properties of composed matrices.
        test_empty_chain_raises: Tests that passing an empty transformation chain
        results in a value error.
        test_specific_known_passive_composition: Tests specific composition of
        passive transformations with a known result.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_basic_passive_composition(self):
        """Tests basic passive transformation composition.

        :return: None
        """
        # Create translation and rotation transformations
        translation = np.array([1.0, 2.0, 3.0])
        angles = np.array([30.0, 0.0, 0.0])

        T1 = ps.transformations.generate_trans_T(translation, True)
        T2 = ps.transformations.generate_rot_T(angles, True, True, "xyz")

        # Test composition
        T_composed = ps.transformations.compose_T_pas(T1, T2)

        # Test that result has correct structure
        self.assertEqual(T_composed.shape, (4, 4))

        # Test that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_composed[3, :], expected_bottom)

    def test_composition_order(self):
        """Tests that composition order matters for non-commuting transformations.

        :return: None
        """
        translation = np.array([1.0, 0.0, 0.0])
        angles = np.array([0.0, 0.0, 90.0])

        T_trans = ps.transformations.generate_trans_T(translation, True)
        T_rot = ps.transformations.generate_rot_T(angles, True, True, "xyz")

        # Test different composition orders
        T1 = ps.transformations.compose_T_pas(
            T_trans, T_rot
        )  # Translation then rotation
        T2 = ps.transformations.compose_T_pas(
            T_rot, T_trans
        )  # Rotation then translation

        # Should produce different results for non-commuting transformations
        self.assertFalse(np.allclose(T1, T2))

    def test_identity_composition(self):
        """Tests composition with identity matrices.

        :return: None
        """
        identity = np.eye(4)
        translation = np.array([1.0, 2.0, 3.0])
        T_trans = ps.transformations.generate_trans_T(translation, True)

        # Composition with identity should leave transformation unchanged
        T_composed1 = ps.transformations.compose_T_pas(identity, T_trans)
        T_composed2 = ps.transformations.compose_T_pas(T_trans, identity)

        npt.assert_allclose(T_composed1, T_trans, atol=1e-14)
        npt.assert_allclose(T_composed2, T_trans, atol=1e-14)

    def test_inverse_composition(self):
        """Tests composition with inverse transformations.

        :return: None
        """
        translation = np.array([1.0, 2.0, 3.0])
        T_trans = ps.transformations.generate_trans_T(translation, True)
        T_trans_inv = ps.transformations.invert_T_pas(T_trans)

        # Composition with inverse should yield identity
        T_composed1 = ps.transformations.compose_T_pas(T_trans, T_trans_inv)
        T_composed2 = ps.transformations.compose_T_pas(T_trans_inv, T_trans)

        npt.assert_allclose(T_composed1, np.eye(4), atol=1e-14)
        npt.assert_allclose(T_composed2, np.eye(4), atol=1e-14)

    def test_multiple_transformations(self):
        """Tests composition of multiple transformations.

        :return: None
        """
        # Create three different transformations
        T1 = ps.transformations.generate_trans_T(np.array([1.0, 0.0, 0.0]), True)
        T2 = ps.transformations.generate_rot_T(
            np.array([0.0, 0.0, 90.0]), True, True, "xyz"
        )
        T3 = ps.transformations.generate_trans_T(np.array([0.0, 1.0, 0.0]), True)

        # Test composition of multiple matrices
        T_composed = ps.transformations.compose_T_pas(T1, T2, T3)

        # Test that result has correct structure
        self.assertEqual(T_composed.shape, (4, 4))

        # Test determinant (should be 1 for proper transformations)
        det = np.linalg.det(T_composed)
        self.assertAlmostEqual(det, 1.0, places=12)

    def test_rotation_translation_composition(self):
        """Tests composition of rotation and translation transformations.

        :return: None
        """
        # Test vector
        test_vector = np.array([1.0, 0.0, 0.0])

        # Individual transformations
        T_rot = ps.transformations.generate_rot_T(
            np.array([0.0, 0.0, 90.0]), True, True, "xyz"
        )
        T_trans = ps.transformations.generate_trans_T(np.array([1.0, 1.0, 0.0]), True)

        # Composed transformation
        T_composed = ps.transformations.compose_T_pas(T_rot, T_trans)

        # Apply transformations sequentially
        v1 = ps.transformations.apply_T_to_vectors(T_rot, test_vector, True)
        v2 = ps.transformations.apply_T_to_vectors(T_trans, v1, True)

        # Apply composed transformation
        v_composed = ps.transformations.apply_T_to_vectors(
            T_composed, test_vector, True
        )

        # Results should be identical
        npt.assert_allclose(v2, v_composed, atol=1e-14)

    def test_matrix_properties(self):
        """Tests properties of composed transformation matrices.

        :return: None
        """
        # Create various transformations
        T1 = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), True, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), True)
        T3 = ps.transformations.generate_rot_T(
            np.array([0.0, 45.0, 0.0]), True, True, "xyz"
        )

        T_composed = ps.transformations.compose_T_pas(T1, T2, T3)

        # Test output shape
        self.assertEqual(T_composed.shape, (4, 4))

        # Test that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_composed[3, :], expected_bottom)

        # Test determinant is 1
        det = np.linalg.det(T_composed)
        self.assertAlmostEqual(det, 1.0, places=12)

    def test_empty_chain_raises(self):
        """Tests that passing an empty transformation chain results in a value error.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.transformations.compose_T_pas()

    def test_specific_known_passive_composition(self):
        """Tests specific composition of passive transformations with a known result.

        :return: None
        """
        # Goal: c_Wn_Ler, "the position of point c (in wing axes, relative to the
        # leading edge root point)"

        # Given:
        # The position of point c (in geometry axes, relative to the CG point)
        c_G_Cg = [0.5, -1.0, 2.0]

        # Given:
        # The position of the leading edge root point (in geometry axes, relative to
        # the CG point)
        Ler_G_CG = [1.0, 2.0, 0.5]

        # Given:
        # The orientation of wing axes relative to geometry axes using an intrinsic
        # z-y'-x" rotation
        angles_G_to_Wn_izyx = [0.0, 0.0, 90.0]

        T_rot_pas_G_to_Wn = ps.transformations.generate_rot_T(
            angles_G_to_Wn_izyx, passive=True, intrinsic=True, order="zyx"
        )
        T_trans_pas_G_Cg_to_G_Ler = ps.transformations.generate_trans_T(
            Ler_G_CG, passive=True
        )

        T_pas_G_Cg_to_Wn_Ler = ps.transformations.compose_T_pas(
            T_trans_pas_G_Cg_to_G_Ler, T_rot_pas_G_to_Wn
        )

        c_Wn_Ler = ps.transformations.apply_T_to_vectors(
            T_pas_G_Cg_to_Wn_Ler, c_G_Cg, has_point=True
        )

        # Expected value calculated using CAD model
        c_Wn_Ler_expected = np.array([-3.0, 0.5, 1.5])

        npt.assert_allclose(c_Wn_Ler, c_Wn_Ler_expected)


class TestComposeTAct(unittest.TestCase):
    """This class contains methods for testing the compose_T_act function.

    This class contains the following public methods:
        test_basic_active_composition: Tests basic active transformation composition.
        test_composition_vs_manual: Tests composition against manual matrix multiplication.
        test_identity_composition: Tests composition with identity matrices.
        test_inverse_composition: Tests composition with inverse transformations.
        test_multiple_transformations: Tests composition of multiple transformations.
        test_matrix_properties: Tests properties of composed matrices.
        test_empty_chain_raises: Tests that passing an empty transformation chain
        results in a value error.
        test_specific_known_active_composition: Tests specific composition of active
        transformations with a known result.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_basic_active_composition(self):
        """Tests basic active transformation composition.

        :return: None
        """
        # Create translation and rotation transformations
        translation = np.array([1.0, 2.0, 3.0])
        angles = np.array([0.0, 0.0, 45.0])

        T1 = ps.transformations.generate_trans_T(translation, False)
        T2 = ps.transformations.generate_rot_T(angles, False, True, "xyz")

        # Test composition
        T_composed = ps.transformations.compose_T_act(T1, T2)

        # Test that result has correct structure
        self.assertEqual(T_composed.shape, (4, 4))

        # Test that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_composed[3, :], expected_bottom)

    def test_composition_vs_manual(self):
        """Tests composition against manual matrix multiplication.

        :return: None
        """
        # Create transformations
        T1 = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), False, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), False)

        # Test composed matrix
        T_composed = ps.transformations.compose_T_act(T1, T2)

        # Manual composition for active transformations
        T_manual = T2 @ T1

        npt.assert_allclose(T_composed, T_manual, atol=1e-14)

    def test_identity_composition(self):
        """Tests composition with identity matrices.

        :return: None
        """
        identity = np.eye(4)
        rotation = np.array([30.0, 45.0, 60.0])
        T_rot = ps.transformations.generate_rot_T(rotation, False, True, "xyz")

        # Composition with identity should leave transformation unchanged
        T_composed1 = ps.transformations.compose_T_act(identity, T_rot)
        T_composed2 = ps.transformations.compose_T_act(T_rot, identity)

        npt.assert_allclose(T_composed1, T_rot, atol=1e-14)
        npt.assert_allclose(T_composed2, T_rot, atol=1e-14)

    def test_inverse_composition(self):
        """Tests composition with inverse transformations.

        :return: None
        """
        rotation = np.array([30.0, 45.0, 60.0])
        T_rot = ps.transformations.generate_rot_T(rotation, False, True, "xyz")
        T_rot_inv = ps.transformations.invert_T_act(T_rot)

        # Composition with inverse should yield identity
        T_composed1 = ps.transformations.compose_T_act(T_rot, T_rot_inv)
        T_composed2 = ps.transformations.compose_T_act(T_rot_inv, T_rot)

        npt.assert_allclose(T_composed1, np.eye(4), atol=1e-14)
        npt.assert_allclose(T_composed2, np.eye(4), atol=1e-14)

    def test_multiple_transformations(self):
        """Tests composition of multiple transformations.

        :return: None
        """
        # Create four different transformations
        T1 = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), False, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([1.0, 0.0, 0.0]), False)
        T3 = ps.transformations.generate_rot_T(
            np.array([0.0, 45.0, 0.0]), False, True, "xyz"
        )
        T4 = ps.transformations.generate_trans_T(np.array([0.0, 1.0, 0.0]), False)

        # Test composition of multiple matrices
        T_composed = ps.transformations.compose_T_act(T1, T2, T3, T4)

        # Test that result has correct structure
        self.assertEqual(T_composed.shape, (4, 4))

        # Test determinant (should be 1 for proper transformations)
        det = np.linalg.det(T_composed)
        self.assertAlmostEqual(det, 1.0, places=12)

    def test_matrix_properties(self):
        """Tests properties of composed transformation matrices.

        :return: None
        """
        # Create various transformations
        T1 = ps.transformations.generate_reflect_T(
            np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), False
        )
        T2 = ps.transformations.generate_rot_T(
            np.array([0.0, 0.0, 90.0]), False, True, "xyz"
        )
        T3 = ps.transformations.generate_trans_T(np.array([2.0, 3.0, 1.0]), False)

        T_composed = ps.transformations.compose_T_act(T1, T2, T3)

        # Test output shape
        self.assertEqual(T_composed.shape, (4, 4))

        # Test that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_composed[3, :], expected_bottom)

        # Test determinant is -1 (includes reflection)
        det = np.linalg.det(T_composed)
        self.assertAlmostEqual(det, -1.0, places=12)

    def test_empty_chain_raises(self):
        """Tests that passing an empty transformation chain results in a value error.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.transformations.compose_T_act()

    def test_specific_known_active_composition(self):
        """Tests specific composition of active transformations with a known result.

        :return: None
        """
        c_G = [1.0, 2.0, 3.0]

        angles_act_izyx = [0.0, 0.0, 90.0]
        t_G = [10.0, 0.0, 0.0]

        rot_T_act = ps.transformations.generate_rot_T(
            angles_act_izyx, passive=False, intrinsic=True, order="zyx"
        )
        trans_T_act = ps.transformations.generate_trans_T(t_G, passive=False)

        T_act = ps.transformations.compose_T_act(rot_T_act, trans_T_act)

        cPrime_G = ps.transformations.apply_T_to_vectors(T_act, c_G, has_point=True)

        # Expected known value. If this expected value is confusing to you, and you
        # expect it to instead by np.array([-2.0, 11.0, 3.0]), see the note in the
        # docstrings of compose_T_act on world-fixed vs. body-fixed transformations.
        cPrime_G_expected = np.array([8.0, 1.0, 3.0])

        npt.assert_allclose(cPrime_G, cPrime_G_expected)


class TestInvertTPas(unittest.TestCase):
    """This class contains methods for testing the invert_T_pas function.

    This class contains the following public methods:
        test_translation_inversion: Tests inversion of pure translation matrices.
        test_rotation_inversion: Tests inversion of pure rotation matrices.
        test_reflection_inversion: Tests inversion of pure reflection matrices.
        test_combined_transformation_inversion: Tests inversion of combined transformations.
        test_identity_inversion: Tests inversion of identity matrix.
        test_double_inversion: Tests that double inversion returns original.
        test_inversion_properties: Tests mathematical properties of inverted matrices.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_translation_inversion(self):
        """Tests inversion of pure translation matrices.

        :return: None
        """
        translation = np.array([1.0, 2.0, 3.0])
        T_trans = ps.transformations.generate_trans_T(translation, True)
        T_trans_inv = ps.transformations.invert_T_pas(T_trans)

        # Test that T @ T_inv = I
        result = T_trans @ T_trans_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

        # Test that T_inv @ T = I
        result2 = T_trans_inv @ T_trans
        npt.assert_allclose(result2, np.eye(4), atol=1e-14)

    def test_rotation_inversion(self):
        """Tests inversion of pure rotation matrices.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])
        T_rot = ps.transformations.generate_rot_T(angles, True, True, "xyz")
        T_rot_inv = ps.transformations.invert_T_pas(T_rot)

        # Test that T @ T_inv = I
        result = T_rot @ T_rot_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

        # Test that T_inv @ T = I
        result2 = T_rot_inv @ T_rot
        npt.assert_allclose(result2, np.eye(4), atol=1e-14)

    def test_reflection_inversion(self):
        """Tests inversion of pure reflection matrices.

        :return: None
        """
        plane_point = np.array([1.0, 2.0, 3.0])
        plane_normal = np.array([0.0, 1.0, 0.0])
        T_reflect = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, True
        )
        T_reflect_inv = ps.transformations.invert_T_pas(T_reflect)

        # For reflections, the inverse should equal the original
        npt.assert_allclose(T_reflect_inv, T_reflect, atol=1e-14)

        # Test that T @ T_inv = I
        result = T_reflect @ T_reflect_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

    def test_combined_transformation_inversion(self):
        """Tests inversion of combined transformations.

        :return: None
        """
        # Create composed transformation
        T1 = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), True, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), True)
        T_combined = ps.transformations.compose_T_pas(T1, T2)

        # Invert composed transformation
        T_combined_inv = ps.transformations.invert_T_pas(T_combined)

        # Test that T @ T_inv = I
        result = T_combined @ T_combined_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

        # Test that T_inv @ T = I
        result2 = T_combined_inv @ T_combined
        npt.assert_allclose(result2, np.eye(4), atol=1e-14)

    def test_identity_inversion(self):
        """Tests inversion of identity matrix.

        :return: None
        """
        identity = np.eye(4)
        identity_inv = ps.transformations.invert_T_pas(identity)

        # Inverse of identity should be identity
        npt.assert_allclose(identity_inv, identity, atol=1e-14)

    def test_double_inversion(self):
        """Tests that double inversion returns the original matrix.

        :return: None
        """
        # Create various transformations
        transformations = [
            ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), True),
            ps.transformations.generate_rot_T(
                np.array([30.0, 45.0, 60.0]), True, True, "xyz"
            ),
            ps.transformations.generate_reflect_T(
                np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0]), True
            ),
        ]

        for i, T in enumerate(transformations):
            with self.subTest(transformation=i):
                T_inv = ps.transformations.invert_T_pas(T)
                T_inv_inv = ps.transformations.invert_T_pas(T_inv)

                # Double inversion should return original
                npt.assert_allclose(T_inv_inv, T, atol=1e-14)

    def test_inversion_properties(self):
        """Tests mathematical properties of inverted matrices.

        :return: None
        """
        # Create transformation
        T1 = ps.transformations.generate_rot_T(
            np.array([45.0, 30.0, 60.0]), True, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([2.0, 3.0, 1.0]), True)
        T_combined = ps.transformations.compose_T_pas(T1, T2)

        T_inv = ps.transformations.invert_T_pas(T_combined)

        # Test output shape
        self.assertEqual(T_inv.shape, (4, 4))

        # Test that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_inv[3, :], expected_bottom)

        # Test determinant (should be 1/det(original))
        det_original = np.linalg.det(T_combined)
        det_inv = np.linalg.det(T_inv)
        self.assertAlmostEqual(det_inv * det_original, 1.0, places=12)


class TestInvertTAct(unittest.TestCase):
    """This class contains methods for testing the invert_T_act function.

    This class contains the following public methods:
        test_translation_inversion: Tests inversion of pure translation matrices.
        test_rotation_inversion: Tests inversion of pure rotation matrices.
        test_reflection_inversion: Tests inversion of pure reflection matrices.
        test_combined_transformation_inversion: Tests inversion of combined transformations.
        test_identity_inversion: Tests inversion of identity matrix.
        test_double_inversion: Tests that double inversion returns original.
        test_inversion_properties: Tests mathematical properties of inverted matrices.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_translation_inversion(self):
        """Tests inversion of pure translation matrices.

        :return: None
        """
        translation = np.array([1.0, 2.0, 3.0])
        T_trans = ps.transformations.generate_trans_T(translation, False)
        T_trans_inv = ps.transformations.invert_T_act(T_trans)

        # Test that T @ T_inv = I
        result = T_trans @ T_trans_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

        # Test that T_inv @ T = I
        result2 = T_trans_inv @ T_trans
        npt.assert_allclose(result2, np.eye(4), atol=1e-14)

    def test_rotation_inversion(self):
        """Tests inversion of pure rotation matrices.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])
        T_rot = ps.transformations.generate_rot_T(angles, False, True, "xyz")
        T_rot_inv = ps.transformations.invert_T_act(T_rot)

        # Test that T @ T_inv = I
        result = T_rot @ T_rot_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

        # Test that T_inv @ T = I
        result2 = T_rot_inv @ T_rot
        npt.assert_allclose(result2, np.eye(4), atol=1e-14)

    def test_reflection_inversion(self):
        """Tests inversion of pure reflection matrices.

        :return: None
        """
        plane_point = np.array([1.0, 2.0, 3.0])
        plane_normal = np.array([0.0, 1.0, 0.0])
        T_reflect = ps.transformations.generate_reflect_T(
            plane_point, plane_normal, False
        )
        T_reflect_inv = ps.transformations.invert_T_act(T_reflect)

        # For reflections, the inverse should equal the original
        npt.assert_allclose(T_reflect_inv, T_reflect, atol=1e-14)

        # Test that T @ T_inv = I
        result = T_reflect @ T_reflect_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

    def test_combined_transformation_inversion(self):
        """Tests inversion of combined transformations.

        :return: None
        """
        # Create composed transformation
        T1 = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), False, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), False)
        T_combined = ps.transformations.compose_T_act(T1, T2)

        # Invert composed transformation
        T_combined_inv = ps.transformations.invert_T_act(T_combined)

        # Test that T @ T_inv = I
        result = T_combined @ T_combined_inv
        npt.assert_allclose(result, np.eye(4), atol=1e-14)

        # Test that T_inv @ T = I
        result2 = T_combined_inv @ T_combined
        npt.assert_allclose(result2, np.eye(4), atol=1e-14)

    def test_identity_inversion(self):
        """Tests inversion of identity matrix.

        :return: None
        """
        identity = np.eye(4)
        identity_inv = ps.transformations.invert_T_act(identity)

        # Inverse of identity should be identity
        npt.assert_allclose(identity_inv, identity, atol=1e-14)

    def test_double_inversion(self):
        """Tests that double inversion returns the original matrix.

        :return: None
        """
        # Create various transformations
        transformations = [
            ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), False),
            ps.transformations.generate_rot_T(
                np.array([30.0, 45.0, 60.0]), False, True, "xyz"
            ),
            ps.transformations.generate_reflect_T(
                np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0]), False
            ),
        ]

        for i, T in enumerate(transformations):
            with self.subTest(transformation=i):
                T_inv = ps.transformations.invert_T_act(T)
                T_inv_inv = ps.transformations.invert_T_act(T_inv)

                # Double inversion should return original
                npt.assert_allclose(T_inv_inv, T, atol=1e-14)

    def test_inversion_properties(self):
        """Tests mathematical properties of inverted matrices.

        :return: None
        """
        # Create transformation
        T1 = ps.transformations.generate_rot_T(
            np.array([45.0, 30.0, 60.0]), False, True, "xyz"
        )
        T2 = ps.transformations.generate_trans_T(np.array([2.0, 3.0, 1.0]), False)
        T_combined = ps.transformations.compose_T_act(T1, T2)

        T_inv = ps.transformations.invert_T_act(T_combined)

        # Test output shape
        self.assertEqual(T_inv.shape, (4, 4))

        # Test that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_inv[3, :], expected_bottom)

        # Test determinant (should be 1/det(original))
        det_original = np.linalg.det(T_combined)
        det_inv = np.linalg.det(T_inv)
        self.assertAlmostEqual(det_inv * det_original, 1.0, places=12)


class TestApplyTToVectors(unittest.TestCase):
    """This class contains methods for testing the apply_T_to_vectors function.

    This class contains the following public methods:
        test_position_vector_transformation: Tests transformation of position vectors.

        test_direction_vector_transformation: Tests transformation of direction vectors.

        test_input_validation: Tests various input types and edge cases.

        test_transformation_consistency: Tests consistency with manual homogeneous
        operations.

        test_various_transformation_types: Tests with different types of
        transformations.

        test_vector_types: Tests with different vector input types.

        test_single_vector_compatibility: Tests that single vector behavior is as
        expected.

        test_array_of_vectors_transformation: Tests transformation of arrays of vectors.

        test_higher_dimensional_arrays: Tests transformation of higher dimensional
        arrays.

        test_position_vs_direction_arrays: Tests that has_point parameter works
        correctly with arrays.

        test_array_input_validation: Tests validation of array inputs.

        test_array_consistency_with_single_applications: Tests that array
        transformation gives same results as individual applications.

        test_edge_cases_arrays: Tests edge cases with array inputs.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_position_vector_transformation(self):
        """Tests transformation of position vectors (has_point=True).

        :return: None
        """
        # Test translation transformation
        translation = np.array([1.0, 2.0, 3.0])
        T_trans = ps.transformations.generate_trans_T(translation, False)

        position = np.array([5.0, 6.0, 7.0])
        transformed_position = ps.transformations.apply_T_to_vectors(
            T_trans, position, True
        )

        expected = position + translation
        npt.assert_array_equal(transformed_position, expected)

    def test_direction_vector_transformation(self):
        """Tests transformation of direction vectors (has_point=False).

        :return: None
        """
        # Test rotation transformation
        angles = np.array([0.0, 0.0, 90.0])
        T_rot = ps.transformations.generate_rot_T(angles, False, True, "xyz")

        direction = np.array([1.0, 0.0, 0.0])
        transformed_direction = ps.transformations.apply_T_to_vectors(
            T_rot, direction, False
        )

        expected = np.array([0.0, 1.0, 0.0])
        npt.assert_allclose(transformed_direction, expected, atol=1e-14)

    def test_input_validation(self):
        """Tests various input types and edge cases.

        :return: None
        """
        T = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), False, True, "xyz"
        )

        # Test with tuple input
        vector_tuple = (1.0, 2.0, 3.0)
        result_tuple = ps.transformations.apply_T_to_vectors(T, vector_tuple, True)
        self.assertIsInstance(result_tuple, np.ndarray)
        self.assertEqual(len(result_tuple), 3)

        # Test with list input
        vector_list = [1.0, 2.0, 3.0]
        result_list = ps.transformations.apply_T_to_vectors(T, vector_list, True)
        self.assertIsInstance(result_list, np.ndarray)
        self.assertEqual(len(result_list), 3)

        # Test with numpy array input
        vector_array = np.array([1.0, 2.0, 3.0])
        result_array = ps.transformations.apply_T_to_vectors(T, vector_array, True)
        self.assertIsInstance(result_array, np.ndarray)
        self.assertEqual(len(result_array), 3)

    def test_transformation_consistency(self):
        """Tests consistency with manual homogeneous coordinate operations.

        :return: None
        """
        # Create transformation
        T = ps.transformations.generate_rot_T(
            np.array([45.0, 30.0, 60.0]), False, True, "xyz"
        )

        # Test vector
        vector = np.array([1.0, 2.0, 3.0])

        # Using apply_T_to_vectors for position vector
        result_position = ps.transformations.apply_T_to_vectors(T, vector, True)

        # Manual homogeneous approach for position vector
        homog_vector = np.append(vector, 1.0)
        manual_result_homog = T @ homog_vector
        manual_result_position = manual_result_homog[:3]

        npt.assert_allclose(result_position, manual_result_position, atol=1e-14)

        # Using apply_T_to_vectors for direction vector
        result_direction = ps.transformations.apply_T_to_vectors(T, vector, False)

        # Manual homogeneous approach for direction vector
        homog_vector_dir = np.append(vector, 0.0)
        manual_result_homog_dir = T @ homog_vector_dir
        manual_result_direction = manual_result_homog_dir[:3]

        npt.assert_allclose(result_direction, manual_result_direction, atol=1e-14)

    def test_various_transformation_types(self):
        """Tests with different types of transformations.

        :return: None
        """
        test_vector = np.array([1.0, 2.0, 3.0])

        # Translation transformation
        T_trans = ps.transformations.generate_trans_T(np.array([1.0, 1.0, 1.0]), False)
        result_trans = ps.transformations.apply_T_to_vectors(T_trans, test_vector, True)
        self.assertEqual(len(result_trans), 3)

        # Rotation transformation
        T_rot = ps.transformations.generate_rot_T(
            np.array([30.0, 45.0, 60.0]), False, True, "xyz"
        )
        result_rot = ps.transformations.apply_T_to_vectors(T_rot, test_vector, True)
        self.assertEqual(len(result_rot), 3)

        # Reflection transformation
        T_reflect = ps.transformations.generate_reflect_T(
            np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), False
        )
        result_reflect = ps.transformations.apply_T_to_vectors(
            T_reflect, test_vector, True
        )
        self.assertEqual(len(result_reflect), 3)

        # Composed transformation
        T_composed = ps.transformations.compose_T_act(T_trans, T_rot)
        result_composed = ps.transformations.apply_T_to_vectors(
            T_composed, test_vector, True
        )
        self.assertEqual(len(result_composed), 3)

    def test_vector_types(self):
        """Tests with different vector input types and sizes.

        :return: None
        """
        T = ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), False)

        # Test with integer vector
        int_vector = np.array([1, 2, 3])
        result_int = ps.transformations.apply_T_to_vectors(T, int_vector, True)
        self.assertEqual(result_int.dtype, np.float64)

        # Test with float vector
        float_vector = np.array([1.5, 2.5, 3.5])
        result_float = ps.transformations.apply_T_to_vectors(T, float_vector, True)
        self.assertEqual(result_float.dtype, np.float64)

        # Test with zero vector
        zero_vector = np.array([0.0, 0.0, 0.0])
        result_zero = ps.transformations.apply_T_to_vectors(T, zero_vector, True)
        expected_zero = np.array(
            [1.0, 2.0, 3.0]
        )  # Translation applied to zero position
        npt.assert_array_equal(result_zero, expected_zero)

    def test_single_vector_compatibility(self):
        """Tests that single vector behavior is as expected.

        :return: None
        """
        T = ps.transformations.generate_rot_T(
            np.array([45.0, 0.0, 0.0]), False, True, "xyz"
        )

        # Test with single vector as (3,) shape
        single_vector = np.array([1.0, 2.0, 3.0])
        result_single = ps.transformations.apply_T_to_vectors(T, single_vector, True)

        # Should return same shape as input
        self.assertEqual(result_single.shape, (3,))
        self.assertIsInstance(result_single, np.ndarray)

        # Test consistency between has_point=True and has_point=False
        result_position = ps.transformations.apply_T_to_vectors(T, single_vector, True)
        result_direction = ps.transformations.apply_T_to_vectors(
            T, single_vector, False
        )
        self.assertEqual(len(result_position), 3)
        self.assertEqual(len(result_direction), 3)

    def test_array_of_vectors_transformation(self):
        """Tests transformation of arrays of vectors.

        :return: None
        """
        # Create translation transformation
        T = ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), False)

        # Test with 2D array (multiple vectors)
        vectors_2d = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        result_2d = ps.transformations.apply_T_to_vectors(T, vectors_2d, True)

        expected_2d = np.array([[2.0, 2.0, 3.0], [1.0, 3.0, 3.0], [1.0, 2.0, 4.0]])
        npt.assert_array_equal(result_2d, expected_2d)
        self.assertEqual(result_2d.shape, (3, 3))

    def test_higher_dimensional_arrays(self):
        """Tests transformation of higher dimensional arrays.

        :return: None
        """
        T = ps.transformations.generate_rot_T(
            np.array([0.0, 0.0, 90.0]), False, True, "xyz"
        )

        # Test with 3D array (2x2x3 array of vectors)
        vectors_3d = np.array(
            [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 1.0], [1.0, 1.0, 0.0]]]
        )
        result_3d = ps.transformations.apply_T_to_vectors(T, vectors_3d, False)

        # Z-rotation by 90 degrees: [1,0,0] -> [0,1,0], [0,1,0] -> [-1,0,0]
        expected_3d = np.array(
            [[[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], [[0.0, 0.0, 1.0], [-1.0, 1.0, 0.0]]]
        )
        npt.assert_allclose(result_3d, expected_3d, atol=1e-14)
        self.assertEqual(result_3d.shape, (2, 2, 3))

    def test_position_vs_direction_arrays(self):
        """Tests that has_point parameter works correctly with arrays.

        :return: None
        """
        # Translation transformation
        translation = np.array([5.0, 10.0, 15.0])
        T = ps.transformations.generate_trans_T(translation, False)

        vectors = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

        # Position vectors should be translated
        result_positions = ps.transformations.apply_T_to_vectors(T, vectors, True)
        expected_positions = vectors + translation
        npt.assert_array_equal(result_positions, expected_positions)

        # Direction vectors should not be translated (pure translation matrix)
        result_directions = ps.transformations.apply_T_to_vectors(T, vectors, False)
        npt.assert_array_equal(result_directions, vectors)  # No change expected

    def test_array_input_validation(self):
        """Tests validation of array inputs.

        :return: None
        """
        T = ps.transformations.generate_rot_T(
            np.array([30.0, 0.0, 0.0]), False, True, "xyz"
        )

        # Test various valid array shapes
        valid_inputs = [
            np.array([1.0, 2.0, 3.0]),  # Single vector (3,)
            np.array([[1.0, 2.0, 3.0]]),  # Single vector as (1,3)
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),  # Multiple vectors (2,3)
            np.array([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]]),  # 3D array (1,2,3)
        ]

        for i, vectors in enumerate(valid_inputs):
            with self.subTest(i=i, shape=vectors.shape):
                result = ps.transformations.apply_T_to_vectors(T, vectors, True)
                self.assertEqual(result.shape, vectors.shape)
                self.assertEqual(result.shape[-1], 3)

    def test_array_consistency_with_single_applications(self):
        """Tests that array transformation gives same results as individual
        applications.

        :return: None
        """
        # Create a complex transformation
        T = ps.transformations.compose_T_act(
            ps.transformations.generate_rot_T(
                np.array([30.0, 45.0, 0.0]), False, True, "xyz"
            ),
            ps.transformations.generate_trans_T(np.array([1.0, 2.0, 3.0]), False),
        )

        # Test vectors
        test_vectors = np.array(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
        )

        # Transform as array
        array_result = ps.transformations.apply_T_to_vectors(T, test_vectors, True)

        # Transform individually and compare
        individual_results = []
        for vector in test_vectors:
            individual_result = ps.transformations.apply_T_to_vectors(T, vector, True)
            individual_results.append(individual_result)
        individual_results = np.array(individual_results)

        npt.assert_allclose(array_result, individual_results, atol=1e-14)

    def test_edge_cases_arrays(self):
        """Tests edge cases with array inputs.

        :return: None
        """
        T = ps.transformations.generate_rot_T(
            np.array([0.0, 0.0, 0.0]), False, True, "xyz"  # Identity rotation
        )

        # Test with empty-like inputs that should still work
        single_zero = np.array([0.0, 0.0, 0.0])
        result_single_zero = ps.transformations.apply_T_to_vectors(T, single_zero, True)
        npt.assert_array_equal(result_single_zero, single_zero)

        # Test with array of zeros
        zeros_array = np.zeros((3, 3))
        result_zeros_array = ps.transformations.apply_T_to_vectors(T, zeros_array, True)
        npt.assert_array_equal(result_zeros_array, zeros_array)


if __name__ == "__main__":
    unittest.main()
