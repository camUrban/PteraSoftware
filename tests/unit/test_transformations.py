"""This module contains classes to test functions in the transformations module.

This module contains the following classes:
    TestGenerateHomog: This class contains methods for testing the generate_homog
    function.
    TestGenerateR: This class contains methods for testing the generate_R function.
    TestGenerateTRot: This class contains methods for testing the generate_T_rot
    function.
    TestGenerateTTrans: This class contains methods for testing the generate_T_trans
    function.
    TestGenerateTReflect: This class contains methods for testing the
    generate_T_reflect function.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps


class TestGenerateHomog(unittest.TestCase):
    """This class contains methods for testing the generate_homog function.

    This class contains the following public methods:
        test_position_vector_conversion: Tests conversion of position vectors.
        test_direction_vector_conversion: Tests conversion of direction vectors.
        test_input_validation: Tests various input types and edge cases.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_position_vector_conversion(self):
        """Tests conversion of position vectors (has_point=True).

        :return: None
        """
        # Test basic position vector
        position = np.array([1.0, 2.0, 3.0])
        homog_position = ps.transformations.generate_homog(position, True)

        expected = np.array([1.0, 2.0, 3.0, 1.0])
        npt.assert_array_equal(homog_position, expected)

        # Test zero position vector
        zero_position = np.array([0.0, 0.0, 0.0])
        homog_zero = ps.transformations.generate_homog(zero_position, True)

        expected_zero = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(homog_zero, expected_zero)

        # Test negative position vector
        neg_position = np.array([-5.0, -10.0, -15.0])
        homog_neg = ps.transformations.generate_homog(neg_position, True)

        expected_neg = np.array([-5.0, -10.0, -15.0, 1.0])
        npt.assert_array_equal(homog_neg, expected_neg)

    def test_direction_vector_conversion(self):
        """Tests conversion of direction vectors (has_point=False).

        :return: None
        """
        # Test basic direction vector
        direction = np.array([1.0, 2.0, 3.0])
        homog_direction = ps.transformations.generate_homog(direction, False)

        expected = np.array([1.0, 2.0, 3.0, 0.0])
        npt.assert_array_equal(homog_direction, expected)

        # Test unit direction vector
        unit_direction = np.array([1.0, 0.0, 0.0])
        homog_unit = ps.transformations.generate_homog(unit_direction, False)

        expected_unit = np.array([1.0, 0.0, 0.0, 0.0])
        npt.assert_array_equal(homog_unit, expected_unit)

        # Test zero direction vector
        zero_direction = np.array([0.0, 0.0, 0.0])
        homog_zero_dir = ps.transformations.generate_homog(zero_direction, False)

        expected_zero_dir = np.array([0.0, 0.0, 0.0, 0.0])
        npt.assert_array_equal(homog_zero_dir, expected_zero_dir)

    def test_input_validation(self):
        """Tests various input types and edge cases.

        :return: None
        """
        # Test with different numeric types
        float_vector = np.array([1.5, 2.5, 3.5])
        homog_float = ps.transformations.generate_homog(float_vector, True)
        self.assertEqual(homog_float.dtype, np.float64)

        # Test with integer input (should convert to float)
        int_vector = np.array([1, 2, 3])
        homog_int = ps.transformations.generate_homog(int_vector, True)
        self.assertEqual(homog_int.dtype, np.float64)
        expected_int = np.array([1.0, 2.0, 3.0, 1.0])
        npt.assert_array_equal(homog_int, expected_int)


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
                        R = ps.transformations.generate_R(
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
                    R_passive = ps.transformations.generate_R(
                        angles, True, intrinsic, order
                    )
                    R_active = ps.transformations.generate_R(
                        angles, False, intrinsic, order
                    )
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
                    R_intrinsic = ps.transformations.generate_R(
                        angles, passive, True, intrinsic_order
                    )
                    R_extrinsic = ps.transformations.generate_R(
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
                            R = ps.transformations.generate_R(
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

        R_act_x90 = ps.transformations.generate_R(angles_x90, False, True, "123")
        npt.assert_allclose(R_act_x90, R_act_x90_expected, atol=1e-14)

        # Test 90-degree rotation about y-axis (order "123", only second angle)
        angles_y90 = np.array([0.0, 90.0, 0.0])
        R_act_y90_expected = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])

        R_act_y90 = ps.transformations.generate_R(angles_y90, False, True, "123")
        npt.assert_allclose(R_act_y90, R_act_y90_expected, atol=1e-14)

        # Test 90-degree rotation about z-axis (order "123", only third angle)
        angles_z90 = np.array([0.0, 0.0, 90.0])
        R_act_z90_expected = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

        R_act_z90 = ps.transformations.generate_R(angles_z90, False, True, "123")
        npt.assert_allclose(R_act_z90, R_act_z90_expected, atol=1e-14)

    def test_specific_known_passive_rotations(self):
        """Tests specific passive rotations with analytically known results.

        :return: None
        """
        # Test 90-degree rotation about x-axis (order "123", only first angle)
        angles_x90 = np.array([90.0, 0.0, 0.0])
        v_A = np.array([0.0, 1.0, 0.0])
        v_B_expected = np.array([0.0, 0.0, -1.0])

        R_pas_x90 = ps.transformations.generate_R(angles_x90, True, True, "123")
        v_B = R_pas_x90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

        # Test 90-degree rotation about y-axis (order "123", only second angle)
        angles_y90 = np.array([0.0, 90.0, 0.0])
        v_A = np.array([0.0, 0.0, 1.0])
        v_B_expected = np.array([-1.0, 0.0, 0.0])

        R_pas_y90 = ps.transformations.generate_R(angles_y90, True, True, "123")
        v_B = R_pas_y90 @ v_A
        npt.assert_allclose(v_B, v_B_expected, atol=1e-14)

        # Test 90-degree rotation about z-axis (order "123", only third angle)
        angles_z90 = np.array([0.0, 0.0, 90.0])
        v_A = np.array([1.0, 0.0, 0.0])
        v_B_expected = np.array([0.0, -1.0, 0])

        R_pas_z90 = ps.transformations.generate_R(angles_z90, True, True, "123")
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
        R1 = ps.transformations.generate_R(angles1, False, True, "123")
        R2 = ps.transformations.generate_R(angles2, False, True, "123")

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
                        R_large = ps.transformations.generate_R(
                            large_angles, passive, intrinsic, order
                        )
                        R_equivalent = ps.transformations.generate_R(
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
                        R = ps.transformations.generate_R(
                            angles, passive, intrinsic, order
                        )

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
                    R = ps.transformations.generate_R(angles, passive, True, "123")

                    # Should be a valid rotation matrix
                    det = np.linalg.det(R)
                    self.assertAlmostEqual(det, 1.0, places=14)

                    identity_test = R.T @ R
                    npt.assert_allclose(identity_test, np.eye(3), atol=1e-14)


class TestGenerateTRot(unittest.TestCase):
    """This class contains methods for testing the generate_T_rot function.

    This class contains the following public methods:
        test_identity_transformation: Tests transformation with identity rotation.
        test_rotation_preservation: Tests that rotation components are preserved.
        test_with_generate_R_matrices: Tests with matrices from generate_R.
        test_homog_coordinate_transformations: Tests transformation of
        homogeneous coordinates.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def test_identity_transformation(self):
        """Tests transformation with identity rotation matrix.

        :return: None
        """
        # Test with identity rotation
        R_identity = np.eye(3)
        T_identity = ps.transformations.generate_T_rot(R_identity)

        expected_T = np.eye(4)
        npt.assert_array_equal(T_identity, expected_T)

    def test_rotation_preservation(self):
        """Tests that rotation components are correctly preserved.

        :return: None
        """
        # Test with a known rotation matrix
        R_test = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

        T_test = ps.transformations.generate_T_rot(R_test)

        # Check that rotation part is preserved
        npt.assert_array_equal(T_test[:3, :3], R_test)

        # Check that translation part is zero
        expected_translation = np.array([0.0, 0.0, 0.0])
        npt.assert_array_equal(T_test[:3, 3], expected_translation)

        # Check that bottom row is [0, 0, 0, 1]
        expected_bottom = np.array([0.0, 0.0, 0.0, 1.0])
        npt.assert_array_equal(T_test[3, :], expected_bottom)

    def test_with_generate_R_matrices(self):
        """Tests with rotation matrices generated by generate_R.

        :return: None
        """
        angles = np.array([30.0, 45.0, 60.0])

        for passive in [True, False]:
            for intrinsic in [True, False]:
                for order in ["123", "321"]:
                    with self.subTest(
                        passive=passive, intrinsic=intrinsic, order=order
                    ):
                        R = ps.transformations.generate_R(
                            angles, passive, intrinsic, order
                        )
                        T = ps.transformations.generate_T_rot(R)

                        # Test that rotation part matches
                        npt.assert_array_equal(T[:3, :3], R)

                        # Test transformation matrix properties
                        self.assertEqual(T.shape, (4, 4))
                        npt.assert_array_equal(T[3, :], [0, 0, 0, 1])
                        npt.assert_array_equal(T[:3, 3], [0, 0, 0])

    def test_homog_coordinate_transformations(self):
        """Tests transformation of homogeneous coordinates.

        :return: None
        """
        # Test rotation of position vectors
        R_act_z90 = ps.transformations.generate_R(
            np.array([0.0, 0.0, 90.0]), False, True, "123"
        )
        T_rot = ps.transformations.generate_T_rot(R_act_z90)

        # Test position vector [1, 0, 0] -> should become [0, 1, 0]
        pos_vec = np.array([1.0, 0.0, 0.0])
        homog_pos = ps.transformations.generate_homog(pos_vec, True)
        transformed_homog = T_rot @ homog_pos

        expected_transformed = np.array([0.0, 1.0, 0.0, 1.0])
        npt.assert_allclose(transformed_homog, expected_transformed, atol=1e-14)

        # Test direction vector [1, 0, 0] -> should become [0, 1, 0]
        dir_vec = np.array([1.0, 0.0, 0.0])
        homog_dir = ps.transformations.generate_homog(dir_vec, False)
        transformed_homog_dir = T_rot @ homog_dir

        expected_transformed_dir = np.array([0.0, 1.0, 0.0, 0.0])
        npt.assert_allclose(transformed_homog_dir, expected_transformed_dir, atol=1e-14)


class TestGenerateTTrans(unittest.TestCase):
    """This class contains methods for testing the generate_T_trans function.

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

        T_trans_pas = ps.transformations.generate_T_trans(b_A_a, True)
        cHomog_A_a = ps.transformations.generate_homog(c_A_a, True)

        cHomog_A_b = T_trans_pas @ cHomog_A_a
        c_A_b = cHomog_A_b[:3]

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
        T_trans_act = ps.transformations.generate_T_trans(cPrime_A_c, False)
        cHomog_A_a = ps.transformations.generate_homog(c_A_a, True)

        cPrimeHomog_A_a = T_trans_act @ cHomog_A_a
        cPrime_A_a = cPrimeHomog_A_a[:3]

        # For active transformation: new_position = old_position + translation
        expected_cPrime_A_a = c_A_a + cPrime_A_c
        npt.assert_array_equal(cPrime_A_a, expected_cPrime_A_a)

    def test_zero_translation(self):
        """Tests with zero translation vector.

        :return: None
        """
        zero_translation = np.array([0.0, 0.0, 0.0])

        # Test passive with zero translation
        T_passive_zero = ps.transformations.generate_T_trans(zero_translation, True)
        npt.assert_array_equal(T_passive_zero, np.eye(4))

        # Test active with zero translation
        T_active_zero = ps.transformations.generate_T_trans(zero_translation, False)
        npt.assert_array_equal(T_active_zero, np.eye(4))

    def test_transformation_properties(self):
        """Tests properties of transformation matrices.

        :return: None
        """
        translations = np.array([1.0, 2.0, 3.0])

        for passive in [True, False]:
            with self.subTest(passive=passive):
                T = ps.transformations.generate_T_trans(translations, passive)

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
        T_trans_pas = ps.transformations.generate_T_trans(translation, True)
        T_trans_act = ps.transformations.generate_T_trans(translation, False)

        # Test that direction vectors are unaffected by translation
        direction = np.array([1.0, 0.0, 0.0])
        directionHomog = ps.transformations.generate_homog(direction, False)

        passive_transformed_directionHomog = T_trans_pas @ directionHomog
        active_transformed_directionHomog = T_trans_act @ directionHomog

        # Direction vectors should be unchanged by pure translation
        npt.assert_array_equal(passive_transformed_directionHomog[:3], direction)
        npt.assert_array_equal(active_transformed_directionHomog[:3], direction)


class TestGenerateTReflect(unittest.TestCase):
    """This class contains methods for testing the generate_T_reflect function.

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

        T_reflect_act = ps.transformations.generate_T_reflect(
            plane_point, plane_normal, False
        )

        # Point [1, 2, 3] should reflect to [1, 2, -3]
        test_point = np.array([1.0, 2.0, 3.0])
        homog_point = ps.transformations.generate_homog(test_point, True)
        reflected_homog = T_reflect_act @ homog_point
        reflected_point = reflected_homog[:3]

        expected_reflected = np.array([1.0, 2.0, -3.0])
        npt.assert_allclose(reflected_point, expected_reflected, atol=1e-14)

        # Test reflection about xz-plane (y=0, normal=[0,1,0])
        plane_normal_y = np.array([0.0, 1.0, 0.0])
        T_reflect_act_y = ps.transformations.generate_T_reflect(
            plane_point, plane_normal_y, False
        )

        reflected_homog_y = T_reflect_act_y @ homog_point
        reflected_point_y = reflected_homog_y[:3]

        expected_reflected_y = np.array([1.0, -2.0, 3.0])
        npt.assert_allclose(reflected_point_y, expected_reflected_y, atol=1e-14)

    def test_reflection_about_offset_planes(self):
        """Tests reflection about planes not passing through the origin.

        :return: None
        """
        # Test reflection about plane z=2 (normal=[0,0,1], point=[0,0,2])
        plane_point = np.array([0.0, 0.0, 2.0])
        plane_normal = np.array([0.0, 0.0, 1.0])

        T_reflect_act = ps.transformations.generate_T_reflect(
            plane_point, plane_normal, False
        )

        # Point [1, 2, 5] should reflect to [1, 2, -1]
        test_point = np.array([1.0, 2.0, 5.0])
        homog_point = ps.transformations.generate_homog(test_point, True)
        reflected_homog = T_reflect_act @ homog_point
        reflected_point = reflected_homog[:3]

        expected_reflected = np.array([1.0, 2.0, -1.0])
        npt.assert_allclose(reflected_point, expected_reflected, atol=1e-14)

        # Point on the plane should remain unchanged
        plane_test_point = np.array([3.0, 4.0, 2.0])
        homog_plane_point = ps.transformations.generate_homog(plane_test_point, True)
        reflected_plane_homog = T_reflect_act @ homog_plane_point
        reflected_plane_point = reflected_plane_homog[:3]

        npt.assert_allclose(reflected_plane_point, plane_test_point, atol=1e-14)

    def test_zero_length_normal_error(self):
        """Tests error handling for zero-length normal vectors.

        :return: None
        """
        plane_point = np.array([0.0, 0.0, 0.0])
        zero_normal = np.array([0.0, 0.0, 0.0])

        with self.assertRaises(ValueError) as context:
            ps.transformations.generate_T_reflect(plane_point, zero_normal, False)

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
                T = ps.transformations.generate_T_reflect(
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

        T_reflect_act = ps.transformations.generate_T_reflect(
            plane_point, plane_normal, False
        )

        test_point = np.array([5.0, 8.0, -2.0])
        homog_point = ps.transformations.generate_homog(test_point, True)

        # Apply reflection twice
        reflected_once = T_reflect_act @ homog_point
        reflected_twice = T_reflect_act @ reflected_once

        final_point = reflected_twice[:3]
        npt.assert_allclose(final_point, test_point, atol=1e-14)

    def test_normal_vector_normalization(self):
        """Tests that non-unit normal vectors are correctly normalized.

        :return: None
        """
        plane_point = np.array([0.0, 0.0, 0.0])
        # Use non-unit normal vector - should be normalized internally
        plane_normal = np.array([0.0, 0.0, 5.0])  # Should normalize to [0,0,1]

        T_reflect_act = ps.transformations.generate_T_reflect(
            plane_point, plane_normal, False
        )

        test_point = np.array([1.0, 2.0, 3.0])
        homog_point = ps.transformations.generate_homog(test_point, True)
        reflected_homog = T_reflect_act @ homog_point
        reflected_point = reflected_homog[:3]

        # Should give same result as unit normal [0,0,1]
        expected_reflected = np.array([1.0, 2.0, -3.0])
        npt.assert_allclose(reflected_point, expected_reflected, atol=1e-14)

    def test_passive_active_equivalence(self):
        """Tests that passive and active modes produce identical matrices.

        :return: None
        """
        plane_point = np.array([1.0, 2.0, 3.0])
        plane_normal = np.array([1.0, 0.0, 1.0])

        T_reflect_pas = ps.transformations.generate_T_reflect(
            plane_point, plane_normal, True
        )
        T_reflect_act = ps.transformations.generate_T_reflect(
            plane_point, plane_normal, False
        )

        # Matrices should be identical for reflections
        npt.assert_array_equal(T_reflect_pas, T_reflect_act)


if __name__ == "__main__":
    unittest.main()
