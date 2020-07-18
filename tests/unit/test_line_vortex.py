""" This module contains a class to test line vortex objects.

This module contains the following classes:
    TestLineVortex: This class contains methods for testing line vortex objects.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import numpy as np

import main as main
import tests.unit.fixtures.vortex_fixtures


class TestLineVortex(unittest.TestCase):
    """ This class contains methods for testing line vortex objects.

    This class contains the following public methods:
        setUp: This method is automatically called before each testing method to set up the fixtures.
        tearDown: This method is automatically called before each testing method to tear down the fixtures.
        test_class: This method tests the class's instantiation.
        test_calculate_normalized_induced_velocity: This method tests the calculation of normalized induced velocity.
        test_calculate_induced_velocity: This method tests the calculation of induced velocity.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
        """ This method is automatically called before each testing method to set up the fixtures.

        :return: None
        """

        # Create the constructing fixtures.
        self.line_vortex_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_line_vortex_fixture()
        )
        self.origin_fixture = tests.unit.fixtures.vortex_fixtures.make_origin_fixture()
        self.termination_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_termination_fixture()
        )
        self.strength_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_strength_fixture()
        )

    def tearDown(self):
        """ This method is automatically called before each testing method to tear down the fixtures.

        :return: None
        """

        # Delete the constructing fixtures.
        del self.line_vortex_fixture
        del self.origin_fixture
        del self.termination_fixture
        del self.strength_fixture

    def test_class(self):
        """ This method tests the class's instantiation.

        :return:
        """

        # Test that the object is of the right type.
        self.assertIsInstance(self.line_vortex_fixture, main.aerodynamics.LineVortex)

        # Test that the vortex's coordinates were correctly set.
        self.assertTrue(
            np.allclose(self.line_vortex_fixture.origin, self.origin_fixture)
        )
        self.assertTrue(
            np.allclose(self.line_vortex_fixture.termination, self.termination_fixture)
        )

        # Test that the vortex's strength was correctly set.
        self.assertEqual(self.line_vortex_fixture.strength, self.strength_fixture)

        # Test that the vortex's center and vector were correctly calculated.
        self.assertTrue(
            np.allclose(self.line_vortex_fixture.center, np.array([0.5, 0.5, 0.5]))
        )
        self.assertTrue(
            np.allclose(self.line_vortex_fixture.vector, np.array([1, 1, 1]))
        )

    def test_calculate_normalized_induced_velocity(self):
        """ This method tests the calculation of normalized induced velocity.

        :return: None
        """

        # Test the velocity is correctly calculated at the origin.
        point_fixture = self.line_vortex_fixture.origin
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated at the termination.
        point_fixture = self.line_vortex_fixture.termination
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated at the center.
        point_fixture = self.line_vortex_fixture.center
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated in along a line collinear with the vortex, but off to one side.
        point_fixture = (
            self.line_vortex_fixture.center + self.line_vortex_fixture.vector
        )
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated in along a line collinear with the vortex, but off to the other
        # side.
        point_fixture = (
            self.line_vortex_fixture.center - self.line_vortex_fixture.vector
        )
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Create a point fixture which does not have zero induced velocity on it.
        point_fixture = self.line_vortex_fixture.center + np.array([1, 0, 0])

        # Perform the normalized induced velocity calculation. See the calculate_normalized_induced_velocity method in
        # the LineVortex object in the aerodynamics module for references.
        r_1_fixture = point_fixture - self.line_vortex_fixture.origin
        r_2_fixture = point_fixture - self.line_vortex_fixture.termination
        r_0_fixture = r_1_fixture - r_2_fixture
        r_1_length_fixture = np.linalg.norm(r_1_fixture)
        r_2_length_fixture = np.linalg.norm(r_2_fixture)

        r_0_dot_r_1_fixture = np.dot(r_0_fixture, r_1_fixture)
        r_0_dot_r_2_fixture = np.dot(r_0_fixture, r_2_fixture)

        r_1_cross_r_2_fixture = np.cross(r_1_fixture, r_2_fixture)
        r_1_cross_r_2_absolute_magnitude_fixture = (
            r_1_cross_r_2_fixture[0] ** 2
            + r_1_cross_r_2_fixture[1] ** 2
            + r_1_cross_r_2_fixture[2] ** 2
        )

        normalized_induced_velocity = (
            r_1_cross_r_2_fixture
            / (4 * np.pi * r_1_cross_r_2_absolute_magnitude_fixture)
            * (
                r_0_dot_r_1_fixture / r_1_length_fixture
                - r_0_dot_r_2_fixture / r_2_length_fixture
            )
        )

        # Check that the calculate_normalized_induced_velocity method's calculations produce the expected results.
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                normalized_induced_velocity,
            )
        )

    def test_calculate_induced_velocity(self):
        """ This method tests the calculation of induced velocity.

        :return: None
        """

        # Test the velocity is correctly calculated at the origin.
        point_fixture = self.line_vortex_fixture.origin
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated at the termination.
        point_fixture = self.line_vortex_fixture.termination
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated at the center.
        point_fixture = self.line_vortex_fixture.center
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated in along a line collinear with the vortex, but off to one side.
        point_fixture = (
            self.line_vortex_fixture.center + self.line_vortex_fixture.vector
        )
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Test the velocity is correctly calculated in along a line collinear with the vortex, but off to the other
        # side.
        point_fixture = (
            self.line_vortex_fixture.center - self.line_vortex_fixture.vector
        )
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        # Create a point fixture which does not have zero induced velocity on it.
        point_fixture = self.line_vortex_fixture.center + np.array([1, 0, 0])

        # Perform the induced velocity calculation. See the calculate_induced_velocity method in the LineVortex object
        # in the aerodynamics module for references.
        r_1_fixture = self.line_vortex_fixture.origin - point_fixture
        r_2_fixture = self.line_vortex_fixture.termination - point_fixture
        r_0_fixture = r_1_fixture - r_2_fixture
        r_1_length_fixture = np.linalg.norm(r_1_fixture)
        r_2_length_fixture = np.linalg.norm(r_2_fixture)

        r_0_dot_r_1_fixture = np.dot(r_0_fixture, r_1_fixture)
        r_0_dot_r_2_fixture = np.dot(r_0_fixture, r_2_fixture)

        r_1_cross_r_2_fixture = np.cross(r_1_fixture, r_2_fixture)
        r_1_cross_r_2_absolute_magnitude_fixture = (
            r_1_cross_r_2_fixture[0] ** 2
            + r_1_cross_r_2_fixture[1] ** 2
            + r_1_cross_r_2_fixture[2] ** 2
        )

        induced_velocity = (
            r_1_cross_r_2_fixture
            * self.line_vortex_fixture.strength
            / (4 * np.pi * r_1_cross_r_2_absolute_magnitude_fixture)
            * (
                r_0_dot_r_1_fixture / r_1_length_fixture
                - r_0_dot_r_2_fixture / r_2_length_fixture
            )
        )

        # Check that the calculate_induced_velocity method's calculations produce the expected results.
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                induced_velocity,
            )
        )
