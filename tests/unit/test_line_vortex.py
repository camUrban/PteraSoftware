# ToDo: Properly document this module.
"""

"""

import unittest
import numpy as np

import aviansoftwareminimumviableproduct as asmvp
import tests.unit


# ToDo: Properly document this class.
class TestLineVortex(unittest.TestCase):
    """

    """

    # ToDo: Properly document this method.
    def setUp(self):
        """

        :return: 
        """

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

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return: 
        """

        del self.line_vortex_fixture
        del self.origin_fixture
        del self.termination_fixture
        del self.strength_fixture

    # ToDo: Properly document this method.
    def test_class(self):
        """

        :return:
        """

        self.assertIsInstance(self.line_vortex_fixture, asmvp.aerodynamics.LineVortex)

        self.assertTrue(
            np.allclose(self.line_vortex_fixture.origin, self.origin_fixture)
        )
        self.assertTrue(
            np.allclose(self.line_vortex_fixture.termination, self.termination_fixture)
        )
        self.assertEqual(self.line_vortex_fixture.strength, self.strength_fixture)

        self.assertTrue(
            np.allclose(self.line_vortex_fixture.center, np.array([0.5, 0.5, 0.5]))
        )
        self.assertTrue(
            np.allclose(self.line_vortex_fixture.vector, np.array([1, 1, 1]))
        )

    # ToDo: Properly document this method.
    def test_calculate_normalized_induced_velocity(self):
        """

        :return:
        """

        point_fixture = self.line_vortex_fixture.origin
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        point_fixture = self.line_vortex_fixture.termination
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        point_fixture = self.line_vortex_fixture.center
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

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

        point_fixture = self.line_vortex_fixture.center + np.array([1, 0, 0])

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

        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                normalized_induced_velocity,
            )
        )

    # ToDo: Properly document this method.
    def test_calculate_induced_velocity(self):
        """

        :return:
        """

        point_fixture = self.line_vortex_fixture.origin
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        point_fixture = self.line_vortex_fixture.termination
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

        point_fixture = self.line_vortex_fixture.center
        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                0,
            )
        )

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

        point_fixture = self.line_vortex_fixture.center + np.array([1, 0, 0])

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

        self.assertTrue(
            np.allclose(
                self.line_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                induced_velocity,
            )
        )
