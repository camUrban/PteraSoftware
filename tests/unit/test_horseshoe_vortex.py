# ToDo: Properly document this module.
"""

"""

import unittest

import numpy as np

import aviansoftwareminimumviableproduct as asmvp
import tests.unit


# ToDo: Properly document this class.
class TestHorseshoeVortex(unittest.TestCase):
    """

    """

    # ToDo: Properly document this method.
    def setUp(self):
        """
        
        :return: 
        """

        self.horseshoe_vortex_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_horseshoe_vortex_fixture()
        )
        self.origin_fixture = tests.unit.fixtures.vortex_fixtures.make_origin_fixture()
        self.termination_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_termination_fixture()
        )
        self.strength_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_strength_fixture()
        )
        self.infinite_leg_direction_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_infinite_leg_direction_fixture()
        )
        self.infinite_leg_length_fixture = (
            tests.unit.fixtures.vortex_fixtures.make_infinite_leg_length_fixture()
        )

    # ToDo: Properly document this method.
    def tearDown(self):
        """
        
        :return: 
        """

        del self.horseshoe_vortex_fixture
        del self.origin_fixture
        del self.termination_fixture
        del self.strength_fixture
        del self.infinite_leg_direction_fixture
        del self.infinite_leg_length_fixture

    # ToDo: Properly document this method.
    def test_class(self):
        """

        :return:
        """

        self.assertIsInstance(
            self.horseshoe_vortex_fixture, asmvp.aerodynamics.HorseshoeVortex
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.finite_leg, asmvp.aerodynamics.LineVortex
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.left_leg, asmvp.aerodynamics.LineVortex
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.right_leg, asmvp.aerodynamics.LineVortex
        )

        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.finite_leg_origin, self.origin_fixture
            )
        )
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.finite_leg_termination,
                self.termination_fixture,
            )
        )
        self.assertEqual(self.horseshoe_vortex_fixture.strength, self.strength_fixture)
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.infinite_leg_direction,
                self.infinite_leg_direction_fixture,
            )
        )
        self.assertEqual(
            self.horseshoe_vortex_fixture.infinite_leg_length,
            self.infinite_leg_length_fixture,
        )

        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.right_leg_origin,
                self.horseshoe_vortex_fixture.finite_leg_origin
                + self.horseshoe_vortex_fixture.infinite_leg_direction
                * self.horseshoe_vortex_fixture.infinite_leg_length,
            )
        )
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.left_leg_termination,
                self.horseshoe_vortex_fixture.finite_leg_termination
                + self.horseshoe_vortex_fixture.infinite_leg_direction
                * self.horseshoe_vortex_fixture.infinite_leg_length,
            )
        )

    # ToDo: Properly document this method.
    def test_calculate_normalized_induced_velocity(self):
        """

        :return:
        """

        point_fixture = self.horseshoe_vortex_fixture.finite_leg.center
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                self.horseshoe_vortex_fixture.right_leg.calculate_normalized_induced_velocity(
                    point=point_fixture
                )
                + self.horseshoe_vortex_fixture.left_leg.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
            )
        )

        point_fixture = self.horseshoe_vortex_fixture.right_leg.center
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                self.horseshoe_vortex_fixture.finite_leg.calculate_normalized_induced_velocity(
                    point=point_fixture
                )
                + self.horseshoe_vortex_fixture.left_leg.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
            )
        )

        point_fixture = self.horseshoe_vortex_fixture.left_leg.center
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
                self.horseshoe_vortex_fixture.right_leg.calculate_normalized_induced_velocity(
                    point=point_fixture
                )
                + self.horseshoe_vortex_fixture.finite_leg.calculate_normalized_induced_velocity(
                    point=point_fixture
                ),
            )
        )

    # ToDo: Properly document this method.
    def test_calculate_induced_velocity(self):
        """

        :return:
        """

        point_fixture = self.horseshoe_vortex_fixture.finite_leg.center
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                self.horseshoe_vortex_fixture.right_leg.calculate_induced_velocity(
                    point=point_fixture
                )
                + self.horseshoe_vortex_fixture.left_leg.calculate_induced_velocity(
                    point=point_fixture
                ),
            )
        )

        point_fixture = self.horseshoe_vortex_fixture.right_leg.center
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                self.horseshoe_vortex_fixture.finite_leg.calculate_induced_velocity(
                    point=point_fixture
                )
                + self.horseshoe_vortex_fixture.left_leg.calculate_induced_velocity(
                    point=point_fixture
                ),
            )
        )

        point_fixture = self.horseshoe_vortex_fixture.left_leg.center
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.calculate_induced_velocity(
                    point=point_fixture
                ),
                self.horseshoe_vortex_fixture.right_leg.calculate_induced_velocity(
                    point=point_fixture
                )
                + self.horseshoe_vortex_fixture.finite_leg.calculate_induced_velocity(
                    point=point_fixture
                ),
            )
        )

    # ToDo: Properly document this method.
    def test_update_strength(self):
        """

        :return:
        """

        old_strength_fixture = self.horseshoe_vortex_fixture.strength

        new_strength_fixture = old_strength_fixture * 5 + 1

        self.horseshoe_vortex_fixture.update_strength(strength=new_strength_fixture)

        self.assertEqual(self.horseshoe_vortex_fixture.strength, new_strength_fixture)
        self.assertEqual(
            self.horseshoe_vortex_fixture.right_leg.strength, new_strength_fixture
        )
        self.assertEqual(
            self.horseshoe_vortex_fixture.finite_leg.strength, new_strength_fixture
        )
        self.assertEqual(
            self.horseshoe_vortex_fixture.left_leg.strength, new_strength_fixture
        )

        self.horseshoe_vortex_fixture.update_strength(strength=old_strength_fixture)
