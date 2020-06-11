
# ToDo: Properly document this module.
"""

"""

import unittest
import numpy as np

import aviansoftwareminimumviableproduct as asmvp
import tests.unit


# ToDo: Properly document this class.
class TestRingVortex(unittest.TestCase):
    """

    """

    # ToDo: Properly document this method.
    def setUp(self):
        """

        :return: 
        """

        self.ring_vortex_fixture = tests.unit.fixtures.vortex_fixtures.make_ring_vortex_fixture()
        self.strength_fixture = tests.unit.fixtures.vortex_fixtures.make_strength_fixture()
        self.front_left_vertex_fixture = tests.unit.fixtures.vortex_fixtures.make_front_left_vertex_fixture()
        self.front_right_vertex_fixture = tests.unit.fixtures.vortex_fixtures.make_front_right_vertex_fixture()
        self.back_left_vertex_fixture = tests.unit.fixtures.vortex_fixtures.make_back_left_vertex_fixture()
        self.back_right_vertex_fixture = tests.unit.fixtures.vortex_fixtures.make_back_right_vertex_fixture()

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return: 
        """

        del self.ring_vortex_fixture
        del self.strength_fixture
        del self.front_left_vertex_fixture
        del self.front_right_vertex_fixture
        del self.back_left_vertex_fixture
        del self.back_right_vertex_fixture
    
    # ToDo: Properly document this method.
    def test_class(self):
        """

        :return:
        """

        self.assertIsInstance(self.ring_vortex_fixture, asmvp.aerodynamics.RingVortex)
        self.assertIsInstance(self.ring_vortex_fixture.front_leg, asmvp.aerodynamics.LineVortex)
        self.assertIsInstance(self.ring_vortex_fixture.left_leg, asmvp.aerodynamics.LineVortex)
        self.assertIsInstance(self.ring_vortex_fixture.back_leg, asmvp.aerodynamics.LineVortex)
        self.assertIsInstance(self.ring_vortex_fixture.right_leg, asmvp.aerodynamics.LineVortex)

        self.assertTrue(np.allclose(self.ring_vortex_fixture.front_leg.origin, self.front_right_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.front_leg.termination, self.front_left_vertex_fixture))
        
        self.assertTrue(np.allclose(self.ring_vortex_fixture.left_leg.origin, self.front_left_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.left_leg.termination, self.back_left_vertex_fixture))
        
        self.assertTrue(np.allclose(self.ring_vortex_fixture.back_leg.origin, self.back_left_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.back_leg.termination, self.back_right_vertex_fixture))

        self.assertTrue(np.allclose(self.ring_vortex_fixture.right_leg.origin, self.back_right_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.right_leg.termination, self.front_right_vertex_fixture))

        self.assertEqual(self.ring_vortex_fixture.strength, self.strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.front_leg.strength, self.strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.left_leg.strength, self.strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.back_leg.strength, self.strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.right_leg.strength, self.strength_fixture)

    # ToDo: Properly document this method.
    def test_calculate_normalized_induced_velocity(self):
        """

        :return:
        """

        point_fixture = self.ring_vortex_fixture.front_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_normalized_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.left_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.back_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.right_leg.calculate_normalized_induced_velocity(point=point_fixture)
            )
        )
        
        point_fixture = self.ring_vortex_fixture.left_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_normalized_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.front_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.back_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.right_leg.calculate_normalized_induced_velocity(point=point_fixture)
            )
        )
        
        point_fixture = self.ring_vortex_fixture.back_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_normalized_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.front_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.left_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.right_leg.calculate_normalized_induced_velocity(point=point_fixture)
            )
        )

        point_fixture = self.ring_vortex_fixture.right_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_normalized_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.front_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.left_leg.calculate_normalized_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.back_leg.calculate_normalized_induced_velocity(point=point_fixture)
            )
        )

    # ToDo: Properly document this method.
    def test_calculate_induced_velocity(self):
        """

        :return:
        """

        point_fixture = self.ring_vortex_fixture.front_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.left_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.back_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.right_leg.calculate_induced_velocity(point=point_fixture)
            )
        )

        point_fixture = self.ring_vortex_fixture.left_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.front_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.back_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.right_leg.calculate_induced_velocity(point=point_fixture)
            )
        )

        point_fixture = self.ring_vortex_fixture.back_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.front_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.left_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.right_leg.calculate_induced_velocity(point=point_fixture)
            )
        )

        point_fixture = self.ring_vortex_fixture.right_leg.center
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.calculate_induced_velocity(point=point_fixture),
                self.ring_vortex_fixture.front_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.left_leg.calculate_induced_velocity(point=point_fixture)
                + self.ring_vortex_fixture.back_leg.calculate_induced_velocity(point=point_fixture)
            )
        )

    # ToDo: Properly document this method.
    def test_update_strength(self):
        """

        :return:
        """

        old_strength_fixture = self.ring_vortex_fixture.strength

        new_strength_fixture = old_strength_fixture * 5 + 1

        self.ring_vortex_fixture.update_strength(strength=new_strength_fixture)

        self.assertEqual(self.ring_vortex_fixture.strength, new_strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.front_leg.strength, new_strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.left_leg.strength, new_strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.back_leg.strength, new_strength_fixture)
        self.assertEqual(self.ring_vortex_fixture.right_leg.strength, new_strength_fixture)

        self.ring_vortex_fixture.update_strength(strength=old_strength_fixture)

    # ToDo: Properly document this method.
    def test_update_position(self):
        """

        :return:
        """

        old_front_right_vertex_fixture = self.ring_vortex_fixture.front_right_vertex
        old_front_left_vertex_fixture = self.ring_vortex_fixture.front_right_vertex
        old_back_left_vertex_fixture = self.ring_vortex_fixture.back_left_vertex
        old_back_right_vertex_fixture = self.ring_vortex_fixture.back_right_vertex

        new_front_right_vertex_fixture = old_front_right_vertex_fixture + np.array([1, 0, 0])
        new_front_left_vertex_fixture = old_front_left_vertex_fixture + np.array([1, 0, 0])
        new_back_left_vertex_fixture = old_back_left_vertex_fixture + np.array([1, 0, 0])
        new_back_right_vertex_fixture = old_back_right_vertex_fixture + np.array([1, 0, 0])

        self.ring_vortex_fixture.update_position(
            front_right_vertex=new_front_right_vertex_fixture,
            front_left_vertex=new_front_left_vertex_fixture,
            back_left_vertex=new_back_left_vertex_fixture,
            back_right_vertex=new_back_right_vertex_fixture
        )

        self.assertTrue(np.allclose(self.ring_vortex_fixture.front_right_vertex, new_front_right_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.front_left_vertex, new_front_left_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.back_left_vertex, new_back_left_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.back_right_vertex, new_back_right_vertex_fixture))

        self.assertTrue(np.allclose(self.ring_vortex_fixture.front_leg.origin, new_front_right_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.front_leg.termination, new_front_left_vertex_fixture))

        self.assertTrue(np.allclose(self.ring_vortex_fixture.left_leg.origin, new_front_left_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.left_leg.termination, new_back_left_vertex_fixture))

        self.assertTrue(np.allclose(self.ring_vortex_fixture.back_leg.origin, new_back_left_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.back_leg.termination, new_back_right_vertex_fixture))

        self.assertTrue(np.allclose(self.ring_vortex_fixture.right_leg.origin, new_back_right_vertex_fixture))
        self.assertTrue(np.allclose(self.ring_vortex_fixture.right_leg.termination, new_front_right_vertex_fixture))

        self.ring_vortex_fixture.update_position(
            front_right_vertex=old_front_right_vertex_fixture,
            front_left_vertex=old_front_left_vertex_fixture,
            back_left_vertex=old_back_left_vertex_fixture,
            back_right_vertex=old_back_right_vertex_fixture
        )
