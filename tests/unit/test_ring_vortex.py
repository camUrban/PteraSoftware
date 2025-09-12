"""This module contains a class to test ring vortex objects.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    TestRingVortex: This is a class with functions to test ring vortex objects.
"""

import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import vortex_fixtures


class TestRingVortex(unittest.TestCase):
    """This is a class with functions to test ring vortex objects.

    This class contains the following public methods:
        setUp: This method is automatically called before each testing method to set
        up the fixtures.

        tearDown: This method is automatically called before each testing method to
        tear down the fixtures.

        test_class: This method tests the class's instantiation.

        test_update_strength: This method tests the update_strength method.

        test_update_position: This method tests the update_position method.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
        """This method is automatically called before each testing method to set up
        the fixtures.

        :return: None
        """

        # Set up the constructing fixtures.
        self.ring_vortex_fixture = vortex_fixtures.make_ring_vortex_fixture()
        self.strength_fixture = vortex_fixtures.make_strength_fixture()
        self.front_left_vertex_fixture = (
            vortex_fixtures.make_front_left_vertex_fixture()
        )
        self.front_right_vertex_fixture = (
            vortex_fixtures.make_front_right_vertex_fixture()
        )
        self.back_left_vertex_fixture = vortex_fixtures.make_back_left_vertex_fixture()
        self.back_right_vertex_fixture = (
            vortex_fixtures.make_back_right_vertex_fixture()
        )

    def tearDown(self):
        """This method is automatically called before each testing method to tear
        down the fixtures.

        :return: None
        """

        # Delete the constructing fixtures.
        del self.ring_vortex_fixture
        del self.strength_fixture
        del self.front_left_vertex_fixture
        del self.front_right_vertex_fixture
        del self.back_left_vertex_fixture
        del self.back_right_vertex_fixture

    def test_class(self):
        """This method tests the class's instantiation.

        :return: None
        """

        # Test that the objects are all the right type.
        self.assertIsInstance(self.ring_vortex_fixture, ps.aerodynamics.RingVortex)
        self.assertIsInstance(
            self.ring_vortex_fixture.front_leg,
            ps.aerodynamics.LineVortex,
        )
        self.assertIsInstance(
            self.ring_vortex_fixture.left_leg, ps.aerodynamics.LineVortex
        )
        self.assertIsInstance(
            self.ring_vortex_fixture.back_leg, ps.aerodynamics.LineVortex
        )
        self.assertIsInstance(
            self.ring_vortex_fixture.right_leg,
            ps.aerodynamics.LineVortex,
        )

        # Test that the vortex objects' coordinates were correctly set.
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.front_leg.origin,
                self.front_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.front_leg.termination,
                self.front_left_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.left_leg.origin, self.front_left_vertex_fixture
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.left_leg.termination,
                self.back_left_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.back_leg.origin, self.back_left_vertex_fixture
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.back_leg.termination,
                self.back_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.right_leg.origin,
                self.back_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.right_leg.termination,
                self.front_right_vertex_fixture,
            )
        )

        # Test that the vortex objects' strengths were correctly set.
        self.assertEqual(self.ring_vortex_fixture.strength, self.strength_fixture)
        self.assertEqual(
            self.ring_vortex_fixture.front_leg.strength, self.strength_fixture
        )
        self.assertEqual(
            self.ring_vortex_fixture.left_leg.strength, self.strength_fixture
        )
        self.assertEqual(
            self.ring_vortex_fixture.back_leg.strength, self.strength_fixture
        )
        self.assertEqual(
            self.ring_vortex_fixture.right_leg.strength, self.strength_fixture
        )

    def test_update_strength(self):
        """This method tests the update_strength method.

        :return: None
        """

        # Create fixtures to hold the current and new strengths.
        old_strength_fixture = self.ring_vortex_fixture.strength
        new_strength_fixture = old_strength_fixture * 5 + 1

        # Update the ring vortex fixture's strength.
        self.ring_vortex_fixture.update_strength(strength=new_strength_fixture)

        # Test that all the strength's have been updated correctly.
        self.assertEqual(self.ring_vortex_fixture.strength, new_strength_fixture)
        self.assertEqual(
            self.ring_vortex_fixture.front_leg.strength, new_strength_fixture
        )
        self.assertEqual(
            self.ring_vortex_fixture.left_leg.strength, new_strength_fixture
        )
        self.assertEqual(
            self.ring_vortex_fixture.back_leg.strength, new_strength_fixture
        )
        self.assertEqual(
            self.ring_vortex_fixture.right_leg.strength, new_strength_fixture
        )

        # Revert the change.
        self.ring_vortex_fixture.update_strength(strength=old_strength_fixture)

    def test_update_position(self):
        """This method tests the update_position method.

        :return: None
        """

        # Create fixtures to hold the old values of the ring vortex's position.
        old_front_right_vertex_fixture = self.ring_vortex_fixture.front_right_vertex
        old_front_left_vertex_fixture = self.ring_vortex_fixture.front_left_vertex
        old_back_left_vertex_fixture = self.ring_vortex_fixture.back_left_vertex
        old_back_right_vertex_fixture = self.ring_vortex_fixture.back_right_vertex

        # Create fixtures to hold the soon-to-be new values of the ring vortex's
        # position.
        new_front_right_vertex_fixture = old_front_right_vertex_fixture + np.array(
            [1, 0, 0]
        )
        new_front_left_vertex_fixture = old_front_left_vertex_fixture + np.array(
            [1, 0, 0]
        )
        new_back_left_vertex_fixture = old_back_left_vertex_fixture + np.array(
            [1, 0, 0]
        )
        new_back_right_vertex_fixture = old_back_right_vertex_fixture + np.array(
            [1, 0, 0]
        )

        # Update the ring vortex fixture's position.
        self.ring_vortex_fixture.update_position(
            front_right_vertex=new_front_right_vertex_fixture,
            front_left_vertex=new_front_left_vertex_fixture,
            back_left_vertex=new_back_left_vertex_fixture,
            back_right_vertex=new_back_right_vertex_fixture,
        )

        # Test that the position of the ring vortex was correctly updated.
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.front_right_vertex,
                new_front_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.front_left_vertex,
                new_front_left_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.back_left_vertex, new_back_left_vertex_fixture
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.back_right_vertex,
                new_back_right_vertex_fixture,
            )
        )

        # Check that the positions of the child objects have been correctly updated.
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.front_leg.origin,
                new_front_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.front_leg.termination,
                new_front_left_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.left_leg.origin, new_front_left_vertex_fixture
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.left_leg.termination,
                new_back_left_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.back_leg.origin, new_back_left_vertex_fixture
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.back_leg.termination,
                new_back_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.right_leg.origin,
                new_back_right_vertex_fixture,
            )
        )
        self.assertTrue(
            np.allclose(
                self.ring_vortex_fixture.right_leg.termination,
                new_front_right_vertex_fixture,
            )
        )

        # Revert the changes.
        self.ring_vortex_fixture.update_position(
            front_right_vertex=old_front_right_vertex_fixture,
            front_left_vertex=old_front_left_vertex_fixture,
            back_left_vertex=old_back_left_vertex_fixture,
            back_right_vertex=old_back_right_vertex_fixture,
        )
