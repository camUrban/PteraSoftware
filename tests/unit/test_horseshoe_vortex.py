"""This module contains a class to test horseshoe vortex objects.

This module contains the following classes:
    TestHorseshoeVortex: This class contains methods for testing horseshoe vortex
    objects.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import vortex_fixtures


class TestHorseshoeVortex(unittest.TestCase):
    """This class contains methods for testing horseshoe vortex objects.

    This class contains the following public methods:
        setUp: This method is automatically called before each testing method to set
        up the fixtures.

        tearDown: This method is automatically called before each testing method to
        tear down the fixtures.

        test_class: This method tests the class's instantiation.

        test_update_strength: This method tests the update_strength method.

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

        # Get the constructing fixtures.
        self.horseshoe_vortex_fixture = vortex_fixtures.make_horseshoe_vortex_fixture()
        self.origin_fixture = vortex_fixtures.make_origin_fixture()
        self.termination_fixture = vortex_fixtures.make_termination_fixture()
        self.strength_fixture = vortex_fixtures.make_strength_fixture()
        self.infinite_leg_direction_fixture = (
            vortex_fixtures.make_infinite_leg_direction_fixture()
        )
        self.infinite_leg_length_fixture = (
            vortex_fixtures.make_infinite_leg_length_fixture()
        )

    def tearDown(self):
        """This method is automatically called before each testing method to tear
        down the fixtures.

        :return: None
        """

        # Delete the constructing fixtures.
        del self.horseshoe_vortex_fixture
        del self.origin_fixture
        del self.termination_fixture
        del self.strength_fixture
        del self.infinite_leg_direction_fixture
        del self.infinite_leg_length_fixture

    def test_class(self):
        """This method tests the class's instantiation.

        :return: None
        """

        # Test that the objects are all the right type.
        self.assertIsInstance(
            self.horseshoe_vortex_fixture,
            ps.aerodynamics.HorseshoeVortex,
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.finite_leg,
            ps.aerodynamics._LineVortex,
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.left_leg,
            ps.aerodynamics._LineVortex,
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.right_leg,
            ps.aerodynamics._LineVortex,
        )

        # Test that the vortex objects' coordinates were correctly set.
        self.assertTrue(
            np.allclose(self.horseshoe_vortex_fixture.Frhvp_G_Cg, self.origin_fixture)
        )
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.Flhvp_G_Cg,
                self.termination_fixture,
            )
        )

        # Test that the horseshoe vortex object's strength was set correctly.
        self.assertEqual(self.horseshoe_vortex_fixture.strength, self.strength_fixture)

        # Test that other class attributes were correctly set.
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.leftLegVector_G,
                self.infinite_leg_direction_fixture,
            )
        )
        self.assertEqual(
            self.horseshoe_vortex_fixture.left_right_leg_lengths,
            self.infinite_leg_length_fixture,
        )

        # Test that the infinite legs' coordinates are correct.
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.Brhvp_G_Cg,
                self.horseshoe_vortex_fixture.Frhvp_G_Cg
                + self.horseshoe_vortex_fixture.leftLegVector_G
                * self.horseshoe_vortex_fixture.left_right_leg_lengths,
            )
        )
        self.assertTrue(
            np.allclose(
                self.horseshoe_vortex_fixture.Blhvp_G_Cg,
                self.horseshoe_vortex_fixture.Flhvp_G_Cg
                + self.horseshoe_vortex_fixture.leftLegVector_G
                * self.horseshoe_vortex_fixture.left_right_leg_lengths,
            )
        )

    def test_update_strength(self):
        """This method tests the update_strength method.

        :return: None
        """

        # Create fixtures to hold the current and new strengths.
        old_strength_fixture = self.horseshoe_vortex_fixture.strength
        new_strength_fixture = old_strength_fixture * 5 + 1

        # Update the horseshoe vortex fixture's strength.
        self.horseshoe_vortex_fixture.update_strength(strength=new_strength_fixture)

        # Test that all the strength's have been updated correctly.
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

        # Revert the change.
        self.horseshoe_vortex_fixture.update_strength(strength=old_strength_fixture)
