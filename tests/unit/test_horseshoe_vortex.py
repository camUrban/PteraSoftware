""" This module contains a class to test horseshoe vortex objects.

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

import src
import tests.unit.fixtures.vortex_fixtures


class TestHorseshoeVortex(unittest.TestCase):
    """This class contains methods for testing horseshoe vortex objects.

    This class contains the following public methods:
        setUp: This method is automatically called before each testing method to set
        up the fixtures.

        tearDown: This method is automatically called before each testing method to
        tear down the fixtures.

        test_class: This method tests the class's instantiation.

        test_calculate_normalized_induced_velocity: This method tests the calculation
        of normalized induced velocity.

        test_calculate_induced_velocity: This method tests the calculation of induced
        velocity.

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

        # Test that the objects are all of the right type.
        self.assertIsInstance(
            self.horseshoe_vortex_fixture, src.aerodynamics.HorseshoeVortex
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.finite_leg, src.aerodynamics.LineVortex
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.left_leg, src.aerodynamics.LineVortex
        )
        self.assertIsInstance(
            self.horseshoe_vortex_fixture.right_leg, src.aerodynamics.LineVortex
        )

        # Test that the vortex objects' coordinates were correctly set.
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

        # Test that the horseshoe vortex object's strength was set correctly.
        self.assertEqual(self.horseshoe_vortex_fixture.strength, self.strength_fixture)

        # Test that other class attributes were correctly set.
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

        # Test that the infinite legs' coordinates are correct.
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

    def test_calculate_normalized_induced_velocity(self):
        """This method tests the calculation of normalized induced velocity.

        :return: None
        """

        # Test the velocity is correctly calculated at a point on the finite leg.
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

        # Test the velocity is correctly calculated at a point on the right leg.
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

        # Test the velocity is correctly calculated at a point on the left leg.
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

    def test_calculate_induced_velocity(self):
        """This method tests the calculation of induced velocity.

        :return: None
        """

        # Test the velocity is correctly calculated at a point on the finite leg.
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

        # Test the velocity is correctly calculated at a point on the right leg.
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

        # Test the velocity is correctly calculated at a point on the left leg.
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
