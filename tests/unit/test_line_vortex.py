"""This module contains a class to test line vortex objects.

This module contains the following classes:
    TestLineVortex: This class contains methods for testing line vortex objects.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import vortex_fixtures


class TestLineVortex(unittest.TestCase):
    """This class contains methods for testing line vortex objects.

    This class contains the following public methods:
        setUp: This method is automatically called before each testing method to set
        up the fixtures.

        tearDown: This method is automatically called before each testing method to
        tear down the fixtures.

        test_class: This method tests the class's instantiation.

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

        # Create the constructing fixtures.
        self.line_vortex_fixture = vortex_fixtures.make_line_vortex_fixture()
        self.origin_fixture = vortex_fixtures.make_origin_fixture()
        self.termination_fixture = vortex_fixtures.make_termination_fixture()
        self.strength_fixture = vortex_fixtures.make_strength_fixture()

    def tearDown(self):
        """This method is automatically called before each testing method to tear
        down the fixtures.

        :return: None
        """

        # Delete the constructing fixtures.
        del self.line_vortex_fixture
        del self.origin_fixture
        del self.termination_fixture
        del self.strength_fixture

    def test_class(self):
        """This method tests the class's instantiation.

        :return:
        """

        # Test that the object is of the right type.
        self.assertIsInstance(self.line_vortex_fixture, ps.aerodynamics.LineVortex)

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
