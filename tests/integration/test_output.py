"""This module is a testing case for the output module.

Note: Most of the tests in this case do not currently test the output against a
expected output. Instead, they test that the methods to create the output don't throw
any errors.

This module contains the following classes:
    TestOutput: This is a class with functions to test the output module.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""
import unittest

import matplotlib.pyplot as plt

import src.pterasoftware
from tests.integration.fixtures import solver_fixtures


class TestOutput(unittest.TestCase):
    """This is a class with functions to test the output module.

    This class contains the following public methods:
        setUp: This method is automatically called before each testing method to set
        up the fixtures.

        tearDown: This method is automatically called before each testing method to
        tear down the fixtures.

        test_plot_results_versus_time: This method tests the plot_results_versus_time
        method.

        test_animate_does_not_throw: This method tests that the animate method does
        not throw any errors.

        test_draw_does_not_throw: This method tests that the draw method does not
        throw any errors.

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
        self.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

    def tearDown(self):
        """This method is automatically called before each testing method to tear
        down the fixtures.

        :return: None
        """

        # Delete the constructing fixtures.
        del self.unsteady_solver

    def test_plot_results_versus_time(self):
        """This method tests the plot_results_versus_time method.

        :return: None
        """

        # Get the number of matplotlib figures before running the test.
        num_figs_before = plt.gcf().number

        # Call the plot_results_versus_time method on the solver fixture. The testing
        # flag is set to true, so the
        # figures will not be displayed.
        src.pterasoftware.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver, testing=True
        )

        # Get the number of matplotlib figures after running the test.
        num_figs_after = plt.gcf().number

        # Test that the method call produced four more figures.
        self.assertEqual(num_figs_before + 4, num_figs_after)

    def test_animate_does_not_throw(self):
        """This method tests that the animate method does not throw any errors.

        :return: None
        """

        # Call the animate function on the unsteady solver fixture.
        src.pterasoftware.output.animate(
            unsteady_solver=self.unsteady_solver,
            show_delta_pressures=False,
            show_wake_vortices=False,
            keep_file=False,
        )

    def test_draw_does_not_throw(self):
        """This method tests that the draw method does not throw any errors.

        :return: None
        """

        # Call the draw function on the unsteady solver fixture.
        src.pterasoftware.output.draw(
            solver=self.unsteady_solver,
            show_delta_pressures=False,
            show_wake_vortices=False,
            show_streamlines=False,
        )
