# NOTE: I haven't yet started refactoring this module.
"""This module is a testing case for the output module.

Note: Most of the tests in this case do not currently test the output against an
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

import pterasoftware as ps
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

    def test_plot_results_versus_time_does_not_throw(self):
        """This method tests that the plot_results_versus_time method doesn't throw
        any errors.

        :return: None
        """

        # Call the plot_results_versus_time method on the solver fixture. The show
        # flag is set to False, so the figures will not be displayed.
        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver, show=False
        )

    def test_animate_does_not_throw(self):
        """This method tests that the animate method does not throw any errors.

        :return: None
        """

        # Call the animate function on the unsteady solver fixture. The testing flag
        # is true so the animation will start automatically after 1 second.
        ps.output.animate(
            unsteady_solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            save=False,
            testing=True,
        )

    def test_draw_does_not_throw(self):
        """This method tests that the draw method does not throw any errors.

        :return: None
        """

        # Call the draw function on the unsteady solver fixture. The testing flag is
        # set to true, so the plotter will close after 1 second.
        ps.output.draw(
            solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            show_streamlines=False,
            testing=True,
        )
