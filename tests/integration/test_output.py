"""This module is a testing case for the output module.

Note: Most of the tests in this case do not currently test against an expected
result. Instead, they test that the functions in output.py don't throw any errors.
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestOutput(unittest.TestCase):
    """This is a class with functions to test the output module."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests in this class.

        :return: None
        """

        # Set up the constructing fixtures.
        cls.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.unsteady_solver.run()

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
        """This method tests that the animate function does not throw any errors.

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
        """This method tests that the draw function does not throw any errors.

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


class TestPrintResults(unittest.TestCase):
    """Tests for the print_results() function."""

    @classmethod
    def setUpClass(cls):
        """Set up solvers for testing print_results().

        :return: None
        """
        cls.steady_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )
        cls.steady_solver.run()

        cls.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.unsteady_solver.run()

    def test_print_results_steady_solver_displays_reynolds_number(self):
        """Test that Reynolds number is displayed for steady solvers.

        :return: None
        """
        import io
        import sys

        # Capture stdout
        captured_output = io.StringIO()
        sys.stdout = captured_output

        try:
            ps.output.print_results(solver=self.steady_solver)
            output = captured_output.getvalue()

            # Verify Reynolds number is in the output
            self.assertIn("Reynolds Number:", output)
            # Verify scientific notation format (e.g., "1.23e+05" or "1.23e+06")
            self.assertRegex(output, r"Reynolds Number:\s+\d+\.\d{2}e[+-]\d{2}")
        finally:
            sys.stdout = sys.__stdout__

    def test_print_results_steady_solver_runs_without_error(self):
        """Test that print_results() runs without error for steady solver.

        :return: None
        """
        # This test ensures no exceptions are raised
        ps.output.print_results(solver=self.steady_solver)

    def test_print_results_unsteady_solver_runs_without_error(self):
        """Test that print_results() runs without error for unsteady solver.

        :return: None
        """
        # This test ensures no exceptions are raised
        ps.output.print_results(solver=self.unsteady_solver)
