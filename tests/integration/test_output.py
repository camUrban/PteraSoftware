"""This module is a testing case for the output module.

Note: Most of the tests in this case do not currently test against an expected
result. Instead, they test that the functions in output.py don't throw any errors.
"""

import logging
import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestOutput(unittest.TestCase):
    """This is a class with functions to test the output module."""

    unsteady_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests in this class.

        :return: None
        """

        # Set up the constructing fixtures.
        cls.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.unsteady_solver.run(show_progress=False)

    def test_plot_results_versus_time_does_not_throw(self) -> None:
        """This method tests that the plot_results_versus_time method doesn't throw
        any errors.

        :return: None
        """

        # Call the plot_results_versus_time method on the solver fixture. The show flag
        # is set to False, so the figures will not be displayed.
        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver, show=False
        )

    def test_animate_does_not_throw(self) -> None:
        """This method tests that the animate function does not throw any errors.

        :return: None
        """

        # Call the animate function on the unsteady solver fixture. The testing flag is
        # true so the animation will start automatically after 1 second.
        ps.output.animate(
            unsteady_solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            save=False,
            testing=True,
        )

    def test_draw_does_not_throw(self) -> None:
        """This method tests that the draw function does not throw any errors.

        :return: None
        """

        # Call the draw function on the unsteady solver fixture. The testing flag is set
        # to true, so the plotter will close after 1 second.
        ps.output.draw(
            solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            show_streamlines=False,
            testing=True,
        )


class TestLogResults(unittest.TestCase):
    """Tests the log_results() function."""

    steady_solver: (
        ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
    )
    unsteady_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up solvers for testing log_results().

        :return: None
        """
        cls.steady_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )
        cls.steady_solver.run()

        cls.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.unsteady_solver.run(show_progress=False)

    def test_log_results_steady_solver_displays_reynolds_number(self) -> None:
        """Test that Reynolds number is logged for steady solvers.

        :return: None
        """
        # Capture log output using assertLogs context manager.
        with self.assertLogs("pterasoftware.output", level=logging.INFO) as log:
            ps.output.log_results(solver=self.steady_solver)

        output = "\n".join(log.output)

        # Verify Reynolds number is in the output.
        self.assertIn("Reynolds Number:", output)
        # Verify scientific notation format with three significant figures (e.g.,
        # "1.23E+05" or "1.20E+06"; the alternate form keeps trailing zeros).
        self.assertRegex(output, r"Reynolds Number:\s+\d\.\d{2}E[+-]\d{2}")

    def test_log_results_steady_solver_runs_without_error(self) -> None:
        """Test that log_results() runs without error for steady solver.

        :return: None
        """
        # This test ensures no exceptions are raised.
        ps.output.log_results(solver=self.steady_solver)

    def test_log_results_unsteady_solver_runs_without_error(self) -> None:
        """Test that log_results() runs without error for unsteady solver.

        :return: None
        """
        # This test ensures no exceptions are raised.
        ps.output.log_results(solver=self.unsteady_solver)


class TestOutputSurfaceEffect(unittest.TestCase):
    """This is a class with functions to test the output module's surface effect
    visualization, including reflected geometry and the image surface plane."""

    unsteady_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests in this class.

        :return: None
        """
        cls.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_surface_effect_solver()
        )
        cls.unsteady_solver.run(show_progress=False)

    def test_draw_with_surface_effect_does_not_throw(self) -> None:
        """This method tests that the draw function does not throw any errors when
        an image surface is defined.

        :return: None
        """
        ps.output.draw(
            solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            show_streamlines=False,
            testing=True,
        )

    def test_draw_with_surface_effect_and_wake_vortices_does_not_throw(self) -> None:
        """This method tests that the draw function does not throw any errors when
        an image surface is defined and wake vortices are shown.

        :return: None
        """
        ps.output.draw(
            solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=True,
            show_streamlines=False,
            testing=True,
        )

    def test_animate_with_surface_effect_does_not_throw(self) -> None:
        """This method tests that the animate function does not throw any errors
        when an image surface is defined.

        :return: None
        """
        ps.output.animate(
            unsteady_solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            save=False,
            testing=True,
        )

    def test_animate_with_surface_effect_and_wake_vortices_does_not_throw(self) -> None:
        """This method tests that the animate function does not throw any errors
        when an image surface is defined and wake vortices are shown.

        :return: None
        """
        ps.output.animate(
            unsteady_solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=True,
            save=False,
            testing=True,
        )


class TestFreeFlightOutput(unittest.TestCase):
    """This is a class with functions to test the output module's free flight
    visualization, which renders each time step's geometry in Earth axes so the body
    flies through the scene."""

    free_flight_solver: (
        ps.free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests in this class.

        :return: None
        """
        cls.free_flight_solver = solver_fixtures.make_simple_glider_free_flight_solver()
        cls.free_flight_solver.run(show_progress=False)

    def test_draw_does_not_throw(self) -> None:
        """This method tests that the draw function does not throw any errors for a
        free flight solver.

        :return: None
        """
        ps.output.draw(
            solver=self.free_flight_solver,
            scalar_type=None,
            show_wake_vortices=False,
            show_streamlines=False,
            testing=True,
        )

    def test_draw_with_wake_vortices_does_not_throw(self) -> None:
        """This method tests that the draw function does not throw any errors for a
        free flight solver when wake vortices are shown.

        :return: None
        """
        ps.output.draw(
            solver=self.free_flight_solver,
            scalar_type=None,
            show_wake_vortices=True,
            show_streamlines=False,
            testing=True,
        )

    def test_animate_does_not_throw(self) -> None:
        """This method tests that the animate function does not throw any errors for a
        free flight solver.

        :return: None
        """
        ps.output.animate(
            unsteady_solver=self.free_flight_solver,
            scalar_type=None,
            show_wake_vortices=False,
            save=False,
            testing=True,
        )

    def test_animate_with_wake_vortices_does_not_throw(self) -> None:
        """This method tests that the animate function does not throw any errors for a
        free flight solver when wake vortices are shown.

        :return: None
        """
        ps.output.animate(
            unsteady_solver=self.free_flight_solver,
            scalar_type=None,
            show_wake_vortices=True,
            save=False,
            testing=True,
        )

    def test_plot_results_versus_time_does_not_throw(self) -> None:
        """This method tests that the plot_results_versus_time function does not throw
        any errors for a free flight solver, which exercises the state-history plots.

        :return: None
        """
        ps.output.plot_results_versus_time(
            unsteady_solver=self.free_flight_solver, show=False
        )

    def test_log_results_logs_state_history(self) -> None:
        """This method tests that log_results logs the first Airplane's initial and
        final six-degree-of-freedom state for a free flight solver.

        :return: None
        """
        with self.assertLogs("pterasoftware.output", level=logging.INFO) as log:
            ps.output.log_results(solver=self.free_flight_solver)

        output = "\n".join(log.output)

        self.assertIn("The First Airplane's Free Flight State History:", output)
        self.assertIn("Initial State", output)
        self.assertIn("Final State", output)
        self.assertIn(
            "Position (of the First Airplane's CG, in Earth Axes, Relative to the "
            "Earth Origin):",
            output,
        )
        self.assertIn(
            "Orientation (of the First Airplane's Body Axes Relative to Earth Axes, "
            "Intrinsic zy'x\" Sequence):",
            output,
        )
        self.assertIn(
            "Velocity (of the First Airplane's CG, in Earth Axes, Observed from the "
            "Earth Frame):",
            output,
        )
        self.assertIn(
            "Angular Velocity (in the First Airplane's Body Axes, Observed from the "
            "Earth Frame):",
            output,
        )
        self.assertIn("Angle of Attack (alpha):", output)
        self.assertIn("Sideslip Angle (beta):", output)

        # Each vector state quantity is broken into one row per component, labeled with
        # its variable-convention name.
        self.assertIn("cgP1X_E_Eo:", output)
        self.assertIn("cgP1Y_E_Eo:", output)
        self.assertIn("cgP1Z_E_Eo:", output)
        self.assertIn("angleX_E_to_BP1_izyx:", output)
        self.assertIn("angleY_E_to_BP1_izyx:", output)
        self.assertIn("angleZ_E_to_BP1_izyx:", output)
        self.assertIn("vCgP1X_E__E:", output)
        self.assertIn("vCgP1Y_E__E:", output)
        self.assertIn("vCgP1Z_E__E:", output)
        self.assertIn("omegaX_BP1__E:", output)
        self.assertIn("omegaY_BP1__E:", output)
        self.assertIn("omegaZ_BP1__E:", output)
