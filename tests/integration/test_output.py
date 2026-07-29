"""This module is a testing case for the output module.

Note: Most of the tests in this case do not currently test against an expected result.
Instead, they test that the functions in output.py don't throw any errors.
"""

import contextlib
import csv
import logging
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import numpy.testing as npt
import webp

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _output_rendering
from tests.integration.fixtures import solver_fixtures


def _count_saved_webp_frames(path: Path) -> int:
    """Counts the frames a saved WebP animation holds by reading its container.

    The WebPData is held in a local rather than passed inline, since it is what keeps
    the underlying buffer alive while the decoder reads it.

    :param path: The path of the saved WebP animation.
    :return: The number of frames the animation holds.
    """
    with open(path, "rb") as webp_file:
        webp_data = webp.WebPData.from_buffer(webp_file.read())
    decoder = webp.WebPAnimDecoder.new(webp_data)
    return int(decoder.anim_info.frame_count)


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
        """This method tests that the plot_results_versus_time method doesn't throw any
        errors.

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
        """This method tests that the draw function does not throw any errors when an
        image surface is defined.

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
        """This method tests that the draw function does not throw any errors when an
        image surface is defined and wake vortices are shown.

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
        """This method tests that the animate function does not throw any errors when an
        image surface is defined.

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
        """This method tests that the animate function does not throw any errors when an
        image surface is defined and wake vortices are shown.

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
        """This method tests that the draw function does not throw any errors for a free
        flight solver.

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
        """This method tests that the draw function does not throw any errors for a free
        flight solver when wake vortices are shown.

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
            "Orientation (of the First Airplane's Body Axes Relative to Earth Axes "
            "Using an Intrinsic zy'x\" Sequence):",
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

    def test_plot_results_versus_time_writes_the_loads_and_state_csvs(self) -> None:
        """Test that save_csv writes a loads CSV and a state CSV for a free flight
        solver, and that the state CSV matches the solver's OperatingPoints.

        :return: None
        """
        movement = self.free_flight_solver.unsteady_problem.movement
        assert isinstance(
            movement, ps.movements.free_flight_movement.FreeFlightMovement
        )
        operating_points = movement.operating_point_movement.operating_points

        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_path = Path(temporary_directory_name)
            ps.output.plot_results_versus_time(
                unsteady_solver=self.free_flight_solver,
                show=False,
                save_csv=True,
                directory=temporary_path,
            )

            # The CSV export writes no PNGs, so the two CSVs must be the only outputs.
            self.assertEqual(
                {path.name for path in temporary_path.iterdir()},
                {"simple_glider_loads.csv", "simple_glider_state.csv"},
            )

            with open(
                temporary_path / "simple_glider_loads.csv",
                newline="",
                encoding="utf-8",
            ) as csv_file:
                loads_rows = list(csv.reader(csv_file))
            with open(
                temporary_path / "simple_glider_state.csv",
                newline="",
                encoding="utf-8",
            ) as csv_file:
                state_rows = list(csv.reader(csv_file))

        # The loads begin at the solver's first results step while the state begins at
        # time step 0, which is why the two quantities are split across two files.
        self.assertEqual(
            len(loads_rows) - 1,
            self.free_flight_solver.num_steps
            - self.free_flight_solver.first_results_step,
        )
        self.assertEqual(len(state_rows) - 1, len(operating_points))

        # The state CSV holds time plus the fourteen state columns, headed by the text
        # the state figures carry.
        state_header_row = state_rows[0]
        self.assertEqual(len(state_header_row), 15)
        self.assertEqual(state_header_row[0], "Time (s)")
        self.assertEqual(
            state_header_row[1],
            "X Position of the First Airplane's CG in Earth Axes Relative to the "
            "Earth Origin (m)",
        )
        self.assertEqual(state_header_row[13], "Angle of Attack (deg)")
        self.assertEqual(state_header_row[14], "Sideslip Angle (deg)")

        # The position and angle of attack columns must match the solver's
        # OperatingPoints exactly.
        state_data_rows = state_rows[1:]
        npt.assert_array_equal(
            np.array([float(row[1]) for row in state_data_rows], dtype=float),
            np.array(
                [
                    this_operating_point.CgP1_E_Eo[0]
                    for this_operating_point in operating_points
                ],
                dtype=float,
            ),
        )
        npt.assert_array_equal(
            np.array([float(row[13]) for row in state_data_rows], dtype=float),
            np.array(
                [
                    this_operating_point.alpha
                    for this_operating_point in operating_points
                ],
                dtype=float,
            ),
        )


class TestOutputFileWriting(unittest.TestCase):
    """This is a class with functions to test that the output functions write their
    files to the destinations the caller asks for and reject the destinations they
    cannot honor."""

    unsteady_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests in this class.

        :return: None
        """
        cls.unsteady_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.unsteady_solver.run(show_progress=False)

    def setUp(self) -> None:
        """Create a temporary directory to hold this test's output files.

        :return: None
        """
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temporary_path = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        """Remove the temporary directory and any files the test wrote into it.

        :return: None
        """
        self.temporary_directory.cleanup()

    def test_draw_saves_under_the_default_name(self) -> None:
        """Test that draw's default path saves the image as draw.webp in the working
        directory.

        :return: None
        """
        with contextlib.chdir(self.temporary_path):
            ps.output.draw(
                solver=self.unsteady_solver,
                scalar_type=None,
                show_wake_vortices=False,
                show_streamlines=False,
                save=True,
                testing=True,
            )

        saved_path = self.temporary_path / "draw.webp"
        self.assertTrue(saved_path.is_file())
        self.assertGreater(saved_path.stat().st_size, 0)

    def test_draw_saves_to_a_custom_path(self) -> None:
        """Test that draw saves the image to a custom path in a subdirectory.

        :return: None
        """
        renders_directory = self.temporary_path / "renders"
        renders_directory.mkdir()
        saved_path = renders_directory / "static_drawing.webp"

        ps.output.draw(
            solver=self.unsteady_solver,
            scalar_type=None,
            show_wake_vortices=False,
            show_streamlines=False,
            save=True,
            path=saved_path,
            testing=True,
        )

        self.assertTrue(saved_path.is_file())
        self.assertGreater(saved_path.stat().st_size, 0)

    def test_draw_rejects_a_missing_directory(self) -> None:
        """Test that draw rejects a path whose directory does not exist and writes
        nothing.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.output.draw(
                solver=self.unsteady_solver,
                scalar_type=None,
                show_wake_vortices=False,
                show_streamlines=False,
                save=True,
                path=self.temporary_path / "missing" / "draw.webp",
                testing=True,
            )

        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_animate_saves_under_the_default_name_while_dropping_frames(self) -> None:
        """Test that animate's default path saves the animation as animate.webp and that
        the aliasing warning stays silent for a static geometry even when frames are
        dropped.

        No frame count is asserted here, since the encoder merges a frame that is
        identical to its predecessor into a longer-lasting one, and a static scene's
        frames can be identical. The flapping animation test asserts the exact count
        instead, because its frames are guaranteed to differ.

        :return: None
        """
        # Ask for a speed 1.9 times what the maximum frame rate can carry, which
        # resolves to saving every second frame.
        speed = (
            1.9 * _output_rendering._MAX_FRAME_RATE * self.unsteady_solver.delta_time
        )

        with self.assertNoLogs("pterasoftware.output", level=logging.WARNING):
            with contextlib.chdir(self.temporary_path):
                ps.output.animate(
                    unsteady_solver=self.unsteady_solver,
                    scalar_type=None,
                    show_wake_vortices=False,
                    window_size=(320, 240),
                    save=True,
                    speed=speed,
                    testing=True,
                )

        saved_path = self.temporary_path / "animate.webp"
        self.assertTrue(saved_path.is_file())
        self.assertGreater(saved_path.stat().st_size, 0)

    def test_animate_rejects_a_missing_directory(self) -> None:
        """Test that animate rejects a path whose directory does not exist and writes
        nothing.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.output.animate(
                unsteady_solver=self.unsteady_solver,
                scalar_type=None,
                show_wake_vortices=False,
                save=True,
                path=self.temporary_path / "missing" / "animate.webp",
                testing=True,
            )

        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_plot_results_versus_time_saves_under_the_default_names(self) -> None:
        """Test that the default directory and prefix save the four load figures under
        the names derived from the Airplane's name in the working directory.

        :return: None
        """
        with contextlib.chdir(self.temporary_path):
            ps.output.plot_results_versus_time(
                unsteady_solver=self.unsteady_solver,
                show=False,
                save=True,
                resolution_dpi=50.0,
            )

        self.assertEqual(
            {path.name for path in self.temporary_path.iterdir()},
            {
                "symmetric_unsteady_validation_airplane_forces.png",
                "symmetric_unsteady_validation_airplane_force_coefficients.png",
                "symmetric_unsteady_validation_airplane_moments.png",
                "symmetric_unsteady_validation_airplane_moment_coefficients.png",
            },
        )

    def test_plot_results_versus_time_saves_with_a_directory_and_prefix(self) -> None:
        """Test that a custom directory and prefix save the four load figures under
        prefixed names in that directory.

        :return: None
        """
        results_directory = self.temporary_path / "results"
        results_directory.mkdir()

        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver,
            show=False,
            save=True,
            directory=results_directory,
            prefix="run_1",
            resolution_dpi=50.0,
        )

        self.assertEqual(
            {path.name for path in results_directory.iterdir()},
            {
                "run_1_symmetric_unsteady_validation_airplane_forces.png",
                "run_1_symmetric_unsteady_validation_airplane_force_coefficients.png",
                "run_1_symmetric_unsteady_validation_airplane_moments.png",
                "run_1_symmetric_unsteady_validation_airplane_moment_coefficients.png",
            },
        )

        # Nothing may land outside the directory the caller named.
        self.assertEqual(
            [path.name for path in self.temporary_path.iterdir()], ["results"]
        )

    def test_plot_results_versus_time_rejects_a_missing_directory(self) -> None:
        """Test that a directory that does not exist is rejected and nothing is written.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.output.plot_results_versus_time(
                unsteady_solver=self.unsteady_solver,
                show=False,
                save=True,
                directory=self.temporary_path / "missing",
            )

        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_plot_results_versus_time_rejects_a_prefix_with_a_separator(self) -> None:
        """Test that a prefix carrying a path separator is rejected and nothing is
        written.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.output.plot_results_versus_time(
                unsteady_solver=self.unsteady_solver,
                show=False,
                save=True,
                directory=self.temporary_path,
                prefix="runs/1",
            )

        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_plot_results_versus_time_writes_nothing_without_the_save_flags(
        self,
    ) -> None:
        """Test that neither figures nor CSVs are written when save and save_csv are
        both False.

        :return: None
        """
        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver,
            show=False,
            directory=self.temporary_path,
        )

        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_plot_results_versus_time_writes_the_loads_csv(self) -> None:
        """Test that save_csv writes one loads CSV whose headers and values match the
        plotted series.

        :return: None
        """
        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_solver,
            show=False,
            save_csv=True,
            directory=self.temporary_path,
        )

        # The CSV export writes no PNGs, so the loads file must be the only output.
        self.assertEqual(
            [path.name for path in self.temporary_path.iterdir()],
            ["symmetric_unsteady_validation_airplane_loads.csv"],
        )

        with open(
            self.temporary_path / "symmetric_unsteady_validation_airplane_loads.csv",
            newline="",
            encoding="utf-8",
        ) as csv_file:
            rows = list(csv.reader(csv_file))
        header_row = rows[0]
        data_rows = rows[1:]

        # The headers pin the column order to the figure order: time, then the forces,
        # force coefficients, moments, and moment coefficients.
        self.assertEqual(
            header_row,
            [
                "Time (s)",
                "Induced Drag in Wind Axes (N)",
                "Side Force in Wind Axes (N)",
                "Lift in Wind Axes (N)",
                "Induced Drag Coefficient in Wind Axes",
                "Side Force Coefficient in Wind Axes",
                "Lift Coefficient in Wind Axes",
                "Rolling Moment in Wind Axes Relative to the CG (N*m)",
                "Pitching Moment in Wind Axes Relative to the CG (N*m)",
                "Yawing Moment in Wind Axes Relative to the CG (N*m)",
                "Rolling Moment Coefficient in Wind Axes Relative to the CG",
                "Pitching Moment Coefficient in Wind Axes Relative to the CG",
                "Yawing Moment Coefficient in Wind Axes Relative to the CG",
            ],
        )

        # The rows span the results-carrying time steps.
        first_results_step = self.unsteady_solver.first_results_step
        num_steps = self.unsteady_solver.num_steps
        delta_time = self.unsteady_solver.delta_time
        self.assertEqual(len(data_rows), num_steps - first_results_step)

        # Gather the expected values from the solver's own loads, with the wind axes x
        # and z force components negated so they read as induced drag and lift, matching
        # the figures.
        expected_induced_drags: list[float] = []
        expected_lifts: list[float] = []
        expected_rolling_moments: list[float] = []
        for steady_problem in self.unsteady_solver.steady_problems[first_results_step:]:
            first_airplane = steady_problem.airplanes[0]
            assert first_airplane.forces_W is not None
            assert first_airplane.moments_W_Cg is not None
            expected_induced_drags.append(-float(first_airplane.forces_W[0]))
            expected_lifts.append(-float(first_airplane.forces_W[2]))
            expected_rolling_moments.append(float(first_airplane.moments_W_Cg[0]))

        # The values must match exactly.
        npt.assert_array_equal(
            np.array([float(row[0]) for row in data_rows], dtype=float),
            np.linspace(
                delta_time * first_results_step,
                delta_time * (num_steps - 1),
                num_steps - first_results_step,
                endpoint=True,
            ),
        )
        npt.assert_array_equal(
            np.array([float(row[1]) for row in data_rows], dtype=float),
            np.array(expected_induced_drags, dtype=float),
        )
        npt.assert_array_equal(
            np.array([float(row[3]) for row in data_rows], dtype=float),
            np.array(expected_lifts, dtype=float),
        )
        npt.assert_array_equal(
            np.array([float(row[7]) for row in data_rows], dtype=float),
            np.array(expected_rolling_moments, dtype=float),
        )


class TestAnimatePlayback(unittest.TestCase):
    """This is a class with functions to test animate's playback speed parameter over a
    solver with flapping geometry."""

    variable_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests in this class.

        :return: None
        """
        cls.variable_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry()
        )
        cls.variable_solver.run(show_progress=False)

    def setUp(self) -> None:
        """Create a temporary directory to hold this test's output files.

        :return: None
        """
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temporary_path = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        """Remove the temporary directory and any files the test wrote into it.

        :return: None
        """
        self.temporary_directory.cleanup()

    def test_animate_rejects_a_speed_too_slow_to_animate(self) -> None:
        """Test that a speed that would save fewer than one frame per second of playback
        is rejected before anything is written.

        :return: None
        """
        with self.assertRaises(ValueError) as context:
            ps.output.animate(
                unsteady_solver=self.variable_solver,
                scalar_type=None,
                show_wake_vortices=False,
                save=True,
                path=self.temporary_path / "animate.webp",
                speed=0.5 * self.variable_solver.delta_time,
                testing=True,
            )

        self.assertIn("too slow", str(context.exception))
        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_animate_rejects_a_speed_too_fast_to_animate(self) -> None:
        """Test that a speed that would save fewer than two frames is rejected before
        anything is written.

        :return: None
        """
        speed = (
            1.5
            * (self.variable_solver.num_steps - 1)
            * _output_rendering._MAX_FRAME_RATE
            * self.variable_solver.delta_time
        )

        with self.assertRaises(ValueError) as context:
            ps.output.animate(
                unsteady_solver=self.variable_solver,
                scalar_type=None,
                show_wake_vortices=False,
                save=True,
                path=self.temporary_path / "animate.webp",
                speed=speed,
                testing=True,
            )

        self.assertIn("too fast", str(context.exception))
        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_animate_warns_and_drops_frames_when_undersampling_the_flapping(
        self,
    ) -> None:
        """Test that saving with a speed that undersamples the flapping logs the
        aliasing warning and that the saved WebP holds the reduced number of frames.

        :return: None
        """
        delta_time = self.variable_solver.delta_time
        num_steps = self.variable_solver.num_steps
        movement = self.variable_solver.unsteady_problem.movement

        # Find the smallest stride that undersamples the fastest motion, then a speed
        # that resolves to it. The 0.5 margin keeps the quotient that the stride is the
        # ceiling of from landing on the next whole number.
        stride = max(
            2,
            math.floor(
                movement.min_period
                / (_output_rendering._MIN_FRAMES_PER_PERIOD * delta_time)
            )
            + 1,
        )
        self.assertLessEqual(stride, num_steps - 1)
        speed = (stride - 0.5) * _output_rendering._MAX_FRAME_RATE * delta_time

        animations_directory = self.temporary_path / "animations"
        animations_directory.mkdir()
        saved_path = animations_directory / "flapping.webp"

        with self.assertLogs("pterasoftware.output", level=logging.WARNING) as log:
            ps.output.animate(
                unsteady_solver=self.variable_solver,
                scalar_type=None,
                show_wake_vortices=False,
                window_size=(320, 240),
                save=True,
                path=saved_path,
                speed=speed,
                testing=True,
            )
        self.assertIn("frames per cycle", log.output[0])

        self.assertTrue(saved_path.is_file())
        self.assertEqual(
            _count_saved_webp_frames(saved_path), (num_steps - 1) // stride + 1
        )


class TestFormationOutput(unittest.TestCase):
    """This is a class with functions to test plotting the results of a formation
    simulation whose second Airplane has a name carrying a path separator."""

    formation_solver: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests in this class.

        :return: None
        """
        cls.formation_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_formation_solver()
        )
        cls.formation_solver.run(show_progress=False)

    def setUp(self) -> None:
        """Create a temporary directory to hold this test's output files.

        :return: None
        """
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temporary_path = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        """Remove the temporary directory and any files the test wrote into it.

        :return: None
        """
        self.temporary_directory.cleanup()

    def test_plot_results_versus_time_rejects_the_separator_named_airplane(
        self,
    ) -> None:
        """Test that saving figures rejects the second Airplane's name before writing
        any of the first Airplane's figures.

        :return: None
        """
        with self.assertRaises(ValueError) as context:
            ps.output.plot_results_versus_time(
                unsteady_solver=self.formation_solver,
                show=False,
                save=True,
                directory=self.temporary_path,
            )
        self.assertIn("Trail/Airplane", str(context.exception))

        # Every Airplane's name is checked before any figure is drawn, so the well-named
        # lead Airplane's figures must not have been written.
        self.assertEqual(list(self.temporary_path.iterdir()), [])

    def test_plot_results_versus_time_rejects_the_separator_named_airplane_for_csvs(
        self,
    ) -> None:
        """Test that exporting CSVs rejects the second Airplane's name before writing
        the first Airplane's loads CSV.

        :return: None
        """
        with self.assertRaises(ValueError) as context:
            ps.output.plot_results_versus_time(
                unsteady_solver=self.formation_solver,
                show=False,
                save_csv=True,
                directory=self.temporary_path,
            )
        self.assertIn("Trail/Airplane", str(context.exception))

        self.assertEqual(list(self.temporary_path.iterdir()), [])
