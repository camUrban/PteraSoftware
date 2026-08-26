"""This module contains classes to test UnsteadyRingVortexLatticeMethodSolvers."""

import unittest
from typing import Any
from unittest.mock import MagicMock, patch

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestUnsteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test UnsteadyRingVortexLatticeMethodSolvers."""

    def test_initialization_accepts_unsteady_problem(self) -> None:
        """Test that initialization accepts an UnsteadyProblem."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        self.assertIsInstance(
            solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )

    def test_initialization_rejects_coupled_unsteady_problem(self) -> None:
        """Test that initialization on the base solver raises TypeError for a
        _CoupledUnsteadyProblem, while still allowing the coupled subclass to pass one
        through super()."""
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        with self.assertRaises(TypeError):
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                coupled_problem
            )

    def test_initialization_rejects_non_problem_types(self) -> None:
        """Test that initialization raises TypeError for non-problem inputs."""
        invalid_inputs: list[Any] = [
            None,
            "not_a_problem",
            123,
            [1, 2, 3],
            {"key": "value"},
        ]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        invalid
                    )

    def test_initialization_rejects_non_zero_body_rates(self) -> None:
        """Test that initialization raises when any per-step operating point carries a
        non-zero body angular velocity, which the base solver does not model."""
        rotating_problem = (
            problem_fixtures.make_with_body_rates_unsteady_problem_fixture()
        )
        with self.assertRaises(ValueError):
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                rotating_problem
            )

    def test_steady_problems_is_read_only(self) -> None:
        """Test that the steady_problems property cannot be reassigned, since it is a
        read-only view of the UnsteadyProblem's SteadyProblems.

        The property is defined on this base solver and inherited by every coupled
        subclass.
        """
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        with self.assertRaises(AttributeError):
            setattr(solver, "steady_problems", ())

    def test_prescribed_wake_raises_before_run(self) -> None:
        """Test that reading prescribed_wake before run raises a RuntimeError."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        with self.assertRaises(RuntimeError):
            _ = solver.prescribed_wake

    def test_prescribed_wake_is_read_only(self) -> None:
        """Test that the prescribed_wake property cannot be reassigned."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        with self.assertRaises(AttributeError):
            setattr(solver, "prescribed_wake", True)

    def test_force_method_raises_before_run(self) -> None:
        """Test that reading force_method before run raises a RuntimeError."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        with self.assertRaises(RuntimeError):
            _ = solver.force_method

    def test_force_method_is_read_only(self) -> None:
        """Test that the force_method property cannot be reassigned."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        with self.assertRaises(AttributeError):
            setattr(solver, "force_method", "joukowski")

    def test_ran_is_false_before_run(self) -> None:
        """Test that ran is False before run has been called."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        self.assertFalse(solver.ran)

    def test_ran_is_read_only(self) -> None:
        """Test that the ran property cannot be reassigned."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        with self.assertRaises(AttributeError):
            setattr(solver, "ran", True)

    def test_run_rejects_invalid_force_method_strings(self) -> None:
        """Test that run raises ValueError for invalid force_method strings."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        invalid_values = ["invalid", "JOUKOWSKI", "Katz", "both", ""]
        for invalid in invalid_values:
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    solver.run(force_method=invalid, show_progress=False)

    def test_run_rejects_non_string_force_methods(self) -> None:
        """Test that run raises TypeError for non-string force_method values."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        invalid_types: list[Any] = [
            123,
            1.0,
            None,
            True,
            ["joukowski"],
            {"method": "joukowski"},
        ]
        for invalid in invalid_types:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    solver.run(force_method=invalid, show_progress=False)


class TestUnsteadyRingVortexLatticeMethodSolverHookDefaults(unittest.TestCase):
    """Tests for the default implementations of the three solver extension hooks added
    to support coupled subclasses: _initialize_step_vortices,
    _reinitialize_step_arrays_hook, and _update_next_step_hook."""

    def setUp(self) -> None:
        """Set up a fresh solver for each hook-default test."""
        self.solver = solver_fixtures.make_unsteady_ring_solver_fixture()

    def test_initialize_step_vortices_initializes_all_on_step_zero(self) -> None:
        """Test that the default _initialize_step_vortices initializes bound vortices
        for all steps when called with step == 0."""
        with patch.object(
            type(self.solver), "_initialize_panel_vortices", autospec=True
        ) as mock_init:
            self.solver._initialize_step_vortices(0)
        mock_init.assert_called_once_with(self.solver)

    def test_initialize_step_vortices_is_noop_for_later_steps(self) -> None:
        """Test that the default _initialize_step_vortices is a no-op for any step
        greater than 0."""
        later_steps = [1, 2, self.solver.num_steps - 1]
        with patch.object(
            type(self.solver), "_initialize_panel_vortices", autospec=True
        ) as mock_init:
            for step in later_steps:
                self.solver._initialize_step_vortices(step)
        mock_init.assert_not_called()

    def test_reinitialize_step_arrays_hook_default_is_noop(self) -> None:
        """Test that the default _reinitialize_step_arrays_hook is a no-op."""
        self.solver._reinitialize_step_arrays_hook()

    def test_update_next_step_hook_default_is_noop(self) -> None:
        """Test that the default _update_next_step_hook is a no-op for all steps."""
        for step in [0, 1, self.solver.num_steps - 1]:
            with self.subTest(step=step):
                self.solver._update_next_step_hook(step)

    def test_models_body_rates_default_is_false(self) -> None:
        """Test that the base solver does not model body rates by default."""
        self.assertFalse(self.solver._models_body_rates)

    def test_current_omegas_default_is_zero(self) -> None:
        """Test that the default _currentOmegasRad_GP1__E returns a zero vector, so the
        base solver contributes no body-rotation velocity."""
        omegasRad_GP1__E = self.solver._currentOmegasRad_GP1__E()
        np.testing.assert_array_equal(omegasRad_GP1__E, np.zeros(3))

    def test_apply_body_rate_without_rotation_is_noop(self) -> None:
        """Test that _apply_body_rate returns the base velocities unchanged when the
        solver models no body rotation."""
        points = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        base_velocity = np.array([[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]])
        self.assertIs(
            self.solver._apply_body_rate(points, base_velocity), base_velocity
        )


class TestUnsteadyRingVortexLatticeMethodSolverFinalizeLoads(unittest.TestCase):
    """Tests for the trapezoid rule mean and RMS load averaging in _finalize_loads."""

    # The Airplane attributes that _finalize_loads reads at each time step in the
    # averaging window. Every one of them has to be present on the stand-in Airplanes,
    # since _finalize_loads collects all of them on every call.
    _LOAD_ATTRIBUTES = (
        "forces_W",
        "forceCoefficients_W",
        "moments_W_CgP1",
        "momentCoefficients_W_CgP1",
        "forces_G",
        "forceCoefficients_G",
        "moments_G_Cg",
        "momentCoefficients_G_Cg",
        "moments_W_Cg",
        "momentCoefficients_W_Cg",
    )

    def setUp(self) -> None:
        """Set up a fresh solver whose movement is variable, which is the case that
        averages over the final cycle."""
        self.solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        self.assertFalse(self.solver.unsteady_problem.movement.static)

    def _finalize_with_history(self, history: np.ndarray) -> None:
        """Run _finalize_loads over a prescribed single Airplane load history.

        Every load quantity is given the same history, so any of the UnsteadyProblem's
        mean and RMS lists reflects it.

        :param history: An ndarray of floats with shape (3, num_steps) holding the
            Airplane's load vector at each time step in the averaging window.
        :return: None
        """
        num_steps = history.shape[1]
        self.solver.num_steps = num_steps
        self.solver._first_averaging_step = 0

        steady_problems = [
            MagicMock(
                airplanes=[
                    MagicMock(
                        **{name: history[:, step] for name in self._LOAD_ATTRIBUTES}
                    )
                ]
            )
            for step in range(num_steps)
        ]

        with patch.object(
            type(self.solver),
            "_get_steady_problem_at",
            autospec=True,
            side_effect=lambda solver, step: steady_problems[step],
        ):
            self.solver._finalize_loads()

    def test_mean_loads_use_the_trapezoid_rule(self) -> None:
        """Test that the mean loads are the trapezoid rule integral of the load history
        divided by the number of intervals."""
        history = np.array(
            [
                [1.0, 2.0, 3.0, 4.0, 5.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [2.0, 2.0, 2.0, 2.0, 2.0],
            ],
            dtype=float,
        )

        self._finalize_with_history(history)

        # Over five evenly spaced samples the trapezoid rule halves the two endpoints,
        # so the ramp integrates to (0.5 * 1 + 2 + 3 + 4 + 0.5 * 5) / 4 = 3.0 and each
        # constant integrates to itself.
        np.testing.assert_allclose(
            self.solver.unsteady_problem.finalMeanForces_W[0],
            np.array([3.0, 0.0, 2.0], dtype=float),
        )

    def test_rms_loads_use_the_trapezoid_rule(self) -> None:
        """Test that the RMS loads are the square root of the trapezoid rule integral of
        the squared load history divided by the number of intervals."""
        history = np.array(
            [
                [1.0, 2.0, 3.0, 4.0, 5.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [2.0, 2.0, 2.0, 2.0, 2.0],
            ],
            dtype=float,
        )

        self._finalize_with_history(history)

        # The squared ramp integrates to (0.5 * 1 + 4 + 9 + 16 + 0.5 * 25) / 4 = 10.5,
        # so its RMS is sqrt(10.5). A constant's RMS is its own magnitude.
        np.testing.assert_allclose(
            self.solver.unsteady_problem.finalRmsForces_W[0],
            np.array([np.sqrt(10.5), 0.0, 2.0], dtype=float),
        )

    def test_mean_load_over_a_full_cycle_recovers_the_cycle_average(self) -> None:
        """Test that averaging a sinusoid sampled over exactly one cycle recovers the
        sinusoid's offset, since its oscillating part integrates to zero."""
        offset = 3.0
        amplitude = 7.0
        phases_rad = np.linspace(0.0, 2.0 * np.pi, 9, dtype=float)
        oscillation = offset + amplitude * np.sin(phases_rad)
        history = np.vstack([oscillation, oscillation, oscillation])

        self._finalize_with_history(history)

        np.testing.assert_allclose(
            self.solver.unsteady_problem.finalMeanForces_W[0],
            np.full(3, offset, dtype=float),
            atol=1e-12,
        )

    def test_single_sample_window_reports_that_sample(self) -> None:
        """Test that an averaging window holding one time step reports that step's loads
        as the mean and its magnitude as the RMS, rather than dividing by zero intervals
        and reporting NaN."""
        history = np.array([[-2.0], [0.0], [4.0]], dtype=float)

        with self.assertLogs(
            "pterasoftware.unsteady_ring_vortex_lattice_method", level="WARNING"
        ):
            self._finalize_with_history(history)

        np.testing.assert_allclose(
            self.solver.unsteady_problem.finalMeanForces_W[0],
            np.array([-2.0, 0.0, 4.0], dtype=float),
        )
        np.testing.assert_allclose(
            self.solver.unsteady_problem.finalRmsForces_W[0],
            np.array([2.0, 0.0, 4.0], dtype=float),
        )

    def test_single_sample_window_warns(self) -> None:
        """Test that a single time step averaging window warns that the loads are that
        step's values and says how to resolve the cycle."""
        history = np.array([[-2.0], [0.0], [4.0]], dtype=float)

        with self.assertLogs(
            "pterasoftware.unsteady_ring_vortex_lattice_method", level="WARNING"
        ) as caught:
            self._finalize_with_history(history)

        self.assertTrue(
            any("single time step" in message for message in caught.output),
            "The single sample warning should name the collapsed averaging window.",
        )
        self.assertTrue(
            any("delta_time" in message for message in caught.output),
            "The single sample warning should say how to resolve the cycle.",
        )

    def test_multiple_sample_window_does_not_warn(self) -> None:
        """Test that an averaging window holding more than one time step does not warn
        about a collapsed window."""
        history = np.array(
            [[1.0, 2.0], [0.0, 0.0], [2.0, 2.0]],
            dtype=float,
        )

        with self.assertNoLogs(
            "pterasoftware.unsteady_ring_vortex_lattice_method", level="WARNING"
        ):
            self._finalize_with_history(history)
