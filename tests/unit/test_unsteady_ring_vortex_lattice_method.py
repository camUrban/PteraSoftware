"""This module contains classes to test UnsteadyRingVortexLatticeMethodSolvers."""

import unittest
from typing import Any
from unittest.mock import patch

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
