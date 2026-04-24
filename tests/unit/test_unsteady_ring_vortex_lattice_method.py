"""This module contains classes to test UnsteadyRingVortexLatticeMethodSolvers."""

import unittest
from unittest.mock import patch

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestUnsteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test UnsteadyRingVortexLatticeMethodSolvers."""

    def test_initialization_accepts_unsteady_problem(self):
        """Test that initialization accepts an UnsteadyProblem."""
        solver = solver_fixtures.make_unsteady_ring_solver_fixture()
        self.assertIsInstance(
            solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )

    def test_initialization_rejects_coupled_unsteady_problem(self):
        """Test that initialization on the base solver raises TypeError for a
        _CoupledUnsteadyProblem, while still allowing the coupled subclass to
        pass one through super().
        """
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        with self.assertRaises(TypeError):
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                coupled_problem
            )

    def test_initialization_rejects_non_problem_types(self):
        """Test that initialization raises TypeError for non-problem inputs."""
        invalid_inputs = [None, "not_a_problem", 123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
                        invalid
                    )


class TestUnsteadyRingVortexLatticeMethodSolverHookDefaults(unittest.TestCase):
    """Tests for the default implementations of the three solver extension hooks
    added to support coupled subclasses: _initialize_step_vortices,
    _reinitialize_step_arrays_hook, and _pre_shed_hook.
    """

    def setUp(self):
        """Set up a fresh solver for each hook-default test."""
        self.solver = solver_fixtures.make_unsteady_ring_solver_fixture()

    def test_initialize_step_vortices_initializes_all_on_step_zero(self):
        """Test that the default _initialize_step_vortices initializes bound
        vortices for all steps when called with step == 0.
        """
        with patch.object(
            type(self.solver), "_initialize_panel_vortices", autospec=True
        ) as mock_init:
            self.solver._initialize_step_vortices(0)
        mock_init.assert_called_once_with(self.solver)

    def test_initialize_step_vortices_is_noop_for_later_steps(self):
        """Test that the default _initialize_step_vortices is a no-op for any
        step greater than 0.
        """
        later_steps = [1, 2, self.solver.num_steps - 1]
        with patch.object(
            type(self.solver), "_initialize_panel_vortices", autospec=True
        ) as mock_init:
            for step in later_steps:
                self.solver._initialize_step_vortices(step)
        mock_init.assert_not_called()

    def test_reinitialize_step_arrays_hook_default_is_noop(self):
        """Test that the default _reinitialize_step_arrays_hook is a no-op."""
        self.assertIsNone(self.solver._reinitialize_step_arrays_hook())

    def test_pre_shed_hook_default_is_noop(self):
        """Test that the default _pre_shed_hook is a no-op for all steps."""
        for step in [0, 1, self.solver.num_steps - 1]:
            with self.subTest(step=step):
                self.assertIsNone(self.solver._pre_shed_hook(step))


if __name__ == "__main__":
    unittest.main()
