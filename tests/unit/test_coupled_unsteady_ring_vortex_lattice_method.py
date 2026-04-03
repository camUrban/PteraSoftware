"""This module contains classes to test CoupledUnsteadyRingVortexLatticeMethodSolver."""

import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import (
    problem_fixtures,
    solver_fixtures,
)


class TestCoupledSolverInitialization(unittest.TestCase):
    """Tests for CoupledUnsteadyRingVortexLatticeMethodSolver initialization."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()

    def test_returns_correct_type(self):
        """Test that the solver is the expected type."""
        self.assertIsInstance(
            self.solver,
            ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver,
        )

    def test_inherits_from_base_solver(self):
        """Test that the coupled solver inherits from the base solver."""
        self.assertIsInstance(
            self.solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )

    def test_has_coupled_unsteady_problem(self):
        """Test that the solver has the coupled_unsteady_problem attribute."""
        self.assertIsInstance(
            self.solver.coupled_unsteady_problem,
            ps.problems.CoupledUnsteadyProblem,
        )

    def test_unsteady_problem_equals_coupled_unsteady_problem(self):
        """Test that unsteady_problem and coupled_unsteady_problem are the same."""
        self.assertIs(
            self.solver.unsteady_problem,
            self.solver.coupled_unsteady_problem,
        )

    def test_has_not_run(self):
        """Test that the solver has not run after initialization."""
        self.assertFalse(self.solver.ran)

    def test_has_empty_steady_problems_data_storage(self):
        """Test that steady_problems_data_storage starts empty."""
        self.assertEqual(len(self.solver.steady_problems_data_storage), 0)


class TestCoupledSolverParameterValidation(unittest.TestCase):
    """Tests for CoupledUnsteadyRingVortexLatticeMethodSolver parameter validation."""

    def test_rejects_regular_unsteady_problem(self):
        """Test that a regular UnsteadyProblem is rejected."""
        unsteady_problem = problem_fixtures.make_basic_unsteady_problem_fixture()
        with self.assertRaises(TypeError):
            ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
                unsteady_problem
            )

    def test_rejects_none(self):
        """Test that None is rejected."""
        with self.assertRaises(TypeError):
            ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
                None
            )

    def test_rejects_invalid_types(self):
        """Test that various invalid types are rejected."""
        invalid_inputs = ["string", 42, [1, 2, 3]]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
                        invalid
                    )


class TestCoupledSolverGetSteadyProblemAt(unittest.TestCase):
    """Tests for the _get_steady_problem_at override."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()

    def test_dispatches_to_coupled_problem(self):
        """Test that _get_steady_problem_at returns the same object as the coupled
        problem's get_steady_problem."""
        from_solver = self.solver._get_steady_problem_at(0)
        from_problem = self.solver.coupled_unsteady_problem.get_steady_problem(0)
        self.assertIs(from_solver, from_problem)

    def test_returns_steady_problem(self):
        """Test that _get_steady_problem_at returns a SteadyProblem."""
        result = self.solver._get_steady_problem_at(0)
        self.assertIsInstance(result, ps.problems.SteadyProblem)


class TestCoupledSolverRun(unittest.TestCase):
    """Tests for running the coupled solver."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures by running the solver.

        This runs the solver once because each run is expensive.
        """
        cls.solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()
        cls.solver.run(
            prescribed_wake=True,
            calculate_streamlines=False,
            show_progress=False,
        )

    def test_ran_is_true(self):
        """Test that the ran flag is True after running."""
        self.assertTrue(self.solver.ran)

    def test_steady_problems_populated(self):
        """Test that steady_problems is populated after running."""
        self.assertIsInstance(self.solver.steady_problems, tuple)
        self.assertEqual(len(self.solver.steady_problems), self.solver.num_steps)

    def test_steady_problems_are_steady_problem_instances(self):
        """Test that each element of steady_problems is a SteadyProblem."""
        for problem in self.solver.steady_problems:
            self.assertIsInstance(problem, ps.problems.SteadyProblem)

    def test_produces_nonzero_forces(self):
        """Test that the solver produces nonzero force coefficients."""
        for airplane in self.solver.current_airplanes:
            force_norm = np.linalg.norm(airplane.forceCoefficients_W)
            self.assertGreater(force_norm, 0.0)

    def test_coupled_problem_has_all_steps(self):
        """Test that the coupled problem contains the correct number of
        SteadyProblems."""
        self.assertEqual(
            len(self.solver.coupled_unsteady_problem.steady_problems),
            self.solver.num_steps,
        )


class TestCoupledSolverInitializeStepGeometry(unittest.TestCase):
    """Tests for the initialize_step_geometry method."""

    def test_rejects_negative_step(self):
        """Test that a negative step is rejected."""
        solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()
        with self.assertRaises(Exception):
            solver.initialize_step_geometry(-1)

    def test_rejects_step_beyond_range(self):
        """Test that a step beyond num_steps is rejected."""
        solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()
        with self.assertRaises(Exception):
            solver.initialize_step_geometry(solver.num_steps)
