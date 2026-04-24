"""This module contains classes to test CoupledUnsteadyRingVortexLatticeMethodSolvers."""

import unittest

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware._coupled_unsteady_ring_vortex_lattice_method import (
    CoupledUnsteadyRingVortexLatticeMethodSolver,
)
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestCoupledUnsteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test
    CoupledUnsteadyRingVortexLatticeMethodSolvers.
    """

    def setUp(self):
        """Set up a fresh problem and solver for each test."""
        self.solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()
        self.problem = self.solver.unsteady_problem

    def test_initialization_accepts_coupled_unsteady_problem(self):
        """Test that initialization accepts a _CoupledUnsteadyProblem."""
        self.assertIsInstance(self.solver, CoupledUnsteadyRingVortexLatticeMethodSolver)
        self.assertIsInstance(
            self.solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )
        self.assertIs(self.solver.unsteady_problem, self.problem)

    def test_initialization_rejects_plain_unsteady_problem(self):
        """Test that initialization raises TypeError for a plain UnsteadyProblem."""
        plain_unsteady_problem = problem_fixtures.make_basic_unsteady_problem_fixture()
        with self.assertRaises(TypeError):
            CoupledUnsteadyRingVortexLatticeMethodSolver(plain_unsteady_problem)

    def test_initialization_rejects_non_problem_types(self):
        """Test that initialization raises TypeError for non-problem inputs."""
        invalid_inputs = ["not_a_problem", 123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    CoupledUnsteadyRingVortexLatticeMethodSolver(invalid)

    def test_initialization_rejects_none(self):
        """Test that initialization raises TypeError for None."""
        with self.assertRaises(TypeError):
            CoupledUnsteadyRingVortexLatticeMethodSolver(None)

    def test_coupled_problem_property_narrows_unsteady_problem(self):
        """Test that the _coupled_problem property returns the same object as
        unsteady_problem, narrowed to _CoupledUnsteadyProblem.
        """
        self.assertIs(self.solver._coupled_problem, self.solver.unsteady_problem)
        self.assertIsInstance(
            self.solver._coupled_problem, ps.problems._CoupledUnsteadyProblem
        )

    def test_get_steady_problem_at_dispatches_through_coupled_problem(self):
        """Test that _get_steady_problem_at dispatches through
        _CoupledUnsteadyProblem.get_steady_problem rather than the cached tuple.

        The parent solver captures steady_problems as a tuple snapshot at
        construction time. The coupled subclass must dispatch through the
        accessor so appends to the problem's _steady_problems backing list are
        visible at later steps.
        """
        self.assertEqual(len(self.solver.steady_problems), 1)

        next_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        self.problem._steady_problems.append(next_steady_problem)

        self.assertEqual(len(self.solver.steady_problems), 1)
        self.assertIs(self.solver._get_steady_problem_at(1), next_steady_problem)

    def test_inherits_empty_slots(self):
        """Test that the subclass declares __slots__ = () so it does not gain an
        instance __dict__ that would defeat the parent's __slots__.
        """
        self.assertEqual(CoupledUnsteadyRingVortexLatticeMethodSolver.__slots__, ())
        with self.assertRaises(AttributeError):
            self.solver.not_a_real_attribute = 42


if __name__ == "__main__":
    unittest.main()
