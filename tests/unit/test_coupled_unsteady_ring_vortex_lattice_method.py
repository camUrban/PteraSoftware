"""This module contains classes to test CoupledUnsteadyRingVortexLatticeMethodSolvers."""

import unittest
from typing import Any

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

    def setUp(self) -> None:
        """Set up a fresh problem and solver for each test."""
        self.solver = solver_fixtures.make_coupled_unsteady_ring_solver_fixture()
        self.problem = self.solver.unsteady_problem

    def test_initialization_accepts_coupled_unsteady_problem(self) -> None:
        """Test that initialization accepts a _CoupledUnsteadyProblem."""
        self.assertIsInstance(self.solver, CoupledUnsteadyRingVortexLatticeMethodSolver)
        self.assertIsInstance(
            self.solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )
        self.assertIs(self.solver.unsteady_problem, self.problem)

    def test_initialization_rejects_plain_unsteady_problem(self) -> None:
        """Test that initialization raises TypeError for a plain UnsteadyProblem."""
        plain_unsteady_problem: Any = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )
        with self.assertRaises(TypeError):
            CoupledUnsteadyRingVortexLatticeMethodSolver(plain_unsteady_problem)

    def test_initialization_rejects_non_problem_types(self) -> None:
        """Test that initialization raises TypeError for non-problem inputs."""
        invalid_inputs: list[Any] = ["not_a_problem", 123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    CoupledUnsteadyRingVortexLatticeMethodSolver(invalid)

    def test_initialization_rejects_none(self) -> None:
        """Test that initialization raises TypeError for None."""
        none_problem: Any = None
        with self.assertRaises(TypeError):
            CoupledUnsteadyRingVortexLatticeMethodSolver(none_problem)

    def test_coupled_unsteady_problem_property_narrows_unsteady_problem(self) -> None:
        """Test that the _coupled_unsteady_problem property returns the same object as
        unsteady_problem, narrowed to _CoupledUnsteadyProblem.
        """
        self.assertIs(
            self.solver._coupled_unsteady_problem, self.solver.unsteady_problem
        )
        self.assertIsInstance(
            self.solver._coupled_unsteady_problem, ps.problems._CoupledUnsteadyProblem
        )

    def test_get_steady_problem_at_dispatches_through_coupled_unsteady_problem(
        self,
    ) -> None:
        """Test that both steady_problems and _get_steady_problem_at reflect the
        problem's growing _steady_problems list.

        The coupled problem seeds its backing list with step zero and appends each
        later step during the solve. The solver's steady_problems property and its
        _get_steady_problem_at dispatch both read that live list, so an append is
        visible immediately through either path.
        """
        self.assertEqual(len(self.solver.steady_problems), 1)

        next_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        assert isinstance(self.problem, ps.problems._CoupledUnsteadyProblem)
        self.problem._steady_problems.append(next_steady_problem)

        self.assertEqual(len(self.solver.steady_problems), 2)
        self.assertIs(self.solver._get_steady_problem_at(1), next_steady_problem)
        self.assertIs(self.solver.steady_problems[1], next_steady_problem)

    def test_inherits_empty_slots(self) -> None:
        """Test that the subclass declares __slots__ = () so it does not gain an
        instance __dict__ that would defeat the parent's __slots__.
        """
        self.assertEqual(CoupledUnsteadyRingVortexLatticeMethodSolver.__slots__, ())
        with self.assertRaises(AttributeError):
            setattr(self.solver, "not_a_real_attribute", 42)
