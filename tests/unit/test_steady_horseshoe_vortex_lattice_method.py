"""This module contains classes to test SteadyHorseshoeVortexLatticeMethodSolvers."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestSteadyHorseshoeVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test
    SteadyHorseshoeVortexLatticeMethodSolvers."""

    def test_initialization_accepts_steady_problem(self) -> None:
        """Test that initialization accepts a SteadyProblem."""
        solver = solver_fixtures.make_steady_horseshoe_solver_fixture()
        self.assertIsInstance(
            solver,
            ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
        )

    def test_initialization_rejects_non_zero_body_rates(self) -> None:
        """Test that initialization raises when the operating point carries a non-zero
        body angular velocity, which this solver does not model."""
        rotating_problem = (
            problem_fixtures.make_with_body_rates_steady_problem_fixture()
        )
        with self.assertRaises(ValueError):
            ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
                rotating_problem
            )

    def test_ran_is_false_before_run(self) -> None:
        """Test that ran is False before run has been called."""
        solver = solver_fixtures.make_steady_horseshoe_solver_fixture()
        self.assertFalse(solver.ran)

    def test_ran_is_read_only(self) -> None:
        """Test that the ran property cannot be reassigned."""
        solver = solver_fixtures.make_steady_horseshoe_solver_fixture()
        with self.assertRaises(AttributeError):
            setattr(solver, "ran", True)
