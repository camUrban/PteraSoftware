"""This module contains classes to test SteadyRingVortexLatticeMethodSolvers."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestSteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test SteadyRingVortexLatticeMethodSolvers."""

    def test_initialization_accepts_steady_problem(self) -> None:
        """Test that initialization accepts a SteadyProblem."""
        solver = solver_fixtures.make_steady_ring_solver_fixture()
        self.assertIsInstance(
            solver,
            ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
        )

    def test_initialization_rejects_non_zero_body_rates(self) -> None:
        """Test that initialization raises when the operating point carries a non-zero
        body angular velocity, which this solver does not model."""
        rotating_problem = (
            problem_fixtures.make_with_body_rates_steady_problem_fixture()
        )
        with self.assertRaises(ValueError):
            ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
                rotating_problem
            )
