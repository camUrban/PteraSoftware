"""This module contains classes to test the
FreeFlightUnsteadyRingVortexLatticeMethodSolver class."""

import unittest

import pterasoftware as ps
from tests.unit.test_free_flight_unsteady_problem import (
    _make_minimal_free_flight_movement,
)


class TestFreeFlightUnsteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """Tests for the FreeFlightUnsteadyRingVortexLatticeMethodSolver class."""

    def setUp(self) -> None:
        """Create a FreeFlightUnsteadyRingVortexLatticeMethodSolver fixture.

        :return: None
        """
        movement = _make_minimal_free_flight_movement()
        self.problem = ps.problems.FreeFlightUnsteadyProblem(movement=movement)
        self.solver = ps.free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver(
            free_flight_unsteady_problem=self.problem,
        )

    def test_initialization_rejects_non_free_flight_unsteady_problem(self) -> None:
        """FreeFlightUnsteadyRingVortexLatticeMethodSolver must raise TypeError when
        passed a problem that is not a FreeFlightUnsteadyProblem.

        :return: None
        """
        with self.assertRaises(TypeError):
            ps.free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver(
                free_flight_unsteady_problem=object(),  # type: ignore[arg-type]
            )

    def test_initialization_accepts_free_flight_unsteady_problem(self) -> None:
        """FreeFlightUnsteadyRingVortexLatticeMethodSolver must accept a valid
        FreeFlightUnsteadyProblem.

        :return: None
        """
        self.assertIsInstance(
            self.solver,
            ps.free_flight_unsteady_ring_vortex_lattice_method.FreeFlightUnsteadyRingVortexLatticeMethodSolver,
        )

    def test_unsteady_problem_is_stored(self) -> None:
        """The solver must expose the FreeFlightUnsteadyProblem passed at
        construction via its unsteady_problem attribute.

        :return: None
        """
        self.assertIs(self.solver.unsteady_problem, self.problem)
