"""This module contains classes to test SteadyRingVortexLatticeMethodSolvers."""

import unittest
from unittest.mock import patch, sentinel

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestSteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test SteadyRingVortexLatticeMethodSolvers."""

    def test_initialization_accepts_steady_problem(self) -> None:
        """Test that initialization accepts a SteadyProblem."""
        solver = solver_fixtures.make_steady_ring_solver_fixture()
        self.assertIsInstance(
            solver,
            ps.SteadyRingVortexLatticeMethodSolver,
        )

    def test_initialization_rejects_non_zero_body_rates(self) -> None:
        """Test that initialization raises when the operating point carries a non-zero
        body angular velocity, which this solver does not model."""
        rotating_problem = (
            problem_fixtures.make_with_body_rates_steady_problem_fixture()
        )
        with self.assertRaises(ValueError):
            ps.SteadyRingVortexLatticeMethodSolver(steady_problem=rotating_problem)

    def test_ran_is_false_before_run(self) -> None:
        """Test that ran is False before run has been called."""
        solver = solver_fixtures.make_steady_ring_solver_fixture()
        self.assertFalse(solver.ran)

    def test_ran_is_read_only(self) -> None:
        """Test that the ran property cannot be reassigned."""
        solver = solver_fixtures.make_steady_ring_solver_fixture()
        with self.assertRaises(AttributeError):
            setattr(solver, "ran", True)


class TestSteadyRingVortexLatticeMethodSolverDeprecatedInternalMembers(
    unittest.TestCase
):
    """This class contains unit tests for the deprecated aliases of the
    SteadyRingVortexLatticeMethodSolver members that are no longer public."""

    def setUp(self) -> None:
        """Set up test fixtures for the deprecated alias tests."""
        self.solver = solver_fixtures.make_steady_ring_solver_fixture()

    def test_deprecated_properties_warn_and_return_internal_values(self) -> None:
        """Test that reading each deprecated property warns and returns the value of its
        internal slot."""
        for name in (
            "num_trailing_edge_panels",
            "vInf_GP1__E",
            "stackFreestreamWingInfluences__E",
            "panels",
            "stackUnitNormals_GP1",
            "panel_areas",
            "stackCpp_GP1_CgP1",
            "stackRbrv_GP1",
            "stackFbrv_GP1",
            "stackLbrv_GP1",
            "stackBbrv_GP1",
            "panel_is_trailing_edge",
            "panel_is_leading_edge",
            "panel_is_left_edge",
            "panel_is_right_edge",
            "stackSeedPoints_GP1_CgP1",
        ):
            with self.subTest(name=name):
                with self.assertWarns(DeprecationWarning):
                    value = getattr(self.solver, name)
                self.assertIs(value, getattr(self.solver, f"_{name}"))

    def test_calculate_solution_velocity_warns_and_delegates(self) -> None:
        """Test that calling calculate_solution_velocity warns and delegates to
        _calculate_solution_velocity."""
        stackP_GP1_CgP1 = [[0.0, 0.0, 0.0]]
        with patch.object(
            ps.SteadyRingVortexLatticeMethodSolver,
            "_calculate_solution_velocity",
            return_value=sentinel.velocities,
        ) as mock_calculate_solution_velocity:
            with self.assertWarns(DeprecationWarning):
                result = self.solver.calculate_solution_velocity(stackP_GP1_CgP1)
        mock_calculate_solution_velocity.assert_called_once_with(
            stackP_GP1_CgP1, bound_singularity_counts=None
        )
        self.assertIs(result, sentinel.velocities)
