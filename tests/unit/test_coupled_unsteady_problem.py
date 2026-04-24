"""This module contains classes to test _CoupledUnsteadyProblems."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    core_movement_fixtures,
    geometry_fixtures,
    operating_point_fixtures,
    problem_fixtures,
)


class TestCoupledUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test _CoupledUnsteadyProblems."""

    def setUp(self):
        """Set up a fresh _CoupledUnsteadyProblem for each test."""
        self.movement = core_movement_fixtures.make_static_core_movement_fixture()
        self.initial_airplane = geometry_fixtures.make_first_airplane_fixture()
        self.initial_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        self.problem = ps.problems._CoupledUnsteadyProblem(
            movement=self.movement,
            initial_airplanes=[self.initial_airplane],
            initial_operating_point=self.initial_operating_point,
        )

    def test_initialization_valid_parameters(self):
        """Test _CoupledUnsteadyProblem initialization with valid parameters."""
        self.assertIsInstance(self.problem, ps.problems._CoupledUnsteadyProblem)
        self.assertIsInstance(self.problem, ps._core.CoreUnsteadyProblem)

    def test_step_zero_seeded_from_initial_inputs(self):
        """Test that _steady_problems is seeded with one SteadyProblem built from
        initial_airplanes and initial_operating_point.
        """
        self.assertEqual(len(self.problem._steady_problems), 1)
        seed = self.problem._steady_problems[0]
        self.assertIsInstance(seed, ps.problems.SteadyProblem)
        self.assertEqual(seed.airplanes, (self.initial_airplane,))
        self.assertIs(seed.operating_point, self.initial_operating_point)

    def test_only_final_results_forced_false(self):
        """Test that only_final_results is always False for coupled problems."""
        self.assertFalse(self.problem.only_final_results)

    def test_delta_time_from_movement(self):
        """Test that delta_time is taken from the movement."""
        self.assertEqual(self.problem.delta_time, self.movement.delta_time)

    def test_num_steps_from_movement(self):
        """Test that num_steps is taken from the movement."""
        self.assertEqual(self.problem.num_steps, self.movement.num_steps)

    def test_max_wake_rows_from_movement(self):
        """Test that max_wake_rows is taken from the movement."""
        self.assertEqual(self.problem.max_wake_rows, self.movement.max_wake_rows)

    def test_lcm_period_from_movement(self):
        """Test that lcm_period is taken from the movement by asserting the derived
        first_averaging_step matches the static-movement formula (num_steps - 1).
        """
        self.assertEqual(self.movement.lcm_period, 0.0)
        self.assertEqual(self.problem.first_averaging_step, self.movement.num_steps - 1)

    def test_movement_property_returns_core_movement(self):
        """Test that the movement property returns the provided CoreMovement."""
        self.assertIs(self.problem.movement, self.movement)

    def test_steady_problems_property_returns_tuple(self):
        """Test that the steady_problems property returns a tuple view of the
        _steady_problems backing list.
        """
        steady_problems = self.problem.steady_problems
        self.assertIsInstance(steady_problems, tuple)
        self.assertEqual(len(steady_problems), 1)
        self.assertIs(steady_problems[0], self.problem._steady_problems[0])

    def test_steady_problems_property_reflects_appends(self):
        """Test that appends to _steady_problems are reflected in the steady_problems
        tuple view on subsequent access.
        """
        self.assertEqual(len(self.problem.steady_problems), 1)

        next_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        self.problem._steady_problems.append(next_steady_problem)

        self.assertEqual(len(self.problem.steady_problems), 2)
        self.assertIs(self.problem.steady_problems[1], next_steady_problem)

    def test_get_steady_problem_returns_requested_step(self):
        """Test that get_steady_problem returns the SteadyProblem at the given
        step.
        """
        self.assertIs(
            self.problem.get_steady_problem(0),
            self.problem._steady_problems[0],
        )

        next_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        self.problem._steady_problems.append(next_steady_problem)
        self.assertIs(self.problem.get_steady_problem(1), next_steady_problem)

    def test_get_steady_problem_rejects_negative_step(self):
        """Test that get_steady_problem raises for negative step values."""
        with self.assertRaises(ValueError):
            self.problem.get_steady_problem(-1)

    def test_get_steady_problem_rejects_step_beyond_initialized(self):
        """Test that get_steady_problem raises for a step index that has not yet
        been initialized.
        """
        with self.assertRaises(ValueError):
            self.problem.get_steady_problem(1)

    def test_get_steady_problem_dynamic_bounds(self):
        """Test that the valid range of get_steady_problem grows as new
        SteadyProblems are appended to _steady_problems.
        """
        with self.assertRaises(ValueError):
            self.problem.get_steady_problem(1)

        next_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        self.problem._steady_problems.append(next_steady_problem)

        self.assertIs(self.problem.get_steady_problem(1), next_steady_problem)

    def test_initialize_next_problem_raises_not_implemented(self):
        """Test that initialize_next_problem raises NotImplementedError on the
        abstract middle class.
        """
        with self.assertRaises(NotImplementedError):
            self.problem.initialize_next_problem(None)


class TestCoupledUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for _CoupledUnsteadyProblem attribute immutability."""

    def setUp(self):
        """Set up a fresh _CoupledUnsteadyProblem for each immutability test."""
        self.problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()

    def test_immutable_movement_property(self):
        """Test that the movement property is read only."""
        new_movement = core_movement_fixtures.make_static_core_movement_fixture()
        with self.assertRaises(AttributeError):
            self.problem.movement = new_movement

    def test_immutable_steady_problems_property(self):
        """Test that the steady_problems property is read only."""
        with self.assertRaises(AttributeError):
            self.problem.steady_problems = ()

    def test_steady_problems_tuple_immutability(self):
        """Test that the steady_problems tuple cannot be mutated."""
        with self.assertRaises(AttributeError):
            self.problem.steady_problems.append(
                problem_fixtures.make_basic_steady_problem_fixture()
            )

    def test_private_steady_problems_list_is_mutable(self):
        """Test that _steady_problems remains mutable so subclass
        initialize_next_problem overrides can append the next step's SteadyProblem
        during the run loop.
        """
        next_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        self.problem._steady_problems.append(next_steady_problem)

        self.assertEqual(len(self.problem._steady_problems), 2)
        self.assertIs(self.problem._steady_problems[1], next_steady_problem)


if __name__ == "__main__":
    unittest.main()
