"""This module contains classes to test CoupledUnsteadyProblems."""

import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import (
    movement_fixtures,
    problem_fixtures,
)


class TestCoupledUnsteadyProblemInitialization(unittest.TestCase):
    """This is a class with functions to test CoupledUnsteadyProblem initialization."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all CoupledUnsteadyProblem tests."""
        cls.basic_coupled_problem = (
            problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        )
        cls.static_coupled_problem = (
            problem_fixtures.make_static_coupled_unsteady_problem_fixture()
        )

    def test_initialization_returns_correct_type(self):
        """Test that CoupledUnsteadyProblem initialization returns the correct type."""
        self.assertIsInstance(
            self.basic_coupled_problem,
            ps.problems.CoupledUnsteadyProblem,
        )

    def test_is_subclass_of_core(self):
        """Test that CoupledUnsteadyProblem is a subclass of CoreUnsteadyProblem."""
        self.assertIsInstance(
            self.basic_coupled_problem,
            ps._core.CoreUnsteadyProblem,
        )

    def test_movement_property_returns_core_movement(self):
        """Test that the movement property returns a CoreMovement instance."""
        self.assertIsInstance(
            self.basic_coupled_problem.movement,
            ps._core.CoreMovement,
        )

    def test_num_steps_matches_movement(self):
        """Test that num_steps matches the movement's num_steps."""
        movement = movement_fixtures.make_basic_movement_fixture()
        self.assertEqual(self.basic_coupled_problem.num_steps, movement.num_steps)

    def test_delta_time_matches_movement(self):
        """Test that delta_time matches the movement's delta_time."""
        movement = movement_fixtures.make_basic_movement_fixture()
        self.assertAlmostEqual(
            self.basic_coupled_problem.delta_time,
            movement.delta_time,
            places=10,
        )

    def test_initial_steady_problems_count(self):
        """Test that initialization creates exactly one SteadyProblem for step 0."""
        self.assertEqual(len(self.basic_coupled_problem.steady_problems), 1)

    def test_initial_steady_problem_is_steady_problem_type(self):
        """Test that the initial SteadyProblem is a SteadyProblem instance."""
        self.assertIsInstance(
            self.basic_coupled_problem.steady_problems[0],
            ps.problems.SteadyProblem,
        )

    def test_initial_steady_problem_has_airplanes(self):
        """Test that the initial SteadyProblem's Airplanes are populated."""
        initial_problem = self.basic_coupled_problem.steady_problems[0]
        self.assertGreater(len(initial_problem.airplanes), 0)
        for airplane in initial_problem.airplanes:
            self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_initial_steady_problem_has_operating_point(self):
        """Test that the initial SteadyProblem has an OperatingPoint."""
        initial_problem = self.basic_coupled_problem.steady_problems[0]
        self.assertIsInstance(
            initial_problem.operating_point,
            ps.operating_point.OperatingPoint,
        )

    def test_only_final_results_default_false(self):
        """Test that only_final_results defaults to False."""
        self.assertFalse(self.basic_coupled_problem.only_final_results)

    def test_only_final_results_true(self):
        """Test that only_final_results can be set to True."""
        movement = movement_fixtures.make_basic_movement_fixture()
        problem = ps.problems.CoupledUnsteadyProblem(
            movement=movement,
            only_final_results=True,
        )
        self.assertTrue(problem.only_final_results)

    def test_static_movement_first_averaging_step(self):
        """Test that static motion sets first_averaging_step to num_steps - 1."""
        self.assertEqual(
            self.static_coupled_problem.first_averaging_step,
            self.static_coupled_problem.num_steps - 1,
        )


class TestCoupledUnsteadyProblemParameterValidation(unittest.TestCase):
    """Tests for CoupledUnsteadyProblem parameter validation."""

    def test_movement_must_be_core_movement(self):
        """Test that movement parameter must be a CoreMovement or subclass."""
        with self.assertRaises(TypeError):
            ps.problems.CoupledUnsteadyProblem(movement="not_a_movement")

    def test_movement_rejects_none(self):
        """Test that movement parameter rejects None."""
        with self.assertRaises(TypeError):
            ps.problems.CoupledUnsteadyProblem(movement=None)

    def test_movement_rejects_invalid_types(self):
        """Test that movement parameter rejects various invalid types."""
        invalid_movements = [123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_movements:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.CoupledUnsteadyProblem(movement=invalid)


class TestCoupledUnsteadyProblemGetSteadyProblem(unittest.TestCase):
    """Tests for CoupledUnsteadyProblem.get_steady_problem."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.coupled_problem = (
            problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        )

    def test_get_step_zero_returns_steady_problem(self):
        """Test that get_steady_problem(0) returns a SteadyProblem."""
        result = self.coupled_problem.get_steady_problem(0)
        self.assertIsInstance(result, ps.problems.SteadyProblem)

    def test_get_step_zero_returns_initial_problem(self):
        """Test that get_steady_problem(0) returns the problem created at init."""
        result = self.coupled_problem.get_steady_problem(0)
        initial = self.coupled_problem.steady_problems[0]
        self.assertIs(result, initial)

    def test_negative_step_raises_value_error(self):
        """Test that negative step index raises ValueError."""
        with self.assertRaises(ValueError):
            self.coupled_problem.get_steady_problem(-1)

    def test_out_of_range_step_raises_value_error(self):
        """Test that step index beyond initialized problems raises ValueError."""
        with self.assertRaises(ValueError):
            self.coupled_problem.get_steady_problem(1)


class TestCoupledUnsteadyProblemInitializeNextProblem(unittest.TestCase):
    """Tests for CoupledUnsteadyProblem.initialize_next_problem."""

    def test_initialize_next_problem_appends_problem(self):
        """Test that initialize_next_problem adds one SteadyProblem."""
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        initial_count = len(coupled_problem.steady_problems)

        # Create a minimal solver to pass as the argument.
        solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
            coupled_problem
        )

        coupled_problem.initialize_next_problem(solver)
        self.assertEqual(len(coupled_problem.steady_problems), initial_count + 1)

    def test_initialize_next_problem_creates_valid_steady_problem(self):
        """Test that the appended SteadyProblem has the correct structure."""
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
            coupled_problem
        )

        coupled_problem.initialize_next_problem(solver)
        new_problem = coupled_problem.get_steady_problem(1)
        self.assertIsInstance(new_problem, ps.problems.SteadyProblem)
        self.assertGreater(len(new_problem.airplanes), 0)
        self.assertIsInstance(
            new_problem.operating_point,
            ps.operating_point.OperatingPoint,
        )

    def test_initialize_multiple_steps(self):
        """Test that calling initialize_next_problem repeatedly grows the list."""
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
            coupled_problem
        )

        num_steps_to_add = 3
        for _ in range(num_steps_to_add):
            coupled_problem.initialize_next_problem(solver)

        self.assertEqual(len(coupled_problem.steady_problems), 1 + num_steps_to_add)

    def test_get_steady_problem_after_initialize(self):
        """Test that get_steady_problem works for newly added steps."""
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
            coupled_problem
        )

        coupled_problem.initialize_next_problem(solver)
        result = coupled_problem.get_steady_problem(1)
        self.assertIsInstance(result, ps.problems.SteadyProblem)


class TestCoupledUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for CoupledUnsteadyProblem attribute immutability."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.coupled_problem = (
            problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        )

    def test_immutable_movement_property(self):
        """Test that movement property is read only."""
        with self.assertRaises(AttributeError):
            self.coupled_problem.movement = None

    def test_steady_problems_returns_tuple(self):
        """Test that steady_problems returns a tuple, not a mutable list."""
        result = self.coupled_problem.steady_problems
        self.assertIsInstance(result, tuple)

    def test_steady_problems_tuple_prevents_mutation(self):
        """Test that the tuple returned by steady_problems cannot be appended to."""
        result = self.coupled_problem.steady_problems
        with self.assertRaises(AttributeError):
            result.append(None)
