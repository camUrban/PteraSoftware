"""This module contains classes to test CoreUnsteadyProblems."""

import math
import unittest

import numpy as np

import pterasoftware as ps


class TestCoreUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test CoreUnsteadyProblems."""

    def test_initialization_static(self):
        """Test CoreUnsteadyProblem initialization with static parameters."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertFalse(core_unsteady_problem.only_final_results)
        self.assertEqual(core_unsteady_problem.delta_time, 0.01)
        self.assertEqual(core_unsteady_problem.num_steps, 50)
        self.assertIsNone(core_unsteady_problem.max_wake_rows)

    def test_initialization_cyclic(self):
        """Test CoreUnsteadyProblem initialization with cyclic parameters."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=300,
            max_wake_rows=None,
            lcm_period=2.0,
        )
        self.assertFalse(core_unsteady_problem.only_final_results)
        self.assertEqual(core_unsteady_problem.delta_time, 0.01)
        self.assertEqual(core_unsteady_problem.num_steps, 300)

    def test_first_averaging_step_static(self):
        """Test first_averaging_step for static CoreUnsteadyProblem."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        # For static (lcm_period == 0), first_averaging_step == num_steps - 1.
        self.assertEqual(core_unsteady_problem.first_averaging_step, 49)

    def test_first_averaging_step_cyclic(self):
        """Test first_averaging_step for cyclic CoreUnsteadyProblem."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=300,
            max_wake_rows=None,
            lcm_period=2.0,
        )
        # Expected: max(0, floor(300 - 2.0 / 0.01)) = max(0, floor(100)) = 100.
        expected = max(0, math.floor(300 - (2.0 / 0.01)))
        self.assertEqual(core_unsteady_problem.first_averaging_step, expected)

    def test_first_averaging_step_cyclic_period_exceeds_duration(self):
        """Test first_averaging_step when lcm_period exceeds total duration."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=100,
            max_wake_rows=None,
            lcm_period=5.0,
        )
        # Expected: max(0, floor(100 - 5.0 / 0.01)) = max(0, -400) = 0.
        self.assertEqual(core_unsteady_problem.first_averaging_step, 0)

    def test_first_results_step_only_final_results_false(self):
        """Test first_results_step when only_final_results is False."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertEqual(core_unsteady_problem.first_results_step, 0)

    def test_first_results_step_only_final_results_true(self):
        """Test first_results_step when only_final_results is True."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=True,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        # When only_final_results is True, first_results_step ==
        # first_averaging_step.
        self.assertEqual(
            core_unsteady_problem.first_results_step,
            core_unsteady_problem.first_averaging_step,
        )

    def test_only_final_results_accepts_numpy_bool(self):
        """Test that only_final_results accepts numpy bool values."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=np.bool_(True),
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertTrue(core_unsteady_problem.only_final_results)
        self.assertIsInstance(core_unsteady_problem.only_final_results, bool)

    def test_initialization_of_load_lists(self):
        """Test that load lists are initialized as empty."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=None,
            lcm_period=0.0,
        )

        # All 12 load lists should be initialized as empty lists.
        load_list_names = [
            "finalForces_W",
            "finalForceCoefficients_W",
            "finalMoments_W_CgP1",
            "finalMomentCoefficients_W_CgP1",
            "finalMeanForces_W",
            "finalMeanForceCoefficients_W",
            "finalMeanMoments_W_CgP1",
            "finalMeanMomentCoefficients_W_CgP1",
            "finalRmsForces_W",
            "finalRmsForceCoefficients_W",
            "finalRmsMoments_W_CgP1",
            "finalRmsMomentCoefficients_W_CgP1",
        ]
        for name in load_list_names:
            with self.subTest(name=name):
                load_list = getattr(core_unsteady_problem, name)
                self.assertIsInstance(load_list, list)
                self.assertEqual(len(load_list), 0)

    def test_max_wake_rows_none(self):
        """Test that max_wake_rows is None when not set."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertIsNone(core_unsteady_problem.max_wake_rows)

    def test_max_wake_rows_positive_int(self):
        """Test that max_wake_rows stores a positive int correctly."""
        core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=25,
            lcm_period=0.0,
        )
        self.assertEqual(core_unsteady_problem.max_wake_rows, 25)


class TestCoreUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for CoreUnsteadyProblem attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.core_unsteady_problem = ps._core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )

    def test_immutable_only_final_results_property(self):
        """Test that only_final_results property is read only."""
        with self.assertRaises(AttributeError):
            self.core_unsteady_problem.only_final_results = True

    def test_immutable_num_steps_property(self):
        """Test that num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.core_unsteady_problem.num_steps = 100

    def test_immutable_delta_time_property(self):
        """Test that delta_time property is read only."""
        with self.assertRaises(AttributeError):
            self.core_unsteady_problem.delta_time = 0.1

    def test_immutable_first_averaging_step_property(self):
        """Test that first_averaging_step property is read only."""
        with self.assertRaises(AttributeError):
            self.core_unsteady_problem.first_averaging_step = 0

    def test_immutable_first_results_step_property(self):
        """Test that first_results_step property is read only."""
        with self.assertRaises(AttributeError):
            self.core_unsteady_problem.first_results_step = 0

    def test_immutable_max_wake_rows_property(self):
        """Test that max_wake_rows property is read only."""
        with self.assertRaises(AttributeError):
            self.core_unsteady_problem.max_wake_rows = 10

    def test_mutable_load_lists(self):
        """Test that load lists remain mutable for solver population."""
        self.core_unsteady_problem.finalForces_W.append(np.array([1.0, 2.0, 3.0]))
        self.assertEqual(len(self.core_unsteady_problem.finalForces_W), 1)

        # Clean up.
        self.core_unsteady_problem.finalForces_W.pop()


if __name__ == "__main__":
    unittest.main()
