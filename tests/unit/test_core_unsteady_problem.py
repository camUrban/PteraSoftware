"""This module contains classes to test CoreUnsteadyProblems."""

import math
import unittest

import numpy as np
import numpy.testing as npt

from pterasoftware import _core

# The names of CoreUnsteadyProblem's load list properties: the eight final lists, then
# the eight final mean lists, each in the order _record_final_loads takes its stacks.
LOAD_LIST_NAMES = (
    "finalForces_W",
    "finalForceCoefficients_W",
    "finalForces_G",
    "finalForceCoefficients_G",
    "finalMoments_G_Cg",
    "finalMomentCoefficients_G_Cg",
    "finalMoments_W_Cg",
    "finalMomentCoefficients_W_Cg",
    "finalMeanForces_W",
    "finalMeanForceCoefficients_W",
    "finalMeanForces_G",
    "finalMeanForceCoefficients_G",
    "finalMeanMoments_G_Cg",
    "finalMeanMomentCoefficients_G_Cg",
    "finalMeanMoments_W_Cg",
    "finalMeanMomentCoefficients_W_Cg",
)


class TestCoreUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test CoreUnsteadyProblems."""

    def test_initialization_static(self) -> None:
        """Test CoreUnsteadyProblem initialization with static parameters."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
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

    def test_initialization_cyclic(self) -> None:
        """Test CoreUnsteadyProblem initialization with cyclic parameters."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=300,
            max_wake_rows=None,
            lcm_period=2.0,
        )
        self.assertFalse(core_unsteady_problem.only_final_results)
        self.assertEqual(core_unsteady_problem.delta_time, 0.01)
        self.assertEqual(core_unsteady_problem.num_steps, 300)

    def test_first_averaging_step_static(self) -> None:
        """Test first_averaging_step for static CoreUnsteadyProblem."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        # For static (lcm_period == 0), first_averaging_step == num_steps - 1.
        self.assertEqual(core_unsteady_problem.first_averaging_step, 49)

    def test_first_averaging_step_cyclic(self) -> None:
        """Test first_averaging_step for cyclic CoreUnsteadyProblem."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=300,
            max_wake_rows=None,
            lcm_period=2.0,
        )
        # Expected: max(0, floor(300 - 2.0 / 0.01)) = max(0, floor(100)) = 100.
        expected = max(0, math.floor(300 - (2.0 / 0.01)))
        self.assertEqual(core_unsteady_problem.first_averaging_step, expected)

    def test_first_averaging_step_cyclic_period_exceeds_duration(self) -> None:
        """Test first_averaging_step when lcm_period exceeds total duration."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=100,
            max_wake_rows=None,
            lcm_period=5.0,
        )
        # Expected: max(0, floor(100 - 5.0 / 0.01)) = max(0, -400) = 0.
        self.assertEqual(core_unsteady_problem.first_averaging_step, 0)

    def test_first_results_step_only_final_results_false(self) -> None:
        """Test first_results_step when only_final_results is False."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertEqual(core_unsteady_problem.first_results_step, 0)

    def test_first_results_step_only_final_results_true(self) -> None:
        """Test first_results_step when only_final_results is True."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=True,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        # When only_final_results is True, first_results_step == first_averaging_step.
        self.assertEqual(
            core_unsteady_problem.first_results_step,
            core_unsteady_problem.first_averaging_step,
        )

    def test_only_final_results_accepts_numpy_bool(self) -> None:
        """Test that only_final_results accepts numpy bool values."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=np.bool(True),
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertTrue(core_unsteady_problem.only_final_results)
        self.assertIsInstance(core_unsteady_problem.only_final_results, bool)

    def test_initialization_of_load_lists(self) -> None:
        """Test that load lists are initialized as empty."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=None,
            lcm_period=0.0,
        )

        # All 16 load lists should be initialized as empty lists.
        for name in LOAD_LIST_NAMES:
            with self.subTest(name=name):
                load_list = getattr(core_unsteady_problem, name)
                self.assertIsInstance(load_list, list)
                self.assertEqual(len(load_list), 0)

    def test_max_wake_rows_none(self) -> None:
        """Test that max_wake_rows is None when not set."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=None,
            lcm_period=0.0,
        )
        self.assertIsNone(core_unsteady_problem.max_wake_rows)

    def test_max_wake_rows_positive_int(self) -> None:
        """Test that max_wake_rows stores a positive int correctly."""
        core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=10,
            max_wake_rows=25,
            lcm_period=0.0,
        )
        self.assertEqual(core_unsteady_problem.max_wake_rows, 25)


class TestCoreUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for CoreUnsteadyProblem attribute immutability."""

    def setUp(self) -> None:
        """Set up test fixtures for immutability tests."""
        self.core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )

    def test_immutable_only_final_results_property(self) -> None:
        """Test that only_final_results property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.core_unsteady_problem, "only_final_results", True)

    def test_immutable_num_steps_property(self) -> None:
        """Test that num_steps property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.core_unsteady_problem, "num_steps", 100)

    def test_immutable_delta_time_property(self) -> None:
        """Test that delta_time property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.core_unsteady_problem, "delta_time", 0.1)

    def test_immutable_first_averaging_step_property(self) -> None:
        """Test that first_averaging_step property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.core_unsteady_problem, "first_averaging_step", 0)

    def test_immutable_first_results_step_property(self) -> None:
        """Test that first_results_step property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.core_unsteady_problem, "first_results_step", 0)

    def test_immutable_max_wake_rows_property(self) -> None:
        """Test that max_wake_rows property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.core_unsteady_problem, "max_wake_rows", 10)

    def test_immutable_load_list_properties(self) -> None:
        """Test that the load list properties are read only."""
        for name in LOAD_LIST_NAMES:
            with self.subTest(name=name):
                with self.assertRaises(AttributeError):
                    setattr(self.core_unsteady_problem, name, [])

    def test_load_lists_cannot_be_mutated(self) -> None:
        """Test that recorded loads cannot be changed through the properties."""
        stackLoads = np.arange(48, dtype=float).reshape(8, 2, 3)
        self.core_unsteady_problem._record_final_loads(
            cycle_averaged=False,
            stackForces_W=stackLoads[0],
            stackForceCoefficients_W=stackLoads[1],
            stackForces_G=stackLoads[2],
            stackForceCoefficients_G=stackLoads[3],
            stackMoments_G_Cg=stackLoads[4],
            stackMomentCoefficients_G_Cg=stackLoads[5],
            stackMoments_W_Cg=stackLoads[6],
            stackMomentCoefficients_W_Cg=stackLoads[7],
        )

        # Each read returns a new list, so appending to one changes nothing.
        self.core_unsteady_problem.finalForces_W.append(np.zeros(3, dtype=float))
        self.assertEqual(len(self.core_unsteady_problem.finalForces_W), 2)

        # The recorded ndarrays are read only.
        with self.assertRaises(ValueError):
            self.core_unsteady_problem.finalForces_W[0][0] = 999.0

        # The recorded ndarrays share no memory with the ndarrays passed in.
        stackLoads[0, 0, 0] = 999.0
        self.assertEqual(self.core_unsteady_problem.finalForces_W[0][0], 0.0)


class TestCoreUnsteadyProblemRecordFinalLoads(unittest.TestCase):
    """Tests for CoreUnsteadyProblem's _record_final_loads method."""

    def setUp(self) -> None:
        """Set up test fixtures for the _record_final_loads tests."""
        self.core_unsteady_problem = _core.CoreUnsteadyProblem(
            only_final_results=False,
            delta_time=0.01,
            num_steps=50,
            max_wake_rows=None,
            lcm_period=0.0,
        )

        # Make one distinct (2, 3) stack per load quantity, as if for two Airplanes.
        self.stackLoads = np.arange(48, dtype=float).reshape(8, 2, 3)

    def record(self, cycle_averaged: bool) -> None:
        """Records the stacks to the CoreUnsteadyProblem.

        :param cycle_averaged: The value to pass as cycle_averaged.
        :return: None
        """
        self.core_unsteady_problem._record_final_loads(
            cycle_averaged=cycle_averaged,
            stackForces_W=self.stackLoads[0],
            stackForceCoefficients_W=self.stackLoads[1],
            stackForces_G=self.stackLoads[2],
            stackForceCoefficients_G=self.stackLoads[3],
            stackMoments_G_Cg=self.stackLoads[4],
            stackMomentCoefficients_G_Cg=self.stackLoads[5],
            stackMoments_W_Cg=self.stackLoads[6],
            stackMomentCoefficients_W_Cg=self.stackLoads[7],
        )

    def test_records_final_lists(self) -> None:
        """Test that cycle_averaged=False records only to the final lists."""
        self.record(cycle_averaged=False)

        for stack_id, name in enumerate(LOAD_LIST_NAMES[:8]):
            with self.subTest(name=name):
                load_list = getattr(self.core_unsteady_problem, name)
                self.assertEqual(len(load_list), 2)
                npt.assert_array_equal(np.array(load_list), self.stackLoads[stack_id])
        for name in LOAD_LIST_NAMES[8:]:
            with self.subTest(name=name):
                self.assertEqual(getattr(self.core_unsteady_problem, name), [])

    def test_records_final_mean_lists(self) -> None:
        """Test that cycle_averaged=True records only to the final mean lists."""
        self.record(cycle_averaged=True)

        for stack_id, name in enumerate(LOAD_LIST_NAMES[8:]):
            with self.subTest(name=name):
                load_list = getattr(self.core_unsteady_problem, name)
                self.assertEqual(len(load_list), 2)
                npt.assert_array_equal(np.array(load_list), self.stackLoads[stack_id])
        for name in LOAD_LIST_NAMES[:8]:
            with self.subTest(name=name):
                self.assertEqual(getattr(self.core_unsteady_problem, name), [])

    def test_recording_again_replaces(self) -> None:
        """Test that a second call replaces the lists instead of extending them."""
        self.record(cycle_averaged=False)
        self.record(cycle_averaged=False)

        self.assertEqual(len(self.core_unsteady_problem.finalForces_W), 2)
