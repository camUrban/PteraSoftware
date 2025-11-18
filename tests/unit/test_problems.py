"""This module contains a class to test SteadyProblems and UnsteadyProblems."""

import unittest
import math

import pterasoftware as ps
from tests.unit.fixtures import geometry_fixtures
from tests.unit.fixtures import operating_point_fixtures
from tests.unit.fixtures import movement_fixtures
from tests.unit.fixtures import problem_fixtures


class TestSteadyProblem(unittest.TestCase):
    """This is a class with functions to test SteadyProblems."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all SteadyProblem tests."""
        # Create fixtures using the problem_fixtures module.
        cls.basic_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        cls.multi_airplane_steady_problem = (
            problem_fixtures.make_multi_airplane_steady_problem_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test SteadyProblem initialization with valid parameters."""
        # Test that basic SteadyProblem initializes correctly.
        self.assertIsInstance(
            self.basic_steady_problem,
            ps.problems.SteadyProblem,
        )
        self.assertIsInstance(self.basic_steady_problem.airplanes, list)
        self.assertEqual(len(self.basic_steady_problem.airplanes), 1)
        self.assertIsInstance(
            self.basic_steady_problem.airplanes[0],
            ps.geometry.airplane.Airplane,
        )
        self.assertIsInstance(
            self.basic_steady_problem.operating_point,
            ps.operating_point.OperatingPoint,
        )

    def test_initialization_multiple_airplanes(self):
        """Test SteadyProblem initialization with multiple Airplanes."""
        # Test that SteadyProblem with multiple Airplanes initializes correctly.
        self.assertIsInstance(
            self.multi_airplane_steady_problem,
            ps.problems.SteadyProblem,
        )
        self.assertEqual(len(self.multi_airplane_steady_problem.airplanes), 2)
        for airplane in self.multi_airplane_steady_problem.airplanes:
            self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_airplanes_parameter_validation_not_list(self):
        """Test that airplanes parameter must be a list."""
        # Test with single Airplane instead of list.
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=geometry_fixtures.make_basic_airplane_fixture(),
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

        # Test with None.
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=None,
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

        # Test with invalid types.
        invalid_airplanes = ["not_a_list", 123, {"key": "value"}]
        for invalid in invalid_airplanes:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.SteadyProblem(
                        airplanes=invalid,
                        operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
                    )

    def test_airplanes_parameter_validation_empty_list(self):
        """Test that airplanes list must have at least one element."""
        with self.assertRaises(ValueError):
            ps.problems.SteadyProblem(
                airplanes=[],
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

    def test_airplanes_parameter_validation_elements_type(self):
        """Test that all elements in airplanes must be Airplanes."""
        # Test with list containing non-Airplane elements.
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=["not_an_airplane"],
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

        # Test with mixed valid and invalid elements.
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=[
                    geometry_fixtures.make_basic_airplane_fixture(),
                    "not_an_airplane",
                ],
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

    def test_operating_point_parameter_validation(self):
        """Test that operating_point parameter is properly validated."""
        # Test with invalid operating_point type.
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=[geometry_fixtures.make_basic_airplane_fixture()],
                operating_point="not_an_operating_point",
            )

        # Test with None.
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=[geometry_fixtures.make_basic_airplane_fixture()],
                operating_point=None,
            )

        # Test with other invalid types.
        invalid_operating_points = [123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_operating_points:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.SteadyProblem(
                        airplanes=[geometry_fixtures.make_basic_airplane_fixture()],
                        operating_point=invalid,
                    )


class TestUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test UnsteadyProblems."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all UnsteadyProblem tests."""
        # Create fixtures using the problem_fixtures module.
        cls.basic_unsteady_problem = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )
        cls.only_final_results_unsteady_problem = (
            problem_fixtures.make_only_final_results_unsteady_problem_fixture()
        )
        cls.static_unsteady_problem = (
            problem_fixtures.make_static_unsteady_problem_fixture()
        )
        cls.cyclic_unsteady_problem = (
            problem_fixtures.make_cyclic_unsteady_problem_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test UnsteadyProblem initialization with valid parameters."""
        # Test that basic UnsteadyProblem initializes correctly.
        self.assertIsInstance(
            self.basic_unsteady_problem,
            ps.problems.UnsteadyProblem,
        )
        self.assertIsInstance(
            self.basic_unsteady_problem.movement,
            ps.movements.movement.Movement,
        )
        self.assertFalse(self.basic_unsteady_problem.only_final_results)

    def test_initialization_only_final_results_true(self):
        """Test UnsteadyProblem initialization with only_final_results=True."""
        # Test that UnsteadyProblem with only_final_results=True initializes
        # correctly.
        self.assertIsInstance(
            self.only_final_results_unsteady_problem,
            ps.problems.UnsteadyProblem,
        )
        self.assertTrue(self.only_final_results_unsteady_problem.only_final_results)

    def test_only_final_results_parameter_validation(self):
        """Test only_final_results parameter validation."""
        # Test with valid bool values.
        movement = movement_fixtures.make_basic_movement_fixture()

        valid_values = [True, False]
        for value in valid_values:
            with self.subTest(value=value):
                unsteady_problem = ps.problems.UnsteadyProblem(
                    movement=movement,
                    only_final_results=value,
                )
                self.assertEqual(unsteady_problem.only_final_results, value)

    def test_movement_parameter_validation(self):
        """Test that movement parameter is properly validated."""
        # Test with invalid movement type.
        with self.assertRaises(TypeError):
            ps.problems.UnsteadyProblem(
                movement="not_a_movement",
            )

        # Test with None.
        with self.assertRaises(TypeError):
            ps.problems.UnsteadyProblem(
                movement=None,
            )

        # Test with other invalid types.
        invalid_movements = [123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_movements:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.UnsteadyProblem(
                        movement=invalid,
                    )

    def test_num_steps_attribute(self):
        """Test that num_steps is set correctly from Movement."""
        # Test that num_steps matches the Movement's num_steps.
        self.assertEqual(
            self.basic_unsteady_problem.num_steps,
            self.basic_unsteady_problem.movement.num_steps,
        )

    def test_delta_time_attribute(self):
        """Test that delta_time is set correctly from Movement."""
        # Test that delta_time matches the Movement's delta_time.
        self.assertEqual(
            self.basic_unsteady_problem.delta_time,
            self.basic_unsteady_problem.movement.delta_time,
        )

    def test_first_averaging_step_static_movement(self):
        """Test first_averaging_step for static Movement."""
        # For static Movement (max_period = 0), first_averaging_step should be
        # num_steps - 1.
        expected_first_averaging_step = self.static_unsteady_problem.num_steps - 1
        self.assertEqual(
            self.static_unsteady_problem.first_averaging_step,
            expected_first_averaging_step,
        )

    def test_first_averaging_step_cyclic_movement(self):
        """Test first_averaging_step for cyclic Movement."""
        # For cyclic Movement (max_period > 0), first_averaging_step should be
        # calculated based on the max_period.
        movement_max_period = self.cyclic_unsteady_problem.movement.max_period
        expected_first_averaging_step = max(
            0,
            math.floor(
                self.cyclic_unsteady_problem.num_steps
                - (movement_max_period / self.cyclic_unsteady_problem.delta_time)
            ),
        )
        self.assertEqual(
            self.cyclic_unsteady_problem.first_averaging_step,
            expected_first_averaging_step,
        )

    def test_first_results_step_only_final_results_false(self):
        """Test first_results_step when only_final_results is False."""
        # When only_final_results is False, first_results_step should be 0.
        self.assertEqual(self.basic_unsteady_problem.first_results_step, 0)

    def test_first_results_step_only_final_results_true(self):
        """Test first_results_step when only_final_results is True."""
        # When only_final_results is True, first_results_step should equal
        # first_averaging_step.
        self.assertEqual(
            self.only_final_results_unsteady_problem.first_results_step,
            self.only_final_results_unsteady_problem.first_averaging_step,
        )

    def test_initialization_of_load_lists(self):
        """Test that load lists are initialized as empty."""
        # All load lists should be initialized as empty lists.
        self.assertIsInstance(self.basic_unsteady_problem.finalForces_W, list)
        self.assertEqual(len(self.basic_unsteady_problem.finalForces_W), 0)

        self.assertIsInstance(
            self.basic_unsteady_problem.finalForceCoefficients_W, list
        )
        self.assertEqual(len(self.basic_unsteady_problem.finalForceCoefficients_W), 0)

        self.assertIsInstance(self.basic_unsteady_problem.finalMoments_W_CgP1, list)
        self.assertEqual(len(self.basic_unsteady_problem.finalMoments_W_CgP1), 0)

        self.assertIsInstance(
            self.basic_unsteady_problem.finalMomentCoefficients_W_CgP1, list
        )
        self.assertEqual(
            len(self.basic_unsteady_problem.finalMomentCoefficients_W_CgP1), 0
        )

        self.assertIsInstance(self.basic_unsteady_problem.finalMeanForces_W, list)
        self.assertEqual(len(self.basic_unsteady_problem.finalMeanForces_W), 0)

        self.assertIsInstance(
            self.basic_unsteady_problem.finalMeanForceCoefficients_W, list
        )
        self.assertEqual(
            len(self.basic_unsteady_problem.finalMeanForceCoefficients_W), 0
        )

        self.assertIsInstance(self.basic_unsteady_problem.finalMeanMoments_W_CgP1, list)
        self.assertEqual(len(self.basic_unsteady_problem.finalMeanMoments_W_CgP1), 0)

        self.assertIsInstance(
            self.basic_unsteady_problem.finalMeanMomentCoefficients_W_CgP1, list
        )
        self.assertEqual(
            len(self.basic_unsteady_problem.finalMeanMomentCoefficients_W_CgP1), 0
        )

        self.assertIsInstance(self.basic_unsteady_problem.finalRmsForces_W, list)
        self.assertEqual(len(self.basic_unsteady_problem.finalRmsForces_W), 0)

        self.assertIsInstance(
            self.basic_unsteady_problem.finalRmsForceCoefficients_W, list
        )
        self.assertEqual(
            len(self.basic_unsteady_problem.finalRmsForceCoefficients_W), 0
        )

        self.assertIsInstance(self.basic_unsteady_problem.finalRmsMoments_W_CgP1, list)
        self.assertEqual(len(self.basic_unsteady_problem.finalRmsMoments_W_CgP1), 0)

        self.assertIsInstance(
            self.basic_unsteady_problem.finalRmsMomentCoefficients_W_CgP1, list
        )
        self.assertEqual(
            len(self.basic_unsteady_problem.finalRmsMomentCoefficients_W_CgP1), 0
        )

    def test_steady_problems_list_initialization(self):
        """Test that steady_problems list is initialized correctly."""
        # steady_problems list should be initialized with correct length.
        self.assertIsInstance(self.basic_unsteady_problem.steady_problems, list)
        self.assertEqual(
            len(self.basic_unsteady_problem.steady_problems),
            self.basic_unsteady_problem.num_steps,
        )

    def test_steady_problems_list_elements_type(self):
        """Test that all elements in steady_problems are SteadyProblems."""
        # All elements in steady_problems should be SteadyProblems.
        for steady_problem in self.basic_unsteady_problem.steady_problems:
            self.assertIsInstance(steady_problem, ps.problems.SteadyProblem)

    def test_steady_problems_list_airplanes(self):
        """Test that each SteadyProblem has correct Airplanes."""
        # Each SteadyProblem should have the same number of Airplanes as the
        # Movement has AirplaneMovements.
        num_airplanes = len(self.basic_unsteady_problem.movement.airplane_movements)

        for steady_problem in self.basic_unsteady_problem.steady_problems:
            self.assertEqual(len(steady_problem.airplanes), num_airplanes)
            for airplane in steady_problem.airplanes:
                self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_steady_problems_list_operating_points(self):
        """Test that each SteadyProblem has an OperatingPoint."""
        # Each SteadyProblem should have an OperatingPoint.
        for steady_problem in self.basic_unsteady_problem.steady_problems:
            self.assertIsInstance(
                steady_problem.operating_point, ps.operating_point.OperatingPoint
            )


if __name__ == "__main__":
    unittest.main()
