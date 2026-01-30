"""This module contains classes to test SteadyProblems and UnsteadyProblems."""

import math
import unittest

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import (
    geometry_fixtures,
    movement_fixtures,
    operating_point_fixtures,
    problem_fixtures,
)


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
        self.assertIsInstance(self.basic_steady_problem.airplanes, tuple)
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

    def test_reynolds_numbers_returns_correct_tuple_length(self):
        """Test that reynolds_numbers returns a tuple with one element per Airplane."""
        # Single Airplane problem should return tuple with one element.
        self.assertIsInstance(self.basic_steady_problem.reynolds_numbers, tuple)
        self.assertEqual(len(self.basic_steady_problem.reynolds_numbers), 1)

        # Multi Airplane problem should return tuple with two elements.
        self.assertIsInstance(
            self.multi_airplane_steady_problem.reynolds_numbers, tuple
        )
        self.assertEqual(len(self.multi_airplane_steady_problem.reynolds_numbers), 2)

    def test_reynolds_numbers_calculation_accuracy(self):
        """Test that reynolds_numbers calculates Re = (V x L) / nu correctly."""
        # Get the values used in calculation.
        v = self.basic_steady_problem.operating_point.vCg__E
        nu = self.basic_steady_problem.operating_point.nu
        c_ref = self.basic_steady_problem.airplanes[0].c_ref

        # Calculate expected Reynolds number.
        expected_re = (v * c_ref) / nu

        # Check the calculated value matches.
        calculated_re = self.basic_steady_problem.reynolds_numbers[0]
        self.assertAlmostEqual(calculated_re, expected_re, places=6)

    def test_reynolds_numbers_multiple_airplanes(self):
        """Test reynolds_numbers with multiple Airplanes with different c_ref."""
        # Get OperatingPoint values.
        v = self.multi_airplane_steady_problem.operating_point.vCg__E
        nu = self.multi_airplane_steady_problem.operating_point.nu

        # Check each Airplane's Reynolds number.
        for i, airplane in enumerate(self.multi_airplane_steady_problem.airplanes):
            expected_re = (v * airplane.c_ref) / nu
            calculated_re = self.multi_airplane_steady_problem.reynolds_numbers[i]
            self.assertAlmostEqual(calculated_re, expected_re, places=6)


class TestSteadyProblemImmutability(unittest.TestCase):
    """Tests for SteadyProblem attribute immutability."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all immutability tests."""
        cls.basic_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

    def test_immutable_airplanes_property(self):
        """Test that airplanes property is read only."""
        new_airplanes = (geometry_fixtures.make_basic_airplane_fixture(),)
        with self.assertRaises(AttributeError):
            self.basic_steady_problem.airplanes = new_airplanes

    def test_immutable_operating_point_property(self):
        """Test that operating_point property is read only."""
        new_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        with self.assertRaises(AttributeError):
            self.basic_steady_problem.operating_point = new_operating_point

    def test_airplanes_tuple_immutability(self):
        """Test that airplanes tuple cannot be modified via append or other methods."""
        # Tuples don't have append, so attempting to call it raises AttributeError.
        with self.assertRaises(AttributeError):
            self.basic_steady_problem.airplanes.append(
                geometry_fixtures.make_basic_airplane_fixture()
            )

    def test_reynolds_numbers_caching(self):
        """Test that reynolds_numbers returns the same cached object on repeated access."""
        # Access reynolds_numbers twice.
        reynolds_first = self.basic_steady_problem.reynolds_numbers
        reynolds_second = self.basic_steady_problem.reynolds_numbers

        # Should return the same cached tuple object.
        self.assertIs(reynolds_first, reynolds_second)


class TestSteadyProblemPanelCoordinates(unittest.TestCase):
    """Tests for SteadyProblem Panel GP1_CgP1 coordinate population."""

    def test_panel_GP1_CgP1_coordinates_are_read_only(self):
        """Test that Panel GP1_CgP1 coordinate arrays are read only."""
        # Create a fresh SteadyProblem.
        steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

        # Get the first Panel.
        first_airplane = steady_problem.airplanes[0]
        first_wing = first_airplane.wings[0]
        self.assertIsNotNone(first_wing.panels)
        first_panel = first_wing.panels[0, 0]

        # Verify that the arrays are read only.
        with self.assertRaises(ValueError):
            first_panel.Frpp_GP1_CgP1[0] = 999.0

        with self.assertRaises(ValueError):
            first_panel.Flpp_GP1_CgP1[0] = 999.0

        with self.assertRaises(ValueError):
            first_panel.Blpp_GP1_CgP1[0] = 999.0

        with self.assertRaises(ValueError):
            first_panel.Brpp_GP1_CgP1[0] = 999.0

    def test_panel_GP1_CgP1_coordinates_multi_airplane(self):
        """Test that Panel GP1_CgP1 coordinates are populated for multiple Airplanes."""
        # Create a SteadyProblem with multiple Airplanes.
        steady_problem = problem_fixtures.make_multi_airplane_steady_problem_fixture()

        # Check that all Panels in all Airplanes have GP1_CgP1 coordinates set.
        for airplane in steady_problem.airplanes:
            for wing in airplane.wings:
                self.assertIsNotNone(wing.panels)
                for panel in np.ravel(wing.panels):
                    self.assertIsNotNone(panel.Frpp_GP1_CgP1)
                    self.assertIsNotNone(panel.Flpp_GP1_CgP1)
                    self.assertIsNotNone(panel.Blpp_GP1_CgP1)
                    self.assertIsNotNone(panel.Brpp_GP1_CgP1)


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
        cls.multi_airplane_unsteady_problem = (
            problem_fixtures.make_multi_airplane_unsteady_problem_fixture()
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
        # Test with valid bool values. A fresh movement fixture is needed for each
        # iteration because UnsteadyProblem sets attributes on Panels that can only be
        # set once.
        valid_values = [True, False]
        for value in valid_values:
            with self.subTest(value=value):
                movement = movement_fixtures.make_basic_movement_fixture()
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
        # For cyclic Movement (lcm_period > 0), first_averaging_step should be
        # calculated based on the lcm_period.
        movement_lcm_period = self.cyclic_unsteady_problem.movement.lcm_period
        expected_first_averaging_step = max(
            0,
            math.floor(
                self.cyclic_unsteady_problem.num_steps
                - (movement_lcm_period / self.cyclic_unsteady_problem.delta_time)
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

    def test_steady_problems_tuple_initialization(self):
        """Test that steady_problems tuple is initialized correctly."""
        # steady_problems tuple should be initialized with correct length.
        self.assertIsInstance(self.basic_unsteady_problem.steady_problems, tuple)
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

    def test_only_final_results_accepts_numpy_bool(self):
        """Test that only_final_results accepts numpy bool values."""
        # Create a fresh movement fixture for this test.
        movement = movement_fixtures.make_basic_movement_fixture()
        unsteady_problem = ps.problems.UnsteadyProblem(
            movement=movement,
            only_final_results=np.bool_(True),
        )
        self.assertTrue(unsteady_problem.only_final_results)
        self.assertIsInstance(unsteady_problem.only_final_results, bool)

    def test_initialization_multiple_airplanes(self):
        """Test UnsteadyProblem initialization with multiple Airplanes."""
        # Test that UnsteadyProblem with multiple Airplanes initializes correctly.
        self.assertIsInstance(
            self.multi_airplane_unsteady_problem,
            ps.problems.UnsteadyProblem,
        )
        # Verify that each SteadyProblem has multiple Airplanes.
        for steady_problem in self.multi_airplane_unsteady_problem.steady_problems:
            self.assertEqual(len(steady_problem.airplanes), 2)


class TestUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for UnsteadyProblem attribute immutability."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all immutability tests."""
        cls.basic_unsteady_problem = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )

    def test_immutable_movement_property(self):
        """Test that movement property is read only."""
        new_movement = movement_fixtures.make_basic_movement_fixture()
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.movement = new_movement

    def test_immutable_only_final_results_property(self):
        """Test that only_final_results property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.only_final_results = True

    def test_immutable_num_steps_property(self):
        """Test that num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.num_steps = 100

    def test_immutable_delta_time_property(self):
        """Test that delta_time property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.delta_time = 0.1

    def test_immutable_first_averaging_step_property(self):
        """Test that first_averaging_step property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.first_averaging_step = 0

    def test_immutable_first_results_step_property(self):
        """Test that first_results_step property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.first_results_step = 0

    def test_immutable_steady_problems_property(self):
        """Test that steady_problems property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.steady_problems = ()

    def test_steady_problems_tuple_immutability(self):
        """Test that steady_problems tuple cannot be modified via append or other
        methods.
        """
        # Tuples don't have append, so attempting to call it raises AttributeError.
        with self.assertRaises(AttributeError):
            self.basic_unsteady_problem.steady_problems.append(
                problem_fixtures.make_basic_steady_problem_fixture()
            )

    def test_mutable_load_lists(self):
        """Test that load lists remain mutable for solver population."""
        # The load lists should be mutable so the solver can populate them.
        self.basic_unsteady_problem.finalForces_W.append(np.array([1.0, 2.0, 3.0]))
        self.assertEqual(len(self.basic_unsteady_problem.finalForces_W), 1)

        # Clean up.
        self.basic_unsteady_problem.finalForces_W.pop()


if __name__ == "__main__":
    unittest.main()
