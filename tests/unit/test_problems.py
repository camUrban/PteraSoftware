"""This module contains classes to test SteadyProblems and UnsteadyProblems."""

import tempfile
import unittest
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pyvista as pv

import pterasoftware as ps
from tests.unit.fixtures import (
    geometry_fixtures,
    movement_fixtures,
    operating_point_fixtures,
    problem_fixtures,
)


class TestSteadyProblem(unittest.TestCase):
    """This is a class with functions to test SteadyProblems."""

    basic_steady_problem: ps.problems.SteadyProblem
    multi_airplane_steady_problem: ps.problems.SteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all SteadyProblem tests."""
        # Create fixtures using the problem_fixtures module.
        cls.basic_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        cls.multi_airplane_steady_problem = (
            problem_fixtures.make_multi_airplane_steady_problem_fixture()
        )

    def test_initialization_valid_parameters(self) -> None:
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

    def test_initialization_multiple_airplanes(self) -> None:
        """Test SteadyProblem initialization with multiple Airplanes."""
        # Test that SteadyProblem with multiple Airplanes initializes correctly.
        self.assertIsInstance(
            self.multi_airplane_steady_problem,
            ps.problems.SteadyProblem,
        )
        self.assertEqual(len(self.multi_airplane_steady_problem.airplanes), 2)
        for airplane in self.multi_airplane_steady_problem.airplanes:
            self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_airplanes_parameter_validation_not_list(self) -> None:
        """Test that airplanes parameter must be a list."""
        # Test with single Airplane instead of list.
        single_airplane: Any = geometry_fixtures.make_basic_airplane_fixture()
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=single_airplane,
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

        # Test with None.
        none_airplanes: Any = None
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=none_airplanes,
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

        # Test with invalid types.
        invalid_airplanes: list[Any] = ["not_a_list", 123, {"key": "value"}]
        for invalid in invalid_airplanes:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.SteadyProblem(
                        airplanes=invalid,
                        operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
                    )

    def test_airplanes_parameter_validation_empty_list(self) -> None:
        """Test that airplanes list must have at least one element."""
        with self.assertRaises(ValueError):
            ps.problems.SteadyProblem(
                airplanes=[],
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

    def test_airplanes_parameter_validation_elements_type(self) -> None:
        """Test that all elements in airplanes must be Airplanes."""
        # Test with list containing non-Airplane elements.
        non_airplane_elements: list[Any] = ["not_an_airplane"]
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=non_airplane_elements,
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

        # Test with mixed valid and invalid elements.
        mixed_elements: list[Any] = [
            geometry_fixtures.make_basic_airplane_fixture(),
            "not_an_airplane",
        ]
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=mixed_elements,
                operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
            )

    def test_operating_point_parameter_validation(self) -> None:
        """Test that operating_point parameter is properly validated."""
        # Test with invalid operating_point type.
        bad_operating_point: Any = "not_an_operating_point"
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=[geometry_fixtures.make_basic_airplane_fixture()],
                operating_point=bad_operating_point,
            )

        # Test with None.
        none_operating_point: Any = None
        with self.assertRaises(TypeError):
            ps.problems.SteadyProblem(
                airplanes=[geometry_fixtures.make_basic_airplane_fixture()],
                operating_point=none_operating_point,
            )

        # Test with other invalid types.
        invalid_operating_points: list[Any] = [123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_operating_points:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.SteadyProblem(
                        airplanes=[geometry_fixtures.make_basic_airplane_fixture()],
                        operating_point=invalid,
                    )

    def test_reynolds_numbers_returns_correct_tuple_length(self) -> None:
        """Test that reynolds_numbers returns a tuple with one element per Airplane."""
        # Single Airplane problem should return tuple with one element.
        self.assertIsInstance(self.basic_steady_problem.reynolds_numbers, tuple)
        self.assertEqual(len(self.basic_steady_problem.reynolds_numbers), 1)

        # Multi Airplane problem should return tuple with two elements.
        self.assertIsInstance(
            self.multi_airplane_steady_problem.reynolds_numbers, tuple
        )
        self.assertEqual(len(self.multi_airplane_steady_problem.reynolds_numbers), 2)

    def test_reynolds_numbers_calculation_accuracy(self) -> None:
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

    def test_reynolds_numbers_multiple_airplanes(self) -> None:
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

    basic_steady_problem: ps.problems.SteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all immutability tests."""
        cls.basic_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

    def test_immutable_airplanes_property(self) -> None:
        """Test that airplanes property is read only."""
        new_airplanes = (geometry_fixtures.make_basic_airplane_fixture(),)
        with self.assertRaises(AttributeError):
            setattr(self.basic_steady_problem, "airplanes", new_airplanes)

    def test_immutable_operating_point_property(self) -> None:
        """Test that operating_point property is read only."""
        new_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        with self.assertRaises(AttributeError):
            setattr(self.basic_steady_problem, "operating_point", new_operating_point)

    def test_airplanes_tuple_immutability(self) -> None:
        """Test that airplanes tuple cannot be modified via append or other methods."""
        # Tuples don't have append, so attempting to call it raises AttributeError.
        airplanes: Any = self.basic_steady_problem.airplanes
        with self.assertRaises(AttributeError):
            airplanes.append(geometry_fixtures.make_basic_airplane_fixture())

    def test_reynolds_numbers_caching(self) -> None:
        """Test that reynolds_numbers returns the same cached object on repeated
        access."""
        # Access reynolds_numbers twice.
        reynolds_first = self.basic_steady_problem.reynolds_numbers
        reynolds_second = self.basic_steady_problem.reynolds_numbers

        # Should return the same cached tuple object.
        self.assertIs(reynolds_first, reynolds_second)


class TestSteadyProblemPanelCoordinates(unittest.TestCase):
    """Tests for SteadyProblem Panel GP1_CgP1 coordinate population."""

    def test_panel_GP1_CgP1_coordinates_are_read_only(self) -> None:
        """Test that Panel GP1_CgP1 coordinate arrays are read only."""
        # Create a fresh SteadyProblem.
        steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

        # Get the first Panel.
        first_airplane = steady_problem.airplanes[0]
        first_wing = first_airplane.wings[0]
        self.assertIsNotNone(first_wing.panels)
        assert first_wing.panels is not None
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

    def test_panel_GP1_CgP1_coordinates_multi_airplane(self) -> None:
        """Test that Panel GP1_CgP1 coordinates are populated for multiple Airplanes."""
        # Create a SteadyProblem with multiple Airplanes.
        steady_problem = problem_fixtures.make_multi_airplane_steady_problem_fixture()

        # Check that all Panels in all Airplanes have GP1_CgP1 coordinates set.
        for airplane in steady_problem.airplanes:
            for wing in airplane.wings:
                self.assertIsNotNone(wing.panels)
                assert wing.panels is not None
                for panel in np.ravel(wing.panels):
                    self.assertIsNotNone(panel.Frpp_GP1_CgP1)
                    self.assertIsNotNone(panel.Flpp_GP1_CgP1)
                    self.assertIsNotNone(panel.Blpp_GP1_CgP1)
                    self.assertIsNotNone(panel.Brpp_GP1_CgP1)


class TestSteadyProblemDiagram(unittest.TestCase):
    """Tests for SteadyProblem.diagram method."""

    def setUp(self) -> None:
        """Set up test fixtures for diagram tests."""
        self.basic_steady_problem = problem_fixtures.make_basic_steady_problem_fixture()
        self.multi_airplane_steady_problem = (
            problem_fixtures.make_multi_airplane_steady_problem_fixture()
        )

    def test_diagram_shows_the_diagram(self) -> None:
        """Test that diagram shows the diagram once."""
        # Patch the Plotter's show method to avoid blocking on window close.
        with patch.object(pv.Plotter, "show") as mock_show:
            self.basic_steady_problem.diagram()

        mock_show.assert_called_once()

    def test_diagram_accepts_numpy_bools(self) -> None:
        """Test that diagram accepts numpy bools for its flags."""
        with patch.object(pv.Plotter, "show") as mock_show:
            self.basic_steady_problem.diagram(
                show_airplane_axes_and_points=np.bool(True),
                show_wing_axes_and_points=np.bool(True),
                show_wing_cross_section_axes_and_points=np.bool(True),
                show_airfoil_axes_and_points=np.bool(True),
                show_airfoils=np.bool(True),
                show_mcls=np.bool(True),
                show_collocation_points=np.bool(True),
                label_collocation_points=np.bool(True),
                math_labels=np.bool(True),
                save=np.bool(False),
            )

        mock_show.assert_called_once()

    def test_diagram_multi_airplane_steady_problem(self) -> None:
        """Test diagram with a SteadyProblem with multiple Airplanes."""
        with patch.object(pv.Plotter, "show") as mock_show:
            self.multi_airplane_steady_problem.diagram()

        mock_show.assert_called_once()

    def test_diagram_saves_the_diagram(self) -> None:
        """Test that diagram saves a WebP to the given path.

        Saving takes a screenshot of the shown Plotter, so show can't be patched out
        here. Rendering off screen instead keeps show from opening a window and blocking
        until it is closed.
        """
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            saved_path = Path(temporary_directory_name) / "steady_problem.webp"
            with patch.object(pv, "OFF_SCREEN", True):
                self.basic_steady_problem.diagram(save=True, path=saved_path)

            self.assertTrue(saved_path.is_file())
            self.assertGreater(saved_path.stat().st_size, 0)

    def test_diagram_invalid_save_type_raises(self) -> None:
        """Test that diagram raises error for invalid save type."""
        bad_save: Any = "invalid"
        with self.assertRaises(TypeError):
            # noinspection PyTypeChecker
            self.basic_steady_problem.diagram(save=bad_save)


class TestUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test UnsteadyProblems."""

    basic_unsteady_problem: ps.problems.UnsteadyProblem
    only_final_results_unsteady_problem: ps.problems.UnsteadyProblem
    multi_airplane_unsteady_problem: ps.problems.UnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all UnsteadyProblem tests."""
        # Create fixtures using the problem_fixtures module.
        cls.basic_unsteady_problem = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )
        cls.only_final_results_unsteady_problem = (
            problem_fixtures.make_only_final_results_unsteady_problem_fixture()
        )
        cls.multi_airplane_unsteady_problem = (
            problem_fixtures.make_multi_airplane_unsteady_problem_fixture()
        )

    def test_initialization_valid_parameters(self) -> None:
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

    def test_initialization_only_final_results_true(self) -> None:
        """Test UnsteadyProblem initialization with only_final_results=True."""
        # Test that UnsteadyProblem with only_final_results=True initializes correctly.
        self.assertIsInstance(
            self.only_final_results_unsteady_problem,
            ps.problems.UnsteadyProblem,
        )
        self.assertTrue(self.only_final_results_unsteady_problem.only_final_results)

    def test_only_final_results_parameter_validation(self) -> None:
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

    def test_movement_parameter_validation(self) -> None:
        """Test that movement parameter is properly validated."""
        # Test with invalid movement type.
        bad_movement: Any = "not_a_movement"
        with self.assertRaises(TypeError):
            ps.problems.UnsteadyProblem(
                movement=bad_movement,
            )

        # Test with None.
        none_movement: Any = None
        with self.assertRaises(TypeError):
            ps.problems.UnsteadyProblem(
                movement=none_movement,
            )

        # Test with other invalid types.
        invalid_movements: list[Any] = [123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_movements:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    ps.problems.UnsteadyProblem(
                        movement=invalid,
                    )

    def test_num_steps_attribute(self) -> None:
        """Test that num_steps is set correctly from Movement."""
        # Test that num_steps matches the Movement's num_steps.
        self.assertEqual(
            self.basic_unsteady_problem.num_steps,
            self.basic_unsteady_problem.movement.num_steps,
        )

    def test_delta_time_attribute(self) -> None:
        """Test that delta_time is set correctly from Movement."""
        # Test that delta_time matches the Movement's delta_time.
        self.assertEqual(
            self.basic_unsteady_problem.delta_time,
            self.basic_unsteady_problem.movement.delta_time,
        )

    def test_steady_problems_tuple_initialization(self) -> None:
        """Test that steady_problems tuple is initialized correctly."""
        # steady_problems tuple should be initialized with correct length.
        self.assertIsInstance(self.basic_unsteady_problem.steady_problems, tuple)
        self.assertEqual(
            len(self.basic_unsteady_problem.steady_problems),
            self.basic_unsteady_problem.num_steps,
        )

    def test_steady_problems_list_elements_type(self) -> None:
        """Test that all elements in steady_problems are SteadyProblems."""
        # All elements in steady_problems should be SteadyProblems.
        for steady_problem in self.basic_unsteady_problem.steady_problems:
            self.assertIsInstance(steady_problem, ps.problems.SteadyProblem)

    def test_steady_problems_list_airplanes(self) -> None:
        """Test that each SteadyProblem has correct Airplanes."""
        # Each SteadyProblem should have the same number of Airplanes as the Movement
        # has AirplaneMovements.
        num_airplanes = len(self.basic_unsteady_problem.movement.airplane_movements)

        for steady_problem in self.basic_unsteady_problem.steady_problems:
            self.assertEqual(len(steady_problem.airplanes), num_airplanes)
            for airplane in steady_problem.airplanes:
                self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_steady_problems_list_operating_points(self) -> None:
        """Test that each SteadyProblem has an OperatingPoint."""
        # Each SteadyProblem should have an OperatingPoint.
        for steady_problem in self.basic_unsteady_problem.steady_problems:
            self.assertIsInstance(
                steady_problem.operating_point, ps.operating_point.OperatingPoint
            )

    def test_initialization_multiple_airplanes(self) -> None:
        """Test UnsteadyProblem initialization with multiple Airplanes."""
        # Test that UnsteadyProblem with multiple Airplanes initializes correctly.
        self.assertIsInstance(
            self.multi_airplane_unsteady_problem,
            ps.problems.UnsteadyProblem,
        )
        # Verify that each SteadyProblem has multiple Airplanes.
        for steady_problem in self.multi_airplane_unsteady_problem.steady_problems:
            self.assertEqual(len(steady_problem.airplanes), 2)


class TestUnsteadyProblemDiagram(unittest.TestCase):
    """Tests for the UnsteadyProblem diagram method."""

    basic_unsteady_problem: ps.problems.UnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all diagram tests."""
        cls.basic_unsteady_problem = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )

    def test_diagram_shows_the_diagram(self) -> None:
        """Test that diagram shows the diagram once."""
        # Patch the Plotter's show method to avoid blocking on window close.
        with patch.object(pv.Plotter, "show") as mock_show:
            self.basic_unsteady_problem.diagram()

        mock_show.assert_called_once()

    def test_diagram_draws_the_last_time_step_by_default(self) -> None:
        """Test that diagram draws the last time step's SteadyProblem by default."""
        with patch.object(
            ps.problems.SteadyProblem, "diagram", autospec=True
        ) as mock_diagram:
            self.basic_unsteady_problem.diagram()

        mock_diagram.assert_called_once()
        self.assertIs(
            mock_diagram.call_args.args[0],
            self.basic_unsteady_problem.steady_problems[-1],
        )

    def test_diagram_draws_the_given_time_step(self) -> None:
        """Test that diagram draws the SteadyProblem of the time step it is given."""
        with patch.object(
            ps.problems.SteadyProblem, "diagram", autospec=True
        ) as mock_diagram:
            self.basic_unsteady_problem.diagram(step=1)

        mock_diagram.assert_called_once()
        self.assertIs(
            mock_diagram.call_args.args[0],
            self.basic_unsteady_problem.steady_problems[1],
        )

    def test_diagram_rejects_an_out_of_range_step(self) -> None:
        """Test that diagram rejects a step past either end of the time steps."""
        num_steps = self.basic_unsteady_problem.num_steps
        for step in [num_steps, -num_steps - 1]:
            with self.subTest(step=step):
                with self.assertRaises(ValueError):
                    self.basic_unsteady_problem.diagram(step=step)

    def test_diagram_warns_for_a_negative_step_before_solving(self) -> None:
        """Test that diagram warns when a negative step counts back from the last time
        step created so far rather than from the last time step.

        An AeroelasticUnsteadyProblem creates each time step's SteadyProblem while it is
        solved, so before then only the first time step's exists.
        """
        aeroelastic_unsteady_problem = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )
        with patch.object(pv.Plotter, "show"):
            with self.assertLogs("pterasoftware.core", level="WARNING") as context:
                aeroelastic_unsteady_problem.diagram()
        self.assertIn("hasn't been solved", context.output[0])

    def test_diagram_does_not_warn_for_a_nonnegative_step_before_solving(
        self,
    ) -> None:
        """Test that diagram doesn't warn when a nonnegative step names a time step
        whose SteadyProblem has been created."""
        aeroelastic_unsteady_problem = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )
        with patch.object(pv.Plotter, "show"):
            with self.assertNoLogs("pterasoftware.core", level="WARNING"):
                aeroelastic_unsteady_problem.diagram(step=0)

    def test_diagram_rejects_a_step_not_yet_created(self) -> None:
        """Test that diagram rejects a step whose SteadyProblem hasn't been created
        yet."""
        aeroelastic_unsteady_problem = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )
        with self.assertRaises(ValueError):
            aeroelastic_unsteady_problem.diagram(step=1)


class TestUnsteadyProblemImmutability(unittest.TestCase):
    """Tests for UnsteadyProblem attribute immutability."""

    basic_unsteady_problem: ps.problems.UnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all immutability tests."""
        cls.basic_unsteady_problem = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )

    def test_immutable_movement_property(self) -> None:
        """Test that movement property is read only."""
        new_movement = movement_fixtures.make_basic_movement_fixture()
        with self.assertRaises(AttributeError):
            setattr(self.basic_unsteady_problem, "movement", new_movement)

    def test_immutable_steady_problems_property(self) -> None:
        """Test that steady_problems property is read only."""
        with self.assertRaises(AttributeError):
            setattr(self.basic_unsteady_problem, "steady_problems", ())

    def test_steady_problems_tuple_immutability(self) -> None:
        """Test that steady_problems tuple cannot be modified via append or other
        methods."""
        # Tuples don't have append, so attempting to call it raises AttributeError.
        steady_problems: Any = self.basic_unsteady_problem.steady_problems
        with self.assertRaises(AttributeError):
            steady_problems.append(problem_fixtures.make_basic_steady_problem_fixture())
