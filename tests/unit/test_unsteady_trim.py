"""This module contains a class to test the analyze_unsteady_trim function."""

import unittest
from typing import Any

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures


class TestAnalyzeUnsteadyTrim(unittest.TestCase):
    """A class with functions to test analyze_unsteady_trim."""

    problem: ps.problems.UnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all analyze_unsteady_trim tests."""
        cls.problem = problem_fixtures.make_basic_unsteady_problem_fixture()

    def test_problem_validation(self) -> None:
        """Test problem parameter validation."""
        bad_problem: Any = "not a problem"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=bad_problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_multiple_airplane_movements_validation(self) -> None:
        """Test that only one AirplaneMovement is allowed."""
        problem = problem_fixtures.make_multi_airplane_unsteady_problem_fixture()

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_boundsVCg__E_validation(self) -> None:
        """Test boundsVCg__E parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=bad_str,
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=bad_short_tuple,
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=bad_element_tuple,
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(10.0, 1.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(0.0, 10.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_alpha_bounds_validation(self) -> None:
        """Test alpha_bounds parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=bad_str,
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=bad_short_tuple,
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=bad_element_tuple,
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(10.0, -10.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_beta_bounds_validation(self) -> None:
        """Test beta_bounds parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=bad_str,
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=bad_short_tuple,
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=bad_element_tuple,
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(10.0, -10.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_boundsExternalFX_W_validation(self) -> None:
        """Test boundsExternalFX_W parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=bad_str,
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=bad_short_tuple,
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=bad_element_tuple,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(10.0, -10.0),
            )

    def test_objective_cut_off_validation(self) -> None:
        """Test objective_cut_off parameter validation."""
        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=0.0,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=-1.0,
            )

    def test_num_calls_validation(self) -> None:
        """Test num_calls parameter validation."""
        bad_num_calls: Any = 1.5
        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=0,
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=bad_num_calls,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=-1,
            )

    def test_base_operating_point_parameter_bounds_validation(self) -> None:
        """Test that the base operating point values must lie within the supplied
        bounds."""
        base_operating_point = (
            self.problem.movement.operating_point_movement.base_operating_point
        )

        with self.subTest(parameter="vCg__E"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(
                        base_operating_point.vCg__E + 1.0,
                        base_operating_point.vCg__E + 10.0,
                    ),
                    alpha_bounds=(-20.0, 20.0),
                    beta_bounds=(-20.0, 20.0),
                    boundsExternalFX_W=(-1000.0, 1000.0),
                )

        with self.subTest(parameter="alpha"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(1.0, 100.0),
                    alpha_bounds=(
                        base_operating_point.alpha + 1.0,
                        base_operating_point.alpha + 10.0,
                    ),
                    beta_bounds=(-20.0, 20.0),
                    boundsExternalFX_W=(-1000.0, 1000.0),
                )

        with self.subTest(parameter="beta"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(1.0, 100.0),
                    alpha_bounds=(-20.0, 20.0),
                    beta_bounds=(
                        base_operating_point.beta + 1.0,
                        base_operating_point.beta + 10.0,
                    ),
                    boundsExternalFX_W=(-1000.0, 1000.0),
                )

        with self.subTest(parameter="externalFX_W"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(1.0, 100.0),
                    alpha_bounds=(-20.0, 20.0),
                    beta_bounds=(-20.0, 20.0),
                    boundsExternalFX_W=(
                        base_operating_point.externalFX_W + 1.0,
                        base_operating_point.externalFX_W + 10.0,
                    ),
                )


if __name__ == "__main__":
    unittest.main()
