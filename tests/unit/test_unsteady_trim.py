import unittest

import pterasoftware as ps
from tests.unit.fixtures import problem_fixtures


class TestAnalyzeUnsteadyTrim(unittest.TestCase):
    """This class contains unit tests for the analyze_unsteady_trim function."""

    def setUp(self):
        """Set up test fixtures."""
        self.problem = problem_fixtures.make_basic_unsteady_problem_fixture()

    def test_problem_validation(self):
        """Test problem parameter validation."""

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem="not a problem",
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_multiple_airplane_movements_validation(self):
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

    def test_boundsVCg__E_validation(self):
        """Test boundsVCg__E parameter validation."""

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E="invalid",
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0,),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=("a", 10.0),
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

    def test_alpha_bounds_validation(self):
        """Test alpha_bounds parameter validation."""

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds="invalid",
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(1.0,),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=("a", 10.0),
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

    def test_beta_bounds_validation(self):
        """Test beta_bounds parameter validation."""

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds="invalid",
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(1.0,),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=("a", 10.0),
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

    def test_boundsExternalFX_W_validation(self):
        """Test boundsExternalFX_W parameter validation."""

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W="invalid",
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(1.0,),
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=("a", 10.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(10.0, -10.0),
            )

    def test_objective_cut_off_validation(self):
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

    def test_num_calls_validation(self):
        """Test num_calls parameter validation."""

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=0,
            )

    def test_base_operating_point_parameter_bounds_validation(self):
        """Test that the base operating point values must lie within the supplied
        bounds.
        """

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
