"""This module contains a testing case for the steady trim function."""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import airplane_fixtures


class TestSteadyTrimHorseshoeVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the steady trim function on a
    SteadHorseshoeVortexLatticeMethodSolver."""

    def setUp(self) -> None:
        """This method sets up the test.

        :return: None
        """
        self.v_x_ans = 2.9222951743478016
        self.alpha_ans = 1.933469345202583
        self.beta_ans = 0.000
        self.thrust_ans = 0.0884579818006783

        self.ans_corruption = 0.05

        corrupted_v_x = self.v_x_ans * (1 + self.ans_corruption)
        corrupted_alpha = self.alpha_ans * (1 + self.ans_corruption)
        corrupted_beta = self.beta_ans * (1 + self.ans_corruption)
        corrupted_thrust = self.thrust_ans * (1 + self.ans_corruption)

        this_airplane = (
            airplane_fixtures.make_multiple_wing_steady_validation_airplane()
        )
        this_operating_point = ps.operating_point.OperatingPoint(
            vCg__E=corrupted_v_x,
            alpha=corrupted_alpha,
            beta=corrupted_beta,
            externalFX_W=corrupted_thrust,
            g_E=(0.0, 0.0, 9.80665),
        )

        # Create the SteadyProblem.
        self.steady_validation_problem = ps.problems.SteadyProblem(
            airplanes=[this_airplane],
            operating_point=this_operating_point,
        )

    def test_g_E_validation(self) -> None:
        """This method tests that a zero g_E is rejected, since the trim analysis places
        the Airplane's weight along g_E's direction.

        :return: None
        """
        problem = ps.problems.SteadyProblem(
            airplanes=[
                airplane_fixtures.make_multiple_wing_steady_validation_airplane()
            ],
            operating_point=ps.operating_point.OperatingPoint(),
        )

        with self.assertRaisesRegex(ValueError, "g_E must be non-zero"):
            ps.trim.analyze_steady_trim(
                problem=problem,
                solver_type="steady horseshoe vortex lattice method",
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_angle_bounds_validation(self) -> None:
        """This method tests that alpha bounds outside (-180.0, 180.0] and beta bounds
        outside (-90.0, 90.0) are rejected, since the search builds an OperatingPoint
        from every trial and varies alpha and beta independently.

        :return: None
        """
        for bad_alpha_bounds in [(-180.0, 20.0), (-20.0, 180.001)]:
            with self.subTest(alpha_bounds=bad_alpha_bounds):
                with self.assertRaisesRegex(ValueError, "range \\(-180.0, 180.0\\]"):
                    ps.trim.analyze_steady_trim(
                        problem=self.steady_validation_problem,
                        solver_type="steady horseshoe vortex lattice method",
                        boundsVCg__E=(1.0, 100.0),
                        alpha_bounds=bad_alpha_bounds,
                        beta_bounds=(-20.0, 20.0),
                        boundsExternalFX_W=(-1000.0, 1000.0),
                    )

        for bad_beta_bounds in [(-90.0, 20.0), (-20.0, 90.0), (-180.0, 180.0)]:
            with self.subTest(beta_bounds=bad_beta_bounds):
                with self.assertRaisesRegex(ValueError, "range \\(-90.0, 90.0\\)"):
                    ps.trim.analyze_steady_trim(
                        problem=self.steady_validation_problem,
                        solver_type="steady horseshoe vortex lattice method",
                        boundsVCg__E=(1.0, 100.0),
                        alpha_bounds=(-20.0, 20.0),
                        beta_bounds=bad_beta_bounds,
                        boundsExternalFX_W=(-1000.0, 1000.0),
                    )

    def test_base_attitude_validation(self) -> None:
        """This method tests that a base attitude that does not resolve to level flight
        is rejected, since the trials resolve their own attitudes to level flight and
        would otherwise silently discard it.

        :return: None
        """
        problem = ps.problems.SteadyProblem(
            airplanes=[
                airplane_fixtures.make_multiple_wing_steady_validation_airplane()
            ],
            operating_point=ps.operating_point.OperatingPoint(
                angles_E_to_BP1_izyx=(0.0, 0.0, 0.0), g_E=(0.0, 0.0, 9.80665)
            ),
        )

        with self.assertRaisesRegex(ValueError, "must resolve to level flight"):
            ps.trim.analyze_steady_trim(
                problem=problem,
                solver_type="steady horseshoe vortex lattice method",
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_function(self) -> None:
        """This method tests that the function finds a pre-known trim condition.

        :return: None
        """
        ans_range = self.ans_corruption * 2

        v_x_delta = max(abs(self.v_x_ans * ans_range), 0.01)
        alpha_delta = max(abs(self.alpha_ans * ans_range), 0.01)
        beta_delta = max(abs(self.beta_ans * ans_range), 0.01)
        thrust_delta = max(abs(self.thrust_ans * ans_range), 0.01)

        v_x_bounds = (self.v_x_ans - v_x_delta, self.v_x_ans + v_x_delta)
        alpha_bounds = (self.alpha_ans - alpha_delta, self.alpha_ans + alpha_delta)
        beta_bounds = (self.beta_ans - beta_delta, self.beta_ans + beta_delta)
        thrust_bounds = (
            self.thrust_ans - thrust_delta,
            self.thrust_ans + thrust_delta,
        )

        trim_conditions = ps.trim.analyze_steady_trim(
            problem=self.steady_validation_problem,
            solver_type="steady horseshoe vortex lattice method",
            boundsVCg__E=v_x_bounds,
            alpha_bounds=alpha_bounds,
            beta_bounds=beta_bounds,
            boundsExternalFX_W=thrust_bounds,
            objective_cut_off=0.01,
            num_calls=100,
        )

        assert trim_conditions[0] is not None
        assert trim_conditions[1] is not None
        assert trim_conditions[3] is not None
        v_error = (trim_conditions[0] - self.v_x_ans) / self.v_x_ans
        alpha_error = (trim_conditions[1] - self.alpha_ans) / self.alpha_ans
        thrust_error = (trim_conditions[3] - self.thrust_ans) / self.thrust_ans

        allowable_error = 0.1

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(abs(v_error) < allowable_error)
        self.assertTrue(abs(alpha_error) < allowable_error)
        self.assertTrue(abs(thrust_error) < allowable_error)
