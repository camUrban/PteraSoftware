# NOTE: I haven't yet started refactoring this module.
"""This module contains a testing case for the steady trim function.

This module contains the following classes:
    TestSteadyTrimHorseshoeVortexLatticeMethod: This is a class for testing the
    steady trim function on with the horseshoe vortex lattice method solver.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import airplane_fixtures


class TestSteadyTrimHorseshoeVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the steady trim function on with the horseshoe
    vortex lattice method solver.

    This class contains the following public methods:
        setUp: This method sets up the test.

        tearDown: This method tears down the test.

        test_function: This method tests that the function finds a pre-known trim
        condition.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
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
        this_operating_point = ps.operating_point.OperatingPoint(vCg__E=corrupted_v_x,
                                                                 alpha=corrupted_alpha,
                                                                 beta=corrupted_beta,
                                                                 externalFX_W=corrupted_thrust)

        # Create the steady problem.
        self.steady_validation_problem = ps.problems.SteadyProblem(
            airplanes=[this_airplane],
            operating_point=this_operating_point,
        )

        del this_airplane
        del this_operating_point

    def tearDown(self):
        """This method tears down the test.

        :return: None
        """

        del self.steady_validation_problem
        del self.v_x_ans
        del self.alpha_ans
        del self.beta_ans
        del self.thrust_ans
        del self.ans_corruption

    def test_function(self):
        """This method tests that the function finds a pre-known trim condition.

        :return: None
        """

        ans_range = self.ans_corruption * 2

        v_x_delta = max(abs(self.v_x_ans * ans_range), 0.01)
        alpha_delta = max(abs(self.alpha_ans * ans_range), 0.01)
        beta_delta = max(abs(self.beta_ans * ans_range), 0.01)
        thrust_delta = max(abs(self.thrust_ans * ans_range), 0.01)

        v_x_bounds = (min((self.v_x_ans - v_x_delta), 0), self.v_x_ans + v_x_delta)
        alpha_bounds = (self.alpha_ans - alpha_delta, self.alpha_ans + alpha_delta)
        beta_bounds = (self.beta_ans - beta_delta, self.beta_ans + beta_delta)
        thrust_bounds = (
            min((self.thrust_ans - thrust_delta), 0),
            self.thrust_ans + thrust_delta,
        )

        trim_conditions = ps.trim.analyze_steady_trim(
            problem=self.steady_validation_problem,
            velocity_bounds=v_x_bounds,
            alpha_bounds=alpha_bounds,
            beta_bounds=beta_bounds,
            external_thrust_bounds=thrust_bounds,
            objective_cut_off=0.01,
            num_calls=100,
        )

        v_error = (trim_conditions[0] - self.v_x_ans) / self.v_x_ans
        alpha_error = (trim_conditions[1] - self.alpha_ans) / self.alpha_ans
        thrust_error = (trim_conditions[3] - self.thrust_ans) / self.thrust_ans

        allowable_error = 0.1

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(abs(v_error) < allowable_error)
        self.assertTrue(abs(alpha_error) < allowable_error)
        self.assertTrue(abs(thrust_error) < allowable_error)
