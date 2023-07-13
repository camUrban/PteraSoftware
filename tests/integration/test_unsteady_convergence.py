"""This module contains a testing case for the unsteady convergence function.

This module contains the following classes:
    TestUnsteadyConvergence: This is a class for testing the unsteady convergence
    function.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""
import unittest

import pterasoftware as ps
from tests.integration.fixtures import problem_fixtures


class TestUnsteadyConvergence(unittest.TestCase):
    """This is a class for testing the unsteady convergence function.

    This class contains the following public methods:
        setUp: This method sets up the test.

        tearDown: This method tears down the test.

        test_unsteady_convergence: This method tests that the function finds
        pre-known convergence parameters for an unsteady problem.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
        """This method sets up the test.

        :return: None
        """

        # Create the unsteady problem.
        self.unsteady_validation_problem = (
            problem_fixtures.make_unsteady_validation_problem_with_static_geometry()
        )

    def tearDown(self):
        """This method tears down the test.

        :return: None
        """

        del self.unsteady_validation_problem

    def test_unsteady_convergence(self):
        """This method tests that the function finds pre-known convergence parameters
        for an unsteady problem.

        :return: None
        """

        converged_parameters = ps.convergence.analyze_unsteady_convergence(
            ref_problem=self.unsteady_validation_problem,
            prescribed_wake=True,
            free_wake=True,
            num_chords_bounds=(2, 6),
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(2, 6),
            convergence_criteria=5.0,
            coefficient_mask=[True, False, True, False, True, False],
        )

        converged_wake_state = converged_parameters[0]
        converged_num_chords = converged_parameters[1]
        converged_panel_ar = converged_parameters[2]
        converged_num_chordwise = converged_parameters[3]

        # ToDo: Update these expected results.
        wake_state_ans = True
        num_chords_ans = 4
        panel_ar_ans = 4
        num_chordwise_ans = 4

        self.assertTrue(converged_wake_state == wake_state_ans)
        self.assertTrue(abs(converged_num_chords - num_chords_ans) <= 1)
        self.assertTrue(abs(converged_panel_ar - panel_ar_ans) <= 1)
        self.assertTrue(abs(converged_num_chordwise - num_chordwise_ans) <= 1)
