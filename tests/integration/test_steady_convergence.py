# REFACTOR: I haven't yet started refactoring this module.
"""This module contains testing cases for the steady convergence function."""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import problem_fixtures


class TestSteadyConvergence(unittest.TestCase):
    """This is a class for testing the steady convergence function."""

    def setUp(self):
        """This method sets up the test.

        :return: None
        """

        # Create the steady problem.
        self.steady_validation_problem = (
            problem_fixtures.make_steady_validation_problem()
        )

    def test_steady_horseshoe_convergence(self):
        """This method tests that the function finds pre-known convergence parameters
        for a horseshoe vortex lattice method solver.

        :return: None
        """

        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=self.steady_validation_problem,
            solver_type="steady horseshoe vortex lattice method",
            panel_aspect_ratio_bounds=(4, 1),
            num_chordwise_panels_bounds=(3, 10),
            convergence_criteria=1.0,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 4

        self.assertTrue(abs(converged_panel_ar - panel_ar_ans) <= 1)
        self.assertTrue(abs(converged_num_chordwise - num_chordwise_ans) <= 1)

    def test_steady_ring_convergence(self):
        """This method tests that the function finds pre-known convergence parameters
        for a ring vortex lattice method solver.

        :return: None
        """

        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=self.steady_validation_problem,
            solver_type="steady ring vortex lattice method",
            panel_aspect_ratio_bounds=(4, 1),
            num_chordwise_panels_bounds=(3, 10),
            convergence_criteria=1.0,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 5

        self.assertTrue(abs(converged_panel_ar - panel_ar_ans) <= 1)
        self.assertTrue(abs(converged_num_chordwise - num_chordwise_ans) <= 1)
