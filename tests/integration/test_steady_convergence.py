"""This module contains testing cases for the steady convergence function."""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import (
    airplane_fixtures,
    operating_point_fixtures,
    problem_fixtures,
)


class TestSteadyConvergence(unittest.TestCase):
    """This is a class for testing the steady convergence function."""

    def setUp(self):
        """This method sets up the test.

        :return: None
        """
        self.steady_validation_problem = (
            problem_fixtures.make_steady_validation_problem()
        )

    def test_steady_horseshoe_convergence(self):
        """This method tests that the function finds pre-known convergence parameters
        for a SteadyHorseshoeVortexLatticeMethodSolver.

        :return: None
        """
        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=self.steady_validation_problem,
            solver_type="steady horseshoe vortex lattice method",
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 4),
            convergence_criteria=5.0,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)

    def test_rejects_exploded_wing(self):
        """This method tests that the function rejects a SteadyProblem whose Airplane has
        a Wing with a non-trapezoidal spanwise mesh.

        :return: None
        """
        exploded_problem = ps.problems.SteadyProblem(
            airplanes=[airplane_fixtures.make_exploded_validation_airplane()],
            operating_point=operating_point_fixtures.make_validation_operating_point(),
        )

        with self.assertRaises(ValueError):
            ps.convergence.analyze_steady_convergence(
                ref_problem=exploded_problem,
                solver_type="steady ring vortex lattice method",
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                convergence_criteria=5.0,
            )

    def test_steady_ring_convergence(self):
        """This method tests that the function finds pre-known convergence parameters
        for a SteadyRingVortexLatticeMethodSolver.

        :return: None
        """
        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=self.steady_validation_problem,
            solver_type="steady ring vortex lattice method",
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 4),
            convergence_criteria=5.0,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
