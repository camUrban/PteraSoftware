"""This module contains a testing case for the unsteady convergence function."""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import (
    airplane_fixtures,
    operating_point_fixtures,
    problem_fixtures,
)


class TestUnsteadyConvergence(unittest.TestCase):
    """This is a class for testing the unsteady convergence function."""

    def setUp(self):
        """This method sets up the test.

        :return: None
        """
        self.unsteady_validation_problem = (
            problem_fixtures.make_unsteady_validation_problem_with_static_geometry()
        )

    def test_unsteady_convergence(self):
        """This method tests that the function finds pre-known convergence parameters
        for an UnsteadyRingVortexLatticeMethodSolver.

        :return: None
        """
        converged_parameters = ps.convergence.analyze_unsteady_convergence(
            ref_problem=self.unsteady_validation_problem,
            prescribed_wake=True,
            free_wake=True,
            num_chords_bounds=(1, 4),
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 5),
            convergence_criteria=5.0,
            show_solver_progress=False,
        )

        converged_wake_state = converged_parameters[0]
        converged_num_chords = converged_parameters[1]
        converged_panel_ar = converged_parameters[2]
        converged_num_chordwise = converged_parameters[3]

        wake_state_ans = True
        num_chords_ans = 2
        panel_ar_ans = 4
        num_chordwise_ans = 3

        self.assertEqual(converged_wake_state, wake_state_ans)
        self.assertEqual(converged_num_chords, num_chords_ans)
        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)

    def test_rejects_exploded_wing(self):
        """This method tests that the function rejects an UnsteadyProblem whose Airplane
        has a Wing with a non-trapezoidal spanwise mesh.

        :return: None
        """
        exploded_airplane = airplane_fixtures.make_exploded_validation_airplane()
        exploded_wing = exploded_airplane.wings[0]

        wing_cross_section_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wing_cross_section
            )
            for wing_cross_section in exploded_wing.wing_cross_sections
        ]
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=exploded_wing,
            wing_cross_section_movements=wing_cross_section_movements,
        )
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=exploded_airplane,
            wing_movements=[wing_movement],
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=(
                    operating_point_fixtures.make_validation_operating_point()
                )
            )
        )
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            num_chords=1,
        )
        exploded_problem = ps.problems.UnsteadyProblem(movement=movement)

        with self.assertRaises(ValueError):
            ps.convergence.analyze_unsteady_convergence(
                ref_problem=exploded_problem,
                prescribed_wake=True,
                free_wake=True,
                num_chords_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                convergence_criteria=5.0,
                show_solver_progress=False,
            )
