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
            rtol=0.05,
            atol=0.001,
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
        self.assertIsNone(converged_parameters[4])

    def test_rejects_exploded_wing(self):
        """This method tests that the function rejects an UnsteadyProblem whose Airplane
        has an exploded Wing, which carries no edge curves and so cannot be refined.

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
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
            )

    def test_edge_defined_unsteady_convergence(self):
        """This method tests that the function finds pre-known convergence parameters for
        an UnsteadyProblem whose Airplane has an edge-defined Wing with static
        WingCrossSectionMovements.

        :return: None
        """
        edge_defined_unsteady_problem = (
            problem_fixtures.make_edge_defined_unsteady_validation_problem()
        )

        # The Panel aspect ratio and chordwise Panel bounds are kept tight because an
        # edge-defined Wing refined to a fine Panel aspect ratio with many chordwise
        # Panels needs many WingCrossSections, which makes the unsteady solves expensive.
        # These bounds still span enough of each parameter direction to detect
        # convergence.
        converged_parameters = ps.convergence.analyze_unsteady_convergence(
            ref_problem=edge_defined_unsteady_problem,
            prescribed_wake=True,
            free_wake=True,
            num_chords_bounds=(1, 4),
            panel_aspect_ratio_bounds=(4, 3),
            num_chordwise_panels_bounds=(1, 3),
            rtol=0.05,
            atol=0.001,
            show_solver_progress=False,
        )

        converged_wake_state = converged_parameters[0]
        converged_num_chords = converged_parameters[1]
        converged_panel_ar = converged_parameters[2]
        converged_num_chordwise = converged_parameters[3]

        wake_state_ans = True
        num_chords_ans = 3
        panel_ar_ans = 4
        num_chordwise_ans = 1

        self.assertEqual(converged_wake_state, wake_state_ans)
        self.assertEqual(converged_num_chords, num_chords_ans)
        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsNone(converged_parameters[4])

    def test_rejects_non_static_edge_defined_movement(self):
        """This method tests that the function rejects an UnsteadyProblem whose
        edge-defined Wing carries a non-static WingCrossSectionMovement, which resampling
        the Wing cannot preserve.

        A non-static WingCrossSectionMovement makes the Movement variable, so num_cycles
        bounds are supplied.

        :return: None
        """
        edge_defined_non_static_problem = (
            problem_fixtures.make_edge_defined_non_static_unsteady_validation_problem()
        )

        with self.assertRaisesRegex(ValueError, "static"):
            ps.convergence.analyze_unsteady_convergence(
                ref_problem=edge_defined_non_static_problem,
                prescribed_wake=True,
                free_wake=True,
                num_cycles_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
            )

    def test_unsteady_convergence_resolves_solver(self):
        """This method tests that the function returns the converged, run solver for an
        UnsteadyRingVortexLatticeMethodSolver when resolve_converged_solver is True.

        :return: None
        """
        converged_parameters = ps.convergence.analyze_unsteady_convergence(
            ref_problem=self.unsteady_validation_problem,
            prescribed_wake=True,
            free_wake=True,
            num_chords_bounds=(1, 4),
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 5),
            rtol=0.05,
            atol=0.001,
            show_solver_progress=False,
            resolve_converged_solver=True,
        )

        converged_wake_state = converged_parameters[0]
        converged_num_chords = converged_parameters[1]
        converged_panel_ar = converged_parameters[2]
        converged_num_chordwise = converged_parameters[3]
        converged_solver = converged_parameters[4]

        wake_state_ans = True
        num_chords_ans = 2
        panel_ar_ans = 4
        num_chordwise_ans = 3

        self.assertEqual(converged_wake_state, wake_state_ans)
        self.assertEqual(converged_num_chords, num_chords_ans)
        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsInstance(
            converged_solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )
        self.assertGreater(
            len(converged_solver.unsteady_problem.finalForceCoefficients_W), 0
        )
