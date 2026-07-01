"""This module contains a testing case for the unsteady convergence function."""

import logging
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _convergence_cache
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

    def test_variable_geometry_unsteady_convergence(self) -> None:
        """This method tests that a variable-geometry UnsteadyProblem runs through
        convergence to completion, exercising the signed cycle-mean final-coefficient
        path that static-geometry cases never reach.

        The bounds are kept minimal because driving a variable-geometry case to true
        convergence would take too long for the test suite, so this run is not expected
        to converge. Its purpose is to exercise that coefficient gathering end to end
        without error.

        :return: None
        """
        variable_geometry_problem = (
            problem_fixtures.make_unsteady_validation_problem_with_variable_geometry()
        )

        # The coarse bounds that keep this smoke test cheap also make the delta_time
        # optimizer warn about poor temporal resolution. That warning is correct in
        # production but is expected noise here, so quiet this one logger for the run.
        movement_logger = logging.getLogger("pterasoftware.movements.movement")
        previous_level = movement_logger.level
        movement_logger.setLevel(logging.ERROR)
        try:
            converged_parameters = ps.convergence.analyze_unsteady_convergence(
                ref_problem=variable_geometry_problem,
                prescribed_wake=True,
                free_wake=False,
                num_cycles_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 3),
                num_chordwise_panels_bounds=(1, 2),
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
            )
        finally:
            movement_logger.setLevel(previous_level)

        self.assertEqual(converged_parameters, (None, None, None, None, None))

    def test_unsteady_cache_reproduces_converged_parameters(self) -> None:
        """This method tests that a run with a cache path finds the same pre-known
        convergence parameters as an uncached run and writes a populated cache file.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

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
                cache_path=cache_path,
            )

            self.assertTrue(cache_path.exists())
            self.assertGreater(len(_convergence_cache.load_solve_cache(cache_path)), 0)

        self.assertEqual(converged_parameters[0], True)
        self.assertEqual(converged_parameters[1], 2)
        self.assertEqual(converged_parameters[2], 4)
        self.assertEqual(converged_parameters[3], 3)
        self.assertIsNone(converged_parameters[4])

    def test_unsteady_cache_warm_run_skips_solves(self) -> None:
        """This method tests that a second run against a warm cache reuses the stored
        solves and does not run the solver again.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

            cold_parameters = ps.convergence.analyze_unsteady_convergence(
                ref_problem=self.unsteady_validation_problem,
                prescribed_wake=True,
                free_wake=False,
                num_chords_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 3),
                num_chordwise_panels_bounds=(1, 2),
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
                cache_path=cache_path,
            )

            # On the warm run every mesh should be a cache hit, so the solver must
            # never run. Patching run to raise turns any solve into a test failure.
            solver_class = (
                ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
            )
            with mock.patch.object(
                solver_class,
                "run",
                side_effect=AssertionError("The solver ran despite a warm cache."),
            ):
                warm_parameters = ps.convergence.analyze_unsteady_convergence(
                    ref_problem=self.unsteady_validation_problem,
                    prescribed_wake=True,
                    free_wake=False,
                    num_chords_bounds=(1, 2),
                    panel_aspect_ratio_bounds=(4, 3),
                    num_chordwise_panels_bounds=(1, 2),
                    rtol=0.05,
                    atol=0.001,
                    show_solver_progress=False,
                    cache_path=cache_path,
                )

        self.assertEqual(warm_parameters[0], cold_parameters[0])
        self.assertEqual(warm_parameters[1], cold_parameters[1])
        self.assertEqual(warm_parameters[2], cold_parameters[2])
        self.assertEqual(warm_parameters[3], cold_parameters[3])

    def test_unsteady_cache_warm_run_skips_delta_time_optimization(self) -> None:
        """This method tests that a second run against a warm cache reuses each mesh's
        stored delta_time and does not re-run the iterative delta_time optimizer.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

            cold_parameters = ps.convergence.analyze_unsteady_convergence(
                ref_problem=self.unsteady_validation_problem,
                prescribed_wake=True,
                free_wake=False,
                num_chords_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 3),
                num_chordwise_panels_bounds=(1, 2),
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
                cache_path=cache_path,
            )

            # On the warm run every mesh's delta_time should be a cache hit, so the
            # iterative optimizer must never run. Patching it to raise turns any
            # optimization into a test failure.
            with mock.patch.object(
                ps.movements.movement,
                "_optimize_delta_time",
                side_effect=AssertionError(
                    "The delta_time optimizer ran despite a warm cache."
                ),
            ):
                warm_parameters = ps.convergence.analyze_unsteady_convergence(
                    ref_problem=self.unsteady_validation_problem,
                    prescribed_wake=True,
                    free_wake=False,
                    num_chords_bounds=(1, 2),
                    panel_aspect_ratio_bounds=(4, 3),
                    num_chordwise_panels_bounds=(1, 2),
                    rtol=0.05,
                    atol=0.001,
                    show_solver_progress=False,
                    cache_path=cache_path,
                )

        self.assertEqual(warm_parameters[0], cold_parameters[0])
        self.assertEqual(warm_parameters[1], cold_parameters[1])
        self.assertEqual(warm_parameters[2], cold_parameters[2])
        self.assertEqual(warm_parameters[3], cold_parameters[3])

    def test_unsteady_cache_reuses_memos_across_different_bounds(self) -> None:
        """This method tests that memos written by one run are reused by a later run
        whose bounds differ, reproducing the uncached result while never re-running the
        delta_time optimizer for the meshes the earlier run already resolved.

        :return: None
        """
        # An uncached run over the narrower bounds gives the reference result.
        uncached_parameters = ps.convergence.analyze_unsteady_convergence(
            ref_problem=self.unsteady_validation_problem,
            prescribed_wake=True,
            free_wake=False,
            num_chords_bounds=(1, 2),
            panel_aspect_ratio_bounds=(3, 2),
            num_chordwise_panels_bounds=(1, 2),
            rtol=0.05,
            atol=0.001,
            show_solver_progress=False,
        )

        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

            # A cold run over wider Panel aspect ratio bounds populates the cache with
            # every mesh the narrower run below will need. The memos are keyed on absolute
            # mesh values, so they carry over even though the two runs index their sweeps
            # differently.
            ps.convergence.analyze_unsteady_convergence(
                ref_problem=self.unsteady_validation_problem,
                prescribed_wake=True,
                free_wake=False,
                num_chords_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 2),
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
                cache_path=cache_path,
            )

            # The narrower run reuses those cached memos, so its delta_time optimizer must
            # never run despite the different bounds. Patching it to raise turns any
            # optimization into a test failure.
            with mock.patch.object(
                ps.movements.movement,
                "_optimize_delta_time",
                side_effect=AssertionError(
                    "The delta_time optimizer ran despite a warm cache."
                ),
            ):
                warm_parameters = ps.convergence.analyze_unsteady_convergence(
                    ref_problem=self.unsteady_validation_problem,
                    prescribed_wake=True,
                    free_wake=False,
                    num_chords_bounds=(1, 2),
                    panel_aspect_ratio_bounds=(3, 2),
                    num_chordwise_panels_bounds=(1, 2),
                    rtol=0.05,
                    atol=0.001,
                    show_solver_progress=False,
                    cache_path=cache_path,
                )

        # Reusing the cross-bounds memos must reproduce the uncached result exactly.
        self.assertEqual(warm_parameters[0], uncached_parameters[0])
        self.assertEqual(warm_parameters[1], uncached_parameters[1])
        self.assertEqual(warm_parameters[2], uncached_parameters[2])
        self.assertEqual(warm_parameters[3], uncached_parameters[3])

    def test_unsteady_cache_path_without_json_suffix_raises(self) -> None:
        """This method tests that a cache_path not ending in .json raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.convergence.analyze_unsteady_convergence(
                ref_problem=self.unsteady_validation_problem,
                prescribed_wake=True,
                free_wake=False,
                num_chords_bounds=(1, 2),
                panel_aspect_ratio_bounds=(4, 3),
                num_chordwise_panels_bounds=(1, 2),
                rtol=0.05,
                atol=0.001,
                show_solver_progress=False,
                cache_path="cache.txt",
            )

    def test_unsteady_cache_path_directory_raises(self) -> None:
        """This method tests that a cache_path pointing at an existing directory raises
        a ValueError.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"
            cache_path.mkdir()

            with self.assertRaises(ValueError):
                ps.convergence.analyze_unsteady_convergence(
                    ref_problem=self.unsteady_validation_problem,
                    prescribed_wake=True,
                    free_wake=False,
                    num_chords_bounds=(1, 2),
                    panel_aspect_ratio_bounds=(4, 3),
                    num_chordwise_panels_bounds=(1, 2),
                    rtol=0.05,
                    atol=0.001,
                    show_solver_progress=False,
                    cache_path=cache_path,
                )

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
