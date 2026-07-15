"""This module contains testing cases for the steady convergence function."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _convergence_cache, _convergence_meshing
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
            rtol=0.05,
            atol=0.001,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsNone(converged_parameters[2])

    def test_rejects_exploded_wing(self):
        """This method tests that the function rejects a SteadyProblem whose Airplane has
        an exploded Wing, which carries no edge curves and so cannot be refined.

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
                rtol=0.05,
                atol=0.001,
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
            rtol=0.05,
            atol=0.001,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsNone(converged_parameters[2])

    def test_steady_horseshoe_convergence_resolves_solver(self):
        """This method tests that the function returns the converged, run solver for a
        SteadyHorseshoeVortexLatticeMethodSolver when resolve_converged_solver is True.

        :return: None
        """
        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=self.steady_validation_problem,
            solver_type="steady horseshoe vortex lattice method",
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 4),
            rtol=0.05,
            atol=0.001,
            resolve_converged_solver=True,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]
        converged_solver = converged_parameters[2]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsInstance(
            converged_solver,
            ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver,
        )
        self.assertIsNotNone(converged_solver.airplanes[0].forceCoefficients_W)

    def test_steady_ring_convergence_resolves_solver(self):
        """This method tests that the function returns the converged, run solver for a
        SteadyRingVortexLatticeMethodSolver when resolve_converged_solver is True.

        :return: None
        """
        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=self.steady_validation_problem,
            solver_type="steady ring vortex lattice method",
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 4),
            rtol=0.05,
            atol=0.001,
            resolve_converged_solver=True,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]
        converged_solver = converged_parameters[2]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsInstance(
            converged_solver,
            ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver,
        )
        self.assertIsNotNone(converged_solver.airplanes[0].forceCoefficients_W)

    def test_edge_defined_steady_convergence(self):
        """This method tests that the function finds pre-known convergence parameters for
        a SteadyProblem whose Airplane has an edge-defined Wing, refined by resampling
        its stored edge curves.

        :return: None
        """
        edge_defined_steady_problem = (
            problem_fixtures.make_edge_defined_steady_validation_problem()
        )

        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=edge_defined_steady_problem,
            solver_type="steady ring vortex lattice method",
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 4),
            rtol=0.05,
            atol=0.001,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 1

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsNone(converged_parameters[2])

    def test_steady_cache_reproduces_converged_parameters(self) -> None:
        """This method tests that a run with a cache path finds the same pre-known
        convergence parameters as an uncached run and writes a populated cache file.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

            converged_parameters = ps.convergence.analyze_steady_convergence(
                ref_problem=self.steady_validation_problem,
                solver_type="steady ring vortex lattice method",
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                rtol=0.05,
                atol=0.001,
                cache_path=cache_path,
            )

            self.assertTrue(cache_path.exists())
            self.assertGreater(len(_convergence_cache.load_solve_cache(cache_path)), 0)

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 2

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsNone(converged_parameters[2])

    def test_cache_path_without_json_suffix_raises(self) -> None:
        """This method tests that a cache_path not ending in .json raises a ValueError.

        :return: None
        """
        with self.assertRaises(ValueError):
            ps.convergence.analyze_steady_convergence(
                ref_problem=self.steady_validation_problem,
                solver_type="steady ring vortex lattice method",
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                rtol=0.05,
                atol=0.001,
                cache_path="cache.txt",
            )

    def test_cache_path_directory_raises(self) -> None:
        """This method tests that a cache_path pointing at an existing directory raises
        a ValueError.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"
            cache_path.mkdir()

            with self.assertRaises(ValueError):
                ps.convergence.analyze_steady_convergence(
                    ref_problem=self.steady_validation_problem,
                    solver_type="steady ring vortex lattice method",
                    panel_aspect_ratio_bounds=(4, 2),
                    num_chordwise_panels_bounds=(1, 4),
                    rtol=0.05,
                    atol=0.001,
                    cache_path=cache_path,
                )

    def test_steady_cache_warm_run_skips_solves(self) -> None:
        """This method tests that a second run against a warm cache reuses the stored
        solves and does not run the solver again.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

            cold_parameters = ps.convergence.analyze_steady_convergence(
                ref_problem=self.steady_validation_problem,
                solver_type="steady ring vortex lattice method",
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                rtol=0.05,
                atol=0.001,
                cache_path=cache_path,
            )

            # On the warm run every mesh should be a cache hit, so the solver must never
            # run. Patching run to raise turns any solve into a test failure.
            solver_class = (
                ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
            )
            with mock.patch.object(
                solver_class,
                "run",
                side_effect=AssertionError("The solver ran despite a warm cache."),
            ):
                warm_parameters = ps.convergence.analyze_steady_convergence(
                    ref_problem=self.steady_validation_problem,
                    solver_type="steady ring vortex lattice method",
                    panel_aspect_ratio_bounds=(4, 2),
                    num_chordwise_panels_bounds=(1, 4),
                    rtol=0.05,
                    atol=0.001,
                    cache_path=cache_path,
                )

        self.assertEqual(warm_parameters[0], cold_parameters[0])
        self.assertEqual(warm_parameters[1], cold_parameters[1])

    def test_steady_cache_warm_run_skips_panel_resolution(self) -> None:
        """This method tests that a second run against a warm cache reuses each mesh's
        stored spanwise Panel counts and does not re-resolve them.

        :return: None
        """
        with tempfile.TemporaryDirectory() as tmp:
            cache_path = Path(tmp) / "cache.json"

            cold_parameters = ps.convergence.analyze_steady_convergence(
                ref_problem=self.steady_validation_problem,
                solver_type="steady ring vortex lattice method",
                panel_aspect_ratio_bounds=(4, 2),
                num_chordwise_panels_bounds=(1, 4),
                rtol=0.05,
                atol=0.001,
                cache_path=cache_path,
            )

            # On the warm run every mesh's spanwise Panel counts should be cache hits,
            # so the resolver must never recompute one. Patching the computation to
            # raise turns any recomputation into a test failure.
            with mock.patch.object(
                _convergence_meshing,
                "_get_wing_section_num_spanwise_panels",
                side_effect=AssertionError(
                    "The spanwise Panel resolver ran despite a warm cache."
                ),
            ):
                warm_parameters = ps.convergence.analyze_steady_convergence(
                    ref_problem=self.steady_validation_problem,
                    solver_type="steady ring vortex lattice method",
                    panel_aspect_ratio_bounds=(4, 2),
                    num_chordwise_panels_bounds=(1, 4),
                    rtol=0.05,
                    atol=0.001,
                    cache_path=cache_path,
                )

        self.assertEqual(warm_parameters[0], cold_parameters[0])
        self.assertEqual(warm_parameters[1], cold_parameters[1])

    def test_mixed_airplane_steady_convergence(self):
        """This method tests that the function finds pre-known convergence parameters for
        a SteadyProblem whose Airplane holds both a trapezoidal Wing and an edge-defined
        Wing, refining each Wing by its own spanwise mesh.

        :return: None
        """
        mixed_steady_problem = problem_fixtures.make_mixed_steady_validation_problem()

        converged_parameters = ps.convergence.analyze_steady_convergence(
            ref_problem=mixed_steady_problem,
            solver_type="steady ring vortex lattice method",
            panel_aspect_ratio_bounds=(4, 2),
            num_chordwise_panels_bounds=(1, 4),
            rtol=0.05,
            atol=0.001,
        )

        converged_panel_ar = converged_parameters[0]
        converged_num_chordwise = converged_parameters[1]

        panel_ar_ans = 4
        num_chordwise_ans = 1

        self.assertEqual(converged_panel_ar, panel_ar_ans)
        self.assertEqual(converged_num_chordwise, num_chordwise_ans)
        self.assertIsNone(converged_parameters[2])
