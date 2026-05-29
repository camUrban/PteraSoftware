"""This module contains classes to test the
AeroelasticUnsteadyRingVortexLatticeMethodSolver class."""

import unittest

import numpy as np

import pterasoftware as ps
from pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method import (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver,
)
from tests.unit.fixtures import problem_fixtures, solver_fixtures


class TestAeroelasticUnsteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This class contains unit tests for the
    AeroelasticUnsteadyRingVortexLatticeMethodSolver class."""

    def setUp(self):
        """Set up test fixtures for solver tests."""
        self.solver = solver_fixtures.make_aeroelastic_unsteady_ring_solver_fixture()

    def test_initialization_accepts_aeroelastic_unsteady_problem(self):
        """Test that initialization succeeds with an AeroelasticUnsteadyProblem."""
        self.assertIsInstance(
            self.solver, AeroelasticUnsteadyRingVortexLatticeMethodSolver
        )

    def test_initialization_rejects_non_aeroelastic_problem(self):
        """Test that a non-AeroelasticUnsteadyProblem raises TypeError."""
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        with self.assertRaises(TypeError):
            AeroelasticUnsteadyRingVortexLatticeMethodSolver(coupled_problem)

    def test_unsteady_problem_is_aeroelastic(self):
        """Test that the stored unsteady_problem is an AeroelasticUnsteadyProblem."""
        self.assertIsInstance(
            self.solver.unsteady_problem,
            ps.problems.AeroelasticUnsteadyProblem,
        )

    def test_aeroelastic_problem_property_narrows_unsteady_problem(self):
        """Test that the _aeroelastic_problem property returns the same object as
        unsteady_problem, narrowed to AeroelasticUnsteadyProblem."""
        self.assertIs(self.solver._aeroelastic_problem, self.solver.unsteady_problem)
        self.assertIsInstance(
            self.solver._aeroelastic_problem,
            ps.problems.AeroelasticUnsteadyProblem,
        )

    def test_num_panels_positive(self):
        """Test that num_panels is a positive integer."""
        self.assertIsInstance(self.solver.num_panels, int)
        self.assertGreater(self.solver.num_panels, 0)

    def test_slep_point_indices_is_ndarray(self):
        """Test that slep_point_indices is a numpy ndarray."""
        self.assertIsInstance(self.solver.slep_point_indices, np.ndarray)

    def test_slep_point_indices_starts_at_zero(self):
        """Test that slep_point_indices[0] is zero (first panel is always the
        first SLEP)."""
        self.assertEqual(self.solver.slep_point_indices[0], 0)

    def test_slep_point_indices_non_decreasing(self):
        """Test that slep_point_indices is non-decreasing."""
        indices = self.solver.slep_point_indices
        self.assertTrue(np.all(np.diff(indices) >= 0))

    def test_slep_point_indices_integer_dtype(self):
        """Test that slep_point_indices has integer dtype."""
        self.assertTrue(np.issubdtype(self.solver.slep_point_indices.dtype, np.integer))

    def test_moments_gp1_slep_initial_empty(self):
        """Test that moments_GP1_Slep is initially empty before any step."""
        self.assertEqual(self.solver.moments_GP1_Slep.size, 0)

    def test_stack_leading_edge_points_initial_empty(self):
        """Test that stack_leading_edge_points is initially empty before any step."""
        self.assertEqual(self.solver.stack_leading_edge_points.size, 0)

    def test_reinitialize_step_arrays_hook_produces_zero_arrays(self):
        """Test that _reinitialize_step_arrays_hook fills all SLEP arrays with
        zeros."""
        self.solver._reinitialize_step_arrays_hook()
        arrays_to_check = [
            self.solver.stackCblvpr_GP1_Slep,
            self.solver.stackCblvpf_GP1_Slep,
            self.solver.stackCblvpl_GP1_Slep,
            self.solver.stackCblvpb_GP1_Slep,
            self.solver.stackCpp_GP1_Slep,
            self.solver.moments_GP1_Slep,
            self.solver.stackFlpp_GP1_CgP1,
        ]
        for arr in arrays_to_check:
            with self.subTest(array=arr):
                self.assertTrue(np.all(arr == 0.0))

    def test_reinitialize_step_arrays_hook_produces_correct_shape(self):
        """Test that _reinitialize_step_arrays_hook creates arrays of shape
        (num_panels, 3)."""
        self.solver._reinitialize_step_arrays_hook()
        expected_shape = (self.solver.num_panels, 3)
        arrays_to_check = [
            self.solver.stackCblvpr_GP1_Slep,
            self.solver.stackCblvpf_GP1_Slep,
            self.solver.stackCblvpl_GP1_Slep,
            self.solver.stackCblvpb_GP1_Slep,
            self.solver.stackCpp_GP1_Slep,
            self.solver.moments_GP1_Slep,
            self.solver.stackFlpp_GP1_CgP1,
        ]
        for arr in arrays_to_check:
            with self.subTest(array=arr):
                self.assertEqual(arr.shape, expected_shape)

    def test_slep_point_indices_length_equals_wcs_count(self):
        """Test that slep_point_indices has one entry per WingCrossSection."""
        first_problem = self.solver.unsteady_problem.steady_problems[0]
        total_wcs = sum(
            len(wing.wing_cross_sections)
            for airplane in first_problem.airplanes
            for wing in airplane.wings
        )
        self.assertEqual(len(self.solver.slep_point_indices), total_wcs)
