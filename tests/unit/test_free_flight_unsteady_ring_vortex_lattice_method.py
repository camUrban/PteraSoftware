"""This module contains classes to test the
FreeFlightUnsteadyRingVortexLatticeMethodSolver class."""

import unittest

import numpy as np

import pterasoftware as ps
from pterasoftware.free_flight_unsteady_ring_vortex_lattice_method import (
    FreeFlightUnsteadyRingVortexLatticeMethodSolver,
)
from tests.unit.fixtures import (
    operating_point_fixtures,
    problem_fixtures,
    solver_fixtures,
)


class TestFreeFlightUnsteadyRingVortexLatticeMethodSolver(unittest.TestCase):
    """This is a class with functions to test
    FreeFlightUnsteadyRingVortexLatticeMethodSolvers.
    """

    def setUp(self):
        """Set up a fresh solver for each test."""
        self.solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()

    def test_initialization_accepts_free_flight_unsteady_problem(self):
        """Test that initialization accepts a FreeFlightUnsteadyProblem."""
        self.assertIsInstance(
            self.solver, FreeFlightUnsteadyRingVortexLatticeMethodSolver
        )
        self.assertIsInstance(
            self.solver,
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
        )
        self.assertIsInstance(
            self.solver.unsteady_problem, ps.problems.FreeFlightUnsteadyProblem
        )

    def test_initialization_rejects_non_free_flight_problem(self):
        """Test that initialization raises TypeError for a coupled problem that is not a
        FreeFlightUnsteadyProblem.
        """
        coupled_problem = problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        with self.assertRaises(TypeError):
            FreeFlightUnsteadyRingVortexLatticeMethodSolver(coupled_problem)

    def test_initialization_rejects_non_problem_types(self):
        """Test that initialization raises TypeError for non-problem inputs."""
        invalid_inputs = [None, "not_a_problem", 123, [1, 2, 3], {"key": "value"}]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    FreeFlightUnsteadyRingVortexLatticeMethodSolver(invalid)

    def test_free_flight_unsteady_problem_property_narrows_unsteady_problem(self):
        """Test that the _free_flight_unsteady_problem property returns the same object
        as unsteady_problem, narrowed to FreeFlightUnsteadyProblem.
        """
        self.assertIs(
            self.solver._free_flight_unsteady_problem, self.solver.unsteady_problem
        )
        self.assertIsInstance(
            self.solver._free_flight_unsteady_problem,
            ps.problems.FreeFlightUnsteadyProblem,
        )

    def test_models_body_rates_is_true(self):
        """Test that the free-flight solver declares that it models body rates, so the
        inherited constructor permits a non zero omegas_BP1__E.
        """
        self.assertTrue(self.solver._models_body_rates)

    def test_permits_non_zero_body_rates(self):
        """Test that the free-flight solver constructs from a problem whose initial
        OperatingPoint carries a non zero body angular velocity, which the base solver
        would reject.
        """
        rotating_problem = (
            problem_fixtures.make_with_body_rates_free_flight_unsteady_problem_fixture()
        )
        solver = FreeFlightUnsteadyRingVortexLatticeMethodSolver(rotating_problem)
        self.assertIsInstance(solver, FreeFlightUnsteadyRingVortexLatticeMethodSolver)

    def test_current_omegas_rad_without_rotation_is_zero(self):
        """Test that _currentOmegasRad_GP1__E returns a zero vector when the current
        OperatingPoint carries no body rotation.
        """
        omegasRad_GP1__E = self.solver._currentOmegasRad_GP1__E()
        np.testing.assert_array_equal(omegasRad_GP1__E, np.zeros(3))

    def test_current_omegas_rad_transforms_body_rate(self):
        """Test that _currentOmegasRad_GP1__E transforms the current OperatingPoint's
        body rate from the body axes in degrees per second to the geometry axes in
        radians per second.

        The body-to-geometry transformation negates the x and z components, and the
        result is converted from degrees per second to radians per second.
        """
        self.solver.current_operating_point = (
            operating_point_fixtures.make_with_body_rates_operating_point_fixture()
        )

        omegasRad_GP1__E = self.solver._currentOmegasRad_GP1__E()

        # The fixture's body rate is (0, 0, 1) degrees per second in body axes, which
        # becomes (0, 0, -1) degrees per second in geometry axes, then is converted to
        # radians per second.
        expected_omegasRad_GP1__E = np.deg2rad(np.array([0.0, 0.0, -1.0]))
        np.testing.assert_allclose(omegasRad_GP1__E, expected_omegasRad_GP1__E)

    def test_current_omegas_rad_negates_x_and_z_preserves_y(self):
        """Test that _currentOmegasRad_GP1__E negates the x and z components and
        preserves the y component when transforming the body rate from the body axes to
        the geometry axes.

        Using a body rate with three distinct non zero components pins down the full
        body to geometry flip, which a body rate aligned with a single axis cannot.
        """
        self.solver.current_operating_point = (
            operating_point_fixtures.make_with_full_body_rates_operating_point_fixture()
        )

        omegasRad_GP1__E = self.solver._currentOmegasRad_GP1__E()

        # The fixture's body rate is (1, 2, 3) degrees per second in body axes, which
        # becomes (-1, 2, -3) degrees per second in geometry axes, then is converted to
        # radians per second.
        expected_omegasRad_GP1__E = np.deg2rad(np.array([-1.0, 2.0, -3.0]))
        np.testing.assert_allclose(omegasRad_GP1__E, expected_omegasRad_GP1__E)

    def test_declares_substep_slots(self):
        """Test that the subclass declares exactly the strongly coupled sub-iteration's
        transient slots and so does not gain an instance __dict__ that would defeat the
        parent's __slots__.
        """
        self.assertEqual(
            FreeFlightUnsteadyRingVortexLatticeMethodSolver.__slots__,
            (
                "_substep_next_step",
                "_substep_next_steady_problem",
                "_substep_next_operating_point",
                "_substep_stackVIndGridWrvp_GP1__E",
                "_substep_gamma_n",
                "_substep_gamma_n_minus_1",
            ),
        )
        with self.assertRaises(AttributeError):
            self.solver.not_a_real_attribute = 42


if __name__ == "__main__":
    unittest.main()
