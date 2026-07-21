"""This module contains classes to test the
FreeFlightUnsteadyRingVortexLatticeMethodSolver class."""

import unittest
from typing import Any
from unittest.mock import MagicMock, patch

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
    FreeFlightUnsteadyRingVortexLatticeMethodSolvers."""

    def setUp(self) -> None:
        """Set up a fresh solver for each test."""
        self.solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()

    def test_initialization_accepts_free_flight_unsteady_problem(self) -> None:
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

    def test_initialization_rejects_non_free_flight_problem(self) -> None:
        """Test that initialization raises TypeError for a coupled problem that is not a
        FreeFlightUnsteadyProblem."""
        coupled_problem: Any = (
            problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
        )
        with self.assertRaises(TypeError):
            FreeFlightUnsteadyRingVortexLatticeMethodSolver(coupled_problem)

    def test_initialization_rejects_non_problem_types(self) -> None:
        """Test that initialization raises TypeError for non-problem inputs."""
        invalid_inputs: list[Any] = [
            None,
            "not_a_problem",
            123,
            [1, 2, 3],
            {"key": "value"},
        ]
        for invalid in invalid_inputs:
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    FreeFlightUnsteadyRingVortexLatticeMethodSolver(invalid)

    def test_free_flight_unsteady_problem_property_narrows_unsteady_problem(
        self,
    ) -> None:
        """Test that the _free_flight_unsteady_problem property returns the same object
        as unsteady_problem, narrowed to FreeFlightUnsteadyProblem."""
        self.assertIs(
            self.solver._free_flight_unsteady_problem, self.solver.unsteady_problem
        )
        self.assertIsInstance(
            self.solver._free_flight_unsteady_problem,
            ps.problems.FreeFlightUnsteadyProblem,
        )

    def test_models_body_rates_is_true(self) -> None:
        """Test that the free-flight solver declares that it models body rates, so the
        inherited constructor permits a non zero omegas_BP1__E."""
        self.assertTrue(self.solver._models_body_rates)

    def test_permits_non_zero_body_rates(self) -> None:
        """Test that the free-flight solver constructs from a problem whose initial
        OperatingPoint carries a non zero body angular velocity, which the base solver
        would reject."""
        rotating_problem = (
            problem_fixtures.make_with_body_rates_free_flight_unsteady_problem_fixture()
        )
        solver = FreeFlightUnsteadyRingVortexLatticeMethodSolver(rotating_problem)
        self.assertIsInstance(solver, FreeFlightUnsteadyRingVortexLatticeMethodSolver)

    def test_current_omegas_without_rotation_is_zero(self) -> None:
        """Test that _currentOmegasRad_GP1__E returns a zero vector when the current
        OperatingPoint carries no body rotation."""
        omegasRad_GP1__E = self.solver._currentOmegasRad_GP1__E()
        np.testing.assert_array_equal(omegasRad_GP1__E, np.zeros(3))

    def test_current_omegas_transforms_body_rate(self) -> None:
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

    def test_current_omegas_negates_x_and_z_preserves_y(self) -> None:
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

    def test_declares_substep_slots(self) -> None:
        """Test that the subclass declares exactly the strongly coupled sub-iteration's
        transient slots and so does not gain an instance __dict__ that would defeat the
        parent's __slots__."""
        self.assertEqual(
            FreeFlightUnsteadyRingVortexLatticeMethodSolver.__slots__,
            (
                "_substep_next_step",
                "_substep_next_steady_problem",
                "_substep_next_operating_point",
                "_substepStackVIndGridWrvp_GP1__E",
                "_substep_gamma_n",
                "_substep_gamma_n_minus_1",
            ),
        )
        with self.assertRaises(AttributeError):
            setattr(self.solver, "not_a_real_attribute", 42)


class TestFreeFlightSolverSubstepDispatch(unittest.TestCase):
    """Tests for the substep-aware geometry and operating point dispatch overrides."""

    def setUp(self) -> None:
        """Set up a fresh solver for each dispatch test."""
        self.solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()

    def test_get_steady_problem_at_dispatches_on_substep_state(self) -> None:
        """Test that _get_steady_problem_at returns the transient next-step
        SteadyProblem only for the next step during a sub-iteration, and otherwise
        defers to the committed accessor."""
        committed = self.solver._get_steady_problem_at(0)
        transient = MagicMock()

        self.solver._substep_next_step = 1
        self.solver._substep_next_steady_problem = transient

        # The transient problem is returned for the matching next step only.
        self.assertIs(self.solver._get_steady_problem_at(1), transient)
        # A non-matching step falls through to the committed accessor.
        self.assertIs(self.solver._get_steady_problem_at(0), committed)

        # With no transient problem set, the committed accessor is used.
        self.solver._substep_next_steady_problem = None
        self.assertIs(self.solver._get_steady_problem_at(0), committed)

    def test_operating_point_at_dispatches_on_substep_state(self) -> None:
        """Test that _operating_point_at returns the trial OperatingPoint only for the
        next step during a sub-iteration, and otherwise defers to the committed
        accessor."""
        committed_operating_point = self.solver._operating_point_at(0)
        trial_operating_point = MagicMock()

        self.solver._substep_next_step = 1
        self.solver._substep_next_operating_point = trial_operating_point

        # The trial operating point is returned for the matching next step only.
        self.assertIs(self.solver._operating_point_at(1), trial_operating_point)
        # A non-matching step falls through to the committed accessor.
        self.assertIs(self.solver._operating_point_at(0), committed_operating_point)

        # With no trial operating point set, the committed accessor is used.
        self.solver._substep_next_operating_point = None
        self.assertIs(self.solver._operating_point_at(0), committed_operating_point)


class TestFreeFlightSolverSubstepLifecycle(unittest.TestCase):
    """Tests for the transient working state freeze_substep and restore_substep
    manage."""

    def setUp(self) -> None:
        """Set up a fresh solver for each lifecycle test."""
        self.solver = solver_fixtures.make_free_flight_unsteady_ring_solver_fixture()

    def test_freeze_substep_snapshots_strengths_as_copies(self) -> None:
        """Test that freeze_substep records the next step, stores the transient
        SteadyProblem, and snapshots the current and previous bound vortex strengths as
        independent copies.

        The snapshots must be copies, because the trials overwrite the solver's strength
        arrays while evaluating the next step, and a view would carry that mutation into
        the snapshot.
        """
        self.solver._current_step = 0
        self.solver._current_bound_vortex_strengths = np.array(
            [1.0, 2.0, 3.0], dtype=float
        )
        self.solver._last_bound_vortex_strengths = np.array(
            [4.0, 5.0, 6.0], dtype=float
        )
        next_steady_problem = MagicMock()

        self.solver.freeze_substep(next_steady_problem)

        self.assertEqual(self.solver._substep_next_step, 1)
        self.assertIs(self.solver._substep_next_steady_problem, next_steady_problem)
        self.assertIsNone(self.solver._substep_next_operating_point)
        np.testing.assert_array_equal(
            self.solver._substep_gamma_n, np.array([1.0, 2.0, 3.0])
        )
        np.testing.assert_array_equal(
            self.solver._substep_gamma_n_minus_1, np.array([4.0, 5.0, 6.0])
        )

        # Mutating the source strengths must not change the snapshots.
        self.solver._current_bound_vortex_strengths[0] = 999.0
        self.solver._last_bound_vortex_strengths[0] = 999.0
        assert self.solver._substep_gamma_n is not None
        assert self.solver._substep_gamma_n_minus_1 is not None
        self.assertEqual(self.solver._substep_gamma_n[0], 1.0)
        self.assertEqual(self.solver._substep_gamma_n_minus_1[0], 4.0)

    def test_freeze_substep_skips_induced_precompute_for_prescribed_wake(self) -> None:
        """Test that freeze_substep leaves the frozen wake induced velocities unset for
        a prescribed wake, which has no induced transport to precompute."""
        self.solver._prescribed_wake = True
        self.solver._current_step = 0
        self.solver._current_bound_vortex_strengths = np.ones(3, dtype=float)
        self.solver._last_bound_vortex_strengths = np.zeros(3, dtype=float)

        self.solver.freeze_substep(MagicMock())

        self.assertIsNone(self.solver._substepStackVIndGridWrvp_GP1__E)

    def test_restore_substep_clears_transient_state(self) -> None:
        """Test that restore_substep re-evaluates the current step and clears every
        transient sub-iteration slot.

        A leftover non-None slot would redirect the next step's geometry and wake reads
        to a stale scratch copy, so the clearing is a correctness contract. The re-
        evaluation of the current step's aerodynamics is stubbed, since reconstructing
        it is the inherited solver's behavior, not this method's.
        """
        self.solver._substep_next_step = 1
        self.solver._substep_next_steady_problem = MagicMock()
        self.solver._substep_next_operating_point = MagicMock()
        self.solver._substepStackVIndGridWrvp_GP1__E = [[np.zeros(3, dtype=float)]]
        self.solver._substep_gamma_n = np.ones(3, dtype=float)
        self.solver._substep_gamma_n_minus_1 = np.zeros(3, dtype=float)

        with patch.object(
            ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
            "_evaluate_step_aerodynamics",
        ) as mock_evaluate_step_aerodynamics:
            self.solver.restore_substep(step=0)

        mock_evaluate_step_aerodynamics.assert_called_once()
        self.assertIsNone(self.solver._substep_next_step)
        self.assertIsNone(self.solver._substep_next_steady_problem)
        self.assertIsNone(self.solver._substep_next_operating_point)
        self.assertIsNone(self.solver._substepStackVIndGridWrvp_GP1__E)
        self.assertIsNone(self.solver._substep_gamma_n)
        self.assertIsNone(self.solver._substep_gamma_n_minus_1)
