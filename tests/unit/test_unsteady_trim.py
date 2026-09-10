"""This module contains a class to test the analyze_unsteady_trim function."""

import unittest
from typing import Any
from unittest.mock import patch

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import movement_fixtures, problem_fixtures


class TestAnalyzeUnsteadyTrim(unittest.TestCase):
    """A class with functions to test analyze_unsteady_trim."""

    problem: ps.problems.UnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all analyze_unsteady_trim tests."""
        cls.problem = problem_fixtures.make_basic_unsteady_problem_fixture()

    def test_problem_validation(self) -> None:
        """Test problem parameter validation."""
        bad_problem: Any = "not a problem"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=bad_problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_multiple_airplane_movements_validation(self) -> None:
        """Test that only one AirplaneMovement is allowed."""
        problem = problem_fixtures.make_multi_airplane_unsteady_problem_fixture()

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_boundsVCg__E_validation(self) -> None:
        """Test boundsVCg__E parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=bad_str,
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=bad_short_tuple,
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=bad_element_tuple,
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(10.0, 1.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(0.0, 10.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_alpha_bounds_validation(self) -> None:
        """Test alpha_bounds parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=bad_str,
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=bad_short_tuple,
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=bad_element_tuple,
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(10.0, -10.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        # The bounds must lie within the range OperatingPoint accepts for alpha.
        for bad_alpha_bounds in [(-180.0, 20.0), (-20.0, 180.001)]:
            with self.subTest(alpha_bounds=bad_alpha_bounds):
                with self.assertRaisesRegex(ValueError, "range \\(-180.0, 180.0\\]"):
                    ps.trim.analyze_unsteady_trim(
                        problem=self.problem,
                        boundsVCg__E=(1.0, 100.0),
                        alpha_bounds=bad_alpha_bounds,
                        beta_bounds=(-20.0, 20.0),
                        boundsExternalFX_W=(-1000.0, 1000.0),
                    )

    def test_beta_bounds_validation(self) -> None:
        """Test beta_bounds parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=bad_str,
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=bad_short_tuple,
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=bad_element_tuple,
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(10.0, -10.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

        # The bounds must exclude the poles at beta = +/-90.0, where OperatingPoint only
        # accepts alpha = 0.0, since the search varies alpha and beta independently.
        for bad_beta_bounds in [(-90.0, 20.0), (-20.0, 90.0), (-180.0, 180.0)]:
            with self.subTest(beta_bounds=bad_beta_bounds):
                with self.assertRaisesRegex(ValueError, "range \\(-90.0, 90.0\\)"):
                    ps.trim.analyze_unsteady_trim(
                        problem=self.problem,
                        boundsVCg__E=(1.0, 100.0),
                        alpha_bounds=(-20.0, 20.0),
                        beta_bounds=bad_beta_bounds,
                        boundsExternalFX_W=(-1000.0, 1000.0),
                    )

    def test_boundsExternalFX_W_validation(self) -> None:
        """Test boundsExternalFX_W parameter validation."""
        bad_str: Any = "invalid"
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=bad_str,
            )

        bad_short_tuple: Any = (1.0,)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=bad_short_tuple,
            )

        bad_element_tuple: Any = ("a", 10.0)
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=bad_element_tuple,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(10.0, -10.0),
            )

    def test_objective_cut_off_validation(self) -> None:
        """Test objective_cut_off parameter validation."""
        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=0.0,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=-1.0,
            )

    def test_num_calls_validation(self) -> None:
        """Test num_calls parameter validation."""
        bad_num_calls: Any = 1.5
        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=0,
            )

        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=bad_num_calls,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                num_calls=-1,
            )

    def test_force_method_validation(self) -> None:
        """Test force_method parameter validation."""
        bad_type: Any = 42
        with self.assertRaises(TypeError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                force_method=bad_type,
            )

        with self.assertRaises(ValueError):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                force_method="invalid",
            )

    def test_base_operating_point_parameter_bounds_validation(self) -> None:
        """Test that the base operating point values must lie within the supplied
        bounds."""
        base_operating_point = (
            self.problem.movement.operating_point_movement.base_operating_point
        )

        with self.subTest(parameter="vCg__E"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(
                        base_operating_point.vCg__E + 1.0,
                        base_operating_point.vCg__E + 10.0,
                    ),
                    alpha_bounds=(-20.0, 20.0),
                    beta_bounds=(-20.0, 20.0),
                    boundsExternalFX_W=(-1000.0, 1000.0),
                )

        with self.subTest(parameter="alpha"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(1.0, 100.0),
                    alpha_bounds=(
                        base_operating_point.alpha + 1.0,
                        base_operating_point.alpha + 10.0,
                    ),
                    beta_bounds=(-20.0, 20.0),
                    boundsExternalFX_W=(-1000.0, 1000.0),
                )

        with self.subTest(parameter="beta"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(1.0, 100.0),
                    alpha_bounds=(-20.0, 20.0),
                    beta_bounds=(
                        base_operating_point.beta + 1.0,
                        base_operating_point.beta + 10.0,
                    ),
                    boundsExternalFX_W=(-1000.0, 1000.0),
                )

        with self.subTest(parameter="externalFX_W"):
            with self.assertRaises(ValueError):
                ps.trim.analyze_unsteady_trim(
                    problem=self.problem,
                    boundsVCg__E=(1.0, 100.0),
                    alpha_bounds=(-20.0, 20.0),
                    beta_bounds=(-20.0, 20.0),
                    boundsExternalFX_W=(
                        base_operating_point.externalFX_W + 1.0,
                        base_operating_point.externalFX_W + 10.0,
                    ),
                )

    def test_requires_unresolved_duration(self) -> None:
        """Test that trim requires a step-size-independent duration."""
        movement = (
            movement_fixtures.make_non_static_movement_with_explicit_num_steps_fixture()
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)

        with self.assertRaisesRegex(
            ValueError, "must define its duration with num_cycles or num_chords"
        ):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_requires_unresolved_wake_length(self) -> None:
        """Test that trim requires a step-size-independent wake length."""
        reference_movement = self.problem.movement
        movement = ps.movements.movement.Movement(
            airplane_movements=list(reference_movement.airplane_movements),
            operating_point_movement=reference_movement.operating_point_movement,
            delta_time=reference_movement.delta_time,
            num_cycles=reference_movement.num_cycles,
            max_wake_rows=2,
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)

        with self.assertRaisesRegex(
            ValueError, "must define its maximum wake length with max_wake_chords"
        ):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_requires_speed_bounds_above_movement_amplitude(self) -> None:
        """Test that each trial's oscillating speed remains positive."""
        reference_movement = self.problem.movement
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=(
                    reference_movement.operating_point_movement.base_operating_point
                ),
                ampVCg__E=1.0,
                periodVCg__E=2.0,
            )
        )
        movement = ps.movements.movement.Movement(
            airplane_movements=list(reference_movement.airplane_movements),
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_cycles=2,
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)

        with self.assertRaisesRegex(ValueError, "must be greater than.*ampVCg__E"):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_g_E_validation(self) -> None:
        """Test that a zero base g_E is rejected, since the trim analysis places the
        Airplane's weight along g_E's direction."""
        with self.assertRaisesRegex(ValueError, "g_E must be non-zero"):
            ps.trim.analyze_unsteady_trim(
                problem=self.problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_base_attitude_validation(self) -> None:
        """Test that a base attitude that does not resolve to level flight is rejected,
        since the trials resolve their own attitudes to level flight and would otherwise
        silently discard it."""
        reference_movement = movement_fixtures.make_static_movement_fixture()
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=ps.operating_point.OperatingPoint(
                    angles_E_to_BP1_izyx=(0.0, 0.0, 0.0), g_E=(0.0, 0.0, 9.80665)
                )
            )
        )
        movement = ps.movements.movement.Movement(
            airplane_movements=list(reference_movement.airplane_movements),
            operating_point_movement=operating_point_movement,
            num_chords=reference_movement.num_chords,
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)

        with self.assertRaisesRegex(ValueError, "must resolve to level flight"):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
            )

    def test_trials_are_level_flight_with_weight_along_g_E(self) -> None:
        """Test that each trial resolves its attitude to level flight and balances the
        Airplane's weight along g_E's direction, rotated into that trial's wind axes.

        The base OperatingPoint has a g_E that is not along Earth +z, so the weight is
        not along the wind z axis even in level flight. The stubbed solver records each
        trial's Earth axes to wind axes matrix and returns exactly the force
        coefficients that cancel the external thrust and the weight placed along g_E, so
        the objective is zero at the initial guess only if the trim analysis places the
        weight the same way. With a tiny cutoff and a single allowed call, the analysis
        then returns the initial guess, and anything else returns Nones.
        """
        reference_movement = movement_fixtures.make_static_movement_fixture()
        base_operating_point = ps.operating_point.OperatingPoint(
            vCg__E=10.0,
            alpha=5.0,
            beta=3.0,
            externalFX_W=7.0,
            g_E=(1.0, -2.0, 4.0),
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point
            )
        )
        movement = ps.movements.movement.Movement(
            airplane_movements=list(reference_movement.airplane_movements),
            operating_point_movement=operating_point_movement,
            num_chords=reference_movement.num_chords,
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)
        base_airplane = movement.airplane_movements[0].base_airplane
        self.assertGreater(base_airplane.weight, 0.0)
        trial_T_pas_E_CgP1_to_W_CgP1s: list[np.ndarray] = []

        class SolverStub:
            """A solver stub that records each trial's Earth axes to wind axes matrix
            and whose loads exactly cancel the external thrust and the weight placed
            along g_E."""

            def __init__(self, unsteady_problem: ps.problems.UnsteadyProblem) -> None:
                self.unsteady_problem = unsteady_problem

            def run(self, **_: Any) -> None:
                trial_operating_point = (
                    self.unsteady_problem.movement.operating_point_movement.base_operating_point
                )
                trial_T_pas_E_CgP1_to_W_CgP1s.append(
                    trial_operating_point.T_pas_E_CgP1_to_W_CgP1
                )
                g_E = trial_operating_point.g_E
                weightForce_E = base_airplane.weight * g_E / np.linalg.norm(g_E)
                R_pas_E_to_W = trial_operating_point.T_pas_E_CgP1_to_W_CgP1[:3, :3]
                externalForces_W = (
                    np.array([trial_operating_point.externalFX_W, 0.0, 0.0])
                    + R_pas_E_to_W @ weightForce_E
                )
                assert base_airplane.s_ref is not None
                externalForceCoefficients_W = (
                    externalForces_W
                    / trial_operating_point.qInf__E
                    / base_airplane.s_ref
                )
                self.unsteady_problem.finalForceCoefficients_W = [
                    -externalForceCoefficients_W
                ]
                self.unsteady_problem.finalMomentCoefficients_W_Cg = [
                    np.zeros(3, dtype=float)
                ]

        with patch(
            "pterasoftware.trim.unsteady_ring_vortex_lattice_method."
            "UnsteadyRingVortexLatticeMethodSolver",
            SolverStub,
        ):
            trim_conditions = ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=1.0e-9,
                num_calls=1,
                show_solver_progress=False,
            )

        self.assertEqual(trim_conditions, (10.0, 5.0, 3.0, 7.0))
        self.assertEqual(len(trial_T_pas_E_CgP1_to_W_CgP1s), 1)
        npt.assert_allclose(trial_T_pas_E_CgP1_to_W_CgP1s[0], np.eye(4), atol=1e-12)

    def test_static_trial_uses_final_load_coefficients(self) -> None:
        """Test that a static trial uses its final-time-step loads."""
        reference_movement = movement_fixtures.make_static_movement_fixture()
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=ps.operating_point.OperatingPoint(
                    g_E=(0.0, 0.0, 9.80665)
                )
            )
        )
        movement = ps.movements.movement.Movement(
            airplane_movements=list(reference_movement.airplane_movements),
            operating_point_movement=operating_point_movement,
            num_chords=reference_movement.num_chords,
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)

        class SolverStub:
            """A solver stub that populates only static-movement loads."""

            def __init__(self, unsteady_problem: ps.problems.UnsteadyProblem) -> None:
                self.unsteady_problem = unsteady_problem

            def run(self, **_: Any) -> None:
                self.unsteady_problem.finalForceCoefficients_W = [
                    np.zeros(3, dtype=float)
                ]
                self.unsteady_problem.finalMomentCoefficients_W_Cg = [
                    np.zeros(3, dtype=float)
                ]

        with patch(
            "pterasoftware.trim.unsteady_ring_vortex_lattice_method."
            "UnsteadyRingVortexLatticeMethodSolver",
            SolverStub,
        ):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(1.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=1_000_000.0,
                show_solver_progress=False,
            )

    def test_trial_preserves_movement_parameters(self) -> None:
        """Test that trim trials preserve unresolved movement parameters."""
        reference_movement = self.problem.movement
        reference_operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=ps.operating_point.OperatingPoint(
                    g_E=(0.0, 0.0, 9.80665)
                ),
                ampVCg__E=1.0,
                periodVCg__E=2.0,
                spacingVCg__E="uniform",
                phaseVCg__E=45.0,
            )
        )
        movement = ps.movements.movement.Movement(
            airplane_movements=list(reference_movement.airplane_movements),
            operating_point_movement=reference_operating_point_movement,
            delta_time=0.1,
            num_cycles=2,
            max_wake_cycles=1,
        )
        problem = ps.problems.UnsteadyProblem(movement=movement)
        trial_movements: list[ps.movements.movement.Movement] = []

        class SolverStub:
            """A solver stub that records the generated trial Movement."""

            def __init__(self, unsteady_problem: ps.problems.UnsteadyProblem) -> None:
                self.unsteady_problem = unsteady_problem
                trial_movements.append(unsteady_problem.movement)

            def run(self, **_: Any) -> None:
                self.unsteady_problem.finalMeanForceCoefficients_W = [
                    np.zeros(3, dtype=float)
                ]
                self.unsteady_problem.finalMeanMomentCoefficients_W_Cg = [
                    np.zeros(3, dtype=float)
                ]

        with patch(
            "pterasoftware.trim.unsteady_ring_vortex_lattice_method."
            "UnsteadyRingVortexLatticeMethodSolver",
            SolverStub,
        ):
            ps.trim.analyze_unsteady_trim(
                problem=problem,
                boundsVCg__E=(2.0, 100.0),
                alpha_bounds=(-20.0, 20.0),
                beta_bounds=(-20.0, 20.0),
                boundsExternalFX_W=(-1000.0, 1000.0),
                objective_cut_off=1_000_000.0,
                show_solver_progress=False,
            )

        trial_movement = trial_movements[0]
        trial_operating_point_movement = trial_movement.operating_point_movement
        self.assertEqual(
            trial_operating_point_movement.ampVCg__E,
            reference_operating_point_movement.ampVCg__E,
        )
        self.assertEqual(
            trial_operating_point_movement.periodVCg__E,
            reference_operating_point_movement.periodVCg__E,
        )
        self.assertEqual(
            trial_operating_point_movement.spacingVCg__E,
            reference_operating_point_movement.spacingVCg__E,
        )
        self.assertEqual(
            trial_operating_point_movement.phaseVCg__E,
            reference_operating_point_movement.phaseVCg__E,
        )
        self.assertEqual(trial_movement.num_cycles, movement.num_cycles)
        self.assertEqual(trial_movement.num_chords, movement.num_chords)
        self.assertEqual(trial_movement.max_wake_cycles, movement.max_wake_cycles)
        self.assertEqual(trial_movement.max_wake_chords, movement.max_wake_chords)
        self.assertEqual(trial_movement.static, movement.static)


if __name__ == "__main__":
    unittest.main()
