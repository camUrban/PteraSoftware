"""This module tests the FreeFlightUnsteadyRingVortexLatticeMethodSolver by running two
coupled free flight cases and checking that the resulting motion is physically sensible.

The first case is the simple glider in an unpowered glide. The glider is statically
stable in pitch and yaw (its geometry, center of gravity, and inertia were tuned for this
and verified in XFLR5), so a successful coupled run produces a bounded, damped glide: the
angle of attack oscillates gently about its trimmed value rather than diverging, lift
balances weight, the glider descends and moves forward, the flight stays laterally
symmetric, and total mechanical energy decreases as drag dissipates it.

The second case is a flapping-wing airframe whose main wing flaps symmetrically. Flapping
produces large, oscillatory aerodynamic loads with significant added mass, the regime
that an explicit, loosely coupled handoff renders unstable and diverges. The strongly
coupled, Aitken-relaxed sub-iteration drives the loads and the body state to within-step
consistency each time step instead, so the run stays finite and bounded. This case
validates that the coupling completes a flapping free flight without diverging, that the
flapping genuinely excites large oscillatory loads, that the sub-iteration converges
within its cap every step, and that symmetric flapping keeps the flight laterally
symmetric.
"""

import logging
import unittest

import numpy as np

from tests.integration.fixtures import solver_fixtures


class TestFreeFlightUnsteadyRingVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the
    FreeFlightUnsteadyRingVortexLatticeMethodSolver on a statically stable glider."""

    @classmethod
    def setUpClass(cls):
        """This method sets up the test by running the simple glider's free flight
        simulation once and extracting its time histories.

        :return: None
        """
        cls.solver = solver_fixtures.make_simple_glider_free_flight_solver()
        cls.solver.run(prescribed_wake=True, show_progress=False)

        problem = cls.solver.unsteady_problem
        operating_points = problem.movement.operating_point_movement.operating_points

        # Extract per-time-step time histories from the dynamically populated
        # OperatingPoints and the problem's recorded aerodynamic load history.
        cls.speeds = np.array(
            [operating_point.vCg__E for operating_point in operating_points]
        )
        cls.alphas = np.array(
            [operating_point.alpha for operating_point in operating_points]
        )
        cls.betas = np.array(
            [operating_point.beta for operating_point in operating_points]
        )
        cls.positions_E_Eo = np.array(
            [operating_point.CgP1_E_Eo for operating_point in operating_points]
        )
        cls.omegas_BP1__E = np.array(
            [operating_point.omegas_BP1__E for operating_point in operating_points]
        )

        # In Ptera Software's wind axes, lift is the negative z-component of the wind-axes
        # force.
        forces_W = np.array(problem.forces_W)
        cls.lifts = -forces_W[:, 2]
        cls.side_forces = forces_W[:, 1]

        cls.weight = cls.solver.current_airplanes[0].weight
        cls.g_E = operating_points[0].g_E
        cls.mass = problem.mass
        cls.inertia_matrix = problem.I_BP1_CgP1

    def test_run_produces_finite_full_length_history(self):
        """This method tests that the coupled run completed for every time step and
        produced finite state and load histories (no divergence or NaNs).

        :return: None
        """
        num_steps = self.solver.unsteady_problem.num_steps

        self.assertEqual(len(self.speeds), num_steps)
        self.assertEqual(len(self.lifts), num_steps)

        self.assertTrue(np.all(np.isfinite(self.speeds)))
        self.assertTrue(np.all(np.isfinite(self.alphas)))
        self.assertTrue(np.all(np.isfinite(self.positions_E_Eo)))
        self.assertTrue(np.all(np.isfinite(self.lifts)))

    def test_lift_approximately_balances_weight(self):
        """This method tests that the mean lift over the glide is close to the glider's
        weight, as it should be near the trimmed condition.

        :return: None
        """
        mean_lift = float(np.mean(self.lifts))
        lift_to_weight_ratio = mean_lift / self.weight

        self.assertGreater(lift_to_weight_ratio, 0.85)
        self.assertLess(lift_to_weight_ratio, 1.15)

    def test_glider_descends(self):
        """This method tests that the glider loses altitude over the glide. In Earth
        axes, the +z direction points down (gravity acts along +z), so a descending
        glider's z-coordinate increases.

        :return: None
        """
        initial_z_E = self.positions_E_Eo[0, 2]
        final_z_E = self.positions_E_Eo[-1, 2]

        self.assertGreater(final_z_E, initial_z_E)

    def test_speed_remains_bounded(self):
        """This method tests that the glide speed stays within a reasonable band of its
        initial value rather than running away.

        :return: None
        """
        initial_speed = self.speeds[0]

        self.assertLess(float(np.max(self.speeds)), 1.2 * initial_speed)
        self.assertGreater(float(np.min(self.speeds)), 0.8 * initial_speed)

    def test_flight_stays_laterally_symmetric(self):
        """This method tests that the glide stays in the longitudinal plane. The glider
        is laterally symmetric and starts with zero sideslip, so the sideslip angle, the
        lateral side force, and the lateral position should all stay negligible.

        :return: None
        """
        self.assertLess(float(np.max(np.abs(self.betas))), 1.0)
        self.assertLess(float(np.max(np.abs(self.side_forces))), 1.0)
        self.assertLess(float(np.max(np.abs(self.positions_E_Eo[:, 1]))), 1.0e-3)

    def test_angle_of_attack_stays_bounded(self):
        """This method tests that the angle of attack stays in a narrow band about its
        trimmed value rather than diverging. This is the signature of the glider's static
        pitch stability: a stable airframe makes a small restoring oscillation, whereas
        an unstable one would run away toward the angle limits.

        :return: None
        """
        trimmed_alpha = self.alphas[0]

        self.assertLess(float(np.max(self.alphas)), trimmed_alpha + 3.0)
        self.assertGreater(float(np.min(self.alphas)), trimmed_alpha - 3.0)

    def test_total_mechanical_energy_decreases(self):
        """This method tests that the glider's total mechanical energy at the end of the
        glide is less than at the start. With no thrust, drag is the only non-conservative
        load, so it must dissipate mechanical energy over the glide.

        :return: None
        """
        # Translational kinetic energy.
        kinetic_energies = 0.5 * self.mass * self.speeds**2

        # Gravitational potential energy. With g_E along +z (down), a descending glider
        # (increasing z) loses potential energy.
        potential_energies = -self.mass * (self.positions_E_Eo @ self.g_E)

        # Rotational kinetic energy, using the body-axes angular velocity (converted from
        # degrees per second to radians per second) and inertia matrix.
        omegasRad_BP1__E = np.deg2rad(self.omegas_BP1__E)
        rotational_kinetic_energies = np.array(
            [
                0.5 * omegaRad @ self.inertia_matrix @ omegaRad
                for omegaRad in omegasRad_BP1__E
            ]
        )

        total_energies = (
            kinetic_energies + potential_energies + rotational_kinetic_energies
        )

        self.assertLess(total_energies[-1], total_energies[0])


class TestFreeFlightUnsteadyRingVortexLatticeMethodFlapping(unittest.TestCase):
    """This is a class for testing the
    FreeFlightUnsteadyRingVortexLatticeMethodSolver on a flapping-wing airframe whose
    strongly coupled sub-iteration must remain stable under large oscillatory loads."""

    _warning_records: list[logging.LogRecord]

    @classmethod
    def setUpClass(cls):
        """This method sets up the test by running the flapping-wing free flight
        simulation once and extracting its time histories.

        It attaches a handler to the solver's logger for the duration of the run so the
        test can later assert that the strongly coupled sub-iteration never hit its
        iteration cap (the solver logs a warning on a capped, non-converged step).

        :return: None
        """
        # Collect the WARNING-level records the solver emits during the run so the test
        # can check for the sub-iteration cap warning.
        cls._warning_records = []

        class _WarningCollector(logging.Handler):
            def emit(self, record: logging.LogRecord) -> None:
                cls._warning_records.append(record)

        solver_logger = logging.getLogger("pterasoftware.problems")
        collector = _WarningCollector(level=logging.WARNING)
        solver_logger.addHandler(collector)

        try:
            cls.solver = solver_fixtures.make_flapping_free_flight_solver()
            cls.solver.run(prescribed_wake=True, show_progress=False)
        finally:
            solver_logger.removeHandler(collector)

        problem = cls.solver.unsteady_problem
        operating_points = problem.movement.operating_point_movement.operating_points

        # Extract per-time-step time histories from the dynamically populated
        # OperatingPoints and the problem's recorded aerodynamic load history.
        cls.betas = np.array(
            [operating_point.beta for operating_point in operating_points]
        )
        cls.positions_E_Eo = np.array(
            [operating_point.CgP1_E_Eo for operating_point in operating_points]
        )

        # In Ptera Software's wind axes, lift is the negative z-component of the wind-axes
        # force and the lateral side force is the y-component.
        forces_W = np.array(problem.forces_W)
        cls.forces_W = forces_W
        cls.lifts = -forces_W[:, 2]
        cls.side_forces = forces_W[:, 1]

        cls.weight = cls.solver.current_airplanes[0].weight

    def test_run_produces_finite_full_length_history(self):
        """This method tests that the coupled run completed for every time step and
        produced finite state and load histories. Under loose coupling the flapping case
        diverges; a finite, full-length history is the signature that the strongly
        coupled sub-iteration kept the integration stable.

        :return: None
        """
        num_steps = self.solver.unsteady_problem.num_steps

        self.assertEqual(len(self.forces_W), num_steps)

        self.assertTrue(np.all(np.isfinite(self.forces_W)))
        self.assertTrue(np.all(np.isfinite(self.betas)))
        self.assertTrue(np.all(np.isfinite(self.positions_E_Eo)))

    def test_sub_iteration_converges_within_cap(self):
        """This method tests that the strongly coupled sub-iteration converged within its
        iteration cap on every free flight time step. The solver logs a warning when a
        step reaches the cap without converging and accepts the capped iterate anyway, so
        the absence of that warning confirms every step reached the convergence tolerance.

        :return: None
        """
        cap_warnings = [
            record for record in self._warning_records if "cap" in record.getMessage()
        ]

        self.assertEqual(cap_warnings, [])

    def test_flapping_excites_oscillatory_loads(self):
        """This method tests that the flapping genuinely drives large, oscillatory
        aerodynamic loads rather than settling to a near-constant value. The peak-to-peak
        lift over the run should exceed the airframe's weight, confirming that the case
        exercises the coupling in the large-load regime it is meant to stress.

        :return: None
        """
        peak_to_peak_lift = float(np.max(self.lifts) - np.min(self.lifts))

        self.assertGreater(peak_to_peak_lift, self.weight)

    def test_flight_stays_laterally_symmetric(self):
        """This method tests that the flight stays in the longitudinal plane. The main
        wing flaps symmetrically and the airframe starts with zero sideslip, so the
        sideslip angle, the lateral side force, and the lateral position should all stay
        negligible despite the large flapping loads.

        :return: None
        """
        self.assertLess(float(np.max(np.abs(self.betas))), 1.0)
        self.assertLess(float(np.max(np.abs(self.side_forces))), 1.0)
        self.assertLess(float(np.max(np.abs(self.positions_E_Eo[:, 1]))), 1.0e-3)
