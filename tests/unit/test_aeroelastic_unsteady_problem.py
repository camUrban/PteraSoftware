"""This module contains classes to test the AeroelasticUnsteadyProblem class."""

import unittest
from unittest.mock import MagicMock, PropertyMock, patch

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import (
    aeroelastic_wing_movement_fixtures,
    movement_fixtures,
    problem_fixtures,
)


class TestAeroelasticUnsteadyProblem(unittest.TestCase):
    """This class contains unit tests for the AeroelasticUnsteadyProblem class."""

    def setUp(self):
        """Set up test fixtures for AeroelasticUnsteadyProblem tests."""
        self.problem = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )

    def test_initialization_accepts_aeroelastic_movement(self):
        """Test that AeroelasticUnsteadyProblem accepts an AeroelasticMovement."""
        self.assertIsInstance(self.problem, ps.problems.AeroelasticUnsteadyProblem)

    def test_initialization_rejects_non_aeroelastic_movement(self):
        """Test that a non-AeroelasticMovement raises TypeError."""
        basic_movement = movement_fixtures.make_basic_movement_fixture()
        with self.assertRaises(TypeError):
            ps.problems.AeroelasticUnsteadyProblem(
                movement=basic_movement,
                wing_density=0.01,
                spring_constant_rad=10.0,
                damping_constant_rad=0.5,
            )

    def test_num_steps(self):
        """Test that num_steps matches the movement's num_steps."""
        self.assertEqual(self.problem.num_steps, 3)

    def test_delta_time(self):
        """Test that delta_time matches the movement's delta_time."""
        self.assertAlmostEqual(self.problem.delta_time, 0.1)

    def test_initial_steady_problems_count(self):
        """Test that exactly one SteadyProblem exists at initialization (step 0)."""
        self.assertEqual(len(self.problem.steady_problems), 1)

    def test_wing_density_stored(self):
        """Test that wing_density is stored correctly."""
        self.assertAlmostEqual(self.problem.wing_density, 0.01)

    def test_spring_constant_rad_stored(self):
        """Test that spring_constant_rad is stored correctly."""
        self.assertAlmostEqual(self.problem.spring_constant_rad, 10.0)

    def test_damping_constant_rad_stored(self):
        """Test that damping_constant_rad is stored correctly."""
        self.assertAlmostEqual(self.problem.damping_constant_rad, 0.5)

    def test_calculate_mass_matrix_shape(self):
        """Test that calculate_mass_matrix returns an array with shape
        (num_chordwise_panels, num_spanwise_panels, 3)."""
        wing = self.problem.steady_problems[0].airplanes[0].wings[0]
        mass_matrix = self.problem.calculate_mass_matrix(wing)
        self.assertEqual(len(mass_matrix.shape), 3)
        self.assertEqual(mass_matrix.shape[0], wing.num_chordwise_panels)
        self.assertEqual(mass_matrix.shape[1], wing.num_spanwise_panels)
        self.assertEqual(mass_matrix.shape[2], 3)

    def test_calculate_mass_matrix_non_negative(self):
        """Test that all entries in the mass matrix are non-negative."""
        wing = self.problem.steady_problems[0].airplanes[0].wings[0]
        mass_matrix = self.problem.calculate_mass_matrix(wing)
        self.assertTrue(np.all(mass_matrix >= 0.0))

    def test_calculate_mass_matrix_components_equal(self):
        """Test that all three spatial components of the mass matrix are equal."""
        wing = self.problem.steady_problems[0].airplanes[0].wings[0]
        mass_matrix = self.problem.calculate_mass_matrix(wing)
        np.testing.assert_array_equal(mass_matrix[:, :, 0], mass_matrix[:, :, 1])
        np.testing.assert_array_equal(mass_matrix[:, :, 1], mass_matrix[:, :, 2])

    def test_calculate_spring_moments_accumulates_state_across_steps(self):
        """Test that a strip's spring-damper ODE is re-seeded from its own state at
        the end of the previous time step, so a constant aerodynamic moment drives
        the strip's deformation angle y component to accumulate across steps toward
        the static equilibrium M / k.

        Uses a single-strip Wing in a slow, overdamped regime where the one-step
        response from rest is a small fraction of M / k. If a strip's state were
        discarded between steps (a restart from rest every step), the deformation
        angle y component would freeze at that one-step response instead of
        converging.
        """
        problem = ps.problems.AeroelasticUnsteadyProblem(
            movement=movement_fixtures.make_basic_aeroelastic_movement_fixture(),
            wing_density=0.01,
            spring_constant_rad=10.0,
            damping_constant_rad=20.0,
        )
        wing = problem.steady_problems[0].airplanes[0].wings[0]
        num_spanwise_panels = wing.num_spanwise_panels
        num_chordwise_panels = wing.num_chordwise_panels
        assert num_spanwise_panels == 1
        static_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_static_aeroelastic_wing_movement_fixture()
        )

        # Choose the per-entry mass so the strip's ODE inertia, I = (1 / 2) * mass *
        # L^2 with mass summed over the strip's mass-matrix entries, is exactly 1.0.
        # With spring_constant_rad = 10.0 and damping_constant_rad = 20.0, the ODE is
        # then overdamped (c^2 > 4 * k * I) with a slowest decay rate of about 0.5 per
        # second, so it is far from settled within one 0.1 s time step.
        L = (wing.wing_cross_sections[0].chord + wing.wing_cross_sections[1].chord) / 2
        num_mass_entries = num_chordwise_panels * num_spanwise_panels * 3
        mass_matrix = np.full(
            (num_chordwise_panels, num_spanwise_panels, 3),
            2.0 / (num_mass_entries * L**2),
        )

        # Apply a constant aerodynamic moment y component (in the first Airplane's
        # geometry axes, relative to the strip's leading edge point), summing to 2.0
        # N*m over the strip's chordwise panels, so the static deformation angle y
        # component is M / k = 0.2 rad.
        aeroMoments_GP1_Slep = np.zeros(
            (num_chordwise_panels, num_spanwise_panels, 3), dtype=float
        )
        aeroMoments_GP1_Slep[:, 0, 1] = 2.0 / num_chordwise_panels
        static_theta_rad = 2.0 / problem.spring_constant_rad

        # March the structural solve, replicating _record_structural_state's
        # post-discard recording after each time step.
        theta_history_rad = []
        for step in range(problem.step_discards + 1, problem.step_discards + 61):
            (
                newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz,
                newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz,
            ) = problem.calculate_spring_moments(
                num_spanwise_panels=num_spanwise_panels,
                wing=wing,
                mass_matrix=mass_matrix,
                aeroMoments_GP1_Slep=aeroMoments_GP1_Slep,
                step=step,
                wing_idx=0,
                wing_movement=static_wing_movement,
            )
            problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[0].append(
                np.array(newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz, dtype=float)
            )
            problem._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[0].append(
                np.array(
                    newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz, dtype=float
                )
            )
            theta_history_rad.append(
                float(newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[1])
            )

        # The regime check: the one-step response must be well below the static
        # equilibrium, or freezing and converging would be indistinguishable.
        self.assertLess(theta_history_rad[0], 0.5 * static_theta_rad)
        # The deformation angle y component must accumulate monotonically, as an
        # overdamped rise from rest.
        for earlier, later in zip(theta_history_rad[:10], theta_history_rad[1:11]):
            self.assertGreater(later, earlier)
        # The deformation angle y component must approach the static equilibrium.
        self.assertAlmostEqual(
            theta_history_rad[-1], static_theta_rad, delta=0.1 * static_theta_rad
        )

    def test_record_structural_state_records_new_state_after_discard_window(self):
        """Test that _record_structural_state appends the new state to both time
        series as fresh arrays when the time step is past the discard window."""
        problem = self.problem
        wing_idx = 0
        num_wing_cross_sections = problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[
            wing_idx
        ][-1].shape[0]
        newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz = np.full(
            num_wing_cross_sections, 0.25, dtype=float
        )
        newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = np.full(
            num_wing_cross_sections, -0.5, dtype=float
        )

        problem._record_structural_state(
            step=problem.step_discards + 1,
            newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz=(
                newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz
            ),
            newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz=(
                newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz
            ),
            wing_idx=wing_idx,
        )

        deformationAnglesYRad_Wcsp_to_Wcs_ixyz = (
            problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = (
            problem._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        self.assertEqual(len(deformationAnglesYRad_Wcsp_to_Wcs_ixyz), 2)
        self.assertEqual(len(deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz), 2)
        np.testing.assert_array_equal(
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz[-1],
            newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz,
        )
        np.testing.assert_array_equal(
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[-1],
            newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz,
        )
        # The recorded entries are fresh arrays, not aliases of the caller's arrays.
        self.assertIsNot(
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz[-1],
            newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz,
        )
        self.assertIsNot(
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[-1],
            newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz,
        )

    def test_record_structural_state_holds_previous_state_during_discard_window(self):
        """Test that _record_structural_state ignores the new state during the discard
        window, appending a copy of the previous entry to both time series."""
        problem = self.problem
        wing_idx = 0
        num_wing_cross_sections = problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[
            wing_idx
        ][-1].shape[0]
        newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz = np.full(
            num_wing_cross_sections, 0.25, dtype=float
        )
        newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = np.full(
            num_wing_cross_sections, -0.5, dtype=float
        )

        problem._record_structural_state(
            step=problem.step_discards,
            newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz=(
                newDeformationAnglesYRad_Wcsp_to_Wcs_ixyz
            ),
            newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz=(
                newDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz
            ),
            wing_idx=wing_idx,
        )

        deformationAnglesYRad_Wcsp_to_Wcs_ixyz = (
            problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = (
            problem._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        self.assertEqual(len(deformationAnglesYRad_Wcsp_to_Wcs_ixyz), 2)
        self.assertEqual(len(deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz), 2)
        # Both series hold the seed (zero) state, not the passed values.
        np.testing.assert_array_equal(
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz[-1],
            np.zeros(num_wing_cross_sections, dtype=float),
        )
        np.testing.assert_array_equal(
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[-1],
            np.zeros(num_wing_cross_sections, dtype=float),
        )
        # The held entries are copies, not aliases of the previous entries.
        self.assertIsNot(
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz[-1],
            deformationAnglesYRad_Wcsp_to_Wcs_ixyz[-2],
        )
        self.assertIsNot(
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[-1],
            deformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[-2],
        )

    def test_spring_numerical_ode_zero_initial_zero_forces_theta_stays_zero(self):
        """Test that with zero initial conditions and zero external forces, theta
        remains zero."""
        t = np.array([0.0, 0.05])
        zero_moment_func = lambda time: 0.0
        theta_rad, theta_derivative_rad = self.problem.spring_numerical_ode(
            t,
            spring_constant_rad=10.0,
            damping_constant_rad=0.5,
            I=1.0,
            theta0_rad=0.0,
            theta_derivative0_rad=0.0,
            aero_moment=0.0,
            inertial_moment_func=zero_moment_func,
        )
        self.assertAlmostEqual(theta_rad, 0.0, places=6)
        self.assertAlmostEqual(theta_derivative_rad, 0.0, places=6)

    def test_spring_numerical_ode_returns_floats(self):
        """Test that spring_numerical_ode returns Python floats."""
        t = np.array([0.0, 0.1])
        zero_moment_func = lambda time: 0.0
        theta_rad, theta_derivative_rad = self.problem.spring_numerical_ode(
            t,
            spring_constant_rad=10.0,
            damping_constant_rad=0.5,
            I=1.0,
            theta0_rad=0.0,
            theta_derivative0_rad=0.0,
            aero_moment=0.0,
            inertial_moment_func=zero_moment_func,
        )
        self.assertIsInstance(theta_rad, float)
        self.assertIsInstance(theta_derivative_rad, float)

    def test_spring_numerical_ode_spring_restores_toward_zero(self):
        """Test that a positive initial displacement decays toward zero with
        positive spring constant and damping."""
        t = np.array([0.0, 0.5, 1.0])
        zero_moment_func = lambda time: 0.0
        theta0_rad = 1.0
        theta_rad, theta_derivative_rad = self.problem.spring_numerical_ode(
            t,
            spring_constant_rad=10.0,
            damping_constant_rad=2.0,
            I=1.0,
            theta0_rad=theta0_rad,
            theta_derivative0_rad=0.0,
            aero_moment=0.0,
            inertial_moment_func=zero_moment_func,
        )
        # With damping and spring, displacement should decrease in magnitude.
        self.assertLess(abs(theta_rad), abs(theta0_rad))

    def test_generate_inertial_moment_function_returns_callable(self):
        """Test that generate_inertial_moment_function returns a callable."""
        moment_func = self.problem.generate_inertial_moment_function(span_I=1.0)
        self.assertTrue(callable(moment_func))

    def test_generate_inertial_moment_function_returns_float_at_zero(self):
        """Test that the moment function returned by generate_inertial_moment_function
        returns a numeric value when evaluated at time zero."""
        moment_func = self.problem.generate_inertial_moment_function(span_I=1.0)
        result = moment_func(0.0)
        self.assertIsInstance(result, (float, np.floating))

    def test_generate_inertial_moment_function_scales_with_span_i(self):
        """Test that the moment function scales linearly with span_I."""
        moment_func_1 = self.problem.generate_inertial_moment_function(span_I=1.0)
        moment_func_2 = self.problem.generate_inertial_moment_function(span_I=2.0)
        t_eval = 0.25
        self.assertAlmostEqual(
            moment_func_2(t_eval), 2.0 * moment_func_1(t_eval), places=10
        )

    def test_generate_inertial_moment_function_sine_matches_analytic_value(self):
        """Test that the sinusoidal-spacing moment equals -I * b^2 * sin(b * t + h) * A
        with the amplitude converted from degrees to radians so the moment is in SI
        N*m."""
        wing_movement = self.problem.wing_movement
        amp_deg = wing_movement.ampAngles_Gs_to_Wn_ixyz[0]
        period = wing_movement.periodAngles_Gs_to_Wn_ixyz[0]
        phase_deg = wing_movement.phaseAngles_Gs_to_Wn_ixyz[0]

        span_I = 3.0
        t_eval = 0.3
        moment_func = self.problem.generate_inertial_moment_function(span_I=span_I)

        b_rad = 2.0 * np.pi / period
        h_rad = np.deg2rad(phase_deg)
        expected = (
            -1.0
            * b_rad**2
            * np.sin(b_rad * t_eval + h_rad)
            * np.deg2rad(amp_deg)
            * span_I
        )
        self.assertAlmostEqual(moment_func(t_eval), expected, places=10)

    def test_generate_inertial_moment_function_static_motion_returns_zero(self):
        """Test that a motion-off wing movement (zero flapping amplitude and period)
        produces an identically zero, finite inertial moment at every time.

        A wing with no prescribed flapping applies no inertial moment, so the
        returned function must evaluate to exactly 0.0 rather than to NaN from the
        2 * pi / period frequency computation with a zero period.
        """
        static_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_static_aeroelastic_wing_movement_fixture()
        )
        moment_func = self.problem.generate_inertial_moment_function(
            span_I=1.0, wing_movement=static_wing_movement
        )
        for t_eval in (0.0, 0.1, 0.5):
            result = moment_func(t_eval)
            self.assertTrue(
                np.isfinite(result),
                msg=f"The inertial moment at time {t_eval} is {result}, which is "
                f"not finite.",
            )
            self.assertAlmostEqual(float(result), 0.0, places=12)

    def test_generate_inertial_moment_function_uniform_spacing_raises(self):
        """Test that generate_inertial_moment_function raises ValueError when the
        wing motion spacing is "uniform" (sawtooth), which is not differentiable."""
        wing_movement = self.problem.wing_movement
        with patch.object(
            type(wing_movement),
            "spacingAngles_Gs_to_Wn_ixyz",
            new_callable=PropertyMock,
            return_value=("uniform", "sine", "sine"),
        ):
            with self.assertRaises(ValueError):
                self.problem.generate_inertial_moment_function(span_I=1.0)

    def test_generate_inertial_moment_function_callable_spacing_with_derivative(self):
        """Test that generate_inertial_moment_function uses the AeroelasticWingMovement's
        second derivative when the spacing is a custom callable."""
        wing_movement = self.problem.wing_movement
        with (
            patch.object(
                type(wing_movement),
                "spacingAngles_Gs_to_Wn_ixyz",
                new_callable=PropertyMock,
                return_value=(lambda t: np.sin(t), "sine", "sine"),
            ),
            patch.object(
                type(wing_movement),
                "spacingAnglesSecondDerivative_Gs_to_Wn_ixyz",
                new_callable=PropertyMock,
                return_value=(lambda t: -np.sin(t), None, None),
            ),
        ):
            moment_func = self.problem.generate_inertial_moment_function(span_I=2.0)
            self.assertTrue(callable(moment_func))
            result = moment_func(0.5)
            # The amplitude is stored in degrees but the moment (N*m) must be in SI
            # units, so the second derivative is scaled by np.deg2rad(amp).
            amp = wing_movement.ampAngles_Gs_to_Wn_ixyz[0]
            expected = np.deg2rad(amp) * -np.sin(0.5) * 2.0
            self.assertAlmostEqual(result, expected, places=10)


class TestRecordNullStepForWing(unittest.TestCase):
    """This class contains unit tests for the _record_null_step_for_wing method and
    the standard WingMovement code path in calculate_wing_deformation."""

    def setUp(self):
        """Set up a fresh AeroelasticUnsteadyProblem for each test."""
        self.problem_aero = (
            problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
        )
        self.problem_std = (
            problem_fixtures.make_aeroelastic_unsteady_problem_with_standard_wing_fixture()
        )

    def test_record_null_step_for_wing_appends_one_entry_to_each_time_series(self):
        """Test that _record_null_step_for_wing appends exactly one entry to each
        per-wing time series list."""
        wing_idx = 0
        problem = self.problem_aero
        num_previous_states = len(
            problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx]
        )
        # The two time series are seeded together and recorded together, so their
        # lengths are always equal.
        self.assertEqual(
            len(
                problem._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx]
            ),
            num_previous_states,
        )

        problem._record_null_step_for_wing(wing_idx)

        self.assertEqual(
            len(problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx]),
            num_previous_states + 1,
        )
        self.assertEqual(
            len(
                problem._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx]
            ),
            num_previous_states + 1,
        )

    def test_record_null_step_for_wing_entry_shape_and_values(self):
        """Test that _record_null_step_for_wing appends zero-valued entries, with one
        element per WingCrossSection, to listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz and
        _listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz."""
        wing_idx = 0
        num_wing_cross_sections = (
            self.problem_aero.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx][
                -1
            ].shape[0]
        )

        self.problem_aero._record_null_step_for_wing(wing_idx)

        lastDeformationAnglesYRad_Wcsp_to_Wcs_ixyz = (
            self.problem_aero.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx][-1]
        )
        lastDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz = (
            self.problem_aero._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[
                wing_idx
            ][-1]
        )
        expected_shape = (num_wing_cross_sections,)

        self.assertEqual(
            lastDeformationAnglesYRad_Wcsp_to_Wcs_ixyz.shape, expected_shape
        )
        self.assertEqual(
            lastDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz.shape, expected_shape
        )
        np.testing.assert_array_equal(
            lastDeformationAnglesYRad_Wcsp_to_Wcs_ixyz, np.zeros(expected_shape)
        )
        np.testing.assert_array_equal(
            lastDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz,
            np.zeros(expected_shape),
        )

    def test_calculate_wing_deformation_returns_none_for_standard_wing_movement(self):
        """Test that calculate_wing_deformation returns None for a Wing backed by a
        standard WingMovement (the else branch)."""
        mock_solver = MagicMock()

        results = self.problem_std.calculate_wing_deformation(
            solver=mock_solver, step=0
        )

        self.assertEqual(len(results), 1)
        self.assertIsNone(results[0])

    def test_calculate_wing_deformation_appends_entries_for_standard_wing_movement(
        self,
    ):
        """Test that calculate_wing_deformation appends to the time series lists when a
        Wing is backed by a standard WingMovement."""
        mock_solver = MagicMock()
        wing_idx = 0

        self.problem_std.calculate_wing_deformation(solver=mock_solver, step=0)

        problem = self.problem_std
        # The two structural state time series are seeded with an initial-state entry
        # at construction, so one recorded step brings them to two entries.
        self.assertEqual(
            len(problem.listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz[wing_idx]), 2
        )
        self.assertEqual(
            len(
                problem._listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz[wing_idx]
            ),
            2,
        )
