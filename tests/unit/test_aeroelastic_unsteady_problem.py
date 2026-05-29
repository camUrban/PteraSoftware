"""This module contains classes to test the AeroelasticUnsteadyProblem class."""

import unittest
from unittest.mock import MagicMock, PropertyMock, patch

import numpy as np

import pterasoftware as ps
from tests.unit.fixtures import movement_fixtures, problem_fixtures


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
                spring_constant=10.0,
                damping_constant=0.5,
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

    def test_spring_constant_stored(self):
        """Test that spring_constant is stored correctly."""
        self.assertAlmostEqual(self.problem.spring_constant, 10.0)

    def test_damping_constant_stored(self):
        """Test that damping_constant is stored correctly."""
        self.assertAlmostEqual(self.problem.damping_constant, 0.5)

    def test_calculate_mass_matrix_shape(self):
        """Test that calculate_mass_matrix returns an array with shape
        (num_chordwise, num_spanwise, 3)."""
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

    def test_calculate_wing_panel_accelerations_no_positions_returns_scalar_zero(self):
        """Test that calculate_wing_panel_accelerations returns zeros(1) when no
        positions are stored yet."""
        accel = self.problem.calculate_wing_panel_accelerations()
        np.testing.assert_array_equal(accel, np.zeros(1))

    def test_calculate_wing_panel_accelerations_one_position_returns_zeros_like(self):
        """Test that calculate_wing_panel_accelerations returns zeros_like the
        position when only one position is stored."""
        dummy_pos = np.ones((2, 1, 3))
        self.problem.positions_per_wing[0].append(dummy_pos)
        accel = self.problem.calculate_wing_panel_accelerations()
        np.testing.assert_array_equal(accel, np.zeros_like(dummy_pos))

    def test_spring_numerical_ode_zero_initial_zero_forces_theta_stays_zero(self):
        """Test that with zero initial conditions and zero external forces, theta
        remains zero."""
        t = np.array([0.0, 0.05])
        zero_torque_func = lambda time: 0.0
        theta, omega = self.problem.spring_numerical_ode(
            t,
            k=10.0,
            c=0.5,
            I=1.0,
            theta0=0.0,
            omega0=0.0,
            aero_torque=0.0,
            inertial_torque_func=zero_torque_func,
        )
        self.assertAlmostEqual(theta, 0.0, places=6)
        self.assertAlmostEqual(omega, 0.0, places=6)

    def test_spring_numerical_ode_returns_floats(self):
        """Test that spring_numerical_ode returns Python floats."""
        t = np.array([0.0, 0.1])
        zero_torque_func = lambda time: 0.0
        theta, omega = self.problem.spring_numerical_ode(
            t,
            k=10.0,
            c=0.5,
            I=1.0,
            theta0=0.0,
            omega0=0.0,
            aero_torque=0.0,
            inertial_torque_func=zero_torque_func,
        )
        self.assertIsInstance(theta, float)
        self.assertIsInstance(omega, float)

    def test_spring_numerical_ode_spring_restores_toward_zero(self):
        """Test that a positive initial displacement decays toward zero with
        positive spring constant and damping."""
        t = np.array([0.0, 0.5, 1.0])
        zero_torque_func = lambda time: 0.0
        theta_init = 1.0
        theta, omega = self.problem.spring_numerical_ode(
            t,
            k=10.0,
            c=2.0,
            I=1.0,
            theta0=theta_init,
            omega0=0.0,
            aero_torque=0.0,
            inertial_torque_func=zero_torque_func,
        )
        # With damping and spring, displacement should decrease in magnitude.
        self.assertLess(abs(theta), abs(theta_init))

    def test_generate_inertial_torque_function_returns_callable(self):
        """Test that generate_inertial_torque_function returns a callable."""
        torque_func = self.problem.generate_inertial_torque_function(span_I=1.0)
        self.assertTrue(callable(torque_func))

    def test_generate_inertial_torque_function_returns_float_at_zero(self):
        """Test that the torque function returned by generate_inertial_torque_function
        returns a numeric value when evaluated at time zero."""
        torque_func = self.problem.generate_inertial_torque_function(span_I=1.0)
        result = torque_func(0.0)
        self.assertIsInstance(result, (float, np.floating))

    def test_generate_inertial_torque_function_scales_with_span_i(self):
        """Test that the torque function scales linearly with span_I."""
        torque_func_1 = self.problem.generate_inertial_torque_function(span_I=1.0)
        torque_func_2 = self.problem.generate_inertial_torque_function(span_I=2.0)
        t_eval = 0.25
        self.assertAlmostEqual(
            torque_func_2(t_eval), 2.0 * torque_func_1(t_eval), places=10
        )

    def test_generate_inertial_torque_function_uniform_spacing_raises(self):
        """Test that generate_inertial_torque_function raises ValueError when the
        wing motion spacing is "uniform" (sawtooth), which is not differentiable."""
        wing_movement = self.problem.wing_movement
        with patch.object(
            type(wing_movement),
            "spacingAngles_Gs_to_Wn_ixyz",
            new_callable=PropertyMock,
            return_value=("uniform", "sine", "sine"),
        ):
            with self.assertRaises(ValueError):
                self.problem.generate_inertial_torque_function(span_I=1.0)

    def test_generate_inertial_torque_function_callable_spacing_with_derivative(self):
        """Test that generate_inertial_torque_function uses the wing movement's second
        derivative when the spacing is a custom callable."""
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
            torque_func = self.problem.generate_inertial_torque_function(span_I=2.0)
            self.assertTrue(callable(torque_func))
            result = torque_func(0.5)
            self.assertAlmostEqual(result, -np.sin(0.5) * 2.0, places=10)

    def test_plot_aeroelastic_results_calls_plot_flap_cycle_curves_four_times(self):
        """Test that _plot_aeroelastic_results calls plot_flap_cycle_curves exactly
        four times with the correct titles."""
        self.problem.per_step_inertial_per_wing[0] = [np.zeros((1, 1, 3))]
        self.problem.per_step_aero_per_wing[0] = [np.zeros((1, 1, 3))]
        self.problem.net_data_per_wing[0] = [np.zeros((2, 3))]
        self.problem.flap_points_per_wing[0] = [np.zeros((1, 1, 3))]

        with patch.object(type(self.problem), "plot_flap_cycle_curves") as mock_plot:
            self.problem._plot_aeroelastic_results()

        self.assertEqual(mock_plot.call_count, 4)
        titles = [call.args[1] for call in mock_plot.call_args_list]
        self.assertIn("Net Deformation", titles)
        self.assertIn("Per Step Inertial Moments", titles)
        self.assertIn("Per Step Aero Moments", titles)
        self.assertIn("Flap Points Z", titles)

    def test_plot_flap_cycle_curves_executes_without_error(self):
        """Test that plot_flap_cycle_curves runs without error when matplotlib
        calls are mocked."""
        with patch("pterasoftware.problems.plt") as mock_plt:
            self.problem.plot_flap_cycle_curves(
                data=[[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
                title="Test Plot",
            )

        mock_plt.figure.assert_called_once()
        mock_plt.show.assert_called_once()
        mock_plt.savefig.assert_called_once_with("Test_Plot.png")

    def test_plot_flap_cycle_curves_with_flap_cycle_overlay(self):
        """Test that plot_flap_cycle_curves plots the flap_cycle overlay when
        provided."""
        with patch("pterasoftware.problems.plt") as mock_plt:
            self.problem.plot_flap_cycle_curves(
                data=[[0.0, 1.0, 2.0]],
                title="Test Plot",
                flap_cycle=[0.5, 1.5, 2.5],
            )

        # One call for the data curve, one for the flap_cycle overlay.
        self.assertEqual(mock_plt.plot.call_count, 2)

    def test_calculate_wing_deformation_calls_plot_at_final_step(self):
        """Test that calculate_wing_deformation calls _plot_aeroelastic_results when
        plot_flap_cycle=True and step equals num_steps-1."""
        problem = problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture(
            plot_flap_cycle=True
        )

        mock_solver = MagicMock()
        dummy_moments = np.zeros((1, 1, 3))
        dummy_thetas = np.zeros(
            problem.steady_problems[0].airplanes[0].wings[0].num_spanwise_panels + 1
        )
        dummy_omegas = np.zeros_like(dummy_thetas)
        dummy_deformation = np.zeros(
            (
                problem.steady_problems[0].airplanes[0].wings[0].num_spanwise_panels
                + 1,
                3,
            )
        )

        problem_type = type(problem)
        with (
            patch.object(
                problem_type, "_extract_aero_moments", return_value=dummy_moments
            ),
            patch.object(
                problem_type, "_calculate_inertial_moments", return_value=dummy_moments
            ),
            patch.object(
                problem_type,
                "calculate_spring_moments",
                return_value=(dummy_thetas, dummy_omegas),
            ),
            patch.object(
                problem_type,
                "_build_deformation_vector",
                return_value=dummy_deformation,
            ),
            patch.object(problem_type, "_apply_moment_updates"),
            patch.object(problem_type, "_plot_aeroelastic_results") as mock_plot,
        ):
            final_step = problem.num_steps - 1
            problem.calculate_wing_deformation(mock_solver, step=final_step)

        mock_plot.assert_called_once()


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

    def test_record_null_step_for_wing_appends_one_entry_to_each_history_list(self):
        """Test that _record_null_step_for_wing appends exactly one entry to each
        per-wing history list."""
        wing_idx = 0
        wing = self.problem_aero.steady_problems[0].airplanes[0].wings[0]
        before = len(self.problem_aero.per_step_inertial_per_wing[wing_idx])

        self.problem_aero._record_null_step_for_wing(wing_idx, wing, step=0)

        self.assertEqual(
            len(self.problem_aero.per_step_inertial_per_wing[wing_idx]), before + 1
        )
        self.assertEqual(
            len(self.problem_aero.per_step_aero_per_wing[wing_idx]), before + 1
        )
        self.assertEqual(len(self.problem_aero.net_data_per_wing[wing_idx]), before + 1)
        self.assertEqual(
            len(self.problem_aero.angular_velocity_data_per_wing[wing_idx]), before + 1
        )
        self.assertEqual(
            len(self.problem_aero.flap_points_per_wing[wing_idx]), before + 1
        )

    def test_record_null_step_for_wing_inertial_and_aero_shape_and_values(self):
        """Test that _record_null_step_for_wing appends zero moment arrays with shape
        (num_chordwise_panels, num_spanwise_panels, 3)."""
        wing_idx = 0
        wing = self.problem_aero.steady_problems[0].airplanes[0].wings[0]

        self.problem_aero._record_null_step_for_wing(wing_idx, wing, step=0)

        expected_shape = (wing.num_chordwise_panels, wing.num_spanwise_panels, 3)
        appended_inertial = self.problem_aero.per_step_inertial_per_wing[wing_idx][-1]
        appended_aero = self.problem_aero.per_step_aero_per_wing[wing_idx][-1]

        self.assertEqual(appended_inertial.shape, expected_shape)
        self.assertEqual(appended_aero.shape, expected_shape)
        np.testing.assert_array_equal(appended_inertial, np.zeros(expected_shape))
        np.testing.assert_array_equal(appended_aero, np.zeros(expected_shape))

    def test_record_null_step_for_wing_net_data_shape_and_values(self):
        """Test that _record_null_step_for_wing appends a zero net_data array with
        shape (num_deformation_rows, 3)."""
        wing_idx = 0
        wing = self.problem_aero.steady_problems[0].airplanes[0].wings[0]
        num_deformation_rows = self.problem_aero.net_deformation_per_wing[
            wing_idx
        ].shape[0]

        self.problem_aero._record_null_step_for_wing(wing_idx, wing, step=0)

        appended_net = self.problem_aero.net_data_per_wing[wing_idx][-1]
        appended_ang = self.problem_aero.angular_velocity_data_per_wing[wing_idx][-1]
        expected_shape = (num_deformation_rows, 3)

        self.assertEqual(appended_net.shape, expected_shape)
        self.assertEqual(appended_ang.shape, expected_shape)
        np.testing.assert_array_equal(appended_net, np.zeros(expected_shape))
        np.testing.assert_array_equal(appended_ang, np.zeros(expected_shape))

    def test_record_null_step_for_wing_flap_points_shape_and_values(self):
        """Test that _record_null_step_for_wing appends a zero flap_points array with
        shape (num_chordwise_panels, num_spanwise_panels, 3)."""
        wing_idx = 0
        wing = self.problem_aero.steady_problems[0].airplanes[0].wings[0]

        self.problem_aero._record_null_step_for_wing(wing_idx, wing, step=0)

        appended_flap = self.problem_aero.flap_points_per_wing[wing_idx][-1]
        expected_shape = (wing.num_chordwise_panels, wing.num_spanwise_panels, 3)

        self.assertEqual(appended_flap.shape, expected_shape)
        np.testing.assert_array_equal(appended_flap, np.zeros(expected_shape))

    def test_calculate_wing_deformation_returns_none_for_standard_wing_movement(self):
        """Test that calculate_wing_deformation returns None for a wing backed by a
        standard WingMovement (the else branch)."""
        mock_solver = MagicMock()

        results = self.problem_std.calculate_wing_deformation(
            solver=mock_solver, step=0
        )

        self.assertEqual(len(results), 1)
        self.assertIsNone(results[0])

    def test_calculate_wing_deformation_appends_history_for_standard_wing_movement(
        self,
    ):
        """Test that calculate_wing_deformation populates history lists when a wing is
        backed by a standard WingMovement."""
        mock_solver = MagicMock()
        wing_idx = 0

        self.problem_std.calculate_wing_deformation(solver=mock_solver, step=0)

        self.assertEqual(len(self.problem_std.per_step_inertial_per_wing[wing_idx]), 1)
        self.assertEqual(len(self.problem_std.per_step_aero_per_wing[wing_idx]), 1)
        self.assertEqual(len(self.problem_std.net_data_per_wing[wing_idx]), 1)
        self.assertEqual(
            len(self.problem_std.angular_velocity_data_per_wing[wing_idx]), 1
        )
        self.assertEqual(len(self.problem_std.flap_points_per_wing[wing_idx]), 1)
