"""This module contains classes to test the AeroelasticUnsteadyProblem class."""

import unittest

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
        self.problem.positions.append(dummy_pos)
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
