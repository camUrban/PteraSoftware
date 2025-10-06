"""This module contains a class to test OperatingPointMovements."""

import unittest
import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps

from tests.unit.fixtures import operating_point_movement_fixtures
from tests.unit.fixtures import operating_point_fixtures


class TestOperatingPointMovement(unittest.TestCase):
    """This is a class with functions to test OperatingPointMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all OperatingPointMovement tests."""
        cls.static_op_movement = (
            operating_point_movement_fixtures.make_static_operating_point_movement_fixture()
        )
        cls.sine_spacing_op_movement = (
            operating_point_movement_fixtures.make_sine_spacing_operating_point_movement_fixture()
        )
        cls.uniform_spacing_op_movement = (
            operating_point_movement_fixtures.make_uniform_spacing_operating_point_movement_fixture()
        )
        cls.phase_offset_op_movement = (
            operating_point_movement_fixtures.make_phase_offset_operating_point_movement_fixture()
        )
        cls.custom_spacing_op_movement = (
            operating_point_movement_fixtures.make_custom_spacing_operating_point_movement_fixture()
        )
        cls.basic_op_movement = (
            operating_point_movement_fixtures.make_basic_operating_point_movement_fixture()
        )
        cls.large_amplitude_op_movement = (
            operating_point_movement_fixtures.make_large_amplitude_operating_point_movement_fixture()
        )
        cls.long_period_op_movement = (
            operating_point_movement_fixtures.make_long_period_operating_point_movement_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test OperatingPointMovement initialization with valid parameters."""
        # Test that basic OperatingPointMovement initializes correctly.
        op_movement = self.basic_op_movement
        self.assertIsInstance(
            op_movement,
            ps.movements.operating_point_movement.OperatingPointMovement,
        )
        self.assertIsInstance(
            op_movement.base_operating_point,
            ps.operating_point.OperatingPoint,
        )
        self.assertEqual(op_movement.ampVCg__E, 5.0)
        self.assertEqual(op_movement.periodVCg__E, 2.0)
        self.assertEqual(op_movement.spacingVCg__E, "sine")
        self.assertEqual(op_movement.phaseVCg__E, 0.0)

    def test_initialization_default_parameters(self):
        """Test OperatingPointMovement initialization with default parameters."""
        base_operating_point = operating_point_fixtures.make_basic_operating_point_fixture()
        op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_operating_point
        )

        self.assertEqual(op_movement.ampVCg__E, 0.0)
        self.assertEqual(op_movement.periodVCg__E, 0.0)
        self.assertEqual(op_movement.spacingVCg__E, "sine")
        self.assertEqual(op_movement.phaseVCg__E, 0.0)

    def test_base_operating_point_validation(self):
        """Test that base_operating_point parameter is properly validated."""
        # Test with invalid type.
        with self.assertRaises(TypeError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point="not an operating point"
            )

        # Test with None.
        with self.assertRaises(TypeError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=None
            )

        # Test with valid OperatingPoint works.
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()
        op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op
        )
        self.assertEqual(op_movement.base_operating_point, base_op)

    def test_ampVCg__E_validation(self):
        """Test ampVCg__E parameter validation."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test valid non-negative values.
        valid_amps = [0.0, 1.0, 5.0, 100.0]
        for amp in valid_amps:
            with self.subTest(amp=amp):
                op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                    base_operating_point=base_op, ampVCg__E=amp
                )
                self.assertEqual(op_movement.ampVCg__E, amp)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op, ampVCg__E=-1.0
            )

        # Test invalid types raise error.
        with self.assertRaises((TypeError, ValueError)):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op, ampVCg__E="invalid"
            )

    def test_periodVCg__E_validation(self):
        """Test periodVCg__E parameter validation."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test valid non-negative values.
        valid_periods = [0.0, 1.0, 5.0, 100.0]
        for period in valid_periods:
            with self.subTest(period=period):
                amp = 1.0 if period > 0 else 0.0
                op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                    base_operating_point=base_op,
                    ampVCg__E=amp,
                    periodVCg__E=period,
                )
                self.assertEqual(op_movement.periodVCg__E, period)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=-1.0,
            )

    def test_spacingVCg__E_validation(self):
        """Test spacingVCg__E parameter validation."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test valid string values.
        valid_spacings = ["sine", "uniform"]
        for spacing in valid_spacings:
            with self.subTest(spacing=spacing):
                op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                    base_operating_point=base_op, spacingVCg__E=spacing
                )
                self.assertEqual(op_movement.spacingVCg__E, spacing)

        # Test invalid string raises error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op, spacingVCg__E="invalid"
            )

        # Test callable is accepted.
        def custom_func(x):
            return np.sin(x)

        op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op, spacingVCg__E=custom_func
        )
        self.assertTrue(callable(op_movement.spacingVCg__E))

        # Test non-callable, non-string raises error.
        with self.assertRaises(TypeError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op, spacingVCg__E=123
            )

    def test_phaseVCg__E_validation(self):
        """Test phaseVCg__E parameter validation."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test valid phase values within range (-180.0, 180.0].
        valid_phases = [0.0, 90.0, 180.0, -90.0, -179.9]
        for phase in valid_phases:
            with self.subTest(phase=phase):
                amp = 1.0 if phase != 0 else 0.0
                period = 1.0 if phase != 0 else 0.0
                op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                    base_operating_point=base_op,
                    ampVCg__E=amp,
                    periodVCg__E=period,
                    phaseVCg__E=phase,
                )
                self.assertEqual(op_movement.phaseVCg__E, phase)

        # Test phase > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                phaseVCg__E=180.1,
            )

        # Test phase <= -180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                phaseVCg__E=-180.0,
            )

    def test_amp_period_relationship(self):
        """Test that if ampVCg__E is 0, periodVCg__E must be 0."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test amp=0 with period=0 works.
        op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op,
            ampVCg__E=0.0,
            periodVCg__E=0.0,
        )
        self.assertIsNotNone(op_movement)

        # Test amp=0 with period!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=0.0,
                periodVCg__E=1.0,
            )

    def test_amp_phase_relationship(self):
        """Test that if ampVCg__E is 0, phaseVCg__E must be 0."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test amp=0 with phase=0 works.
        op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op,
            ampVCg__E=0.0,
            periodVCg__E=0.0,
            phaseVCg__E=0.0,
        )
        self.assertIsNotNone(op_movement)

        # Test amp=0 with phase!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=0.0,
                periodVCg__E=0.0,
                phaseVCg__E=45.0,
            )

    def test_spacing_sine_produces_sinusoidal_motion(self):
        """Test that sine spacing actually produces sinusoidal motion for vCg__E."""
        num_steps = 100
        delta_time = 0.01
        operating_points = self.sine_spacing_op_movement.generate_operating_points(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract vCg__E values from generated OperatingPoints.
        vCg_values = np.array([op.vCg__E for op in operating_points], dtype=float)

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_vCg = 100.0 + 10.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated values match the expected sine wave.
        npt.assert_allclose(vCg_values, expected_vCg, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_produces_triangular_wave(self):
        """Test that uniform spacing actually produces triangular wave motion for vCg__E."""
        num_steps = 100
        delta_time = 0.01
        operating_points = self.uniform_spacing_op_movement.generate_operating_points(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract vCg__E values from generated OperatingPoints.
        vCg_values = np.array([op.vCg__E for op in operating_points], dtype=float)

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_vCg = 100.0 + 10.0 * signal.sawtooth(
            2 * np.pi * times / 1.0 + np.pi / 2, 0.5
        )

        # Assert that the generated values match the expected triangular wave.
        npt.assert_allclose(vCg_values, expected_vCg, rtol=1e-10, atol=1e-14)

    def test_max_period_static_movement(self):
        """Test that max_period returns 0.0 for static movement."""
        op_movement = self.static_op_movement
        self.assertEqual(op_movement.max_period, 0.0)

    def test_max_period_with_movement(self):
        """Test that max_period returns correct period for movement."""
        op_movement = self.basic_op_movement
        # periodVCg__E is 2.0, so max should be 2.0.
        self.assertEqual(op_movement.max_period, 2.0)

    def test_max_period_long_period(self):
        """Test that max_period returns correct value for long period movement."""
        op_movement = self.long_period_op_movement
        # periodVCg__E is 10.0, so max should be 10.0.
        self.assertEqual(op_movement.max_period, 10.0)

    def test_generate_operating_points_parameter_validation(self):
        """Test that generate_operating_points validates num_steps and delta_time."""
        op_movement = self.basic_op_movement

        # Test invalid num_steps.
        with self.assertRaises((ValueError, TypeError)):
            op_movement.generate_operating_points(num_steps=0, delta_time=0.01)

        with self.assertRaises((ValueError, TypeError)):
            op_movement.generate_operating_points(num_steps=-1, delta_time=0.01)

        with self.assertRaises(TypeError):
            op_movement.generate_operating_points(
                num_steps="invalid", delta_time=0.01
            )

        # Test invalid delta_time.
        with self.assertRaises((ValueError, TypeError)):
            op_movement.generate_operating_points(num_steps=10, delta_time=0.0)

        with self.assertRaises((ValueError, TypeError)):
            op_movement.generate_operating_points(num_steps=10, delta_time=-0.01)

        with self.assertRaises(TypeError):
            op_movement.generate_operating_points(num_steps=10, delta_time="invalid")

    def test_generate_operating_points_returns_correct_length(self):
        """Test that generate_operating_points returns list of correct length."""
        op_movement = self.basic_op_movement

        test_num_steps = [1, 5, 10, 50, 100]
        for num_steps in test_num_steps:
            with self.subTest(num_steps=num_steps):
                operating_points = op_movement.generate_operating_points(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(operating_points), num_steps)

    def test_generate_operating_points_returns_correct_types(self):
        """Test that generate_operating_points returns OperatingPoints."""
        op_movement = self.basic_op_movement
        operating_points = op_movement.generate_operating_points(
            num_steps=10, delta_time=0.01
        )

        # Verify all elements are OperatingPoints.
        for op in operating_points:
            self.assertIsInstance(op, ps.operating_point.OperatingPoint)

    def test_generate_operating_points_preserves_non_changing_attributes(self):
        """Test that generate_operating_points preserves non-changing attributes."""
        op_movement = self.basic_op_movement
        base_op = op_movement.base_operating_point

        operating_points = op_movement.generate_operating_points(
            num_steps=10, delta_time=0.01
        )

        # Check that non-changing attributes are preserved.
        for op in operating_points:
            self.assertEqual(op.rho, base_op.rho)
            self.assertEqual(op.alpha, base_op.alpha)
            self.assertEqual(op.beta, base_op.beta)
            self.assertEqual(op.externalFX_W, base_op.externalFX_W)
            self.assertEqual(op.nu, base_op.nu)

    def test_generate_operating_points_static_movement(self):
        """Test that static movement produces constant vCg__E."""
        op_movement = self.static_op_movement
        base_op = op_movement.base_operating_point

        operating_points = op_movement.generate_operating_points(
            num_steps=50, delta_time=0.01
        )

        # All OperatingPoints should have same vCg__E.
        for op in operating_points:
            self.assertEqual(op.vCg__E, base_op.vCg__E)

    def test_generate_operating_points_different_num_steps(self):
        """Test generate_operating_points with various num_steps values."""
        op_movement = self.basic_op_movement

        num_steps_list = [1, 10, 25, 100, 200]
        for num_steps in num_steps_list:
            with self.subTest(num_steps=num_steps):
                operating_points = op_movement.generate_operating_points(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(operating_points), num_steps)

    def test_generate_operating_points_different_delta_time(self):
        """Test generate_operating_points with various delta_time values."""
        op_movement = self.basic_op_movement

        delta_time_list = [0.001, 0.01, 0.1, 1.0]
        num_steps = 50
        for delta_time in delta_time_list:
            with self.subTest(delta_time=delta_time):
                operating_points = op_movement.generate_operating_points(
                    num_steps=num_steps, delta_time=delta_time
                )
                self.assertEqual(len(operating_points), num_steps)

    def test_phase_offset_shifts_initial_value(self):
        """Test that phase shifts initial vCg__E correctly."""
        op_movement = self.phase_offset_op_movement
        operating_points = op_movement.generate_operating_points(
            num_steps=100, delta_time=0.01
        )

        # Extract vCg__E values.
        vCg_values = np.array([op.vCg__E for op in operating_points], dtype=float)

        # Verify that phase offset causes non-zero initial value different from base.
        # With 90 degree phase offset, the first value should not be at the base.
        self.assertFalse(
            np.isclose(vCg_values[0], op_movement.base_operating_point.vCg__E, atol=1e-10)
        )

    def test_custom_spacing_function_works(self):
        """Test that custom spacing function works for vCg__E."""
        op_movement = self.custom_spacing_op_movement
        operating_points = op_movement.generate_operating_points(
            num_steps=100, delta_time=0.01
        )

        # Extract vCg__E values.
        vCg_values = np.array([op.vCg__E for op in operating_points], dtype=float)

        # Verify that values vary (not constant).
        self.assertFalse(np.allclose(vCg_values, vCg_values[0]))

        # Verify that values are within expected range.
        # For custom_harmonic with amp=15.0 and base=100.0, values should be roughly
        # in [85.0, 115.0].
        self.assertTrue(np.all(vCg_values >= 80.0))
        self.assertTrue(np.all(vCg_values <= 120.0))

    def test_custom_function_validation_invalid_start_value(self):
        """Test that custom function with invalid start value raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function that doesn't start at 0.
        def invalid_nonzero_start(x):
            return np.sin(x) + 1.0

        # Should raise error during generation.
        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_nonzero_start,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_invalid_end_value(self):
        """Test that custom function with invalid end value raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function that doesn't return to 0 at 2*pi.
        def invalid_nonzero_end(x):
            return np.sin(x) + 0.1

        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_nonzero_end,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_invalid_mean(self):
        """Test that custom function with invalid mean raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function with non-zero mean.
        def invalid_nonzero_mean(x):
            return np.sin(x) + 0.5

        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_nonzero_mean,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_invalid_amplitude(self):
        """Test that custom function with invalid amplitude raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function with wrong amplitude.
        def invalid_wrong_amplitude(x):
            return 2.0 * np.sin(x)

        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_wrong_amplitude,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_not_periodic(self):
        """Test that custom function that is not periodic raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function that is not periodic.
        def invalid_not_periodic(x):
            return np.tanh(x)

        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_not_periodic,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_returns_non_finite(self):
        """Test that custom function returning NaN or Inf raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function that returns NaN.
        def invalid_non_finite(x):
            return np.where(x < np.pi, np.sin(x), np.nan)

        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_non_finite,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_wrong_shape(self):
        """Test that custom function returning wrong shape raises error."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Define invalid custom function that returns wrong shape.
        def invalid_wrong_shape(x):
            return np.sin(x)[: len(x) // 2]

        with self.assertRaises(ValueError):
            op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_op,
                ampVCg__E=1.0,
                periodVCg__E=1.0,
                spacingVCg__E=invalid_wrong_shape,
            )
            op_movement.generate_operating_points(num_steps=10, delta_time=0.01)

    def test_unsafe_amplitude_causes_error(self):
        """Test that amplitude too high for base vCg__E causes error during generation."""
        # Use low-speed operating point with vCg__E = 10.0.
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Create OperatingPointMovement with amplitude that will drive vCg__E negative.
        op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op,
            ampVCg__E=15.0,
            periodVCg__E=1.0,
            spacingVCg__E="sine",
            phaseVCg__E=-90.0,
        )

        # Generating OperatingPoints should raise ValueError when vCg__E goes negative.
        with self.assertRaises(ValueError) as context:
            op_movement.generate_operating_points(num_steps=100, delta_time=0.01)

        # Verify the error message is about vCg__E validation.
        self.assertIn("vCg__E", str(context.exception))

    def test_boundary_phase_values(self):
        """Test phase at boundary values (-179.9, 0.0, and 180.0)."""
        base_op = operating_point_fixtures.make_basic_operating_point_fixture()

        # Test phase = 0.0 works.
        op_movement1 = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op,
            ampVCg__E=0.0,
            periodVCg__E=0.0,
            phaseVCg__E=0.0,
        )
        self.assertEqual(op_movement1.phaseVCg__E, 0.0)

        # Test phase = 180.0 works (upper boundary, inclusive).
        op_movement2 = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op,
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            phaseVCg__E=180.0,
        )
        self.assertEqual(op_movement2.phaseVCg__E, 180.0)

        # Test phase = -179.9 works (near lower boundary).
        op_movement3 = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=base_op,
            ampVCg__E=1.0,
            periodVCg__E=1.0,
            phaseVCg__E=-179.9,
        )
        self.assertEqual(op_movement3.phaseVCg__E, -179.9)

    def test_large_amplitude_movement(self):
        """Test OperatingPointMovement with large amplitude."""
        op_movement = self.large_amplitude_op_movement
        operating_points = op_movement.generate_operating_points(
            num_steps=100, delta_time=0.01
        )

        # Verify that all OperatingPoints have valid vCg__E values.
        for op in operating_points:
            self.assertGreater(op.vCg__E, 0.0)

        # Extract vCg__E values.
        vCg_values = np.array([op.vCg__E for op in operating_points], dtype=float)

        # Verify that values vary significantly.
        self.assertGreater(np.max(vCg_values) - np.min(vCg_values), 50.0)


if __name__ == "__main__":
    unittest.main()