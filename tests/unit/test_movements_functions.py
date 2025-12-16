"""This module contains a class to test movements functions."""

import unittest

import numpy as np
import numpy.testing as npt
from scipy import signal

# noinspection PyProtectedMember
from pterasoftware.movements import _functions
from tests.unit.fixtures import movements_functions_fixtures


class TestMovementsFunctions(unittest.TestCase):
    """This is a class with functions to test movements functions."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all movements function tests."""
        # Parameter fixtures.
        (
            cls.scalar_amps,
            cls.scalar_periods,
            cls.scalar_phases,
            cls.scalar_bases,
        ) = movements_functions_fixtures.make_scalar_parameters_fixture()

        (
            cls.array_amps,
            cls.array_periods,
            cls.array_phases,
            cls.array_bases,
        ) = movements_functions_fixtures.make_array_parameters_fixture()

        (
            cls.static_amps,
            cls.static_periods,
            cls.static_phases,
            cls.static_bases,
        ) = movements_functions_fixtures.make_static_parameters_fixture()

        (
            cls.mixed_static_amps,
            cls.mixed_static_periods,
            cls.mixed_static_phases,
            cls.mixed_static_bases,
        ) = movements_functions_fixtures.make_mixed_static_parameters_fixture()

        (
            cls.phase_offset_amps,
            cls.phase_offset_periods,
            cls.phase_offset_phases,
            cls.phase_offset_bases,
        ) = movements_functions_fixtures.make_phase_offset_parameters_fixture()

        (
            cls.large_amplitude_amps,
            cls.large_amplitude_periods,
            cls.large_amplitude_phases,
            cls.large_amplitude_bases,
        ) = movements_functions_fixtures.make_large_amplitude_parameters_fixture()

        (
            cls.small_period_amps,
            cls.small_period_periods,
            cls.small_period_phases,
            cls.small_period_bases,
        ) = movements_functions_fixtures.make_small_period_parameters_fixture()

        (
            cls.negative_phase_amps,
            cls.negative_phase_periods,
            cls.negative_phase_phases,
            cls.negative_phase_bases,
        ) = movements_functions_fixtures.make_negative_phase_parameters_fixture()

        (
            cls.max_phase_amps,
            cls.max_phase_periods,
            cls.max_phase_phases,
            cls.max_phase_bases,
        ) = movements_functions_fixtures.make_max_phase_parameters_fixture()

        (
            cls.min_phase_amps,
            cls.min_phase_periods,
            cls.min_phase_phases,
            cls.min_phase_bases,
        ) = movements_functions_fixtures.make_min_phase_parameters_fixture()

        # Time step fixtures.
        cls.num_steps, cls.delta_time = (
            movements_functions_fixtures.make_num_steps_and_delta_time_fixture()
        )
        cls.small_num_steps = (
            movements_functions_fixtures.make_small_num_steps_fixture()
        )
        cls.large_num_steps = (
            movements_functions_fixtures.make_large_num_steps_fixture()
        )

        # Custom function fixtures.
        cls.valid_custom_sine = staticmethod(
            movements_functions_fixtures.make_valid_custom_sine_function_fixture()
        )
        cls.valid_custom_triangle = staticmethod(
            movements_functions_fixtures.make_valid_custom_triangle_function_fixture()
        )
        cls.valid_custom_harmonic = staticmethod(
            movements_functions_fixtures.make_valid_custom_harmonic_function_fixture()
        )
        cls.invalid_custom_wrong_start = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_wrong_start_fixture()
        )
        cls.invalid_custom_wrong_end = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_wrong_end_fixture()
        )
        cls.invalid_custom_wrong_amplitude = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_wrong_amplitude_fixture()
        )
        cls.invalid_custom_not_periodic = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_not_periodic_fixture()
        )
        cls.invalid_custom_returns_nan = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_returns_nan_fixture()
        )
        cls.invalid_custom_returns_inf = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_returns_inf_fixture()
        )
        cls.invalid_custom_wrong_shape = staticmethod(
            movements_functions_fixtures.make_invalid_custom_function_wrong_shape_fixture()
        )

    def test_oscillating_sinspaces_scalar_parameters(self):
        """Test oscillating_sinspaces with scalar parameters."""
        result = _functions.oscillating_sinspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        self.assertEqual(result.shape, (self.num_steps,))

        # Verify output values match expected sine wave.
        times = np.linspace(
            0, self.num_steps * self.delta_time, self.num_steps, endpoint=False
        )
        expected = (
            self.scalar_amps
            * np.sin(
                2 * np.pi * times / self.scalar_periods + np.deg2rad(self.scalar_phases)
            )
            + self.scalar_bases
        )

        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_sinspaces_array_parameters(self):
        """Test oscillating_sinspaces with array parameters."""
        result = _functions.oscillating_sinspaces(
            amps=self.array_amps,
            periods=self.array_periods,
            phases=self.array_phases,
            bases=self.array_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        expected_shape = (3, self.num_steps)
        self.assertEqual(result.shape, expected_shape)

        # Verify output values for each element.
        times = np.linspace(
            0, self.num_steps * self.delta_time, self.num_steps, endpoint=False
        )
        for i in range(3):
            expected = (
                self.array_amps[i]
                * np.sin(
                    2 * np.pi * times / self.array_periods[i]
                    + np.deg2rad(self.array_phases[i])
                )
                + self.array_bases[i]
            )
            npt.assert_allclose(result[i, :], expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_sinspaces_static_parameters(self):
        """Test oscillating_sinspaces with static parameters."""
        result = _functions.oscillating_sinspaces(
            amps=self.static_amps,
            periods=self.static_periods,
            phases=self.static_phases,
            bases=self.static_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output is constant and equal to base.
        expected = np.full(self.num_steps, self.static_bases, dtype=float)
        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_sinspaces_mixed_static_parameters(self):
        """Test oscillating_sinspaces with mixed static and dynamic parameters."""
        result = _functions.oscillating_sinspaces(
            amps=self.mixed_static_amps,
            periods=self.mixed_static_periods,
            phases=self.mixed_static_phases,
            bases=self.mixed_static_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        expected_shape = (3, self.num_steps)
        self.assertEqual(result.shape, expected_shape)

        # Verify that static element (index 1) is constant.
        expected_static = np.full(
            self.num_steps, self.mixed_static_bases[1], dtype=float
        )
        npt.assert_allclose(result[1, :], expected_static, rtol=1e-10, atol=1e-14)

        # Verify that dynamic elements are oscillating.
        self.assertFalse(np.allclose(result[0, :], result[0, 0]))
        self.assertFalse(np.allclose(result[2, :], result[2, 0]))

    def test_oscillating_sinspaces_phase_offset(self):
        """Test oscillating_sinspaces with phase offset."""
        result = _functions.oscillating_sinspaces(
            amps=self.phase_offset_amps,
            periods=self.phase_offset_periods,
            phases=self.phase_offset_phases,
            bases=self.phase_offset_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify that phase offset shifts the waveform.
        # At t=0, sine with 90 degree phase should equal the amplitude.
        npt.assert_allclose(result[0], self.phase_offset_amps, rtol=1e-10, atol=1e-14)

    def test_oscillating_sinspaces_large_amplitude(self):
        """Test oscillating_sinspaces with large amplitude."""
        result = _functions.oscillating_sinspaces(
            amps=self.large_amplitude_amps,
            periods=self.large_amplitude_periods,
            phases=self.large_amplitude_phases,
            bases=self.large_amplitude_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output range.
        expected_min = self.large_amplitude_bases - self.large_amplitude_amps
        expected_max = self.large_amplitude_bases + self.large_amplitude_amps
        self.assertGreaterEqual(result.min(), expected_min - 1e-10)
        self.assertLessEqual(result.max(), expected_max + 1e-10)

    def test_oscillating_sinspaces_small_period(self):
        """Test oscillating_sinspaces with small period."""
        result = _functions.oscillating_sinspaces(
            amps=self.small_period_amps,
            periods=self.small_period_periods,
            phases=self.small_period_phases,
            bases=self.small_period_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify that small period results in multiple oscillations.
        # Count zero crossings as a proxy for oscillations.
        zero_crossings = np.sum(np.diff(np.sign(result - self.small_period_bases)) != 0)
        self.assertGreater(zero_crossings, 10)

    def test_oscillating_sinspaces_small_num_steps(self):
        """Test oscillating_sinspaces with small num_steps."""
        result = _functions.oscillating_sinspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.small_num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        self.assertEqual(result.shape, (self.small_num_steps,))

    def test_oscillating_sinspaces_large_num_steps(self):
        """Test oscillating_sinspaces with large num_steps."""
        result = _functions.oscillating_sinspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.large_num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        self.assertEqual(result.shape, (self.large_num_steps,))

    def test_oscillating_sinspaces_negative_phase(self):
        """Test oscillating_sinspaces with negative phase."""
        result = _functions.oscillating_sinspaces(
            amps=self.negative_phase_amps,
            periods=self.negative_phase_periods,
            phases=self.negative_phase_phases,
            bases=self.negative_phase_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify that negative phase shifts the waveform backward.
        # At t=0, sine with -90 degree phase should equal negative amplitude.
        npt.assert_allclose(
            result[0], -self.negative_phase_amps, rtol=1e-10, atol=1e-14
        )

    def test_oscillating_sinspaces_max_phase(self):
        """Test oscillating_sinspaces with maximum phase."""
        result = _functions.oscillating_sinspaces(
            amps=self.max_phase_amps,
            periods=self.max_phase_periods,
            phases=self.max_phase_phases,
            bases=self.max_phase_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify that max phase (180 degrees) inverts the waveform.
        # At t=0, sine with 180 degree phase should be approximately 0.
        npt.assert_allclose(result[0], 0.0, rtol=1e-10, atol=1e-14)

    def test_oscillating_sinspaces_min_phase(self):
        """Test oscillating_sinspaces with minimum phase."""
        result = _functions.oscillating_sinspaces(
            amps=self.min_phase_amps,
            periods=self.min_phase_periods,
            phases=self.min_phase_phases,
            bases=self.min_phase_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify result is computed without error.
        self.assertEqual(result.shape, (self.num_steps,))

    def test_oscillating_linspaces_scalar_parameters(self):
        """Test oscillating_linspaces with scalar parameters."""
        result = _functions.oscillating_linspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        self.assertEqual(result.shape, (self.num_steps,))

        # Verify output values match expected triangular wave.
        times = np.linspace(
            0, self.num_steps * self.delta_time, self.num_steps, endpoint=False
        )
        expected = (
            self.scalar_amps
            * signal.sawtooth(
                2 * np.pi * times / self.scalar_periods
                + (np.pi / 2)
                + np.deg2rad(self.scalar_phases),
                0.5,
            )
            + self.scalar_bases
        )

        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_linspaces_array_parameters(self):
        """Test oscillating_linspaces with array parameters."""
        result = _functions.oscillating_linspaces(
            amps=self.array_amps,
            periods=self.array_periods,
            phases=self.array_phases,
            bases=self.array_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        expected_shape = (3, self.num_steps)
        self.assertEqual(result.shape, expected_shape)

        # Verify output values for each element.
        times = np.linspace(
            0, self.num_steps * self.delta_time, self.num_steps, endpoint=False
        )
        for i in range(3):
            expected = (
                self.array_amps[i]
                * signal.sawtooth(
                    2 * np.pi * times / self.array_periods[i]
                    + (np.pi / 2)
                    + np.deg2rad(self.array_phases[i]),
                    0.5,
                )
                + self.array_bases[i]
            )
            npt.assert_allclose(result[i, :], expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_linspaces_static_parameters(self):
        """Test oscillating_linspaces with static parameters."""
        result = _functions.oscillating_linspaces(
            amps=self.static_amps,
            periods=self.static_periods,
            phases=self.static_phases,
            bases=self.static_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output is constant and equal to base.
        expected = np.full(self.num_steps, self.static_bases, dtype=float)
        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_linspaces_mixed_static_parameters(self):
        """Test oscillating_linspaces with mixed static and dynamic parameters."""
        result = _functions.oscillating_linspaces(
            amps=self.mixed_static_amps,
            periods=self.mixed_static_periods,
            phases=self.mixed_static_phases,
            bases=self.mixed_static_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output shape.
        expected_shape = (3, self.num_steps)
        self.assertEqual(result.shape, expected_shape)

        # Verify that static element (index 1) is constant.
        expected_static = np.full(
            self.num_steps, self.mixed_static_bases[1], dtype=float
        )
        npt.assert_allclose(result[1, :], expected_static, rtol=1e-10, atol=1e-14)

        # Verify that dynamic elements are oscillating.
        self.assertFalse(np.allclose(result[0, :], result[0, 0]))
        self.assertFalse(np.allclose(result[2, :], result[2, 0]))

    def test_oscillating_linspaces_phase_offset(self):
        """Test oscillating_linspaces with phase offset."""
        result = _functions.oscillating_linspaces(
            amps=self.phase_offset_amps,
            periods=self.phase_offset_periods,
            phases=self.phase_offset_phases,
            bases=self.phase_offset_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify that phase offset shifts the waveform.
        # At t=0, triangular wave with 90 degree phase should be near maximum.
        self.assertGreater(result[0], 0.5 * self.phase_offset_amps)

    def test_oscillating_linspaces_large_amplitude(self):
        """Test oscillating_linspaces with large amplitude."""
        result = _functions.oscillating_linspaces(
            amps=self.large_amplitude_amps,
            periods=self.large_amplitude_periods,
            phases=self.large_amplitude_phases,
            bases=self.large_amplitude_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify output range.
        expected_min = self.large_amplitude_bases - self.large_amplitude_amps
        expected_max = self.large_amplitude_bases + self.large_amplitude_amps
        self.assertGreaterEqual(result.min(), expected_min - 1e-10)
        self.assertLessEqual(result.max(), expected_max + 1e-10)

    def test_oscillating_customspaces_scalar_parameters_with_sine(self):
        """Test oscillating_customspaces with scalar parameters and valid sine
        function."""
        result = _functions.oscillating_customspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            bases=self.scalar_bases,
            phases=self.scalar_phases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
            custom_function=self.valid_custom_sine,
        )

        # Verify output shape.
        self.assertEqual(result.shape, (self.num_steps,))

        # Verify result matches oscillating_sinspaces.
        expected = _functions.oscillating_sinspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )
        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_customspaces_array_parameters_with_triangle(self):
        """Test oscillating_customspaces with array parameters and valid triangle
        function."""
        result = _functions.oscillating_customspaces(
            amps=self.array_amps,
            periods=self.array_periods,
            bases=self.array_bases,
            phases=self.array_phases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
            custom_function=self.valid_custom_triangle,
        )

        # Verify output shape.
        expected_shape = (3, self.num_steps)
        self.assertEqual(result.shape, expected_shape)

    def test_oscillating_customspaces_with_harmonic(self):
        """Test oscillating_customspaces with valid harmonic function."""
        result = _functions.oscillating_customspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            bases=self.scalar_bases,
            phases=self.scalar_phases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
            custom_function=self.valid_custom_harmonic,
        )

        # Verify output shape.
        self.assertEqual(result.shape, (self.num_steps,))

    def test_oscillating_customspaces_static_parameters(self):
        """Test oscillating_customspaces with static parameters."""
        result = _functions.oscillating_customspaces(
            amps=self.static_amps,
            periods=self.static_periods,
            bases=self.static_bases,
            phases=self.static_phases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
            custom_function=self.valid_custom_sine,
        )

        # Verify output is constant and equal to base.
        expected = np.full(self.num_steps, self.static_bases, dtype=float)
        npt.assert_allclose(result, expected, rtol=1e-10, atol=1e-14)

    def test_oscillating_customspaces_invalid_function_wrong_start(self):
        """Test oscillating_customspaces with invalid function that doesn't start at 0."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_wrong_start,
            )

        self.assertIn("must start at 0", str(context.exception))

    def test_oscillating_customspaces_invalid_function_wrong_end(self):
        """Test oscillating_customspaces with invalid function that doesn't return to 0
        after one period."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_wrong_end,
            )

        self.assertIn("must return to 0 after one period", str(context.exception))

    def test_oscillating_customspaces_invalid_function_wrong_amplitude(self):
        """Test oscillating_customspaces with invalid function with amplitude not equal
        to 1."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_wrong_amplitude,
            )

        self.assertIn("must have amplitude of 1", str(context.exception))

    def test_oscillating_customspaces_invalid_function_not_periodic(self):
        """Test oscillating_customspaces with invalid function that is not periodic."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_not_periodic,
            )

        self.assertIn("must be periodic", str(context.exception))

    def test_oscillating_customspaces_invalid_function_returns_nan(self):
        """Test oscillating_customspaces with invalid function that returns NaN."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_returns_nan,
            )

        self.assertIn("finite values only", str(context.exception))

    def test_oscillating_customspaces_invalid_function_returns_inf(self):
        """Test oscillating_customspaces with invalid function that returns Inf."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_returns_inf,
            )

        self.assertIn("finite values only", str(context.exception))

    def test_oscillating_customspaces_invalid_function_wrong_shape(self):
        """Test oscillating_customspaces with invalid function that returns wrong
        shape."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_customspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                bases=self.scalar_bases,
                phases=self.scalar_phases,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
                custom_function=self.invalid_custom_wrong_shape,
            )

        self.assertIn("same shape", str(context.exception))

    def test_oscillating_sinspaces_and_linspaces_different_outputs(self):
        """Test that oscillating_sinspaces and oscillating_linspaces produce different
        outputs for the same parameters."""
        sinspaces_result = _functions.oscillating_sinspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        linspaces_result = _functions.oscillating_linspaces(
            amps=self.scalar_amps,
            periods=self.scalar_periods,
            phases=self.scalar_phases,
            bases=self.scalar_bases,
            num_steps=self.num_steps,
            delta_time=self.delta_time,
        )

        # Verify that the outputs are different (except at specific points).
        self.assertFalse(np.allclose(sinspaces_result, linspaces_result))

    def test_oscillating_functions_validation_invalid_amps_periods_mismatch(self):
        """Test that oscillating functions raise error when amps is zero but periods is
        not."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_sinspaces(
                amps=0.0,
                periods=1.0,
                phases=0.0,
                bases=0.0,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
            )

        self.assertIn("If an element in amps is 0.0", str(context.exception))

    def test_oscillating_functions_validation_invalid_phase_for_static(self):
        """Test that oscillating functions raise error when amps and periods are zero
        but phases is not."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_sinspaces(
                amps=0.0,
                periods=0.0,
                phases=90.0,
                bases=0.0,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
            )

        self.assertIn(
            "corresponding element in phases must also be 0.0", str(context.exception)
        )

    def test_oscillating_functions_validation_invalid_num_steps(self):
        """Test that oscillating functions raise error with invalid num_steps."""
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            _functions.oscillating_sinspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                phases=self.scalar_phases,
                bases=self.scalar_bases,
                num_steps=0,
                delta_time=self.delta_time,
            )

    def test_oscillating_functions_validation_invalid_delta_time(self):
        """Test that oscillating functions raise error with invalid delta_time."""
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            _functions.oscillating_sinspaces(
                amps=self.scalar_amps,
                periods=self.scalar_periods,
                phases=self.scalar_phases,
                bases=self.scalar_bases,
                num_steps=self.num_steps,
                delta_time=0.0,
            )

    def test_oscillating_functions_validation_mismatched_array_shapes(self):
        """Test that oscillating functions raise error with mismatched array shapes."""
        with self.assertRaises(ValueError) as context:
            _functions.oscillating_sinspaces(
                amps=np.array([1.0, 2.0], dtype=float),
                periods=np.array([1.0, 2.0, 3.0], dtype=float),
                phases=0.0,
                bases=0.0,
                num_steps=self.num_steps,
                delta_time=self.delta_time,
            )

        self.assertIn("must have the same shape", str(context.exception))


if __name__ == "__main__":
    unittest.main()
