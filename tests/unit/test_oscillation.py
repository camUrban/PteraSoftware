"""This module contains a class to test oscillation functions."""

import unittest

import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _oscillation
from tests.unit.fixtures import oscillation_fixtures


class TestOscillation(unittest.TestCase):
    """This is a class with functions to test oscillation functions."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all oscillation function tests."""
        # Parameter fixtures.
        (
            cls.static_amp,
            cls.static_period,
            cls.static_phase,
            cls.static_base,
        ) = oscillation_fixtures.make_static_parameters_fixture()

        (
            cls.phase_offset_amp,
            cls.phase_offset_period,
            cls.phase_offset_phase,
            cls.phase_offset_base,
        ) = oscillation_fixtures.make_phase_offset_parameters_fixture()

        (
            cls.max_phase_amp,
            cls.max_phase_period,
            cls.max_phase_phase,
            cls.max_phase_base,
        ) = oscillation_fixtures.make_max_phase_parameters_fixture()

        # Time fixture.
        cls.time = oscillation_fixtures.make_time_fixture()

        # Custom function fixtures.
        cls.valid_custom_sine = staticmethod(
            oscillation_fixtures.make_valid_custom_sine_function_fixture()
        )

    def test_oscillating_sin_at_time_static_parameters(self):
        """Test oscillating_sin_at_time with static parameters."""
        result = _oscillation.oscillating_sin_at_time(
            amp=self.static_amp,
            period=self.static_period,
            phase=self.static_phase,
            base=self.static_base,
            time=self.time,
        )

        # Verify output is equal to base.
        npt.assert_allclose(result, self.static_base, rtol=1e-10, atol=1e-14)

    def test_oscillating_sin_at_time_phase_offset(self):
        """Test oscillating_sin_at_time with phase offset."""
        result = _oscillation.oscillating_sin_at_time(
            amp=self.phase_offset_amp,
            period=self.phase_offset_period,
            phase=self.phase_offset_phase,
            base=self.phase_offset_base,
            time=0.0,
        )

        # Verify that phase offset shifts the waveform.
        # At t=0.0, sine with 90.0 degree phase should equal the amplitude.
        npt.assert_allclose(result, self.phase_offset_amp, rtol=1e-10, atol=1e-14)

    def test_oscillating_sin_at_time_max_phase(self):
        """Test oscillating_sin_at_time with maximum phase."""
        result = _oscillation.oscillating_sin_at_time(
            amp=self.max_phase_amp,
            period=self.max_phase_period,
            phase=self.max_phase_phase,
            base=self.max_phase_base,
            time=0.0,
        )

        # Verify that max phase (180.0 degrees) inverts the waveform.
        # At t=0.0, sine with 180.0 degree phase should be approximately 0.0.
        npt.assert_allclose(result, 0.0, rtol=1e-10, atol=1e-14)

    def test_oscillating_lin_at_time_static_parameters(self):
        """Test oscillating_lin_at_time with static parameters."""
        result = _oscillation.oscillating_lin_at_time(
            amp=self.static_amp,
            period=self.static_period,
            phase=self.static_phase,
            base=self.static_base,
            time=self.time,
        )

        # Verify output is equal to base.
        npt.assert_allclose(result, self.static_base, rtol=1e-10, atol=1e-14)

    def test_oscillating_lin_at_time_phase_offset(self):
        """Test oscillating_lin_at_time with phase offset."""
        result = _oscillation.oscillating_lin_at_time(
            amp=self.phase_offset_amp,
            period=self.phase_offset_period,
            phase=self.phase_offset_phase,
            base=self.phase_offset_base,
            time=0.0,
        )

        # Verify that phase offset shifts the waveform.
        # At t=0.0, triangular wave with 90.0 degree phase should be near maximum.
        self.assertGreater(result, 0.5 * self.phase_offset_amp)

    def test_oscillating_custom_at_time_static_parameters(self):
        """Test oscillating_custom_at_time with static parameters."""
        result = _oscillation.oscillating_custom_at_time(
            amp=self.static_amp,
            period=self.static_period,
            base=self.static_base,
            phase=self.static_phase,
            time=self.time,
            custom_function=self.valid_custom_sine,
        )

        # Verify output is equal to base.
        npt.assert_allclose(result, self.static_base, rtol=1e-10, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
