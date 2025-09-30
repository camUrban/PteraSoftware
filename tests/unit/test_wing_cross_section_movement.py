"""This module contains a class to test WingCrossSectionMovements.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    TestWingCrossSectionMovement: This is a class with functions to test
    WingCrossSectionMovements.
"""

import unittest
import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps
from tests.unit.fixtures import movement_fixtures


class TestWingCrossSectionMovement(unittest.TestCase):
    """This is a class with functions to test WingCrossSectionMovements."""

    def setUp(self):
        """Set up test fixtures for WingCrossSectionMovement tests."""
        self.sine_spacing_Lp_movement = (
            movement_fixtures.make_sine_spacing_Lp_wing_cross_section_movement_fixture()
        )
        self.uniform_spacing_Lp_movement = (
            movement_fixtures.make_uniform_spacing_Lp_wing_cross_section_movement_fixture()
        )
        self.mixed_spacing_Lp_movement = (
            movement_fixtures.make_mixed_spacing_Lp_wing_cross_section_movement_fixture()
        )
        self.sine_spacing_angles_movement = (
            movement_fixtures.make_sine_spacing_angles_wing_cross_section_movement_fixture()
        )
        self.uniform_spacing_angles_movement = (
            movement_fixtures.make_uniform_spacing_angles_wing_cross_section_movement_fixture()
        )
        self.mixed_spacing_angles_movement = (
            movement_fixtures.make_mixed_spacing_angles_wing_cross_section_movement_fixture()
        )

    def tearDown(self):
        """Clean up test fixtures."""
        del self.sine_spacing_Lp_movement
        del self.uniform_spacing_Lp_movement
        del self.mixed_spacing_Lp_movement
        del self.sine_spacing_angles_movement
        del self.uniform_spacing_angles_movement
        del self.mixed_spacing_angles_movement

    def test_spacing_sine_for_Lp_Wcsp_Lpp(self):
        """Test that sine spacing actually produces sinusoidal motion for
        Lp_Wcsp_Lpp."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = self.sine_spacing_Lp_movement.generate_wing_cross_sections(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated WingCrossSections.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 1.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected sine wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_Lp_Wcsp_Lpp(self):
        """Test that uniform spacing actually produces triangular wave motion for
        Lp_Wcsp_Lpp."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.uniform_spacing_Lp_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract x-positions from generated WingCrossSections.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 1.0 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)

        # Assert that the generated positions match the expected triangular wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_Lp_Wcsp_Lpp(self):
        """Test that mixed spacing types work correctly for Lp_Wcsp_Lpp."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = self.mixed_spacing_Lp_movement.generate_wing_cross_sections(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions from generated WingCrossSections.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])
        y_positions = np.array([wcs.Lp_Wcsp_Lpp[1] for wcs in wing_cross_sections])
        z_positions = np.array([wcs.Lp_Wcsp_Lpp[2] for wcs in wing_cross_sections])

        # Calculate expected values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.5 + 1.0 * np.sin(2 * np.pi * times / 1.0)
        expected_y = 2.0 + 1.5 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)
        expected_z = 0.2 + 0.5 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected values.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_spacing_sine_for_angles_Wcsp_to_Wcs_izyx(self):
        """Test that sine spacing actually produces sinusoidal motion for
        angles_Wcsp_to_Wcs_izyx."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.sine_spacing_angles_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract angles from generated WingCrossSections.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_izyx[0] for wcs in wing_cross_sections]
        )

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_angles = 10.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated angles match the expected sine wave.
        npt.assert_allclose(angles_z, expected_angles, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_angles_Wcsp_to_Wcs_izyx(self):
        """Test that uniform spacing actually produces triangular wave motion for
        angles_Wcsp_to_Wcs_izyx."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.uniform_spacing_angles_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract angles from generated WingCrossSections.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_izyx[0] for wcs in wing_cross_sections]
        )

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_angles = 10.0 * signal.sawtooth(
            2 * np.pi * times / 1.0 + np.pi / 2, 0.5
        )

        # Assert that the generated angles match the expected triangular wave.
        npt.assert_allclose(angles_z, expected_angles, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_angles_Wcsp_to_Wcs_izyx(self):
        """Test that mixed spacing types work correctly for angles_Wcsp_to_Wcs_izyx."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.mixed_spacing_angles_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract angles from generated WingCrossSections.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_izyx[0] for wcs in wing_cross_sections]
        )
        angles_y = np.array(
            [wcs.angles_Wcsp_to_Wcs_izyx[1] for wcs in wing_cross_sections]
        )
        angles_x = np.array(
            [wcs.angles_Wcsp_to_Wcs_izyx[2] for wcs in wing_cross_sections]
        )

        # Calculate expected values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_angles_z = 10.0 * np.sin(2 * np.pi * times / 1.0)
        expected_angles_y = 20.0 * signal.sawtooth(
            2 * np.pi * times / 1.0 + np.pi / 2, 0.5
        )
        expected_angles_x = 5.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated angles match the expected values.
        npt.assert_allclose(angles_z, expected_angles_z, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(angles_y, expected_angles_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(angles_x, expected_angles_x, rtol=1e-10, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
