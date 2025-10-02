"""This module contains a class to test WingMovements."""

import unittest
import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps

from tests.unit.fixtures import geometry_fixtures
from tests.unit.fixtures import wing_cross_section_movement_fixtures
from tests.unit.fixtures import wing_movement_fixtures


class TestWingMovement(unittest.TestCase):
    """This is a class with functions to test WingMovements."""

    def setUp(self):
        """Set up test fixtures for WingMovement tests."""
        # Spacing test fixtures for prelimLer_G_Cg.
        self.sine_spacing_prelimLer_wing_movement = (
            wing_movement_fixtures.make_sine_spacing_prelimLer_wing_movement_fixture()
        )
        self.uniform_spacing_prelimLer_wing_movement = (
            wing_movement_fixtures.make_uniform_spacing_prelimLer_wing_movement_fixture()
        )
        self.mixed_spacing_prelimLer_wing_movement = (
            wing_movement_fixtures.make_mixed_spacing_prelimLer_wing_movement_fixture()
        )

        # Spacing test fixtures for angles_G_to_prelimWn_izyx.
        self.sine_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_sine_spacing_angles_wing_movement_fixture()
        )
        self.uniform_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_uniform_spacing_angles_wing_movement_fixture()
        )
        self.mixed_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_mixed_spacing_angles_wing_movement_fixture()
        )

        # Additional test fixtures.
        self.static_wing_movement = (
            wing_movement_fixtures.make_static_wing_movement_fixture()
        )
        self.basic_wing_movement = (
            wing_movement_fixtures.make_basic_wing_movement_fixture()
        )
        self.prelimLer_only_wing_movement = (
            wing_movement_fixtures.make_prelimLer_only_wing_movement_fixture()
        )
        self.angles_only_wing_movement = (
            wing_movement_fixtures.make_angles_only_wing_movement_fixture()
        )
        self.phase_offset_prelimLer_wing_movement = (
            wing_movement_fixtures.make_phase_offset_prelimLer_wing_movement_fixture()
        )
        self.phase_offset_angles_wing_movement = (
            wing_movement_fixtures.make_phase_offset_angles_wing_movement_fixture()
        )
        self.multiple_periods_wing_movement = (
            wing_movement_fixtures.make_multiple_periods_wing_movement_fixture()
        )
        self.custom_spacing_prelimLer_wing_movement = (
            wing_movement_fixtures.make_custom_spacing_prelimLer_wing_movement_fixture()
        )
        self.custom_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_custom_spacing_angles_wing_movement_fixture()
        )
        self.mixed_custom_and_standard_spacing_wing_movement = (
            wing_movement_fixtures.make_mixed_custom_and_standard_spacing_wing_movement_fixture()
        )

    def test_spacing_sine_for_prelimLer_G_Cg(self):
        """Test that sine spacing actually produces sinusoidal motion for
        prelimLer_G_Cg."""
        num_steps = 100
        delta_time = 0.01
        wings = self.sine_spacing_prelimLer_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated Wings.
        x_positions = np.array([wing.prelimLer_G_Cg[0] for wing in wings])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.2 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected sine wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_prelimLer_G_Cg(self):
        """Test that uniform spacing actually produces triangular wave motion for
        prelimLer_G_Cg."""
        num_steps = 100
        delta_time = 0.01
        wings = self.uniform_spacing_prelimLer_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated Wings.
        x_positions = np.array([wing.prelimLer_G_Cg[0] for wing in wings])

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.2 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)

        # Assert that the generated positions match the expected triangular wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_prelimLer_G_Cg(self):
        """Test that mixed spacing types work correctly for prelimLer_G_Cg."""
        num_steps = 100
        delta_time = 0.01
        wings = self.mixed_spacing_prelimLer_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions from generated Wings.
        x_positions = np.array([wing.prelimLer_G_Cg[0] for wing in wings])
        y_positions = np.array([wing.prelimLer_G_Cg[1] for wing in wings])
        z_positions = np.array([wing.prelimLer_G_Cg[2] for wing in wings])

        # Calculate expected values for each dimension.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.2 * np.sin(2 * np.pi * times / 1.0)
        expected_y = 0.15 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)
        expected_z = 0.1 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected values.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_spacing_sine_for_angles_G_to_prelimWn_izyx(self):
        """Test that sine spacing actually produces sinusoidal motion for
        angles_G_to_prelimWn_izyx."""
        num_steps = 100
        delta_time = 0.01
        wings = self.sine_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-angles from generated Wings.
        x_angles = np.array([wing.angles_G_to_prelimWn_izyx[0] for wing in wings])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated angles match the expected sine wave.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_angles_G_to_prelimWn_izyx(self):
        """Test that uniform spacing actually produces triangular wave motion for
        angles_G_to_prelimWn_izyx."""
        num_steps = 100
        delta_time = 0.01
        wings = self.uniform_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-angles from generated Wings.
        x_angles = np.array([wing.angles_G_to_prelimWn_izyx[0] for wing in wings])

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)

        # Assert that the generated angles match the expected triangular wave.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_angles_G_to_prelimWn_izyx(self):
        """Test that mixed spacing types work correctly for
        angles_G_to_prelimWn_izyx."""
        num_steps = 100
        delta_time = 0.01
        wings = self.mixed_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract angles from generated Wings.
        x_angles = np.array([wing.angles_G_to_prelimWn_izyx[0] for wing in wings])
        y_angles = np.array([wing.angles_G_to_prelimWn_izyx[1] for wing in wings])
        z_angles = np.array([wing.angles_G_to_prelimWn_izyx[2] for wing in wings])

        # Calculate expected values for each dimension.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * np.sin(2 * np.pi * times / 1.0)
        expected_y = 15.0 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)
        expected_z = 8.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated angles match the expected values.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_angles, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_angles, expected_z, rtol=1e-10, atol=1e-14)

    def test_static_wing_movement_produces_constant_wings(self):
        """Test that static WingMovement produces Wings with constant parameters."""
        num_steps = 50
        delta_time = 0.02
        wings = self.static_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract parameters from all Wings.
        prelimLers_G_Cg = np.array([wing.prelimLer_G_Cg for wing in wings])
        angles_G_to_prelimWn_izyxs = np.array(
            [wing.angles_G_to_prelimWn_izyx for wing in wings]
        )

        # Assert that all Wings have the same parameters.
        npt.assert_allclose(
            prelimLers_G_Cg,
            np.tile(wings[0].prelimLer_G_Cg, (num_steps, 1)),
            rtol=1e-10,
            atol=1e-14,
        )
        npt.assert_allclose(
            angles_G_to_prelimWn_izyxs,
            np.tile(wings[0].angles_G_to_prelimWn_izyx, (num_steps, 1)),
            rtol=1e-10,
            atol=1e-14,
        )

    def test_generate_wings_returns_correct_number(self):
        """Test that generate_wings returns the correct number of Wings."""
        num_steps = 75
        delta_time = 0.015
        wings = self.basic_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        self.assertEqual(len(wings), num_steps)

    def test_generate_wings_preserves_wing_properties(self):
        """Test that generate_wings preserves non-changing Wing properties."""
        num_steps = 30
        delta_time = 0.02
        wings = self.basic_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Check that all Wings have the same non-changing properties.
        base_wing = self.basic_wing_movement.base_wing
        for wing in wings:
            self.assertEqual(wing.name, base_wing.name)
            self.assertEqual(wing.symmetric, base_wing.symmetric)
            self.assertEqual(wing.mirror_only, base_wing.mirror_only)
            npt.assert_array_equal(
                wing.symmetry_normal_Wn, base_wing.symmetry_normal_Wn
            )
            npt.assert_array_equal(
                wing.symmetry_point_Wn_Ler, base_wing.symmetry_point_Wn_Ler
            )
            self.assertEqual(wing.num_chordwise_panels, base_wing.num_chordwise_panels)
            self.assertEqual(wing.chordwise_spacing, base_wing.chordwise_spacing)
            self.assertEqual(
                len(wing.wing_cross_sections), len(base_wing.wing_cross_sections)
            )

    def test_phase_offset_prelimLer_produces_shifted_motion(self):
        """Test that phase offset for prelimLer_G_Cg produces phase-shifted
        motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.phase_offset_prelimLer_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions from generated Wings.
        x_positions = np.array([wing.prelimLer_G_Cg[0] for wing in wings])
        y_positions = np.array([wing.prelimLer_G_Cg[1] for wing in wings])
        z_positions = np.array([wing.prelimLer_G_Cg[2] for wing in wings])

        # Calculate expected phase-shifted sine waves.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.1 * np.sin(2 * np.pi * times / 1.0 + np.deg2rad(90.0))
        expected_y = 0.08 * np.sin(2 * np.pi * times / 1.0 + np.deg2rad(-45.0))
        expected_z = 0.06 * np.sin(2 * np.pi * times / 1.0 + np.deg2rad(60.0))

        # Assert that the generated positions match the expected phase-shifted waves.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_phase_offset_angles_produces_shifted_motion(self):
        """Test that phase offset for angles_G_to_prelimWn_izyx produces
        phase-shifted motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.phase_offset_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract angles from generated Wings.
        x_angles = np.array([wing.angles_G_to_prelimWn_izyx[0] for wing in wings])
        y_angles = np.array([wing.angles_G_to_prelimWn_izyx[1] for wing in wings])
        z_angles = np.array([wing.angles_G_to_prelimWn_izyx[2] for wing in wings])

        # Calculate expected phase-shifted sine waves.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * np.sin(2 * np.pi * times / 1.0 + np.deg2rad(45.0))
        expected_y = 12.0 * np.sin(2 * np.pi * times / 1.0 + np.deg2rad(90.0))
        expected_z = 8.0 * np.sin(2 * np.pi * times / 1.0 + np.deg2rad(-30.0))

        # Assert that the generated angles match the expected phase-shifted waves.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_angles, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_angles, expected_z, rtol=1e-10, atol=1e-14)

    def test_max_period_static_movement(self):
        """Test that max_period returns 0.0 for static WingMovement."""
        max_period = self.static_wing_movement.max_period
        self.assertEqual(max_period, 0.0)

    def test_max_period_prelimLer_only_movement(self):
        """Test that max_period correctly identifies the maximum period for
        prelimLer-only WingMovement."""
        max_period = self.prelimLer_only_wing_movement.max_period
        self.assertEqual(max_period, 1.5)

    def test_max_period_angles_only_movement(self):
        """Test that max_period correctly identifies the maximum period for
        angles-only WingMovement."""
        max_period = self.angles_only_wing_movement.max_period
        self.assertEqual(max_period, 1.5)

    def test_max_period_multiple_periods_movement(self):
        """Test that max_period correctly identifies the maximum period when
        different dimensions have different periods."""
        max_period = self.multiple_periods_wing_movement.max_period

        # The maximum should be from either the WingMovement's own motion or from
        # WingCrossSectionMovements.
        expected_max = max(3.0, 2.5, 2.0)
        self.assertEqual(max_period, expected_max)

    def test_initialization_with_valid_parameters(self):
        """Test WingMovement initialization with valid parameters."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampPrelimLer_G_Cg=(0.1, 0.05, 0.02),
            periodPrelimLer_G_Cg=(1.0, 1.0, 1.0),
            spacingPrelimLer_G_Cg=("sine", "uniform", "sine"),
            phasePrelimLer_G_Cg=(0.0, 45.0, -30.0),
            ampAngles_G_to_prelimWn_izyx=(5.0, 3.0, 2.0),
            periodAngles_G_to_prelimWn_izyx=(1.0, 1.0, 1.0),
            spacingAngles_G_to_prelimWn_izyx=("uniform", "sine", "uniform"),
            phaseAngles_G_to_prelimWn_izyx=(30.0, 0.0, -45.0),
        )

        self.assertIsInstance(wing_movement, ps.movements.wing_movement.WingMovement)
        self.assertEqual(wing_movement.base_wing, base_wing)
        self.assertEqual(
            len(wing_movement.wing_cross_section_movements),
            len(base_wing.wing_cross_sections),
        )
        npt.assert_array_equal(
            wing_movement.ampPrelimLer_G_Cg, np.array([0.1, 0.05, 0.02])
        )
        npt.assert_array_equal(
            wing_movement.periodPrelimLer_G_Cg, np.array([1.0, 1.0, 1.0])
        )
        self.assertEqual(
            wing_movement.spacingPrelimLer_G_Cg, ("sine", "uniform", "sine")
        )
        npt.assert_array_equal(
            wing_movement.phasePrelimLer_G_Cg, np.array([0.0, 45.0, -30.0])
        )
        npt.assert_array_equal(
            wing_movement.ampAngles_G_to_prelimWn_izyx, np.array([5.0, 3.0, 2.0])
        )
        npt.assert_array_equal(
            wing_movement.periodAngles_G_to_prelimWn_izyx, np.array([1.0, 1.0, 1.0])
        )
        self.assertEqual(
            wing_movement.spacingAngles_G_to_prelimWn_izyx,
            ("uniform", "sine", "uniform"),
        )
        npt.assert_array_equal(
            wing_movement.phaseAngles_G_to_prelimWn_izyx, np.array([30.0, 0.0, -45.0])
        )

        del base_wing
        del wcs_movements
        del wing_movement

    def test_initialization_invalid_base_wing(self):
        """Test that WingMovement initialization fails with invalid base_wing."""
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
        ]

        with self.assertRaises(TypeError):
            ps.movements.wing_movement.WingMovement(
                base_wing="not_a_wing",
                wing_cross_section_movements=wcs_movements,
            )

        del wcs_movements

    def test_initialization_invalid_wing_cross_section_movements_type(self):
        """Test that WingMovement initialization fails with invalid
        wing_cross_section_movements type."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()

        with self.assertRaises(TypeError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements="not_a_list",
            )

        del base_wing

    def test_initialization_invalid_wing_cross_section_movements_length(self):
        """Test that WingMovement initialization fails when
        wing_cross_section_movements length doesn't match base_wing."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
        ]

        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
            )

        del base_wing
        del wcs_movements

    def test_initialization_ampPrelimLer_G_Cg_validation(self):
        """Test ampPrelimLer_G_Cg parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with negative amplitude.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampPrelimLer_G_Cg=(-0.1, 0.0, 0.0),
                periodPrelimLer_G_Cg=(1.0, 0.0, 0.0),
            )

        del base_wing
        del wcs_movements

    def test_initialization_periodPrelimLer_G_Cg_validation(self):
        """Test periodPrelimLer_G_Cg parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with zero amplitude but non-zero period.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                periodPrelimLer_G_Cg=(1.0, 0.0, 0.0),
            )

        del base_wing
        del wcs_movements

    def test_initialization_phasePrelimLer_G_Cg_validation(self):
        """Test phasePrelimLer_G_Cg parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with phase out of valid range.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampPrelimLer_G_Cg=(0.1, 0.0, 0.0),
                periodPrelimLer_G_Cg=(1.0, 0.0, 0.0),
                phasePrelimLer_G_Cg=(181.0, 0.0, 0.0),
            )

        # Test with zero amplitude but non-zero phase.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                periodPrelimLer_G_Cg=(0.0, 0.0, 0.0),
                phasePrelimLer_G_Cg=(45.0, 0.0, 0.0),
            )

        del base_wing
        del wcs_movements

    def test_initialization_ampAngles_G_to_prelimWn_izyx_validation(self):
        """Test ampAngles_G_to_prelimWn_izyx parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with amplitude >= 180 degrees.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_G_to_prelimWn_izyx=(180.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
            )

        # Test with negative amplitude.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_G_to_prelimWn_izyx=(-10.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
            )

        del base_wing
        del wcs_movements

    def test_initialization_periodAngles_G_to_prelimWn_izyx_validation(self):
        """Test periodAngles_G_to_prelimWn_izyx parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with zero amplitude but non-zero period.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
            )

        del base_wing
        del wcs_movements

    def test_initialization_phaseAngles_G_to_prelimWn_izyx_validation(self):
        """Test phaseAngles_G_to_prelimWn_izyx parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with phase out of valid range.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_G_to_prelimWn_izyx=(10.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(1.0, 0.0, 0.0),
                phaseAngles_G_to_prelimWn_izyx=(181.0, 0.0, 0.0),
            )

        # Test with zero amplitude but non-zero phase.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                periodAngles_G_to_prelimWn_izyx=(0.0, 0.0, 0.0),
                phaseAngles_G_to_prelimWn_izyx=(45.0, 0.0, 0.0),
            )

        del base_wing
        del wcs_movements

    def test_custom_spacing_prelimLer_produces_expected_motion(self):
        """Test that custom spacing function for prelimLer_G_Cg produces
        expected motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.custom_spacing_prelimLer_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated Wings.
        x_positions = np.array([wing.prelimLer_G_Cg[0] for wing in wings])

        # Calculate expected custom harmonic wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        x_rad = 2 * np.pi * times / 1.0
        expected_x = (
            0.15
            * (3.0 / (2.0 * np.sqrt(2.0)))
            * (np.sin(x_rad) + (1.0 / 3.0) * np.sin(3.0 * x_rad))
        )

        # Assert that the generated positions match the expected custom wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_custom_spacing_angles_produces_expected_motion(self):
        """Test that custom spacing function for angles_G_to_prelimWn_izyx
        produces expected motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.custom_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-angles from generated Wings.
        x_angles = np.array([wing.angles_G_to_prelimWn_izyx[0] for wing in wings])

        # Calculate expected custom harmonic wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        x_rad = 2 * np.pi * times / 1.0
        expected_x = (
            10.0
            * (3.0 / (2.0 * np.sqrt(2.0)))
            * (np.sin(x_rad) + (1.0 / 3.0) * np.sin(3.0 * x_rad))
        )

        # Assert that the generated angles match the expected custom wave.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)
