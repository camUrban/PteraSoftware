"""This module contains classes to test WingMovements."""

import unittest

import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps
from tests.unit.fixtures import (
    geometry_fixtures,
    wing_cross_section_movement_fixtures,
    wing_movement_fixtures,
)


class TestWingMovement(unittest.TestCase):
    """This is a class with functions to test WingMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all WingMovement tests."""
        # Spacing test fixtures for Ler_Gs_Cgs.
        cls.sine_spacing_Ler_wing_movement = (
            wing_movement_fixtures.make_sine_spacing_Ler_wing_movement_fixture()
        )
        cls.uniform_spacing_Ler_wing_movement = (
            wing_movement_fixtures.make_uniform_spacing_Ler_wing_movement_fixture()
        )
        cls.mixed_spacing_Ler_wing_movement = (
            wing_movement_fixtures.make_mixed_spacing_Ler_wing_movement_fixture()
        )

        # Spacing test fixtures for angles_Gs_to_Wn_ixyz.
        cls.sine_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_sine_spacing_angles_wing_movement_fixture()
        )
        cls.uniform_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_uniform_spacing_angles_wing_movement_fixture()
        )
        cls.mixed_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_mixed_spacing_angles_wing_movement_fixture()
        )

        # Additional test fixtures.
        cls.static_wing_movement = (
            wing_movement_fixtures.make_static_wing_movement_fixture()
        )
        cls.basic_wing_movement = (
            wing_movement_fixtures.make_basic_wing_movement_fixture()
        )
        cls.Ler_only_wing_movement = (
            wing_movement_fixtures.make_Ler_only_wing_movement_fixture()
        )
        cls.angles_only_wing_movement = (
            wing_movement_fixtures.make_angles_only_wing_movement_fixture()
        )
        cls.phase_offset_Ler_wing_movement = (
            wing_movement_fixtures.make_phase_offset_Ler_wing_movement_fixture()
        )
        cls.phase_offset_angles_wing_movement = (
            wing_movement_fixtures.make_phase_offset_angles_wing_movement_fixture()
        )
        cls.multiple_periods_wing_movement = (
            wing_movement_fixtures.make_multiple_periods_wing_movement_fixture()
        )
        cls.custom_spacing_Ler_wing_movement = (
            wing_movement_fixtures.make_custom_spacing_Ler_wing_movement_fixture()
        )
        cls.custom_spacing_angles_wing_movement = (
            wing_movement_fixtures.make_custom_spacing_angles_wing_movement_fixture()
        )
        cls.mixed_custom_and_standard_spacing_wing_movement = (
            wing_movement_fixtures.make_mixed_custom_and_standard_spacing_wing_movement_fixture()
        )
        cls.rotation_point_offset_wing_movement = (
            wing_movement_fixtures.make_rotation_point_offset_wing_movement_fixture()
        )

    def test_spacing_sine_for_Ler_Gs_Cgs(self):
        """Test that sine spacing actually produces sinusoidal motion for
        Ler_Gs_Cgs."""
        num_steps = 100
        delta_time = 0.01
        wings = self.sine_spacing_Ler_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated Wings.
        x_positions = np.array([wing.Ler_Gs_Cgs[0] for wing in wings])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.2 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected sine wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_Ler_Gs_Cgs(self):
        """Test that uniform spacing actually produces triangular wave motion for
        Ler_Gs_Cgs."""
        num_steps = 100
        delta_time = 0.01
        wings = self.uniform_spacing_Ler_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated Wings.
        x_positions = np.array([wing.Ler_Gs_Cgs[0] for wing in wings])

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.2 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)

        # Assert that the generated positions match the expected triangular wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_Ler_Gs_Cgs(self):
        """Test that mixed spacing types work correctly for Ler_Gs_Cgs."""
        num_steps = 100
        delta_time = 0.01
        wings = self.mixed_spacing_Ler_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions from generated Wings.
        x_positions = np.array([wing.Ler_Gs_Cgs[0] for wing in wings])
        y_positions = np.array([wing.Ler_Gs_Cgs[1] for wing in wings])
        z_positions = np.array([wing.Ler_Gs_Cgs[2] for wing in wings])

        # Calculate expected values for each dimension.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.2 * np.sin(2 * np.pi * times / 1.0)
        expected_y = 0.15 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)
        expected_z = 0.1 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected values.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_spacing_sine_for_angles_Gs_to_Wn_ixyz(self):
        """Test that sine spacing actually produces sinusoidal motion for
        angles_Gs_to_Wn_ixyz."""
        num_steps = 100
        delta_time = 0.01
        wings = self.sine_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-angles from generated Wings.
        x_angles = np.array([wing.angles_Gs_to_Wn_ixyz[0] for wing in wings])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated angles match the expected sine wave.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_angles_Gs_to_Wn_ixyz(self):
        """Test that uniform spacing actually produces triangular wave motion for
        angles_Gs_to_Wn_ixyz."""
        num_steps = 100
        delta_time = 0.01
        wings = self.uniform_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-angles from generated Wings.
        x_angles = np.array([wing.angles_Gs_to_Wn_ixyz[0] for wing in wings])

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)

        # Assert that the generated angles match the expected triangular wave.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_angles_Gs_to_Wn_ixyz(self):
        """Test that mixed spacing types work correctly for
        angles_Gs_to_Wn_ixyz."""
        num_steps = 100
        delta_time = 0.01
        wings = self.mixed_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract angles from generated Wings.
        x_angles = np.array([wing.angles_Gs_to_Wn_ixyz[0] for wing in wings])
        y_angles = np.array([wing.angles_Gs_to_Wn_ixyz[1] for wing in wings])
        z_angles = np.array([wing.angles_Gs_to_Wn_ixyz[2] for wing in wings])

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
        Lers_G_Cg = np.array([wing.Ler_Gs_Cgs for wing in wings])
        angles_Gs_to_Wn_ixyzs = np.array([wing.angles_Gs_to_Wn_ixyz for wing in wings])

        # Assert that all Wings have the same parameters.
        npt.assert_allclose(
            Lers_G_Cg,
            np.tile(wings[0].Ler_Gs_Cgs, (num_steps, 1)),
            rtol=1e-10,
            atol=1e-14,
        )
        npt.assert_allclose(
            angles_Gs_to_Wn_ixyzs,
            np.tile(wings[0].angles_Gs_to_Wn_ixyz, (num_steps, 1)),
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
            npt.assert_array_equal(wing.symmetryNormal_G, base_wing.symmetryNormal_G)
            npt.assert_array_equal(
                wing.symmetryPoint_G_Cg, base_wing.symmetryPoint_G_Cg
            )
            self.assertEqual(wing.num_chordwise_panels, base_wing.num_chordwise_panels)
            self.assertEqual(wing.chordwise_spacing, base_wing.chordwise_spacing)
            self.assertEqual(
                len(wing.wing_cross_sections), len(base_wing.wing_cross_sections)
            )

    def test_phase_offset_Ler_produces_shifted_motion(self):
        """Test that phase offset for Ler_Gs_Cgs produces phase-shifted
        motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.phase_offset_Ler_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions from generated Wings.
        x_positions = np.array([wing.Ler_Gs_Cgs[0] for wing in wings])
        y_positions = np.array([wing.Ler_Gs_Cgs[1] for wing in wings])
        z_positions = np.array([wing.Ler_Gs_Cgs[2] for wing in wings])

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
        """Test that phase offset for angles_Gs_to_Wn_ixyz produces
        phase-shifted motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.phase_offset_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract angles from generated Wings.
        x_angles = np.array([wing.angles_Gs_to_Wn_ixyz[0] for wing in wings])
        y_angles = np.array([wing.angles_Gs_to_Wn_ixyz[1] for wing in wings])
        z_angles = np.array([wing.angles_Gs_to_Wn_ixyz[2] for wing in wings])

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

    def test_max_period_Ler_only_movement(self):
        """Test that max_period correctly identifies the maximum period for
        Ler-only WingMovement."""
        max_period = self.Ler_only_wing_movement.max_period
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
            ampLer_Gs_Cgs=(0.1, 0.05, 0.02),
            periodLer_Gs_Cgs=(1.0, 1.0, 1.0),
            spacingLer_Gs_Cgs=("sine", "uniform", "sine"),
            phaseLer_Gs_Cgs=(0.0, 45.0, -30.0),
            ampAngles_Gs_to_Wn_ixyz=(5.0, 3.0, 2.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 1.0, 1.0),
            spacingAngles_Gs_to_Wn_ixyz=("uniform", "sine", "uniform"),
            phaseAngles_Gs_to_Wn_ixyz=(30.0, 0.0, -45.0),
        )

        self.assertIsInstance(wing_movement, ps.movements.wing_movement.WingMovement)
        self.assertEqual(wing_movement.base_wing, base_wing)
        self.assertEqual(
            len(wing_movement.wing_cross_section_movements),
            len(base_wing.wing_cross_sections),
        )
        npt.assert_array_equal(wing_movement.ampLer_Gs_Cgs, np.array([0.1, 0.05, 0.02]))
        npt.assert_array_equal(
            wing_movement.periodLer_Gs_Cgs, np.array([1.0, 1.0, 1.0])
        )
        self.assertEqual(wing_movement.spacingLer_Gs_Cgs, ("sine", "uniform", "sine"))
        npt.assert_array_equal(
            wing_movement.phaseLer_Gs_Cgs, np.array([0.0, 45.0, -30.0])
        )
        npt.assert_array_equal(
            wing_movement.ampAngles_Gs_to_Wn_ixyz, np.array([5.0, 3.0, 2.0])
        )
        npt.assert_array_equal(
            wing_movement.periodAngles_Gs_to_Wn_ixyz, np.array([1.0, 1.0, 1.0])
        )
        self.assertEqual(
            wing_movement.spacingAngles_Gs_to_Wn_ixyz,
            ("uniform", "sine", "uniform"),
        )
        npt.assert_array_equal(
            wing_movement.phaseAngles_Gs_to_Wn_ixyz, np.array([30.0, 0.0, -45.0])
        )

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

    def test_initialization_invalid_wing_cross_section_movements_type(self):
        """Test that WingMovement initialization fails with invalid
        wing_cross_section_movements type."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()

        with self.assertRaises(TypeError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements="not_a_list",
            )

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

    def test_initialization_ampLer_Gs_Cgs_validation(self):
        """Test ampLer_Gs_Cgs parameter validation."""
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
                ampLer_Gs_Cgs=(-0.1, 0.0, 0.0),
                periodLer_Gs_Cgs=(1.0, 0.0, 0.0),
            )

    def test_initialization_periodLer_Gs_Cgs_validation(self):
        """Test periodLer_Gs_Cgs parameter validation."""
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
                ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
                periodLer_Gs_Cgs=(1.0, 0.0, 0.0),
            )

    def test_initialization_phaseLer_Gs_Cgs_validation(self):
        """Test phaseLer_Gs_Cgs parameter validation."""
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
                ampLer_Gs_Cgs=(0.1, 0.0, 0.0),
                periodLer_Gs_Cgs=(1.0, 0.0, 0.0),
                phaseLer_Gs_Cgs=(181.0, 0.0, 0.0),
            )

        # Test with zero amplitude but non-zero phase.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
                periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
                phaseLer_Gs_Cgs=(45.0, 0.0, 0.0),
            )

    def test_initialization_ampAngles_Gs_to_Wn_ixyz_validation(self):
        """Test ampAngles_Gs_to_Wn_ixyz parameter validation."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        # Test with amplitude > 180 degrees.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_Gs_to_Wn_ixyz=(180.1, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            )

        # Test with negative amplitude.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_Gs_to_Wn_ixyz=(-10.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            )

    def test_initialization_periodAngles_Gs_to_Wn_ixyz_validation(self):
        """Test periodAngles_Gs_to_Wn_ixyz parameter validation."""
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
                ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            )

    def test_initialization_phaseAngles_Gs_to_Wn_ixyz_validation(self):
        """Test phaseAngles_Gs_to_Wn_ixyz parameter validation."""
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
                ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
                phaseAngles_Gs_to_Wn_ixyz=(181.0, 0.0, 0.0),
            )

        # Test with zero amplitude but non-zero phase.
        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                ampAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                phaseAngles_Gs_to_Wn_ixyz=(45.0, 0.0, 0.0),
            )

    def test_custom_spacing_Ler_produces_expected_motion(self):
        """Test that custom spacing function for Ler_Gs_Cgs produces
        expected motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.custom_spacing_Ler_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-positions from generated Wings.
        x_positions = np.array([wing.Ler_Gs_Cgs[0] for wing in wings])

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
        """Test that custom spacing function for angles_Gs_to_Wn_ixyz
        produces expected motion."""
        num_steps = 100
        delta_time = 0.01
        wings = self.custom_spacing_angles_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract x-angles from generated Wings.
        x_angles = np.array([wing.angles_Gs_to_Wn_ixyz[0] for wing in wings])

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

    def test_rotation_point_offset_zero_matches_default(self):
        """Test that zero rotation point offset produces identical results to default
        behavior."""
        num_steps = 50
        delta_time = 0.02

        # Create two WingMovements: one with explicit zero offset, one without.
        base_wing = geometry_fixtures.make_origin_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
            wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
        ]

        movement_default = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        )

        movement_zero_offset = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=[
                wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
                wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
            ],
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            rotationPointOffset_Gs_Ler=(0.0, 0.0, 0.0),
        )

        wings_default = movement_default.generate_wings(num_steps, delta_time)
        wings_zero = movement_zero_offset.generate_wings(num_steps, delta_time)

        for i in range(num_steps):
            npt.assert_allclose(
                wings_default[i].Ler_Gs_Cgs,
                wings_zero[i].Ler_Gs_Cgs,
                rtol=1e-10,
                atol=1e-14,
            )
            npt.assert_allclose(
                wings_default[i].angles_Gs_to_Wn_ixyz,
                wings_zero[i].angles_Gs_to_Wn_ixyz,
                rtol=1e-10,
                atol=1e-14,
            )

    def test_rotation_point_offset_produces_position_adjustment(self):
        """Test that non zero rotation point offset causes position changes when
        angles oscillate."""
        num_steps = 100
        delta_time = 0.01
        wings = self.rotation_point_offset_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # With rotation about a point offset in y (0.5m in y direction from Ler),
        # z positions should vary when rotating about x axis.
        # The offset is perpendicular to the rotation axis, so position changes occur.
        z_positions = np.array([wing.Ler_Gs_Cgs[2] for wing in wings])
        y_positions = np.array([wing.Ler_Gs_Cgs[1] for wing in wings])

        # Verify that z positions are not all zero (they should oscillate).
        self.assertFalse(np.allclose(z_positions, 0.0))

        # For rotation about x axis with offset P = (0, 0.5, 0), the position
        # adjustment is (I - R) @ P where R is the rotation matrix about x.
        # The active rotation matrix for angle theta about x is:
        # R = [[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]]
        # So (I - R) @ [0, 0.5, 0] = [0, 0.5*(1 - cos(theta)), -0.5*sin(theta)]
        # Thus y_adj = 0.5*(1 - cos(theta)) and z_adj = -0.5*sin(theta)
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        angles_x_rad = np.deg2rad(10.0 * np.sin(2 * np.pi * times / 1.0))
        expected_y = 0.5 * (1.0 - np.cos(angles_x_rad))
        expected_z = -0.5 * np.sin(angles_x_rad)

        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_rotation_point_offset_preserves_angles(self):
        """Test that rotation point offset does not affect the Wing angles."""
        num_steps = 50
        delta_time = 0.02
        wings = self.rotation_point_offset_wing_movement.generate_wings(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract angles from generated Wings.
        x_angles = np.array([wing.angles_Gs_to_Wn_ixyz[0] for wing in wings])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 10.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that angles are unaffected by rotation point offset.
        npt.assert_allclose(x_angles, expected_x, rtol=1e-10, atol=1e-14)

    def test_rotation_point_offset_initialization(self):
        """Test that rotationPointOffset_Gs_Ler is correctly initialized."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
            rotationPointOffset_Gs_Ler=(0.25, 0.1, -0.05),
        )

        npt.assert_array_equal(
            wing_movement.rotationPointOffset_Gs_Ler, np.array([0.25, 0.1, -0.05])
        )

    def test_rotation_point_offset_validation_wrong_size(self):
        """Test that invalid rotationPointOffset_Gs_Ler size raises error."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        with self.assertRaises(ValueError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                rotationPointOffset_Gs_Ler=(0.1, 0.2),
            )

    def test_rotation_point_offset_validation_non_numeric(self):
        """Test that non numeric rotationPointOffset_Gs_Ler raises error."""
        base_wing = geometry_fixtures.make_type_1_wing_fixture()
        wcs_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
            for _ in base_wing.wing_cross_sections
        ]

        with self.assertRaises(TypeError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wcs_movements,
                rotationPointOffset_Gs_Ler=("a", "b", "c"),
            )

    def test_all_periods_static_movement(self):
        """Test that all_periods returns empty tuple for static WingMovement."""
        wing_movement = self.static_wing_movement
        self.assertEqual(wing_movement.all_periods, ())

    def test_all_periods_Ler_only(self):
        """Test that all_periods returns correct periods for Ler only movement."""
        wing_movement = self.Ler_only_wing_movement
        # periodLer_Gs_Cgs is (1.5, 1.5, 1.5), all non zero.
        # periodAngles_Gs_to_Wn_ixyz is (0.0, 0.0, 0.0).
        # WingCrossSectionMovements are static (all zeros).
        # Should return tuple with three 1.5 values.
        self.assertEqual(wing_movement.all_periods, (1.5, 1.5, 1.5))

    def test_all_periods_angles_only(self):
        """Test that all_periods returns correct periods for angles only movement."""
        wing_movement = self.angles_only_wing_movement
        # periodLer_Gs_Cgs is (0.0, 0.0, 0.0).
        # periodAngles_Gs_to_Wn_ixyz is (1.5, 1.5, 1.5), all non zero.
        # WingCrossSectionMovements are static (all zeros).
        # Should return tuple with three 1.5 values.
        self.assertEqual(wing_movement.all_periods, (1.5, 1.5, 1.5))

    def test_all_periods_mixed(self):
        """Test that all_periods returns all non zero periods for mixed movement."""
        wing_movement = self.multiple_periods_wing_movement
        # periodLer_Gs_Cgs is (1.0, 2.0, 3.0).
        # periodAngles_Gs_to_Wn_ixyz is (0.5, 1.5, 2.5).
        # WingCrossSectionMovements include one with multiple periods.
        # Should return tuple with all non zero values from WingCrossSectionMovements first,
        # then WingMovement's own periods.
        all_periods = wing_movement.all_periods

        # Verify WingMovement's own periods are included.
        self.assertIn(1.0, all_periods)
        self.assertIn(2.0, all_periods)
        self.assertIn(3.0, all_periods)
        self.assertIn(0.5, all_periods)
        self.assertIn(1.5, all_periods)
        self.assertIn(2.5, all_periods)

    def test_all_periods_contains_duplicates(self):
        """Test that all_periods contains duplicate periods if they appear multiple
        times.
        """
        wing_movement = self.basic_wing_movement
        all_periods = wing_movement.all_periods

        # Both periodLer_Gs_Cgs and periodAngles_Gs_to_Wn_ixyz are (2.0, 2.0, 2.0).
        # This contributes six 2.0 values. Plus WingCrossSectionMovement periods.
        # Count how many times 2.0 appears (should be at least 6 from WingMovement).
        count_2_0 = all_periods.count(2.0)
        self.assertGreaterEqual(count_2_0, 6)

    def test_all_periods_partial_movement(self):
        """Test all_periods with only some dimensions having non zero periods."""
        wing_movement = self.sine_spacing_Ler_wing_movement
        # periodLer_Gs_Cgs is (1.0, 0.0, 0.0), only first element is non zero.
        # periodAngles_Gs_to_Wn_ixyz is (0.0, 0.0, 0.0).
        # WingCrossSectionMovements are static.
        # Should return tuple with one 1.0 value.
        self.assertEqual(wing_movement.all_periods, (1.0,))


class TestWingMovementImmutability(unittest.TestCase):
    """Tests for WingMovement attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.wing_movement = wing_movement_fixtures.make_basic_wing_movement_fixture()

    def test_immutable_base_wing_property(self):
        """Test that base_wing property is read only."""
        new_wing = geometry_fixtures.make_type_1_wing_fixture()
        with self.assertRaises(AttributeError):
            self.wing_movement.base_wing = new_wing

    def test_immutable_wing_cross_section_movements_property(self):
        """Test that wing_cross_section_movements property is read only."""
        new_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
        ]
        with self.assertRaises(AttributeError):
            self.wing_movement.wing_cross_section_movements = new_movements

    def test_immutable_ampLer_Gs_Cgs_property(self):
        """Test that ampLer_Gs_Cgs property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.ampLer_Gs_Cgs = np.array([1.0, 2.0, 3.0])

    def test_immutable_ampLer_Gs_Cgs_array_read_only(self):
        """Test that ampLer_Gs_Cgs array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.ampLer_Gs_Cgs[0] = 999.0

    def test_immutable_periodLer_Gs_Cgs_property(self):
        """Test that periodLer_Gs_Cgs property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.periodLer_Gs_Cgs = np.array([1.0, 2.0, 3.0])

    def test_immutable_periodLer_Gs_Cgs_array_read_only(self):
        """Test that periodLer_Gs_Cgs array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.periodLer_Gs_Cgs[0] = 999.0

    def test_immutable_spacingLer_Gs_Cgs_property(self):
        """Test that spacingLer_Gs_Cgs property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.spacingLer_Gs_Cgs = ("uniform", "uniform", "uniform")

    def test_immutable_phaseLer_Gs_Cgs_property(self):
        """Test that phaseLer_Gs_Cgs property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.phaseLer_Gs_Cgs = np.array([45.0, 45.0, 45.0])

    def test_immutable_phaseLer_Gs_Cgs_array_read_only(self):
        """Test that phaseLer_Gs_Cgs array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.phaseLer_Gs_Cgs[0] = 999.0

    def test_immutable_ampAngles_Gs_to_Wn_ixyz_property(self):
        """Test that ampAngles_Gs_to_Wn_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.ampAngles_Gs_to_Wn_ixyz = np.array([1.0, 2.0, 3.0])

    def test_immutable_ampAngles_Gs_to_Wn_ixyz_array_read_only(self):
        """Test that ampAngles_Gs_to_Wn_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.ampAngles_Gs_to_Wn_ixyz[0] = 999.0

    def test_immutable_periodAngles_Gs_to_Wn_ixyz_property(self):
        """Test that periodAngles_Gs_to_Wn_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.periodAngles_Gs_to_Wn_ixyz = np.array([1.0, 2.0, 3.0])

    def test_immutable_periodAngles_Gs_to_Wn_ixyz_array_read_only(self):
        """Test that periodAngles_Gs_to_Wn_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.periodAngles_Gs_to_Wn_ixyz[0] = 999.0

    def test_immutable_spacingAngles_Gs_to_Wn_ixyz_property(self):
        """Test that spacingAngles_Gs_to_Wn_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.spacingAngles_Gs_to_Wn_ixyz = (
                "uniform",
                "uniform",
                "uniform",
            )

    def test_immutable_phaseAngles_Gs_to_Wn_ixyz_property(self):
        """Test that phaseAngles_Gs_to_Wn_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.phaseAngles_Gs_to_Wn_ixyz = np.array([45.0, 45.0, 45.0])

    def test_immutable_phaseAngles_Gs_to_Wn_ixyz_array_read_only(self):
        """Test that phaseAngles_Gs_to_Wn_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.phaseAngles_Gs_to_Wn_ixyz[0] = 999.0

    def test_immutable_rotationPointOffset_Gs_Ler_property(self):
        """Test that rotationPointOffset_Gs_Ler property is read only."""
        with self.assertRaises(AttributeError):
            self.wing_movement.rotationPointOffset_Gs_Ler = np.array([1.0, 2.0, 3.0])

    def test_immutable_rotationPointOffset_Gs_Ler_array_read_only(self):
        """Test that rotationPointOffset_Gs_Ler array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.wing_movement.rotationPointOffset_Gs_Ler[0] = 999.0


class TestWingMovementCaching(unittest.TestCase):
    """Tests for WingMovement caching behavior."""

    def setUp(self):
        """Set up test fixtures for caching tests."""
        self.wing_movement = wing_movement_fixtures.make_basic_wing_movement_fixture()

    def test_all_periods_caching_returns_same_object(self):
        """Test that repeated access to all_periods returns the same cached object."""
        all_periods_1 = self.wing_movement.all_periods
        all_periods_2 = self.wing_movement.all_periods
        self.assertIs(all_periods_1, all_periods_2)

    def test_max_period_caching_returns_same_value(self):
        """Test that repeated access to max_period returns the same cached value."""
        max_period_1 = self.wing_movement.max_period
        max_period_2 = self.wing_movement.max_period
        # Since floats are immutable, we check equality rather than identity.
        self.assertEqual(max_period_1, max_period_2)


class TestWingMovementDeepcopy(unittest.TestCase):
    """Tests for WingMovement deepcopy behavior."""

    def setUp(self):
        """Set up test fixtures for deepcopy tests."""
        self.wing_movement = wing_movement_fixtures.make_basic_wing_movement_fixture()

    def test_deepcopy_returns_new_instance(self):
        """Test that deepcopy returns a new WingMovement instance."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.movements.wing_movement.WingMovement)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_attribute_values(self):
        """Test that deepcopy preserves all attribute values."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        # Check numpy array attributes.
        npt.assert_array_equal(copied.ampLer_Gs_Cgs, original.ampLer_Gs_Cgs)
        npt.assert_array_equal(copied.periodLer_Gs_Cgs, original.periodLer_Gs_Cgs)
        npt.assert_array_equal(copied.phaseLer_Gs_Cgs, original.phaseLer_Gs_Cgs)
        npt.assert_array_equal(
            copied.ampAngles_Gs_to_Wn_ixyz, original.ampAngles_Gs_to_Wn_ixyz
        )
        npt.assert_array_equal(
            copied.periodAngles_Gs_to_Wn_ixyz, original.periodAngles_Gs_to_Wn_ixyz
        )
        npt.assert_array_equal(
            copied.phaseAngles_Gs_to_Wn_ixyz, original.phaseAngles_Gs_to_Wn_ixyz
        )
        npt.assert_array_equal(
            copied.rotationPointOffset_Gs_Ler, original.rotationPointOffset_Gs_Ler
        )

        # Check tuple attributes.
        self.assertEqual(copied.spacingLer_Gs_Cgs, original.spacingLer_Gs_Cgs)
        self.assertEqual(
            copied.spacingAngles_Gs_to_Wn_ixyz,
            original.spacingAngles_Gs_to_Wn_ixyz,
        )

    def test_deepcopy_numpy_arrays_are_independent(self):
        """Test that deepcopied numpy arrays are independent objects."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        # Verify arrays are different objects.
        self.assertIsNot(copied.ampLer_Gs_Cgs, original.ampLer_Gs_Cgs)
        self.assertIsNot(copied.periodLer_Gs_Cgs, original.periodLer_Gs_Cgs)
        self.assertIsNot(copied.phaseLer_Gs_Cgs, original.phaseLer_Gs_Cgs)
        self.assertIsNot(
            copied.ampAngles_Gs_to_Wn_ixyz, original.ampAngles_Gs_to_Wn_ixyz
        )
        self.assertIsNot(
            copied.periodAngles_Gs_to_Wn_ixyz, original.periodAngles_Gs_to_Wn_ixyz
        )
        self.assertIsNot(
            copied.phaseAngles_Gs_to_Wn_ixyz, original.phaseAngles_Gs_to_Wn_ixyz
        )
        self.assertIsNot(
            copied.rotationPointOffset_Gs_Ler, original.rotationPointOffset_Gs_Ler
        )

    def test_deepcopy_numpy_arrays_cannot_be_modified_in_place(self):
        """Test that deepcopied numpy arrays raise ValueError on in place modification."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        # Verify that attempting to modify copied arrays raises ValueError.
        with self.assertRaises(ValueError):
            copied.ampLer_Gs_Cgs[0] = 999.0

        with self.assertRaises(ValueError):
            copied.periodLer_Gs_Cgs[0] = 999.0

        with self.assertRaises(ValueError):
            copied.phaseLer_Gs_Cgs[0] = 999.0

        with self.assertRaises(ValueError):
            copied.ampAngles_Gs_to_Wn_ixyz[0] = 999.0

        with self.assertRaises(ValueError):
            copied.periodAngles_Gs_to_Wn_ixyz[0] = 999.0

        with self.assertRaises(ValueError):
            copied.phaseAngles_Gs_to_Wn_ixyz[0] = 999.0

        with self.assertRaises(ValueError):
            copied.rotationPointOffset_Gs_Ler[0] = 999.0

    def test_deepcopy_base_wing_is_independent(self):
        """Test that deepcopied base_wing is an independent object."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        # Verify base_wing is a different object.
        self.assertIsNot(copied.base_wing, original.base_wing)

        # Verify attributes are equal.
        self.assertEqual(copied.base_wing.name, original.base_wing.name)
        npt.assert_array_equal(
            copied.base_wing.Ler_Gs_Cgs, original.base_wing.Ler_Gs_Cgs
        )
        npt.assert_array_equal(
            copied.base_wing.angles_Gs_to_Wn_ixyz,
            original.base_wing.angles_Gs_to_Wn_ixyz,
        )

    def test_deepcopy_wing_cross_section_movements_are_independent(self):
        """Test that deepcopied wing_cross_section_movements are independent objects."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        # Verify wing_cross_section_movements tuple is different.
        self.assertIsNot(
            copied.wing_cross_section_movements, original.wing_cross_section_movements
        )

        # Verify each WingCrossSectionMovement is a different object.
        for i in range(len(original.wing_cross_section_movements)):
            self.assertIsNot(
                copied.wing_cross_section_movements[i],
                original.wing_cross_section_movements[i],
            )

    def test_deepcopy_resets_caches_to_none(self):
        """Test that deepcopy resets cached derived properties to None."""
        import copy

        original = self.wing_movement

        # Access cached properties to populate caches.
        _ = original.all_periods
        _ = original.max_period

        # Verify original caches are populated.
        self.assertIsNotNone(original._all_periods)
        self.assertIsNotNone(original._max_period)

        # Deepcopy the object.
        copied = copy.deepcopy(original)

        # Verify copied caches are reset to None.
        self.assertIsNone(copied._all_periods)
        self.assertIsNone(copied._max_period)

    def test_deepcopy_cached_properties_can_be_recomputed(self):
        """Test that cached properties work correctly after deepcopy."""
        import copy

        original = self.wing_movement

        # Get original cached values.
        original_all_periods = original.all_periods
        original_max_period = original.max_period

        # Deepcopy the object.
        copied = copy.deepcopy(original)

        # Verify cached properties can be computed and match original.
        self.assertEqual(copied.all_periods, original_all_periods)
        self.assertEqual(copied.max_period, original_max_period)

    def test_deepcopy_generate_wings_produces_same_results(self):
        """Test that generate_wings produces same results after deepcopy."""
        import copy

        original = self.wing_movement
        copied = copy.deepcopy(original)

        num_steps = 50
        delta_time = 0.01

        original_wings = original.generate_wings(
            num_steps=num_steps, delta_time=delta_time
        )
        copied_wings = copied.generate_wings(num_steps=num_steps, delta_time=delta_time)

        # Verify same number of Wings.
        self.assertEqual(len(copied_wings), len(original_wings))

        # Verify each Wing has matching attributes.
        for original_wing, copied_wing in zip(original_wings, copied_wings):
            npt.assert_array_equal(copied_wing.Ler_Gs_Cgs, original_wing.Ler_Gs_Cgs)
            npt.assert_array_equal(
                copied_wing.angles_Gs_to_Wn_ixyz, original_wing.angles_Gs_to_Wn_ixyz
            )
            self.assertEqual(copied_wing.name, original_wing.name)

    def test_deepcopy_handles_memo_correctly(self):
        """Test that deepcopy handles the memo dict correctly for circular references."""
        import copy

        original = self.wing_movement
        memo = {}

        # First deepcopy.
        copied1 = copy.deepcopy(original, memo)

        # Verify original is in memo.
        self.assertIn(id(original), memo)
        self.assertIs(memo[id(original)], copied1)

        # Second deepcopy with same memo should return same object.
        copied2 = copy.deepcopy(original, memo)
        self.assertIs(copied1, copied2)


if __name__ == "__main__":
    unittest.main()
