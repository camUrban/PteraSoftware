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
