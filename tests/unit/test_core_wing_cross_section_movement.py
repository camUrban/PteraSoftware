"""This module contains classes to test CoreWingCrossSectionMovements."""

import copy
import unittest

import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps
from tests.unit.fixtures import core_wing_cross_section_movement_fixtures


class TestCoreWingCrossSectionMovement(unittest.TestCase):
    """This is a class with functions to test CoreWingCrossSectionMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all CoreWingCrossSectionMovement tests."""
        # Spacing test fixtures.
        cls.sine_spacing_Lp_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_sine_spacing_Lp_core_wing_cross_section_movement_fixture()
        )
        cls.uniform_spacing_Lp_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_uniform_spacing_Lp_core_wing_cross_section_movement_fixture()
        )
        cls.mixed_spacing_Lp_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_mixed_spacing_Lp_core_wing_cross_section_movement_fixture()
        )
        cls.sine_spacing_angles_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_sine_spacing_angles_core_wing_cross_section_movement_fixture()
        )
        cls.uniform_spacing_angles_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_uniform_spacing_angles_core_wing_cross_section_movement_fixture()
        )
        cls.mixed_spacing_angles_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_mixed_spacing_angles_core_wing_cross_section_movement_fixture()
        )

        # Additional test fixtures.
        cls.static_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture()
        )
        cls.basic_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_basic_core_wing_cross_section_movement_fixture()
        )
        cls.Lp_only_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_Lp_only_core_wing_cross_section_movement_fixture()
        )
        cls.angles_only_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_angles_only_core_wing_cross_section_movement_fixture()
        )
        cls.phase_offset_Lp_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_phase_offset_Lp_core_wing_cross_section_movement_fixture()
        )
        cls.phase_offset_angles_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_phase_offset_angles_core_wing_cross_section_movement_fixture()
        )
        cls.multiple_periods_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_multiple_periods_core_wing_cross_section_movement_fixture()
        )
        cls.custom_spacing_Lp_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_custom_spacing_Lp_core_wing_cross_section_movement_fixture()
        )
        cls.custom_spacing_angles_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_custom_spacing_angles_core_wing_cross_section_movement_fixture()
        )
        cls.mixed_custom_and_standard_spacing_core_wcs_movement = (
            core_wing_cross_section_movement_fixtures.make_mixed_custom_and_standard_spacing_core_wing_cross_section_movement_fixture()
        )

    def test_spacing_sine_for_Lp_Wcsp_Lpp(self):
        """Test that sine spacing actually produces sinusoidal motion for
        Lp_Wcsp_Lpp."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.sine_spacing_Lp_core_wcs_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
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
            self.uniform_spacing_Lp_core_wcs_movement.generate_wing_cross_sections(
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
        wing_cross_sections = (
            self.mixed_spacing_Lp_core_wcs_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract positions from generated WingCrossSections.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])
        y_positions = np.array([wcs.Lp_Wcsp_Lpp[1] for wcs in wing_cross_sections])
        z_positions = np.array([wcs.Lp_Wcsp_Lpp[2] for wcs in wing_cross_sections])

        # Calculate expected values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.5 + 1.0 * np.sin(2 * np.pi * times / 1.0)
        expected_y = 2.0 + 1.5 * signal.sawtooth(
            2 * np.pi * times / 1.0 + np.pi / 2, 0.5
        )
        expected_z = 0.2 + 0.5 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected values.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_spacing_sine_for_angles_Wcsp_to_Wcs_ixyz(self):
        """Test that sine spacing actually produces sinusoidal motion for
        angles_Wcsp_to_Wcs_ixyz."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.sine_spacing_angles_core_wcs_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract angles from generated WingCrossSections.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[0] for wcs in wing_cross_sections]
        )

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_angles = 10.0 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated angles match the expected sine wave.
        npt.assert_allclose(angles_z, expected_angles, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_angles_Wcsp_to_Wcs_ixyz(self):
        """Test that uniform spacing actually produces triangular wave motion for
        angles_Wcsp_to_Wcs_ixyz."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.uniform_spacing_angles_core_wcs_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract angles from generated WingCrossSections.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[0] for wcs in wing_cross_sections]
        )

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_angles = 10.0 * signal.sawtooth(
            2 * np.pi * times / 1.0 + np.pi / 2, 0.5
        )

        # Assert that the generated angles match the expected triangular wave.
        npt.assert_allclose(angles_z, expected_angles, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_angles_Wcsp_to_Wcs_ixyz(self):
        """Test that mixed spacing types work correctly for angles_Wcsp_to_Wcs_ixyz."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.mixed_spacing_angles_core_wcs_movement.generate_wing_cross_sections(
                num_steps=num_steps,
                delta_time=delta_time,
            )
        )

        # Extract angles from generated WingCrossSections.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[0] for wcs in wing_cross_sections]
        )
        angles_y = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[1] for wcs in wing_cross_sections]
        )
        angles_x = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[2] for wcs in wing_cross_sections]
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

    def test_initialization_valid_parameters(self):
        """Test CoreWingCrossSectionMovement initialization with valid parameters."""
        # Test that basic CoreWingCrossSectionMovement initializes correctly.
        core_wcs_movement = self.basic_core_wcs_movement
        self.assertIsInstance(
            core_wcs_movement,
            ps._core.CoreWingCrossSectionMovement,
        )
        self.assertIsInstance(
            core_wcs_movement.base_wing_cross_section,
            ps.geometry.wing_cross_section.WingCrossSection,
        )
        npt.assert_array_equal(
            core_wcs_movement.ampLp_Wcsp_Lpp, np.array([0.4, 0.3, 0.15])
        )
        npt.assert_array_equal(
            core_wcs_movement.periodLp_Wcsp_Lpp, np.array([2.0, 2.0, 2.0])
        )
        self.assertEqual(core_wcs_movement.spacingLp_Wcsp_Lpp, ("sine", "sine", "sine"))
        npt.assert_array_equal(
            core_wcs_movement.phaseLp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0])
        )
        npt.assert_array_equal(
            core_wcs_movement.ampAngles_Wcsp_to_Wcs_ixyz, np.array([15.0, 10.0, 5.0])
        )
        npt.assert_array_equal(
            core_wcs_movement.periodAngles_Wcsp_to_Wcs_ixyz, np.array([2.0, 2.0, 2.0])
        )
        self.assertEqual(
            core_wcs_movement.spacingAngles_Wcsp_to_Wcs_ixyz, ("sine", "sine", "sine")
        )
        npt.assert_array_equal(
            core_wcs_movement.phaseAngles_Wcsp_to_Wcs_ixyz, np.array([0.0, 0.0, 0.0])
        )

    def test_base_wing_cross_section_validation(self):
        """Test that base_wing_cross_section parameter validation works correctly."""
        from tests.unit.fixtures import geometry_fixtures

        # Test non-WingCrossSection raises error.
        with self.assertRaises(TypeError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section="not a wing cross section"
            )

        # Test None raises error.
        with self.assertRaises(TypeError):
            ps._core.CoreWingCrossSectionMovement(base_wing_cross_section=None)

        # Test valid WingCrossSection works.
        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs
        )
        self.assertEqual(core_wcs_movement.base_wing_cross_section, base_wcs)

    def test_ampLp_Wcsp_Lpp_validation(self):
        """Test ampLp_Wcsp_Lpp parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid values.
        valid_amps = [
            (0.0, 0.0, 0.0),
            (1.0, 2.0, 3.0),
            [0.5, 1.5, 2.5],
            np.array([0.1, 0.2, 0.3]),
        ]
        for amp in valid_amps:
            with self.subTest(amp=amp):
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs, ampLp_Wcsp_Lpp=amp
                )
                npt.assert_array_equal(core_wcs_movement.ampLp_Wcsp_Lpp, amp)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs, ampLp_Wcsp_Lpp=(-1.0, 0.0, 0.0)
            )

        # Test invalid types raise error.
        # noinspection PyTypeChecker
        with self.assertRaises((TypeError, ValueError)):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs, ampLp_Wcsp_Lpp="invalid"
            )

    def test_periodLp_Wcsp_Lpp_validation(self):
        """Test periodLp_Wcsp_Lpp parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid values.
        valid_periods = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0), [0.5, 1.5, 2.5]]
        for period in valid_periods:
            with self.subTest(period=period):
                # Need matching amps for non-zero periods.
                amp = tuple(1.0 if p > 0 else 0.0 for p in period)
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=amp,
                    periodLp_Wcsp_Lpp=period,
                )
                npt.assert_array_equal(core_wcs_movement.periodLp_Wcsp_Lpp, period)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                periodLp_Wcsp_Lpp=(-1.0, 1.0, 1.0),
            )

    def test_spacingLp_Wcsp_Lpp_validation(self):
        """Test spacingLp_Wcsp_Lpp parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid string values.
        valid_spacings = [
            ("sine", "sine", "sine"),
            ("uniform", "uniform", "uniform"),
            ("sine", "uniform", "sine"),
        ]
        for spacing in valid_spacings:
            with self.subTest(spacing=spacing):
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs, spacingLp_Wcsp_Lpp=spacing
                )
                self.assertEqual(core_wcs_movement.spacingLp_Wcsp_Lpp, spacing)

        # Test invalid string raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                spacingLp_Wcsp_Lpp=("invalid", "sine", "sine"),
            )

    def test_phaseLp_Wcsp_Lpp_validation(self):
        """Test phaseLp_Wcsp_Lpp parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid phase values within range (-180.0, 180.0].
        valid_phases = [
            (0.0, 0.0, 0.0),
            (90.0, 180.0, -90.0),
            (179.9, 0.0, -179.9),
        ]
        for phase in valid_phases:
            with self.subTest(phase=phase):
                # Need non-zero amps for non-zero phases.
                amp = tuple(1.0 if p != 0 else 0.0 for p in phase)
                period = tuple(1.0 if p != 0 else 0.0 for p in phase)
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=amp,
                    periodLp_Wcsp_Lpp=period,
                    phaseLp_Wcsp_Lpp=phase,
                )
                npt.assert_array_equal(core_wcs_movement.phaseLp_Wcsp_Lpp, phase)

        # Test phase > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                phaseLp_Wcsp_Lpp=(180.1, 0.0, 0.0),
            )

        # Test phase <= -180.0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                phaseLp_Wcsp_Lpp=(-180.0, 0.0, 0.0),
            )

    def test_ampAngles_Wcsp_to_Wcs_ixyz_validation(self):
        """Test ampAngles_Wcsp_to_Wcs_ixyz parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid amplitude values within range [0.0, 180.0].
        valid_amps = [
            (0.0, 0.0, 0.0),
            (45.0, 90.0, 135.0),
            (179.9, 0.0, 90.0),
            (180.0, 0.0, 0.0),
        ]
        for amp in valid_amps:
            with self.subTest(amp=amp):
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs, ampAngles_Wcsp_to_Wcs_ixyz=amp
                )
                npt.assert_array_equal(
                    core_wcs_movement.ampAngles_Wcsp_to_Wcs_ixyz, amp
                )

        # Test amplitude > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(180.1, 0.0, 0.0),
            )

        # Test negative amplitude raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(-1.0, 0.0, 0.0),
            )

    def test_periodAngles_Wcsp_to_Wcs_ixyz_validation(self):
        """Test periodAngles_Wcsp_to_Wcs_ixyz parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid periods.
        valid_periods = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0)]
        for period in valid_periods:
            with self.subTest(period=period):
                amp = tuple(10.0 if p > 0 else 0.0 for p in period)
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampAngles_Wcsp_to_Wcs_ixyz=amp,
                    periodAngles_Wcsp_to_Wcs_ixyz=period,
                )
                npt.assert_array_equal(
                    core_wcs_movement.periodAngles_Wcsp_to_Wcs_ixyz, period
                )

        # Test negative period raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 10.0, 10.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(-1.0, 1.0, 1.0),
            )

    def test_spacingAngles_Wcsp_to_Wcs_ixyz_validation(self):
        """Test spacingAngles_Wcsp_to_Wcs_ixyz parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid string values.
        valid_spacings = [
            ("sine", "sine", "sine"),
            ("uniform", "uniform", "uniform"),
            ("sine", "uniform", "sine"),
        ]
        for spacing in valid_spacings:
            with self.subTest(spacing=spacing):
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    spacingAngles_Wcsp_to_Wcs_ixyz=spacing,
                )
                self.assertEqual(
                    core_wcs_movement.spacingAngles_Wcsp_to_Wcs_ixyz, spacing
                )

        # Test invalid string raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                spacingAngles_Wcsp_to_Wcs_ixyz=("invalid", "sine", "sine"),
            )

    def test_phaseAngles_Wcsp_to_Wcs_ixyz_validation(self):
        """Test phaseAngles_Wcsp_to_Wcs_ixyz parameter validation."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test valid phase values within range (-180.0, 180.0].
        valid_phases = [(0.0, 0.0, 0.0), (90.0, 180.0, -90.0), (179.9, 0.0, -179.9)]
        for phase in valid_phases:
            with self.subTest(phase=phase):
                amp = tuple(10.0 if p != 0 else 0.0 for p in phase)
                period = tuple(1.0 if p != 0 else 0.0 for p in phase)
                core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampAngles_Wcsp_to_Wcs_ixyz=amp,
                    periodAngles_Wcsp_to_Wcs_ixyz=period,
                    phaseAngles_Wcsp_to_Wcs_ixyz=phase,
                )
                npt.assert_array_equal(
                    core_wcs_movement.phaseAngles_Wcsp_to_Wcs_ixyz, phase
                )

        # Test phase > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(10.0, 10.0, 10.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 1.0, 1.0),
                phaseAngles_Wcsp_to_Wcs_ixyz=(180.1, 0.0, 0.0),
            )

    def test_amp_period_relationship_Lp(self):
        """Test that if ampLp_Wcsp_Lpp element is 0, corresponding period must be 0."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amp=0 with period=0 works.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
        )
        self.assertIsNotNone(core_wcs_movement)

        # Test amp=0 with period!=0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 1.0, 0.0),
            )

    def test_amp_phase_relationship_Lp(self):
        """Test that if ampLp_Wcsp_Lpp element is 0, corresponding phase must be 0."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amp=0 with phase=0 works.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
            phaseLp_Wcsp_Lpp=(0.0, -90.0, 0.0),
        )
        self.assertIsNotNone(core_wcs_movement)

        # Test amp=0 with phase!=0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                phaseLp_Wcsp_Lpp=(45.0, -90.0, 0.0),
            )

    def test_amp_period_relationship_angles(self):
        """Test that if ampAngles_Wcsp_to_Wcs_ixyz element is 0, corresponding period must be 0."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amp=0 with period=0 works.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 1.0, 0.0),
        )
        self.assertIsNotNone(core_wcs_movement)

        # Test amp=0 with period!=0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 1.0, 0.0),
            )

    def test_amp_phase_relationship_angles(self):
        """Test that if ampAngles_Wcsp_to_Wcs_ixyz element is 0, corresponding phase must be 0."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amp=0 with phase=0 works.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 1.0, 0.0),
            phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, -90.0, 0.0),
        )
        self.assertIsNotNone(core_wcs_movement)

        # Test amp=0 with phase!=0 raises error.
        with self.assertRaises(ValueError):
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 1.0, 0.0),
                phaseAngles_Wcsp_to_Wcs_ixyz=(45.0, -90.0, 0.0),
            )

    def test_max_period_static_movement(self):
        """Test that max_period returns 0.0 for static movement."""
        core_wcs_movement = self.static_core_wcs_movement
        self.assertEqual(core_wcs_movement.max_period, 0.0)

    def test_max_period_Lp_only(self):
        """Test that max_period returns correct period for Lp-only movement."""
        core_wcs_movement = self.Lp_only_core_wcs_movement
        # periodLp_Wcsp_Lpp is (1.5, 1.5, 1.5), so max should be 1.5.
        self.assertEqual(core_wcs_movement.max_period, 1.5)

    def test_max_period_angles_only(self):
        """Test that max_period returns correct period for angles-only movement."""
        core_wcs_movement = self.angles_only_core_wcs_movement
        # periodAngles_Wcsp_to_Wcs_ixyz is (1.5, 1.5, 1.5), so max should be 1.5.
        self.assertEqual(core_wcs_movement.max_period, 1.5)

    def test_max_period_mixed(self):
        """Test that max_period returns maximum of all periods for mixed movement."""
        core_wcs_movement = self.multiple_periods_core_wcs_movement
        # periodLp_Wcsp_Lpp is (1.0, 2.0, 3.0).
        # periodAngles_Wcsp_to_Wcs_ixyz is (0.5, 1.5, 2.5).
        # Maximum should be 3.0.
        self.assertEqual(core_wcs_movement.max_period, 3.0)

    def test_max_period_multiple_dimensions(self):
        """Test max_period with multiple dimensions having different periods."""
        core_wcs_movement = self.basic_core_wcs_movement
        # Both periodLp_Wcsp_Lpp and periodAngles_Wcsp_to_Wcs_ixyz are (2.0, 2.0, 2.0).
        # Maximum should be 2.0.
        self.assertEqual(core_wcs_movement.max_period, 2.0)

    def test_all_periods_static_movement(self):
        """Test that all_periods returns empty tuple for static movement."""
        wing_cross_section_movement = self.static_core_wcs_movement
        self.assertEqual(wing_cross_section_movement.all_periods, ())

    def test_all_periods_Lp_only(self):
        """Test that all_periods returns correct periods for Lp only movement."""
        wing_cross_section_movement = self.Lp_only_core_wcs_movement
        # periodLp_Wcsp_Lpp is (1.5, 1.5, 1.5), all non zero.
        # periodAngles_Wcsp_to_Wcs_ixyz is (0.0, 0.0, 0.0).
        # Should return tuple with three 1.5 values.
        self.assertEqual(wing_cross_section_movement.all_periods, (1.5, 1.5, 1.5))

    def test_all_periods_angles_only(self):
        """Test that all_periods returns correct periods for angles only movement."""
        wing_cross_section_movement = self.angles_only_core_wcs_movement
        # periodLp_Wcsp_Lpp is (0.0, 0.0, 0.0).
        # periodAngles_Wcsp_to_Wcs_ixyz is (1.5, 1.5, 1.5), all non zero.
        # Should return tuple with three 1.5 values.
        self.assertEqual(wing_cross_section_movement.all_periods, (1.5, 1.5, 1.5))

    def test_all_periods_mixed(self):
        """Test that all_periods returns all non zero periods for mixed movement."""
        wing_cross_section_movement = self.multiple_periods_core_wcs_movement
        # periodLp_Wcsp_Lpp is (1.0, 2.0, 3.0).
        # periodAngles_Wcsp_to_Wcs_ixyz is (0.5, 1.5, 2.5).
        # Should return tuple with all six values.
        expected = (1.0, 2.0, 3.0, 0.5, 1.5, 2.5)
        self.assertEqual(wing_cross_section_movement.all_periods, expected)

    def test_all_periods_contains_duplicates(self):
        """Test that all_periods contains duplicate periods if they appear multiple
        times.
        """
        wing_cross_section_movement = self.basic_core_wcs_movement
        # Both periodLp_Wcsp_Lpp and periodAngles_Wcsp_to_Wcs_ixyz are (2.0, 2.0, 2.0).
        # Should return tuple with six 2.0 values (not deduplicated).
        expected = (2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        self.assertEqual(wing_cross_section_movement.all_periods, expected)

    def test_all_periods_partial_movement(self):
        """Test all_periods with only some dimensions having non zero periods."""
        wing_cross_section_movement = self.sine_spacing_Lp_core_wcs_movement
        # periodLp_Wcsp_Lpp is (1.0, 0.0, 0.0), only first element is non zero.
        # periodAngles_Wcsp_to_Wcs_ixyz is (0.0, 0.0, 0.0).
        # Should return tuple with one 1.0 value.
        self.assertEqual(wing_cross_section_movement.all_periods, (1.0,))

    def test_generate_wing_cross_sections_parameter_validation(self):
        """Test that generate_wing_cross_sections validates num_steps and delta_time."""
        core_wcs_movement = self.basic_core_wcs_movement

        # Test invalid num_steps.
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            core_wcs_movement.generate_wing_cross_sections(num_steps=0, delta_time=0.01)

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=-1, delta_time=0.01
            )

        with self.assertRaises(TypeError):
            core_wcs_movement.generate_wing_cross_sections(
                num_steps="invalid", delta_time=0.01
            )

        # Test invalid delta_time.
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            core_wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.0)

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=-0.01
            )

        with self.assertRaises(TypeError):
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time="invalid"
            )

    def test_generate_wing_cross_sections_returns_correct_length(self):
        """Test that generate_wing_cross_sections returns list of correct length."""
        core_wcs_movement = self.basic_core_wcs_movement

        test_num_steps = [1, 5, 10, 50, 100]
        for num_steps in test_num_steps:
            with self.subTest(num_steps=num_steps):
                wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(wing_cross_sections), num_steps)

    def test_generate_wing_cross_sections_returns_correct_types(self):
        """Test that generate_wing_cross_sections returns WingCrossSections."""
        core_wcs_movement = self.basic_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=10, delta_time=0.01
        )

        # Verify all elements are WingCrossSections.
        for wcs in wing_cross_sections:
            self.assertIsInstance(wcs, ps.geometry.wing_cross_section.WingCrossSection)

    def test_generate_wing_cross_sections_preserves_non_changing_attributes(self):
        """Test that generate_wing_cross_sections preserves non-changing attributes."""
        core_wcs_movement = self.basic_core_wcs_movement
        base_wcs = core_wcs_movement.base_wing_cross_section

        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=10, delta_time=0.01
        )

        # Check that non-changing attributes are preserved.
        for wcs in wing_cross_sections:
            self.assertEqual(wcs.airfoil, base_wcs.airfoil)
            self.assertEqual(wcs.chord, base_wcs.chord)
            self.assertEqual(wcs.num_spanwise_panels, base_wcs.num_spanwise_panels)
            self.assertEqual(
                wcs.control_surface_symmetry_type,
                base_wcs.control_surface_symmetry_type,
            )
            self.assertEqual(
                wcs.control_surface_hinge_point, base_wcs.control_surface_hinge_point
            )
            self.assertEqual(
                wcs.control_surface_deflection, base_wcs.control_surface_deflection
            )
            self.assertEqual(wcs.spanwise_spacing, base_wcs.spanwise_spacing)

    def test_generate_wing_cross_sections_static_movement(self):
        """Test that static movement produces constant positions and angles."""
        core_wcs_movement = self.static_core_wcs_movement
        base_wcs = core_wcs_movement.base_wing_cross_section

        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=50, delta_time=0.01
        )

        # All WingCrossSections should have same Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_ixyz.
        for wcs in wing_cross_sections:
            npt.assert_array_equal(wcs.Lp_Wcsp_Lpp, base_wcs.Lp_Wcsp_Lpp)
            npt.assert_array_equal(
                wcs.angles_Wcsp_to_Wcs_ixyz, base_wcs.angles_Wcsp_to_Wcs_ixyz
            )

    def test_phase_offset_Lp(self):
        """Test that phase shifts initial position correctly for Lp_Wcsp_Lpp."""
        core_wcs_movement = self.phase_offset_Lp_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=100, delta_time=0.01
        )

        # Extract positions.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])
        y_positions = np.array([wcs.Lp_Wcsp_Lpp[1] for wcs in wing_cross_sections])
        z_positions = np.array([wcs.Lp_Wcsp_Lpp[2] for wcs in wing_cross_sections])

        # Verify that phase offset causes non-zero initial values.
        # With phase offsets, the first values should not all be at the base position.
        self.assertFalse(np.allclose(x_positions[0], 0.0, atol=1e-10))
        self.assertFalse(np.allclose(y_positions[0], 0.0, atol=1e-10))
        self.assertFalse(np.allclose(z_positions[0], 0.0, atol=1e-10))

    def test_phase_offset_angles(self):
        """Test that phase shifts initial angles correctly for angles_Wcsp_to_Wcs_ixyz."""
        core_wcs_movement = self.phase_offset_angles_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=100, delta_time=0.01
        )

        # Extract angles.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[0] for wcs in wing_cross_sections]
        )
        angles_y = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[1] for wcs in wing_cross_sections]
        )
        angles_x = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[2] for wcs in wing_cross_sections]
        )

        # Verify that phase offset causes non-zero initial values.
        # With phase offsets, the first values should not all be at the base angles.
        self.assertFalse(np.allclose(angles_z[0], 0.0, atol=1e-10))
        self.assertFalse(np.allclose(angles_y[0], 0.0, atol=1e-10))
        self.assertFalse(np.allclose(angles_x[0], 0.0, atol=1e-10))

    def test_single_dimension_movement_Lp(self):
        """Test that only one dimension of Lp_Wcsp_Lpp moves."""
        core_wcs_movement = self.sine_spacing_Lp_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=50, delta_time=0.01
        )

        # Extract positions.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])
        y_positions = np.array([wcs.Lp_Wcsp_Lpp[1] for wcs in wing_cross_sections])
        z_positions = np.array([wcs.Lp_Wcsp_Lpp[2] for wcs in wing_cross_sections])

        # Only x should vary, y and z should be constant.
        self.assertFalse(np.allclose(x_positions, x_positions[0]))
        npt.assert_array_equal(y_positions, y_positions[0])
        npt.assert_array_equal(z_positions, z_positions[0])

    def test_single_dimension_movement_angles(self):
        """Test that only one dimension of angles_Wcsp_to_Wcs_ixyz moves."""
        core_wcs_movement = self.sine_spacing_angles_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=50, delta_time=0.01
        )

        # Extract angles.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[0] for wcs in wing_cross_sections]
        )
        angles_y = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[1] for wcs in wing_cross_sections]
        )
        angles_x = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[2] for wcs in wing_cross_sections]
        )

        # Only z should vary, y and x should be constant.
        self.assertFalse(np.allclose(angles_z, angles_z[0]))
        npt.assert_array_equal(angles_y, angles_y[0])
        npt.assert_array_equal(angles_x, angles_x[0])

    def test_boundary_amplitude_angles(self):
        """Test amplitude at boundary value (180.0 degrees)."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amplitude at 180.0 works.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampAngles_Wcsp_to_Wcs_ixyz=(180.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
        )
        self.assertEqual(core_wcs_movement.ampAngles_Wcsp_to_Wcs_ixyz[0], 180.0)

    def test_boundary_phase_values(self):
        """Test phase at boundary values (-179.9, 0.0, and 180.0)."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test phase = 0.0 works.
        core_wcs_movement1 = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )
        self.assertEqual(core_wcs_movement1.phaseLp_Wcsp_Lpp[0], 0.0)

        # Test phase = 180.0 works (upper boundary, inclusive).
        core_wcs_movement2 = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            phaseLp_Wcsp_Lpp=(180.0, 0.0, 0.0),
        )
        self.assertEqual(core_wcs_movement2.phaseLp_Wcsp_Lpp[0], 180.0)

        # Test phase = -179.9 works (near lower boundary).
        core_wcs_movement3 = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
            phaseLp_Wcsp_Lpp=(-179.9, 0.0, 0.0),
        )
        self.assertEqual(core_wcs_movement3.phaseLp_Wcsp_Lpp[0], -179.9)

    def test_custom_spacing_function_Lp(self):
        """Test that custom spacing function works for Lp_Wcsp_Lpp."""
        core_wcs_movement = self.custom_spacing_Lp_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=100, delta_time=0.01
        )

        # Extract x-positions.
        x_positions = np.array([wcs.Lp_Wcsp_Lpp[0] for wcs in wing_cross_sections])

        # Verify that values vary (not constant).
        self.assertFalse(np.allclose(x_positions, x_positions[0]))

        # Verify that values are within expected range.
        # For custom_harmonic with amp=1.0, values should be in [-1.0, 1.0].
        self.assertTrue(np.all(x_positions >= -1.1))
        self.assertTrue(np.all(x_positions <= 1.1))

    def test_custom_spacing_function_angles(self):
        """Test that custom spacing function works for angles_Wcsp_to_Wcs_ixyz."""
        core_wcs_movement = self.custom_spacing_angles_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=100, delta_time=0.01
        )

        # Extract z-angles.
        angles_z = np.array(
            [wcs.angles_Wcsp_to_Wcs_ixyz[0] for wcs in wing_cross_sections]
        )

        # Verify that values vary (not constant).
        self.assertFalse(np.allclose(angles_z, angles_z[0]))

        # Verify that values are within expected range.
        # For custom_triangle with amp=10.0, values should be in [-10.0, 10.0].
        self.assertTrue(np.all(angles_z >= -11.0))
        self.assertTrue(np.all(angles_z <= 11.0))

    def test_custom_spacing_function_mixed_with_standard(self):
        """Test that custom and standard spacing functions can be mixed."""
        core_wcs_movement = self.mixed_custom_and_standard_spacing_core_wcs_movement
        wing_cross_sections = core_wcs_movement.generate_wing_cross_sections(
            num_steps=100, delta_time=0.01
        )

        # Verify that WingCrossSections are generated successfully.
        self.assertEqual(len(wing_cross_sections), 100)
        for wcs in wing_cross_sections:
            self.assertIsInstance(wcs, ps.geometry.wing_cross_section.WingCrossSection)

    def test_custom_function_validation_invalid_start_value(self):
        """Test that custom function with invalid start value raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that doesn't start at 0.
        def invalid_nonzero_start(x):
            return np.sin(x) + 1.0

        # Should raise error during initialization or generation.
        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_nonzero_start, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_custom_function_validation_invalid_end_value(self):
        """Test that custom function with invalid end value raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that doesn't return to 0 at 2*pi.
        def invalid_nonzero_end(x):
            return np.sin(x) + 0.1

        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_nonzero_end, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_custom_function_validation_invalid_mean(self):
        """Test that custom function with invalid mean raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function with non-zero mean.
        def invalid_nonzero_mean(x):
            return np.sin(x) + 0.5

        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_nonzero_mean, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_custom_function_validation_invalid_amplitude(self):
        """Test that custom function with invalid amplitude raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function with wrong amplitude.
        def invalid_wrong_amplitude(x):
            return 2.0 * np.sin(x)

        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_wrong_amplitude, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_custom_function_validation_not_periodic(self):
        """Test that custom function that is not periodic raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that is not periodic.
        def invalid_not_periodic(x):
            return np.tanh(x)

        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_not_periodic, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_custom_function_validation_returns_non_finite(self):
        """Test that custom function returning NaN or Inf raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that returns NaN.
        def invalid_non_finite(x):
            return np.where(x < np.pi, np.sin(x), np.nan)

        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_non_finite, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_custom_function_validation_wrong_shape(self):
        """Test that custom function returning wrong shape raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that returns wrong shape.
        def invalid_wrong_shape(x):
            return np.sin(x)[: len(x) // 2]

        with self.assertRaises(ValueError):
            core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=(invalid_wrong_shape, "sine", "sine"),
            )
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time=0.01
            )

    def test_unsafe_amplitude_causes_error_Lp(self):
        """Test that amplitude too high for base Lp value causes error during generation."""
        from tests.unit.fixtures import geometry_fixtures

        # Use root fixture with Lp_Wcsp_Lpp = [0.0, 0.0, 0.0].
        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Create CoreWingCrossSectionMovement with amplitude that will drive the second element in
        # Lp_Wcsp_Lpp negative, which is never allowed by WingCrossSection.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
            periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
            spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
            phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        )

        # Generating WingCrossSections should raise ValueError when Lp goes negative.
        with self.assertRaises(ValueError) as context:
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=100, delta_time=0.01
            )

        # Verify the error message is about Lp_Wcsp_Lpp validation.
        self.assertIn("Lp_Wcsp_Lpp", str(context.exception))

    def test_unsafe_amplitude_causes_error_angles(self):
        """Test that amplitude too high for base angle value causes error during generation."""
        from tests.unit.fixtures import geometry_fixtures

        # Use root fixture with angles = [0.0, 0.0, 0.0].
        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Create CoreWingCrossSectionMovement with amplitude that will drive angles out of valid range.
        # Valid range for angles is (-180, 180], so amplitude 181 with base 0 will exceed.
        core_wcs_movement = ps._core.CoreWingCrossSectionMovement(
            base_wing_cross_section=base_wcs,
            ampAngles_Wcsp_to_Wcs_ixyz=(179.0, 0.0, 0.0),
            periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
            phaseAngles_Wcsp_to_Wcs_ixyz=(90.0, 0.0, 0.0),
        )

        # Generating WingCrossSections should raise ValueError when angles exceed range.
        with self.assertRaises(ValueError) as context:
            core_wcs_movement.generate_wing_cross_sections(
                num_steps=100, delta_time=0.01
            )

        # Verify the error message is about angles_Wcsp_to_Wcs_ixyz validation.
        self.assertIn("angles_Wcsp_to_Wcs_ixyz", str(context.exception))


class TestCoreWingCrossSectionMovementImmutability(unittest.TestCase):
    """Tests for CoreWingCrossSectionMovement attribute immutability."""

    def setUp(self):
        """Set up test fixtures for immutability tests."""
        self.core_wing_cross_section_movement = (
            core_wing_cross_section_movement_fixtures.make_basic_core_wing_cross_section_movement_fixture()
        )

    def test_immutable_base_wing_cross_section_property(self):
        """Test that base_wing_cross_section property is read only."""
        from tests.unit.fixtures import geometry_fixtures

        new_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.base_wing_cross_section = (
                new_wing_cross_section
            )

    def test_immutable_ampLp_Wcsp_Lpp_property(self):
        """Test that ampLp_Wcsp_Lpp property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.ampLp_Wcsp_Lpp = np.array(
                [1.0, 2.0, 3.0]
            )

    def test_immutable_ampLp_Wcsp_Lpp_array_read_only(self):
        """Test that ampLp_Wcsp_Lpp array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.core_wing_cross_section_movement.ampLp_Wcsp_Lpp[0] = 999.0

    def test_immutable_periodLp_Wcsp_Lpp_property(self):
        """Test that periodLp_Wcsp_Lpp property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.periodLp_Wcsp_Lpp = np.array(
                [1.0, 2.0, 3.0]
            )

    def test_immutable_periodLp_Wcsp_Lpp_array_read_only(self):
        """Test that periodLp_Wcsp_Lpp array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.core_wing_cross_section_movement.periodLp_Wcsp_Lpp[0] = 999.0

    def test_immutable_spacingLp_Wcsp_Lpp_property(self):
        """Test that spacingLp_Wcsp_Lpp property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.spacingLp_Wcsp_Lpp = (
                "uniform",
                "uniform",
                "uniform",
            )

    def test_immutable_phaseLp_Wcsp_Lpp_property(self):
        """Test that phaseLp_Wcsp_Lpp property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.phaseLp_Wcsp_Lpp = np.array(
                [45.0, 45.0, 45.0]
            )

    def test_immutable_phaseLp_Wcsp_Lpp_array_read_only(self):
        """Test that phaseLp_Wcsp_Lpp array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.core_wing_cross_section_movement.phaseLp_Wcsp_Lpp[0] = 999.0

    def test_immutable_ampAngles_Wcsp_to_Wcs_ixyz_property(self):
        """Test that ampAngles_Wcsp_to_Wcs_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz = np.array(
                [1.0, 2.0, 3.0]
            )

    def test_immutable_ampAngles_Wcsp_to_Wcs_ixyz_array_read_only(self):
        """Test that ampAngles_Wcsp_to_Wcs_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.core_wing_cross_section_movement.ampAngles_Wcsp_to_Wcs_ixyz[0] = 999.0

    def test_immutable_periodAngles_Wcsp_to_Wcs_ixyz_property(self):
        """Test that periodAngles_Wcsp_to_Wcs_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.periodAngles_Wcsp_to_Wcs_ixyz = (
                np.array([1.0, 2.0, 3.0])
            )

    def test_immutable_periodAngles_Wcsp_to_Wcs_ixyz_array_read_only(self):
        """Test that periodAngles_Wcsp_to_Wcs_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.core_wing_cross_section_movement.periodAngles_Wcsp_to_Wcs_ixyz[0] = (
                999.0
            )

    def test_immutable_spacingAngles_Wcsp_to_Wcs_ixyz_property(self):
        """Test that spacingAngles_Wcsp_to_Wcs_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.spacingAngles_Wcsp_to_Wcs_ixyz = (
                "uniform",
                "uniform",
                "uniform",
            )

    def test_immutable_phaseAngles_Wcsp_to_Wcs_ixyz_property(self):
        """Test that phaseAngles_Wcsp_to_Wcs_ixyz property is read only."""
        with self.assertRaises(AttributeError):
            self.core_wing_cross_section_movement.phaseAngles_Wcsp_to_Wcs_ixyz = (
                np.array([45.0, 45.0, 45.0])
            )

    def test_immutable_phaseAngles_Wcsp_to_Wcs_ixyz_array_read_only(self):
        """Test that phaseAngles_Wcsp_to_Wcs_ixyz array cannot be modified in place."""
        with self.assertRaises(ValueError):
            self.core_wing_cross_section_movement.phaseAngles_Wcsp_to_Wcs_ixyz[0] = (
                999.0
            )


class TestCoreWingCrossSectionMovementCaching(unittest.TestCase):
    """Tests for CoreWingCrossSectionMovement caching behavior."""

    def setUp(self):
        """Set up test fixtures for caching tests."""
        self.core_wing_cross_section_movement = (
            core_wing_cross_section_movement_fixtures.make_basic_core_wing_cross_section_movement_fixture()
        )

    def test_all_periods_caching_returns_same_object(self):
        """Test that repeated access to all_periods returns the same cached object."""
        all_periods_1 = self.core_wing_cross_section_movement.all_periods
        all_periods_2 = self.core_wing_cross_section_movement.all_periods
        self.assertIs(all_periods_1, all_periods_2)

    def test_max_period_caching_returns_same_value(self):
        """Test that repeated access to max_period returns the same cached value."""
        max_period_1 = self.core_wing_cross_section_movement.max_period
        max_period_2 = self.core_wing_cross_section_movement.max_period
        # Since floats are immutable, we check equality rather than identity.
        self.assertEqual(max_period_1, max_period_2)


class TestCoreWingCrossSectionMovementDeepcopy(unittest.TestCase):
    """Tests for CoreWingCrossSectionMovement deepcopy behavior."""

    def setUp(self):
        """Set up test fixtures for deepcopy tests."""
        self.core_wing_cross_section_movement = (
            core_wing_cross_section_movement_fixtures.make_basic_core_wing_cross_section_movement_fixture()
        )

    def test_deepcopy_returns_new_instance(self):
        """Test that deepcopy returns a new CoreWingCrossSectionMovement instance."""
        original = self.core_wing_cross_section_movement
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps._core.CoreWingCrossSectionMovement)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_attribute_values(self):
        """Test that deepcopy preserves all attribute values."""
        original = self.core_wing_cross_section_movement
        copied = copy.deepcopy(original)

        # Check numpy array attributes.
        npt.assert_array_equal(copied.ampLp_Wcsp_Lpp, original.ampLp_Wcsp_Lpp)
        npt.assert_array_equal(copied.periodLp_Wcsp_Lpp, original.periodLp_Wcsp_Lpp)
        npt.assert_array_equal(copied.phaseLp_Wcsp_Lpp, original.phaseLp_Wcsp_Lpp)
        npt.assert_array_equal(
            copied.ampAngles_Wcsp_to_Wcs_ixyz, original.ampAngles_Wcsp_to_Wcs_ixyz
        )
        npt.assert_array_equal(
            copied.periodAngles_Wcsp_to_Wcs_ixyz, original.periodAngles_Wcsp_to_Wcs_ixyz
        )
        npt.assert_array_equal(
            copied.phaseAngles_Wcsp_to_Wcs_ixyz, original.phaseAngles_Wcsp_to_Wcs_ixyz
        )

        # Check tuple attributes.
        self.assertEqual(copied.spacingLp_Wcsp_Lpp, original.spacingLp_Wcsp_Lpp)
        self.assertEqual(
            copied.spacingAngles_Wcsp_to_Wcs_ixyz,
            original.spacingAngles_Wcsp_to_Wcs_ixyz,
        )

    def test_deepcopy_numpy_arrays_are_independent(self):
        """Test that deepcopied numpy arrays are independent objects."""
        original = self.core_wing_cross_section_movement
        copied = copy.deepcopy(original)

        # Verify arrays are different objects.
        self.assertIsNot(copied.ampLp_Wcsp_Lpp, original.ampLp_Wcsp_Lpp)
        self.assertIsNot(copied.periodLp_Wcsp_Lpp, original.periodLp_Wcsp_Lpp)
        self.assertIsNot(copied.phaseLp_Wcsp_Lpp, original.phaseLp_Wcsp_Lpp)
        self.assertIsNot(
            copied.ampAngles_Wcsp_to_Wcs_ixyz, original.ampAngles_Wcsp_to_Wcs_ixyz
        )
        self.assertIsNot(
            copied.periodAngles_Wcsp_to_Wcs_ixyz, original.periodAngles_Wcsp_to_Wcs_ixyz
        )
        self.assertIsNot(
            copied.phaseAngles_Wcsp_to_Wcs_ixyz, original.phaseAngles_Wcsp_to_Wcs_ixyz
        )

    def test_deepcopy_numpy_arrays_cannot_be_modified_in_place(self):
        """Test that deepcopied numpy arrays raise ValueError on in place modification."""
        original = self.core_wing_cross_section_movement
        copied = copy.deepcopy(original)

        # Verify that attempting to modify copied arrays raises ValueError.
        with self.assertRaises(ValueError):
            copied.ampLp_Wcsp_Lpp[0] = 999.0

        with self.assertRaises(ValueError):
            copied.periodLp_Wcsp_Lpp[0] = 999.0

        with self.assertRaises(ValueError):
            copied.phaseLp_Wcsp_Lpp[0] = 999.0

        with self.assertRaises(ValueError):
            copied.ampAngles_Wcsp_to_Wcs_ixyz[0] = 999.0

        with self.assertRaises(ValueError):
            copied.periodAngles_Wcsp_to_Wcs_ixyz[0] = 999.0

        with self.assertRaises(ValueError):
            copied.phaseAngles_Wcsp_to_Wcs_ixyz[0] = 999.0

    def test_deepcopy_base_wing_cross_section_is_independent(self):
        """Test that deepcopied base_wing_cross_section is an independent object."""
        original = self.core_wing_cross_section_movement
        copied = copy.deepcopy(original)

        # Verify base_wing_cross_section is a different object.
        self.assertIsNot(
            copied.base_wing_cross_section, original.base_wing_cross_section
        )

        # Verify attributes are equal.
        self.assertEqual(
            copied.base_wing_cross_section.chord,
            original.base_wing_cross_section.chord,
        )
        npt.assert_array_equal(
            copied.base_wing_cross_section.Lp_Wcsp_Lpp,
            original.base_wing_cross_section.Lp_Wcsp_Lpp,
        )
        npt.assert_array_equal(
            copied.base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            original.base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
        )

    def test_deepcopy_resets_caches_to_none(self):
        """Test that deepcopy resets cached derived properties to None."""
        original = self.core_wing_cross_section_movement

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
        original = self.core_wing_cross_section_movement

        # Get original cached values.
        original_all_periods = original.all_periods
        original_max_period = original.max_period

        # Deepcopy the object.
        copied = copy.deepcopy(original)

        # Verify cached properties can be computed and match original.
        self.assertEqual(copied.all_periods, original_all_periods)
        self.assertEqual(copied.max_period, original_max_period)

    def test_deepcopy_generate_wing_cross_sections_produces_same_results(self):
        """Test that generate_wing_cross_sections produces same results after deepcopy."""
        original = self.core_wing_cross_section_movement
        copied = copy.deepcopy(original)

        num_steps = 50
        delta_time = 0.01

        original_wcs_list = original.generate_wing_cross_sections(
            num_steps=num_steps, delta_time=delta_time
        )
        copied_wcs_list = copied.generate_wing_cross_sections(
            num_steps=num_steps, delta_time=delta_time
        )

        # Verify same number of WingCrossSections.
        self.assertEqual(len(copied_wcs_list), len(original_wcs_list))

        # Verify each WingCrossSection has matching attributes.
        for original_wcs, copied_wcs in zip(original_wcs_list, copied_wcs_list):
            npt.assert_array_equal(copied_wcs.Lp_Wcsp_Lpp, original_wcs.Lp_Wcsp_Lpp)
            npt.assert_array_equal(
                copied_wcs.angles_Wcsp_to_Wcs_ixyz, original_wcs.angles_Wcsp_to_Wcs_ixyz
            )
            self.assertEqual(copied_wcs.chord, original_wcs.chord)

    def test_deepcopy_handles_memo_correctly(self):
        """Test that deepcopy handles the memo dict correctly for circular references."""
        original = self.core_wing_cross_section_movement
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
