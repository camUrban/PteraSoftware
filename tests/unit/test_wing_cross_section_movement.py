"""This module contains a class to test WingCrossSectionMovements."""

import unittest
import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps

from tests.unit.fixtures import wing_cross_section_movement_fixtures


class TestWingCrossSectionMovement(unittest.TestCase):
    """This is a class with functions to test WingCrossSectionMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all WingCrossSectionMovement tests."""
        # Spacing test fixtures.
        cls.sine_spacing_Lp_wcs_movement = (
            wing_cross_section_movement_fixtures.make_sine_spacing_Lp_wing_cross_section_movement_fixture()
        )
        cls.uniform_spacing_Lp_wcs_movement = (
            wing_cross_section_movement_fixtures.make_uniform_spacing_Lp_wing_cross_section_movement_fixture()
        )
        cls.mixed_spacing_Lp_wcs_movement = (
            wing_cross_section_movement_fixtures.make_mixed_spacing_Lp_wing_cross_section_movement_fixture()
        )
        cls.sine_spacing_angles_wcs_movement = (
            wing_cross_section_movement_fixtures.make_sine_spacing_angles_wing_cross_section_movement_fixture()
        )
        cls.uniform_spacing_angles_wcs_movement = (
            wing_cross_section_movement_fixtures.make_uniform_spacing_angles_wing_cross_section_movement_fixture()
        )
        cls.mixed_spacing_angles_wcs_movement = (
            wing_cross_section_movement_fixtures.make_mixed_spacing_angles_wing_cross_section_movement_fixture()
        )

        # Additional test fixtures.
        cls.static_wcs_movement = (
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture()
        )
        cls.basic_wcs_movement = (
            wing_cross_section_movement_fixtures.make_basic_wing_cross_section_movement_fixture()
        )
        cls.Lp_only_wcs_movement = (
            wing_cross_section_movement_fixtures.make_Lp_only_wing_cross_section_movement_fixture()
        )
        cls.angles_only_wcs_movement = (
            wing_cross_section_movement_fixtures.make_angles_only_wing_cross_section_movement_fixture()
        )
        cls.phase_offset_Lp_wcs_movement = (
            wing_cross_section_movement_fixtures.make_phase_offset_Lp_wing_cross_section_movement_fixture()
        )
        cls.phase_offset_angles_wcs_movement = (
            wing_cross_section_movement_fixtures.make_phase_offset_angles_wing_cross_section_movement_fixture()
        )
        cls.multiple_periods_wcs_movement = (
            wing_cross_section_movement_fixtures.make_multiple_periods_wing_cross_section_movement_fixture()
        )
        cls.custom_spacing_Lp_wcs_movement = (
            wing_cross_section_movement_fixtures.make_custom_spacing_Lp_wing_cross_section_movement_fixture()
        )
        cls.custom_spacing_angles_wcs_movement = (
            wing_cross_section_movement_fixtures.make_custom_spacing_angles_wing_cross_section_movement_fixture()
        )
        cls.mixed_custom_and_standard_spacing_wcs_movement = (
            wing_cross_section_movement_fixtures.make_mixed_custom_and_standard_spacing_wing_cross_section_movement_fixture()
        )

    def test_spacing_sine_for_Lp_Wcsp_Lpp(self):
        """Test that sine spacing actually produces sinusoidal motion for
        Lp_Wcsp_Lpp."""
        num_steps = 100
        delta_time = 0.01
        wing_cross_sections = (
            self.sine_spacing_Lp_wcs_movement.generate_wing_cross_sections(
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
            self.uniform_spacing_Lp_wcs_movement.generate_wing_cross_sections(
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
            self.mixed_spacing_Lp_wcs_movement.generate_wing_cross_sections(
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
            self.sine_spacing_angles_wcs_movement.generate_wing_cross_sections(
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
            self.uniform_spacing_angles_wcs_movement.generate_wing_cross_sections(
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
            self.mixed_spacing_angles_wcs_movement.generate_wing_cross_sections(
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
        """Test WingCrossSectionMovement initialization with valid parameters."""
        # Test that basic WingCrossSectionMovement initializes correctly.
        wcs_movement = self.basic_wcs_movement
        self.assertIsInstance(
            wcs_movement,
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement,
        )
        self.assertIsInstance(
            wcs_movement.base_wing_cross_section,
            ps.geometry.wing_cross_section.WingCrossSection,
        )
        npt.assert_array_equal(wcs_movement.ampLp_Wcsp_Lpp, np.array([0.4, 0.3, 0.15]))
        npt.assert_array_equal(
            wcs_movement.periodLp_Wcsp_Lpp, np.array([2.0, 2.0, 2.0])
        )
        self.assertEqual(wcs_movement.spacingLp_Wcsp_Lpp, ("sine", "sine", "sine"))
        npt.assert_array_equal(wcs_movement.phaseLp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0]))
        npt.assert_array_equal(
            wcs_movement.ampAngles_Wcsp_to_Wcs_ixyz, np.array([15.0, 10.0, 5.0])
        )
        npt.assert_array_equal(
            wcs_movement.periodAngles_Wcsp_to_Wcs_ixyz, np.array([2.0, 2.0, 2.0])
        )
        self.assertEqual(
            wcs_movement.spacingAngles_Wcsp_to_Wcs_ixyz, ("sine", "sine", "sine")
        )
        npt.assert_array_equal(
            wcs_movement.phaseAngles_Wcsp_to_Wcs_ixyz, np.array([0.0, 0.0, 0.0])
        )

    def test_base_wing_cross_section_validation(self):
        """Test that base_wing_cross_section parameter validation works correctly."""
        from tests.unit.fixtures import geometry_fixtures

        # Test non-WingCrossSection raises error.
        with self.assertRaises(TypeError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section="not a wing cross section"
            )

        # Test None raises error.
        with self.assertRaises(TypeError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=None
            )

        # Test valid WingCrossSection works.
        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs
            )
        )
        self.assertEqual(wcs_movement.base_wing_cross_section, base_wcs)

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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs, ampLp_Wcsp_Lpp=amp
                    )
                )
                npt.assert_array_equal(wcs_movement.ampLp_Wcsp_Lpp, amp)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs, ampLp_Wcsp_Lpp=(-1.0, 0.0, 0.0)
            )

        # Test invalid types raise error.
        # noinspection PyTypeChecker
        with self.assertRaises((TypeError, ValueError)):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs,
                        ampLp_Wcsp_Lpp=amp,
                        periodLp_Wcsp_Lpp=period,
                    )
                )
                npt.assert_array_equal(wcs_movement.periodLp_Wcsp_Lpp, period)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs, spacingLp_Wcsp_Lpp=spacing
                    )
                )
                self.assertEqual(wcs_movement.spacingLp_Wcsp_Lpp, spacing)

        # Test invalid string raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs,
                        ampLp_Wcsp_Lpp=amp,
                        periodLp_Wcsp_Lpp=period,
                        phaseLp_Wcsp_Lpp=phase,
                    )
                )
                npt.assert_array_equal(wcs_movement.phaseLp_Wcsp_Lpp, phase)

        # Test phase > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                periodLp_Wcsp_Lpp=(1.0, 1.0, 1.0),
                phaseLp_Wcsp_Lpp=(180.1, 0.0, 0.0),
            )

        # Test phase <= -180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs, ampAngles_Wcsp_to_Wcs_ixyz=amp
                    )
                )
                npt.assert_array_equal(wcs_movement.ampAngles_Wcsp_to_Wcs_ixyz, amp)

        # Test amplitude > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(180.1, 0.0, 0.0),
            )

        # Test negative amplitude raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs,
                        ampAngles_Wcsp_to_Wcs_ixyz=amp,
                        periodAngles_Wcsp_to_Wcs_ixyz=period,
                    )
                )
                npt.assert_array_equal(
                    wcs_movement.periodAngles_Wcsp_to_Wcs_ixyz, period
                )

        # Test negative period raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs,
                        spacingAngles_Wcsp_to_Wcs_ixyz=spacing,
                    )
                )
                self.assertEqual(wcs_movement.spacingAngles_Wcsp_to_Wcs_ixyz, spacing)

        # Test invalid string raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
                wcs_movement = (
                    ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                        base_wing_cross_section=base_wcs,
                        ampAngles_Wcsp_to_Wcs_ixyz=amp,
                        periodAngles_Wcsp_to_Wcs_ixyz=period,
                        phaseAngles_Wcsp_to_Wcs_ixyz=phase,
                    )
                )
                npt.assert_array_equal(wcs_movement.phaseAngles_Wcsp_to_Wcs_ixyz, phase)

        # Test phase > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
            )
        )
        self.assertIsNotNone(wcs_movement)

        # Test amp=0 with period!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 1.0, 0.0),
            )

    def test_amp_phase_relationship_Lp(self):
        """Test that if ampLp_Wcsp_Lpp element is 0, corresponding phase must be 0."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amp=0 with phase=0 works.
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                phaseLp_Wcsp_Lpp=(0.0, -90.0, 0.0),
            )
        )
        self.assertIsNotNone(wcs_movement)

        # Test amp=0 with phase!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
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
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 1.0, 0.0),
            )
        )
        self.assertIsNotNone(wcs_movement)

        # Test amp=0 with period!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 1.0, 0.0),
            )

    def test_amp_phase_relationship_angles(self):
        """Test that if ampAngles_Wcsp_to_Wcs_ixyz element is 0, corresponding phase must be 0."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test amp=0 with phase=0 works.
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 1.0, 0.0),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, -90.0, 0.0),
            )
        )
        self.assertIsNotNone(wcs_movement)

        # Test amp=0 with phase!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 10.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 1.0, 0.0),
                phaseAngles_Wcsp_to_Wcs_ixyz=(45.0, -90.0, 0.0),
            )

    def test_max_period_static_movement(self):
        """Test that max_period returns 0.0 for static movement."""
        wcs_movement = self.static_wcs_movement
        self.assertEqual(wcs_movement.max_period, 0.0)

    def test_max_period_Lp_only(self):
        """Test that max_period returns correct period for Lp-only movement."""
        wcs_movement = self.Lp_only_wcs_movement
        # periodLp_Wcsp_Lpp is (1.5, 1.5, 1.5), so max should be 1.5.
        self.assertEqual(wcs_movement.max_period, 1.5)

    def test_max_period_angles_only(self):
        """Test that max_period returns correct period for angles-only movement."""
        wcs_movement = self.angles_only_wcs_movement
        # periodAngles_Wcsp_to_Wcs_ixyz is (1.5, 1.5, 1.5), so max should be 1.5.
        self.assertEqual(wcs_movement.max_period, 1.5)

    def test_max_period_mixed(self):
        """Test that max_period returns maximum of all periods for mixed movement."""
        wcs_movement = self.multiple_periods_wcs_movement
        # periodLp_Wcsp_Lpp is (1.0, 2.0, 3.0).
        # periodAngles_Wcsp_to_Wcs_ixyz is (0.5, 1.5, 2.5).
        # Maximum should be 3.0.
        self.assertEqual(wcs_movement.max_period, 3.0)

    def test_max_period_multiple_dimensions(self):
        """Test max_period with multiple dimensions having different periods."""
        wcs_movement = self.basic_wcs_movement
        # Both periodLp_Wcsp_Lpp and periodAngles_Wcsp_to_Wcs_ixyz are (2.0, 2.0, 2.0).
        # Maximum should be 2.0.
        self.assertEqual(wcs_movement.max_period, 2.0)

    def test_generate_wing_cross_sections_parameter_validation(self):
        """Test that generate_wing_cross_sections validates num_steps and delta_time."""
        wcs_movement = self.basic_wcs_movement

        # Test invalid num_steps.
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            wcs_movement.generate_wing_cross_sections(num_steps=0, delta_time=0.01)

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            wcs_movement.generate_wing_cross_sections(num_steps=-1, delta_time=0.01)

        with self.assertRaises(TypeError):
            wcs_movement.generate_wing_cross_sections(
                num_steps="invalid", delta_time=0.01
            )

        # Test invalid delta_time.
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.0)

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=-0.01)

        with self.assertRaises(TypeError):
            wcs_movement.generate_wing_cross_sections(
                num_steps=10, delta_time="invalid"
            )

    def test_generate_wing_cross_sections_returns_correct_length(self):
        """Test that generate_wing_cross_sections returns list of correct length."""
        wcs_movement = self.basic_wcs_movement

        test_num_steps = [1, 5, 10, 50, 100]
        for num_steps in test_num_steps:
            with self.subTest(num_steps=num_steps):
                wing_cross_sections = wcs_movement.generate_wing_cross_sections(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(wing_cross_sections), num_steps)

    def test_generate_wing_cross_sections_returns_correct_types(self):
        """Test that generate_wing_cross_sections returns WingCrossSections."""
        wcs_movement = self.basic_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
            num_steps=10, delta_time=0.01
        )

        # Verify all elements are WingCrossSections.
        for wcs in wing_cross_sections:
            self.assertIsInstance(wcs, ps.geometry.wing_cross_section.WingCrossSection)

    def test_generate_wing_cross_sections_preserves_non_changing_attributes(self):
        """Test that generate_wing_cross_sections preserves non-changing attributes."""
        wcs_movement = self.basic_wcs_movement
        base_wcs = wcs_movement.base_wing_cross_section

        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = self.static_wcs_movement
        base_wcs = wcs_movement.base_wing_cross_section

        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
            num_steps=50, delta_time=0.01
        )

        # All WingCrossSections should have same Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_ixyz.
        for wcs in wing_cross_sections:
            npt.assert_array_equal(wcs.Lp_Wcsp_Lpp, base_wcs.Lp_Wcsp_Lpp)
            npt.assert_array_equal(
                wcs.angles_Wcsp_to_Wcs_ixyz, base_wcs.angles_Wcsp_to_Wcs_ixyz
            )

    def test_generate_wing_cross_sections_different_num_steps(self):
        """Test generate_wing_cross_sections with various num_steps values."""
        wcs_movement = self.basic_wcs_movement

        num_steps_list = [1, 10, 25, 100, 200]
        for num_steps in num_steps_list:
            with self.subTest(num_steps=num_steps):
                wing_cross_sections = wcs_movement.generate_wing_cross_sections(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(wing_cross_sections), num_steps)

    def test_generate_wing_cross_sections_different_delta_time(self):
        """Test generate_wing_cross_sections with various delta_time values."""
        wcs_movement = self.basic_wcs_movement

        delta_time_list = [0.001, 0.01, 0.1, 1.0]
        num_steps = 50
        for delta_time in delta_time_list:
            with self.subTest(delta_time=delta_time):
                wing_cross_sections = wcs_movement.generate_wing_cross_sections(
                    num_steps=num_steps, delta_time=delta_time
                )
                self.assertEqual(len(wing_cross_sections), num_steps)

    def test_phase_offset_Lp(self):
        """Test that phase shifts initial position correctly for Lp_Wcsp_Lpp."""
        wcs_movement = self.phase_offset_Lp_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = self.phase_offset_angles_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = self.sine_spacing_Lp_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = self.sine_spacing_angles_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(180.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
            )
        )
        self.assertEqual(wcs_movement.ampAngles_Wcsp_to_Wcs_ixyz[0], 180.0)

    def test_boundary_phase_values(self):
        """Test phase at boundary values (-179.9, 0.0, and 180.0)."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Test phase = 0.0 works.
        wcs_movement1 = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            )
        )
        self.assertEqual(wcs_movement1.phaseLp_Wcsp_Lpp[0], 0.0)

        # Test phase = 180.0 works (upper boundary, inclusive).
        wcs_movement2 = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                phaseLp_Wcsp_Lpp=(180.0, 0.0, 0.0),
            )
        )
        self.assertEqual(wcs_movement2.phaseLp_Wcsp_Lpp[0], 180.0)

        # Test phase = -179.9 works (near lower boundary).
        wcs_movement3 = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                phaseLp_Wcsp_Lpp=(-179.9, 0.0, 0.0),
            )
        )
        self.assertEqual(wcs_movement3.phaseLp_Wcsp_Lpp[0], -179.9)

    def test_custom_spacing_function_Lp(self):
        """Test that custom spacing function works for Lp_Wcsp_Lpp."""
        wcs_movement = self.custom_spacing_Lp_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = self.custom_spacing_angles_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
        wcs_movement = self.mixed_custom_and_standard_spacing_wcs_movement
        wing_cross_sections = wcs_movement.generate_wing_cross_sections(
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
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_nonzero_start, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_invalid_end_value(self):
        """Test that custom function with invalid end value raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that doesn't return to 0 at 2*pi.
        def invalid_nonzero_end(x):
            return np.sin(x) + 0.1

        with self.assertRaises(ValueError):
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_nonzero_end, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_invalid_mean(self):
        """Test that custom function with invalid mean raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function with non-zero mean.
        def invalid_nonzero_mean(x):
            return np.sin(x) + 0.5

        with self.assertRaises(ValueError):
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_nonzero_mean, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_invalid_amplitude(self):
        """Test that custom function with invalid amplitude raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function with wrong amplitude.
        def invalid_wrong_amplitude(x):
            return 2.0 * np.sin(x)

        with self.assertRaises(ValueError):
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_wrong_amplitude, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_not_periodic(self):
        """Test that custom function that is not periodic raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that is not periodic.
        def invalid_not_periodic(x):
            return np.tanh(x)

        with self.assertRaises(ValueError):
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_not_periodic, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_returns_non_finite(self):
        """Test that custom function returning NaN or Inf raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that returns NaN.
        def invalid_non_finite(x):
            return np.where(x < np.pi, np.sin(x), np.nan)

        with self.assertRaises(ValueError):
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_non_finite, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_custom_function_validation_wrong_shape(self):
        """Test that custom function returning wrong shape raises error."""
        from tests.unit.fixtures import geometry_fixtures

        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Define invalid custom function that returns wrong shape.
        def invalid_wrong_shape(x):
            return np.sin(x)[: len(x) // 2]

        with self.assertRaises(ValueError):
            wcs_movement = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=base_wcs,
                    ampLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    periodLp_Wcsp_Lpp=(1.0, 0.0, 0.0),
                    spacingLp_Wcsp_Lpp=(invalid_wrong_shape, "sine", "sine"),
                )
            )
            wcs_movement.generate_wing_cross_sections(num_steps=10, delta_time=0.01)

    def test_unsafe_amplitude_causes_error_Lp(self):
        """Test that amplitude too high for base Lp value causes error during generation."""
        from tests.unit.fixtures import geometry_fixtures

        # Use root fixture with Lp_Wcsp_Lpp = [0.0, 0.0, 0.0].
        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Create WingCrossSectionMovement with amplitude that will drive the second element in
        # Lp_Wcsp_Lpp negative, which is never allowed by WingCrossSection.
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            )
        )

        # Generating WingCrossSections should raise ValueError when Lp goes negative.
        with self.assertRaises(ValueError) as context:
            wcs_movement.generate_wing_cross_sections(num_steps=100, delta_time=0.01)

        # Verify the error message is about Lp_Wcsp_Lpp validation.
        self.assertIn("Lp_Wcsp_Lpp", str(context.exception))

    def test_unsafe_amplitude_causes_error_angles(self):
        """Test that amplitude too high for base angle value causes error during generation."""
        from tests.unit.fixtures import geometry_fixtures

        # Use root fixture with angles = [0.0, 0.0, 0.0].
        base_wcs = geometry_fixtures.make_root_wing_cross_section_fixture()

        # Create WingCrossSectionMovement with amplitude that will drive angles out of valid range.
        # Valid range for angles is (-180, 180], so amplitude 181 with base 0 will exceed.
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wcs,
                ampAngles_Wcsp_to_Wcs_ixyz=(179.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(1.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(90.0, 0.0, 0.0),
            )
        )

        # Generating WingCrossSections should raise ValueError when angles exceed range.
        with self.assertRaises(ValueError) as context:
            wcs_movement.generate_wing_cross_sections(num_steps=100, delta_time=0.01)

        # Verify the error message is about angles_Wcsp_to_Wcs_ixyz validation.
        self.assertIn("angles_Wcsp_to_Wcs_ixyz", str(context.exception))


if __name__ == "__main__":
    unittest.main()
