"""This module contains a class to test AirplaneMovements."""

import unittest
import numpy as np
import numpy.testing as npt
from scipy import signal

import pterasoftware as ps

from tests.unit.fixtures import airplane_movement_fixtures
from tests.unit.fixtures import geometry_fixtures
from tests.unit.fixtures import wing_movement_fixtures


class TestAirplaneMovement(unittest.TestCase):
    """This is a class with functions to test AirplaneMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all AirplaneMovement tests."""
        # Spacing test fixtures.
        cls.sine_spacing_Cg_airplane_movement = (
            airplane_movement_fixtures.make_sine_spacing_Cg_airplane_movement_fixture()
        )
        cls.uniform_spacing_Cg_airplane_movement = (
            airplane_movement_fixtures.make_uniform_spacing_Cg_airplane_movement_fixture()
        )
        cls.mixed_spacing_Cg_airplane_movement = (
            airplane_movement_fixtures.make_mixed_spacing_Cg_airplane_movement_fixture()
        )

        # Additional test fixtures.
        cls.static_airplane_movement = (
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        )
        cls.basic_airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        cls.Cg_airplane_movement = (
            airplane_movement_fixtures.make_Cg_airplane_movement_fixture()
        )
        cls.phase_offset_Cg_airplane_movement = (
            airplane_movement_fixtures.make_phase_offset_Cg_airplane_movement_fixture()
        )
        cls.multiple_periods_airplane_movement = (
            airplane_movement_fixtures.make_multiple_periods_airplane_movement_fixture()
        )
        cls.custom_spacing_Cg_airplane_movement = (
            airplane_movement_fixtures.make_custom_spacing_Cg_airplane_movement_fixture()
        )
        cls.mixed_custom_and_standard_spacing_airplane_movement = (
            airplane_movement_fixtures.make_mixed_custom_and_standard_spacing_airplane_movement_fixture()
        )

    def test_spacing_sine_for_Cg_GP1_CgP1(self):
        """Test that sine spacing actually produces sinusoidal motion for Cg_GP1_CgP1."""
        num_steps = 10
        delta_time = 0.01
        airplanes = self.sine_spacing_Cg_airplane_movement.generate_airplanes(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions (in Earth axes, relative to the simulation starting
        # point) from generated Airplanes.
        x_positions = np.array([airplane.Cg_GP1_CgP1[0] for airplane in airplanes])

        # Calculate expected sine wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.1 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected sine wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_uniform_for_Cg_GP1_CgP1(self):
        """Test that uniform spacing actually produces triangular wave motion for
        Cg_GP1_CgP1."""
        num_steps = 10
        delta_time = 0.01
        airplanes = self.uniform_spacing_Cg_airplane_movement.generate_airplanes(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions (in Earth axes, relative to the simulation starting
        # point) from generated Airplanes.
        x_positions = np.array([airplane.Cg_GP1_CgP1[0] for airplane in airplanes])

        # Calculate expected triangular wave values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.1 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)

        # Assert that the generated positions match the expected triangular wave.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)

    def test_spacing_mixed_for_Cg_GP1_CgP1(self):
        """Test that mixed spacing types work correctly for Cg_GP1_CgP1."""
        num_steps = 10
        delta_time = 0.01
        airplanes = self.mixed_spacing_Cg_airplane_movement.generate_airplanes(
            num_steps=num_steps,
            delta_time=delta_time,
        )

        # Extract positions (in Earth axes, relative to the simulation starting
        # point) from generated Airplanes.
        x_positions = np.array([airplane.Cg_GP1_CgP1[0] for airplane in airplanes])
        y_positions = np.array([airplane.Cg_GP1_CgP1[1] for airplane in airplanes])
        z_positions = np.array([airplane.Cg_GP1_CgP1[2] for airplane in airplanes])

        # Calculate expected values.
        times = np.linspace(0, num_steps * delta_time, num_steps, endpoint=False)
        expected_x = 0.1 * np.sin(2 * np.pi * times / 1.0)
        expected_y = 0.08 * signal.sawtooth(2 * np.pi * times / 1.0 + np.pi / 2, 0.5)
        expected_z = 0.06 * np.sin(2 * np.pi * times / 1.0)

        # Assert that the generated positions match the expected values.
        npt.assert_allclose(x_positions, expected_x, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(y_positions, expected_y, rtol=1e-10, atol=1e-14)
        npt.assert_allclose(z_positions, expected_z, rtol=1e-10, atol=1e-14)

    def test_initialization_valid_parameters(self):
        """Test AirplaneMovement initialization with valid parameters."""
        # Test that basic AirplaneMovement initializes correctly.
        airplane_movement = self.basic_airplane_movement
        self.assertIsInstance(
            airplane_movement, ps.movements.airplane_movement.AirplaneMovement
        )
        self.assertIsInstance(
            airplane_movement.base_airplane, ps.geometry.airplane.Airplane
        )
        self.assertIsInstance(airplane_movement.wing_movements, list)
        self.assertEqual(len(airplane_movement.wing_movements), 1)
        self.assertIsInstance(
            airplane_movement.wing_movements[0],
            ps.movements.wing_movement.WingMovement,
        )
        npt.assert_array_equal(
            airplane_movement.ampCg_GP1_CgP1, np.array([0.0, 0.0, 0.0])
        )
        npt.assert_array_equal(
            airplane_movement.periodCg_GP1_CgP1, np.array([0.0, 0.0, 0.0])
        )
        self.assertEqual(airplane_movement.spacingCg_GP1_CgP1, ("sine", "sine", "sine"))
        npt.assert_array_equal(
            airplane_movement.phaseCg_GP1_CgP1, np.array([0.0, 0.0, 0.0])
        )

    def test_base_airplane_validation(self):
        """Test that base_airplane parameter validation works correctly."""
        # Test non-Airplane raises error.
        with self.assertRaises(TypeError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane="not an airplane",
                wing_movements=[
                    wing_movement_fixtures.make_static_wing_movement_fixture()
                ],
            )

        # Test None raises error.
        with self.assertRaises(TypeError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=None,
                wing_movements=[
                    wing_movement_fixtures.make_static_wing_movement_fixture()
                ],
            )

        # Test valid Airplane works.
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane, wing_movements=wing_movements
        )
        self.assertEqual(airplane_movement.base_airplane, base_airplane)

    def test_wing_movements_validation(self):
        """Test that wing_movements parameter validation works correctly."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()

        # Test non-list raises error.
        with self.assertRaises(TypeError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane, wing_movements="not a list"
            )

        # Test empty list raises error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane, wing_movements=[]
            )

        # Test wrong number of WingMovements raises error.
        base_airplane_multi_wing = geometry_fixtures.make_multi_wing_airplane_fixture()
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane_multi_wing,
                wing_movements=[
                    wing_movement_fixtures.make_static_wing_movement_fixture()
                ],
            )

        # Test non-WingMovement element raises error.
        with self.assertRaises(TypeError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane, wing_movements=["not a wing movement"]
            )

        # Test valid list works.
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane, wing_movements=wing_movements
        )
        self.assertEqual(airplane_movement.wing_movements, wing_movements)

    def test_ampCg_GP1_CgP1_validation(self):
        """Test ampCg_GP1_CgP1 parameter validation."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

        # Test valid values.
        valid_amps = [
            (0.0, 0.0, 0.0),
            (1.0, 2.0, 3.0),
            [0.5, 1.5, 2.5],
            np.array([0.1, 0.2, 0.3]),
        ]
        for amp in valid_amps:
            with self.subTest(amp=amp):
                airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
                    base_airplane=base_airplane,
                    wing_movements=wing_movements,
                    ampCg_GP1_CgP1=amp,
                )
                npt.assert_array_equal(airplane_movement.ampCg_GP1_CgP1, amp)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1=(-1.0, 0.0, 0.0),
            )

        # Test invalid types raise error.
        # noinspection PyTypeChecker
        with self.assertRaises((TypeError, ValueError)):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1="invalid",
            )

    def test_periodCg_GP1_CgP1_validation(self):
        """Test periodCg_GP1_CgP1 parameter validation."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

        # Test valid values.
        valid_periods = [(0.0, 0.0, 0.0), (1.0, 2.0, 3.0), [0.5, 1.5, 2.5]]
        for period in valid_periods:
            with self.subTest(period=period):
                # Need matching amps for non-zero periods.
                amp = tuple(1.0 if p > 0 else 0.0 for p in period)
                airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
                    base_airplane=base_airplane,
                    wing_movements=wing_movements,
                    ampCg_GP1_CgP1=amp,
                    periodCg_GP1_CgP1=period,
                )
                npt.assert_array_equal(airplane_movement.periodCg_GP1_CgP1, period)

        # Test negative values raise error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1=(1.0, 1.0, 1.0),
                periodCg_GP1_CgP1=(-1.0, 1.0, 1.0),
            )

    def test_spacingCg_GP1_CgP1_validation(self):
        """Test spacingCg_GP1_CgP1 parameter validation."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

        # Test valid string values.
        valid_spacings = [
            ("sine", "sine", "sine"),
            ("uniform", "uniform", "uniform"),
            ("sine", "uniform", "sine"),
        ]
        for spacing in valid_spacings:
            with self.subTest(spacing=spacing):
                airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
                    base_airplane=base_airplane,
                    wing_movements=wing_movements,
                    spacingCg_GP1_CgP1=spacing,
                )
                self.assertEqual(airplane_movement.spacingCg_GP1_CgP1, spacing)

        # Test invalid string raises error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                spacingCg_GP1_CgP1=("invalid", "sine", "sine"),
            )

    def test_phaseCg_GP1_CgP1_validation(self):
        """Test phaseCg_GP1_CgP1 parameter validation."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

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
                airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
                    base_airplane=base_airplane,
                    wing_movements=wing_movements,
                    ampCg_GP1_CgP1=amp,
                    periodCg_GP1_CgP1=period,
                    phaseCg_GP1_CgP1=phase,
                )
                npt.assert_array_equal(airplane_movement.phaseCg_GP1_CgP1, phase)

        # Test phase > 180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1=(1.0, 1.0, 1.0),
                periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
                phaseCg_GP1_CgP1=(180.1, 0.0, 0.0),
            )

        # Test phase <= -180.0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1=(1.0, 1.0, 1.0),
                periodCg_GP1_CgP1=(1.0, 1.0, 1.0),
                phaseCg_GP1_CgP1=(-180.0, 0.0, 0.0),
            )

    def test_amp_period_relationship_Cg(self):
        """Test that if ampCg_GP1_CgP1 element is 0, corresponding period must be 0."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

        # Test amp=0 with period=0 works.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.0, 1.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 1.0, 0.0),
        )
        self.assertIsNotNone(airplane_movement)

        # Test amp=0 with period!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1=(0.0, 1.0, 0.0),
                periodCg_GP1_CgP1=(1.0, 1.0, 0.0),
            )

    def test_amp_phase_relationship_Cg(self):
        """Test that if ampCg_GP1_CgP1 element is 0, corresponding phase must be 0."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

        # Test amp=0 with phase=0 works.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(0.0, 1.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 1.0, 0.0),
            phaseCg_GP1_CgP1=(0.0, -90.0, 0.0),
        )
        self.assertIsNotNone(airplane_movement)

        # Test amp=0 with phase!=0 raises error.
        with self.assertRaises(ValueError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
                ampCg_GP1_CgP1=(0.0, 1.0, 0.0),
                periodCg_GP1_CgP1=(0.0, 1.0, 0.0),
                phaseCg_GP1_CgP1=(45.0, -90.0, 0.0),
            )

    def test_max_period_static_movement(self):
        """Test that max_period returns 0.0 for static movement."""
        airplane_movement = self.static_airplane_movement
        self.assertEqual(airplane_movement.max_period, 0.0)

    def test_max_period_Cg(self):
        """Test that max_period returns correct period."""
        airplane_movement = self.Cg_airplane_movement
        # periodCg_GP1_CgP1 is (1.5, 1.5, 1.5), so max should be 1.5.
        self.assertEqual(airplane_movement.max_period, 1.5)

    def test_generate_airplanes_parameter_validation(self):
        """Test that generate_airplanes validates num_steps and delta_time."""
        airplane_movement = self.basic_airplane_movement

        # Test invalid num_steps.
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            airplane_movement.generate_airplanes(num_steps=0, delta_time=0.01)

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            airplane_movement.generate_airplanes(num_steps=-1, delta_time=0.01)

        with self.assertRaises(TypeError):
            airplane_movement.generate_airplanes(num_steps="invalid", delta_time=0.01)

        # Test invalid delta_time.
        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            airplane_movement.generate_airplanes(num_steps=10, delta_time=0.0)

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            airplane_movement.generate_airplanes(num_steps=10, delta_time=-0.01)

        with self.assertRaises(TypeError):
            airplane_movement.generate_airplanes(num_steps=10, delta_time="invalid")

    def test_generate_airplanes_returns_correct_length(self):
        """Test that generate_airplanes returns list of correct length."""
        airplane_movement = self.basic_airplane_movement

        test_num_steps = [1, 5, 10, 50, 100]
        for num_steps in test_num_steps:
            with self.subTest(num_steps=num_steps):
                airplanes = airplane_movement.generate_airplanes(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(airplanes), num_steps)

    def test_generate_airplanes_returns_correct_types(self):
        """Test that generate_airplanes returns Airplanes."""
        airplane_movement = self.basic_airplane_movement
        airplanes = airplane_movement.generate_airplanes(num_steps=10, delta_time=0.01)

        # Verify all elements are Airplanes.
        for airplane in airplanes:
            self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_generate_airplanes_preserves_non_changing_attributes(self):
        """Test that generate_airplanes preserves non-changing attributes."""
        airplane_movement = self.basic_airplane_movement
        base_airplane = airplane_movement.base_airplane

        airplanes = airplane_movement.generate_airplanes(num_steps=10, delta_time=0.01)

        # Check that non-changing attributes are preserved. Note: s_ref, c_ref,
        # and b_ref are NOT included because they are calculated from the Wings,
        # which change due to WingMovement or WingCrossSectionMovement.
        for airplane in airplanes:
            self.assertEqual(airplane.name, base_airplane.name)
            self.assertEqual(airplane.weight, base_airplane.weight)

    def test_generate_airplanes_static_movement(self):
        """Test that static movement produces constant positions and angles."""
        airplane_movement = self.static_airplane_movement
        base_airplane = airplane_movement.base_airplane

        airplanes = airplane_movement.generate_airplanes(num_steps=50, delta_time=0.01)

        # All Airplanes should have same Cg_GP1_CgP1.
        for airplane in airplanes:
            npt.assert_array_equal(airplane.Cg_GP1_CgP1, base_airplane.Cg_GP1_CgP1)

    def test_generate_airplanes_different_num_steps(self):
        """Test generate_airplanes with various num_steps values."""
        airplane_movement = self.basic_airplane_movement

        num_steps_list = [1, 10, 25, 100, 200]
        for num_steps in num_steps_list:
            with self.subTest(num_steps=num_steps):
                airplanes = airplane_movement.generate_airplanes(
                    num_steps=num_steps, delta_time=0.01
                )
                self.assertEqual(len(airplanes), num_steps)

    def test_generate_airplanes_different_delta_time(self):
        """Test generate_airplanes with various delta_time values."""
        airplane_movement = self.basic_airplane_movement

        delta_time_list = [0.001, 0.01, 0.1, 1.0]
        num_steps = 10
        for delta_time in delta_time_list:
            with self.subTest(delta_time=delta_time):
                airplanes = airplane_movement.generate_airplanes(
                    num_steps=num_steps, delta_time=delta_time
                )
                self.assertEqual(len(airplanes), num_steps)

    def test_phase_offset_Cg(self):
        """Test that phase shifts initial position correctly for Cg_GP1_CgP1."""
        airplane_movement = self.phase_offset_Cg_airplane_movement
        airplanes = airplane_movement.generate_airplanes(num_steps=100, delta_time=0.01)

        # Extract positions (in Earth axes, relative to the simulation starting
        # point).
        x_positions = np.array([airplane.Cg_GP1_CgP1[0] for airplane in airplanes])
        y_positions = np.array([airplane.Cg_GP1_CgP1[1] for airplane in airplanes])
        z_positions = np.array([airplane.Cg_GP1_CgP1[2] for airplane in airplanes])

        # Verify that phase offset causes non-zero initial values.
        # With phase offsets, the first values should not all be at the base position.
        self.assertFalse(np.allclose(x_positions[0], 0.0, atol=1e-10))
        self.assertFalse(np.allclose(y_positions[0], 0.0, atol=1e-10))
        self.assertFalse(np.allclose(z_positions[0], 0.0, atol=1e-10))

    def test_single_dimension_movement_Cg(self):
        """Test that only one dimension of Cg_GP1_CgP1 moves."""
        airplane_movement = self.sine_spacing_Cg_airplane_movement
        airplanes = airplane_movement.generate_airplanes(num_steps=50, delta_time=0.01)

        # Extract positions (in Earth axes, relative to the simulation starting
        # point).
        x_positions = np.array([airplane.Cg_GP1_CgP1[0] for airplane in airplanes])
        y_positions = np.array([airplane.Cg_GP1_CgP1[1] for airplane in airplanes])
        z_positions = np.array([airplane.Cg_GP1_CgP1[2] for airplane in airplanes])

        # Only x should vary, y and z should be constant.
        self.assertFalse(np.allclose(x_positions, x_positions[0]))
        npt.assert_array_equal(y_positions, y_positions[0])
        npt.assert_array_equal(z_positions, z_positions[0])

    def test_boundary_phase_values(self):
        """Test phase at boundary values (-179.9, 0.0, and 180.0)."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]

        # Test phase = 0.0 works.
        airplane_movement1 = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(1.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
        self.assertEqual(airplane_movement1.phaseCg_GP1_CgP1[0], 0.0)

        # Test phase = 180.0 works (upper boundary, inclusive).
        airplane_movement2 = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(1.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
            phaseCg_GP1_CgP1=(180.0, 0.0, 0.0),
        )
        self.assertEqual(airplane_movement2.phaseCg_GP1_CgP1[0], 180.0)

        # Test phase = -179.9 works (near lower boundary).
        airplane_movement3 = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=(1.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(1.0, 0.0, 0.0),
            phaseCg_GP1_CgP1=(-179.9, 0.0, 0.0),
        )
        self.assertEqual(airplane_movement3.phaseCg_GP1_CgP1[0], -179.9)

    def test_custom_spacing_function_Cg(self):
        """Test that custom spacing function works for Cg_GP1_CgP1."""
        airplane_movement = self.custom_spacing_Cg_airplane_movement
        airplanes = airplane_movement.generate_airplanes(num_steps=100, delta_time=0.01)

        # Extract x-positions (in Earth axes, relative to the simulation starting
        # point).
        x_positions = np.array([airplane.Cg_GP1_CgP1[0] for airplane in airplanes])

        # Verify that values vary (not constant).
        self.assertFalse(np.allclose(x_positions, x_positions[0]))

        # Verify that values are within expected range.
        # For custom_harmonic with amp=0.08, values should be in [-0.08, 0.08].
        self.assertTrue(np.all(x_positions >= -0.09))
        self.assertTrue(np.all(x_positions <= 0.09))

    def test_custom_spacing_function_mixed_with_standard(self):
        """Test that custom and standard spacing functions can be mixed."""
        airplane_movement = self.mixed_custom_and_standard_spacing_airplane_movement
        airplanes = airplane_movement.generate_airplanes(num_steps=100, delta_time=0.01)

        # Verify that Airplanes are generated successfully.
        self.assertEqual(len(airplanes), 100)
        for airplane in airplanes:
            self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)


if __name__ == "__main__":
    unittest.main()
