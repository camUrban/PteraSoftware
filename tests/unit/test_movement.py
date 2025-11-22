"""This module contains a class to test Movements."""

import unittest

import pterasoftware as ps

from tests.unit.fixtures import movement_fixtures
from tests.unit.fixtures import airplane_movement_fixtures
from tests.unit.fixtures import operating_point_fixtures
from tests.unit.fixtures import geometry_fixtures


class TestMovement(unittest.TestCase):
    """This is a class with functions to test Movements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all Movement tests."""
        cls.static_movement = movement_fixtures.make_static_movement_fixture()
        cls.basic_movement = movement_fixtures.make_basic_movement_fixture()
        cls.static_movement_with_explicit_num_steps = (
            movement_fixtures.make_static_movement_with_explicit_num_steps_fixture()
        )
        cls.non_static_movement_with_explicit_num_steps = (
            movement_fixtures.make_non_static_movement_with_explicit_num_steps_fixture()
        )
        cls.movement_with_custom_delta_time = (
            movement_fixtures.make_movement_with_custom_delta_time_fixture()
        )
        cls.movement_with_multiple_airplanes = (
            movement_fixtures.make_movement_with_multiple_airplanes_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test Movement initialization with valid parameters."""
        movement = self.basic_movement
        self.assertIsInstance(movement, ps.movements.movement.Movement)
        self.assertIsInstance(movement.airplane_movements, list)
        self.assertEqual(len(movement.airplane_movements), 1)
        self.assertIsInstance(
            movement.airplane_movements[0],
            ps.movements.airplane_movement.AirplaneMovement,
        )
        self.assertIsInstance(
            movement.operating_point_movement,
            ps.movements.operating_point_movement.OperatingPointMovement,
        )
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)
        self.assertEqual(movement.num_cycles, 1)
        self.assertIsNone(movement.num_chords)
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

    def test_airplane_movements_validation_not_list(self):
        """Test that airplane_movements must be a list."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.Movement(
                airplane_movements="not a list",
                operating_point_movement=operating_point_movement,
                num_chords=10,
            )

    def test_airplane_movements_validation_empty_list(self):
        """Test that airplane_movements must have at least one element."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=[],
                operating_point_movement=operating_point_movement,
                num_chords=10,
            )

    def test_airplane_movements_validation_invalid_element_type(self):
        """Test that all elements in airplane_movements must be AirplaneMovements."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.Movement(
                airplane_movements=["not an airplane movement"],
                operating_point_movement=operating_point_movement,
                num_chords=10,
            )

    def test_operating_point_movement_validation(self):
        """Test that operating_point_movement must be an OperatingPointMovement."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]

        with self.assertRaises(TypeError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement="not an operating point movement",
                num_chords=10,
            )

    def test_delta_time_validation_positive(self):
        """Test that delta_time must be positive if provided."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with positive delta_time works.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.01,
            num_chords=10,
        )
        self.assertEqual(movement.delta_time, 0.01)

    def test_delta_time_validation_zero(self):
        """Test that delta_time cannot be zero."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=0.0,
                num_chords=10,
            )

    def test_delta_time_validation_negative(self):
        """Test that delta_time cannot be negative."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=-0.01,
                num_chords=10,
            )

    def test_static_movement_requires_num_chords(self):
        """Test that static Movement with num_steps=None requires num_chords."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_chords is None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_chords=None,
            )

    def test_static_movement_cannot_have_num_cycles(self):
        """Test that static Movement with num_steps=None cannot have num_cycles."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_cycles is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=3,
            )

    def test_non_static_movement_requires_num_cycles(self):
        """Test that non-static Movement with num_steps=None requires num_cycles."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_cycles is None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=None,
            )

    def test_non_static_movement_cannot_have_num_chords(self):
        """Test that non-static Movement with num_steps=None cannot have num_chords."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_chords is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=3,
                num_chords=10,
            )

    def test_num_steps_overrides_num_cycles_and_num_chords(self):
        """Test that when num_steps is set, num_cycles and num_chords must be None."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is set and num_cycles is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_steps=100,
                num_cycles=3,
            )

        # Should raise error when num_steps is set and num_chords is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_steps=100,
                num_chords=10,
            )

    def test_num_cycles_validation(self):
        """Test num_cycles parameter validation."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with valid positive integer.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=5,
        )
        self.assertEqual(movement.num_cycles, 5)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        num_cycles=invalid_value,
                    )

    def test_num_chords_validation(self):
        """Test num_chords parameter validation."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with valid positive integer.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_chords=15,
        )
        self.assertEqual(movement.num_chords, 15)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "ten"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        num_chords=invalid_value,
                    )

    def test_num_steps_validation(self):
        """Test num_steps parameter validation."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with valid positive integer.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=200,
        )
        self.assertEqual(movement.num_steps, 200)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "hundred"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        num_steps=invalid_value,
                    )

    def test_static_property_for_static_movement(self):
        """Test that static property returns True for static Movement."""
        movement = self.static_movement
        self.assertTrue(movement.static)

    def test_static_property_for_non_static_movement(self):
        """Test that static property returns False for non-static Movement."""
        movement = self.basic_movement
        self.assertFalse(movement.static)

    def test_max_period_for_static_movement(self):
        """Test that max_period returns 0.0 for static Movement."""
        movement = self.static_movement
        self.assertEqual(movement.max_period, 0.0)

    def test_max_period_for_non_static_movement(self):
        """Test that max_period returns correct value for non-static Movement."""
        movement = self.basic_movement
        # The basic_movement has period of 2.0 for all motion.
        self.assertEqual(movement.max_period, 2.0)

    def test_airplanes_generation(self):
        """Test that airplanes are generated correctly."""
        movement = self.basic_movement

        # Check that airplanes attribute is a list of lists.
        self.assertIsInstance(movement.airplanes, list)
        self.assertEqual(len(movement.airplanes), len(movement.airplane_movements))

        # Check that each element is a list of Airplanes.
        for airplane_list in movement.airplanes:
            self.assertIsInstance(airplane_list, list)
            self.assertEqual(len(airplane_list), movement.num_steps)
            for airplane in airplane_list:
                self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_operating_points_generation(self):
        """Test that operating_points are generated correctly."""
        movement = self.basic_movement

        # Check that operating_points attribute is a list.
        self.assertIsInstance(movement.operating_points, list)
        self.assertEqual(len(movement.operating_points), movement.num_steps)

        # Check that each element is an OperatingPoint.
        for operating_point in movement.operating_points:
            self.assertIsInstance(operating_point, ps.operating_point.OperatingPoint)

    def test_delta_time_automatic_calculation(self):
        """Test that delta_time is automatically calculated when not provided."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=3,
        )

        # Check that delta_time was calculated and is positive.
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)

    def test_num_steps_automatic_calculation_for_static(self):
        """Test that num_steps is automatically calculated for static Movement."""
        movement = self.static_movement

        # Check that num_steps was calculated and is positive.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

        # Check that it's based on num_chords.
        self.assertIsNotNone(movement.num_chords)

    def test_num_steps_automatic_calculation_for_non_static(self):
        """Test that num_steps is automatically calculated for non-static Movement."""
        movement = self.basic_movement

        # Check that num_steps was calculated and is positive.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

        # Check that it's based on num_cycles.
        self.assertIsNotNone(movement.num_cycles)

    def test_explicit_num_steps_for_static(self):
        """Test that explicit num_steps works for static Movement."""
        movement = self.static_movement_with_explicit_num_steps
        self.assertEqual(movement.num_steps, 5)
        self.assertIsNone(movement.num_cycles)
        self.assertIsNone(movement.num_chords)

    def test_explicit_num_steps_for_non_static(self):
        """Test that explicit num_steps works for non-static Movement."""
        movement = self.non_static_movement_with_explicit_num_steps
        self.assertEqual(movement.num_steps, 10)
        self.assertIsNone(movement.num_cycles)
        self.assertIsNone(movement.num_chords)

    def test_custom_delta_time(self):
        """Test that custom delta_time is used correctly."""
        movement = self.movement_with_custom_delta_time
        self.assertEqual(movement.delta_time, 0.05)

    def test_multiple_airplanes(self):
        """Test Movement with multiple AirplaneMovements."""
        movement = self.movement_with_multiple_airplanes

        # Check that we have multiple airplane_movements.
        self.assertEqual(len(movement.airplane_movements), 2)

        # Check that airplanes list has correct structure.
        self.assertEqual(len(movement.airplanes), 2)
        for airplane_list in movement.airplanes:
            self.assertEqual(len(airplane_list), movement.num_steps)

    def test_max_period_with_multiple_airplanes(self):
        """Test that max_period returns maximum across all AirplaneMovements."""
        movement = self.movement_with_multiple_airplanes

        # The movement has one static (period 0.0) and one with period 2.0.
        # Should return 2.0.
        self.assertEqual(movement.max_period, 2.0)

    def test_delta_time_averaging_with_multiple_airplanes(self):
        """Test that delta_time is averaged across multiple Airplanes when auto-calculated."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture(),
            airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_chords=10,
        )

        # Check that delta_time was calculated.
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)

    def test_num_steps_calculation_uses_ceil_for_static(self):
        """Test that num_steps calculation uses math.ceil for static Movement."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_chords=5,
        )

        # num_steps should be at least num_chords * c_ref / (delta_time * vCg__E).
        # Since we use math.ceil, it should be an integer.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

    def test_num_steps_calculation_uses_ceil_for_non_static(self):
        """Test that num_steps calculation uses math.ceil for non-static Movement."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_cycles=3,
        )

        # num_steps should be at least num_cycles * max_period / delta_time.
        # Since we use math.ceil, it should be an integer.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

    def test_type_4_to_5_transition_raises_error(self):
        """Test that a Wing transitioning from type 4 to type 5 symmetry raises
        an error."""
        # Create a type 4 Wing (symmetric=True, coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_4_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a WingMovement with rotation that will cause the symmetry plane
        # to become non-coincident (type 4->5 transition).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        )

        # Create an AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Attempting to create a Movement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.Movement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                num_cycles=1,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_type_3_to_2_transition_raises_error(self):
        """Test that a Wing transitioning from type 3 to type 2 symmetry raises
        an error."""
        # Create a type 3 Wing (mirror_only=True, non-coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_3_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a WingMovement with rotation that will cause the symmetry plane
        # to become coincident (type 3->2 transition).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
            ampLer_Gs_Cgs=(0.0, 0.5, 0.0),
            periodLer_Gs_Cgs=(0.0, 1.0, 0.0),
        )

        # Create an AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Attempting to create a Movement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.Movement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                delta_time=0.25,
                num_steps=5,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_type_2_to_3_transition_raises_error(self):
        """Test that a Wing transitioning from type 2 to type 3 symmetry raises
        an error."""
        # Create a type 2 Wing (mirror_only=True, coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_2_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a WingMovement with rotation that will cause the symmetry plane
        # to become non-coincident (type 2->3 transition).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        )

        # Create an AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Attempting to create a Movement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.Movement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                num_cycles=1,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_static_movement_with_symmetric_wing_succeeds(self):
        """Test that a Movement with a symmetric Wing but no motion succeeds."""
        # Create a type 4 Wing.
        base_wing = geometry_fixtures.make_type_4_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create static WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a static WingMovement (no rotation or translation).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
        )

        # Create a static AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Creating a Movement should succeed without raising an error.
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            num_chords=3,
        )

        # Verify the Movement was created successfully.
        self.assertIsInstance(movement, ps.movements.movement.Movement)


if __name__ == "__main__":
    unittest.main()
