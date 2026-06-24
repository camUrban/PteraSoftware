"""This module contains classes to test FreeFlightMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    airplane_movement_fixtures,
    free_flight_airplane_movement_fixtures,
    free_flight_movement_fixtures,
    free_flight_operating_point_movement_fixtures,
    geometry_fixtures,
    operating_point_fixtures,
)


class TestFreeFlightMovement(unittest.TestCase):
    """This is a class with functions to test FreeFlightMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all FreeFlightMovement tests."""
        cls.basic_free_flight_movement = (
            free_flight_movement_fixtures.make_basic_free_flight_movement_fixture()
        )
        cls.static_free_flight_movement = (
            free_flight_movement_fixtures.make_static_free_flight_movement_fixture()
        )
        cls.free_flight_movement_with_multiple_airplanes = (
            free_flight_movement_fixtures.make_free_flight_movement_with_multiple_airplanes_fixture()
        )

    def test_is_subclass_of_core_movement(self):
        """Test that FreeFlightMovement is a subclass of CoreMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.free_flight_movement.FreeFlightMovement,
                ps._core.CoreMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that FreeFlightMovement instantiation returns a FreeFlightMovement."""
        self.assertIsInstance(
            self.basic_free_flight_movement,
            ps.movements.free_flight_movement.FreeFlightMovement,
        )

    def test_rejects_non_free_flight_airplane_movement_children(self):
        """Test that FreeFlightMovement rejects AirplaneMovement instances that are not
        FreeFlightAirplaneMovements.
        """
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.free_flight_movement.FreeFlightMovement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=0.1,
                prescribed_num_steps=2,
                free_num_steps=3,
            )

    def test_rejects_non_free_flight_operating_point_movement(self):
        """Test that FreeFlightMovement rejects OperatingPointMovement instances that
        are not FreeFlightOperatingPointMovements.
        """
        airplane_movements = [
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.free_flight_movement.FreeFlightMovement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=0.1,
                prescribed_num_steps=2,
                free_num_steps=3,
            )

    def test_prescribed_num_steps_validation(self):
        """Test prescribed_num_steps parameter validation."""
        airplane_movements = [
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        ]
        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        invalid_values = [0, -5, 2.5, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.free_flight_movement.FreeFlightMovement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        delta_time=0.1,
                        prescribed_num_steps=invalid_value,
                        free_num_steps=3,
                    )

    def test_free_num_steps_validation(self):
        """Test free_num_steps parameter validation."""
        airplane_movements = [
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        ]
        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        invalid_values = [0, -5, 2.5, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.free_flight_movement.FreeFlightMovement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        delta_time=0.1,
                        prescribed_num_steps=2,
                        free_num_steps=invalid_value,
                    )

    def test_delta_time_validation(self):
        """Test that delta_time must be a positive number."""
        airplane_movements = [
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        ]
        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        invalid_values = [0.0, -0.01, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.free_flight_movement.FreeFlightMovement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        delta_time=invalid_value,
                        prescribed_num_steps=2,
                        free_num_steps=3,
                    )

    def test_num_steps_equals_sum_of_prescribed_and_free(self):
        """Test that num_steps equals prescribed_num_steps plus free_num_steps."""
        self.assertEqual(self.basic_free_flight_movement.num_steps, 5)
        self.assertEqual(
            self.basic_free_flight_movement.num_steps,
            self.basic_free_flight_movement.prescribed_num_steps
            + self.basic_free_flight_movement.free_num_steps,
        )

    def test_prescribed_num_steps_property(self):
        """Test that prescribed_num_steps is stored and accessible via the property."""
        self.assertEqual(self.basic_free_flight_movement.prescribed_num_steps, 2)

    def test_free_num_steps_property(self):
        """Test that free_num_steps is stored and accessible via the property."""
        self.assertEqual(self.basic_free_flight_movement.free_num_steps, 3)

    def test_airplane_movements_returns_tuple_of_free_flight_airplane_movements(self):
        """Test that airplane_movements returns a tuple of
        FreeFlightAirplaneMovements.
        """
        airplane_movements = self.basic_free_flight_movement.airplane_movements
        self.assertIsInstance(airplane_movements, tuple)
        for airplane_movement in airplane_movements:
            self.assertIsInstance(
                airplane_movement,
                ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement,
            )

    def test_operating_point_movement_returns_free_flight_operating_point_movement(
        self,
    ):
        """Test that operating_point_movement returns a
        FreeFlightOperatingPointMovement.
        """
        self.assertIsInstance(
            self.basic_free_flight_movement.operating_point_movement,
            ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement,
        )

    def test_airplanes_structure(self):
        """Test that airplanes is a tuple of tuples of Airplanes with the correct
        dimensions.
        """
        airplanes = self.basic_free_flight_movement.airplanes
        self.assertIsInstance(airplanes, tuple)

        # The first index runs over the FreeFlightAirplaneMovements.
        self.assertEqual(
            len(airplanes),
            len(self.basic_free_flight_movement.airplane_movements),
        )

        # The second index runs over the time steps.
        for airplane_list in airplanes:
            self.assertIsInstance(airplane_list, tuple)
            self.assertEqual(
                len(airplane_list), self.basic_free_flight_movement.num_steps
            )
            for airplane in airplane_list:
                self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_airplanes_structure_with_multiple_airplanes(self):
        """Test that airplanes has one inner tuple per FreeFlightAirplaneMovement."""
        airplanes = self.free_flight_movement_with_multiple_airplanes.airplanes
        self.assertEqual(len(airplanes), 2)
        for airplane_list in airplanes:
            self.assertEqual(
                len(airplane_list),
                self.free_flight_movement_with_multiple_airplanes.num_steps,
            )

    def test_static_property_for_static_free_flight_movement(self):
        """Test that static returns True for a static FreeFlightMovement."""
        self.assertTrue(self.static_free_flight_movement.static)

    def test_static_property_for_non_static_free_flight_movement(self):
        """Test that static returns False for a non static FreeFlightMovement."""
        self.assertFalse(self.basic_free_flight_movement.static)

    def test_max_period_for_static_free_flight_movement(self):
        """Test that max_period returns 0.0 for a static FreeFlightMovement."""
        self.assertEqual(self.static_free_flight_movement.max_period, 0.0)

    def test_max_period_for_non_static_free_flight_movement(self):
        """Test that max_period returns the correct value for a non static
        FreeFlightMovement.
        """
        # The basic FreeFlightMovement's prescribed motion has a period of 2.0.
        self.assertEqual(self.basic_free_flight_movement.max_period, 2.0)

    def test_min_period_for_non_static_free_flight_movement(self):
        """Test that min_period returns the correct value for a non static
        FreeFlightMovement.
        """
        # The basic FreeFlightMovement's prescribed motion has a period of 2.0.
        self.assertEqual(self.basic_free_flight_movement.min_period, 2.0)

    def test_lcm_period_for_non_static_free_flight_movement(self):
        """Test that lcm_period returns the correct value for a non static
        FreeFlightMovement.
        """
        # The basic FreeFlightMovement's prescribed motion has a period of 2.0, so the
        # LCM of its identical periods is that period.
        self.assertEqual(self.basic_free_flight_movement.lcm_period, 2.0)

    def test_max_wake_rows_default_none(self):
        """Test that max_wake_rows defaults to None."""
        self.assertIsNone(self.basic_free_flight_movement.max_wake_rows)

    def test_max_wake_rows_stored_and_accessible(self):
        """Test that max_wake_rows is stored and accessible via the property."""
        airplane_movements = [
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        ]
        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        free_flight_movement = ps.movements.free_flight_movement.FreeFlightMovement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            prescribed_num_steps=2,
            free_num_steps=3,
            max_wake_rows=5,
        )
        self.assertEqual(free_flight_movement.max_wake_rows, 5)

    def test_max_wake_rows_validation(self):
        """Test that max_wake_rows must be a positive int."""
        airplane_movements = [
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        ]
        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        invalid_values = [0, -1, 2.5, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.free_flight_movement.FreeFlightMovement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        delta_time=0.1,
                        prescribed_num_steps=2,
                        free_num_steps=3,
                        max_wake_rows=invalid_value,
                    )

    def test_symmetry_change_raises_error(self):
        """Test that a Wing transitioning between symmetry types raises an error."""
        # Create a type 4 Wing (symmetric=True, coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_4_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry. The
        # base Airplane is the first in a simulation, so its Cg_GP1_CgP1 is all zeros.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create FreeFlightWingCrossSectionMovements using the actual WingCrossSections.
        wing_cross_section_movements = [
            ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
                base_wing_cross_section=wing_cross_section,
            )
            for wing_cross_section in processed_wing.wing_cross_sections
        ]

        # Create a FreeFlightWingMovement with rotation that will cause the symmetry
        # plane to become non coincident (a type 4 to type 5 transition).
        wing_movement = ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wing_cross_section_movements,
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        )

        airplane_movement = (
            ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=[wing_movement],
            )
        )

        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        # Attempting to create a FreeFlightMovement should raise a ValueError. With a
        # delta_time of 0.25 and a period of 1.0, the second time step lands at the peak
        # of the rotation oscillation, so the symmetry plane is non coincident there.
        with self.assertRaises(ValueError) as context:
            ps.movements.free_flight_movement.FreeFlightMovement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                delta_time=0.25,
                prescribed_num_steps=1,
                free_num_steps=1,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_static_motion_with_symmetric_wing_succeeds(self):
        """Test that a FreeFlightMovement with a symmetric Wing but no motion
        succeeds.
        """
        # Create a type 4 Wing (symmetric=True, coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_4_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create static FreeFlightWingCrossSectionMovements using the actual
        # WingCrossSections.
        wing_cross_section_movements = [
            ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
                base_wing_cross_section=wing_cross_section,
            )
            for wing_cross_section in processed_wing.wing_cross_sections
        ]

        # Create a static FreeFlightWingMovement (no rotation or translation).
        wing_movement = ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wing_cross_section_movements,
        )

        airplane_movement = (
            ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=[wing_movement],
            )
        )

        operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )

        # Creating a FreeFlightMovement should succeed without raising an error.
        free_flight_movement = ps.movements.free_flight_movement.FreeFlightMovement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            prescribed_num_steps=2,
            free_num_steps=3,
        )

        self.assertIsInstance(
            free_flight_movement,
            ps.movements.free_flight_movement.FreeFlightMovement,
        )


class TestFreeFlightMovementImmutability(unittest.TestCase):
    """Tests for FreeFlightMovement attribute immutability."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all immutability tests."""
        cls.basic_free_flight_movement = (
            free_flight_movement_fixtures.make_basic_free_flight_movement_fixture()
        )

    def test_immutable_prescribed_num_steps_property(self):
        """Test that prescribed_num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_free_flight_movement.prescribed_num_steps = 100

    def test_immutable_free_num_steps_property(self):
        """Test that free_num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_free_flight_movement.free_num_steps = 100

    def test_immutable_airplanes_property(self):
        """Test that airplanes property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_free_flight_movement.airplanes = ()

    def test_immutable_airplane_movements_property(self):
        """Test that airplane_movements property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_free_flight_movement.airplane_movements = ()
