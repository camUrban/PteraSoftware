"""This module contains classes to test FreeFlightAirplaneMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    free_flight_airplane_movement_fixtures,
    free_flight_wing_movement_fixtures,
    geometry_fixtures,
    wing_movement_fixtures,
)


class TestFreeFlightAirplaneMovement(unittest.TestCase):
    """This is a class with functions to test FreeFlightAirplaneMovements."""

    def test_is_subclass_of_core(self):
        """Test that FreeFlightAirplaneMovement is a subclass of
        CoreAirplaneMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement,
                ps._core.CoreAirplaneMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that FreeFlightAirplaneMovement instantiation returns a
        FreeFlightAirplaneMovement.
        """
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [
            free_flight_wing_movement_fixtures.make_basic_free_flight_wing_movement_fixture()
        ]
        free_flight_airplane_movement = (
            ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
            )
        )
        self.assertIsInstance(
            free_flight_airplane_movement,
            ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement,
        )

    def test_rejects_non_free_flight_wing_movement_children(self):
        """Test that FreeFlightAirplaneMovement rejects WingMovement instances that
        are not FreeFlightWingMovements.
        """
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]
        with self.assertRaises(TypeError):
            ps.movements.free_flight_airplane_movement.FreeFlightAirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
            )

    def test_wing_movements_property_returns_tuple(self):
        """Test that the wing_movements property returns a tuple of the
        FreeFlightAirplaneMovement's wing movements.
        """
        free_flight_airplane_movement = (
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        )
        self.assertIsInstance(free_flight_airplane_movement.wing_movements, tuple)
        self.assertEqual(len(free_flight_airplane_movement.wing_movements), 1)

    def test_generate_airplanes_returns_airplanes(self):
        """Test that generate_airplanes returns Airplanes when called through the
        public class.
        """
        free_flight_airplane_movement = (
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        )
        airplanes = free_flight_airplane_movement.generate_airplanes(
            num_steps=5, delta_time=0.01
        )
        self.assertEqual(len(airplanes), 5)
        for airplane in airplanes:
            self.assertIsInstance(
                airplane,
                ps.geometry.airplane.Airplane,
            )

    def test_generate_airplane_at_time_step_returns_airplane(self):
        """Test that generate_airplane_at_time_step returns an Airplane."""
        free_flight_airplane_movement = (
            free_flight_airplane_movement_fixtures.make_basic_free_flight_airplane_movement_fixture()
        )
        airplane = free_flight_airplane_movement.generate_airplane_at_time_step(
            step=2, delta_time=0.01
        )
        self.assertIsInstance(
            airplane,
            ps.geometry.airplane.Airplane,
        )
