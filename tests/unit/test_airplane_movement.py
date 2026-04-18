"""This module contains classes to test AirplaneMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    airplane_movement_fixtures,
    geometry_fixtures,
    wing_movement_fixtures,
)


class TestAirplaneMovement(unittest.TestCase):
    """This is a class with functions to test AirplaneMovements."""

    def test_is_subclass_of_core(self):
        """Test that AirplaneMovement is a subclass of CoreAirplaneMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.airplane_movement.AirplaneMovement,
                ps._core.CoreAirplaneMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that AirplaneMovement instantiation returns an AirplaneMovement."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [wing_movement_fixtures.make_static_wing_movement_fixture()]
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
        )
        self.assertIsInstance(
            airplane_movement,
            ps.movements.airplane_movement.AirplaneMovement,
        )

    def test_rejects_core_wing_movement_children(self):
        """Test that AirplaneMovement rejects CoreWingMovement instances."""
        from tests.unit.fixtures import core_wing_movement_fixtures

        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements = [
            core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
        ]
        with self.assertRaises(TypeError):
            ps.movements.airplane_movement.AirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
            )

    def test_generate_airplanes_returns_airplanes(self):
        """Test that generate_airplanes returns Airplanes when called through the
        public class.
        """
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        airplanes = airplane_movement.generate_airplanes(num_steps=5, delta_time=0.01)
        self.assertEqual(len(airplanes), 5)
        for airplane in airplanes:
            self.assertIsInstance(
                airplane,
                ps.geometry.airplane.Airplane,
            )


if __name__ == "__main__":
    unittest.main()
