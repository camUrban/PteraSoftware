"""This module contains classes to test AeroelasticOperatingPointMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    aeroelastic_operating_point_movement_fixtures,
    operating_point_fixtures,
)


class TestAeroelasticOperatingPointMovement(unittest.TestCase):
    """This is a class with functions to test AeroelasticOperatingPointMovements."""

    def test_is_subclass_of_core(self):
        """Test that AeroelasticOperatingPointMovement is a subclass of
        CoreOperatingPointMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement,
                ps._core.CoreOperatingPointMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that AeroelasticOperatingPointMovement instantiation returns an
        AeroelasticOperatingPointMovement.
        """
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        aeroelastic_operating_point_movement = ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
            base_operating_point=base_operating_point,
        )
        self.assertIsInstance(
            aeroelastic_operating_point_movement,
            ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement,
        )

    def test_generate_operating_points_returns_operating_points(self):
        """Test that generate_operating_points returns OperatingPoints when called
        through the public class.
        """
        aeroelastic_operating_point_movement = (
            aeroelastic_operating_point_movement_fixtures.make_sine_spacing_aeroelastic_operating_point_movement_fixture()
        )
        operating_points = (
            aeroelastic_operating_point_movement.generate_operating_points(
                num_steps=5, delta_time=0.01
            )
        )
        self.assertEqual(len(operating_points), 5)
        for operating_point in operating_points:
            self.assertIsInstance(
                operating_point,
                ps.operating_point.OperatingPoint,
            )


if __name__ == "__main__":
    unittest.main()
