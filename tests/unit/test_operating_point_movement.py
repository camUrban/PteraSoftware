"""This module contains classes to test OperatingPointMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    operating_point_fixtures,
    operating_point_movement_fixtures,
)


class TestOperatingPointMovement(unittest.TestCase):
    """This is a class with functions to test OperatingPointMovements."""

    def test_is_subclass_of_core(self):
        """Test that OperatingPointMovement is a subclass of
        CoreOperatingPointMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.operating_point_movement.OperatingPointMovement,
                ps._core.CoreOperatingPointMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that OperatingPointMovement instantiation returns an
        OperatingPointMovement.
        """
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point,
            )
        )
        self.assertIsInstance(
            operating_point_movement,
            ps.movements.operating_point_movement.OperatingPointMovement,
        )

    def test_generate_operating_points_returns_operating_points(self):
        """Test that generate_operating_points returns OperatingPoints when called
        through the public class.
        """
        operating_point_movement = (
            operating_point_movement_fixtures.make_sine_spacing_operating_point_movement_fixture()
        )
        operating_points = operating_point_movement.generate_operating_points(
            num_steps=5, delta_time=0.01
        )
        self.assertEqual(len(operating_points), 5)
        for operating_point in operating_points:
            self.assertIsInstance(
                operating_point,
                ps.operating_point.OperatingPoint,
            )


if __name__ == "__main__":
    unittest.main()
