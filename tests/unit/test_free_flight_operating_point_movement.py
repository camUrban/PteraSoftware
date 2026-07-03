"""This module contains classes to test FreeFlightOperatingPointMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    free_flight_operating_point_movement_fixtures,
    operating_point_fixtures,
)


class TestFreeFlightOperatingPointMovement(unittest.TestCase):
    """This is a class with functions to test FreeFlightOperatingPointMovements."""

    def test_is_subclass_of_core(self):
        """Test that FreeFlightOperatingPointMovement is a subclass of
        CoreOperatingPointMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement,
                ps._core.CoreOperatingPointMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that FreeFlightOperatingPointMovement instantiation returns a
        FreeFlightOperatingPointMovement.
        """
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        free_flight_operating_point_movement = ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
            base_operating_point=base_operating_point,
        )
        self.assertIsInstance(
            free_flight_operating_point_movement,
            ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement,
        )

    def test_operating_points_initialized_with_base_operating_point(self):
        """Test that operating_points starts as a list holding only the base
        OperatingPoint.
        """
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        free_flight_operating_point_movement = ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
            base_operating_point=base_operating_point,
        )
        self.assertEqual(len(free_flight_operating_point_movement.operating_points), 1)
        self.assertIs(
            free_flight_operating_point_movement.operating_points[0],
            base_operating_point,
        )

    def test_operating_points_is_mutable(self):
        """Test that operating_points supports the solver appending new
        OperatingPoints.
        """
        free_flight_operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )
        new_operating_point = (
            operating_point_fixtures.make_high_speed_operating_point_fixture()
        )
        free_flight_operating_point_movement.operating_points.append(
            new_operating_point
        )
        self.assertEqual(len(free_flight_operating_point_movement.operating_points), 2)
        self.assertIs(
            free_flight_operating_point_movement.operating_points[1],
            new_operating_point,
        )

    def test_generate_operating_points_returns_operating_points(self):
        """Test that generate_operating_points returns OperatingPoints when called
        through the public class.
        """
        free_flight_operating_point_movement = (
            free_flight_operating_point_movement_fixtures.make_basic_free_flight_operating_point_movement_fixture()
        )
        operating_points = (
            free_flight_operating_point_movement.generate_operating_points(
                num_steps=5, delta_time=0.01
            )
        )
        self.assertEqual(len(operating_points), 5)
        for operating_point in operating_points:
            self.assertIsInstance(
                operating_point,
                ps.operating_point.OperatingPoint,
            )
