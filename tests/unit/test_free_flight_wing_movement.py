"""This module contains classes to test FreeFlightWingMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    free_flight_wing_movement_fixtures,
    geometry_fixtures,
    wing_cross_section_movement_fixtures,
)


class TestFreeFlightWingMovement(unittest.TestCase):
    """This is a class with functions to test FreeFlightWingMovements."""

    def test_is_subclass_of_core(self):
        """Test that FreeFlightWingMovement is a subclass of CoreWingMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.free_flight_wing_movement.FreeFlightWingMovement,
                ps._core.CoreWingMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that FreeFlightWingMovement instantiation returns a
        FreeFlightWingMovement.
        """
        base_wing = geometry_fixtures.make_origin_wing_fixture()
        wing_cross_section_movements = [
            ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
                base_wing_cross_section=wing_cross_section
            )
            for wing_cross_section in base_wing.wing_cross_sections
        ]
        free_flight_wing_movement = (
            ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wing_cross_section_movements,
            )
        )
        self.assertIsInstance(
            free_flight_wing_movement,
            ps.movements.free_flight_wing_movement.FreeFlightWingMovement,
        )

    def test_rejects_non_free_flight_wing_cross_section_movement_children(self):
        """Test that FreeFlightWingMovement rejects WingCrossSectionMovement instances
        that are not FreeFlightWingCrossSectionMovements.
        """
        base_wing = geometry_fixtures.make_origin_wing_fixture()
        wing_cross_section_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
            wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
        ]
        with self.assertRaises(TypeError):
            ps.movements.free_flight_wing_movement.FreeFlightWingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wing_cross_section_movements,
            )

    def test_generate_wings_returns_wings(self):
        """Test that generate_wings returns Wings when called through the public
        class.
        """
        free_flight_wing_movement = (
            free_flight_wing_movement_fixtures.make_basic_free_flight_wing_movement_fixture()
        )
        wings = free_flight_wing_movement.generate_wings(num_steps=5, delta_time=0.01)
        self.assertEqual(len(wings), 5)
        for wing in wings:
            self.assertIsInstance(
                wing,
                ps.geometry.wing.Wing,
            )

    def test_generate_wing_at_time_step_returns_wing(self):
        """Test that generate_wing_at_time_step returns a Wing."""
        free_flight_wing_movement = (
            free_flight_wing_movement_fixtures.make_basic_free_flight_wing_movement_fixture()
        )
        wing = free_flight_wing_movement.generate_wing_at_time_step(
            step=2, delta_time=0.01
        )
        self.assertIsInstance(
            wing,
            ps.geometry.wing.Wing,
        )
