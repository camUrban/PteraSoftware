"""This module contains classes to test FreeFlightWingCrossSectionMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    free_flight_wing_cross_section_movement_fixtures,
    geometry_fixtures,
)


class TestFreeFlightWingCrossSectionMovement(unittest.TestCase):
    """This is a class with functions to test FreeFlightWingCrossSectionMovements."""

    def test_is_subclass_of_core(self):
        """Test that FreeFlightWingCrossSectionMovement is a subclass of
        CoreWingCrossSectionMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement,
                ps._core.CoreWingCrossSectionMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that FreeFlightWingCrossSectionMovement instantiation returns a
        FreeFlightWingCrossSectionMovement.
        """
        base_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        free_flight_wing_cross_section_movement = ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
        )
        self.assertIsInstance(
            free_flight_wing_cross_section_movement,
            ps.movements.free_flight_wing_cross_section_movement.FreeFlightWingCrossSectionMovement,
        )

    def test_generate_wing_cross_sections_returns_wing_cross_sections(self):
        """Test that generate_wing_cross_sections returns WingCrossSections when called
        through the public class.
        """
        free_flight_wing_cross_section_movement = (
            free_flight_wing_cross_section_movement_fixtures.make_basic_free_flight_wing_cross_section_movement_fixture()
        )
        wing_cross_sections = (
            free_flight_wing_cross_section_movement.generate_wing_cross_sections(
                num_steps=5, delta_time=0.01
            )
        )
        self.assertEqual(len(wing_cross_sections), 5)
        for wing_cross_section in wing_cross_sections:
            self.assertIsInstance(
                wing_cross_section,
                ps.geometry.wing_cross_section.WingCrossSection,
            )

    def test_generate_wing_cross_section_at_time_step_returns_wing_cross_section(self):
        """Test that generate_wing_cross_section_at_time_step returns a
        WingCrossSection.
        """
        free_flight_wing_cross_section_movement = (
            free_flight_wing_cross_section_movement_fixtures.make_basic_free_flight_wing_cross_section_movement_fixture()
        )
        wing_cross_section = free_flight_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=2, delta_time=0.01
        )
        self.assertIsInstance(
            wing_cross_section,
            ps.geometry.wing_cross_section.WingCrossSection,
        )
