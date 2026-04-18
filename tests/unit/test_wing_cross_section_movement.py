"""This module contains classes to test WingCrossSectionMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    geometry_fixtures,
    wing_cross_section_movement_fixtures,
)


class TestWingCrossSectionMovement(unittest.TestCase):
    """This is a class with functions to test WingCrossSectionMovements."""

    def test_is_subclass_of_core(self):
        """Test that WingCrossSectionMovement is a subclass of
        CoreWingCrossSectionMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement,
                ps._core.CoreWingCrossSectionMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that WingCrossSectionMovement instantiation returns a
        WingCrossSectionMovement.
        """
        base_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        wing_cross_section_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_cross_section,
            )
        )
        self.assertIsInstance(
            wing_cross_section_movement,
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement,
        )

    def test_generate_wing_cross_sections_returns_wing_cross_sections(self):
        """Test that generate_wing_cross_sections returns WingCrossSections when
        called through the public class.
        """
        wing_cross_section_movement = (
            wing_cross_section_movement_fixtures.make_basic_wing_cross_section_movement_fixture()
        )
        wing_cross_sections = wing_cross_section_movement.generate_wing_cross_sections(
            num_steps=5, delta_time=0.01
        )
        self.assertEqual(len(wing_cross_sections), 5)
        for wing_cross_section in wing_cross_sections:
            self.assertIsInstance(
                wing_cross_section,
                ps.geometry.wing_cross_section.WingCrossSection,
            )


if __name__ == "__main__":
    unittest.main()
