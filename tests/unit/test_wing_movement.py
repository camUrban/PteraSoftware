"""This module contains classes to test WingMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    geometry_fixtures,
    wing_cross_section_movement_fixtures,
    wing_movement_fixtures,
)


class TestWingMovement(unittest.TestCase):
    """This is a class with functions to test WingMovements."""

    def test_is_subclass_of_core(self):
        """Test that WingMovement is a subclass of CoreWingMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.wing_movement.WingMovement,
                ps._core.CoreWingMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that WingMovement instantiation returns a WingMovement."""
        base_wing = geometry_fixtures.make_origin_wing_fixture()
        wing_cross_section_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
            wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
        ]
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wing_cross_section_movements,
        )
        self.assertIsInstance(
            wing_movement,
            ps.movements.wing_movement.WingMovement,
        )

    def test_rejects_core_wing_cross_section_movement_children(self):
        """Test that WingMovement rejects CoreWingCrossSectionMovement instances."""
        from tests.unit.fixtures import core_wing_cross_section_movement_fixtures

        base_wing = geometry_fixtures.make_origin_wing_fixture()
        wing_cross_section_movements = [
            core_wing_cross_section_movement_fixtures.make_static_core_wing_cross_section_movement_fixture(),
            core_wing_cross_section_movement_fixtures.make_static_tip_core_wing_cross_section_movement_fixture(),
        ]
        with self.assertRaises(TypeError):
            ps.movements.wing_movement.WingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wing_cross_section_movements,
            )

    def test_generate_wings_returns_wings(self):
        """Test that generate_wings returns Wings when called through the public
        class.
        """
        wing_movement = wing_movement_fixtures.make_basic_wing_movement_fixture()
        wings = wing_movement.generate_wings(num_steps=5, delta_time=0.01)
        self.assertEqual(len(wings), 5)
        for wing in wings:
            self.assertIsInstance(
                wing,
                ps.geometry.wing.Wing,
            )


if __name__ == "__main__":
    unittest.main()
