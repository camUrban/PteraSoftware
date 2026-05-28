"""This module contains classes to test the AeroelasticAirplaneMovement class."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import movement_fixtures


class TestAeroelasticAirplaneMovementStandardWing(unittest.TestCase):
    """This class contains unit tests for AeroelasticAirplaneMovement with a standard
    WingMovement child, exercising the else branch in generate_airplane_at_time_step."""

    def setUp(self):
        """Set up an AeroelasticAirplaneMovement backed by a standard WingMovement."""
        movement = (
            movement_fixtures.make_aeroelastic_movement_with_standard_wing_fixture()
        )
        self.airplane_movement = movement.airplane_movements[0]
        self.delta_time = movement.delta_time

    def test_generate_airplane_at_time_step_standard_wing_returns_airplane(self):
        """Test that generate_airplane_at_time_step returns an Airplane when the wing
        is backed by a standard WingMovement (the else branch)."""
        result = self.airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=self.delta_time
        )

        self.assertIsInstance(result, ps.geometry.airplane.Airplane)

    def test_generate_airplane_at_time_step_standard_wing_returns_one_wing(self):
        """Test that the Airplane returned via the standard WingMovement path contains
        exactly one Wing."""
        result = self.airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=self.delta_time
        )

        self.assertEqual(len(result.wings), 1)

    def test_generate_airplane_at_time_step_standard_wing_wing_has_panels(self):
        """Test that the Wing produced via the standard WingMovement path has its
        panels populated (not None)."""
        result = self.airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=self.delta_time
        )

        self.assertIsNotNone(result.wings[0].panels)

    def test_generate_airplane_at_time_step_standard_wing_deformation_ignored(self):
        """Test that passing wing_deformation_angles_ixyz=None produces the same
        Airplane as calling without deformation for a standard WingMovement."""
        result_no_kwarg = self.airplane_movement.generate_airplane_at_time_step(
            step=1, delta_time=self.delta_time
        )
        result_none_kwarg = self.airplane_movement.generate_airplane_at_time_step(
            step=1,
            delta_time=self.delta_time,
            wing_deformation_angles_ixyz=None,
        )

        self.assertEqual(len(result_no_kwarg.wings), len(result_none_kwarg.wings))
        self.assertIsNotNone(result_none_kwarg.wings[0].panels)
