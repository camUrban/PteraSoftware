"""This module contains classes to test Movements and related functions."""

import unittest
from unittest.mock import patch

import pterasoftware as ps
from tests.unit.fixtures import (
    airplane_movement_fixtures,
    geometry_fixtures,
    movement_fixtures,
    operating_point_fixtures,
)


class TestMovement(unittest.TestCase):
    """This is a class with functions to test Movements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all Movement tests."""
        cls.static_movement = movement_fixtures.make_static_movement_fixture()
        cls.basic_movement = movement_fixtures.make_basic_movement_fixture()
        cls.static_movement_with_explicit_num_steps = (
            movement_fixtures.make_static_movement_with_explicit_num_steps_fixture()
        )
        cls.non_static_movement_with_explicit_num_steps = (
            movement_fixtures.make_non_static_movement_with_explicit_num_steps_fixture()
        )
        cls.movement_with_custom_delta_time = (
            movement_fixtures.make_movement_with_custom_delta_time_fixture()
        )
        cls.movement_with_multiple_airplanes = (
            movement_fixtures.make_movement_with_multiple_airplanes_fixture()
        )

    def test_airplane_movements_validation_not_list(self):
        """Test that airplane_movements must be a list."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.Movement(
                airplane_movements="not a list",
                operating_point_movement=operating_point_movement,
                num_chords=10,
            )

    def test_airplane_movements_validation_empty_list(self):
        """Test that airplane_movements must have at least one element."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=[],
                operating_point_movement=operating_point_movement,
                num_chords=10,
            )

    def test_airplane_movements_validation_invalid_element_type(self):
        """Test that all elements in airplane_movements must be AirplaneMovements."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.Movement(
                airplane_movements=["not an airplane movement"],
                operating_point_movement=operating_point_movement,
                num_chords=10,
            )

    def test_operating_point_movement_validation(self):
        """Test that operating_point_movement must be an OperatingPointMovement."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]

        with self.assertRaises(TypeError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement="not an operating point movement",
                num_chords=10,
            )

    def test_delta_time_validation_positive(self):
        """Test that delta_time must be positive if provided."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with positive delta_time works. Use num_steps=1 to speed up the test;
        # we only need to verify delta_time is accepted, not generate many airplanes.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.01,
            num_steps=1,
        )
        self.assertEqual(movement.delta_time, 0.01)

    def test_delta_time_validation_zero(self):
        """Test that delta_time cannot be zero."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=0.0,
                num_chords=10,
            )

    def test_delta_time_validation_negative(self):
        """Test that delta_time cannot be negative."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # noinspection PyTypeChecker
        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=-0.01,
                num_chords=10,
            )

    def test_static_movement_requires_num_chords(self):
        """Test that static Movement with num_steps=None requires num_chords."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_chords is None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_chords=None,
            )

    def test_static_movement_cannot_have_num_cycles(self):
        """Test that static Movement with num_steps=None cannot have num_cycles."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_cycles is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=3,
            )

    def test_non_static_movement_requires_num_cycles(self):
        """Test that non static Movement with num_steps=None requires num_cycles."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_cycles is None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=None,
            )

    def test_non_static_movement_cannot_have_num_chords(self):
        """Test that non static Movement with num_steps=None cannot have num_chords."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is None and num_chords is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_cycles=3,
                num_chords=10,
            )

    def test_num_steps_overrides_num_cycles_and_num_chords(self):
        """Test that when num_steps is set, num_cycles and num_chords must be None."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Should raise error when num_steps is set and num_cycles is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_steps=100,
                num_cycles=3,
            )

        # Should raise error when num_steps is set and num_chords is not None.
        with self.assertRaises(ValueError):
            ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                num_steps=100,
                num_chords=10,
            )

    def test_num_cycles_validation(self):
        """Test num_cycles parameter validation."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with valid positive integer. Use num_cycles=1 to speed up the test;
        # the validation logic doesn't depend on the specific value.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=1,
        )
        self.assertEqual(movement.num_cycles, 1)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        num_cycles=invalid_value,
                    )

    def test_num_chords_validation(self):
        """Test num_chords parameter validation."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with valid positive integer. Use num_chords=1 to speed up the test;
        # the validation logic doesn't depend on the specific value.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_chords=1,
        )
        self.assertEqual(movement.num_chords, 1)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "ten"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        num_chords=invalid_value,
                    )

    def test_num_steps_validation(self):
        """Test num_steps parameter validation."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Test with valid positive integer. Use num_steps=1 to speed up the test;
        # the validation logic doesn't depend on the specific value.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=1,
        )
        self.assertEqual(movement.num_steps, 1)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "hundred"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                # noinspection PyTypeChecker
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        num_steps=invalid_value,
                    )

    def test_static_property_for_static_movement(self):
        """Test that static property returns True for static Movement."""
        movement = self.static_movement
        self.assertTrue(movement.static)

    def test_static_property_for_non_static_movement(self):
        """Test that static property returns False for non static Movement."""
        movement = self.basic_movement
        self.assertFalse(movement.static)

    def test_max_period_for_static_movement(self):
        """Test that max_period returns 0.0 for static Movement."""
        movement = self.static_movement
        self.assertEqual(movement.max_period, 0.0)

    def test_max_period_for_non_static_movement(self):
        """Test that max_period returns correct value for non static Movement."""
        movement = self.basic_movement
        # The basic_movement has period of 2.0 for all motion.
        self.assertEqual(movement.max_period, 2.0)

    def test_lcm_period_for_static_movement(self):
        """Test that lcm_period returns 0.0 for static Movement."""
        movement = self.static_movement
        self.assertEqual(movement.lcm_period, 0.0)

    def test_lcm_period_for_single_period_movement(self):
        """Test that lcm_period returns correct value when all periods are the same."""
        movement = self.basic_movement
        # The basic_movement has period of 2.0 for all motion.
        # LCM of identical periods should equal that period.
        self.assertEqual(movement.lcm_period, 2.0)

    def test_lcm_period_with_multiple_wings_same_airplane(self):
        """Test that lcm_period collects all periods, not just max from each
        AirplaneMovement.

        This test creates a single Airplane with two Wings having different periods
        (3.0 s and 4.0 s). The correct LCM is 12.0 s. If the implementation only uses
        max_period from the AirplaneMovement, lcm_period would incorrectly return 4.0 s
        instead of 12.0 s.
        """
        # Create two Wings for the same Airplane.
        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing_1, base_wing_2],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Wing_1: tip WingCrossSectionMovement has period 3.0 s.
        wcs_movements_wing_1 = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_1 = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing_1,
            wing_cross_section_movements=wcs_movements_wing_1,
        )

        # Wing_2: tip WingCrossSectionMovement has period 4.0 s.
        wcs_movements_wing_2 = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(4.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_2 = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing_2,
            wing_cross_section_movements=wcs_movements_wing_2,
        )

        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement_1, wing_movement_2],
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use num_steps=1 instead of num_cycles=1 to speed up this test. The lcm_period
        # property is calculated from the Movement parameters (periods), not from the
        # generated Airplanes, so we only need to generate one Airplane to test the
        # period calculation logic.
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The max_period should be 4.0 (the max of 3.0 and 4.0).
        self.assertEqual(movement.max_period, 4.0)

        # The lcm_period should be LCM(3.0, 4.0) = 12.0, Not 4.0. This test will Fail if
        # lcm_period only uses max_period from each AirplaneMovement instead of
        # collecting all individual periods.
        self.assertEqual(movement.lcm_period, 12.0)

    def test_lcm_period_with_multiple_cross_sections_same_wing(self):
        """Test that lcm_period collects all periods from WingCrossSectionMovements.

        This test creates a single Wing with three WingCrossSections having different
        periods (root static, middle 3.0 s, tip 4.0 s). The correct LCM is 12.0 s. If
        the implementation only uses max_period from each WingMovement, lcm_period
        would incorrectly return 4.0 s instead of 12.0 s.
        """
        # Create a Wing with three WingCrossSections.
        test_airfoil = ps.geometry.airfoil.Airfoil(name="naca2412")

        root_wcs = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil,
            num_spanwise_panels=4,
            chord=2.0,
            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )

        middle_wcs = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil,
            num_spanwise_panels=4,
            chord=1.5,
            Lp_Wcsp_Lpp=(0.0, 1.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )

        tip_wcs = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil,
            num_spanwise_panels=None,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 1.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )

        base_wing = ps.geometry.wing.Wing(
            wing_cross_sections=[root_wcs, middle_wcs, tip_wcs],
            name="Test Wing",
        )

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Root WingCrossSectionMovement must be static.
        # Middle has period 3.0 s, tip has period 4.0 s.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing.wing_cross_sections[2],
                periodLp_Wcsp_Lpp=(4.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wcs_movements,
        )

        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use num_steps=1 instead of num_cycles=1 to speed up this test. The lcm_period
        # property is calculated from the Movement parameters (periods), not from the
        # generated Airplanes, so we only need to generate one Airplane to test the
        # period calculation logic.
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The max_period should be 4.0 (the max of 3.0 and 4.0).
        self.assertEqual(movement.max_period, 4.0)

        # The lcm_period should be LCM(3.0, 4.0) = 12.0, not 4.0. This test will fail if
        # lcm_period only uses max_period from each WingMovement instead of collecting
        # all individual periods from WingCrossSectionMovements.
        self.assertEqual(movement.lcm_period, 12.0)

    def test_lcm_period_with_multiple_airplanes(self):
        """Test that lcm_period calculates LCM correctly with multiple periods."""
        # Create AirplaneMovements with different periods

        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_airplane_1 = ps.geometry.airplane.Airplane(
            wings=[base_wing_1],
            name="Test Airplane 1",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Make WingCrossSectionMovements for the first Airplane's Wing's root and
        # tip WingCrossSections. The root WingCrossSectionMovement must be static.
        # The tip WingCrossSectionMovement will have a period of 2.0 s.
        wcs_movements_1 = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(2.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_1 = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing_1,
            wing_cross_section_movements=wcs_movements_1,
        )

        airplane_movement_1 = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane_1,
            wing_movements=[wing_movement_1],
        )

        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_airplane_2 = ps.geometry.airplane.Airplane(
            wings=[base_wing_2],
            name="Test Airplane 2",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Make WingCrossSectionMovements for the second Airplane's Wing's root and
        # tip WingCrossSections. The root WingCrossSectionMovement must be static.
        # The tip WingCrossSectionMovement will have a period of 3.0 s.
        wcs_movements_2 = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_2 = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing_2,
            wing_cross_section_movements=wcs_movements_2,
        )

        airplane_movement_2 = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane_2,
            wing_movements=[wing_movement_2],
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use num_steps=1 instead of num_cycles=1 to speed up this test. The lcm_period
        # property is calculated from the Movement parameters (periods), not from the
        # generated Airplanes, so we only need to generate one Airplane to test the
        # period calculation logic.
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement_1, airplane_movement_2],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The LCM of 2.0 and 3.0 should be 6.0.
        self.assertEqual(movement.lcm_period, 6.0)

        # The max_period should still be 3.0.
        self.assertEqual(movement.max_period, 3.0)

    def test_delta_time_automatic_calculation(self):
        """Test that delta_time is automatically calculated when not provided."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use num_cycles=1 to speed up the test while still testing auto-calculation.
        # The auto-calculation logic doesn't depend on the specific value of num_cycles.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=1,
        )

        # Check that delta_time was calculated and is positive.
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)

    def test_num_steps_automatic_calculation_for_static(self):
        """Test that num_steps is automatically calculated for static Movement."""
        movement = self.static_movement

        # Check that num_steps was calculated and is positive.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

        # Check that it's based on num_chords.
        self.assertIsNotNone(movement.num_chords)

    def test_num_steps_automatic_calculation_for_non_static(self):
        """Test that num_steps is automatically calculated for non static Movement."""
        movement = self.basic_movement

        # Check that num_steps was calculated and is positive.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

        # Check that it's based on num_cycles.
        self.assertIsNotNone(movement.num_cycles)

    def test_explicit_num_steps_for_static(self):
        """Test that explicit num_steps works for static Movement."""
        movement = self.static_movement_with_explicit_num_steps
        self.assertEqual(movement.num_steps, 5)
        self.assertIsNone(movement.num_cycles)
        self.assertIsNone(movement.num_chords)

    def test_explicit_num_steps_for_non_static(self):
        """Test that explicit num_steps works for non static Movement."""
        movement = self.non_static_movement_with_explicit_num_steps
        self.assertEqual(movement.num_steps, 10)
        self.assertIsNone(movement.num_cycles)
        self.assertIsNone(movement.num_chords)

    def test_custom_delta_time(self):
        """Test that custom delta_time is used correctly."""
        movement = self.movement_with_custom_delta_time
        self.assertEqual(movement.delta_time, 0.05)

    def test_multiple_airplanes(self):
        """Test Movement with multiple AirplaneMovements."""
        movement = self.movement_with_multiple_airplanes

        # Check that we have multiple airplane_movements.
        self.assertEqual(len(movement.airplane_movements), 2)

        # Check that airplanes list has correct structure.
        self.assertEqual(len(movement.airplanes), 2)
        for airplane_list in movement.airplanes:
            self.assertEqual(len(airplane_list), movement.num_steps)

    def test_max_period_with_multiple_airplanes(self):
        """Test that max_period returns maximum across all AirplaneMovements."""
        movement = self.movement_with_multiple_airplanes

        # The movement has one static (period 0.0) and one with period 2.0.
        # Should return 2.0.
        self.assertEqual(movement.max_period, 2.0)

    def test_delta_time_averaging_with_multiple_airplanes(self):
        """Test that delta_time is averaged across multiple Airplanes when auto-calculated."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture(),
            airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use num_chords=1 to speed up the test while still testing auto-calculation.
        # The averaging logic doesn't depend on the specific value of num_chords.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_chords=1,
        )

        # Check that delta_time was calculated.
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)

    def test_num_steps_calculation_uses_ceil_for_static(self):
        """Test that num_steps calculation uses math.ceil for static Movement."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_chords=5,
        )

        # num_steps should be at least num_chords * c_ref / (delta_time * vCg__E).
        # Since we use math.ceil, it should be an integer.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

    def test_num_steps_calculation_uses_ceil_for_non_static(self):
        """Test that num_steps calculation uses math.ceil for non static Movement."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_cycles=3,
        )

        # num_steps should be at least num_cycles * max_period / delta_time.
        # Since we use math.ceil, it should be an integer.
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

    def test_type_4_to_5_transition_raises_error(self):
        """Test that a Wing transitioning from type 4 to type 5 symmetry raises
        an error."""
        # Create a type 4 Wing (symmetric=True, coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_4_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a WingMovement with rotation that will cause the symmetry plane
        # to become non-coincident (type 4->5 transition).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        )

        # Create an AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Attempting to create a Movement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.Movement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                num_cycles=1,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_type_3_to_2_transition_raises_error(self):
        """Test that a Wing transitioning from type 3 to type 2 symmetry raises
        an error."""
        # Create a type 3 Wing (mirror_only=True, non-coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_3_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a WingMovement with rotation that will cause the symmetry plane
        # to become coincident (type 3->2 transition).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
            ampLer_Gs_Cgs=(0.0, 0.5, 0.0),
            periodLer_Gs_Cgs=(0.0, 1.0, 0.0),
        )

        # Create an AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Attempting to create a Movement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.Movement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                delta_time=0.25,
                num_steps=5,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_type_2_to_3_transition_raises_error(self):
        """Test that a Wing transitioning from type 2 to type 3 symmetry raises
        an error."""
        # Create a type 2 Wing (mirror_only=True, coincident symmetry plane).
        base_wing = geometry_fixtures.make_type_2_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a WingMovement with rotation that will cause the symmetry plane
        # to become non-coincident (type 2->3 transition).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        )

        # Create an AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Attempting to create a Movement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.Movement(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                num_cycles=1,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_static_movement_with_symmetric_wing_succeeds(self):
        """Test that a Movement with a symmetric Wing but no motion succeeds."""
        # Create a type 4 Wing.
        base_wing = geometry_fixtures.make_type_4_wing_fixture()

        # Create an Airplane with the base Wing first, so it processes symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Now reference the Wing from the Airplane (after symmetry processing).
        processed_wing = base_airplane.wings[0]

        # Create static WingCrossSectionMovements using the actual WingCrossSections.
        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
            for wcs in processed_wing.wing_cross_sections
        ]

        # Create a static WingMovement (no rotation or translation).
        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=processed_wing,
            wing_cross_section_movements=wcs_movements,
        )

        # Create a static AirplaneMovement.
        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        # Create an OperatingPointMovement.
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Creating a Movement should succeed without raising an error.
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            num_chords=3,
        )

        # Verify the Movement was created successfully.
        self.assertIsInstance(movement, ps.movements.movement.Movement)

    def test_delta_time_invalid_string_raises_error(self):
        """Test that delta_time with invalid string raises ValueError."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        invalid_strings = [
            "invalid",
            "auto",
            "OPTIMIZE",
            "Optimize",
            "",
            "analytically_optimize",
        ]
        for invalid_string in invalid_strings:
            with self.subTest(invalid_string=invalid_string):
                with self.assertRaises(ValueError):
                    ps.movements.movement.Movement(
                        airplane_movements=airplane_movements,
                        operating_point_movement=operating_point_movement,
                        delta_time=invalid_string,
                        num_cycles=1,
                    )

    def test_delta_time_optimize_for_static_movement(self):
        """Test that delta_time='optimize' works for static Movement."""
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time="optimize",
            num_chords=3,
        )

        # Verify the Movement was created and delta_time is a positive float.
        self.assertIsInstance(movement, ps.movements.movement.Movement)
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)
        self.assertTrue(movement.static)

    def test_delta_time_optimize_calls_both_optimizers(self):
        """Test that delta_time='optimize' calls both analytical and iterative
        optimizers, using analytical result as initial guess for iterative.

        This test uses mocking to avoid running the expensive optimizations. The actual
        optimization behavior is tested in TestOptimizeDeltaTime and
        TestAnalyticallyOptimizeDeltaTime.
        """
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Mock both optimization functions.
        fake_analytical_delta_time = 0.0111111111
        fake_iterative_delta_time = 0.0123456789

        with (
            patch(
                "pterasoftware.movements.movement._analytically_optimize_delta_time"
            ) as mock_analytical,
            patch(
                "pterasoftware.movements.movement._optimize_delta_time"
            ) as mock_iterative,
        ):
            mock_analytical.return_value = fake_analytical_delta_time
            mock_iterative.return_value = fake_iterative_delta_time

            # Use num_steps=1 to speed up the test. The optimizers are mocked, so the
            # only time spent is generating airplanes after getting the delta_time.
            movement = ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time="optimize",
                num_steps=1,
            )

            # Verify the analytical optimizer was called exactly once.
            mock_analytical.assert_called_once()

            # Verify the iterative optimizer was called exactly once with the
            # analytical result as initial_delta_time.
            mock_iterative.assert_called_once()
            call_kwargs = mock_iterative.call_args[1]
            self.assertEqual(
                call_kwargs["initial_delta_time"], fake_analytical_delta_time
            )

            # Verify the Movement used the iterative optimizer's return value.
            self.assertEqual(movement.delta_time, fake_iterative_delta_time)

    def test_min_period_for_static_movement(self):
        """Test that min_period returns 0.0 for static Movement."""
        movement = self.static_movement
        self.assertEqual(movement.min_period, 0.0)

    def test_min_period_for_single_period_movement(self):
        """Test that min_period returns correct value when all periods are the same."""
        movement = self.basic_movement
        # The basic_movement has period of 2.0 for all motion.
        self.assertEqual(movement.min_period, 2.0)

    def test_min_period_with_multiple_periods(self):
        """Test that min_period returns the smallest non zero period."""
        # Reuse the multi airplane fixture which has periods 2.0 and 0.0.
        # The min non zero period should be 2.0.
        movement = self.movement_with_multiple_airplanes
        self.assertEqual(movement.min_period, 2.0)

    def test_min_period_with_different_periods(self):
        """Test that min_period returns the smallest period across Wings with
        different periods."""
        # Create an Airplane with two Wings having periods 3.0 and 4.0.
        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing_1, base_wing_2],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        wing_cross_section_movements_wing_1 = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_1 = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing_1,
            wing_cross_section_movements=wing_cross_section_movements_wing_1,
        )

        wing_cross_section_movements_wing_2 = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(4.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_2 = ps.movements.wing_movement.WingMovement(
            base_wing=base_wing_2,
            wing_cross_section_movements=wing_cross_section_movements_wing_2,
        )

        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement_1, wing_movement_2],
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The min non zero period should be 3.0 (the smaller of 3.0 and 4.0).
        self.assertEqual(movement.min_period, 3.0)

    def test_delta_time_none_calls_analytical_optimizer(self):
        """Test that delta_time=None (default) calls the analytical optimizer.

        This test uses mocking to avoid running the actual optimization. The actual
        optimization behavior is tested in TestAnalyticallyOptimizeDeltaTime.
        """
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Mock _analytically_optimize_delta_time to return a known value instantly.
        fake_optimized_delta_time = 0.0987654321

        with patch(
            "pterasoftware.movements.movement._analytically_optimize_delta_time"
        ) as mock_optimize:
            mock_optimize.return_value = fake_optimized_delta_time

            # Use num_steps=1 to speed up the test. The optimizer is mocked, so the
            # only time spent is generating Airplanes after getting the delta_time.
            movement = ps.movements.movement.Movement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=None,
                num_steps=1,
            )

            # Verify the optimizer was called exactly once.
            mock_optimize.assert_called_once()

            # Verify the Movement used the optimizer's return value.
            self.assertEqual(movement.delta_time, fake_optimized_delta_time)


class TestAnalyticallyOptimizeDeltaTime(unittest.TestCase):
    """This is a class with functions to test the _analytically_optimize_delta_time
    function."""

    def test_returns_positive_float(self):
        """Test that _analytically_optimize_delta_time returns a positive float."""
        from pterasoftware.movements.movement import (
            _analytically_optimize_delta_time,
        )

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _analytically_optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

    def test_returns_initial_for_static_movement(self):
        """Test that _analytically_optimize_delta_time returns initial_delta_time
        for static Movement."""
        from pterasoftware.movements.movement import (
            _analytically_optimize_delta_time,
        )

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _analytically_optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        # For static movement, should return the initial estimate.
        self.assertEqual(optimized_delta_time, initial_delta_time)

    def test_result_is_reasonable(self):
        """Test that _analytically_optimize_delta_time produces a reasonable result.

        The result should be within a reasonable range of the initial estimate (within
        two orders of magnitude).
        """
        from pterasoftware.movements.movement import (
            _analytically_optimize_delta_time,
        )

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _analytically_optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        # The result should be within two orders of magnitude of the initial estimate.
        self.assertGreater(optimized_delta_time, initial_delta_time / 100.0)
        self.assertLess(optimized_delta_time, initial_delta_time * 100.0)


class TestComputeWakeAreaMismatch(unittest.TestCase):
    """This is a class with functions to test the _compute_wake_area_mismatch
    function."""

    def test_returns_non_negative_value(self):
        """Test that _compute_wake_area_mismatch returns a non-negative value."""
        from pterasoftware.movements.movement import _compute_wake_area_mismatch

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Calculate the initial delta_time estimate.
        c_ref = airplane_movements[0].base_airplane.c_ref
        assert c_ref is not None
        delta_time = (
            c_ref
            / airplane_movements[0].base_airplane.wings[0].num_chordwise_panels
            / operating_point_movement.base_operating_point.vCg__E
        )

        mismatch = _compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
        )

        self.assertIsInstance(mismatch, float)
        self.assertGreaterEqual(mismatch, 0.0)

    def test_returns_zero_for_static_single_step(self):
        """Test that _compute_wake_area_mismatch returns 0.0 when no comparisons
        are made."""
        from pterasoftware.movements.movement import _compute_wake_area_mismatch

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use a delta_time that results in only 1 step for static movement.
        # With max_period = 0, num_steps will be 1, so step > 0 never runs.
        delta_time = 0.01

        mismatch = _compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
        )

        # With only 1 step, no comparisons are made, so mismatch should be 0.0.
        self.assertEqual(mismatch, 0.0)

    def test_does_not_mutate_original_movements(self):
        """Test that _compute_wake_area_mismatch does not mutate original objects."""
        from pterasoftware.movements.movement import _compute_wake_area_mismatch

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Store reference to original base_airplane.
        original_base_airplane = airplane_movements[0].base_airplane

        delta_time = 0.01

        _compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
        )

        # Verify the original base_airplane reference is unchanged.
        self.assertIs(
            airplane_movements[0].base_airplane,
            original_base_airplane,
        )


class TestOptimizeDeltaTimeStatic(unittest.TestCase):
    """This is a class with functions to test the _optimize_delta_time_static function."""

    def test_returns_positive_float(self):
        """Test that _optimize_delta_time_static returns a positive float."""
        from pterasoftware.movements.movement import _optimize_delta_time_static

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _optimize_delta_time_static(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
            mismatch_cutoff=0.01,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

    def test_early_termination_with_acceptable_initial(self):
        """Test that _optimize_delta_time_static terminates early if initial mismatch
        is below cutoff."""
        from pterasoftware.movements.movement import _optimize_delta_time_static

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        # Use a high cutoff that the initial should satisfy.
        optimized_delta_time = _optimize_delta_time_static(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
            mismatch_cutoff=1.0,
        )

        # Should return the initial since it should be below the high cutoff.
        self.assertEqual(optimized_delta_time, initial_delta_time)


class TestOptimizeDeltaTimeNonStatic(unittest.TestCase):
    """This is a class with functions to test the _optimize_delta_time_non_static
    function."""

    def test_returns_positive_float(self):
        """Test that _optimize_delta_time_non_static returns a positive float."""
        from pterasoftware.movements.movement import (
            _lcm_multiple,
            _optimize_delta_time_non_static,
        )

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Calculate the LCM period.
        all_periods = []
        for airplane_movement in airplane_movements:
            all_periods.extend(airplane_movement.all_periods)
        lcm_period = _lcm_multiple(all_periods)

        # Use a larger initial_delta_time to reduce the brute force search range.
        initial_delta_time = 0.1

        optimized_delta_time = _optimize_delta_time_non_static(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
            lcm_period=lcm_period,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

    def test_result_divides_lcm_period_evenly(self):
        """Test that _optimize_delta_time_non_static result divides LCM period evenly."""
        from pterasoftware.movements.movement import (
            _lcm_multiple,
            _optimize_delta_time_non_static,
        )

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Calculate the LCM period.
        all_periods = []
        for airplane_movement in airplane_movements:
            all_periods.extend(airplane_movement.all_periods)
        lcm_period = _lcm_multiple(all_periods)

        initial_delta_time = 0.1

        optimized_delta_time = _optimize_delta_time_non_static(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
            lcm_period=lcm_period,
        )

        # The result should divide the LCM period evenly into an integer number of
        # steps.
        num_steps = lcm_period / optimized_delta_time
        self.assertAlmostEqual(num_steps, round(num_steps), places=10)


class TestOptimizeDeltaTime(unittest.TestCase):
    """This is a class with functions to test the _optimize_delta_time function."""

    def test_returns_positive_float_within_bounds(self):
        """Test that _optimize_delta_time returns a positive float within expected
        bounds."""
        from pterasoftware.movements.movement import _optimize_delta_time

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use a larger initial_delta_time to reduce the brute force search range.
        # The fixture has a period of 2.0 s, so 0.1 s gives ~20 steps per period,
        # and the search range is ~10 to ~40 steps (~30 evaluations).
        initial_delta_time = 0.1

        # For non static movements, brute force search is used (mismatch_cutoff is
        # ignored).
        optimized_delta_time = _optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

        # The optimization searches num_steps from 0.5x to 2x the initial estimate.
        # Calculate the exact bounds using the same equations as the function.
        lcm_period = 2.0  # Fixture has a period of 2.0 s.
        initial_num_steps = lcm_period / initial_delta_time
        min_num_steps = max(1, int(initial_num_steps / 2))
        max_num_steps = int(initial_num_steps * 2) + 1
        min_delta_time = lcm_period / max_num_steps
        max_delta_time = lcm_period / min_num_steps
        self.assertGreaterEqual(optimized_delta_time, min_delta_time)
        self.assertLessEqual(optimized_delta_time, max_delta_time)

    def test_works_with_static_movement(self):
        """Test that _optimize_delta_time works with static AirplaneMovement."""
        from pterasoftware.movements.movement import _optimize_delta_time

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

    def test_dispatches_to_static_for_static_movement(self):
        """Test that _optimize_delta_time dispatches to _optimize_delta_time_static
        for static movements."""
        from pterasoftware.movements.movement import _optimize_delta_time

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        # Mock _optimize_delta_time_static to verify it's called.
        with patch(
            "pterasoftware.movements.movement._optimize_delta_time_static"
        ) as mock_static:
            mock_static.return_value = 0.012
            optimized_delta_time = _optimize_delta_time(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
            )

            # Verify _optimize_delta_time_static was called.
            mock_static.assert_called_once()
            self.assertEqual(optimized_delta_time, 0.012)

    def test_dispatches_to_non_static_for_non_static_movement(self):
        """Test that _optimize_delta_time dispatches to _optimize_delta_time_non_static
        for non static movements."""
        from pterasoftware.movements.movement import _optimize_delta_time

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.1

        # Mock _optimize_delta_time_non_static to verify it's called.
        with patch(
            "pterasoftware.movements.movement._optimize_delta_time_non_static"
        ) as mock_non_static:
            mock_non_static.return_value = 0.05
            optimized_delta_time = _optimize_delta_time(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
            )

            # Verify _optimize_delta_time_non_static was called.
            mock_non_static.assert_called_once()
            self.assertEqual(optimized_delta_time, 0.05)


class TestMovementGeneratedAttributes(unittest.TestCase):
    """Tests for Movement's generated airplanes and operating_points attributes."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all tests."""
        cls.static_movement = movement_fixtures.make_static_movement_fixture()
        cls.basic_movement = movement_fixtures.make_basic_movement_fixture()
        cls.movement_with_multiple_airplanes = (
            movement_fixtures.make_movement_with_multiple_airplanes_fixture()
        )

    def test_airplanes_correct_shape(self):
        """Test that airplanes has the correct shape (num_airplane_movements x num_steps)."""
        movement = self.basic_movement

        # Check outer dimension matches number of AirplaneMovements.
        self.assertEqual(len(movement.airplanes), len(movement.airplane_movements))

        # Check each inner dimension matches num_steps.
        for airplane_list in movement.airplanes:
            self.assertEqual(len(airplane_list), movement.num_steps)

    def test_airplanes_correct_types(self):
        """Test that airplanes contains Airplane objects."""
        movement = self.basic_movement

        for airplane_list in movement.airplanes:
            for airplane in airplane_list:
                self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_operating_points_correct_length(self):
        """Test that operating_points has correct length (num_steps)."""
        movement = self.basic_movement
        self.assertEqual(len(movement.operating_points), movement.num_steps)

    def test_operating_points_correct_types(self):
        """Test that operating_points contains OperatingPoint objects."""
        movement = self.basic_movement

        for operating_point in movement.operating_points:
            self.assertIsInstance(operating_point, ps.operating_point.OperatingPoint)

    def test_airplanes_with_multiple_airplane_movements(self):
        """Test that airplanes correctly handles multiple AirplaneMovements."""
        movement = self.movement_with_multiple_airplanes

        # Should have 2 airplane lists (one per AirplaneMovement).
        self.assertEqual(len(movement.airplanes), 2)

        # Each list should have num_steps Airplanes.
        for airplane_list in movement.airplanes:
            self.assertEqual(len(airplane_list), movement.num_steps)

    def test_airplanes_static_movement_are_consistent(self):
        """Test that static Movement generates consistent Airplanes across time steps."""
        movement = self.static_movement

        # For static movement, all Airplanes should have the same Wing positions.
        first_airplane = movement.airplanes[0][0]
        for airplane in movement.airplanes[0][1:]:
            # Check that Wing positions are the same.
            for wing_id, wing in enumerate(airplane.wings):
                first_wing = first_airplane.wings[wing_id]
                # Compare Ler_Gs_Cgs (Wing positions relative to Airplane CG).
                import numpy.testing as npt

                npt.assert_array_almost_equal(wing.Ler_Gs_Cgs, first_wing.Ler_Gs_Cgs)


class TestMovementImmutability(unittest.TestCase):
    """Tests for Movement attribute immutability."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all immutability tests."""
        cls.static_movement = movement_fixtures.make_static_movement_fixture()
        cls.basic_movement = movement_fixtures.make_basic_movement_fixture()

    def test_immutable_airplane_movements_property(self):
        """Test that airplane_movements property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_movement.airplane_movements = []

    def test_immutable_airplane_movements_tuple(self):
        """Test that airplane_movements returns a tuple (immutable sequence)."""
        airplane_movements = self.basic_movement.airplane_movements
        self.assertIsInstance(airplane_movements, tuple)

    def test_immutable_operating_point_movement_property(self):
        """Test that operating_point_movement property is read only."""
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )
        with self.assertRaises(AttributeError):
            self.basic_movement.operating_point_movement = operating_point_movement

    def test_immutable_delta_time_property(self):
        """Test that delta_time property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_movement.delta_time = 0.05

    def test_immutable_num_cycles_property(self):
        """Test that num_cycles property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_movement.num_cycles = 5

    def test_immutable_num_chords_property(self):
        """Test that num_chords property is read only."""
        with self.assertRaises(AttributeError):
            self.static_movement.num_chords = 10

    def test_immutable_num_steps_property(self):
        """Test that num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_movement.num_steps = 100

    def test_immutable_airplanes_property(self):
        """Test that airplanes property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_movement.airplanes = ()

    def test_immutable_airplanes_tuple_of_tuples(self):
        """Test that airplanes returns a tuple of tuples (immutable structure)."""
        airplanes = self.basic_movement.airplanes
        self.assertIsInstance(airplanes, tuple)
        for airplane_list in airplanes:
            self.assertIsInstance(airplane_list, tuple)

    def test_immutable_operating_points_property(self):
        """Test that operating_points property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_movement.operating_points = ()

    def test_immutable_operating_points_tuple(self):
        """Test that operating_points returns a tuple (immutable sequence)."""
        operating_points = self.basic_movement.operating_points
        self.assertIsInstance(operating_points, tuple)


class TestMovementCaching(unittest.TestCase):
    """Tests for Movement caching behavior."""

    def test_lcm_period_cache_is_populated_after_access(self):
        """Test that _lcm_period cache is populated after first access."""
        # Create a fresh Movement to test cache population.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=1,
        )

        # Cache should be None initially.
        self.assertIsNone(movement._lcm_period)

        # Access the property.
        _ = movement.lcm_period

        # Cache should now be populated.
        self.assertIsNotNone(movement._lcm_period)

    def test_max_period_cache_is_populated_after_init(self):
        """Test that _max_period cache is populated during __init__ because the static
        property is accessed during __init__ to determine num_steps calculation, and
        static depends on max_period.
        """
        # Create a fresh Movement.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=1,
        )

        # max_period is accessed during __init__ via the static property
        # (which checks self.max_period == 0), so the cache should already be populated.
        self.assertIsNotNone(movement._max_period)

    def test_static_cache_is_populated_after_init(self):
        """Test that _static cache is populated during __init__ because it is accessed
        to determine num_steps calculation.
        """
        # Create a fresh Movement.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=1,
        )

        # static is accessed during __init__ to determine num_steps calculation,
        # so the cache should already be populated.
        self.assertIsNotNone(movement._static)

    def test_min_period_cache_is_populated_after_access(self):
        """Test that _min_period cache is populated after first access."""
        # Create a fresh Movement to test cache population.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=1,
        )

        # Cache should be None initially.
        self.assertIsNone(movement._min_period)

        # Access the property.
        _ = movement.min_period

        # Cache should now be populated.
        self.assertIsNotNone(movement._min_period)


class TestMovementDeepcopy(unittest.TestCase):
    """Tests for Movement deepcopy behavior."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all deepcopy tests."""
        cls.static_movement = movement_fixtures.make_static_movement_fixture()
        cls.basic_movement = movement_fixtures.make_basic_movement_fixture()

    def test_deepcopy_returns_new_instance(self):
        """Test that deepcopy returns a new Movement instance."""
        import copy

        original = self.basic_movement
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.movements.movement.Movement)
        self.assertIsNot(original, copied)

    def test_deepcopy_preserves_attribute_values(self):
        """Test that deepcopy preserves all attribute values."""
        import copy

        original = self.basic_movement
        copied = copy.deepcopy(original)

        # Check scalar attributes.
        self.assertEqual(copied.delta_time, original.delta_time)
        self.assertEqual(copied.num_cycles, original.num_cycles)
        self.assertEqual(copied.num_chords, original.num_chords)
        self.assertEqual(copied.num_steps, original.num_steps)

        # Check derived properties.
        self.assertEqual(copied.static, original.static)
        self.assertEqual(copied.max_period, original.max_period)
        self.assertEqual(copied.min_period, original.min_period)
        self.assertEqual(copied.lcm_period, original.lcm_period)

    def test_deepcopy_airplane_movements_are_independent(self):
        """Test that deepcopied airplane_movements are independent objects."""
        import copy

        original = self.basic_movement
        copied = copy.deepcopy(original)

        # Verify airplane_movements tuples are different objects.
        self.assertIsNot(copied.airplane_movements, original.airplane_movements)

        # Verify each AirplaneMovement is a different object.
        for orig_am, copied_am in zip(
            original.airplane_movements, copied.airplane_movements
        ):
            self.assertIsNot(copied_am, orig_am)

    def test_deepcopy_operating_point_movement_is_independent(self):
        """Test that deepcopied operating_point_movement is an independent object."""
        import copy

        original = self.basic_movement
        copied = copy.deepcopy(original)

        # Verify operating_point_movement is a different object.
        self.assertIsNot(
            copied.operating_point_movement, original.operating_point_movement
        )

    def test_deepcopy_airplanes_are_independent(self):
        """Test that deepcopied airplanes are independent objects."""
        import copy

        original = self.basic_movement
        copied = copy.deepcopy(original)

        # Verify airplanes tuples are different objects.
        self.assertIsNot(copied.airplanes, original.airplanes)

        # Verify each airplane list is a different tuple.
        for orig_airplane_list, copied_airplane_list in zip(
            original.airplanes, copied.airplanes
        ):
            self.assertIsNot(copied_airplane_list, orig_airplane_list)

    def test_deepcopy_operating_points_are_independent(self):
        """Test that deepcopied operating_points are independent objects."""
        import copy

        original = self.basic_movement
        copied = copy.deepcopy(original)

        # Verify operating_points tuple is a different object.
        self.assertIsNot(copied.operating_points, original.operating_points)

        # Verify each OperatingPoint is a different object.
        for orig_op, copied_op in zip(
            original.operating_points, copied.operating_points
        ):
            self.assertIsNot(copied_op, orig_op)

    def test_deepcopy_static_movement(self):
        """Test that deepcopy works correctly for static Movement."""
        import copy

        original = self.static_movement
        copied = copy.deepcopy(original)

        self.assertIsInstance(copied, ps.movements.movement.Movement)
        self.assertIsNot(original, copied)
        self.assertTrue(copied.static)
        self.assertEqual(copied.max_period, 0.0)


class TestLcmFunctions(unittest.TestCase):
    """Tests for _lcm and _lcm_multiple module level functions."""

    def test_lcm_of_two_positive_numbers(self):
        """Test _lcm returns correct LCM for two positive numbers."""
        result = ps.movements.movement._lcm(2.0, 3.0)
        self.assertEqual(result, 6.0)

    def test_lcm_of_same_numbers(self):
        """Test _lcm returns the number when both inputs are the same."""
        result = ps.movements.movement._lcm(4.0, 4.0)
        self.assertEqual(result, 4.0)

    def test_lcm_with_first_zero(self):
        """Test _lcm returns 0.0 when first input is zero."""
        result = ps.movements.movement._lcm(0.0, 5.0)
        self.assertEqual(result, 0.0)

    def test_lcm_with_second_zero(self):
        """Test _lcm returns 0.0 when second input is zero."""
        result = ps.movements.movement._lcm(5.0, 0.0)
        self.assertEqual(result, 0.0)

    def test_lcm_with_both_zero(self):
        """Test _lcm returns 0.0 when both inputs are zero."""
        result = ps.movements.movement._lcm(0.0, 0.0)
        self.assertEqual(result, 0.0)

    def test_lcm_of_one_and_any_number(self):
        """Test _lcm of 1.0 and any number returns that number."""
        result = ps.movements.movement._lcm(1.0, 7.0)
        self.assertEqual(result, 7.0)

        result = ps.movements.movement._lcm(7.0, 1.0)
        self.assertEqual(result, 7.0)

    def test_lcm_of_multiples(self):
        """Test _lcm of a number and its multiple returns the larger number."""
        result = ps.movements.movement._lcm(3.0, 9.0)
        self.assertEqual(result, 9.0)

        result = ps.movements.movement._lcm(9.0, 3.0)
        self.assertEqual(result, 9.0)

    def test_lcm_multiple_empty_list(self):
        """Test _lcm_multiple returns 0.0 for empty list."""
        result = ps.movements.movement._lcm_multiple([])
        self.assertEqual(result, 0.0)

    def test_lcm_multiple_all_zeros(self):
        """Test _lcm_multiple returns 0.0 when all periods are zero."""
        result = ps.movements.movement._lcm_multiple([0.0, 0.0, 0.0])
        self.assertEqual(result, 0.0)

    def test_lcm_multiple_single_nonzero(self):
        """Test _lcm_multiple returns the value for a single non zero period."""
        result = ps.movements.movement._lcm_multiple([5.0])
        self.assertEqual(result, 5.0)

    def test_lcm_multiple_single_zero(self):
        """Test _lcm_multiple returns 0.0 for a single zero period."""
        result = ps.movements.movement._lcm_multiple([0.0])
        self.assertEqual(result, 0.0)

    def test_lcm_multiple_mixed_with_zeros(self):
        """Test _lcm_multiple correctly ignores zeros in the list."""
        result = ps.movements.movement._lcm_multiple([0.0, 2.0, 0.0, 3.0, 0.0])
        self.assertEqual(result, 6.0)

    def test_lcm_multiple_three_periods(self):
        """Test _lcm_multiple returns correct LCM for three periods."""
        result = ps.movements.movement._lcm_multiple([2.0, 3.0, 4.0])
        self.assertEqual(result, 12.0)

    def test_lcm_multiple_many_same_periods(self):
        """Test _lcm_multiple returns the period when all are the same."""
        result = ps.movements.movement._lcm_multiple([5.0, 5.0, 5.0, 5.0])
        self.assertEqual(result, 5.0)

    def test_lcm_multiple_coprime_periods(self):
        """Test _lcm_multiple of coprime numbers returns their product."""
        # 2, 3, and 5 are coprime, so LCM = 2 * 3 * 5 = 30.
        result = ps.movements.movement._lcm_multiple([2.0, 3.0, 5.0])
        self.assertEqual(result, 30.0)

    def test_lcm_non_integer_periods(self):
        """Test _lcm returns correct LCM for non integer periods."""
        # LCM(1.5, 2.5) = 7.5 (both divide 7.5 evenly: 7.5/1.5=5, 7.5/2.5=3).
        result = ps.movements.movement._lcm(1.5, 2.5)
        self.assertAlmostEqual(result, 7.5, places=6)

    def test_lcm_multiple_non_integer_periods(self):
        """Test _lcm_multiple returns correct LCM for non integer periods."""
        # LCM(1.5, 2.0, 2.5) = 30.0.
        # 30.0/1.5=20, 30.0/2.0=15, 30.0/2.5=12.
        result = ps.movements.movement._lcm_multiple([1.5, 2.0, 2.5])
        self.assertAlmostEqual(result, 30.0, places=6)

    def test_lcm_small_periods(self):
        """Test _lcm handles small periods correctly without precision issues."""
        # LCM(0.001, 0.002) = 0.002.
        result = ps.movements.movement._lcm(0.001, 0.002)
        self.assertAlmostEqual(result, 0.002, places=9)

    def test_lcm_multiple_small_periods(self):
        """Test _lcm_multiple handles small periods correctly."""
        # LCM(0.01, 0.02, 0.03) = 0.06.
        result = ps.movements.movement._lcm_multiple([0.01, 0.02, 0.03])
        self.assertAlmostEqual(result, 0.06, places=9)


class TestAnalyticallyOptimizeDeltaTimeEdgeCases(unittest.TestCase):
    """Tests for edge cases in the _analytically_optimize_delta_time function."""

    def test_with_multiple_airplanes(self):
        """Test _analytically_optimize_delta_time works with multiple Airplanes."""
        from pterasoftware.movements.movement import (
            _analytically_optimize_delta_time,
        )

        # Create two AirplaneMovements with different motion.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
            airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _analytically_optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

    def test_with_multiple_wings_per_airplane(self):
        """Test _analytically_optimize_delta_time works with multiple Wings per
        Airplane."""
        from pterasoftware.movements.movement import (
            _analytically_optimize_delta_time,
        )

        # Create a multi wing Airplane with Cg_GP1_CgP1 at origin.
        # Use simple tapered wings.
        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing_1, base_wing_2],
            name="Multi Wing Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Create WingMovements for each Wing in the processed Airplane.
        wing_movements = []
        for wing in base_airplane.wings:
            wcs_movements = [
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=wcs,
                    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                )
                for wcs in wing.wing_cross_sections
            ]
            # Set non zero amplitude on tip for motion.
            wcs_movements[-1] = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=wing.wing_cross_sections[-1],
                    periodLp_Wcsp_Lpp=(2.0, 0.0, 0.0),
                    ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
                )
            )
            wing_movements.append(
                ps.movements.wing_movement.WingMovement(
                    base_wing=wing,
                    wing_cross_section_movements=wcs_movements,
                )
            )

        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        optimized_delta_time = _analytically_optimize_delta_time(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)

    def test_coarse_temporal_resolution_warning(self):
        """Test that _analytically_optimize_delta_time warns when steps per min period
        is less than 20."""
        import logging

        from pterasoftware.movements.movement import (
            _analytically_optimize_delta_time,
        )

        # Create an AirplaneMovement with very fast motion (short period) and few
        # chordwise panels. This should result in fewer than 20 steps per min period.
        # The basic fixture has period 2.0 s. With few chordwise panels and fast motion,
        # the trailing edge panels are large relative to the motion displacement.
        # Create a Wing with only 2 chordwise panels to make trailing edge panels large.
        base_airplane = geometry_fixtures.make_2_chordwise_panels_airplane_fixture()
        wing = base_airplane.wings[0]

        wcs_movements = [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            )
            for wcs in wing.wing_cross_sections
        ]
        # Add fast motion with short period to tip.
        wcs_movements[-1] = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wing.wing_cross_sections[-1],
                periodLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.01, 0.0, 0.0),
            )
        )

        wing_movement = ps.movements.wing_movement.WingMovement(
            base_wing=wing,
            wing_cross_section_movements=wcs_movements,
        )

        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01

        # Capture warnings.
        with self.assertLogs(
            "pterasoftware.movements.movement", level=logging.WARNING
        ) as log_context:
            _analytically_optimize_delta_time(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
            )

        # Verify the warning was issued.
        warning_found = any(
            "time steps per minimum period" in msg for msg in log_context.output
        )
        self.assertTrue(
            warning_found,
            "Expected warning about time steps per minimum period not found.",
        )


class TestComputeWakeAreaMismatchEdgeCases(unittest.TestCase):
    """Tests for edge cases in the _compute_wake_area_mismatch function."""

    def test_with_multiple_airplanes(self):
        """Test _compute_wake_area_mismatch works with multiple Airplanes."""
        from pterasoftware.movements.movement import _compute_wake_area_mismatch

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture(),
            airplane_movement_fixtures.make_static_airplane_movement_fixture(),
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Calculate a reasonable delta_time.
        c_ref = airplane_movements[0].base_airplane.c_ref
        assert c_ref is not None
        delta_time = (
            c_ref
            / airplane_movements[0].base_airplane.wings[0].num_chordwise_panels
            / operating_point_movement.base_operating_point.vCg__E
        )

        mismatch = _compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
        )

        self.assertIsInstance(mismatch, float)
        self.assertGreaterEqual(mismatch, 0.0)

    def test_with_multiple_wings_per_airplane(self):
        """Test _compute_wake_area_mismatch works with multiple Wings per Airplane."""
        from pterasoftware.movements.movement import _compute_wake_area_mismatch

        # Create a multi wing Airplane with Cg_GP1_CgP1 at origin.
        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing_1, base_wing_2],
            name="Multi Wing Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Create WingMovements for each Wing.
        wing_movements = []
        for wing in base_airplane.wings:
            wcs_movements = [
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=wcs,
                    periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                    ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                )
                for wcs in wing.wing_cross_sections
            ]
            wcs_movements[-1] = (
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=wing.wing_cross_sections[-1],
                    periodLp_Wcsp_Lpp=(2.0, 0.0, 0.0),
                    ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
                )
            )
            wing_movements.append(
                ps.movements.wing_movement.WingMovement(
                    base_wing=wing,
                    wing_cross_section_movements=wcs_movements,
                )
            )

        airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
        )

        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        c_ref = base_airplane.c_ref
        assert c_ref is not None
        delta_time = (
            c_ref
            / base_airplane.wings[0].num_chordwise_panels
            / operating_point_movement.base_operating_point.vCg__E
        )

        mismatch = _compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
        )

        self.assertIsInstance(mismatch, float)
        self.assertGreaterEqual(mismatch, 0.0)

    def test_with_non_static_movement_multiple_steps(self):
        """Test _compute_wake_area_mismatch computes correctly over multiple time
        steps for non static movement."""
        from pterasoftware.movements.movement import _compute_wake_area_mismatch

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # Use a delta_time that results in multiple time steps over the max period.
        # The basic fixture has period 2.0 s. With delta_time=0.5, we get 4 steps.
        delta_time = 0.5

        mismatch = _compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
        )

        # The mismatch should be computed (not zero) since we have multiple steps.
        self.assertIsInstance(mismatch, float)
        self.assertGreaterEqual(mismatch, 0.0)


class TestOptimizeDeltaTimeNonStaticWarnings(unittest.TestCase):
    """Tests for warnings in _optimize_delta_time_non_static."""

    def test_warns_when_at_lower_bound(self):
        """Test that _optimize_delta_time_non_static warns when optimum is at lower
        bound.

        This test uses mocking to avoid running the expensive optimization. We mock
        _compute_wake_area_mismatch to return decreasing values as num_steps increases,
        forcing the best value to be at min_num_steps (lower bound).
        """
        import logging

        from pterasoftware.movements.movement import _optimize_delta_time_non_static

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        lcm_period = 2.0  # Fixture has period 2.0.
        initial_delta_time = 0.1  # Results in ~20 steps, search range ~10 to ~41.

        # Mock _compute_wake_area_mismatch to return values that decrease with
        # increasing delta_time (fewer steps), so the best is at min_num_steps
        # (lower bound = largest delta_time).
        def mock_mismatch(dt, am, opm):
            # Mismatch decreases as delta_time increases (fewer steps).
            return 1.0 / dt

        with (
            patch(
                "pterasoftware.movements.movement._compute_wake_area_mismatch",
                side_effect=mock_mismatch,
            ),
            self.assertLogs(
                "pterasoftware.movements.movement", level=logging.WARNING
            ) as log_context,
        ):
            _optimize_delta_time_non_static(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
                lcm_period=lcm_period,
            )

        # Check that lower bound warning was issued.
        lower_bound_warning_found = any(
            "lower bound" in msg.lower() for msg in log_context.output
        )
        self.assertTrue(
            lower_bound_warning_found,
            "Expected warning about lower bound not found.",
        )

    def test_warns_when_at_upper_bound(self):
        """Test that _optimize_delta_time_non_static warns when optimum is at upper
        bound.

        This test uses mocking to avoid running the expensive optimization. We mock
        _compute_wake_area_mismatch to return increasing values as num_steps increases,
        forcing the best value to be at max_num_steps (upper bound).
        """
        import logging

        from pterasoftware.movements.movement import _optimize_delta_time_non_static

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        lcm_period = 2.0  # Fixture has period 2.0.
        initial_delta_time = 0.1  # Results in ~20 steps, search range ~10 to ~41.

        # Mock _compute_wake_area_mismatch to return values that decrease with
        # decreasing delta_time (more steps), so the best is at max_num_steps
        # (upper bound = smallest delta_time).
        def mock_mismatch(dt, am, opm):
            # Mismatch decreases as delta_time decreases (more steps).
            return dt * 10.0

        with (
            patch(
                "pterasoftware.movements.movement._compute_wake_area_mismatch",
                side_effect=mock_mismatch,
            ),
            self.assertLogs(
                "pterasoftware.movements.movement", level=logging.WARNING
            ) as log_context,
        ):
            _optimize_delta_time_non_static(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
                lcm_period=lcm_period,
            )

        # Check that upper bound warning was issued.
        upper_bound_warning_found = any(
            "upper bound" in msg.lower() for msg in log_context.output
        )
        self.assertTrue(
            upper_bound_warning_found,
            "Expected warning about upper bound not found.",
        )


class TestOptimizeDeltaTimeStaticWarnings(unittest.TestCase):
    """Tests for warning logic in _optimize_delta_time_static.

    Note: scipy's bounded optimizer uses xatol=0.001 tolerance, so it may not converge
    exactly to bounds. These tests verify the warning logic by mocking the optimizer
    to return values exactly at the bounds.
    """

    def test_warning_logic_for_lower_bound(self):
        """Test that _optimize_delta_time_static warning logic triggers correctly
        for lower bound.

        This test directly verifies the warning logic by patching the optimizer to
        return a value at the lower bound.
        """
        import logging

        import scipy.optimize as sp_opt

        from pterasoftware.movements.movement import _optimize_delta_time_static

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01
        lower_bound = initial_delta_time / 2.0

        # Mock the optimizer to return a value exactly at the lower bound.
        class MockResult:
            success = True
            x = lower_bound

        with (
            patch(
                "pterasoftware.movements.movement._compute_wake_area_mismatch",
                return_value=0.5,  # Any value above cutoff.
            ),
            patch.object(
                sp_opt,
                "minimize_scalar",
                return_value=MockResult(),
            ),
            self.assertLogs(
                "pterasoftware.movements.movement", level=logging.WARNING
            ) as log_context,
        ):
            _optimize_delta_time_static(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
                mismatch_cutoff=0.0,
            )

        # Check that lower bound warning was issued.
        lower_bound_warning_found = any(
            "lower bound" in msg.lower() for msg in log_context.output
        )
        self.assertTrue(
            lower_bound_warning_found,
            "Expected warning about lower bound not found.",
        )

    def test_warning_logic_for_upper_bound(self):
        """Test that _optimize_delta_time_static warning logic triggers correctly
        for upper bound.

        This test directly verifies the warning logic by patching the optimizer to
        return a value at the upper bound.
        """
        import logging

        import scipy.optimize as sp_opt

        from pterasoftware.movements.movement import _optimize_delta_time_static

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.01
        upper_bound = initial_delta_time * 2.0

        # Mock the optimizer to return a value exactly at the upper bound.
        class MockResult:
            success = True
            x = upper_bound

        with (
            patch(
                "pterasoftware.movements.movement._compute_wake_area_mismatch",
                return_value=0.5,  # Any value above cutoff.
            ),
            patch.object(
                sp_opt,
                "minimize_scalar",
                return_value=MockResult(),
            ),
            self.assertLogs(
                "pterasoftware.movements.movement", level=logging.WARNING
            ) as log_context,
        ):
            _optimize_delta_time_static(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                initial_delta_time=initial_delta_time,
                mismatch_cutoff=0.0,
            )

        # Check that upper bound warning was issued.
        upper_bound_warning_found = any(
            "upper bound" in msg.lower() for msg in log_context.output
        )
        self.assertTrue(
            upper_bound_warning_found,
            "Expected warning about upper bound not found.",
        )


class TestMovementWithOperatingPointMovementPeriod(unittest.TestCase):
    """Tests for Movement when OperatingPointMovement has non zero period."""

    def test_lcm_period_includes_operating_point_movement_period(self):
        """Test that lcm_period includes the OperatingPointMovement period."""
        # Create a static AirplaneMovement.
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]

        # Create an OperatingPointMovement with a non zero period.
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point,
                ampVCg__E=1.0,
                periodVCg__E=3.0,
            )
        )

        # Create the Movement with explicit num_steps to avoid auto calculation.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The lcm_period should be 3.0 (from OperatingPointMovement).
        # Since the AirplaneMovement is static (period 0.0), the only period is 3.0.
        self.assertEqual(movement.lcm_period, 3.0)

    def test_lcm_period_combines_airplane_and_operating_point_periods(self):
        """Test that lcm_period combines periods from AirplaneMovement and
        OperatingPointMovement."""
        # Create a non static AirplaneMovement with period 2.0.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]

        # Create an OperatingPointMovement with period 3.0.
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point,
                ampVCg__E=1.0,
                periodVCg__E=3.0,
            )
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The lcm_period should be LCM(2.0, 3.0) = 6.0.
        self.assertEqual(movement.lcm_period, 6.0)

    def test_min_period_includes_operating_point_movement_period(self):
        """Test that min_period includes the OperatingPointMovement period."""
        # Create a non static AirplaneMovement with period 2.0.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]

        # Create an OperatingPointMovement with a shorter period (1.5).
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point,
                ampVCg__E=1.0,
                periodVCg__E=1.5,
            )
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The min_period should be 1.5 (from OperatingPointMovement, smaller than 2.0).
        self.assertEqual(movement.min_period, 1.5)

    def test_max_period_includes_operating_point_movement_period(self):
        """Test that max_period includes the OperatingPointMovement period."""
        # Create a non static AirplaneMovement with period 2.0.
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]

        # Create an OperatingPointMovement with a longer period (5.0).
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = (
            ps.movements.operating_point_movement.OperatingPointMovement(
                base_operating_point=base_operating_point,
                ampVCg__E=1.0,
                periodVCg__E=5.0,
            )
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The max_period should be 5.0 (from OperatingPointMovement, larger than 2.0).
        self.assertEqual(movement.max_period, 5.0)


if __name__ == "__main__":
    unittest.main()
