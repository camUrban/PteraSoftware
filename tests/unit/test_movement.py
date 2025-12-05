"""This module contains classes to test Movements and CoupledMovements."""

import unittest

import numpy as np

import pterasoftware as ps

from tests.unit.fixtures import movement_fixtures
from tests.unit.fixtures import airplane_movement_fixtures
from tests.unit.fixtures import operating_point_fixtures
from tests.unit.fixtures import geometry_fixtures


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

    def test_initialization_valid_parameters(self):
        """Test Movement initialization with valid parameters."""
        movement = self.basic_movement
        self.assertIsInstance(movement, ps.movements.movement.Movement)
        self.assertIsInstance(movement.airplane_movements, list)
        self.assertEqual(len(movement.airplane_movements), 1)
        self.assertIsInstance(
            movement.airplane_movements[0],
            ps.movements.airplane_movement.AirplaneMovement,
        )
        self.assertIsInstance(
            movement.operating_point_movement,
            ps.movements.operating_point_movement.OperatingPointMovement,
        )
        self.assertIsInstance(movement.delta_time, float)
        self.assertGreater(movement.delta_time, 0.0)
        self.assertEqual(movement.num_cycles, 1)
        self.assertIsNone(movement.num_chords)
        self.assertIsInstance(movement.num_steps, int)
        self.assertGreater(movement.num_steps, 0)

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

        # Test with positive delta_time works.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.01,
            num_chords=10,
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
        """Test that non-static Movement with num_steps=None requires num_cycles."""
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
        """Test that non-static Movement with num_steps=None cannot have num_chords."""
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

        # Test with valid positive integer.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=5,
        )
        self.assertEqual(movement.num_cycles, 5)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "three"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
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

        # Test with valid positive integer.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_chords=15,
        )
        self.assertEqual(movement.num_chords, 15)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "ten"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
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

        # Test with valid positive integer.
        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=200,
        )
        self.assertEqual(movement.num_steps, 200)

        # Test with invalid values.
        invalid_values = [0, -5, 2.5, "hundred"]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
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
        """Test that static property returns False for non-static Movement."""
        movement = self.basic_movement
        self.assertFalse(movement.static)

    def test_max_period_for_static_movement(self):
        """Test that max_period returns 0.0 for static Movement."""
        movement = self.static_movement
        self.assertEqual(movement.max_period, 0.0)

    def test_max_period_for_non_static_movement(self):
        """Test that max_period returns correct value for non-static Movement."""
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

        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            num_cycles=1,
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

        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            num_cycles=1,
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

        # Create Movement with both AirplaneMovements (periods 2.0 and 3.0).
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement_1, airplane_movement_2],
            operating_point_movement=operating_point_movement,
            num_cycles=1,
        )

        # The LCM of 2.0 and 3.0 should be 6.0.
        self.assertEqual(movement.lcm_period, 6.0)

        # The max_period should still be 3.0.
        self.assertEqual(movement.max_period, 3.0)

    def test_airplanes_generation(self):
        """Test that airplanes are generated correctly."""
        movement = self.basic_movement

        # Check that airplanes attribute is a list of lists.
        self.assertIsInstance(movement.airplanes, list)
        self.assertEqual(len(movement.airplanes), len(movement.airplane_movements))

        # Check that each element is a list of Airplanes.
        for airplane_list in movement.airplanes:
            self.assertIsInstance(airplane_list, list)
            self.assertEqual(len(airplane_list), movement.num_steps)
            for airplane in airplane_list:
                self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_operating_points_generation(self):
        """Test that operating_points are generated correctly."""
        movement = self.basic_movement

        # Check that operating_points attribute is a list.
        self.assertIsInstance(movement.operating_points, list)
        self.assertEqual(len(movement.operating_points), movement.num_steps)

        # Check that each element is an OperatingPoint.
        for operating_point in movement.operating_points:
            self.assertIsInstance(operating_point, ps.operating_point.OperatingPoint)

    def test_delta_time_automatic_calculation(self):
        """Test that delta_time is automatically calculated when not provided."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=3,
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
        """Test that num_steps is automatically calculated for non-static Movement."""
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
        """Test that explicit num_steps works for non-static Movement."""
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

        movement = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_chords=10,
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
        """Test that num_steps calculation uses math.ceil for non-static Movement."""
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

        invalid_strings = ["invalid", "auto", "OPTIMIZE", "Optimize", ""]
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

    def test_delta_time_optimize_within_expected_bounds(self):
        """Test that optimized delta_time is within the expected search bounds."""
        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        # First, create a Movement with automatic delta_time to get the initial
        # estimate.
        movement_auto = ps.movements.movement.Movement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_cycles=1,
        )
        initial_estimate = movement_auto.delta_time

        # Now create a Movement with optimized delta_time.
        # Need fresh fixtures since the previous ones may have been modified.
        airplane_movements_opt = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement_opt = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )
        movement_optimized = ps.movements.movement.Movement(
            airplane_movements=airplane_movements_opt,
            operating_point_movement=operating_point_movement_opt,
            delta_time="optimize",
            num_cycles=1,
        )

        # Verify the optimized delta_time is within the search bounds.
        # The optimization searches within [initial / 10, initial * 2].
        self.assertGreaterEqual(movement_optimized.delta_time, initial_estimate / 10)
        self.assertLessEqual(movement_optimized.delta_time, initial_estimate * 2)


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


class TestOptimizeDeltaTime(unittest.TestCase):
    """This is a class with functions to test the _optimize_delta_time function."""

    def test_result_within_bounds(self):
        """Test that _optimize_delta_time returns a value within expected bounds."""
        from pterasoftware.movements.movement import _optimize_delta_time

        airplane_movements = [
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.03

        optimized_delta_time = _optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        # The optimization searches within [initial / 10, initial * 2].
        self.assertGreaterEqual(optimized_delta_time, initial_delta_time / 10)
        self.assertLessEqual(optimized_delta_time, initial_delta_time * 2)

    def test_works_with_static_movement(self):
        """Test that _optimize_delta_time works with static AirplaneMovement."""
        from pterasoftware.movements.movement import _optimize_delta_time

        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        initial_delta_time = 0.02

        optimized_delta_time = _optimize_delta_time(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            initial_delta_time=initial_delta_time,
        )

        self.assertIsInstance(optimized_delta_time, float)
        self.assertGreater(optimized_delta_time, 0.0)


class TestCoupledMovement(unittest.TestCase):
    """This is a class with functions to test CoupledMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all CoupledMovement tests."""
        cls.basic_coupled_movement = (
            movement_fixtures.make_basic_coupled_movement_fixture()
        )
        cls.static_coupled_movement = (
            movement_fixtures.make_static_coupled_movement_fixture()
        )
        cls.coupled_movement_with_angular_speed = (
            movement_fixtures.make_coupled_movement_with_angular_speed_fixture()
        )
        cls.coupled_movement_with_attitude_angles = (
            movement_fixtures.make_coupled_movement_with_attitude_angles_fixture()
        )
        cls.full_coupled_movement = (
            movement_fixtures.make_full_coupled_movement_fixture()
        )
        cls.coupled_movement_with_custom_delta_time = (
            movement_fixtures.make_coupled_movement_with_custom_delta_time_fixture()
        )
        cls.coupled_movement_with_more_prescribed_steps = (
            movement_fixtures.make_coupled_movement_with_more_prescribed_steps_fixture()
        )
        cls.coupled_movement_with_more_free_steps = (
            movement_fixtures.make_coupled_movement_with_more_free_steps_fixture()
        )

    def test_initialization_valid_parameters(self):
        """Test CoupledMovement initialization with valid parameters."""
        coupled_movement = self.basic_coupled_movement
        self.assertIsInstance(coupled_movement, ps.movements.movement.CoupledMovement)
        self.assertIsInstance(
            coupled_movement.airplane_movement,
            ps.movements.airplane_movement.AirplaneMovement,
        )
        self.assertIsInstance(coupled_movement.coupled_operating_points, list)
        self.assertEqual(len(coupled_movement.coupled_operating_points), 1)
        self.assertIsInstance(
            coupled_movement.coupled_operating_points[0],
            ps.operating_point.CoupledOperatingPoint,
        )
        self.assertIsInstance(coupled_movement.delta_time, float)
        self.assertGreater(coupled_movement.delta_time, 0.0)
        self.assertIsInstance(coupled_movement.prescribed_num_steps, int)
        self.assertGreater(coupled_movement.prescribed_num_steps, 0)
        self.assertIsInstance(coupled_movement.free_num_steps, int)
        self.assertGreater(coupled_movement.free_num_steps, 0)
        self.assertIsInstance(coupled_movement.num_steps, int)
        self.assertEqual(
            coupled_movement.num_steps,
            coupled_movement.prescribed_num_steps + coupled_movement.free_num_steps,
        )

    def test_airplane_movement_validation_not_airplane_movement(self):
        """Test that airplane_movement must be an AirplaneMovement."""
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.CoupledMovement(
                airplane_movement="not an airplane movement",
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.1,
                prescribed_num_steps=5,
                free_num_steps=5,
            )

    def test_initial_coupled_operating_point_validation_not_coupled_operating_point(
        self,
    ):
        """Test that initial_coupled_operating_point must be a
        CoupledOperatingPoint."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point="not a coupled operating point",
                delta_time=0.1,
                prescribed_num_steps=5,
                free_num_steps=5,
            )

    def test_initial_coupled_operating_point_validation_regular_operating_point(self):
        """Test that initial_coupled_operating_point cannot be a regular
        OperatingPoint."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        regular_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=regular_operating_point,
                delta_time=0.1,
                prescribed_num_steps=5,
                free_num_steps=5,
            )

    def test_delta_time_validation_positive(self):
        """Test that delta_time must be positive."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        # Test with positive delta_time works.
        coupled_movement = ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.01,
            prescribed_num_steps=5,
            free_num_steps=5,
        )
        self.assertEqual(coupled_movement.delta_time, 0.01)

    def test_delta_time_validation_zero(self):
        """Test that delta_time cannot be zero."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.0,
                prescribed_num_steps=5,
                free_num_steps=5,
            )

    def test_delta_time_validation_negative(self):
        """Test that delta_time cannot be negative."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=-0.01,
                prescribed_num_steps=5,
                free_num_steps=5,
            )

    def test_delta_time_converted_to_float(self):
        """Test that delta_time is converted to float if provided as int."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        coupled_movement = ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=1,
            prescribed_num_steps=5,
            free_num_steps=5,
        )
        self.assertIsInstance(coupled_movement.delta_time, float)
        self.assertEqual(coupled_movement.delta_time, 1.0)

    def test_prescribed_num_steps_validation_positive(self):
        """Test that prescribed_num_steps must be a positive integer."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        # Test with valid positive integer.
        coupled_movement = ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.1,
            prescribed_num_steps=10,
            free_num_steps=5,
        )
        self.assertEqual(coupled_movement.prescribed_num_steps, 10)

    def test_prescribed_num_steps_validation_zero(self):
        """Test that prescribed_num_steps cannot be zero."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.1,
                prescribed_num_steps=0,
                free_num_steps=5,
            )

    def test_prescribed_num_steps_validation_negative(self):
        """Test that prescribed_num_steps cannot be negative."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.1,
                prescribed_num_steps=-5,
                free_num_steps=5,
            )

    def test_prescribed_num_steps_validation_not_integer(self):
        """Test that prescribed_num_steps must be an integer."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        invalid_values = [5.5, "five", 3.0]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):
                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.CoupledMovement(
                        airplane_movement=airplane_movement,
                        initial_coupled_operating_point=initial_coupled_operating_point,
                        delta_time=0.1,
                        prescribed_num_steps=invalid_value,
                        free_num_steps=5,
                    )

    def test_free_num_steps_validation_positive(self):
        """Test that free_num_steps must be a positive integer."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        # Test with valid positive integer.
        coupled_movement = ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.1,
            prescribed_num_steps=5,
            free_num_steps=10,
        )
        self.assertEqual(coupled_movement.free_num_steps, 10)

    def test_free_num_steps_validation_zero(self):
        """Test that free_num_steps cannot be zero."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.1,
                prescribed_num_steps=5,
                free_num_steps=0,
            )

    def test_free_num_steps_validation_negative(self):
        """Test that free_num_steps cannot be negative."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        with self.assertRaises((ValueError, TypeError)):
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.1,
                prescribed_num_steps=5,
                free_num_steps=-5,
            )

    def test_free_num_steps_validation_not_integer(self):
        """Test that free_num_steps must be an integer."""
        airplane_movement = (
            airplane_movement_fixtures.make_basic_airplane_movement_fixture()
        )
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        invalid_values = [5.5, "five", 3.0]
        for invalid_value in invalid_values:
            with self.subTest(invalid_value=invalid_value):

                with self.assertRaises((ValueError, TypeError)):
                    ps.movements.movement.CoupledMovement(
                        airplane_movement=airplane_movement,
                        initial_coupled_operating_point=initial_coupled_operating_point,
                        delta_time=0.1,
                        prescribed_num_steps=5,
                        free_num_steps=invalid_value,
                    )

    def test_num_steps_calculation(self):
        """Test that num_steps is correctly calculated from prescribed and free steps."""
        coupled_movement = self.basic_coupled_movement
        expected_num_steps = (
            coupled_movement.prescribed_num_steps + coupled_movement.free_num_steps
        )
        self.assertEqual(coupled_movement.num_steps, expected_num_steps)

    def test_num_steps_with_more_prescribed_steps(self):
        """Test num_steps when prescribed steps outnumber free steps."""
        coupled_movement = self.coupled_movement_with_more_prescribed_steps
        self.assertEqual(coupled_movement.prescribed_num_steps, 15)
        self.assertEqual(coupled_movement.free_num_steps, 3)
        self.assertEqual(coupled_movement.num_steps, 18)

    def test_num_steps_with_more_free_steps(self):
        """Test num_steps when free steps outnumber prescribed steps."""
        coupled_movement = self.coupled_movement_with_more_free_steps
        self.assertEqual(coupled_movement.prescribed_num_steps, 3)
        self.assertEqual(coupled_movement.free_num_steps, 15)
        self.assertEqual(coupled_movement.num_steps, 18)

    def test_airplanes_generation(self):
        """Test that airplanes are generated correctly."""
        coupled_movement = self.basic_coupled_movement

        # Check that airplanes attribute is a list of Airplanes.
        self.assertIsInstance(coupled_movement.airplanes, list)
        self.assertEqual(len(coupled_movement.airplanes), coupled_movement.num_steps)

        # Check that each element is an Airplane.
        for airplane in coupled_movement.airplanes:
            self.assertIsInstance(airplane, ps.geometry.airplane.Airplane)

    def test_coupled_operating_points_initialization(self):
        """Test that coupled_operating_points is initialized with the initial point."""
        coupled_movement = self.basic_coupled_movement

        # Check that coupled_operating_points starts with just the initial point.
        self.assertIsInstance(coupled_movement.coupled_operating_points, list)
        self.assertEqual(len(coupled_movement.coupled_operating_points), 1)
        self.assertIsInstance(
            coupled_movement.coupled_operating_points[0],
            ps.operating_point.CoupledOperatingPoint,
        )

    def test_static_property_for_static_coupled_movement(self):
        """Test that static property returns True for static CoupledMovement."""
        coupled_movement = self.static_coupled_movement
        self.assertTrue(coupled_movement.static)

    def test_static_property_for_non_static_coupled_movement(self):
        """Test that static property returns False for non-static CoupledMovement."""
        coupled_movement = self.basic_coupled_movement
        self.assertFalse(coupled_movement.static)

    def test_max_period_for_static_coupled_movement(self):
        """Test that max_period returns 0.0 for static CoupledMovement."""
        coupled_movement = self.static_coupled_movement
        self.assertEqual(coupled_movement.max_period, 0.0)

    def test_max_period_for_non_static_coupled_movement(self):
        """Test that max_period returns correct value for non-static CoupledMovement."""
        coupled_movement = self.basic_coupled_movement
        # The basic_coupled_movement uses basic_airplane_movement with period 2.0.
        self.assertEqual(coupled_movement.max_period, 2.0)

    def test_max_period_returns_airplane_movement_max_period(self):
        """Test that max_period returns the AirplaneMovement's max_period."""
        coupled_movement = self.full_coupled_movement
        expected_max_period = coupled_movement.airplane_movement.max_period
        self.assertEqual(coupled_movement.max_period, expected_max_period)

    def test_custom_delta_time(self):
        """Test that custom delta_time is used correctly."""
        coupled_movement = self.coupled_movement_with_custom_delta_time
        self.assertEqual(coupled_movement.delta_time, 0.02)

    def test_coupled_movement_with_angular_speed(self):
        """Test CoupledMovement with non zero angular speed in initial operating
        point."""
        coupled_movement = self.coupled_movement_with_angular_speed
        initial_op = coupled_movement.coupled_operating_points[0]

        # Check that the initial CoupledOperatingPoint has non zero angular speed.
        self.assertIsInstance(initial_op, ps.operating_point.CoupledOperatingPoint)
        # The fixture uses (0.5, 0.3, 0.1) for angular speed.
        np.testing.assert_array_equal(initial_op.omegas_BP1__E, (0.5, 0.3, 0.1))

    def test_coupled_movement_with_attitude_angles(self):
        """Test CoupledMovement with non zero attitude angles in initial operating
        point."""
        coupled_movement = self.coupled_movement_with_attitude_angles
        initial_op = coupled_movement.coupled_operating_points[0]

        # Check that the initial CoupledOperatingPoint has non zero attitude angles.
        self.assertIsInstance(initial_op, ps.operating_point.CoupledOperatingPoint)
        # The fixture uses (15.0, 10.0, 5.0) for attitude angles.
        np.testing.assert_array_equal(
            initial_op.angles_E_to_BP1_izyx, (15.0, 10.0, 5.0)
        )

    def test_full_coupled_movement(self):
        """Test CoupledMovement with all parameters set to non zero values."""
        coupled_movement = self.full_coupled_movement
        initial_op = coupled_movement.coupled_operating_points[0]

        # Check delta_time.
        self.assertEqual(coupled_movement.delta_time, 0.05)

        # Check prescribed and free steps.
        self.assertEqual(coupled_movement.prescribed_num_steps, 10)
        self.assertEqual(coupled_movement.free_num_steps, 15)
        self.assertEqual(coupled_movement.num_steps, 25)

        # Check that the initial CoupledOperatingPoint has expected values.
        self.assertEqual(initial_op.vCg__E, 15.0)
        self.assertEqual(initial_op.alpha, 10.0)
        self.assertEqual(initial_op.beta, 5.0)
        self.assertEqual(initial_op.externalFX_W, 25.0)

    def test_type_4_to_5_transition_raises_error(self):
        """Test that a Wing transitioning from type 4 to type 5 symmetry raises
        an error in CoupledMovement."""
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

        # Create a CoupledOperatingPoint.
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        # Attempting to create a CoupledMovement should raise a ValueError.
        with self.assertRaises(ValueError) as context:
            ps.movements.movement.CoupledMovement(
                airplane_movement=airplane_movement,
                initial_coupled_operating_point=initial_coupled_operating_point,
                delta_time=0.1,
                prescribed_num_steps=5,
                free_num_steps=5,
            )

        # Verify the error message mentions symmetry.
        self.assertIn("symmetry", str(context.exception).lower())

    def test_static_coupled_movement_with_symmetric_wing_succeeds(self):
        """Test that a CoupledMovement with a symmetric Wing but no motion succeeds."""
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

        # Create a CoupledOperatingPoint.
        initial_coupled_operating_point = (
            operating_point_fixtures.make_basic_coupled_operating_point_fixture()
        )

        # Creating a CoupledMovement should succeed without raising an error.
        coupled_movement = ps.movements.movement.CoupledMovement(
            airplane_movement=airplane_movement,
            initial_coupled_operating_point=initial_coupled_operating_point,
            delta_time=0.1,
            prescribed_num_steps=3,
            free_num_steps=3,
        )

        # Verify the CoupledMovement was created successfully.
        self.assertIsInstance(coupled_movement, ps.movements.movement.CoupledMovement)


if __name__ == "__main__":
    unittest.main()
