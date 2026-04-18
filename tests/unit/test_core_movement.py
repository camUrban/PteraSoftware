"""This module contains classes to test CoreMovements."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import (
    core_airplane_movement_fixtures,
    core_movement_fixtures,
    core_operating_point_movement_fixtures,
    geometry_fixtures,
    operating_point_fixtures,
)


class TestCoreMovement(unittest.TestCase):
    """This is a class with functions to test CoreMovements."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all CoreMovement tests."""
        cls.static_core_movement = (
            core_movement_fixtures.make_static_core_movement_fixture()
        )
        cls.basic_core_movement = (
            core_movement_fixtures.make_basic_core_movement_fixture()
        )

    def test_static_property_for_static_core_movement(self):
        """Test that static property returns True for static CoreMovement."""
        self.assertTrue(self.static_core_movement.static)

    def test_static_property_for_non_static_core_movement(self):
        """Test that static property returns False for non static CoreMovement."""
        self.assertFalse(self.basic_core_movement.static)

    def test_max_period_for_static_core_movement(self):
        """Test that max_period returns 0.0 for static CoreMovement."""
        self.assertEqual(self.static_core_movement.max_period, 0.0)

    def test_max_period_for_non_static_core_movement(self):
        """Test that max_period returns correct value for non static CoreMovement."""
        # The basic_core_movement has period of 2.0 for all motion.
        self.assertEqual(self.basic_core_movement.max_period, 2.0)

    def test_lcm_period_for_static_core_movement(self):
        """Test that lcm_period returns 0.0 for static CoreMovement."""
        self.assertEqual(self.static_core_movement.lcm_period, 0.0)

    def test_lcm_period_for_single_period_core_movement(self):
        """Test that lcm_period returns correct value when all periods are the same."""
        # The basic_core_movement has period of 2.0 for all motion.
        # LCM of identical periods should equal that period.
        self.assertEqual(self.basic_core_movement.lcm_period, 2.0)

    def test_lcm_period_with_multiple_wings_same_airplane(self):
        """Test that lcm_period collects all periods, not just max from each
        CoreAirplaneMovement.

        This test creates a single Airplane with two Wings having different periods
        (3.0 s and 4.0 s). The correct LCM is 12.0 s. If the implementation only uses
        max_period from the CoreAirplaneMovement, lcm_period would incorrectly return
        4.0 s instead of 12.0 s.
        """
        # Create two Wings for the same Airplane.
        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing_1, base_wing_2],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Wing_1: tip CoreWingCrossSectionMovement has period 3.0 s.
        wing_cross_section_movements_wing_1 = [
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_1 = ps._core.CoreWingMovement(
            base_wing=base_wing_1,
            wing_cross_section_movements=wing_cross_section_movements_wing_1,
        )

        # Wing_2: tip CoreWingCrossSectionMovement has period 4.0 s.
        wing_cross_section_movements_wing_2 = [
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(4.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_2 = ps._core.CoreWingMovement(
            base_wing=base_wing_2,
            wing_cross_section_movements=wing_cross_section_movements_wing_2,
        )

        airplane_movement = ps._core.CoreAirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement_1, wing_movement_2],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        operating_point_movement = (
            core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The max_period should be 4.0 (the max of 3.0 and 4.0).
        self.assertEqual(core_movement.max_period, 4.0)

        # The lcm_period should be LCM(3.0, 4.0) = 12.0, not 4.0. This test will fail
        # if lcm_period only uses max_period from each CoreAirplaneMovement instead of
        # collecting all individual periods.
        self.assertEqual(core_movement.lcm_period, 12.0)

    def test_lcm_period_with_multiple_cross_sections_same_wing(self):
        """Test that lcm_period collects all periods from
        CoreWingCrossSectionMovements.

        This test creates a single Wing with three WingCrossSections having different
        periods (root static, middle 3.0 s, tip 4.0 s). The correct LCM is 12.0 s. If
        the implementation only uses max_period from each CoreWingMovement, lcm_period
        would incorrectly return 4.0 s instead of 12.0 s.
        """
        # Create a Wing with three WingCrossSections.
        test_airfoil = ps.geometry.airfoil.Airfoil(name="naca2412")

        root_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil,
            num_spanwise_panels=4,
            chord=2.0,
            Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )

        middle_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil,
            num_spanwise_panels=4,
            chord=1.5,
            Lp_Wcsp_Lpp=(0.0, 1.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )

        tip_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            airfoil=test_airfoil,
            num_spanwise_panels=None,
            chord=1.0,
            Lp_Wcsp_Lpp=(0.0, 1.5, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        )

        base_wing = ps.geometry.wing.Wing(
            wing_cross_sections=[
                root_wing_cross_section,
                middle_wing_cross_section,
                tip_wing_cross_section,
            ],
            name="Test Wing",
        )

        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Test Airplane",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Root CoreWingCrossSectionMovement must be static.
        # Middle has period 3.0 s, tip has period 4.0 s.
        wing_cross_section_movements = [
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing.wing_cross_sections[2],
                periodLp_Wcsp_Lpp=(4.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement = ps._core.CoreWingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wing_cross_section_movements,
        )

        airplane_movement = ps._core.CoreAirplaneMovement(
            base_airplane=base_airplane,
            wing_movements=[wing_movement],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        operating_point_movement = (
            core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The max_period should be 4.0 (the max of 3.0 and 4.0).
        self.assertEqual(core_movement.max_period, 4.0)

        # The lcm_period should be LCM(3.0, 4.0) = 12.0, not 4.0. This test will fail
        # if lcm_period only uses max_period from each CoreWingMovement instead of
        # collecting all individual periods from CoreWingCrossSectionMovements.
        self.assertEqual(core_movement.lcm_period, 12.0)

    def test_lcm_period_with_multiple_airplanes(self):
        """Test that lcm_period calculates LCM correctly with multiple periods."""
        # Create CoreAirplaneMovements with different periods.

        base_wing_1 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_airplane_1 = ps.geometry.airplane.Airplane(
            wings=[base_wing_1],
            name="Test Airplane 1",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Make CoreWingCrossSectionMovements for the first Airplane's Wing's root and
        # tip WingCrossSections. The root CoreWingCrossSectionMovement must be static.
        # The tip CoreWingCrossSectionMovement will have a period of 2.0 s.
        wing_cross_section_movements_1 = [
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_1.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(2.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_1 = ps._core.CoreWingMovement(
            base_wing=base_wing_1,
            wing_cross_section_movements=wing_cross_section_movements_1,
        )

        airplane_movement_1 = ps._core.CoreAirplaneMovement(
            base_airplane=base_airplane_1,
            wing_movements=[wing_movement_1],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        base_wing_2 = geometry_fixtures.make_simple_tapered_wing_fixture()
        base_airplane_2 = ps.geometry.airplane.Airplane(
            wings=[base_wing_2],
            name="Test Airplane 2",
            Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        # Make CoreWingCrossSectionMovements for the second Airplane's Wing's root and
        # tip WingCrossSections. The root CoreWingCrossSectionMovement must be static.
        # The tip CoreWingCrossSectionMovement will have a period of 3.0 s.
        wing_cross_section_movements_2 = [
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[0],
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
            ),
            ps._core.CoreWingCrossSectionMovement(
                base_wing_cross_section=base_wing_2.wing_cross_sections[1],
                periodLp_Wcsp_Lpp=(3.0, 0.0, 0.0),
                ampLp_Wcsp_Lpp=(0.1, 0.0, 0.0),
            ),
        ]

        wing_movement_2 = ps._core.CoreWingMovement(
            base_wing=base_wing_2,
            wing_cross_section_movements=wing_cross_section_movements_2,
        )

        airplane_movement_2 = ps._core.CoreAirplaneMovement(
            base_airplane=base_airplane_2,
            wing_movements=[wing_movement_2],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )

        operating_point_movement = (
            core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=[airplane_movement_1, airplane_movement_2],
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The LCM of 2.0 and 3.0 should be 6.0.
        self.assertEqual(core_movement.lcm_period, 6.0)

        # The max_period should still be 3.0.
        self.assertEqual(core_movement.max_period, 3.0)

    def test_max_period_with_multiple_airplanes(self):
        """Test that max_period returns maximum across all CoreAirplaneMovements."""
        # Create a static and a basic CoreAirplaneMovement.
        static_airplane_movement = (
            core_airplane_movement_fixtures.make_static_core_airplane_movement_fixture()
        )
        basic_airplane_movement = (
            core_airplane_movement_fixtures.make_basic_core_airplane_movement_fixture()
        )

        operating_point_movement = (
            core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=[static_airplane_movement, basic_airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=0.01,
            num_steps=1,
        )

        # The CoreMovement has one static (period 0.0) and one with period 2.0.
        # Should return 2.0.
        self.assertEqual(core_movement.max_period, 2.0)


class TestCoreMovementImmutability(unittest.TestCase):
    """Tests for CoreMovement attribute immutability."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures once for all immutability tests."""
        cls.basic_core_movement = (
            core_movement_fixtures.make_basic_core_movement_fixture()
        )

    def test_immutable_airplane_movements_property(self):
        """Test that airplane_movements property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_core_movement.airplane_movements = []

    def test_immutable_airplane_movements_tuple(self):
        """Test that airplane_movements returns a tuple (immutable sequence)."""
        airplane_movements = self.basic_core_movement.airplane_movements
        self.assertIsInstance(airplane_movements, tuple)

    def test_immutable_operating_point_movement_property(self):
        """Test that operating_point_movement property is read only."""
        operating_point_movement = (
            core_operating_point_movement_fixtures.make_static_core_operating_point_movement_fixture()
        )
        with self.assertRaises(AttributeError):
            self.basic_core_movement.operating_point_movement = operating_point_movement

    def test_immutable_delta_time_property(self):
        """Test that delta_time property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_core_movement.delta_time = 0.05

    def test_immutable_num_steps_property(self):
        """Test that num_steps property is read only."""
        with self.assertRaises(AttributeError):
            self.basic_core_movement.num_steps = 100


class TestCoreMovementWithOperatingPointMovementPeriod(unittest.TestCase):
    """Tests for CoreMovement when CoreOperatingPointMovement has non zero period."""

    def test_lcm_period_includes_operating_point_movement_period(self):
        """Test that lcm_period includes the CoreOperatingPointMovement period."""
        # Create a static CoreAirplaneMovement.
        airplane_movements = [
            core_airplane_movement_fixtures.make_static_core_airplane_movement_fixture()
        ]

        # Create a CoreOperatingPointMovement with a non zero period.
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = ps._core.CoreOperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=1.0,
            periodVCg__E=3.0,
        )

        # Create the CoreMovement with explicit num_steps to avoid auto calculation.
        core_movement = ps._core.CoreMovement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The lcm_period should be 3.0 (from CoreOperatingPointMovement).
        # Since the CoreAirplaneMovement is static (period 0.0), the only period is 3.0.
        self.assertEqual(core_movement.lcm_period, 3.0)

    def test_lcm_period_combines_airplane_and_operating_point_periods(self):
        """Test that lcm_period combines periods from CoreAirplaneMovement and
        CoreOperatingPointMovement."""
        # Create a non static CoreAirplaneMovement with period 2.0.
        airplane_movements = [
            core_airplane_movement_fixtures.make_basic_core_airplane_movement_fixture()
        ]

        # Create a CoreOperatingPointMovement with period 3.0.
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = ps._core.CoreOperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=1.0,
            periodVCg__E=3.0,
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The lcm_period should be LCM(2.0, 3.0) = 6.0.
        self.assertEqual(core_movement.lcm_period, 6.0)

    def test_min_period_includes_operating_point_movement_period(self):
        """Test that min_period includes the CoreOperatingPointMovement period."""
        # Create a non static CoreAirplaneMovement with period 2.0.
        airplane_movements = [
            core_airplane_movement_fixtures.make_basic_core_airplane_movement_fixture()
        ]

        # Create a CoreOperatingPointMovement with a shorter period (1.5).
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = ps._core.CoreOperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=1.0,
            periodVCg__E=1.5,
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The min_period should be 1.5 (from CoreOperatingPointMovement, smaller than
        # 2.0).
        self.assertEqual(core_movement.min_period, 1.5)

    def test_max_period_includes_operating_point_movement_period(self):
        """Test that max_period includes the CoreOperatingPointMovement period."""
        # Create a non static CoreAirplaneMovement with period 2.0.
        airplane_movements = [
            core_airplane_movement_fixtures.make_basic_core_airplane_movement_fixture()
        ]

        # Create a CoreOperatingPointMovement with a longer period (5.0).
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
        operating_point_movement = ps._core.CoreOperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=1.0,
            periodVCg__E=5.0,
        )

        core_movement = ps._core.CoreMovement(
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            delta_time=0.1,
            num_steps=1,
        )

        # The max_period should be 5.0 (from CoreOperatingPointMovement, larger than
        # 2.0).
        self.assertEqual(core_movement.max_period, 5.0)


if __name__ == "__main__":
    unittest.main()
