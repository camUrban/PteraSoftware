"""This module contains classes to test AeroelasticMovements."""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import (
    aeroelastic_airplane_movement_fixtures,
    aeroelastic_operating_point_movement_fixtures,
    airplane_movement_fixtures,
    movement_fixtures,
    operating_point_fixtures,
)


class TestAeroelasticMovement(unittest.TestCase):
    """This is a class with functions to test AeroelasticMovements."""

    def test_is_subclass_of_core(self):
        """Test that AeroelasticMovement is a subclass of CoreMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.aeroelastic_movement.AeroelasticMovement,
                ps._core.CoreMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that AeroelasticMovement instantiation returns an
        AeroelasticMovement.
        """
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )
        self.assertIsInstance(
            aeroelastic_movement,
            ps.movements.aeroelastic_movement.AeroelasticMovement,
        )

    def test_rejects_non_aeroelastic_airplane_movement(self):
        """Test that AeroelasticMovement rejects a CoreAirplaneMovement that is not an
        AeroelasticAirplaneMovement.
        """
        # A standard AirplaneMovement is a CoreAirplaneMovement but not an
        # AeroelasticAirplaneMovement, so it must be rejected.
        airplane_movements = [
            airplane_movement_fixtures.make_static_airplane_movement_fixture()
        ]
        operating_point_movement = (
            aeroelastic_operating_point_movement_fixtures.make_sine_spacing_aeroelastic_operating_point_movement_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.aeroelastic_movement.AeroelasticMovement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=0.1,
                num_steps=1,
            )

    def test_rejects_non_aeroelastic_operating_point_movement(self):
        """Test that AeroelasticMovement rejects a CoreOperatingPointMovement that is
        not an AeroelasticOperatingPointMovement.
        """
        # A standard OperatingPointMovement is a CoreOperatingPointMovement but not an
        # AeroelasticOperatingPointMovement, so it must be rejected.
        airplane_movements = [
            aeroelastic_airplane_movement_fixtures.make_static_aeroelastic_airplane_movement_fixture()
        ]
        operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point_fixtures.make_basic_operating_point_fixture()
        )

        with self.assertRaises(TypeError):
            ps.movements.aeroelastic_movement.AeroelasticMovement(
                airplane_movements=airplane_movements,
                operating_point_movement=operating_point_movement,
                delta_time=0.1,
                num_steps=1,
            )

    def test_airplane_movements_property_returns_tuple(self):
        """Test that the airplane_movements property returns a tuple of the
        AeroelasticMovement's AeroelasticAirplaneMovements.
        """
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )
        self.assertIsInstance(aeroelastic_movement.airplane_movements, tuple)
        self.assertEqual(len(aeroelastic_movement.airplane_movements), 1)
        for airplane_movement in aeroelastic_movement.airplane_movements:
            self.assertIsInstance(
                airplane_movement,
                ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement,
            )

    def test_operating_point_movement_property_returns_aeroelastic(self):
        """Test that the operating_point_movement property returns an
        AeroelasticOperatingPointMovement.
        """
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )
        self.assertIsInstance(
            aeroelastic_movement.operating_point_movement,
            ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement,
        )

    def test_operating_points_property_returns_tuple(self):
        """Test that the operating_points property returns a tuple of the
        OperatingPoints pre generated for every time step.
        """
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )
        self.assertIsInstance(aeroelastic_movement.operating_points, tuple)
        # The OperatingPoints are prescribed and generated upfront, one per time step.
        self.assertEqual(
            len(aeroelastic_movement.operating_points),
            aeroelastic_movement.num_steps,
        )
        for operating_point in aeroelastic_movement.operating_points:
            self.assertIsInstance(
                operating_point,
                ps.operating_point.OperatingPoint,
            )

    def test_generate_airplane_at_time_step_returns_airplane(self):
        """Test that generate_airplane_at_time_step returns an Airplane."""
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )
        airplane = aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0, step=2
        )
        self.assertIsInstance(
            airplane,
            ps.geometry.airplane.Airplane,
        )

    def test_generate_airplane_at_time_step_delegates_with_internal_delta_time(self):
        """Test that generate_airplane_at_time_step delegates to the selected
        AeroelasticAirplaneMovement using the AeroelasticMovement's own delta_time.
        """
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )

        # The AeroelasticMovement does not take a delta_time argument here; it must use
        # its stored delta_time when delegating to the AeroelasticAirplaneMovement.
        result = aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0, step=2
        )
        expected = aeroelastic_movement.airplane_movements[
            0
        ].generate_airplane_at_time_step(2, aeroelastic_movement.delta_time)

        npt.assert_array_equal(result.Cg_GP1_CgP1, expected.Cg_GP1_CgP1)
        result_wing_cross_sections = result.wings[0].wing_cross_sections
        expected_wing_cross_sections = expected.wings[0].wing_cross_sections
        for index in range(len(expected_wing_cross_sections)):
            npt.assert_array_equal(
                result_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                expected_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
            )

    def test_generate_airplane_at_time_step_applies_deformation(self):
        """Test that generate_airplane_at_time_step threads per Wing deformation down
        to the AeroelasticAirplaneMovement, so each WingCrossSection's angles differ
        from the prescribed angles by exactly that row's deformation.
        """
        aeroelastic_movement = (
            movement_fixtures.make_basic_aeroelastic_movement_fixture()
        )
        # The root (index 0) deformation must stay zero, since the clamped root
        # WingCrossSection's angles_Wcsp_to_Wcs_ixyz must remain (0, 0, 0).
        wing_deformation_angles_ixyz = [
            np.array([[0.0, 0.0, 0.0], [3.0, -2.0, 1.0]], dtype=float)
        ]

        prescribed_airplane = aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0, step=1
        )
        deformed_airplane = aeroelastic_movement.generate_airplane_at_time_step(
            airplane_movement_index=0,
            step=1,
            wing_deformation_angles_ixyz=wing_deformation_angles_ixyz,
        )

        prescribed_wing_cross_sections = prescribed_airplane.wings[
            0
        ].wing_cross_sections
        deformed_wing_cross_sections = deformed_airplane.wings[0].wing_cross_sections
        for index in range(len(deformed_wing_cross_sections)):
            npt.assert_allclose(
                deformed_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz
                - prescribed_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                wing_deformation_angles_ixyz[0][index],
                rtol=1e-10,
                atol=1e-14,
            )


if __name__ == "__main__":
    unittest.main()
