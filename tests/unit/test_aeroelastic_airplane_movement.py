"""This module contains classes to test the AeroelasticAirplaneMovement class."""

import unittest
from typing import Any

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import (
    aeroelastic_airplane_movement_fixtures,
    core_wing_movement_fixtures,
    geometry_fixtures,
    movement_fixtures,
)


class TestAeroelasticAirplaneMovement(unittest.TestCase):
    """This is a class with functions to test AeroelasticAirplaneMovements."""

    def test_is_subclass_of_core(self) -> None:
        """Test that AeroelasticAirplaneMovement is a subclass of
        CoreAirplaneMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement,
                ps._core.CoreAirplaneMovement,
            )
        )

    def test_instantiation_returns_correct_type(self) -> None:
        """Test that AeroelasticAirplaneMovement instantiation returns an
        AeroelasticAirplaneMovement."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_static_aeroelastic_airplane_movement_fixture()
        )
        self.assertIsInstance(
            aeroelastic_airplane_movement,
            ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement,
        )

    def test_rejects_core_wing_movement_children(self) -> None:
        """Test that AeroelasticAirplaneMovement rejects CoreWingMovement instances that
        are neither AeroelasticWingMovements nor WingMovements."""
        base_airplane = geometry_fixtures.make_first_airplane_fixture()
        wing_movements: Any = [
            core_wing_movement_fixtures.make_static_core_wing_movement_fixture()
        ]
        with self.assertRaises(TypeError):
            ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=wing_movements,
            )

    def test_rejects_type_4_symmetric_base_wing_meshed_after_wing_movement(
        self,
    ) -> None:
        """Test that AeroelasticAirplaneMovement rejects an AeroelasticWingMovement
        whose base Wing has type 4 symmetry, including when that base Wing was still
        unmeshed when the AeroelasticWingMovement was constructed."""
        # Build the AeroelasticWingMovement around the raw, unmeshed Wing, whose
        # symmetry type is still None.
        base_wing = geometry_fixtures.make_type_4_wing_fixture()
        wing_cross_section_movements = [
            ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
                base_wing_cross_section=wing_cross_section
            )
            for wing_cross_section in base_wing.wing_cross_sections
        ]
        wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=base_wing,
            wing_cross_section_movements=wing_cross_section_movements,
        )

        # Constructing the base Airplane meshes the same Wing object in place,
        # determining its type 4 symmetry.
        base_airplane = ps.geometry.airplane.Airplane(
            wings=[base_wing],
            name="Type 4 Test Airplane",
            Cg_GP1_CgP1=[0.0, 0.0, 0.0],
            weight=1500.0,
        )
        self.assertEqual(base_wing.symmetry_type, 4)

        with self.assertRaises(ValueError):
            ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
                base_airplane=base_airplane,
                wing_movements=[wing_movement],
            )

    def test_wing_movements_property_returns_tuple(self) -> None:
        """Test that the wing_movements property returns a tuple of the
        AeroelasticAirplaneMovement's wing movements."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_static_aeroelastic_airplane_movement_fixture()
        )
        self.assertIsInstance(aeroelastic_airplane_movement.wing_movements, tuple)
        self.assertEqual(len(aeroelastic_airplane_movement.wing_movements), 1)

    def test_generate_airplanes_returns_airplanes(self) -> None:
        """Test that generate_airplanes returns Airplanes when called through the public
        class."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_basic_aeroelastic_airplane_movement_fixture()
        )
        airplanes = aeroelastic_airplane_movement.generate_airplanes(
            num_steps=5, delta_time=0.01
        )
        self.assertEqual(len(airplanes), 5)
        for airplane in airplanes:
            self.assertIsInstance(
                airplane,
                ps.geometry.airplane.Airplane,
            )

    def test_generate_airplane_at_time_step_returns_airplane(self) -> None:
        """Test that generate_airplane_at_time_step returns an Airplane when the Wing is
        backed by an AeroelasticWingMovement."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_basic_aeroelastic_airplane_movement_fixture()
        )
        airplane = aeroelastic_airplane_movement.generate_airplane_at_time_step(
            step=2, delta_time=0.01
        )
        self.assertIsInstance(
            airplane,
            ps.geometry.airplane.Airplane,
        )


class TestAeroelasticAirplaneMovementStandardWing(unittest.TestCase):
    """This class contains unit tests for AeroelasticAirplaneMovement with a standard
    WingMovement child, exercising the else branch in generate_airplane_at_time_step."""

    def setUp(self) -> None:
        """Set up an AeroelasticAirplaneMovement backed by a standard WingMovement."""
        movement = (
            movement_fixtures.make_aeroelastic_movement_with_standard_wing_fixture()
        )
        self.airplane_movement = movement.airplane_movements[0]
        self.delta_time = movement.delta_time

    def test_generate_airplane_at_time_step_standard_wing_returns_airplane(
        self,
    ) -> None:
        """Test that generate_airplane_at_time_step returns an Airplane when the wing is
        backed by a standard WingMovement (the else branch)."""
        result = self.airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=self.delta_time
        )

        self.assertIsInstance(result, ps.geometry.airplane.Airplane)

    def test_generate_airplane_at_time_step_standard_wing_returns_one_wing(
        self,
    ) -> None:
        """Test that the Airplane returned via the standard WingMovement path contains
        exactly one Wing."""
        result = self.airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=self.delta_time
        )

        self.assertEqual(len(result.wings), 1)

    def test_generate_airplane_at_time_step_standard_wing_wing_has_panels(self) -> None:
        """Test that the Wing produced via the standard WingMovement path has its panels
        populated (not None)."""
        result = self.airplane_movement.generate_airplane_at_time_step(
            step=0, delta_time=self.delta_time
        )

        self.assertIsNotNone(result.wings[0].panels)

    def test_generate_airplane_at_time_step_standard_wing_deformation_ignored(
        self,
    ) -> None:
        """Test that passing deformationAngles_Wcsp_to_Wcs_ixyz=None produces the same
        Airplane as calling without deformation for a standard WingMovement."""
        result_no_kwarg = self.airplane_movement.generate_airplane_at_time_step(
            step=1, delta_time=self.delta_time
        )
        result_none_kwarg = self.airplane_movement.generate_airplane_at_time_step(
            step=1,
            delta_time=self.delta_time,
            deformationAngles_Wcsp_to_Wcs_ixyz=None,
        )

        self.assertEqual(len(result_no_kwarg.wings), len(result_none_kwarg.wings))
        self.assertIsNotNone(result_none_kwarg.wings[0].panels)


class TestAeroelasticAirplaneMovementDeformation(unittest.TestCase):
    """This is a class with functions to test the structural deformation behavior of
    AeroelasticAirplaneMovement.generate_airplane_at_time_step."""

    def test_default_deformation_matches_explicit_none(self) -> None:
        """Test that omitting deformationAngles_Wcsp_to_Wcs_ixyz produces the same
        result as explicitly passing None."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_basic_aeroelastic_airplane_movement_fixture()
        )

        default_airplane = aeroelastic_airplane_movement.generate_airplane_at_time_step(
            step=2, delta_time=0.01
        )
        explicit_none_airplane = (
            aeroelastic_airplane_movement.generate_airplane_at_time_step(
                step=2, delta_time=0.01, deformationAngles_Wcsp_to_Wcs_ixyz=None
            )
        )

        npt.assert_array_equal(
            default_airplane.Cg_GP1_CgP1,
            explicit_none_airplane.Cg_GP1_CgP1,
        )
        default_wing_cross_sections = default_airplane.wings[0].wing_cross_sections
        explicit_none_wing_cross_sections = explicit_none_airplane.wings[
            0
        ].wing_cross_sections
        for index in range(len(default_wing_cross_sections)):
            npt.assert_array_equal(
                default_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                explicit_none_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
            )

    def test_deformation_adds_to_prescribed_angles(self) -> None:
        """Test that each row of a Wing's deformation is added to the corresponding
        WingCrossSection's prescribed angles_Wcsp_to_Wcs_ixyz."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_static_aeroelastic_airplane_movement_fixture()
        )
        base_airplane = aeroelastic_airplane_movement.base_airplane
        # The root (index 0) deformation must stay zero, since the clamped root
        # WingCrossSection's angles_Wcsp_to_Wcs_ixyz must remain (0, 0, 0).
        deformationAngles_Wcsp_to_Wcs_ixyz: list[np.ndarray | None] = [
            np.array([[0.0, 0.0, 0.0], [3.0, -2.0, 1.0]], dtype=float)
        ]

        airplane = aeroelastic_airplane_movement.generate_airplane_at_time_step(
            step=1,
            delta_time=0.01,
            deformationAngles_Wcsp_to_Wcs_ixyz=deformationAngles_Wcsp_to_Wcs_ixyz,
        )

        # With static movement, the prescribed angles equal the base angles, so each
        # WingCrossSection's result should be the base angles plus that row's
        # deformation.
        base_wing_cross_sections = base_airplane.wings[0].wing_cross_sections
        wing_cross_sections = airplane.wings[0].wing_cross_sections
        deformation_angles = deformationAngles_Wcsp_to_Wcs_ixyz[0]
        assert deformation_angles is not None
        for index in range(len(base_wing_cross_sections)):
            expected_angles = (
                base_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz
                + deformation_angles[index]
            )
            npt.assert_allclose(
                wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                expected_angles,
                rtol=1e-10,
                atol=1e-14,
            )

    def test_zero_deformation_matches_no_deformation(self) -> None:
        """Test that a zero deformation produces the same result as no deformation."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_basic_aeroelastic_airplane_movement_fixture()
        )

        no_deformation_airplane = (
            aeroelastic_airplane_movement.generate_airplane_at_time_step(
                step=2, delta_time=0.01
            )
        )
        zero_deformation_airplane = (
            aeroelastic_airplane_movement.generate_airplane_at_time_step(
                step=2,
                delta_time=0.01,
                deformationAngles_Wcsp_to_Wcs_ixyz=[np.zeros((2, 3), dtype=float)],
            )
        )

        no_deformation_wing_cross_sections = no_deformation_airplane.wings[
            0
        ].wing_cross_sections
        zero_deformation_wing_cross_sections = zero_deformation_airplane.wings[
            0
        ].wing_cross_sections
        for index in range(len(no_deformation_wing_cross_sections)):
            npt.assert_allclose(
                zero_deformation_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                no_deformation_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                rtol=1e-10,
                atol=1e-14,
            )

    def test_deformation_adds_on_top_of_oscillation(self) -> None:
        """Test that deformation is added on top of the oscillating prescribed
        WingCrossSection angles."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_basic_aeroelastic_airplane_movement_fixture()
        )
        # The root (index 0) deformation must stay zero, since the clamped root
        # WingCrossSection's angles_Wcsp_to_Wcs_ixyz must remain (0, 0, 0).
        deformationAngles_Wcsp_to_Wcs_ixyz: list[np.ndarray | None] = [
            np.array([[0.0, 0.0, 0.0], [-1.0, 4.0, -3.0]], dtype=float)
        ]

        prescribed_airplane = (
            aeroelastic_airplane_movement.generate_airplane_at_time_step(
                step=3, delta_time=0.01
            )
        )
        deformed_airplane = (
            aeroelastic_airplane_movement.generate_airplane_at_time_step(
                step=3,
                delta_time=0.01,
                deformationAngles_Wcsp_to_Wcs_ixyz=deformationAngles_Wcsp_to_Wcs_ixyz,
            )
        )

        # Each WingCrossSection's deformed angles should differ from its prescribed
        # angles by exactly that row's deformation.
        prescribed_wing_cross_sections = prescribed_airplane.wings[
            0
        ].wing_cross_sections
        deformed_wing_cross_sections = deformed_airplane.wings[0].wing_cross_sections
        deformation_angles = deformationAngles_Wcsp_to_Wcs_ixyz[0]
        assert deformation_angles is not None
        for index in range(len(deformed_wing_cross_sections)):
            npt.assert_allclose(
                deformed_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz
                - prescribed_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                deformation_angles[index],
                rtol=1e-10,
                atol=1e-14,
            )


class TestAeroelasticAirplaneMovementMixedWings(unittest.TestCase):
    """This is a class with functions to test that AeroelasticAirplaneMovement threads
    per Wing deformation only to its AeroelasticWingMovement children, leaving its
    standard WingMovement children unaffected."""

    def test_deformation_threaded_to_aeroelastic_wing_only(self) -> None:
        """Test that, with one AeroelasticWingMovement child and one WingMovement child,
        deformation is applied to the aeroelastic Wing but the standard Wing matches its
        base."""
        aeroelastic_airplane_movement = (
            aeroelastic_airplane_movement_fixtures.make_mixed_wing_aeroelastic_airplane_movement_fixture()
        )
        base_airplane = aeroelastic_airplane_movement.base_airplane
        # Deform only the aeroelastic Wing's tip WingCrossSection. The standard Wing's
        # entry is ignored, so None is passed for it.
        deformationAngles_Wcsp_to_Wcs_ixyz = [
            np.array([[0.0, 0.0, 0.0], [3.0, -2.0, 1.0]], dtype=float),
            None,
        ]

        airplane = aeroelastic_airplane_movement.generate_airplane_at_time_step(
            step=1,
            delta_time=0.01,
            deformationAngles_Wcsp_to_Wcs_ixyz=deformationAngles_Wcsp_to_Wcs_ixyz,
        )

        # The aeroelastic Wing (index 0) receives the deformation on top of its static
        # base angles.
        base_aeroelastic_wing_cross_sections = base_airplane.wings[
            0
        ].wing_cross_sections
        aeroelastic_wing_cross_sections = airplane.wings[0].wing_cross_sections
        deformation_angles = deformationAngles_Wcsp_to_Wcs_ixyz[0]
        assert deformation_angles is not None
        for index in range(len(base_aeroelastic_wing_cross_sections)):
            expected_angles = (
                base_aeroelastic_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz
                + deformation_angles[index]
            )
            npt.assert_allclose(
                aeroelastic_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                expected_angles,
                rtol=1e-10,
                atol=1e-14,
            )

        # The standard Wing (index 1) is advanced without deformation, so it matches its
        # static base.
        base_standard_wing_cross_sections = base_airplane.wings[1].wing_cross_sections
        standard_wing_cross_sections = airplane.wings[1].wing_cross_sections
        for index in range(len(base_standard_wing_cross_sections)):
            npt.assert_allclose(
                standard_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                base_standard_wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                rtol=1e-10,
                atol=1e-14,
            )
