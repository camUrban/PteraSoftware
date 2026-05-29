"""This module contains classes to test AeroelasticWingCrossSectionMovements."""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import (
    aeroelastic_wing_cross_section_movement_fixtures,
    geometry_fixtures,
)


class TestAeroelasticWingCrossSectionMovement(unittest.TestCase):
    """This is a class with functions to test AeroelasticWingCrossSectionMovements."""

    def test_is_subclass_of_core(self):
        """Test that AeroelasticWingCrossSectionMovement is a subclass of
        CoreWingCrossSectionMovement.
        """
        self.assertTrue(
            issubclass(
                ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement,
                ps._core.CoreWingCrossSectionMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that AeroelasticWingCrossSectionMovement instantiation returns an
        AeroelasticWingCrossSectionMovement.
        """
        base_wing_cross_section = (
            geometry_fixtures.make_root_wing_cross_section_fixture()
        )
        aeroelastic_wing_cross_section_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=base_wing_cross_section,
        )
        self.assertIsInstance(
            aeroelastic_wing_cross_section_movement,
            ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement,
        )

    def test_generate_wing_cross_sections_returns_wing_cross_sections(self):
        """Test that generate_wing_cross_sections returns WingCrossSections when called
        through the public class.
        """
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_basic_aeroelastic_wing_cross_section_movement_fixture()
        )
        wing_cross_sections = (
            aeroelastic_wing_cross_section_movement.generate_wing_cross_sections(
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
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_basic_aeroelastic_wing_cross_section_movement_fixture()
        )
        wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=2, delta_time=0.01
        )
        self.assertIsInstance(
            wing_cross_section,
            ps.geometry.wing_cross_section.WingCrossSection,
        )

    def test_no_deformation_matches_static_base(self):
        """Test that, with no prescribed movement and no deformation, the generated
        WingCrossSection matches the base WingCrossSection.
        """
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_static_aeroelastic_wing_cross_section_movement_fixture()
        )
        base_wing_cross_section = (
            aeroelastic_wing_cross_section_movement.base_wing_cross_section
        )

        wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=3, delta_time=0.01
        )

        npt.assert_allclose(
            wing_cross_section.Lp_Wcsp_Lpp,
            base_wing_cross_section.Lp_Wcsp_Lpp,
            rtol=1e-10,
            atol=1e-14,
        )
        npt.assert_allclose(
            wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            rtol=1e-10,
            atol=1e-14,
        )

    def test_default_deformation_matches_explicit_none(self):
        """Test that omitting deformation_angles_ixyz produces the same result as
        explicitly passing None.
        """
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_basic_aeroelastic_wing_cross_section_movement_fixture()
        )

        default_wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=2, delta_time=0.01
        )
        explicit_none_wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=2, delta_time=0.01, deformation_angles_ixyz=None
        )

        npt.assert_array_equal(
            default_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            explicit_none_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
        )
        npt.assert_array_equal(
            default_wing_cross_section.Lp_Wcsp_Lpp,
            explicit_none_wing_cross_section.Lp_Wcsp_Lpp,
        )

    def test_deformation_adds_to_prescribed_angles(self):
        """Test that deformation_angles_ixyz is added to the prescribed
        angles_Wcsp_to_Wcs_ixyz.
        """
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_static_aeroelastic_wing_cross_section_movement_fixture()
        )
        base_wing_cross_section = (
            aeroelastic_wing_cross_section_movement.base_wing_cross_section
        )
        deformation_angles_ixyz = np.array([5.0, 5.0, 5.0], dtype=float)

        wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=1,
            delta_time=0.01,
            deformation_angles_ixyz=deformation_angles_ixyz,
        )

        # With static movement, the prescribed angles equal the base angles, so the
        # result should be the base angles plus the deformation.
        expected_angles = (
            base_wing_cross_section.angles_Wcsp_to_Wcs_ixyz + deformation_angles_ixyz
        )
        npt.assert_allclose(
            wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            expected_angles,
            rtol=1e-10,
            atol=1e-14,
        )

    def test_deformation_leaves_Lp_unchanged(self):
        """Test that deformation_angles_ixyz does not affect the generated
        Lp_Wcsp_Lpp.
        """
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_static_aeroelastic_wing_cross_section_movement_fixture()
        )
        base_wing_cross_section = (
            aeroelastic_wing_cross_section_movement.base_wing_cross_section
        )
        deformation_angles_ixyz = np.array([5.0, 5.0, 5.0], dtype=float)

        wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=1,
            delta_time=0.01,
            deformation_angles_ixyz=deformation_angles_ixyz,
        )

        npt.assert_allclose(
            wing_cross_section.Lp_Wcsp_Lpp,
            base_wing_cross_section.Lp_Wcsp_Lpp,
            rtol=1e-10,
            atol=1e-14,
        )

    def test_zero_deformation_matches_no_deformation(self):
        """Test that a zero deformation produces the same result as no deformation."""
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_basic_aeroelastic_wing_cross_section_movement_fixture()
        )

        no_deformation_wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=2, delta_time=0.01
        )
        zero_deformation_wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=2,
            delta_time=0.01,
            deformation_angles_ixyz=np.zeros(3, dtype=float),
        )

        npt.assert_allclose(
            zero_deformation_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            no_deformation_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            rtol=1e-10,
            atol=1e-14,
        )

    def test_deformation_adds_on_top_of_oscillation(self):
        """Test that deformation is added on top of the oscillating prescribed
        angles.
        """
        aeroelastic_wing_cross_section_movement = (
            aeroelastic_wing_cross_section_movement_fixtures.make_basic_aeroelastic_wing_cross_section_movement_fixture()
        )
        deformation_angles_ixyz = np.array([3.0, -2.0, 1.0], dtype=float)

        prescribed_wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=3, delta_time=0.01
        )
        deformed_wing_cross_section = aeroelastic_wing_cross_section_movement.generate_wing_cross_section_at_time_step(
            step=3,
            delta_time=0.01,
            deformation_angles_ixyz=deformation_angles_ixyz,
        )

        # The deformed angles should differ from the prescribed angles by exactly the
        # deformation.
        npt.assert_allclose(
            deformed_wing_cross_section.angles_Wcsp_to_Wcs_ixyz
            - prescribed_wing_cross_section.angles_Wcsp_to_Wcs_ixyz,
            deformation_angles_ixyz,
            rtol=1e-10,
            atol=1e-14,
        )


if __name__ == "__main__":
    unittest.main()
