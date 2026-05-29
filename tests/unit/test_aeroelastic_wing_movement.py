"""This module contains classes to test AeroelasticWingMovements."""

import copy
import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import (
    aeroelastic_wing_movement_fixtures,
    geometry_fixtures,
    wing_cross_section_movement_fixtures,
)


def _make_base_wing_and_wing_cross_section_movements():
    """Build a base Wing and matching AeroelasticWingCrossSectionMovements.

    The AeroelasticWingCrossSectionMovements are built from the base Wing's own
    WingCrossSections so that the resulting AeroelasticWingMovement is valid.

    :return: A 2-tuple of (base_wing, wing_cross_section_movements) suitable for
        constructing an AeroelasticWingMovement.
    """
    base_wing = geometry_fixtures.make_origin_wing_fixture()
    wing_cross_section_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=wing_cross_section
        )
        for wing_cross_section in base_wing.wing_cross_sections
    ]
    return base_wing, wing_cross_section_movements


class TestAeroelasticWingMovement(unittest.TestCase):
    """This is a class with functions to test AeroelasticWingMovements."""

    def test_is_subclass_of_core(self):
        """Test that AeroelasticWingMovement is a subclass of CoreWingMovement."""
        self.assertTrue(
            issubclass(
                ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement,
                ps._core.CoreWingMovement,
            )
        )

    def test_instantiation_returns_correct_type(self):
        """Test that AeroelasticWingMovement instantiation returns an
        AeroelasticWingMovement.
        """
        base_wing, wing_cross_section_movements = (
            _make_base_wing_and_wing_cross_section_movements()
        )
        aeroelastic_wing_movement = (
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wing_cross_section_movements,
            )
        )
        self.assertIsInstance(
            aeroelastic_wing_movement,
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement,
        )

    def test_rejects_non_aeroelastic_wing_cross_section_movement_children(self):
        """Test that AeroelasticWingMovement rejects WingCrossSectionMovement instances
        that are not AeroelasticWingCrossSectionMovements.
        """
        base_wing = geometry_fixtures.make_origin_wing_fixture()
        wing_cross_section_movements = [
            wing_cross_section_movement_fixtures.make_static_wing_cross_section_movement_fixture(),
            wing_cross_section_movement_fixtures.make_static_tip_wing_cross_section_movement_fixture(),
        ]
        with self.assertRaises(TypeError):
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wing_cross_section_movements,
            )

    def test_generate_wings_returns_wings(self):
        """Test that generate_wings returns Wings when called through the public
        class.
        """
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )
        wings = aeroelastic_wing_movement.generate_wings(num_steps=5, delta_time=0.01)
        self.assertEqual(len(wings), 5)
        for wing in wings:
            self.assertIsInstance(
                wing,
                ps.geometry.wing.Wing,
            )

    def test_generate_wing_at_time_step_returns_wing(self):
        """Test that generate_wing_at_time_step returns a Wing."""
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )
        wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=2, delta_time=0.01
        )
        self.assertIsInstance(
            wing,
            ps.geometry.wing.Wing,
        )


class TestAeroelasticWingMovementDeformation(unittest.TestCase):
    """This is a class with functions to test the structural deformation behavior of
    AeroelasticWingMovement.generate_wing_at_time_step.
    """

    def test_no_deformation_matches_static_base(self):
        """Test that, with no prescribed movement and no deformation, the generated
        Wing matches the base Wing.
        """
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_static_aeroelastic_wing_movement_fixture()
        )
        base_wing = aeroelastic_wing_movement.base_wing

        wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=3, delta_time=0.01
        )

        npt.assert_allclose(
            wing.Ler_Gs_Cgs,
            base_wing.Ler_Gs_Cgs,
            rtol=1e-10,
            atol=1e-14,
        )
        npt.assert_allclose(
            wing.angles_Gs_to_Wn_ixyz,
            base_wing.angles_Gs_to_Wn_ixyz,
            rtol=1e-10,
            atol=1e-14,
        )
        for index in range(len(base_wing.wing_cross_sections)):
            npt.assert_allclose(
                wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                base_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                rtol=1e-10,
                atol=1e-14,
            )

    def test_default_deformation_matches_explicit_none(self):
        """Test that omitting deformation_angles_ixyz produces the same result as
        explicitly passing None.
        """
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )

        default_wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=2, delta_time=0.01
        )
        explicit_none_wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=2, delta_time=0.01, deformation_angles_ixyz=None
        )

        npt.assert_array_equal(
            default_wing.Ler_Gs_Cgs,
            explicit_none_wing.Ler_Gs_Cgs,
        )
        for index in range(len(default_wing.wing_cross_sections)):
            npt.assert_array_equal(
                default_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                explicit_none_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
            )

    def test_deformation_adds_to_prescribed_angles(self):
        """Test that each row of deformation_angles_ixyz is added to the corresponding
        WingCrossSection's prescribed angles_Wcsp_to_Wcs_ixyz.
        """
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_static_aeroelastic_wing_movement_fixture()
        )
        base_wing = aeroelastic_wing_movement.base_wing
        # The root (index 0) deformation must stay zero, since the clamped root
        # WingCrossSection's angles_Wcsp_to_Wcs_ixyz must remain (0, 0, 0).
        deformation_angles_ixyz = np.array(
            [[0.0, 0.0, 0.0], [3.0, -2.0, 1.0]], dtype=float
        )

        wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=1,
            delta_time=0.01,
            deformation_angles_ixyz=deformation_angles_ixyz,
        )

        # With static movement, the prescribed angles equal the base angles, so each
        # WingCrossSection's result should be the base angles plus that row's
        # deformation.
        for index in range(len(base_wing.wing_cross_sections)):
            expected_angles = (
                base_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz
                + deformation_angles_ixyz[index]
            )
            npt.assert_allclose(
                wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                expected_angles,
                rtol=1e-10,
                atol=1e-14,
            )

    def test_zero_deformation_matches_no_deformation(self):
        """Test that a zero deformation produces the same result as no deformation."""
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )

        no_deformation_wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=2, delta_time=0.01
        )
        zero_deformation_wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=2,
            delta_time=0.01,
            deformation_angles_ixyz=np.zeros((2, 3), dtype=float),
        )

        for index in range(len(no_deformation_wing.wing_cross_sections)):
            npt.assert_allclose(
                zero_deformation_wing.wing_cross_sections[
                    index
                ].angles_Wcsp_to_Wcs_ixyz,
                no_deformation_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                rtol=1e-10,
                atol=1e-14,
            )

    def test_deformation_adds_on_top_of_oscillation(self):
        """Test that deformation is added on top of the oscillating prescribed
        WingCrossSection angles.
        """
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )
        # The root (index 0) deformation must stay zero, since the clamped root
        # WingCrossSection's angles_Wcsp_to_Wcs_ixyz must remain (0, 0, 0).
        deformation_angles_ixyz = np.array(
            [[0.0, 0.0, 0.0], [-1.0, 4.0, -3.0]], dtype=float
        )

        prescribed_wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=3, delta_time=0.01
        )
        deformed_wing = aeroelastic_wing_movement.generate_wing_at_time_step(
            step=3,
            delta_time=0.01,
            deformation_angles_ixyz=deformation_angles_ixyz,
        )

        # Each WingCrossSection's deformed angles should differ from its prescribed
        # angles by exactly that row's deformation.
        for index in range(len(deformed_wing.wing_cross_sections)):
            npt.assert_allclose(
                deformed_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz
                - prescribed_wing.wing_cross_sections[index].angles_Wcsp_to_Wcs_ixyz,
                deformation_angles_ixyz[index],
                rtol=1e-10,
                atol=1e-14,
            )


class TestAeroelasticWingMovementRotationPointOffset(unittest.TestCase):
    """This is a class with functions to test that a non zero
    rotationPointOffset_Gs_Ler shifts the generated Wing's Ler_Gs_Cgs during angular
    motion.
    """

    def test_offset_has_no_effect_with_zero_angular_motion(self):
        """Test that, at a time step with zero angular motion, a rotation point offset
        leaves the generated Wing's Ler_Gs_Cgs unchanged relative to the no offset
        case.
        """
        basic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )
        offset_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_rotation_offset_aeroelastic_wing_movement_fixture()
        )

        # At step 0 the sine spaced angular motion is zero, so the rotation matrix is
        # the identity and the offset adjustment vanishes.
        basic_wing = basic_wing_movement.generate_wing_at_time_step(
            step=0, delta_time=0.1
        )
        offset_wing = offset_wing_movement.generate_wing_at_time_step(
            step=0, delta_time=0.1
        )

        npt.assert_allclose(
            offset_wing.Ler_Gs_Cgs,
            basic_wing.Ler_Gs_Cgs,
            rtol=1e-10,
            atol=1e-14,
        )

    def test_offset_shifts_ler_with_nonzero_angular_motion(self):
        """Test that, at a time step with non zero angular motion, a rotation point
        offset shifts the generated Wing's Ler_Gs_Cgs but leaves its
        angles_Gs_to_Wn_ixyz unchanged relative to the no offset case.
        """
        basic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )
        offset_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_rotation_offset_aeroelastic_wing_movement_fixture()
        )

        # At step 1 the sine spaced angular motion is non zero, so the offset
        # adjustment shifts Ler_Gs_Cgs.
        basic_wing = basic_wing_movement.generate_wing_at_time_step(
            step=1, delta_time=0.1
        )
        offset_wing = offset_wing_movement.generate_wing_at_time_step(
            step=1, delta_time=0.1
        )

        self.assertFalse(np.allclose(offset_wing.Ler_Gs_Cgs, basic_wing.Ler_Gs_Cgs))
        npt.assert_allclose(
            offset_wing.angles_Gs_to_Wn_ixyz,
            basic_wing.angles_Gs_to_Wn_ixyz,
            rtol=1e-10,
            atol=1e-14,
        )


class TestAeroelasticWingMovementSecondDerivativeValidation(unittest.TestCase):
    """This is a class with functions to test the
    spacingAnglesSecondDerivative_Gs_to_Wn_ixyz validation logic in
    AeroelasticWingMovement.__init__.
    """

    def setUp(self):
        """Set up a shared base Wing and AeroelasticWingCrossSectionMovements for each
        test.
        """
        self.base_wing, self.wing_cross_section_movements = (
            _make_base_wing_and_wing_cross_section_movements()
        )

    def _make_wing_movement(
        self,
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
    ):
        """Construct an AeroelasticWingMovement with the given angular spacing and second
        derivative arguments.

        :param spacingAnglesSecondDerivative_Gs_to_Wn_ixyz: The value to pass for the
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz parameter.
        :param spacingAngles_Gs_to_Wn_ixyz: The value to pass for the
            spacingAngles_Gs_to_Wn_ixyz parameter. The default is ("sine", "sine",
            "sine").
        :return: The constructed AeroelasticWingMovement.
        """
        return ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=self.base_wing,
            wing_cross_section_movements=self.wing_cross_section_movements,
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=spacingAngles_Gs_to_Wn_ixyz,
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
        )

    def test_non_sequence_raises_value_error(self):
        """Test that passing a non sequence value for
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz raises ValueError.
        """
        with self.assertRaises(ValueError):
            self._make_wing_movement(
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz="bad_value"
            )

    def test_wrong_length_two_elements_raises_value_error(self):
        """Test that passing a sequence with 2 elements raises ValueError."""
        with self.assertRaises(ValueError):
            self._make_wing_movement(
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[None, None]
            )

    def test_wrong_length_four_elements_raises_value_error(self):
        """Test that passing a sequence with 4 elements raises ValueError."""
        with self.assertRaises(ValueError):
            self._make_wing_movement(
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[None, None, None, None]
            )

    def test_non_callable_element_raises_type_error(self):
        """Test that passing a sequence whose element is neither callable nor None
        raises TypeError.
        """
        with self.assertRaises(TypeError):
            self._make_wing_movement(
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[42, None, None]
            )

    def test_valid_list_of_none_values_accepted(self):
        """Test that a valid 3-element list of None values is accepted and stored."""
        aeroelastic_wing_movement = self._make_wing_movement(
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[None, None, None]
        )

        self.assertEqual(
            aeroelastic_wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            (None, None, None),
        )

    def test_valid_callable_element_accepted_and_stored(self):
        """Test that a valid list with a callable at index 0 is accepted, converted to
        a tuple, and stored correctly.
        """

        def custom_spacing(t):
            return 0.0

        def deriv_func(t):
            return -1.0 * t

        aeroelastic_wing_movement = self._make_wing_movement(
            spacingAngles_Gs_to_Wn_ixyz=(custom_spacing, "sine", "sine"),
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[deriv_func, None, None],
        )

        result = aeroelastic_wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz
        self.assertIsInstance(result, tuple)
        self.assertIs(result[0], deriv_func)
        self.assertIsNone(result[1])
        self.assertIsNone(result[2])

    def test_property_returns_none_defaults_when_none_passed(self):
        """Test that the property returns (None, None, None) when
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=None is passed (the default).
        """
        aeroelastic_wing_movement = self._make_wing_movement(
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=None
        )

        self.assertEqual(
            aeroelastic_wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            (None, None, None),
        )

    def test_callable_spacing_no_derivative_raises(self):
        """Test that a custom callable spacingAngles_Gs_to_Wn_ixyz component with no
        matching second derivative raises ValueError.
        """

        def custom_spacing(t):
            return 0.0

        with self.assertRaises(ValueError):
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
                base_wing=self.base_wing,
                wing_cross_section_movements=self.wing_cross_section_movements,
                ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
                spacingAngles_Gs_to_Wn_ixyz=(custom_spacing, "sine", "sine"),
                phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            )

    def test_non_callable_spacing_with_derivative_raises(self):
        """Test that a non-None spacingAnglesSecondDerivative_Gs_to_Wn_ixyz component
        paired with a named (non callable) spacingAngles_Gs_to_Wn_ixyz component raises
        ValueError.
        """

        def deriv_func(t):
            return -1.0 * t

        with self.assertRaises(ValueError):
            self._make_wing_movement(
                spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[deriv_func, None, None],
            )


class TestAeroelasticWingMovementDeepCopy(unittest.TestCase):
    """This is a class with functions to test deep copying AeroelasticWingMovements."""

    def test_deepcopy_returns_independent_aeroelastic_wing_movement(self):
        """Test that deep copying returns an independent AeroelasticWingMovement."""
        aeroelastic_wing_movement = (
            aeroelastic_wing_movement_fixtures.make_basic_aeroelastic_wing_movement_fixture()
        )

        copied_aeroelastic_wing_movement = copy.deepcopy(aeroelastic_wing_movement)

        self.assertIsInstance(
            copied_aeroelastic_wing_movement,
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement,
        )
        self.assertIsNot(copied_aeroelastic_wing_movement, aeroelastic_wing_movement)
        self.assertIsNot(
            copied_aeroelastic_wing_movement.base_wing,
            aeroelastic_wing_movement.base_wing,
        )

    def test_deepcopy_preserves_second_derivative(self):
        """Test that deep copying preserves the
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz tuple.
        """
        base_wing, wing_cross_section_movements = (
            _make_base_wing_and_wing_cross_section_movements()
        )

        def custom_spacing(t):
            return 0.0

        def deriv_func(t):
            return -1.0 * t

        aeroelastic_wing_movement = (
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
                base_wing=base_wing,
                wing_cross_section_movements=wing_cross_section_movements,
                ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
                spacingAngles_Gs_to_Wn_ixyz=(custom_spacing, "sine", "sine"),
                phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[deriv_func, None, None],
            )
        )

        copied_aeroelastic_wing_movement = copy.deepcopy(aeroelastic_wing_movement)

        self.assertEqual(
            copied_aeroelastic_wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            aeroelastic_wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
        )
        self.assertIs(
            copied_aeroelastic_wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz[
                0
            ],
            deriv_func,
        )


if __name__ == "__main__":
    unittest.main()
