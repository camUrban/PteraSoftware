"""This module contains classes to test the AeroelasticWingMovement class."""

import unittest

import pterasoftware as ps
from tests.unit.fixtures import movement_fixtures


def _make_base_wing_and_wcs_movements():
    """Build a minimal base Wing and matching AeroelasticWingCrossSectionMovements.

    :return: A 2-tuple of (base_wing, wcs_movements) suitable for constructing an
        AeroelasticWingMovement.
    """
    movement = movement_fixtures.make_basic_aeroelastic_movement_fixture()
    base_wing = movement.airplane_movements[0].wing_movements[0].base_wing
    wcs_movements = [
        ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
            base_wing_cross_section=wcs
        )
        for wcs in base_wing.wing_cross_sections
    ]
    return base_wing, wcs_movements


class TestAeroelasticWingMovementSecondDerivativeValidation(unittest.TestCase):
    """This class contains unit tests for the spacingAnglesSecondDerivative_Gs_to_Wn_ixyz
    validation logic in AeroelasticWingMovement.__init__."""

    def setUp(self):
        """Set up shared base Wing and WingCrossSectionMovements for each test."""
        self.base_wing, self.wcs_movements = _make_base_wing_and_wcs_movements()

    def _make_wing_movement(self, spacingAnglesSecondDerivative_Gs_to_Wn_ixyz):
        """Helper that constructs an AeroelasticWingMovement with the given derivative
        argument and sine spacing.

        :param spacingAnglesSecondDerivative_Gs_to_Wn_ixyz: The value to pass for the
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz parameter.
        :return: The constructed AeroelasticWingMovement.
        """
        return ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=self.base_wing,
            wing_cross_section_movements=self.wcs_movements,
            ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
        )

    def test_non_sequence_raises_value_error(self):
        """Test that passing a non-sequence value for
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz raises ValueError."""
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
        raises TypeError."""
        with self.assertRaises(TypeError):
            self._make_wing_movement(
                spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[42, None, None]
            )

    def test_valid_list_of_none_values_accepted(self):
        """Test that a valid 3-element list of None values is accepted and stored."""
        wing_movement = self._make_wing_movement(
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[None, None, None]
        )

        self.assertIsNotNone(wing_movement)
        self.assertEqual(
            wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            (None, None, None),
        )

    def test_valid_callable_element_accepted_and_stored(self):
        """Test that a valid list with a callable at index 0 is accepted, converted to
        a tuple, and stored correctly."""

        def deriv_func(t):
            return -1.0 * t

        wing_movement = self._make_wing_movement(
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=[deriv_func, None, None]
        )

        result = wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz
        self.assertIsInstance(result, tuple)
        self.assertIs(result[0], deriv_func)
        self.assertIsNone(result[1])
        self.assertIsNone(result[2])

    def test_property_returns_none_defaults_when_none_passed(self):
        """Test that the property returns (None, None, None) when
        spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=None is passed (the default)."""
        wing_movement = self._make_wing_movement(
            spacingAnglesSecondDerivative_Gs_to_Wn_ixyz=None
        )

        self.assertEqual(
            wing_movement.spacingAnglesSecondDerivative_Gs_to_Wn_ixyz,
            (None, None, None),
        )

    def test_callable_spacing_no_derivative_raises(self):
        """Test that a custom callable spacingAngles_Gs_to_Wn_ixyz component with no
        matching second derivative raises ValueError."""

        def custom_spacing(t):
            return 0.0

        with self.assertRaises(ValueError):
            ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
                base_wing=self.base_wing,
                wing_cross_section_movements=self.wcs_movements,
                ampAngles_Gs_to_Wn_ixyz=(10.0, 0.0, 0.0),
                periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
                spacingAngles_Gs_to_Wn_ixyz=(custom_spacing, "sine", "sine"),
                phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
            )
