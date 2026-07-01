"""This module contains classes to test the convergence analysis functions."""

import unittest

import numpy as np

from pterasoftware import convergence


class TestConvergedParameterId(unittest.TestCase):
    """This class contains methods for testing
    convergence._converged_parameter_id.
    """

    def test_single_returns_this_id(self) -> None:
        """Test that a single tested value returns this iteration's own index."""
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=0, single=True, converged=False
            ),
            0,
        )

    def test_single_takes_precedence_over_converged(self) -> None:
        """Test that a single tested value returns this iteration's own index even
        when the converged flag is set.
        """
        self.assertEqual(
            convergence._converged_parameter_id(this_id=3, single=True, converged=True),
            3,
        )

    def test_converged_returns_coarser_id(self) -> None:
        """Test that a converged iteration returns the incrementally coarser
        index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=4, single=False, converged=True
            ),
            3,
        )

    def test_saturated_returns_this_id(self) -> None:
        """Test that an iteration that passed without converging returns this
        iteration's own index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=5, single=False, converged=False
            ),
            5,
        )

    def test_converged_id_depends_on_this_id(self) -> None:
        """Test that the coarser index tracks this iteration's index rather than a
        fixed value.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=7, single=False, converged=True
            ),
            6,
        )

    def test_returns_int(self) -> None:
        """Test that the converged index is returned as an int."""
        result = convergence._converged_parameter_id(
            this_id=2, single=False, converged=True
        )
        self.assertIsInstance(result, int)


class TestValidateCoefficientMask(unittest.TestCase):
    """This class contains methods for testing
    convergence._validate_coefficient_mask.
    """

    def test_none_returns_all_true(self) -> None:
        """Test that None returns a (6,) mask of all True."""
        result = convergence._validate_coefficient_mask(None)
        self.assertTrue(np.array_equal(result, np.ones(6, dtype=bool)))

    def test_tuple_is_returned_as_bool_array(self) -> None:
        """Test that a valid tuple is returned as an equivalent (6,) bool ndarray."""
        result = convergence._validate_coefficient_mask(
            (True, False, True, False, True, False)
        )
        self.assertTrue(
            np.array_equal(result, np.array([1, 0, 1, 0, 1, 0], dtype=bool))
        )

    def test_non_tuple_raises_type_error(self) -> None:
        """Test that a non-tuple, non-None mask raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_coefficient_mask([True] * 6)  # type: ignore[arg-type]

    def test_wrong_length_raises_value_error(self) -> None:
        """Test that a mask without exactly six elements raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_coefficient_mask(
                (True, True, True)  # type: ignore[arg-type]
            )

    def test_non_bool_element_raises_type_error(self) -> None:
        """Test that a mask with a non-bool element raises a TypeError."""
        with self.assertRaises(TypeError):
            convergence._validate_coefficient_mask(
                (True, 1, True, True, True, True)  # type: ignore[arg-type]
            )

    def test_all_false_raises_value_error(self) -> None:
        """Test that a mask with no True element raises a ValueError."""
        with self.assertRaises(ValueError):
            convergence._validate_coefficient_mask((False,) * 6)


class TestCheckCoefficientConvergence(unittest.TestCase):
    """This class contains methods for testing
    convergence._check_coefficient_convergence.
    """

    def setUp(self) -> None:
        """Set up a full mask and the tolerances used across the tests.

        :return: None
        """
        self.mask = np.ones(6, dtype=bool)
        self.rtol = 0.05
        self.atol = 0.001

    def test_identical_coefficients_converge(self) -> None:
        """Test that identical coefficients converge with a perfect metric."""
        these = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        converged, metric, _ = convergence._check_coefficient_convergence(
            these, these.copy(), self.rtol, self.atol, self.mask
        )
        self.assertTrue(converged)
        self.assertEqual(metric, 100.0)

    def test_large_relative_change_does_not_converge(self) -> None:
        """Test that a coefficient changing by more than the relative tolerance does not
        converge and is reported as the limiting coefficient.
        """
        coarser = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        these = coarser.copy()
        these[0, 2] = 2.5
        converged, _, limiting_id = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertFalse(converged)
        self.assertEqual(limiting_id, 2)

    def test_masked_out_coefficient_is_ignored(self) -> None:
        """Test that masking out the only offending coefficient makes the check
        converge.
        """
        coarser = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        these = coarser.copy()
        these[0, 2] = 2.5
        mask = np.array([True, True, False, True, True, True], dtype=bool)
        converged, _, _ = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, mask
        )
        self.assertTrue(converged)

    def test_absolute_tolerance_floors_near_zero(self) -> None:
        """Test that a coefficient near zero converges via the absolute tolerance floor
        even though its relative change is large.
        """
        these = np.zeros((1, 6), dtype=float)
        coarser = np.zeros((1, 6), dtype=float)
        coarser[0, 0] = 0.5 * self.atol
        converged, _, _ = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertTrue(converged)

    def test_all_airplanes_must_converge(self) -> None:
        """Test that the check fails when any one Airplane has an unconverged
        coefficient.
        """
        coarser = np.array(
            [
                [1.0, 0.0, 2.0, 0.1, 0.0, 0.05],
                [1.0, 0.0, 2.0, 0.1, 0.0, 0.05],
            ]
        )
        these = coarser.copy()
        these[1, 0] = 2.0
        converged, _, limiting_id = convergence._check_coefficient_convergence(
            these, coarser, self.rtol, self.atol, self.mask
        )
        self.assertFalse(converged)
        self.assertEqual(limiting_id, 0)

    def test_returns_bool_float_int(self) -> None:
        """Test that the result is a bool, a float, and an int."""
        these = np.array([[1.0, 0.0, 2.0, 0.1, 0.0, 0.05]])
        converged, metric, limiting_id = convergence._check_coefficient_convergence(
            these, these.copy(), self.rtol, self.atol, self.mask
        )
        self.assertIsInstance(converged, bool)
        self.assertIsInstance(metric, float)
        self.assertIsInstance(limiting_id, int)
