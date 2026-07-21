"""This module contains classes to test the fixed point relaxation functions."""

import unittest

import numpy as np

# noinspection PyProtectedMember
from pterasoftware import _fixed_point_relaxation


class TestWeightedNorm(unittest.TestCase):
    """This class contains methods for testing _fixed_point_relaxation.weighted_norm."""

    def test_unit_weights_equals_euclidean_norm(self) -> None:
        """Test that unit weights reduce the weighted norm to the Euclidean norm."""
        weights = np.array([1.0, 1.0, 1.0], dtype=float)
        vector = np.array([3.0, 4.0, 0.0], dtype=float)
        self.assertAlmostEqual(
            _fixed_point_relaxation.weighted_norm(weights, vector), 5.0, places=12
        )

    def test_weights_scale_each_component(self) -> None:
        """Test that the weights scale each component before the norm is taken."""
        weights = np.array([2.0, 1.0, 1.0], dtype=float)
        vector = np.array([3.0, 4.0, 0.0], dtype=float)
        # The weighted vector is [6, 4, 0], whose norm is sqrt(52).
        self.assertAlmostEqual(
            _fixed_point_relaxation.weighted_norm(weights, vector),
            np.sqrt(52.0),
            places=12,
        )

    def test_zero_vector_has_zero_norm(self) -> None:
        """Test that the zero vector has a zero weighted norm."""
        weights = np.array([2.0, 5.0, 0.5], dtype=float)
        vector = np.zeros(3, dtype=float)
        self.assertAlmostEqual(
            _fixed_point_relaxation.weighted_norm(weights, vector), 0.0, places=12
        )

    def test_returns_python_float(self) -> None:
        """Test that the weighted norm is returned as a Python float."""
        weights = np.array([1.0, 1.0], dtype=float)
        vector = np.array([1.0, 0.0], dtype=float)
        result = _fixed_point_relaxation.weighted_norm(weights, vector)
        self.assertIsInstance(result, float)


class TestIsConverged(unittest.TestCase):
    """This class contains methods for testing _fixed_point_relaxation.is_converged."""

    def test_zero_residual_is_converged(self) -> None:
        """Test that a zero residual is always converged."""
        weights = np.array([1.0, 1.0, 1.0], dtype=float)
        residual = np.zeros(3, dtype=float)
        increment = np.array([1.0, 2.0, 3.0], dtype=float)
        self.assertTrue(
            _fixed_point_relaxation.is_converged(
                weights, residual, increment, 1e-6, 1e-10
            )
        )

    def test_large_residual_is_not_converged(self) -> None:
        """Test that a large residual against a zero increment is not converged."""
        weights = np.array([1.0, 1.0, 1.0], dtype=float)
        residual = np.array([1.0, 1.0, 1.0], dtype=float)
        increment = np.zeros(3, dtype=float)
        self.assertFalse(
            _fixed_point_relaxation.is_converged(
                weights, residual, increment, 1e-6, 1e-10
            )
        )

    def test_absolute_floor_accepts_near_trim_step(self) -> None:
        """Test that the absolute floor governs when the increment vanishes.

        With a zero increment the relative term is zero, so a residual below the
        absolute tolerance is accepted and one above it is rejected.
        """
        weights = np.array([1.0, 1.0, 1.0], dtype=float)
        increment = np.zeros(3, dtype=float)
        below_floor = np.array([5e-11, 0.0, 0.0], dtype=float)
        above_floor = np.array([2e-10, 0.0, 0.0], dtype=float)
        self.assertTrue(
            _fixed_point_relaxation.is_converged(
                weights, below_floor, increment, 1e-6, 1e-10
            )
        )
        self.assertFalse(
            _fixed_point_relaxation.is_converged(
                weights, above_floor, increment, 1e-6, 1e-10
            )
        )

    def test_relative_term_governs_ordinary_step(self) -> None:
        """Test that the relative term governs when the increment is appreciable."""
        weights = np.array([1.0, 1.0, 1.0], dtype=float)
        increment = np.array([0.1, 0.0, 0.0], dtype=float)
        converged_residual = np.array([5e-8, 0.0, 0.0], dtype=float)
        diverged_residual = np.array([5e-7, 0.0, 0.0], dtype=float)
        self.assertTrue(
            _fixed_point_relaxation.is_converged(
                weights, converged_residual, increment, 1e-6, 1e-10
            )
        )
        self.assertFalse(
            _fixed_point_relaxation.is_converged(
                weights, diverged_residual, increment, 1e-6, 1e-10
            )
        )

    def test_boundary_is_inclusive(self) -> None:
        """Test that a residual exactly on the tolerance is accepted."""
        weights = np.array([1.0, 1.0, 1.0], dtype=float)
        increment = np.array([4.0, 0.0, 0.0], dtype=float)
        # With a relative tolerance of 0.5 and a zero absolute tolerance, the bound is
        # 0.5 * 4 = 2, and a residual of norm exactly 2 must be accepted.
        residual = np.array([2.0, 0.0, 0.0], dtype=float)
        self.assertTrue(
            _fixed_point_relaxation.is_converged(weights, residual, increment, 0.5, 0.0)
        )

    def test_weighting_can_change_the_verdict(self) -> None:
        """Test that the weighting changes the convergence verdict.

        The same residual and increment are unconverged under unit weights but converged
        once the dominant residual component is down weighted, confirming the weights
        enter both sides of the test.
        """
        residual = np.array([0.0, 0.0, 1.0], dtype=float)
        increment = np.array([1.0, 0.0, 0.0], dtype=float)
        unit_weights = np.array([1.0, 1.0, 1.0], dtype=float)
        down_weights = np.array([1.0, 1.0, 0.1], dtype=float)
        self.assertFalse(
            _fixed_point_relaxation.is_converged(
                unit_weights, residual, increment, 0.5, 0.0
            )
        )
        self.assertTrue(
            _fixed_point_relaxation.is_converged(
                down_weights, residual, increment, 0.5, 0.0
            )
        )

    def test_returns_python_bool(self) -> None:
        """Test that the convergence verdict is returned as a Python bool."""
        weights = np.array([1.0, 1.0], dtype=float)
        residual = np.zeros(2, dtype=float)
        increment = np.ones(2, dtype=float)
        result = _fixed_point_relaxation.is_converged(
            weights, residual, increment, 1e-6, 1e-10
        )
        self.assertIsInstance(result, bool)


class TestAitkenRelaxationFactor(unittest.TestCase):
    """This class contains methods for testing
    _fixed_point_relaxation.aitken_relaxation_factor."""

    def test_matches_hand_computed_value(self) -> None:
        """Test the factor against a hand computed unit weight case.

        With previous residual [1, 0], residual [0.5, 0], and previous factor 0.5, the
        residual change is [-0.5, 0], so the factor is -0.5 * (-0.5 / 0.25) = 1.0.
        """
        weights = np.array([1.0, 1.0], dtype=float)
        previous_residual = np.array([1.0, 0.0], dtype=float)
        residual = np.array([0.5, 0.0], dtype=float)
        factor = _fixed_point_relaxation.aitken_relaxation_factor(
            weights, residual, previous_residual, 0.5, 0.5, 1e-20
        )
        self.assertAlmostEqual(factor, 1.0, places=12)

    def test_weighting_changes_the_factor(self) -> None:
        """Test that the weighting enters the Aitken inner products.

        With weights [1, 2], previous residual [2, 1], residual [1, -1], and previous
        factor 0.5, the weighted residual change is [-1, -4] with squared norm 17 and the
        weighted inner product is -10, so the factor is -0.5 * (-10 / 17) = 5 / 17.
        """
        weights = np.array([1.0, 2.0], dtype=float)
        previous_residual = np.array([2.0, 1.0], dtype=float)
        residual = np.array([1.0, -1.0], dtype=float)
        factor = _fixed_point_relaxation.aitken_relaxation_factor(
            weights, residual, previous_residual, 0.5, 0.5, 1e-20
        )
        self.assertAlmostEqual(factor, 5.0 / 17.0, places=12)

    def test_guard_reverts_on_unchanged_residual(self) -> None:
        """Test that an unchanged residual reverts the factor to the initial factor.

        When the residual does not change between sub-iterations the denominator is
        zero, so the guard returns the initial factor rather than dividing by zero.
        """
        weights = np.array([1.0, 1.0], dtype=float)
        previous_residual = np.array([1.0, 2.0], dtype=float)
        residual = np.array([1.0, 2.0], dtype=float)
        factor = _fixed_point_relaxation.aitken_relaxation_factor(
            weights, residual, previous_residual, 0.9, 0.5, 1e-20
        )
        self.assertEqual(factor, 0.5)

    def test_guard_reverts_on_negligible_residual_change(self) -> None:
        """Test that a negligible residual change reverts to the initial factor.

        A residual change tiny relative to the residual collapses the denominator below
        the divergence tolerance, so the guard returns the initial factor instead of the
        untrustworthy extrapolation.
        """
        weights = np.array([1.0, 1.0], dtype=float)
        previous_residual = np.array([1.0, 0.0], dtype=float)
        residual = np.array([1.0 + 1e-12, 0.0], dtype=float)
        factor = _fixed_point_relaxation.aitken_relaxation_factor(
            weights, residual, previous_residual, 0.3, 0.7, 1e-20
        )
        self.assertEqual(factor, 0.7)

    def test_no_guard_when_denominator_is_healthy(self) -> None:
        """Test that a healthy denominator does not trigger the guard.

        With a residual change well above the divergence tolerance the factor is the
        Aitken value, distinct from the initial factor.
        """
        weights = np.array([1.0, 1.0], dtype=float)
        previous_residual = np.array([1.0, 0.0], dtype=float)
        residual = np.array([0.5, 0.0], dtype=float)
        factor = _fixed_point_relaxation.aitken_relaxation_factor(
            weights, residual, previous_residual, 0.5, 0.25, 1e-20
        )
        self.assertNotAlmostEqual(factor, 0.25, places=12)
        self.assertAlmostEqual(factor, 1.0, places=12)

    def test_returns_python_float(self) -> None:
        """Test that the relaxation factor is returned as a Python float."""
        weights = np.array([1.0, 1.0], dtype=float)
        previous_residual = np.array([1.0, 0.0], dtype=float)
        residual = np.array([0.5, 0.0], dtype=float)
        factor = _fixed_point_relaxation.aitken_relaxation_factor(
            weights, residual, previous_residual, 0.5, 0.5, 1e-20
        )
        self.assertIsInstance(factor, float)
