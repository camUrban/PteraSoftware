"""Contains functions for Aitken-relaxed fixed-point sub-iteration."""

from __future__ import annotations

import numpy as np


def weighted_norm(weights: np.ndarray, vector: np.ndarray) -> float:
    """Returns the norm of a vector scaled by a diagonal weighting.

    The weighting is the diagonal of a diagonal matrix D, so this returns the Euclidean
    norm of D times the vector. It nondimensionalizes a residual or increment whose
    components carry different units before the norm collapses them to a scalar.

    :param weights: A (N,) ndarray of floats representing the diagonal of the weighting
        matrix D.
    :param vector: A (N,) ndarray of floats representing the vector to take the weighted
        norm of.
    :return: A float representing the Euclidean norm of the weighted vector.
    """
    return float(np.linalg.norm(weights * vector))


def is_converged(
    weights: np.ndarray,
    residual: np.ndarray,
    increment: np.ndarray,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> bool:
    """Returns whether a sub-iteration's residual passes the mixed convergence test.

    **Notes:**

    The test is the standard mixed local-error form: the weighted residual norm must not
    exceed a relative term scaled by the weighted increment norm plus an absolute floor.
    The weighting nondimensionalizes both norms, so the two tolerances are dimensionless
    constants. The relative term governs ordinary sub-iterations, while the absolute
    floor governs near-trim sub-iterations where the increment vanishes.

    :param weights: A (N,) ndarray of floats representing the diagonal of the weighting
        matrix D.
    :param residual: A (N,) ndarray of floats representing the sub-iteration's residual,
        the difference between the dynamics-propagated state and the trial state.
    :param increment: A (N,) ndarray of floats representing the increment of the trial
        state over the snapshot state at the start of the time step.
    :param relative_tolerance: A float representing the relative convergence tolerance,
        scaling the weighted increment norm.
    :param absolute_tolerance: A float representing the absolute convergence tolerance,
        the floor on the weighted residual norm.
    :return: A bool that is True when the weighted residual norm satisfies the mixed
        convergence bound, and False otherwise.
    """
    residual_norm = weighted_norm(weights, residual)
    increment_norm = weighted_norm(weights, increment)

    return bool(
        residual_norm <= relative_tolerance * increment_norm + absolute_tolerance
    )


def aitken_relaxation_factor(
    weights: np.ndarray,
    residual: np.ndarray,
    previous_residual: np.ndarray,
    previous_factor: float,
    initial_factor: float,
    divergence_tolerance: float,
) -> float:
    """Returns the Aitken delta-squared relaxation factor for a sub-iteration.

    **Notes:**

    The factor is the Aitken delta-squared update computed in the weighted inner
    product, from the residual and relaxation factor of the previous sub-iteration. It
    is defined for the second sub-iteration onward; the first sub-iteration uses the
    initial factor directly, so this function is not called there.

    The denominator is the squared weighted norm of the change in residual between
    consecutive sub-iterations. When that change collapses relative to the residual the
    extrapolation is untrustworthy, so the factor reverts to the initial factor whenever
    the denominator falls to or below the divergence tolerance times the squared
    weighted norm of the previous residual.

    :param weights: A (N,) ndarray of floats representing the diagonal of the weighting
        matrix D.
    :param residual: A (N,) ndarray of floats representing the current sub-iteration's
        residual.
    :param previous_residual: A (N,) ndarray of floats representing the previous sub-
        iteration's residual.
    :param previous_factor: A float representing the previous sub-iteration's relaxation
        factor.
    :param initial_factor: A float representing the initial relaxation factor, returned
        as the fallback when the denominator guard triggers.
    :param divergence_tolerance: A float representing the relative threshold below which
        the denominator is treated as collapsed and the factor reverts to the initial
        factor.
    :return: A float representing the relaxation factor to apply to the current sub-
        iteration's residual.
    """
    weighted_previous_residual = weights * previous_residual
    weighted_residual_change = weights * (residual - previous_residual)

    denominator = float(np.dot(weighted_residual_change, weighted_residual_change))

    previous_residual_squared_norm = float(
        np.dot(weighted_previous_residual, weighted_previous_residual)
    )

    # When consecutive residuals barely differ the denominator collapses and the
    # extrapolation is untrustworthy, so revert to the initial factor.
    if denominator <= divergence_tolerance * previous_residual_squared_norm:
        return initial_factor

    numerator = float(np.dot(weighted_previous_residual, weighted_residual_change))

    return -previous_factor * numerator / denominator
