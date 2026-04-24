"""Benchmarks correctness and speedup of the parallel Biot-Savart line vortex kernels.

Compares each parallel Numba kernel against a serial reference to verify matching
numerical results and measure the parallel speedup across a range of problem sizes.
Covers both the collapsed kernel (output shape (N,3)) and the expanded kernel (output
shape (N,M,3)).

The expanded kernel allocates an output array of shape (num_points, num_vortices, 3),
which can grow to several gigabytes for the larger _SIZES entries. Comment out the
largest entries if memory becomes a constraint.
"""

import math
import sys
import timeit

import numba
import numpy as np
from numba import njit

# noinspection PyProtectedMember
from pterasoftware import _aerodynamics_functions

# noinspection PyProtectedMember
from pterasoftware._aerodynamics_functions import (
    _eps,
    _four_lamb,
    _four_pi,
    _squire,
    _tol,
)

# noinspection PyProtectedMember
_collapsed_parallel = _aerodynamics_functions._collapsed_velocities_from_line_vortices

# noinspection PyProtectedMember
_expanded_parallel = _aerodynamics_functions._expanded_velocities_from_line_vortices


@njit(cache=True, fastmath=False)
def _collapsed_serial_reference(
    stackP_GP1_CgP1: np.ndarray,
    stackSlvp_GP1_CgP1: np.ndarray,
    stackElvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of LineVortices and
    finds the cumulative induced velocity at every point.

    This is a local copy of the package's original serial collapsed Biot-Savart kernel,
    kept in the benchmark so the parallel kernel can be timed against a same algorithm
    serial baseline. The body must be kept in sync with the package kernel if the
    algorithm changes. The correctness check in this benchmark will catch any numerical
    drift immediately because it compares against the parallel kernel for bit equality.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackSlvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M LineVortices' starting vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param stackElvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M LineVortices' ending vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        LineVortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        LineVortices. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across calls.
    :param ages: For bound LineVortices, this must be None. For LineVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M LineVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,3) ndarray of floats for the cumulative induced velocity at each of
        the N points (in the first Airplane's geometry axes, observed from the Earth
        frame). The units are in meters per second.
    """
    num_vortices = stackSlvp_GP1_CgP1.shape[0]
    num_points = stackP_GP1_CgP1.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in
    # the first Airplane's geometry axes, observed from the Earth frame).
    stackVInd_GP1__E = np.zeros((num_points, 3))

    # If the user didn't specify any ages, set the age of each LineVortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
        Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]

        # The r0_GP1 vector goes from the LineVortex's start point to its end point (in
        # the first Airplane's geometry axes).
        r0X_GP1 = Elvp_GP1_CgP1[0] - Slvp_GP1_CgP1[0]
        r0Y_GP1 = Elvp_GP1_CgP1[1] - Slvp_GP1_CgP1[1]
        r0Z_GP1 = Elvp_GP1_CgP1[2] - Slvp_GP1_CgP1[2]

        # Find r0_GP1's length.
        r0 = math.sqrt(r0X_GP1**2.0 + r0Y_GP1**2.0 + r0Z_GP1**2.0)

        # Skip degenerate filaments where the start and end points coincide.
        if r0 < _eps:
            singularity_counts[0] += 1
            continue

        strength = strengths[vortex_id]
        age = ages[vortex_id]
        r_c0 = r_c0s[vortex_id]

        # Pre compute r0 * _tol outside the inner loop.
        r0_times_tol = r0 * _tol

        # Calculate the radius of the LineVortex's core squared. The initial core radius
        # ensures nonzero regularization even for bound vortices with zero age.
        r_c_sq = r_c0**2.0 + _four_lamb * (nu + _squire * abs(strength)) * age

        c_1 = strength / _four_pi
        c_2 = r0**2.0 * r_c_sq

        for point_id in range(num_points):
            P_GP1_CgP1 = stackP_GP1_CgP1[point_id]

            # The r1_GP1 vector goes from P_GP1_CgP1 to the LineVortex's start point (in
            # the first Airplane's geometry axes).
            r1X_GP1 = Slvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r1Y_GP1 = Slvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r1Z_GP1 = Slvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r2_GP1 vector goes from P_GP1_CgP1 to the LineVortex's end point (in
            # the first Airplane's geometry axes).
            r2X_GP1 = Elvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r2Y_GP1 = Elvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r2Z_GP1 = Elvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r3_GP1 vector is the cross product of r1_GP1 and r2_GP1 (in the first
            # Airplane's geometry axes).
            r3X_GP1 = r1Y_GP1 * r2Z_GP1 - r1Z_GP1 * r2Y_GP1
            r3Y_GP1 = r1Z_GP1 * r2X_GP1 - r1X_GP1 * r2Z_GP1
            r3Z_GP1 = r1X_GP1 * r2Y_GP1 - r1Y_GP1 * r2X_GP1

            # Find the lengths of r1_GP1 and r2_GP1.
            r1 = math.sqrt(r1X_GP1**2.0 + r1Y_GP1**2.0 + r1Z_GP1**2.0)
            r2 = math.sqrt(r2X_GP1**2.0 + r2Y_GP1**2.0 + r2Z_GP1**2.0)

            # Check for singularities using scale invariant criteria. The vertex
            # proximity checks (r1/r0 and r2/r0 but refactored below to use
            # multiplication instead of slower division) guard 1/r singularities.
            if r1 < r0_times_tol:
                singularity_counts[1] += 1
                continue
            if r2 < r0_times_tol:
                singularity_counts[2] += 1
                continue

            # Cache squared length of r3_GP1 as it is used in the c_4 calculation.
            r3_sq = r3X_GP1**2.0 + r3Y_GP1**2.0 + r3Z_GP1**2.0

            # Find the length of r3_GP1.
            r3 = math.sqrt(r3_sq)

            # Cache r1 * r2 as it is used for the collinearity check and twice in the
            # c_4 calculation.
            r1_times_r2 = r1 * r2

            c_3 = r1X_GP1 * r2X_GP1 + r1Y_GP1 * r2Y_GP1 + r1Z_GP1 * r2Z_GP1

            # The collinearity check (r3/(r1*r2) = |sin(theta)| but with the same
            # multiplication instead of division refactor) guards catastrophic
            # cancellation in 1-cos(theta).
            if r3 < (_tol * r1_times_r2):
                # Collinearity can indicate one of two things. If the point is collinear
                # and between the filament's vertices, it is a true singularity (the
                # Biot-Savart equation diverges), so we exclude the contribution as it
                # is the influence of the filament on itself. If the point is collinear
                # and off to one side of the filament, it isn't a true singularity, as
                # the Biot-Savart equation (if calculated with infinite precision)
                # correctly returns zero induced velocity. However, we still run into
                # the catastrophic cancellation issue, so we again manually return zero
                # induced velocity contribution. These two situations are distinguished
                # by the sign of the c_3 (the dot product of r1 and r2).
                if c_3 < 0.0:
                    singularity_counts[3] += 1
                continue

            c_4 = c_1 * (r1 + r2) * (r1_times_r2 - c_3) / (r1_times_r2 * (r3_sq + c_2))
            stackVInd_GP1__E[point_id, 0] += c_4 * r3X_GP1
            stackVInd_GP1__E[point_id, 1] += c_4 * r3Y_GP1
            stackVInd_GP1__E[point_id, 2] += c_4 * r3Z_GP1
    return stackVInd_GP1__E


_collapsed_serial = _collapsed_serial_reference


@njit(cache=True, fastmath=False)
def _expanded_serial_reference(
    stackP_GP1_CgP1: np.ndarray,
    stackSlvp_GP1_CgP1: np.ndarray,
    stackElvp_GP1_CgP1: np.ndarray,
    strengths: np.ndarray,
    r_c0s: np.ndarray,
    singularity_counts: np.ndarray,
    ages: np.ndarray | None = None,
    nu: float = 0.0,
) -> np.ndarray:
    """Takes in a group of points and the attributes of a group of LineVortices and
    finds the induced velocity at every point due to each LineVortex.

    This is a local copy of the package's original serial expanded Biot-Savart kernel,
    kept in the benchmark so the parallel kernel can be timed against a same algorithm
    serial baseline. The body must be kept in sync with the package kernel if the
    algorithm changes. The correctness check in this benchmark will catch any numerical
    drift immediately because it compares against the parallel kernel for bit equality.

    :param stackP_GP1_CgP1: A (N,3) ndarray of floats representing the positions of N
        points (in the first Airplane's geometry axes, relative to the first Airplane's
        CG). The units are in meters.
    :param stackSlvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of M
        LineVortices' starting vertices (in the first Airplane's geometry axes, relative
        to the first Airplane's CG). The units are in meters.
    :param stackElvp_GP1_CgP1: A (M,3) ndarray of floats representing the positions of
        the M LineVortices' ending vertices (in the first Airplane's geometry axes,
        relative to the first Airplane's CG). The units are in meters.
    :param strengths: A (M,) ndarray of floats representing the strengths of the M
        LineVortices. The units are in meters squared per second.
    :param r_c0s: A (M,) ndarray of floats representing the initial core radii of the M
        LineVortices. The units are in meters.
    :param singularity_counts: A (4,) ndarray of int64 representing the cumulative
        counts of singularity events. Index mapping: [0] degenerate filament, [1] vertex
        start proximity, [2] vertex end proximity, [3] collinearity. Counts are
        incremented in place and accumulate across calls.
    :param ages: For bound LineVortices, this must be None. For LineVortices that have
        been shed into the wake, it must be a (M,) ndarray of floats representing the
        ages of the M LineVortices in seconds. The default is None.
    :param nu: A non negative float representing the kinematic viscosity of the fluid.
        The units are in meters squared per second. The default is 0.0.
    :return: A (N,M,3) ndarray of floats for the induced velocity at each of the N
        points (in the first Airplane's geometry axes, observed from the Earth frame)
        due to each of the M LineVortices. The units are in meters per second.
    """
    num_vortices = stackSlvp_GP1_CgP1.shape[0]
    num_points = stackP_GP1_CgP1.shape[0]

    # Initialize an empty array, which we will fill with the induced velocities (in the
    # first Airplane's geometry axes, observed from the Earth frame).
    gridVInd_GP1__E = np.zeros((num_points, num_vortices, 3))

    # If the user didn't specify any ages, set the age of each LineVortex to 0.0
    # seconds.
    if ages is None:
        ages = np.zeros(num_vortices)

    for vortex_id in range(num_vortices):
        Slvp_GP1_CgP1 = stackSlvp_GP1_CgP1[vortex_id]
        Elvp_GP1_CgP1 = stackElvp_GP1_CgP1[vortex_id]

        # The r0_GP1 vector goes from the LineVortex's start point to its end point (in
        # the first Airplane's geometry axes).
        r0X_GP1 = Elvp_GP1_CgP1[0] - Slvp_GP1_CgP1[0]
        r0Y_GP1 = Elvp_GP1_CgP1[1] - Slvp_GP1_CgP1[1]
        r0Z_GP1 = Elvp_GP1_CgP1[2] - Slvp_GP1_CgP1[2]

        # Find r0_GP1's length.
        r0 = math.sqrt(r0X_GP1**2.0 + r0Y_GP1**2.0 + r0Z_GP1**2.0)

        # Skip degenerate filaments where the start and end points coincide.
        if r0 < _eps:
            singularity_counts[0] += 1
            continue

        strength = strengths[vortex_id]
        age = ages[vortex_id]
        r_c0 = r_c0s[vortex_id]

        # Pre compute r0 * _tol outside the inner loop.
        r0_times_tol = r0 * _tol

        # Calculate the radius of the LineVortex's core squared. The initial core radius
        # ensures nonzero regularization even for bound vortices with zero age.
        r_c_sq = r_c0**2.0 + _four_lamb * (nu + _squire * abs(strength)) * age

        c_1 = strength / _four_pi
        c_2 = r0**2.0 * r_c_sq

        for point_id in range(num_points):
            P_GP1_CgP1 = stackP_GP1_CgP1[point_id]

            # The r1_GP1 vector goes from P_GP1_CgP1 to the LineVortex's start point (in
            # the first Airplane's geometry axes).
            r1X_GP1 = Slvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r1Y_GP1 = Slvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r1Z_GP1 = Slvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r2_GP1 vector goes from P_GP1_CgP1 to the LineVortex's end point (in
            # the first Airplane's geometry axes).
            r2X_GP1 = Elvp_GP1_CgP1[0] - P_GP1_CgP1[0]
            r2Y_GP1 = Elvp_GP1_CgP1[1] - P_GP1_CgP1[1]
            r2Z_GP1 = Elvp_GP1_CgP1[2] - P_GP1_CgP1[2]

            # The r3_GP1 vector is the cross product of r1_GP1 and r2_GP1 (in the first
            # Airplane's geometry axes).
            r3X_GP1 = r1Y_GP1 * r2Z_GP1 - r1Z_GP1 * r2Y_GP1
            r3Y_GP1 = r1Z_GP1 * r2X_GP1 - r1X_GP1 * r2Z_GP1
            r3Z_GP1 = r1X_GP1 * r2Y_GP1 - r1Y_GP1 * r2X_GP1

            # Find the lengths of r1_GP1 and r2_GP1.
            r1 = math.sqrt(r1X_GP1**2.0 + r1Y_GP1**2.0 + r1Z_GP1**2.0)
            r2 = math.sqrt(r2X_GP1**2.0 + r2Y_GP1**2.0 + r2Z_GP1**2.0)

            # Check for singularities using scale invariant criteria. The vertex
            # proximity checks (r1/r0 and r2/r0 but refactored below to use
            # multiplication instead of slower division) guard 1/r singularities.
            if r1 < r0_times_tol:
                singularity_counts[1] += 1
                continue
            if r2 < r0_times_tol:
                singularity_counts[2] += 1
                continue

            # Cache squared length of r3_GP1 as it is used in the c_4 calculation.
            r3_sq = r3X_GP1**2.0 + r3Y_GP1**2.0 + r3Z_GP1**2.0

            # Find the length of r3_GP1.
            r3 = math.sqrt(r3_sq)

            # Cache r1 * r2 as it is used for the collinearity check and twice in the
            # c_4 calculation.
            r1_times_r2 = r1 * r2

            c_3 = r1X_GP1 * r2X_GP1 + r1Y_GP1 * r2Y_GP1 + r1Z_GP1 * r2Z_GP1

            # The collinearity check (r3/(r1*r2) = |sin(theta)| but with the same
            # multiplication instead of division refactor) guards catastrophic
            # cancellation in 1-cos(theta).
            if r3 < (_tol * r1_times_r2):
                # Collinearity can indicate one of two things. If the point is collinear
                # and between the filament's vertices, it is a true singularity (the
                # Biot-Savart equation diverges), so we exclude the contribution as it
                # is the influence of the filament on itself. If the point is collinear
                # and off to one side of the filament, it isn't a true singularity, as
                # the Biot-Savart equation (if calculated with infinite precision)
                # correctly returns zero induced velocity. However, we still run into
                # the catastrophic cancellation issue, so we again manually return zero
                # induced velocity contribution. These two situations are distinguished
                # by the sign of the c_3 (the dot product of r1 and r2).
                if c_3 < 0.0:
                    singularity_counts[3] += 1
                continue

            c_4 = c_1 * (r1 + r2) * (r1_times_r2 - c_3) / (r1_times_r2 * (r3_sq + c_2))
            gridVInd_GP1__E[point_id, vortex_id, 0] = c_4 * r3X_GP1
            gridVInd_GP1__E[point_id, vortex_id, 1] = c_4 * r3Y_GP1
            gridVInd_GP1__E[point_id, vortex_id, 2] = c_4 * r3Z_GP1
    return gridVInd_GP1__E


_expanded_serial = _expanded_serial_reference

_SIZES = [
    ("PD Small", 1000, 500),
    ("PD Medium", 2000, 1000),
    ("PD Large", 5000, 2000),
    ("PD XLarge", 10000, 5000),
    ("PD Huge", 20000, 10000),
    ("VD Small", 500, 1000),
    ("VD Medium", 1000, 2000),
    ("VD Large", 2000, 5000),
    ("VD XLarge", 5000, 10000),
    ("VD Huge", 10000, 20000),
]


def _make_inputs(num_points: int, num_vortices: int) -> tuple:
    """Creates a deterministic set of inputs for a given problem size."""
    np.random.seed(42)
    stackP_GP1_CgP1 = np.random.randn(num_points, 3).astype(np.float64)
    stackSlvp_GP1_CgP1 = np.random.randn(num_vortices, 3).astype(np.float64)
    stackElvp_GP1_CgP1 = (
        stackSlvp_GP1_CgP1 + np.random.randn(num_vortices, 3).astype(np.float64) * 0.5
    )
    strengths = np.random.randn(num_vortices).astype(np.float64)
    r_c0s = np.abs(np.random.randn(num_vortices).astype(np.float64)) + 0.01
    return (
        stackP_GP1_CgP1,
        stackSlvp_GP1_CgP1,
        stackElvp_GP1_CgP1,
        strengths,
        r_c0s,
    )


def _correctness_table(label: str, serial, parallel) -> bool:
    """Verifies that one (serial, parallel) kernel pair produces matching results.

    :param label: A short label printed as a sub header above the table (for example,
        "Collapsed" or "Expanded").
    :param serial: The serial kernel to use as the baseline.
    :param parallel: The parallel kernel to compare against the baseline.
    :return: True if every problem size passes within tolerance; False otherwise.
    """
    print(label)
    print("-" * 70)
    print(
        f"{'Size':<10} {'Points':<8} {'Vortices':<10} "
        f"{'Max Error':>11} {'Results':>9} {'Counts':>8} {'Status':>7}"
    )
    print("-" * 70)

    all_pass = True

    for name, num_points, num_vortices in _SIZES:
        inputs = _make_inputs(num_points, num_vortices)

        counts_s = np.zeros(4, dtype=np.int64)
        result_s = serial(*inputs, counts_s)

        counts_p = np.zeros(4, dtype=np.int64)
        result_p = parallel(*inputs, counts_p)

        max_error = float(np.max(np.abs(result_s - result_p)))

        # Each output element is written by exactly one thread (the collapsed kernel
        # accumulates into a per point row; the expanded kernel writes a single value
        # per (point, vortex) cell), so results should be bit identical to the serial
        # baseline.
        results_match = np.array_equal(result_s, result_p)
        counts_match = np.array_equal(counts_s, counts_p)

        status = "Pass" if (results_match and counts_match) else "Fail"
        print(
            f"{name:<10} {num_points:<8} {num_vortices:<10} "
            f"{max_error:>11.2e} {str(results_match):>9} {str(counts_match):>8} "
            f"{status:>7}"
        )

        if not (results_match and counts_match):
            all_pass = False

    print()
    return all_pass


def bench_correctness() -> dict[str, bool]:
    """Verifies that every parallel kernel matches its serial reference.

    :return: A dict mapping each kernel label to True if it passes for every problem
        size and False otherwise.
    """
    print("=" * 70)
    print("CORRECTNESS")
    print("=" * 70)
    print()

    results = {
        "Collapsed": _correctness_table(
            "Collapsed", _collapsed_serial, _collapsed_parallel
        ),
        "Expanded": _correctness_table(
            "Expanded", _expanded_serial, _expanded_parallel
        ),
    }

    print('*PD = "Point-Dominated", VD = "Vortex-Dominated"')
    print()
    return results


def _performance_table(label: str, serial, parallel) -> float:
    """Measures performance and speedup for one (serial, parallel) kernel pair.

    :param label: A short label printed as a sub header above the table (for example,
        "Collapsed" or "Expanded").
    :param serial: The serial kernel to time as the baseline.
    :param parallel: The parallel kernel to time and compare against the baseline.
    :return: The geometric mean speedup (serial / parallel) across all sizes.
    """
    print(label)
    print("-" * 70)
    print(
        f"{'Size':<10} {'Points':<8} {'Vortices':<10} "
        f"{'Serial (ms)':>12} {'Parallel (ms)':>14} {'Speedup':>10}"
    )
    print("-" * 70)

    speedups = []

    for name, num_points, num_vortices in _SIZES:
        inputs = _make_inputs(num_points, num_vortices)

        # Warm up each kernel once so JIT compilation is not attributed to the first
        # timed call inside autorange.
        serial(*inputs, np.zeros(4, dtype=np.int64))
        parallel(*inputs, np.zeros(4, dtype=np.int64))

        # Use autorange to pick a per-sample call count large enough that clock noise
        # is a negligible fraction, then take the minimum of 5 samples because timing
        # noise (OS jitter, GC, context switches) can only inflate a sample, never
        # shorten it.
        counts_s = np.zeros(4, dtype=np.int64)
        timer_s = timeit.Timer(lambda: serial(*inputs, counts_s))
        number_s, _ = timer_s.autorange()
        time_s = min(timer_s.repeat(repeat=5, number=number_s)) / number_s

        counts_p = np.zeros(4, dtype=np.int64)
        timer_p = timeit.Timer(lambda: parallel(*inputs, counts_p))
        number_p, _ = timer_p.autorange()
        time_p = min(timer_p.repeat(repeat=5, number=number_p)) / number_p

        speedup = time_s / time_p if time_p > 0 else 0.0
        speedups.append(speedup)

        print(
            f"{name:<10} {num_points:<8} {num_vortices:<10} "
            f"{time_s * 1000:>12.2f} {time_p * 1000:>14.2f} {speedup:>9.2f}x"
        )

    # Use the geometric mean because arithmetic mean of ratios overweights
    # high speedups: 0.5x and 2.0x should average to 1.0x (break-even), not
    # 1.25x.
    gmean_speedup = float(np.exp(np.mean(np.log(speedups))))
    max_speedup = float(np.max(speedups))

    print()
    print(f"  Geometric mean speedup: {gmean_speedup:.2f}x")
    print(f"  Maximum speedup:        {max_speedup:.2f}x")
    print()

    return gmean_speedup


def bench_performance() -> dict[str, float]:
    """Measures performance and speedup for every (serial, parallel) kernel pair.

    :return: A dict mapping each kernel label to its geometric mean speedup across all
        sizes.
    """
    print("=" * 70)
    print("PERFORMANCE")
    print("=" * 70)
    print()

    speedups = {
        "Collapsed": _performance_table(
            "Collapsed", _collapsed_serial, _collapsed_parallel
        ),
        "Expanded": _performance_table(
            "Expanded", _expanded_serial, _expanded_parallel
        ),
    }

    print('*PD = "Point-Dominated", VD = "Vortex-Dominated"')
    print()
    return speedups


if __name__ == "__main__":
    print()
    print("Parallel Biot-Savart Kernel Benchmark")
    print(f"Python {sys.version}")
    print(f"Numba threads: {numba.get_num_threads()}")
    print(f"Numba threading layer: {numba.threading_layer()}")
    print()

    correct_summary = bench_correctness()
    gmean_speedup_summary = bench_performance()

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for label in correct_summary:
        correct_str = "Pass" if correct_summary[label] else "Fail"
        print(f"{label} correctness:            {correct_str}")
        print(f"{label} geometric mean speedup: {gmean_speedup_summary[label]:.2f}x")
