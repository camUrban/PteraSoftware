"""Benchmarks correctness and speedup of the parallel Biot-Savart line vortex kernel.

Compares the parallel Numba kernel against the serial kernel to verify matching
numerical results and measure the parallel speedup across a range of problem sizes.
"""

import sys
import timeit

import numba
import numpy as np

# noinspection PyProtectedMember
from pterasoftware import _aerodynamics_functions

# noinspection PyProtectedMember
_serial = _aerodynamics_functions._collapsed_velocities_from_line_vortices

# noinspection PyProtectedMember
_parallel = _aerodynamics_functions._collapsed_velocities_from_line_vortices_parallel

_SIZES = [
    ("PD Small", 100, 50),
    ("PD Medium", 200, 100),
    ("PD Large", 500, 200),
    ("PD XLarge", 1000, 500),
    ("PD Huge", 2000, 1000),
    ("VD Small", 50, 100),
    ("VD Medium", 100, 200),
    ("VD Large", 200, 500),
    ("VD XLarge", 500, 1000),
    ("VD Huge", 1000, 2000),
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


def bench_correctness() -> bool:
    """Verifies that the parallel and serial kernels produce matching results.

    :return: True if every problem size passes within tolerance; False otherwise.
    """
    print("=" * 70)
    print("CORRECTNESS")
    print("=" * 70)
    print()

    print(
        f"{'Size':<10} {'Points':<8} {'Vortices':<10} "
        f"{'Max Error':>11} {'Results':>9} {'Counts':>8} {'Status':>7}"
    )
    print("-" * 70)

    all_pass = True

    for name, num_points, num_vortices in _SIZES:
        inputs = _make_inputs(num_points, num_vortices)

        counts_s = np.zeros(4, dtype=np.int64)
        result_s = _serial(*inputs, counts_s)

        counts_p = np.zeros(4, dtype=np.int64)
        result_p = _parallel(*inputs, counts_p)

        max_error = float(np.max(np.abs(result_s - result_p)))

        # Parallelizing should preserve sum order, so results should be bit identical.
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
    print('*PD = "Point-Dominated", VD = "Vortex-Dominated"')
    print()
    return all_pass


def bench_performance() -> float:
    """Measures performance and speedup across problem sizes.

    :return: The geometric mean speedup (serial / parallel) across all sizes.
    """
    print("=" * 70)
    print("PERFORMANCE")
    print("=" * 70)
    print()

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
        _serial(*inputs, np.zeros(4, dtype=np.int64))
        _parallel(*inputs, np.zeros(4, dtype=np.int64))

        # Use autorange to pick a per-sample call count large enough that clock noise
        # is a negligible fraction, then take the minimum of 5 samples because timing
        # noise (OS jitter, GC, context switches) can only inflate a sample, never
        # shorten it.
        counts_s = np.zeros(4, dtype=np.int64)
        timer_s = timeit.Timer(lambda: _serial(*inputs, counts_s))
        number_s, _ = timer_s.autorange()
        time_s = min(timer_s.repeat(repeat=5, number=number_s)) / number_s

        counts_p = np.zeros(4, dtype=np.int64)
        timer_p = timeit.Timer(lambda: _parallel(*inputs, counts_p))
        number_p, _ = timer_p.autorange()
        time_p = min(timer_p.repeat(repeat=5, number=number_p)) / number_p

        speedup = time_s / time_p if time_p > 0 else 0.0
        speedups.append(speedup)

        print(
            f"{name:<10} {num_points:<8} {num_vortices:<10} "
            f"{time_s * 1000:>12.2f} {time_p * 1000:>14.2f} {speedup:>9.2f}x"
        )

    print()
    print('*PD = "Point-Dominated", VD = "Vortex-Dominated"')

    # Use the geometric mean because arithmetic mean of ratios overweights
    # high speedups: 0.5x and 2.0x should average to 1.0x (break-even), not
    # 1.25x.
    gmean_speedup = float(np.exp(np.mean(np.log(speedups))))
    max_speedup = float(np.max(speedups))

    print()
    print(f"Geometric mean speedup: {gmean_speedup:.2f}x")
    print(f"Maximum speedup:        {max_speedup:.2f}x")
    print()

    return gmean_speedup


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
    print(f"Correctness:            {'Pass' if correct_summary else 'Fail'}")
    print(f"Geometric mean speedup: {gmean_speedup_summary:.2f}x")
