"""Verifies the closed-form wake area mismatch implementation against a reference.

The production implementation in
pterasoftware.movements.movement._compute_wake_area_mismatch evaluates the metric
analytically from panel corners, skipping solver instantiation and the geometry
initialization that the solver runs at each time step. This benchmark preserves the
original full-solver implementation as a reference and uses it to:

- Verify the production implementation matches the reference (within floating
  point tolerance) at multiple delta_time values across a Strouhal sweep.
- Verify the iterative optimizer selects the same delta_time with either
  implementation.
- Report the end-to-end wall-time speedup the production implementation delivers.

The reference body must be kept in sync with
UnsteadyRingVortexLatticeMethodSolver.initialize_step_geometry, which defines the
bound and wake RingVortex geometries the closed-form implementation reproduces.

The geometry is based on examples/unsteady_ring_vortex_lattice_method_solver_variable.py
and is reduced to a single mirrored flapping main wing (no v-tail) so the benchmark
runs in a few minutes rather than tens of minutes. The kinematics match the example:
15 deg flap amplitude, 10 m/s freestream. The Strouhal number is varied by changing
the flap period; freestream velocity, flap amplitude, span, and chord are held fixed.
"""

import copy
import math
import timeit

import numpy as np

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _vortices, problems, unsteady_ring_vortex_lattice_method
from pterasoftware.movements import airplane_movement as airplane_movement_mod
from pterasoftware.movements import (
    operating_point_movement as operating_point_movement_mod,
)

# noinspection PyProtectedMember
from pterasoftware.movements.movement import (
    Movement,
    _compute_wake_area_mismatch,
    _optimize_delta_time,
)

# Geometry constants drawn from the variable-geometry flapping example.
ROOT_CHORD = 1.75
TIP_CHORD = 1.5
SEMI_SPAN = 6.0
ROOT_Y = 0.5
NUM_SPANWISE_PANELS = 8
NUM_CHORDWISE_PANELS = 6
N_AIRFOIL_POINTS = 400

# Kinematics. Strouhal is varied by changing the flap period.
FLAP_AMP_DEG = 15.0
VELOCITY = 10.0

# Strouhal sweep. The wide range stresses the optimizer at both ends:
# low St -> long period, large num_steps_per_lcm, many evaluations,
# high St -> short period, few num_steps_per_lcm, few evaluations.
STROUHAL_VALUES = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8]

# Multiplicative offsets from the analytical seed delta_time at which to assert
# the production and reference implementations agree. Spans the same 0.5x to 2x
# bracket the iterative search uses, plus a couple of intermediate points.
DT_MULTIPLIERS = [0.5, 0.75, 1.0, 1.5, 2.0]


def _strouhal_to_period(strouhal: float) -> float:
    """Converts a Strouhal number to the flap period for the benchmark geometry.

    Uses St = f * A / U with A = peak-to-peak vertical tip excursion =
    2 * SEMI_SPAN * sin(FLAP_AMP_DEG).

    :param strouhal: The dimensionless Strouhal number.
    :return: The flap period in seconds.
    """
    peak_to_peak = 2.0 * SEMI_SPAN * math.sin(np.deg2rad(FLAP_AMP_DEG))
    frequency = strouhal * VELOCITY / peak_to_peak
    return 1.0 / frequency


def _make_movement_pieces(
    period: float,
) -> tuple[
    airplane_movement_mod.AirplaneMovement,
    operating_point_movement_mod.OperatingPointMovement,
]:
    """Builds the AirplaneMovement and OperatingPointMovement for a given flap period.

    The caller is responsible for constructing the Movement, which triggers the
    analytical optimizer.

    :param period: The flap period in seconds.
    :return: A tuple of (AirplaneMovement, OperatingPointMovement) configured for
        the benchmark geometry at the given flap period.
    """
    airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        num_spanwise_panels=NUM_SPANWISE_PANELS,
                        chord=ROOT_CHORD,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        spanwise_spacing="cosine",
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            n_points_per_side=N_AIRFOIL_POINTS,
                        ),
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        num_spanwise_panels=None,
                        chord=TIP_CHORD,
                        Lp_Wcsp_Lpp=(0.75, SEMI_SPAN, 1.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        spanwise_spacing=None,
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            n_points_per_side=N_AIRFOIL_POINTS,
                        ),
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, ROOT_Y, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=NUM_CHORDWISE_PANELS,
                chordwise_spacing="uniform",
            ),
        ],
        name="Bench Airplane",
    )

    main_root_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=airplane.wings[0].wing_cross_sections[0],
        )
    )
    main_tip_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=airplane.wings[0].wing_cross_sections[1],
        )
    )
    reflected_root_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=airplane.wings[1].wing_cross_sections[0],
        )
    )
    reflected_tip_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=airplane.wings[1].wing_cross_sections[1],
        )
    )

    main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=[main_root_movement, main_tip_movement],
        ampAngles_Gs_to_Wn_ixyz=(FLAP_AMP_DEG, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(period, 0.0, 0.0),
    )
    reflected_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[1],
        wing_cross_section_movements=[reflected_root_movement, reflected_tip_movement],
        ampAngles_Gs_to_Wn_ixyz=(FLAP_AMP_DEG, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(period, 0.0, 0.0),
    )

    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[main_wing_movement, reflected_wing_movement],
    )

    operating_point = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=VELOCITY,
        alpha=1.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )
    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point,
        )
    )

    return airplane_movement, operating_point_movement


def _reference_compute_wake_area_mismatch(
    delta_time: float,
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    num_steps: int | None = None,
) -> float:
    """Full-solver reference implementation of the wake area mismatch metric.

    Creates a temporary Problem and solver, steps through some number of time steps
    (geometry only, no aerodynamic solve), and computes the average area mismatch at
    each step. Used to verify the production closed-form implementation in
    pterasoftware.movements.movement._compute_wake_area_mismatch.

    Must be kept in sync with
    UnsteadyRingVortexLatticeMethodSolver.initialize_step_geometry: that is what
    defines the bound and wake RingVortex geometries the closed-form implementation
    must reproduce.

    :param delta_time: The delta_time value to test. It must be a positive float. Its
        units are in seconds.
    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param num_steps: The number of time steps to simulate. If None (the default),
        derived from ceil(max_period / delta_time). Callers that already know the
        intended integer step count should pass it explicitly to avoid floating point
        rounding in the inferred step count, which can otherwise tip ceil(N + epsilon)
        to N + 1 when delta_time is the result of a quantized division.
    :return: The average area mismatch. The absolute percent error between the area of
        shed wake RingVortices and the area of their parent bound RingVortices (at the
        time step where they were shed). Averaged across all time steps and all pairs
        of child and parent RingVortices. A lower value indicates better matching.
    """
    airplane_movements_copy = copy.deepcopy(airplane_movements)
    operating_point_movement_copy = copy.deepcopy(operating_point_movement)

    if num_steps is None:
        max_airplane_movement_period = 0.0
        for airplane_movement in airplane_movements_copy:
            max_airplane_movement_period = max(
                airplane_movement.max_period, max_airplane_movement_period
            )

        max_period = max(
            max_airplane_movement_period, operating_point_movement_copy.max_period
        )

        # Calculate the number of time steps to traverse the max period (or just a
        # single time step if there is no movement).
        num_steps = 1
        if max_period > 0.0:
            num_steps = math.ceil(max_period / delta_time)

    temp_movement = Movement(
        airplane_movements=airplane_movements_copy,
        operating_point_movement=operating_point_movement_copy,
        delta_time=delta_time,
        num_steps=num_steps,
    )

    temp_problem = problems.UnsteadyProblem(movement=temp_movement)
    temp_solver = (
        unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            temp_problem
        )
    )

    total_mismatch = 0.0
    num_comparisons = 0

    for step in range(num_steps):
        temp_solver.initialize_step_geometry(step)

        # At time steps after the first, compare wake first row RingVortex areas
        # (current time step) to bound trailing edge RingVortex areas (previous time
        # step).
        if step > 0:
            current_airplanes = temp_solver.steady_problems[step].airplanes
            previous_airplanes = temp_solver.steady_problems[step - 1].airplanes

            for airplane_id, airplane in enumerate(current_airplanes):
                previous_airplane = previous_airplanes[airplane_id]

                for wing_id, wing in enumerate(airplane.wings):
                    previous_wing = previous_airplane.wings[wing_id]

                    wake_ring_vortices = wing.wake_ring_vortices

                    assert wake_ring_vortices is not None

                    num_spanwise = wake_ring_vortices.shape[1]

                    previous_panels = previous_wing.panels
                    if previous_panels is None:
                        continue

                    num_chordwise_panels = previous_wing.num_chordwise_panels
                    trailing_edge_chordwise_index = num_chordwise_panels - 1

                    for spanwise_id in range(num_spanwise):
                        wake_rv: _vortices.ring_vortex.RingVortex = wake_ring_vortices[
                            0, spanwise_id
                        ]
                        wake_area = wake_rv.area

                        trailing_edge_panel = previous_panels[
                            trailing_edge_chordwise_index, spanwise_id
                        ]
                        _bound_rv = trailing_edge_panel.ring_vortex

                        assert _bound_rv is not None
                        bound_rv: _vortices.ring_vortex.RingVortex = _bound_rv

                        bound_area = bound_rv.area

                        epsilon = 1e-12
                        if abs(bound_area) > epsilon:
                            total_mismatch += abs(wake_area - bound_area) / bound_area
                            num_comparisons += 1

    if num_comparisons == 0:
        return 0.0

    return total_mismatch / num_comparisons


def _reference_optimize_delta_time_non_static(
    airplane_movements: list[airplane_movement_mod.AirplaneMovement],
    operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
    initial_delta_time: float,
    lcm_period: float,
) -> float:
    """Reference brute force optimizer for non static Movements.

    Mirrors the brute force walk that the production optimizer used before the cached
    non static evaluator landed: searches integer num_steps_per_lcm_cycle values from
    0.5x to 2x the initial estimate, calling _reference_compute_wake_area_mismatch
    once per candidate with the candidate's integer step count passed explicitly.

    The bench owns this driver directly rather than asking production to expose a kill
    switch that forces the optimizer down a slow path. That keeps production free of
    benchmark scaffolding and puts the responsibility of staging the slow vs fast
    comparison on the bench, where it belongs.

    :param airplane_movements: The AirplaneMovements defining the motion.
    :param operating_point_movement: The OperatingPointMovement.
    :param initial_delta_time: The seed estimate, typically the analytical optimizer's
        result. It must be a positive float. Its units are in seconds.
    :param lcm_period: The LCM of all motion periods. It must be a positive float. Its
        units are in seconds.
    :return: The optimized delta_time. Its units are in seconds.
    """
    initial_num_steps = lcm_period / initial_delta_time
    min_num_steps = max(1, int(initial_num_steps / 2))
    max_num_steps = int(initial_num_steps * 2) + 1

    best_num_steps = min_num_steps
    best_mismatch = float("inf")
    for num_steps in range(min_num_steps, max_num_steps + 1):
        delta_time = lcm_period / num_steps
        # Pass num_steps explicitly so the reference function does not re derive it
        # via ceil(max_period / delta_time). The optimizer knows the integer step
        # count it intends; the float round trip can otherwise tip the inferred
        # step count to num_steps + 1.
        mismatch = _reference_compute_wake_area_mismatch(
            delta_time=delta_time,
            airplane_movements=airplane_movements,
            operating_point_movement=operating_point_movement,
            num_steps=num_steps,
        )
        if mismatch < best_mismatch:
            best_mismatch = mismatch
            best_num_steps = num_steps

    return lcm_period / best_num_steps


def _relative_difference(production_value: float, reference_value: float) -> float:
    """Returns the relative difference between two values.

    Falls back to the absolute difference if reference_value is 0.

    :param production_value: The production implementation's value.
    :param reference_value: The reference implementation's value.
    :return: |production_value - reference_value| / |reference_value|, or
        |production_value - reference_value| if reference_value is 0.
    """
    if reference_value == 0.0:
        return abs(production_value - reference_value)
    return abs(production_value - reference_value) / abs(reference_value)


def _check_one(strouhal: float) -> dict[str, float]:
    """Runs the agreement check and timing comparison at one Strouhal number.

    :param strouhal: The dimensionless Strouhal number.
    :return: A dict mapping result names to floats. Includes the input Strouhal
        number, the flap period, the analytical seed delta_time, the maximum
        relative mismatch difference observed across DT_MULTIPLIERS, both
        optimizer wall times, the speedup, both optimizer delta_times, and the
        relative difference between them.
    """
    period = _strouhal_to_period(strouhal)
    airplane_movement, operating_point_movement = _make_movement_pieces(period)

    # Build Movement with delta_time=None so the analytical optimizer produces the
    # seed used for both subsequent iterative calls. Both modes start from the same
    # seed, isolating the iterative phase from any analytical phase variation.
    movement = Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        delta_time=None,
        num_cycles=1,
    )
    seed_dt = movement.delta_time
    lcm_period = movement.lcm_period

    print(
        f"--- St = {strouhal:.2f} "
        f"(period = {period:.4f} s, seed_dt = {seed_dt:.6f} s) ---"
    )

    # Mismatch agreement at fixed delta_time candidates.
    max_rel_diff = 0.0
    for multiplier in DT_MULTIPLIERS:
        dt = seed_dt * multiplier
        mismatch_production = _compute_wake_area_mismatch(
            delta_time=dt,
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
        )
        mismatch_reference = _reference_compute_wake_area_mismatch(
            delta_time=dt,
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
        )
        rel_diff = _relative_difference(mismatch_production, mismatch_reference)
        max_rel_diff = max(max_rel_diff, rel_diff)
        print(
            f"  dt = {dt:.6f} s  |  production = {mismatch_production:.10f}  |  "
            f"reference = {mismatch_reference:.10f}  |  rel_diff = {rel_diff:.2e}"
        )
    print(f"  max relative mismatch difference: {max_rel_diff:.2e}")

    # End-to-end optimizer comparison. The production optimizer runs unmodified.
    # The reference run uses the bench's own brute force driver, which calls the
    # reference (full solver) mismatch function directly per candidate. This keeps
    # the bench self contained instead of patching production internals to force
    # the slow path. Each side is timed with timeit.Timer.autorange to pick the
    # per sample call count and timer.repeat with min over 3 repeats to suppress
    # OS jitter, mirroring bench_parallel_biot_savart's pattern. The per call
    # cost is well above 0.2 s for both sides at every Strouhal, so autorange
    # picks number=1 and the total cost is roughly 4 calls per timer (one for
    # autorange's discovery probe plus three repeats).
    print("  running production optimizer (3 repeats with autorange)...", flush=True)
    production_results: list[float] = []
    production_timer = timeit.Timer(
        lambda: production_results.append(
            _optimize_delta_time(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                initial_delta_time=seed_dt,
            )
        )
    )
    production_number, _ = production_timer.autorange()
    production_time = (
        min(production_timer.repeat(repeat=3, number=production_number))
        / production_number
    )
    production_dt = production_results[-1]
    print(
        f"  production min wall: {production_time:.2f}s "
        f"(number={production_number})",
        flush=True,
    )

    print("  running reference optimizer (3 repeats with autorange)...", flush=True)
    reference_results: list[float] = []
    reference_timer = timeit.Timer(
        lambda: reference_results.append(
            _reference_optimize_delta_time_non_static(
                airplane_movements=[airplane_movement],
                operating_point_movement=operating_point_movement,
                initial_delta_time=seed_dt,
                lcm_period=lcm_period,
            )
        )
    )
    reference_number, _ = reference_timer.autorange()
    reference_time = (
        min(reference_timer.repeat(repeat=3, number=reference_number))
        / reference_number
    )
    reference_dt = reference_results[-1]
    print(
        f"  reference min wall: {reference_time:.2f}s " f"(number={reference_number})",
        flush=True,
    )

    rel_diff_optim = _relative_difference(production_dt, reference_dt)
    speedup = reference_time / production_time if production_time > 0 else float("inf")

    print(
        f"  optimizer wall:  production = {production_time:.2f}s  |  "
        f"reference = {reference_time:.2f}s  |  speedup = {speedup:.2f}x"
    )
    print(
        f"  optimizer dt:    production = {production_dt:.6f} s  |  "
        f"reference = {reference_dt:.6f} s  |  rel_diff = {rel_diff_optim:.2e}"
    )
    print()

    return {
        "strouhal": strouhal,
        "period_s": period,
        "seed_dt_s": seed_dt,
        "max_rel_diff_mismatch": max_rel_diff,
        "production_time_s": production_time,
        "reference_time_s": reference_time,
        "speedup": speedup,
        "production_dt_s": production_dt,
        "reference_dt_s": reference_dt,
        "rel_diff_optim": rel_diff_optim,
    }


def _print_summary(results: list[dict[str, float]]) -> None:
    """Prints a single summary table across the Strouhal sweep.

    :param results: A list of per Strouhal result dicts produced by _check_one.
    :return: None
    """
    print("=" * 78)
    print("Summary")
    print("=" * 78)
    header = (
        f"{'St':>5}  {'mismatch diff':>14}  {'prod wall':>10}  "
        f"{'ref wall':>10}  {'speedup':>8}  {'optim diff':>11}"
    )
    print(header)
    print("-" * len(header))
    for result in results:
        print(
            f"{result['strouhal']:>5.2f}  {result['max_rel_diff_mismatch']:>14.2e}  "
            f"{result['production_time_s']:>9.2f}s  "
            f"{result['reference_time_s']:>9.2f}s  "
            f"{result['speedup']:>7.2f}x  {result['rel_diff_optim']:>11.2e}"
        )

    total_production = sum(result["production_time_s"] for result in results)
    total_reference = sum(result["reference_time_s"] for result in results)

    # Use the geometric mean because arithmetic mean of ratios overweights high
    # speedups: 0.5x and 2.0x should average to 1.0x (break even), not 1.25x.
    speedups = [result["speedup"] for result in results]
    gmean_speedup = float(np.exp(np.mean(np.log(speedups))))
    max_speedup = float(np.max(speedups))

    worst_mismatch_diff = max(result["max_rel_diff_mismatch"] for result in results)
    worst_optim_diff = max(result["rel_diff_optim"] for result in results)

    print()
    print(
        f"Total optimizer wall:  production = {total_production:.2f}s  |  "
        f"reference = {total_reference:.2f}s"
    )
    print(f"Geometric mean speedup: {gmean_speedup:.2f}x")
    print(f"Maximum speedup:        {max_speedup:.2f}x")
    print(
        f"Worst mismatch relative difference across all checks: {worst_mismatch_diff:.2e}"
    )
    print(
        f"Worst optimizer dt relative difference across all St: {worst_optim_diff:.2e}"
    )


if __name__ == "__main__":
    ps.set_up_logging()

    print("Production vs reference wake area mismatch comparison")
    print("=" * 78)
    print(
        f"Geometry: single mirrored flapping main wing, "
        f"{NUM_SPANWISE_PANELS} x {NUM_CHORDWISE_PANELS} panels per half wing"
    )
    print(f"Kinematics: amp = {FLAP_AMP_DEG} deg, vCg__E = {VELOCITY} m/s")
    print(f"Strouhal sweep: {STROUHAL_VALUES}")
    print(f"delta_time multipliers for agreement check: {DT_MULTIPLIERS}")
    print()

    all_results: list[dict[str, float]] = []
    for st in STROUHAL_VALUES:
        all_results.append(_check_one(st))

    _print_summary(all_results)
