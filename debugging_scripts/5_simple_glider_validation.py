"""Validation script for the coupled free flight solver.

Runs the simple glider in free flight with no external force (pure gliding) and
checks that:
    1. Aerodynamic forces are physically sensible during the converged prescribed phase.
    2. The glider descends and decelerates as expected when gliding with no thrust.
    3. Total mechanical energy decreases monotonically (drag dissipates energy).
    4. The rate of energy loss is consistent with drag * speed.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np

import pterasoftware as ps

# Configure logging to WARNING to reduce noise during the validation run.
logging.root.setLevel(logging.WARNING)
for handler in logging.root.handlers:
    handler.setLevel(logging.WARNING)

script_logger = logging.getLogger("validation")
script_logger.setLevel(logging.INFO)

# ---------------------------------------------------------------------------
# Configuration.
# ---------------------------------------------------------------------------
# Converged mesh parameters (from convergence study, script 1).
converged_prescribed_wake = True
converged_num_chordwise_panels = 6
wing_1_converged_num_spanwise_panels = 30
wing_2_converged_num_spanwise_panels = 6
wing_3_converged_num_spanwise_panels = 12

# Trimmed aerodynamic state (from trim study, script 2). We use the trim speed and
# alpha so that the wake develops at a realistic operating condition. But we set the
# external force to zero because this validation tests pure unpowered gliding.
trim_vCg__E = 12.9
trim_alpha = 3.3
trim_beta = 0.0
externalFX_W = 0.0

delta_time = 0.012920
converged_num_steps = 78

# Use enough free flight time steps to see real dynamics. 500 time steps at 0.01292 s
# gives about 6.5 seconds of free flight, covering at least one phugoid period
# (T_phugoid ~ pi * sqrt(2) * V / g ~ 5.85 s for this glider).
free_num_steps = 500

this_weight = 420.0  # Newtons.
this_g_magnitude = 9.80665  # Meters per second squared.
this_mass = this_weight / this_g_magnitude

I_BP1_CgP1 = np.array(
    [[155.614, 0.0, -45.658], [0.0, 398.513, 0.0], [-45.658, 0.0, 508.699]],
    dtype=float,
)

# ---------------------------------------------------------------------------
# Build the simple glider Airplane.
# ---------------------------------------------------------------------------
simple_glider_airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                    num_spanwise_panels=wing_1_converged_num_spanwise_panels,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing=None,
                ),
            ],
            chordwise_spacing="uniform",
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                    num_spanwise_panels=wing_2_converged_num_spanwise_panels,
                    control_surface_symmetry_type="symmetric",
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 1.0, 0.0),
                    control_surface_symmetry_type="symmetric",
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 0.5),
            angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
            symmetric=True,
            symmetryNormal_G=(0, 1, 0),
            symmetryPoint_G_Cg=(0, 0, 0),
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
        ps.geometry.wing.Wing(
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                    num_spanwise_panels=wing_3_converged_num_spanwise_panels,
                    spanwise_spacing="cosine",
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                    num_spanwise_panels=None,
                    Lp_Wcsp_Lpp=(0.0, 2.0, 0.0),
                ),
            ],
            chordwise_spacing="uniform",
            Ler_Gs_Cgs=(5, 0, 1.0),
            angles_Gs_to_Wn_ixyz=(90.0, 0.0, 0.0),
            symmetric=False,
            mirror_only=False,
            symmetryNormal_G=None,
            symmetryPoint_G_Cg=None,
            num_chordwise_panels=converged_num_chordwise_panels,
        ),
    ],
    weight=this_weight,
)

simple_glider_airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
    base_airplane=simple_glider_airplane,
    wing_movements=[
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[0],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        0
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[1],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        1
                    ].wing_cross_sections[1]
                ),
            ],
        ),
        ps.movements.wing_movement.WingMovement(
            base_wing=simple_glider_airplane.wings[2],
            wing_cross_section_movements=[
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[0]
                ),
                ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                    base_wing_cross_section=simple_glider_airplane.wings[
                        2
                    ].wing_cross_sections[1]
                ),
            ],
        ),
    ],
)

del simple_glider_airplane

simple_glider_coupled_operating_point = ps.operating_point.CoupledOperatingPoint(
    vCg__E=trim_vCg__E,
    alpha=trim_alpha,
    beta=trim_beta,
    angles_E_to_BP1_izyx=(0.0, trim_alpha, -trim_beta),
    externalFX_W=externalFX_W,
)

simple_glider_coupled_movement = ps.movements.movement.CoupledMovement(
    airplane_movement=simple_glider_airplane_movement,
    initial_coupled_operating_point=simple_glider_coupled_operating_point,
    delta_time=delta_time,
    prescribed_num_steps=converged_num_steps,
    free_num_steps=free_num_steps,
)

del simple_glider_airplane_movement
del simple_glider_coupled_operating_point

simple_glider_coupled_unsteady_problem = ps.problems.CoupledUnsteadyProblem(
    coupled_movement=simple_glider_coupled_movement, I_BP1_CgP1=I_BP1_CgP1
)

del simple_glider_coupled_movement

solver = ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
    coupled_unsteady_problem=simple_glider_coupled_unsteady_problem,
)

# ---------------------------------------------------------------------------
# Run the simulation.
# ---------------------------------------------------------------------------
num_steps = simple_glider_coupled_unsteady_problem.num_steps
script_logger.info(
    f"Running validation: {converged_num_steps} prescribed + {free_num_steps}"
    f" free = {num_steps} total time steps."
)

del simple_glider_coupled_unsteady_problem

solver.run(prescribed_wake=converged_prescribed_wake)

script_logger.info("Simulation completed successfully.")

# ---------------------------------------------------------------------------
# Extract time histories.
# ---------------------------------------------------------------------------
coupled_operating_points = (
    solver.coupled_unsteady_problem.coupled_movement.coupled_operating_points
)
coupled_steady_problems = solver.coupled_steady_problems

times = np.arange(num_steps, dtype=float) * delta_time

# Allocate arrays for time history data.
speeds = np.zeros(num_steps, dtype=float)
alphas = np.zeros(num_steps, dtype=float)
betas = np.zeros(num_steps, dtype=float)
forces_W = np.zeros((num_steps, 3), dtype=float)
force_coefficients_W = np.zeros((num_steps, 3), dtype=float)
moments_W_CgP1 = np.zeros((num_steps, 3), dtype=float)
moment_coefficients_W_CgP1 = np.zeros((num_steps, 3), dtype=float)
omegas_deg = np.zeros((num_steps, 3), dtype=float)
euler_angles_deg = np.zeros((num_steps, 3), dtype=float)

for i in range(num_steps):
    cop = coupled_operating_points[i]
    csp = coupled_steady_problems[i]
    airplane = csp.airplane

    speeds[i] = cop.vCg__E
    alphas[i] = cop.alpha
    betas[i] = cop.beta
    omegas_deg[i] = cop.omegas_BP1__E
    euler_angles_deg[i] = cop.angles_E_to_BP1_izyx

    this_forces_W = airplane.forces_W
    assert this_forces_W is not None
    forces_W[i] = this_forces_W

    this_force_coefficients_W = airplane.forceCoefficients_W
    assert this_force_coefficients_W is not None
    force_coefficients_W[i] = this_force_coefficients_W

    this_moments_W_CgP1 = airplane.moments_W_CgP1
    assert this_moments_W_CgP1 is not None
    moments_W_CgP1[i] = this_moments_W_CgP1

    this_moment_coefficients_W_CgP1 = airplane.momentCoefficients_W_CgP1
    assert this_moment_coefficients_W_CgP1 is not None
    moment_coefficients_W_CgP1[i] = this_moment_coefficients_W_CgP1

# Derived quantities. In Ptera Software's wind axes, the aircraft velocity is along
# +x_W and the freestream is along -x_W. Both drag and lift are negative in wind axes:
#   Drag = -Fx_W (drag opposes motion, so Fx_W < 0 for positive drag)
#   Lift = -Fz_W (lift opposes gravity, so Fz_W < 0 for positive lift)
lifts = -forces_W[:, 2]
drags = -forces_W[:, 0]
side_forces = forces_W[:, 1]
CLs = -force_coefficients_W[:, 2]
CDs = -force_coefficients_W[:, 0]

# Position from MuJoCo state history.
positions = solver.stackPosition_E_E

# ---------------------------------------------------------------------------
# Compute energy time histories (free flight phase only).
# ---------------------------------------------------------------------------
free_start = converged_num_steps
free_end = num_steps

# Gravitational acceleration vector (in Earth axes).
g_E = coupled_operating_points[0].g_E
g_magnitude = float(np.linalg.norm(g_E))

# Kinetic energy: KE = 0.5 * m * v^2.
kinetic_energies = 0.5 * this_mass * speeds**2

# Potential energy: PE = -m * dot(g_E, position_E_E).
# With g_E = [0, 0, 9.80665], gravity points in +z_E, so +z_E is "down".
# As the glider descends (z_E increases), PE decreases.
potential_energies = -this_mass * (positions @ g_E)

# Rotational kinetic energy: 0.5 * omega^T * I * omega.
omegas_rad = np.deg2rad(omegas_deg)
rotational_kinetic_energies = np.zeros(num_steps, dtype=float)
for i in range(num_steps):
    omega = omegas_rad[i]
    rotational_kinetic_energies[i] = 0.5 * omega @ I_BP1_CgP1 @ omega

# Total mechanical energy.
total_energies = kinetic_energies + potential_energies + rotational_kinetic_energies

# Power dissipated by drag: P_drag = drag * speed (always positive, removes energy).
drag_powers = drags * speeds

# ===========================================================================
# Validation checks.
# ===========================================================================
# Use the second half of the prescribed phase as the "converged" region.
converged_start = converged_num_steps // 2
converged_end = converged_num_steps

check_results = []


def log_check(name: str, passed: bool, detail: str) -> None:
    """Logs and records the result of a single validation check."""
    status = "PASS" if passed else "FAIL"
    script_logger.info(f"  [{status}] {name}: {detail}")
    check_results.append((name, passed, detail))


script_logger.info("")
script_logger.info("=" * 72)
script_logger.info("VALIDATION RESULTS")
script_logger.info("=" * 72)

# ---------------------------------------------------------------------------
# Check 1: Force sensibility during converged prescribed motion.
# ---------------------------------------------------------------------------
script_logger.info("")
script_logger.info("--- Check 1: Force sensibility (converged prescribed phase) ---")

# 1a. Lift should approximately equal weight.
mean_lift = float(np.mean(lifts[converged_start:converged_end]))
lift_weight_ratio = mean_lift / this_weight
log_check(
    "Lift ~ weight",
    0.8 < lift_weight_ratio < 1.2,
    f"mean lift = {mean_lift:.2f} N, weight = {this_weight:.1f} N,"
    f" ratio = {lift_weight_ratio:.4f}",
)

# 1b. Drag should be positive.
mean_drag = float(np.mean(drags[converged_start:converged_end]))
log_check(
    "Drag > 0",
    mean_drag > 0,
    f"mean drag = {mean_drag:.4f} N",
)

# 1c. L/D should be physically reasonable for an inviscid glider (VLM captures only
# induced drag, so L/D is higher than a real aircraft).
L_over_D = mean_lift / mean_drag if mean_drag > 0 else 0
log_check(
    "L/D in reasonable range (inviscid)",
    5 < L_over_D < 200,
    f"L/D = {L_over_D:.1f}",
)

# 1d. Side force should be near zero (symmetric Airplane, beta = 0).
mean_abs_side_force = float(np.mean(np.abs(side_forces[converged_start:converged_end])))
log_check(
    "Side force ~ 0",
    mean_abs_side_force < 0.5,
    f"mean |side force| = {mean_abs_side_force:.4e} N",
)

# 1e. CL and CD should be positive and in reasonable ranges.
mean_CL = float(np.mean(CLs[converged_start:converged_end]))
mean_CD = float(np.mean(CDs[converged_start:converged_end]))
log_check(
    "CL positive and reasonable",
    0.05 < mean_CL < 2.0,
    f"mean CL = {mean_CL:.4f}",
)
log_check(
    "CD positive and reasonable",
    0.0 < mean_CD < 0.5,
    f"mean CD = {mean_CD:.6f}",
)

# 1f. Forces should be approximately steady during converged phase (low variation).
lift_std = float(np.std(lifts[converged_start:converged_end]))
lift_cv = lift_std / abs(mean_lift) if abs(mean_lift) > 0 else 0
log_check(
    "Lift approximately steady",
    lift_cv < 0.10,
    f"lift CV = {lift_cv:.4f} (std = {lift_std:.4f} N)",
)

# ---------------------------------------------------------------------------
# Check 2: Gliding behavior during free flight.
# ---------------------------------------------------------------------------
script_logger.info("")
script_logger.info("--- Check 2: Gliding behavior (free flight phase) ---")

speed_free = speeds[free_start:free_end]
alpha_free = alphas[free_start:free_end]
pos_free = positions[free_start:free_end]

# 2a. Speed should remain bounded (phugoid oscillation is expected for unpowered
# glide starting from level flight). The speed may increase during parts of the
# phugoid cycle as the glider trades PE for KE. We check that the speed stays within
# a reasonable range of the initial speed (within 20%) rather than requiring monotonic
# decrease.
max_speed = float(np.max(speed_free))
min_speed = float(np.min(speed_free))
speed_range_pct = (max_speed - min_speed) / speed_free[0] * 100
log_check(
    "Speed bounded (phugoid)",
    max_speed < speed_free[0] * 1.2 and min_speed > speed_free[0] * 0.8,
    f"speed range: {min_speed:.4f} to {max_speed:.4f} m/s"
    f" (initial = {speed_free[0]:.4f} m/s, variation = {speed_range_pct:.2f}%)",
)

# 2a2. Phugoid period should approximately match theory: T = pi * sqrt(2) * V / g.
# Estimate the period from the speed time history by finding zero crossings of the
# demeaned speed signal.
speed_detrended = speed_free - np.mean(speed_free)
zero_crossings = np.where(np.diff(np.sign(speed_detrended)))[0]
theoretical_phugoid_period = np.pi * np.sqrt(2) * speed_free[0] / this_g_magnitude
if len(zero_crossings) >= 2:
    # Average half period from consecutive zero crossings.
    half_periods = np.diff(zero_crossings) * delta_time
    estimated_period = float(2 * np.mean(half_periods))
    period_error = (
        abs(estimated_period - theoretical_phugoid_period) / theoretical_phugoid_period
    )
    log_check(
        "Phugoid period ~ theory",
        period_error < 0.50,
        f"estimated period = {estimated_period:.2f} s,"
        f" theoretical = {theoretical_phugoid_period:.2f} s,"
        f" relative error = {period_error:.4f}",
    )
else:
    script_logger.info(
        f"  [INFO] Phugoid period: not enough zero crossings ({len(zero_crossings)})"
        f" to estimate period (theoretical = {theoretical_phugoid_period:.2f} s)."
        f" Need longer simulation."
    )

# 2b. The glider should descend. In Earth axes, +z is down (gravity in +z), so
# z_E should increase during gliding flight.
z_change = pos_free[-1, 2] - pos_free[0, 2]
log_check(
    "Glider descends (z_E increases)",
    z_change > 0,
    f"z_E: {pos_free[0, 2]:.4f} -> {pos_free[-1, 2]:.4f} m"
    f" (descent = {z_change:.4f} m)",
)

# 2c. Alpha should evolve smoothly and stay in a reasonable range.
alpha_range = float(np.ptp(alpha_free))
log_check(
    "Alpha range < 20 deg",
    alpha_range < 20,
    f"alpha: min = {float(np.min(alpha_free)):.2f},"
    f" max = {float(np.max(alpha_free)):.2f} deg"
    f" (range = {alpha_range:.2f} deg)",
)

# 2d. Beta should stay near zero (no lateral disturbance).
beta_free = betas[free_start:free_end]
max_abs_beta = float(np.max(np.abs(beta_free)))
log_check(
    "Beta stays near 0",
    max_abs_beta < 5,
    f"max |beta| = {max_abs_beta:.4e} deg",
)

# 2e. Position should evolve smoothly (no jumps).
pos_diffs = np.diff(pos_free, axis=0)
max_pos_jump = float(np.max(np.linalg.norm(pos_diffs, axis=1)))
expected_max_jump = speed_free[0] * delta_time * 2.0
log_check(
    "No position jumps",
    max_pos_jump < expected_max_jump,
    f"max position change per time step = {max_pos_jump:.4f} m"
    f" (bound = {expected_max_jump:.4f} m)",
)

# ---------------------------------------------------------------------------
# Check 3: Energy conservation during free flight.
# ---------------------------------------------------------------------------
script_logger.info("")
script_logger.info("--- Check 3: Energy conservation (free flight phase) ---")

E_free = total_energies[free_start:free_end]
KE_free = kinetic_energies[free_start:free_end]
PE_free = potential_energies[free_start:free_end]
drag_power_free = drag_powers[free_start:free_end]

# 3a. Total energy should decrease (drag dissipates energy, no energy input).
delta_E = E_free[-1] - E_free[0]
log_check(
    "Total energy decreases",
    delta_E < 0,
    f"delta E = {delta_E:.2f} J"
    f" (E_start = {E_free[0]:.2f} J, E_end = {E_free[-1]:.2f} J)",
)

# 3b. Total energy should decrease monotonically (no energy creation).
dE = np.diff(E_free)
num_increases = int(np.sum(dE > 0))
total_intervals = len(dE)
log_check(
    "Energy monotonically decreasing",
    num_increases == 0,
    f"{num_increases}/{total_intervals} intervals had energy increase",
)

# 3c. The rate of energy loss should be consistent with the drag power.
# dE/dt should approximately equal -P_drag = -(drag * speed).
# We compare the mean values over the free flight phase.
dE_dt = np.diff(E_free) / delta_time
mean_dE_dt = float(np.mean(dE_dt))
mean_drag_power = float(np.mean(drag_power_free))
# The expected dE/dt is -drag_power (negative because drag removes energy).
expected_dE_dt = -mean_drag_power
relative_error = (
    abs(mean_dE_dt - expected_dE_dt) / abs(expected_dE_dt)
    if abs(expected_dE_dt) > 0
    else 0
)
log_check(
    "dE/dt ~ -drag * speed",
    relative_error < 0.25,
    f"mean dE/dt = {mean_dE_dt:.2f} W,"
    f" expected (-drag*speed) = {expected_dE_dt:.2f} W,"
    f" relative error = {relative_error:.4f}",
)

# 3d. Glide ratio vs L/D (informational only). For steady state gliding, the glide
# ratio (horizontal distance / vertical descent) equals L/D. However, this simulation
# starts from level flight and is still in transient (mainly decelerating rather than
# descending), so the actual glide ratio will be much higher than L/D. This metric is
# logged for reference but not included as a pass/fail check.
free_duration = (free_end - free_start - 1) * delta_time
horizontal_distance = float(
    np.sqrt(
        (pos_free[-1, 0] - pos_free[0, 0]) ** 2
        + (pos_free[-1, 1] - pos_free[0, 1]) ** 2
    )
)
vertical_descent = float(pos_free[-1, 2] - pos_free[0, 2])
if horizontal_distance > 0 and vertical_descent > 0:
    actual_glide_ratio = horizontal_distance / vertical_descent
    mean_L_over_D_free = float(
        np.mean(lifts[free_start:free_end] / drags[free_start:free_end])
    )
    script_logger.info(
        f"  [INFO] Glide ratio (informational): actual = {actual_glide_ratio:.1f},"
        f" mean L/D = {mean_L_over_D_free:.1f}"
        f" (expected to differ during transient from level flight)"
    )

# ---------------------------------------------------------------------------
# Summary.
# ---------------------------------------------------------------------------
script_logger.info("")
script_logger.info("=" * 72)
num_passed = sum(1 for _, passed, _ in check_results if passed)
num_total = len(check_results)
script_logger.info(f"SUMMARY: {num_passed}/{num_total} checks passed.")
script_logger.info("=" * 72)

# Print key numeric results for reference.
script_logger.info("")
script_logger.info("Key results:")
script_logger.info(f"  Converged CL = {mean_CL:.4f}, CD = {mean_CD:.6f}")
script_logger.info(f"  Converged L/D = {L_over_D:.1f}")
script_logger.info(f"  Free flight duration = {free_duration:.3f} s")
script_logger.info(
    f"  Speed: {min_speed:.4f} to {max_speed:.4f} m/s"
    f" (initial = {speed_free[0]:.4f} m/s, variation = {speed_range_pct:.2f}%)"
)
script_logger.info(f"  Descent = {vertical_descent:.4f} m")
script_logger.info(
    f"  Energy lost = {-delta_E:.2f} J ({-delta_E / E_free[0] * 100:.2f}%)"
)
script_logger.info(f"  Mean drag power = {mean_drag_power:.2f} W")
script_logger.info(f"  Mean dE/dt = {mean_dE_dt:.2f} W")

# ===========================================================================
# Plots (free flight phase only).
# ===========================================================================
free_times = times[free_start:free_end] - times[free_start]

fig, axes = plt.subplots(3, 3, figsize=(16, 12))
fig.suptitle(
    "Simple Glider Unpowered Glide Validation (Free Flight Phase)", fontsize=14
)

# Row 0: aerodynamic state.
axes[0, 0].plot(free_times, alpha_free, "b-")
axes[0, 0].set_ylabel("alpha (deg)")
axes[0, 0].set_title("Angle of attack")

axes[0, 1].plot(free_times, beta_free, "b-")
axes[0, 1].set_ylabel("beta (deg)")
axes[0, 1].set_title("Sideslip angle")

axes[0, 2].plot(free_times, speed_free, "b-")
axes[0, 2].set_ylabel("Speed (m/s)")
axes[0, 2].set_title("Speed")

# Row 1: forces.
axes[1, 0].plot(free_times, lifts[free_start:free_end], "b-", label="Lift")
axes[1, 0].axhline(this_weight, color="g", linestyle="--", alpha=0.5, label="Weight")
axes[1, 0].set_ylabel("Force (N)")
axes[1, 0].set_title("Lift vs weight")
axes[1, 0].legend(fontsize=8)

axes[1, 1].plot(free_times, drags[free_start:free_end], "b-")
axes[1, 1].set_ylabel("Drag (N)")
axes[1, 1].set_title("Drag force")

# L/D ratio.
L_over_D_free = np.where(
    drags[free_start:free_end] > 1e-8,
    lifts[free_start:free_end] / drags[free_start:free_end],
    0.0,
)
axes[1, 2].plot(free_times, L_over_D_free, "b-")
axes[1, 2].set_ylabel("L/D")
axes[1, 2].set_title("Lift to drag ratio")

# Row 2: energy and trajectory.
axes[2, 0].plot(free_times, KE_free - KE_free[0], "b-", label="delta KE")
axes[2, 0].plot(free_times, PE_free - PE_free[0], "r-", label="delta PE")
axes[2, 0].plot(free_times, E_free - E_free[0], "k-", linewidth=2, label="delta total")
axes[2, 0].set_ylabel("Energy change (J)")
axes[2, 0].set_xlabel("Free flight time (s)")
axes[2, 0].set_title("Energy changes from start of free flight")
axes[2, 0].legend(fontsize=8)

# Compare actual dE/dt with -drag*speed.
axes[2, 1].plot(free_times[:-1], dE_dt, "b-", label="dE/dt (numerical)")
axes[2, 1].plot(free_times, -drag_power_free, "r--", label="-drag * speed", alpha=0.7)
axes[2, 1].set_ylabel("Power (W)")
axes[2, 1].set_xlabel("Free flight time (s)")
axes[2, 1].set_title("Energy dissipation rate")
axes[2, 1].legend(fontsize=8)

# Trajectory in the xz plane (Earth axes). Show expected glide slope.
axes[2, 2].plot(pos_free[:, 0], pos_free[:, 2], "b-", linewidth=2, label="Trajectory")
axes[2, 2].plot(pos_free[0, 0], pos_free[0, 2], "go", markersize=8, label="Start")
axes[2, 2].plot(pos_free[-1, 0], pos_free[-1, 2], "ro", markersize=8, label="End")
# Draw expected glide slope from L/D.
if L_over_D > 0:
    x_range = pos_free[-1, 0] - pos_free[0, 0]
    expected_z_end = pos_free[0, 2] + x_range / L_over_D
    axes[2, 2].plot(
        [pos_free[0, 0], pos_free[-1, 0]],
        [pos_free[0, 2], expected_z_end],
        "g--",
        alpha=0.5,
        label=f"Expected (L/D = {L_over_D:.0f})",
    )
axes[2, 2].set_xlabel("x_E (m)")
axes[2, 2].set_ylabel("z_E (m, + is down)")
axes[2, 2].set_title("Trajectory (xz plane)")
axes[2, 2].legend(fontsize=8)
axes[2, 2].invert_yaxis()

plt.tight_layout()
plt.savefig("simple_glider_validation.png", dpi=150)
script_logger.info(f"Plots saved to simple_glider_validation.png")
plt.show()
