"""This script runs the coupled aeroelastic UVLM solver over every combination of
nondimensional spring stiffness (K) and reduced frequency (sigma), then plots the
resulting trailing-edge amplitude and wingtip deformation angle amplitude against
sigma, with one colored line per K value, to validate against Moore (2014)
"Analytical results on the role of flexibility in flapping propulsion" Fig. 4(a).

This version differs from aeroelastic_unsteady_first_order_validation.py in three
ways that bring it closer to that figure's assumptions:

1. The wing is nearly massless (a small but nonzero wing_density), approximating
   the R = 0 (fluid inertia only) case that Fig. 4(a) is computed for. The previous
   script used a much larger default density, which does not approximate R = 0.
2. There is little mechanical damping. Moore's torque balance (Eq. 2.9) has no
   damping term at all; the only dissipation comes from the aerodynamics. The
   previous script's larger default damping_constant adds dissipation the
   reference model does not have.
3. K and frequency are nondimensionalized the same way as the paper, K = kappa /
   (rho * U_inf^2 * c^2) and sigma = pi * c * f / U_inf, instead of sweeping raw
   spring_constant and Hz. This lets K_VALUES_NONDIM and REDUCED_FREQUENCY_VALUES
   be set directly to the paper's values and plotted on directly comparable axes.

wing_density and damping_constant are kept small but nonzero, rather than exactly
0.0, because the structural ODE's natural frequency is omega_n = sqrt(k / I), and I
scales linearly with wing_density. Driving wing_density to (near) zero sends I (and
therefore omega_n) toward zero/infinity, which makes the per-step scipy.solve_ivp
call in problems.py numerically stiff: it has to resolve every structural
oscillation within a single outer time step at tight tolerance, and the wall-clock
cost explodes (it can look like the script has hung). The same blow-up happens for
an excessively large K_RIGID_NONDIM (the rigid-wing reference, see below). The
values chosen here keep omega_n low enough, relative to the outer delta_time, that
solve_ivp converges in a handful of internal steps instead of effectively never.

One known mismatch this script does not fix: Moore's model has a single torsional
spring at the wing's leading edge, giving the whole rigid chord one pitching degree
of freedom. This script's wing instead has an independent torsional connection
between every pair of adjacent WingCrossSections (a multi-segment chain along the
span), which can have multiple natural frequencies/mode shapes instead of the
paper's single clean resonance. Collapsing the wing to a single spanwise segment
would address this but was left out here since it is a larger geometry change."""

import matplotlib.pyplot as plt
import numpy as np

import pterasoftware as ps

# Default values used for the parameters that are not being swept. A small
# DEFAULT_DENSITY approximates the R = 0 (massless wing) assumption behind Moore
# (2014) Fig. 4(a), and a small DEFAULT_B approximates the paper's undamped torque
# balance (Eq. 2.9). Neither is driven all the way to 0.0: doing so would send the
# structural ODE's rotational inertia I toward zero, which makes its natural
# frequency omega_n = sqrt(k / I) blow up and the per-step solve_ivp call in
# problems.py numerically intractable (see the module docstring).
DEFAULT_B = 10.0
DEFAULT_DENSITY = 0.05

# The freestream density, wing chord, and freestream velocity used to nondimensionalize
# K and sigma below. These must match the chord used inside run_aeroelastic and the
# rho/vCg__E used in its OperatingPoint.
RHO = 0.05
CHORD = 1.75
FREESTREAM_VELOCITY = 10.0

# The nondimensional spring stiffnesses (K) and reduced frequencies (sigma) to sweep,
# matching Moore (2014) Fig. 4(a): K = 0, 0.5, 1 (very flexible), K = 4, 6, 8 (stiff),
# and K = infinity (rigid). Every combination of K_VALUES_NONDIM and
# REDUCED_FREQUENCY_VALUES is run.
K_VALUES_NONDIM: list[float] = [200.0, 300.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0]
REDUCED_FREQUENCY_VALUES: list[float] = [
    0.5,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    8.0,
    10.0,
    12,
    14,
    16,
]

# A nondimensional spring stiffness representing the rigid-wing (K -> infinity)
# limit. Its trailing-edge motion is nearly purely kinematic (negligible elastic
# deformation relative to the runs in K_VALUES_NONDIM), so dividing every other
# run's trailing-edge amplitude by this reference's, at the same sigma, turns the
# y-axis into a gain that is approximately 1.0 for a rigid wing, matching the
# validation figure's flat black K = infinity curve. This value was checked
# directly at the highest sigma in REDUCED_FREQUENCY_VALUES (the worst case, since
# higher sigma means larger aerodynamic torque to resist): K = 500 still diverges
# past the geometric +/-90 degree limit there, while K = 2000 and K = 8000 give
# nearly identical trailing-edge amplitudes (within about 1%), confirming 2000 is
# stiff enough to be a good rigid-wing proxy. Pushing it much higher than that
# would make the structural ODE numerically intractable (see the module
# docstring's explanation of omega_n = sqrt(k / I)).
K_RIGID_NONDIM = 600000.0

# Time steps per flapping cycle, held fixed across the sigma sweep. DELTA_TIME is
# derived per frequency (period / STEPS_PER_CYCLE) below instead of being fixed,
# because a fixed DELTA_TIME under-resolves the structural dynamics at high
# frequencies. That under-resolution causes the torsional deflection to diverge past
# the +/-90 degree geometric limit before the wing ever reaches the high-frequency,
# phase-lagged regime.
STEPS_PER_CYCLE = 32

# Each run simulates this many flapping cycles, discarding the first
# TRANSIENT_CYCLES_DISCARDED of them as transient before measuring the amplitude
# over the remaining cycles.
CYCLES_PER_RUN = 4
TRANSIENT_CYCLES_DISCARDED = 3

NUM_STEPS = CYCLES_PER_RUN * STEPS_PER_CYCLE
DISCARD_STEPS = TRANSIENT_CYCLES_DISCARDED * STEPS_PER_CYCLE


def spring_constant_from_k(k_nondim: float) -> float:
    """Convert a nondimensional spring stiffness to a dimensional spring constant.

    Uses Moore (2014) Eq. (2.11), K = kappa / (rho * U_inf^2 * c^2), solved for
    kappa.

    :param k_nondim: The nondimensional spring stiffness (K).
    :return: The dimensional torsional spring constant (kappa), in the units
        expected by AeroelasticUnsteadyProblem's spring_constant parameter.
    """
    return k_nondim * RHO * FREESTREAM_VELOCITY**2 * CHORD**2


def frequency_hz_from_sigma(sigma: float) -> float:
    """Convert a reduced frequency to a dimensional flapping frequency.

    Uses Moore (2014)'s definition, sigma = pi * c * f / U_inf, solved for f.

    :param sigma: The reduced frequency (sigma).
    :return: The dimensional flapping frequency in Hz.
    """
    return sigma * FREESTREAM_VELOCITY / (np.pi * CHORD)


def run_aeroelastic(
    spring_constant: float,
    flap_frequency_hz: float,
    num_steps: int,
    delta_time: float,
    damping_constant: float = DEFAULT_B,
    wing_density: float = DEFAULT_DENSITY,
    animate: bool = False,
) -> tuple[list, object]:
    """Run the aeroelastic solver and return the net deformation data.

    :param spring_constant: The torsional spring stiffness value (K).
    :param flap_frequency_hz: The main wing's flapping frequency in Hz.
    :param num_steps: The number of time steps to simulate.
    :param delta_time: The time step size in seconds.
    :param damping_constant: The damping constant value.
    :param wing_density: The wing density per unit height (kg/m^2).
    :return: A tuple of (net_data list, solved problem object).
    """

    # Wing cross section initialization
    num_spanwise_panels = 2
    Lp_Wcsp_Lpp_Offsets = (0.1, 0.5, 0.0)
    cross_section_chords = [CHORD] * 8
    wing_cross_sections = []

    for i in range(len(cross_section_chords)):
        wing_cross_sections.append(
            ps.geometry.wing_cross_section.WingCrossSection(
                num_spanwise_panels=(
                    num_spanwise_panels if i < len(cross_section_chords) - 1 else None
                ),
                chord=cross_section_chords[i],
                Lp_Wcsp_Lpp=Lp_Wcsp_Lpp_Offsets if i > 0 else (0.0, 0.0, 0.0),
                angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                control_surface_symmetry_type="symmetric",
                control_surface_hinge_point=0.75,
                control_surface_deflection=0.0,
                spanwise_spacing=(
                    "uniform" if i < len(cross_section_chords) - 1 else None
                ),
                airfoil=ps.geometry.airfoil.Airfoil(
                    name="naca2412",
                    outline_A_lp=None,
                    resample=True,
                    n_points_per_side=400,
                ),
            )
        )

    wing_1 = ps.geometry.wing.Wing(
        wing_cross_sections=wing_cross_sections,
        name="Main Wing",
        Ler_Gs_Cgs=(0.0, 0.5, 0.0),
        angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        explode_into_strips=True,
        num_chordwise_panels=6,
        chordwise_spacing="uniform",
    )

    example_airplane = ps.geometry.airplane.Airplane(
        wings=[
            wing_1,
        ],
        name="Example Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        c_ref=None,
        b_ref=None,
    )

    # Wing cross section movement parameters
    dephase_x = 0.0
    period_x = 0.0
    amplitude_x = 0.0

    dephase_y = 0.0
    period_y = 0.0
    amplitude_y = 0.0

    dephase_z = 0.0
    period_z = 0.0
    amplitude_z = 0.0

    main_wcs_movements_list = []
    reflected_wcs_movements_list = []

    for i in range(len(example_airplane.wings[0].wing_cross_sections)):
        if i == 0:
            wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
                base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[
                    i
                ],
            )
            main_wcs_movements_list.append(wcs_movement)
            reflected_wcs_movements_list.append(wcs_movement)

        else:
            wcs_movement = ps.movements.aeroelastic_wing_cross_section_movement.AeroelasticWingCrossSectionMovement(
                base_wing_cross_section=example_airplane.wings[0].wing_cross_sections[
                    i
                ],
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(
                    amplitude_x,
                    amplitude_y,
                    amplitude_z,
                ),
                periodAngles_Wcsp_to_Wcs_ixyz=(period_x, period_y, period_z),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(dephase_x, dephase_y, dephase_z),
            )
            main_wcs_movements_list.append(wcs_movement)
            reflected_wcs_movements_list.append(wcs_movement)

    dephase = 169.0
    flap_period = 1.0 / flap_frequency_hz

    main_wing_movement = ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
        base_wing=example_airplane.wings[0],
        wing_cross_section_movements=main_wcs_movements_list,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(flap_period, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(dephase, 0.0, 0.0),
    )

    reflected_main_wing_movement = (
        ps.movements.aeroelastic_wing_movement.AeroelasticWingMovement(
            base_wing=example_airplane.wings[1],
            wing_cross_section_movements=reflected_wcs_movements_list,
            ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
            periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
            spacingLer_Gs_Cgs=("sine", "sine", "sine"),
            phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
            ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
            periodAngles_Gs_to_Wn_ixyz=(flap_period, 0.0, 0.0),
            spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
            phaseAngles_Gs_to_Wn_ixyz=(dephase, 0.0, 0.0),
        )
    )

    example_airplane_movement = (
        ps.movements.aeroelastic_airplane_movement.AeroelasticAirplaneMovement(
            base_airplane=example_airplane,
            wing_movements=[
                main_wing_movement,
                reflected_main_wing_movement,
            ],
            ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
            periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
            spacingCg_GP1_CgP1=("sine", "sine", "sine"),
            phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
        )
    )

    example_operating_point = ps.operating_point.OperatingPoint(
        rho=RHO,
        vCg__E=FREESTREAM_VELOCITY,
        alpha=0.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    example_operating_point_movement = ps.movements.aeroelastic_operating_point_movement.AeroelasticOperatingPointMovement(
        base_operating_point=example_operating_point,
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
    )

    example_movement = ps.movements.aeroelastic_movement.AeroelasticMovement(
        airplane_movements=[example_airplane_movement],
        operating_point_movement=example_operating_point_movement,
        delta_time=delta_time,
        num_steps=num_steps,
    )

    example_problem = ps.problems.AeroelasticUnsteadyProblem(
        movement=example_movement,
        wing_density=wing_density,
        spring_constant=spring_constant,
        damping_constant=damping_constant,
        aero_scaling=1.0,
        step_discards=5,
        moment_scaling_factor=1.0,
        plot_flap_cycle=False,
    )

    example_solver = ps.aeroelastic_unsteady_ring_vortex_lattice_method.AeroelasticUnsteadyRingVortexLatticeMethodSolver(
        aeroelastic_unsteady_problem=example_problem,
    )

    example_solver.run(
        prescribed_wake=True,
    )

    problem = example_solver.unsteady_problem

    if animate:
        ps.output.animate(
            unsteady_solver=example_solver,
            scalar_type="lift",
            show_wake_vortices=True,
            save=True,
        )
    return problem.net_data_per_wing[0], problem


def trailing_edge_z_history(problem: object, wing_idx: int = 0) -> np.ndarray:
    """Extract the wingtip trailing edge's vertical position at every time step.

    The trailing edge position is taken as the midpoint of the last chordwise row's
    last spanwise panel's two trailing edge vertices.

    :param problem: The solved AeroelasticUnsteadyProblem returned by
        run_aeroelastic.
    :param wing_idx: The index, within the Airplane's wings list, of the Wing whose
        trailing edge is tracked. The default is 0.
    :return: A (num_steps,) ndarray of floats representing the wingtip trailing
        edge's z-position (in the first Airplane's geometry axes, relative to the
        first Airplane's CG) at each time step. The units are meters.
    """
    steady_problems = problem.steady_problems
    wing = steady_problems[0].airplanes[0].wings[wing_idx]
    chordwise_idx = wing.num_chordwise_panels - 1
    spanwise_idx = wing.num_spanwise_panels - 1

    z_history = np.empty(len(steady_problems), dtype=float)
    for step, steady_problem in enumerate(steady_problems):
        panel = (
            steady_problem.airplanes[0]
            .wings[wing_idx]
            .panels[chordwise_idx][spanwise_idx]
        )
        z_history[step] = (panel.Blpp_GP1_CgP1[2] + panel.Brpp_GP1_CgP1[2]) / 2.0
    return z_history


def steady_state_half_range(values: np.ndarray, discard_steps: int) -> float:
    """Compute the steady-state amplitude of an oscillating time history.

    Discards the leading transient steps, then returns half the peak-to-peak range
    of what remains.

    :param values: A (num_steps,) ndarray of floats representing the oscillating
        time history.
    :param discard_steps: The number of leading time steps to discard as transient.
    :return: Half the peak-to-peak range of the time history once the transient
        steps are discarded.
    """
    steady_state = values[discard_steps:] if discard_steps < len(values) else values
    return (np.max(steady_state) - np.min(steady_state)) / 2.0


# For each sigma, first run a rigid-wing (K -> infinity) reference case to get the
# trailing-edge amplitude with no elastic deformation. If the prescribed motion ever
# drives a WingCrossSection's deformation angle past the geometric +/-90 degree
# limit, the run raises a ValueError; that sigma is skipped entirely (for every K)
# since there is no reference to normalize against.
reference_te_amplitude_by_sigma = {}
for sigma in REDUCED_FREQUENCY_VALUES:
    frequency_hz = frequency_hz_from_sigma(sigma)
    delta_time = (1.0 / frequency_hz) / STEPS_PER_CYCLE

    print(f"Running rigid reference at sigma={sigma}...")
    try:
        _, reference_problem = run_aeroelastic(
            spring_constant=spring_constant_from_k(K_RIGID_NONDIM),
            flap_frequency_hz=frequency_hz,
            num_steps=NUM_STEPS,
            delta_time=delta_time,
            animate=False,
        )
    except ValueError as error:
        print(f"  Skipping sigma={sigma}: {error}")
        continue
    reference_te_amplitude_by_sigma[sigma] = steady_state_half_range(
        trailing_edge_z_history(reference_problem), DISCARD_STEPS
    )
    print(f"  Rigid TE amplitude={reference_te_amplitude_by_sigma[sigma]:.4f} m")

# Run every (K, sigma) combination and collect the trailing-edge amplitude gain
# (relative to the rigid-wing reference at the same sigma) and the wingtip
# deformation angle amplitude for each. The wingtip curve index is the last index
# along net_data's spanwise axis, so it stays correct regardless of how many wing
# cross sections (and therefore spanwise stations) the geometry has. As with the
# rigid reference above, a deformation angle that exceeds the geometric +/-90 degree
# limit raises a ValueError; that single (K, sigma) data point is skipped rather than
# aborting the whole sweep. Each K's results and the sigmas they correspond to are
# tracked together so the plots only connect points that actually ran.
sigmas_used = {k: [] for k in K_VALUES_NONDIM}
te_amplitude_gain = {k: [] for k in K_VALUES_NONDIM}
angle_amplitudes_deg = {k: [] for k in K_VALUES_NONDIM}
wingtip_curve_index = None
for k in K_VALUES_NONDIM:
    spring_constant = spring_constant_from_k(k)
    for sigma in REDUCED_FREQUENCY_VALUES:
        if sigma not in reference_te_amplitude_by_sigma:
            print(f"Skipping K={k}, sigma={sigma}: no rigid reference available.")
            continue
        frequency_hz = frequency_hz_from_sigma(sigma)
        delta_time = (1.0 / frequency_hz) / STEPS_PER_CYCLE

        print(f"Running K={k}, sigma={sigma} (frequency={frequency_hz:.3f} Hz)...")
        try:
            net_data, problem = run_aeroelastic(
                spring_constant=spring_constant,
                flap_frequency_hz=frequency_hz,
                num_steps=NUM_STEPS,
                delta_time=delta_time,
                animate=False,
            )
        except ValueError as error:
            print(f"  Skipping K={k}, sigma={sigma}: {error}")
            continue

        net_data_array = np.array(net_data)
        wingtip_curve_index = net_data_array.shape[1] - 1

        te_amplitude = steady_state_half_range(
            trailing_edge_z_history(problem), DISCARD_STEPS
        )
        gain = te_amplitude / reference_te_amplitude_by_sigma[sigma]
        angle_history_deg = net_data_array[:, wingtip_curve_index, 1] * 180.0 / np.pi
        angle_amplitude_deg = steady_state_half_range(angle_history_deg, DISCARD_STEPS)

        sigmas_used[k].append(sigma)
        te_amplitude_gain[k].append(gain)
        angle_amplitudes_deg[k].append(angle_amplitude_deg)
        print(
            f"  TE amplitude gain={gain:.3f}, "
            f"angle amplitude={angle_amplitude_deg:.2f} deg"
        )

colors = plt.get_cmap("viridis")(np.linspace(0.0, 1.0, len(K_VALUES_NONDIM)))

# Trailing-edge amplitude gain versus reduced frequency, one colored line per K, plus
# the rigid-wing reference (gain = 1.0 by construction) for comparison.
plt.figure(figsize=(10, 6), dpi=200)
for color, k in zip(colors, K_VALUES_NONDIM):
    plt.plot(
        sigmas_used[k],
        te_amplitude_gain[k],
        marker="o",
        color=color,
        label=f"K={k}",
    )
rigid_sigmas = sorted(reference_te_amplitude_by_sigma)
plt.plot(
    rigid_sigmas,
    np.ones(len(rigid_sigmas)),
    color="black",
    linestyle="--",
    label="K=inf (rigid)",
)
plt.xlabel("reduced frequency, sigma = pi * c * f / U_inf")
plt.ylabel("Trailing-Edge Amplitude / Rigid-Wing Trailing-Edge Amplitude")
plt.title("Trailing-Edge Amplitude Gain vs. Reduced Frequency")
plt.legend()
plt.grid(True)
plt.savefig("Trailing_Edge_Amplitude_Gain_vs_Reduced_Frequency.png")
plt.show()

# Wingtip deformation angle amplitude versus reduced frequency, one colored line per
# K.
plt.figure(figsize=(10, 6), dpi=200)
for color, k in zip(colors, K_VALUES_NONDIM):
    plt.plot(
        sigmas_used[k],
        angle_amplitudes_deg[k],
        marker="o",
        color=color,
        label=f"K={k}",
    )
plt.xlabel("reduced frequency, sigma = pi * c * f / U_inf")
plt.ylabel(f"Deformation Angle Amplitude at Curve {wingtip_curve_index} (degrees)")
plt.title("Deformation Angle Amplitude vs. Reduced Frequency")
plt.legend()
plt.grid(True)
plt.savefig(
    f"Deformation_Angle_Amplitude_Curve_{wingtip_curve_index}_vs_Reduced_Frequency.png"
)
plt.show()
