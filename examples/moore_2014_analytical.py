"""Analytical solutions from Moore (2014) for an aeroelastic flapping wing.

Reference: Moore, M. N. J. (2015). "Analytical results on the role of flexibility
in flapping propulsion." Journal of Fluid Mechanics, 757, 599-612.

Implements the exact solution for the kinematics of a thin flapping wing with a
torsional spring at the leading edge, valid for any wing-to-fluid inertia ratio
R >= 0. The paper's Fig. 3 uses R=0 (massless); this script generalises to R > 0.

Key dimensionless parameters:

    sigma = pi * f * c / U_inf            (reduced frequency, = pi*f*c/U)
    U     = 2 * U_inf / (c * f)           (dimensionless freestream, = 2*pi/sigma)
    K     = kappa / (rho * U_inf^2 * c^2) (dimensionless spring stiffness)
    R     = (b * rho_s) / (c * rho)       (wing-to-fluid inertia ratio)
    C(sigma)                               (Theodorsen function)

where kappa is the spring constant [N*m/rad], b is the wing thickness [m],
rho_s is the solid wing density [kg/m^3], rho is the fluid density [kg/m^3],
and c is the chord [m].

Mapping from Ptera Software aeroelastic parameters:

    K = spring_constant / (rho * vCg__E^2 * c^2)
    R = wing_density / (rho * c)        (wing_density = b*rho_s, [kg/m^2])
    sigma = pi * flap_frequency * c / vCg__E

Outputs (PNG files in the working directory):

    moore_2014_fig2_equiv.png  -- Wing shape snapshots (Fig. 2 equiv., R > 0)
    moore_2014_fig4a_equiv.png -- Trailing-edge amplitude vs sigma (Fig. 4a equiv.)
    moore_2014_fig36_equiv.png -- C_T and C_P vs sigma (Fig. 3/6 equiv.)
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.special import hankel2

# =============================================================================
# Part 1: Symbolic verification that D_0 and D_1 (Eqs. 2.14-2.15) follow from
# Eq. (2.10) under the time-harmonic assumption theta ~ beta_1 * e^{2*pi*j*t}.
# =============================================================================


def _derive_D_symbolically():
    """Derive D_0, D_1 symbolically from Eq. (2.10) and verify Eqs. (2.14-2.15).

    Substitutes a_0, a_1, a_2 (Eqs. 2.6-2.7) and the torque expressions
    (Eqs. 2.12-2.13) into the torque balance, then linearises in the small
    amplitude limit (theta_0 ~ beta_1). The result is a linear algebraic
    equation D_0*beta_0 + D_1*beta_1 = 0.

    The coefficient of beta_0 in the residual equals (pi^2/6)*D_0 and the
    coefficient of beta_1 equals (pi^2/6)*D_1, which matches the paper exactly.

    :return: A tuple (D0, D1, residual_D0, residual_D1) where the residuals
        should simplify to zero.
    """
    pi = sp.pi
    j = sp.I  # imaginary unit

    U, K, R = sp.symbols("U K R", positive=True)
    C = sp.Symbol("C")  # complex Theodorsen function (left symbolic)
    b0, b1 = sp.symbols("beta_0 beta_1")

    # Eqs. (2.6)-(2.7): flow coefficients in terms of kinematics
    a0 = -2 * pi * j * U * C * b0 + 2 * pi * j * U * (1 - C) * b1 - 2 * U**2 * C * b1
    a1 = 2 * pi**2 * b0 - 4 * pi * j * U * b1
    a2 = pi**2 * b1

    # Fluid torque (Eq. 2.13) and inertial torque (Eq. 2.12)
    N_f = sp.Rational(1, 2) * pi * (a0 + 2 * a1 + a2)
    N_i = 8 * pi**2 * R * (b0 - 2 * b1)

    # Eq. (2.10) with time-harmonic theta = beta_1 * exp(2*pi*j*t) so that
    # d^2 theta/dt^2 = -4*pi^2 * beta_1 (the e^{2*pi*j*t} factor cancels).
    # Rearrange to 0 = rhs_b0 * beta_0 + rhs_b1 * beta_1.
    inertia_lhs = sp.Rational(16, 3) * R * (-4 * pi**2) * b1
    rhs = sp.expand(-4 * K * U**2 * b1 + N_f + N_i - inertia_lhs)

    # Extract linear coefficients by substitution (expression is linear in b0, b1)
    coeff_b0 = rhs.subs({b0: 1, b1: 0})
    coeff_b1 = rhs.subs({b0: 0, b1: 1})

    # Paper's D_0 and D_1 are related by coeff = (pi^2 / 6) * D
    scale = pi**2 / sp.Integer(6)
    D0_paper = 48 * R + 12 * pi - 6 * j * U * C
    D1_paper = (
        32 * R
        + 3 * pi
        - 24 * K * U**2 / pi**2
        - 18 * j * U
        - 6 * j * U * C
        - 6 * U**2 * C / pi
    )

    res_D0 = sp.expand(coeff_b0 - scale * D0_paper)
    res_D1 = sp.expand(coeff_b1 - scale * D1_paper)
    return D0_paper, D1_paper, res_D0, res_D1


print("=" * 60)
print("Part 1: Symbolic verification of D_0 and D_1 (Eqs. 2.14-2.15)")
print("=" * 60)
D0_sym, D1_sym, res_D0, res_D1 = _derive_D_symbolically()
print(f"D_0 = {D0_sym}")
print(f"D_1 = {D1_sym}")
print(f"Residual D_0 (should be 0): {res_D0}")
print(f"Residual D_1 (should be 0): {res_D1}")
print()


# =============================================================================
# Part 2: Numerical evaluation functions
# =============================================================================


def theodorsen(sigma):
    """Theodorsen function C(sigma) via Hankel functions of the second kind.

    Equivalent to K_1(j*sigma) / (K_0(j*sigma) + K_1(j*sigma)) using modified
    Bessel functions, which simplifies to H1_2(sigma) / (H1_2 + i*H0_2).

    :param sigma: Reduced frequency sigma = pi*f*c/U_inf (array or scalar,
        must be positive).
    :return: Complex Theodorsen function values, same shape as sigma.
    """
    h0 = hankel2(0, sigma)
    h1 = hankel2(1, sigma)
    return h1 / (h1 + 1j * h0)


def kinematics(sigma, K, R, eps=1.0):
    """Exact solution for beta_0 and beta_1 from Eq. (2.16).

    Valid for any R >= 0. In the rigid limit K -> inf, beta_1 -> 0 and
    beta_0 -> 2*eps (pure heaving, no pitching).

    :param sigma: Reduced frequency (array or scalar, must be positive).
    :param K: Dimensionless spring stiffness K = kappa/(rho*U_inf^2*c^2).
    :param R: Wing/fluid inertia ratio R = wing_density/(rho*c).
    :param eps: Dimensionless heaving amplitude A/c. Cancels in normalised
        outputs (C_T, C_P, trailing-edge amplitude ratio). Default 1.0.
    :return: Tuple (beta_0, beta_1) of complex arrays.
    """
    U = 2.0 * np.pi / sigma
    C = theodorsen(sigma)

    D0 = 48.0 * R + 12.0 * np.pi - 6.0j * U * C
    D1 = (
        32.0 * R
        + 3.0 * np.pi
        - 24.0 * K * U**2 / np.pi**2
        - 18.0j * U
        - 6.0j * U * C
        - 6.0 * U**2 * C / np.pi
    )
    denom = 2.0 * D0 + D1
    return 2.0 * eps * D1 / denom, -2.0 * eps * D0 / denom


def _a_coeffs(sigma, beta_0, beta_1):
    """Compute flow coefficients a_0, a_1, a_2 from Eqs. (2.6)-(2.7).

    :param sigma: Reduced frequency.
    :param beta_0: Complex kinematic amplitude beta_0.
    :param beta_1: Complex kinematic amplitude beta_1.
    :return: Tuple (a_0, a_1, a_2).
    """
    U = 2.0 * np.pi / sigma
    C = theodorsen(sigma)
    a_0 = (
        -2.0 * np.pi * 1j * U * C * beta_0
        + 2.0 * np.pi * 1j * U * (1.0 - C) * beta_1
        - 2.0 * U**2 * C * beta_1
    )
    a_1 = 2.0 * np.pi**2 * beta_0 - 4.0 * np.pi * 1j * U * beta_1
    a_2 = np.pi**2 * beta_1
    return a_0, a_1, a_2


def _thrust_bar(a_0, beta_0, beta_1, sigma):
    """Dimensionless time-averaged thrust from Eq. (B3).

    :param a_0: Complex flow coefficient a_0.
    :param beta_0: Complex kinematic amplitude beta_0.
    :param beta_1: Complex kinematic amplitude beta_1.
    :param sigma: Reduced frequency.
    :return: Dimensionless time-averaged thrust T_bar.
    """
    U = 2.0 * np.pi / sigma
    a0r, a0i = np.real(a_0), np.imag(a_0)
    b0r, b0i = np.real(beta_0), np.imag(beta_0)
    b1r, b1i = np.real(beta_1), np.imag(beta_1)
    return (
        np.pi / (4.0 * U**2) * (a0r**2 + a0i**2)
        + np.pi / 2.0 * (a0r * b1r + a0i * b1i)
        + np.pi**3 * (b0r * b1r + b0i * b1i)
    )


def _power_bar(a_0, beta_0, beta_1, sigma):
    """Dimensionless time-averaged input power from Eq. (B5).

    :param a_0: Complex flow coefficient a_0.
    :param beta_0: Complex kinematic amplitude beta_0.
    :param beta_1: Complex kinematic amplitude beta_1.
    :param sigma: Reduced frequency.
    :return: Dimensionless time-averaged power P_bar.
    """
    U = 2.0 * np.pi / sigma
    a0r, a0i = np.real(a_0), np.imag(a_0)
    b0r, b0i = np.real(beta_0), np.imag(beta_0)
    b1r, b1i = np.real(beta_1), np.imag(beta_1)
    return np.pi**2 / 2.0 * (
        a0r * (b0r - b1i) - a0i * (b0i - b1r)
    ) + 2.0 * np.pi**3 * U * (b0r * b1r + b0i * b1i)


def performance(sigma_arr, K, R, eps=1.0):
    """Thrust and power coefficients C_T and C_P over a frequency sweep.

    :param sigma_arr: Array of reduced frequencies (must be positive).
    :param K: Dimensionless spring stiffness.
    :param R: Wing/fluid inertia ratio.
    :param eps: Dimensionless heaving amplitude (default 1.0, cancels in output).
    :return: Tuple (C_T, C_P) arrays.
    """
    b0, b1 = kinematics(sigma_arr, K, R, eps)
    a0, _a1, _a2 = _a_coeffs(sigma_arr, b0, b1)
    scale = 4.0 * np.pi**3 * eps**2
    return (
        _thrust_bar(a0, b0, b1, sigma_arr) / scale,
        _power_bar(a0, b0, b1, sigma_arr) / scale,
    )


def trailing_edge_amplitude(sigma_arr, K, R):
    """Trailing-edge amplitude ratio |h_1| / eps from Eq. (3.3).

    h_1 = (D_1 - 2*D_0) / (D_1 + 2*D_0) = beta_0/2 + beta_1 (at x = +1).
    For the rigid wing (K -> inf), |h_1| -> 1.0 identically.

    :param sigma_arr: Array of reduced frequencies.
    :param K: Dimensionless spring stiffness.
    :param R: Wing/fluid inertia ratio.
    :return: Array of |h_1| values (normalised by eps).
    """
    b0, b1 = kinematics(sigma_arr, K, R, eps=1.0)
    return np.abs(b0 / 2.0 + b1)


def sigma_resonant(K, R):
    """Approximate resonant reduced frequency from Eq. (3.4), valid for K >> 1.

    :param K: Dimensionless spring stiffness.
    :param R: Wing/fluid inertia ratio.
    :return: Approximate sigma_r.
    """
    return np.sqrt(96.0 * K / (128.0 * R + 27.0 * np.pi))


# =============================================================================
# Part 3: Plot configuration
# =============================================================================

# Use a very large K to approximate the rigid (K -> inf) limit numerically.
K_RIGID = 1.0e8

K_VALUES = [0.0, 0.5, 1.0, 4.0, 6.0, 8.0, K_RIGID]
K_LABELS = {
    0.0: "K=0",
    0.5: "K=0.5",
    1.0: "K=1",
    4.0: "K=4",
    6.0: "K=6",
    8.0: "K=8",
    K_RIGID: r"K=$\infty$ (rigid)",
}
K_COLORS = {
    0.0: "gold",
    0.5: "orange",
    1.0: "darkorange",
    4.0: "red",
    6.0: "firebrick",
    8.0: "darkred",
    K_RIGID: "black",
}

R_PLOT = 1.0  # wing/fluid inertia ratio used in the multi-K plots
SIGMA = np.linspace(0.05, 8.0, 1000, dtype=float)

# Stiff K values for which the resonant-frequency dashed lines are meaningful
K_STIFF = {4.0, 6.0, 8.0}


# =============================================================================
# Part 4: Figure 2 equivalent -- wing shape snapshots for R > 0
# =============================================================================

print("Generating Figure 2 equivalent (wing shapes)...")

FIG2_K = 4.0
FIG2_R = R_PLOT
FIG2_N_SNAPS = 20

sigma_r_fig2 = sigma_resonant(FIG2_K, FIG2_R)
FIG2_SIGMAS = [0.5, round(float(sigma_r_fig2), 2), 5.0]
FIG2_LABELS = [
    rf"$\sigma = {FIG2_SIGMAS[0]}$  (low frequency)",
    rf"$\sigma = {FIG2_SIGMAS[1]:.2f}$  (near resonance, Eq. 3.4 estimate)",
    rf"$\sigma = {FIG2_SIGMAS[2]}$  (high frequency)",
]

# eps for visualization -- small enough that pitching is visible but wing stays thin
FIG2_EPS = 0.15
# Chord length in dimensionless data units (used in unnormalized mode).
FIG2_CHORD = 1.2
# True  -> all segments the same pixel length; only direction (pitch angle) varies.
# False -> segment length reflects true TE displacement in data coordinates.
FIG2_NORMALIZE_SEGMENTS = True
# Target segment length in screen pixels (only used when FIG2_NORMALIZE_SEGMENTS).
_PIXEL_LEN = 110.0

# Pre-compute kinematics for all sigma values so we can find the global y range
# before drawing anything.  This ensures all three panels share the same y scale.
_t_dense = np.linspace(0.0, 1.0, 400, dtype=float)
_y_extremes: list[float] = []
_fig2_kin: dict = {}
for _sv in FIG2_SIGMAS:
    _b0, _b1 = kinematics(float(_sv), FIG2_K, FIG2_R, eps=FIG2_EPS)
    _fig2_kin[_sv] = (_b0, _b1)
    for _t in _t_dense:
        _ph = np.exp(2.0j * np.pi * _t)
        _y_extremes.append(float(np.real((_b0 / 2.0 - _b1) * _ph)))
        _y_extremes.append(float(np.real((_b0 / 2.0 + _b1) * _ph)))

_y_abs_max = max(abs(y) for y in _y_extremes)
# Shared y limits with a small margin; same on every subplot.
_y_lim = (-_y_abs_max * 1.25, _y_abs_max * 1.25)

t_snaps = np.linspace(0.0, 1.0, FIG2_N_SNAPS, endpoint=False, dtype=float)
snap_colors = plt.cm.rainbow(np.linspace(0.0, 1.0, FIG2_N_SNAPS, dtype=float))

fig2, axes2 = plt.subplots(3, 1, figsize=(12, 4), dpi=150)
fig2.suptitle(
    rf"Wing swimming through space (K={FIG2_K}, R={FIG2_R})  --  "
    r"Fig. 2 equivalent for $R > 0$  [Moore 2014]",
    fontsize=10,
)

# Pass 1: configure axes (limits, labels) so tight_layout can finalize positions.
for ax, sigma_val, label in zip(axes2, FIG2_SIGMAS, FIG2_LABELS):
    U_val = 2.0 * np.pi / float(sigma_val)
    _pad = U_val / FIG2_N_SNAPS * 0.4
    ax.set_xlim(-_pad, U_val + _pad)
    ax.set_ylim(_y_lim)
    ax.axhline(0.0, color="gray", linewidth=0.5, linestyle="--")
    ax.set_ylabel(r"$h$")
    ax.set_title(label, fontsize=9)

axes2[-1].set_xlabel(
    r"$U \cdot t$   (distance traveled, leading-edge dot marks each snapshot)"
)
plt.tight_layout()
# Force a canvas draw so ax.transData reflects the final axis positions (needed
# for the pixel-space normalization pass below).
fig2.canvas.draw()

# Pass 2: draw segments.
for ax, sigma_val in zip(axes2, FIG2_SIGMAS):
    b0, b1 = _fig2_kin[sigma_val]
    U_val = 2.0 * np.pi / float(sigma_val)
    _trans = ax.transData
    _inv = _trans.inverted()

    for t_val, color in zip(t_snaps, snap_colors):
        phasor = np.exp(2.0j * np.pi * t_val)

        x_le = U_val * t_val
        y_le = float(np.real((b0 / 2.0 - b1) * phasor))
        # Pitch angle: theta ~= Re(beta_1 * phasor), small-angle convention.
        pitch = float(np.real(b1 * phasor))

        if FIG2_NORMALIZE_SEGMENTS:
            # Normalize to _PIXEL_LEN so every segment looks the same size on
            # screen regardless of which panel's x-axis range it sits in.
            pix_le = _trans.transform((x_le, y_le))
            pix_tip = _trans.transform((x_le + 1.0, y_le + pitch))
            pix_dir = pix_tip - pix_le
            pix_norm = float(np.hypot(pix_dir[0], pix_dir[1]))
            if pix_norm > 1e-10:
                pix_te = pix_le + pix_dir * (_PIXEL_LEN / pix_norm)
            else:
                pix_te = pix_le + np.array([_PIXEL_LEN, 0.0])
            x_te, y_te = _inv.transform(pix_te)
        else:
            # Unnormalized: segment length reflects true TE displacement.
            x_te = x_le + FIG2_CHORD
            y_te = y_le + FIG2_CHORD * pitch

        ax.plot(
            [x_le, x_te],
            [y_le, y_te],
            color=color,
            linewidth=1.5,
            alpha=0.9,
            solid_capstyle="round",
        )
        ax.plot(x_le, y_le, "o", color=color, markersize=3, zorder=5)

plt.savefig("moore_2014_fig2_equiv.png")
print("  Saved: moore_2014_fig2_equiv.png")
plt.show()


# =============================================================================
# Part 5: Figure 4a equivalent -- trailing-edge amplitude vs reduced frequency
# =============================================================================

print("Generating Figure 4a equivalent (trailing-edge amplitude)...")

fig4, ax4 = plt.subplots(figsize=(8, 5), dpi=150)

for K in K_VALUES:
    h1 = trailing_edge_amplitude(SIGMA, K, R_PLOT)
    ax4.plot(SIGMA, h1, color=K_COLORS[K], label=K_LABELS[K], linewidth=1.8)
    if K in K_STIFF:
        sr = sigma_resonant(K, R_PLOT)
        ax4.axvline(sr, color=K_COLORS[K], linewidth=0.8, linestyle="--", alpha=0.5)

ax4.axhline(1.0, color="gray", linewidth=0.6, linestyle=":", label="_nolegend_")
ax4.set_xlabel(r"Reduced frequency  $\sigma = \pi f c / U_\infty$")
ax4.set_ylabel(r"Trailing-edge amplitude  $|h_1| / \varepsilon$")
ax4.set_title(
    rf"Trailing-edge amplitude vs. reduced frequency  (R = {R_PLOT})"
    "\n"
    r"Fig. 4a equivalent from Moore (2014) for $R > 0$"
    "  --  dashed verticals from Eq. (3.4)"
)
ax4.legend(fontsize=9, loc="upper right")
ax4.set_xlim(0.0, 8.0)
ax4.set_ylim(0.0, 3.5)
ax4.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"moore_2014_fig4a_R{R_PLOT}.png")
print(f"  Saved: moore_2014_fig4a_R{R_PLOT}.png")
plt.show()


# =============================================================================
# Part 6: Figure 3/6 equivalent -- C_T and C_P vs reduced frequency for R > 0
# =============================================================================

print(f"Generating Figure 3/6 equivalent (C_T, C_P vs sigma, R={R_PLOT})...")

fig36, (ax_ct, ax_cp) = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

for K in K_VALUES:
    C_T, C_P = performance(SIGMA, K, R_PLOT)
    ax_ct.plot(SIGMA, C_T, color=K_COLORS[K], label=K_LABELS[K], linewidth=1.8)
    ax_cp.plot(SIGMA, C_P, color=K_COLORS[K], label=K_LABELS[K], linewidth=1.8)

for ax in (ax_ct, ax_cp):
    ax.axhline(0.0, color="gray", linewidth=0.5, linestyle="--")
    ax.set_xlabel(r"Reduced frequency  $\sigma$")
    ax.set_xlim(0.0, 5.0)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

ax_ct.set_ylabel(r"Thrust coefficient  $C_T$")
ax_ct.set_title(rf"Thrust coefficient  (R = {R_PLOT})")
ax_cp.set_ylabel(r"Power coefficient  $C_P$")
ax_cp.set_title(rf"Power coefficient  (R = {R_PLOT})")
fig36.suptitle(
    rf"Propulsive performance for R = {R_PLOT}  --  " "Fig. 3/6 equivalent [Moore 2014]"
)
plt.tight_layout()
plt.savefig(f"moore_2014_fig36_R{R_PLOT}.png")
print(f"  Saved: moore_2014_fig36_R{R_PLOT}.png")
plt.show()


# =============================================================================
# Part 7: Ptera Software parameter mapping
# =============================================================================

print()
print("=" * 60)
print("Part 7: Ptera Software parameter mapping to Moore (2014)")
print("=" * 60)

# Parameters from aeroelastic_unsteady_first_order_validation.py
rho_fluid = 0.1  # kg/m^3 (fluid density)
U_inf = 10.0  # m/s (vCg__E)
chord = 1.75  # m
wing_density_ptera = 6.0  # kg/m^2  (b * rho_s)
spring_constant_ptera = 20000.0  # N*m/rad

K_ptera = spring_constant_ptera / (rho_fluid * U_inf**2 * chord**2)
R_ptera = wing_density_ptera / (rho_fluid * chord)
sigma_r_ptera = sigma_resonant(K_ptera, R_ptera)
f_r_ptera = sigma_r_ptera * U_inf / (np.pi * chord)

print(
    f"rho_fluid        = {rho_fluid} kg/m^3\n"
    f"U_inf            = {U_inf} m/s\n"
    f"chord            = {chord} m\n"
    f"wing_density     = {wing_density_ptera} kg/m^2\n"
    f"spring_constant  = {spring_constant_ptera} N*m/rad\n"
)
print(f"Dimensionless K  = {K_ptera:.1f}")
print(f"Dimensionless R  = {R_ptera:.1f}")
print(f"sigma_r (Eq. 3.4) = {sigma_r_ptera:.3f}")
print(f"Resonant frequency = {f_r_ptera:.2f} Hz")
print()
print("Reduced frequency vs. flapping frequency (sigma = pi*f*c/U_inf):")
for f_hz in [1.0, 3.0, 5.0, 6.0, 8.0, 10.0]:
    sigma_hz = np.pi * f_hz * chord / U_inf
    print(f"  f = {f_hz:5.1f} Hz  ->  sigma = {sigma_hz:.3f}")
