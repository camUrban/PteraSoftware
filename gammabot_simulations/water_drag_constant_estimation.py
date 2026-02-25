"""Estimate water drag constant from GammaBot coasting data.

Fits the analytical solution to the drag ODE (m*a = -0.5*rho*Cd*A*v^2) to
experimental displacement vs. time data. The initial velocity v_0 is computed
numerically from the first two data points, and the lumped drag parameter
c_1 = 0.5*rho*Cd*A/m is found via nonlinear least squares fitting.

The ODE is: a = -c_1 * v^2 (c_1 > 0)
The analytical solution is: x(t) = ln(1 + c_1*v_0*t) / c_1

Usage:
    python water_drag_constant_estimation.py [--no-show]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit


def format_coefficient(value: float, sig_figs: int = 5) -> str:
    """Format a coefficient to specified significant figures.

    :param value: The coefficient value to format.
    :param sig_figs: Number of significant figures.
    :return: Formatted string representation.
    """
    if value == 0:
        return "0"
    return f"{value:.{sig_figs}g}"


def displacement_model(
    t: NDArray[np.float64], c1: float, v0: float
) -> NDArray[np.float64]:
    """Analytical solution for displacement under quadratic drag.

    Solves the ODE a = -c_1 * v^2 with x(0) = 0 and v(0) = v_0, where c_1 > 0.
    Solution: x(t) = ln(1 + c_1*v_0*t) / c_1

    :param t: Time array.
    :param c1: Lumped drag parameter (0.5*rho*Cd*A/m), positive.
    :param v0: Initial velocity.
    :return: Displacement array.
    """
    return np.log(1 + c1 * v0 * t) / c1


def make_model_with_fixed_v0(
    v0: float,
) -> callable:
    """Create a displacement model function with v0 fixed.

    :param v0: Fixed initial velocity value.
    :return: Model function that only takes t and c1 as arguments.
    """

    def model(t: NDArray[np.float64], c1: float) -> NDArray[np.float64]:
        return displacement_model(t, c1, v0)

    return model


def compute_r_squared(
    y_data: NDArray[np.float64], y_pred: NDArray[np.float64]
) -> float:
    """Compute R^2 (coefficient of determination).

    :param y_data: Observed data.
    :param y_pred: Predicted values.
    :return: R^2 value.
    """
    ss_res = np.sum((y_data - y_pred) ** 2)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot != 0 else 0.0


def main() -> None:
    """Load coasting data and fit the drag model."""
    parser = argparse.ArgumentParser(
        description="Estimate water drag constant from coasting data.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display the plot interactively",
    )
    args = parser.parse_args()

    # Load experimental data from CSV
    script_dir = Path(__file__).parent
    csv_file = script_dir / "experimental_coasting.csv"
    data = np.genfromtxt(csv_file, delimiter=",", skip_header=1)
    time = data[:, 0].astype(np.float64)
    displacement = data[:, 1].astype(np.float64)

    # Compute initial velocity from first 20 data points using linear fit (in mm/s)
    # This reduces noise sensitivity compared to using just 2 points
    n_points = 20
    coeffs = np.polyfit(time[:n_points], displacement[:n_points], 1)
    v0 = coeffs[0]  # Slope is the initial velocity

    # Create model function with fixed v0
    model = make_model_with_fixed_v0(v0)

    # Initial guess for c1 (must be positive)
    c1_initial = 0.1

    # Fit the model to displacement data
    popt, pcov = curve_fit(
        model,
        time,
        displacement,
        p0=[c1_initial],
        maxfev=10000,
    )
    c1 = popt[0]

    # Compute fitted displacement and R^2
    displacement_fit = model(time, c1)
    r_squared = compute_r_squared(displacement, displacement_fit)

    # Convert results to meters for display
    # v_0: mm/s -> m/s (divide by 1000)
    # c_1: 1/mm -> 1/m (multiply by 1000)
    v0_m = v0 / 1000
    c1_m = c1 * 1000

    # Build equation string
    equation = "x(t) = ln(1 + c_1*v_0*t) / c_1"

    # Print results
    print("Water Drag Constant Estimation")
    print("=" * 55)
    print()
    print("Model: a = -c_1*v^2, where c_1 = 0.5*rho*Cd*A/m")
    print(f"Solution: {equation}")
    print()
    print("Fitted Parameters:")
    print(f"  v_0 = {format_coefficient(v0_m)} m/s (from data)")
    print(f"  c_1 = {format_coefficient(c1_m)} 1/m")
    print(f"  R^2 = {format_coefficient(r_squared)}")
    print()

    # Save fit results to text file
    fit_file = script_dir / "water_drag_fit.txt"
    fit_file.write_text(
        f"Water Drag Constant Estimation\n"
        f"===============================\n"
        f"\n"
        f"Model: a = -c_1*v^2, where c_1 = 0.5*rho*Cd*A/m\n"
        f"Solution: {equation}\n"
        f"\n"
        f"Fitted Parameters:\n"
        f"  v_0 = {format_coefficient(v0_m)} m/s (from data)\n"
        f"  c_1 = {format_coefficient(c1_m)} 1/m\n"
        f"  R^2 = {format_coefficient(r_squared)}\n"
        f"\n"
        f"Where:\n"
        f"  c_1 = 0.5*rho*Cd*A/m (drag parameter)\n"
        f"  v_0 = initial velocity\n"
        f"  x(t) = displacement (m)\n"
        f"  t = time (s)\n"
    )
    print(f"Fit saved to: {fit_file}")
    print()

    # Generate smooth fitted curve
    t_fit = np.linspace(time[0], time[-1], 200)
    x_fit = model(t_fit, c1)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot experimental data
    ax.plot(
        time,
        displacement,
        "^",
        color="tab:green",
        markersize=8,
        label="Experimental",
    )

    # Plot fitted curve
    ax.plot(
        t_fit,
        x_fit,
        "-",
        color="tab:blue",
        linewidth=2,
        label="Fitted Model",
    )

    ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Displacement (mm)")
    ax.set_title("GammaBot Coasting - Water Drag Model Fit")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")

    # Add equations and fitted parameters annotation
    eq_display = (
        r"$a = -c_1 v^2$" + "\n"
        r"$x(t) = \frac{\ln(1 + c_1 v_0 t)}{c_1}$" + "\n\n"
        f"$v_0$ = {format_coefficient(v0_m)} m/s\n"
        f"$c_1$ = {format_coefficient(c1_m)} 1/m\n"
        f"$R^2$ = {format_coefficient(r_squared)}"
    )
    ax.text(
        0.95,
        0.05,
        eq_display,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )

    plt.tight_layout()

    # Save plot
    filename = script_dir / "water_drag_estimation.png"
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {filename}")

    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
