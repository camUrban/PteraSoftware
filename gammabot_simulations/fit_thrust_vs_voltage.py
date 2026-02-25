"""Fit linear models to experimental thrust vs voltage data for GammaBot wings.

Performs linear regression on experimental thrust data for both left and right wings,
prints fit equations with R^2 values, and plots the results.

Usage:
    python fit_thrust_vs_voltage.py [--no-show]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path for local imports
sys.path.insert(0, str(Path(__file__).parent))

import experimental_thrust


def format_coefficient(value: float, sig_figs: int = 5) -> str:
    """Format a coefficient to specified significant figures.

    :param value: The coefficient value to format.
    :param sig_figs: Number of significant figures.
    :return: Formatted string representation.
    """
    if value == 0:
        return "0"
    return f"{value:.{sig_figs}g}"


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Perform linear regression and compute R^2.

    :param x: Independent variable array.
    :param y: Dependent variable array.
    :return: Tuple of (slope, intercept, r_squared).
    """
    coeffs = np.polyfit(x, y, 1)
    slope = float(coeffs[0])
    intercept = float(coeffs[1])

    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = float(1 - ss_res / ss_tot)

    return slope, intercept, r_squared


def main() -> None:
    """Perform linear fits and plot results."""
    parser = argparse.ArgumentParser(
        description="Fit linear models to experimental thrust vs voltage data.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display the plot interactively",
    )
    args = parser.parse_args()

    exp_data = experimental_thrust.THRUST_MN

    left_voltages = np.array(sorted(exp_data["left"].keys()), dtype=np.float64)
    left_thrust = np.array(
        [exp_data["left"][v] for v in left_voltages], dtype=np.float64
    )

    right_voltages = np.array(sorted(exp_data["right"].keys()), dtype=np.float64)
    right_thrust = np.array(
        [exp_data["right"][v] for v in right_voltages], dtype=np.float64
    )

    left_slope, left_intercept, left_r2 = linear_fit(left_voltages, left_thrust)
    right_slope, right_intercept, right_r2 = linear_fit(right_voltages, right_thrust)

    left_eq = (
        f"y = {format_coefficient(left_slope)}x + {format_coefficient(left_intercept)}"
    )
    right_eq = f"y = {format_coefficient(right_slope)}x + {format_coefficient(right_intercept)}"

    print("Linear Fit Results (Thrust [mN] vs Voltage [V])")
    print("=" * 55)
    print()
    print("Left Wing:")
    print(f"  {left_eq}")
    print(f"  R^2 = {format_coefficient(left_r2)}")
    print()
    print("Right Wing:")
    print(f"  {right_eq}")
    print(f"  R^2 = {format_coefficient(right_r2)}")
    print()

    script_dir = Path(__file__).parent

    left_fit_file = script_dir / "left_wing_thrust_fit.txt"
    left_fit_file.write_text(
        f"Left Wing Thrust vs Voltage Linear Fit\n"
        f"=======================================\n"
        f"Equation: {left_eq}\n"
        f"R^2 = {format_coefficient(left_r2)}\n"
        f"\n"
        f"Where:\n"
        f"  y = Thrust (mN)\n"
        f"  x = Voltage Amplitude (V)\n"
    )
    print(f"Left wing fit saved to: {left_fit_file}")

    right_fit_file = script_dir / "right_wing_thrust_fit.txt"
    right_fit_file.write_text(
        f"Right Wing Thrust vs Voltage Linear Fit\n"
        f"========================================\n"
        f"Equation: {right_eq}\n"
        f"R^2 = {format_coefficient(right_r2)}\n"
        f"\n"
        f"Where:\n"
        f"  y = Thrust (mN)\n"
        f"  x = Voltage Amplitude (V)\n"
    )
    print(f"Right wing fit saved to: {right_fit_file}")
    print()

    x_fit = np.linspace(145, 185, 100)
    left_fit = left_slope * x_fit + left_intercept
    right_fit = right_slope * x_fit + right_intercept

    all_thrust = (
        list(left_thrust) + list(right_thrust) + list(left_fit) + list(right_fit)
    )
    y_min = min(min(all_thrust), 0)
    y_max = max(max(all_thrust), 0)
    y_margin = (y_max - y_min) * 0.1 if y_max != y_min else 0.1
    y_limits = (y_min - y_margin, y_max + y_margin)

    x_min = 145
    x_max = 185

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    ax_left.plot(
        left_voltages,
        left_thrust,
        "^",
        color="tab:green",
        markersize=8,
        label="Experimental",
    )
    ax_left.plot(
        x_fit,
        left_fit,
        "-",
        color="tab:blue",
        linewidth=2,
        label="Linear Fit",
    )

    ax_right.plot(
        right_voltages,
        right_thrust,
        "^",
        color="tab:green",
        markersize=8,
        label="Experimental",
    )
    ax_right.plot(
        x_fit,
        right_fit,
        "-",
        color="tab:blue",
        linewidth=2,
        label="Linear Fit",
    )

    ax_left.axhline(y=0, color="gray", linestyle="--", linewidth=0.8)
    ax_left.set_xlabel("Voltage Amplitude (V)")
    ax_left.set_ylabel("Average Thrust (mN)")
    ax_left.set_title("Left Wing")
    ax_left.grid(True, alpha=0.3)
    ax_left.set_ylim(y_limits)
    ax_left.set_xlim(x_min, x_max)
    ax_left.legend(loc="upper left")

    left_eq_display = f"{left_eq}\n$R^2$ = {format_coefficient(left_r2)}"
    ax_left.text(
        0.95,
        0.05,
        left_eq_display,
        transform=ax_left.transAxes,
        fontsize=10,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )

    ax_right.axhline(y=0, color="gray", linestyle="--", linewidth=0.8)
    ax_right.set_xlabel("Voltage Amplitude (V)")
    ax_right.set_title("Right Wing")
    ax_right.grid(True, alpha=0.3)
    ax_right.set_xlim(x_min, x_max)
    ax_right.legend(loc="upper left")

    right_eq_display = f"{right_eq}\n$R^2$ = {format_coefficient(right_r2)}"
    ax_right.text(
        0.95,
        0.05,
        right_eq_display,
        transform=ax_right.transAxes,
        fontsize=10,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )

    fig.suptitle(
        "Experimental Thrust vs Voltage - Linear Fit\n" "GammaBot Single Wing Flapping",
        fontsize=12,
    )

    plt.tight_layout()

    script_dir = Path(__file__).parent
    filename = script_dir / "thrust_linear_fit.png"
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {filename}")

    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
