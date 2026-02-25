"""Plot induced thrust vs voltage for single-wing flapping simulations.

Automatically plots all available mesh types (Coarse and/or Fine) from a run.

Usage:     python plot_single_wing_thrust.py <run_id> [--no-show] [--no-experimental]
python plot_single_wing_thrust.py --latest [--no-show] [--no-experimental]

Examples:     python plot_single_wing_thrust.py 2026-01-18_17-04-43     python
plot_single_wing_thrust.py --latest     python plot_single_wing_thrust.py --latest --no-
show     python plot_single_wing_thrust.py --latest --no-experimental
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt

# Add parent directory to path for local imports
sys.path.insert(0, str(Path(__file__).parent))

import experimental_thrust


@dataclass
class SingleWingResult:
    """Result data for a single-wing flapping simulation."""

    config_name: str
    voltage: int
    flapping_wing: str  # "left" or "right"
    drag_force: float  # induced drag force in Newtons (negative = thrust)
    mesh_type: str  # "Coarse" or "Fine"


def save_analysis_scripts_to_snapshot(run_dir: Path) -> None:
    """Save this script and experimental_thrust.py to the codebase snapshot.

    :param run_dir: Path to the run directory containing codebase_snapshot/.
    """
    snapshot_dir = run_dir / "codebase_snapshot" / "gammabot_simulations"
    if not snapshot_dir.exists():
        print(f"Warning: Snapshot directory not found: {snapshot_dir}")
        return

    script_dir = Path(__file__).parent
    files_to_copy = ["plot_single_wing_thrust.py", "experimental_thrust.py"]

    for filename in files_to_copy:
        src = script_dir / filename
        dst = snapshot_dir / filename
        if src.exists():
            shutil.copy2(src, dst)
            print(f"Saved {filename} to snapshot")
        else:
            print(f"Warning: {filename} not found")


def find_latest_run(results_dir: Path) -> str | None:
    """Find the most recent run based on timestamp directory names.

    :param results_dir: Path to the simulation_results directory.
    :return: The name of the latest run directory, or None if no runs found.
    """
    if not results_dir.exists():
        return None

    # Get all directories that look like timestamps (YYYY-MM-DD_HH-MM-SS)
    timestamp_pattern = re.compile(r"^\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}$")
    timestamp_dirs = [
        d.name
        for d in results_dir.iterdir()
        if d.is_dir() and timestamp_pattern.match(d.name)
    ]

    if not timestamp_dirs:
        # Fall back to all directories sorted alphabetically
        all_dirs = [d.name for d in results_dir.iterdir() if d.is_dir()]
        return max(all_dirs) if all_dirs else None

    # Sort and return the latest (timestamps sort correctly as strings)
    return max(timestamp_dirs)


def parse_voltage_from_config(config_name: str) -> tuple[int, int]:
    """Extract left and right voltage values from a config name.

    :param config_name: Configuration name (e.g., "L150V_R0V_170Hz").
    :return: A tuple of (left_voltage, right_voltage).
    :raises ValueError: If the config name doesn't match the expected pattern.
    """
    pattern = r"L(\d+)V_R(\d+)V_170Hz"
    match = re.match(pattern, config_name)
    if not match:
        raise ValueError(f"Config name doesn't match expected pattern: {config_name}")
    return int(match.group(1)), int(match.group(2))


def is_single_wing_config(config_name: str) -> tuple[bool, str | None, int | None]:
    """Check if a config represents single-wing flapping.

    :param config_name: Configuration name to check.
    :return: A tuple of (is_single_wing, flapping_wing, voltage) where flapping_wing is
        "left" or "right" and voltage is the active wing's voltage.
    """
    try:
        left_v, right_v = parse_voltage_from_config(config_name)
    except ValueError:
        return False, None, None

    if left_v == 0 and right_v > 0:
        return True, "right", right_v
    elif right_v == 0 and left_v > 0:
        return True, "left", left_v
    else:
        return False, None, None


def extract_drag_force_from_results(results_file: Path) -> float | None:
    """Extract the induced drag force from a print_results.txt file.

    :param results_file: Path to the print_results.txt file.
    :return: The Drag force value in Newtons, or None if not found.
    """
    if not results_file.exists():
        return None

    text = results_file.read_text()

    # Look for the Drag value in the results (now in scientific notation)
    # Format: "FX_W:      1.234e-05 N      Drag:             -1.234e-05 N"
    pattern = r"Drag:\s+([-\d.eE+]+)\s*N"
    match = re.search(pattern, text)
    if match:
        return float(match.group(1))
    return None


def collect_all_single_wing_results(run_dir: Path) -> list[SingleWingResult]:
    """Collect results from all single-wing flapping simulations in a run.

    Collects results from all mesh types (Coarse, Fine, Extra-Fine) if available.

    :param run_dir: Path to the run directory (e.g., simulation_results/run_id/).
    :return: A list of SingleWingResult objects with mesh_type populated.
    """
    results: list[SingleWingResult] = []

    if not run_dir.exists():
        print(f"Run directory not found: {run_dir}")
        return results

    mesh_types = ["Coarse", "Fine", "Extra-Fine"]

    # Iterate over config directories
    for config_dir in run_dir.iterdir():
        if not config_dir.is_dir():
            continue

        config_name = config_dir.name

        # Skip non-config directories (like codebase_snapshot)
        is_single, flapping_wing, voltage = is_single_wing_config(config_name)
        if not is_single or flapping_wing is None or voltage is None:
            continue

        # Look for results in both mesh types
        for mesh_type in mesh_types:
            results_file = config_dir / mesh_type / "print_results.txt"
            drag_force = extract_drag_force_from_results(results_file)

            if drag_force is not None:
                results.append(
                    SingleWingResult(
                        config_name=config_name,
                        voltage=voltage,
                        flapping_wing=flapping_wing,
                        drag_force=drag_force,
                        mesh_type=mesh_type,
                    )
                )

    return results


def plot_thrust_vs_voltage(
    results: list[SingleWingResult],
    run_id: str,
    show: bool = True,
    save: bool = True,
    output_dir: Path | None = None,
    show_experimental: bool = True,
) -> None:
    """Plot average induced thrust vs voltage for single-wing simulations.

    Plots all available mesh types (Coarse and/or Fine) with different markers. Always
    shows all experimental data points regardless of which simulation results are
    available.

    :param results: A list of SingleWingResult objects (may include multiple mesh
        types).
    :param run_id: The run identifier (used in title).
    :param show: Whether to display the plot interactively.
    :param save: Whether to save the plot to a file.
    :param output_dir: Directory to save the plot. If None, uses current directory.
    :param show_experimental: Whether to overlay experimental data on the plot.
    """
    # Get experimental data (always fetch to determine x-axis range)
    exp_data = experimental_thrust.THRUST_MN
    exp_left_voltages = sorted(exp_data["left"].keys())
    exp_left_thrust_mN = [exp_data["left"][v] for v in exp_left_voltages]
    exp_right_voltages = sorted(exp_data["right"].keys())
    exp_right_thrust_mN = [exp_data["right"][v] for v in exp_right_voltages]

    # Determine which mesh types are present in results
    mesh_types_present = sorted(set(r.mesh_type for r in results))

    # Style configuration for each mesh type
    mesh_styles = {
        "Coarse": {"marker": "o", "linestyle": "-", "color": "tab:blue"},
        "Fine": {"marker": "s", "linestyle": "-", "color": "tab:orange"},
        "Extra-Fine": {"marker": "D", "linestyle": "-", "color": "tab:red"},
    }

    # Collect all thrust values to determine y-axis limits
    all_thrust: list[float] = []

    # Add experimental data to thrust range calculation
    if show_experimental:
        all_thrust.extend(exp_left_thrust_mN)
        all_thrust.extend(exp_right_thrust_mN)

    # Add simulation data to thrust range calculation
    for r in results:
        thrust_mN = -r.drag_force * 1000
        all_thrust.append(thrust_mN)

    # Determine y-axis limits
    if all_thrust:
        y_min = min(min(all_thrust), 0)
        y_max = max(max(all_thrust), 0)
        y_margin = (y_max - y_min) * 0.1 if y_max != y_min else 0.1
        y_limits = (y_min - y_margin, y_max + y_margin)
    else:
        y_limits = (-1, 1)

    # Determine x-axis limits based on experimental data (always show full range)
    all_voltages = exp_left_voltages + exp_right_voltages
    for r in results:
        all_voltages.append(r.voltage)
    if all_voltages:
        x_min = min(all_voltages) - 5
        x_max = max(all_voltages) + 5
    else:
        x_min, x_max = 145, 185

    # Create figure with two subplots
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Plot experimental data first (so simulation data appears on top)
    if show_experimental:
        if exp_left_voltages:
            ax_left.plot(
                exp_left_voltages,
                exp_left_thrust_mN,
                "^--",
                color="tab:green",
                linewidth=2,
                markersize=8,
                label="Experimental",
            )
        if exp_right_voltages:
            ax_right.plot(
                exp_right_voltages,
                exp_right_thrust_mN,
                "^--",
                color="tab:green",
                linewidth=2,
                markersize=8,
                label="Experimental",
            )

    # Plot simulation data for each mesh type
    for mesh_type in mesh_types_present:
        style = mesh_styles.get(
            mesh_type, {"marker": "x", "linestyle": ":", "color": "tab:gray"}
        )

        # Filter results for this mesh type
        mesh_results = [r for r in results if r.mesh_type == mesh_type]

        # Left wing results
        left_results = [r for r in mesh_results if r.flapping_wing == "left"]
        left_results.sort(key=lambda r: r.voltage)
        if left_results:
            left_voltages = [r.voltage for r in left_results]
            left_thrust_mN = [-r.drag_force * 1000 for r in left_results]
            ax_left.plot(
                left_voltages,
                left_thrust_mN,
                marker=style["marker"],
                linestyle=style["linestyle"],
                color=style["color"],
                linewidth=2,
                markersize=8,
                label=f"Simulation ({mesh_type})",
            )

        # Right wing results
        right_results = [r for r in mesh_results if r.flapping_wing == "right"]
        right_results.sort(key=lambda r: r.voltage)
        if right_results:
            right_voltages = [r.voltage for r in right_results]
            right_thrust_mN = [-r.drag_force * 1000 for r in right_results]
            ax_right.plot(
                right_voltages,
                right_thrust_mN,
                marker=style["marker"],
                linestyle=style["linestyle"],
                color=style["color"],
                linewidth=2,
                markersize=8,
                label=f"Simulation ({mesh_type})",
            )

    # Configure left wing subplot
    ax_left.axhline(y=0, color="gray", linestyle="--", linewidth=0.8)
    ax_left.set_xlabel("Voltage Amplitude (V)")
    ax_left.set_ylabel("Average Induced Thrust (mN)")
    ax_left.set_title("Left Wing Flapping Only")
    ax_left.grid(True, alpha=0.3)
    ax_left.set_ylim(y_limits)
    ax_left.set_xlim(x_min, x_max)
    ax_left.legend(loc="best")

    # Configure right wing subplot
    ax_right.axhline(y=0, color="gray", linestyle="--", linewidth=0.8)
    ax_right.set_xlabel("Voltage Amplitude (V)")
    ax_right.set_title("Right Wing Flapping Only")
    ax_right.grid(True, alpha=0.3)
    ax_right.set_xlim(x_min, x_max)
    ax_right.legend(loc="best")

    # Build title showing which mesh types are included
    if mesh_types_present:
        mesh_label = " + ".join(mesh_types_present)
    else:
        mesh_label = "No simulation data"

    fig.suptitle(
        f"Induced Thrust vs Voltage - Single Wing Flapping\n"
        f"Run: {run_id} ({mesh_label})",
        fontsize=12,
    )

    plt.tight_layout()

    if save:
        save_dir = output_dir if output_dir else Path.cwd()
        filename = save_dir / "single_wing_thrust.png"
        fig.savefig(filename, dpi=150, bbox_inches="tight")
        print(f"Plot saved to: {filename}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    """Parse arguments and generate the plot."""
    parser = argparse.ArgumentParser(
        description="Plot induced thrust vs voltage for single-wing simulations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "run_id",
        nargs="?",
        default=None,
        help="Run identifier (tag or timestamp, e.g., '2026-01-18_17-04-43')",
    )
    parser.add_argument(
        "--latest",
        action="store_true",
        help="Use the most recent run (based on timestamp)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display the plot interactively",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Don't save the plot to a file",
    )
    parser.add_argument(
        "--no-experimental",
        action="store_true",
        help="Don't overlay experimental data on the plot",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.run_id and not args.latest:
        parser.error("Either run_id or --latest is required")

    # Construct path to run directory
    script_dir = Path(__file__).parent
    results_dir = script_dir / "simulation_results"

    if args.latest:
        run_id = find_latest_run(results_dir)
        if run_id is None:
            print("Error: No runs found in simulation_results/")
            return
        print(f"Using latest run: {run_id}")
    else:
        run_id = args.run_id

    run_dir = results_dir / run_id

    if not run_dir.exists():
        print(f"Error: Run directory not found: {run_dir}")
        print("Available runs:")
        results_dir = script_dir / "simulation_results"
        if results_dir.exists():
            for d in sorted(results_dir.iterdir()):
                if d.is_dir():
                    print(f"  {d.name}")
        return

    # Collect results from all mesh types
    print(f"Collecting results from: {run_dir}")
    results = collect_all_single_wing_results(run_dir)

    # Report what was found
    mesh_types_found = sorted(set(r.mesh_type for r in results))
    if mesh_types_found:
        print(f"Mesh types found: {', '.join(mesh_types_found)}")
    else:
        print("No simulation results found (will show experimental data only)")

    if results:
        print(f"Found {len(results)} single-wing simulation results:")
        for r in sorted(
            results, key=lambda x: (x.mesh_type, x.flapping_wing, x.voltage)
        ):
            thrust_mN = -r.drag_force * 1000  # Convert to mN
            print(
                f"  [{r.mesh_type}] {r.config_name}: {r.flapping_wing} wing, "
                f"{r.voltage}V, thrust={thrust_mN:.4f} mN"
            )

    # Generate plot (works even with no simulation results if experimental is enabled)
    output_dir = run_dir
    plot_thrust_vs_voltage(
        results=results,
        run_id=run_id,
        show=not args.no_show,
        save=not args.no_save,
        output_dir=output_dir,
        show_experimental=not args.no_experimental,
    )

    # Save analysis scripts to the codebase snapshot
    save_analysis_scripts_to_snapshot(run_dir)


if __name__ == "__main__":
    main()
