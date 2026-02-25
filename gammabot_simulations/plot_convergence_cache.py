"""Plot convergence of a coefficient or load from the convergence cache.

Creates a single column of subplots (one per panel aspect ratio) showing how the
selected quantity varies with num_chordwise_panels refinement. All wake lengths are
overlaid on each subplot for direct comparison. Prescribed wake uses solid lines with
circle markers, free wake uses dashed lines with triangle markers.

The number of subplots and lines per subplot is determined entirely by what parameter
combinations are present in the cache file.

Usage:
    python plot_convergence_cache.py
"""

from __future__ import annotations

import contextlib
import json
import math
import sys
from pathlib import Path
from typing import Generator

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

COEFFICIENT_NAMES = ("cFX", "cFY", "cFZ", "cMX", "cMY", "cMZ")
LOAD_NAMES = ("FX", "FY", "FZ", "MX", "MY", "MZ")
_N_TO_MGF = 1.0 / 9.80665e-6

_LOAD_UNITS_MGF = ("mgf", "mgf", "mgf", "N*m", "N*m", "N*m")
_LOAD_SCALES_MGF = (_N_TO_MGF, _N_TO_MGF, _N_TO_MGF, 1.0, 1.0, 1.0)
_LOAD_UNITS_N = ("N", "N", "N", "N*m", "N*m", "N*m")
_LOAD_SCALES_N = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

# Perceptually uniform sequential colormap for wake length lines.
# Only the middle 75% of the colormap is sampled to avoid the extreme light/dark ends.
WAKE_LENGTH_COLORMAP = "viridis"
_CMAP_V_MIN = 0.125
_CMAP_V_MAX = 0.975

# IQR multiplier for outlier detection when computing y-axis limits.
_IQR_MULTIPLIER = 3.0

# Experimental thrust in mgf (displayed only when plotting FX).
_EXPERIMENTAL_THRUST_MGF = 87.0

# ── Configuration ──────────────────────────────────────────────────────────────
CACHE_PATH = Path(__file__).parent / "convergence_cache" / "L170V_R180V_170Hz.json"
QUANTITY = "FX"
FORCE_UNIT_MGF = True  # True for mgf (whole numbers), False for N (scientific notation)
# ───────────────────────────────────────────────────────────────────────────────


@contextlib.contextmanager
def _lock_cache_file(cache_path: Path) -> Generator[None, None, None]:
    """Acquire an exclusive file lock for safe concurrent cache access.

    Uses a .lock file adjacent to the cache file. On Windows uses msvcrt locking, on
    Unix uses fcntl file locking.

    :param cache_path: Path to the cache file to lock.
    """
    lock_path = cache_path.with_suffix(cache_path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_fh = open(lock_path, "w")
    try:
        if sys.platform == "win32":
            import msvcrt

            # Write a byte so there is content to lock.
            lock_fh.write(" ")
            lock_fh.flush()
            lock_fh.seek(0)
            # LK_LOCK retries for approximately 10 seconds before raising OSError.
            msvcrt.locking(lock_fh.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl

            fcntl.flock(lock_fh, fcntl.LOCK_EX)
        yield
    finally:
        try:
            if sys.platform == "win32":
                import msvcrt

                lock_fh.seek(0)
                msvcrt.locking(lock_fh.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(lock_fh, fcntl.LOCK_UN)
        except OSError:
            pass
        lock_fh.close()


def _load_cache(cache_path: Path) -> dict:
    """Load the convergence cache JSON file with file locking.

    :param cache_path: Path to the cache JSON file.
    :return: Parsed JSON data as a dictionary.
    """
    with _lock_cache_file(cache_path):
        with open(cache_path, "r") as f:
            return json.load(f)


def _parse_cache(data: dict) -> list[dict]:
    """Parse cache entries into a list of structured records.

    :param data: Raw cache dictionary (will be modified to remove _optimized_dt).
    :return: List of dictionaries with parsed fields.
    """
    data.pop("_optimized_dt", None)

    records = []
    for key, value in data.items():
        parts = key.split(",")
        records.append(
            {
                "prescribed_wake": parts[0].strip() == "True",
                "wake_length": int(parts[1].strip()),
                "panel_ar": int(parts[2].strip()),
                "num_chordwise": int(parts[3].strip()),
                "coefficients": value["coefficients"],
                "loads": value["loads"],
                "time": value["time"],
            }
        )
    return records


def _resolve_quantity(
    quantity_name: str,
    force_unit_mgf: bool = True,
) -> tuple[str, int, str, str, float]:
    """Resolve a quantity name to a data key, index, display name, unit, and scale.

    :param quantity_name: Coefficient name (e.g. "cFX") or load name (e.g. "FX").
    :param force_unit_mgf: If True, forces use mgf; if False, forces use N.
    :return: Tuple of (data_key, index, display_name, unit, scale).
    """
    for i, name in enumerate(COEFFICIENT_NAMES):
        if quantity_name == name:
            return "coefficients", i, name, "", 1.0

    load_units = _LOAD_UNITS_MGF if force_unit_mgf else _LOAD_UNITS_N
    load_scales = _LOAD_SCALES_MGF if force_unit_mgf else _LOAD_SCALES_N

    for i, name in enumerate(LOAD_NAMES):
        if quantity_name == name:
            return "loads", i, name, load_units[i], load_scales[i]

    all_names = ", ".join(COEFFICIENT_NAMES) + ", " + ", ".join(LOAD_NAMES)
    raise ValueError(f"Unknown quantity '{quantity_name}'. Valid names: {all_names}.")


def _filter_records(
    records: list[dict],
    panel_ar: int,
    wake_length: int,
    prescribed_wake: bool,
) -> tuple[list[int], list[dict]]:
    """Filter and sort records for a specific parameter combination.

    :param records: All parsed cache records.
    :param panel_ar: Panel aspect ratio to filter by.
    :param wake_length: Wake length to filter by.
    :param prescribed_wake: Whether to filter for prescribed (True) or free (False) wake.
    :return: Tuple of (num_chordwise values, filtered records), both sorted by
        num_chordwise.
    """
    filtered = [
        r
        for r in records
        if r["panel_ar"] == panel_ar
        and r["wake_length"] == wake_length
        and r["prescribed_wake"] == prescribed_wake
    ]
    filtered.sort(key=lambda r: r["num_chordwise"])
    x = [r["num_chordwise"] for r in filtered]
    return x, filtered


def _plot_wake_lines(
    ax: plt.Axes,
    records: list[dict],
    panel_ar: int,
    wake_length: int,
    data_key: str,
    value_index: int,
    scale: float = 1.0,
    color: str | tuple = "tab:blue",
    label_prefix: str = "",
) -> None:
    """Plot prescribed and free wake lines on a single axes.

    :param ax: Matplotlib axes to plot on.
    :param records: All parsed cache records.
    :param panel_ar: Panel aspect ratio to filter by.
    :param wake_length: Wake length to filter by.
    :param data_key: Record key to read values from ("coefficients" or "loads").
    :param value_index: Index within the data array to plot.
    :param scale: Multiplicative factor applied to raw values before plotting.
    :param color: Color for both lines. Prescribed uses solid, free uses dashed.
    :param label_prefix: Prefix for legend labels (e.g. "WL=2, ").
    """
    # Prescribed wake.
    x_vals, filtered = _filter_records(records, panel_ar, wake_length, True)
    if filtered:
        y_vals = [r[data_key][value_index] * scale for r in filtered]
        ax.plot(
            x_vals,
            y_vals,
            "o-",
            color=color,
            markersize=5,
            label=f"{label_prefix}Prescribed",
        )

    # Free wake.
    x_vals, filtered = _filter_records(records, panel_ar, wake_length, False)
    if filtered:
        y_vals = [r[data_key][value_index] * scale for r in filtered]
        ax.plot(
            x_vals,
            y_vals,
            "^--",
            color=color,
            markersize=5,
            label=f"{label_prefix}Free",
        )


def _compute_y_limits(
    records: list[dict],
    data_key: str,
    value_index: int,
    scale: float,
) -> tuple[float, float]:
    """Compute robust y-axis limits that always include 0 and exclude outliers.

    Uses the interquartile range (IQR) method: values beyond Q1 - 1.5*IQR or
    Q3 + 1.5*IQR are treated as outliers and excluded from the limits.

    :param records: All parsed cache records.
    :param data_key: Record key to read values from ("coefficients" or "loads").
    :param value_index: Index within the data array to plot.
    :param scale: Multiplicative factor applied to raw values.
    :return: Tuple of (y_min, y_max) with a small padding margin.
    """
    all_y = np.array([r[data_key][value_index] * scale for r in records])

    q1 = np.percentile(all_y, 25)
    q3 = np.percentile(all_y, 75)
    iqr = q3 - q1
    lower_fence = q1 - _IQR_MULTIPLIER * iqr
    upper_fence = q3 + _IQR_MULTIPLIER * iqr

    inliers = all_y[(all_y >= lower_fence) & (all_y <= upper_fence)]
    if len(inliers) == 0:
        inliers = all_y

    y_min = min(0.0, float(np.min(inliers)))
    y_max = max(0.0, float(np.max(inliers)))

    # Add a small margin so points don't sit on the axis edge.
    margin = (y_max - y_min) * 0.05
    return y_min - margin, y_max + margin


def plot_convergence(
    cache_path: Path,
    data_key: str,
    value_index: int,
    display_name: str,
    unit: str,
    scale: float,
) -> None:
    """Create the convergence plot with one subplot per panel aspect ratio.

    All wake lengths are overlaid on each subplot with distinct colors sampled from a
    perceptually uniform sequential colormap. The number of subplots and lines is
    determined entirely by the cache contents. All subplots share the same y-axis scale
    (always including 0), with extreme outliers excluded from the limits.

    :param cache_path: Path to the convergence cache JSON file.
    :param data_key: Record key to read values from ("coefficients" or "loads").
    :param value_index: Index within the data array to plot (0 through 5).
    :param display_name: Human readable name for the plotted quantity.
    :param unit: Unit string to append to y axis labels (e.g. "mgf", "N*m").
    :param scale: Multiplicative factor applied to raw values before plotting.
    """
    data = _load_cache(cache_path)
    records = _parse_cache(data)

    if not records:
        print("No simulation records found in cache.")
        return

    # Discover parameter values from the cache.
    panel_ars = sorted(set(r["panel_ar"] for r in records), reverse=True)
    wake_lengths = sorted(set(r["wake_length"] for r in records))
    all_chordwise = sorted(set(r["num_chordwise"] for r in records))

    # Build a colormap that samples only the middle 75% of the full range.
    full_cmap = plt.get_cmap(WAKE_LENGTH_COLORMAP)
    trimmed_colors = full_cmap(np.linspace(_CMAP_V_MIN, _CMAP_V_MAX, 256))
    cmap = mcolors.LinearSegmentedColormap.from_list("trimmed", trimmed_colors)
    wl_norm = mcolors.Normalize(vmin=min(wake_lengths), vmax=max(wake_lengths))

    # Compute shared y-axis limits (includes 0, excludes outliers).
    y_min, y_max = _compute_y_limits(records, data_key, value_index, scale)

    num_rows = len(panel_ars)
    y_label = f"{display_name} ({unit})" if unit else display_name

    fig, axes = plt.subplots(
        num_rows,
        1,
        figsize=(14, 4 * num_rows),
        squeeze=False,
        sharey=True,
    )
    fig.suptitle(
        f"Convergence of {display_name}  \u2014  {cache_path.stem}",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )

    # Experimental thrust line (only for FX).
    show_experimental = display_name == "FX"
    if show_experimental:
        exp_value = (
            _EXPERIMENTAL_THRUST_MGF
            if unit == "mgf"
            else _EXPERIMENTAL_THRUST_MGF * 9.80665e-6
        )

    for row, panel_ar in enumerate(panel_ars):
        ax = axes[row, 0]

        # Plot experimental band and line first so they sit behind the data.
        if show_experimental:
            ax.axhspan(
                exp_value * 0.9,
                exp_value * 1.1,
                color="red",
                alpha=0.1,
                label=f"\u00b110% range",
            )
            ax.axhline(
                exp_value,
                color="red",
                linewidth=1.5,
                linestyle="-",
                label=f"Experimental ({_EXPERIMENTAL_THRUST_MGF:.0f} mgf)",
            )

        for wake_length in wake_lengths:
            color = cmap(wl_norm(wake_length))
            _plot_wake_lines(
                ax,
                records,
                panel_ar,
                wake_length,
                data_key,
                value_index,
                scale=scale,
                color=color,
                label_prefix=f"WL={wake_length}, ",
            )

        ax.set_title(f"Panel AR = {panel_ar}", fontsize=11)
        ax.set_xlabel("Num Chordwise Panels")
        ax.set_ylabel(y_label)
        ax.grid(True, alpha=0.3)
        ax.set_xticks(all_chordwise)

    # Apply shared y-axis limits (includes 0, outliers excluded).
    axes[0, 0].set_ylim(y_min, y_max)

    # Format y-axis ticks based on the unit.
    if unit in ("N", "N*m"):
        for ax in axes.flat:
            ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2e"))
    elif unit == "mgf":
        for ax in axes.flat:
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.0f}"))

    # Simulation legend: placed below the figure, 4 items per column.
    all_handles, all_labels = axes[0, 0].get_legend_handles_labels()
    if show_experimental:
        # The first 2 handles are the band and line; the rest are simulation data.
        exp_handles = all_handles[:2]
        exp_labels = all_labels[:2]
        sim_handles = all_handles[2:]
        sim_labels = all_labels[2:]
    else:
        sim_handles = all_handles
        sim_labels = all_labels

    sim_ncol = max(1, math.ceil(len(sim_handles) / 4))
    sim_legend = fig.legend(
        sim_handles,
        sim_labels,
        loc="lower center",
        ncol=sim_ncol,
        fontsize=8,
        bbox_to_anchor=(0.5, 0.0),
    )

    # Experimental legend: separate, placed just above the simulation legend.
    if show_experimental:
        fig.legend(
            exp_handles,
            exp_labels,
            loc="lower center",
            ncol=2,
            fontsize=8,
            bbox_to_anchor=(0.5, 0.06),
        )
        fig.add_artist(sim_legend)

    fig.tight_layout(rect=[0.05, 0.14, 0.95, 0.96])

    save_path = cache_path.parent.parent / "convergence_plot.png"
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot to {save_path}")

    plt.show()


if __name__ == "__main__":
    data_key, value_index, display_name, unit, scale = _resolve_quantity(
        QUANTITY, force_unit_mgf=FORCE_UNIT_MGF
    )
    print(f"Plotting {display_name} from {CACHE_PATH.name}...")
    plot_convergence(CACHE_PATH, data_key, value_index, display_name, unit, scale)
