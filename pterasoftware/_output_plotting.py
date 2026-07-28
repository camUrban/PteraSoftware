"""Contains the functions that draw the Matplotlib figures for the visualizations."""

from __future__ import annotations

import matplotlib.legend_handler
import matplotlib.pyplot as plt
import numpy as np

from . import _output_rendering, _transformations
from . import operating_point as operating_point_mod

# Define the colors and line widths used by the results plots. The text color matches
# the one the rendered visualizations use, so the two kinds of output look related.
_FIGURE_BACKGROUND_COLOR = "None"
_TEXT_COLOR_NORMALIZED: tuple[float, float, float] = (
    _output_rendering.TEXT_COLOR[0] / 255,
    _output_rendering.TEXT_COLOR[1] / 255,
    _output_rendering.TEXT_COLOR[2] / 255,
)

# Lines are drawn from thickest to thinnest so that all remain visible even when they
# overlap.
_MAX_LINE_WIDTH = 3.5
_MIN_LINE_WIDTH = 1.5
_LEGEND_LINE_WIDTH = (_MAX_LINE_WIDTH + _MIN_LINE_WIDTH) / 2


def get_operating_point_velocity(
    this_operating_point: operating_point_mod.OperatingPoint,
) -> np.ndarray:
    """Returns the first Airplane's CG velocity (in Earth axes, observed from the Earth
    frame) for a free flight OperatingPoint.

    The CG velocity is the negative of the freestream velocity, since the freestream (a
    still airmass) is entirely due to the first Airplane's motion. The OperatingPoint
    stores the freestream velocity in the first Airplane's geometry axes, so it is
    rotated into Earth axes here.

    :param this_operating_point: The OperatingPoint whose CG velocity will be returned.
    :return: A (3,) ndarray of floats representing the first Airplane's CG velocity (in
        Earth axes, observed from the Earth frame) in meters per second.
    """
    vInf_E__E = _transformations.apply_T_to_vectors(
        this_operating_point.T_pas_GP1_CgP1_to_E_CgP1,
        this_operating_point.vInf_GP1__E,
        is_position=False,
    )
    return -vInf_E__E


def plot_time_history(
    times: np.ndarray,
    series: list[np.ndarray],
    labels: list[str],
    colors: list[str],
    title: str,
    subtitle: str,
    y_label: str,
    save: bool,
    save_name: str,
) -> None:
    """Plots one time-history figure, which is a set of series that share a y-axis and
    are plotted against time.

    Every figure that plot_results_versus_time produces is drawn through this function,
    both the per-Airplane load figures and the free flight state-history figures, so all
    of them share one styling implementation.

    :param times: A (num_steps,) ndarray of floats representing the time, in seconds, at
        each time step.
    :param series: A list of (num_steps,) ndarrays of floats, one per line to plot.
    :param labels: A list of the legend labels, one per series.
    :param colors: A list of the line colors, one per series.
    :param title: The figure's title.
    :param subtitle: A smaller line below the title describing the axes, points, and
        frames of the plotted quantity. Pass an empty string to omit.
    :param y_label: The figure's y-axis label.
    :param save: Set this to True to save the figure as a PNG.
    :param save_name: The file name to save the figure under if save is True.
    :return: None
    """
    figure, axes = plt.subplots()

    # Remove the plot's top and right spines.
    axes.spines.right.set_visible(False)
    axes.spines.top.set_visible(False)

    # Format the plot's spine and label colors.
    axes.spines.bottom.set_color(_TEXT_COLOR_NORMALIZED)
    axes.spines.left.set_color(_TEXT_COLOR_NORMALIZED)
    axes.xaxis.label.set_color(_TEXT_COLOR_NORMALIZED)
    axes.yaxis.label.set_color(_TEXT_COLOR_NORMALIZED)

    # Format the plot's tick colors.
    axes.tick_params(axis="x", colors=_TEXT_COLOR_NORMALIZED)
    axes.tick_params(axis="y", colors=_TEXT_COLOR_NORMALIZED)

    # Format the plot's background colors.
    figure.patch.set_facecolor(_FIGURE_BACKGROUND_COLOR)
    axes.set_facecolor(_FIGURE_BACKGROUND_COLOR)

    # Populate the plot. Lines are drawn from thickest to thinnest so that all remain
    # visible even when the curves overlap.
    num_series = len(series)
    widths = np.linspace(_MAX_LINE_WIDTH, _MIN_LINE_WIDTH, num_series)
    for series_id, (this_series, label, color) in enumerate(
        zip(series, labels, colors)
    ):
        axes.plot(
            times,
            this_series,
            label=label,
            color=color,
            linewidth=widths[series_id],
            solid_capstyle="butt",
        )

    # Name the plot's axis labels, title, and subtitle.
    axes.set_xlabel("Time (s)", color=_TEXT_COLOR_NORMALIZED)
    axes.set_ylabel(y_label, color=_TEXT_COLOR_NORMALIZED)
    figure.suptitle(title, color=_TEXT_COLOR_NORMALIZED)
    if subtitle:
        axes.set_title(subtitle, color=_TEXT_COLOR_NORMALIZED, fontsize="small")

    # Format the plot's legend.
    axes.legend(
        facecolor=_FIGURE_BACKGROUND_COLOR,
        edgecolor=_FIGURE_BACKGROUND_COLOR,
        labelcolor=_TEXT_COLOR_NORMALIZED,
        handler_map={
            plt.Line2D: matplotlib.legend_handler.HandlerLine2D(
                update_func=lambda h, orig: (
                    h.update_from(orig),
                    h.set_linewidth(_LEGEND_LINE_WIDTH),
                )
            )
        },
    )

    # Save the figure as a PNG if the user wants to do so.
    if save:
        figure.savefig(save_name, dpi=300)
