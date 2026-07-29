"""Contains the functions that draw the Matplotlib figures for the visualizations."""

from __future__ import annotations

import csv
from pathlib import Path

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

# The fraction of the data's span added as padding on each side of the y axis. It is
# three times matplotlib's default so the legend, which sits inside the axes at its
# best-effort position, usually has empty space to land in.
_Y_AXIS_MARGIN = 0.15


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


def csv_headers(labels: list[str], subtitle: str, y_label: str) -> list[str]:
    """Composes one figure's CSV column headers out of the text that figure already
    carries.

    A header is a transformed form of three pieces of a figure: the quantity from the
    legend label, the axes, point, and frame from the subtitle, and the unit from the y
    axis label. Deriving them rather than writing them out a second time is what keeps a
    column and the figure beside it describing the same quantity the same way.

    :param labels: The figure's legend labels, one per series.
    :param subtitle: The figure's subtitle, parenthesized as the figure carries it. Pass
        an empty string for a figure that has none.
    :param y_label: The figure's y axis label.
    :return: A list of column headers, one per label.
    """
    # A y axis label pairs a quantity with an optional unit, as in "Position (m)" or
    # "Force Coefficient".
    if y_label.endswith(")"):
        quantity, _, unit = y_label.rpartition(" (")
        unit = "(" + unit
    else:
        quantity = y_label
        unit = ""

    # A header spells a moment's unit "N*m" where the axis label spells it "N m". The
    # axis label sits beneath a figure that names the quantity, so the looser spelling
    # reads fine there, while a header stands alone.
    unit = unit.replace("N m", "N*m")

    # A subtitle is parenthesized and comma separated, and a header is neither.
    context = subtitle.strip("()").replace(",", "")

    headers = []
    for label in labels:
        # A figure that labels its series by component alone names the quantity in its
        # title. A column has no title, so the y axis label's quantity stands in.
        if label.endswith(" Component"):
            label = label.replace("Component", quantity)
        headers.append(" ".join(piece for piece in (label, context, unit) if piece))
    return headers


def write_time_history_csv(
    times: np.ndarray,
    headers: list[str],
    columns: list[np.ndarray],
    save_path: Path,
) -> None:
    """Writes one time history to a CSV file, with time as the first column.

    The columns arrive already carrying the signs and the selection that the figures
    plot, so a row reads the same way the corresponding figure does. The headers are
    composed by the caller alongside the figures, since they are drawn from the same
    legend labels, titles, subtitles, and y axis labels that the figures use.

    :param times: A (num_steps,) ndarray of floats representing the time, in seconds, at
        each time step.
    :param headers: The column headers, one per column, not counting time.
    :param columns: A list of (num_steps,) ndarrays of floats, one per column, not
        counting time.
    :param save_path: The fully resolved file path to write the CSV to.
    :return: None
    """
    if len(headers) != len(columns):
        raise ValueError("headers and columns must be the same length.")

    # newline="" is what the csv module requires so that it controls the line endings
    # itself rather than having the text layer translate them a second time. The module
    # then writes RFC 4180's CRLF unless told otherwise, so the terminator is set to LF
    # to match every other text file this project writes.
    with open(save_path, "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file, lineterminator="\n")
        writer.writerow(["Time (s)"] + headers)
        for step, time_value in enumerate(times):
            writer.writerow([time_value] + [column[step] for column in columns])


def plot_time_history(
    times: np.ndarray,
    series: list[np.ndarray],
    labels: list[str],
    colors: list[str],
    title: str,
    subtitle: str,
    y_label: str,
    figure_size_in: tuple[float, float],
    save: bool,
    save_path: Path,
    resolution_dpi: float,
) -> None:
    """Plots one time-history figure, which is a set of series that share a y axis and
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
    :param y_label: The figure's y axis label.
    :param figure_size_in: The figure's width and height in inches.
    :param save: Set this to True to save the figure as a PNG.
    :param save_path: The fully resolved file path to save the figure to if save is
        True. The caller composes it, so this function neither knows nor decides how the
        figures are named.
    :param resolution_dpi: The dots per inch at which to save the PNG if save is True.
    :return: None
    """
    figure, axes = plt.subplots(figsize=figure_size_in, layout="constrained")

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

    # Pad the y axis beyond matplotlib's default so the legend usually has empty space
    # to land in.
    axes.margins(y=_Y_AXIS_MARGIN)

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

    # The subtitle centers over the axes, but the title, being a figure-level artist,
    # centers over the figure, whose midpoint sits left of the axes' midpoint because
    # constrained layout widens the left margin to fit the y axis text. One layout pass
    # finds the axes' final position, and the title is then re-centered over it so the
    # two stay aligned. The layout engine leaves the title's x alone on later draws.
    figure.draw_without_rendering()
    axes_position = axes.get_position()
    figure.suptitle(
        title,
        color=_TEXT_COLOR_NORMALIZED,
        x=(axes_position.x0 + axes_position.x1) / 2,
    )

    # Save the figure as a PNG if the user wants to do so.
    if save:
        figure.savefig(save_path, dpi=resolution_dpi)
