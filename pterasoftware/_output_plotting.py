"""Contains the functions that draw the Matplotlib figures for the visualizations."""

from __future__ import annotations

import base64
import csv
import io
import xml.etree.ElementTree
from pathlib import Path

import fontTools.subset
import fontTools.ttLib
import matplotlib
import matplotlib.legend_handler
import matplotlib.pyplot as plt
import matplotlib.text
import numpy as np

from . import _fonts
from . import _operating_point as operating_point_mod
from . import _output_rendering, _transformations

# Define the file formats the results plots can be saved in.
VALID_FILE_FORMATS = ("png", "svg", "pdf")

# Define the namespace that the elements of an SVG are qualified with.
SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"

# Define the colors and line widths used by the results plots. The text color matches
# the one the rendered visualizations use, so the two kinds of output look related.
FIGURE_BACKGROUND_COLOR = "None"
TEXT_COLOR_NORMALIZED: tuple[float, float, float] = (
    _output_rendering.TEXT_COLOR[0] / 255,
    _output_rendering.TEXT_COLOR[1] / 255,
    _output_rendering.TEXT_COLOR[2] / 255,
)

# Lines are drawn from thickest to thinnest so that all remain visible even when they
# overlap. The widths are spread evenly about a middle width, reaching this fraction of
# it above and below, and the legend draws every line at the middle width.
LINE_WIDTH = 2.5
LINE_WIDTH_SPREAD = 0.4

# The fraction of the data's span added as padding on each side of the y axis. It is
# three times matplotlib's default so the legend, which sits inside the axes at its
# best-effort position, usually has empty space to land in.
Y_AXIS_MARGIN = 0.15


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


def embed_font_in_svg(svg: str) -> str:
    """Returns an SVG with the vendored font embedded in it, subset to the characters
    that its text uses.

    Matplotlib writes an SVG's text as text elements that name their font without
    carrying it, so a viewer would otherwise draw them in whatever font it has installed
    under that name. The font is embedded as a base64 encoded @font-face rule, so the
    text keeps its typeface in any viewer that supports such rules while staying
    selectable. Only the glyphs the text uses are kept, which keeps the file small.

    :param svg: The SVG, as Matplotlib writes it with svg.fonttype set to "none".
    :return: The SVG with the subset font embedded.
    """
    root = xml.etree.ElementTree.fromstring(svg)
    used_text = "".join(
        "".join(text_element.itertext())
        for text_element in root.iter(SVG_NAMESPACE + "text")
    )

    # The font file carries an FFTM table, which is FontForge's record of when the font
    # was built. The subsetter does not know how to subset it, and it warns before
    # dropping it, so it is dropped outright instead.
    options = fontTools.subset.Options()
    options.drop_tables += ["FFTM"]
    subsetter = fontTools.subset.Subsetter(options)
    subsetter.populate(text=used_text)
    font = fontTools.ttLib.TTFont(_fonts.FONT_PATH)
    subsetter.subset(font)
    font_buffer = io.BytesIO()
    font.save(font_buffer)
    font_data = base64.b64encode(font_buffer.getvalue()).decode("ascii")

    # Matplotlib always opens an SVG's definitions with a style element of its own, so
    # the font's rule is placed in a style element just ahead of it.
    font_style = (
        f'<style type="text/css">@font-face {{font-family: "{_fonts.FONT_FAMILY}"; '
        f'src: url(data:font/ttf;base64,{font_data}) format("truetype")}}</style>'
    )
    if "<defs>" not in svg:
        raise ValueError("svg must have a defs element to embed the font in.")
    return svg.replace("<defs>", "<defs>\n  " + font_style, 1)


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
    show_titles: bool = True,
    font_size: float | None = None,
    text_color: (
        tuple[float, float, float] | tuple[float, float, float, float] | None
    ) = None,
    line_width: float | None = None,
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
    :param save: Set this to True to save the figure.
    :param save_path: The fully resolved file path to save the figure to if save is
        True. Its suffix picks the file format, and it must be one of the formats in
        VALID_FILE_FORMATS. The caller composes it, so this function neither knows nor
        decides how the figures are named.
    :param resolution_dpi: The dots per inch at which to save the figure if save is
        True. It only affects a PNG, since the vector formats have no resolution.
    :param show_titles: Set this to False to omit the title and subtitle. The default is
        True.
    :param font_size: The size, in points, of every piece of text. Pass None to size
        each piece of text by Matplotlib's defaults, with the subtitle smaller than the
        rest. The default is None.
    :param text_color: The RGB or RGBA color, with components from 0.0 to 1.0, of the
        text, the axis spines, and the ticks. Pass None to use the color the rendered
        visualizations' text uses. The default is None.
    :param line_width: The middle line width, in points. The lines' widths are spread
        evenly from LINE_WIDTH_SPREAD times it above this width to the same amount
        below, and the legend draws every line at this width. Pass None to use
        LINE_WIDTH. The default is None.
    :return: None
    """
    if text_color is None:
        text_color = TEXT_COLOR_NORMALIZED
    if line_width is None:
        line_width = LINE_WIDTH

    figure, axes = plt.subplots(figsize=figure_size_in, layout="constrained")

    # Remove the plot's top and right spines.
    axes.spines.right.set_visible(False)
    axes.spines.top.set_visible(False)

    # Format the plot's spine and label colors.
    axes.spines.bottom.set_color(text_color)
    axes.spines.left.set_color(text_color)
    axes.xaxis.label.set_color(text_color)
    axes.yaxis.label.set_color(text_color)

    # Format the plot's tick colors.
    axes.tick_params(axis="x", colors=text_color)
    axes.tick_params(axis="y", colors=text_color)

    # Format the plot's background colors.
    figure.patch.set_facecolor(FIGURE_BACKGROUND_COLOR)
    axes.set_facecolor(FIGURE_BACKGROUND_COLOR)

    # Populate the plot. Lines are drawn from thickest to thinnest so that all remain
    # visible even when the curves overlap.
    num_series = len(series)
    widths = line_width * np.linspace(
        1.0 + LINE_WIDTH_SPREAD, 1.0 - LINE_WIDTH_SPREAD, num_series
    )
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
    axes.margins(y=Y_AXIS_MARGIN)

    # Name the plot's axis labels, title, and subtitle.
    axes.set_xlabel("Time (s)", color=text_color)
    axes.set_ylabel(y_label, color=text_color)
    if show_titles:
        title_text = figure.suptitle(title, color=text_color)
        if subtitle:
            axes.set_title(subtitle, color=text_color, fontsize="small")

    # Format the plot's legend.
    axes.legend(
        facecolor=FIGURE_BACKGROUND_COLOR,
        edgecolor=FIGURE_BACKGROUND_COLOR,
        labelcolor=text_color,
        handler_map={
            plt.Line2D: matplotlib.legend_handler.HandlerLine2D(
                update_func=lambda h, orig: (
                    h.update_from(orig),
                    h.set_linewidth(line_width),
                )
            )
        },
    )

    # Set every piece of text in the vendored font. The font is selected by its file
    # rather than by its family name, so a different font installed under the same name
    # cannot stand in for it. The family name is set too, since an SVG names its font by
    # family. A tick that a later draw adds copies its text properties from the first
    # tick, so it inherits the font, and the size below, as well.
    for text in figure.findobj(matplotlib.text.Text):
        font_properties = text.get_fontproperties().copy()
        font_properties.set_family(_fonts.FONT_FAMILY)
        font_properties.set_file(_fonts.FONT_PATH)
        if font_size is not None:
            font_properties.set_size(font_size)
        text.set_fontproperties(font_properties)

    # The subtitle centers over the axes, but the title, being a figure-level artist,
    # centers over the figure, whose midpoint sits left of the axes' midpoint because
    # constrained layout widens the left margin to fit the y axis text. One layout pass
    # finds the axes' final position, and the title is then re-centered over it so the
    # two stay aligned. The layout engine leaves the title's x alone on later draws.
    # Re-titling resets the title's size to Matplotlib's default unless one is passed,
    # so the size it already has is passed along.
    if show_titles:
        figure.draw_without_rendering()
        axes_position = axes.get_position()
        figure.suptitle(
            title,
            color=text_color,
            x=(axes_position.x0 + axes_position.x1) / 2,
            fontsize=title_text.get_fontsize(),
        )

    # Save the figure if the user wants to do so, in the format its path's suffix names.
    # The two settings below only affect the vector formats. A PDF embeds the font as
    # TrueType rather than as Type 3, and an SVG writes its text as text rather than as
    # paths, so that it stays selectable once the font is embedded in it.
    if save:
        with matplotlib.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none"}):
            if save_path.suffix == ".svg":
                svg_buffer = io.BytesIO()
                figure.savefig(svg_buffer, format="svg", dpi=resolution_dpi)
                svg = embed_font_in_svg(svg_buffer.getvalue().decode("utf-8"))
                save_path.write_bytes(svg.encode("utf-8"))
            else:
                figure.savefig(save_path, dpi=resolution_dpi)
