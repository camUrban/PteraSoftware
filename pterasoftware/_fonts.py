"""Contains the fonts used by the visualization functions, and the function that embeds
them in an SVG."""

import base64
import importlib.resources
import io
import xml.etree.ElementTree
from pathlib import Path

import fontTools.subset
import fontTools.ttLib
import matplotlib.font_manager

from . import _logging

_logger = _logging.get_logger("output")

# Use Liberation Sans for every piece of text in the visualizations other than the
# diagrams' axes and point labels, so the rendered scenes, the plots, and the vector
# exports all share one typeface. It is metric compatible with Arial and Helvetica. The
# font file is vendored in the _font_data directory so that the text looks the same on
# every machine, and so that it can be embedded in the saved files. See that directory's
# LIBERATION_FONTS_LICENSE.md. Matplotlib and VTK both load a font from a path on disk,
# so the resource is resolved to one here. The family name must match the one the font
# file declares.
FONT_FAMILY = "Liberation Sans"
FONT_PATH = Path(
    str(
        importlib.resources.files("pterasoftware")
        .joinpath("_font_data")
        .joinpath("LiberationSans-Regular.ttf")
    )
)

# Use Liberation Mono, from the same release as Liberation Sans, for the diagrams' axes
# and point labels, so they read as the variable names they are written like. It is
# metric compatible with Courier New, and it is licensed and loaded the same way as
# Liberation Sans.
MONO_FONT_FAMILY = "Liberation Mono"
MONO_FONT_PATH = Path(
    str(
        importlib.resources.files("pterasoftware")
        .joinpath("_font_data")
        .joinpath("LiberationMono-Regular.ttf")
    )
)

# Define the namespace that the elements of an SVG are qualified with.
_SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"


def _parse_svg_style(style: str) -> dict[str, str]:
    """Returns the properties set by an SVG element's style attribute.

    :param style: The style attribute's value, such as "font-size: 14px; font-family:
        'STIXGeneral'".
    :return: A dict mapping each property's name to its value, both stripped of
        surrounding whitespace.
    """
    properties: dict[str, str] = {}
    for declaration in style.split(";"):
        name, separator, value = declaration.partition(":")
        if separator:
            properties[name.strip()] = value.strip()
    return properties


def _find_font_path(family: str, weight: int, style: str) -> Path | None:
    """Returns the file of a font face, or None if no file is known for it.

    The vendored fonts' regular faces are always taken from their vendored files, so a
    different font installed under the same family name cannot stand in for them. Any
    other face, such as the STIX faces that Matplotlib writes math in, is taken from the
    fonts Matplotlib knows, which include the ones it ships.

    :param family: The face's family name, such as "STIXGeneral".
    :param weight: The face's numeric weight, such as 400 for regular or 700 for bold.
    :param style: The face's style, such as "normal" or "italic".
    :return: The path of the face's font file, or None if no file is known for it.
    """
    if weight == 400 and style == "normal":
        if family == FONT_FAMILY:
            return FONT_PATH
        if family == MONO_FONT_FAMILY:
            return MONO_FONT_PATH
    for font_entry in matplotlib.font_manager.fontManager.ttflist:
        if (
            font_entry.name == family
            and font_entry.style == style
            and font_entry.weight == weight
        ):
            return Path(font_entry.fname)
    return None


def embed_fonts_in_svg(svg: str) -> str:
    """Returns an SVG with the fonts its text uses embedded in it, each subset to the
    characters written in it.

    Matplotlib writes an SVG's text as text elements that name their font without
    carrying it, so a viewer would otherwise draw them in whatever font it has installed
    under that name. Each face the text uses, which is a family at a weight and a style,
    is embedded as a base64 encoded @font-face rule, so the text keeps its typeface in
    any viewer that supports such rules while staying selectable. Text written as math
    is split into tspan elements, each naming its own face, such as the STIX faces, and
    those faces are embedded too. Only the glyphs each face draws are kept, which keeps
    the file small. A face whose font file isn't known is left out with a warning, so
    its text falls back to whatever the viewer has installed.

    :param svg: The SVG, as Matplotlib writes it with svg.fonttype set to "none".
    :return: The SVG with the subset fonts embedded.
    """
    if "<defs>" not in svg:
        raise ValueError("svg must have a defs element to embed the fonts in.")

    # Gather the characters each face draws. A text element written as math holds its
    # characters in tspan elements, each of which sets its own style over the text
    # element's, while any other text element holds its characters directly.
    used_text_by_face: dict[tuple[str, int, str], str] = {}
    root = xml.etree.ElementTree.fromstring(svg)
    for text_element in root.iter(_SVG_NAMESPACE + "text"):
        text_style = _parse_svg_style(text_element.get("style", ""))
        tspans = list(text_element.iter(_SVG_NAMESPACE + "tspan"))
        if tspans:
            runs = [
                (
                    {**text_style, **_parse_svg_style(tspan.get("style", ""))},
                    tspan.text or "",
                )
                for tspan in tspans
            ]
        else:
            runs = [(text_style, "".join(text_element.itertext()))]
        for run_style, characters in runs:
            # A font-family property lists fallbacks after the family it names, and a
            # font-weight property can be a number or a name, such as "bold".
            family = run_style.get("font-family", "").split(",")[0].strip().strip("'\"")
            weight_text = run_style.get("font-weight", "400")
            weight = (
                int(weight_text)
                if weight_text.isdigit()
                else matplotlib.font_manager.weight_dict[weight_text]
            )
            face = (family, weight, run_style.get("font-style", "normal"))
            used_text_by_face[face] = used_text_by_face.get(face, "") + characters

    font_rules: list[str] = []
    for (family, weight, style), used_text in used_text_by_face.items():
        font_path = _find_font_path(family, weight, style)
        if font_path is None:
            _logger.warning(
                _logging.indent()
                + 'No font file is known for the "%s" face at weight %d and style '
                '"%s", so it is not embedded in the SVG.',
                family,
                weight,
                style,
            )
            continue

        # The vendored font files carry an FFTM table, which is FontForge's record of
        # when the font was built. The subsetter does not know how to subset it, and it
        # warns before dropping it, so it is dropped outright instead.
        options = fontTools.subset.Options()
        options.drop_tables += ["FFTM"]
        subsetter = fontTools.subset.Subsetter(options)
        subsetter.populate(text=used_text)
        font = fontTools.ttLib.TTFont(font_path)
        subsetter.subset(font)
        font_buffer = io.BytesIO()
        font.save(font_buffer)
        font_data = base64.b64encode(font_buffer.getvalue()).decode("ascii")
        font_rules.append(
            f'@font-face {{font-family: "{family}"; font-weight: {weight}; '
            f"font-style: {style}; "
            f'src: url(data:font/ttf;base64,{font_data}) format("truetype")}}'
        )

    # Matplotlib always opens an SVG's definitions with a style element of its own, so
    # the fonts' rules are placed in a style element just ahead of it.
    font_style = '<style type="text/css">' + "\n".join(font_rules) + "</style>"
    return svg.replace("<defs>", "<defs>\n  " + font_style, 1)
