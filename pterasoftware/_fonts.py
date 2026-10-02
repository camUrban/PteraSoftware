"""Contains the fonts used by the visualization functions, and the function that embeds
one in an SVG."""

import base64
import importlib.resources
import io
import xml.etree.ElementTree
from pathlib import Path

import fontTools.subset
import fontTools.ttLib

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
        for text_element in root.iter(_SVG_NAMESPACE + "text")
    )

    # The font file carries an FFTM table, which is FontForge's record of when the font
    # was built. The subsetter does not know how to subset it, and it warns before
    # dropping it, so it is dropped outright instead.
    options = fontTools.subset.Options()
    options.drop_tables += ["FFTM"]
    subsetter = fontTools.subset.Subsetter(options)
    subsetter.populate(text=used_text)
    font = fontTools.ttLib.TTFont(FONT_PATH)
    subsetter.subset(font)
    font_buffer = io.BytesIO()
    font.save(font_buffer)
    font_data = base64.b64encode(font_buffer.getvalue()).decode("ascii")

    # Matplotlib always opens an SVG's definitions with a style element of its own, so
    # the font's rule is placed in a style element just ahead of it.
    font_style = (
        f'<style type="text/css">@font-face {{font-family: "{FONT_FAMILY}"; '
        f'src: url(data:font/ttf;base64,{font_data}) format("truetype")}}</style>'
    )
    if "<defs>" not in svg:
        raise ValueError("svg must have a defs element to embed the font in.")
    return svg.replace("<defs>", "<defs>\n  " + font_style, 1)
