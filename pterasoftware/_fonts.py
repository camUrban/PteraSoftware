"""Contains the font used by the visualization functions."""

import importlib.resources
from pathlib import Path

# Use Liberation Sans for every piece of text in the visualizations, so the rendered
# scenes, the plots, and the vector exports all share one typeface. It is metric
# compatible with Arial and Helvetica. The font file is vendored in the _font_data
# directory so that the text looks the same on every machine, and so that it can be
# embedded in the saved files. See that directory's LIBERATION_FONTS_LICENSE.md.
# Matplotlib and VTK both load a font from a path on disk, so the resource is resolved
# to one here. The family name must match the one the font file declares.
FONT_FAMILY = "Liberation Sans"
FONT_PATH = Path(
    str(
        importlib.resources.files("pterasoftware")
        .joinpath("_font_data")
        .joinpath("LiberationSans-Regular.ttf")
    )
)
