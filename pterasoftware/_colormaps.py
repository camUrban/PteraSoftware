"""Contains the color maps and color palettes used by the visualization functions."""

import importlib.resources

import matplotlib.colors
import numpy as np

_data_dir = importlib.resources.files("pterasoftware").joinpath("_colormap_data")


def _load_colormap(name: str) -> matplotlib.colors.ListedColormap:
    """Returns a ListedColormap built from one of the color map data files in the
    _colormap_data directory.

    :param name: The name of the color map. It must match a "<name>_rgb.txt" file in the
        _colormap_data directory.
    :return: A ListedColormap built from the named data file's colors.
    """
    with _data_dir.joinpath(name + "_rgb.txt").open("r") as rgb_file:
        rgb = np.loadtxt(rgb_file)
    return matplotlib.colors.ListedColormap(rgb, name=name)


def _load_hex_palette(name: str) -> list[str]:
    """Returns a list of hex color strings built from one of the palette data files in
    the _colormap_data directory.

    :param name: The name of the palette. It must match a "<name>_hex.txt" file in the
        _colormap_data directory.
    :return: A list of hex color strings, in the order they appear in the data file.
    """
    palette_text = _data_dir.joinpath(name + "_hex.txt").read_text()
    return [line.strip() for line in palette_text.splitlines() if line.strip()]


# Use cmocean's "speed" and "delta" color maps. Their colors are vendored in the
# _colormap_data directory so that Ptera Software does not depend on the cmocean
# package. See that directory's CMOCEAN_LICENSE.md.
sequential_color_map = _load_colormap("speed")
diverging_color_map = _load_colormap("delta")

# For the figure lines, use the "Prism" qualitative color map from
# carto.com/carto-colors. Its colors are vendored in the _colormap_data directory. See
# that directory's CARTOCOLORS_LICENSE.md.
prism = _load_hex_palette("prism")
