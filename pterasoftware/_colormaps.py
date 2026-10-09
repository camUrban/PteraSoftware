"""Contains the color maps used by the visualization functions."""

import importlib.resources

import matplotlib.colors
import numpy as np

DATA_DIR = importlib.resources.files("pterasoftware").joinpath("_colormap_data")


def load_colormap(name: str) -> matplotlib.colors.ListedColormap:
    """Returns a ListedColormap built from one of the color map data files in the
    _colormap_data directory.

    :param name: The name of the color map. It must match a "<name>_rgb.txt" file in the
        _colormap_data directory.
    :return: A ListedColormap built from the named data file's colors.
    """
    with DATA_DIR.joinpath(name + "_rgb.txt").open("r") as rgb_file:
        rgb = np.loadtxt(rgb_file)
    return matplotlib.colors.ListedColormap(rgb, name=name)


# Use cmocean's "speed" and "delta" color maps. Their colors are vendored in the
# _colormap_data directory so that Ptera Software does not depend on the cmocean
# package. See that directory's CMOCEAN_LICENSE.md.
SEQUENTIAL_COLOR_MAP = load_colormap("speed")
DIVERGING_COLOR_MAP = load_colormap("delta")
