"""Tests for the color maps module."""

import re
import unittest

import matplotlib.colors
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _colormaps


class TestLoadColormap(unittest.TestCase):
    """Tests for the _load_colormap function."""

    def test_raises_for_an_unknown_color_map(self):
        """_load_colormap should raise FileNotFoundError for an unknown color map."""
        with self.assertRaises(FileNotFoundError):
            _colormaps._load_colormap("not_a_color_map")


class TestLoadHexPalette(unittest.TestCase):
    """Tests for the _load_hex_palette function."""

    def test_raises_for_an_unknown_palette(self):
        """_load_hex_palette should raise FileNotFoundError for an unknown palette."""
        with self.assertRaises(FileNotFoundError):
            _colormaps._load_hex_palette("not_a_palette")


class TestSequentialColorMap(unittest.TestCase):
    """Tests for the sequential color map.

    The expected values pin the vendored copy of cmocean's "speed" color map. Ptera
    Software no longer depends on cmocean, so nothing upstream will catch it if this
    color map's data file drifts.
    """

    def test_is_a_listed_color_map_with_256_colors(self):
        """The sequential color map should be a ListedColormap with 256 colors."""
        self.assertIsInstance(
            _colormaps.sequential_color_map, matplotlib.colors.ListedColormap
        )
        self.assertEqual(_colormaps.sequential_color_map.N, 256)

    def test_has_the_expected_first_color(self):
        """The sequential color map's first color should be pale yellow."""
        npt.assert_allclose(
            _colormaps.sequential_color_map(0.0),
            (0.9996253193176977, 0.9913711226010461, 0.8041012438578545, 1.0),
        )

    def test_has_the_expected_last_color(self):
        """The sequential color map's last color should be dark green."""
        npt.assert_allclose(
            _colormaps.sequential_color_map(1.0),
            (0.09053276383981979, 0.13733860758438335, 0.07325761429945674, 1.0),
        )


class TestDivergingColorMap(unittest.TestCase):
    """Tests for the diverging color map.

    The expected values pin the vendored copy of cmocean's "delta" color map. Ptera
    Software no longer depends on cmocean, so nothing upstream will catch it if this
    color map's data file drifts.
    """

    def test_is_a_listed_color_map_with_512_colors(self):
        """The diverging color map should be a ListedColormap with 512 colors."""
        self.assertIsInstance(
            _colormaps.diverging_color_map, matplotlib.colors.ListedColormap
        )
        self.assertEqual(_colormaps.diverging_color_map.N, 512)

    def test_has_the_expected_first_color(self):
        """The diverging color map's first color should be dark blue."""
        npt.assert_allclose(
            _colormaps.diverging_color_map(0.0),
            (0.0659773860137986, 0.12386004993819841, 0.24948115997128678, 1.0),
        )

    def test_has_the_expected_last_color(self):
        """The diverging color map's last color should be dark green."""
        npt.assert_allclose(
            _colormaps.diverging_color_map(1.0),
            (0.09053276383981979, 0.13733860758438335, 0.07325761429945674, 1.0),
        )


class TestPrism(unittest.TestCase):
    """Tests for the Prism qualitative color palette.

    The expected values pin the vendored copy of CARTOColors' "Prism" palette.
    """

    def test_has_12_colors(self):
        """The Prism palette should have 12 colors."""
        self.assertEqual(len(_colormaps.prism), 12)

    def test_colors_are_uppercase_hex_strings(self):
        """The Prism palette's colors should all be uppercase hex strings."""
        for color in _colormaps.prism:
            self.assertRegex(color, re.compile(r"^#[0-9A-F]{6}$"))

    def test_has_the_expected_first_color(self):
        """The Prism palette's first color should be purple."""
        self.assertEqual(_colormaps.prism[0], "#5F4690")

    def test_has_the_expected_last_color(self):
        """The Prism palette's last color should be gray."""
        self.assertEqual(_colormaps.prism[-1], "#666666")
