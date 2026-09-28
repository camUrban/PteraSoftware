"""Tests for the fonts module."""

import unittest

import fontTools.ttLib

# noinspection PyProtectedMember
from pterasoftware import _fonts


class TestFont(unittest.TestCase):
    """Tests for the vendored font.

    The expected values pin the vendored copy of Liberation Sans, so a replaced font
    file cannot silently change which font the visualizations name.
    """

    def test_font_path_is_a_file(self) -> None:
        """The font path should point at the vendored font file."""
        self.assertTrue(_fonts.FONT_PATH.is_file())

    def test_font_family_matches_the_font_files_family_name(self) -> None:
        """The font family should be the family name the font file declares."""
        font = fontTools.ttLib.TTFont(_fonts.FONT_PATH)
        self.assertEqual(font["name"].getDebugName(1), _fonts.FONT_FAMILY)

    def test_font_is_the_regular_style(self) -> None:
        """The font file should hold the regular style, which all the text uses."""
        font = fontTools.ttLib.TTFont(_fonts.FONT_PATH)
        self.assertEqual(font["name"].getDebugName(2), "Regular")

    def test_font_is_the_vendored_version(self) -> None:
        """The font file should be the vendored release of Liberation Sans."""
        font = fontTools.ttLib.TTFont(_fonts.FONT_PATH)
        self.assertEqual(font["name"].getDebugName(5), "Version 2.1.5")
