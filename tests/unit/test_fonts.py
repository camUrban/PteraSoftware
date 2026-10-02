"""Tests for the fonts module."""

import base64
import io
import re
import unittest

import fontTools.ttLib

# noinspection PyProtectedMember
from pterasoftware import _fonts
from tests.unit.fixtures import fonts_fixtures


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


class TestMonoFont(unittest.TestCase):
    """Tests for the vendored monospaced font.

    The expected values pin the vendored copy of Liberation Mono, so a replaced font
    file cannot silently change which font the diagrams' labels name.
    """

    def test_mono_font_path_is_a_file(self) -> None:
        """The monospaced font path should point at the vendored font file."""
        self.assertTrue(_fonts.MONO_FONT_PATH.is_file())

    def test_mono_font_family_matches_the_font_files_family_name(self) -> None:
        """The monospaced font family should be the family name the font file
        declares."""
        font = fontTools.ttLib.TTFont(_fonts.MONO_FONT_PATH)
        self.assertEqual(font["name"].getDebugName(1), _fonts.MONO_FONT_FAMILY)

    def test_mono_font_is_the_regular_style(self) -> None:
        """The monospaced font file should hold the regular style, which the labels
        use."""
        font = fontTools.ttLib.TTFont(_fonts.MONO_FONT_PATH)
        self.assertEqual(font["name"].getDebugName(2), "Regular")

    def test_mono_font_is_the_vendored_version(self) -> None:
        """The monospaced font file should be the vendored release of Liberation
        Mono."""
        font = fontTools.ttLib.TTFont(_fonts.MONO_FONT_PATH)
        self.assertEqual(font["name"].getDebugName(5), "Version 2.1.5")


class TestEmbedFontsInSvg(unittest.TestCase):
    """This class contains methods for testing _fonts.embed_fonts_in_svg."""

    def test_embeds_one_font_face_rule_ahead_of_matplotlibs_style(self) -> None:
        """Test that the font is embedded once, in a style element placed first in the
        definitions."""
        svg = _fonts.embed_fonts_in_svg(fonts_fixtures.make_svg_fixture())
        self.assertEqual(svg.count("@font-face"), 1)
        self.assertLess(svg.index("@font-face"), svg.index("*{stroke-linejoin"))

    def test_names_the_font_by_its_family(self) -> None:
        """Test that the embedded font takes the family name the text elements use."""
        svg = _fonts.embed_fonts_in_svg(fonts_fixtures.make_svg_fixture())
        self.assertIn(f'font-family: "{_fonts.FONT_FAMILY}"', svg)

    def test_subsets_the_font_to_the_characters_the_text_uses(self) -> None:
        """Test that the embedded font maps every character the text uses and drops a
        character it does not."""
        svg = _fonts.embed_fonts_in_svg(fonts_fixtures.make_svg_fixture())
        font_data = re.search(r"base64,([A-Za-z0-9+/=]+)", svg)
        self.assertIsNotNone(font_data)
        assert font_data is not None
        font = fontTools.ttLib.TTFont(io.BytesIO(base64.b64decode(font_data.group(1))))
        character_map = font.getBestCmap()
        for character in "LiftDrag":
            self.assertIn(ord(character), character_map)
        self.assertNotIn(ord("Z"), character_map)

    def test_raises_for_an_svg_without_definitions(self) -> None:
        """Test that an SVG with no defs element is rejected rather than returned
        without its font."""
        svg = fonts_fixtures.make_svg_fixture()
        svg = svg.replace(" <defs>\n", "").replace(" </defs>\n", "")
        with self.assertRaises(ValueError):
            _fonts.embed_fonts_in_svg(svg)

    def test_embeds_each_face_the_text_uses(self) -> None:
        """Test that each face the text uses, including the faces of math written in
        tspan elements, is embedded once, with its family, weight, and style."""
        svg = _fonts.embed_fonts_in_svg(fonts_fixtures.make_math_svg_fixture())
        self.assertEqual(svg.count("@font-face"), 3)
        for rule in (
            'font-family: "STIXGeneral"; font-weight: 700; font-style: italic;',
            'font-family: "STIXGeneral"; font-weight: 400; font-style: normal;',
            f'font-family: "{_fonts.MONO_FONT_FAMILY}"; font-weight: 400; '
            "font-style: normal;",
        ):
            with self.subTest(rule=rule):
                self.assertIn(rule, svg)

    def test_subsets_each_face_to_the_characters_written_in_it(self) -> None:
        """Test that each embedded face maps the characters written in it and drops
        those written only in another face."""
        svg = _fonts.embed_fonts_in_svg(fonts_fixtures.make_math_svg_fixture())
        expected_characters = {
            "font-weight: 700; font-style: italic;": ("x", "A"),
            "font-weight: 400; font-style: normal; src": ("A/", "x"),
        }
        for descriptors, (kept, dropped) in expected_characters.items():
            with self.subTest(descriptors=descriptors):
                font_data = re.search(
                    '"STIXGeneral"; ' + descriptors + r".*?base64,([A-Za-z0-9+/=]+)",
                    svg,
                )
                self.assertIsNotNone(font_data)
                assert font_data is not None
                font = fontTools.ttLib.TTFont(
                    io.BytesIO(base64.b64decode(font_data.group(1)))
                )
                character_map = font.getBestCmap()
                for character in kept:
                    self.assertIn(ord(character), character_map)
                for character in dropped:
                    self.assertNotIn(ord(character), character_map)

    def test_leaves_out_a_face_with_no_known_font_file_with_a_warning(self) -> None:
        """Test that a face no font file is known for isn't embedded, and that a warning
        says so."""
        with self.assertLogs("pterasoftware.output", level="WARNING") as logs:
            svg = _fonts.embed_fonts_in_svg(
                fonts_fixtures.make_unknown_font_svg_fixture()
            )
        self.assertNotIn("@font-face", svg)
        self.assertIn("No Such Font Family", logs.output[0])
