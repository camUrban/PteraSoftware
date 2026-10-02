"""This module contains functions to create fixtures for the fonts tests."""

# noinspection PyProtectedMember
from pterasoftware import _fonts


def make_svg_fixture() -> str:
    """Makes a fixture that is a minimal SVG, shaped as Matplotlib writes one with
    svg.fonttype set to "none".

    Its two text elements hold "Lift" and "Drag", so a test can tell which characters
    the SVG's text uses.

    :return: The SVG's text.
    """
    return (
        '<?xml version="1.0" encoding="utf-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" version="1.1">\n'
        " <defs>\n"
        '  <style type="text/css">*{stroke-linejoin: round}</style>\n'
        " </defs>\n"
        f" <text style=\"font-family: '{_fonts.FONT_FAMILY}'\">Lift</text>\n"
        f" <text style=\"font-family: '{_fonts.FONT_FAMILY}'\">Drag</text>\n"
        "</svg>\n"
    )
