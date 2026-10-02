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


def make_math_svg_fixture() -> str:
    """Makes a fixture that is a minimal SVG, shaped as Matplotlib writes one with
    svg.fonttype set to "none", holding one text element of math and one of Liberation
    Mono.

    The math's text element holds its characters in tspan elements, each naming its own
    face. Its "x" is in the bold italic STIXGeneral face and its "A" and "/" are in the
    regular one, so the SVG uses three faces in all.

    :return: The SVG's text.
    """
    return (
        '<?xml version="1.0" encoding="utf-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" version="1.1">\n'
        " <defs>\n"
        '  <style type="text/css">*{stroke-linejoin: round}</style>\n'
        " </defs>\n"
        " <text>\n"
        '  <tspan style="font-style: italic; font-weight: 700; font-size: 14px; '
        "font-family: 'STIXGeneral'\">x</tspan>\n"
        "  <tspan style=\"font-size: 10px; font-family: 'STIXGeneral'\">A</tspan>\n"
        "  <tspan style=\"font-size: 14px; font-family: 'STIXGeneral'\">/</tspan>\n"
        " </text>\n"
        f" <text style=\"font-family: '{_fonts.MONO_FONT_FAMILY}'\">Cg</text>\n"
        "</svg>\n"
    )


def make_unknown_font_svg_fixture() -> str:
    """Makes a fixture that is a minimal SVG whose one text element names a font family
    that no font file is known for.

    :return: The SVG's text.
    """
    return (
        '<?xml version="1.0" encoding="utf-8" standalone="no"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" version="1.1">\n'
        " <defs>\n"
        '  <style type="text/css">*{stroke-linejoin: round}</style>\n'
        " </defs>\n"
        " <text style=\"font-family: 'No Such Font Family'\">Lift</text>\n"
        "</svg>\n"
    )
