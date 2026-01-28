"""This module contains functions to create Panels for use in tests."""

import numpy as np

# noinspection PyProtectedMember
from pterasoftware import _panel


def make_basic_panel_fixture():
    """This method makes a fixture that is a Panel for testing purposes.

    The Panel is a simple rectangular panel with 1.0 m chord and 0.5 m span, lying
    flat in the xy plane.

    :return: A Panel configured for testing.
    """
    basic_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    return basic_panel_fixture


def make_leading_edge_panel_fixture():
    """This method makes a fixture that is a Panel at the leading edge.

    :return: A Panel at the leading edge configured for testing.
    """
    leading_edge_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([0.25, 0.0, 0.0]),
        Brpp_G_Cg=np.array([0.25, 0.5, 0.0]),
        is_leading_edge=True,
        is_trailing_edge=False,
    )

    return leading_edge_panel_fixture


def make_trailing_edge_panel_fixture():
    """This method makes a fixture that is a Panel at the trailing edge.

    :return: A Panel at the trailing edge configured for testing.
    """
    trailing_edge_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.75, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.75, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=True,
    )

    return trailing_edge_panel_fixture


def make_leading_and_trailing_edge_panel_fixture():
    """This method makes a fixture that is a Panel at both leading and trailing edges.

    This represents a Panel on a single chordwise panel Wing.

    :return: A Panel at both leading and trailing edges configured for testing.
    """
    leading_and_trailing_edge_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=True,
        is_trailing_edge=True,
    )

    return leading_and_trailing_edge_panel_fixture


def make_square_panel_fixture():
    """This method makes a fixture that is a square Panel.

    The Panel is 1.0 m x 1.0 m, lying flat in the xy plane.

    :return: A square Panel configured for testing.
    """
    square_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 1.0, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 1.0, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    return square_panel_fixture


def make_high_aspect_ratio_panel_fixture():
    """This method makes a fixture that is a high aspect ratio Panel.

    The Panel has span 4.0 m and chord 1.0 m.

    :return: A high aspect ratio Panel configured for testing.
    """
    high_aspect_ratio_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 4.0, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 4.0, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    return high_aspect_ratio_panel_fixture


def make_twisted_panel_fixture():
    """This method makes a fixture that is a twisted (non planar) Panel.

    :return: A twisted Panel configured for testing.
    """
    twisted_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 1.0, 0.5]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([2.0, 0.0, -0.5]),
        Brpp_G_Cg=np.array([2.0, 1.0, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    return twisted_panel_fixture


def make_small_panel_fixture():
    """This method makes a fixture that is a very small Panel.

    The Panel is 0.01 m x 0.01 m.

    :return: A very small Panel configured for testing.
    """
    small_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.00, 0.01, 0.0]),
        Flpp_G_Cg=np.array([0.00, 0.00, 0.0]),
        Blpp_G_Cg=np.array([0.01, 0.00, 0.0]),
        Brpp_G_Cg=np.array([0.01, 0.01, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    return small_panel_fixture


def make_panel_with_set_once_attributes_fixture():
    """This method makes a fixture that is a Panel with set once attributes populated.

    The Panel has is_right_edge, is_left_edge, local_chordwise_position, and
    local_spanwise_position set. It does NOT have GP1 positions set.

    :return: A Panel with set once mesh attributes populated.
    """
    panel_with_set_once_attributes_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    # Set the mesh position attributes.
    panel_with_set_once_attributes_fixture.is_right_edge = True
    panel_with_set_once_attributes_fixture.is_left_edge = False
    panel_with_set_once_attributes_fixture.local_chordwise_position = 2
    panel_with_set_once_attributes_fixture.local_spanwise_position = 5

    return panel_with_set_once_attributes_fixture


def make_panel_with_gp1_positions_fixture():
    """This method makes a fixture that is a Panel with GP1 positions set.

    The GP1 positions are offset from the local positions by [10.0, 20.0, 5.0].

    :return: A Panel with GP1 positions set.
    """
    offset = np.array([10.0, 20.0, 5.0])

    panel_with_gp1_positions_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    # Set the GP1 positions with an offset.
    panel_with_gp1_positions_fixture.Frpp_GP1_CgP1 = np.array([0.0, 0.5, 0.0]) + offset
    panel_with_gp1_positions_fixture.Flpp_GP1_CgP1 = np.array([0.0, 0.0, 0.0]) + offset
    panel_with_gp1_positions_fixture.Blpp_GP1_CgP1 = np.array([1.0, 0.0, 0.0]) + offset
    panel_with_gp1_positions_fixture.Brpp_GP1_CgP1 = np.array([1.0, 0.5, 0.0]) + offset

    return panel_with_gp1_positions_fixture


def make_fully_configured_panel_fixture():
    """This method makes a fixture that is a Panel with all set once attributes set.

    This includes mesh position attributes and GP1 positions.

    :return: A Panel with all set once attributes populated.
    """
    offset = np.array([10.0, 20.0, 5.0])

    fully_configured_panel_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=True,
        is_trailing_edge=False,
    )

    # Set mesh position attributes.
    fully_configured_panel_fixture.is_right_edge = False
    fully_configured_panel_fixture.is_left_edge = True
    fully_configured_panel_fixture.local_chordwise_position = 0
    fully_configured_panel_fixture.local_spanwise_position = 3

    # Set GP1 positions.
    fully_configured_panel_fixture.Frpp_GP1_CgP1 = np.array([0.0, 0.5, 0.0]) + offset
    fully_configured_panel_fixture.Flpp_GP1_CgP1 = np.array([0.0, 0.0, 0.0]) + offset
    fully_configured_panel_fixture.Blpp_GP1_CgP1 = np.array([1.0, 0.0, 0.0]) + offset
    fully_configured_panel_fixture.Brpp_GP1_CgP1 = np.array([1.0, 0.5, 0.0]) + offset

    return fully_configured_panel_fixture


def make_panel_with_cached_properties_fixture():
    """This method makes a fixture that is a Panel with cached properties accessed.

    Accessing the cached properties causes them to be computed and stored.

    :return: A Panel with cached properties already computed.
    """
    panel_with_cached_properties_fixture = _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0]),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0]),
        is_leading_edge=False,
        is_trailing_edge=False,
    )

    # Access cached properties to populate them.
    _ = panel_with_cached_properties_fixture.rightLeg_G
    _ = panel_with_cached_properties_fixture.frontLeg_G
    _ = panel_with_cached_properties_fixture.leftLeg_G
    _ = panel_with_cached_properties_fixture.backLeg_G
    _ = panel_with_cached_properties_fixture.Frbvp_G_Cg
    _ = panel_with_cached_properties_fixture.Flbvp_G_Cg
    _ = panel_with_cached_properties_fixture.Cpp_G_Cg
    _ = panel_with_cached_properties_fixture.unitNormal_G
    _ = panel_with_cached_properties_fixture.area
    _ = panel_with_cached_properties_fixture.aspect_ratio

    return panel_with_cached_properties_fixture
