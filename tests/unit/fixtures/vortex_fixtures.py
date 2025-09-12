"""This module contains functions to create vortex objects for use in tests.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_origin_fixture: This method makes a fixture that holds the coordinates of a
    position.

    make_termination_fixture: This method makes a fixture that holds the coordinates
    of a position.

    make_strength_fixture: This method makes a fixture that holds a value.

    make_line_vortex_fixture: This method makes a fixture that is a line vortex object.

    make_infinite_leg_direction_fixture: This method makes a direction fixture.

    make_infinite_leg_length_fixture: This method makes a fixture that holds a value.

    make_horseshoe_vortex_fixture: This method makes a fixture that is a horseshoe
    vortex object.

    make_front_left_vertex_fixture: This method makes a fixture that holds the
    coordinates of a position.

    make_front_right_vertex_fixture: This method makes a fixture that holds the
    coordinates of a position.

    make_back_left_vertex_fixture: This method makes a fixture that holds the
    coordinates of a position.

    make_back_right_vertex_fixture: This method makes a fixture that holds the
    coordinates of a position.

    make_ring_vortex_fixture: This method makes a fixture that is a ring vortex
    object.
"""

import numpy as np

import pterasoftware as ps


def make_origin_fixture():
    """This method makes a fixture that holds the coordinates of a position.

    :return termination: (,3) array of floats
        This holds the origin's coordinates.
    """
    origin = np.zeros(3)

    return origin


def make_termination_fixture():
    """This method makes a fixture that holds the coordinates of a position.

    :return termination: (,3) array of floats
        This holds the termination's coordinates.
    """
    termination = np.ones(3)

    return termination


def make_strength_fixture():
    """This method makes a fixture that holds a value.

    :return strength: float
        This is the strength of the vortex. Its units are meters squared per second.
    """
    strength = 10.0

    return strength


def make_line_vortex_fixture():
    """This method makes a fixture that is a line vortex object.

    :return line_vortex_fixture: LineVortex
        This is the line vortex object.
    """
    # Initialize the origin, termination, and strength.
    origin_fixture = make_origin_fixture()
    termination_fixture = make_termination_fixture()
    strength_fixture = make_strength_fixture()

    # Create the line vortex object.
    line_vortex_fixture = ps.aerodynamics.LineVortex(
        origin=origin_fixture,
        termination=termination_fixture,
        strength=strength_fixture,
    )

    # Delete the constructing fixtures.
    del origin_fixture
    del termination_fixture
    del strength_fixture

    # Return the vortex fixture.
    return line_vortex_fixture


def make_infinite_leg_direction_fixture():
    """This method makes a direction fixture.

    :return direction: (,3) array of floats
        This holds the direction's components.
    """
    direction = np.array([1.0, 0.0, 0.0])

    return direction


def make_infinite_leg_length_fixture():
    """This method makes a fixture that holds a value.

    :return length: float
        This is the length fixture.
    """
    length = 20.0

    return length


def make_horseshoe_vortex_fixture():
    """This method makes a fixture that is a horseshoe vortex object.

    :return horseshoe_vortex_fixture: HorseshoeVortex
        This is the horseshoe vortex object.
    """
    # Initialize the constructing fixtures.
    origin_fixture = make_origin_fixture()
    termination_fixture = make_termination_fixture()
    strength_fixture = make_strength_fixture()
    infinite_leg_direction_fixture = make_infinite_leg_direction_fixture()
    infinite_leg_length_fixture = make_infinite_leg_length_fixture()

    # Create the horseshoe vortex object.
    horseshoe_vortex_fixture = ps.aerodynamics.HorseshoeVortex(
        finite_leg_origin=origin_fixture,
        finite_leg_termination=termination_fixture,
        strength=strength_fixture,
        infinite_leg_direction=infinite_leg_direction_fixture,
        infinite_leg_length=infinite_leg_length_fixture,
    )

    # Delete the constructing fixtures.
    del origin_fixture
    del termination_fixture
    del strength_fixture
    del infinite_leg_direction_fixture
    del infinite_leg_length_fixture

    # Return the horseshoe vortex fixture.
    return horseshoe_vortex_fixture


def make_front_left_vertex_fixture():
    """This method makes a fixture that holds the coordinates of a position.

    :return Flpp_G_Cg: (,3) array of floats
        This holds the fixture's coordinates.
    """
    front_left_vertex = np.zeros(3)

    return front_left_vertex


def make_front_right_vertex_fixture():
    """This method makes a fixture that holds the coordinates of a position.

    :return Frpp_G_Cg: (,3) array of floats
        This holds the fixture's coordinates.
    """
    front_right_vertex = np.ones(3)

    return front_right_vertex


def make_back_left_vertex_fixture():
    """This method makes a fixture that holds the coordinates of a position.

    :return Blpp_G_Cg: (,3) array of floats
        This holds the fixture's coordinates.
    """
    back_left_vertex = np.array([1.0, 1.0, 0.0])

    return back_left_vertex


def make_back_right_vertex_fixture():
    """This method makes a fixture that holds the coordinates of a position.

    :return Brpp_G_Cg: (,3) array of floats
        This holds the fixture's coordinates.
    """
    back_right_vertex = np.array([1.0, 1.0, 0.0])

    return back_right_vertex


def make_ring_vortex_fixture():
    """This method makes a fixture that is a ring vortex object.

    :return ring_vortex_fixture: RingVortex
        This is the ring vortex object.
    """
    # Initialize the constructing fixtures.
    front_left_vertex_fixture = make_front_left_vertex_fixture()
    front_right_vertex_fixture = make_front_right_vertex_fixture()
    back_left_vertex_fixture = make_back_left_vertex_fixture()
    back_right_vertex_fixture = make_back_right_vertex_fixture()
    strength_fixture = make_strength_fixture()

    # Create the ring vortex object.
    ring_vortex_fixture = ps.aerodynamics.RingVortex(
        front_left_vertex=front_left_vertex_fixture,
        front_right_vertex=front_right_vertex_fixture,
        back_left_vertex=back_left_vertex_fixture,
        back_right_vertex=back_right_vertex_fixture,
        strength=strength_fixture,
    )

    # Delete the constructing fixtures.
    del front_left_vertex_fixture
    del front_right_vertex_fixture
    del back_left_vertex_fixture
    del back_right_vertex_fixture
    del strength_fixture

    # Return the ring vortex fixture.
    return ring_vortex_fixture
