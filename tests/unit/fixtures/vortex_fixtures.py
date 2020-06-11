
# ToDo: Properly document this module.
"""

"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly document this method.
def make_origin_fixture():
    """

    :return:
    """

    return np.zeros(3)


# ToDo: Properly document this method.
def make_termination_fixture():
    """

    :return:
    """

    return np.ones(3)


# ToDo: Properly document this method.
def make_strength_fixture():
    """

    :return:
    """

    return 10


# ToDo: Properly document this method.
def make_line_vortex_fixture():
    """

    :return:
    """

    origin_fixture = make_origin_fixture()
    termination_fixture = make_termination_fixture()
    strength_fixture = make_strength_fixture()

    line_vortex_fixture = asmvp.aerodynamics.LineVortex(
        origin=origin_fixture,
        termination=termination_fixture,
        strength=strength_fixture
    )

    del origin_fixture
    del termination_fixture
    del strength_fixture

    return line_vortex_fixture


# ToDo: Properly document this method.
def make_infinite_leg_direction_fixture():
    """

    :return:
    """

    return np.array([1, 0, 0])


# ToDo: Properly document this method.
def make_infinite_leg_length_fixture():
    """

    :return:
    """

    return 20


# ToDo: Properly document this method.
def make_horseshoe_vortex_fixture():
    """

    :return:
    """

    origin_fixture = make_origin_fixture()
    termination_fixture = make_termination_fixture()
    strength_fixture = make_strength_fixture()
    infinite_leg_direction_fixture = make_infinite_leg_direction_fixture()
    infinite_leg_length_fixture = make_infinite_leg_length_fixture()

    horseshoe_vortex_fixture = asmvp.aerodynamics.HorseshoeVortex(
        finite_leg_origin=origin_fixture,
        finite_leg_termination=termination_fixture,
        strength=strength_fixture,
        infinite_leg_direction=infinite_leg_direction_fixture,
        infinite_leg_length=infinite_leg_length_fixture
    )

    del origin_fixture
    del termination_fixture
    del strength_fixture
    del infinite_leg_direction_fixture
    del infinite_leg_length_fixture

    return horseshoe_vortex_fixture


# ToDo: Properly document this method.
def make_front_left_vertex_fixture():
    """

    :return:
    """

    return np.zeros(3)


# ToDo: Properly document this method.
def make_front_right_vertex_fixture():
    """

    :return:
    """

    return np.ones(3)


# ToDo: Properly document this method.
def make_back_left_vertex_fixture():
    """

    :return:
    """

    return np.array([1, 1, 0])


# ToDo: Properly document this method.
def make_back_right_vertex_fixture():
    """

    :return:
    """

    return np.array([1, 0, 0])


# ToDo: Properly document this method.
def make_ring_vortex_fixture():
    """

    :return:
    """

    front_left_vertex_fixture = make_front_left_vertex_fixture()
    front_right_vertex_fixture = make_front_right_vertex_fixture()
    back_left_vertex_fixture = make_back_left_vertex_fixture()
    back_right_vertex_fixture = make_back_right_vertex_fixture()
    strength_fixture = make_strength_fixture()

    ring_vortex_fixture = asmvp.aerodynamics.RingVortex(
        front_left_vertex=front_left_vertex_fixture,
        front_right_vertex=front_right_vertex_fixture,
        back_left_vertex=back_left_vertex_fixture,
        back_right_vertex=back_right_vertex_fixture,
        strength=strength_fixture
    )

    del front_left_vertex_fixture
    del front_right_vertex_fixture
    del back_left_vertex_fixture
    del back_right_vertex_fixture
    del strength_fixture

    return ring_vortex_fixture
