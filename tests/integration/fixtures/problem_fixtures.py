"""This module creates problem objects to be used as fixtures."""

import pterasoftware as ps
from tests.integration.fixtures import (
    airplane_fixtures,
    movement_fixtures,
    operating_point_fixtures,
)


def make_steady_validation_problem():
    """This function creates a SteadyProblem to be used as a fixture.

    :return steady_validation_problem: SteadyProblem
        This is the SteadyProblem fixture.
    """
    steady_validation_airplane = airplane_fixtures.make_steady_validation_airplane()
    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    steady_validation_problem = ps.problems.SteadyProblem(
        airplanes=[steady_validation_airplane],
        operating_point=steady_validation_operating_point,
    )

    return steady_validation_problem


def make_steady_multiple_wing_validation_problem():
    """This function creates a SteadyProblem with multi-wing geometry to be used as a
    fixture.

    :return steady_validation_problem: SteadyProblem
        This is the SteadyProblem fixture.
    """
    steady_validation_airplane = (
        airplane_fixtures.make_multiple_wing_steady_validation_airplane()
    )
    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    steady_validation_problem = ps.problems.SteadyProblem(
        airplanes=[steady_validation_airplane],
        operating_point=steady_validation_operating_point,
    )

    return steady_validation_problem


def make_unsteady_validation_problem_with_static_geometry():
    """This function creates an UnsteadyProblem with static geometry to be used as a
    fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the UnsteadyProblem fixture.
    """
    unsteady_validation_movement = movement_fixtures.make_static_validation_movement()

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_unsteady_validation_problem_with_variable_geometry():
    """This function creates an UnsteadyProblem with variable geometry to be used as
    a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the UnsteadyProblem fixture.
    """
    unsteady_validation_movement = movement_fixtures.make_variable_validation_movement()

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_unsteady_validation_problem_with_multiple_wing_static_geometry():
    """This function creates an UnsteadyProblem with multi-wing, static geometry to
    be used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the UnsteadyProblem fixture.
    """
    unsteady_validation_movement = (
        movement_fixtures.make_multiple_wing_static_validation_movement()
    )

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_unsteady_validation_problem_with_multiple_wing_variable_geometry():
    """This function creates an UnsteadyProblem with multi-wing, variable geometry to
    be used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the UnsteadyProblem fixture.
    """
    unsteady_validation_movement = (
        movement_fixtures.make_multiple_wing_variable_validation_movement()
    )

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_surface_effect_steady_problem():
    """This function creates a SteadyProblem with an image surface for surface effect
    testing.

    :return surface_effect_steady_problem: SteadyProblem
        This is the SteadyProblem fixture.
    """
    surface_effect_airplane = airplane_fixtures.make_surface_effect_airplane()
    surface_effect_operating_point = (
        operating_point_fixtures.make_surface_effect_operating_point()
    )

    surface_effect_steady_problem = ps.problems.SteadyProblem(
        airplanes=[surface_effect_airplane],
        operating_point=surface_effect_operating_point,
    )

    return surface_effect_steady_problem


def make_surface_effect_free_air_steady_problem():
    """This function creates a SteadyProblem without an image surface, for use as a
    free-air baseline in surface effect validation tests.

    :return free_air_steady_problem: SteadyProblem
        This is the SteadyProblem fixture.
    """
    surface_effect_airplane = airplane_fixtures.make_surface_effect_airplane()
    free_air_operating_point = (
        operating_point_fixtures.make_surface_effect_free_air_operating_point()
    )

    free_air_steady_problem = ps.problems.SteadyProblem(
        airplanes=[surface_effect_airplane],
        operating_point=free_air_operating_point,
    )

    return free_air_steady_problem


def make_surface_effect_unsteady_problem():
    """This function creates an UnsteadyProblem with an image surface for surface
    effect testing.

    :return surface_effect_unsteady_problem: UnsteadyProblem
        This is the UnsteadyProblem fixture.
    """
    surface_effect_movement = movement_fixtures.make_surface_effect_static_movement()

    surface_effect_unsteady_problem = ps.problems.UnsteadyProblem(
        movement=surface_effect_movement
    )

    return surface_effect_unsteady_problem


def make_surface_effect_free_air_unsteady_problem():
    """This function creates an UnsteadyProblem without an image surface, for use as
    a free-air baseline in surface effect validation tests.

    :return free_air_unsteady_problem: UnsteadyProblem
        This is the UnsteadyProblem fixture.
    """
    free_air_movement = movement_fixtures.make_surface_effect_free_air_static_movement()

    free_air_unsteady_problem = ps.problems.UnsteadyProblem(movement=free_air_movement)

    return free_air_unsteady_problem
