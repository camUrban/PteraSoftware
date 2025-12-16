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
