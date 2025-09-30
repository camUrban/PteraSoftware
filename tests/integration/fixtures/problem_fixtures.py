# NOTE: I haven't yet started refactoring this module.
"""This module creates problem objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_steady_validation_problem: This function creates a steady problem object to
    be used as a fixture.

    make_steady_multiple_wing_validation_problem: This function creates a steady
    problem object with multi-wing geometry to be used as a fixture.

    make_unsteady_validation_problem_with_static_geometry: This function creates an
    unsteady problem object with static geometry to be used as a fixture.

    make_unsteady_validation_problem_with_variable_geometry: This function creates an
    unsteady problem object with variable geometry to be used as a fixture.

    make_unsteady_validation_problem_with_multiple_wing_static_geometry: This
    function creates an unsteady problem object with multi-wing, static geometry to
    be used as a fixture.

    make_unsteady_validation_problem_with_multiple_wing_variable_geometry: This
    function creates an unsteady problem object with multi-wing, variable geometry to
    be used as a fixture.
"""

import pterasoftware as ps
from tests.integration.fixtures import airplane_fixtures
from tests.integration.fixtures import movement_fixtures
from tests.integration.fixtures import operating_point_fixtures


def make_steady_validation_problem():
    """This function creates a steady problem object to be used as a fixture.

    :return steady_validation_problem: SteadyProblem
        This is the problem fixture.
    """

    # Create the constructing fixtures.
    steady_validation_airplane = airplane_fixtures.make_steady_validation_airplane()
    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    # Create the problem fixture.
    steady_validation_problem = ps.problems.SteadyProblem(
        airplanes=[steady_validation_airplane],
        operating_point=steady_validation_operating_point,
    )

    # Delete the constructing fixtures.
    del steady_validation_airplane
    del steady_validation_operating_point

    # Return the problem fixture.
    return steady_validation_problem


def make_steady_multiple_wing_validation_problem():
    """This function creates a steady problem object with multi-wing geometry to be
    used as a fixture.

    :return steady_validation_problem: SteadyProblem
        This is the problem fixture.
    """

    # Create the constructing fixtures.
    steady_validation_airplane = (
        airplane_fixtures.make_multiple_wing_steady_validation_airplane()
    )
    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    # Create the problem fixture.
    steady_validation_problem = ps.problems.SteadyProblem(
        airplanes=[steady_validation_airplane],
        operating_point=steady_validation_operating_point,
    )

    # Delete the constructing fixtures.
    del steady_validation_airplane
    del steady_validation_operating_point

    # Return the problem fixture.
    return steady_validation_problem


def make_unsteady_validation_problem_with_static_geometry():
    """This function creates an unsteady problem object with static geometry to be
    used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the problem fixture.
    """

    # Create the constructing fixture.
    unsteady_validation_movement = movement_fixtures.make_static_validation_movement()

    # Create the problem fixture.
    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    # Delete the constructing fixture.
    del unsteady_validation_movement

    # Return the problem fixture.
    return unsteady_validation_problem


def make_unsteady_validation_problem_with_variable_geometry():
    """This function creates an unsteady problem object with variable geometry to be
    used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the problem fixture.
    """

    # Create the constructing fixture.
    unsteady_validation_movement = movement_fixtures.make_variable_validation_movement()

    # Create the problem fixture.
    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    # Delete the constructing fixture.
    del unsteady_validation_movement

    # Return the problem fixture.
    return unsteady_validation_problem


def make_unsteady_validation_problem_with_multiple_wing_static_geometry():
    """This function creates an unsteady problem object with multi-wing,
    static geometry to be used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the problem fixture.
    """

    # Create the constructing fixture.
    unsteady_validation_movement = (
        movement_fixtures.make_multiple_wing_static_validation_movement()
    )

    # Create the problem fixture.
    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    # Delete the constructing fixture.
    del unsteady_validation_movement

    # Return the problem fixture.
    return unsteady_validation_problem


def make_unsteady_validation_problem_with_multiple_wing_variable_geometry():
    """This function creates an unsteady problem object with multi-wing, variable
    geometry to be used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem
        This is the problem fixture.
    """

    # Create the constructing fixture.
    unsteady_validation_movement = (
        movement_fixtures.make_multiple_wing_variable_validation_movement()
    )

    # Create the problem fixture.
    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    # Delete the constructing fixture.
    del unsteady_validation_movement

    # Return the problem fixture.
    return unsteady_validation_problem
