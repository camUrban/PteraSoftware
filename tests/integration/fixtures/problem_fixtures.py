# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp
import tests as tests


# ToDo: Properly document this method.
def make_steady_validation_problem():
    """

    :return:
    """

    steady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_steady_validation_airplane()
    )
    steady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )

    steady_validation_problem = asmvp.problems.SteadyProblem(
        airplane=steady_validation_airplane,
        operating_point=steady_validation_operating_point,
    )

    del steady_validation_airplane
    del steady_validation_operating_point

    return steady_validation_problem


# ToDo: Properly document this method.
def make_steady_multiple_wing_validation_problem():
    """

    :return:
    """

    steady_validation_airplane = (
        tests.integration.fixtures.airplane_fixtures.make_multiple_wing_steady_validation_airplane()
    )
    steady_validation_operating_point = (
        tests.integration.fixtures.operating_point_fixtures.make_validation_operating_point()
    )

    steady_validation_problem = asmvp.problems.SteadyProblem(
        airplane=steady_validation_airplane,
        operating_point=steady_validation_operating_point,
    )

    del steady_validation_airplane
    del steady_validation_operating_point

    return steady_validation_problem


# ToDo: Properly document this method.
def make_unsteady_validation_problem_with_static_geometry():
    """

    :return:
    """

    unsteady_validation_movement = (
        tests.integration.fixtures.movement_fixtures.make_static_validation_movement()
    )

    unsteady_validation_problem = asmvp.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    del unsteady_validation_movement

    return unsteady_validation_problem


# ToDo: Properly document this method.
def make_unsteady_validation_problem_with_variable_geometry():
    """

    :return:
    """

    unsteady_validation_movement = (
        tests.integration.fixtures.movement_fixtures.make_variable_validation_movement()
    )

    unsteady_validation_problem = asmvp.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    del unsteady_validation_movement

    return unsteady_validation_problem
