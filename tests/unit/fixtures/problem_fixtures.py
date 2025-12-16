"""This module contains functions to create problem objects for use in tests."""

import pterasoftware as ps

from . import geometry_fixtures, movement_fixtures, operating_point_fixtures


def make_basic_steady_problem_fixture():
    """This method makes a fixture that is a SteadyProblem for general testing.

    :return basic_steady_problem_fixture: SteadyProblem
        This is the SteadyProblem configured for general testing.
    """
    # Create an Airplane.
    first_airplane = geometry_fixtures.make_first_airplane_fixture()

    # Create a basic OperatingPoint.
    basic_operating_point = (
        operating_point_fixtures.make_basic_operating_point_fixture()
    )

    # Create the SteadyProblem.
    basic_steady_problem_fixture = ps.problems.SteadyProblem(
        airplanes=[first_airplane],
        operating_point=basic_operating_point,
    )

    return basic_steady_problem_fixture


def make_multi_airplane_steady_problem_fixture():
    """This method makes a fixture that is a SteadyProblem with multiple Airplanes.

    :return multi_airplane_steady_problem_fixture: SteadyProblem
        This is the SteadyProblem with multiple Airplanes.
    """
    # Create multiple Airplanes.
    airplane1 = geometry_fixtures.make_first_airplane_fixture()
    airplane2 = geometry_fixtures.make_basic_airplane_fixture()

    # Create an OperatingPoint.
    operating_point = operating_point_fixtures.make_basic_operating_point_fixture()

    # Create the SteadyProblem with multiple Airplanes.
    multi_airplane_steady_problem_fixture = ps.problems.SteadyProblem(
        airplanes=[airplane1, airplane2],
        operating_point=operating_point,
    )

    return multi_airplane_steady_problem_fixture


def make_basic_unsteady_problem_fixture():
    """This method makes a fixture that is an UnsteadyProblem for general testing.

    :return basic_unsteady_problem_fixture: UnsteadyProblem
        This is the UnsteadyProblem configured for general testing.
    """
    # Create a basic Movement.
    basic_movement = movement_fixtures.make_basic_movement_fixture()

    # Create the UnsteadyProblem.
    basic_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=basic_movement,
        only_final_results=False,
    )

    return basic_unsteady_problem_fixture


def make_only_final_results_unsteady_problem_fixture():
    """This method makes a fixture that is an UnsteadyProblem with only_final_results
    set to True.

    :return only_final_results_unsteady_problem_fixture: UnsteadyProblem
        This is the UnsteadyProblem with only_final_results set to True.
    """
    # Create a basic Movement.
    basic_movement = movement_fixtures.make_basic_movement_fixture()

    # Create the UnsteadyProblem with only_final_results=True.
    only_final_results_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=basic_movement,
        only_final_results=True,
    )

    return only_final_results_unsteady_problem_fixture


def make_static_unsteady_problem_fixture():
    """This method makes a fixture that is an UnsteadyProblem with static Movement.

    :return static_unsteady_problem_fixture: UnsteadyProblem
        This is the UnsteadyProblem with static Movement.
    """
    # Create a static Movement.
    static_movement = movement_fixtures.make_static_movement_fixture()

    # Create the UnsteadyProblem with static Movement.
    static_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=static_movement,
        only_final_results=False,
    )

    return static_unsteady_problem_fixture


def make_cyclic_unsteady_problem_fixture():
    """This method makes a fixture that is an UnsteadyProblem with cyclic Movement.

    :return cyclic_unsteady_problem_fixture: UnsteadyProblem
        This is the UnsteadyProblem with cyclic Movement.
    """
    # Create a cyclic Movement.
    cyclic_movement = movement_fixtures.make_cyclic_movement_fixture()

    # Create the UnsteadyProblem with cyclic Movement.
    cyclic_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=cyclic_movement,
        only_final_results=False,
    )

    return cyclic_unsteady_problem_fixture
