"""This module contains functions to create problem objects for use in tests."""

import pterasoftware as ps

from . import (
    core_movement_fixtures,
    geometry_fixtures,
    movement_fixtures,
    operating_point_fixtures,
)


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


def make_multi_airplane_unsteady_problem_fixture():
    """This method makes a fixture that is an UnsteadyProblem with multiple Airplanes.

    :return multi_airplane_unsteady_problem_fixture: UnsteadyProblem
        This is the UnsteadyProblem with multiple Airplanes.
    """
    # Create a Movement with multiple AirplaneMovements.
    multi_airplane_movement = (
        movement_fixtures.make_movement_with_multiple_airplanes_fixture()
    )

    # Create the UnsteadyProblem with multiple Airplanes.
    multi_airplane_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=multi_airplane_movement,
        only_final_results=False,
    )

    return multi_airplane_unsteady_problem_fixture


def make_basic_coupled_unsteady_problem_fixture():
    """This method makes a fixture that is a _CoupledUnsteadyProblem for general testing.

    :return basic_coupled_unsteady_problem_fixture: _CoupledUnsteadyProblem
        This is the _CoupledUnsteadyProblem configured for general testing.
    """
    # SteadyProblem sets GP1_CgP1 attributes on each Panel exactly once, so a fresh
    # Airplane is required for every _CoupledUnsteadyProblem instance.
    basic_coupled_unsteady_problem_fixture = ps.problems._CoupledUnsteadyProblem(
        movement=core_movement_fixtures.make_static_core_movement_fixture(),
        initial_airplanes=[geometry_fixtures.make_first_airplane_fixture()],
        initial_operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
    )

    return basic_coupled_unsteady_problem_fixture
