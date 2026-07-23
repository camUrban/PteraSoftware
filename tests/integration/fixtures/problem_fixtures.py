"""This module creates problem objects to be used as fixtures."""

import numpy as np

import pterasoftware as ps
from tests.integration.fixtures import (
    airplane_fixtures,
    movement_fixtures,
    operating_point_fixtures,
)


def make_steady_validation_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem to be used as a fixture.

    :return steady_validation_problem: SteadyProblem This is the SteadyProblem fixture.
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


def make_steady_formation_validation_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem with two identical, widely separated
    Airplanes to be used as a fixture.

    :return steady_formation_validation_problem: SteadyProblem This is the SteadyProblem
        fixture.
    """
    first_airplane = airplane_fixtures.make_steady_validation_airplane()

    # Offset the second Airplane far enough from the first that their aerodynamic
    # interaction is negligible.
    second_airplane = first_airplane.deep_copy_with_Cg_GP1_CgP1(
        np.array([10.0, 500.0, -20.0])
    )

    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    steady_formation_validation_problem = ps.problems.SteadyProblem(
        airplanes=[first_airplane, second_airplane],
        operating_point=steady_validation_operating_point,
    )

    return steady_formation_validation_problem


def make_steady_multiple_wing_validation_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem with multi-wing geometry to be used as a
    fixture.

    :return steady_validation_problem: SteadyProblem This is the SteadyProblem fixture.
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


def make_edge_defined_steady_validation_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem with an edge-defined Airplane to be used as
    a fixture.

    :return edge_defined_steady_validation_problem: SteadyProblem This is the
        SteadyProblem fixture.
    """
    edge_defined_validation_airplane = (
        airplane_fixtures.make_edge_defined_validation_airplane()
    )
    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    edge_defined_steady_validation_problem = ps.problems.SteadyProblem(
        airplanes=[edge_defined_validation_airplane],
        operating_point=steady_validation_operating_point,
    )

    return edge_defined_steady_validation_problem


def make_mixed_steady_validation_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem with an Airplane holding both a trapezoidal
    Wing and an edge-defined Wing, to be used as a fixture.

    :return mixed_steady_validation_problem: SteadyProblem This is the SteadyProblem
        fixture.
    """
    mixed_validation_airplane = airplane_fixtures.make_mixed_validation_airplane()
    steady_validation_operating_point = (
        operating_point_fixtures.make_validation_operating_point()
    )

    mixed_steady_validation_problem = ps.problems.SteadyProblem(
        airplanes=[mixed_validation_airplane],
        operating_point=steady_validation_operating_point,
    )

    return mixed_steady_validation_problem


def make_unsteady_validation_problem_with_static_geometry() -> (
    ps.problems.UnsteadyProblem
):
    """This function creates an UnsteadyProblem with static geometry to be used as a
    fixture.

    :return unsteady_validation_problem: UnsteadyProblem This is the UnsteadyProblem
        fixture.
    """
    unsteady_validation_movement = movement_fixtures.make_static_validation_movement()

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_edge_defined_unsteady_validation_problem() -> ps.problems.UnsteadyProblem:
    """This function creates an UnsteadyProblem with an edge-defined Airplane and static
    WingCrossSectionMovements, to be used as a fixture.

    :return edge_defined_unsteady_validation_problem: UnsteadyProblem This is the
        UnsteadyProblem fixture.
    """
    edge_defined_validation_movement = (
        movement_fixtures.make_edge_defined_static_validation_movement()
    )

    edge_defined_unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=edge_defined_validation_movement
    )

    return edge_defined_unsteady_validation_problem


def make_edge_defined_non_static_unsteady_validation_problem() -> (
    ps.problems.UnsteadyProblem
):
    """This function creates an UnsteadyProblem with an edge-defined Airplane whose
    WingCrossSectionMovements are not all static, to be used as a fixture for testing
    that edge-defined convergence rejects it.

    :return edge_defined_non_static_unsteady_validation_problem: UnsteadyProblem This is
        the UnsteadyProblem fixture.
    """
    edge_defined_non_static_validation_movement = (
        movement_fixtures.make_edge_defined_non_static_validation_movement()
    )

    edge_defined_non_static_unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=edge_defined_non_static_validation_movement
    )

    return edge_defined_non_static_unsteady_validation_problem


def make_unsteady_validation_problem_with_variable_geometry() -> (
    ps.problems.UnsteadyProblem
):
    """This function creates an UnsteadyProblem with variable geometry to be used as a
    fixture.

    :return unsteady_validation_problem: UnsteadyProblem This is the UnsteadyProblem
        fixture.
    """
    unsteady_validation_movement = movement_fixtures.make_variable_validation_movement()

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_unsteady_validation_problem_with_multiple_wing_static_geometry() -> (
    ps.problems.UnsteadyProblem
):
    """This function creates an UnsteadyProblem with multi-wing, static geometry to be
    used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem This is the UnsteadyProblem
        fixture.
    """
    unsteady_validation_movement = (
        movement_fixtures.make_multiple_wing_static_validation_movement()
    )

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_unsteady_validation_problem_with_multiple_wing_variable_geometry() -> (
    ps.problems.UnsteadyProblem
):
    """This function creates an UnsteadyProblem with multi-wing, variable geometry to be
    used as a fixture.

    :return unsteady_validation_problem: UnsteadyProblem This is the UnsteadyProblem
        fixture.
    """
    unsteady_validation_movement = (
        movement_fixtures.make_multiple_wing_variable_validation_movement()
    )

    unsteady_validation_problem = ps.problems.UnsteadyProblem(
        movement=unsteady_validation_movement
    )

    return unsteady_validation_problem


def make_surface_effect_steady_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem with an image surface for surface effect
    testing.

    :return surface_effect_steady_problem: SteadyProblem This is the SteadyProblem
        fixture.
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


def make_surface_effect_free_air_steady_problem() -> ps.problems.SteadyProblem:
    """This function creates a SteadyProblem without an image surface, for use as a
    free-air baseline in surface effect validation tests.

    :return free_air_steady_problem: SteadyProblem This is the SteadyProblem fixture.
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


def make_surface_effect_unsteady_problem() -> ps.problems.UnsteadyProblem:
    """This function creates an UnsteadyProblem with an image surface for surface effect
    testing.

    :return surface_effect_unsteady_problem: UnsteadyProblem This is the UnsteadyProblem
        fixture.
    """
    surface_effect_movement = movement_fixtures.make_surface_effect_static_movement()

    surface_effect_unsteady_problem = ps.problems.UnsteadyProblem(
        movement=surface_effect_movement
    )

    return surface_effect_unsteady_problem


def make_surface_effect_free_air_unsteady_problem() -> ps.problems.UnsteadyProblem:
    """This function creates an UnsteadyProblem without an image surface, for use as a
    free-air baseline in surface effect validation tests.

    :return free_air_unsteady_problem: UnsteadyProblem This is the UnsteadyProblem
        fixture.
    """
    free_air_movement = movement_fixtures.make_surface_effect_free_air_static_movement()

    free_air_unsteady_problem = ps.problems.UnsteadyProblem(movement=free_air_movement)

    return free_air_unsteady_problem


def make_simple_glider_free_flight_problem() -> ps.problems.FreeFlightUnsteadyProblem:
    """This function creates the simple glider's FreeFlightUnsteadyProblem to be used as
    a fixture.

    The inertia matrix is the one tuned alongside the glider's planform geometry for
    static stability and verified in XFLR5. It is expressed in the first Airplane's body
    axes relative to the first Airplane's center of gravity, which is at the geometry
    origin. The off-diagonal terms are the body-axes products of inertia. No external
    loads are applied, so the glider flies an unpowered glide driven only by its
    aerodynamics, gravity, and inertia.

    :return simple_glider_free_flight_problem: FreeFlightUnsteadyProblem This is the
        simple glider FreeFlightUnsteadyProblem fixture.
    """
    simple_glider_free_flight_movement = (
        movement_fixtures.make_simple_glider_free_flight_movement()
    )

    # Derive the mass from the glider's weight and gravitational field so that the
    # weight == mass * |g_E| consistency holds (the glider's weight and standard gravity
    # are set on the airplane and operating-point fixtures).
    base_airplane = simple_glider_free_flight_movement.airplane_movements[
        0
    ].base_airplane
    base_g_E = (
        simple_glider_free_flight_movement.operating_point_movement.base_operating_point.g_E
    )
    mass = base_airplane.weight / float(np.linalg.norm(base_g_E))

    I_BP1_CgP1 = np.array(
        [
            [155.614, 0.0, -45.658],
            [0.0, 398.513, 0.0],
            [-45.658, 0.0, 508.699],
        ],
        dtype=float,
    )

    simple_glider_free_flight_problem = ps.problems.FreeFlightUnsteadyProblem(
        movement=simple_glider_free_flight_movement,
        mass=mass,
        I_BP1_CgP1=I_BP1_CgP1,
        external_loads_fn=None,
    )

    return simple_glider_free_flight_problem


def make_flapping_free_flight_problem() -> ps.problems.FreeFlightUnsteadyProblem:
    """This function creates the flapping-wing FreeFlightUnsteadyProblem to be used as a
    fixture.

    The mass properties are reused from the simple glider; this airframe is a coupling
    demonstration and is not separately tuned or verified for static stability. The
    inertia matrix is expressed in the first Airplane's body axes relative to the first
    Airplane's center of gravity, which is at the geometry origin. The off-diagonal
    terms are the body-axes products of inertia. No external loads are applied, so the
    airplane is driven only by its flapping aerodynamics, gravity, and inertia.

    :return flapping_free_flight_problem: FreeFlightUnsteadyProblem This is the
        flapping-wing FreeFlightUnsteadyProblem fixture.
    """
    flapping_free_flight_movement = (
        movement_fixtures.make_flapping_free_flight_movement()
    )

    # Derive the mass from the airplane's weight and gravitational field so that the
    # weight == mass * |g_E| consistency holds (the weight and standard gravity are set
    # on the airplane and operating-point fixtures).
    base_airplane = flapping_free_flight_movement.airplane_movements[0].base_airplane
    base_g_E = (
        flapping_free_flight_movement.operating_point_movement.base_operating_point.g_E
    )
    mass = base_airplane.weight / float(np.linalg.norm(base_g_E))

    I_BP1_CgP1 = np.array(
        [
            [155.614, 0.0, -45.658],
            [0.0, 398.513, 0.0],
            [-45.658, 0.0, 508.699],
        ],
        dtype=float,
    )

    flapping_free_flight_problem = ps.problems.FreeFlightUnsteadyProblem(
        movement=flapping_free_flight_movement,
        mass=mass,
        I_BP1_CgP1=I_BP1_CgP1,
        external_loads_fn=None,
    )

    return flapping_free_flight_problem
