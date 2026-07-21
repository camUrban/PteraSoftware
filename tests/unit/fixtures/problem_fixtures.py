"""This module contains functions to create problem objects for use in tests."""

from collections.abc import Callable

import numpy as np

import pterasoftware as ps

from . import (
    airplane_movement_fixtures,
    core_movement_fixtures,
    geometry_fixtures,
    movement_fixtures,
    operating_point_fixtures,
)


def make_basic_aeroelastic_unsteady_problem_fixture() -> (
    ps.problems.AeroelasticUnsteadyProblem
):
    """This method makes a fixture that is an AeroelasticUnsteadyProblem for testing.

    :return basic_aeroelastic_unsteady_problem_fixture: AeroelasticUnsteadyProblem This
        is the AeroelasticUnsteadyProblem configured for general testing.
    """
    # Create the AeroelasticMovement.
    aeroelastic_movement = movement_fixtures.make_basic_aeroelastic_movement_fixture()

    # Create and return the AeroelasticUnsteadyProblem.
    basic_aeroelastic_unsteady_problem_fixture = ps.problems.AeroelasticUnsteadyProblem(
        movement=aeroelastic_movement,
        wing_density=0.01,
        spring_constant_rad=10.0,
        damping_constant_rad=0.5,
    )

    return basic_aeroelastic_unsteady_problem_fixture


def make_aeroelastic_unsteady_problem_with_standard_wing_fixture() -> (
    ps.problems.AeroelasticUnsteadyProblem
):
    """This method makes a fixture that is an AeroelasticUnsteadyProblem whose Wing is
    backed by a standard WingMovement rather than an AeroelasticWingMovement.

    :return aeroelastic_unsteady_problem_with_standard_wing_fixture:
        AeroelasticUnsteadyProblem This is the AeroelasticUnsteadyProblem configured for
        testing the non-aeroelastic Wing code path.
    """
    # Create the AeroelasticMovement with a standard WingMovement child.
    movement = movement_fixtures.make_aeroelastic_movement_with_standard_wing_fixture()

    # Create and return the AeroelasticUnsteadyProblem.
    aeroelastic_unsteady_problem_with_standard_wing_fixture = (
        ps.problems.AeroelasticUnsteadyProblem(
            movement=movement,
            wing_density=0.01,
            spring_constant_rad=10.0,
            damping_constant_rad=0.5,
        )
    )

    return aeroelastic_unsteady_problem_with_standard_wing_fixture


def make_basic_steady_problem_fixture() -> ps.problems.SteadyProblem:
    """This method makes a fixture that is a SteadyProblem for general testing.

    :return basic_steady_problem_fixture: SteadyProblem This is the SteadyProblem
        configured for general testing.
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


def make_multi_airplane_steady_problem_fixture() -> ps.problems.SteadyProblem:
    """This method makes a fixture that is a SteadyProblem with multiple Airplanes.

    :return multi_airplane_steady_problem_fixture: SteadyProblem This is the
        SteadyProblem with multiple Airplanes.
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


def make_basic_unsteady_problem_fixture() -> ps.problems.UnsteadyProblem:
    """This method makes a fixture that is an UnsteadyProblem for general testing.

    :return basic_unsteady_problem_fixture: UnsteadyProblem This is the UnsteadyProblem
        configured for general testing.
    """
    # Create a basic Movement.
    basic_movement = movement_fixtures.make_basic_movement_fixture()

    # Create the UnsteadyProblem.
    basic_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=basic_movement,
        only_final_results=False,
    )

    return basic_unsteady_problem_fixture


def make_only_final_results_unsteady_problem_fixture() -> ps.problems.UnsteadyProblem:
    """This method makes a fixture that is an UnsteadyProblem with only_final_results
    set to True.

    :return only_final_results_unsteady_problem_fixture: UnsteadyProblem This is the
        UnsteadyProblem with only_final_results set to True.
    """
    # Create a basic Movement.
    basic_movement = movement_fixtures.make_basic_movement_fixture()

    # Create the UnsteadyProblem with only_final_results=True.
    only_final_results_unsteady_problem_fixture = ps.problems.UnsteadyProblem(
        movement=basic_movement,
        only_final_results=True,
    )

    return only_final_results_unsteady_problem_fixture


def make_multi_airplane_unsteady_problem_fixture() -> ps.problems.UnsteadyProblem:
    """This method makes a fixture that is an UnsteadyProblem with multiple Airplanes.

    :return multi_airplane_unsteady_problem_fixture: UnsteadyProblem This is the
        UnsteadyProblem with multiple Airplanes.
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


def make_with_body_rates_steady_problem_fixture() -> ps.problems.SteadyProblem:
    """This method makes a fixture that is a SteadyProblem with a non zero omegas_BP1__E
    for testing that the steady solvers reject body rotation.

    :return with_body_rates_steady_problem_fixture: SteadyProblem This is the
        SteadyProblem whose OperatingPoint has a non zero body angular velocity.
    """
    return ps.problems.SteadyProblem(
        airplanes=[geometry_fixtures.make_first_airplane_fixture()],
        operating_point=operating_point_fixtures.make_with_body_rates_operating_point_fixture(),
    )


def make_with_body_rates_unsteady_problem_fixture() -> ps.problems.UnsteadyProblem:
    """This method makes a fixture that is an UnsteadyProblem whose Movement carries a
    non zero omegas_BP1__E on its base OperatingPoint, for testing that the unsteady
    solver rejects body rotation.

    :return with_body_rates_unsteady_problem_fixture: UnsteadyProblem This is the
        UnsteadyProblem whose generated per-step OperatingPoints all carry a non zero
        body angular velocity.
    """
    airplane_movement = (
        airplane_movement_fixtures.make_basic_airplane_movement_fixture()
    )
    operating_point_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point_fixtures.make_with_body_rates_operating_point_fixture()
    )
    movement = ps.movements.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_cycles=1,
    )
    return ps.problems.UnsteadyProblem(
        movement=movement,
        only_final_results=False,
    )


def make_basic_free_flight_unsteady_problem_fixture(
    base_operating_point: ps.operating_point.OperatingPoint | None = None,
    external_loads_fn: (
        Callable[
            [
                ps.operating_point.OperatingPoint,
                ps.geometry.airplane.Airplane,
            ],
            tuple[np.ndarray, np.ndarray],
        ]
        | None
    ) = None,
    mujoco_assets: dict[str, bytes] | None = None,
) -> ps.problems.FreeFlightUnsteadyProblem:
    """This method makes a fixture that is a FreeFlightUnsteadyProblem for general
    testing.

    :param base_operating_point: OperatingPoint or None The base OperatingPoint for the
        FreeFlightOperatingPointMovement. If None, a basic OperatingPoint with no body
        rotation is used. The default is None.
    :param external_loads_fn: Callable or None A callable that computes additional
        forces and moments to apply to the Airplane during the simulation. If None, no
        additional loads are applied. The default is None.
    :param mujoco_assets: dict or None A dict mapping virtual filenames to their binary
        contents for the MuJoCo model. If None, no extra assets are provided. The
        default is None.
    :return basic_free_flight_unsteady_problem_fixture: FreeFlightUnsteadyProblem This
        is the FreeFlightUnsteadyProblem configured for general testing.
    """
    # Build the FreeFlightMovement using a static motion configuration.
    base_airplane = geometry_fixtures.make_first_airplane_fixture()
    base_wing = base_airplane.wings[0]

    # Create WingCrossSectionMovements (one per WingCrossSection in the Wing).
    wing_cross_section_movements = []
    for wing_cross_section in base_wing.wing_cross_sections:
        wing_cross_section_movements.append(
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wing_cross_section,
            )
        )

    # Create the WingMovement.
    wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=base_wing,
        wing_cross_section_movements=wing_cross_section_movements,
    )

    # Create the AirplaneMovement.
    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=base_airplane,
        wing_movements=[wing_movement],
    )

    # Create the FreeFlightOperatingPointMovement. Free flight requires the Airplane's
    # weight, the mass, and the gravitational field to agree (weight == mass * |g_E|).
    # The shared OperatingPoint fixtures default to no gravity, so rebuild the
    # OperatingPoint with standard gravity here (preserving every other field) and
    # derive the matching mass from the Airplane's weight.
    if base_operating_point is None:
        base_operating_point = (
            operating_point_fixtures.make_basic_operating_point_fixture()
        )
    base_operating_point = ps.operating_point.OperatingPoint(
        rho=base_operating_point.rho,
        vCg__E=base_operating_point.vCg__E,
        alpha=base_operating_point.alpha,
        beta=base_operating_point.beta,
        angles_E_to_BP1_izyx=base_operating_point.angles_E_to_BP1_izyx,
        CgP1_E_Eo=base_operating_point.CgP1_E_Eo,
        surfaceNormal_E=base_operating_point.surfaceNormal_E,
        surfacePoint_E_Eo=base_operating_point.surfacePoint_E_Eo,
        externalFX_W=base_operating_point.externalFX_W,
        nu=base_operating_point.nu,
        g_E=(0.0, 0.0, 9.80665),
        omegas_BP1__E=base_operating_point.omegas_BP1__E,
    )
    mass = base_airplane.weight / 9.80665
    op_movement = ps.movements.free_flight_operating_point_movement.FreeFlightOperatingPointMovement(
        base_operating_point=base_operating_point,
    )

    # Create the FreeFlightMovement.
    movement = ps.movements.free_flight_movement.FreeFlightMovement(
        airplane_movements=[airplane_movement],
        operating_point_movement=op_movement,
        delta_time=0.01,
        prescribed_num_steps=3,
        free_num_steps=2,
    )

    # Create the FreeFlightUnsteadyProblem. The initial geometry and operating point are
    # derived from the movement at the first time step.
    basic_free_flight_unsteady_problem_fixture = ps.problems.FreeFlightUnsteadyProblem(
        movement=movement,
        mass=mass,
        I_BP1_CgP1=np.diag([1.0, 1.0, 1.0]),
        external_loads_fn=external_loads_fn,
        mujoco_assets=mujoco_assets,
    )

    return basic_free_flight_unsteady_problem_fixture


def make_with_body_rates_free_flight_unsteady_problem_fixture() -> (
    ps.problems.FreeFlightUnsteadyProblem
):
    """This method makes a fixture that is a FreeFlightUnsteadyProblem whose base
    OperatingPoint carries a non zero omegas_BP1__E, for testing that the free-flight
    solver permits body rotation.

    :return with_body_rates_free_flight_unsteady_problem_fixture:
        FreeFlightUnsteadyProblem This is the FreeFlightUnsteadyProblem whose initial
        OperatingPoint carries a non zero body angular velocity.
    """
    return make_basic_free_flight_unsteady_problem_fixture(
        base_operating_point=operating_point_fixtures.make_with_body_rates_operating_point_fixture()
    )


def make_basic_coupled_unsteady_problem_fixture() -> (
    ps.problems._CoupledUnsteadyProblem
):
    """This method makes a fixture that is a _CoupledUnsteadyProblem for general
    testing.

    :return basic_coupled_unsteady_problem_fixture: _CoupledUnsteadyProblem This is the
        _CoupledUnsteadyProblem configured for general testing.
    """
    # SteadyProblem sets GP1_CgP1 attributes on each Panel exactly once, so a fresh
    # Airplane is required for every _CoupledUnsteadyProblem instance.
    # noinspection PyProtectedMember
    basic_coupled_unsteady_problem_fixture = ps.problems._CoupledUnsteadyProblem(
        movement=core_movement_fixtures.make_static_core_movement_fixture(),
        initial_airplanes=[geometry_fixtures.make_first_airplane_fixture()],
        initial_operating_point=operating_point_fixtures.make_basic_operating_point_fixture(),
    )

    return basic_coupled_unsteady_problem_fixture
