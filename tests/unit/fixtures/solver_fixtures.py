"""This module contains functions to create solver objects for use in unit tests."""

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware._coupled_unsteady_ring_vortex_lattice_method import (
    CoupledUnsteadyRingVortexLatticeMethodSolver,
)
from pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method import (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver,
)
from pterasoftware.free_flight_unsteady_ring_vortex_lattice_method import (
    FreeFlightUnsteadyRingVortexLatticeMethodSolver,
)

from . import problem_fixtures


def make_aeroelastic_unsteady_ring_solver_fixture() -> (
    AeroelasticUnsteadyRingVortexLatticeMethodSolver
):
    """This method makes a fixture that is an
    AeroelasticUnsteadyRingVortexLatticeMethodSolver for general testing.

    :return solver: AeroelasticUnsteadyRingVortexLatticeMethodSolver
        This is the AeroelasticUnsteadyRingVortexLatticeMethodSolver fixture.
    """
    aeroelastic_problem = (
        problem_fixtures.make_basic_aeroelastic_unsteady_problem_fixture()
    )

    solver = AeroelasticUnsteadyRingVortexLatticeMethodSolver(aeroelastic_problem)

    return solver


def make_steady_horseshoe_solver_fixture() -> (
    ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver
):
    """This method makes a fixture that is a
    SteadyHorseshoeVortexLatticeMethodSolver for general testing.

    :return solver: SteadyHorseshoeVortexLatticeMethodSolver
        This is the SteadyHorseshoeVortexLatticeMethodSolver fixture.
    """
    steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

    solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_problem
    )

    return solver


def make_steady_ring_solver_fixture() -> (
    ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
):
    """This method makes a fixture that is a SteadyRingVortexLatticeMethodSolver for
    general testing.

    :return solver: SteadyRingVortexLatticeMethodSolver
        This is the SteadyRingVortexLatticeMethodSolver fixture.
    """
    steady_problem = problem_fixtures.make_basic_steady_problem_fixture()

    solver = ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
        steady_problem
    )

    return solver


def make_unsteady_ring_solver_fixture() -> (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
):
    """This method makes a fixture that is an
    UnsteadyRingVortexLatticeMethodSolver for general testing.

    :return solver: UnsteadyRingVortexLatticeMethodSolver
        This is the UnsteadyRingVortexLatticeMethodSolver fixture.
    """
    unsteady_problem = problem_fixtures.make_basic_unsteady_problem_fixture()

    solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_problem
        )
    )

    return solver


def make_coupled_unsteady_ring_solver_fixture() -> (
    CoupledUnsteadyRingVortexLatticeMethodSolver
):
    """This method makes a fixture that is a
    CoupledUnsteadyRingVortexLatticeMethodSolver for general testing.

    :return solver: CoupledUnsteadyRingVortexLatticeMethodSolver
        This is the CoupledUnsteadyRingVortexLatticeMethodSolver fixture.
    """
    coupled_unsteady_problem = (
        problem_fixtures.make_basic_coupled_unsteady_problem_fixture()
    )

    solver = CoupledUnsteadyRingVortexLatticeMethodSolver(coupled_unsteady_problem)

    return solver


def make_free_flight_unsteady_ring_solver_fixture() -> (
    FreeFlightUnsteadyRingVortexLatticeMethodSolver
):
    """This method makes a fixture that is a
    FreeFlightUnsteadyRingVortexLatticeMethodSolver for general testing.

    :return solver: FreeFlightUnsteadyRingVortexLatticeMethodSolver
        This is the FreeFlightUnsteadyRingVortexLatticeMethodSolver fixture.
    """
    free_flight_unsteady_problem = (
        problem_fixtures.make_basic_free_flight_unsteady_problem_fixture()
    )

    solver = FreeFlightUnsteadyRingVortexLatticeMethodSolver(
        free_flight_unsteady_problem
    )

    return solver
