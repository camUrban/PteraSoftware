"""This module contains functions to create solver objects for use in unit tests."""

import pterasoftware as ps

from . import problem_fixtures


def make_steady_horseshoe_solver_fixture():
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


def make_steady_ring_solver_fixture():
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


def make_unsteady_ring_solver_fixture():
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
