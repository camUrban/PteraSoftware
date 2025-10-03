"""This module creates solver objects to be used as fixtures."""

import pterasoftware as ps
from tests.integration.fixtures import problem_fixtures


def make_steady_horseshoe_vortex_lattice_method_validation_solver():
    """This function creates a SteadyHorseshoeVortexLatticeMethodSolver to be used as
    a fixture.

    :return steady_horseshoe_vortex_lattice_method_validation_solver:
    SteadyHorseshoeVortexLatticeMethodSolver
        This is the SteadyHorseshoeVortexLatticeMethodSolver fixture.
    """
    steady_validation_problem = problem_fixtures.make_steady_validation_problem()

    steady_horseshoe_vortex_lattice_method_validation_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_validation_problem
    )

    return steady_horseshoe_vortex_lattice_method_validation_solver


def make_steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver():
    """This function creates a SteadyHorseshoeVortexLatticeMethodSolver with
    multi-wing geometry to be used as a fixture.

    :return steady_horseshoe_vortex_lattice_method_validation_solver:
    SteadyHorseshoeVortexLatticeMethodSolver
        This is the SteadyHorseshoeVortexLatticeMethodSolver fixture.
    """
    steady_validation_problem = (
        problem_fixtures.make_steady_multiple_wing_validation_problem()
    )

    steady_horseshoe_vortex_lattice_method_validation_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_validation_problem
    )

    return steady_horseshoe_vortex_lattice_method_validation_solver


def make_steady_ring_vortex_lattice_method_validation_solver():
    """This function creates a SteadyRingVortexLatticeMethodSolver to be used as a
    fixture.

    :return steady_ring_vortex_lattice_method_validation_solver:
    SteadyRingVortexLatticeMethodSolver
        This is the SteadyRingVortexLatticeMethodSolver fixture.
    """
    steady_validation_problem = problem_fixtures.make_steady_validation_problem()

    steady_ring_vortex_lattice_method_validation_solver = (
        ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
            steady_validation_problem
        )
    )

    return steady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry():
    """This function creates a UnsteadyRingVortexLatticeMethodSolver with static
    geometry to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the UnsteadyRingVortexLatticeMethodSolver fixture.
    """
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_static_geometry()
    )

    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    return unsteady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry():
    """This function creates a UnsteadyRingVortexLatticeMethodSolver with variable
    geometry to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the UnsteadyRingVortexLatticeMethodSolver fixture.
    """
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_variable_geometry()
    )

    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    return unsteady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_static_geometry():
    """This function creates a UnsteadyRingVortexLatticeMethodSolver with multi-wing,
    static geometry to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the UnsteadyRingVortexLatticeMethodSolver fixture.
    """
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_multiple_wing_static_geometry()
    )

    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    return unsteady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_variable_geometry():
    """This function creates a UnsteadyRingVortexLatticeMethodSolver with multi-wing
    variable geometry to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the UnsteadyRingVortexLatticeMethodSolver fixture.
    """
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_multiple_wing_variable_geometry()
    )

    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    return unsteady_ring_vortex_lattice_method_validation_solver
