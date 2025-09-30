# NOTE: I haven't yet started refactoring this module.
"""This module creates solver objects to be used as fixtures.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    make_steady_horseshoe_vortex_lattice_method_validation_solver: This function
    creates a solver object using the horseshoe vortex lattice method to be used as a
    fixture.

    make_steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver: This
    function creates a solver object with multi-wing geometry using the horseshoe
    vortex lattice method to be used as a fixture.

    make_steady_ring_vortex_lattice_method_validation_solver: This function creates a
    solver object using the ring vortex lattice method to be used as a fixture.

    make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry:
    This function creates a solver object with static geometry using the unsteady
    ring vortex lattice method to be used as a fixture.

    make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry:
    This function creates a solver object with variable geometry using the unsteady
    ring vortex lattice method to be used as a fixture.

    make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_static_geometry:
    This function creates a solver object with multi-wing, static geometry using the
    unsteady ring vortex lattice method to be used as a fixture.

    make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_variable_geometry:
    This function creates a solver object with multi-wing, variable geometry using
    the unsteady ring vortex lattice method to be used as a fixture.
"""

import pterasoftware as ps
from tests.integration.fixtures import problem_fixtures


def make_steady_horseshoe_vortex_lattice_method_validation_solver():
    """This function creates a solver object using the horseshoe vortex lattice
    method to be used as a fixture.

    :return steady_horseshoe_vortex_lattice_method_validation_solver:
    SteadyHorseshoeVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    steady_validation_problem = problem_fixtures.make_steady_validation_problem()

    # Use the problem fixture to create the solver fixture.
    steady_horseshoe_vortex_lattice_method_validation_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_validation_problem
    )

    # Delete the constructing fixture.
    del steady_validation_problem

    # Return the solver.
    return steady_horseshoe_vortex_lattice_method_validation_solver


def make_steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver():
    """This function creates a solver object with multi-wing geometry using the
    horseshoe vortex lattice method to be used as a fixture.

    :return steady_horseshoe_vortex_lattice_method_validation_solver:
    SteadyHorseshoeVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    steady_validation_problem = (
        problem_fixtures.make_steady_multiple_wing_validation_problem()
    )

    # Use the problem fixture to create the solver fixture.
    steady_horseshoe_vortex_lattice_method_validation_solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_validation_problem
    )

    # Delete the constructing fixture.
    del steady_validation_problem

    # Return the solver.
    return steady_horseshoe_vortex_lattice_method_validation_solver


def make_steady_ring_vortex_lattice_method_validation_solver():
    """This function creates a solver object using the ring vortex lattice method to
    be used as a fixture.

    :return steady_ring_vortex_lattice_method_validation_solver:
    SteadyRingVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    steady_validation_problem = problem_fixtures.make_steady_validation_problem()

    # Use the problem fixture to create the solver fixture.
    steady_ring_vortex_lattice_method_validation_solver = (
        ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
            steady_validation_problem
        )
    )

    # Delete the constructing fixture.
    del steady_validation_problem

    # Return the solver.
    return steady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry():
    """This function creates a solver object with static geometry using the unsteady
    ring vortex lattice method to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_static_geometry()
    )

    # Use the problem fixture to create the solver fixture.
    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    # Delete the constructing fixture.
    del unsteady_validation_problem

    # Return the solver.
    return unsteady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry():
    """This function creates a solver object with variable geometry using the
    unsteady ring vortex lattice method to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_variable_geometry()
    )

    # Use the problem fixture to create the solver fixture.
    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    # Delete the constructing fixture.
    del unsteady_validation_problem

    # Return the solver.
    return unsteady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_static_geometry():
    """This function creates a solver object with multi-wing, static geometry using
    the unsteady ring vortex lattice method to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_multiple_wing_static_geometry()
    )

    # Use the problem fixture to create the solver fixture.
    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    # Delete the constructing fixture.
    del unsteady_validation_problem

    # Return the solver.
    return unsteady_ring_vortex_lattice_method_validation_solver


def make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_variable_geometry():
    """This function creates a solver object with multi-wing variable geometry using
    the unsteady ring vortex lattice method to be used as a fixture.

    :return unsteady_ring_vortex_lattice_method_validation_solver:
    UnsteadyRingVortexLatticeMethodSolver
        This is the solver fixture.
    """
    # Create the problem fixture.
    unsteady_validation_problem = (
        problem_fixtures.make_unsteady_validation_problem_with_multiple_wing_variable_geometry()
    )

    # Use the problem fixture to create the solver fixture.
    unsteady_ring_vortex_lattice_method_validation_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_validation_problem
        )
    )

    # Delete the constructing fixture.
    del unsteady_validation_problem

    # Return the solver.
    return unsteady_ring_vortex_lattice_method_validation_solver
