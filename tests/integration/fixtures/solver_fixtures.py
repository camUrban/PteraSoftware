# ToDo: Properly document this module.
"""

"""

import aviansoftwareminimumviableproduct as asmvp
import tests.integration


# ToDo: Properly document this method.
def make_steady_horseshoe_vortex_lattice_method_validation_solver():
    """

    :return: 
    """

    steady_validation_problem = (
        tests.integration.fixtures.problem_fixtures.make_steady_validation_problem()
    )

    steady_horseshoe_vortex_lattice_method_validation_solver = asmvp.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
        steady_validation_problem
    )

    del steady_validation_problem

    return steady_horseshoe_vortex_lattice_method_validation_solver


# ToDo: Properly document this method.
def make_steady_ring_vortex_lattice_method_validation_solver():
    """

    :return: 
    """

    steady_validation_problem = (
        tests.integration.fixtures.problem_fixtures.make_steady_validation_problem()
    )

    steady_ring_vortex_lattice_method_validation_solver = asmvp.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver(
        steady_validation_problem
    )

    del steady_validation_problem

    return steady_ring_vortex_lattice_method_validation_solver


# ToDo: Properly document this method.
def make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry():
    """

    :return:
    """

    unsteady_validation_problem = (
        tests.integration.fixtures.problem_fixtures.make_unsteady_validation_problem_with_static_geometry()
    )

    unsteady_ring_vortex_lattice_method_validation_solver = asmvp.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_validation_problem
    )

    del unsteady_validation_problem

    return unsteady_ring_vortex_lattice_method_validation_solver


# ToDo: Properly document this method.
def make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry():
    """

    :return:
    """

    unsteady_validation_problem = (
        tests.integration.fixtures.problem_fixtures.make_unsteady_validation_problem_with_variable_geometry()
    )

    unsteady_ring_vortex_lattice_method_validation_solver = asmvp.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
        unsteady_validation_problem
    )

    del unsteady_validation_problem

    return unsteady_ring_vortex_lattice_method_validation_solver
