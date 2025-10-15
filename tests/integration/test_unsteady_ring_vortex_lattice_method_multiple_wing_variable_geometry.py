"""This is a testing case for the UnsteadyRingVortexLatticeMethodSolver with
variable, multi-wing geometry.

Note: This case does not currently test the solver's output against an expected
output. Instead, it just tests that the solver doesn't throw an error.
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestUnsteadyRingVortexLatticeMethodMultipleWingVariableGeometry(
    unittest.TestCase
):
    """This is a class for testing the UnsteadyRingVortexLatticeMethodSolver on
    multi-wing, variable geometry."""

    def setUp(self):
        """This method sets up the test.

        :return: None
        """
        self.unsteady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_variable_geometry()
        )

    def test_method_does_not_throw(self):
        """This method tests that the UnsteadyRingVortexLatticeMethodSolver doesn't
        throw an error when it runs or when the animate function is called using it.

        :return: None
        """
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(
            prescribed_wake=True,
            logging_level="Critical",
        )

        ps.output.animate(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            scalar_type="induced drag",
            show_wake_vortices=True,
            save=False,
            testing=True,
        )
