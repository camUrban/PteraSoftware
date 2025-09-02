"""This is a testing case for the unsteady ring vortex lattice method solver with
variable, multi-wing geometry.

Note: This case does not currently test the solver's output against an expected
output. Instead, it just tests that the solver doesn't throw an error.

This module contains the following classes:
    TestUnsteadyRingVortexLatticeMethodMultipleWingVariableGeometry: This is a class
    for testing the unsteady ring vortex lattice method solver on variable,
    multi-wing geometry.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestUnsteadyRingVortexLatticeMethodMultipleWingVariableGeometry(
    unittest.TestCase
):
    """This is a class for testing the unsteady ring vortex lattice method solver on
    multi-wing, variable geometry.

    This class contains the following public methods:
        setUp: This method sets up the test.

        tearDown: This method tears down the test.

        test_method_does_not_throw: This method tests that the solver does not throw
        any errors.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
        """This method sets up the test.

        :return: None
        """

        # Create the unsteady method solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_variable_geometry()
        )

    def tearDown(self):
        """This method tears down the test.

        :return: None
        """

        del self.unsteady_ring_vortex_lattice_method_validation_solver

    def test_method_does_not_throw(self):
        """This method tests that the solver does not throw any errors.

        :return: None
        """

        # Run the solver.
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
