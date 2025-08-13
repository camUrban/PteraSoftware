"""This is a testing case for the unsteady ring vortex lattice method solver with
variable geometry.

Note: This case does not currently test the solver's output against an expected
output. Instead, it just tests that the solver doesn't throw an error.

This module contains the following classes:
    TestUnsteadyRingVortexLatticeMethodVariableGeometry: This is a class for testing
    the unsteady ring vortex lattice method solver on variable geometry.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestUnsteadyRingVortexLatticeMethodVariableGeometry(unittest.TestCase):
    """This is a class for testing the unsteady ring vortex lattice method solver on
    variable geometry.

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
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry()
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
        # Run the unsteady solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(
            prescribed_wake=True,
        )

        ps.output.animate(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            scalar_type="lift",
            show_wake_vortices=True,
            save=False,
            testing=True,
        )

        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show=False,
            save=False,
        )
