# ToDo: Properly document this script.
"""

"""

import unittest

import tests.integration.fixtures.solver_fixtures


# ToDo: Properly document this class.
class TestUnsteadyRingVortexLatticeMethodMultipleWingVariableGeometry(
    unittest.TestCase
):
    """ This is a class for testing the unsteady ring vortex lattice method solver on multi-wing, variable geometry.

    This class contains the following public methods:
        setUp: This method sets up the test.
        tearDown: This method tears down the test.
        test_method: This method tests the solver's output.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
        """ This method sets up the test.

        :return: None
        """

        # Create the unsteady method solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver = (
            tests.integration.fixtures.solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_multiple_wing_variable_geometry()
        )

    def tearDown(self):
        """ This method tears down the test.

        :return: None
        """

        del self.unsteady_ring_vortex_lattice_method_validation_solver

    # ToDo: Properly document this method.
    def test_method(self):
        """ This method tests the solver's output.

        :return: None
        """

        # Run the solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(
            verbose=True, prescribed_wake=True,
        )

        import pterasoftware as ps

        ps.output.animate(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show_delta_pressures=True,
            show_wake_vortices=True,
        )
