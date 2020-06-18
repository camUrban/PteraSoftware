
# ToDo: Properly document this module.

import unittest
import tests.integration


# ToDo: Properly document this class.
class TestUnsteadyRingVortexLatticeMethod(unittest.TestCase):
    """

    """

    # ToDo: Properly document this method.
    def setUp(self):
        """

        :return:
        """

        self.unsteady_ring_vortex_lattice_method_validation_solver = (
            tests
            .integration
            .fixtures
            .solver_fixtures
            .make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry()
        )

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return:
        """

        del self.unsteady_ring_vortex_lattice_method_validation_solver

    # ToDo: Properly document this method.
    def test_method_does_not_throw(self):
        """

        :return:
        """

        self.unsteady_ring_vortex_lattice_method_validation_solver.run(verbose=True)
