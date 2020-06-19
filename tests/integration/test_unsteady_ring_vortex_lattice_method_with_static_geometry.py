
# ToDo: Properly document this module.
"""This is a testing case for the unsteady ring vortex lattice method solver.

    Based on an equivalent XFLR5 testing case, the expected output for this case is:
        CL:     0.588
        CDi:    0.011
        Cl:     -0.000
        Cm:     -0.197
        Cn:     -0.000

    Note: The expected output was created using XFLR5's inviscid VLM2 analysis type, which is a ring vortex lattice
    method solver. The geometry in this case is static. Therefore the results of this unsteady solver should converge to
    be close to XFLR5's static result.
"""

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
            tests.integration.fixtures.solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

    # ToDo: Properly document this method.
    def tearDown(self):
        """

        :return:
        """

        del self.unsteady_ring_vortex_lattice_method_validation_solver

    # ToDo: Properly document this method.
    def test_method(self):
        """

        :return:
        """

        # Run the solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(verbose=True)

        CDi_expected = 0.011
        CDi_error = abs(self.unsteady_ring_vortex_lattice_method_validation_solver.CDi - CDi_expected) / CDi_expected

        CL_expected = 0.588
        CL_error = abs(self.unsteady_ring_vortex_lattice_method_validation_solver.CL - CL_expected) / CL_expected

        Cl_expected = -0.000
        Cl_error = abs(self.unsteady_ring_vortex_lattice_method_validation_solver.Cl - Cl_expected)

        Cm_expected = -0.197
        Cm_error = abs(self.unsteady_ring_vortex_lattice_method_validation_solver.Cm - Cm_expected) / Cm_expected

        Cn_expected = -0.000
        Cn_error = abs(self.unsteady_ring_vortex_lattice_method_validation_solver.Cn - Cn_expected)

        allowable_error = 0.25

        self.assertTrue(CDi_error < allowable_error)
        self.assertTrue(CL_error < allowable_error)
        self.assertTrue(Cl_error < allowable_error)
        self.assertTrue(Cm_error < allowable_error)
        self.assertTrue(Cn_error < allowable_error)
