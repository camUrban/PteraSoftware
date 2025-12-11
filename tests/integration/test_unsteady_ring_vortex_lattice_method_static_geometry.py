"""This is a testing case for the UnsteadyRingVortexLatticeMethodSolver with static
geometry.

Based on an equivalent XFLR5 testing case, the expected output for this case is:
    CL:     0.485
    CDi:    0.015
    Cm:     -0.166

Note: The expected output was created using XFLR5's inviscid VLM2 analysis type,
which is a ring vortex lattice method solver. The geometry in this case is static.
Therefore, the results of this unsteady solver should converge to be close to XFLR5's
static result.
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestUnsteadyRingVortexLatticeMethodStaticGeometry(unittest.TestCase):
    """This is a class for testing the UnsteadyRingVortexLatticeMethodSolver on
    static geometry."""

    @classmethod
    def setUpClass(cls):
        """This method sets up the test.

        :return: None
        """
        cls.unsteady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

    def test_method(self):
        """This method tests the UnsteadyRingVortexLatticeMethodSolver's output. It
        also tests that the solver doesn't throw an error when the animate and
        plot_results_versus_time functions are called using it.

        :return: None
        """
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(
            prescribed_wake=True,
        )

        this_solver = self.unsteady_ring_vortex_lattice_method_validation_solver
        this_airplane = this_solver.current_airplanes[0]

        # Calculate the percent errors of the output.
        c_di_expected = 0.015
        c_di_calculated = -this_airplane.forceCoefficients_W[0]
        c_di_error = abs(c_di_calculated - c_di_expected) / c_di_expected

        c_l_expected = 0.485
        c_l_calculated = -this_airplane.forceCoefficients_W[2]
        c_l_error = abs(c_l_calculated - c_l_expected) / c_l_expected

        c_m_expected = -0.166
        c_m_calculated = this_airplane.momentCoefficients_W_CgP1[1]
        c_m_error = abs(c_m_calculated - c_m_expected) / c_m_expected

        # Set the allowable percent error.
        allowable_error = 0.10

        ps.output.animate(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show_wake_vortices=True,
            scalar_type="lift",
            save=False,
            testing=True,
        )

        ps.output.plot_results_versus_time(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show=False,
            save=False,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(abs(c_di_error) < allowable_error)
        self.assertTrue(abs(c_l_error) < allowable_error)
        self.assertTrue(abs(c_m_error) < allowable_error)
