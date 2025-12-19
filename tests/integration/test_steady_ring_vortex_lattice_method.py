"""This module is a testing case for the SteadyRingVortexLatticeMethodSolver.

Based on an identical XFLR5 VLM2 testing case, the expected output for this case is:
    CL:     0.784
    CDi:    0.019
    Cm:     -0.678
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestSteadyRingVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the SteadyRingVortexLatticeMethodSolver."""

    @classmethod
    def setUpClass(cls):
        """This method sets up the test.

        :return: None
        """
        cls.steady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )

    def test_method(self):
        """This method tests the SteadyRingVortexLatticeMethodSolver's output. It
        also tests that the solver doesn't throw an error when the draw function is
        called using it.

        :return: None
        """
        self.steady_ring_vortex_lattice_method_validation_solver.run()

        # Calculate the percent errors of the output.
        c_di_expected = 0.019
        c_di_calculated = -(
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].forceCoefficients_W[0]
        )
        c_di_error = abs((c_di_calculated - c_di_expected) / c_di_expected)

        c_l_expected = 0.784
        c_l_calculated = -(
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].forceCoefficients_W[2]
        )
        c_l_error = abs((c_l_calculated - c_l_expected) / c_l_expected)

        c_m_expected = -0.678
        c_m_calculated = (
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].momentCoefficients_W_CgP1[1]
        )
        c_m_error = abs((c_m_calculated - c_m_expected) / c_m_expected)

        # Set the allowable percent error.
        allowable_error = 0.10

        ps.output.draw(
            solver=self.steady_ring_vortex_lattice_method_validation_solver,
            show_wake_vortices=False,
            show_streamlines=True,
            scalar_type="lift",
            testing=True,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(c_di_error < allowable_error)
        self.assertTrue(c_l_error < allowable_error)
        self.assertTrue(c_m_error < allowable_error)

    def test_method_katz(self):
        """This method tests the SteadyRingVortexLatticeMethodSolver using the Katz
        force calculation method. It verifies that the Katz method produces physically
        reasonable results.

        :return: None
        """
        self.steady_ring_vortex_lattice_method_validation_solver.run(
            force_method="Katz"
        )

        # Get the force coefficients calculated using the Katz method.
        c_di_katz = -(
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].forceCoefficients_W[0]
        )
        c_l_katz = -(
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].forceCoefficients_W[2]
        )
        c_m_katz = self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
            0
        ].momentCoefficients_W_CgP1[1]

        # Verify that the Katz method produces physically reasonable results.
        # Lift should be positive for positive angle of attack.
        self.assertGreater(c_l_katz, 0.0)
        # Induced drag should be positive.
        self.assertGreater(c_di_katz, 0.0)
        # Moment should be negative (nose down) for this configuration.
        self.assertLess(c_m_katz, 0.0)

        # Verify that the Katz method produces results of similar magnitude to the
        # Joukowski method.
        c_l_expected = 0.784
        c_di_expected = 0.019
        c_m_expected = -0.678

        # Allow larger tolerance since methods differ slightly.
        allowable_error = 0.30

        c_l_error = abs((c_l_katz - c_l_expected) / c_l_expected)
        c_di_error = abs((c_di_katz - c_di_expected) / c_di_expected)
        c_m_error = abs((c_m_katz - c_m_expected) / c_m_expected)

        self.assertTrue(c_l_error < allowable_error)
        self.assertTrue(c_di_error < allowable_error)
        self.assertTrue(c_m_error < allowable_error)
