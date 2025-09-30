# NOTE: I haven't yet started refactoring this module.
"""This module is a testing case for the steady ring vortex lattice method solver.

Based on an identical XFLR5 VLM2 testing case, the expected output for this case is:
    CL:     0.784
    CDi:    0.019
    Cm:     -0.678

Note: The expected output was created using XFLR5's inviscid VLM2 analysis type,
which is a ring vortex lattice method solver.

This module contains the following classes:
    TestSteadyRingVortexLatticeMethod: This is a class for testing the steady ring
    vortex lattice method solver.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestSteadyRingVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the steady ring vortex lattice method solver.

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
        """This method sets up the test.

        :return: None
        """

        # Create the steady method solver.
        self.steady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )

    def tearDown(self):
        """This method tears down the test.

        :return: None
        """

        del self.steady_ring_vortex_lattice_method_validation_solver

    def test_method(self):
        """This method tests the solver's output.

        :return: None
        """

        # Run the solver.
        self.steady_ring_vortex_lattice_method_validation_solver.run()

        # Calculate the percent errors of the output.
        c_di_expected = 0.019
        c_di_calculated = (
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].total_near_field_force_coefficients_wind_axes[0]
        )
        c_di_error = abs((c_di_calculated - c_di_expected) / c_di_expected)

        c_l_expected = 0.784
        c_l_calculated = (
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].total_near_field_force_coefficients_wind_axes[2]
        )
        c_l_error = abs((c_l_calculated - c_l_expected) / c_l_expected)

        c_m_expected = -0.678
        c_m_calculated = (
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[
                0
            ].total_near_field_moment_coefficients_wind_axes[1]
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
