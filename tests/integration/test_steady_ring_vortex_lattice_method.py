""" This module is a testing case for the steady ring vortex lattice method solver.

Based on an identical XFLR5 testing case, the expected output for this case is:
    CL:     0.788
    CDi:    0.019
    Cm:     -0.687

Note: The expected output was created using XFLR5's inviscid VLM2 analysis type, which is a ring vortex lattice method
solver.

This module contains the following classes:
    TestRingHorseshoeVortexLatticeMethod: This is a class for testing the steady ring vortex lattice method solver.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import tests.integration.fixtures.solver_fixtures

import pterasoftware as ps


class TestSteadyRingVortexLatticeMethod(unittest.TestCase):
    """ This is a class for testing the steady ring vortex lattice method solver.

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

        # Create the steady method solver.
        self.steady_ring_vortex_lattice_method_validation_solver = (
            tests.integration.fixtures.solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )

    def tearDown(self):
        """ This method tears down the test.

        :return: None
        """

        del self.steady_ring_vortex_lattice_method_validation_solver

    def test_method(self):
        """ This method tests the solver's output.

        :return: None
        """

        # Run the solver.
        self.steady_ring_vortex_lattice_method_validation_solver.run(verbose=True)

        # Calculate the percent errors of the output.
        CDi_expected = 0.019
        CDi_calculated = self.steady_ring_vortex_lattice_method_validation_solver.airplane.total_near_field_force_coefficients_wind_axes[
            0
        ]
        CDi_error = abs(CDi_calculated - CDi_expected) / CDi_expected

        CL_expected = 0.788
        CL_calculated = self.steady_ring_vortex_lattice_method_validation_solver.airplane.total_near_field_force_coefficients_wind_axes[
            2
        ]
        CL_error = abs(CL_calculated - CL_expected) / CL_expected

        Cm_expected = -0.687
        Cm_calculated = self.steady_ring_vortex_lattice_method_validation_solver.airplane.total_near_field_force_coefficients_wind_axes[
            1
        ]
        Cm_error = abs(Cm_calculated - Cm_expected) / Cm_expected

        # Set the allowable percent error.
        allowable_error = 0.10

        ps.output.draw(
            solver=self.steady_ring_vortex_lattice_method_validation_solver,
            show_wake_vortices=False,
            show_streamlines=True,
            show_delta_pressures=True,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(CDi_error < allowable_error)
        self.assertTrue(CL_error < allowable_error)
        self.assertTrue(Cm_error < allowable_error)
