"""This module is a testing case for the steady horseshoe vortex lattice method solver.

Based on an identical XFLR5 testing case, the expected output for the single-wing case
is:
    CL:     0.790
    CDi:    0.019
    Cm:     -0.690

Based on an identical XFLR5 testing case, the expected output for the multi-wing case
is:
    CL:     0.524
    CDi:    0.007
    Cm:     -0.350

Note: The expected output was created using XFLR5's inviscid VLM1 analysis type,
which is a horseshoe vortex lattice method solver.

This module contains the following classes:
    TestSteadyHorseshoeVortexLatticeMethod: This is a class for testing the steady
    horseshoe vortex lattice method solver.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""
import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestSteadyHorseshoeVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the steady horseshoe vortex lattice method solver.

    This class contains the following public methods:
        setUp: This method sets up the test.

        tearDown: This method tears down the test.

        test_method: This method tests the solver's output.

        test_method_multiple_wings: This method tests the solver's output with
        multi-wing geometry.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def setUp(self):
        """This method sets up the test.

        :return: None
        """

        # Create the steady method solvers.
        self.steady_horseshoe_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_steady_horseshoe_vortex_lattice_method_validation_solver()
        )
        self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver()
        )

    def tearDown(self):
        """This method tears down the test.

        :return: None
        """

        del self.steady_horseshoe_vortex_lattice_method_validation_solver
        del self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver

    def test_method(self):
        """This method tests the solver's output.

        :return: None
        """
        # Run the solver.
        self.steady_horseshoe_vortex_lattice_method_validation_solver.run()

        # Calculate the percent errors of the output.
        c_di_expected = 0.019
        c_di_calculated = (
            self.steady_horseshoe_vortex_lattice_method_validation_solver.airplanes[
                0
            ].total_near_field_force_coefficients_wind_axes[0]
        )
        c_di_error = abs(c_di_calculated - c_di_expected) / c_di_expected

        c_l_expected = 0.790
        c_l_calculated = (
            self.steady_horseshoe_vortex_lattice_method_validation_solver.airplanes[
                0
            ].total_near_field_force_coefficients_wind_axes[2]
        )
        c_l_error = abs(c_l_calculated - c_l_expected) / c_l_expected

        c_m_expected = -0.690
        c_m_calculated = (
            self.steady_horseshoe_vortex_lattice_method_validation_solver.airplanes[
                0
            ].total_near_field_moment_coefficients_wind_axes[1]
        )
        c_m_error = abs(c_m_calculated - c_m_expected) / c_m_expected

        # Set the allowable percent error.
        allowable_error = 0.10

        ps.output.draw(
            solver=self.steady_horseshoe_vortex_lattice_method_validation_solver,
            show_wake_vortices=False,
            show_streamlines=True,
            show_delta_pressures=True,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(abs(c_di_error) < allowable_error)
        self.assertTrue(abs(c_l_error) < allowable_error)
        self.assertTrue(abs(c_m_error) < allowable_error)

    def test_method_multiple_wings(self):
        """This method tests the solver's output with multi-wing geometry.

        :return: None
        """

        # Run the solver.
        self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver.run()

        # Calculate the percent errors of the output.
        c_di_expected = 0.007
        c_di_calculated = self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver.airplanes[
            0
        ].total_near_field_force_coefficients_wind_axes[
            0
        ]
        c_di_error = abs(c_di_calculated - c_di_expected) / c_di_expected

        c_l_expected = 0.524
        c_l_calculated = self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver.airplanes[
            0
        ].total_near_field_force_coefficients_wind_axes[
            2
        ]
        c_l_error = abs(c_l_calculated - c_l_expected) / c_l_expected

        c_m_expected = -0.350
        c_m_calculated = self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver.airplanes[
            0
        ].total_near_field_moment_coefficients_wind_axes[
            1
        ]
        c_m_error = abs(c_m_calculated - c_m_expected) / c_m_expected

        # Set the allowable percent error.
        allowable_error = 0.10

        ps.output.draw(
            solver=self.steady_multiple_wing_horseshoe_vortex_lattice_method_validation_solver,
            show_delta_pressures=True,
            show_streamlines=True,
            show_wake_vortices=False,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(abs(c_di_error) < allowable_error)
        self.assertTrue(abs(c_l_error) < allowable_error)
        self.assertTrue(abs(c_m_error) < allowable_error)
