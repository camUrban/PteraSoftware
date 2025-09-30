# NOTE: I haven't yet started refactoring this module.
"""This is a testing case for the unsteady ring vortex lattice method solver with
static geometry.

Based on an equivalent XFLR5 testing case, the expected output for this case is:
    CL:     0.485
    CDi:    0.015
    Cm:     -0.166

Note: The expected output was created using XFLR5's inviscid VLM2 analysis type,
which is a ring vortex lattice method solver. The geometry in this case is static.
Therefore, the results of this unsteady solver should converge to be close to XFLR5's
static result.

This module contains the following classes:
    TestUnsteadyRingVortexLatticeMethodStaticGeometry: This is a class for testing
    the unsteady ring vortex lattice method solver on static geometry.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestUnsteadyRingVortexLatticeMethodStaticGeometry(unittest.TestCase):
    """This is a class for testing the unsteady ring vortex lattice method solver on
    static geometry.

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

        # Create the unsteady method solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

    def tearDown(self):
        """This method tears down the test.

        :return: None
        """

        del self.unsteady_ring_vortex_lattice_method_validation_solver

    def test_method(self):
        """This method tests the solver's output.

        :return: None
        """

        # Run the solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(
            prescribed_wake=True,
            logging_level="Critical",
        )

        this_solver = self.unsteady_ring_vortex_lattice_method_validation_solver
        this_airplane = this_solver.current_airplanes[0]

        # Calculate the percent errors of the output.
        c_di_expected = 0.015
        c_di_calculated = this_airplane.total_near_field_force_coefficients_wind_axes[0]
        c_di_error = abs(c_di_calculated - c_di_expected) / c_di_expected

        c_l_expected = 0.485
        c_l_calculated = this_airplane.total_near_field_force_coefficients_wind_axes[2]
        c_l_error = abs(c_l_calculated - c_l_expected) / c_l_expected

        c_m_expected = -0.166
        c_m_calculated = this_airplane.total_near_field_moment_coefficients_wind_axes[1]
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
