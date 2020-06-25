""" This is a testing case for the unsteady ring vortex lattice method solver with static geometry.

Based on an equivalent XFLR5 testing case, the expected output for this case is:
    CL:     0.588
    CDi:    0.011
    Cm:     -0.197

Note: The expected output was created using XFLR5's inviscid VLM2 analysis type, which is a ring vortex lattice method
solver. The geometry in this case is static. Therefore the results of this unsteady solver should converge to be close
to XFLR5's static result.

This module contains the following classes:
    TestUnsteadyRingVortexLatticeMethodStaticGeometry: This is a class for testing the unsteady ring vortex lattice
                                                       method solver on static geometry.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import unittest

import tests.integration


class TestUnsteadyRingVortexLatticeMethodStaticGeometry(unittest.TestCase):
    """ This is a class for testing the unsteady ring vortex lattice method solver on static geometry.

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
            tests.integration.fixtures.solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

    def tearDown(self):
        """ This method tears down the test.

        :return: None
        """

        del self.unsteady_ring_vortex_lattice_method_validation_solver

    def test_method(self):
        """ This method tests the solver's output.

        :return: None
        """

        # Run the solver.
        self.unsteady_ring_vortex_lattice_method_validation_solver.run(verbose=True)

        # Calculate the percent errors of the output.
        CDi_expected = 0.011
        CDi_calculated = self.unsteady_ring_vortex_lattice_method_validation_solver.current_airplane.total_near_field_force_coefficients_wind_axes[
            0
        ]
        CDi_error = abs(CDi_calculated - CDi_expected) / CDi_expected

        CL_expected = 0.588
        CL_calculated = self.unsteady_ring_vortex_lattice_method_validation_solver.current_airplane.total_near_field_force_coefficients_wind_axes[
            2
        ]
        CL_error = abs(CL_calculated - CL_expected) / CL_expected

        Cm_expected = -0.197
        Cm_calculated = self.unsteady_ring_vortex_lattice_method_validation_solver.current_airplane.total_near_field_force_coefficients_wind_axes[
            1
        ]
        Cm_error = abs(Cm_calculated - Cm_expected) / Cm_expected

        # Set the allowable percent error.
        allowable_error = 0.25

        import aviansoftwareminimumviableproduct as asmvp

        asmvp.output.animate(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show_wake_vortices=True,
            show_delta_pressures=True,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(CDi_error < allowable_error)
        self.assertTrue(CL_error < allowable_error)
        self.assertTrue(Cm_error < allowable_error)
