"""This module tests the SteadyHorseshoeVortexLatticeMethodSolver's surface effect
implementation by comparing results with and without an image surface.

Ground effect (flying close to a solid surface) produces well-known aerodynamic trends:
the image vortex system reduces downwash on the real airplane, which increases the
effective lift and decreases the induced drag. This test validates that the method of
images implementation produces these physically correct trends.
"""

import unittest

from tests.integration.fixtures import solver_fixtures


class TestSteadyHorseshoeVortexLatticeMethodSurfaceEffect(unittest.TestCase):
    """This is a class for testing the SteadyHorseshoeVortexLatticeMethodSolver's
    surface effect implementation."""

    @classmethod
    def setUpClass(cls):
        """This method sets up the test by creating and running both solvers.

        :return: None
        """
        cls.surface_effect_solver = (
            solver_fixtures.make_steady_horseshoe_vortex_lattice_method_surface_effect_solver()
        )
        cls.free_air_solver = (
            solver_fixtures.make_steady_horseshoe_vortex_lattice_method_free_air_solver()
        )

        cls.surface_effect_solver.run()
        cls.free_air_solver.run()

    def test_surface_effect_increases_lift(self):
        """This method tests that ground effect increases the lift coefficient
        magnitude.

        In ground effect, the image vortex system reduces downwash, which increases
        the effective angle of attack and therefore the lift. The lift coefficient
        (negative of the z-component of forceCoefficients_W) should be larger in
        magnitude with the surface effect than without.

        :return: None
        """
        c_l_surface_effect = -self.surface_effect_solver.airplanes[
            0
        ].forceCoefficients_W[2]
        c_l_free_air = -self.free_air_solver.airplanes[0].forceCoefficients_W[2]

        self.assertGreater(c_l_surface_effect, c_l_free_air)

    def test_surface_effect_decreases_induced_drag(self):
        """This method tests that ground effect decreases the induced drag coefficient
        magnitude.

        In ground effect, the reduced downwash tilts the lift vector forward, which
        decreases the induced drag. The induced drag coefficient (negative of the
        x-component of forceCoefficients_W) should be smaller in magnitude with the
        surface effect than without.

        :return: None
        """
        c_di_surface_effect = -self.surface_effect_solver.airplanes[
            0
        ].forceCoefficients_W[0]
        c_di_free_air = -self.free_air_solver.airplanes[0].forceCoefficients_W[0]

        self.assertLess(c_di_surface_effect, c_di_free_air)
