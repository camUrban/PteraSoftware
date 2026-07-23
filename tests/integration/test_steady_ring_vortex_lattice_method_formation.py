"""This module is a testing case for the SteadyRingVortexLatticeMethodSolver with a two
Airplane formation.

The formation holds two identical Airplanes, the second offset far enough from the first
that their aerodynamic interaction is negligible. This makes the second Airplane's loads
relative to its own CG physically comparable to the first Airplane's loads.
"""

import unittest

import numpy as np
import numpy.testing as npt

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _transformations
from tests.integration.fixtures import solver_fixtures


class TestSteadyRingVortexLatticeMethodFormation(unittest.TestCase):
    """This is a class for testing the SteadyRingVortexLatticeMethodSolver on a two
    Airplane formation."""

    steady_ring_vortex_lattice_method_formation_solver: (
        ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """This method sets up the test.

        :return: None
        """
        cls.steady_ring_vortex_lattice_method_formation_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_formation_solver()
        )
        cls.steady_ring_vortex_lattice_method_formation_solver.run()

    def test_far_apart_airplanes_report_matching_own_axes_loads(self) -> None:
        """This method tests that the two far apart Airplanes report matching loads in
        their own geometry axes and relative to their own CGs.

        :return: None
        """
        first_airplane, second_airplane = (
            self.steady_ring_vortex_lattice_method_formation_solver.airplanes
        )

        for attribute_name in ("forces_G", "moments_G_Cg", "moments_W_Cg"):
            with self.subTest(attribute_name=attribute_name):
                first_loads = getattr(first_airplane, attribute_name)
                second_loads = getattr(second_airplane, attribute_name)
                assert first_loads is not None
                assert second_loads is not None

                # Compare loosely, with an absolute floor scaled to the load vector's
                # magnitude so that near zero components (e.g., the symmetric Airplanes'
                # side forces and roll moments) don't fail the relative comparison.
                npt.assert_allclose(
                    second_loads,
                    first_loads,
                    rtol=0.01,
                    atol=0.01 * float(np.linalg.norm(first_loads)),
                )

    def test_second_airplane_moment_variants_differ_by_re_referencing_term(
        self,
    ) -> None:
        """This method tests that the second Airplane's moment variants differ by
        exactly the predicted re referencing term.

        This is an identity that holds regardless of the Airplanes' separation: the
        moment (in wind axes, relative to the first Airplane's CG) minus the moment (in
        wind axes, relative to the second Airplane's own CG) equals the wind axes
        rotation of the cross product of Cg_GP1_CgP1 with the total force.

        :return: None
        """
        second_airplane = (
            self.steady_ring_vortex_lattice_method_formation_solver.airplanes[1]
        )
        operating_point = (
            self.steady_ring_vortex_lattice_method_formation_solver.operating_point
        )

        assert second_airplane.forces_G is not None
        expectedDifference_W = _transformations.apply_T_to_vectors(
            operating_point.T_pas_GP1_CgP1_to_W_CgP1,
            np.cross(second_airplane.Cg_GP1_CgP1, second_airplane.forces_G),
            is_position=False,
        )

        assert second_airplane.moments_W_CgP1 is not None
        assert second_airplane.moments_W_Cg is not None
        npt.assert_allclose(
            second_airplane.moments_W_CgP1 - second_airplane.moments_W_Cg,
            expectedDifference_W,
            atol=1e-6,
        )

    def test_first_airplane_moment_variants_are_identical(self) -> None:
        """This method tests that the first Airplane's moment variants relative to its
        own CG are identical to those relative to the first Airplane's CG.

        :return: None
        """
        first_airplane = (
            self.steady_ring_vortex_lattice_method_formation_solver.airplanes[0]
        )

        assert first_airplane.moments_W_Cg is not None
        assert first_airplane.moments_W_CgP1 is not None
        npt.assert_array_equal(
            first_airplane.moments_W_Cg, first_airplane.moments_W_CgP1
        )
        assert first_airplane.momentCoefficients_W_Cg is not None
        assert first_airplane.momentCoefficients_W_CgP1 is not None
        npt.assert_array_equal(
            first_airplane.momentCoefficients_W_Cg,
            first_airplane.momentCoefficients_W_CgP1,
        )
