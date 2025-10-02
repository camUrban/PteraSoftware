"""This module contains a class to test aerodynamics functions."""

import unittest
import numpy as np
import numpy.testing as npt

import pterasoftware as ps
from tests.unit.fixtures import aerodynamics_functions_fixtures


class TestAerodynamicsFunctions(unittest.TestCase):
    """This is a class with functions to test aerodynamics functions."""

    def setUp(self):
        """Set up test fixtures for aerodynamics function tests."""
        # Evaluation point fixtures.
        self.single_point = aerodynamics_functions_fixtures.make_single_point_fixture()
        self.grid_of_points = (
            aerodynamics_functions_fixtures.make_grid_of_points_fixture()
        )
        self.line_of_points = (
            aerodynamics_functions_fixtures.make_line_of_points_fixture()
        )
        self.random_points = (
            aerodynamics_functions_fixtures.make_random_points_fixture()
        )

        # Create fixtures for ndarrays of RingVortices.
        (
            self.simple_ring_Brrvp,
            self.simple_ring_Frrvp,
            self.simple_ring_Flrvp,
            self.simple_ring_Blrvp,
            self.simple_ring_strengths,
        ) = aerodynamics_functions_fixtures.make_simple_ring_vortex_arrays_fixture()

        (
            self.multiple_ring_Brrvp,
            self.multiple_ring_Frrvp,
            self.multiple_ring_Flrvp,
            self.multiple_ring_Blrvp,
            self.multiple_ring_strengths,
        ) = aerodynamics_functions_fixtures.make_multiple_ring_vortex_arrays_fixture()

        # Create fixtures for ndarrays of HorseshoeVortices.
        (
            self.simple_horseshoe_Brhvp,
            self.simple_horseshoe_Frhvp,
            self.simple_horseshoe_Flhvp,
            self.simple_horseshoe_Blhvp,
            self.simple_horseshoe_strengths,
        ) = (
            aerodynamics_functions_fixtures.make_simple_horseshoe_vortex_arrays_fixture()
        )

        (
            self.multiple_horseshoe_Brhvp,
            self.multiple_horseshoe_Frhvp,
            self.multiple_horseshoe_Flhvp,
            self.multiple_horseshoe_Blhvp,
            self.multiple_horseshoe_strengths,
        ) = (
            aerodynamics_functions_fixtures.make_multiple_horseshoe_vortex_arrays_fixture()
        )

        # Create age and viscosity fixtures.
        self.ages = aerodynamics_functions_fixtures.make_ages_fixture()
        self.zero_ages = aerodynamics_functions_fixtures.make_zero_ages_fixture()
        self.kinematic_viscosity = (
            aerodynamics_functions_fixtures.make_kinematic_viscosity_fixture()
        )

    def tearDown(self):
        """Clean up test fixtures."""
        del self.single_point
        del self.grid_of_points
        del self.line_of_points
        del self.random_points
        del self.simple_ring_Brrvp
        del self.simple_ring_Frrvp
        del self.simple_ring_Flrvp
        del self.simple_ring_Blrvp
        del self.simple_ring_strengths
        del self.multiple_ring_Brrvp
        del self.multiple_ring_Frrvp
        del self.multiple_ring_Flrvp
        del self.multiple_ring_Blrvp
        del self.multiple_ring_strengths
        del self.simple_horseshoe_Brhvp
        del self.simple_horseshoe_Frhvp
        del self.simple_horseshoe_Flhvp
        del self.simple_horseshoe_Blhvp
        del self.simple_horseshoe_strengths
        del self.multiple_horseshoe_Brhvp
        del self.multiple_horseshoe_Frhvp
        del self.multiple_horseshoe_Flhvp
        del self.multiple_horseshoe_Blhvp
        del self.multiple_horseshoe_strengths
        del self.ages
        del self.zero_ages
        del self.kinematic_viscosity

    def test_collapsed_velocities_from_ring_vortices_single_point(self):
        """Test collapsed_velocities_from_ring_vortices with single evaluation point."""
        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.simple_ring_Brrvp,
            stackFrrvp_G_Cg=self.simple_ring_Frrvp,
            stackFlrvp_G_Cg=self.simple_ring_Flrvp,
            stackBlrvp_G_Cg=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify the output shape.
        self.assertEqual(velocities.shape, (1, 3))

        # Verify the output is not all zeros (unless vortex has zero strength).
        if np.all(self.simple_ring_strengths != 0):
            self.assertFalse(np.allclose(velocities, 0.0))

    def test_collapsed_velocities_from_ring_vortices_multiple_points(self):
        """Test collapsed_velocities_from_ring_vortices with multiple evaluation
        points."""
        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.grid_of_points,
            stackBrrvp_G_Cg=self.simple_ring_Brrvp,
            stackFrrvp_G_Cg=self.simple_ring_Frrvp,
            stackFlrvp_G_Cg=self.simple_ring_Flrvp,
            stackBlrvp_G_Cg=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (25, 3))

    def test_collapsed_velocities_from_ring_vortices_multiple_vortices(self):
        """Test collapsed_velocities_from_ring_vortices with multiple RingVortices."""
        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.multiple_ring_Brrvp,
            stackFrrvp_G_Cg=self.multiple_ring_Frrvp,
            stackFlrvp_G_Cg=self.multiple_ring_Flrvp,
            stackBlrvp_G_Cg=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

    def test_collapsed_velocities_from_ring_vortices_with_ages(self):
        """Test collapsed_velocities_from_ring_vortices with age parameters."""
        # Call the function with ages.
        velocities_with_ages = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.multiple_ring_Brrvp,
            stackFrrvp_G_Cg=self.multiple_ring_Frrvp,
            stackFlrvp_G_Cg=self.multiple_ring_Flrvp,
            stackBlrvp_G_Cg=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            ages=self.ages,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities_with_ages.shape, (1, 3))

    def test_collapsed_velocities_from_ring_vortices_zero_strength(self):
        """Test collapsed_velocities_from_ring_vortices with zero strength
        RingVortices."""
        # Create zero strength array.
        zero_strengths = np.zeros_like(self.simple_ring_strengths)

        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.simple_ring_Brrvp,
            stackFrrvp_G_Cg=self.simple_ring_Frrvp,
            stackFlrvp_G_Cg=self.simple_ring_Flrvp,
            stackBlrvp_G_Cg=self.simple_ring_Blrvp,
            strengths=zero_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output is zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 3), dtype=float))

    def test_collapsed_velocities_from_ring_vortices_chordwise_segments_single_point(
        self,
    ):
        """Test collapsed_velocities_from_ring_vortices_chordwise_segments with single
        evaluation point."""
        # Call the function.
        velocities = (
            ps.aerodynamics.collapsed_velocities_from_ring_vortices_chordwise_segments(
                stackP_G_Cg=self.single_point,
                stackBrrvp_G_Cg=self.simple_ring_Brrvp,
                stackFrrvp_G_Cg=self.simple_ring_Frrvp,
                stackFlrvp_G_Cg=self.simple_ring_Flrvp,
                stackBlrvp_G_Cg=self.simple_ring_Blrvp,
                strengths=self.simple_ring_strengths,
                ages=None,
                nu=self.kinematic_viscosity,
            )
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

    def test_collapsed_velocities_from_ring_vortices_chordwise_segments_multiple_points(
        self,
    ):
        """Test collapsed_velocities_from_ring_vortices_chordwise_segments with multiple
        evaluation points."""
        # Call the function.
        velocities = (
            ps.aerodynamics.collapsed_velocities_from_ring_vortices_chordwise_segments(
                stackP_G_Cg=self.line_of_points,
                stackBrrvp_G_Cg=self.simple_ring_Brrvp,
                stackFrrvp_G_Cg=self.simple_ring_Frrvp,
                stackFlrvp_G_Cg=self.simple_ring_Flrvp,
                stackBlrvp_G_Cg=self.simple_ring_Blrvp,
                strengths=self.simple_ring_strengths,
                ages=None,
                nu=self.kinematic_viscosity,
            )
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (10, 3))

    def test_expanded_velocities_from_ring_vortices_single_point(self):
        """Test expanded_velocities_from_ring_vortices with single evaluation point."""
        # Call the function.
        velocities = ps.aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.simple_ring_Brrvp,
            stackFrrvp_G_Cg=self.simple_ring_Frrvp,
            stackFlrvp_G_Cg=self.simple_ring_Flrvp,
            stackBlrvp_G_Cg=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape (1 evaluation point, 1 vortex).
        self.assertEqual(velocities.shape, (1, 1, 3))

    def test_expanded_velocities_from_ring_vortices_multiple_vortices(self):
        """Test expanded_velocities_from_ring_vortices with multiple RingVortices."""
        # Call the function.
        velocities = ps.aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.multiple_ring_Brrvp,
            stackFrrvp_G_Cg=self.multiple_ring_Frrvp,
            stackFlrvp_G_Cg=self.multiple_ring_Flrvp,
            stackBlrvp_G_Cg=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape (1 evaluation point, 3 vortices).
        self.assertEqual(velocities.shape, (1, 3, 3))

    def test_expanded_velocities_from_ring_vortices_multiple_points(self):
        """Test expanded_velocities_from_ring_vortices with multiple evaluation
        points."""
        # Call the function.
        velocities = ps.aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.line_of_points,
            stackBrrvp_G_Cg=self.simple_ring_Brrvp,
            stackFrrvp_G_Cg=self.simple_ring_Frrvp,
            stackFlrvp_G_Cg=self.simple_ring_Flrvp,
            stackBlrvp_G_Cg=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape (10 evaluation points, 1 vortex).
        self.assertEqual(velocities.shape, (10, 1, 3))

    def test_expanded_and_collapsed_ring_vortices_consistency(self):
        """Test that expanded and collapsed RingVortex functions are consistent."""
        # Get expanded velocities.
        expanded_velocities = ps.aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.multiple_ring_Brrvp,
            stackFrrvp_G_Cg=self.multiple_ring_Frrvp,
            stackFlrvp_G_Cg=self.multiple_ring_Flrvp,
            stackBlrvp_G_Cg=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Get collapsed velocities.
        collapsed_velocities = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.single_point,
            stackBrrvp_G_Cg=self.multiple_ring_Brrvp,
            stackFrrvp_G_Cg=self.multiple_ring_Frrvp,
            stackFlrvp_G_Cg=self.multiple_ring_Flrvp,
            stackBlrvp_G_Cg=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Sum expanded velocities along vortex axis.
        summed_expanded = np.sum(expanded_velocities, axis=1)

        # Verify that summed expanded equals collapsed.
        npt.assert_array_almost_equal(summed_expanded, collapsed_velocities, decimal=10)

    def test_collapsed_velocities_from_horseshoe_vortices_single_point(self):
        """Test collapsed_velocities_from_horseshoe_vortices with single evaluation
        point."""
        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.single_point,
            stackBrhvp_G_Cg=self.simple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.simple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.simple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.simple_horseshoe_Blhvp,
            strengths=self.simple_horseshoe_strengths,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

        # Verify output is not all zeros (unless vortex has zero strength).
        if np.all(self.simple_horseshoe_strengths != 0):
            self.assertFalse(np.allclose(velocities, 0.0))

    def test_collapsed_velocities_from_horseshoe_vortices_multiple_points(self):
        """Test collapsed_velocities_from_horseshoe_vortices with multiple evaluation
        points."""
        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.grid_of_points,
            stackBrhvp_G_Cg=self.simple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.simple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.simple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.simple_horseshoe_Blhvp,
            strengths=self.simple_horseshoe_strengths,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (25, 3))

    def test_collapsed_velocities_from_horseshoe_vortices_multiple_vortices(self):
        """Test collapsed_velocities_from_horseshoe_vortices with multiple
        HorseshoeVortices."""
        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.single_point,
            stackBrhvp_G_Cg=self.multiple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.multiple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.multiple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.multiple_horseshoe_Blhvp,
            strengths=self.multiple_horseshoe_strengths,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

    def test_collapsed_velocities_from_horseshoe_vortices_zero_strength(self):
        """Test collapsed_velocities_from_horseshoe_vortices with zero strength
        HorseshoeVortices."""
        # Create zero strength array.
        zero_strengths = np.zeros_like(self.simple_horseshoe_strengths)

        # Call the function.
        velocities = ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.single_point,
            stackBrhvp_G_Cg=self.simple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.simple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.simple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.simple_horseshoe_Blhvp,
            strengths=zero_strengths,
        )

        # Verify output is zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 3), dtype=float))

    def test_expanded_velocities_from_horseshoe_vortices_single_point(self):
        """Test expanded_velocities_from_horseshoe_vortices with single evaluation
        point."""
        # Call the function.
        velocities = ps.aerodynamics.expanded_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.single_point,
            stackBrhvp_G_Cg=self.simple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.simple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.simple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.simple_horseshoe_Blhvp,
            strengths=self.simple_horseshoe_strengths,
        )

        # Verify output shape (1 evaluation point, 1 vortex).
        self.assertEqual(velocities.shape, (1, 1, 3))

    def test_expanded_velocities_from_horseshoe_vortices_multiple_vortices(self):
        """Test expanded_velocities_from_horseshoe_vortices with multiple
        HorseshoeVortices."""
        # Call the function.
        velocities = ps.aerodynamics.expanded_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.single_point,
            stackBrhvp_G_Cg=self.multiple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.multiple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.multiple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.multiple_horseshoe_Blhvp,
            strengths=self.multiple_horseshoe_strengths,
        )

        # Verify output shape (1 evaluation point, 2 vortices).
        self.assertEqual(velocities.shape, (1, 2, 3))

    def test_expanded_velocities_from_horseshoe_vortices_multiple_points(self):
        """Test expanded_velocities_from_horseshoe_vortices with multiple evaluation
        points."""
        # Call the function.
        velocities = ps.aerodynamics.expanded_velocities_from_horseshoe_vortices(
            stackP_G_Cg=self.line_of_points,
            stackBrhvp_G_Cg=self.simple_horseshoe_Brhvp,
            stackFrhvp_G_Cg=self.simple_horseshoe_Frhvp,
            stackFlhvp_G_Cg=self.simple_horseshoe_Flhvp,
            stackBlhvp_G_Cg=self.simple_horseshoe_Blhvp,
            strengths=self.simple_horseshoe_strengths,
        )

        # Verify the output shape (10 evaluation points, 1 vortex).
        self.assertEqual(velocities.shape, (10, 1, 3))

    def test_expanded_and_collapsed_horseshoe_vortices_consistency(self):
        """Test that expanded and collapsed HorseshoeVortex functions are
        consistent."""
        # Get the expanded velocities.
        expanded_velocities = (
            ps.aerodynamics.expanded_velocities_from_horseshoe_vortices(
                stackP_G_Cg=self.single_point,
                stackBrhvp_G_Cg=self.multiple_horseshoe_Brhvp,
                stackFrhvp_G_Cg=self.multiple_horseshoe_Frhvp,
                stackFlhvp_G_Cg=self.multiple_horseshoe_Flhvp,
                stackBlhvp_G_Cg=self.multiple_horseshoe_Blhvp,
                strengths=self.multiple_horseshoe_strengths,
            )
        )

        # Get collapsed velocities.
        collapsed_velocities = (
            ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
                stackP_G_Cg=self.single_point,
                stackBrhvp_G_Cg=self.multiple_horseshoe_Brhvp,
                stackFrhvp_G_Cg=self.multiple_horseshoe_Frhvp,
                stackFlhvp_G_Cg=self.multiple_horseshoe_Flhvp,
                stackBlhvp_G_Cg=self.multiple_horseshoe_Blhvp,
                strengths=self.multiple_horseshoe_strengths,
            )
        )

        # Sum expanded velocities along vortex axis.
        summed_expanded = np.sum(expanded_velocities, axis=1)

        # Verify that summed expanded equals collapsed.
        npt.assert_array_almost_equal(summed_expanded, collapsed_velocities, decimal=10)

    def test_velocity_functions_with_random_points(self):
        """Test velocity functions with random evaluation points."""
        # Test RingVortex function.
        ring_velocities = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=self.random_points,
            stackBrrvp_G_Cg=self.simple_ring_Brrvp,
            stackFrrvp_G_Cg=self.simple_ring_Frrvp,
            stackFlrvp_G_Cg=self.simple_ring_Flrvp,
            stackBlrvp_G_Cg=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(ring_velocities.shape, (20, 3))

        # Test HorseshoeVortex function.
        horseshoe_velocities = (
            ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
                stackP_G_Cg=self.random_points,
                stackBrhvp_G_Cg=self.simple_horseshoe_Brhvp,
                stackFrhvp_G_Cg=self.simple_horseshoe_Frhvp,
                stackFlhvp_G_Cg=self.simple_horseshoe_Flhvp,
                stackBlhvp_G_Cg=self.simple_horseshoe_Blhvp,
                strengths=self.simple_horseshoe_strengths,
            )
        )

        # Verify output shape.
        self.assertEqual(horseshoe_velocities.shape, (20, 3))

    @staticmethod
    def ref_calculate_biot_savart_velocity(S_A_a, E_A_a, P_A_a, gamma):
        """Calculate induced velocity using Biot-Savart formula for a conceptual line
        vortex.

        This is a reference implementation to validate the aerodynamics functions.

        Formula:
            v_A__I = (gamma/(4*pi)) * (r3_A / |r3_A|^2) * [r0_A · (r1Hat_A - r2Hat_A)]

        Where:
            r1_A = P_A_a - S_A_a
            r2_A = P_A_a - E_A_a
            r0_A = E_A_a - S_A_a
            r3_A = r1_A x r2_A
            r1Hat_A = r1_A/|r1_A|
            r2Hat_A = r2_A/|r2_A|

        :param S_A_a: (3,) ndarray of floats
            Start point of the conceptual line vortex (in A axes, relative to point
            a) in meters.
        :param E_A_a: (3,) ndarray of floats
            End point of the conceptual line vortex (in A axes, relative to point a)
            in meters.
        :param P_A_a: (3,) ndarray of floats
            Evaluation point (in A axes, relative to point a) in meters.
        :param gamma: float
            Conceptual line vortex strength in meters squared per second.
        :return v_A__I: (3,) ndarray of floats
            Induced velocity (in A axes, observed from an inertial frame) in meters per
            second.
        """
        # Get machine epsilon for singularity checks.
        eps = np.finfo(float).eps

        # Calculate vectors.
        r1_A = P_A_a - S_A_a
        r2_A = P_A_a - E_A_a
        r0_A = E_A_a - S_A_a

        # Calculate cross product vector.
        r3_A = np.cross(r1_A, r2_A)

        # Find vector lengths
        r1_norm = np.linalg.norm(r1_A)
        r2_norm = np.linalg.norm(r2_A)
        r3_norm = np.linalg.norm(r3_A)

        # Find cross product vector's length squared
        r3_norm_sq = r3_norm**2

        # Check for singularities (point on or very close to vortex).
        if r1_norm < eps or r2_norm < eps or r3_norm_sq < eps:
            return np.zeros(3, dtype=float)

        # Calculate unit vectors.
        r1Hat_A = r1_A / r1_norm
        r2Hat_A = r2_A / r2_norm

        # Biot-Savart formula.
        v_A__I = (
            (gamma / (4.0 * np.pi))
            * (r3_A / r3_norm_sq)
            * np.dot(r0_A, r1Hat_A - r2Hat_A)
        )

        return v_A__I

    def test_ref_biot_savart_single_line_vortex_along_axis(self):
        """Test reference Biot-Savart implementation for a conceptual line vortex
        along z-axis with a perpendicular evaluation point."""
        A_A_a = np.array([0.0, 0.0, 0.0], dtype=float)
        B_A_a = np.array([0.0, 0.0, 1.0], dtype=float)
        gamma = 2.0

        P_A_a = np.array([1.0, 0.0, 0.5], dtype=float)

        # Calculate the expected velocity using the reference implementation.
        v_A__I = self.ref_calculate_biot_savart_velocity(A_A_a, B_A_a, P_A_a, gamma)

        # Verify the velocity is in y-direction (due to right-hand rule). For a
        # conceptual line vortex along +z and evaluation point at +x, the velocity
        # should be in +y.
        self.assertGreater(v_A__I[1], 0.0)
        self.assertAlmostEqual(v_A__I[0], 0.0, places=10)
        self.assertAlmostEqual(v_A__I[2], 0.0, places=10)

    def test_ref_biot_savart_symmetric_configuration(self):
        """Test reference Biot-Savart implementation for symmetric vortex
        configuration."""
        # Two line vortices forming a V-shape, symmetric about xz-plane (in A axes).
        # Vortex 1: from (1, 1, 0) to origin (in A axes, relative to the a point).
        # Vortex 2: from origin to (1, -1, 0) (in A axes, relative to the a point).
        S1_A_a = np.array([1.0, 1.0, 0.0], dtype=float)
        E1_A_a = np.array([0.0, 0.0, 0.0], dtype=float)
        S2_A_a = np.array([0.0, 0.0, 0.0], dtype=float)
        E2_A_a = np.array([1.0, -1.0, 0.0], dtype=float)

        # Evaluation point on xz-plane at (0.5, 0, 0) (in A axes, relative to the a
        # point).
        P_A_a = np.array([0.5, 0.0, 0.0], dtype=float)

        # Same strength for both vortices.
        gamma = 1.0

        # Calculate velocities (in A axes, observed from an inertial frame).
        v1_A__I = self.ref_calculate_biot_savart_velocity(S1_A_a, E1_A_a, P_A_a, gamma)
        v2_A__I = self.ref_calculate_biot_savart_velocity(S2_A_a, E2_A_a, P_A_a, gamma)

        # The total velocity (in A axes, observed from an inertial frame) should have
        # a zero y-component (due to symmetry) and should have a non-zero magnitude.
        v_A__I = v1_A__I + v2_A__I
        self.assertAlmostEqual(v_A__I[1], 0.0, places=10)
        self.assertGreater(np.linalg.norm(v_A__I), 0.0)

    def test_ref_biot_savart_point_on_vortex_returns_zero(self):
        """Test that the reference Biot-Savart implementation returns zero velocity
        for an evaluation point on a conceptual line vortex."""
        # Line vortex from origin to (1, 0, 0).
        A = np.array([0.0, 0.0, 0.0], dtype=float)
        B = np.array([1.0, 0.0, 0.0], dtype=float)

        # Evaluation point on the vortex.
        P = np.array([0.5, 0.0, 0.0], dtype=float)

        # Vortex strength.
        gamma = 1.0

        # Calculate velocity - should be zero due to singularity.
        velocity = self.ref_calculate_biot_savart_velocity(A, B, P, gamma)

        # Verify velocity is zero.
        npt.assert_array_almost_equal(velocity, np.zeros(3, dtype=float), decimal=10)

    def test_ref_biot_savart_point_at_vortex_endpoint_returns_zero(self):
        """Test that, using the reference Biot-Savart implementation, an evaluation
        point at a conceptual line vortex endpoint returns zero velocity."""
        # Line vortex from origin to (1, 0, 0).
        A = np.array([0.0, 0.0, 0.0], dtype=float)
        B = np.array([1.0, 0.0, 0.0], dtype=float)

        # Evaluation point at start point.
        P_start = A.copy()

        # Evaluation point at end point.
        P_end = B.copy()

        # Vortex strength.
        gamma = 1.0

        # Calculate velocities - both should be zero.
        velocity_start = self.ref_calculate_biot_savart_velocity(A, B, P_start, gamma)
        velocity_end = self.ref_calculate_biot_savart_velocity(A, B, P_end, gamma)

        # Verify velocities are zero.
        npt.assert_array_almost_equal(
            velocity_start, np.zeros(3, dtype=float), decimal=10
        )
        npt.assert_array_almost_equal(
            velocity_end, np.zeros(3, dtype=float), decimal=10
        )

    def test_ref_biot_savart_far_field_falloff(self):
        """Test that, using the reference Biot-Savart implementation, induced
        velocity decreases with distance from a conceptual line vortex."""
        # Line vortex from (0, -0.5, 0) to (0, 0.5, 0).
        A = np.array([0.0, -0.5, 0.0], dtype=float)
        B = np.array([0.0, 0.5, 0.0], dtype=float)

        # Vortex strength.
        gamma = 1.0

        # Evaluation points at increasing distances along x-axis.
        distances = [1.0, 2.0, 5.0, 10.0]
        velocities = []

        for dist in distances:
            P = np.array([dist, 0.0, 0.0], dtype=float)
            v = self.ref_calculate_biot_savart_velocity(A, B, P, gamma)
            velocities.append(np.linalg.norm(v))

        # Verify that velocity magnitude decreases with distance.
        for i in range(len(velocities) - 1):
            self.assertGreater(velocities[i], velocities[i + 1])

    def test_ref_biot_savart_opposite_strength_cancellation(self):
        """Test that, using the reference Biot-Savart implementation, two conceptual
        line vortices with opposite strengths produce opposite velocities."""
        # Line vortex from origin to (1, 0, 0).
        A = np.array([0.0, 0.0, 0.0], dtype=float)
        B = np.array([1.0, 0.0, 0.0], dtype=float)

        # Evaluation point above vortex.
        P = np.array([0.5, 0.0, 1.0], dtype=float)

        # Calculate velocity with positive strength.
        v_positive = self.ref_calculate_biot_savart_velocity(A, B, P, 1.0)

        # Calculate velocity with negative strength.
        v_negative = self.ref_calculate_biot_savart_velocity(A, B, P, -1.0)

        # Verify velocities are opposite.
        npt.assert_array_almost_equal(v_positive, -v_negative, decimal=10)

    def test_ring_vortex_decomposition_against_ref_biot_savart(self):
        """Test that a RingVortex induces the same velocity as four conceptual line
        vortices."""
        # Create a unit square in the xy-plane (in geometry axes, relative to the CG).
        Br_G_Cg = np.array([1.0, 1.0, 0.0], dtype=float)
        Fr_G_Cg = np.array([0.0, 1.0, 0.0], dtype=float)
        Fl_G_Cg = np.array([0.0, 0.0, 0.0], dtype=float)
        Bl_G_Cg = np.array([1.0, 0.0, 0.0], dtype=float)

        # Create an evaluation point above the square's center (in geometry axes,
        # relative to the CG).
        P_G_Cg = np.array([0.5, 0.5, 1.0], dtype=float)

        # Calculate velocity induced (in geometry axes, observed from the Earth frame)
        # by each conceptual line vortex using the reference implementation.
        vRight_G__E = self.ref_calculate_biot_savart_velocity(
            Br_G_Cg, Fr_G_Cg, P_G_Cg, 1.0
        )
        vFront_G__E = self.ref_calculate_biot_savart_velocity(
            Fr_G_Cg, Fl_G_Cg, P_G_Cg, 1.0
        )
        vLeft_G__E = self.ref_calculate_biot_savart_velocity(
            Fl_G_Cg, Bl_G_Cg, P_G_Cg, 1.0
        )
        vBack_G__E = self.ref_calculate_biot_savart_velocity(
            Bl_G_Cg, Br_G_Cg, P_G_Cg, 1.0
        )

        # Sum to get total velocity induced by all the conceptual line vortices (in
        # geometry axes, observed from the Earth frame).
        vLines_G__E = vRight_G__E + vFront_G__E + vLeft_G__E + vBack_G__E

        stackP_G_Cg = P_G_Cg.reshape(1, 3)
        stackBrrvp_G_Cg = Br_G_Cg.reshape(1, 3)
        stackFrrvp_G_Cg = Fr_G_Cg.reshape(1, 3)
        stackFlrvp_G_Cg = Fl_G_Cg.reshape(1, 3)
        stackBlrvp_G_Cg = Bl_G_Cg.reshape(1, 3)
        strengths = np.array([1.0], dtype=float)

        # Calculate velocity induced by an equivalent RingVortex (in geometry axes,
        # observed from the Earth frame) using the
        # collapsed_velocities_from_ring_vortices function.
        vRing_G__E = ps.aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackBrrvp_G_Cg=stackBrrvp_G_Cg,
            stackFrrvp_G_Cg=stackFrrvp_G_Cg,
            stackFlrvp_G_Cg=stackFlrvp_G_Cg,
            stackBlrvp_G_Cg=stackBlrvp_G_Cg,
            strengths=strengths,
            ages=None,
            nu=0.0,
        )[0]

        npt.assert_array_almost_equal(vRing_G__E, vLines_G__E, decimal=10)

    def test_horseshoe_vortex_decomposition_against_ref_biot_savart(self):
        """Test that a HorseshoeVortex induces the same velocity as three conceptual
        line vortices."""
        # Create horseshoe vortex with finite leg along y-axis and infinite legs along
        # x-axis.
        Fr = np.array([0.0, 0.5, 0.0], dtype=float)
        Fl = np.array([0.0, -0.5, 0.0], dtype=float)
        Br = np.array([20.0, 0.5, 0.0], dtype=float)
        Bl = np.array([20.0, -0.5, 0.0], dtype=float)

        # Evaluation point in front of the horseshoe.
        P = np.array([-1.0, 0.0, 0.0], dtype=float)

        # Vortex strength.
        gamma = 1.0

        # Calculate velocity induced by each line segment.
        v_right = self.ref_calculate_biot_savart_velocity(Br, Fr, P, gamma)
        v_finite = self.ref_calculate_biot_savart_velocity(Fr, Fl, P, gamma)
        v_left = self.ref_calculate_biot_savart_velocity(Fl, Bl, P, gamma)

        # Sum to get total velocity.
        expected_total = v_right + v_finite + v_left

        # Calculate velocity using HorseshoeVortex function.
        stackP_G_Cg = P.reshape(1, 3)
        stackBrhvp_G_Cg = Br.reshape(1, 3)
        stackFrhvp_G_Cg = Fr.reshape(1, 3)
        stackFlhvp_G_Cg = Fl.reshape(1, 3)
        stackBlhvp_G_Cg = Bl.reshape(1, 3)
        strengths = np.array([gamma], dtype=float)

        computed_total = ps.aerodynamics.collapsed_velocities_from_horseshoe_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackBrhvp_G_Cg=stackBrhvp_G_Cg,
            stackFrhvp_G_Cg=stackFrhvp_G_Cg,
            stackFlhvp_G_Cg=stackFlhvp_G_Cg,
            stackBlhvp_G_Cg=stackBlhvp_G_Cg,
            strengths=strengths,
        )[0]

        # Verify that computed velocity matches expected velocity.
        npt.assert_array_almost_equal(computed_total, expected_total, decimal=10)


if __name__ == "__main__":
    unittest.main()
