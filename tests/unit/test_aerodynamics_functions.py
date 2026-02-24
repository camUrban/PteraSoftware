"""This module contains a class to test aerodynamics functions."""

import math
import unittest

import numpy as np
import numpy.testing as npt

# noinspection PyProtectedMember
from pterasoftware import _aerodynamics_functions
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

        # Create initial core radius fixtures.
        self.simple_ring_rc0s = aerodynamics_functions_fixtures.make_rc0s_fixture(1)
        self.multiple_ring_rc0s = aerodynamics_functions_fixtures.make_rc0s_fixture(3)
        self.simple_horseshoe_rc0s = aerodynamics_functions_fixtures.make_rc0s_fixture(
            1
        )
        self.multiple_horseshoe_rc0s = (
            aerodynamics_functions_fixtures.make_rc0s_fixture(2)
        )

        # Create age and viscosity fixtures.
        self.ages = aerodynamics_functions_fixtures.make_ages_fixture()
        self.zero_ages = aerodynamics_functions_fixtures.make_zero_ages_fixture()
        self.kinematic_viscosity = (
            aerodynamics_functions_fixtures.make_kinematic_viscosity_fixture()
        )

    def test_collapsed_velocities_from_ring_vortices_single_point(self):
        """Test collapsed_velocities_from_ring_vortices with single evaluation point."""
        # Call the function.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.single_point,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
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
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.grid_of_points,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (25, 3))

    def test_collapsed_velocities_from_ring_vortices_multiple_vortices(self):
        """Test collapsed_velocities_from_ring_vortices with multiple RingVortices."""
        # Call the function.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.single_point,
            stackBrrvp_GP1_CgP1=self.multiple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.multiple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.multiple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            r_c0s=self.multiple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

    def test_collapsed_velocities_from_ring_vortices_with_ages(self):
        """Test collapsed_velocities_from_ring_vortices with age parameters."""
        # Call the function with ages.
        velocities_with_ages = (
            _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrrvp_GP1_CgP1=self.multiple_ring_Brrvp,
                stackFrrvp_GP1_CgP1=self.multiple_ring_Frrvp,
                stackFlrvp_GP1_CgP1=self.multiple_ring_Flrvp,
                stackBlrvp_GP1_CgP1=self.multiple_ring_Blrvp,
                strengths=self.multiple_ring_strengths,
                r_c0s=self.multiple_ring_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
                ages=self.ages,
                nu=self.kinematic_viscosity,
            )
        )

        # Verify output shape.
        self.assertEqual(velocities_with_ages.shape, (1, 3))

    def test_collapsed_velocities_from_ring_vortices_zero_strength(self):
        """Test collapsed_velocities_from_ring_vortices with zero strength
        RingVortices."""
        # Create zero strength array.
        zero_strengths = np.zeros_like(self.simple_ring_strengths)

        # Call the function.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.single_point,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=zero_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
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
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices_chordwise_segments(
            stackP_GP1_CgP1=self.single_point,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

    def test_collapsed_velocities_from_ring_vortices_chordwise_segments_multiple_points(
        self,
    ):
        """Test collapsed_velocities_from_ring_vortices_chordwise_segments with multiple
        evaluation points."""
        # Call the function.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices_chordwise_segments(
            stackP_GP1_CgP1=self.line_of_points,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (10, 3))

    def test_expanded_velocities_from_ring_vortices_single_point(self):
        """Test expanded_velocities_from_ring_vortices with single evaluation point."""
        # Call the function.
        velocities = _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.single_point,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape (1 evaluation point, 1 vortex).
        self.assertEqual(velocities.shape, (1, 1, 3))

    def test_expanded_velocities_from_ring_vortices_multiple_vortices(self):
        """Test expanded_velocities_from_ring_vortices with multiple RingVortices."""
        # Call the function.
        velocities = _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.single_point,
            stackBrrvp_GP1_CgP1=self.multiple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.multiple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.multiple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.multiple_ring_Blrvp,
            strengths=self.multiple_ring_strengths,
            r_c0s=self.multiple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape (1 evaluation point, 3 vortices).
        self.assertEqual(velocities.shape, (1, 3, 3))

    def test_expanded_velocities_from_ring_vortices_multiple_points(self):
        """Test expanded_velocities_from_ring_vortices with multiple evaluation
        points."""
        # Call the function.
        velocities = _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.line_of_points,
            stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
            strengths=self.simple_ring_strengths,
            r_c0s=self.simple_ring_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=None,
            nu=self.kinematic_viscosity,
        )

        # Verify output shape (10 evaluation points, 1 vortex).
        self.assertEqual(velocities.shape, (10, 1, 3))

    def test_expanded_and_collapsed_ring_vortices_consistency(self):
        """Test that expanded and collapsed RingVortex functions are consistent."""
        # Get expanded velocities.
        expanded_velocities = (
            _aerodynamics_functions.expanded_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrrvp_GP1_CgP1=self.multiple_ring_Brrvp,
                stackFrrvp_GP1_CgP1=self.multiple_ring_Frrvp,
                stackFlrvp_GP1_CgP1=self.multiple_ring_Flrvp,
                stackBlrvp_GP1_CgP1=self.multiple_ring_Blrvp,
                strengths=self.multiple_ring_strengths,
                r_c0s=self.multiple_ring_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
                ages=None,
                nu=self.kinematic_viscosity,
            )
        )

        # Get collapsed velocities.
        collapsed_velocities = (
            _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrrvp_GP1_CgP1=self.multiple_ring_Brrvp,
                stackFrrvp_GP1_CgP1=self.multiple_ring_Frrvp,
                stackFlrvp_GP1_CgP1=self.multiple_ring_Flrvp,
                stackBlrvp_GP1_CgP1=self.multiple_ring_Blrvp,
                strengths=self.multiple_ring_strengths,
                r_c0s=self.multiple_ring_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
                ages=None,
                nu=self.kinematic_viscosity,
            )
        )

        # Sum expanded velocities along vortex axis.
        summed_expanded = np.sum(expanded_velocities, axis=1)

        # Verify that summed expanded equals collapsed.
        npt.assert_array_almost_equal(summed_expanded, collapsed_velocities, decimal=10)

    def test_collapsed_velocities_from_horseshoe_vortices_single_point(self):
        """Test collapsed_velocities_from_horseshoe_vortices with single evaluation
        point."""
        # Call the function.
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.simple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.simple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.simple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.simple_horseshoe_Blhvp,
                strengths=self.simple_horseshoe_strengths,
                r_c0s=self.simple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
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
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.grid_of_points,
                stackBrhvp_GP1_CgP1=self.simple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.simple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.simple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.simple_horseshoe_Blhvp,
                strengths=self.simple_horseshoe_strengths,
                r_c0s=self.simple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (25, 3))

    def test_collapsed_velocities_from_horseshoe_vortices_multiple_vortices(self):
        """Test collapsed_velocities_from_horseshoe_vortices with multiple
        HorseshoeVortices."""
        # Call the function.
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.multiple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.multiple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.multiple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.multiple_horseshoe_Blhvp,
                strengths=self.multiple_horseshoe_strengths,
                r_c0s=self.multiple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify output shape.
        self.assertEqual(velocities.shape, (1, 3))

    def test_collapsed_velocities_from_horseshoe_vortices_zero_strength(self):
        """Test collapsed_velocities_from_horseshoe_vortices with zero strength
        HorseshoeVortices."""
        # Create zero strength array.
        zero_strengths = np.zeros_like(self.simple_horseshoe_strengths)

        # Call the function.
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.simple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.simple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.simple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.simple_horseshoe_Blhvp,
                strengths=zero_strengths,
                r_c0s=self.simple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify output is zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 3), dtype=float))

    def test_expanded_velocities_from_horseshoe_vortices_single_point(self):
        """Test expanded_velocities_from_horseshoe_vortices with single evaluation
        point."""
        # Call the function.
        velocities = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.simple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.simple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.simple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.simple_horseshoe_Blhvp,
                strengths=self.simple_horseshoe_strengths,
                r_c0s=self.simple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify output shape (1 evaluation point, 1 vortex).
        self.assertEqual(velocities.shape, (1, 1, 3))

    def test_expanded_velocities_from_horseshoe_vortices_multiple_vortices(self):
        """Test expanded_velocities_from_horseshoe_vortices with multiple
        HorseshoeVortices."""
        # Call the function.
        velocities = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.multiple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.multiple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.multiple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.multiple_horseshoe_Blhvp,
                strengths=self.multiple_horseshoe_strengths,
                r_c0s=self.multiple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify output shape (1 evaluation point, 2 vortices).
        self.assertEqual(velocities.shape, (1, 2, 3))

    def test_expanded_velocities_from_horseshoe_vortices_multiple_points(self):
        """Test expanded_velocities_from_horseshoe_vortices with multiple evaluation
        points."""
        # Call the function.
        velocities = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.line_of_points,
                stackBrhvp_GP1_CgP1=self.simple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.simple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.simple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.simple_horseshoe_Blhvp,
                strengths=self.simple_horseshoe_strengths,
                r_c0s=self.simple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify the output shape (10 evaluation points, 1 vortex).
        self.assertEqual(velocities.shape, (10, 1, 3))

    def test_expanded_and_collapsed_horseshoe_vortices_consistency(self):
        """Test that expanded and collapsed HorseshoeVortex functions are
        consistent."""
        # Get the expanded velocities.
        expanded_velocities = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.multiple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.multiple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.multiple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.multiple_horseshoe_Blhvp,
                strengths=self.multiple_horseshoe_strengths,
                r_c0s=self.multiple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Get collapsed velocities.
        collapsed_velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.single_point,
                stackBrhvp_GP1_CgP1=self.multiple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.multiple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.multiple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.multiple_horseshoe_Blhvp,
                strengths=self.multiple_horseshoe_strengths,
                r_c0s=self.multiple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Sum expanded velocities along vortex axis.
        summed_expanded = np.sum(expanded_velocities, axis=1)

        # Verify that summed expanded equals collapsed.
        npt.assert_array_almost_equal(summed_expanded, collapsed_velocities, decimal=10)

    def test_velocity_functions_with_random_points(self):
        """Test velocity functions with random evaluation points."""
        # Test RingVortex function.
        ring_velocities = (
            _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.random_points,
                stackBrrvp_GP1_CgP1=self.simple_ring_Brrvp,
                stackFrrvp_GP1_CgP1=self.simple_ring_Frrvp,
                stackFlrvp_GP1_CgP1=self.simple_ring_Flrvp,
                stackBlrvp_GP1_CgP1=self.simple_ring_Blrvp,
                strengths=self.simple_ring_strengths,
                r_c0s=self.simple_ring_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
                ages=None,
                nu=self.kinematic_viscosity,
            )
        )

        # Verify output shape.
        self.assertEqual(ring_velocities.shape, (20, 3))

        # Test HorseshoeVortex function.
        horseshoe_velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=self.random_points,
                stackBrhvp_GP1_CgP1=self.simple_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.simple_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.simple_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.simple_horseshoe_Blhvp,
                strengths=self.simple_horseshoe_strengths,
                r_c0s=self.simple_horseshoe_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify output shape.
        self.assertEqual(horseshoe_velocities.shape, (20, 3))

    @staticmethod
    def ref_calculate_biot_savart_velocity(S_A_a, E_A_a, P_A_a, gamma):
        """Calculate induced velocity using Biot-Savart formula for a conceptual,
        coreless line vortex.

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
        # Get machine epsilon for absolute singularity check and set tolerance for
        # relative checks.
        eps = np.finfo(float).eps
        tol = 1.0e-10

        # Calculate vectors.
        r1_A = P_A_a - S_A_a
        r2_A = P_A_a - E_A_a
        r0_A = E_A_a - S_A_a

        # Calculate cross product vector.
        r3_A = np.cross(r1_A, r2_A)

        # Find vector lengths
        r0_norm = np.linalg.norm(r0_A)
        r1_norm = np.linalg.norm(r1_A)
        r2_norm = np.linalg.norm(r2_A)
        r3_norm = np.linalg.norm(r3_A)

        # Find cross product vector's length squared
        r3_norm_sq = r3_norm**2

        # Perform absolute check against a degenerate line vortex (near zero length).
        if r0_norm < eps:
            return np.zeros(3, dtype=float)

        # Check for relative singularities.
        if (
            r1_norm / r0_norm < tol
            or r2_norm / r0_norm < tol
            or r3_norm / (r1_norm * r2_norm) < tol
        ):
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
        self.assertGreater(float(v_A__I[1]), 0.0)
        self.assertAlmostEqual(float(v_A__I[0]), 0.0, places=10)
        self.assertAlmostEqual(float(v_A__I[2]), 0.0, places=10)

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

        stackP_GP1_CgP1 = P_G_Cg.reshape(1, 3)
        stackBrrvp_GP1_CgP1 = Br_G_Cg.reshape(1, 3)
        stackFrrvp_GP1_CgP1 = Fr_G_Cg.reshape(1, 3)
        stackFlrvp_GP1_CgP1 = Fl_G_Cg.reshape(1, 3)
        stackBlrvp_GP1_CgP1 = Bl_G_Cg.reshape(1, 3)
        strengths = np.array([1.0], dtype=float)

        # Use 0 for the initial core radius as the reference implementation is a
        # coreless model.
        r_c0s = np.zeros(1, dtype=float)

        # Calculate velocity induced by an equivalent RingVortex (in geometry axes,
        # observed from the Earth frame) using the
        # collapsed_velocities_from_ring_vortices function.
        vRing_G__E = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackBrrvp_GP1_CgP1=stackBrrvp_GP1_CgP1,
            stackFrrvp_GP1_CgP1=stackFrrvp_GP1_CgP1,
            stackFlrvp_GP1_CgP1=stackFlrvp_GP1_CgP1,
            stackBlrvp_GP1_CgP1=stackBlrvp_GP1_CgP1,
            strengths=strengths,
            r_c0s=r_c0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
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
        stackP_GP1_CgP1 = P.reshape(1, 3)
        stackBrhvp_GP1_CgP1 = Br.reshape(1, 3)
        stackFrhvp_GP1_CgP1 = Fr.reshape(1, 3)
        stackFlhvp_GP1_CgP1 = Fl.reshape(1, 3)
        stackBlhvp_GP1_CgP1 = Bl.reshape(1, 3)
        strengths = np.array([gamma], dtype=float)

        # Use 0 for the initial core radius as the reference implementation is a
        # coreless model.
        r_c0s = np.zeros(1, dtype=float)

        computed_total = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrhvp_GP1_CgP1=stackBrhvp_GP1_CgP1,
                stackFrhvp_GP1_CgP1=stackFrhvp_GP1_CgP1,
                stackFlhvp_GP1_CgP1=stackFlhvp_GP1_CgP1,
                stackBlhvp_GP1_CgP1=stackBlhvp_GP1_CgP1,
                strengths=strengths,
                r_c0s=r_c0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )[0]
        )

        # Verify that computed velocity matches expected velocity.
        npt.assert_array_almost_equal(computed_total, expected_total, decimal=10)


class TestSingularityGuards(unittest.TestCase):
    """This is a class with functions to test the scale invariant singularity
    guards in the Biot-Savart kernels."""

    def setUp(self):
        """Set up fixtures for singularity guard tests."""
        # Degenerate RingVortex fixture (all corners at the origin).
        (
            self.degenerate_ring_Brrvp,
            self.degenerate_ring_Frrvp,
            self.degenerate_ring_Flrvp,
            self.degenerate_ring_Blrvp,
            self.degenerate_ring_strengths,
        ) = aerodynamics_functions_fixtures.make_degenerate_ring_vortex_arrays_fixture()

        # Degenerate HorseshoeVortex fixture (all points at the origin).
        (
            self.degenerate_horseshoe_Brhvp,
            self.degenerate_horseshoe_Frhvp,
            self.degenerate_horseshoe_Flhvp,
            self.degenerate_horseshoe_Blhvp,
            self.degenerate_horseshoe_strengths,
        ) = (
            aerodynamics_functions_fixtures.make_degenerate_horseshoe_vortex_arrays_fixture()
        )

        # Simple (non degenerate) RingVortex fixture for vertex proximity and
        # collinearity tests. Corners: Br=(1,0.5,0), Fr=(0,0.5,0),
        # Fl=(0,-0.5,0), Bl=(1,-0.5,0).
        (
            self.ring_Brrvp,
            self.ring_Frrvp,
            self.ring_Flrvp,
            self.ring_Blrvp,
            self.ring_strengths,
        ) = aerodynamics_functions_fixtures.make_simple_ring_vortex_arrays_fixture()

        # Simple (non degenerate) HorseshoeVortex fixture for vertex proximity
        # and collinearity tests. Corners: Br=(20,0.5,0), Fr=(0,0.5,0),
        # Fl=(0,-0.5,0), Bl=(20,-0.5,0).
        (
            self.horseshoe_Brhvp,
            self.horseshoe_Frhvp,
            self.horseshoe_Flhvp,
            self.horseshoe_Blhvp,
            self.horseshoe_strengths,
        ) = (
            aerodynamics_functions_fixtures.make_simple_horseshoe_vortex_arrays_fixture()
        )

        # Zero initial core radii for coreless comparison with the reference
        # Biot-Savart implementation.
        self.zero_rc0s = np.zeros(1, dtype=float)

    # ---- Degenerate filament tests (r0 < eps) ----

    def test_collapsed_ring_vortex_degenerate_filament_returns_zero(self):
        """Test that collapsed_velocities_from_ring_vortices returns zero
        velocity for a degenerate RingVortex where all corners coincide."""
        # Evaluate at a point away from the degenerate vortex.
        point = np.array([[1.0, 0.0, 0.0]], dtype=float)

        # Call the function.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.degenerate_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.degenerate_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.degenerate_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.degenerate_ring_Blrvp,
            strengths=self.degenerate_ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Verify all velocities are zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 3), dtype=float))

    def test_expanded_ring_vortex_degenerate_filament_returns_zero(self):
        """Test that expanded_velocities_from_ring_vortices returns zero
        velocity for a degenerate RingVortex where all corners coincide."""
        # Evaluate at a point away from the degenerate vortex.
        point = np.array([[1.0, 0.0, 0.0]], dtype=float)

        # Call the function.
        velocities = _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.degenerate_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.degenerate_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.degenerate_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.degenerate_ring_Blrvp,
            strengths=self.degenerate_ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Verify all velocities are zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 1, 3), dtype=float))

    def test_collapsed_horseshoe_vortex_degenerate_filament_returns_zero(self):
        """Test that collapsed_velocities_from_horseshoe_vortices returns zero
        velocity for a degenerate HorseshoeVortex where all points coincide."""
        # Evaluate at a point away from the degenerate vortex.
        point = np.array([[1.0, 0.0, 0.0]], dtype=float)

        # Call the function.
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=point,
                stackBrhvp_GP1_CgP1=self.degenerate_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.degenerate_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.degenerate_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.degenerate_horseshoe_Blhvp,
                strengths=self.degenerate_horseshoe_strengths,
                r_c0s=self.zero_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify all velocities are zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 3), dtype=float))

    def test_expanded_horseshoe_vortex_degenerate_filament_returns_zero(self):
        """Test that expanded_velocities_from_horseshoe_vortices returns zero
        velocity for a degenerate HorseshoeVortex where all points coincide."""
        # Evaluate at a point away from the degenerate vortex.
        point = np.array([[1.0, 0.0, 0.0]], dtype=float)

        # Call the function.
        velocities = (
            _aerodynamics_functions.expanded_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=point,
                stackBrhvp_GP1_CgP1=self.degenerate_horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.degenerate_horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.degenerate_horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.degenerate_horseshoe_Blhvp,
                strengths=self.degenerate_horseshoe_strengths,
                r_c0s=self.zero_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Verify all velocities are zero.
        npt.assert_array_almost_equal(velocities, np.zeros((1, 1, 3), dtype=float))

    # ---- Vertex proximity tests (r1/r0 < tol or r2/r0 < tol) ----

    def test_collapsed_ring_vortex_point_at_corner_matches_non_singular_legs(self):
        """Test that collapsed_velocities_from_ring_vortices returns the
        contribution from only the non singular legs when the evaluation point
        is at a RingVortex corner.

        The evaluation point is placed at the front right corner. The right leg
        (Br to Fr) hits the r2/r0 < tol guard, and the front leg (Fr to Fl)
        hits the r1/r0 < tol guard. The left and back legs are non singular.
        """
        # Place the evaluation point at the front right corner.
        point = self.ring_Frrvp.copy()
        gamma = float(self.ring_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Compute expected velocity from the two non singular legs using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        v_left = ref(self.ring_Flrvp[0], self.ring_Blrvp[0], P, gamma)
        v_back = ref(self.ring_Blrvp[0], self.ring_Brrvp[0], P, gamma)
        expected = v_left + v_back

        # Verify the kernel result matches the non singular legs' contribution.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)

    def test_expanded_ring_vortex_point_at_corner_matches_non_singular_legs(self):
        """Test that expanded_velocities_from_ring_vortices returns the
        contribution from only the non singular legs when the evaluation point
        is at a RingVortex corner.

        The evaluation point is placed at the front right corner. The right leg
        (Br to Fr) hits the r2/r0 < tol guard, and the front leg (Fr to Fl)
        hits the r1/r0 < tol guard. The left and back legs are non singular.
        """
        # Place the evaluation point at the front right corner.
        point = self.ring_Frrvp.copy()
        gamma = float(self.ring_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Compute expected velocity from the two non singular legs using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        v_left = ref(self.ring_Flrvp[0], self.ring_Blrvp[0], P, gamma)
        v_back = ref(self.ring_Blrvp[0], self.ring_Brrvp[0], P, gamma)
        expected = v_left + v_back

        # Verify the kernel result matches the non singular legs' contribution.
        npt.assert_array_almost_equal(velocities[0, 0], expected, decimal=10)

    def test_collapsed_horseshoe_vortex_point_at_corner_matches_non_singular_legs(
        self,
    ):
        """Test that collapsed_velocities_from_horseshoe_vortices returns the
        contribution from only the non singular leg when the evaluation point is
        at a HorseshoeVortex corner.

        The evaluation point is placed at the front right corner. The right leg
        (Br to Fr) hits the r2/r0 < tol guard, and the finite leg (Fr to Fl)
        hits the r1/r0 < tol guard. Only the left leg is non singular.
        """
        # Place the evaluation point at the front right corner.
        point = self.horseshoe_Frhvp.copy()
        gamma = float(self.horseshoe_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=point,
                stackBrhvp_GP1_CgP1=self.horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.horseshoe_Blhvp,
                strengths=self.horseshoe_strengths,
                r_c0s=self.zero_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Compute expected velocity from the one non singular leg using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        expected = ref(self.horseshoe_Flhvp[0], self.horseshoe_Blhvp[0], P, gamma)

        # Verify the kernel result matches the non singular leg's contribution.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)

    # ---- Collinearity tests (r3/(r1*r2) < tol) ----

    def test_collapsed_ring_vortex_point_on_leg_extension_matches_non_singular_legs(
        self,
    ):
        """Test that collapsed_velocities_from_ring_vortices returns the
        contribution from only the non singular legs when the evaluation point
        is on the extension of a RingVortex leg.

        The evaluation point is placed at (5, 0.5, 0), which is collinear with
        the right leg (Br to Fr, both at y=0.5). The right leg hits the
        r3/(r1*r2) < tol guard. The other three legs are non singular.
        """
        # Place the evaluation point on the extension of the right leg.
        point = np.array([[5.0, 0.5, 0.0]], dtype=float)
        gamma = float(self.ring_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Compute expected velocity from the three non singular legs using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        v_front = ref(self.ring_Frrvp[0], self.ring_Flrvp[0], P, gamma)
        v_left = ref(self.ring_Flrvp[0], self.ring_Blrvp[0], P, gamma)
        v_back = ref(self.ring_Blrvp[0], self.ring_Brrvp[0], P, gamma)
        expected = v_front + v_left + v_back

        # Verify the kernel result matches the non singular legs' contribution.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)

    def test_expanded_ring_vortex_point_on_leg_extension_matches_non_singular_legs(
        self,
    ):
        """Test that expanded_velocities_from_ring_vortices returns the
        contribution from only the non singular legs when the evaluation point
        is on the extension of a RingVortex leg.

        The evaluation point is placed at (5, 0.5, 0), which is collinear with
        the right leg (Br to Fr, both at y=0.5). The right leg hits the
        r3/(r1*r2) < tol guard. The other three legs are non singular.
        """
        # Place the evaluation point on the extension of the right leg.
        point = np.array([[5.0, 0.5, 0.0]], dtype=float)
        gamma = float(self.ring_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Compute expected velocity from the three non singular legs using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        v_front = ref(self.ring_Frrvp[0], self.ring_Flrvp[0], P, gamma)
        v_left = ref(self.ring_Flrvp[0], self.ring_Blrvp[0], P, gamma)
        v_back = ref(self.ring_Blrvp[0], self.ring_Brrvp[0], P, gamma)
        expected = v_front + v_left + v_back

        # Verify the kernel result matches the non singular legs' contribution.
        npt.assert_array_almost_equal(velocities[0, 0], expected, decimal=10)

    def test_collapsed_ring_vortex_point_on_leg_midpoint_matches_non_singular_legs(
        self,
    ):
        """Test that collapsed_velocities_from_ring_vortices returns the
        contribution from only the non singular legs when the evaluation point
        is at the midpoint of a RingVortex leg.

        The evaluation point is placed at (0.5, 0.5, 0), which is the midpoint
        of the right leg (Br to Fr). This triggers the r3/(r1*r2) < tol guard
        for the right leg due to collinearity. The other three legs are non
        singular.
        """
        # Place the evaluation point at the midpoint of the right leg.
        point = np.array([[0.5, 0.5, 0.0]], dtype=float)
        gamma = float(self.ring_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
        )

        # Compute expected velocity from the three non singular legs using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        v_front = ref(self.ring_Frrvp[0], self.ring_Flrvp[0], P, gamma)
        v_left = ref(self.ring_Flrvp[0], self.ring_Blrvp[0], P, gamma)
        v_back = ref(self.ring_Blrvp[0], self.ring_Brrvp[0], P, gamma)
        expected = v_front + v_left + v_back

        # Verify the kernel result matches the non singular legs' contribution.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)

    def test_collapsed_horseshoe_vortex_point_on_finite_leg_matches_non_singular_legs(
        self,
    ):
        """Test that collapsed_velocities_from_horseshoe_vortices returns the
        contribution from only the non singular legs when the evaluation point
        is on the HorseshoeVortex's finite leg.

        The evaluation point is placed at (0, 0, 0), which is the midpoint of
        the finite leg (Fr to Fl, both at x=0). This triggers the
        r3/(r1*r2) < tol guard for the finite leg due to collinearity. The
        right and left legs are non singular.
        """
        # Place the evaluation point at the midpoint of the finite leg.
        point = np.array([[0.0, 0.0, 0.0]], dtype=float)
        gamma = float(self.horseshoe_strengths[0])

        # Call the function with zero core radius for coreless comparison.
        velocities = (
            _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
                stackP_GP1_CgP1=point,
                stackBrhvp_GP1_CgP1=self.horseshoe_Brhvp,
                stackFrhvp_GP1_CgP1=self.horseshoe_Frhvp,
                stackFlhvp_GP1_CgP1=self.horseshoe_Flhvp,
                stackBlhvp_GP1_CgP1=self.horseshoe_Blhvp,
                strengths=self.horseshoe_strengths,
                r_c0s=self.zero_rc0s,
                singularity_counts=np.zeros(4, dtype=np.int64),
            )
        )

        # Compute expected velocity from the two non singular legs using the
        # reference Biot-Savart implementation.
        ref = TestAerodynamicsFunctions.ref_calculate_biot_savart_velocity
        P = point[0]
        v_right = ref(self.horseshoe_Brhvp[0], self.horseshoe_Frhvp[0], P, gamma)
        v_left = ref(self.horseshoe_Flhvp[0], self.horseshoe_Blhvp[0], P, gamma)
        expected = v_right + v_left

        # Verify the kernel result matches the non singular legs' contribution.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)


class TestCoreRadiusFormula(unittest.TestCase):
    """This is a class with functions to test the Ramasamy-Leishman core radius
    formula and its effect on induced velocities."""

    def setUp(self):
        """Set up fixtures for core radius formula tests."""
        # Simple RingVortex fixture.
        (
            self.ring_Brrvp,
            self.ring_Frrvp,
            self.ring_Flrvp,
            self.ring_Blrvp,
            self.ring_strengths,
        ) = aerodynamics_functions_fixtures.make_simple_ring_vortex_arrays_fixture()

        # Evaluation point above the RingVortex center, away from any
        # singularity.
        self.center_point = np.array([[0.5, 0.0, 1.0]], dtype=float)

    @staticmethod
    def ref_calculate_regularized_biot_savart_velocity(S_A_a, E_A_a, P_A_a, gamma, r_c):
        """Calculate induced velocity using the regularized Biot-Savart formula
        for a conceptual line vortex with a finite core radius.

        Extends the coreless Biot-Savart formula with Lamb-Oseen core
        regularization.

        r0_A = E_A_a - S_A_a
        r1_A = P_A_a - S_A_a
        r2_A = P_A_a - E_A_a
        r3_A = r1_A x r2_A

        r1 = |r1_A|
        r2 = |r2_A|
        r3 = |r3_A|

        c_1 = (gamma / (4 * pi))

        c_2_num = (r1 + r2) * (r1 * r2 - r1_A . r2_A)
        c_2_den = r1 * r2 * (r3^2 + |r2_A - r1_A|^2 * r_c^2)

        v_A__I = c_1 * (c_2_num / c_2_den) * r3_A

        :param S_A_a: A (3,) ndarray of floats representing the start point of
            the line vortex (in A axes, relative to point a) in meters.
        :param E_A_a: A (3,) ndarray of floats representing the end point of
            the line vortex (in A axes, relative to point a) in meters.
        :param P_A_a: A (3,) ndarray of floats representing the evaluation
            point (in A axes, relative to point a) in meters.
        :param gamma: A float representing the line vortex strength in meters
            squared per second.
        :param r_c: A non negative float representing the core radius in
            meters.
        :return v_A__I: A (3,) ndarray of floats representing the induced
            velocity (in A axes, observed from an inertial frame) in meters per
            second.
        """
        eps = np.finfo(float).eps
        tol = 1.0e-10

        r1_A = P_A_a - S_A_a
        r2_A = P_A_a - E_A_a
        r0_A = E_A_a - S_A_a
        r3_A = np.cross(r1_A, r2_A)

        r0 = np.linalg.norm(r0_A)
        r1 = np.linalg.norm(r1_A)
        r2 = np.linalg.norm(r2_A)
        r3 = np.linalg.norm(r3_A)

        if r0 < eps:
            return np.zeros(3, dtype=float)

        if r1 / r0 < tol or r2 / r0 < tol or r3 / (r1 * r2) < tol:
            return np.zeros(3, dtype=float)

        c_1 = gamma / (4.0 * math.pi)

        c_2_num = (r1 + r2) * (r1 * r2 - np.dot(r1_A, r2_A))
        c_2_den = r1 * r2 * (r3**2.0 + np.linalg.norm(r2_A - r1_A) ** 2.0 * r_c**2.0)

        v_A__I = c_1 * (c_2_num / c_2_den) * r3_A

        return v_A__I

    def _call_collapsed_ring(self, rc0s, ages=None, nu=0.0):
        """Helper to call collapsed_velocities_from_ring_vortices with the
        simple RingVortex fixture and the center evaluation point.
        """
        return _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.center_point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=rc0s,
            singularity_counts=np.zeros(4, dtype=np.int64),
            ages=ages,
            nu=nu,
        )

    def _compute_reference_ring_velocity(self, r_c):
        """Helper to compute the expected RingVortex velocity by summing all
        four legs' contributions using the regularized reference
        implementation.
        """
        ref = self.ref_calculate_regularized_biot_savart_velocity
        P = self.center_point[0]
        gamma = float(self.ring_strengths[0])
        Br = self.ring_Brrvp[0]
        Fr = self.ring_Frrvp[0]
        Fl = self.ring_Flrvp[0]
        Bl = self.ring_Blrvp[0]

        v_right = ref(Br, Fr, P, gamma, r_c)
        v_front = ref(Fr, Fl, P, gamma, r_c)
        v_left = ref(Fl, Bl, P, gamma, r_c)
        v_back = ref(Bl, Br, P, gamma, r_c)

        return v_right + v_front + v_left + v_back

    def test_collapsed_ring_vortex_velocity_decreases_with_increasing_rc0(self):
        """Test that increasing the initial core radius decreases the induced
        velocity magnitude for a RingVortex."""
        rc0_values = [0.0, 0.01, 0.03, 0.1, 0.5]
        magnitudes = []

        for rc0 in rc0_values:
            rc0s = np.array([rc0], dtype=float)
            velocities = self._call_collapsed_ring(rc0s)
            magnitudes.append(np.linalg.norm(velocities[0]))

        # Verify that velocity magnitude strictly decreases.
        for i in range(len(magnitudes) - 1):
            self.assertGreater(magnitudes[i], magnitudes[i + 1])

    def test_collapsed_ring_vortex_velocity_decreases_with_increasing_age(self):
        """Test that increasing the vortex age decreases the induced velocity
        magnitude for a RingVortex due to core radius growth."""
        rc0s = np.array([0.03], dtype=float)
        nu = 1.5e-5
        age_values = [0.0, 0.1, 0.5, 1.0, 5.0]
        magnitudes = []

        for age in age_values:
            ages = np.array([age], dtype=float)
            velocities = self._call_collapsed_ring(rc0s, ages=ages, nu=nu)
            magnitudes.append(np.linalg.norm(velocities[0]))

        # Verify that velocity magnitude strictly decreases.
        for i in range(len(magnitudes) - 1):
            self.assertGreater(magnitudes[i], magnitudes[i + 1])

    def test_collapsed_ring_vortex_velocity_decreases_with_increasing_viscosity(self):
        """Test that increasing the kinematic viscosity decreases the induced
        velocity magnitude for an aged RingVortex due to faster core radius
        growth."""
        rc0s = np.array([0.03], dtype=float)
        ages = np.array([1.0], dtype=float)
        nu_values = [0.0, 1.0e-5, 1.0e-4, 1.0e-3]
        magnitudes = []

        for nu in nu_values:
            velocities = self._call_collapsed_ring(rc0s, ages=ages, nu=nu)
            magnitudes.append(np.linalg.norm(velocities[0]))

        # Verify that velocity magnitude strictly decreases.
        for i in range(len(magnitudes) - 1):
            self.assertGreater(magnitudes[i], magnitudes[i + 1])

    def test_collapsed_ring_vortex_large_rc0_suppresses_velocity(self):
        """Test that a very large initial core radius suppresses the induced
        velocity to near zero."""
        rc0s = np.array([1000.0], dtype=float)

        # Call the function.
        velocities = self._call_collapsed_ring(rc0s)

        # Verify that the velocity magnitude is extremely small.
        magnitude = np.linalg.norm(velocities[0])
        self.assertLess(magnitude, 1.0e-6)

    def test_collapsed_ring_vortex_regularized_matches_reference(self):
        """Test that the regularized kernel output matches the regularized
        reference Biot-Savart implementation for a bound RingVortex with a
        nonzero initial core radius.

        With ages=None and nu=0.0, the core radius equals r_c0.
        """
        r_c0 = 0.03
        rc0s = np.array([r_c0], dtype=float)

        # For bound vortices (ages=None), r_c = r_c0.
        r_c = r_c0

        # Call the kernel.
        velocities = self._call_collapsed_ring(rc0s)

        # Compute expected velocity from all four legs using the regularized
        # reference implementation.
        expected = self._compute_reference_ring_velocity(r_c)

        # Verify the kernel result matches the reference.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)

    def test_collapsed_ring_vortex_aged_core_radius_matches_reference(self):
        """Test that the regularized kernel output matches the regularized
        reference Biot-Savart implementation for an aged RingVortex.

        With age=1.0, nu=1.5e-5, and strength=1.0, the core radius grows from
        r_c0 according to the Ramasamy-Leishman formula.
        """
        r_c0 = 0.03
        age = 1.0
        nu = 1.5e-5
        gamma = float(self.ring_strengths[0])
        rc0s = np.array([r_c0], dtype=float)
        ages = np.array([age], dtype=float)

        # Manually compute r_c using the Ramasamy-Leishman formula with the
        # module's physical constants.
        lamb = _aerodynamics_functions._lamb
        squire = _aerodynamics_functions._squire
        r_c = np.sqrt(r_c0**2 + 4.0 * lamb * (nu + squire * abs(gamma)) * age)

        # Call the kernel.
        velocities = self._call_collapsed_ring(rc0s, ages=ages, nu=nu)

        # Compute expected velocity from all four legs using the regularized
        # reference implementation.
        expected = self._compute_reference_ring_velocity(r_c)

        # Verify the kernel result matches the reference.
        npt.assert_array_almost_equal(velocities[0], expected, decimal=10)


class TestSingularityCounters(unittest.TestCase):
    """This is a class with functions to test that singularity counters increment
    correctly for known singular configurations."""

    def setUp(self):
        """Set up fixtures for singularity counter tests."""
        # Simple (non degenerate) RingVortex fixture.
        (
            self.ring_Brrvp,
            self.ring_Frrvp,
            self.ring_Flrvp,
            self.ring_Blrvp,
            self.ring_strengths,
        ) = aerodynamics_functions_fixtures.make_simple_ring_vortex_arrays_fixture()

        # Degenerate RingVortex fixture (all corners at the origin).
        (
            self.degenerate_ring_Brrvp,
            self.degenerate_ring_Frrvp,
            self.degenerate_ring_Flrvp,
            self.degenerate_ring_Blrvp,
            self.degenerate_ring_strengths,
        ) = aerodynamics_functions_fixtures.make_degenerate_ring_vortex_arrays_fixture()

        # Simple (non degenerate) HorseshoeVortex fixture.
        (
            self.horseshoe_Brhvp,
            self.horseshoe_Frhvp,
            self.horseshoe_Flhvp,
            self.horseshoe_Blhvp,
            self.horseshoe_strengths,
        ) = (
            aerodynamics_functions_fixtures.make_simple_horseshoe_vortex_arrays_fixture()
        )

        # Zero initial core radii for coreless tests.
        self.zero_rc0s = np.zeros(1, dtype=float)

    def test_degenerate_filament_increments_counter_zero(self):
        """Test that a degenerate RingVortex (all corners at origin) increments
        singularity_counts[0] (degenerate filament)."""
        singularity_counts = np.zeros(4, dtype=np.int64)
        point = np.array([[1.0, 0.0, 0.0]], dtype=float)

        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.degenerate_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.degenerate_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.degenerate_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.degenerate_ring_Blrvp,
            strengths=self.degenerate_ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # All four legs are degenerate filaments.
        self.assertGreater(singularity_counts[0], 0)

    def test_vertex_proximity_increments_counter_one_and_two(self):
        """Test that placing the evaluation point at a RingVortex corner
        increments singularity_counts[1] (vertex start proximity) and
        singularity_counts[2] (vertex end proximity)."""
        singularity_counts = np.zeros(4, dtype=np.int64)

        # Place the evaluation point at the front right corner.
        point = self.ring_Frrvp.copy()

        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # The right leg (Br to Fr) hits the r2/r0 guard (counter 2) because
        # the evaluation point is at Fr. The front leg (Fr to Fl) hits the
        # r1/r0 guard (counter 1) because the evaluation point is at Fr.
        self.assertGreater(singularity_counts[1], 0)
        self.assertGreater(singularity_counts[2], 0)

    def test_on_filament_collinearity_increments_counter_three(self):
        """Test that placing the evaluation point on the a RingVortex leg increments
        singularity_counts[3] (on-filament collinearity)."""
        singularity_counts = np.zeros(4, dtype=np.int64)

        # Place the evaluation point on the right leg.
        point = np.array([[0.5, 0.5, 0.0]], dtype=float)

        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # The right leg (Br to Fr) is collinear with the evaluation point.
        self.assertGreater(singularity_counts[3], 0)

    def test_off_filament_collinearity_does_not_increment_counter(self):
        """Test that placing the evaluation point on the extension of a RingVortex leg
        doesn't increment singularity_counts[3] (off-filament collinearity)."""
        singularity_counts = np.zeros(4, dtype=np.int64)

        # Place the evaluation point on the extension of the right leg.
        point = np.array([[5.0, 0.5, 0.0]], dtype=float)

        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # The point is collinear with the right leg (Br to Fr) but lies off the
        # filament (c_3 > 0), so the counter is not incremented.
        self.assertEqual(singularity_counts[3], 0)

    def test_non_singular_configuration_has_zero_counts(self):
        """Test that a non singular configuration produces zero singularity
        counts."""
        singularity_counts = np.zeros(4, dtype=np.int64)

        # Place the evaluation point well above the RingVortex.
        point = np.array([[0.5, 0.0, 5.0]], dtype=float)

        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # No singularity events should have occurred.
        self.assertEqual(singularity_counts.sum(), 0)

    def test_counts_accumulate_across_calls(self):
        """Test that singularity counts accumulate when the same array is
        passed to multiple calls."""
        singularity_counts = np.zeros(4, dtype=np.int64)
        point = np.array([[1.0, 0.0, 0.0]], dtype=float)

        # First call: degenerate RingVortex.
        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.degenerate_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.degenerate_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.degenerate_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.degenerate_ring_Blrvp,
            strengths=self.degenerate_ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        first_total = singularity_counts.sum()
        self.assertGreater(first_total, 0)

        # Second call: same degenerate RingVortex.
        _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.degenerate_ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.degenerate_ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.degenerate_ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.degenerate_ring_Blrvp,
            strengths=self.degenerate_ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # The total should have at least doubled.
        self.assertGreaterEqual(singularity_counts.sum(), 2 * first_total)

    def test_horseshoe_vortex_counter_increments(self):
        """Test that singularity counters work correctly for HorseshoeVortex
        wrapper functions."""
        singularity_counts = np.zeros(4, dtype=np.int64)

        # Place the evaluation point at the front right corner.
        point = self.horseshoe_Frhvp.copy()

        _aerodynamics_functions.collapsed_velocities_from_horseshoe_vortices(
            stackP_GP1_CgP1=point,
            stackBrhvp_GP1_CgP1=self.horseshoe_Brhvp,
            stackFrhvp_GP1_CgP1=self.horseshoe_Frhvp,
            stackFlhvp_GP1_CgP1=self.horseshoe_Flhvp,
            stackBlhvp_GP1_CgP1=self.horseshoe_Blhvp,
            strengths=self.horseshoe_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # At least one singularity event should have been recorded.
        self.assertGreater(singularity_counts.sum(), 0)

    def test_expanded_ring_vortex_counter_increments(self):
        """Test that singularity counters work correctly for the expanded
        RingVortex wrapper."""
        singularity_counts = np.zeros(4, dtype=np.int64)

        # Place the evaluation point at the front right corner.
        point = self.ring_Frrvp.copy()

        _aerodynamics_functions.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=point,
            stackBrrvp_GP1_CgP1=self.ring_Brrvp,
            stackFrrvp_GP1_CgP1=self.ring_Frrvp,
            stackFlrvp_GP1_CgP1=self.ring_Flrvp,
            stackBlrvp_GP1_CgP1=self.ring_Blrvp,
            strengths=self.ring_strengths,
            r_c0s=self.zero_rc0s,
            singularity_counts=singularity_counts,
        )

        # At least one singularity event should have been recorded.
        self.assertGreater(singularity_counts.sum(), 0)


class TestLogSingularityCounts(unittest.TestCase):
    """This is a class with functions to test the _log_singularity_counts
    helper function."""

    def test_zero_counts_does_not_log(self):
        """Test that _log_singularity_counts does not log when all counts are
        zero."""
        import logging

        from pterasoftware._functions import log_unexpected_singularity_counts

        logger = logging.getLogger("test_zero_counts")
        singularity_counts = np.zeros(4, dtype=np.int64)

        with self.assertLogs(logger, level="DEBUG") as cm:
            logger.debug("sentinel")
            log_unexpected_singularity_counts(
                logger, logging.WARNING, "test_context", singularity_counts
            )

        # Only the sentinel message should have been logged.
        self.assertEqual(len(cm.output), 1)
        self.assertIn("sentinel", cm.output[0])

    def test_nonzero_counts_logs_message(self):
        """Test that _log_singularity_counts logs a message when counts are
        nonzero."""
        import logging

        from pterasoftware._functions import log_unexpected_singularity_counts

        logger = logging.getLogger("test_nonzero_counts")
        singularity_counts = np.array([2, 0, 3, 0], dtype=np.int64)

        with self.assertLogs(logger, level="WARNING") as cm:
            log_unexpected_singularity_counts(
                logger, logging.WARNING, "test_context", singularity_counts
            )

        # Verify the log message contains the context and the total count.
        self.assertEqual(len(cm.output), 1)
        self.assertIn("test_context", cm.output[0])
        self.assertIn("5 singularity skip(s)", cm.output[0])
        self.assertIn("degenerate filament=2", cm.output[0])
        self.assertIn("vertex end proximity=3", cm.output[0])

    def test_log_level_is_respected(self):
        """Test that _log_singularity_counts uses the specified logging
        level."""
        import logging

        from pterasoftware._functions import log_unexpected_singularity_counts

        logger = logging.getLogger("test_log_level")
        singularity_counts = np.array([1, 0, 0, 0], dtype=np.int64)

        with self.assertLogs(logger, level="ERROR") as cm:
            log_unexpected_singularity_counts(
                logger, logging.ERROR, "error_context", singularity_counts
            )

        self.assertEqual(len(cm.output), 1)
        self.assertIn("ERROR", cm.output[0])


if __name__ == "__main__":
    unittest.main()
