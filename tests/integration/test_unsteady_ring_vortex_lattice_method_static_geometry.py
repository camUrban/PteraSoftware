"""This is a testing case for the UnsteadyRingVortexLatticeMethodSolver with static
geometry.

Based on an equivalent XFLR5 testing case, the expected output for this case is:     CL:
0.485     CDi:    0.015     Cm:     -0.166

The expected output was created using XFLR5's inviscid VLM2 analysis type, which is a
ring vortex lattice method solver. The geometry in this case is static. Therefore, the
results of this unsteady solver should converge to be close to XFLR5's static result.
"""

import unittest
from unittest.mock import patch

import numpy as np
import pyvista as pv

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestUnsteadyRingVortexLatticeMethodStaticGeometry(unittest.TestCase):
    """This is a class for testing the UnsteadyRingVortexLatticeMethodSolver on static
    geometry."""

    unsteady_ring_vortex_lattice_method_validation_solver: (
        ps.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """This method sets up the test.

        The solver is run once here and shared by every test, since the run dominates
        this module's time and none of the tests change the solver.

        :return: None
        """
        cls.unsteady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.unsteady_ring_vortex_lattice_method_validation_solver.run(
            prescribed_wake=True,
            show_progress=False,
        )

    def test_method(self) -> None:
        """This method tests the UnsteadyRingVortexLatticeMethodSolver's output.

        It also tests that the solver doesn't throw an error when the animate and
        plot_results_versus_time functions are called using it.

        :return: None
        """
        this_solver = self.unsteady_ring_vortex_lattice_method_validation_solver
        this_airplane = this_solver.current_airplanes[0]
        assert this_airplane.forceCoefficients_W is not None
        assert this_airplane.momentCoefficients_W_CgP1 is not None

        # Calculate the percent errors of the output.
        c_di_expected = 0.015
        c_di_calculated = -this_airplane.forceCoefficients_W[0]
        c_di_error = abs(c_di_calculated - c_di_expected) / c_di_expected

        c_l_expected = 0.485
        c_l_calculated = -this_airplane.forceCoefficients_W[2]
        c_l_error = abs(c_l_calculated - c_l_expected) / c_l_expected

        c_m_expected = -0.166
        c_m_calculated = this_airplane.momentCoefficients_W_CgP1[1]
        c_m_error = abs(c_m_calculated - c_m_expected) / c_m_expected

        # Set the allowable percent error.
        allowable_error = 0.10

        ps.animate(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show_wake_vortices=True,
            scalar_type="lift",
            save=False,
            testing=True,
        )

        ps.plot_results_versus_time(
            unsteady_solver=self.unsteady_ring_vortex_lattice_method_validation_solver,
            show=False,
            save=False,
        )

        # Assert that the percent errors are less than the allowable error.
        self.assertTrue(abs(c_di_error) < allowable_error)
        self.assertTrue(abs(c_l_error) < allowable_error)
        self.assertTrue(abs(c_m_error) < allowable_error)

    def test_diagram_shows_the_diagram(self) -> None:
        """Test that diagram shows the last time step's diagram once.

        :return: None
        """
        # Patch the Plotter's show method to avoid blocking on window close.
        with patch.object(pv.Plotter, "show") as mock_show:
            self.unsteady_ring_vortex_lattice_method_validation_solver.diagram()

        mock_show.assert_called_once()

    def test_diagram_shows_the_simplified_vortices(self) -> None:
        """Test that diagram shows the diagram once with the vortices simplified.

        Each simplified vortex leg adds its own arrow tip, so this draws an early time
        step, whose short wake keeps the test fast while still simplifying both bound
        and wake ring vortices.

        :return: None
        """
        this_solver = self.unsteady_ring_vortex_lattice_method_validation_solver
        self.assertGreater(len(this_solver.listStackFrwrvp_GP1_CgP1[2]), 0)

        with patch.object(pv.Plotter, "show") as mock_show:
            this_solver.diagram(step=2, simplify_vortices=True)

        mock_show.assert_called_once()

    def test_diagram_shows_the_first_time_step(self) -> None:
        """Test that diagram shows the first time step's diagram once.

        No wake ring vortices have been shed at the first time step, so this draws the
        bound ring vortices with an empty stack of wake ring vortices.

        :return: None
        """
        with patch.object(pv.Plotter, "show") as mock_show:
            self.unsteady_ring_vortex_lattice_method_validation_solver.diagram(step=0)

        mock_show.assert_called_once()

    def test_diagram_accepts_numpy_bools(self) -> None:
        """Test that diagram accepts numpy bools for its flags.

        :return: None
        """
        with patch.object(pv.Plotter, "show") as mock_show:
            self.unsteady_ring_vortex_lattice_method_validation_solver.diagram(
                show_airplane_axes_and_points=np.bool(True),
                show_wing_axes_and_points=np.bool(True),
                show_wing_cross_section_axes_and_points=np.bool(True),
                show_airfoil_axes_and_points=np.bool(True),
                show_airfoils=np.bool(True),
                show_mcls=np.bool(True),
                show_collocation_points=np.bool(True),
                label_collocation_points=np.bool(True),
                simplify_vortices=np.bool(False),
                math_labels=np.bool(True),
                save=np.bool(False),
            )

        mock_show.assert_called_once()

    def test_diagram_rejects_an_out_of_range_step(self) -> None:
        """Test that diagram rejects a step past either end of the time steps.

        :return: None
        """
        num_steps = self.unsteady_ring_vortex_lattice_method_validation_solver.num_steps
        for step in [num_steps, -num_steps - 1]:
            with self.subTest(step=step):
                with self.assertRaises(ValueError):
                    self.unsteady_ring_vortex_lattice_method_validation_solver.diagram(
                        step=step
                    )

    def test_diagram_raises_before_the_solver_runs(self) -> None:
        """Test that diagram raises a RuntimeError for a solver that hasn't run, since
        the vortices it draws are placed during the run.

        :return: None
        """
        unrun_solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        with patch.object(pv.Plotter, "show") as mock_show:
            with self.assertRaises(RuntimeError):
                unrun_solver.diagram()

        mock_show.assert_not_called()
