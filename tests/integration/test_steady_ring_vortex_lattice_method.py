"""This module is a testing case for the SteadyRingVortexLatticeMethodSolver.

Based on an identical XFLR5 VLM2 testing case, the expected output for this case is: CL:
0.784     CDi:    0.019     Cm:     -0.678
"""

import unittest
from unittest.mock import patch

import numpy as np
import pyvista as pv

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestSteadyRingVortexLatticeMethod(unittest.TestCase):
    """This is a class for testing the SteadyRingVortexLatticeMethodSolver."""

    steady_ring_vortex_lattice_method_validation_solver: (
        ps.steady_ring_vortex_lattice_method.SteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """This method sets up the test.

        The solver is run once here and shared by every test, none of which change it.

        :return: None
        """
        cls.steady_ring_vortex_lattice_method_validation_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )
        cls.steady_ring_vortex_lattice_method_validation_solver.run()

    def test_method(self) -> None:
        """This method tests the SteadyRingVortexLatticeMethodSolver's output.

        It also tests that the solver doesn't throw an error when the draw function is
        called using it.

        :return: None
        """
        this_airplane = (
            self.steady_ring_vortex_lattice_method_validation_solver.airplanes[0]
        )
        assert this_airplane.forceCoefficients_W is not None
        assert this_airplane.momentCoefficients_W_CgP1 is not None

        # Calculate the percent errors of the output.
        c_di_expected = 0.019
        c_di_calculated = -this_airplane.forceCoefficients_W[0]
        c_di_error = abs((c_di_calculated - c_di_expected) / c_di_expected)

        c_l_expected = 0.784
        c_l_calculated = -this_airplane.forceCoefficients_W[2]
        c_l_error = abs((c_l_calculated - c_l_expected) / c_l_expected)

        c_m_expected = -0.678
        c_m_calculated = this_airplane.momentCoefficients_W_CgP1[1]
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

    def test_diagram_shows_the_diagram(self) -> None:
        """Test that diagram shows the diagram once.

        :return: None
        """
        # Patch the Plotter's show method to avoid blocking on window close.
        with patch.object(pv.Plotter, "show") as mock_show:
            self.steady_ring_vortex_lattice_method_validation_solver.diagram()

        mock_show.assert_called_once()

    def test_diagram_shows_the_simplified_vortices(self) -> None:
        """Test that diagram shows the diagram once with the vortices simplified.

        :return: None
        """
        with patch.object(pv.Plotter, "show") as mock_show:
            self.steady_ring_vortex_lattice_method_validation_solver.diagram(
                simplify_vortices=True
            )

        mock_show.assert_called_once()

    def test_diagram_accepts_numpy_bools(self) -> None:
        """Test that diagram accepts numpy bools for its flags.

        :return: None
        """
        with patch.object(pv.Plotter, "show") as mock_show:
            self.steady_ring_vortex_lattice_method_validation_solver.diagram(
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

    def test_diagram_raises_before_the_solver_runs(self) -> None:
        """Test that diagram raises a RuntimeError for a solver that hasn't run, since
        the vortices it draws are placed during the run.

        :return: None
        """
        unrun_solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )
        with patch.object(pv.Plotter, "show") as mock_show:
            with self.assertRaises(RuntimeError):
                unrun_solver.diagram()

        mock_show.assert_not_called()
