"""This module tests that output functions accept deserialized solver objects.

Note: These tests do not verify specific output values. They verify that the output
functions in output.py accept solver objects that have been serialized and deserialized
via save and load without throwing any errors.
"""

import tempfile
import unittest
from pathlib import Path

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestSteadySolverSerializationOutput(unittest.TestCase):
    """This class contains methods for testing that output functions accept a
    deserialized SteadyRingVortexLatticeMethodSolver.
    """

    @classmethod
    def setUpClass(cls):
        """Solve a steady ring solver, save it, and load it back.

        :return: None
        """
        solver = (
            solver_fixtures.make_steady_ring_vortex_lattice_method_validation_solver()
        )
        solver.run()

        cls.temporary_directory = tempfile.TemporaryDirectory()
        path = Path(cls.temporary_directory.name) / "solver.json"
        ps.save(path, solver)
        cls.loaded_solver = ps.load(path)

    @classmethod
    def tearDownClass(cls):
        """Clean up the temporary directory.

        :return: None
        """
        cls.temporary_directory.cleanup()

    def test_draw_with_scalar_does_not_throw(self):
        """Tests that the draw function accepts a deserialized steady solver with a
        scalar type that requires Panel force data.

        :return: None
        """
        ps.output.draw(
            solver=self.loaded_solver,
            scalar_type="lift",
            show_wake_vortices=False,
            show_streamlines=True,
            testing=True,
        )

    def test_log_results_does_not_throw(self):
        """Tests that the log_results function accepts a deserialized steady solver.

        :return: None
        """
        ps.output.log_results(solver=self.loaded_solver)


class TestUnsteadySolverSerializationOutput(unittest.TestCase):
    """This class contains methods for testing that output functions accept a
    deserialized UnsteadyRingVortexLatticeMethodSolver.
    """

    @classmethod
    def setUpClass(cls):
        """Solve an unsteady solver, save it, and load it back.

        :return: None
        """
        solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        solver.run(show_progress=False)

        cls.temporary_directory = tempfile.TemporaryDirectory()
        path = Path(cls.temporary_directory.name) / "solver.json"
        ps.save(path, solver)
        cls.loaded_solver = ps.load(path)

    @classmethod
    def tearDownClass(cls):
        """Clean up the temporary directory.

        :return: None
        """
        cls.temporary_directory.cleanup()

    def test_draw_with_scalar_does_not_throw(self):
        """Tests that the draw function accepts a deserialized unsteady solver with a
        scalar type that requires Panel force data.

        :return: None
        """
        ps.output.draw(
            solver=self.loaded_solver,
            scalar_type="lift",
            show_wake_vortices=False,
            show_streamlines=False,
            testing=True,
        )

    def test_animate_with_scalar_does_not_throw(self):
        """Tests that the animate function accepts a deserialized unsteady solver with
        a scalar type that requires Panel force data.

        :return: None
        """
        ps.output.animate(
            unsteady_solver=self.loaded_solver,
            scalar_type="lift",
            show_wake_vortices=True,
            save=False,
            testing=True,
        )

    def test_plot_results_versus_time_does_not_throw(self):
        """Tests that the plot_results_versus_time function accepts a deserialized
        unsteady solver.

        :return: None
        """
        ps.output.plot_results_versus_time(
            unsteady_solver=self.loaded_solver, show=False
        )

    def test_log_results_does_not_throw(self):
        """Tests that the log_results function accepts a deserialized unsteady solver.

        :return: None
        """
        ps.output.log_results(solver=self.loaded_solver)
