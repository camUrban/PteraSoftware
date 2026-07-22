"""This module tests that the Joukowski and Katz force methods produce comparable
results for the UnsteadyRingVortexLatticeMethodSolver.

The tests verify that both methods complete without errors on static and variable
geometry cases, and that they produce forces within a reasonable tolerance of each
other.
"""

import unittest

import numpy as np

import pterasoftware as ps
from tests.integration.fixtures import solver_fixtures


class TestForceMethodsComparison(unittest.TestCase):
    """This is a class for comparing the Joukowski and Katz force calculation
    methods."""

    static_solver_joukowski: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )
    static_solver_katz: (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
    )

    @classmethod
    def setUpClass(cls) -> None:
        """Set up test fixtures once for all tests."""
        # Create solver for static geometry case.
        cls.static_solver_joukowski = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )
        cls.static_solver_katz = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_static_geometry()
        )

        # Run both solvers.
        cls.static_solver_joukowski.run(
            prescribed_wake=True,
            show_progress=False,
            force_method="joukowski",
        )
        cls.static_solver_katz.run(
            prescribed_wake=True,
            show_progress=False,
            force_method="katz",
        )

    def test_static_geometry_methods_produce_similar_lift(self) -> None:
        """Test that both methods produce lift within 25% of each other."""
        joukowski_airplane = self.static_solver_joukowski.current_airplanes[0]
        katz_airplane = self.static_solver_katz.current_airplanes[0]
        assert joukowski_airplane.forceCoefficients_W is not None
        assert katz_airplane.forceCoefficients_W is not None

        c_l_joukowski = -joukowski_airplane.forceCoefficients_W[2]
        c_l_katz = -katz_airplane.forceCoefficients_W[2]

        # Calculate relative difference.
        avg_c_l = (abs(c_l_joukowski) + abs(c_l_katz)) / 2
        if avg_c_l > 0:
            relative_diff = abs(c_l_joukowski - c_l_katz) / avg_c_l
        else:
            relative_diff = 0.0

        # Both methods compute lift from the same bound circulation distribution, and
        # the Katz method's chordwise pressure term telescopes over each chordwise strip
        # to the same Kutta-Joukowski strip total that the Joukowski method integrates,
        # so the two methods should agree closely on lift. This tolerance also guards
        # the Panel-length denominator in the chordwise vorticity density: dividing by
        # center-to-center distances instead double counts the leading edge Panels'
        # circulation and pushes this difference to about 24%.
        allowable_difference = 0.05
        self.assertLess(
            relative_diff,
            allowable_difference,
            f"Lift coefficients differ by {relative_diff:.1%}: "
            f"Joukowski={c_l_joukowski:.4f}, Katz={c_l_katz:.4f}",
        )

    def test_static_geometry_methods_produce_similar_drag(self) -> None:
        """Test that both methods produce drag within 100% of each other."""
        joukowski_airplane = self.static_solver_joukowski.current_airplanes[0]
        katz_airplane = self.static_solver_katz.current_airplanes[0]
        assert joukowski_airplane.forceCoefficients_W is not None
        assert katz_airplane.forceCoefficients_W is not None

        c_di_joukowski = -joukowski_airplane.forceCoefficients_W[0]
        c_di_katz = -katz_airplane.forceCoefficients_W[0]

        # Calculate relative difference.
        avg_c_di = (abs(c_di_joukowski) + abs(c_di_katz)) / 2
        if avg_c_di > 0:
            relative_diff = abs(c_di_joukowski - c_di_katz) / avg_c_di
        else:
            relative_diff = 0.0

        # The two methods agree on this case's drag to within about 7%. Drag is still
        # the most sensitive quantity to the force calculation method, so allow a 25%
        # difference. This tolerance guards the two pieces of the Katz method's induced
        # drag that pressure integration cannot check: omitting the trailing edge back
        # legs from the drag downwash drives the Katz drag negative (a difference of
        # about 200%), and building the decomposition axes from the full solution
        # velocity instead of the kinematic velocity double counts the induced drag (a
        # difference of about 69%).
        allowable_difference = 0.25
        self.assertLess(
            relative_diff,
            allowable_difference,
            f"Drag coefficients differ by {relative_diff:.1%}: "
            f"Joukowski={c_di_joukowski:.4f}, Katz={c_di_katz:.4f}",
        )

    def test_static_geometry_methods_produce_similar_drag_transients(self) -> None:
        """Test that both methods produce similar drag transients after the impulsive
        start."""
        # Collect each solver's drag coefficient history over all time steps.
        joukowski_history = []
        for steady_problem in self.static_solver_joukowski.steady_problems:
            airplane = steady_problem.airplanes[0]
            assert airplane.forceCoefficients_W is not None
            joukowski_history.append(-airplane.forceCoefficients_W[0])
        katz_history = []
        for steady_problem in self.static_solver_katz.steady_problems:
            airplane = steady_problem.airplanes[0]
            assert airplane.forceCoefficients_W is not None
            katz_history.append(-airplane.forceCoefficients_W[0])

        # Isolate each method's transient by subtracting its own converged value, and
        # drop the first time step, whose drag is dominated by the impulsive start.
        # Comparing the transients rather than the raw histories removes the steady
        # state offset between the methods, which the drag comparison test above already
        # covers.
        joukowski_transient = np.array(joukowski_history)
        joukowski_transient = (joukowski_transient - joukowski_transient[-1])[1:]
        katz_transient = np.array(katz_history)
        katz_transient = (katz_transient - katz_transient[-1])[1:]

        # Calculate the RMS difference between the transients, normalized by the RMS of
        # the Joukowski transient.
        relative_diff = float(
            np.sqrt(np.mean((katz_transient - joukowski_transient) ** 2))
            / np.sqrt(np.mean(joukowski_transient**2))
        )

        # The transient decay is governed by the unsteady term in each method's load
        # calculation, so this comparison guards the sign of the unsteady term in the
        # Katz method's induced drag. With the correct sign, this difference is about
        # 0.9, reflecting the Katz method's larger transient amplitude. With the wrong
        # sign, the Katz transient rises sharply after the impulsive start instead of
        # decaying, and this difference grows to about 5.4.
        allowable_difference = 2.0
        self.assertLess(
            relative_diff,
            allowable_difference,
            f"Drag transients differ by {relative_diff:.1%} of the Joukowski "
            f"transient's RMS, which suggests a sign error in the unsteady term of "
            f"the Katz method's induced drag calculation.",
        )

    def test_static_geometry_methods_produce_similar_moment(self) -> None:
        """Test that both methods produce moment within 50% of each other."""
        joukowski_airplane = self.static_solver_joukowski.current_airplanes[0]
        katz_airplane = self.static_solver_katz.current_airplanes[0]
        assert joukowski_airplane.momentCoefficients_W_CgP1 is not None
        assert katz_airplane.momentCoefficients_W_CgP1 is not None

        c_m_joukowski = joukowski_airplane.momentCoefficients_W_CgP1[1]
        c_m_katz = katz_airplane.momentCoefficients_W_CgP1[1]

        # Calculate relative difference.
        avg_c_m = (abs(c_m_joukowski) + abs(c_m_katz)) / 2
        if avg_c_m > 0:
            relative_diff = abs(c_m_joukowski - c_m_katz) / avg_c_m
        else:
            relative_diff = 0.0

        # Allow 50% difference for moment.
        allowable_difference = 0.50
        self.assertLess(
            relative_diff,
            allowable_difference,
            f"Moment coefficients differ by {relative_diff:.1%}: "
            f"Joukowski={c_m_joukowski:.4f}, Katz={c_m_katz:.4f}",
        )


class TestForceMethodsVariableGeometry(unittest.TestCase):
    """This is a class for testing that both force methods work with variable
    geometry."""

    def test_variable_geometry_joukowski_completes(self) -> None:
        """Test that Joukowski method completes on variable geometry."""
        solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry()
        )
        solver.run(
            prescribed_wake=True,
            show_progress=False,
            force_method="joukowski",
        )
        self.assertTrue(solver.ran)

    def test_variable_geometry_katz_completes(self) -> None:
        """Test that Katz method completes on variable geometry."""
        solver = (
            solver_fixtures.make_unsteady_ring_vortex_lattice_method_validation_solver_with_variable_geometry()
        )
        solver.run(
            prescribed_wake=True,
            show_progress=False,
            force_method="katz",
        )
        self.assertTrue(solver.ran)
