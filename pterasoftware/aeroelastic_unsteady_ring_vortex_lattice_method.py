"""Contains the AeroelasticUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

AeroelasticUnsteadyRingVortexLatticeMethodSolver: A subclass of
CoupledUnsteadyRingVortexLatticeMethodSolver that solves AeroelasticUnsteadyProblems,
extending the coupled solver with Strip Leading Edge Point (SLEP) functionality for
computing aerodynamic moments about the strip leading edge so that wing deformations can
be coupled with aerodynamic loads.

**Contains the following functions:**

None
"""

from __future__ import annotations

from typing import cast

import numpy as np

from . import _functions, problems
from ._coupled_unsteady_ring_vortex_lattice_method import (
    CoupledUnsteadyRingVortexLatticeMethodSolver,
)


class AeroelasticUnsteadyRingVortexLatticeMethodSolver(
    CoupledUnsteadyRingVortexLatticeMethodSolver
):
    """A subclass of CoupledUnsteadyRingVortexLatticeMethodSolver that adds SLEP (Strip
    Leading Edge Point) functionality for aeroelastic simulations.

    This solver extends the coupled solver with moment calculations about each panel's
    strip leading edge point, which is important for analyzing wing loading and
    deformation characteristics relative to the wing root.

    **Key additions over parent CoupledUnsteadyRingVortexLatticeMethodSolver:**
    initializes and maintains SLEP index mapping and position arrays, overrides
    ``_reinitialize_step_arrays_hook`` to reset SLEP arrays each step, overrides
    ``_load_calculation_moment_processing_hook`` to compute SLEP moments, and computes
    bound vortex positions relative to strip leading edge points.
    """

    __slots__ = (
        "slep_point_indices",
        "stackCblvpr_GP1_Slep",
        "stackCblvpf_GP1_Slep",
        "stackCblvpl_GP1_Slep",
        "stackCblvpb_GP1_Slep",
        "stackCpp_GP1_Slep",
        "stackFlpp_GP1_CgP1",
        "moments_GP1_Slep",
        "stack_leading_edge_points",
    )

    def __init__(
        self,
        aeroelastic_unsteady_problem: problems.AeroelasticUnsteadyProblem,
    ) -> None:
        """Initialize the solver for an AeroelasticUnsteadyProblem.

        Sets up the solver infrastructure and initializes SLEP (Strip Leading Edge
        Point) related attributes.

        :param aeroelastic_unsteady_problem: The AeroelasticUnsteadyProblem to be
            solved.
        :return: None
        """
        if not isinstance(
            aeroelastic_unsteady_problem, problems.AeroelasticUnsteadyProblem
        ):
            raise TypeError(
                "aeroelastic_unsteady_problem must be an " "AeroelasticUnsteadyProblem."
            )

        super().__init__(aeroelastic_unsteady_problem)

        first_steady_problem: problems.SteadyProblem = self._get_steady_problem_at(0)

        # Initialize SLEP (Strip Leading Edge Point) information. For each airplane
        # and wing, we track the panel index where each new spanwise strip begins.
        # This allows efficient computation of moments about the strip leading edge
        # (wing root to tip).
        panel_count = 0
        slep_point_indices_list: list[int] = []
        for airplane in first_steady_problem.airplanes:
            for wing in airplane.wings:
                for wing_cross_section in wing.wing_cross_sections:
                    # Record the first panel index for this wing cross-section
                    # (start of strip).
                    slep_point_indices_list.append(panel_count)
                    if wing_cross_section.num_spanwise_panels is not None:
                        panel_count += wing_cross_section.num_spanwise_panels
        self.slep_point_indices: np.ndarray = np.array(
            slep_point_indices_list, dtype=int
        )

        # The current time step's center bound LineVortex points for the right,
        # front, left, and back legs (in the first Airplane's geometry axes,
        # relative to the local strip leading edge point).
        self.stackCblvpr_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpf_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpl_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpb_GP1_Slep: np.ndarray = np.empty(0, dtype=float)

        # The colocation panel points and the front left panel point (in the first
        # Airplane's geometry axes, relative to the local strip leading edge point
        # and the first Airplane's CG respectively).
        self.stackCpp_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackFlpp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.moments_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stack_leading_edge_points: np.ndarray = np.empty(0, dtype=float)

    @property
    def _aeroelastic_problem(self) -> problems.AeroelasticUnsteadyProblem:
        """The solver's AeroelasticUnsteadyProblem, narrowed from the inherited
        unsteady_problem.

        The inherited unsteady_problem slot is typed as the base CoreUnsteadyProblem so
        the parent solver can hold any coupled problem. This solver's constructor only
        accepts an AeroelasticUnsteadyProblem, so the cast here is safe.

        :return: This solver's AeroelasticUnsteadyProblem.
        """
        return cast(problems.AeroelasticUnsteadyProblem, self.unsteady_problem)

    def _reinitialize_step_arrays_hook(self) -> None:
        """Reinitialize SLEP arrays at the start of each time step.

        :return: None
        """
        self.stackCblvpr_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpf_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpl_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCblvpb_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.stackCpp_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.moments_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.stackFlpp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

    def _load_calculation_moment_processing_hook(
        self,
        rightLegForces_GP1,
        frontLegForces_GP1,
        leftLegForces_GP1,
        backLegForces_GP1,
        unsteady_forces_GP1,
    ) -> np.ndarray:
        """Override parent to compute moments about both center-of-gravity and SLEP.

        This hook extends the parent class's moment calculation by additionally
        computing moments about each panel's Strip Leading Edge Point (SLEP). This is
        used for analyzing wing loading and deformation characteristics relative to the
        wing root.

        The method first calls the parent's implementation to get CG-based moments, then
        updates bound vortex positions relative to SLEP points, recalculates all moment
        contributions in the SLEP frame, and stores the SLEP moments in
        self.moments_GP1_Slep.

        :return: moments_GP1_CgP1, a (N,3) ndarray of floats representing the moments
            (in the first Airplane's geometry axes, relative to the first Airplane's CG)
            on every Panel at the current time step. SLEP moments are stored separately
            in self.moments_GP1_Slep.
        """
        moments_GP1_CgP1 = super()._load_calculation_moment_processing_hook(
            rightLegForces_GP1,
            frontLegForces_GP1,
            leftLegForces_GP1,
            backLegForces_GP1,
            unsteady_forces_GP1,
        )

        self._update_bound_vortex_positions_relative_to_slep_points()

        rightLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self.stackCblvpr_GP1_Slep, rightLegForces_GP1
        )
        frontLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self.stackCblvpf_GP1_Slep, frontLegForces_GP1
        )
        leftLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self.stackCblvpl_GP1_Slep, leftLegForces_GP1
        )
        backLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self.stackCblvpb_GP1_Slep, backLegForces_GP1
        )

        # The unsteady moment is calculated at the collocation point because the
        # unsteady force acts on the bound RingVortex, whose center is at the
        # collocation point, not at the Panel's centroid.
        unsteady_moments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self.stackCpp_GP1_Slep, unsteady_forces_GP1
        )

        self.moments_GP1_Slep = (
            rightLegMoments_GP1_Slep
            + frontLegMoments_GP1_Slep
            + leftLegMoments_GP1_Slep
            + backLegMoments_GP1_Slep
            + unsteady_moments_GP1_Slep
        )

        return moments_GP1_CgP1

    def _update_bound_vortex_positions_relative_to_slep_points(self) -> None:
        """Transform bound RingVortex leg center positions from CG-relative to SLEP-
        relative.

        For each panel, this method: 1. Gets the front-left panel point (leading edge)
        from each panel 2. Maps panels to their corresponding strip's leading edge point
        using slep_point_indices 3. Subtracts the SLEP position from all vortex leg
        center positions 4. Subtracts the SLEP position from collocation points

        This prepares positions for computing moments about the strip leading edge,
        which is important for analyzing local wing loading and deformations.

        :return: None
        """
        for panel_num, panel in enumerate(self.panels):
            self.stackFlpp_GP1_CgP1[panel_num] = panel.Flpp_GP1_CgP1
        slep_points = self.stackFlpp_GP1_CgP1[self.slep_point_indices]
        slep_map = (
            np.searchsorted(
                self.slep_point_indices, np.arange(self.num_panels), side="right"
            )
            - 1
        )
        self.stack_leading_edge_points = np.array([slep_points[i] for i in slep_map])
        self.stackCblvpr_GP1_Slep = (
            self.stackCblvpr_GP1_CgP1 - self.stack_leading_edge_points
        )
        self.stackCblvpf_GP1_Slep = (
            self.stackCblvpf_GP1_CgP1 - self.stack_leading_edge_points
        )
        self.stackCblvpl_GP1_Slep = (
            self.stackCblvpl_GP1_CgP1 - self.stack_leading_edge_points
        )
        self.stackCblvpb_GP1_Slep = (
            self.stackCblvpb_GP1_CgP1 - self.stack_leading_edge_points
        )

        # Find the collocation point positions relative to the SLEP points.
        self.stackCpp_GP1_Slep = self.stackCpp_GP1_CgP1 - self.stack_leading_edge_points
