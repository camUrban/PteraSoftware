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
        "_slep_outboard_is_left",
        "stackCblvpr_GP1_Slep",
        "stackCblvpf_GP1_Slep",
        "stackCblvpl_GP1_Slep",
        "stackCblvpb_GP1_Slep",
        "stackCpp_GP1_Slep",
        "moments_GP1_Slep",
        "stackSlep_GP1_CgP1",
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

        # Initialize SLEP (Strip Leading Edge Point) information. The SLEP is the
        # leading edge point of the strip's outboard bounding WingCrossSection, the
        # same WingCrossSection that receives the strip's twist, so each strip's
        # moments and its deformation share one reference. The solver's flat panel
        # stack is ordered chord-major within each wing (all spanwise positions of
        # the first chordwise row, then all spanwise positions of the second
        # chordwise row, and so on), with the wings stacked in order. For each
        # panel, we record the flat index of the panel at the same wing and
        # spanwise position in the first chordwise row, because that panel's
        # outboard front point is the strip's leading edge point. Which grid corner
        # is outboard depends on the meshing direction: the front-right point on a
        # root-to-tip grid, and the front-left point on a mirror-meshed
        # (tip-to-root) grid.
        panel_count = 0
        slep_point_indices_list: list[int] = []
        slep_outboard_is_left_list: list[bool] = []
        for airplane in first_steady_problem.airplanes:
            for wing in airplane.wings:
                num_spanwise_panels = wing.num_spanwise_panels
                assert num_spanwise_panels is not None
                for _ in range(wing.num_chordwise_panels):
                    for spanwise_position in range(num_spanwise_panels):
                        slep_point_indices_list.append(panel_count + spanwise_position)
                num_wing_panels = wing.num_chordwise_panels * num_spanwise_panels
                slep_outboard_is_left_list.extend([wing.mirror_only] * num_wing_panels)
                panel_count += num_wing_panels
        self.slep_point_indices: np.ndarray = np.array(
            slep_point_indices_list, dtype=int
        )
        self._slep_outboard_is_left: np.ndarray = np.array(
            slep_outboard_is_left_list, dtype=bool
        )

        # The current time step's center bound LineVortex points for the right,
        # front, left, and back legs (in the first Airplane's geometry axes,
        # relative to the local strip leading edge point).
        self.stackCblvpr_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpf_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpl_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpb_GP1_Slep: np.ndarray = np.empty(0, dtype=float)

        # The colocation panel points (in the first Airplane's geometry axes,
        # relative to the local strip leading edge point), and each panel's own
        # strip's leading edge point (in the first Airplane's geometry axes,
        # relative to the first Airplane's CG).
        self.stackCpp_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.moments_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackSlep_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

    @property
    def _aeroelastic_unsteady_problem(self) -> problems.AeroelasticUnsteadyProblem:
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
        self.stackSlep_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

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

        Gathers the outboard front point from each panel (the front-right point on a
        root-to-tip grid, and the front-left point on a mirror-meshed grid), maps each
        panel to its strip's leading edge point using slep_point_indices, and subtracts
        that SLEP position from the vortex leg center positions and the collocation
        points.

        This prepares positions for computing moments about the strip leading edge,
        which is important for analyzing local wing loading and deformations.

        :return: None
        """
        # Stage each panel's outboard front point. Only the first-chordwise-row
        # entries are consumed by the gather below, since every strip's SLEP is a
        # first-row panel's corner.
        outboardFrontPoints_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        for panel_num, panel in enumerate(self.panels):
            if self._slep_outboard_is_left[panel_num]:
                outboardFrontPoints_GP1_CgP1[panel_num] = panel.Flpp_GP1_CgP1
            else:
                outboardFrontPoints_GP1_CgP1[panel_num] = panel.Frpp_GP1_CgP1
        self.stackSlep_GP1_CgP1 = outboardFrontPoints_GP1_CgP1[self.slep_point_indices]
        self.stackCblvpr_GP1_Slep = self.stackCblvpr_GP1_CgP1 - self.stackSlep_GP1_CgP1
        self.stackCblvpf_GP1_Slep = self.stackCblvpf_GP1_CgP1 - self.stackSlep_GP1_CgP1
        self.stackCblvpl_GP1_Slep = self.stackCblvpl_GP1_CgP1 - self.stackSlep_GP1_CgP1
        self.stackCblvpb_GP1_Slep = self.stackCblvpb_GP1_CgP1 - self.stackSlep_GP1_CgP1

        # Find the collocation point positions relative to the SLEP points.
        self.stackCpp_GP1_Slep = self.stackCpp_GP1_CgP1 - self.stackSlep_GP1_CgP1
