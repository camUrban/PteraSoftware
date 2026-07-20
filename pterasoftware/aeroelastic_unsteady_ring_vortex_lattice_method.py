"""Contains the AeroelasticUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

AeroelasticUnsteadyRingVortexLatticeMethodSolver: A subclass of
CoupledUnsteadyRingVortexLatticeMethodSolver that solves AeroelasticUnsteadyProblems,
extending the coupled solver with Strip Leading Edge Point (SLEP) functionality for
computing the aerodynamic moments (in the first Airplane's geometry axes, relative to
the strip leading edge points) so that wing deformations can be coupled with aerodynamic
loads.

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

    This solver extends the coupled solver with calculations of the moments (in the
    first Airplane's geometry axes) about each panel's strip leading edge point, which
    is important for analyzing wing loading and deformation characteristics relative to
    the wing root.

    **Key additions over parent CoupledUnsteadyRingVortexLatticeMethodSolver:**
    initializes and maintains SLEP index mapping and position arrays, overrides
    ``_reinitialize_step_arrays_hook`` to reset SLEP arrays each step, overrides
    ``_load_calculation_moment_processing_hook`` to compute the moments (in the first
    Airplane's geometry axes) about the strip leading edge points, and computes the
    bound vortex positions (in the first Airplane's geometry axes) relative to the strip
    leading edge points.

    **Structural coupling output:** moments_GP1_Slep is the solver's one public SLEP
    attribute: a (num_panels, 3) ndarray of floats representing the moments (in the
    first Airplane's geometry axes, each relative to its panel's strip leading edge
    point) on every Panel at the current time step, which the
    AeroelasticUnsteadyProblem's structural solve reads as its aerodynamic forcing. The
    SLEP index mapping and the relative-position arrays that produce it are private
    working state.
    """

    __slots__ = (
        "_slep_point_indices",
        "_slep_outboard_is_left",
        "_stackCblvpr_GP1_Slep",
        "_stackCblvpf_GP1_Slep",
        "_stackCblvpl_GP1_Slep",
        "_stackCblvpb_GP1_Slep",
        "_stackCpp_GP1_Slep",
        "moments_GP1_Slep",
        "_stackSlep_GP1_CgP1",
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

        # Initialize SLEP (strip leading edge point) information. The SLEP is the
        # leading edge point of the strip's outboard bounding WingCrossSection, the same
        # WingCrossSection that receives the strip's torsional deformation, so each
        # strip's moments and its deformation share one reference. The solver's flat
        # Panel stack is ordered chord-major within each Wing (all spanwise positions of
        # the first chordwise row, then all spanwise positions of the second chordwise
        # row, and so on), with the Wings stacked in order. For each Panel, we record
        # the flat index of the Panel at the same Wing and spanwise position in the
        # first chordwise row, because that Panel's outboard front point is the strip's
        # leading edge point. Which grid corner is outboard depends on the meshing
        # direction: the front-right point on a root-to-tip grid, and the front-left
        # point on a mirror-meshed (tip-to-root) grid.
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
        self._slep_point_indices: np.ndarray = np.array(
            slep_point_indices_list, dtype=int
        )
        self._slep_outboard_is_left: np.ndarray = np.array(
            slep_outboard_is_left_list, dtype=bool
        )

        # These are the current time step's center bound line vortex points for the
        # right, front, left, and back legs (in the first Airplane's geometry axes,
        # relative to the local strip leading edge point).
        self._stackCblvpr_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self._stackCblvpf_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self._stackCblvpl_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self._stackCblvpb_GP1_Slep: np.ndarray = np.empty(0, dtype=float)

        # These are the collocation Panel points (in the first Airplane's geometry axes,
        # relative to the local strip leading edge point), and each Panel's own strip's
        # leading edge point (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG).
        self._stackCpp_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.moments_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self._stackSlep_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

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
        self._stackCblvpr_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self._stackCblvpf_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self._stackCblvpl_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self._stackCblvpb_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self._stackCpp_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self.moments_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self._stackSlep_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

    def _load_calculation_moment_processing_hook(
        self,
        rightLegForces_GP1: np.ndarray,
        frontLegForces_GP1: np.ndarray,
        leftLegForces_GP1: np.ndarray,
        backLegForces_GP1: np.ndarray,
        unsteady_forces_GP1: np.ndarray,
    ) -> np.ndarray:
        """Override parent to compute moments about both the first Airplane's CG and the
        strip leading edge points.

        This hook extends the parent class's moment calculation by additionally
        computing the moments (in the first Airplane's geometry axes) about each panel's
        strip leading edge point (SLEP). This is used for analyzing wing loading and
        deformation characteristics relative to the wing root.

        The method first calls the parent's implementation to get the moments (in the
        first Airplane's geometry axes, relative to the first Airplane's CG), then
        updates the bound vortex positions to be relative to the strip leading edge
        points, recalculates all moment contributions (in the first Airplane's geometry
        axes, each relative to its panel's strip leading edge point), and stores them in
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
            self._stackCblvpr_GP1_Slep, rightLegForces_GP1
        )
        frontLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self._stackCblvpf_GP1_Slep, frontLegForces_GP1
        )
        leftLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self._stackCblvpl_GP1_Slep, leftLegForces_GP1
        )
        backLegMoments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self._stackCblvpb_GP1_Slep, backLegForces_GP1
        )

        # The unsteady moment is calculated at the collocation point because the
        # unsteady force acts on the bound ring vortex, whose center is at the
        # collocation point, not at the Panel's centroid.
        unsteady_moments_GP1_Slep = _functions.numba_1d_explicit_cross(
            self._stackCpp_GP1_Slep, unsteady_forces_GP1
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
        """Transform the bound RingVortex leg center positions (in the first Airplane's
        geometry axes) from relative to the first Airplane's CG to relative to each
        panel's strip leading edge point.

        Gathers the outboard front point from each panel (the front-right point on a
        root-to-tip grid, and the front-left point on a mirror-meshed grid), maps each
        panel to its strip's leading edge point using _slep_point_indices, and subtracts
        that SLEP position from the vortex leg center positions and the collocation
        points.

        This prepares positions for computing the moments (in the first Airplane's
        geometry axes) about the strip leading edge points, which is important for
        analyzing local wing loading and deformations.

        :return: None
        """
        # Stage each Panel's outboard front point. Only the first-chordwise-row entries
        # are consumed by the gather below, since every strip's SLEP is a first-row
        # Panel's corner.
        outboardFrontPoints_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
        for panel_num, panel in enumerate(self.panels):
            if self._slep_outboard_is_left[panel_num]:
                outboardFrontPoints_GP1_CgP1[panel_num] = panel.Flpp_GP1_CgP1
            else:
                outboardFrontPoints_GP1_CgP1[panel_num] = panel.Frpp_GP1_CgP1
        self._stackSlep_GP1_CgP1 = outboardFrontPoints_GP1_CgP1[
            self._slep_point_indices
        ]
        self._stackCblvpr_GP1_Slep = (
            self.stackCblvpr_GP1_CgP1 - self._stackSlep_GP1_CgP1
        )
        self._stackCblvpf_GP1_Slep = (
            self.stackCblvpf_GP1_CgP1 - self._stackSlep_GP1_CgP1
        )
        self._stackCblvpl_GP1_Slep = (
            self.stackCblvpl_GP1_CgP1 - self._stackSlep_GP1_CgP1
        )
        self._stackCblvpb_GP1_Slep = (
            self.stackCblvpb_GP1_CgP1 - self._stackSlep_GP1_CgP1
        )

        # Find the collocation point positions relative to the SLEP points.
        self._stackCpp_GP1_Slep = self.stackCpp_GP1_CgP1 - self._stackSlep_GP1_CgP1
