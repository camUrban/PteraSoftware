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
    initializes and maintains the SLEP index mapping and position arrays, overrides
    ``_reinitialize_step_arrays_hook`` to reset the SLEP arrays each step, and overrides
    ``_process_panel_loads_hook`` to compute the moments (in the first Airplane's
    geometry axes) about the strip leading edge points from the per Panel total loads.

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

        # These are the moments (in the first Airplane's geometry axes, each relative to
        # its panel's strip leading edge point) on every Panel, and each Panel's own
        # strip's leading edge point (in the first Airplane's geometry axes, relative to
        # the first Airplane's CG).
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
        self.moments_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
        self._stackSlep_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

    def _process_panel_loads_hook(
        self,
        forces_GP1: np.ndarray,
        moments_GP1_CgP1: np.ndarray,
    ) -> None:
        """Override parent to additionally compute the moments about the strip leading
        edge points.

        This hook extends the parent class's load processing by computing the moments
        (in the first Airplane's geometry axes, each relative to its panel's strip
        leading edge point) on every Panel and storing them in self.moments_GP1_Slep,
        which the AeroelasticUnsteadyProblem's structural solve reads as its aerodynamic
        forcing.

        Each Panel's moment about its strip leading edge point follows from the moment
        transport theorem: it is the Panel's moment about the first Airplane's CG minus
        the cross product of the strip leading edge point's position (relative to that
        CG) with the Panel's total force. Because this consumes only the per Panel total
        loads, it holds for every force method regardless of where each method applies
        its component forces.

        :param forces_GP1: A (N,3) ndarray of floats representing the forces (in the
            first Airplane's geometry axes) on every Panel at the current time step. The
            units are in Newtons.
        :param moments_GP1_CgP1: A (N,3) ndarray of floats representing the moments (in
            the first Airplane's geometry axes, relative to the first Airplane's CG) on
            every Panel at the current time step. The units are in Newton-meters.
        :return: None
        """
        self._populate_slep_positions()

        self.moments_GP1_Slep = moments_GP1_CgP1 - _functions.numba_1d_explicit_cross(
            self._stackSlep_GP1_CgP1, forces_GP1
        )

        super()._process_panel_loads_hook(forces_GP1, moments_GP1_CgP1)

    def _populate_slep_positions(self) -> None:
        """Populate each Panel's strip leading edge point position (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        Gathers the outboard front point from each panel (the front-right point on a
        root-to-tip grid, and the front-left point on a mirror-meshed grid) and maps
        each panel to its strip's leading edge point using _slep_point_indices.

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
