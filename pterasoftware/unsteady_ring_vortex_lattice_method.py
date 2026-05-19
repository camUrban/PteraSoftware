"""Contains the UnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

UnsteadyRingVortexLatticeMethodSolver: A class used to solve UnsteadyProblems with the
unsteady ring vortex lattice method.

**Contains the following functions:**

None
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import cast

import numpy as np
from tqdm import tqdm

from . import (
    _aerodynamics_functions,
    _core,
    _functions,
    _logging,
    _panel,
    _parameter_validation,
    _transformations,
    geometry,
    operating_point,
    problems,
)

_logger = _logging.get_logger("unsteady_ring_vortex_lattice_method")


# REFACTOR: Add unit tests for trapezoid-rule-based averages for the mean and RMS loads
#  and load coefficients.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
class UnsteadyRingVortexLatticeMethodSolver:
    """A class used to solve UnsteadyProblems with the unsteady ring vortex lattice
    method.

    **Contains the following methods:**

    run: Runs the solver on the UnsteadyProblem.

    initialize_step_geometry: Initializes geometry for a specific step without solving.

    calculate_solution_velocity: Finds the fluid velocity (in the first Airplane's
    geometry axes, observed from the Earth frame) at one or more points (in the first
    Airplane's geometry axes, relative to the first Airplane's CG) due to the freestream
    velocity and the induced velocity from every ring vortex.
    """

    __slots__ = (
        "unsteady_problem",
        "_max_wake_rows",
        "num_steps",
        "delta_time",
        "first_results_step",
        "_first_averaging_step",
        "_current_step",
        "_prescribed_wake",
        "steady_problems",
        "current_airplanes",
        "current_operating_point",
        "num_airplanes",
        "num_panels",
        "_currentVInf_GP1__E",
        "_currentStackFreestreamWingInfluences__E",
        "_currentGridWingWingInfluences__E",
        "_currentStackWakeWingInfluences__E",
        "_current_bound_vortex_strengths",
        "_last_bound_vortex_strengths",
        "panels",
        "stackUnitNormals_GP1",
        "panel_areas",
        "stackCpp_GP1_CgP1",
        "_stackLastCpp_GP1_CgP1",
        "stackBrbrvp_GP1_CgP1",
        "stackFrbrvp_GP1_CgP1",
        "stackFlbrvp_GP1_CgP1",
        "stackBlbrvp_GP1_CgP1",
        "_lastStackBrbrvp_GP1_CgP1",
        "_lastStackFrbrvp_GP1_CgP1",
        "_lastStackFlbrvp_GP1_CgP1",
        "_lastStackBlbrvp_GP1_CgP1",
        "_listStackBrbrvp_GP1_CgP1",
        "_listStackFrbrvp_GP1_CgP1",
        "_listStackFlbrvp_GP1_CgP1",
        "_listStackBlbrvp_GP1_CgP1",
        "_per_wing_panel_offsets",
        "_per_wing_num_chordwise_panels",
        "_per_wing_num_spanwise_panels",
        "_per_wing_spanwise_cumsum",
        "stackCblvpr_GP1_CgP1",
        "stackCblvpf_GP1_CgP1",
        "stackCblvpl_GP1_CgP1",
        "stackCblvpb_GP1_CgP1",
        "_lastStackCblvpr_GP1_CgP1",
        "_lastStackCblvpf_GP1_CgP1",
        "_lastStackCblvpl_GP1_CgP1",
        "_lastStackCblvpb_GP1_CgP1",
        "stackRbrv_GP1",
        "stackFbrv_GP1",
        "stackLbrv_GP1",
        "stackBbrv_GP1",
        "panel_is_trailing_edge",
        "panel_is_leading_edge",
        "panel_is_left_edge",
        "panel_is_right_edge",
        "_current_wake_vortex_strengths",
        "_current_wake_vortex_ages",
        "_currentStackBrwrvp_GP1_CgP1",
        "_currentStackFrwrvp_GP1_CgP1",
        "_currentStackFlwrvp_GP1_CgP1",
        "_currentStackBlwrvp_GP1_CgP1",
        "list_num_wake_vortices",
        "_list_wake_vortex_strengths",
        "listStackBrwrvp_GP1_CgP1",
        "listStackFrwrvp_GP1_CgP1",
        "listStackFlwrvp_GP1_CgP1",
        "listStackBlwrvp_GP1_CgP1",
        "_currentStackBoundRc0s",
        "_currentStackWakeRc0s",
        "stackSeedPoints_GP1_CgP1",
        "gridStreamlinePoints_GP1_CgP1",
        "ran",
    )

    def __init__(self, unsteady_problem: _core.CoreUnsteadyProblem) -> None:
        """The initialization method.

        :param unsteady_problem: The UnsteadyProblem to be solved.
        :return: None
        """
        # Guard direct instantiation of the base solver against coupled problems while
        # allowing subclasses to pass their own CoreUnsteadyProblem variants via super().
        if type(self) is UnsteadyRingVortexLatticeMethodSolver and not isinstance(
            unsteady_problem, problems.UnsteadyProblem
        ):
            raise TypeError("unsteady_problem must be an UnsteadyProblem.")
        self.unsteady_problem = unsteady_problem

        self._max_wake_rows = self.unsteady_problem.max_wake_rows
        self.num_steps = self.unsteady_problem.num_steps
        self.delta_time = self.unsteady_problem.delta_time
        self.first_results_step = self.unsteady_problem.first_results_step
        self._first_averaging_step = self.unsteady_problem.first_averaging_step
        self._current_step: int = 0
        self._prescribed_wake: bool = True

        self.steady_problems = self.unsteady_problem.steady_problems

        first_steady_problem: problems.SteadyProblem = self._get_steady_problem_at(0)

        self.current_airplanes: tuple[geometry.airplane.Airplane, ...] = ()
        self.current_operating_point: operating_point.OperatingPoint = (
            first_steady_problem.operating_point
        )
        self.num_airplanes: int = len(first_steady_problem.airplanes)

        num_panels = 0
        for airplane in first_steady_problem.airplanes:
            num_panels += airplane.num_panels
        self.num_panels: int = num_panels

        # Initialize attributes to hold aerodynamic data that pertain to the simulation.
        self._currentVInf_GP1__E: np.ndarray = (
            first_steady_problem.operating_point.vInf_GP1__E
        )
        self._currentStackFreestreamWingInfluences__E: np.ndarray = np.empty(
            0, dtype=float
        )
        self._currentGridWingWingInfluences__E: np.ndarray = np.empty(0, dtype=float)
        self._currentStackWakeWingInfluences__E: np.ndarray = np.empty(0, dtype=float)
        # Initialized to ones so that initialize_step_geometry (which can run
        # before any strength solve) finds a sensible placeholder when shedding
        # the next step's wake. run() overwrites this each step before solving.
        self._current_bound_vortex_strengths: np.ndarray = np.ones(
            self.num_panels, dtype=float
        )
        # _last_bound_vortex_strengths starts as zeros so step 0 sees no previous
        # step contribution. At each step's _calculate_vortex_strengths the just
        # solved current strengths are captured here for use by the next step.
        self._last_bound_vortex_strengths: np.ndarray = np.zeros(
            self.num_panels, dtype=float
        )

        # Initialize attributes to hold geometric data that pertain to this
        # UnsteadyProblem.
        self.panels: np.ndarray = np.empty(0, dtype=object)
        self.stackUnitNormals_GP1: np.ndarray = np.empty(0, dtype=float)
        self.panel_areas: np.ndarray = np.empty(0, dtype=float)

        # The current and last time step's collocation panel points (in the first
        # Airplane's geometry axes, relative to the first Airplane's CG).
        self.stackCpp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._stackLastCpp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

        # The current and last time step's back right, front right, front left,
        # and back left bound ring vortex points (in the first Airplane's geometry
        # axes, relative to the first Airplane's CG).
        self.stackBrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackFrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackFlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackBlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackBrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackFrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackFlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackBlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

        # Per step bound ring vortex corner stacks. Filled by
        # _initialize_panel_vortices_at and read by _collapse_geometry,
        # _populate_next_airplanes_wake_vortex_points, and the last step block
        # of the next step's _collapse_geometry. Pre-allocated in run().
        self._listStackBrbrvp_GP1_CgP1: list[np.ndarray] = []
        self._listStackFrbrvp_GP1_CgP1: list[np.ndarray] = []
        self._listStackFlbrvp_GP1_CgP1: list[np.ndarray] = []
        self._listStackBlbrvp_GP1_CgP1: list[np.ndarray] = []

        # Per (airplane, wing) flat panel offsets and shapes. Built immediately
        # below so the bound vortex pipeline and the wake convection loop can
        # index into the per step list arrays without re-deriving offsets.
        self._per_wing_panel_offsets: list[list[int]] = []
        self._per_wing_num_chordwise_panels: list[list[int]] = []
        self._per_wing_num_spanwise_panels: list[list[int]] = []
        self._per_wing_spanwise_cumsum: list[list[int]] = []
        running_offset = 0
        cumulative_spanwise = 0
        for airplane in first_steady_problem.airplanes:
            airplane_offsets: list[int] = []
            airplane_chordwise: list[int] = []
            airplane_spanwise: list[int] = []
            airplane_cumsum: list[int] = []
            for wing in airplane.wings:
                _wing_num_spanwise_panels = wing.num_spanwise_panels
                assert _wing_num_spanwise_panels is not None

                airplane_offsets.append(running_offset)
                airplane_chordwise.append(wing.num_chordwise_panels)
                airplane_spanwise.append(_wing_num_spanwise_panels)
                airplane_cumsum.append(cumulative_spanwise)
                running_offset += wing.num_chordwise_panels * _wing_num_spanwise_panels
                cumulative_spanwise += _wing_num_spanwise_panels
            self._per_wing_panel_offsets.append(airplane_offsets)
            self._per_wing_num_chordwise_panels.append(airplane_chordwise)
            self._per_wing_num_spanwise_panels.append(airplane_spanwise)
            self._per_wing_spanwise_cumsum.append(airplane_cumsum)
        total_spanwise_panels = cumulative_spanwise

        # Pre-allocate the per step bound ring vortex corner stacks. These are
        # filled by _initialize_panel_vortices_at(step) and read by the bound
        # vortex pipeline thereafter.
        self._listStackBrbrvp_GP1_CgP1 = [
            np.zeros((self.num_panels, 3), dtype=float) for _ in range(self.num_steps)
        ]
        self._listStackFrbrvp_GP1_CgP1 = [
            np.zeros((self.num_panels, 3), dtype=float) for _ in range(self.num_steps)
        ]
        self._listStackFlbrvp_GP1_CgP1 = [
            np.zeros((self.num_panels, 3), dtype=float) for _ in range(self.num_steps)
        ]
        self._listStackBlbrvp_GP1_CgP1 = [
            np.zeros((self.num_panels, 3), dtype=float) for _ in range(self.num_steps)
        ]

        # Pre-allocate the per step wake stacks. The number of wake ring vortices
        # at step S is min(S, max_wake_rows) (or S if no truncation) times the
        # total spanwise panel count summed over all wings.
        wake_sizes_per_step = []
        for step in range(self.num_steps):
            num_chordwise_wake_rows = step
            if self._max_wake_rows is not None:
                num_chordwise_wake_rows = min(step, self._max_wake_rows)
            wake_sizes_per_step.append(num_chordwise_wake_rows * total_spanwise_panels)

        self.list_num_wake_vortices: list[int] = list(wake_sizes_per_step)
        self._list_wake_vortex_strengths: list[np.ndarray] = [
            np.zeros(n, dtype=float) for n in wake_sizes_per_step
        ]
        self.listStackBrwrvp_GP1_CgP1: list[np.ndarray] = [
            np.zeros((n, 3), dtype=float) for n in wake_sizes_per_step
        ]
        self.listStackFrwrvp_GP1_CgP1: list[np.ndarray] = [
            np.zeros((n, 3), dtype=float) for n in wake_sizes_per_step
        ]
        self.listStackFlwrvp_GP1_CgP1: list[np.ndarray] = [
            np.zeros((n, 3), dtype=float) for n in wake_sizes_per_step
        ]
        self.listStackBlwrvp_GP1_CgP1: list[np.ndarray] = [
            np.zeros((n, 3), dtype=float) for n in wake_sizes_per_step
        ]

        # The current and last time step's center bound line vortex points for the
        # right, front, left, and back legs (in the first Airplane's geometry axes,
        # relative to the first Airplane's CG).
        self.stackCblvpr_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpf_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpl_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpb_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackCblvpr_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackCblvpf_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackCblvpl_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackCblvpb_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

        # Right, front, left, and back bound ring vortex vectors (in the first
        # Airplane's geometry axes).
        self.stackRbrv_GP1: np.ndarray = np.empty(0, dtype=float)
        self.stackFbrv_GP1: np.ndarray = np.empty(0, dtype=float)
        self.stackLbrv_GP1: np.ndarray = np.empty(0, dtype=float)
        self.stackBbrv_GP1: np.ndarray = np.empty(0, dtype=float)

        # Initialize variables to hold aerodynamic data that pertains details about
        # each Panel's location on its Wing.
        self.panel_is_trailing_edge: np.ndarray = np.empty(0, dtype=bool)
        self.panel_is_leading_edge: np.ndarray = np.empty(0, dtype=bool)
        self.panel_is_left_edge: np.ndarray = np.empty(0, dtype=bool)
        self.panel_is_right_edge: np.ndarray = np.empty(0, dtype=bool)

        # Initialize variables to hold aerodynamic data that pertains to the wake at
        # the current time step.
        self._current_wake_vortex_strengths: np.ndarray = np.empty(0, dtype=float)
        self._current_wake_vortex_ages: np.ndarray = np.empty(0, dtype=float)

        # The current time step's back right, front right, front left, and back left
        # wake ring vortex points (in the first Airplane's geometry axes, relative to
        # the first Airplane's CG).
        self._currentStackBrwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._currentStackFrwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._currentStackFlwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._currentStackBlwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

        # The list attributes above (list_num_wake_vortices,
        # _list_wake_vortex_strengths, listStack{Br,Fr,Fl,Bl}wrvp_GP1_CgP1)
        # were pre-allocated above this block.

        self._currentStackBoundRc0s: np.ndarray = np.empty(0, dtype=float)
        self._currentStackWakeRc0s: np.ndarray = np.empty(0, dtype=float)

        self.stackSeedPoints_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.gridStreamlinePoints_GP1_CgP1: np.ndarray = np.empty((0, 3), dtype=float)

        self.ran = False

    def run(
        self,
        prescribed_wake: bool | np.bool_ = True,
        calculate_streamlines: bool | np.bool_ = True,
        show_progress: bool | np.bool_ = True,
    ) -> None:
        """Runs the solver on the UnsteadyProblem.

        :param prescribed_wake: Set this to True to solve using a prescribed wake model.
            Set to False to use a free-wake, which may be more accurate but will make
            the fun method significantly slower. Can be a bool or a numpy bool and will
            be converted internally to a bool. The default is True.
        :param calculate_streamlines: Set this to True to calculate streamlines
            emanating from the back of the wing after running the solver. It can be a
            bool or a numpy bool and will be converted internally to a bool. The default
            is True.
        :param show_progress: Set this to True to show the TQDM progress bar. For
            showing the progress bar and displaying log statements, set up logging using
            the setup_logging function. It can be a bool or a numpy bool and will be
            converted internally to a bool. The default is True.
        :return: None
        """
        self._prescribed_wake = _parameter_validation.boolLike_return_bool(
            prescribed_wake, "prescribed_wake"
        )
        calculate_streamlines = _parameter_validation.boolLike_return_bool(
            calculate_streamlines, "calculate_streamlines"
        )
        show_progress = _parameter_validation.boolLike_return_bool(
            show_progress, "show_progress"
        )

        # The per step list arrays for both bound and wake state were
        # pre-allocated in __init__. Recompute the total panel count (used by
        # the progress bar weighting) here.
        num_wing_panels = self.num_panels

        # The following loop attempts to predict how much time each time step will
        # take, relative to the other time steps. This data will be used to generate
        # estimates of how much longer a simulation will take, and create a smoothly
        # advancing progress bar.

        # Initialize list that will hold the approximate, relative times. This has
        # one more element than the number of time steps, because I will also use the
        # progress bar during the simulation initialization.
        approx_times = np.zeros(self.num_steps + 1, dtype=float)
        for step in range(1, self.num_steps):
            # Calculate the total number of ring vortices analyzed during this step.
            num_wake_ring_vortices = self.list_num_wake_vortices[step]
            num_ring_vortices = num_wing_panels + num_wake_ring_vortices

            # The following constant multipliers were determined empirically. Thus
            # far, they seem to provide for adequately smooth progress bar updating.
            if step == 1:
                approx_times[step] = num_ring_vortices * 70
            elif step == 2:
                approx_times[step] = num_ring_vortices * 30
            else:
                approx_times[step] = num_ring_vortices * 3

        approx_partial_time = np.sum(approx_times)
        approx_times[0] = round(approx_partial_time / 100)
        approx_total_time = np.sum(approx_times)

        with tqdm(
            total=approx_total_time,
            unit="",
            unit_scale=True,
            ncols=100,
            desc="Simulating",
            disable=not show_progress,
            bar_format="{desc}:{percentage:3.0f}% |{bar}| Elapsed: {elapsed}, "
            "Remaining: {remaining}",
        ) as bar:
            # Update the progress bar based on the initialization step's predicted
            # approximate, relative computing time.
            bar.update(n=float(approx_times[0]))

            # Iterate through the time steps.
            for step in range(self.num_steps):

                # Save attributes to hold the current step, Airplanes,
                # and OperatingPoint, and freestream velocity (in the first
                # Airplane's geometry axes, observed from the Earth frame).
                self._current_step = step

                # Initialize this step's bound ring vortices. The default does an
                # upfront init for all steps on step 0 and is a no-op thereafter;
                # coupled subclasses override this hook to init one step at a time.
                self._initialize_step_vortices(step)
                current_problem: problems.SteadyProblem = self._get_steady_problem_at(
                    self._current_step
                )
                self.current_airplanes = current_problem.airplanes
                self.current_operating_point = current_problem.operating_point
                self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E
                _logger.debug(
                    "Beginning time step "
                    + str(self._current_step)
                    + " out of "
                    + str(self.num_steps - 1)
                    + "."
                )

                # TODO: I think these steps are redundant, at least during the first
                #  time step. Consider dropping them.
                # Initialize attributes to hold aerodynamic data that pertain to the
                # simulation at this time step.
                self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E
                self._currentStackFreestreamWingInfluences__E = np.zeros(
                    self.num_panels, dtype=float
                )
                self._currentGridWingWingInfluences__E = np.zeros(
                    (self.num_panels, self.num_panels), dtype=float
                )
                self._currentStackWakeWingInfluences__E = np.zeros(
                    self.num_panels, dtype=float
                )
                self._current_bound_vortex_strengths = np.ones(
                    self.num_panels, dtype=float
                )
                # _last_bound_vortex_strengths is left alone here. At step 0 it is
                # the zeros allocated in __init__. At step > 0 it holds the
                # previous step's solved strengths, captured at the end of the
                # previous step's _calculate_vortex_strengths.

                # Initialize attributes to hold geometric data that pertain to this
                # UnsteadyProblem.
                self.panels = np.empty(self.num_panels, dtype=object)
                self.stackUnitNormals_GP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.panel_areas = np.zeros(self.num_panels, dtype=float)

                self.stackCpp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self._stackLastCpp_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )

                self.stackBrbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFrbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFlbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackBlbrvp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackBrbrvp_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )
                self._lastStackFrbrvp_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )
                self._lastStackFlbrvp_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )
                self._lastStackBlbrvp_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )

                self.stackCblvpr_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpf_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpl_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpb_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackCblvpr_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )
                self._lastStackCblvpf_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )
                self._lastStackCblvpl_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )
                self._lastStackCblvpb_GP1_CgP1 = np.zeros(
                    (self.num_panels, 3), dtype=float
                )

                self.stackRbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackLbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)
                self.stackBbrv_GP1 = np.zeros((self.num_panels, 3), dtype=float)

                # Initialize variables to hold details about each Panel's location on
                # its Wing.
                self.panel_is_trailing_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_leading_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_left_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_right_edge = np.zeros(self.num_panels, dtype=bool)

                # Hook: subclasses may reinitialize step-specific arrays here.
                self._reinitialize_step_arrays_hook()

                # Get the pre-allocated (but still all zero) arrays of wake
                # information that are associated with this time step.
                self._current_wake_vortex_strengths = self._list_wake_vortex_strengths[
                    step
                ]
                self._currentStackBrwrvp_GP1_CgP1 = self.listStackBrwrvp_GP1_CgP1[step]
                self._currentStackFrwrvp_GP1_CgP1 = self.listStackFrwrvp_GP1_CgP1[step]
                self._currentStackFlwrvp_GP1_CgP1 = self.listStackFlwrvp_GP1_CgP1[step]
                self._currentStackBlwrvp_GP1_CgP1 = self.listStackBlwrvp_GP1_CgP1[step]

                self._currentStackBoundRc0s = np.zeros(self.num_panels, dtype=float)
                num_wake_vortices = self.list_num_wake_vortices[step]
                self._current_wake_vortex_ages = np.zeros(
                    num_wake_vortices, dtype=float
                )
                self._currentStackWakeRc0s = np.zeros(num_wake_vortices, dtype=float)

                self.stackSeedPoints_GP1_CgP1 = np.zeros((0, 3), dtype=float)

                # Collapse the geometry matrices into 1D ndarrays of attributes.
                _logger.debug("Collapsing the geometry.")
                self._collapse_geometry()

                # Find the matrix of Wing Wing influence coefficients associated with
                # the Airplanes' geometries at this time step.
                _logger.debug("Calculating the Wing Wing influences.")
                self._calculate_wing_wing_influences()

                # Find the normal velocity (in the first Airplane's geometry axes,
                # observed from the Earth frame) at every collocation point due
                # solely to the freestream.
                _logger.debug("Calculating the freestream Wing influences.")
                self._calculate_freestream_wing_influences()

                # Find the normal velocity (in the first Airplane's geometry axes,
                # observed from the Earth frame) at every collocation point due
                # solely to the wake ring vortices.
                _logger.debug("Calculating the wake Wing influences.")
                self._calculate_wake_wing_influences()

                # Solve for each bound ring vortex's strength.
                _logger.debug("Calculating bound ring vortex strengths.")
                self._calculate_vortex_strengths()

                # Solve for the forces (in the first Airplane's geometry axes) and
                # moments (in the first Airplane's geometry axes, relative to the
                # first Airplane's CG) on each Panel.
                if self._current_step >= self.first_results_step:
                    _logger.debug("Calculating forces and moments.")
                    self._calculate_loads()

                # Hook: subclasses may inject work between load calculation and wake
                # shedding (e.g. coupled problems update the next step's geometry
                # from this step's solver results).
                self._pre_shed_hook(step)

                # Shed ring vortices into the wake.
                _logger.debug("Shedding ring vortices into the wake.")
                self._populate_next_airplanes_wake()

                # Snapshot this step's solved bound ring vortex strengths so the
                # next step can use them as the "last" strengths in its
                # _calculate_loads (for the unsteady force term and the back leg
                # effective strength). The copy is needed because the next step
                # will rebind _current_bound_vortex_strengths to a new array.
                self._last_bound_vortex_strengths = (
                    self._current_bound_vortex_strengths.copy()
                )

                # Update the progress bar based on this time step's predicted
                # approximate, relative computing time.
                bar.update(n=float(approx_times[step + 1]))

            _logger.debug("Calculating averaged or final forces and moments.")
            self._finalize_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges, if requested.
        if calculate_streamlines:
            _logger.debug("Calculating streamlines.")
            _functions.calculate_streamlines(self)

        # Mark that the solver has run.
        self.ran = True

    def initialize_step_geometry(self, step: int) -> None:
        """Initializes geometry for a specific step without solving.

        Sets up bound ring vortices and wake ring vortices for the specified time step,
        but does not solve the aerodynamic system. Use this for geometry only analysis
        like delta_time optimization.

        This method must be called sequentially for each step starting from 0, as wake
        vortices at step N depend on the geometry from step N - 1.

        :param step: The time step to initialize geometry for. It is zero indexed. It
            must be a non negative int and be less than the total number of steps.
        :return: None
        """
        step = _parameter_validation.int_in_range_return_int(
            step, "step", 0, True, self.num_steps, False
        )

        # Initialize bound ring vortices. The base solver's hook does an upfront init
        # for all steps when step is 0 and is a no-op otherwise; coupled subclasses
        # override to initialize only the specified step.
        self._initialize_step_vortices(step)

        # Set the current step and related state.
        self._current_step = step
        current_problem: problems.SteadyProblem = self._get_steady_problem_at(step)
        self.current_airplanes = current_problem.airplanes
        self.current_operating_point = current_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E

        # Populate the wake for the next step (if not the last step).
        if step < self.num_steps - 1:
            self._populate_next_airplanes_wake_vortex_points()
            self._populate_next_airplanes_wake_vortices()

    def _initialize_step_vortices(self, step: int) -> None:
        """Initializes this time step's bound ring vortices.

        The default implementation initializes bound ring vortices for all time steps
        upfront on step 0 and is a no-op on subsequent steps. Coupled subclasses
        override this to initialize only the given step, since their geometry is
        determined dynamically from the solver's results at the previous step.

        :param step: The time step to initialize.
        :return: None
        """
        if step == 0:
            _logger.debug("Initializing all Airplanes' bound ring vortices.")
            self._initialize_panel_vortices()

    def _reinitialize_step_arrays_hook(self) -> None:
        """Hook for subclasses to reinitialize step specific arrays.

        Called once per time step in run(), after the standard per step arrays are
        reinitialized and before the wake arrays are retrieved. The default
        implementation is a no op. Subclasses may override this to zero out or
        reallocate feature specific arrays at the start of each step.

        :return: None
        """

    def _pre_shed_hook(self, step: int) -> None:
        """Hook for subclasses to inject work between load calculation and wake shed.

        Called once per time step in run(), after this step's loads have been calculated
        and before wake ring vortices are shed. The default implementation is a no op.
        Coupled subclasses override this to update the next time step's geometry from
        the current step's solver results.

        :param step: The current time step.
        :return: None
        """

    def _initialize_panel_vortices(self) -> None:
        """Calculates the locations of the bound ring vortex vertices for all time
        steps, and then initializes them.

        Every Panel has a ring vortex, which is a quadrangle whose front leg is a line
        vortex at the Panel's quarter chord. The left and right legs are line vortices
        running along the Panel's left and right legs. If the Panel is not along the
        trailing edge, they extend backwards and meet the back line vortex, at the rear
        Panel's quarter chord. Otherwise, they extend backwards and meet the back line
        vortex one quarter chord back from the Panel's back leg.

        :return: None
        """
        for step in range(self.num_steps):
            self._initialize_panel_vortices_at(step)

    def _collapse_geometry(self) -> None:
        """Collapses the bound vortex and wake state into 1D ndarrays for the current
        time step.

        Bound ring vortex corner positions for this and the previous step are read from
        the per step list arrays populated by _initialize_panel_vortices_at. Wake ring
        vortex corner positions, strengths, and ages are already stored in the per step
        list arrays; this method only needs to wire the per Panel scalars and last step
        leg derivatives.

        :return: None
        """
        step = self._current_step
        currentStackFr = self._listStackFrbrvp_GP1_CgP1[step]
        currentStackFl = self._listStackFlbrvp_GP1_CgP1[step]
        currentStackBl = self._listStackBlbrvp_GP1_CgP1[step]
        currentStackBr = self._listStackBrbrvp_GP1_CgP1[step]

        # Iterate through the current time step's Airplanes' Wings.
        global_panel_position = 0
        global_wake_ring_vortex_position = 0
        for airplane in self.current_airplanes:
            for wing in airplane.wings:
                _standard_mean_chord = wing.standard_mean_chord
                assert _standard_mean_chord is not None
                wing_r_c0 = 0.03 * _standard_mean_chord

                _panels = wing.panels
                assert _panels is not None

                _num_spanwise_panels = wing.num_spanwise_panels
                assert _num_spanwise_panels is not None

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(_panels)

                # Iterate through the 1D ndarray of this Wing's Panels and write
                # the per Panel scalars plus the corner-derived bound leg arrays.
                panel: _panel.Panel
                for panel in panels:
                    _functions.update_ring_vortex_solvers_panel_attributes(
                        ring_vortex_solver=self,
                        global_panel_position=global_panel_position,
                        panel=panel,
                        Frrvp_GP1_CgP1=currentStackFr[global_panel_position],
                        Flrvp_GP1_CgP1=currentStackFl[global_panel_position],
                        Blrvp_GP1_CgP1=currentStackBl[global_panel_position],
                        Brrvp_GP1_CgP1=currentStackBr[global_panel_position],
                    )
                    self._currentStackBoundRc0s[global_panel_position] = wing_r_c0
                    global_panel_position += 1

                # Set the wake characteristic core radius for every wake
                # ring vortex contributed by this Wing at this step. The wake
                # corner positions and strengths are stored in the per step list arrays
                # aliased to self._currentStack* in run() and populated by
                # _populate_next_airplanes_wake_vortices. Ages are derived from row
                # position because each retained wake row is one delta_time older than
                # the row before it.
                num_chordwise_wake_rows = step
                if self._max_wake_rows is not None:
                    num_chordwise_wake_rows = min(step, self._max_wake_rows)
                num_wing_wake_vortices = num_chordwise_wake_rows * _num_spanwise_panels
                if num_wing_wake_vortices > 0:
                    block_start = global_wake_ring_vortex_position
                    block_end = block_start + num_wing_wake_vortices

                    # The initial core radius is constant across this Wing's
                    # wake block, so it fills in a single slice.
                    self._currentStackWakeRc0s[block_start:block_end] = wing_r_c0

                    # Each chordwise wake row is one delta_time older than the
                    # row shed after it, so row index c (0-based, newest
                    # first) has age (c + 1) * delta_time. Repeat each row's
                    # age across that row's spanwise wake vortices to match
                    # the block layout used by
                    # _populate_next_airplanes_wake_vortices.
                    row_ages = (
                        np.arange(1, num_chordwise_wake_rows + 1, dtype=float)
                        * self.delta_time
                    )
                    self._current_wake_vortex_ages[block_start:block_end] = np.repeat(
                        row_ages, _num_spanwise_panels
                    )

                    global_wake_ring_vortex_position += num_wing_wake_vortices

        if self._current_step > 0:
            last_step = self._current_step - 1
            lastStackFr = self._listStackFrbrvp_GP1_CgP1[last_step]
            lastStackFl = self._listStackFlbrvp_GP1_CgP1[last_step]
            lastStackBl = self._listStackBlbrvp_GP1_CgP1[last_step]
            lastStackBr = self._listStackBrbrvp_GP1_CgP1[last_step]

            # Bound corner stacks for the last step.
            self._lastStackFrbrvp_GP1_CgP1[:] = lastStackFr
            self._lastStackFlbrvp_GP1_CgP1[:] = lastStackFl
            self._lastStackBlbrvp_GP1_CgP1[:] = lastStackBl
            self._lastStackBrbrvp_GP1_CgP1[:] = lastStackBr

            # Last step bound leg center points, derived inline from corners.
            # Right leg: back right -> front right. Front leg: front right ->
            # front left. Left leg: front left -> back left. Back leg: back
            # left -> back right.
            self._lastStackCblvpr_GP1_CgP1[:] = 0.5 * (lastStackBr + lastStackFr)
            self._lastStackCblvpf_GP1_CgP1[:] = 0.5 * (lastStackFr + lastStackFl)
            self._lastStackCblvpl_GP1_CgP1[:] = 0.5 * (lastStackFl + lastStackBl)
            self._lastStackCblvpb_GP1_CgP1[:] = 0.5 * (lastStackBl + lastStackBr)

            # Last step Panel collocation points. Panel topology is invariant,
            # but Panel positions move when the geometry is unsteady, so these
            # need to be re-read from the previous step's Airplanes.
            global_panel_position = 0
            last_problem = self._get_steady_problem_at(last_step)
            for last_airplane in last_problem.airplanes:
                for last_wing in last_airplane.wings:
                    _last_panels = last_wing.panels
                    assert _last_panels is not None
                    for last_panel in np.ravel(_last_panels):
                        self._stackLastCpp_GP1_CgP1[global_panel_position, :] = (
                            last_panel.Cpp_GP1_CgP1
                        )
                        global_panel_position += 1

    def _calculate_wing_wing_influences(self) -> None:
        """Finds the current time step's SteadyProblem's 2D ndarray of Wing Wing
        influence coefficients (observed from the Earth frame).

        When an image surface is defined on the OperatingPoint, the influence
        coefficients also include the contributions from image bound ring vortices
        reflected across that surface.

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at each Panel's
        # collocation point by each bound ring vortex. The answer is normalized
        # because the solver's list of bound ring vortex strengths was initialized to
        # all be 1.0. This will be updated once the correct strengths are calculated.
        singularity_counts = np.zeros(4, dtype=np.int64)
        gridNormVIndCpp_GP1_E = (
            _aerodynamics_functions.expanded_velocities_from_ring_vortices(
                stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                strengths=self._current_bound_vortex_strengths,
                r_c0s=self._currentStackBoundRc0s,
                singularity_counts=singularity_counts,
                ages=None,
                nu=self.current_operating_point.nu,
            )
        )

        # Add the image contribution if an image surface is defined.
        surfaceReflect_T_act_GP1_CgP1 = (
            self.current_operating_point.surfaceReflect_T_act_GP1_CgP1
        )
        if surfaceReflect_T_act_GP1_CgP1 is not None:
            stackReflectedCpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                self.stackCpp_GP1_CgP1,
                has_point=True,
            )
            gridImageVIndCpp_GP1__E = (
                _aerodynamics_functions.expanded_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=stackReflectedCpp_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                    strengths=self._current_bound_vortex_strengths,
                    r_c0s=self._currentStackBoundRc0s,
                    singularity_counts=singularity_counts,
                    ages=None,
                    nu=self.current_operating_point.nu,
                )
            )
            gridNormVIndCpp_GP1_E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                gridImageVIndCpp_GP1__E,
                has_point=False,
            )

        unexpected_singularity_counts = np.copy(singularity_counts)

        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.ERROR,
            "_calculate_wing_wing_influences",
            unexpected_singularity_counts,
        )

        # Take the batch dot product of the normalized induced velocities (in the
        # first Airplane's geometry axes, observed from the Earth frame) with each
        # Panel's unit normal direction (in the first Airplane's geometry axes). This
        # is now the 2D ndarray of Wing Wing influence coefficients (observed from
        # the Earth frame).
        self._currentGridWingWingInfluences__E = np.einsum(
            "...k,...k->...",
            gridNormVIndCpp_GP1_E,
            np.expand_dims(self.stackUnitNormals_GP1, axis=1),
        )

    def _calculate_freestream_wing_influences(self) -> None:
        """Finds the 1D ndarray of freestream Wing influence coefficients (observed from
        the Earth frame) at the current time step.

        **Notes:**

        This method also includes the influence coefficients due to motion defined in
        Movement (observed from the Earth frame) at every collocation point.

        :return: None
        """
        # Find the normal components of the freestream only Wing influence
        # coefficients (observed from the Earth frame) at each Panel's collocation
        # point by taking a batch dot product.
        currentStackFreestreamOnlyWingInfluences__E = np.einsum(
            "ij,j->i",
            self.stackUnitNormals_GP1,
            self._currentVInf_GP1__E,
        )

        # Get the current apparent velocities at each Panel's collocation point due
        # to any motion defined in Movement (in the first Airplane's geometry axes,
        # observed from the Earth frame).
        currentStackMovementV_GP1_E = (
            self._calculate_current_movement_velocities_at_collocation_points()
        )

        # Get the current motion influence coefficients at each Panel's collocation
        # point (observed from the Earth frame) by taking a batch dot product.
        currentStackMovementInfluences__E = np.einsum(
            "ij,ij->i",
            self.stackUnitNormals_GP1,
            currentStackMovementV_GP1_E,
        )

        # Calculate the total current freestream Wing influence coefficients by
        # summing the freestream-only influence coefficients and the motion influence
        # coefficients (all observed from the Earth frame).
        self._currentStackFreestreamWingInfluences__E = (
            currentStackFreestreamOnlyWingInfluences__E
            + currentStackMovementInfluences__E
        )

    def _calculate_wake_wing_influences(self) -> None:
        """Finds the 1D ndarray of the wake Wing influence coefficients (observed from
        the Earth frame) at the current time step.

        When an image surface is defined on the OperatingPoint, the influence
        coefficients also include the contributions from image wake ring vortices
        reflected across that surface.

        **Notes:**

        If the current time step is the first time step, no wake has been shed, so this
        method will return zero for all the wake Wing influence coefficients (observed
        from the Earth frame).

        :return: None
        """
        if self._current_step > 0:
            # Get the velocities (in the first Airplane's geometry axes, observed
            # from the Earth frame) induced by the wake ring vortices at each Panel's
            # collocation point.
            singularity_counts = np.zeros(4, dtype=np.int64)
            currentStackWakeV_GP1_E = (
                _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self._currentStackBrwrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self._currentStackFrwrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self._currentStackFlwrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self._currentStackBlwrvp_GP1_CgP1,
                    strengths=self._current_wake_vortex_strengths,
                    r_c0s=self._currentStackWakeRc0s,
                    singularity_counts=singularity_counts,
                    ages=self._current_wake_vortex_ages,
                    nu=self.current_operating_point.nu,
                )
            )

            # Add the image contribution if an image surface is defined.
            surfaceReflect_T_act_GP1_CgP1 = (
                self.current_operating_point.surfaceReflect_T_act_GP1_CgP1
            )
            if surfaceReflect_T_act_GP1_CgP1 is not None:
                stackReflectedCpp_GP1_CgP1 = _transformations.apply_T_to_vectors(
                    surfaceReflect_T_act_GP1_CgP1,
                    self.stackCpp_GP1_CgP1,
                    has_point=True,
                )
                currentStackImageWakeV_GP1_E = (
                    _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                        stackP_GP1_CgP1=stackReflectedCpp_GP1_CgP1,
                        stackBrrvp_GP1_CgP1=self._currentStackBrwrvp_GP1_CgP1,
                        stackFrrvp_GP1_CgP1=self._currentStackFrwrvp_GP1_CgP1,
                        stackFlrvp_GP1_CgP1=self._currentStackFlwrvp_GP1_CgP1,
                        stackBlrvp_GP1_CgP1=self._currentStackBlwrvp_GP1_CgP1,
                        strengths=self._current_wake_vortex_strengths,
                        r_c0s=self._currentStackWakeRc0s,
                        singularity_counts=singularity_counts,
                        ages=self._current_wake_vortex_ages,
                        nu=self.current_operating_point.nu,
                    )
                )
                currentStackWakeV_GP1_E += _transformations.apply_T_to_vectors(
                    surfaceReflect_T_act_GP1_CgP1,
                    currentStackImageWakeV_GP1_E,
                    has_point=False,
                )

            unexpected_singularity_counts = np.copy(singularity_counts)

            _functions.log_unexpected_singularity_counts(
                _logger,
                logging.INFO,
                "_calculate_wake_wing_influences",
                unexpected_singularity_counts,
            )

            # Get the current wake Wing influence coefficients (observed from the
            # Earth frame) by taking a batch dot product with each Panel's normal
            # vector (in the first Airplane's geometry axes).
            self._currentStackWakeWingInfluences__E = np.einsum(
                "ij,ij->i", currentStackWakeV_GP1_E, self.stackUnitNormals_GP1
            )

        else:
            # If this is the first time step, set all the current Wake-wing influence
            # coefficients to 0.0 (observed from the Earth frame) because no wake
            # ring vortices have been shed.
            self._currentStackWakeWingInfluences__E = np.zeros(
                self.num_panels, dtype=float
            )

    def _calculate_vortex_strengths(self) -> None:
        """Solves for the strength of each Panel's bound ring vortex.

        :return: None
        """
        self._current_bound_vortex_strengths = np.linalg.solve(
            self._currentGridWingWingInfluences__E,
            -self._currentStackWakeWingInfluences__E
            - self._currentStackFreestreamWingInfluences__E,
        )

    def calculate_solution_velocity(
        self,
        stackP_GP1_CgP1: np.ndarray | Sequence[Sequence[float | int]],
        bound_singularity_counts: np.ndarray | None = None,
        wake_singularity_counts: np.ndarray | None = None,
    ) -> np.ndarray:
        """Finds the fluid velocity (in the first Airplane's geometry axes, observed
        from the Earth frame) at one or more points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG) due to the freestream velocity and
        the induced velocity from every ring vortex.

        When an image surface is defined on the OperatingPoint, the returned velocity
        also includes the induced velocity from image bound and wake ring vortices
        reflected across that surface.

        **Notes:**

        This method assumes that the correct strengths for the ring vortices have
        already been calculated and set.

        This method also does not include the velocity due to the Movement's motion at
        any of the points provided, as it has no way of knowing if any of the points lie
        on panels.

        :param stackP_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (N,3) representing the positions of the evaluation points (in the
            first Airplane's geometry axes, relative to the first Airplane's CG). Can be
            a tuple, list, or ndarray. Values are converted to floats internally. The
            units are in meters.
        :param bound_singularity_counts: An optional (4,) ndarray of int64 for
            accumulating singularity event counts from bound ring vortices. If None,
            counts are discarded.
        :param wake_singularity_counts: An optional (4,) ndarray of int64 for
            accumulating singularity event counts from wake ring vortices. If None,
            counts are discarded.
        :return: A (N,3) ndarray of floats representing the velocity (in the first
            Airplane's geometry axes, observed from the Earth frame) at each evaluation
            point due to the summed effects of the freestream velocity and the induced
            velocity from every ring vortex. The units are in meters per second.
        """
        stackP_GP1_CgP1 = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_GP1_CgP1, "stackP_GP1_CgP1"
            )
        )

        if bound_singularity_counts is None:
            bound_singularity_counts = np.zeros(4, dtype=np.int64)
        if wake_singularity_counts is None:
            wake_singularity_counts = np.zeros(4, dtype=np.int64)

        stackBoundRingVInd_GP1_E = (
            _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                strengths=self._current_bound_vortex_strengths,
                r_c0s=self._currentStackBoundRc0s,
                singularity_counts=bound_singularity_counts,
                ages=None,
                nu=self.current_operating_point.nu,
            )
        )
        stackWakeRingVInd_GP1_E = (
            _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self._currentStackBrwrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self._currentStackFrwrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self._currentStackFlwrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self._currentStackBlwrvp_GP1_CgP1,
                strengths=self._current_wake_vortex_strengths,
                r_c0s=self._currentStackWakeRc0s,
                singularity_counts=wake_singularity_counts,
                ages=self._current_wake_vortex_ages,
                nu=self.current_operating_point.nu,
            )
        )

        # Add the image contributions if an image surface is defined.
        surfaceReflect_T_act_GP1_CgP1 = (
            self.current_operating_point.surfaceReflect_T_act_GP1_CgP1
        )
        if surfaceReflect_T_act_GP1_CgP1 is not None:
            stackReflectedP_GP1_CgP1 = _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackP_GP1_CgP1,
                has_point=True,
            )
            stackImageBoundRingVInd_GP1_E = (
                _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=stackReflectedP_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                    strengths=self._current_bound_vortex_strengths,
                    r_c0s=self._currentStackBoundRc0s,
                    singularity_counts=bound_singularity_counts,
                    ages=None,
                    nu=self.current_operating_point.nu,
                )
            )
            stackBoundRingVInd_GP1_E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackImageBoundRingVInd_GP1_E,
                has_point=False,
            )
            stackImageWakeRingVInd_GP1_E = (
                _aerodynamics_functions.collapsed_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=stackReflectedP_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self._currentStackBrwrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self._currentStackFrwrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self._currentStackFlwrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self._currentStackBlwrvp_GP1_CgP1,
                    strengths=self._current_wake_vortex_strengths,
                    r_c0s=self._currentStackWakeRc0s,
                    singularity_counts=wake_singularity_counts,
                    ages=self._current_wake_vortex_ages,
                    nu=self.current_operating_point.nu,
                )
            )
            stackWakeRingVInd_GP1_E += _transformations.apply_T_to_vectors(
                surfaceReflect_T_act_GP1_CgP1,
                stackImageWakeRingVInd_GP1_E,
                has_point=False,
            )

        return cast(
            np.ndarray,
            stackBoundRingVInd_GP1_E
            + stackWakeRingVInd_GP1_E
            + self._currentVInf_GP1__E,
        )

    def _calculate_loads(self) -> None:
        """Calculates the forces (in the first Airplane's geometry axes) and moments (in
        the first Airplane's geometry axes, relative to the first Airplane's CG) on
        every Panel at the current time step.

        **Notes:**

        This method assumes that the correct strengths for the ring vortices have
        already been calculated and set.

        This method used to accidentally double-count the load on each Panel due to the
        left and right line vortex legs. Additionally, it didn't include contributions
        to the load on each Panel from their back line vortex legs. Thankfully, these
        issues only introduced small errors in most typical simulations. They have both
        now been fixed by (1) using a 1/2 factor for each "effective" vortex strength
        shared between two Panels, and (2) including the effects each Panel's back line
        vortex with its own effective strength.

        :return: None
        """
        # Initialize a variable to hold the global Panel position as we iterate
        # through them.
        global_panel_position = 0

        # Initialize three 1D ndarrays to hold the effective strength of the Panels'
        # ring vortices' line vortices.
        effective_right_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_front_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_left_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_back_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)

        # Iterate through the Airplanes' Wings. Within a Wing, Panels are laid
        # out in row major (chordwise outer, spanwise inner) order, so
        # neighbouring Panel strengths can be found at fixed offsets from the
        # current global Panel position: +1 right, -1 left,
        # +num_spanwise_panels back, -num_spanwise_panels front.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:
                _panels = wing.panels
                assert _panels is not None

                num_spanwise = wing.num_spanwise_panels
                assert num_spanwise is not None

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(_panels)

                # Iterate through this Wing's 1D ndarray of Panels.
                panel: _panel.Panel
                for panel in panels:
                    this_strength = self._current_bound_vortex_strengths[
                        global_panel_position
                    ]

                    if panel.is_right_edge:
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            this_strength
                        )
                    else:
                        # Set the effective right line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex
                        # strength and that of the Panel to the right.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._current_bound_vortex_strengths[
                                global_panel_position + 1
                            ]
                        ) / 2

                    if panel.is_leading_edge:
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            this_strength
                        )
                    else:
                        # Set the effective front line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex
                        # strength and that of the Panel in front of it.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._current_bound_vortex_strengths[
                                global_panel_position - num_spanwise
                            ]
                        ) / 2

                    if panel.is_left_edge:
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            this_strength
                        )
                    else:
                        # Set the effective left line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex
                        # strength and that of the Panel to the left.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._current_bound_vortex_strengths[
                                global_panel_position - 1
                            ]
                        ) / 2

                    if panel.is_trailing_edge:
                        if self._current_step == 0:
                            # No wake ring vortex exists yet to cancel out this
                            # Panel's back line vortex contribution, so the
                            # effective back strength is the full ring vortex
                            # strength.
                            effective_back_line_vortex_strengths[
                                global_panel_position
                            ] = this_strength
                        else:
                            # The Panel's back line vortex is partially cancelled
                            # by the front line vortex of the wake ring vortex
                            # immediately to its rear. That wake ring vortex
                            # carries the strength this Panel's bound ring vortex
                            # had last time step, so the effective back strength
                            # is the difference between this step and the last.
                            effective_back_line_vortex_strengths[
                                global_panel_position
                            ] = (
                                this_strength
                                - self._last_bound_vortex_strengths[
                                    global_panel_position
                                ]
                            )
                    else:
                        # Set the effective back line vortex strength to 1/2 the
                        # difference between this Panel's bound ring vortex
                        # strength and that of the Panel to the back.
                        effective_back_line_vortex_strengths[global_panel_position] = (
                            this_strength
                            - self._current_bound_vortex_strengths[
                                global_panel_position + num_spanwise
                            ]
                        ) / 2

                    global_panel_position += 1

        # Calculate the velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at the center of every Panels' ring vortex's right
        # line vortex, front line vortex, left line vortex, and back line vortex.
        bound_singularity_counts = np.zeros(4, dtype=np.int64)
        wake_singularity_counts = np.zeros(4, dtype=np.int64)
        stackVelocityRightLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(
                stackP_GP1_CgP1=self.stackCblvpr_GP1_CgP1,
                bound_singularity_counts=bound_singularity_counts,
                wake_singularity_counts=wake_singularity_counts,
            )
            + self._calculate_current_movement_velocities_at_right_leg_centers()
        )
        stackVelocityFrontLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(
                stackP_GP1_CgP1=self.stackCblvpf_GP1_CgP1,
                bound_singularity_counts=bound_singularity_counts,
                wake_singularity_counts=wake_singularity_counts,
            )
            + self._calculate_current_movement_velocities_at_front_leg_centers()
        )
        stackVelocityLeftLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(
                stackP_GP1_CgP1=self.stackCblvpl_GP1_CgP1,
                bound_singularity_counts=bound_singularity_counts,
                wake_singularity_counts=wake_singularity_counts,
            )
            + self._calculate_current_movement_velocities_at_left_leg_centers()
        )
        stackVelocityBackLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(
                stackP_GP1_CgP1=self.stackCblvpb_GP1_CgP1,
                bound_singularity_counts=bound_singularity_counts,
                wake_singularity_counts=wake_singularity_counts,
            )
            + self._calculate_current_movement_velocities_at_back_leg_centers()
        )

        unexpected_bound_singularity_counts = np.copy(bound_singularity_counts)
        unexpected_wake_singularity_counts = np.copy(wake_singularity_counts)

        # Subtract the expected structural collinearity before logging. For each Wing
        # with C chordwise and S spanwise Panels, the four leg center evaluations
        # produce (8 * C * S - 2 * C - 2 * S) bound collinearity singularities from
        # ring vortex self and adjacent shared edge pairs. When there is a wake (time
        # step > 0), each trailing edge Panel's back leg center is also collinear with
        # and on-filament for the first wake row's front leg, adding S wake collinearity
        # singularities per Wing.
        expected_bound_collinearity = 0
        expected_wake_collinearity = 0
        for airplane in self.current_airplanes:
            for wing in airplane.wings:
                num_chordwise = wing.num_chordwise_panels
                num_spanwise = wing.num_spanwise_panels
                assert num_spanwise is not None
                n = num_chordwise * num_spanwise
                expected_bound_collinearity += (
                    8 * n - 2 * num_chordwise - 2 * num_spanwise
                )
                if self._current_step > 0:
                    expected_wake_collinearity += num_spanwise
        unexpected_bound_singularity_counts[3] -= expected_bound_collinearity
        unexpected_wake_singularity_counts[3] -= expected_wake_collinearity
        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.ERROR,
            "_calculate_loads (bound)",
            unexpected_bound_singularity_counts,
        )
        _functions.log_unexpected_singularity_counts(
            _logger,
            logging.INFO,
            "_calculate_loads (wake)",
            unexpected_wake_singularity_counts,
        )

        # Using the effective line vortex strengths and the Kutta-Joukowski theorem,
        # find the forces (in the first Airplane's geometry axes) on the Panels'
        # ring vortex's right line vortex, front line vortex, left line vortex, and back
        # line vortex using the effective vortex strengths.
        rightLegForces_GP1 = (
            self.current_operating_point.rho
            * np.expand_dims(effective_right_line_vortex_strengths, axis=1)
            * _functions.numba_1d_explicit_cross(
                stackVelocityRightLineVortexCenters_GP1__E, self.stackRbrv_GP1
            )
        )
        frontLegForces_GP1 = (
            self.current_operating_point.rho
            * np.expand_dims(effective_front_line_vortex_strengths, axis=1)
            * _functions.numba_1d_explicit_cross(
                stackVelocityFrontLineVortexCenters_GP1__E, self.stackFbrv_GP1
            )
        )
        leftLegForces_GP1 = (
            self.current_operating_point.rho
            * np.expand_dims(effective_left_line_vortex_strengths, axis=1)
            * _functions.numba_1d_explicit_cross(
                stackVelocityLeftLineVortexCenters_GP1__E, self.stackLbrv_GP1
            )
        )
        backLegForces_GP1 = (
            self.current_operating_point.rho
            * np.expand_dims(effective_back_line_vortex_strengths, axis=1)
            * np.cross(
                stackVelocityBackLineVortexCenters_GP1__E,
                self.stackBbrv_GP1,
                axis=-1,
            )
        )

        # The unsteady force calculation below includes a negative sign to account for a
        # sign convention mismatch between Ptera Software and the reference literature.
        # Ptera Software defines ring vortices with counter-clockwise (CCW) vertex
        # ordering, while the references use clockwise (CW) ordering. Both define panel
        # normals as pointing upward. This convention difference only affects the
        # unsteady force term because it depends on both vortex strength and the normal
        # vector. When converting from CCW to CW, the strength changes sign but the
        # normal vector does not, requiring a sign correction. In contrast, steady
        # Kutta-Joukowski forces depend on the strength and the line vortex vectors. Both
        # have flipped signs, causing the negatives to cancel. See issue #27:
        # https://github.com/camUrban/PteraSoftware/issues/27

        # Calculate the unsteady component of the force on each Panel (in geometry
        # axes), which is derived from the unsteady Bernoulli equation.
        unsteady_forces_GP1 = -(
            self.current_operating_point.rho
            * np.expand_dims(
                (
                    self._current_bound_vortex_strengths
                    - self._last_bound_vortex_strengths
                ),
                axis=1,
            )
            * np.expand_dims(self.panel_areas, axis=1)
            * self.stackUnitNormals_GP1
            / self.delta_time
        )

        forces_GP1 = (
            rightLegForces_GP1
            + frontLegForces_GP1
            + leftLegForces_GP1
            + backLegForces_GP1
            + unsteady_forces_GP1
        )

        # Find the moments (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) on the Panels' ring vortex's right line vortex,
        # front line vortex, left line vortex, and back line vortex.
        rightLegMoments_GP1_CgP1 = _functions.numba_1d_explicit_cross(
            self.stackCblvpr_GP1_CgP1, rightLegForces_GP1
        )
        frontLegMoments_GP1_CgP1 = _functions.numba_1d_explicit_cross(
            self.stackCblvpf_GP1_CgP1, frontLegForces_GP1
        )
        leftLegMoments_GP1_CgP1 = _functions.numba_1d_explicit_cross(
            self.stackCblvpl_GP1_CgP1, leftLegForces_GP1
        )
        backLegMoments_GP1_CgP1 = _functions.numba_1d_explicit_cross(
            self.stackCblvpb_GP1_CgP1, backLegForces_GP1
        )

        # The unsteady moment is calculated at the collocation point because the
        # unsteady force acts on the bound ring vortex, whose center is at the
        # collocation point, not at the Panel's centroid.

        # Find the moments (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) due to the unsteady component of the force on each Panel.
        unsteady_moments_GP1_CgP1 = _functions.numba_1d_explicit_cross(
            self.stackCpp_GP1_CgP1, unsteady_forces_GP1
        )

        moments_GP1_CgP1 = (
            rightLegMoments_GP1_CgP1
            + frontLegMoments_GP1_CgP1
            + leftLegMoments_GP1_CgP1
            + backLegMoments_GP1_CgP1
            + unsteady_moments_GP1_CgP1
        )

        # TODO: Transform forces_GP1 and moments_GP1_CgP1 to each Airplane's local
        #  geometry axes before passing to process_solver_loads.
        _functions.process_solver_loads(self, forces_GP1, moments_GP1_CgP1)

    def _populate_next_airplanes_wake(self) -> None:
        """Updates the next time step's Airplanes' wakes.

        :return: None
        """
        # Populate the locations of the next time step's Airplanes' wake ring vortex
        # points.
        self._populate_next_airplanes_wake_vortex_points()

        # Populate the locations of the next time step's Airplanes' wake ring vortices.
        self._populate_next_airplanes_wake_vortices()

    def _populate_next_airplanes_wake_vortex_points(self) -> None:
        """Populates the locations of the next time step's Airplanes' wake ring vortex
        points.

        :return: None
        """
        # Check that this isn't the last time step.
        if self._current_step < self.num_steps - 1:
            bound_singularity_counts = np.zeros(4, dtype=np.int64)
            wake_singularity_counts = np.zeros(4, dtype=np.int64)

            # Get the next time step's Airplanes.
            next_problem: problems.SteadyProblem = self._get_steady_problem_at(
                self._current_step + 1
            )
            next_airplanes = next_problem.airplanes

            # Iterate through this time step's Airplanes' successor objects.
            for airplane_id, next_airplane in enumerate(next_airplanes):

                # Iterate through the next Airplane's Wings.
                for wing_id, next_wing in enumerate(next_airplane.wings):

                    # Get the Wings at this position from the current Airplane.
                    this_airplane = self.current_airplanes[airplane_id]
                    this_wing = this_airplane.wings[wing_id]

                    # Check if this is the first time step.
                    if self._current_step == 0:

                        # Get the current Wing's number of chordwise and spanwise
                        # panels.
                        num_spanwise_panels = this_wing.num_spanwise_panels
                        assert num_spanwise_panels is not None

                        num_chordwise_panels = this_wing.num_chordwise_panels

                        # Set the chordwise position to be at the trailing edge.
                        chordwise_panel_id = num_chordwise_panels - 1

                        # Initialize a ndarray to hold the points of the new row of
                        # wake ring vortices (in the first Airplane's geometry axes,
                        # relative to the first Airplane's CG).
                        newRowWrvp_GP1_CgP1 = np.zeros(
                            (1, num_spanwise_panels + 1, 3), dtype=float
                        )

                        next_step = self._current_step + 1
                        next_stackBl = self._listStackBlbrvp_GP1_CgP1[next_step]
                        next_stackBr = self._listStackBrbrvp_GP1_CgP1[next_step]
                        wing_panel_offset = self._per_wing_panel_offsets[airplane_id][
                            wing_id
                        ]
                        te_panel_base = (
                            wing_panel_offset + chordwise_panel_id * num_spanwise_panels
                        )

                        # Iterate through the spanwise Panel positions.
                        for spanwise_panel_id in range(num_spanwise_panels):
                            te_global_idx = te_panel_base + spanwise_panel_id

                            # The position of the new front left wake ring vortex's
                            # point is the next time step's Panel's bound
                            # ring vortex's back left point.
                            newRowWrvp_GP1_CgP1[0, spanwise_panel_id] = next_stackBl[
                                te_global_idx
                            ]

                            # If the Panel is at the right edge of the Wing, add
                            # its back right bound ring vortex point to the row
                            # of new wake ring vortex points.
                            if spanwise_panel_id == (num_spanwise_panels - 1):
                                newRowWrvp_GP1_CgP1[0, spanwise_panel_id + 1] = (
                                    next_stackBr[te_global_idx]
                                )

                        # Set the next time step's Wing's grid of wake ring vortex
                        # points to a copy of the row of new wake ring vortex points.
                        # This is correct because it is currently the first time step.
                        next_wing.gridWrvp_GP1_CgP1 = np.copy(newRowWrvp_GP1_CgP1)

                        # If the wake is prescribed, the velocity at every point is
                        # the freestream velocity (in the first Airplane's geometry
                        # axes, observed from the Earth frame). Otherwise, batch one
                        # solution-velocity call across all of the row's points.
                        if self._prescribed_wake:
                            vRowWrvp_GP1__E = self._currentVInf_GP1__E
                        else:
                            stackRowWrvp_GP1_CgP1 = next_wing.gridWrvp_GP1_CgP1.reshape(
                                -1, 3
                            )
                            stackVRowWrvp_GP1__E = self.calculate_solution_velocity(
                                stackRowWrvp_GP1_CgP1,
                                bound_singularity_counts=bound_singularity_counts,
                                wake_singularity_counts=wake_singularity_counts,
                            )
                            vRowWrvp_GP1__E = stackVRowWrvp_GP1__E.reshape(
                                next_wing.gridWrvp_GP1_CgP1.shape
                            )

                        # Build the second new row by advecting the first row.
                        secondNewRowWrvp_GP1_CgP1 = (
                            next_wing.gridWrvp_GP1_CgP1
                            + vRowWrvp_GP1__E * self.delta_time
                        )

                        # Update the next time step's Wing's grid of wake ring vortex
                        # points by vertically stacking the new second row below it.
                        next_wing.gridWrvp_GP1_CgP1 = np.vstack(
                            (
                                next_wing.gridWrvp_GP1_CgP1,
                                secondNewRowWrvp_GP1_CgP1,
                            )
                        )

                    # If this isn't the first time step, then do this.
                    else:
                        _thisGridWrvp_GP1_CgP1 = this_wing.gridWrvp_GP1_CgP1
                        assert _thisGridWrvp_GP1_CgP1 is not None

                        # Set the next time step's Wing's grid of wake ring vortex
                        # points to a copy of this time step's Wing's grid of wake
                        # ring vortex points.
                        next_wing.gridWrvp_GP1_CgP1 = np.copy(_thisGridWrvp_GP1_CgP1)

                        # If the wake is prescribed, the velocity at every point is
                        # the freestream velocity (in the first Airplane's geometry
                        # axes, observed from the Earth frame). Otherwise, batch one
                        # solution-velocity call across the entire aged grid.
                        if self._prescribed_wake:
                            vGridWrvp_GP1__E = self._currentVInf_GP1__E
                        else:
                            stackGridWrvp_GP1_CgP1 = (
                                next_wing.gridWrvp_GP1_CgP1.reshape(-1, 3)
                            )
                            stackVGridWrvp_GP1__E = self.calculate_solution_velocity(
                                stackGridWrvp_GP1_CgP1,
                                bound_singularity_counts=bound_singularity_counts,
                                wake_singularity_counts=wake_singularity_counts,
                            )
                            vGridWrvp_GP1__E = stackVGridWrvp_GP1__E.reshape(
                                next_wing.gridWrvp_GP1_CgP1.shape
                            )

                        # Advect the entire aged grid in one vector add.
                        next_wing.gridWrvp_GP1_CgP1 = (
                            next_wing.gridWrvp_GP1_CgP1
                            + vGridWrvp_GP1__E * self.delta_time
                        )

                        # Find the chordwise position of the Wing's trailing edge.
                        chordwise_panel_id = this_wing.num_chordwise_panels - 1

                        _num_spanwise_panels = this_wing.num_spanwise_panels
                        assert _num_spanwise_panels is not None

                        # Initialize a new ndarray to hold the new row of wake
                        # ring vortex vertices.
                        newRowWrvp_GP1_CgP1 = np.zeros(
                            (1, _num_spanwise_panels + 1, 3), dtype=float
                        )

                        next_step = self._current_step + 1
                        next_stackBl = self._listStackBlbrvp_GP1_CgP1[next_step]
                        next_stackBr = self._listStackBrbrvp_GP1_CgP1[next_step]
                        wing_panel_offset = self._per_wing_panel_offsets[airplane_id][
                            wing_id
                        ]
                        te_panel_base = (
                            wing_panel_offset
                            + chordwise_panel_id * _num_spanwise_panels
                        )

                        # Iterate spanwise through the trailing edge Panels.
                        for spanwise_panel_id in range(_num_spanwise_panels):
                            te_global_idx = te_panel_base + spanwise_panel_id

                            # Add the Panel's back left bound ring vortex point to
                            # the grid of new wake ring vortex points.
                            newRowWrvp_GP1_CgP1[0, spanwise_panel_id] = next_stackBl[
                                te_global_idx
                            ]

                            # If the Panel is at the right edge of the Wing, add
                            # its back right bound ring vortex point to the grid
                            # of new wake ring vortex vertices.
                            if spanwise_panel_id == (_num_spanwise_panels - 1):
                                newRowWrvp_GP1_CgP1[0, spanwise_panel_id + 1] = (
                                    next_stackBr[te_global_idx]
                                )

                        # Stack the new row of wake ring vortex points above the
                        # Wing's grid of wake ring vortex points.
                        next_wing.gridWrvp_GP1_CgP1 = np.vstack(
                            (
                                newRowWrvp_GP1_CgP1,
                                next_wing.gridWrvp_GP1_CgP1,
                            )
                        )

                        # If wake truncation is enabled, discard the oldest (most
                        # downstream) rows of wake points. The point grid has one more
                        # row than the vortex grid because points form vertices and
                        # vortices form cells.
                        if (
                            self._max_wake_rows is not None
                            and next_wing.gridWrvp_GP1_CgP1.shape[0]
                            > self._max_wake_rows + 1
                        ):
                            next_wing.gridWrvp_GP1_CgP1 = next_wing.gridWrvp_GP1_CgP1[
                                : self._max_wake_rows + 1
                            ]

            unexpected_bound_singularity_counts = np.copy(bound_singularity_counts)
            unexpected_wake_singularity_counts = np.copy(wake_singularity_counts)

            _functions.log_unexpected_singularity_counts(
                _logger,
                logging.DEBUG,
                "_populate_next_airplanes_wake_vortex_points (bound)",
                unexpected_bound_singularity_counts,
            )
            _functions.log_unexpected_singularity_counts(
                _logger,
                logging.DEBUG,
                "_populate_next_airplanes_wake_vortex_points (wake)",
                unexpected_wake_singularity_counts,
            )

    def _populate_next_airplanes_wake_vortices(self) -> None:
        """Populates the next time step's wake corner stacks and strengths.

        Each Wing's wake POINT grid (Wing.gridWrvp_GP1_CgP1) was already advected by
        _populate_next_airplanes_wake_vortex_points. This method derives the next step's
        per wake-cell corner positions and strengths and writes them directly into the
        per step list arrays. Strengths come from this step's solved bound ring vortex
        strengths (for the new front row) and from this step's wake (for inherited
        rows). The oldest row of this step's wake is dropped when truncation is in
        effect.

        :return: None
        """
        if self._current_step >= self.num_steps - 1:
            return

        this_step = self._current_step
        next_step = self._current_step + 1

        next_num_chordwise_rows = next_step
        if self._max_wake_rows is not None:
            next_num_chordwise_rows = min(next_step, self._max_wake_rows)

        this_num_chordwise_rows = this_step
        if self._max_wake_rows is not None:
            this_num_chordwise_rows = min(this_step, self._max_wake_rows)

        # Output buffers for the next step.
        next_strengths = self._list_wake_vortex_strengths[next_step]
        next_stackFr = self.listStackFrwrvp_GP1_CgP1[next_step]
        next_stackFl = self.listStackFlwrvp_GP1_CgP1[next_step]
        next_stackBl = self.listStackBlwrvp_GP1_CgP1[next_step]
        next_stackBr = self.listStackBrwrvp_GP1_CgP1[next_step]

        # This step's wake snapshots, used as the source of inherited rows.
        # When this_num_chordwise_rows is 0 (i.e. step 0 -> 1), they are zero
        # length and never indexed.
        this_strengths = self._list_wake_vortex_strengths[this_step]

        next_problem = self._get_steady_problem_at(next_step)
        for airplane_id, next_airplane in enumerate(next_problem.airplanes):
            for wing_id, next_wing in enumerate(next_airplane.wings):
                wing_num_spanwise = self._per_wing_num_spanwise_panels[airplane_id][
                    wing_id
                ]
                cumulative_prior_spanwise = self._per_wing_spanwise_cumsum[airplane_id][
                    wing_id
                ]

                # Per wing flat block bases inside this step's and the next
                # step's global wake stacks.
                next_wing_wake_base = (
                    next_num_chordwise_rows * cumulative_prior_spanwise
                )
                this_wing_wake_base = (
                    this_num_chordwise_rows * cumulative_prior_spanwise
                )

                # Trailing edge bound panel global positions for this Wing.
                te_chordwise = (
                    self._per_wing_num_chordwise_panels[airplane_id][wing_id] - 1
                )
                te_panel_base = (
                    self._per_wing_panel_offsets[airplane_id][wing_id]
                    + te_chordwise * wing_num_spanwise
                )

                nextGridWrvp_GP1_CgP1 = next_wing.gridWrvp_GP1_CgP1
                assert nextGridWrvp_GP1_CgP1 is not None
                assert nextGridWrvp_GP1_CgP1.shape[0] == next_num_chordwise_rows + 1
                assert nextGridWrvp_GP1_CgP1.shape[1] == wing_num_spanwise + 1

                # Front row (c == 0): newly shed.
                front_start = next_wing_wake_base
                front_end = front_start + wing_num_spanwise
                next_stackFl[front_start:front_end] = nextGridWrvp_GP1_CgP1[
                    0, :wing_num_spanwise
                ]
                next_stackFr[front_start:front_end] = nextGridWrvp_GP1_CgP1[
                    0, 1 : wing_num_spanwise + 1
                ]
                next_stackBl[front_start:front_end] = nextGridWrvp_GP1_CgP1[
                    1, :wing_num_spanwise
                ]
                next_stackBr[front_start:front_end] = nextGridWrvp_GP1_CgP1[
                    1, 1 : wing_num_spanwise + 1
                ]
                next_strengths[front_start:front_end] = (
                    self._current_bound_vortex_strengths[
                        te_panel_base : te_panel_base + wing_num_spanwise
                    ]
                )

                # Inherited rows (c >= 1): aged versions of this step's rows
                # 0 .. inherited_rows - 1. When wake truncation drops the
                # oldest row, the highest-c row of this step is excluded.
                for c in range(1, next_num_chordwise_rows):
                    old_c = c - 1
                    new_start = next_wing_wake_base + c * wing_num_spanwise
                    new_end = new_start + wing_num_spanwise
                    old_start = this_wing_wake_base + old_c * wing_num_spanwise
                    old_end = old_start + wing_num_spanwise

                    next_stackFl[new_start:new_end] = nextGridWrvp_GP1_CgP1[
                        c, :wing_num_spanwise
                    ]
                    next_stackFr[new_start:new_end] = nextGridWrvp_GP1_CgP1[
                        c, 1 : wing_num_spanwise + 1
                    ]
                    next_stackBl[new_start:new_end] = nextGridWrvp_GP1_CgP1[
                        c + 1, :wing_num_spanwise
                    ]
                    next_stackBr[new_start:new_end] = nextGridWrvp_GP1_CgP1[
                        c + 1, 1 : wing_num_spanwise + 1
                    ]
                    next_strengths[new_start:new_end] = this_strengths[
                        old_start:old_end
                    ]

    def _calculate_current_movement_velocities_at_collocation_points(
        self,
    ) -> np.ndarray:
        """Finds the apparent velocities (in the first Airplane's geometry axes,
        observed from the Earth frame) at each Panel's collocation point due to any
        motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at each
            Panel's collocation point due to any motion defined in Movement. If the
            current time step is the first time step, these velocities will all be all
            zeros. Its units are in meters per second.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return cast(
            np.ndarray,
            -(self.stackCpp_GP1_CgP1 - self._stackLastCpp_GP1_CgP1) / self.delta_time,
        )

    def _calculate_current_movement_velocities_at_right_leg_centers(self) -> np.ndarray:
        """Finds the apparent velocities (in the first Airplane's geometry axes,
        observed from the Earth frame) at the center point of each bound ring vortex's
        right leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound ring vortex's right leg due to any motion defined in
            Movement. If the current time step is the first time step, these velocities
            will all be all zeros. Its units are in meters per second.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return cast(
            np.ndarray,
            -(self.stackCblvpr_GP1_CgP1 - self._lastStackCblvpr_GP1_CgP1)
            / self.delta_time,
        )

    def _calculate_current_movement_velocities_at_front_leg_centers(self) -> np.ndarray:
        """Finds the apparent velocities (in the first Airplane's geometry axes,
        observed from the Earth frame) at the center point of each bound ring vortex's
        front leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound ring vortex's front leg due to any motion defined in
            Movement. If the current time step is the first time step, these velocities
            will all be all zeros. Its units are in meters per second.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return cast(
            np.ndarray,
            -(self.stackCblvpf_GP1_CgP1 - self._lastStackCblvpf_GP1_CgP1)
            / self.delta_time,
        )

    def _calculate_current_movement_velocities_at_left_leg_centers(self) -> np.ndarray:
        """Finds the apparent velocities (in the first Airplane's geometry axes,
        observed from the Earth frame) at the center point of each bound ring vortex's
        left leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound ring vortex's left leg due to any motion defined in
            Movement. If the current time step is the first time step, these velocities
            will all be all zeros. Its units are in meters per second.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return cast(
            np.ndarray,
            -(self.stackCblvpl_GP1_CgP1 - self._lastStackCblvpl_GP1_CgP1)
            / self.delta_time,
        )

    def _calculate_current_movement_velocities_at_back_leg_centers(self) -> np.ndarray:
        """Finds the apparent velocities (in the first Airplane's geometry axes,
        observed from the Earth frame) at the center point of each bound ring vortex's
        back leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound ring vortex's back leg due to any motion defined in
            Movement. If the current time step is the first time step, these velocities
            will all be all zeros. Its units are in meters per second.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return cast(
            np.ndarray,
            -(self.stackCblvpb_GP1_CgP1 - self._lastStackCblvpb_GP1_CgP1)
            / self.delta_time,
        )

    def _finalize_loads(self) -> None:
        """For cases with static geometry, finds the final loads and load coefficients
        for each of the SteadyProblem's Airplanes. For cases with variable geometry,
        finds the final cycle-averaged and cycle-root-mean-squared loads and load
        coefficients for each of the SteadyProblem's Airplanes.

        :return: None
        """
        # Get this solver's time step characteristics. Note that the first time step
        # ( time step 0), occurs at 0 seconds.
        num_steps_to_average = self.num_steps - self._first_averaging_step

        # Determine if this SteadyProblem's geometry is static or variable.
        static = self.unsteady_problem.movement.static

        # Initialize ndarrays to hold each Airplane's loads and load coefficients at
        # each of the time steps that calculated the loads.
        forces_W = np.zeros((self.num_airplanes, 3, num_steps_to_average), dtype=float)
        force_coefficients_W = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average), dtype=float
        )
        moments_W_CgP1 = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average), dtype=float
        )
        moment_coefficients_W_CgP1 = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average), dtype=float
        )

        # Initialize a variable to track position in the loads ndarrays.
        results_step = 0

        # Iterate through the time steps with loads and add the loads to their
        # respective ndarrays.
        for step in range(self._first_averaging_step, self.num_steps):

            # Get the Airplanes from the SteadyProblem at this time step.
            this_steady_problem: problems.SteadyProblem = self._get_steady_problem_at(
                step
            )
            these_airplanes = this_steady_problem.airplanes

            # Iterate through this time step's Airplanes.
            for airplane_id, airplane in enumerate(these_airplanes):
                forces_W[airplane_id, :, results_step] = airplane.forces_W
                force_coefficients_W[airplane_id, :, results_step] = (
                    airplane.forceCoefficients_W
                )
                moments_W_CgP1[airplane_id, :, results_step] = airplane.moments_W_CgP1
                moment_coefficients_W_CgP1[airplane_id, :, results_step] = (
                    airplane.momentCoefficients_W_CgP1
                )

            results_step += 1

        # For each Airplane, calculate and then save the final or cycle-averaged and
        # RMS loads and load coefficients. For variable geometry cases, use the
        # trapezoidal rule to compute the time-averaged mean and RMS over the final
        # cycle.
        first_problem: problems.SteadyProblem = self._get_steady_problem_at(0)
        for airplane_id, airplane in enumerate(first_problem.airplanes):
            if static:
                self.unsteady_problem.finalForces_W.append(forces_W[airplane_id, :, -1])
                self.unsteady_problem.finalForceCoefficients_W.append(
                    force_coefficients_W[airplane_id, :, -1]
                )
                self.unsteady_problem.finalMoments_W_CgP1.append(
                    moments_W_CgP1[airplane_id, :, -1]
                )
                self.unsteady_problem.finalMomentCoefficients_W_CgP1.append(
                    moment_coefficients_W_CgP1[airplane_id, :, -1]
                )
            else:
                # The number of intervals for the trapezoidal rule is one less
                # than the number of samples.
                num_intervals = num_steps_to_average - 1

                self.unsteady_problem.finalMeanForces_W.append(
                    np.trapezoid(forces_W[airplane_id], axis=-1) / num_intervals
                )
                self.unsteady_problem.finalMeanForceCoefficients_W.append(
                    np.trapezoid(force_coefficients_W[airplane_id], axis=-1)
                    / num_intervals
                )
                self.unsteady_problem.finalMeanMoments_W_CgP1.append(
                    np.trapezoid(moments_W_CgP1[airplane_id], axis=-1) / num_intervals
                )
                self.unsteady_problem.finalMeanMomentCoefficients_W_CgP1.append(
                    np.trapezoid(moment_coefficients_W_CgP1[airplane_id], axis=-1)
                    / num_intervals
                )

                self.unsteady_problem.finalRmsForces_W.append(
                    np.sqrt(
                        np.trapezoid(
                            np.square(forces_W[airplane_id]),
                            axis=-1,
                        )
                        / num_intervals
                    )
                )
                self.unsteady_problem.finalRmsForceCoefficients_W.append(
                    np.sqrt(
                        np.trapezoid(
                            np.square(force_coefficients_W[airplane_id]),
                            axis=-1,
                        )
                        / num_intervals
                    )
                )
                self.unsteady_problem.finalRmsMoments_W_CgP1.append(
                    np.sqrt(
                        np.trapezoid(
                            np.square(moments_W_CgP1[airplane_id]),
                            axis=-1,
                        )
                        / num_intervals
                    )
                )
                self.unsteady_problem.finalRmsMomentCoefficients_W_CgP1.append(
                    np.sqrt(
                        np.trapezoid(
                            np.square(moment_coefficients_W_CgP1[airplane_id]),
                            axis=-1,
                        )
                        / num_intervals
                    )
                )

    def _get_steady_problem_at(self, step: int) -> problems.SteadyProblem:
        """Gets the SteadyProblem at a given time step.

        Dynamic dispatch is used with _CoreUnsteadyProblems to provide different ways of
        accessing SteadyProblems based on the solver type without added code
        duplication. However, other methods must behave the same way regardless of
        solver type.

        :param step: The time step of the desired SteadyProblem.
        :return: The SteadyProblem at the given time step.
        """
        return self.steady_problems[step]

    def _initialize_panel_vortices_at(self, step: int) -> None:
        """Calculates the bound ring vortex corner positions at a given time step and
        stores them in the per step list arrays.

        :param step: The time step at which to initialize the Panels' bound ring vortex
            corner positions.
        :return: None
        """
        steady_problem = self._get_steady_problem_at(step)

        # Find the freestream velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at this time step.
        this_operating_point = steady_problem.operating_point
        vInf_GP1__E = this_operating_point.vInf_GP1__E

        stackFr = self._listStackFrbrvp_GP1_CgP1[step]
        stackFl = self._listStackFlbrvp_GP1_CgP1[step]
        stackBl = self._listStackBlbrvp_GP1_CgP1[step]
        stackBr = self._listStackBrbrvp_GP1_CgP1[step]

        # Iterate through this SteadyProblem's Airplanes' Wings.
        global_panel_position = 0
        for airplane_id, airplane in enumerate(steady_problem.airplanes):
            for wing_id, wing in enumerate(airplane.wings):
                _num_spanwise_panels = wing.num_spanwise_panels
                assert _num_spanwise_panels is not None

                _panels = wing.panels
                assert _panels is not None

                # Iterate through the Wing's chordwise and spanwise positions.
                for chordwise_position in range(wing.num_chordwise_panels):
                    for spanwise_position in range(_num_spanwise_panels):
                        # Pull the Panel out of the Wing's 2D ndarray of Panels.
                        panel: _panel.Panel = _panels[
                            chordwise_position, spanwise_position
                        ]

                        _Flbvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1
                        assert _Flbvp_GP1_CgP1 is not None

                        _Frbvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1
                        assert _Frbvp_GP1_CgP1 is not None

                        # Front bound ring vortex corner points coincide with the
                        # Panel's front bound vortex points.
                        Flrvp_GP1_CgP1 = _Flbvp_GP1_CgP1
                        Frrvp_GP1_CgP1 = _Frbvp_GP1_CgP1

                        # Define the location of the back left and back right
                        # bound ring vortex points based on whether the Panel is
                        # along the trailing edge.
                        if not panel.is_trailing_edge:
                            next_chordwise_panel: _panel.Panel = _panels[
                                chordwise_position + 1, spanwise_position
                            ]

                            _nextFlbvp_GP1_CgP1 = next_chordwise_panel.Flbvp_GP1_CgP1
                            assert _nextFlbvp_GP1_CgP1 is not None

                            _nextFrbvp_GP1_CgP1 = next_chordwise_panel.Frbvp_GP1_CgP1
                            assert _nextFrbvp_GP1_CgP1 is not None

                            Blrvp_GP1_CgP1 = _nextFlbvp_GP1_CgP1
                            Brrvp_GP1_CgP1 = _nextFrbvp_GP1_CgP1
                        else:
                            # As these vertices are directly behind the trailing
                            # edge, they are spaced back from their Panel's vertex
                            # by one quarter of the distance traveled by the
                            # trailing edge during a time step. This is to more
                            # accurately predict drag. More information can be
                            # found on pages 37-39 of "Modeling of aerodynamic
                            # forces in flapping flight with the Unsteady Vortex
                            # Lattice Method" by Thomas Lambert.
                            if step == 0:
                                _Blpp_GP1_CgP1 = panel.Blpp_GP1_CgP1
                                assert _Blpp_GP1_CgP1 is not None

                                _Brpp_GP1_CgP1 = panel.Brpp_GP1_CgP1
                                assert _Brpp_GP1_CgP1 is not None

                                Blrvp_GP1_CgP1 = (
                                    _Blpp_GP1_CgP1
                                    + vInf_GP1__E * self.delta_time * 0.25
                                )
                                Brrvp_GP1_CgP1 = (
                                    _Brpp_GP1_CgP1
                                    + vInf_GP1__E * self.delta_time * 0.25
                                )
                            else:
                                last_steady_problem = self._get_steady_problem_at(
                                    step - 1
                                )
                                last_airplane = last_steady_problem.airplanes[
                                    airplane_id
                                ]
                                last_wing = last_airplane.wings[wing_id]

                                _last_panels = last_wing.panels
                                assert _last_panels is not None

                                last_panel: _panel.Panel = _last_panels[
                                    chordwise_position, spanwise_position
                                ]

                                _thisBlpp_GP1_CgP1 = panel.Blpp_GP1_CgP1
                                assert _thisBlpp_GP1_CgP1 is not None

                                _lastBlpp_GP1_CgP1 = last_panel.Blpp_GP1_CgP1
                                assert _lastBlpp_GP1_CgP1 is not None

                                # Subtract (thisBlpp_GP1_CgP1 - lastBlpp_GP1_CgP1)
                                # / self.delta_time from vInf_GP1__E to get the
                                # apparent fluid velocity due to motion (observed
                                # in the Earth frame, in the first Airplane's
                                # geometry axes). This is the vector pointing
                                # opposite the velocity from motion.
                                Blrvp_GP1_CgP1 = (
                                    _thisBlpp_GP1_CgP1
                                    + (
                                        vInf_GP1__E
                                        - (_thisBlpp_GP1_CgP1 - _lastBlpp_GP1_CgP1)
                                        / self.delta_time
                                    )
                                    * self.delta_time
                                    * 0.25
                                )

                                _thisBrpp_GP1_CgP1 = panel.Brpp_GP1_CgP1
                                assert _thisBrpp_GP1_CgP1 is not None

                                _lastBrpp_GP1_CgP1 = last_panel.Brpp_GP1_CgP1
                                assert _lastBrpp_GP1_CgP1 is not None

                                Brrvp_GP1_CgP1 = (
                                    _thisBrpp_GP1_CgP1
                                    + (
                                        vInf_GP1__E
                                        - (_thisBrpp_GP1_CgP1 - _lastBrpp_GP1_CgP1)
                                        / self.delta_time
                                    )
                                    * self.delta_time
                                    * 0.25
                                )

                        stackFr[global_panel_position, :] = Frrvp_GP1_CgP1
                        stackFl[global_panel_position, :] = Flrvp_GP1_CgP1
                        stackBl[global_panel_position, :] = Blrvp_GP1_CgP1
                        stackBr[global_panel_position, :] = Brrvp_GP1_CgP1

                        global_panel_position += 1
