"""Contains the CoupledUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

CoupledUnsteadyRingVortexLatticeMethodSolver: A subclass of
UnsteadyRingVortexLatticeMethodSolver that solves CoupledUnsteadyProblems using the
unsteady ring vortex lattice method. This solver handles step-by-step geometry
initialization and supports two-way coupling between the aerodynamic solver and external
models.

**Contains the following functions:**

None
"""

from __future__ import annotations

import logging

import numpy as np
from tqdm import tqdm

from . import (
    _functions,
    _logging,
    _parameter_validation,
    geometry,
    problems,
)
from .unsteady_ring_vortex_lattice_method import UnsteadyRingVortexLatticeMethodSolver

_logger = _logging.get_logger("coupled_unsteady_ring_vortex_lattice_method")


class CoupledUnsteadyRingVortexLatticeMethodSolver(
    UnsteadyRingVortexLatticeMethodSolver
):
    """A subclass of UnsteadyRingVortexLatticeMethodSolver that solves
    CoupledUnsteadyProblems.

    This solver handles CoupledUnsteadyProblems where geometry is initialized and
    updated on a per-step basis (step-by-step), rather than being fully precomputed as
    in the parent class. At each time step, the solver generates the current step's
    geometry from the CoupledUnsteadyProblem, solves the aerodynamic system, then calls
    the problem's initialize_next_problem method to generate the next step's geometry.
    This enables two-way coupling with external models (structural deformation, rigid
    body dynamics) that modify geometry based on the current aerodynamic solution.

    **Key differences from parent UnsteadyRingVortexLatticeMethodSolver:**

    - Inherits core aerodynamic solver logic from parent (wall-built inheritance). -
    Overrides _get_steady_problem_at() to dynamically retrieve problems from the
    CoupledUnsteadyProblem. - Initializes bound RingVortices per-step inside the loop,
    rather than all-at-once   before the loop. - Calls
    CoupledUnsteadyProblem.initialize_next_problem() between time steps to   allow
    external models to update geometry.

    **Inherited methods (used directly from parent):**

    calculate_solution_velocity: Finds the fluid velocity at one or more points due to
    freestream and induced velocity from every RingVortex.

    All movement velocity calculation methods, aerodynamic influence calculation
    methods, vortex strength calculation methods, load calculation methods, and wake
    population methods.

    **Custom methods:**

    run: Runs the solver on the CoupledUnsteadyProblem with per-step geometry
    initialization and coupling hooks.

    initialize_step_geometry: Initializes geometry for a specific step without solving.
    """

    def __init__(
        self,
        coupled_unsteady_problem: problems.CoupledUnsteadyProblem,
    ) -> None:
        """The initialization method.

        :param coupled_unsteady_problem: The CoupledUnsteadyProblem to be solved.
            SteadyProblems are retrieved dynamically from this problem during iteration
            via _get_steady_problem_at().
        :return: None
        """
        if not isinstance(coupled_unsteady_problem, problems.CoupledUnsteadyProblem):
            raise TypeError(
                "coupled_unsteady_problem must be a CoupledUnsteadyProblem."
            )

        # self.coupled_unsteady_problem must be defined before the call to
        # super().__init__() because the parent class's __init__ method calls
        # _get_steady_problem_at(), which this class overrides to dispatch through
        # self.coupled_unsteady_problem.
        self.coupled_unsteady_problem = coupled_unsteady_problem
        super().__init__(coupled_unsteady_problem)

        # Store computed SteadyProblems for each time step. After the solve completes,
        # this list is converted to a tuple and assigned to self.steady_problems.
        self.steady_problems_data_storage: list[problems.SteadyProblem] = []

    def run(
        self,
        prescribed_wake: bool | np.bool_ = True,
        calculate_streamlines: bool | np.bool_ = True,
        show_progress: bool | np.bool_ = True,
    ) -> None:
        """Runs the solver on the CoupledUnsteadyProblem.

        Unlike the parent class, which precomputes all geometry before the main loop,
        this method initializes geometry per-step and calls the CoupledUnsteadyProblem's
        initialize_next_problem method between steps to allow external models to update
        the geometry.

        :param prescribed_wake: Set this to True to solve using a prescribed wake model.
            Set to False to use a free-wake, which may be more accurate but will make
            the run method significantly slower. Can be a bool or a numpy bool and will
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

        # Cache the wing geometry from the initial step. Unlike the parent class (which
        # precomputes all steps), coupled geometry is only known for step 0 at this
        # point. The number of panels is assumed constant across all steps.
        this_problem: problems.SteadyProblem = self._get_steady_problem_at(0)
        these_airplanes = this_problem.airplanes
        num_wing_panels = 0
        these_wings: list[tuple[geometry.wing.Wing, ...]] = []
        for airplane in these_airplanes:
            these_wings.append(airplane.wings)
            num_wing_panels += airplane.num_panels

        # Iterate through the Wings to get the total number of spanwise Panels.
        this_num_spanwise_panels = 0
        for this_wing_set in these_wings:
            for this_wing in this_wing_set:
                _this_wing_num_spanwise_panels = this_wing.num_spanwise_panels
                assert _this_wing_num_spanwise_panels is not None

                this_num_spanwise_panels += _this_wing_num_spanwise_panels

        # Pre-allocate wake arrays for all time steps using the initial geometry's
        # spanwise panel count (which is constant across steps).
        for step in range(self.num_steps):
            this_num_chordwise_wake_rows = step
            if self._max_wake_rows is not None:
                this_num_chordwise_wake_rows = min(step, self._max_wake_rows)
            this_num_wake_ring_vortices = (
                this_num_chordwise_wake_rows * this_num_spanwise_panels
            )

            this_wake_ring_vortex_strengths = np.zeros(
                this_num_wake_ring_vortices, dtype=float
            )
            this_wake_ring_vortex_ages = np.zeros(
                this_num_wake_ring_vortices, dtype=float
            )
            thisStackBrwrvp_GP1_CgP1 = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            thisStackFrwrvp_GP1_CgP1 = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            thisStackFlwrvp_GP1_CgP1 = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            thisStackBlwrvp_GP1_CgP1 = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            this_wake_rc0s = np.zeros(this_num_wake_ring_vortices, dtype=float)

            self.list_num_wake_vortices.append(this_num_wake_ring_vortices)
            self._list_wake_vortex_strengths.append(this_wake_ring_vortex_strengths)
            self._list_wake_vortex_ages.append(this_wake_ring_vortex_ages)
            self.listStackBrwrvp_GP1_CgP1.append(thisStackBrwrvp_GP1_CgP1)
            self.listStackFrwrvp_GP1_CgP1.append(thisStackFrwrvp_GP1_CgP1)
            self.listStackFlwrvp_GP1_CgP1.append(thisStackFlwrvp_GP1_CgP1)
            self.listStackBlwrvp_GP1_CgP1.append(thisStackBlwrvp_GP1_CgP1)
            self._list_wake_rc0s.append(this_wake_rc0s)

        # Estimate progress bar timing using the initial geometry's panel count.
        approx_times = np.zeros(self.num_steps + 1, dtype=float)
        for step in range(self.num_steps):
            if step != 0:
                num_wing_ring_vortices = num_wing_panels
                num_wake_ring_vortices = self.list_num_wake_vortices[step]
                num_ring_vortices = num_wing_ring_vortices + num_wake_ring_vortices

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
            bar.update(n=float(approx_times[0]))

            for step in range(self.num_steps):

                self._current_step = step
                current_problem: problems.SteadyProblem = self._get_steady_problem_at(
                    self._current_step
                )

                # Initialize the current step's bound RingVortices per-step (not
                # all-at-once like the parent).
                _logger.debug(f"Initializing step {step}'s RingVortices.")
                self._initialize_panel_vortices_at(step)
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

                # Initialize per-step aerodynamic and geometric arrays.
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
                self._last_bound_vortex_strengths = np.zeros(
                    self.num_panels, dtype=float
                )

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

                self.panel_is_trailing_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_leading_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_left_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_right_edge = np.zeros(self.num_panels, dtype=bool)

                self._current_wake_vortex_strengths = self._list_wake_vortex_strengths[
                    step
                ]
                self._current_wake_vortex_ages = self._list_wake_vortex_ages[step]
                self._currentStackBrwrvp_GP1_CgP1 = self.listStackBrwrvp_GP1_CgP1[step]
                self._currentStackFrwrvp_GP1_CgP1 = self.listStackFrwrvp_GP1_CgP1[step]
                self._currentStackFlwrvp_GP1_CgP1 = self.listStackFlwrvp_GP1_CgP1[step]
                self._currentStackBlwrvp_GP1_CgP1 = self.listStackBlwrvp_GP1_CgP1[step]
                self._currentStackBoundRc0s = np.zeros(self.num_panels, dtype=float)
                self._currentStackWakeRc0s = self._list_wake_rc0s[step]
                self.stackSeedPoints_GP1_CgP1 = np.zeros((0, 3), dtype=float)

                # Collapse the geometry matrices into 1D ndarrays of attributes.
                _logger.debug("Collapsing the geometry.")
                self._collapse_geometry()

                _logger.debug("Calculating the Wing Wing influences.")
                self._calculate_wing_wing_influences()

                _logger.debug("Calculating the freestream Wing influences.")
                self._calculate_freestream_wing_influences()

                _logger.debug("Calculating the wake Wing influences.")
                self._calculate_wake_wing_influences()

                _logger.debug("Calculating bound RingVortex strengths.")
                self._calculate_vortex_strengths()

                if self._current_step >= self.first_results_step:
                    _logger.debug("Calculating forces and moments.")
                    self._calculate_loads()

                # Coupling: generate the next step's geometry from the
                # CoupledUnsteadyProblem, then initialize its bound RingVortices
                # before shedding wake.
                if self._current_step < self.num_steps - 1:
                    self.coupled_unsteady_problem.initialize_next_problem(self)
                    self._initialize_panel_vortices_at(step + 1)

                    _logger.debug("Shedding RingVortices into the wake.")
                    self._populate_next_airplanes_wake()

                self.steady_problems_data_storage.append(
                    self._get_steady_problem_at(step)
                )
                bar.update(n=float(approx_times[step + 1]))

            _logger.debug("Calculating averaged or final forces and moments.")
            self._finalize_loads()

        if calculate_streamlines:
            _logger.debug("Calculating streamlines.")
            _functions.calculate_streamlines(self)

        self.steady_problems = tuple(self.steady_problems_data_storage)
        self.ran = True

    def initialize_step_geometry(self, step: int) -> None:
        """Initializes geometry for a specific step without solving.

        Sets up bound RingVortices and wake RingVortices for the specified time step,
        but does not solve the aerodynamic system. Use this for geometry-only analysis
        like delta_time optimization.

        This method must be called sequentially for each step starting from 0, as wake
        vortices at step N depend on the geometry from step N - 1.

        :param step: The time step to initialize geometry for. It is zero indexed. It
            must be a non-negative int and be less than the total number of steps.
        :return: None
        """
        step = _parameter_validation.int_in_range_return_int(
            step, "step", 0, True, self.num_steps, False
        )

        if step == 0:
            self._initialize_panel_vortices_at(0)

        self._current_step = step
        current_problem: problems.SteadyProblem = self._get_steady_problem_at(step)
        self.current_airplanes = current_problem.airplanes
        self.current_operating_point = current_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E

        if step < self.num_steps - 1:
            self._populate_next_airplanes_wake_vortex_points()
            self._populate_next_airplanes_wake_vortices()

    def _get_steady_problem_at(self, step: int) -> problems.SteadyProblem:
        """Gets the SteadyProblem at a given time step via the CoupledUnsteadyProblem.

        This override is the key abstraction point that enables the coupled solver. The
        parent retrieves from a pre-computed tuple; this retrieves from the
        CoupledUnsteadyProblem's dynamically growing list.

        :param step: The time step of the desired SteadyProblem. It must be between 0
            and the number of problems initialized so far, exclusive.
        :return: The SteadyProblem at the given time step.
        """
        return self.coupled_unsteady_problem.get_steady_problem(step)
