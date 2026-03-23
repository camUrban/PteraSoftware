"""Contains the CoupledUnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

CoupledUnsteadyRingVortexLatticeMethodSolver: A class used to solve CoupledUnsteadyProblems with
the unsteady ring vortex lattice method.

**Contains the following functions:**

None
"""

from __future__ import annotations

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

_logger = _logging.get_logger("unsteady_ring_vortex_lattice_method")


# TEST: Consider adding unit tests for this function.
# TEST: Assess how comprehensive this function's integration tests are and update or
#  extend them if needed.
class CoupledUnsteadyRingVortexLatticeMethodSolver(UnsteadyRingVortexLatticeMethodSolver):
    """A class used to solve UnsteadyProblems with the unsteady ring vortex lattice
    method.

    **Contains the following methods:**

    run: Runs the solver on the UnsteadyProblem.

    initialize_step_geometry: Initializes geometry for a specific step without solving.

    calculate_solution_velocity: Finds the fluid velocity (in the first Airplane's
    geometry axes, observed from the Earth frame) at one or more points (in the first
    Airplane's geometry axes, relative to the first Airplane's CG) due to the freestream
    velocity and the induced velocity from every RingVortex.
    """

    def __init__(
        self, coupled_unsteady_problem: problems.CoupledUnsteadyProblem
    ) -> None:
        """The initialization method.

        :param unsteady_problem: The UnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(coupled_unsteady_problem, problems.CoupledUnsteadyProblem):
            raise TypeError(
                "coupled_unsteady_problem must be a CoupledUnsteadyProblem."
            )
        # self.coupled_unsteady_problem must be defined before the call to super().__init__()
        # because the parent class's __init__ method calls methods that rely on self.coupled_unsteady_problem being defined.
        self.coupled_unsteady_problem = coupled_unsteady_problem
        super().__init__(problems.UnsteadyProblem(coupled_unsteady_problem.movement))

        self.num_steps = self.coupled_unsteady_problem.num_steps
        self.delta_time = self.coupled_unsteady_problem.delta_time
        self.first_results_step = self.coupled_unsteady_problem.first_results_step
        self._first_averaging_step = self.coupled_unsteady_problem.first_averaging_step

        self.steady_problems = []

        first_steady_problem: problems.SteadyProblem = (
            self.get_steady_problem_at(0)
        )

        # number of panels overide and strip leading edge point initialization
        num_panels = 0
        panel_count = 0
        self.slep_point_indices = []
        for airplane in first_steady_problem.airplanes:
            num_panels += airplane.num_panels
            for wing in airplane.wings:
                for wing_cross_section in wing.wing_cross_sections:
                    self.slep_point_indices.append(panel_count)
                    if wing_cross_section.num_spanwise_panels is not None:
                        panel_count += wing_cross_section.num_spanwise_panels
        self.slep_point_indices = np.array(self.slep_point_indices, dtype=int)
        self.num_panels: int = num_panels

        # The current time step's center bound LineVortex points for the right,
        # front, left, and back legs (in the first Airplane's geometry axes,
        # relative to the local strip leading edge point).
        self.stackCblvpr_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpf_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpl_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stackCblvpb_GP1_Slep: np.ndarray = np.empty(0, dtype=float)

        # The colocation panel points and the front left panel point (in the first Airplane's
        # geometry axes, relative to the local strip leading edge point and the first
        # Airplane's CG respectively).
        self.stackCpp_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        # Leading edge of the panel points
        self.stack_Flpp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.moments_GP1_Slep: np.ndarray = np.empty(0, dtype=float)
        self.stack_leading_edge_points: np.ndarray = np.empty(0, dtype=float)

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

        # Loop through this time step's Airplanes to create a list of their Wings.
        # Here we calculate all of our values from our first ariplane to start our main run loop
        this_problem: problems.SteadyProblem = (
            self.coupled_unsteady_problem.get_steady_problem(0)
        )
        these_airplanes = this_problem.airplanes
        num_wing_panels = 0
        these_wings: list[list[geometry.wing.Wing]] = []
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

        for step in range(self.num_steps):
            # The number of wake RingVortices is the time step number multiplied by
            # the number of spanwise Panels. This works because the first time step
            # number is 0.
            this_num_wake_ring_vortices = step * this_num_spanwise_panels

            # Allocate the ndarrays for this time step.
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

            # Append this time step's ndarrays to the lists of ndarrays.
            self.list_num_wake_vortices.append(this_num_wake_ring_vortices)
            self._list_wake_vortex_strengths.append(this_wake_ring_vortex_strengths)
            self._list_wake_vortex_ages.append(this_wake_ring_vortex_ages)
            self.listStackBrwrvp_GP1_CgP1.append(thisStackBrwrvp_GP1_CgP1)
            self.listStackFrwrvp_GP1_CgP1.append(thisStackFrwrvp_GP1_CgP1)
            self.listStackFlwrvp_GP1_CgP1.append(thisStackFlwrvp_GP1_CgP1)
            self.listStackBlwrvp_GP1_CgP1.append(thisStackBlwrvp_GP1_CgP1)

        # The following loop attempts to predict how much time each time step will
        # take, relative to the other time steps. This data will be used to generate
        # estimates of how much longer a simulation will take, and create a smoothly
        # advancing progress bar.

        # Initialize list that will hold the approximate, relative times. This has
        # one more element than the number of time steps, because I will also use the
        # progress bar during the simulation initialization.
        approx_times = np.zeros(self.num_steps + 1, dtype=float)
        if step != 0:
            # Calculate the total number of RingVortices analyzed during this step.
            num_wing_ring_vortices = num_wing_panels
            num_wake_ring_vortices = self.list_num_wake_vortices[step]
            num_ring_vortices = num_wing_ring_vortices + num_wake_ring_vortices

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
            # Initialize all the Airplanes' bound RingVortices.
            _logger.debug("Initializing all Airplanes' bound RingVortices.")

            # Update the progress bar based on the initialization step's predicted
            # approximate, relative computing time.
            bar.update(n=float(approx_times[0]))

            # Iterate through the time steps.
            for step in range(self.num_steps):

                # Save attributes to hold the current step, Airplanes,
                # and OperatingPoint, and freestream velocity (in the first
                # Airplane's geometry axes, observed from the Earth frame).
                self._current_step = step
                current_problem: problems.SteadyProblem = (
                    self.coupled_unsteady_problem.get_steady_problem(self._current_step)
                )

                # Initialize all the current step's bound RingVortices.
                _logger.debug(f"Initializing step {step}'s RingVortices")
                self._initialize_panel_vortex(current_problem, step)
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
                self._last_bound_vortex_strengths = np.zeros(
                    self.num_panels, dtype=float
                )

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

                self.stackCblvpr_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpf_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpl_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpb_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCpp_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
                self.moments_GP1_Slep = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFlpp_GP1_CgP1 = np.zeros((self.num_panels, 3), dtype=float)

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

                # Get the pre-allocated (but still all zero) arrays of wake
                # information that are associated with this time step.
                self._current_wake_vortex_strengths = self._list_wake_vortex_strengths[
                    step
                ]
                self._current_wake_vortex_ages = self._list_wake_vortex_ages[step]
                self._currentStackBrwrvp_GP1_CgP1 = self.listStackBrwrvp_GP1_CgP1[step]
                self._currentStackFrwrvp_GP1_CgP1 = self.listStackFrwrvp_GP1_CgP1[step]
                self._currentStackFlwrvp_GP1_CgP1 = self.listStackFlwrvp_GP1_CgP1[step]
                self._currentStackBlwrvp_GP1_CgP1 = self.listStackBlwrvp_GP1_CgP1[step]

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
                # solely to the wake RingVortices.
                _logger.debug("Calculating the wake Wing influences.")
                self._calculate_wake_wing_influences()

                # Solve for each bound RingVortex's strength.
                _logger.debug("Calculating bound RingVortex strengths.")
                self._calculate_vortex_strengths()

                # Solve for the forces (in the first Airplane's geometry axes) and
                # moments (in the first Airplane's geometry axes, relative to the
                # first Airplane's CG) on each Panel.
                if self._current_step >= self.first_results_step:
                    _logger.debug("Calculating forces and moments.")
                    self._calculate_loads()

                # Shed RingVortices into the wake.
                # Check if the current time step is not the last step.
                if self._current_step < self.num_steps - 1:
                    self.coupled_unsteady_problem.initialize_next_problem(self)
                    self._initialize_panel_vortex(
                        self.coupled_unsteady_problem.get_steady_problem(step + 1),
                        step + 1,
                    )
                    # Shed RingVortices into the wake.
                    _logger.debug("Shedding RingVortices into the wake.")
                    self._populate_next_airplanes_wake()

                # Update the progress bar based on this time step's predicted
                # approximate, relative computing time.
                self.steady_problems.append(
                    self.coupled_unsteady_problem.get_steady_problem(step)
                )
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

        Sets up bound RingVortices and wake RingVortices for the specified time step,
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

        # Initialize bound RingVortices for all steps on the first call.
        if step == 0:
            self._initialize_panel_vortex(
                self.coupled_unsteady_problem.get_steady_problem(0), 0
            )

        # Set the current step and related state.
        self._current_step = step
        current_problem: problems.SteadyProblem = self.get_steady_problem_at(step)
        self.current_airplanes = current_problem.airplanes
        self.current_operating_point = current_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E

        # Populate the wake for the next step (if not the last step).
        if step < self.num_steps - 1:
            self._populate_next_airplanes_wake_vortex_points()
            self._populate_next_airplanes_wake_vortices()

    def _load_calculation_moment_processing_hook(
        self,
        rightLegForces_GP1,
        frontLegForces_GP1,
        leftLegForces_GP1,
        backLegForces_GP1,
        unsteady_forces_GP1,
    ) -> np.ndarray:
        """Extends the parent implementation to also compute moments
        about the SLEP point, stored in self.moments_GP1_Slep.

        :return: moments_GP1_CgP1, a (N,3) ndarray of floats representing the moments
        (in the first Airplane's geometry axes, relative to the first Airplane's CG)
        on every Panel at the current time step."""
        # Find the moments (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) on the Panels' RingVortex's right LineVortex,
        # front LineVortex, left LineVortex, and back LineVortex.
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

        # Find the moments (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG) due to the unsteady component of the force on each Panel.
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
        """Updates the bound RingVortex position variables to be relative to the
        Airplane's SLEP points.

        :return: None
        """
        # Find the bound RingVortex leg center positions relative to the SLEP points.
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

    def get_steady_problem_at(self, step: int) -> problems.SteadyProblem:
        """Gets the SteadyProblem at a given time step. This is used for dynamic dispatch
        with coupled unsteady problem as we want to have a different way of getting the steady
        problem based on the solver type, but we want functions to work the same way regardless
        of the solver type so that we don't need ot duplicate functionality across solvers.

        :param step: An int representing the time step of the desired SteadyProblem. It
            must be between 0 and num_steps - 1, inclusive.
        :return: The SteadyProblem at the given time step.
        """
        if step < 0 or step >= self.num_steps:
            raise ValueError(
                f"Step must be between 0 and {self.num_steps - 1}, inclusive."
            )
        return self.coupled_unsteady_problem.get_steady_problem(step)
