"""Contains the UnsteadyRingVortexLatticeMethodSolver class.

**Contains the following classes:**

UnsteadyRingVortexLatticeMethodSolver: A class used to solve UnsteadyProblems with the
unsteady ring vortex lattice method.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

import numpy as np
from tqdm import tqdm

from . import geometry
from . import movements
from . import _aerodynamics
from . import _functions
from . import _logging
from . import _panel
from . import _parameter_validation
from . import operating_point
from . import problems

_logger = _logging.get_logger("unsteady_ring_vortex_lattice_method")


# TEST: Consider adding unit tests for this function.
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
    velocity and the induced velocity from every RingVortex.
    """

    def __init__(self, unsteady_problem: problems.UnsteadyProblem) -> None:
        """The initialization method.

        :param unsteady_problem: The UnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(unsteady_problem, problems.UnsteadyProblem):
            raise TypeError("unsteady_problem must be an UnsteadyProblem.")
        self.unsteady_problem = unsteady_problem

        self.num_steps = self.unsteady_problem.num_steps
        self.delta_time = self.unsteady_problem.delta_time
        self.first_results_step = self.unsteady_problem.first_results_step
        self._first_averaging_step = self.unsteady_problem.first_averaging_step
        self._current_step: int = 0
        self._prescribed_wake: bool = True

        self.steady_problems = self.unsteady_problem.steady_problems

        first_steady_problem: problems.SteadyProblem = self.steady_problems[0]

        self.current_airplanes: list[geometry.airplane.Airplane] = []
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
        self._current_bound_vortex_strengths: np.ndarray = np.empty(0, dtype=float)
        self._last_bound_vortex_strengths: np.ndarray = np.empty(0, dtype=float)

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
        # and back left bound RingVortex points (in the first Airplane's geometry
        # axes, relative to the first Airplane's CG).
        self.stackBrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackFrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackFlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.stackBlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackBrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackFrbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackFlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._lastStackBlbrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

        # The current and last time step's center bound LineVortex points for the
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

        # Right, front, left, and back bound RingVortex vectors (in the first
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
        # wake RingVortex points (in the first Airplane's geometry axes, relative to
        # the first Airplane's CG).
        self._currentStackBrwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._currentStackFrwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._currentStackFlwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self._currentStackBlwrvp_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)

        # Initialize lists to store aerodynamic data about the wake at each time
        # step. These attributes are used by the output module to animate the wake.
        self.list_num_wake_vortices: list[int] = []
        # TODO: Determine if these private attributes are needed and if not
        #  delete them.
        self._list_wake_vortex_strengths: list[np.ndarray] = []
        self._list_wake_vortex_ages: list[np.ndarray] = []
        self.listStackBrwrvp_GP1_CgP1: list[np.ndarray] = []
        self.listStackFrwrvp_GP1_CgP1: list[np.ndarray] = []
        self.listStackFlwrvp_GP1_CgP1: list[np.ndarray] = []
        self.listStackBlwrvp_GP1_CgP1: list[np.ndarray] = []

        self.stackSeedPoints_GP1_CgP1: np.ndarray = np.empty(0, dtype=float)
        self.gridStreamlinePoints_GP1_CgP1: np.ndarray = np.empty((0, 3), dtype=float)

        self.ran = False

    def run(
        self,
        prescribed_wake: bool | np.bool_ = True,
        calculate_streamlines: bool | np.bool_ = True,
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
        :return: None
        """
        self._prescribed_wake = _parameter_validation.boolLike_return_bool(
            prescribed_wake, "prescribed_wake"
        )
        calculate_streamlines = _parameter_validation.boolLike_return_bool(
            calculate_streamlines, "calculate_streamlines"
        )

        # The following loop iterates through the time steps to populate currently
        # empty attributes with lists of pre-allocated arrays. During the simulation,
        # these arrays will be filled with data that describe the wake. Using this
        # method eliminates the need for computationally expensive on-the-fly
        # allocation and object copying.
        for step in range(self.num_steps):
            this_problem: problems.SteadyProblem = self.steady_problems[step]
            these_airplanes = this_problem.airplanes

            # Loop through this time step's Airplanes to create a list of their Wings.
            these_wings: list[list[geometry.wing.Wing]] = []
            for airplane in these_airplanes:
                these_wings.append(airplane.wings)

            # Iterate through the Wings to get the total number of spanwise Panels.
            this_num_spanwise_panels = 0
            for this_wing_set in these_wings:
                for this_wing in this_wing_set:
                    _this_wing_num_spanwise_panels = this_wing.num_spanwise_panels
                    assert _this_wing_num_spanwise_panels is not None

                    this_num_spanwise_panels += _this_wing_num_spanwise_panels

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
        for step in range(1, self.num_steps):
            this_problem = self.steady_problems[step]
            these_airplanes = this_problem.airplanes

            # Iterate through this time step's Airplanes to get the total number of
            # Wing Panels.
            num_wing_panels = 0
            for airplane in these_airplanes:
                num_wing_panels += airplane.num_panels

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
            disable=False,
            bar_format="{desc}:{percentage:3.0f}% |{bar}| Elapsed: {elapsed}, "
            "Remaining: {remaining}",
        ) as bar:
            # Initialize all the Airplanes' bound RingVortices.
            _logger.info("Initializing all Airplanes' bound RingVortices.")
            self._initialize_panel_vortices()

            # Update the progress bar based on the initialization step's predicted
            # approximate, relative computing time.
            bar.update(n=float(approx_times[0]))

            # Iterate through the time steps.
            for step in range(self.num_steps):

                # Save attributes to hold the current step, Airplanes,
                # and OperatingPoint, and freestream velocity (in the first
                # Airplane's geometry axes, observed from the Earth frame).
                self._current_step = step
                current_problem: problems.SteadyProblem = self.steady_problems[
                    self._current_step
                ]
                self.current_airplanes = current_problem.airplanes
                self.current_operating_point = current_problem.operating_point
                self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E
                _logger.info(
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
                _logger.info("Collapsing the geometry.")
                self._collapse_geometry()

                # Find the matrix of Wing Wing influence coefficients associated with
                # the Airplanes' geometries at this time step.
                _logger.info("Calculating the Wing Wing influences.")
                self._calculate_wing_wing_influences()

                # Find the normal velocity (in the first Airplane's geometry axes,
                # observed from the Earth frame) at every collocation point due
                # solely to the freestream.
                _logger.info("Calculating the freestream Wing influences.")
                self._calculate_freestream_wing_influences()

                # Find the normal velocity (in the first Airplane's geometry axes,
                # observed from the Earth frame) at every collocation point due
                # solely to the wake RingVortices.
                _logger.info("Calculating the wake Wing influences.")
                self._calculate_wake_wing_influences()

                # Solve for each bound RingVortex's strength.
                _logger.info("Calculating bound RingVortex strengths.")
                self._calculate_vortex_strengths()

                # Solve for the forces (in the first Airplane's geometry axes) and
                # moments (in the first Airplane's geometry axes, relative to the
                # first Airplane's CG) on each Panel.
                if self._current_step >= self.first_results_step:
                    _logger.info("Calculating forces and moments.")
                    self._calculate_loads()

                # Shed RingVortices into the wake.
                _logger.info("Shedding RingVortices into the wake.")
                self._populate_next_airplanes_wake()

                # Update the progress bar based on this time step's predicted
                # approximate, relative computing time.
                bar.update(n=float(approx_times[step + 1]))

            _logger.info("Calculating averaged or final forces and moments.")
            self._finalize_loads()

        # Solve for the location of the streamlines coming off the Wings' trailing
        # edges, if requested.
        if calculate_streamlines:
            _logger.info("Calculating streamlines.")
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
            self._initialize_panel_vortices()

        # Set the current step and related state.
        self._current_step = step
        current_problem: problems.SteadyProblem = self.steady_problems[step]
        self.current_airplanes = current_problem.airplanes
        self.current_operating_point = current_problem.operating_point
        self._currentVInf_GP1__E = self.current_operating_point.vInf_GP1__E

        # Populate the wake for the next step (if not the last step).
        if step < self.num_steps - 1:
            self._populate_next_airplanes_wake_vortex_points()
            self._populate_next_airplanes_wake_vortices()

    def _initialize_panel_vortices(self) -> None:
        """Calculates the locations of the bound RingVortex vertices, and then
        initializes them.

        Every Panel has a RingVortex, which is a quadrangle whose front leg is a
        LineVortex at the Panel's quarter chord. The left and right legs are
        LineVortices running along the Panel's left and right legs. If the Panel is not
        along the trailing edge, they extend backwards and meet the back LineVortex, at
        the rear Panel's quarter chord. Otherwise, they extend backwards and meet the
        back LineVortex one quarter chord back from the Panel's back leg.

        :return: None
        """
        for steady_problem_id, steady_problem in enumerate(self.steady_problems):
            # Find the freestream velocity (in the first Airplane's geometry axes,
            # observed from the Earth frame) at this time step.
            this_operating_point = steady_problem.operating_point
            vInf_GP1__E = this_operating_point.vInf_GP1__E

            # Iterate through this SteadyProblem's Airplanes' Wings.
            for airplane_id, airplane in enumerate(steady_problem.airplanes):
                for wing_id, wing in enumerate(airplane.wings):
                    _num_spanwise_panels = wing.num_spanwise_panels
                    assert _num_spanwise_panels is not None

                    # Iterate through the Wing's chordwise and spanwise positions.
                    for chordwise_position in range(wing.num_chordwise_panels):
                        for spanwise_position in range(_num_spanwise_panels):
                            _panels = wing.panels
                            assert _panels is not None

                            # Pull the Panel out of the Wing's 2D ndarray of Panels.
                            panel: _panel.Panel = _panels[
                                chordwise_position, spanwise_position
                            ]

                            _Flbvp_GP1_CgP1 = panel.Flbvp_GP1_CgP1
                            assert _Flbvp_GP1_CgP1 is not None

                            _Frbvp_GP1_CgP1 = panel.Frbvp_GP1_CgP1
                            assert _Frbvp_GP1_CgP1 is not None

                            # Find the location of this Panel's front left and
                            # front right RingVortex points (in the first Airplane's
                            # geometry axes, relative to the first Airplane's CG).
                            Flrvp_GP1_CgP1 = _Flbvp_GP1_CgP1
                            Frrvp_GP1_CgP1 = _Frbvp_GP1_CgP1

                            # Define the location of the back left and back right
                            # RingVortex points based on whether the Panel is along
                            # the trailing edge or not.
                            if not panel.is_trailing_edge:
                                next_chordwise_panel: _panel.Panel = _panels[
                                    chordwise_position + 1, spanwise_position
                                ]

                                _nextFlbvp_GP1_CgP1 = (
                                    next_chordwise_panel.Flbvp_GP1_CgP1
                                )
                                assert _nextFlbvp_GP1_CgP1 is not None

                                _nextFrbvp_GP1_CgP1 = (
                                    next_chordwise_panel.Frbvp_GP1_CgP1
                                )
                                assert _nextFrbvp_GP1_CgP1 is not None

                                Blrvp_GP1_CgP1 = _nextFlbvp_GP1_CgP1
                                Brrvp_GP1_CgP1 = _nextFrbvp_GP1_CgP1
                            else:
                                # As these vertices are directly behind the trailing
                                # edge, they are spaced back from their Panel's
                                # vertex by one quarter of the distance traveled by
                                # the trailing edge during a time step. This is to
                                # more accurately predict drag. More information can
                                # be found on pages 37-39 of "Modeling of aerodynamic
                                # forces in flapping flight with the Unsteady Vortex
                                # Lattice Method" by Thomas Lambert.
                                if steady_problem_id == 0:
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
                                    last_steady_problem = self.steady_problems[
                                        steady_problem_id - 1
                                    ]
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

                                    # We subtract (thisBlpp_GP1_CgP1 -
                                    # lastBlpp_GP1_CgP1) / self.delta_time from
                                    # vInf_GP1__E, because we want the apparent fluid
                                    # velocity due to motion (observed in the Earth
                                    # frame, in the first Airplane's geometry axes).
                                    # This is the vector pointing opposite the
                                    # velocity from motion.
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

                                    # The comment from above about apparent fluid
                                    # velocity due to motion applies here as well.
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

                            # Initialize the Panel's RingVortex.
                            panel.ring_vortex = _aerodynamics.RingVortex(
                                Flrvp_GP1_CgP1=Flrvp_GP1_CgP1,
                                Frrvp_GP1_CgP1=Frrvp_GP1_CgP1,
                                Blrvp_GP1_CgP1=Blrvp_GP1_CgP1,
                                Brrvp_GP1_CgP1=Brrvp_GP1_CgP1,
                                strength=1.0,
                            )

    def _collapse_geometry(self) -> None:
        """Converts attributes of the UnsteadyProblem's geometry into 1D ndarrays.

        This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize variables to hold the global position of the Panel and the wake
        # RingVortex as we iterate through them.
        global_panel_position = 0
        global_wake_ring_vortex_position = 0

        # Iterate through the current time step's Airplanes' Wings.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:
                _panels = wing.panels
                assert _panels is not None

                _wake_ring_vortices = wing.wake_ring_vortices
                assert _wake_ring_vortices is not None

                # Convert this Wing's 2D ndarray of Panels and wake RingVortices into
                # 1D ndarrays.
                panels = np.ravel(_panels)
                wake_ring_vortices = np.ravel(_wake_ring_vortices)

                # Iterate through the 1D ndarray of this Wing's Panels.
                panel: _panel.Panel
                for panel in panels:
                    # Update the solver's list of attributes with this Panel's
                    # attributes.
                    _functions.update_ring_vortex_solvers_panel_attributes(
                        ring_vortex_solver=self,
                        global_panel_position=global_panel_position,
                        panel=panel,
                    )

                    # Increment the global Panel position variable.
                    global_panel_position += 1

                # Iterate through the 1D ndarray of this Wing's wake RingVortices.
                wake_ring_vortex: _aerodynamics.RingVortex
                for wake_ring_vortex in wake_ring_vortices:
                    # Update the solver's list of attributes with this wake
                    # RingVortex's attributes.
                    self._current_wake_vortex_strengths[
                        global_wake_ring_vortex_position
                    ] = wake_ring_vortex.strength
                    self._current_wake_vortex_ages[global_wake_ring_vortex_position] = (
                        wake_ring_vortex.age
                    )
                    self._currentStackFrwrvp_GP1_CgP1[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Frrvp_GP1_CgP1
                    self._currentStackFlwrvp_GP1_CgP1[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Flrvp_GP1_CgP1
                    self._currentStackBlwrvp_GP1_CgP1[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Blrvp_GP1_CgP1
                    self._currentStackBrwrvp_GP1_CgP1[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Brrvp_GP1_CgP1

                    # Increment the global wake RingVortex position variable.
                    global_wake_ring_vortex_position += 1

        if self._current_step > 0:

            # Reset the global Panel position variable.
            global_panel_position = 0

            last_problem = self.steady_problems[self._current_step - 1]
            last_airplanes = last_problem.airplanes

            # Iterate through the last time step's Airplanes' Wings.
            for last_airplane in last_airplanes:
                for last_wing in last_airplane.wings:
                    _last_panels = last_wing.panels
                    assert _last_panels is not None

                    # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                    last_panels = np.ravel(_last_panels)

                    # Iterate through the 1D ndarray of this Wing's Panels.
                    last_panel: _panel.Panel
                    for last_panel in last_panels:
                        # Update the solver's list of attributes with this Panel's
                        # attributes.
                        self._stackLastCpp_GP1_CgP1[global_panel_position, :] = (
                            last_panel.Cpp_GP1_CgP1
                        )

                        last_ring_vortex = last_panel.ring_vortex
                        assert last_ring_vortex is not None

                        self._last_bound_vortex_strengths[global_panel_position] = (
                            last_ring_vortex.strength
                        )

                        # TODO: Test if we can replace the calls to LineVortex
                        #  attributes with calls to RingVortex attributes.
                        self._lastStackBrbrvp_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.right_leg.Slvp_GP1_CgP1
                        )
                        self._lastStackFrbrvp_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.right_leg.Elvp_GP1_CgP1
                        )
                        self._lastStackFlbrvp_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.left_leg.Slvp_GP1_CgP1
                        )
                        self._lastStackBlbrvp_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.left_leg.Elvp_GP1_CgP1
                        )
                        self._lastStackCblvpr_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.right_leg.Clvp_GP1_CgP1
                        )
                        self._lastStackCblvpf_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.front_leg.Clvp_GP1_CgP1
                        )
                        self._lastStackCblvpl_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.left_leg.Clvp_GP1_CgP1
                        )
                        self._lastStackCblvpb_GP1_CgP1[global_panel_position, :] = (
                            last_ring_vortex.back_leg.Clvp_GP1_CgP1
                        )

                        # Increment the global Panel position variable.
                        global_panel_position += 1

    def _calculate_wing_wing_influences(self) -> None:
        """Finds the current time step's SteadyProblem's 2D ndarray of Wing Wing
        influence coefficients (observed from the Earth frame).

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in the first Airplane's
        # geometry axes, observed from the Earth frame) induced at each Panel's
        # collocation point by each bound RingVortex. The answer is normalized
        # because the solver's list of bound RingVortex strengths was initialized to
        # all be 1.0. This will be updated once the correct strengths are calculated.
        gridNormVIndCpp_GP1_E = _aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
            stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
            stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
            stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
            stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
            strengths=self._current_bound_vortex_strengths,
            ages=None,
            nu=self.current_operating_point.nu,
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

        **Notes:**

        If the current time step is the first time step, no wake has been shed, so this
        method will return zero for all the wake Wing influence coefficients (observed
        from the Earth frame).

        :return: None
        """
        if self._current_step > 0:
            # Get the velocities (in the first Airplane's geometry axes, observed
            # from the Earth frame) induced by the wake RingVortices at each Panel's
            # collocation point.
            currentStackWakeV_GP1_E = (
                _aerodynamics.collapsed_velocities_from_ring_vortices(
                    stackP_GP1_CgP1=self.stackCpp_GP1_CgP1,
                    stackBrrvp_GP1_CgP1=self._currentStackBrwrvp_GP1_CgP1,
                    stackFrrvp_GP1_CgP1=self._currentStackFrwrvp_GP1_CgP1,
                    stackFlrvp_GP1_CgP1=self._currentStackFlwrvp_GP1_CgP1,
                    stackBlrvp_GP1_CgP1=self._currentStackBlwrvp_GP1_CgP1,
                    strengths=self._current_wake_vortex_strengths,
                    ages=self._current_wake_vortex_ages,
                    nu=self.current_operating_point.nu,
                )
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
            # RingVortices have been shed.
            self._currentStackWakeWingInfluences__E = np.zeros(
                self.num_panels, dtype=float
            )

    def _calculate_vortex_strengths(self) -> None:
        """Solves for the strength of each Panel's bound RingVortex.

        :return: None
        """
        self._current_bound_vortex_strengths = np.linalg.solve(
            self._currentGridWingWingInfluences__E,
            -self._currentStackWakeWingInfluences__E
            - self._currentStackFreestreamWingInfluences__E,
        )

        # Update the bound RingVortices' strengths.
        _panels = self.panels
        assert _panels is not None
        for panel_num in range(_panels.size):
            panel: _panel.Panel = _panels[panel_num]

            this_ring_vortex = panel.ring_vortex
            assert this_ring_vortex is not None

            this_ring_vortex.update_strength(
                self._current_bound_vortex_strengths[panel_num]
            )

    def calculate_solution_velocity(
        self, stackP_GP1_CgP1: np.ndarray | Sequence[Sequence[float | int]]
    ) -> np.ndarray:
        """Finds the fluid velocity (in the first Airplane's geometry axes, observed
        from the Earth frame) at one or more points (in the first Airplane's geometry
        axes, relative to the first Airplane's CG) due to the freestream velocity and
        the induced velocity from every RingVortex.

        **Notes:**

        This method assumes that the correct strengths for the RingVortices have already
        been calculated and set.

        This method also does not include the velocity due to the Movement's motion at
        any of the points provided, as it has no way of knowing if any of the points lie
        on panels.

        :param stackP_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (N,3) representing the positions of the evaluation points (in the
            first Airplane's geometry axes, relative to the first Airplane's CG). Can be
            a tuple, list, or ndarray. Values are converted to floats internally. The
            units are in meters.
        :return: A (N,3) ndarray of floats representing the velocity (in the first
            Airplane's geometry axes, observed from the Earth frame) at each evaluation
            point due to the summed effects of the freestream velocity and the induced
            velocity from every RingVortex and HorseshoeVortex. The units are in meters
            per second.
        """
        stackP_GP1_CgP1 = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_GP1_CgP1, "stackP_GP1_CgP1"
            )
        )

        stackBoundRingVInd_GP1_E = (
            _aerodynamics.collapsed_velocities_from_ring_vortices(
                stackP_GP1_CgP1=stackP_GP1_CgP1,
                stackBrrvp_GP1_CgP1=self.stackBrbrvp_GP1_CgP1,
                stackFrrvp_GP1_CgP1=self.stackFrbrvp_GP1_CgP1,
                stackFlrvp_GP1_CgP1=self.stackFlbrvp_GP1_CgP1,
                stackBlrvp_GP1_CgP1=self.stackBlbrvp_GP1_CgP1,
                strengths=self._current_bound_vortex_strengths,
                ages=None,
                nu=self.current_operating_point.nu,
            )
        )
        stackWakeRingVInd_GP1_E = _aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_GP1_CgP1=stackP_GP1_CgP1,
            stackBrrvp_GP1_CgP1=self._currentStackBrwrvp_GP1_CgP1,
            stackFrrvp_GP1_CgP1=self._currentStackFrwrvp_GP1_CgP1,
            stackFlrvp_GP1_CgP1=self._currentStackFlwrvp_GP1_CgP1,
            stackBlrvp_GP1_CgP1=self._currentStackBlwrvp_GP1_CgP1,
            strengths=self._current_wake_vortex_strengths,
            ages=self._current_wake_vortex_ages,
            nu=self.current_operating_point.nu,
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

        This method assumes that the correct strengths for the RingVortices and
        HorseshoeVortices have already been calculated and set.

        This method used to accidentally double-count the load on each Panel due to the
        left and right LineVortex legs. Additionally, it didn't include contributions to
        the load on each Panel from their back LineVortex legs. Thankfully, these issues
        only introduced small errors in most typical simulations. They have both now
        been fixed by (1) using a 1/2 factor for each "effective" vortex strength shared
        between two Panels, and (2) including the effects each Panel's back LineVortex
        with its own effective strength.

        :return: None
        """
        # Initialize a variable to hold the global Panel position as we iterate
        # through them.
        global_panel_position = 0

        # Initialize three 1D ndarrays to hold the effective strength of the Panels'
        # RingVortices' LineVortices.
        effective_right_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_front_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_left_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)
        effective_back_line_vortex_strengths = np.zeros(self.num_panels, dtype=float)

        # Iterate through the Airplanes' Wings.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:
                _panels = wing.panels
                assert _panels is not None

                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(_panels)

                # Iterate through this Wing's 1D ndarray of Panels.
                panel: _panel.Panel
                for panel in panels:
                    _local_chordwise_position = panel.local_chordwise_position
                    assert _local_chordwise_position is not None

                    _local_spanwise_position = panel.local_spanwise_position
                    assert _local_spanwise_position is not None

                    if panel.is_right_edge:
                        # Set the effective right LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_right: _panel.Panel = _panels[
                            _local_chordwise_position,
                            _local_spanwise_position + 1,
                        ]

                        ring_vortex_to_right = panel_to_right.ring_vortex
                        assert ring_vortex_to_right is not None

                        # Set the effective right LineVortex strength to 1/2 the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the right.
                        effective_right_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - ring_vortex_to_right.strength
                        ) / 2

                    if panel.is_leading_edge:
                        # Set the effective front LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_front: _panel.Panel = _panels[
                            _local_chordwise_position - 1,
                            _local_spanwise_position,
                        ]

                        ring_vortex_to_front = panel_to_front.ring_vortex
                        assert ring_vortex_to_front is not None

                        # Set the effective front LineVortex strength to 1/2 the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel in front of it.
                        effective_front_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - ring_vortex_to_front.strength
                        ) / 2

                    if panel.is_left_edge:
                        # Set the effective left LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_left: _panel.Panel = _panels[
                            _local_chordwise_position,
                            _local_spanwise_position - 1,
                        ]

                        ring_vortex_to_left = panel_to_left.ring_vortex
                        assert ring_vortex_to_left is not None

                        # Set the effective left LineVortex strength to 1/2 the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the left.
                        effective_left_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - ring_vortex_to_left.strength
                        ) / 2

                    if panel.is_trailing_edge:
                        if self._current_step == 0:
                            # Set the effective back LineVortex strength to this
                            # Panel's RingVortex's strength, as, for the first time
                            # step, there isn't a wake RingVortex to cancel it out.
                            effective_back_line_vortex_strengths[
                                global_panel_position
                            ] = self._current_bound_vortex_strengths[
                                global_panel_position
                            ]
                        else:
                            # Set the effective back LineVortex strength to the
                            # difference between this Panel's RingVortex's strength and
                            # its strength at the last time step. This models the effect
                            # of the Panel's back LineVortex being partially cancelled
                            # out by the front LineVortex of the wake RingVortex
                            # immediately to this Panel's rear. This works because that
                            # wake RingVortex has the same strength this Panel's
                            # RingVortex had last time step.
                            effective_back_line_vortex_strengths[
                                global_panel_position
                            ] = (
                                self._current_bound_vortex_strengths[
                                    global_panel_position
                                ]
                                - self._last_bound_vortex_strengths[
                                    global_panel_position
                                ]
                            )
                    else:
                        panel_to_back: _panel.Panel = _panels[
                            _local_chordwise_position + 1,
                            _local_spanwise_position,
                        ]

                        _ring_vortex_to_back = panel_to_back.ring_vortex
                        assert _ring_vortex_to_back is not None

                        # Set the effective back LineVortex strength to 1/2 the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the back.
                        effective_back_line_vortex_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - _ring_vortex_to_back.strength
                        ) / 2

                    # Increment the global Panel position variable.
                    global_panel_position += 1

        # Calculate the velocity (in the first Airplane's geometry axes, observed
        # from the Earth frame) at the center of every Panels' RingVortex's right
        # LineVortex, front LineVortex, left LineVortex, and back LineVortex.
        stackVelocityRightLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(stackP_GP1_CgP1=self.stackCblvpr_GP1_CgP1)
            + self._calculate_current_movement_velocities_at_right_leg_centers()
        )
        stackVelocityFrontLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(stackP_GP1_CgP1=self.stackCblvpf_GP1_CgP1)
            + self._calculate_current_movement_velocities_at_front_leg_centers()
        )
        stackVelocityLeftLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(stackP_GP1_CgP1=self.stackCblvpl_GP1_CgP1)
            + self._calculate_current_movement_velocities_at_left_leg_centers()
        )
        stackVelocityBackLineVortexCenters_GP1__E = (
            self.calculate_solution_velocity(stackP_GP1_CgP1=self.stackCblvpb_GP1_CgP1)
            + self._calculate_current_movement_velocities_at_back_leg_centers()
        )

        # Using the effective LineVortex strengths and the Kutta-Joukowski theorem,
        # find the forces (in the first Airplane's geometry axes) on the Panels'
        # RingVortex's right LineVortex, front LineVortex, left LineVortex, and back
        # LineVortex using the effective vortex strengths.
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
        # Ptera Software defines RingVortices with counter-clockwise (CCW) vertex
        # ordering, while the references use clockwise (CW) ordering. Both define panel
        # normals as pointing upward. This convention difference only affects the
        # unsteady force term because it depends on both vortex strength and the normal
        # vector. When converting from CCW to CW, the strength changes sign but the
        # normal vector does not, requiring a sign correction. In contrast, steady
        # Kutta-Joukowski forces depend on the strength and the LineVortex vectors. Both
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
        # first Airplane's CG) on the Panels' RingVortex's right LineVortex,
        # front LineVortex, left LineVortex, and back LineVortex.
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
        # unsteady force acts on the bound RingVortex, whose center is at the
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
        # Populate the locations of the next time step's Airplanes' wake RingVortex
        # points.
        self._populate_next_airplanes_wake_vortex_points()

        # Populate the locations of the next time step's Airplanes' wake RingVortices.
        self._populate_next_airplanes_wake_vortices()

    def _populate_next_airplanes_wake_vortex_points(self) -> None:
        """Populates the locations of the next time step's Airplanes' wake RingVortex
        points.

        **Notes:**

        This method is not vectorized but its loops only consume 1.1% of the runtime, so
        I have kept it as is for increased readability.

        :return: None
        """
        # Check that this isn't the last time step.
        if self._current_step < self.num_steps - 1:

            # Get the next time step's Airplanes.
            next_problem: problems.SteadyProblem = self.steady_problems[
                self._current_step + 1
            ]
            next_airplanes = next_problem.airplanes

            # Get the current Airplanes' combined number of Wings.
            num_wings = 0
            for airplane in self.current_airplanes:
                num_wings += len(airplane.wings)

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
                        # wake RingVortices (in the first Airplane's geometry axes,
                        # relative to the first Airplane's CG).
                        newRowWrvp_GP1_CgP1 = np.zeros(
                            (1, num_spanwise_panels + 1, 3), dtype=float
                        )

                        # Iterate through the spanwise Panel positions.
                        for spanwise_panel_id in range(num_spanwise_panels):
                            _next_panels = next_wing.panels
                            assert _next_panels is not None

                            # Get the next time step's Wing's Panel at this location.
                            next_panel: _panel.Panel = _next_panels[
                                chordwise_panel_id, spanwise_panel_id
                            ]

                            # The position of the new front left wake RingVortex's
                            # point is the next time step's Panel's bound
                            # RingVortex's back left point.
                            next_ring_vortex = next_panel.ring_vortex
                            assert next_ring_vortex is not None

                            newFlwrvp_GP1_CgP1 = next_ring_vortex.Blrvp_GP1_CgP1

                            # Add this to the row of new wake RingVortex points.
                            newRowWrvp_GP1_CgP1[0, spanwise_panel_id] = (
                                newFlwrvp_GP1_CgP1
                            )

                            # If the Panel is at the right edge of the Wing, add its
                            # back right bound RingVortex point to the row of new
                            # wake RingVortex points.
                            if spanwise_panel_id == (num_spanwise_panels - 1):
                                newRowWrvp_GP1_CgP1[0, spanwise_panel_id + 1] = (
                                    next_ring_vortex.Brrvp_GP1_CgP1
                                )

                        # Set the next time step's Wing's grid of wake RingVortex
                        # points to a copy of the row of new wake RingVortex points.
                        # This is correct because it is currently the first time step.
                        next_wing.gridWrvp_GP1_CgP1 = np.copy(newRowWrvp_GP1_CgP1)

                        # Initialize variables to hold the number of spanwise wake
                        # RingVortex points.
                        num_spanwise_points = num_spanwise_panels + 1

                        # Initialize a new ndarray to hold the second new row of wake
                        # RingVortex points (in the first Airplane's geometry axes,
                        # relative to the first Airplane's CG).
                        secondNewRowWrvp_GP1_CgP1 = np.zeros(
                            (1, num_spanwise_panels + 1, 3), dtype=float
                        )

                        # Iterate through the spanwise points.
                        for spanwise_point_id in range(num_spanwise_points):
                            # Get the corresponding point from the first row.
                            Wrvp_GP1_CgP1 = next_wing.gridWrvp_GP1_CgP1[
                                0, spanwise_point_id
                            ]
                            assert Wrvp_GP1_CgP1 is not None

                            # If the wake is prescribed, set the velocity at this
                            # point to the freestream velocity (in the first
                            # Airplane's geometry axes, observed from the Earth
                            # frame). Otherwise, set the velocity to the solution
                            # velocity at this point (in the first Airplane's
                            # geometry axes, observed from the Earth frame).
                            if self._prescribed_wake:
                                vWrvp_GP1__E = self._currentVInf_GP1__E
                            else:
                                vWrvp_GP1__E = self.calculate_solution_velocity(
                                    np.expand_dims(Wrvp_GP1_CgP1, axis=0)
                                )

                            # Update the second new row with the interpolated
                            # position of the first point.
                            secondNewRowWrvp_GP1_CgP1[0, spanwise_point_id] = (
                                Wrvp_GP1_CgP1 + vWrvp_GP1__E * self.delta_time
                            )

                        # Update the next time step's Wing's grid of wake RingVortex
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

                        # Set the next time step's Wing's grid of wake RingVortex
                        # points to a copy of this time step's Wing's grid of wake
                        # RingVortex points.
                        next_wing.gridWrvp_GP1_CgP1 = np.copy(_thisGridWrvp_GP1_CgP1)

                        # Get the number of chordwise and spanwise points.
                        num_chordwise_points = next_wing.gridWrvp_GP1_CgP1.shape[0]
                        num_spanwise_points = next_wing.gridWrvp_GP1_CgP1.shape[1]

                        # Iterate through the chordwise and spanwise point positions.
                        for chordwise_point_id in range(num_chordwise_points):
                            for spanwise_point_id in range(num_spanwise_points):
                                # Get the wake RingVortex point at this position.
                                Wrvp_GP1_CgP1 = next_wing.gridWrvp_GP1_CgP1[
                                    chordwise_point_id,
                                    spanwise_point_id,
                                ]

                                # If the wake is prescribed, set the velocity at this
                                # point to the freestream velocity (in the first
                                # Airplane's geometry axes, observed from the Earth
                                # frame). Otherwise, set the velocity to the solution
                                # velocity at this point (in the first Airplane's
                                # geometry axes, observed from the Earth frame).
                                if self._prescribed_wake:
                                    vWrvp_GP1__E = self._currentVInf_GP1__E
                                else:
                                    vWrvp_GP1__E = np.squeeze(
                                        self.calculate_solution_velocity(
                                            np.expand_dims(Wrvp_GP1_CgP1, axis=0)
                                        )
                                    )

                                # Update this point with its interpolated position.
                                next_wing.gridWrvp_GP1_CgP1[
                                    chordwise_point_id, spanwise_point_id
                                ] += (vWrvp_GP1__E * self.delta_time)

                        # Find the chordwise position of the Wing's trailing edge.
                        chordwise_panel_id = this_wing.num_chordwise_panels - 1

                        _num_spanwise_panels = this_wing.num_spanwise_panels
                        assert _num_spanwise_panels is not None

                        # Initialize a new ndarray to hold the new row of wake
                        # RingVortex vertices.
                        newRowWrvp_GP1_CgP1 = np.zeros(
                            (1, _num_spanwise_panels + 1, 3), dtype=float
                        )

                        # Iterate spanwise through the trailing edge Panels.
                        for spanwise_panel_id in range(_num_spanwise_panels):
                            _next_panels = next_wing.panels
                            assert _next_panels is not None

                            # Get the Panel at this location on the next time step's
                            # Airplane's Wing.
                            this_next_panel: _panel.Panel = _next_panels[
                                chordwise_panel_id, spanwise_panel_id
                            ]

                            # Add the Panel's back left bound RingVortex point to the
                            # grid of new wake RingVortex points.
                            next_ring_vortex = this_next_panel.ring_vortex
                            assert next_ring_vortex is not None

                            newRowWrvp_GP1_CgP1[0, spanwise_panel_id] = (
                                next_ring_vortex.Blrvp_GP1_CgP1
                            )

                            # If the Panel is at the right edge of the Wing, add its
                            # back right bound RingVortex point to the grid of new
                            # wake RingVortex vertices.
                            if spanwise_panel_id == (_num_spanwise_panels - 1):
                                newRowWrvp_GP1_CgP1[0, spanwise_panel_id + 1] = (
                                    next_ring_vortex.Brrvp_GP1_CgP1
                                )

                        # Stack the new row of wake RingVortex points above the
                        # Wing's grid of wake RingVortex points.
                        next_wing.gridWrvp_GP1_CgP1 = np.vstack(
                            (
                                newRowWrvp_GP1_CgP1,
                                next_wing.gridWrvp_GP1_CgP1,
                            )
                        )

    def _populate_next_airplanes_wake_vortices(self) -> None:
        """Populates the locations and strengths of the next time step's wake
        RingVortices.

        **Notes:**

        This method is not vectorized but its loops only consume 0.4% of the runtime, so
        I have kept it as is for increased readability.

        :return: None
        """
        # Check if the current time step is not the last step.
        if self._current_step < self.num_steps - 1:

            # Get the next time step's Airplanes.
            next_problem = self.steady_problems[self._current_step + 1]
            next_airplanes = next_problem.airplanes

            # Iterate through the next time step's Airplanes.
            for airplane_id, next_airplane in enumerate(next_airplanes):

                # For a given Airplane in the next time step, iterate through its
                # predecessor's Wings.
                for wing_id, this_wing in enumerate(
                    self.current_airplanes[airplane_id].wings
                ):
                    next_wing = next_airplane.wings[wing_id]

                    # Get the next time step's Wing's grid of wake RingVortex points.
                    nextGridWrvp_GP1_CgP1 = next_wing.gridWrvp_GP1_CgP1
                    assert nextGridWrvp_GP1_CgP1 is not None

                    # Find the number of chordwise and spanwise points in the next
                    # Wing's grid of wake RingVortex points.
                    num_chordwise_points = nextGridWrvp_GP1_CgP1.shape[0]
                    num_spanwise_points = nextGridWrvp_GP1_CgP1.shape[1]

                    this_wing_wake_ring_vortices = (
                        self.current_airplanes[airplane_id]
                        .wings[wing_id]
                        .wake_ring_vortices
                    )
                    assert this_wing_wake_ring_vortices is not None

                    # Initialize a new ndarray to hold the new row of wake RingVortices.
                    new_row_of_wake_ring_vortices = np.empty(
                        (1, num_spanwise_points - 1), dtype=object
                    )

                    # Create a new ndarray by stacking the new row of wake
                    # RingVortices on top of the current Wing's grid of wake
                    # RingVortices and assign it to the next time step's Wing.
                    next_wing.wake_ring_vortices = np.vstack(
                        (new_row_of_wake_ring_vortices, this_wing_wake_ring_vortices)
                    )

                    # Iterate through the wake RingVortex point positions.
                    for chordwise_point_id in range(num_chordwise_points):
                        for spanwise_point_id in range(num_spanwise_points):
                            # Set bools to determine if this point is on the right
                            # and/or trailing edge of the wake.
                            has_point_to_right = (
                                spanwise_point_id + 1
                            ) < num_spanwise_points
                            has_point_behind = (
                                chordwise_point_id + 1
                            ) < num_chordwise_points

                            if has_point_to_right and has_point_behind:
                                # If this point isn't on the right or trailing edge
                                # of the wake, get the four points that will be
                                # associated with the corresponding RingVortex at
                                # this position (in the first Airplane's geometry
                                # axes, relative to the first Airplane's CG),
                                # for the next time step.
                                Flwrvp_GP1_CgP1 = nextGridWrvp_GP1_CgP1[
                                    chordwise_point_id, spanwise_point_id
                                ]
                                Frwrvp_GP1_CgP1 = nextGridWrvp_GP1_CgP1[
                                    chordwise_point_id,
                                    spanwise_point_id + 1,
                                ]
                                Blwrvp_GP1_CgP1 = nextGridWrvp_GP1_CgP1[
                                    chordwise_point_id + 1,
                                    spanwise_point_id,
                                ]
                                Brwrvp_GP1_CgP1 = nextGridWrvp_GP1_CgP1[
                                    chordwise_point_id + 1,
                                    spanwise_point_id + 1,
                                ]

                                if chordwise_point_id > 0:
                                    # If this isn't the front of the wake, update the
                                    # position of the wake RingVortex at this
                                    # location for the next time step.
                                    next_wake_ring_vortices = (
                                        next_wing.wake_ring_vortices
                                    )
                                    assert next_wake_ring_vortices is not None
                                    next_wake_ring_vortex = cast(
                                        _aerodynamics.RingVortex,
                                        next_wake_ring_vortices[
                                            chordwise_point_id, spanwise_point_id
                                        ],
                                    )

                                    next_wake_ring_vortex.update_position(
                                        Flrvp_GP1_CgP1=Flwrvp_GP1_CgP1,
                                        Frrvp_GP1_CgP1=Frwrvp_GP1_CgP1,
                                        Blrvp_GP1_CgP1=Blwrvp_GP1_CgP1,
                                        Brrvp_GP1_CgP1=Brwrvp_GP1_CgP1,
                                    )

                                    # Also, update the age of the wake RingVortex at
                                    # this position for the next time step.
                                    if self._current_step == 0:
                                        next_wake_ring_vortex.age = self.delta_time
                                    else:
                                        next_wake_ring_vortex.age += self.delta_time

                                if chordwise_point_id == 0:
                                    _panels = this_wing.panels
                                    assert _panels is not None

                                    # If this position corresponds to the front of
                                    # the wake, get the strength from the Panel's
                                    # bound RingVortex.
                                    this_panel: _panel.Panel = _panels[
                                        this_wing.num_chordwise_panels - 1,
                                        spanwise_point_id,
                                    ]

                                    this_ring_vortex = this_panel.ring_vortex
                                    assert this_ring_vortex is not None

                                    this_strength_copy = this_ring_vortex.strength

                                    # Then, for the next time step, make a new wake
                                    # RingVortex at this position in the wake,
                                    # with that bound RingVortex's strength, and add
                                    # it to the grid of the next time step's wake
                                    # RingVortices.
                                    next_wing.wake_ring_vortices[
                                        chordwise_point_id,
                                        spanwise_point_id,
                                    ] = _aerodynamics.RingVortex(
                                        Flrvp_GP1_CgP1=Flwrvp_GP1_CgP1,
                                        Frrvp_GP1_CgP1=Frwrvp_GP1_CgP1,
                                        Blrvp_GP1_CgP1=Blwrvp_GP1_CgP1,
                                        Brrvp_GP1_CgP1=Brwrvp_GP1_CgP1,
                                        strength=this_strength_copy,
                                    )

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
        observed from the Earth frame) at the center point of each bound RingVortex's
        right leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound RingVortex's right leg due to any motion defined in
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
        observed from the Earth frame) at the center point of each bound RingVortex's
        front leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound RingVortex's front leg due to any motion defined in
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
        observed from the Earth frame) at the center point of each bound RingVortex's
        left leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound RingVortex's left leg due to any motion defined in
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
        observed from the Earth frame) at the center point of each bound RingVortex's
        back leg due to any motion defined in Movement at the current time step.

        **Notes:**

        At each point, any apparent velocity due to Movement is opposite the motion due
        to Movement.

        :return: A (M, 3) ndarray of floats representing the apparent velocity (in the
            first Airplane's geometry axes, observed from the Earth frame) at the center
            point of each bound RingVortex's back leg due to any motion defined in
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
        this_movement: movements.movement.Movement = self.unsteady_problem.movement
        static = this_movement.static

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
            this_steady_problem: problems.SteadyProblem = self.steady_problems[step]
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
        # RMS loads and load coefficients.
        first_problem: problems.SteadyProblem = self.steady_problems[0]
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
                self.unsteady_problem.finalMeanForces_W.append(
                    np.mean(forces_W[airplane_id], axis=-1)
                )
                self.unsteady_problem.finalMeanForceCoefficients_W.append(
                    np.mean(force_coefficients_W[airplane_id], axis=-1)
                )
                self.unsteady_problem.finalMeanMoments_W_CgP1.append(
                    np.mean(moments_W_CgP1[airplane_id], axis=-1)
                )
                self.unsteady_problem.finalMeanMomentCoefficients_W_CgP1.append(
                    np.mean(moment_coefficients_W_CgP1[airplane_id], axis=-1)
                )

                self.unsteady_problem.finalRmsForces_W.append(
                    np.sqrt(
                        np.mean(
                            np.square(forces_W[airplane_id]),
                            axis=-1,
                        )
                    )
                )
                self.unsteady_problem.finalRmsForceCoefficients_W.append(
                    np.sqrt(
                        np.mean(
                            np.square(force_coefficients_W[airplane_id]),
                            axis=-1,
                        )
                    )
                )
                self.unsteady_problem.finalRmsMoments_W_CgP1.append(
                    np.sqrt(
                        np.mean(
                            np.square(moments_W_CgP1[airplane_id]),
                            axis=-1,
                        )
                    )
                )
                self.unsteady_problem.finalRmsMomentCoefficients_W_CgP1.append(
                    np.sqrt(
                        np.mean(
                            np.square(moment_coefficients_W_CgP1[airplane_id]),
                            axis=-1,
                        )
                    )
                )
