"""This module contains the class definition of this package's unsteady ring vortex
lattice solver.

This module contains the following classes:
    UnsteadyRingVortexLatticeMethodSolver: This is an aerodynamics solver that uses
    an unsteady ring vortex lattice method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import logging

import numpy as np
from tqdm import tqdm

from . import _aerodynamics
from . import _functions
from . import _parameter_validation
from . import _panel
from . import geometry
from . import problems


class UnsteadyRingVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses an unsteady ring vortex lattice method.

    This class contains the following public methods:

        run: This method runs the solver on the UnsteadyProblem.

        calculate_solution_velocity: This function takes in a group of points (in
        geometry axes, relative to the CG). At every point, it finds the fluid
        velocity (in geometry axes, observed from the Earth frame) at that point due
        to the freestream velocity and the induced velocity from every RingVortex.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, unsteady_problem):
        """This is the initialization method.

        :param unsteady_problem: UnsteadyProblem
            This is the UnsteadyProblem to be solved.
        :return: None
        """
        if not isinstance(unsteady_problem, problems.UnsteadyProblem):
            raise TypeError("unsteady_problem must be an UnsteadyProblem.")
        self.unsteady_problem = unsteady_problem

        self.num_steps = self.unsteady_problem.num_steps
        self.delta_time = self.unsteady_problem.delta_time
        self.first_results_step = self.unsteady_problem.first_results_step
        self._first_averaging_step = self.unsteady_problem.first_averaging_step
        self._current_step = None

        self.steady_problems = self.unsteady_problem.steady_problems

        self.current_airplanes = None
        self.current_operating_point = None
        self.num_airplanes = len(self.steady_problems[0].airplanes)
        num_panels = 0
        for airplane in self.steady_problems[0].airplanes:
            num_panels += airplane.num_panels
        self.num_panels = num_panels

        # Initialize attributes to hold aerodynamic data that pertain to the
        # simulation.
        self._currentVInf_G__E = None
        self._currentStackFreestreamWingInfluences__E = None
        self._currentGridWingWingInfluences__E = None
        self._currentStackWakeWingInfluences__E = None
        self._current_bound_vortex_strengths = None
        self._last_bound_vortex_strengths = None

        # Initialize attributes to hold geometric data that pertain to this
        # UnsteadyProblem.
        self.panels = None
        self.stackUnitNormals_G = None
        self.panel_areas = None

        # The current and last time step's collocation panel points (in geometry
        # axes, relative to the CG)
        self.stackCpp_G_Cg = None
        self._stackLastCpp_G_Cg = None

        # The current and last time step's back-right, front-right, front-left,
        # and back-left bound ring vortex points (in geometry axes, relative to
        # the CG).
        self.stackBrbrvp_G_Cg = None
        self.stackFrbrvp_G_Cg = None
        self.stackFlbrvp_G_Cg = None
        self.stackBlbrvp_G_Cg = None
        self._lastStackBrbrvp_G_Cg = None
        self._lastStackFrbrvp_G_Cg = None
        self._lastStackFlbrvp_G_Cg = None
        self._lastStackBlbrvp_G_Cg = None

        # The current and last time step's center bound line vortex points for
        # the right, front, left, and back legs (in geometry axes, relative to
        # the CG).
        self.stackCblvpr_G_Cg = None
        self.stackCblvpf_G_Cg = None
        self.stackCblvpl_G_Cg = None
        self.stackCblvpb_G_Cg = None
        self._lastStackCblvpr_G_Cg = None
        self._lastStackCblvpf_G_Cg = None
        self._lastStackCblvpl_G_Cg = None
        self._lastStackCblvpb_G_Cg = None

        # Right, front, left, and back bound ring vortex vectors (in geometry
        # axes).
        self.stackRbrv_G = None
        self.stackFbrv_G = None
        self.stackLbrv_G = None
        self.stackBbrv_G = None

        # Initialize variables to hold aerodynamic data that pertains details
        # about each Panel's location on its Wing.
        self.panel_is_trailing_edge = None
        self.panel_is_leading_edge = None
        self.panel_is_left_edge = None
        self.panel_is_right_edge = None

        # Initialize variables to hold aerodynamic data that pertains to the
        # wake at the current time step.
        self._current_wake_vortex_strengths = None
        self._current_wake_vortex_ages = None

        # The current time step's back-right, front-right, front-left, and
        # back-left wake ring vortex points (in geometry axes, relative to the
        # CG).
        self._currentStackBrwrvp_G_Cg = None
        self._currentStackFrwrvp_G_Cg = None
        self._currentStackFlwrvp_G_Cg = None
        self._currentStackBlwrvp_G_Cg = None

        # Initialize lists to store aerodynamic data about the wake at each time
        # step. These attributes are used by the output module to animate the
        # wake.
        self.list_num_wake_vortices = []
        # TODO: Determine if these private attributes are needed and if not
        #  delete them.
        self._list_wake_vortex_strengths = []
        self._list_wake_vortex_ages = []
        self.listStackBrwrvp_G_Cg = []
        self.listStackFrwrvp_G_Cg = []
        self.listStackFlwrvp_G_Cg = []
        self.listStackBlwrvp_G_Cg = []

        self.stackSeedPoints_G_Cg = None
        self.stackStreamlinePoints_G_Cg = None

    def run(
        self,
        logging_level="Warning",
        prescribed_wake=True,
        calculate_streamlines=True,
    ):
        """This method runs the solver on the UnsteadyProblem.

        :param logging_level: str, optional

            This parameter determines the detail of information that the solver's
            logger will output while running. The options are, in order of detail and
            severity, "Debug", "Info", "Warning", "Error", "Critical". The default
            value is "Warning".

        :param prescribed_wake: boolLike, optional

            This parameter determines if the solver uses a prescribed wake model. If
            False it will use a free-wake, which may be more accurate but will make
            the solver significantly slower. The default is True. It can be a boolean
            or a NumPy boolean and will be converted internally to a boolean.

        :param calculate_streamlines: boolLike, optional

            This parameter determines if the solver calculates streamlines emanating
            from the back of the wing after running the solver. It can be a boolean
            or a NumPy boolean and will be converted internally to a boolean.

        :return: None
        """
        logging_level = _parameter_validation.string_return_string(
            logging_level, "logging_level"
        )
        logging_level_value = _functions.convert_logging_level_name_to_value(
            logging_level
        )
        logging.basicConfig(level=logging_level_value)

        prescribed_wake = _parameter_validation.boolLike_return_bool(
            prescribed_wake, "prescribed_wake"
        )
        calculate_streamlines = _parameter_validation.boolLike_return_bool(
            calculate_streamlines, "calculate_streamlines"
        )

        # The following loop iterates through the time steps to populate
        # currently empty attributes with lists of pre-allocated arrays. During
        # the simulation, these arrays will be filled with data that describe
        # the wake. Using this method eliminates the need for computationally
        # expensive on-the-fly allocation and object copying.
        for step in range(self.num_steps):
            this_problem = self.steady_problems[step]
            these_airplanes = this_problem.airplanes

            # Loop through this time step's Airplanes to create a list of their Wings.
            these_wings = []
            for airplane in these_airplanes:
                these_wings.append(airplane.wings)

            # Iterate through the Wings to get the total number of spanwise Panels.
            this_num_spanwise_panels = 0
            for this_wing_set in these_wings:
                for this_wing in this_wing_set:
                    this_num_spanwise_panels += this_wing.num_spanwise_panels

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
            thisStackBrwrvp_G_Cg = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            thisStackFrwrvp_G_Cg = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            thisStackFlwrvp_G_Cg = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )
            thisStackBlwrvp_G_Cg = np.zeros(
                (this_num_wake_ring_vortices, 3), dtype=float
            )

            # Append this time step's ndarrays to the lists of ndarrays.
            self.list_num_wake_vortices.append(this_num_wake_ring_vortices)
            self._list_wake_vortex_strengths.append(this_wake_ring_vortex_strengths)
            self._list_wake_vortex_ages.append(this_wake_ring_vortex_ages)
            self.listStackBrwrvp_G_Cg.append(thisStackBrwrvp_G_Cg)
            self.listStackFrwrvp_G_Cg.append(thisStackFrwrvp_G_Cg)
            self.listStackFlwrvp_G_Cg.append(thisStackFlwrvp_G_Cg)
            self.listStackBlwrvp_G_Cg.append(thisStackBlwrvp_G_Cg)

        # The following loop attempts to predict how much time each time step will take,
        # relative to the other time steps. This data will be used to generate estimates
        # of how much longer a simulation will take, and create a smoothly advancing
        # progress bar.

        # Initialize list that will hold the approximate, relative times. This has
        # one more element than the number of time steps, because I will also use the
        # progress bar during the simulation initialization.
        approx_times = np.zeros(self.num_steps + 1, dtype=float)
        for step in range(1, self.num_steps):
            this_problem = self.steady_problems[step]
            these_airplanes = this_problem.airplanes

            # Iterate through this time step's Airplanes to get the total number of Wing
            # Panels.
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

        # Unless the logging level is at or above Warning, run the simulation with a
        # progress bar.
        with tqdm(
            total=approx_total_time,
            unit="",
            unit_scale=True,
            ncols=100,
            desc="Simulating",
            disable=logging_level_value != logging.WARNING,
            bar_format="{desc}:{percentage:3.0f}% |{bar}| Elapsed: {elapsed}, "
            "Remaining: {remaining}",
        ) as bar:
            # Initialize all the Airplanes' bound RingVortices.
            logging.info("Initializing all Airplanes' bound RingVortices.")
            self._initialize_panel_vortices()

            # Update the progress bar based on the initialization step's predicted
            # approximate, relative computing time.
            bar.update(n=float(approx_times[0]))

            # Iterate through the time steps.
            for step in range(self.num_steps):

                # Save attributes to hold the current step, Airplanes, and
                # OperatingPoint, and freestream velocity (in geometry axes, observed
                # from the Earth frame).
                self._current_step = step
                self.current_airplanes = self.steady_problems[
                    self._current_step
                ].airplanes
                self.current_operating_point = self.steady_problems[
                    self._current_step
                ].operating_point
                self._currentVInf_G__E = self.current_operating_point.vInf_G__E
                logging.info(
                    "Beginning time step "
                    + str(self._current_step)
                    + " out of "
                    + str(self.num_steps - 1)
                    + "."
                )

                # TODO: I think these steps are redundant, at least during the
                #  first time step. Consider dropping them.
                # Initialize attributes to hold aerodynamic data that pertain
                # to the simulation at this time step.
                self._currentVInf_G__E = self.current_operating_point.vInf_G__E
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
                self.stackUnitNormals_G = np.zeros((self.num_panels, 3), dtype=float)
                self.panel_areas = np.zeros(self.num_panels, dtype=float)

                self.stackCpp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._stackLastCpp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)

                self.stackBrbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFrbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFlbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self.stackBlbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackBrbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackFrbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackFlbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackBlbrvp_G_Cg = np.zeros((self.num_panels, 3), dtype=float)

                self.stackCblvpr_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpf_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpl_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self.stackCblvpb_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackCblvpr_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackCblvpf_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackCblvpl_G_Cg = np.zeros((self.num_panels, 3), dtype=float)
                self._lastStackCblvpb_G_Cg = np.zeros((self.num_panels, 3), dtype=float)

                self.stackRbrv_G = np.zeros((self.num_panels, 3), dtype=float)
                self.stackFbrv_G = np.zeros((self.num_panels, 3), dtype=float)
                self.stackLbrv_G = np.zeros((self.num_panels, 3), dtype=float)
                self.stackBbrv_G = np.zeros((self.num_panels, 3), dtype=float)

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
                self._currentStackBrwrvp_G_Cg = self.listStackBrwrvp_G_Cg[step]
                self._currentStackFrwrvp_G_Cg = self.listStackFrwrvp_G_Cg[step]
                self._currentStackFlwrvp_G_Cg = self.listStackFlwrvp_G_Cg[step]
                self._currentStackBlwrvp_G_Cg = self.listStackBlwrvp_G_Cg[step]

                self.stackSeedPoints_G_Cg = np.zeros((0, 3), dtype=float)

                # Collapse the geometry matrices into 1D ndarrays of attributes.
                logging.info("Collapsing the geometry.")
                self._collapse_geometry()

                # Find the matrix of Wing-Wing influence coefficients associated with
                # the Airplanes' geometries at this time step.
                logging.info("Calculating the Wing-Wing influences.")
                self._calculate_wing_wing_influences()

                # Find the normal velocity (in geometry axes, observed from the
                # Earth frame) at every collocation point due solely to the
                # freestream.
                logging.info("Calculating the freestream-Wing influences.")
                self._calculate_freestream_wing_influences()

                # Find the normal velocity (in geometry axes, observed from the
                # Earth frame) at every collocation point due solely to the
                # wake RingVortices.
                logging.info("Calculating the wake-Wing influences.")
                self._calculate_wake_wing_influences()

                # Solve for each bound RingVortex's strength.
                logging.info("Calculating bound RingVortex strengths.")
                self._calculate_vortex_strengths()

                # Solve for the forces (in geometry axes) and moments (in geometry axes,
                # relative to the CG) on each Panel.
                if self._current_step >= self.first_results_step:
                    logging.info("Calculating forces and moments.")
                    self._calculate_loads()

                # Shed RingVortices into the wake.
                logging.info("Shedding RingVortices into the wake.")
                self._populate_next_airplanes_wake(prescribed_wake=prescribed_wake)

                # Update the progress bar based on this time step's predicted
                # approximate, relative computing time.
                bar.update(n=float(approx_times[step + 1]))

            logging.info("Calculating averaged or final forces and moments.")
            self._finalize_loads()

        # Solve for the location of the streamlines coming off the Wings'
        # trailing edges, if requested.
        if calculate_streamlines:
            logging.info("Calculating streamlines.")
            _functions.calculate_streamlines(self)

    def _initialize_panel_vortices(self):
        """This method calculates the locations of the Airplanes' bound RingVortices'
        points, and then initializes the bound RingVortices.

        Every Panel has a RingVortex, which is a quadrangle whose front leg is a
        LineVortex at the Panel's quarter chord. The left and right legs are
        LineVortices running along the Panel's left and right legs. If the Panel is
        not along the trailing edge, they extend backwards and meet the back
        LineVortex, at the rear Panel's quarter chord. Otherwise, they extend
        backwards and meet the back LineVortex one quarter chord back from the
        Panel's back leg.

        :return: None
        """
        for steady_problem in self.steady_problems:
            # Find the freestream velocity (in geometry axes, observed from the Earth
            # frame) at this time step.
            vInf_G__E = steady_problem.operating_point.vInf_G__E

            # Iterate through this SteadyProblem's Airplanes' Wings.
            for airplane in steady_problem.airplanes:
                for wing in airplane.wings:

                    # Iterate through the Wing's chordwise and spanwise positions.
                    for chordwise_position in range(wing.num_chordwise_panels):
                        for spanwise_position in range(wing.num_spanwise_panels):
                            # Pull the panel object out of the Wing's 2D ndarray of Panels.
                            panel = wing.panels[chordwise_position, spanwise_position]

                            # Find the location of this Panel's front-left and
                            # front-right RingVortex points (in geometry axes,
                            # relative to the CG).
                            Flrvp_G_Cg = panel.Flbvp_G_Cg
                            Frrvp_G_Cg = panel.Frbvp_G_Cg

                            # Define the location of the back-left and back-right
                            # RingVortex points based on whether the Panel is along the
                            # trailing edge or not.
                            if not panel.is_trailing_edge:
                                next_chordwise_panel = wing.panels[
                                    chordwise_position + 1, spanwise_position
                                ]
                                Blrvp_G_Cg = next_chordwise_panel.Flbvp_G_Cg
                                Brrvp_G_Cg = next_chordwise_panel.Frbvp_G_Cg
                            else:
                                # As these vertices are directly behind the trailing
                                # edge, they are spaced back from their Panel's
                                # vertex by one quarter the distance traveled during
                                # a time step. This is to more accurately predict
                                # drag. More information can be found on pages 37-39
                                # of "Modeling of aerodynamic forces in flapping
                                # flight with the Unsteady Vortex Lattice Method" by
                                # Thomas Lambert.

                                # FIXME: I think this might be a bug. It looks like we
                                #  are already spacing the points back from the Panels'
                                #  rear points by a quarter-chord of each Panel. This
                                #  is spacing even further beyond that. Double check if
                                #  that is actually correct.

                                Blrvp_G_Cg = (
                                    Flrvp_G_Cg
                                    + (panel.Blpp_G_Cg - panel.Flpp_G_Cg)
                                    + vInf_G__E * self.delta_time * 0.25
                                )
                                Brrvp_G_Cg = (
                                    Frrvp_G_Cg
                                    + (panel.Brpp_G_Cg - panel.Frpp_G_Cg)
                                    + vInf_G__E * self.delta_time * 0.25
                                )

                            # Initialize the Panel's RingVortex.
                            panel.ring_vortex = _aerodynamics.RingVortex(
                                Flrvp_G_Cg=Flrvp_G_Cg,
                                Frrvp_G_Cg=Frrvp_G_Cg,
                                Blrvp_G_Cg=Blrvp_G_Cg,
                                Brrvp_G_Cg=Brrvp_G_Cg,
                                strength=None,
                            )

    def _collapse_geometry(self):
        """This method converts attributes of the UnsteadyProblem's geometry into 1D
        ndarrays. This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize variables to hold the global position of the Panel and the
        # wake RingVortex as we iterate through them.
        global_panel_position = 0
        global_wake_ring_vortex_position = 0

        # Iterate through the current time step's Airplanes' Wings.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:

                # Convert this Wing's 2D ndarray of Panels and wake RingVortices into
                # 1D ndarrays.
                panels = np.ravel(wing.panels)
                wake_ring_vortices = np.ravel(wing.wake_ring_vortices)

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
                    # Update the solver's list of attributes with this wake RingVortex's
                    # attributes.
                    self._current_wake_vortex_strengths[
                        global_wake_ring_vortex_position
                    ] = wake_ring_vortex.strength
                    self._current_wake_vortex_ages[global_wake_ring_vortex_position] = (
                        wake_ring_vortex.age
                    )
                    self._currentStackFrwrvp_G_Cg[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Frrvp_G_Cg
                    self._currentStackFlwrvp_G_Cg[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Flrvp_G_Cg
                    self._currentStackBlwrvp_G_Cg[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Blrvp_G_Cg
                    self._currentStackBrwrvp_G_Cg[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Brrvp_G_Cg

                    # Increment the global wake RingVortex position variable.
                    global_wake_ring_vortex_position += 1

        if self._current_step > 0:

            # Reset the global Panel position variable.
            global_panel_position = 0

            last_airplanes = self.steady_problems[self._current_step - 1].airplanes

            # Iterate through the last time step's Airplanes' Wings.
            for last_airplane in last_airplanes:
                for wing in last_airplane.wings:

                    # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                    panels = np.ravel(wing.panels)

                    # Iterate through the 1D ndarray of this Wing's Panels.
                    panel: _panel.Panel
                    for panel in panels:
                        # Update the solver's list of attributes with this Panel's
                        # attributes.
                        self._stackLastCpp_G_Cg[global_panel_position, :] = (
                            panel.Cpp_G_Cg
                        )
                        self._last_bound_vortex_strengths[global_panel_position] = (
                            panel.ring_vortex.strength
                        )
                        # TODO: Test if we can replace the calls to LineVortex
                        #  attributes with calls to RingVortex attributes.
                        self._lastStackBrbrvp_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.right_leg.Slvp_G_Cg
                        )
                        self._lastStackFrbrvp_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.right_leg.Elvp_G_Cg
                        )
                        self._lastStackFlbrvp_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.left_leg.Slvp_G_Cg
                        )
                        self._lastStackBlbrvp_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.left_leg.Elvp_G_Cg
                        )
                        self._lastStackCblvpr_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.right_leg.Clvp_G_Cg
                        )
                        self._lastStackCblvpf_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.front_leg.Clvp_G_Cg
                        )
                        self._lastStackCblvpl_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.left_leg.Clvp_G_Cg
                        )
                        self._lastStackCblvpb_G_Cg[global_panel_position, :] = (
                            panel.ring_vortex.back_leg.Clvp_G_Cg
                        )

                        # Increment the global Panel position variable.
                        global_panel_position += 1

    def _calculate_wing_wing_influences(self):
        """This method finds the 2d ndarray of Wing-Wing influence coefficients (
        observed from the Earth frame).

        :return: None
        """
        # Find the 2D ndarray of normalized velocities (in geometry axes, observed
        # from the Earth frame) induced at each Panel's collocation point by each bound
        # RingVortex. The answer is normalized because the solver's list of bound
        # RingVortex strengths was initialized to all be 1.0. This will be updated once
        # the correct strengths are calculated.
        gridNormVIndCpp_G__E = _aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.stackCpp_G_Cg,
            stackBrrvp_G_Cg=self.stackBrbrvp_G_Cg,
            stackFrrvp_G_Cg=self.stackFrbrvp_G_Cg,
            stackFlrvp_G_Cg=self.stackFlbrvp_G_Cg,
            stackBlrvp_G_Cg=self.stackBlbrvp_G_Cg,
            strengths=self._current_bound_vortex_strengths,
            ages=None,
            nu=self.current_operating_point.nu,
        )

        # Take the batch dot product of the normalized induced velocities (in
        # geometry axes, observed from the Earth frame) with each Panel's unit
        # normal direction (in geometry axes). This is now the 2D ndarray of
        # Wing-Wing influence coefficients (observed from the Earth frame).
        self._currentGridWingWingInfluences__E = np.einsum(
            "...k,...k->...",
            gridNormVIndCpp_G__E,
            np.expand_dims(self.stackUnitNormals_G, axis=1),
        )

    def _calculate_freestream_wing_influences(self):
        """This method finds the 1D ndarray of freestream-Wing influence coefficients
        (observed from the Earth frame).

        Note: This method also includes the influence coefficients due to motion
        defined in Movement (observed from the Earth frame) at every collocation point.

        :return: None
        """
        # Find the normal components of the freestream-only-Wing influence
        # coefficients (observed from the Earth frame) at each Panel's
        # collocation point by taking a batch dot product.
        currentStackFreestreamOnlyWingInfluences__E = np.einsum(
            "ij,j->i",
            self.stackUnitNormals_G,
            self._currentVInf_G__E,
        )

        # Get the current apparent velocities at each Panel's collocation point
        # due to any motion defined in Movement (in geometry axes, observed from
        # the Earth frame).
        currentStackMovementV_G__E = (
            self._calculate_current_movement_velocities_at_collocation_points()
        )

        # Get the current motion influence coefficients at each Panel's
        # collocation point (observed from the Earth frame) by taking a batch dot
        # product.
        currentStackMovementInfluences__E = np.einsum(
            "ij,ij->i",
            self.stackUnitNormals_G,
            currentStackMovementV_G__E,
        )

        # Calculate the total current freestream-Wing influence coefficients by
        # summing the freestream-only influence coefficients and the motion
        # influence coefficients (all observed from the Earth frame).
        self._currentStackFreestreamWingInfluences__E = (
            currentStackFreestreamOnlyWingInfluences__E
            + currentStackMovementInfluences__E
        )

    def _calculate_wake_wing_influences(self):
        """This method finds the 1D ndarray of the wake-Wing influence coefficients (
        observed from the Earth frame) associated with the UnsteadyProblem at the
        current time step.

        Note: If the current time step is the first time step, no wake has been shed,
        so this method will return zero for all the wake-Wing influence coefficients
        (observed from the Earth frame).

        :return: None
        """
        if self._current_step > 0:
            # Get the velocities (in geometry axes, observed from the Earth frame)
            # induced by the wake RingVortices at each Panel's collocation point.
            currentStackWakeV_G__E = (
                _aerodynamics.collapsed_velocities_from_ring_vortices(
                    stackP_G_Cg=self.stackCpp_G_Cg,
                    stackBrrvp_G_Cg=self._currentStackBrwrvp_G_Cg,
                    stackFrrvp_G_Cg=self._currentStackFrwrvp_G_Cg,
                    stackFlrvp_G_Cg=self._currentStackFlwrvp_G_Cg,
                    stackBlrvp_G_Cg=self._currentStackBlwrvp_G_Cg,
                    strengths=self._current_wake_vortex_strengths,
                    ages=self._current_wake_vortex_ages,
                    nu=self.current_operating_point.nu,
                )
            )

            # Get the current wake-Wing influence coefficients (observed from
            # the Earth frame) by taking a batch dot product with each Panel's
            # normal vector (in geometry axes).
            self._currentStackWakeWingInfluences__E = np.einsum(
                "ij,ij->i", currentStackWakeV_G__E, self.stackUnitNormals_G
            )

        else:
            # If this is the first time step, set all the current Wake-wing
            # influence coefficients to 0.0 (observed from the Earth frame)
            # because no wake RingVortices have been shed.
            self._currentStackWakeWingInfluences__E = np.zeros(
                self.num_panels, dtype=float
            )

    def _calculate_vortex_strengths(self):
        """Solve for the strength of each Panel's bound RingVortex.

        :return: None
        """
        self._current_bound_vortex_strengths = np.linalg.solve(
            self._currentGridWingWingInfluences__E,
            -self._currentStackWakeWingInfluences__E
            - self._currentStackFreestreamWingInfluences__E,
        )

        # Update the bound RingVortices' strengths.
        for panel_num in range(self.panels.size):
            panel = self.panels[panel_num]

            panel.ring_vortex.update_strength(
                self._current_bound_vortex_strengths[panel_num]
            )

    def calculate_solution_velocity(self, stackP_G_Cg):
        """This function takes in a group of points (in geometry axes, relative to
        the CG). At every point, it finds the fluid velocity (in geometry axes,
        observed from the Earth frame) at that point due to the freestream velocity
        and the induced velocity from every RingVortex.

        Note: This method assumes that the correct strengths for the RingVortices and
        HorseshoeVortices have already been calculated and set. This method also does
        not include the velocity due to the Movement's motion at any of the points
        provided, as it has no way of knowing if any of the points lie on panels.

        :param stackP_G_Cg: (N,3) array-like of numbers

            Positions of the evaluation points (in geometry axes, relative to the
            CG). Can be any array-like object (tuple, list, or ndarray) with size (N,
            3) that has numeric elements (int or float). Values are converted to
            floats internally. The units are in meters.

        :return: (N,3) ndarray of floats

            The velocity (in geometry axes, observed from the Earth frame) at every
            evaluation point due to the summed effects of the freestream velocity and
            the induced velocity from every RingVortex. The units are in meters per
            second.
        """
        stackP_G_Cg = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                stackP_G_Cg, "stackP_G_Cg"
            )
        )

        stackBoundRingVInd_G__E = _aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackBrrvp_G_Cg=self.stackBrbrvp_G_Cg,
            stackFrrvp_G_Cg=self.stackFrbrvp_G_Cg,
            stackFlrvp_G_Cg=self.stackFlbrvp_G_Cg,
            stackBlrvp_G_Cg=self.stackBlbrvp_G_Cg,
            strengths=self._current_bound_vortex_strengths,
            ages=None,
            nu=self.current_operating_point.nu,
        )
        stackWakeRingVInd_G__E = _aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=stackP_G_Cg,
            stackBrrvp_G_Cg=self._currentStackBrwrvp_G_Cg,
            stackFrrvp_G_Cg=self._currentStackFrwrvp_G_Cg,
            stackFlrvp_G_Cg=self._currentStackFlwrvp_G_Cg,
            stackBlrvp_G_Cg=self._currentStackBlwrvp_G_Cg,
            strengths=self._current_wake_vortex_strengths,
            ages=self._current_wake_vortex_ages,
            nu=self.current_operating_point.nu,
        )

        return stackBoundRingVInd_G__E + stackWakeRingVInd_G__E + self._currentVInf_G__E

    def _calculate_loads(self):
        """Calculate the forces (in geometry axes) and moments (in geometry axes,
        relative to the CG) on every Panel.

        Citation: This method uses logic described on pages 9-11 of "Modeling of
        aerodynamic forces in flapping flight with the Unsteady Vortex Lattice
        Method" by Thomas Lambert.

        Note: This method assumes that the correct strengths for the RingVortices and
        HorseshoeVortices have already been calculated and set.

        :return: None
        """
        # Initialize a variable to hold the global Panel position as we iterate
        # through them.
        global_panel_position = 0

        # Initialize three 1D ndarrays to hold the effective strength of the Panels'
        # RingVortices' LineVortices.
        effective_right_vortex_line_strengths = np.zeros(self.num_panels, dtype=float)
        effective_front_vortex_line_strengths = np.zeros(self.num_panels, dtype=float)
        effective_left_vortex_line_strengths = np.zeros(self.num_panels, dtype=float)

        # Iterate through the Airplanes' Wings.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:
                # Convert this Wing's 2D ndarray of Panels into a 1D ndarray.
                panels = np.ravel(wing.panels)

                # Iterate through this Wing's 1D ndarray of Panels.
                for panel in panels:

                    # FIXME: After rereading pages 9-10 of "Modeling of
                    #  aerodynamic forces in flapping flight with the Unsteady
                    #  Vortex Lattice Method" by Thomas Lambert, I think our
                    #  implementation here is critically wrong. Consider we have
                    #  a wing with a (1,2) ndarray of Panels. Let's call them
                    #  Panel A and Panel B. With our current method, we calculate
                    #  the force on Panel A's right LineVortex as though it had a
                    #  strength of Gamma_a - Gamma_b, and the force on Panel B's
                    #  left LineVortex as Gamma_b - Gamma_a. I think these forces
                    #  will precisely cancel-out!

                    if panel.is_right_edge:
                        # Set the effective right LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_right_vortex_line_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_right = wing.panels[
                            panel.local_chordwise_position,
                            panel.local_spanwise_position + 1,
                        ]

                        # Set the effective right LineVortex strength to the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the right.
                        effective_right_vortex_line_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - panel_to_right.ring_vortex.strength
                        )

                    if panel.is_leading_edge:
                        # Set the effective front LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_front_vortex_line_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_front = wing.panels[
                            panel.local_chordwise_position - 1,
                            panel.local_spanwise_position,
                        ]

                        # Set the effective front LineVortex strength to the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel in front of it.
                        effective_front_vortex_line_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - panel_to_front.ring_vortex.strength
                        )

                    if panel.is_left_edge:
                        # Set the effective left LineVortex strength to this Panel's
                        # RingVortex's strength.
                        effective_left_vortex_line_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                        )
                    else:
                        panel_to_left = wing.panels[
                            panel.local_chordwise_position,
                            panel.local_spanwise_position - 1,
                        ]

                        # Set the effective left LineVortex strength to the
                        # difference between this Panel's RingVortex's strength,
                        # and the RingVortex's strength of the Panel to the left.
                        effective_left_vortex_line_strengths[global_panel_position] = (
                            self._current_bound_vortex_strengths[global_panel_position]
                            - panel_to_left.ring_vortex.strength
                        )

                    # Increment the global Panel position variable.
                    global_panel_position += 1

        # Calculate the velocity (in geometry axes, observed from the Earth frame) at
        # the center of every Panels' RingVortex's right LineVortex,
        # front LineVortex, and left LineVortex.
        stackVelocityRightLineVortexCenters_G__E = (
            self.calculate_solution_velocity(stackP_G_Cg=self.stackCblvpr_G_Cg)
            + self._calculate_current_movement_velocities_at_right_leg_centers()
        )
        stackVelocityFrontLineVortexCenters_G__E = (
            self.calculate_solution_velocity(stackP_G_Cg=self.stackCblvpf_G_Cg)
            + self._calculate_current_movement_velocities_at_front_leg_centers()
        )
        stackVelocityLeftLineVortexCenters_G__E = (
            self.calculate_solution_velocity(stackP_G_Cg=self.stackCblvpl_G_Cg)
            + self._calculate_current_movement_velocities_at_left_leg_centers()
        )

        # Using the effective LineVortex strengths and the Kutta-Joukowski theorem,
        # find the forces (in geometry axes) on the Panels' RingVortex's right
        # LineVortex, front LineVortex, and left LineVortex using the effective
        # vortex strengths. Also calculate the unsteady component of the force on each
        # Panel, which is derived from the unsteady Bernoulli equation.
        rightLegForces_G = (
            self.current_operating_point.rho
            * np.expand_dims(effective_right_vortex_line_strengths, axis=1)
            * _functions.numba_1d_explicit_cross(
                stackVelocityRightLineVortexCenters_G__E,
                self.stackRbrv_G,
            )
        )
        frontLegForces_G = (
            self.current_operating_point.rho
            * np.expand_dims(effective_front_vortex_line_strengths, axis=1)
            * _functions.numba_1d_explicit_cross(
                stackVelocityFrontLineVortexCenters_G__E,
                self.stackFbrv_G,
            )
        )
        leftLegForces_G = (
            self.current_operating_point.rho
            * np.expand_dims(effective_left_vortex_line_strengths, axis=1)
            * _functions.numba_1d_explicit_cross(
                stackVelocityLeftLineVortexCenters_G__E,
                self.stackLbrv_G,
            )
        )

        # FIXME: The unsteady force calculation may be wrong. In short, it looks
        #  like we are getting a result that's -1 times what we expect. However,
        #  all the reference equations support our implementation. More testing
        #  and investigation is needed. See
        #  https://github.com/camUrban/PteraSoftware/issues/27 for more details
        #  and updates.

        # Calculate the unsteady component of the force on each Panel (in
        # geometry axes), which is derived from the unsteady Bernoulli equation.
        unsteady_forces_G = (
            self.current_operating_point.rho
            * np.expand_dims(
                (
                    self._current_bound_vortex_strengths
                    - self._last_bound_vortex_strengths
                ),
                axis=1,
            )
            * np.expand_dims(self.panel_areas, axis=1)
            * self.stackUnitNormals_G
            / self.delta_time
        )

        forces_G = (
            rightLegForces_G + frontLegForces_G + leftLegForces_G + unsteady_forces_G
        )

        # Find the moments (in geometry axes, relative to the CG) on the Panels'
        # RingVortex's right LineVortex, front LineVortex, and left LineVortex.
        rightLegMoments_G_Cg = _functions.numba_1d_explicit_cross(
            self.stackCblvpr_G_Cg,
            rightLegForces_G,
        )
        frontLegMoments_G_Cg = _functions.numba_1d_explicit_cross(
            self.stackCblvpf_G_Cg,
            frontLegForces_G,
        )
        leftLegMoments_G_Cg = _functions.numba_1d_explicit_cross(
            self.stackCblvpl_G_Cg,
            leftLegForces_G,
        )

        # FIXME: The unsteady moment calculation may also be wrong. Right now,
        #  we are using the position of the collocation points (in geometry axes,
        #  relative to the CG) to calculate the moment, but shouldn't we instead
        #  use the  to the locations of the Panels' centroids?

        # Find the moments (in geometry axes, relative to the CG) due to the
        # unsteady component of the force on each Panel.
        unsteady_moments_G_Cg = _functions.numba_1d_explicit_cross(
            self.stackCpp_G_Cg,
            unsteady_forces_G,
        )

        moments_G_Cg = (
            rightLegMoments_G_Cg
            + frontLegMoments_G_Cg
            + leftLegMoments_G_Cg
            + unsteady_moments_G_Cg
        )

        _functions.process_solver_loads(self, forces_G, moments_G_Cg)

    def _populate_next_airplanes_wake(self, prescribed_wake=True):
        """This method updates the next time step's Airplanes' wakes.

        :param prescribed_wake: Bool, optional
            This parameter determines if the solver uses a prescribed wake model. If
            false it will use a free-wake, which may be more accurate but will make
            the solver significantly slower. The default is True.
        :return: None
        """
        # Populate the locations of the next time step's Airplanes' wake RingVortex
        # points.
        self._populate_next_airplanes_wake_vortex_points(
            prescribed_wake=prescribed_wake
        )

        # Populate the locations of the next time step's Airplanes' wake RingVortices.
        self._populate_next_airplanes_wake_vortices()

    def _populate_next_airplanes_wake_vortex_points(self, prescribed_wake=True):
        """This method populates the locations of the next time step's Airplanes'
        wake RingVortex points.

        This method is not vectorized but its loops only consume 1.1% of the runtime,
        so I have kept it as is for increased readability.

        :param prescribed_wake: Bool, optional
            This parameter determines if the solver uses a prescribed wake model. If
            false it will use a free-wake, which may be more accurate but will make
            the solver significantly slower. The default is True.
        :return: None
        """
        # Check that this isn't the last time step.
        if self._current_step < self.num_steps - 1:

            # Get the next time step's Airplanes.
            next_airplanes = self.steady_problems[self._current_step + 1].airplanes

            # Get the current Airplanes' combined number of Wings.
            num_wings = 0
            airplane: geometry.airplane.Airplane
            for airplane in self.current_airplanes:
                num_wings += len(airplane.wings)

            # Iterate through this time step's Airplanes' successor objects.
            next_airplane: geometry.airplane.Airplane
            for airplane_id, next_airplane in enumerate(next_airplanes):

                # Iterate through the next Airplane's Wings.
                next_wing: geometry.wing.Wing
                for wing_id, next_wing in enumerate(next_airplane.wings):

                    # Get the Wings at this position from the current Airplane.
                    this_wing = self.current_airplanes[airplane_id].wings[wing_id]

                    # Check if this is the first time step.
                    if self._current_step == 0:

                        # Get the current Wing's number of chordwise and spanwise
                        # panels.
                        num_spanwise_panels = this_wing.num_spanwise_panels
                        num_chordwise_panels = this_wing.num_chordwise_panels

                        # Set the chordwise position to be at the trailing edge.
                        chordwise_panel_id = num_chordwise_panels - 1

                        # Initialize a ndarray to hold the points of the new row of
                        # wake RingVortices (in geometry axes, relative to the CG).
                        newRowWrvp_G_Cg = np.zeros(
                            (1, num_spanwise_panels + 1, 3), dtype=float
                        )

                        # Iterate through the spanwise Panel positions.
                        for spanwise_panel_id in range(num_spanwise_panels):

                            # Get the next time step's Wing's Panel at this location.
                            next_panel = next_wing.panels[
                                chordwise_panel_id, spanwise_panel_id
                            ]

                            # The position of the new front left wake RingVortex's
                            # point is the next time step's Panel's bound
                            # RingVortex's back left point.
                            newFlwrvp_G_Cg = next_panel.ring_vortex.Blrvp_G_Cg

                            # Add this to the row of new wake RingVortex points.
                            newRowWrvp_G_Cg[0, spanwise_panel_id] = newFlwrvp_G_Cg

                            # If the Panel is at the right edge of the Wing, add its
                            # back right bound RingVortex point to the row of new
                            # wake RingVortex points.
                            if spanwise_panel_id == (num_spanwise_panels - 1):
                                newRowWrvp_G_Cg[0, spanwise_panel_id + 1] = (
                                    next_panel.ring_vortex.Brrvp_G_Cg
                                )

                        # Set the next time step's Wing's grid of wake RingVortex
                        # points to a copy of the row of new wake RingVortex points.
                        # This is correct because it is currently the first time step.
                        next_wing.gridWrvp_G_Cg = np.copy(newRowWrvp_G_Cg)

                        # Initialize variables to hold the number of spanwise wake
                        # RingVortex points.
                        num_spanwise_points = num_spanwise_panels + 1

                        # Initialize a new ndarray to hold the second new row of wake
                        # RingVortex points (in geometry axes, relative to the CG).
                        secondNewRowWrvp_G_Cg = np.zeros(
                            (1, num_spanwise_panels + 1, 3), dtype=float
                        )

                        # Iterate through the spanwise points.
                        for spanwise_point_id in range(num_spanwise_points):
                            # Get the corresponding point from the first row.
                            Wrvp_G_Cg = next_wing.gridWrvp_G_Cg[0, spanwise_point_id]

                            # If the wake is prescribed, set the velocity at this
                            # point to the freestream velocity (in geometry axes,
                            # observed from the Earth frame). Otherwise, set the
                            # velocity to the solution velocity at this point (in
                            # geometry axes, observed from the Earth frame).
                            if prescribed_wake:
                                vWrvp_G__E = self._currentVInf_G__E
                            else:
                                vWrvp_G__E = self.calculate_solution_velocity(
                                    np.expand_dims(Wrvp_G_Cg, axis=0)
                                )

                            # Update the second new row with the interpolated
                            # position of the first point.
                            secondNewRowWrvp_G_Cg[0, spanwise_point_id] = (
                                Wrvp_G_Cg + vWrvp_G__E * self.delta_time
                            )

                        # Update the next time step's Wing's grid of wake RingVortex
                        # points by vertically stacking the new second row below it.
                        next_wing.gridWrvp_G_Cg = np.vstack(
                            (
                                next_wing.gridWrvp_G_Cg,
                                secondNewRowWrvp_G_Cg,
                            )
                        )

                    # If this isn't the first time step, then do this.
                    else:
                        # Set the next time step's Wing's grid of wake RingVortex
                        # points to a copy of this time step's Wing's grid of wake
                        # RingVortex points.
                        next_wing.gridWrvp_G_Cg = np.copy(this_wing.gridWrvp_G_Cg)

                        # Get the number of chordwise and spanwise points.
                        num_chordwise_points = next_wing.gridWrvp_G_Cg.shape[0]
                        num_spanwise_points = next_wing.gridWrvp_G_Cg.shape[1]

                        # Iterate through the chordwise and spanwise point positions.
                        for chordwise_point_id in range(num_chordwise_points):
                            for spanwise_point_id in range(num_spanwise_points):
                                # Get the wake RingVortex point at this position.
                                Wrvp_G_Cg = next_wing.gridWrvp_G_Cg[
                                    chordwise_point_id,
                                    spanwise_point_id,
                                ]

                                # If the wake is prescribed, set the velocity at this
                                # point to the freestream velocity (in geometry axes,
                                # observed from the Earth frame). Otherwise, set the
                                # velocity to the solution velocity at this point (in
                                # geometry axes, observed from the Earth frame).
                                if prescribed_wake:
                                    vWrvp_G__E = self._currentVInf_G__E
                                else:
                                    vWrvp_G__E = np.squeeze(
                                        self.calculate_solution_velocity(
                                            np.expand_dims(Wrvp_G_Cg, axis=0)
                                        )
                                    )

                                # Update this point with its interpolated position.
                                next_wing.gridWrvp_G_Cg[
                                    chordwise_point_id, spanwise_point_id
                                ] += (vWrvp_G__E * self.delta_time)

                        # Find the chordwise position of the Wing's trailing edge.
                        chordwise_panel_id = this_wing.num_chordwise_panels - 1

                        # Initialize a new ndarray to hold the new row of wake
                        # RingVortex vertices.
                        newRowWrvp_G_Cg = np.zeros(
                            (1, this_wing.num_spanwise_panels + 1, 3), dtype=float
                        )

                        # Iterate spanwise through the trailing edge Panels.
                        for spanwise_panel_id in range(this_wing.num_spanwise_panels):
                            # Get the Panel at this location on the next time step's
                            # Airplane's Wing.
                            next_panel = next_wing.panels[
                                chordwise_panel_id, spanwise_panel_id
                            ]

                            # Add the Panel's back left bound RingVortex point to the
                            # grid of new wake RingVortex points.
                            newRowWrvp_G_Cg[0, spanwise_panel_id] = (
                                next_panel.ring_vortex.Blrvp_G_Cg
                            )

                            # If the Panel is at the right edge of the Wing, add its
                            # back right bound RingVortex point to the grid of new
                            # wake RingVortex vertices.
                            if spanwise_panel_id == (this_wing.num_spanwise_panels - 1):
                                newRowWrvp_G_Cg[0, spanwise_panel_id + 1] = (
                                    next_panel.ring_vortex.Brrvp_G_Cg
                                )

                        # Stack the new row of wake RingVortex points above the
                        # Wing's grid of wake RingVortex points.
                        next_wing.gridWrvp_G_Cg = np.vstack(
                            (
                                newRowWrvp_G_Cg,
                                next_wing.gridWrvp_G_Cg,
                            )
                        )

    def _populate_next_airplanes_wake_vortices(self):
        """This method populates the locations and strengths of the next time step's
        wake RingVortices.

        This method is not vectorized but its loops only consume 0.4% of the runtime,
        so I have kept it as is for increased readability.

        :return: None
        """
        # Check if the current time step is not the last step.
        if self._current_step < self.num_steps - 1:

            # Get the next time step's Airplanes.
            next_airplanes = self.steady_problems[self._current_step + 1].airplanes

            # Iterate through the next time step's Airplanes.
            next_airplane: geometry.airplane.Airplane
            for airplane_id, next_airplane in enumerate(next_airplanes):

                # For a given Airplane in the next time step, iterate through its
                # predecessor's Wings.
                this_wing: geometry.wing.Wing
                for wing_id, this_wing in enumerate(
                    self.current_airplanes[airplane_id].wings
                ):
                    next_wing = next_airplane.wings[wing_id]

                    # Get the next time step's Wing's grid of wake RingVortex points.
                    nextGridWrvp_G_Cg = next_wing.gridWrvp_G_Cg

                    # Find the number of chordwise and spanwise points in the next
                    # Wing's grid of wake RingVortex points.
                    num_chordwise_points = nextGridWrvp_G_Cg.shape[0]
                    num_spanwise_points = nextGridWrvp_G_Cg.shape[1]

                    this_wing_wake_ring_vortices = (
                        self.current_airplanes[airplane_id]
                        .wings[wing_id]
                        .wake_ring_vortices
                    )

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
                            # Set booleans to determine if this point is on the
                            # right and/or trailing edge of the wake.
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
                                # this position (in geometry axes, relative to the
                                # CG), for the next time step.
                                Flwrvp_G_Cg = nextGridWrvp_G_Cg[
                                    chordwise_point_id, spanwise_point_id
                                ]
                                Frwrvp_G_Cg = nextGridWrvp_G_Cg[
                                    chordwise_point_id,
                                    spanwise_point_id + 1,
                                ]
                                Blwrvp_G_Cg = nextGridWrvp_G_Cg[
                                    chordwise_point_id + 1,
                                    spanwise_point_id,
                                ]
                                Brwrvp_G_Cg = nextGridWrvp_G_Cg[
                                    chordwise_point_id + 1,
                                    spanwise_point_id + 1,
                                ]

                                if chordwise_point_id > 0:
                                    # If this isn't the front of the wake, update the
                                    # position of the wake RingVortex at this
                                    # location for the next time step.
                                    next_wake_ring_vortex = (
                                        next_wing.wake_ring_vortices[
                                            chordwise_point_id,
                                            spanwise_point_id,
                                        ]
                                    )
                                    assert isinstance(
                                        next_wake_ring_vortex, _aerodynamics.RingVortex
                                    )
                                    next_wake_ring_vortex.update_position(
                                        Flrvp_G_Cg=Flwrvp_G_Cg,
                                        Frrvp_G_Cg=Frwrvp_G_Cg,
                                        Blrvp_G_Cg=Blwrvp_G_Cg,
                                        Brrvp_G_Cg=Brwrvp_G_Cg,
                                    )

                                    # Also, update the age of the wake RingVortex at
                                    # this position for the next time step.
                                    if self._current_step == 0:
                                        next_wake_ring_vortex.age = self.delta_time
                                    else:
                                        next_wake_ring_vortex.age += self.delta_time

                                if chordwise_point_id == 0:
                                    # If this position corresponds to the front of
                                    # the wake, get the strength from the Panel's
                                    # bound RingVortex.
                                    this_strength_copy = this_wing.panels[
                                        this_wing.num_chordwise_panels - 1,
                                        spanwise_point_id,
                                    ].ring_vortex.strength

                                    # Then, for the next time step, make a new wake
                                    # RingVortex at this position in the wake,
                                    # with that bound RingVortex's strength, and add
                                    # it to the grid of the next time step's wake
                                    # RingVortices.
                                    next_wing.wake_ring_vortices[
                                        chordwise_point_id,
                                        spanwise_point_id,
                                    ] = _aerodynamics.RingVortex(
                                        Flrvp_G_Cg=Flwrvp_G_Cg,
                                        Frrvp_G_Cg=Frwrvp_G_Cg,
                                        Blrvp_G_Cg=Blwrvp_G_Cg,
                                        Brrvp_G_Cg=Brwrvp_G_Cg,
                                        strength=this_strength_copy,
                                    )

    def _calculate_current_movement_velocities_at_collocation_points(self):
        """Get the current apparent velocities (in geometry axes, observed from the
        Earth frame) at each Panel's collocation point due to any motion defined in
        Movement.

        Note: At each point, any apparent velocity due to Movement is opposite the
        motion due to Movement.

        :return: (M, 3) ndarray of floats

            This is a ndarray containing the apparent velocity (in geometry axes,
            observed from the Earth frame) at each Panel's collocation point due to
            any motion defined in Movement. Its units are in meters per second. If
            the current time step is the first time step, these velocities will all
            be all zeros.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return -(self.stackCpp_G_Cg - self._stackLastCpp_G_Cg) / self.delta_time

    def _calculate_current_movement_velocities_at_right_leg_centers(self):
        """Get the current apparent velocities (in geometry axes, observed from the
        Earth frame) at the center point of each bound RingVortex's right leg due to
        any motion defined in Movement.

        Note: At each point, any apparent velocity due to Movement is opposite the
        motion due to Movement.

        :return: (M, 3) ndarray of floats

            This is a ndarray containing the apparent velocity (in geometry axes,
            observed from the Earth frame) at the center point of each bound
            RingVortex's right leg due to any motion defined in Movement. Its units
            are in meters per second. If the current time step is the first time
            step, these velocities will all be all zeros.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return -(self.stackCblvpr_G_Cg - self._lastStackCblvpr_G_Cg) / self.delta_time

    def _calculate_current_movement_velocities_at_front_leg_centers(self):
        """Get the current apparent velocities (in geometry axes, observed from the
        Earth frame) at the center point of each bound RingVortex's front leg due to
        any motion defined in Movement.

        Note: At each point, any apparent velocity due to Movement is opposite the
        motion due to Movement.

        :return: (M, 3) ndarray of floats

            This is a ndarray containing the apparent velocity (in geometry axes,
            observed from the Earth frame) at the center point of each bound
            RingVortex's front leg due to any motion defined in Movement. Its units
            are in meters per second. If the current time step is the first time
            step, these velocities will all be all zeros.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return -(self.stackCblvpf_G_Cg - self._lastStackCblvpf_G_Cg) / self.delta_time

    def _calculate_current_movement_velocities_at_left_leg_centers(self):
        """Get the current apparent velocities (in geometry axes, observed from the
        Earth frame) at the center point of each bound RingVortex's left leg due to
        any motion defined in Movement.

        Note: At each point, any apparent velocity due to Movement is opposite the
        motion due to Movement.

        :return: (M, 3) ndarray of floats

            This is a ndarray containing the apparent velocity (in geometry axes,
            observed from the Earth frame) at the center point of each bound
            RingVortex's left leg due to any motion defined in Movement. Its units
            are in meters per second. If the current time step is the first time
            step, these velocities will all be all zeros.
        """
        # Check if this is the current time step. If so, return all zeros.
        if self._current_step < 1:
            return np.zeros((self.num_panels, 3), dtype=float)

        return -(self.stackCblvpl_G_Cg - self._lastStackCblvpl_G_Cg) / self.delta_time

    def _finalize_loads(self):
        """For cases with static geometry, this function finds the final loads and
        load coefficients for each of the SteadyProblem's Airplanes. For cases with
        variable geometry, it finds the final cycle-averaged and
        cycle-root-mean-squared loads and load coefficients for each of the
        SteadyProblem's Airplanes.

        :return: None
        """
        # Get this solver's time step characteristics. Note that the first time step (
        # time step 0), occurs at 0 seconds.
        num_steps_to_average = self.num_steps - self._first_averaging_step

        # Determine if this SteadyProblem's geometry is static or variable.
        static = self.unsteady_problem.movement.static

        # Initialize ndarrays to hold each Airplane's loads and load coefficients at
        # each of the time steps that calculated the loads.
        forces_W = np.zeros((self.num_airplanes, 3, num_steps_to_average), dtype=float)
        force_coefficients_W = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average), dtype=float
        )
        moments_W_Cg = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average), dtype=float
        )
        moment_coefficients_W_Cg = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average), dtype=float
        )

        # Initialize a variable to track position in the loads ndarrays.
        results_step = 0

        # Iterate through the time steps with loads and add the loads to their
        # respective ndarrays.
        for step in range(self._first_averaging_step, self.num_steps):

            # Get the Airplanes from the SteadyProblem at this time step.
            these_airplanes = self.steady_problems[step].airplanes

            # Iterate through this time step's Airplanes.
            airplane: geometry.airplane.Airplane
            for airplane_id, airplane in enumerate(these_airplanes):
                forces_W[airplane_id, :, results_step] = airplane.forces_W
                force_coefficients_W[airplane_id, :, results_step] = (
                    airplane.forceCoefficients_W
                )
                moments_W_Cg[airplane_id, :, results_step] = airplane.moments_W_Cg
                moment_coefficients_W_Cg[airplane_id, :, results_step] = (
                    airplane.momentCoefficients_W_Cg
                )

            results_step += 1

        # For each Airplane, calculate and then save the final or cycle-averaged and
        # RMS loads and load coefficients.
        airplane: geometry.airplane.Airplane
        for airplane_id, airplane in enumerate(self.steady_problems[0].airplanes):
            if static:
                self.unsteady_problem.finalForces_W.append(forces_W[airplane_id, :, -1])
                self.unsteady_problem.finalForceCoefficients_W.append(
                    force_coefficients_W[airplane_id, :, -1]
                )
                self.unsteady_problem.finalMoments_W_Cg.append(
                    moments_W_Cg[airplane_id, :, -1]
                )
                self.unsteady_problem.finalMomentCoefficients_W_Cg.append(
                    moment_coefficients_W_Cg[airplane_id, :, -1]
                )
            else:
                self.unsteady_problem.finalMeanForces_W.append(
                    np.mean(forces_W[airplane_id], axis=-1)
                )
                self.unsteady_problem.finalMeanForceCoefficients_W.append(
                    np.mean(force_coefficients_W[airplane_id], axis=-1)
                )
                self.unsteady_problem.finalMeanMoments_W_Cg.append(
                    np.mean(moments_W_Cg[airplane_id], axis=-1)
                )
                self.unsteady_problem.finalMeanMomentCoefficients_W_Cg.append(
                    np.mean(moment_coefficients_W_Cg[airplane_id], axis=-1)
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
                self.unsteady_problem.finalRmsMoments_W_Cg.append(
                    np.sqrt(
                        np.mean(
                            np.square(moments_W_Cg[airplane_id]),
                            axis=-1,
                        )
                    )
                )
                self.unsteady_problem.finalRmsMomentCoefficients_W_Cg.append(
                    np.sqrt(
                        np.mean(
                            np.square(moment_coefficients_W_Cg[airplane_id]),
                            axis=-1,
                        )
                    )
                )
