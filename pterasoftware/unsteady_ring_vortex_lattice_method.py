# NOTE: I haven't yet started refactoring this module.
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

from . import aerodynamics
from . import functions


# NOTE: I haven't yet started refactoring this class.
class UnsteadyRingVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses an unsteady ring vortex lattice method.

    This class contains the following public methods:
        run: This method runs the solver on the unsteady problem.

        initialize_panel_vortices: This method calculates the locations of the
        airplanes' bound vortex vertices, and then initializes their panels' bound
        vortices.

        collapse_geometry: This method converts attributes of the problem's geometry
        into 1D ndarrays. This facilitates vectorization, which speeds up the solver.

        calculate_wing_wing_influences: This method finds the matrix of wing-wing
        influence coefficients associated with the airplanes' geometries.

        calculate_freestream_wing_influences: This method finds the vector of
        freestream-wing influences associated with the problem at this time step.

        calculate_wake_wing_influences: This method finds the vector of the wake-wing
        influences associated with the problem at this time step.

        calculate_vortex_strengths: This method solves for each panel's vortex
        strength.

        calculate_solution_velocity: This function takes in a group of points. At
        every point, it finds the induced velocity due to every vortex and the
        freestream velocity.

        calculate_forces_and_moments: This method finds the forces and moments.

        populate_next_airplanes_wake: This method updates the next time step's
        airplanes' wakes.

        populate_next_airplanes_wake_vortex_vertices: This method populates the
        locations of the next airplanes' wake vortex vertices.

        populate_next_airplanes_wake_vortices: This method populates the locations of
        the next airplanes' wake vortices.

        calculate_current_flapping_velocities_at_collocation_points: This method
        finds the apparent flow velocity due to flapping at the centers of the
        current airplanes' collocation points.

        calculate_current_flapping_velocities_at_right_leg_centers: This method finds
        the apparent flow velocity due to flapping at the centers of the current
        airplanes' bound ring vortices' right legs.

        calculate_current_flapping_velocities_at_front_leg_centers: This method finds
        the apparent flow velocity due to flapping at the centers of the current
        airplanes' bound ring vortices' front legs.

        calculate_current_flapping_velocities_at_left_leg_centers: This method finds
        the apparent flow velocity due to flapping at the centers of the current
        airplanes' bound ring vortices' left legs.

        finalize_forces_and_moments: For cases with static geometry,
        this function finds the final loads and load coefficients for each of the
        problem's airplanes. For cases with variable geometry, it finds the final
        cycle-averaged and cycle-root-mean-squared loads and load coefficients.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I haven't yet started refactoring this method.
    def __init__(self, unsteady_problem):
        """This is the initialization method.

        :param unsteady_problem: UnsteadyProblem
            This is the unsteady problem to be solved.
        :return: None
        """
        self.unsteady_problem = unsteady_problem
        self.num_steps = self.unsteady_problem.num_steps
        self.delta_time = self.unsteady_problem.delta_time
        self.steady_problems = self.unsteady_problem.steady_problems
        self.first_results_step = self.unsteady_problem.first_results_step
        self.first_averaging_step = self.unsteady_problem.first_averaging_step
        self.num_airplanes = len(self.steady_problems[0].airplanes)

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.current_step = None
        self.current_airplanes = None
        self.current_operating_point = None
        self.current_freestream_velocity_G = None
        self.current_wing_wing_influences = None
        self.current_freestream_wing_influences = None
        self.current_wake_wing_influences = None
        self.current_vortex_strengths = None
        self.streamline_points = None

        # Initialize attributes to hold geometric data that pertains to this problem.
        self.panels = None
        self.panel_normal_directions = None
        self.panel_areas = None
        self.panel_collocation_points = None
        self.panel_back_right_vortex_vertices = None
        self.panel_front_right_vortex_vertices = None
        self.panel_front_left_vortex_vertices = None
        self.panel_back_left_vortex_vertices = None
        self.panel_right_vortex_centers = None
        self.panel_right_vortex_vectors = None
        self.panel_front_vortex_centers = None
        self.panel_front_vortex_vectors = None
        self.panel_left_vortex_centers = None
        self.panel_left_vortex_vectors = None
        self.panel_back_vortex_centers = None
        self.panel_back_vortex_vectors = None
        self.panel_moment_references = None
        self.seed_points = None

        # Initialize variables to hold aerodynamic data that pertains details about
        # this panel's location on its wing.
        self.panel_is_trailing_edge = None
        self.panel_is_leading_edge = None
        self.panel_is_left_edge = None
        self.panel_is_right_edge = None

        # Initialize variables to hold aerodynamic data that pertains to this
        # problem's last time step.
        self.last_panel_collocation_points = None
        self.last_panel_vortex_strengths = None
        self.last_panel_back_right_vortex_vertices = None
        self.last_panel_front_right_vortex_vertices = None
        self.last_panel_front_left_vortex_vertices = None
        self.last_panel_back_left_vortex_vertices = None
        self.last_panel_right_vortex_centers = None
        self.last_panel_front_vortex_centers = None
        self.last_panel_left_vortex_centers = None
        self.last_panel_back_vortex_centers = None

        # Initialize lists to store aerodynamic data about the wake at each time
        # step. These attributes are used by several functions in the output module
        # to animate the wake.
        self.num_wake_ring_vortices_list = []
        self.wake_ring_vortex_strengths_list = []
        self.wake_ring_vortex_ages_list = []
        self.wake_ring_vortex_front_right_vertices_list = []
        self.wake_ring_vortex_front_left_vertices_list = []
        self.wake_ring_vortex_back_left_vertices_list = []
        self.wake_ring_vortex_back_right_vertices_list = []

        # Initialize variables to hold aerodynamic data that pertains to the wake at
        # the current time step.
        self.current_wake_ring_vortex_strengths = None
        self.current_wake_ring_vortex_front_right_vertices = None
        self.current_wake_ring_vortex_front_left_vertices = None
        self.current_wake_ring_vortex_back_left_vertices = None
        self.current_wake_ring_vortex_back_right_vertices = None
        self.current_wake_ring_vortex_ages = None
        self.num_panels = None

    # NOTE: I haven't yet started refactoring this method.
    def run(
        self,
        logging_level="Warning",
        prescribed_wake=True,
        calculate_streamlines=True,
    ):
        """This method runs the solver on the unsteady problem.

        :param logging_level: str, optional
            This parameter determines the detail of information that the solver's
            logger will output while running. The options are, in order of detail and
            severity, "Debug", "Info", "Warning", "Error", "Critical". The default
            value is "Warning".
        :param prescribed_wake: Bool, optional
            This parameter determines if the solver uses a prescribed wake model. If
            false it will use a free-wake, which may be more accurate but will make
            the solver significantly slower. The default is True.
        :param calculate_streamlines: Bool, optional
            This parameter determines if the solver uses calculates streamlines
            emanating from the back of the wing after running the solver.
        :return: None
        """
        # Configure the problem's logger.
        logging_level_value = functions.convert_logging_level_name_to_value(
            logging_level
        )
        logging.basicConfig(level=logging_level_value)

        # The following loop iterates through the steps to populate currently empty
        # attributes with lists of pre-allocated arrays. During the simulation,
        # these arrays will be filled with data that describe the wake. Using this
        # method eliminates the need for computationally expensive on-the-fly
        # allocation and object copying.
        for step in range(self.num_steps):
            this_problem = self.steady_problems[step]
            these_airplanes = this_problem.airplanes

            # Loop through this step's airplanes to create a list of their wings.
            these_wings = []
            for airplane in these_airplanes:
                these_wings.append(airplane.wings)

            # Iterate through the wings to get the total number of spanwise panels.
            this_num_spanwise_panels = 0
            for this_wing_set in these_wings:
                for this_wing in this_wing_set:
                    this_num_spanwise_panels += this_wing.num_spanwise_panels

            # The number of wake vortices is the step number multiplied by the number
            # of spanwise panels. This works because the first step number is zero.
            this_num_wake_ring_vortices = step * this_num_spanwise_panels

            # Allocate the arrays for this step.
            this_wake_ring_vortex_strengths = np.zeros(this_num_wake_ring_vortices)
            this_wake_ring_vortex_ages = np.zeros(this_num_wake_ring_vortices)
            this_wake_ring_vortex_front_right_vertices = np.zeros(
                (this_num_wake_ring_vortices, 3)
            )
            this_wake_ring_vortex_front_left_vertices = np.zeros(
                (this_num_wake_ring_vortices, 3)
            )
            this_wake_ring_vortex_back_left_vertices = np.zeros(
                (this_num_wake_ring_vortices, 3)
            )
            this_wake_ring_vortex_back_right_vertices = np.zeros(
                (this_num_wake_ring_vortices, 3)
            )

            # Append this step's arrays to the list of arrays.
            self.num_wake_ring_vortices_list.append(this_num_wake_ring_vortices)
            self.wake_ring_vortex_strengths_list.append(this_wake_ring_vortex_strengths)
            self.wake_ring_vortex_ages_list.append(this_wake_ring_vortex_ages)
            self.wake_ring_vortex_front_right_vertices_list.append(
                this_wake_ring_vortex_front_right_vertices
            )
            self.wake_ring_vortex_front_left_vertices_list.append(
                this_wake_ring_vortex_front_left_vertices
            )
            self.wake_ring_vortex_back_left_vertices_list.append(
                this_wake_ring_vortex_back_left_vertices
            )
            self.wake_ring_vortex_back_right_vertices_list.append(
                this_wake_ring_vortex_back_right_vertices
            )

        # The following loop attempts to predict how much time each step will take,
        # relative to the other steps. This data will be used to generate estimates
        # of how much longer a simulation will take, and create a smoothly advancing
        # progress bar.

        # Initialize list that will hold the approximate, relative times. This has
        # one more element than the number of steps, because I will also use the
        # progress bar during the problem initialization.
        approx_times = np.zeros(self.num_steps + 1)
        for step in range(1, self.num_steps):
            this_problem = self.steady_problems[step]
            these_airplanes = this_problem.airplanes

            # Iterate through this step's airplanes to get the total number of wing
            # panels.
            num_wing_panels = 0
            for airplane in these_airplanes:
                num_wing_panels += airplane.num_panels

            # Calculate the total number of vortices analyzed during this step.
            num_wing_ring_vortices = num_wing_panels
            num_wake_ring_vortices = self.num_wake_ring_vortices_list[step]
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

            # Initialize all the airplanes' panels' vortices.
            logging.info("Initializing all airplanes' panel vortices.")
            self.initialize_panel_vortices()

            # Update the progress bar based on the initialization step's predicted
            # approximate, relative computing time.
            bar.update(n=approx_times[0])

            # Iterate through the time steps.
            for step in range(self.num_steps):

                # Save attributes to hold the current step, airplane, and operating
                # point.
                self.current_step = step
                self.current_airplanes = self.steady_problems[
                    self.current_step
                ].airplanes
                self.current_operating_point = self.steady_problems[
                    self.current_step
                ].operating_point
                self.current_freestream_velocity_G = (
                    self.current_operating_point.vInf_G__E
                )
                logging.info(
                    "Beginning time step "
                    + str(self.current_step)
                    + " out of "
                    + str(self.num_steps - 1)
                    + "."
                )

                # Calculate the number of panels for this time step.
                self.num_panels = 0
                for airplane in self.current_airplanes:
                    self.num_panels += airplane.num_panels

                # Initialize attributes to hold aerodynamic data that pertains to this
                # problem.
                self.current_wing_wing_influences = np.zeros(
                    (self.num_panels, self.num_panels)
                )
                self.current_freestream_velocity_G = (
                    self.current_operating_point.vInf_G__E
                )
                self.current_freestream_wing_influences = np.zeros(self.num_panels)
                self.current_wake_wing_influences = np.zeros(self.num_panels)
                self.current_vortex_strengths = np.ones(self.num_panels)

                # Initialize attributes to hold geometric data that pertains to this
                # problem.
                self.panels = np.empty(self.num_panels, dtype=object)
                self.panel_normal_directions = np.zeros((self.num_panels, 3))
                self.panel_areas = np.zeros(self.num_panels)
                self.panel_collocation_points = np.zeros((self.num_panels, 3))
                self.panel_back_right_vortex_vertices = np.zeros((self.num_panels, 3))
                self.panel_front_right_vortex_vertices = np.zeros((self.num_panels, 3))
                self.panel_front_left_vortex_vertices = np.zeros((self.num_panels, 3))
                self.panel_back_left_vortex_vertices = np.zeros((self.num_panels, 3))
                self.panel_right_vortex_centers = np.zeros((self.num_panels, 3))
                self.panel_right_vortex_vectors = np.zeros((self.num_panels, 3))
                self.panel_front_vortex_centers = np.zeros((self.num_panels, 3))
                self.panel_front_vortex_vectors = np.zeros((self.num_panels, 3))
                self.panel_left_vortex_centers = np.zeros((self.num_panels, 3))
                self.panel_left_vortex_vectors = np.zeros((self.num_panels, 3))
                self.panel_back_vortex_centers = np.zeros((self.num_panels, 3))
                self.panel_back_vortex_vectors = np.zeros((self.num_panels, 3))
                self.seed_points = np.zeros((0, 3))
                self.panel_moment_references = np.zeros((self.num_panels, 3))

                # Initialize variables to hold details about each panel's location on
                # its wing.
                self.panel_is_trailing_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_leading_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_left_edge = np.zeros(self.num_panels, dtype=bool)
                self.panel_is_right_edge = np.zeros(self.num_panels, dtype=bool)

                # Initialize variables to hold details about the last time step's
                # airplanes' panels.
                self.last_panel_collocation_points = np.zeros((self.num_panels, 3))
                self.last_panel_vortex_strengths = np.zeros(self.num_panels)
                self.last_panel_back_right_vortex_vertices = np.zeros(
                    (self.num_panels, 3)
                )
                self.last_panel_front_right_vortex_vertices = np.zeros(
                    (self.num_panels, 3)
                )
                self.last_panel_front_left_vortex_vertices = np.zeros(
                    (self.num_panels, 3)
                )
                self.last_panel_back_left_vortex_vertices = np.zeros(
                    (self.num_panels, 3)
                )
                self.last_panel_right_vortex_centers = np.zeros((self.num_panels, 3))
                self.last_panel_front_vortex_centers = np.zeros((self.num_panels, 3))
                self.last_panel_left_vortex_centers = np.zeros((self.num_panels, 3))
                self.last_panel_back_vortex_centers = np.zeros((self.num_panels, 3))

                # Get the pre-allocated (but still all zero) arrays of wake
                # information that are associated with this time step.
                self.current_wake_ring_vortex_strengths = (
                    self.wake_ring_vortex_strengths_list[step]
                )
                self.current_wake_ring_vortex_ages = self.wake_ring_vortex_ages_list[
                    step
                ]
                self.current_wake_ring_vortex_front_right_vertices = (
                    self.wake_ring_vortex_front_right_vertices_list[step]
                )
                self.current_wake_ring_vortex_front_left_vertices = (
                    self.wake_ring_vortex_front_left_vertices_list[step]
                )
                self.current_wake_ring_vortex_back_left_vertices = (
                    self.wake_ring_vortex_back_left_vertices_list[step]
                )
                self.current_wake_ring_vortex_back_right_vertices = (
                    self.wake_ring_vortex_back_right_vertices_list[step]
                )

                # Collapse this problem's geometry matrices into 1D arrays of
                # attributes.
                logging.info("Collapsing the geometry.")
                self.collapse_geometry()

                # Find the matrix of wing-wing influence coefficients associated with
                # the airplanes' geometries.
                logging.info("Calculating the wing-wing influences.")
                self.calculate_wing_wing_influences()

                # Find the vector of freestream-wing influence coefficients associated
                # with this problem.
                logging.info("Calculating the freestream-wing influences.")
                self.calculate_freestream_wing_influences()

                # Find the vector of wake-wing influence coefficients associated with
                # this problem.
                logging.info("Calculating the wake-wing influences.")
                self.calculate_wake_wing_influences()

                # Solve for each panel's vortex strength.
                logging.info("Calculating vortex strengths.")
                self.calculate_vortex_strengths()

                # Solve for the forces and moments on each panel.
                if self.current_step >= self.first_results_step:
                    logging.info("Calculating forces and moments.")
                    self.calculate_forces_and_moments()

                # Solve for the forces and moments on each panel.
                logging.info("Shedding wake vortices.")
                self.populate_next_airplanes_wake(prescribed_wake=prescribed_wake)

                # Update the progress bar based on this step's predicted approximate,
                # relative computing time.
                bar.update(n=approx_times[step + 1])

            logging.info("Calculating averaged or final forces and moments.")
            self.finalize_forces_and_moments()

        # Solve for the location of the streamlines if requested.
        if calculate_streamlines:
            logging.info("Calculating streamlines.")
            functions.calculate_streamlines(self)

    # NOTE: I haven't yet started refactoring this method.
    def initialize_panel_vortices(self):
        """This method calculates the locations of the airplanes' bound vortex
        vertices, and then initializes their panels' bound vortices.

        Every panel has a ring vortex, which is a quadrangle whose front vortex leg
        is at the panel's quarter chord. The left and right vortex legs run along the
        panel's left and right legs. If the panel is not along the trailing edge,
        they extend backwards and meet the back vortex leg at a length of one quarter
        of the rear panel's chord back from the rear panel's front leg. Otherwise,
        they extend back backwards and meet the back vortex leg at a length of one
        quarter of the current panel's chord back from the current panel's back leg.

        :return: None
        """
        # Iterate through all the steady problem objects.
        for steady_problem in self.steady_problems:

            # Get the freestream velocity at this time step's problem.
            this_freestream_velocity_G = steady_problem.operating_point.vInf_G__E

            # Iterate through this problem's airplanes' wings.
            for airplane in steady_problem.airplanes:
                for wing in airplane.wings:

                    # Iterate through the wing's chordwise and spanwise positions.
                    for chordwise_position in range(wing.num_chordwise_panels):
                        for spanwise_position in range(wing.num_spanwise_panels):

                            # Get the panel object from the wing's list of panels.
                            panel = wing.panels[chordwise_position, spanwise_position]

                            # Find the location of the panel's front left and right
                            # vortex vertices.
                            front_left_vortex_vertex = panel.Flbvp_G_Cg
                            front_right_vortex_vertex = panel.Frbvp_G_Cg

                            # Define the back left and right vortex vertices based on
                            # whether the panel is along the trailing edge or not.
                            if not panel.is_trailing_edge:
                                next_chordwise_panel = wing.panels[
                                    chordwise_position + 1, spanwise_position
                                ]
                                back_left_vortex_vertex = (
                                    next_chordwise_panel.Flbvp_G_Cg
                                )
                                back_right_vortex_vertex = (
                                    next_chordwise_panel.Frbvp_G_Cg
                                )
                            else:
                                # As these vertices are directly behind the trailing
                                # edge, they are spaced back from their panel's
                                # vertex by one quarter the distance traveled during
                                # a time step. This is to more accurately predict
                                # drag. More information can be found on pages 37-39
                                # of "Modeling of aerodynamic forces in flapping
                                # flight with the Unsteady Vortex Lattice Method" by
                                # Thomas Lambert.
                                back_left_vortex_vertex = (
                                    front_left_vortex_vertex
                                    + (panel.Blpp_G_Cg - panel.Flpp_G_Cg)
                                    + this_freestream_velocity_G
                                    * self.delta_time
                                    * 0.25
                                )
                                back_right_vortex_vertex = (
                                    front_right_vortex_vertex
                                    + (panel.Brpp_G_Cg - panel.Frpp_G_Cg)
                                    + this_freestream_velocity_G
                                    * self.delta_time
                                    * 0.25
                                )

                            # Initialize the panel's ring vortex.
                            panel.ring_vortex = aerodynamics.RingVortex(
                                Flrvp_G_Cg=front_left_vortex_vertex,
                                Frrvp_G_Cg=front_right_vortex_vertex,
                                Blrvp_G_Cg=back_left_vortex_vertex,
                                Brrvp_G_Cg=back_right_vortex_vertex,
                                strength=None,
                            )

    # NOTE: I haven't yet started refactoring this method.
    def collapse_geometry(self):
        """This method converts attributes of the problem's geometry into 1D
        ndarrays. This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize a variable to hold the global position of the panel as we
        # iterate through them.
        global_panel_position = 0
        global_wake_ring_vortex_position = 0

        # Iterate through the current step's airplanes' wings.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:

                # Convert this wing's 2D array of panels into a 1D array.
                panels = np.ravel(wing.panels)
                wake_ring_vortices = np.ravel(wing.wake_ring_vortices)

                # Iterate through the 1D array of this wing's panels.
                for panel in panels:
                    # Update the solver's list of attributes with this panel's
                    # attributes.
                    functions.update_ring_vortex_solvers_panel_attributes(
                        solver=self,
                        global_panel_position=global_panel_position,
                        panel=panel,
                        airplane=airplane,
                    )

                    # Increment the global panel position.
                    global_panel_position += 1

                for wake_ring_vortex in wake_ring_vortices:
                    self.current_wake_ring_vortex_strengths[
                        global_wake_ring_vortex_position
                    ] = wake_ring_vortex.strength
                    self.current_wake_ring_vortex_ages[
                        global_wake_ring_vortex_position
                    ] = wake_ring_vortex.age
                    self.current_wake_ring_vortex_front_right_vertices[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Frpp_G_Cg
                    self.current_wake_ring_vortex_front_left_vertices[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Flpp_G_Cg
                    self.current_wake_ring_vortex_back_left_vertices[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Blpp_G_Cg
                    self.current_wake_ring_vortex_back_right_vertices[
                        global_wake_ring_vortex_position, :
                    ] = wake_ring_vortex.Brpp_G_Cg

                    global_wake_ring_vortex_position += 1

        # Initialize a variable to hold the global position of the panel as we
        # iterate through them.
        global_panel_position = 0

        if self.current_step > 0:

            last_airplanes = self.steady_problems[self.current_step - 1].airplanes

            # Iterate through the last step's airplanes' wings.
            for last_airplane in last_airplanes:
                for wing in last_airplane.wings:

                    # Convert this wing's 2D array of panels into a 1D array.
                    panels = np.ravel(wing.panels)

                    # Iterate through the 1D array of this wing's panels.
                    for panel in panels:
                        # Update the solver's list of attributes with this panel's
                        # attributes.
                        self.last_panel_collocation_points[global_panel_position, :] = (
                            panel.collocation_point
                        )
                        self.last_panel_vortex_strengths[global_panel_position] = (
                            panel.ring_vortex.strength
                        )
                        self.last_panel_back_right_vortex_vertices[
                            global_panel_position, :
                        ] = panel.ring_vortex.rightLeg_G.Slvp_G_Cg
                        self.last_panel_front_right_vortex_vertices[
                            global_panel_position, :
                        ] = panel.ring_vortex.rightLeg_G.Elvp_G_Cg
                        self.last_panel_front_left_vortex_vertices[
                            global_panel_position, :
                        ] = panel.ring_vortex.leftLeg_G.Slvp_G_Cg
                        self.last_panel_back_left_vortex_vertices[
                            global_panel_position, :
                        ] = panel.ring_vortex.leftLeg_G.Elvp_G_Cg
                        self.last_panel_right_vortex_centers[
                            global_panel_position, :
                        ] = panel.ring_vortex.rightLeg_G.Crvp_G_Cg
                        self.last_panel_front_vortex_centers[
                            global_panel_position, :
                        ] = panel.ring_vortex.frontLeg_G.Crvp_G_Cg
                        self.last_panel_left_vortex_centers[
                            global_panel_position, :
                        ] = panel.ring_vortex.leftLeg_G.Crvp_G_Cg
                        self.last_panel_back_vortex_centers[
                            global_panel_position, :
                        ] = panel.ring_vortex.backLeg_G.Crvp_G_Cg

                        # Increment the global panel position.
                        global_panel_position += 1

    # NOTE: I haven't yet started refactoring this method.
    def calculate_wing_wing_influences(self):
        """This method finds the matrix of wing-wing influence coefficients
        associated with the airplanes' geometries.

        :return: None
        """
        # Find the matrix of normalized velocities induced at every panel's
        # collocation point by every panel's ring vortex. The answer is normalized
        # because the solver's vortex strength list was initialized to all ones. This
        # will be updated once the correct vortex strengths are calculated.
        total_influences = aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.panel_collocation_points,
            stackBrrvp_G_Cg=self.panel_back_right_vortex_vertices,
            stackFrrvp_G_Cg=self.panel_front_right_vortex_vertices,
            stackFlrvp_G_Cg=self.panel_front_left_vortex_vertices,
            stackBlrvp_G_Cg=self.panel_back_left_vortex_vertices,
            strengths=self.current_vortex_strengths,
            ages=None,
            nu=self.current_operating_point.nu,
        )

        # Take the batch dot product of the normalized velocities with each panel's
        # normal direction. This is now the problem's matrix of wing-wing influence
        # coefficients.
        self.current_wing_wing_influences = np.einsum(
            "...k,...k->...",
            total_influences,
            np.expand_dims(self.panel_normal_directions, axis=1),
        )

    # NOTE: I haven't yet started refactoring this method.
    def calculate_freestream_wing_influences(self):
        """This method finds the vector of freestream-wing influence coefficients
        associated with this problem.

        Note: This method also includes the influence due to flapping at every
        collocation point in the freestream-wing influence vector.

        :return: None
        """
        # Find the normal components of the freestream velocity on every panel by
        # taking a batch dot product.
        freestream_influences = np.einsum(
            "ij,j->i",
            self.panel_normal_directions,
            self.current_freestream_velocity_G,
        )

        # Get the current flapping velocities at every collocation point.
        current_flapping_velocities_at_collocation_points = (
            self.calculate_current_flapping_velocities_at_collocation_points()
        )

        # Find the normal components of every panel's flapping velocities at their
        # collocation points by taking a batch dot product.
        flapping_influences = np.einsum(
            "ij,ij->i",
            self.panel_normal_directions,
            current_flapping_velocities_at_collocation_points,
        )

        # Calculate the total current freestream-wing influences by summing the
        # freestream influences and the flapping influences.
        self.current_freestream_wing_influences = (
            freestream_influences + flapping_influences
        )

    # NOTE: I haven't yet started refactoring this method.
    def calculate_wake_wing_influences(self):
        """This method finds the vector of the wake-wing influences associated with
        the problem at this time step.

        Note: If the current time step is the first time step, no wake has yet been
        shed, and this method will set the current wake-wing influence vector to all
        zeros.

        :return: None
        """
        # Check if this time step is not the first time step.
        if self.current_step > 0:

            # Get the wake induced velocities. This is a (M x 3) array with the x, y,
            # and z components of the velocity induced by the entire wake at each of
            # the M panels.
            velocities_from_wake = aerodynamics.collapsed_velocities_from_ring_vortices(
                stackP_G_Cg=self.panel_collocation_points,
                stackBrrvp_G_Cg=self.current_wake_ring_vortex_back_right_vertices,
                stackFrrvp_G_Cg=self.current_wake_ring_vortex_front_right_vertices,
                stackFlrvp_G_Cg=self.current_wake_ring_vortex_front_left_vertices,
                stackBlrvp_G_Cg=self.current_wake_ring_vortex_back_left_vertices,
                strengths=self.current_wake_ring_vortex_strengths,
                ages=self.current_wake_ring_vortex_ages,
                nu=self.current_operating_point.nu,
            )

            # Set the current wake-wing influences to the normal component of the
            # wake induced velocities at each panel.
            self.current_wake_wing_influences = np.einsum(
                "ij,ij->i", velocities_from_wake, self.panel_normal_directions
            )

        else:

            # If this is the first time step, set the current wake-wing influences to
            # zero everywhere, as there is no wake yet.
            self.current_wake_wing_influences = np.zeros(self.num_panels)

    # NOTE: I haven't yet started refactoring this method.
    def calculate_vortex_strengths(self):
        """This method solves for each panel's vortex strength.

        :return: None
        """
        # Solve for the strength of each panel's vortex.
        self.current_vortex_strengths = np.linalg.solve(
            self.current_wing_wing_influences,
            -self.current_wake_wing_influences
            - self.current_freestream_wing_influences,
        )

        # Iterate through the panels and update their vortex strengths.
        for panel_num in range(self.panels.size):
            # Get the panel at this location.
            panel = self.panels[panel_num]

            # Update this panel's ring vortex strength.
            panel.ring_vortex.update_strength(self.current_vortex_strengths[panel_num])

    # NOTE: I haven't yet started refactoring this method.
    def calculate_solution_velocity(self, points):
        """This function takes in a group of points. At every point, it finds the
        induced velocity due to every vortex and the freestream velocity.

        Note: The velocity calculated by this method is in geometry axes. Also,
        this method assumes that the correct vortex strengths have already been
        calculated. This method also does not include the velocity due to flapping at
        any of the points provided, as it has no way of knowing if any of the points
        lie on panels.

        This method uses vectorization, and therefore is much faster for batch
        operations than using the vortex objects' class methods for calculating
        induced velocity.

        :param points: 2D array of floats
            This variable is an array of shape (N x 3), where N is the number of
            points. Each row contains the x, y, and z float coordinates of that
            point's position in meters.
        :return: 2D array of floats
            The output is the summed effects from every vortex, and from the
            freestream on a given point. The result will be of shape (N x 3),
            where each row identifies the velocity at a point. The results units are
            meters per second.
        """
        # Find the vector of velocities induced at every point by every panel's ring
        # vortex. The effect of every ring vortex on each point will be summed.
        velocities_from_wings = aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=points,
            stackBrrvp_G_Cg=self.panel_back_right_vortex_vertices,
            stackFrrvp_G_Cg=self.panel_front_right_vortex_vertices,
            stackFlrvp_G_Cg=self.panel_front_left_vortex_vertices,
            stackBlrvp_G_Cg=self.panel_back_left_vortex_vertices,
            strengths=self.current_vortex_strengths,
            ages=None,
            nu=self.current_operating_point.nu,
        )

        # Find the vector of velocities induced at every point by every wake ring
        # vortex. The effect of every wake ring vortex on each point will be summed.
        velocities_from_wake = aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=points,
            stackBrrvp_G_Cg=self.current_wake_ring_vortex_back_right_vertices,
            stackFrrvp_G_Cg=self.current_wake_ring_vortex_front_right_vertices,
            stackFlrvp_G_Cg=self.current_wake_ring_vortex_front_left_vertices,
            stackBlrvp_G_Cg=self.current_wake_ring_vortex_back_left_vertices,
            strengths=self.current_wake_ring_vortex_strengths,
            ages=self.current_wake_ring_vortex_ages,
            nu=self.current_operating_point.nu,
        )

        # Find the total influence of the vortices, which is the sum of the influence
        # due to the bound ring vortices and the wake ring vortices.
        total_vortex_velocities = velocities_from_wings + velocities_from_wake

        # Calculate and return the sum of the velocities induced by the vortices and
        # freestream at every point.
        return total_vortex_velocities + self.current_freestream_velocity_G

    # NOTE: I haven't yet started refactoring this method.
    def calculate_forces_and_moments(self):
        """This method finds the forces and moments calculated from the near field.

        Citation: This method uses logic described on pages 9-11 of "Modeling of
        aerodynamic forces in flapping flight with the Unsteady Vortex Lattice
        Method" by Thomas Lambert.

        Note: The forces and moments calculated are in geometry axes. The moment is
        about the airplanes' reference points, which should be at their centers of
        gravity. The units are Newtons and Newton-meters.

        :return: None
        """
        # Initialize a variable to hold the global panel position as the panels are
        # iterated through.
        global_panel_position = 0

        # Initialize three lists of variables, which will hold the effective strength
        # of the line vortices comprising each panel's ring vortex.
        effective_right_vortex_line_strengths = np.zeros(self.num_panels)
        effective_front_vortex_line_strengths = np.zeros(self.num_panels)
        effective_left_vortex_line_strengths = np.zeros(self.num_panels)

        # Iterate through each of this step's airplanes' wings.
        for airplane in self.current_airplanes:
            for wing in airplane.wings:

                # Convert this wing's 2D array of panels into a 1D array.
                panels = np.ravel(wing.panels)

                # Iterate through this wing's 1D array panels.
                for panel in panels:

                    # Check if this panel is on its wing's right edge.
                    if panel.is_right_edge:

                        # Change the effective right vortex line strength from zero
                        # to this panel's ring vortex's strength.
                        effective_right_vortex_line_strengths[global_panel_position] = (
                            self.current_vortex_strengths[global_panel_position]
                        )

                    else:

                        # Get the panel directly to the right of this panel.
                        panel_to_right = wing.panels[
                            panel.local_chordwise_position,
                            panel.local_spanwise_position + 1,
                        ]

                        # Change the effective right vortex line strength from zero
                        # to the difference between this panel's ring vortex's
                        # strength, and the ring vortex strength of the panel to the
                        # right of it.
                        effective_right_vortex_line_strengths[global_panel_position] = (
                            self.current_vortex_strengths[global_panel_position]
                            - panel_to_right.ring_vortex.strength
                        )

                    # Check if this panel is on its wing's leading edge.
                    if panel.is_leading_edge:

                        # Change the effective front vortex line strength from zero
                        # to this panel's ring vortex's strength.
                        effective_front_vortex_line_strengths[global_panel_position] = (
                            self.current_vortex_strengths[global_panel_position]
                        )
                    else:

                        # Get the panel directly in front of this panel.
                        panel_to_front = wing.panels[
                            panel.local_chordwise_position - 1,
                            panel.local_spanwise_position,
                        ]

                        # Change the effective front vortex line strength from zero
                        # to the difference between this panel's ring vortex's
                        # strength, and the ring vortex strength of the panel in
                        # front of it.
                        effective_front_vortex_line_strengths[global_panel_position] = (
                            self.current_vortex_strengths[global_panel_position]
                            - panel_to_front.ring_vortex.strength
                        )

                    # Check if this panel is on its wing's left edge.
                    if panel.is_left_edge:

                        # Change the effective left vortex line strength from zero to
                        # this panel's ring vortex's strength.
                        effective_left_vortex_line_strengths[global_panel_position] = (
                            self.current_vortex_strengths[global_panel_position]
                        )
                    else:

                        # Get the panel directly to the left of this panel.
                        panel_to_left = wing.panels[
                            panel.local_chordwise_position,
                            panel.local_spanwise_position - 1,
                        ]

                        # Change the effective left vortex line strength from zero to
                        # the difference between this panel's ring vortex's strength,
                        # and the ring vortex strength of the panel to the left of it.
                        effective_left_vortex_line_strengths[global_panel_position] = (
                            self.current_vortex_strengths[global_panel_position]
                            - panel_to_left.ring_vortex.strength
                        )

                    # Increment the global panel position.
                    global_panel_position += 1

        # Calculate the solution velocities at the centers of the panel's front leg,
        # left leg, and right leg.
        velocities_at_ring_vortex_front_leg_centers = (
            self.calculate_solution_velocity(points=self.panel_front_vortex_centers)
            + self.calculate_current_flapping_velocities_at_front_leg_centers()
        )
        velocities_at_ring_vortex_left_leg_centers = (
            self.calculate_solution_velocity(points=self.panel_left_vortex_centers)
            + self.calculate_current_flapping_velocities_at_left_leg_centers()
        )
        velocities_at_ring_vortex_right_leg_centers = (
            self.calculate_solution_velocity(points=self.panel_right_vortex_centers)
            + self.calculate_current_flapping_velocities_at_right_leg_centers()
        )

        # Using the effective line vortex strengths, and the Kutta-Joukowski theorem
        # to find the force in geometry axes on the front leg, left leg,
        # and right leg. Also calculate the unsteady component of the force on each
        # panel, which is derived from the unsteady Bernoulli equation.
        forces_on_ring_vortex_right_legs_G = (
            self.current_operating_point.density
            * np.expand_dims(effective_right_vortex_line_strengths, axis=1)
            * functions.numba_1d_explicit_cross(
                velocities_at_ring_vortex_right_leg_centers,
                self.panel_right_vortex_vectors,
            )
        )
        forces_on_ring_vortex_front_legs_G = (
            self.current_operating_point.density
            * np.expand_dims(effective_front_vortex_line_strengths, axis=1)
            * functions.numba_1d_explicit_cross(
                velocities_at_ring_vortex_front_leg_centers,
                self.panel_front_vortex_vectors,
            )
        )
        forces_on_ring_vortex_left_legs_G = (
            self.current_operating_point.density
            * np.expand_dims(effective_left_vortex_line_strengths, axis=1)
            * functions.numba_1d_explicit_cross(
                velocities_at_ring_vortex_left_leg_centers,
                self.panel_left_vortex_vectors,
            )
        )
        unsteady_forces_G = (
            self.current_operating_point.density
            * np.expand_dims(
                (self.current_vortex_strengths - self.last_panel_vortex_strengths),
                axis=1,
            )
            * np.expand_dims(self.panel_areas, axis=1)
            * self.panel_normal_directions
            / self.delta_time
        )

        # Sum the forces on the legs, and the unsteady force, to calculate the total
        # force, in geometry axes, on each panel.
        forces_G = (
            forces_on_ring_vortex_front_legs_G
            + forces_on_ring_vortex_left_legs_G
            + forces_on_ring_vortex_right_legs_G
            + unsteady_forces_G
        )

        # Find the moment in geometry axes on the front leg, left leg,
        # and right leg. Also find the moment on each panel due to the unsteady force.
        moments_on_ring_vortex_front_legs_G = functions.numba_1d_explicit_cross(
            self.panel_front_vortex_centers - self.panel_moment_references,
            forces_on_ring_vortex_front_legs_G,
        )
        moments_on_ring_vortex_left_legs_G = functions.numba_1d_explicit_cross(
            self.panel_left_vortex_centers - self.panel_moment_references,
            forces_on_ring_vortex_left_legs_G,
        )
        moments_on_ring_vortex_right_legs_G = functions.numba_1d_explicit_cross(
            self.panel_right_vortex_centers - self.panel_moment_references,
            forces_on_ring_vortex_right_legs_G,
        )
        unsteady_moments_G_Cg = functions.numba_1d_explicit_cross(
            self.panel_collocation_points - self.panel_moment_references,
            unsteady_forces_G,
        )

        # Sum the moments on the legs, and the unsteady moment, to calculate the
        # total moment, in geometry axes, on each panel.
        moments_G_Cg = (
            moments_on_ring_vortex_front_legs_G
            + moments_on_ring_vortex_left_legs_G
            + moments_on_ring_vortex_right_legs_G
            + unsteady_moments_G_Cg
        )

        functions.process_unsteady_solver_loads(
            unsteady_solver=self, forces_G=forces_G, moments_G_Cg=moments_G_Cg
        )

    # NOTE: I haven't yet started refactoring this method.
    def populate_next_airplanes_wake(self, prescribed_wake=True):
        """This method updates the next time step's airplanes' wakes.

        :param prescribed_wake: Bool, optional
            This parameter determines if the solver uses a prescribed wake model. If
            false it will use a free-wake, which may be more accurate but will make
            the solver significantly slower. The default is True.
        :return: None
        """
        # Populate the locations of the next step's airplanes' wake vortex vertices:
        self.populate_next_airplanes_wake_vortex_vertices(
            prescribed_wake=prescribed_wake
        )

        # Populate the locations of the next step's airplanes' wake vortices.
        self.populate_next_airplanes_wake_vortices()

    # NOTE: I haven't yet started refactoring this method.
    def populate_next_airplanes_wake_vortex_vertices(self, prescribed_wake=True):
        """This method populates the locations of the next airplanes' wake vortex
        vertices.

        This method is not vectorized but its loops only consume 1.1% of the runtime,
        so I have kept it as is for increased readability.

        :param prescribed_wake: Bool, optional
            This parameter determines if the solver uses a prescribed wake model. If
            false it will use a free-wake, which may be more accurate but will make
            the solver significantly slower. The default is True.
        :return: None
        """
        # Check if this is not the last step.
        if self.current_step < self.num_steps - 1:

            # Get the next time step's airplane objects.
            next_airplanes = self.steady_problems[self.current_step + 1].airplanes

            # Get the current airplanes' combined number of wings.
            num_wings = 0
            for airplane in self.current_airplanes:
                num_wings += len(airplane.wings)

            # Iterate through this time step's airplanes' successor objects.
            for airplane_id, next_airplane in enumerate(next_airplanes):

                # Iterate through the next airplane's wings.
                for wing_id, next_wing in enumerate(next_airplane.wings):

                    # Get the wing objects at this position from the current airplane.
                    this_wing = self.current_airplanes[airplane_id].wings[wing_id]

                    # Check if this is the first step.
                    if self.current_step == 0:

                        # Get the current wing's number of chordwise and spanwise
                        # panels.
                        num_spanwise_panels = this_wing.num_spanwise_panels
                        num_chordwise_panels = this_wing.num_chordwise_panels

                        # Set the chordwise position to be at the trailing edge.
                        chordwise_position = num_chordwise_panels - 1

                        # Initialize a matrix to hold the vertices of the new row of
                        # wake ring vortices.
                        first_row_of_wake_ring_vortex_vertices = np.zeros(
                            (1, num_spanwise_panels + 1, 3)
                        )

                        # Iterate through the spanwise panel positions.
                        for spanwise_position in range(num_spanwise_panels):

                            # Get the next wing's panel object at this location.
                            next_panel = next_wing.panels[
                                chordwise_position, spanwise_position
                            ]

                            # The position of the next front left wake ring vortex
                            # vertex is the next panel's ring vortex's back left
                            # vertex.
                            next_front_left_vertex = next_panel.ring_vortex.Blpp_G_Cg

                            # Add this to the new row of wake ring vortex vertices.
                            first_row_of_wake_ring_vortex_vertices[
                                0, spanwise_position
                            ] = next_front_left_vertex

                            # Check if this panel is on the right edge of the wing.
                            if spanwise_position == (num_spanwise_panels - 1):
                                # The position of the next front right wake ring
                                # vortex vertex is the next panel's ring vortex's
                                # back right vertex.
                                next_front_right_vertex = (
                                    next_panel.ring_vortex.Brpp_G_Cg
                                )

                                # Add this to the new row of wake ring vortex vertices.
                                first_row_of_wake_ring_vortex_vertices[
                                    0, spanwise_position + 1
                                ] = next_front_right_vertex

                        # Set the next wing's matrix of wake ring vortex vertices to
                        # a copy of the row of new wake ring vortex vertices. This is
                        # correct because this is the first time step.
                        next_wing.wake_ring_vortex_vertices = np.copy(
                            first_row_of_wake_ring_vortex_vertices
                        )

                        # Initialize variables to hold the number of spanwise vertices.
                        num_spanwise_vertices = num_spanwise_panels + 1

                        # Initialize a new matrix to hold the second row of wake ring
                        # vortex vertices.
                        second_row_of_wake_ring_vortex_vertices = np.zeros(
                            (1, num_spanwise_panels + 1, 3)
                        )

                        # Iterate through the spanwise vertex positions.
                        for spanwise_vertex_position in range(num_spanwise_vertices):

                            # Get the corresponding vertex from the first row.
                            wake_ring_vortex_vertex = (
                                next_wing.wake_ring_vortex_vertices[
                                    0, spanwise_vertex_position
                                ]
                            )

                            if prescribed_wake:

                                # If the wake is prescribed, set the velocity at this
                                # vertex to the freestream velocity.
                                velocity_at_first_row_wake_ring_vortex_vertex = (
                                    self.current_freestream_velocity_G
                                )
                            else:

                                # If the wake is not prescribed, set the velocity at
                                # this vertex to the solution velocity at this point.
                                velocity_at_first_row_wake_ring_vortex_vertex = (
                                    self.calculate_solution_velocity(
                                        np.expand_dims(wake_ring_vortex_vertex, axis=0)
                                    )
                                )

                            # Update the second row with the interpolated position of
                            # the first vertex.
                            second_row_of_wake_ring_vortex_vertices[
                                0, spanwise_vertex_position
                            ] = (
                                wake_ring_vortex_vertex
                                + velocity_at_first_row_wake_ring_vortex_vertex
                                * self.delta_time
                            )

                        # Update the wing's wake ring vortex vertex matrix by
                        # vertically stacking the second row below it.
                        next_wing.wake_ring_vortex_vertices = np.vstack(
                            (
                                next_wing.wake_ring_vortex_vertices,
                                second_row_of_wake_ring_vortex_vertices,
                            )
                        )

                    # If this isn't the first step, then do this.
                    else:

                        # Set the next wing's wake ring vortex vertex matrix to a
                        # copy of this wing's wake ring vortex vertex matrix.
                        next_wing.wake_ring_vortex_vertices = np.copy(
                            this_wing.wake_ring_vortex_vertices
                        )

                        # Get the number of chordwise and spanwise vertices.
                        num_chordwise_vertices = (
                            next_wing.wake_ring_vortex_vertices.shape[0]
                        )
                        num_spanwise_vertices = (
                            next_wing.wake_ring_vortex_vertices.shape[1]
                        )

                        # Iterate through the chordwise and spanwise vertex positions.
                        for chordwise_vertex_position in range(num_chordwise_vertices):
                            for spanwise_vertex_position in range(
                                num_spanwise_vertices
                            ):

                                # Get the wake ring vortex vertex at this position.
                                wake_ring_vortex_vertex = (
                                    next_wing.wake_ring_vortex_vertices[
                                        chordwise_vertex_position,
                                        spanwise_vertex_position,
                                    ]
                                )

                                if prescribed_wake:

                                    # If the wake is prescribed, set the velocity at
                                    # this vertex to the freestream velocity.
                                    velocity_at_first_row_wake_vortex_vertex = (
                                        self.current_freestream_velocity_G
                                    )
                                else:

                                    # If the wake is not prescribed, set the velocity
                                    # at this vertex to the solution velocity at this
                                    # point.
                                    velocity_at_first_row_wake_vortex_vertex = (
                                        np.squeeze(
                                            self.calculate_solution_velocity(
                                                np.expand_dims(
                                                    wake_ring_vortex_vertex, axis=0
                                                )
                                            )
                                        )
                                    )

                                # Update the vertex at this point with its
                                # interpolated position.
                                next_wing.wake_ring_vortex_vertices[
                                    chordwise_vertex_position, spanwise_vertex_position
                                ] += (
                                    velocity_at_first_row_wake_vortex_vertex
                                    * self.delta_time
                                )

                        # Set the chordwise position to the trailing edge.
                        chordwise_position = this_wing.num_chordwise_panels - 1

                        # Initialize a new matrix to hold the new first row of wake
                        # ring vortex vertices.
                        first_row_of_wake_ring_vortex_vertices = np.zeros(
                            (1, this_wing.num_spanwise_panels + 1, 3)
                        )

                        # Iterate spanwise through the trailing edge panels.
                        for spanwise_position in range(this_wing.num_spanwise_panels):

                            # Get the panel object at this location on the next
                            # airplane's wing object.
                            next_panel = next_wing.panels[
                                chordwise_position, spanwise_position
                            ]

                            # Add the panel object's back left ring vortex vertex to
                            # the matrix of new wake ring vortex vertices.
                            first_row_of_wake_ring_vortex_vertices[
                                0, spanwise_position
                            ] = next_panel.ring_vortex.Blpp_G_Cg

                            if spanwise_position == (this_wing.num_spanwise_panels - 1):
                                # If the panel object is at the right edge of the
                                # wing, add its back right ring vortex vertex to the
                                # matrix of new wake ring vortex vertices.
                                first_row_of_wake_ring_vortex_vertices[
                                    0, spanwise_position + 1
                                ] = next_panel.ring_vortex.Brpp_G_Cg

                        # Stack the new first row of wake ring vortex vertices above
                        # the wing's matrix of wake ring vortex vertices.
                        next_wing.wake_ring_vortex_vertices = np.vstack(
                            (
                                first_row_of_wake_ring_vortex_vertices,
                                next_wing.wake_ring_vortex_vertices,
                            )
                        )

    # NOTE: I haven't yet started refactoring this method.
    def populate_next_airplanes_wake_vortices(self):
        """This method populates the locations of the next airplanes' wake vortices.

        This method is not vectorized but its loops only consume 0.4% of the runtime,
        so I have kept it as is for increased readability.

        :return: None
        """

        # Check if the current step is not the last step.
        if self.current_step < self.num_steps - 1:

            # Get the next step's airplane objects.
            next_airplanes = self.steady_problems[self.current_step + 1].airplanes

            for airplane_id, next_airplane in enumerate(next_airplanes):

                # Iterate through the copy of the current airplane's wing positions.
                for wing_id, this_wing in enumerate(
                    self.current_airplanes[airplane_id].wings
                ):
                    next_wing = next_airplane.wings[wing_id]

                    # Get the next wing's matrix of wake ring vortex vertices.
                    next_wing_wake_ring_vortex_vertices = (
                        next_wing.wake_ring_vortex_vertices
                    )

                    this_wing_wake_ring_vortices = (
                        self.current_airplanes[airplane_id]
                        .wings[wing_id]
                        .wake_ring_vortices
                    )

                    # Find the number of chordwise and spanwise vertices in the next
                    # wing's matrix of wake ring vortex vertices.
                    num_chordwise_vertices = next_wing_wake_ring_vortex_vertices.shape[
                        0
                    ]
                    num_spanwise_vertices = next_wing_wake_ring_vortex_vertices.shape[1]

                    # Initialize a new matrix to hold the new row of wake ring
                    # vortices.
                    new_row_of_wake_ring_vortices = np.empty(
                        (1, num_spanwise_vertices - 1), dtype=object
                    )

                    # Stack the new matrix on top of the copy of this wing's matrix
                    # and assign it to the next wing.
                    next_wing.wake_ring_vortices = np.vstack(
                        (new_row_of_wake_ring_vortices, this_wing_wake_ring_vortices)
                    )

                    # Iterate through the vertex positions.
                    for chordwise_vertex_position in range(num_chordwise_vertices):
                        for spanwise_vertex_position in range(num_spanwise_vertices):

                            # Set booleans to determine if this vertex is on the
                            # right and/or trailing edge of the wake.
                            has_right_vertex = (
                                spanwise_vertex_position + 1
                            ) < num_spanwise_vertices
                            has_back_vertex = (
                                chordwise_vertex_position + 1
                            ) < num_chordwise_vertices

                            if has_right_vertex and has_back_vertex:

                                # If this position is not on the right or trailing
                                # edge of the wake, get the four vertices that will
                                # be associated with the corresponding ring vortex at
                                # this position.
                                front_left_vertex = next_wing_wake_ring_vortex_vertices[
                                    chordwise_vertex_position, spanwise_vertex_position
                                ]
                                front_right_vertex = (
                                    next_wing_wake_ring_vortex_vertices[
                                        chordwise_vertex_position,
                                        spanwise_vertex_position + 1,
                                    ]
                                )
                                back_left_vertex = next_wing_wake_ring_vortex_vertices[
                                    chordwise_vertex_position + 1,
                                    spanwise_vertex_position,
                                ]
                                back_right_vertex = next_wing_wake_ring_vortex_vertices[
                                    chordwise_vertex_position + 1,
                                    spanwise_vertex_position + 1,
                                ]

                                if chordwise_vertex_position > 0:

                                    # If this isn't the front of the wake, update the
                                    # position of the ring vortex at this location.
                                    this_wake_ring_vortex: aerodynamics.RingVortex = (
                                        next_wing.wake_ring_vortices[
                                            chordwise_vertex_position,
                                            spanwise_vertex_position,
                                        ]
                                    )
                                    this_wake_ring_vortex.update_position(
                                        Flrvp_G_Cg=front_left_vertex,
                                        Frrvp_G_Cg=front_right_vertex,
                                        Blrvp_G_Cg=back_left_vertex,
                                        Brrvp_G_Cg=back_right_vertex,
                                    )

                                    # Also, update the age of this RingVortex.
                                    if self.current_step == 0:
                                        this_wake_ring_vortex.age = self.delta_time
                                    else:
                                        this_wake_ring_vortex.age += self.delta_time

                                if chordwise_vertex_position == 0:
                                    # If this is the front of the wake, get the
                                    # vortex strength from the wing panel's ring
                                    # vortex direction in front of it.
                                    this_strength_copy = this_wing.panels[
                                        this_wing.num_chordwise_panels - 1,
                                        spanwise_vertex_position,
                                    ].ring_vortex.strength

                                    # Then, make a new ring vortex at this location,
                                    # with the panel's ring vortex's strength,
                                    # and add it to the matrix of ring vortices.
                                    next_wing.wake_ring_vortices[
                                        chordwise_vertex_position,
                                        spanwise_vertex_position,
                                    ] = aerodynamics.RingVortex(
                                        Flrvp_G_Cg=front_left_vertex,
                                        Frrvp_G_Cg=front_right_vertex,
                                        Blrvp_G_Cg=back_left_vertex,
                                        Brrvp_G_Cg=back_right_vertex,
                                        strength=this_strength_copy,
                                    )

    # NOTE: I haven't yet started refactoring this method.
    def calculate_current_flapping_velocities_at_collocation_points(self):
        """This method finds the apparent flow velocity due to flapping at the
        centers of the current airplanes' collocation points.

        Note: The apparent flow velocity due to flapping is opposite the direction of
        the wing's motion.

        :return: size (M x 3) array of floats, where M is the current airplanes'
        number of panels
            This is an array containing the x, y, and z components of the velocity
            due to flapping at every one of the current airplanes' collocation
            points. Its units are in meters per second. If the current time step is
            the first time step, all the flapping velocities will be zero.
        """
        # Check if the current step is the first step.
        if self.current_step < 1:
            # Set the flapping velocities to be zero for all points. Then, return the
            # flapping velocities.
            return np.zeros((self.num_panels, 3))

        # Get the current airplane's collocation points, and the last airplane's
        # collocation points.
        these_collocations = self.panel_collocation_points
        last_collocations = self.last_panel_collocation_points

        # Calculate and return the flapping velocities.
        return -(these_collocations - last_collocations) / self.delta_time

    # NOTE: I haven't yet started refactoring this method.
    def calculate_current_flapping_velocities_at_right_leg_centers(self):
        """This method finds the apparent flow velocity due to flapping at the
        centers of the current airplanes' bound ring vortices' right legs.

        Note: The apparent flow velocity due to flapping is opposite the direction of
        the wing's motion.

        :return: size (M x 3) array of floats, where M is the current airplanes'
        number of panels
            This is an array containing the x, y, and z components of the velocity
            due to flapping at every one of the current airplanes' bound vortices'
            right legs' centers. Its units are in meters per second. If the current
            time step is the first time step, all the flapping velocities will be zero.
        """
        # Check if the current step is the first step.
        if self.current_step < 1:
            # Set the flapping velocities to be zero for all points. Then, return the
            # flapping velocities.
            return np.zeros((self.num_panels, 3))

        # Get the current airplanes' bound vortices' right legs' centers, and the
        # last airplanes' bound vortices' right legs' centers.
        these_right_leg_centers = self.panel_right_vortex_centers
        last_right_leg_centers = self.last_panel_right_vortex_centers

        # Calculate and return the flapping velocities.
        return -(these_right_leg_centers - last_right_leg_centers) / self.delta_time

    # NOTE: I haven't yet started refactoring this method.
    def calculate_current_flapping_velocities_at_front_leg_centers(self):
        """This method finds the apparent flow velocity due to flapping at the
        centers of the current airplanes' bound ring vortices' front legs.

        Note: The apparent flow velocity due to flapping is opposite the direction of
        the wing's motion.

        :return: size (M x 3) array of floats, where M is the current airplanes'
        number of panels
            This is an array containing the x, y, and z components of the velocity
            due to flapping at every one of the current airplanes' bound vortices'
            front legs' centers. Its units are in meters per second. If the current
            time step is the first time step, all the flapping velocities will be zero.
        """
        # Check if the current step is the first step.
        if self.current_step < 1:
            # Set the flapping velocities to be zero for all points. Then, return the
            # flapping velocities.
            return np.zeros((self.num_panels, 3))

        # Get the current airplanes' bound vortices' front legs' centers, and the
        # last airplanes' bound vortices' front legs' centers.
        these_front_leg_centers = self.panel_front_vortex_centers
        last_front_leg_centers = self.last_panel_front_vortex_centers

        # Calculate and return the flapping velocities.
        return -(these_front_leg_centers - last_front_leg_centers) / self.delta_time

    # NOTE: I haven't yet started refactoring this method.
    def calculate_current_flapping_velocities_at_left_leg_centers(self):
        """This method finds the apparent flow velocity due to flapping at the
        centers of the current airplanes' bound ring vortices' left legs.

        Note: The apparent flow velocity due to flapping is opposite the direction of
        the wing's motion.

        :return: size (M x 3) array of floats, where M is the current airplanes'
        number of panels
            This is an array containing the x, y, and z components of the velocity
            due to flapping at every one of the current airplanes' bound vortices'
            left legs' centers. Its units are in meters per second. If the current
            time step is the first time step, all the flapping velocities will be zero.
        """
        # Check if the current step is the first step.
        if self.current_step < 1:
            # Set the flapping velocities to be zero for all points. Then, return the
            # flapping velocities.
            return np.zeros((self.num_panels, 3))

        # Get the current airplanes' bound vortices' left legs' centers, and the last
        # airplanes' bound vortices' left legs' centers.
        these_left_leg_centers = self.panel_left_vortex_centers
        last_left_leg_centers = self.last_panel_left_vortex_centers

        # Calculate and return the flapping velocities.
        return -(these_left_leg_centers - last_left_leg_centers) / self.delta_time

    # NOTE: I haven't yet started refactoring this method.
    def finalize_forces_and_moments(self):
        """For cases with static geometry, this function finds the final loads and
        load coefficients for each of the problem's airplanes. For cases with
        variable geometry, it finds the final cycle-averaged and
        cycle-root-mean-squared loads and load coefficients.

        :return: None
        """
        # Get this solver's time step characteristics. Note that the first time step (
        # time step 0), occurs at 0 seconds.
        num_steps_to_average = self.num_steps - self.first_averaging_step

        # Determine if this problem's geometry is static or variable.
        is_static = self.unsteady_problem.max_period == 0

        # Initialize matrices to hold the forces, moments, and coefficients at each of
        # the time steps that has results.
        total_forces_W = np.zeros((self.num_airplanes, 3, num_steps_to_average))
        total_force_coefficients_W = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average)
        )
        total_moments_W_Cg = np.zeros((self.num_airplanes, 3, num_steps_to_average))
        total_moment_coefficients_W_Cg = np.zeros(
            (self.num_airplanes, 3, num_steps_to_average)
        )

        # Initialize a variable to track position in the results arrays.
        results_step = 0

        # Iterate through the time steps with results and add the results to their
        # respective matrices.
        for step in range(self.first_averaging_step, self.num_steps):

            # Get the airplanes from the problem at this step.
            these_airplanes = self.steady_problems[step].airplanes

            # Iterate through this step's airplanes.
            for airplane_id, airplane in enumerate(these_airplanes):
                total_forces_W[airplane_id, :, results_step] = airplane.total_force_W
                total_force_coefficients_W[airplane_id, :, results_step] = (
                    airplane.total_force_coefficients_W
                )
                total_moments_W_Cg[airplane_id, :, results_step] = (
                    airplane.total_moment_W
                )
                total_moment_coefficients_W_Cg[airplane_id, :, results_step] = (
                    airplane.total_moment_coefficients_W_Cg
                )

            results_step += 1

        # For each airplane object, calculate and then save the final or
        # cycle-averaged forces, moments, force coefficients, and moment coefficients.
        for airplane_id, airplane in enumerate(self.steady_problems[0].airplanes):

            if is_static:
                self.unsteady_problem.final_forces_W.append(
                    total_forces_W[airplane_id, :, -1]
                )
                self.unsteady_problem.final_force_coefficients_W.append(
                    total_force_coefficients_W[airplane_id, :, -1]
                )
                self.unsteady_problem.final_moments_W_Cg.append(
                    total_moments_W_Cg[airplane_id, :, -1]
                )
                self.unsteady_problem.final_moment_coefficients_W_Cg.append(
                    total_moment_coefficients_W_Cg[airplane_id, :, -1]
                )
            else:
                mean_forces = np.mean(total_forces_W[airplane_id], axis=-1)
                mean_force_coefficients = np.mean(
                    total_force_coefficients_W[airplane_id], axis=-1
                )
                mean_moments = np.mean(total_moments_W_Cg[airplane_id], axis=-1)
                mean_moment_coefficients = np.mean(
                    total_moment_coefficients_W_Cg[airplane_id], axis=-1
                )

                rms_forces = np.sqrt(
                    np.mean(
                        np.square(total_forces_W[airplane_id]),
                        axis=-1,
                    )
                )
                rms_force_coefficients = np.sqrt(
                    np.mean(
                        np.square(total_force_coefficients_W[airplane_id]),
                        axis=-1,
                    )
                )
                rms_moments = np.sqrt(
                    np.mean(
                        np.square(total_moments_W_Cg[airplane_id]),
                        axis=-1,
                    )
                )
                rms_moment_coefficients = np.sqrt(
                    np.mean(
                        np.square(total_moment_coefficients_W_Cg[airplane_id]),
                        axis=-1,
                    )
                )

                self.unsteady_problem.final_mean_forces_W.append(mean_forces)
                self.unsteady_problem.final_mean_force_coefficients_W.append(
                    mean_force_coefficients
                )
                self.unsteady_problem.final_mean_moments_W_Cg.append(mean_moments)
                self.unsteady_problem.final_mean_moment_coefficients_W_Cg.append(
                    mean_moment_coefficients
                )

                self.unsteady_problem.final_rms_forces_W.append(rms_forces)
                self.unsteady_problem.final_rms_force_coefficients_W.append(
                    rms_force_coefficients
                )
                self.unsteady_problem.final_rms_moments_W_Cg.append(rms_moments)
                self.unsteady_problem.final_rms_moment_coefficients_W_Cg.append(
                    rms_moment_coefficients
                )
