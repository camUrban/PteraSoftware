"""This module contains the class definition of this package's steady ring vortex
lattice solver.

This module contains the following classes:
    SteadyRingVortexLatticeMethodSolver: This is an aerodynamics solver that uses a
    steady ring vortex lattice method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import logging

from . import aerodynamics
from . import functions


class SteadyRingVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses a steady ring vortex lattice method.

    Citation:
        Adapted from:         aerodynamics.vlm3.py in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    This class contains the following public methods:
        run: Run the solver on the steady problem.

        initialize_panel_vortices: This method calculates the locations of the vortex
        vertices, and then initializes the panels' vortices.

        collapse_geometry: This method converts attributes of the problem's geometry
        into 1D ndarrays. This facilitates vectorization, which speeds up the solver.

        calculate_wing_wing_influences: This method finds the matrix of wing-wing
        influence coefficients associated with this airplane's geometry.

        calculate_vortex_strengths: This method solves for each panel's vortex
        strength.

        calculate_solution_velocity: This function takes in a group of points. At
        every point, it finds the induced velocity due to every vortex and the
        freestream velocity.

        calculate_near_field_forces_and_moments: This method finds the forces and
        moments calculated from the near field.

        calculate_streamlines: This method calculates the location of the streamlines
        coming off the back of the wings.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, steady_problem):
        """This is the initialization method.

        :param steady_problem: SteadyProblem
            This is the steady problem to be solved.
        """
        # Initialize this solution's attributes.
        self.steady_problem = steady_problem
        self.airplanes = self.steady_problem.airplanes
        self.operating_point = self.steady_problem.operating_point
        self.num_airplanes = len(self.airplanes)

        # Calculate the total number of panels for all of this problem's airplanes.
        self.num_panels = 0
        for airplane in self.airplanes:
            self.num_panels += airplane.num_panels

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.wing_wing_influences = np.zeros((self.num_panels, self.num_panels))
        self.freestream_velocity = (
            self.operating_point.calculate_freestream_velocity_geometry_axes()
        )
        self.freestream_wing_influences = np.zeros(self.num_panels)
        self.vortex_strengths = np.ones(self.num_panels)
        self.panel_normal_directions = np.zeros((self.num_panels, 3))
        self.panel_areas = np.zeros(self.num_panels)
        self.panel_collocation_points = np.zeros((self.num_panels, 3))
        self.panel_back_right_vortex_vertices = np.zeros((self.num_panels, 3))
        self.panel_front_right_vortex_vertices = np.zeros((self.num_panels, 3))
        self.panel_front_left_vortex_vertices = np.zeros((self.num_panels, 3))
        self.panel_back_left_vortex_vertices = np.zeros((self.num_panels, 3))
        self.panels = np.empty(self.num_panels, dtype=object)
        self.panel_right_vortex_centers = np.zeros((self.num_panels, 3))
        self.panel_right_vortex_vectors = np.zeros((self.num_panels, 3))
        self.panel_front_vortex_centers = np.zeros((self.num_panels, 3))
        self.panel_front_vortex_vectors = np.zeros((self.num_panels, 3))
        self.panel_left_vortex_centers = np.zeros((self.num_panels, 3))
        self.panel_left_vortex_vectors = np.zeros((self.num_panels, 3))
        self.panel_back_vortex_centers = np.zeros((self.num_panels, 3))
        self.panel_back_vortex_vectors = np.zeros((self.num_panels, 3))
        self.panel_moment_references = np.zeros((self.num_panels, 3))
        self.seed_points = np.empty((0, 3))
        self.streamline_points = None

        # Initialize variables that will hold data which characterizes this panel's
        # horseshoe vortex. If the panel does not have a horseshoe vortex,
        # these values will not be updated. However, this does not adversely affect
        # the results, because the default horseshoe vortex strength is zero. The
        # default coordinates will also be updated by the collapse geometry method
        # for panels that have a horseshoe vortex.
        self.horseshoe_vortex_back_right_vertex = np.zeros((self.num_panels, 3))
        self.horseshoe_vortex_front_right_vertex = np.zeros((self.num_panels, 3))
        self.horseshoe_vortex_front_left_vertex = np.zeros((self.num_panels, 3))
        self.horseshoe_vortex_back_left_vertex = np.zeros((self.num_panels, 3))
        self.horseshoe_vortex_strengths = np.zeros(self.num_panels)

        # Initialize variables to hold details about this panel's location on its wing.
        self.panel_is_trailing_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_leading_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_left_edge = np.zeros(self.num_panels, dtype=bool)
        self.panel_is_right_edge = np.zeros(self.num_panels, dtype=bool)

    def run(self, logging_level="Warning"):
        """Run the solver on the steady problem.

        :param logging_level: str, optional
            This parameter determines the detail of information that the solver's
            logger will output while running. The options are, in order of detail and
            severity, "Debug", "Info", "Warning", "Error", "Critical". The default
            value is "Warning".
        :return: None
        """
        # Configure the problem's logger.
        logging_level_value = functions.convert_logging_level_name_to_value(
            logging_level
        )
        logging.basicConfig(level=logging_level_value)

        # Initialize this problem's panels to have vortices congruent with this
        # solver type.
        logging.info("Initializing panel vortices.")
        self.initialize_panel_vortices()

        # Collapse this problem's geometry matrices into 1D ndarrays of attributes.
        logging.info("Collapsing geometry.")
        self.collapse_geometry()

        # Find the matrix of wing-wing influence coefficients associated with this
        # current_airplane's geometry.
        logging.info("Calculating the wing-wing influences.")
        self.calculate_wing_wing_influences()

        # Find the vector of freestream-wing influence coefficients associated with
        # this problem.
        logging.info("Calculating the freestream-wing influences.")
        functions.calculate_steady_freestream_wing_influences(steady_solver=self)

        # Solve for each panel's vortex strength.
        logging.info("Calculating vortex strengths.")
        self.calculate_vortex_strengths()

        # Solve for the near field forces and moments on each panel.
        logging.info("Calculating near field forces.")
        self.calculate_near_field_forces_and_moments()

        # Solve for the location of the streamlines coming off the back of the wings.
        logging.info("Calculating streamlines.")
        functions.calculate_streamlines(self)

    def initialize_panel_vortices(self):
        """This method calculates the locations of the vortex vertices, and then
        initializes the panels' vortices.

        Every panel has a ring vortex, which is a quadrangle whose front vortex leg
        is at the panel's quarter chord. The left and right vortex legs run along the
        panel's left and right legs. If the panel is not along the trailing edge,
        they extend backwards and meet the back vortex leg at the rear panel's
        quarter chord. Otherwise, they extend back backwards and meet the back vortex
        leg one quarter chord back from the current panel's back leg.

        Panels that are at the trailing edge of a wing have a horseshoe vortex in
        addition to their ring vortex. The horseshoe vortex's finite leg runs along
        the ring vortex's back leg but in the opposite direction. Its infinite legs
        point backwards in the direction of the freestream. The ring vortex and
        horseshoe vortex have the same strength, so the back leg of the effects of
        the ring vortex's back leg and the horseshoe vortex's finite leg cancel each
        other.

        :return: None
        """
        # Find the freestream direction in geometry axes.
        freestream_direction = (
            self.operating_point.calculate_freestream_direction_geometry_axes()
        )

        # Iterate through each airplane's wings.
        for airplane in self.airplanes:
            for wing in airplane.wings:

                # Find a suitable length for the "infinite" legs of the horseshoe
                # vortices on this wing. At twenty-times the wing's span, these legs
                # are essentially infinite.
                infinite_leg_length = wing.span * 20

                # Iterate through the wing's chordwise and spanwise positions.
                for chordwise_position in range(wing.num_chordwise_panels):
                    for spanwise_position in range(wing.num_spanwise_panels):

                        # Pull the panel object out of the wing's list of panels.
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
                            back_left_vortex_vertex = next_chordwise_panel.Flbvp_G_Cg
                            back_right_vortex_vertex = next_chordwise_panel.Frbvp_G_Cg
                        else:
                            back_left_vortex_vertex = front_left_vortex_vertex + (
                                panel.Blpp_G_Cg - panel.Flpp_G_Cg
                            )
                            back_right_vortex_vertex = front_right_vortex_vertex + (
                                panel.Brpp_G_Cg - panel.Frpp_G_Cg
                            )

                            # If the panel is along the trailing edge, initialize its
                            # horseshoe vortex.
                            panel.horseshoe_vortex = aerodynamics.HorseshoeVortex(
                                Frhvp_G_Cg=back_right_vortex_vertex,
                                Flhvp_G_Cg=back_left_vortex_vertex,
                                leftLegVector_G=freestream_direction,
                                left_right_leg_lengths=infinite_leg_length,
                                strength=None,
                            )

                        # Initialize the panel's ring vortex.
                        panel.ring_vortex = aerodynamics.RingVortex(
                            Flrvp_G_Cg=front_left_vortex_vertex,
                            Frrvp_G_Cg=front_right_vortex_vertex,
                            Blrvp_G_Cg=back_left_vortex_vertex,
                            Brrvp_G_Cg=back_right_vortex_vertex,
                            strength=None,
                        )

    def collapse_geometry(self):
        """This method converts attributes of the problem's geometry into 1D
        ndarrays. This facilitates vectorization, which speeds up the solver.

        :return: None
        """
        # Initialize a variable to hold the global position of the panel as we
        # iterate through them.
        global_panel_position = 0

        # Iterate through each airplane's wings.
        for airplane in self.airplanes:
            for wing in airplane.wings:

                # Convert this wing's 2D array of panels into a 1D array.
                panels = np.ravel(wing.panels)

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

                    if panel.is_trailing_edge:
                        # Also, update the attribute lists horseshoe vortex
                        # attributes at this position with this panel's horseshoe
                        # vortex attributes
                        self.horseshoe_vortex_back_right_vertex[
                            global_panel_position
                        ] = panel.horseshoe_vortex.rightLeg_G.Slvp_G_Cg
                        self.horseshoe_vortex_front_right_vertex[
                            global_panel_position
                        ] = panel.horseshoe_vortex.rightLeg_G.Elvp_G_Cg
                        self.horseshoe_vortex_front_left_vertex[
                            global_panel_position
                        ] = panel.horseshoe_vortex.leftLeg_G.Slvp_G_Cg
                        self.horseshoe_vortex_back_left_vertex[
                            global_panel_position
                        ] = panel.horseshoe_vortex.leftLeg_G.Elvp_G_Cg

                        # Set the horseshoe vortex strength at this position to 1.0.
                        # This will be updated after the correct vortex strengths are
                        # calculated.
                        self.horseshoe_vortex_strengths[global_panel_position] = 1.0

                    # Increment the global panel position.
                    global_panel_position += 1

    def calculate_wing_wing_influences(self):
        """This method finds the matrix of wing-wing influence coefficients
        associated with the airplanes' geometry.

        :return: None
        """
        # Find the matrix of normalized velocities induced at every panel's
        # collocation point by every panel's ring vortex. The answer is normalized
        # because the solver's vortex strength list was initialized to all ones. This
        # will be updated once the correct vortex strengths' are calculated.
        ring_vortex_influences = aerodynamics.expanded_velocities_from_ring_vortices(
            stackP_G_Cg=self.panel_collocation_points,
            stackBrrvp_G_Cg=self.panel_back_right_vortex_vertices,
            stackFrrvp_G_Cg=self.panel_front_right_vortex_vertices,
            stackFlrvp_G_Cg=self.panel_front_left_vortex_vertices,
            stackBlrvp_G_Cg=self.panel_back_left_vortex_vertices,
            strengths=self.vortex_strengths,
        )

        # Find the matrix of normalized velocities induced at every panel's
        # collocation point by every panel's horseshoe vortex. The answer is
        # normalized because the solver's horseshoe vortex strength list was
        # initialized to ones for locations which have horseshoe vortices, and zeros
        # everywhere else. The strengths at the positions with horseshoe vortices
        # will be updated once the correct vortex strengths' are calculated. The
        # positions elsewhere will remain zero.
        horseshoe_vortex_influences = (
            aerodynamics.expanded_velocities_from_horseshoe_vortices(
                stackP_G_Cg=self.panel_collocation_points,
                stackBrhvp_G_Cg=self.horseshoe_vortex_back_right_vertex,
                stackFrhvp_G_Cg=self.horseshoe_vortex_front_right_vertex,
                stackFlhvp_G_Cg=self.horseshoe_vortex_front_left_vertex,
                stackBlhvp_G_Cg=self.horseshoe_vortex_back_left_vertex,
                strengths=self.horseshoe_vortex_strengths,
            )
        )

        # Calculate the total normalized influences, which is the sum of the
        # influences by the ring vortices and by the horseshoe vortices.
        total_influences = ring_vortex_influences + horseshoe_vortex_influences

        # Take the batch dot product of the normalized velocities with each panel's
        # normal direction. This is now the problem's matrix of wing-wing influence
        # coefficients.
        self.wing_wing_influences = np.einsum(
            "...k,...k->...",
            total_influences,
            np.expand_dims(self.panel_normal_directions, axis=1),
        )

    def calculate_vortex_strengths(self):
        """Solve for each panel's vortex strengths.

        :return: None
        """
        # Solve for the strength of each panel's vortex.
        self.vortex_strengths = np.linalg.solve(
            self.wing_wing_influences, -self.freestream_wing_influences
        )

        # Iterate through the panels and update their vortex strengths.
        for panel_num in range(self.panels.size):

            # Get the panel at this location.
            panel = self.panels[panel_num]

            # Update this panel's ring vortex strength.
            panel.ring_vortex.update_strength(self.vortex_strengths[panel_num])

            # Check if the panel has a horseshoe vortex.
            if panel.horseshoe_vortex is not None:
                # Update the panel's horseshoe vortex strength.
                panel.horseshoe_vortex.update_strength(self.vortex_strengths[panel_num])
                self.horseshoe_vortex_strengths[panel_num] = self.vortex_strengths[
                    panel_num
                ]

    def calculate_solution_velocity(self, points):
        """This function takes in a group of points. At every point, it finds the
        induced velocity due to every vortex and the freestream velocity.

        Note: The velocity calculated by this method is in geometry axes. Also,
        this method assumes that the correct vortex strengths have already been
        calculated.

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
        # Find the matrix of velocities induced at every point by every panel's ring
        # vortex. The effect of every ring vortex on each point will be summed.
        ring_vortex_influences = aerodynamics.collapsed_velocities_from_ring_vortices(
            stackP_G_Cg=points,
            stackBrrvp_G_Cg=self.panel_back_right_vortex_vertices,
            stackFrrvp_G_Cg=self.panel_front_right_vortex_vertices,
            stackFlrvp_G_Cg=self.panel_front_left_vortex_vertices,
            stackBlrvp_G_Cg=self.panel_back_left_vortex_vertices,
            strengths=self.vortex_strengths,
        )

        # Find the matrix of velocities induced at every panel's collocation point by
        # every panel's horseshoe vortex. This can be called in a batch fashion,
        # because every panel has horseshoe vortex attributes that are defined.
        # However, panels not on the trailing edge have horseshoe vortex strengths to
        # zero, which eliminates their effects. The effect of every horseshoe vortex
        # on each point will be summed.
        horseshoe_vortex_influences = (
            aerodynamics.collapsed_velocities_from_horseshoe_vortices(
                stackP_G_Cg=points,
                stackBrhvp_G_Cg=self.horseshoe_vortex_back_right_vertex,
                stackFrhvp_G_Cg=self.horseshoe_vortex_front_right_vertex,
                stackFlhvp_G_Cg=self.horseshoe_vortex_front_left_vertex,
                stackBlhvp_G_Cg=self.horseshoe_vortex_back_left_vertex,
                strengths=self.horseshoe_vortex_strengths,
            )
        )

        # Find the total influence of the vortices, which is the sum of the influence
        # due to the ring vortices and the horseshoe vortices.
        total_influences = ring_vortex_influences + horseshoe_vortex_influences

        # Calculate and return the freestream velocity added to the velocity induced
        # by the vortices. This is in geometry axes.
        return total_influences + self.freestream_velocity

    def calculate_near_field_forces_and_moments(self):
        """This method finds the forces and moments calculated from the near field.

        Citation: This method uses logic described on pages 9-11 of "Modeling of
        aerodynamic forces in flapping flight with the Unsteady Vortex Lattice
        Method" by Thomas Lambert.

        Note: The forces and moments calculated are in geometry axes. The moment is
        about the airplane's reference point, which should be at the center of
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

        # Iterate through the airplanes' wings.
        for airplane in self.airplanes:
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
                            self.vortex_strengths[global_panel_position]
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
                            self.vortex_strengths[global_panel_position]
                            - panel_to_right.ring_vortex.strength
                        )

                    # Check if this panel is on its wing's leading edge.
                    if panel.is_leading_edge:

                        # Change the effective front vortex line strength from zero
                        # to this panel's ring vortex's strength.
                        effective_front_vortex_line_strengths[global_panel_position] = (
                            self.vortex_strengths[global_panel_position]
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
                            self.vortex_strengths[global_panel_position]
                            - panel_to_front.ring_vortex.strength
                        )

                    # Check if this panel is on its wing's left edge.
                    if panel.is_left_edge:

                        # Change the effective left vortex line strength from zero to
                        # this panel's ring vortex's strength.
                        effective_left_vortex_line_strengths[global_panel_position] = (
                            self.vortex_strengths[global_panel_position]
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
                            self.vortex_strengths[global_panel_position]
                            - panel_to_left.ring_vortex.strength
                        )

                    # Increment the global panel position.
                    global_panel_position += 1

        # Calculate the solution velocities at the centers of the panel's front leg,
        # left leg, and right leg.
        velocities_at_ring_vortex_front_leg_centers = self.calculate_solution_velocity(
            points=self.panel_front_vortex_centers
        )
        velocities_at_ring_vortex_left_leg_centers = self.calculate_solution_velocity(
            points=self.panel_left_vortex_centers
        )
        velocities_at_ring_vortex_right_leg_centers = self.calculate_solution_velocity(
            points=self.panel_right_vortex_centers
        )

        # Using the effective line vortex strengths, and the Kutta-Joukowski theorem
        # to find the near field force in geometry axes on the front leg, left leg,
        # and right leg.
        near_field_forces_on_ring_vortex_right_legs_geometry_axes = (
            self.operating_point.density
            * np.expand_dims(effective_right_vortex_line_strengths, axis=1)
            * np.cross(
                velocities_at_ring_vortex_right_leg_centers,
                self.panel_right_vortex_vectors,
                axis=-1,
            )
        )
        near_field_forces_on_ring_vortex_front_legs_geometry_axes = (
            self.operating_point.density
            * np.expand_dims(effective_front_vortex_line_strengths, axis=1)
            * np.cross(
                velocities_at_ring_vortex_front_leg_centers,
                self.panel_front_vortex_vectors,
                axis=-1,
            )
        )
        near_field_forces_on_ring_vortex_left_legs_geometry_axes = (
            self.operating_point.density
            * np.expand_dims(effective_left_vortex_line_strengths, axis=1)
            * np.cross(
                velocities_at_ring_vortex_left_leg_centers,
                self.panel_left_vortex_vectors,
                axis=-1,
            )
        )

        # Sum the forces on the legs to calculate the total near field force,
        # in geometry axes, on each panel.
        near_field_forces_geometry_axes = (
            near_field_forces_on_ring_vortex_front_legs_geometry_axes
            + near_field_forces_on_ring_vortex_left_legs_geometry_axes
            + near_field_forces_on_ring_vortex_right_legs_geometry_axes
        )

        # Find the near field moment in geometry axes on the front leg, left leg,
        # and right leg.
        near_field_moments_on_ring_vortex_front_legs_geometry_axes = np.cross(
            self.panel_front_vortex_centers - self.panel_moment_references,
            near_field_forces_on_ring_vortex_front_legs_geometry_axes,
            axis=-1,
        )
        near_field_moments_on_ring_vortex_left_legs_geometry_axes = np.cross(
            self.panel_left_vortex_centers - self.panel_moment_references,
            near_field_forces_on_ring_vortex_left_legs_geometry_axes,
            axis=-1,
        )
        near_field_moments_on_ring_vortex_right_legs_geometry_axes = np.cross(
            self.panel_right_vortex_centers - self.panel_moment_references,
            near_field_forces_on_ring_vortex_right_legs_geometry_axes,
            axis=-1,
        )

        # Sum the moments on the legs to calculate the total near field moment,
        # in geometry axes, on each panel.
        near_field_moments_geometry_axes = (
            near_field_moments_on_ring_vortex_front_legs_geometry_axes
            + near_field_moments_on_ring_vortex_left_legs_geometry_axes
            + near_field_moments_on_ring_vortex_right_legs_geometry_axes
        )

        functions.process_steady_solver_forces(
            steady_solver=self,
            near_field_forces_geometry_axes=near_field_forces_geometry_axes,
            near_field_moments_geometry_axes=near_field_moments_geometry_axes,
        )
