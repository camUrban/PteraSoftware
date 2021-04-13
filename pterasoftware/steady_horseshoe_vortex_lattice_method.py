""" This module contains the class definition of this package's steady horseshoe
vortex lattice solver.

This module contains the following classes:
    SteadyHorseshoeVortexLatticeMethodSolver: This is an aerodynamics solver that
    uses a steady horseshoe vortex lattice method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

import pterasoftware as ps


class SteadyHorseshoeVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses a steady horseshoe vortex lattice
    method.

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

        calculate_freestream_wing_influences: Find the normal velocity speed at every
        collocation points without the influence of the vortices.

        calculate_vortex_strengths: Solve for each panels' vortex strengths.

        calculate_near_field_forces_and_moments: Find the the forces and moments
        calculated from the near field.

        calculate_streamlines: Calculates the location of the streamlines coming off
        the back of the wings.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, steady_problem):
        """This is the initialization method.

        :param steady_problem: SteadyProblem
            This is the steady problem to be solved.
        :return: None
        """

        # Initialize this solution's attributes.
        self.airplane = steady_problem.airplane
        self.operating_point = steady_problem.operating_point

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.wing_wing_influences = np.zeros(
            (self.airplane.num_panels, self.airplane.num_panels)
        )
        self.freestream_velocity = (
            self.operating_point.calculate_freestream_velocity_geometry_axes()
        )
        self.freestream_wing_influences = np.zeros(self.airplane.num_panels)
        self.vortex_strengths = np.zeros(self.airplane.num_panels)
        self.panel_normal_directions = np.zeros((self.airplane.num_panels, 3))
        self.panel_areas = np.zeros(self.airplane.num_panels)
        self.panel_collocation_points = np.zeros((self.airplane.num_panels, 3))
        self.panel_vortex_strengths = np.zeros(self.airplane.num_panels)
        self.panel_back_right_vortex_vertices = np.zeros((self.airplane.num_panels, 3))
        self.panel_front_right_vortex_vertices = np.zeros((self.airplane.num_panels, 3))
        self.panel_front_left_vortex_vertices = np.zeros((self.airplane.num_panels, 3))
        self.panel_back_left_vortex_vertices = np.zeros((self.airplane.num_panels, 3))
        self.panels = np.empty(self.airplane.num_panels, dtype=object)
        self.panel_bound_vortex_centers = np.zeros((self.airplane.num_panels, 3))
        self.panel_bound_vortex_vectors = np.zeros((self.airplane.num_panels, 3))
        self.seed_points = np.empty((0, 3))
        self.streamline_points = None

    def run(self, verbose=True):
        """Run the solver on the steady problem.

        :param verbose: Bool, optional
            This parameter determines if the solver prints output to the console.
            It's default value is True.
        :return: None
        """

        # Initialize this problem's panels to have vortices congruent with this
        # solver type.
        if verbose:
            print("Initializing panel vortices.")
        self.initialize_panel_vortices()

        # Collapse this problem's geometry matrices into 1D ndarrays of attributes.
        if verbose:
            print("Collapsing geometry.")
        self.collapse_geometry()

        # Find the matrix of aerodynamic influence coefficients associated with this
        # problem's geometry.
        if verbose:
            print("\nCalculating the wing-wing influences.")
        self.calculate_wing_wing_influences()

        # Find the normal freestream speed at every collocation points without vortices.
        if verbose:
            print("\nCalculating the freestream-wing influences.")
        self.calculate_freestream_wing_influences()

        # Solve for each panel's vortex strengths.
        if verbose:
            print("\nCalculating vortex strengths.")
        self.calculate_vortex_strengths()

        # Solve for the near field forces and moments on each panel.
        if verbose:
            print("\nCalculating near field forces.")
        self.calculate_near_field_forces_and_moments()

        # Solve for the location of the streamlines coming off the back of the wings.
        if verbose:
            print("\nCalculating streamlines.")
        self.calculate_streamlines()

        # Print out the total forces.
        if verbose:
            print("\n\nForces in Wind Axes:")
            print(
                "\tInduced Drag:\t\t\t",
                np.round(self.airplane.total_near_field_force_wind_axes[0], 3),
                " N",
            )
            print(
                "\tSide Force:\t\t\t\t",
                np.round(self.airplane.total_near_field_force_wind_axes[1], 3),
                " N",
            )
            print(
                "\tLift:\t\t\t\t\t",
                np.round(self.airplane.total_near_field_force_wind_axes[2], 3),
                " N",
            )

        # Print out the total moments.
        if verbose:
            print("\nMoments in Wind Axes:")
            print(
                "\tRolling Moment:\t\t\t",
                np.round(self.airplane.total_near_field_moment_wind_axes[0], 3),
                " Nm",
            )
            print(
                "\tPitching Moment:\t\t",
                np.round(self.airplane.total_near_field_moment_wind_axes[1], 3),
                " Nm",
            )
            print(
                "\tYawing Moment:\t\t\t",
                np.round(self.airplane.total_near_field_moment_wind_axes[2], 3),
                " Nm",
            )

        # Print out the coefficients.
        if verbose:
            print("\nCoefficients in Wind Axes:")
            print(
                "\tCDi:\t\t\t\t\t",
                np.round(
                    self.airplane.total_near_field_force_coefficients_wind_axes[0], 3
                ),
            )
            print(
                "\tCY:\t\t\t\t\t\t",
                np.round(
                    self.airplane.total_near_field_force_coefficients_wind_axes[1], 3
                ),
            )
            print(
                "\tCL:\t\t\t\t\t\t",
                np.round(
                    self.airplane.total_near_field_force_coefficients_wind_axes[2], 3
                ),
            )
            print(
                "\tCl:\t\t\t\t\t\t",
                np.round(
                    self.airplane.total_near_field_moment_coefficients_wind_axes[0], 3
                ),
            )
            print(
                "\tCm:\t\t\t\t\t\t",
                np.round(
                    self.airplane.total_near_field_moment_coefficients_wind_axes[1], 3
                ),
            )
            print(
                "\tCn:\t\t\t\t\t\t",
                np.round(
                    self.airplane.total_near_field_moment_coefficients_wind_axes[2], 3
                ),
            )

    def initialize_panel_vortices(self):
        """This method calculates the locations of the vortex vertices, and then
        initializes the panels' vortices.

        Every panel has a horseshoe vortex. The vortex's finite leg runs along the
        panel's quarter chord from right to
        left. It's infinite legs points backwards in the positive x direction.

        :return: None
        """

        # Find the freestream direction in geometry axes.
        freestream_direction = (
            self.operating_point.calculate_freestream_direction_geometry_axes()
        )

        # Iterate through the current_airplane's wings.
        for wing in self.airplane.wings:

            # Find a suitable length for the "infinite" legs of the horseshoe
            # vortices on this wing. At twenty-times the
            # wing's span, these legs are essentially infinite.
            infinite_leg_length = wing.span * 20

            # Iterate through the wing's chordwise and spanwise panel positions.
            for chordwise_position in range(wing.num_chordwise_panels):
                for spanwise_position in range(wing.num_spanwise_panels):
                    # Pull the panel object out of the wing's list of panels.
                    panel = wing.panels[chordwise_position, spanwise_position]

                    # Find the location of the panel's front and right vortex vertices.
                    front_left_vortex_vertex = panel.front_left_vortex_vertex
                    front_right_vortex_vertex = panel.front_right_vortex_vertex

                    # Initialize the horseshoe vortex at this panel.
                    panel.horseshoe_vortex = ps.aerodynamics.HorseshoeVortex(
                        finite_leg_origin=front_right_vortex_vertex,
                        finite_leg_termination=front_left_vortex_vertex,
                        strength=None,
                        infinite_leg_direction=freestream_direction,
                        infinite_leg_length=infinite_leg_length,
                    )

    def collapse_geometry(self):
        """This method converts attributes of the problem's geometry into 1D
        ndarrays. This facilitates vectorization,
        which speeds up the solver.

        :return: None
        """

        # Initialize a variable to hold the global position of the panel as we
        # iterate through them.
        global_panel_position = 0

        # Iterate through the airplane's wings.
        for wing in self.airplane.wings:

            # Convert this wing's 2D array of panels into a 1D array.
            panels = np.ravel(wing.panels)

            # Iterate through the 1D array of this wing's panels.
            for panel in panels:

                # Update the solver's list of attributes with this panel's attributes.
                self.panels[global_panel_position] = panel
                self.panel_normal_directions[
                    global_panel_position, :
                ] = panel.normal_direction
                self.panel_areas[global_panel_position] = panel.area
                self.panel_collocation_points[
                    global_panel_position, :
                ] = panel.collocation_point
                self.panel_back_right_vortex_vertices[
                    global_panel_position, :
                ] = panel.horseshoe_vortex.right_leg.origin
                self.panel_front_right_vortex_vertices[
                    global_panel_position, :
                ] = panel.horseshoe_vortex.right_leg.termination
                self.panel_front_left_vortex_vertices[
                    global_panel_position, :
                ] = panel.horseshoe_vortex.left_leg.origin
                self.panel_back_left_vortex_vertices[
                    global_panel_position, :
                ] = panel.horseshoe_vortex.left_leg.termination
                self.panel_bound_vortex_centers[
                    global_panel_position, :
                ] = panel.horseshoe_vortex.finite_leg.center
                self.panel_bound_vortex_vectors[
                    global_panel_position, :
                ] = panel.horseshoe_vortex.finite_leg.vector

                # Check if this panel is on the trailing edge.
                if panel.is_trailing_edge:
                    # If it is, calculate it's streamline seed point and add it to
                    # the solver's array of seed points.
                    self.seed_points = np.vstack(
                        (
                            self.seed_points,
                            panel.back_left_vertex
                            + 0.5 * (panel.back_right_vertex - panel.back_left_vertex),
                        )
                    )

                # Increment the global panel position.
                global_panel_position += 1

    def calculate_wing_wing_influences(self):
        """This method finds the matrix of wing-wing influence coefficients
        associated with this airplane's geometry.

        :return: None
        """

        # Find the matrix of normalized velocities induced at every panel's
        # collocation point by every panel's horseshoe
        # vortex.
        induced_velocities = (
            ps.aerodynamics.calculate_velocity_induced_by_horseshoe_vortices(
                points=self.panel_collocation_points,
                back_right_vortex_vertices=self.panel_back_right_vortex_vertices,
                front_right_vortex_vertices=self.panel_front_right_vortex_vertices,
                front_left_vortex_vertices=self.panel_front_left_vortex_vertices,
                back_left_vortex_vertices=self.panel_back_left_vortex_vertices,
                strengths=np.ones(self.airplane.num_panels),
                collapse=False,
            )
        )

        # Take the batch dot product of the normalized velocities with each panel's
        # normal direction. This is now the
        # problem's matrix of wing-wing influence coefficients.
        self.wing_wing_influences = np.einsum(
            "...k,...k->...",
            induced_velocities,
            np.expand_dims(self.panel_normal_directions, axis=1),
        )

    def calculate_freestream_wing_influences(self):
        """This method finds the vector of freestream-wing influence coefficients
        associated with this problem.

        :return: None
        """

        # Take the batch dot product of the freestream velocity with each panel's
        # normal direction. This is now the
        # problem's 1D array of freestream-wing influence coefficients.
        self.freestream_wing_influences = np.einsum(
            "ij,j->i", self.panel_normal_directions, self.freestream_velocity
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

            # Update this panel's horseshoe vortex strength.
            panel.horseshoe_vortex.update_strength(self.vortex_strengths[panel_num])

    def calculate_near_field_forces_and_moments(self):
        """Find the the forces and moments calculated from the near field.

        Note: The forces and moments calculated are in geometry axes. The moment is
        about the airplane's reference
              point, which should be at the center of gravity. The units are Newtons
              and Newton-meters.

        :return: None
        """

        # Calculate the velocities induced at every panel's bound vortex center.
        induced_velocities = (
            ps.aerodynamics.calculate_velocity_induced_by_horseshoe_vortices(
                points=self.panel_bound_vortex_centers,
                back_right_vortex_vertices=self.panel_back_right_vortex_vertices,
                front_right_vortex_vertices=self.panel_front_right_vortex_vertices,
                front_left_vortex_vertices=self.panel_front_left_vortex_vertices,
                back_left_vortex_vertices=self.panel_back_left_vortex_vertices,
                strengths=self.vortex_strengths,
                collapse=True,
            )
        )

        # Add the freestream velocity to the induced velocities to calculate the
        # total velocity at every panel's bound
        # vortex center.
        total_velocities = induced_velocities + self.freestream_velocity

        # Calculate the near field force, in geometry axes, on each panel's bound
        # vortex.
        near_field_forces_geometry_axes = (
            self.operating_point.density
            * np.expand_dims(self.vortex_strengths, axis=1)
            * np.cross(total_velocities, self.panel_bound_vortex_vectors, axis=-1)
        )

        # Calculate the near field moments, in geometry axes, on each panel's bound
        # vortex.
        near_field_moments_geometry_axes = np.cross(
            self.panel_bound_vortex_centers - self.airplane.xyz_ref,
            near_field_forces_geometry_axes,
            axis=-1,
        )

        # Initialize a variable to hold the global panel position.
        global_panel_position = 0

        # Iterate through this solver's panels.
        for panel in self.panels:
            # Update the force and moment on this panel.
            panel.near_field_force_geometry_axes = near_field_forces_geometry_axes[
                global_panel_position, :
            ]
            panel.near_field_moment_geometry_axes = near_field_moments_geometry_axes[
                global_panel_position, :
            ]

            # Update the pressure on this panel.
            panel.update_pressure()

            # Increment the global panel position.
            global_panel_position += 1

        # Sum up the near field forces and moments on every panel to find the total
        # force and moment on the geometry.
        total_near_field_force_geometry_axes = np.sum(
            near_field_forces_geometry_axes, axis=0
        )
        total_near_field_moment_geometry_axes = np.sum(
            near_field_moments_geometry_axes, axis=0
        )

        # Find the total near field force in wind axes from the rotation matrix and
        # the total near field force in
        # geometry axes.
        self.airplane.total_near_field_force_wind_axes = (
            np.transpose(
                self.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes()
            )
            @ total_near_field_force_geometry_axes
        )

        # Find the total near field moment in wind axes from the rotation matrix and
        # the total near field moment in
        # geometry axes.
        self.airplane.total_near_field_moment_wind_axes = (
            np.transpose(
                self.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes()
            )
            @ total_near_field_moment_geometry_axes
        )

        # Calculate the current_airplane's induced drag coefficient
        induced_drag_coefficient = (
            -self.airplane.total_near_field_force_wind_axes[0]
            / self.operating_point.calculate_dynamic_pressure()
            / self.airplane.s_ref
        )

        # Calculate the current_airplane's side force coefficient.
        side_force_coefficient = (
            self.airplane.total_near_field_force_wind_axes[1]
            / self.operating_point.calculate_dynamic_pressure()
            / self.airplane.s_ref
        )

        # Calculate the current_airplane's lift coefficient.
        lift_coefficient = (
            -self.airplane.total_near_field_force_wind_axes[2]
            / self.operating_point.calculate_dynamic_pressure()
            / self.airplane.s_ref
        )

        # Calculate the current_airplane's rolling moment coefficient.
        rolling_moment_coefficient = (
            self.airplane.total_near_field_moment_wind_axes[0]
            / self.operating_point.calculate_dynamic_pressure()
            / self.airplane.s_ref
            / self.airplane.b_ref
        )

        # Calculate the current_airplane's pitching moment coefficient.
        pitching_moment_coefficient = (
            self.airplane.total_near_field_moment_wind_axes[1]
            / self.operating_point.calculate_dynamic_pressure()
            / self.airplane.s_ref
            / self.airplane.c_ref
        )

        # Calculate the current_airplane's yawing moment coefficient.
        yawing_moment_coefficient = (
            self.airplane.total_near_field_moment_wind_axes[2]
            / self.operating_point.calculate_dynamic_pressure()
            / self.airplane.s_ref
            / self.airplane.b_ref
        )

        self.airplane.total_near_field_force_coefficients_wind_axes = np.array(
            [induced_drag_coefficient, side_force_coefficient, lift_coefficient]
        )
        self.airplane.total_near_field_moment_coefficients_wind_axes = np.array(
            [
                rolling_moment_coefficient,
                pitching_moment_coefficient,
                yawing_moment_coefficient,
            ]
        )

    def calculate_streamlines(self, num_steps=10, delta_time=0.1):
        """Calculates the location of the streamlines coming off the back of the wings.

        :param num_steps: int, optional
            This is the integer number of points along each streamline (not including
            the initial points). It can be
            increased for higher fidelity visuals. The default value is 10.
        :param delta_time: float, optional
            This is the time in seconds between each time current_step It can be
            decreased for higher fidelity visuals
            or to make the streamlines shorter. It's default value is 0.1 seconds.
        :return: None
        """

        # Initialize a array to hold this problem's matrix of streamline points.
        self.streamline_points = np.expand_dims(self.seed_points, axis=0)

        # Iterate through the streamline steps.
        for step in range(num_steps):
            # Get the last row of streamline points.
            last_row_streamline_points = self.streamline_points[-1, :, :]

            # Find the induced velocities at this row of points.
            induced_velocities = (
                ps.aerodynamics.calculate_velocity_induced_by_horseshoe_vortices(
                    points=last_row_streamline_points,
                    back_right_vortex_vertices=self.panel_back_right_vortex_vertices,
                    front_right_vortex_vertices=self.panel_front_right_vortex_vertices,
                    front_left_vortex_vertices=self.panel_front_left_vortex_vertices,
                    back_left_vortex_vertices=self.panel_back_left_vortex_vertices,
                    strengths=self.vortex_strengths,
                    collapse=True,
                )
            )

            # Add the freestream velocity to the induced velocity to get the total
            # velocity at each of the last row of
            # streamline points.
            total_velocities = induced_velocities + self.freestream_velocity

            # Interpolate the positions on a new row of streamline points.
            new_row_streamline_points = (
                last_row_streamline_points + total_velocities * delta_time
            )

            # Stack the new row of streamline points to the bottom of the matrix of
            # streamline points.
            self.streamline_points = np.vstack(
                (
                    self.streamline_points,
                    np.expand_dims(new_row_streamline_points, axis=0),
                )
            )
