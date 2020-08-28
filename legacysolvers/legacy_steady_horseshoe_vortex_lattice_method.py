""" This module contains the class definition of a depreciated steady horseshoe vortex lattice solver.

This module contains the following classes:
    LegacySteadyHorseshoeVortexLatticeMethodSolver: This is an aerodynamics solver that uses a steady horseshoe vortex
                                                    lattice method. It has not been vectorized for speed.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

import pterasoftware as ps


class LegacySteadyHorseshoeVortexLatticeMethodSolver:
    """ This is an aerodynamics solver that uses a steady horseshoe vortex lattice method. It has not been vectorized for
    speed.

    Citation:
        Adapted from:         aerodynamics.vlm3.py in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    This class contains the following public methods:
        run: Run the solver on the steady problem.
        initialize_panel_vortices: This method calculates the locations of the vortex vertices, and then initializes the
                                   panels' vortices.
        calculate_wing_wing_influences: Find the matrix of aerodynamic influence coefficients associated with this
                                        problem's geometry.
        calculate_freestream_wing_influences: Find the normal velocity speed at every collocation point without the
                                              influence of the vortices.
        calculate_vortex_strengths: Solve for each panels' vortex strength.
        calculate_solution_velocity: Find the velocity at a given point due to the freestream and the vortices.
        calculate_near_field_forces_and_moments: Find the the forces and moments calculated from the near field.
        calculate_streamlines: Calculates the location of the streamlines coming off the back of the wings.

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
        self.streamline_points = None

    def run(self, verbose=True):
        """Run the solver on the steady problem.

        :param verbose: Bool, optional
            This parameter determines if the solver prints output to the console. It's default value is True.
        :return: None
        """

        # Initialize this problem's panels to have vortices congruent with this solver type.
        if verbose:
            print("Initializing panel vortices.")
        self.initialize_panel_vortices()

        # Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
        if verbose:
            print("\nCalculating the wing-wing influences.")
        self.calculate_wing_wing_influences()

        # Find the normal freestream speed at every collocation point without vortices.
        if verbose:
            print("\nCalculating the freestream-wing influences.")
        self.calculate_freestream_wing_influences()

        # Solve for each panel's vortex strength.
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
        """This method calculates the locations of the vortex vertices, and then initializes the panels' vortices.

        Every panel has a horseshoe vortex. The vortex's finite leg runs along the panel's quarter chord from right to
        left. It's infinite legs point backwards in the positive x direction.

        :return: None
        """

        # Find the freestream direction in geometry axes.
        freestream_direction = (
            self.operating_point.calculate_freestream_direction_geometry_axes()
        )

        # Iterate through the current_airplane's wings.
        for wing in self.airplane.wings:

            # Find a suitable length for the "infinite" legs of the horseshoe vortices on this wing. At twenty-times the
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

    def calculate_wing_wing_influences(self):
        """This method finds the matrix of wing-wing influence coefficients associated with this current_airplane's
        geometry.

        :return: None
        """

        # Initialize two variables to hold the global positions of the current collocation and vortex panel.
        global_collocation_panel_index = 0
        global_vortex_panel_index = 0

        # Iterate through the current_airplane's wings. This wing contains the panel with the collocation point where
        # the vortex influence is to be calculated.
        for collocation_panel_wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            collocation_panels = np.ravel(collocation_panel_wing.panels)

            # Iterate through the list of panels with the collocation points.
            for collocation_panel in collocation_panels:

                # Iterate through the current_airplane's wings. This wing contains the panel with the vortex whose
                # influence on the collocation point is to be calculated.
                for vortex_panel_wing in self.airplane.wings:

                    # Convert the 2D ndarray of this wing's panels into a 1D list.
                    vortex_panels = np.ravel(vortex_panel_wing.panels)

                    # Iterate through the list of panels with the vortices.
                    for vortex_panel in vortex_panels:

                        # Calculate the velocity induced at this collocation point by this vortex if the vortex's
                        # strength was 1.
                        normalized_induced_velocity_at_collocation_point = vortex_panel.calculate_normalized_induced_velocity(
                            collocation_panel.collocation_point
                        )

                        # Find the normal direction of the panel with the collocation point.
                        collocation_panel_normal_direction = (
                            collocation_panel.normal_direction
                        )

                        # Calculate the normal component of the velocity induced at this collocation point by this
                        # vortex if the vortex's strength was 1.
                        normal_normalized_induced_velocity_at_collocation_point = np.dot(
                            normalized_induced_velocity_at_collocation_point,
                            collocation_panel_normal_direction,
                        )

                        # Add this value to the solver's aerodynamic influence coefficient matrix.
                        self.wing_wing_influences[
                            global_collocation_panel_index, global_vortex_panel_index
                        ] = normal_normalized_induced_velocity_at_collocation_point

                        # At the next loop, we will analyze the effect of the next vortex panel on this collocation
                        # panel, so increment the global vortex panel's index.
                        global_vortex_panel_index += 1

                # At the next loop, we will analyze the effects of all the vortex panels on the next collocation panel,
                # so increment the collocation panel's global index. As we are starting over with the vortex panels, set
                # the vortex panel's global index to zero.
                global_collocation_panel_index += 1
                global_vortex_panel_index = 0

    def calculate_freestream_wing_influences(self):
        """This method finds the vector of freestream-wing influence coefficients associated with this problem.

        :return: None
        """

        # Initialize a variable to hold the global position of the current collocation panel.
        global_collocation_panel_index = 0

        # Iterate through the current_airplane's wings.
        for collocation_panel_wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            collocation_panels = np.ravel(collocation_panel_wing.panels)

            # Iterate through the list of panels with the collocation points.
            for collocation_panel in collocation_panels:

                # Update solver's list of freestream influences.
                self.freestream_wing_influences[
                    global_collocation_panel_index
                ] = np.dot(self.freestream_velocity, collocation_panel.normal_direction)

                # At the next loop, we will analyze the freestream influences on the next collocation panel, so
                # increment the panel's global index.
                global_collocation_panel_index += 1

    def calculate_vortex_strengths(self):
        """Solve for each panel's vortex strength.

        :return: None
        """

        # Solve for the strength of each panel's vortex.
        self.vortex_strengths = np.linalg.solve(
            self.wing_wing_influences, -self.freestream_wing_influences
        )

        # Initialize a variable to hold the global position of the current panel.
        global_panel_index = 0

        # Iterate through the current_airplane's wings.
        for wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            wing_panels = np.ravel(wing.panels)

            # Iterate through this list of panels.
            for panel in wing_panels:

                # Update each panel's vortex strength.
                panel.horseshoe_vortex.update_strength(
                    self.vortex_strengths[global_panel_index]
                )

                # At the next loop, we will update the vortex strength of the next panel, so increment the panel's
                # global index.
                global_panel_index += 1

    def calculate_solution_velocity(self, point):
        """Find the velocity at a given point due to the freestream and the vortices.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the velocity at.
        :return: None
        """

        # Initialize an ndarray to hold the solution velocity.
        velocity_induced_by_vortices = np.zeros(3)

        # Iterate through the current_airplane's wings.
        for wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            wing_panels = np.ravel(wing.panels)

            # Iterate through this list of panels.
            for panel in wing_panels:
                # Add the velocity induced by this panel's vortex at this point to the total induced velocity.
                velocity_induced_by_vortices += panel.calculate_induced_velocity(point)

        # Compute the freestream velocity in geometry axes.
        freestream = self.operating_point.calculate_freestream_velocity_geometry_axes()

        # Return the freestream velocity in geometry axes to the velocity induced by the vortices.
        return velocity_induced_by_vortices + freestream

    def calculate_near_field_forces_and_moments(self):
        """Find the the forces and moments calculated from the near field.

        :return: None
        """

        total_near_field_force_geometry_axes = np.zeros(3)
        total_near_field_moment_geometry_axes = np.zeros(3)

        # Iterate through the current_airplane's wings.
        for wing in self.airplane.wings:

            # Find the number of chordwise and spanwise panels in this wing.
            num_chordwise_panels = wing.num_chordwise_panels
            num_spanwise_panels = wing.num_spanwise_panels

            # Iterate through the chordwise and spanwise locations of the wing's panels.
            for chordwise_location in range(num_chordwise_panels):
                for spanwise_location in range(num_spanwise_panels):
                    # Find the panel at this location.
                    panel = wing.panels[chordwise_location, spanwise_location]

                    # Find the finite leg of this panel's horseshoe vortex.
                    bound_vortex = panel.horseshoe_vortex.finite_leg

                    # Calculate the velocity at the center of this finite leg.
                    velocity_at_vortex_center = self.calculate_solution_velocity(
                        bound_vortex.center
                    )

                    # Calculate the force and moment on this finite leg.
                    panel.near_field_force_geometry_axes = (
                        self.operating_point.density
                        * bound_vortex.strength
                        * np.cross(velocity_at_vortex_center, bound_vortex.vector)
                    )
                    panel.near_field_moment_geometry_axes = np.cross(
                        bound_vortex.center - self.airplane.xyz_ref,
                        panel.near_field_force_geometry_axes,
                    )

                    # Update the pressure on this panel.
                    panel.update_pressure()

                    # Add the near field force and moment on this panel to the total near field force and moment on this
                    # current_airplane.
                    total_near_field_force_geometry_axes += (
                        panel.near_field_force_geometry_axes
                    )
                    total_near_field_moment_geometry_axes += (
                        panel.near_field_moment_geometry_axes
                    )

        # Find the total near field force in wind axes from the rotation matrix and the total near field force in
        # geometry axes.
        self.airplane.total_near_field_force_wind_axes = (
            np.transpose(
                self.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes()
            )
            @ total_near_field_force_geometry_axes
        )

        # Find the total near field moment in wind axes from the rotation matrix and the total near field moment in
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
            This is the integer number of points along each streamline (not including the initial point). It can be
            increased for higher fidelity visuals. The default value is 10.
        :param delta_time: float, optional
            This is the time in seconds between each time current_step It can be decreased for higher fidelity visuals
            or to make the streamlines shorter. It's default value is 0.1 seconds.
        :return: None
        """

        # Initialize the streamline points attribute using the number of steps in this simulation.
        self.streamline_points = np.empty((num_steps + 1, 0, 3))

        # Iterate through the current_airplane's wings.
        for wing in self.airplane.wings:

            # Initialize an ndarray to hold the points along the streamline. It is shape (M x N x 3), where M is the
            # number of points in the streamline (not including the initial point), N is the number of spanwise panels
            # and thus the number of streamlines), and each bucket in this 2D array holds the x, y, and z components of
            # the streamline point's location.
            wing.streamline_points = np.zeros(
                (num_steps + 1, wing.num_spanwise_panels, 3)
            )

            # The chordwise position is along the trailing edge.
            chordwise_position = wing.num_chordwise_panels - 1

            # Increment through the wing's spanwise positions.
            for spanwise_position in range(wing.num_spanwise_panels):

                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]

                # The seed point is at the center of the panel's back leg.
                seed_point = panel.back_left_vertex + 0.5 * (
                    panel.back_right_vertex - panel.back_left_vertex
                )

                # Add the seed point to the array of streamline points.
                wing.streamline_points[0, spanwise_position, :] = seed_point

                # Iterate through the time steps.
                for step in range(num_steps):

                    # Find the last point.
                    last_point = wing.streamline_points[step, spanwise_position, :]

                    # Calculate the location of the new point, and add it to the array of streamline points.
                    wing.streamline_points[step + 1, spanwise_position, :] = (
                        last_point
                        + delta_time * self.calculate_solution_velocity(last_point)
                    )

            # Stack the current wing's streamline point matrix on to the solver's streamline point matrix.
            self.streamline_points = np.hstack(
                (self.streamline_points, wing.streamline_points)
            )
