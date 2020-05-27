
"""This module contains the class definition of this package's unsteady ring vortex lattice solver.

This module contains the following classes:
    UnsteadyRingVortexLatticeMethodSolver: This is an aerodynamics solver that uses an unsteady ring vortex lattice
                                           method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


class UnsteadyRingVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses an unsteady ring vortex lattice method.

    This class contains the following public methods:
        run: Run the solver on the unsteady problem.
        initialize_panel_vortices: This method calculates the locations of the vortex vertices, and then initializes the
                                   panels' vortices.
        set_up_geometry: Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
        set_up_operating_point: Find the normal freestream speed at every collocation point without the influence of the
                                vortices.
        calculate_vortex_strengths: Solve for each panel's vortex strength.
        calculate_solution_velocity: Find the velocity at a given point due to the freestream and the vortices.
        calculate_velocity_influences: Find the velocity at a given point due to the vorticity of every vortex if their
                                       strengths were all set to 1.0 meters squared per second.
        calculate_near_field_forces_and_moments: Find the the forces and moments calculated from the near field.
        calculate_streamlines: Calculates the location of the streamlines coming off the back of the wings.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, unsteady_problem):
        """This is the initialization method.

        :param unsteady_problem: UnsteadyProblem
            This is the unsteady problem to be solved.
        """

        # Initialize this solution's attributes.
        self.unsteady_problem = unsteady_problem
        self.airplane = self.unsteady_problem.airplane
        self.operating_point = self.unsteady_problem.operating_point
        self.movement = self.unsteady_problem.movement

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.aerodynamic_influence_coefficients = np.zeros((self.airplane.num_panels, self.airplane.num_panels))
        self.freestream_velocity = np.zeros(3)
        self.normal_directions = np.zeros((self.airplane.num_panels, 3))
        self.freestream_influences = np.zeros(self.airplane.num_panels)
        self.vortex_strengths = np.zeros(self.airplane.num_panels)
        self.total_near_field_force_wind_axes = np.zeros(3)
        self.total_near_field_moment_wind_axes = np.zeros(3)
        self.CL = None
        self.CDi = None
        self.CY = None
        self.Cl = None
        self.Cm = None
        self.Cn = None

    def run(self):
        """Run the solver on the unsteady problem.

        :return: None
        """

        # Initialize this problem's panels to have vortices congruent with this solver type.
        print("Initializing panel vortices...")
        self.initialize_panel_vortices()
        print("Panel vortices initialized.")

        # Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
        print("\nSetting up geometry...")
        self.set_up_geometry()
        print("Geometry set up.")

        print("\nShedding wake vortices...")

        for i in range(1600):

            print("\nBeginning time step " + str(i+1) + " out of " + str(10) + "...")

            # Find the normal freestream speed at every collocation point without vortices.
            print("Setting up operating point...")
            self.set_up_operating_point()
            print("Operating point set up.")

            # Solve for each panel's vortex strength.
            print("Calculating vortex strengths...")
            self.calculate_vortex_strengths()
            print("Vortex strengths calculated.")

            # Solve for the near field forces and moments on each panel.
            print("Calculating wake rollup...")
            self.calculate_wake_rollup()
            print("Wake rollup calculated.")

            print("Finished time step " + str(i+1) + " out of " + str(10) + ".")

        print("\nWake vortices shed.")

        # Solve for the near field forces and moments on each panel.
        print("\nCalculating near field forces...")
        self.calculate_near_field_forces_and_moments()
        print("Near field forces calculated.")

        # Solve for the location of the streamlines coming off the back of the wings.
        print("\nCalculating streamlines...")
        self.calculate_streamlines()
        print("Streamlines calculated.")

        # Print out the total forces.
        print("\n\nTotal Forces in Wind Axes:")
        print("\tInduced Drag:\t\t\t", np.round(self.total_near_field_force_wind_axes[0], 3), " N")
        print("\tSide Force:\t\t\t\t", np.round(self.total_near_field_force_wind_axes[1], 3), " N")
        print("\tLift:\t\t\t\t\t", np.round(self.total_near_field_force_wind_axes[2], 3), " N")

        # Print out the total moments.
        print("\nTotal Moments in Wind Axes:")
        print("\tRolling Moment:\t\t\t", np.round(self.total_near_field_moment_wind_axes[0], 3), " Nm")
        print("\tPitching Moment:\t\t", np.round(self.total_near_field_moment_wind_axes[1], 3), " Nm")
        print("\tYawing Moment:\t\t\t", np.round(self.total_near_field_moment_wind_axes[2], 3), " Nm")

        # Print out the coefficients.
        print("\nCoefficients in Wind Axes:")
        print("\tCDi:\t\t\t\t\t", np.round(self.CDi, 3))
        print("\tCY:\t\t\t\t\t\t", np.round(self.CY, 3))
        print("\tCL:\t\t\t\t\t\t", np.round(self.CL, 3))
        print("\tCl:\t\t\t\t\t\t", np.round(self.Cl, 3))
        print("\tCm:\t\t\t\t\t\t", np.round(self.Cm, 3))
        print("\tCn:\t\t\t\t\t\t", np.round(self.Cn, 3))

    def initialize_panel_vortices(self):
        """This method calculates the locations of the vortex vertices, and then initializes the panel's vortices.

        Every panel has a ring vortex, which is a quadrangle whose front vortex leg is at the panel's quarter chord.
        The left and right vortex legs run along the panel's left and right legs. If the panel is not along the
        trailing edge, they extend backwards and meet the back vortex leg at a length of one quarter of the rear
        panel's chord back from the rear panel's front leg. Otherwise, they extend back backwards and meet the back
        vortex leg at a length of one quarter of the current panel's chord back from the current panel's back leg.

        :return: None
        """

        # Iterate through the airplane's wings.
        for wing in self.airplane.wings:

            # Iterate through the wing's chordwise and spanwise positions.
            for chordwise_position in range(wing.num_chordwise_panels):
                for spanwise_position in range(wing.num_spanwise_panels):

                    # Pull the panel object out of the wing's list of panels.
                    panel = wing.panels[chordwise_position, spanwise_position]

                    # Find the location of the panel's front left and right vortex vertices.
                    front_left_vortex_vertex = panel.front_left_vortex_vertex
                    front_right_vortex_vertex = panel.front_right_vortex_vertex

                    # Define the back left and right vortex vertices based on whether the panel is along the trailing
                    # edge or not.
                    if not panel.is_trailing_edge:
                        next_chordwise_panel = wing.panels[chordwise_position + 1, spanwise_position]
                        back_left_vortex_vertex = next_chordwise_panel.front_left_vortex_vertex
                        back_right_vortex_vertex = next_chordwise_panel.front_right_vortex_vertex
                    else:
                        back_left_vortex_vertex = front_left_vortex_vertex + (
                                panel.back_left_vertex - panel.front_left_vertex)
                        back_right_vortex_vertex = front_right_vortex_vertex + (
                                panel.back_right_vertex - panel.front_right_vertex)

                    # Initialize the panel's ring vortex.
                    panel.ring_vortex = asmvp.aerodynamics.RingVortex(
                        front_right_vertex=front_right_vortex_vertex,
                        front_left_vertex=front_left_vortex_vertex,
                        back_left_vertex=back_left_vortex_vertex,
                        back_right_vertex=back_right_vortex_vertex,
                        strength=None
                    )

    def set_up_geometry(self):
        """Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.

        :return: None
        """

        # Iterate through the airplane's wings. This wing contains the panel with the collocation point where the
        # vortex influence is to be calculated.
        for collocation_panel_wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            collocation_panels = np.ravel(collocation_panel_wing.panels)

            # Iterate through the list of panels with the collocation points.
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panels):

                # Iterate through the airplane's wings. This wing contains the panel with the vortex whose influence
                # on the collocation point is to be calculated.
                for vortex_panel_wing in self.airplane.wings:

                    # Convert the 2D ndarray of this wing's panels into a 1D list.
                    vortex_panel_wings_panels = np.ravel(vortex_panel_wing.panels)

                    # Iterate through the list of panels with the vortices.
                    for vortex_panel_index, vortex_panel in np.ndenumerate(vortex_panel_wings_panels):
                        # Calculate the velocity induced at this collocation point by this vortex if the vortex's
                        # strength was 1.
                        normalized_induced_velocity_at_collocation_point = (
                            vortex_panel.calculate_normalized_induced_velocity(collocation_panel.collocation_point))

                        # Find the normal direction of the panel with the collocation point.
                        collocation_panel_normal_direction = collocation_panel.normal_direction

                        # Calculate the normal component of the velocity induced at this collocation point by this
                        # vortex if the vortex's strength was 1.
                        normal_normalized_induced_velocity_at_collocation_point = np.dot(
                            normalized_induced_velocity_at_collocation_point, collocation_panel_normal_direction)

                        # Add this value to the solver's aerodynamic influence coefficient matrix.
                        self.aerodynamic_influence_coefficients[collocation_panel_index, vortex_panel_index] = (
                            normal_normalized_induced_velocity_at_collocation_point)

    def set_up_operating_point(self):
        """Find the normal freestream speed at every collocation point without the influence of the vortices.

        :return: None
        """

        # This calculates and updates the direction the wind is going to, in geometry axes coordinates.
        self.freestream_velocity = np.expand_dims(self.operating_point.calculate_freestream_velocity_geometry_axes(), 0)

        # Iterate through the airplane's wings.
        for collocation_panel_wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            collocation_panels = np.ravel(collocation_panel_wing.panels)

            # Iterate through the list of panels with the collocation points.
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panels):
                # Update the solver's list of normal directions.
                self.normal_directions[collocation_panel_index] = collocation_panel.normal_direction

                wake_velocity = np.zeros(3)

                for wake_ring_vortex_panel_wing in self.airplane.wings:

                    wake_ring_vortex_panels = np.ravel(wake_ring_vortex_panel_wing.panels)

                    for wake_ring_vortex_panel in wake_ring_vortex_panels:

                        if wake_ring_vortex_panel.wake_ring_vortices is not None:
                            wake_velocity += wake_ring_vortex_panel.calculate_velocity_induced_by_wake_ring_vortices(
                                collocation_panel.collocation_point)

                # Update solver's list of freestream influences.
                self.freestream_influences[collocation_panel_index] = (
                    np.dot(self.freestream_velocity + wake_velocity, collocation_panel.normal_direction))

    def calculate_vortex_strengths(self):
        """Solve for each panel's vortex strength.

        :return:
        """

        # Solve for the strength of each panel's vortex.
        self.vortex_strengths = np.linalg.solve(self.aerodynamic_influence_coefficients, -self.freestream_influences)

        # Iterate through the airplane's wings.
        for wing in self.airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            wing_panels = np.ravel(wing.panels)

            # Iterate through this list of panels.
            for panel_index, panel in np.ndenumerate(wing_panels):

                # Update each panel's ring vortex strength.
                panel.ring_vortex.update_strength(self.vortex_strengths[panel_index])

                # If the panel has a horseshoe vortex, update its strength.
                if panel.horseshoe_vortex is not None:
                    panel.horseshoe_vortex.update_strength(self.vortex_strengths[panel_index])

    # ToDo: Properly document this method.
    def calculate_solution_velocity(self, point):
        """

        :param point:
        :return:
        """

        velocity_induced_by_vortices = np.zeros(3)

        for wing in self.airplane.wings:

            wing_panels = np.ravel(wing.panels)
            for panel in wing_panels:
                velocity_induced_by_vortices += panel.calculate_induced_velocity(point)
                if panel.wake_ring_vortices is not None:
                    velocity_induced_by_vortices += panel.calculate_velocity_induced_by_wake_ring_vortices(point)

        freestream = self.operating_point.calculate_freestream_velocity_geometry_axes()

        return velocity_induced_by_vortices + freestream

    # ToDo: Properly cite and document this method.
    def calculate_near_field_forces_and_moments(self):
        """

        :return:
        """
        total_near_field_force_geometry_axes = np.zeros(3)
        total_near_field_moment_geometry_axes = np.zeros(3)

        for wing in self.airplane.wings:
            panels = wing.panels
            num_chordwise_panels = wing.num_chordwise_panels
            num_spanwise_panels = wing.num_spanwise_panels
            for chordwise_location in range(num_chordwise_panels):
                for spanwise_location in range(num_spanwise_panels):
                    panel = panels[chordwise_location, spanwise_location]

                    is_right_edge = (spanwise_location + 1 == num_spanwise_panels)
                    is_leading_edge = panel.is_leading_edge
                    is_left_edge = (spanwise_location == 0)

                    if is_right_edge:
                        right_bound_vortex_strength = panel.ring_vortex.strength
                    else:
                        right_bound_vortex_strength = 0
                    if is_leading_edge:
                        front_bound_vortex_strength = panel.ring_vortex.strength
                    else:
                        front_bound_vortex_strength = (
                                panel.ring_vortex.strength
                                - panels[chordwise_location - 1, spanwise_location].ring_vortex.strength
                        )
                    if is_left_edge:
                        left_bound_vortex_strength = panel.ring_vortex.strength
                    else:
                        left_bound_vortex_strength = (
                                panel.ring_vortex.strength
                                - panels[chordwise_location, spanwise_location - 1].ring_vortex.strength
                        )

                    panel.near_field_force_geometry_axes = np.zeros(3)
                    panel.near_field_moment_geometry_axes = np.zeros(3)

                    if right_bound_vortex_strength != 0:
                        velocity_at_right_bound_vortex_center = self.calculate_solution_velocity(
                            panel.ring_vortex.right_leg.center)
                        right_bound_vortex_force_in_geometry_axes = (
                                self.operating_point.density
                                * right_bound_vortex_strength
                                * np.cross(velocity_at_right_bound_vortex_center, panel.ring_vortex.right_leg.vector)
                        )
                        panel.near_field_force_geometry_axes += right_bound_vortex_force_in_geometry_axes
                        panel.near_field_moment_geometry_axes += np.cross(
                            panel.ring_vortex.right_leg.center
                            - self.airplane.xyz_ref, right_bound_vortex_force_in_geometry_axes
                        )
                    if front_bound_vortex_strength != 0:
                        velocity_at_front_bound_vortex_center = self.calculate_solution_velocity(
                            panel.ring_vortex.front_leg.center)
                        front_bound_vortex_force_in_geometry_axes = (
                                self.operating_point.density
                                * front_bound_vortex_strength
                                * np.cross(velocity_at_front_bound_vortex_center, panel.ring_vortex.front_leg.vector)
                        )
                        panel.near_field_force_geometry_axes += front_bound_vortex_force_in_geometry_axes
                        panel.near_field_moment_geometry_axes += np.cross(
                            panel.ring_vortex.front_leg.center
                            - self.airplane.xyz_ref, front_bound_vortex_force_in_geometry_axes
                        )
                    if left_bound_vortex_strength != 0:
                        velocity_at_left_bound_vortex_center = self.calculate_solution_velocity(
                            panel.ring_vortex.left_leg.center)
                        left_bound_vortex_force_in_geometry_axes = (
                                self.operating_point.density
                                * left_bound_vortex_strength
                                * np.cross(velocity_at_left_bound_vortex_center, panel.ring_vortex.left_leg.vector)
                        )
                        panel.near_field_force_geometry_axes += left_bound_vortex_force_in_geometry_axes
                        panel.near_field_moment_geometry_axes += np.cross(
                            panel.ring_vortex.left_leg.center
                            - self.airplane.xyz_ref, left_bound_vortex_force_in_geometry_axes
                        )

                    panel.update_pressure()

                    total_near_field_force_geometry_axes += panel.near_field_force_geometry_axes
                    total_near_field_moment_geometry_axes += panel.near_field_moment_geometry_axes

        self.total_near_field_force_wind_axes = (
                np.transpose(self.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes())
                @ total_near_field_force_geometry_axes
        )
        self.total_near_field_moment_wind_axes = (
                np.transpose(self.operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes())
                @ total_near_field_moment_geometry_axes
        )

        self.CDi = (
                -self.total_near_field_force_wind_axes[0]
                / self.operating_point.calculate_dynamic_pressure()
                / self.airplane.s_ref
        )
        self.CY = (
                self.total_near_field_force_wind_axes[1]
                / self.operating_point.calculate_dynamic_pressure()
                / self.airplane.s_ref
        )
        self.CL = (
                -self.total_near_field_force_wind_axes[2]
                / self.operating_point.calculate_dynamic_pressure()
                / self.airplane.s_ref
        )
        self.Cl = (
                self.total_near_field_moment_wind_axes[0]
                / self.operating_point.calculate_dynamic_pressure()
                / self.airplane.s_ref
                / self.airplane.b_ref
        )
        self.Cm = (
                self.total_near_field_moment_wind_axes[1]
                / self.operating_point.calculate_dynamic_pressure()
                / self.airplane.s_ref
                / self.airplane.c_ref
        )
        self.Cn = (
                self.total_near_field_moment_wind_axes[2]
                / self.operating_point.calculate_dynamic_pressure()
                / self.airplane.s_ref
                / self.airplane.b_ref
        )

    # ToDo: Properly document this method.
    def calculate_streamlines(self):
        """

        :return:
        """

        airplane = self.airplane

        num_steps = 10
        delta_time = 1

        for wing in airplane.wings:
            wing.streamline_points = np.zeros((num_steps + 1, wing.num_spanwise_panels, 3))
            chordwise_position = wing.num_chordwise_panels - 1

            # Increment through the wing's chordwise and spanwise positions.
            for spanwise_position in range(wing.num_spanwise_panels):

                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]
                seed_point = panel.back_left_vertex + 0.5 * (panel.back_right_vertex - panel.back_left_vertex)
                wing.streamline_points[0, spanwise_position, :] = seed_point
                for step in range(num_steps):
                    last_point = wing.streamline_points[step, spanwise_position, :]

                    wing.streamline_points[step + 1, spanwise_position, :] = (
                            last_point
                            + delta_time
                            * self.calculate_solution_velocity(last_point)
                    )

    # ToDo: Properly document and cite this method.
    def calculate_wake_rollup(self):

        delta_time = 0.0125

        # Iterate through the airplane's wings.
        for wing in self.airplane.wings:

            chordwise_position = wing.num_chordwise_panels - 1

            # Iterate through the wing's chordwise and spanwise positions.
            for spanwise_position in range(wing.num_spanwise_panels):

                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]

                wake_ring_vortices = panel.wake_ring_vortices

                for wake_ring_vortex in wake_ring_vortices:
                    front_right_vertex = wake_ring_vortex.front_right_vertex
                    front_left_vertex = wake_ring_vortex.front_left_vertex
                    back_left_vertex = wake_ring_vortex.back_left_vertex
                    back_right_vertex = wake_ring_vortex.back_right_vertex

                    velocity_at_front_right_vertex = self.calculate_solution_velocity(front_right_vertex)
                    velocity_at_front_left_vertex = self.calculate_solution_velocity(front_left_vertex)
                    velocity_at_back_left_vertex = self.calculate_solution_velocity(back_left_vertex)
                    velocity_at_back_right_vertex = self.calculate_solution_velocity(back_right_vertex)

                    front_right_vertex += delta_time * velocity_at_front_right_vertex
                    front_left_vertex += delta_time * velocity_at_front_left_vertex
                    back_left_vertex += delta_time * velocity_at_back_left_vertex
                    back_right_vertex += delta_time * velocity_at_back_right_vertex

                    wake_ring_vortex.update_location(
                        front_right_vertex=front_right_vertex,
                        front_left_vertex=front_left_vertex,
                        back_left_vertex=back_left_vertex,
                        back_right_vertex=back_right_vertex
                    )

                new_wake_ring_vortex_front_right_vertex = panel.ring_vortex.back_right_vertex
                new_wake_ring_vortex_front_left_vertex = panel.ring_vortex.back_left_vertex

                velocity_at_front_right_vertex = self.calculate_solution_velocity(
                    new_wake_ring_vortex_front_right_vertex)
                velocity_at_front_left_vertex = self.calculate_solution_velocity(
                    new_wake_ring_vortex_front_left_vertex)

                new_wake_ring_vortex_back_left_vertex = (new_wake_ring_vortex_front_left_vertex + delta_time
                                                         * velocity_at_front_left_vertex)
                new_wake_ring_vortex_back_right_vertex = (new_wake_ring_vortex_front_right_vertex + delta_time
                                                          * velocity_at_front_right_vertex)

                new_wake_ring_vortex_strength = panel.ring_vortex.strength

                new_wake_ring_vortex = asmvp.aerodynamics.RingVortex(
                        front_right_vertex=new_wake_ring_vortex_front_right_vertex,
                        front_left_vertex=new_wake_ring_vortex_front_left_vertex,
                        back_left_vertex=new_wake_ring_vortex_back_left_vertex,
                        back_right_vertex=new_wake_ring_vortex_back_right_vertex,
                        strength=new_wake_ring_vortex_strength
                )

                panel.wake_ring_vortices = np.hstack((new_wake_ring_vortex, wake_ring_vortices))
