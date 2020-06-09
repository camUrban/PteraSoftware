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
        calculate_wing_wing_influences: Find the matrix of aerodynamic influence coefficients associated with this
                                        problem's geometry.
        calculate_freestream_wing_influences: Find the normal freestream speed at every collocation point without the
                                              influence of the vortices.
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

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.current_step = None
        self.current_airplane = None
        self.current_operating_point = None
        self.current_wing_wing_influences = None
        self.current_freestream_wing_influences = None
        self.current_freestream_velocity = None
        self.current_wake_wing_influences = None
        self.current_vortex_strengths = None
        self.current_total_near_field_force_wind_axes = None
        self.current_total_near_field_moment_wind_axes = None
        self.current_CL = None
        self.current_CDi = None
        self.current_CY = None
        self.current_Cl = None
        self.current_Cm = None
        self.current_Cn = None

    # ToDo: Properly document this method.
    def run(self):
        """Run the solver on the unsteady problem.

        :return: None
        """

        print("Initializing all airplanes' panel vortices...")
        for airplane in self.unsteady_problem.movement.airplanes:
            self.current_airplane = airplane
            self.initialize_panel_vortices()
        print("All airplanes' panel vortices initialized.")

        for step in range(self.unsteady_problem.movement.num_steps):
            self.current_step = step

            print("\nBeginning time current_step "
                  + str(self.current_step)
                  + " out of "
                  + str(self.unsteady_problem.movement.num_steps - 1)
                  + "...")

            self.current_airplane = self.unsteady_problem.movement.airplanes[self.current_step]
            self.current_operating_point = self.unsteady_problem.movement.operating_points[self.current_step]

            # Initialize attributes to hold aerodynamic data that pertains to this problem.
            self.current_wing_wing_influences = None
            self.current_freestream_wing_influences = None
            self.current_freestream_velocity = np.expand_dims(
                self.current_operating_point.calculate_freestream_velocity_geometry_axes(), 0
            )
            self.current_wake_wing_influences = None
            self.current_vortex_strengths = None
            self.current_total_near_field_force_wind_axes = None
            self.current_total_near_field_moment_wind_axes = None
            self.current_CL = None
            self.current_CDi = None
            self.current_CY = None
            self.current_Cl = None
            self.current_Cm = None
            self.current_Cn = None

            # Find the matrix of wing-wing influence coefficients associated with this current_airplane's geometry.
            print("Calculating the wing-wing influences...")
            self.calculate_wing_wing_influences()
            print("Wing-wing influences calculated.")

            # Find the vector of freestream-wing influence coefficients associated with this problem.
            print("Calculating the freestream-wing influences...")
            self.calculate_freestream_wing_influences()
            print("Freestream-wing influences calculated.")

            self.calculate_wake_wing_influences()

            # Solve for each panel's vortex strength.
            print("Calculating vortex strengths...")
            self.calculate_vortex_strengths()
            print("Vortex strengths calculated.")

            # Solve for the near field forces and moments on each panel.
            print("Shedding wake vortices...")
            self.calculate_wake_rollup()
            print("Wake vortices shed.")

            self.debug_wake_vortices()

            print("Finished time current_step "
                  + str(self.current_step)
                  + " out of "
                  + str(self.unsteady_problem.movement.num_steps - 1)
                  + ".")

        # Solve for the near field forces and moments on each panel.
        print("\nCalculating near field forces...")
        self.calculate_near_field_forces_and_moments()
        print("Near field forces calculated.")

        print('Final Airplane\'s wake vortex vertices:')
        asmvp.output.draw(self.current_airplane, show_delta_pressures=True)

        # Print out the total forces.
        print("\n\nTotal Forces in Wind Axes:")
        print("\tInduced Drag:\t\t\t", np.round(self.current_total_near_field_force_wind_axes[0], 3), " N")
        print("\tSide Force:\t\t\t\t", np.round(self.current_total_near_field_force_wind_axes[1], 3), " N")
        print("\tLift:\t\t\t\t\t", np.round(self.current_total_near_field_force_wind_axes[2], 3), " N")

        # Print out the total moments.
        print("\nTotal Moments in Wind Axes:")
        print("\tRolling Moment:\t\t\t", np.round(self.current_total_near_field_moment_wind_axes[0], 3), " Nm")
        print("\tPitching Moment:\t\t", np.round(self.current_total_near_field_moment_wind_axes[1], 3), " Nm")
        print("\tYawing Moment:\t\t\t", np.round(self.current_total_near_field_moment_wind_axes[2], 3), " Nm")

        # Print out the coefficients.
        print("\nCoefficients in Wind Axes:")
        print("\tcurrent_CDi:\t\t\t\t\t", np.round(self.current_CDi, 3))
        print("\tcurrent_CY:\t\t\t\t\t\t", np.round(self.current_CY, 3))
        print("\tcurrent_CL:\t\t\t\t\t\t", np.round(self.current_CL, 3))
        print("\tcurrent_Cl:\t\t\t\t\t\t", np.round(self.current_Cl, 3))
        print("\tcurrent_Cm:\t\t\t\t\t\t", np.round(self.current_Cm, 3))
        print("\tcurrent_Cn:\t\t\t\t\t\t", np.round(self.current_Cn, 3))

    # ToDo: Properly document this method.
    def initialize_panel_vortices(self):
        """This method calculates the locations of the vortex vertices, and then initializes the panel's vortices.

        Every panel has a ring vortex, which is a quadrangle whose front vortex leg is at the panel's quarter chord.
        The left and right vortex legs run along the panel's left and right legs. If the panel is not along the
        trailing edge, they extend backwards and meet the back vortex leg at a length of one quarter of the rear
        panel's chord back from the rear panel's front leg. Otherwise, they extend back backwards and meet the back
        vortex leg at a length of one quarter of the current panel's chord back from the current panel's back leg.

        :return: None
        """

        # Iterate through the current_airplane's wings.
        for wing in self.current_airplane.wings:

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

    # ToDo: Properly document this method.
    def calculate_wing_wing_influences(self):
        """Find the matrix of wing-wing influence coefficients associated with this current_airplane's geometry.

        :return: None
        """

        # Initialize the wing-wing influence coefficient matrix.
        self.current_wing_wing_influences = np.zeros(
            (self.current_airplane.num_panels, self.current_airplane.num_panels)
        )

        # Iterate through the current_airplane's wings. This wing contains the panel with the collocation point where
        # the vortex influence is to be calculated.
        for collocation_panel_wing in self.current_airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            collocation_panels = np.ravel(collocation_panel_wing.panels)

            # Iterate through the list of panels with the collocation points.
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panels):

                # Iterate through the current_airplane's wings. This wing contains the panel with the vortex whose
                # influence on the collocation point is to be calculated.
                for vortex_panel_wing in self.current_airplane.wings:

                    # Convert the 2D ndarray of this wing's panels into a 1D list.
                    vortex_panels = np.ravel(vortex_panel_wing.panels)

                    # Iterate through the list of panels with the vortices.
                    for vortex_panel_index, vortex_panel in np.ndenumerate(vortex_panels):

                        # Calculate the velocity induced at this collocation point by this vortex if the vortex's
                        # strength was 1.
                        normalized_induced_velocity_at_collocation_point = (
                            vortex_panel.calculate_normalized_induced_velocity(collocation_panel.collocation_point)
                        )

                        # Find the normal direction of the panel with the collocation point.
                        collocation_panel_normal_direction = collocation_panel.normal_direction

                        # Calculate the normal component of the velocity induced at this collocation point by this
                        # vortex if the vortex's strength was 1.
                        normal_normalized_induced_velocity_at_collocation_point = np.dot(
                            normalized_induced_velocity_at_collocation_point, collocation_panel_normal_direction
                        )

                        # Add this value to the solver's aerodynamic influence coefficient matrix.
                        self.current_wing_wing_influences[collocation_panel_index, vortex_panel_index] = (
                            normal_normalized_induced_velocity_at_collocation_point
                        )

    # ToDo: Properly document this method.
    def calculate_freestream_wing_influences(self):
        """Find the vector of freestream-wing influence coefficients associated with this problem.

        Note: This contains the influence due to the freestream and due to any geometry unsteady_problem.movement
        relative to the freestream.

        :return: None
        """

        # Initialize the vector of freestream-wing influence coefficients.
        self.current_freestream_wing_influences = np.zeros(self.current_airplane.num_panels)

        for collocation_panel_wing_position in range(len(self.current_airplane.wings)):

            collocation_panel_wing = self.current_airplane.wings[collocation_panel_wing_position]

            for collocation_panel_chordwise_position in range(collocation_panel_wing.num_chordwise_panels):

                for collocation_panel_spanwise_position in range(collocation_panel_wing.num_spanwise_panels):
                    collocation_panel_position = (collocation_panel_chordwise_position
                                                  * collocation_panel_wing.num_chordwise_panels
                                                  + collocation_panel_spanwise_position)

                    collocation_panel = collocation_panel_wing.panels[collocation_panel_chordwise_position,
                                                                      collocation_panel_spanwise_position]

                    collocation_panel_normal_direction = collocation_panel.normal_direction

                    velocity_induced_by_flapping_at_collocation_point = (
                        self.unsteady_problem.movement.get_flapping_velocity_at_point_on_panel(
                            wing_position=collocation_panel_wing_position,
                            panel_chordwise_position=collocation_panel_chordwise_position,
                            panel_spanwise_position=collocation_panel_spanwise_position,
                            point_name='collocation_point',
                            current_step=self.current_step)
                    )

                    self.current_freestream_wing_influences[collocation_panel_position] = (
                        np.dot(
                            (
                                    self.current_freestream_velocity
                                    + velocity_induced_by_flapping_at_collocation_point
                            ),
                            collocation_panel_normal_direction
                        )
                    )

    # ToDo: Properly document this method.
    def calculate_wake_wing_influences(self):

        self.current_wake_wing_influences = np.zeros(self.current_airplane.num_panels)

        for collocation_panel_wing_position in range(len(self.current_airplane.wings)):

            collocation_panel_wing = self.current_airplane.wings[collocation_panel_wing_position]

            for collocation_panel_chordwise_position in range(collocation_panel_wing.num_chordwise_panels):

                for collocation_panel_spanwise_position in range(collocation_panel_wing.num_spanwise_panels):
                    collocation_panel_position = (collocation_panel_chordwise_position
                                                  * collocation_panel_wing.num_chordwise_panels
                                                  + collocation_panel_spanwise_position)

                    collocation_panel = collocation_panel_wing.panels[collocation_panel_chordwise_position,
                                                                      collocation_panel_spanwise_position]

                    collocation_point = collocation_panel.collocation_point

                    collocation_panel_normal_direction = collocation_panel.normal_direction

                    wake_velocity = np.zeros(3)

                    for wake_ring_vortex_wing in self.current_airplane.wings:

                        wake_ring_vortices = np.ravel(wake_ring_vortex_wing.wake_ring_vortices)

                        for wake_ring_vortex in wake_ring_vortices:

                            if wake_ring_vortex is not None:
                                wake_velocity += (
                                    wake_ring_vortex.calculate_induced_velocity(collocation_point)
                                )

                    self.current_wake_wing_influences[collocation_panel_position] = (
                        np.dot(wake_velocity, collocation_panel_normal_direction)
                    )

    def calculate_vortex_strengths(self):
        """Solve for each panel's vortex strength.

        :return:
        """

        # Solve for the strength of each panel's vortex.
        self.current_vortex_strengths = np.linalg.solve(
            self.current_wing_wing_influences,
            - self.current_wake_wing_influences - self.current_freestream_wing_influences
        )

        # Iterate through the current_airplane's wings.
        for wing in self.current_airplane.wings:

            # Convert the 2D ndarray of this wing's panels into a 1D list.
            wing_panels = np.ravel(wing.panels)

            # Iterate through this list of panels.
            for panel_index, panel in np.ndenumerate(wing_panels):

                # Update each panel's ring vortex strength.
                panel.ring_vortex.update_strength(self.current_vortex_strengths[panel_index])

                # If the panel has a horseshoe vortex, update its strength.
                if panel.horseshoe_vortex is not None:
                    panel.horseshoe_vortex.update_strength(self.current_vortex_strengths[panel_index])

    # ToDo: Properly document this method.
    def calculate_solution_velocity(self, point):
        """

        :param point:
        :return:
        """

        velocity_induced_by_vortices = np.zeros(3)

        for wing in self.current_airplane.wings:

            wing_panels = np.ravel(wing.panels)
            for panel in wing_panels:
                velocity_induced_by_vortices += panel.calculate_induced_velocity(point)

            wake_ring_vortices = np.ravel(wing.wake_ring_vortices)
            for wake_ring_vortex in wake_ring_vortices:
                if wake_ring_vortex is not None:
                    velocity_induced_by_vortices += wake_ring_vortex.calculate_induced_velocity(point)

        return velocity_induced_by_vortices + self.current_freestream_velocity

    # ToDo: Properly cite and document this method.
    def calculate_near_field_forces_and_moments(self):
        """

        :return:
        """
        total_near_field_force_geometry_axes = np.zeros(3)
        total_near_field_moment_geometry_axes = np.zeros(3)

        for wing in self.current_airplane.wings:
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
                                self.current_operating_point.density
                                * right_bound_vortex_strength
                                * np.cross(velocity_at_right_bound_vortex_center, panel.ring_vortex.right_leg.vector)
                        )
                        panel.near_field_force_geometry_axes += right_bound_vortex_force_in_geometry_axes
                        panel.near_field_moment_geometry_axes += np.cross(
                            panel.ring_vortex.right_leg.center
                            - self.current_airplane.xyz_ref, right_bound_vortex_force_in_geometry_axes
                        )
                    if front_bound_vortex_strength != 0:
                        velocity_at_front_bound_vortex_center = self.calculate_solution_velocity(
                            panel.ring_vortex.front_leg.center)
                        front_bound_vortex_force_in_geometry_axes = (
                                self.current_operating_point.density
                                * front_bound_vortex_strength
                                * np.cross(velocity_at_front_bound_vortex_center, panel.ring_vortex.front_leg.vector)
                        )
                        panel.near_field_force_geometry_axes += front_bound_vortex_force_in_geometry_axes
                        panel.near_field_moment_geometry_axes += np.cross(
                            panel.ring_vortex.front_leg.center
                            - self.current_airplane.xyz_ref, front_bound_vortex_force_in_geometry_axes
                        )
                    if left_bound_vortex_strength != 0:
                        velocity_at_left_bound_vortex_center = self.calculate_solution_velocity(
                            panel.ring_vortex.left_leg.center)
                        left_bound_vortex_force_in_geometry_axes = (
                                self.current_operating_point.density
                                * left_bound_vortex_strength
                                * np.cross(velocity_at_left_bound_vortex_center, panel.ring_vortex.left_leg.vector)
                        )
                        panel.near_field_force_geometry_axes += left_bound_vortex_force_in_geometry_axes
                        panel.near_field_moment_geometry_axes += np.cross(
                            panel.ring_vortex.left_leg.center
                            - self.current_airplane.xyz_ref, left_bound_vortex_force_in_geometry_axes
                        )

                    panel.update_pressure()

                    total_near_field_force_geometry_axes += panel.near_field_force_geometry_axes
                    total_near_field_moment_geometry_axes += panel.near_field_moment_geometry_axes

        self.current_total_near_field_force_wind_axes = (
                np.transpose(self.current_operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes())
                @ total_near_field_force_geometry_axes
        )
        self.current_total_near_field_moment_wind_axes = (
                np.transpose(self.current_operating_point.calculate_rotation_matrix_wind_axes_to_geometry_axes())
                @ total_near_field_moment_geometry_axes
        )

        self.current_CDi = (
                -self.current_total_near_field_force_wind_axes[0]
                / self.current_operating_point.calculate_dynamic_pressure()
                / self.current_airplane.s_ref
        )
        self.current_CY = (
                self.current_total_near_field_force_wind_axes[1]
                / self.current_operating_point.calculate_dynamic_pressure()
                / self.current_airplane.s_ref
        )
        self.current_CL = (
                -self.current_total_near_field_force_wind_axes[2]
                / self.current_operating_point.calculate_dynamic_pressure()
                / self.current_airplane.s_ref
        )
        self.current_Cl = (
                self.current_total_near_field_moment_wind_axes[0]
                / self.current_operating_point.calculate_dynamic_pressure()
                / self.current_airplane.s_ref
                / self.current_airplane.b_ref
        )
        self.current_Cm = (
                self.current_total_near_field_moment_wind_axes[1]
                / self.current_operating_point.calculate_dynamic_pressure()
                / self.current_airplane.s_ref
                / self.current_airplane.c_ref
        )
        self.current_Cn = (
                self.current_total_near_field_moment_wind_axes[2]
                / self.current_operating_point.calculate_dynamic_pressure()
                / self.current_airplane.s_ref
                / self.current_airplane.b_ref
        )

    # ToDo: Properly document this method.
    def calculate_streamlines(self):
        """

        :return:
        """

        airplane = self.current_airplane

        num_steps = 10
        delta_time = 0.1

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
        """
        
        :return: 
        """

        self.calculate_next_airplanes_wake_vortex_vertices()
        self.calculate_next_airplanes_wake_vortices()

    # ToDo: Properly document this method.
    def calculate_next_airplanes_wake_vortex_vertices(self):

        delta_time = self.unsteady_problem.movement.delta_time
        step = self.current_step
        num_steps = self.unsteady_problem.movement.num_steps
        this_airplane = self.current_airplane

        if step < num_steps - 1:
            next_airplane = self.unsteady_problem.movement.airplanes[step + 1]
            num_wings = len(this_airplane.wings)

            for wing_num in range(num_wings):
                this_wing = this_airplane.wings[wing_num]
                next_wing = next_airplane.wings[wing_num]

                if step == 0:
                    num_spanwise_panels = this_wing.num_spanwise_panels
                    num_chordwise_panels = this_wing.num_chordwise_panels
                    chordwise_position = num_chordwise_panels - 1
                    new_row_of_wake_ring_vertex_vertices = np.zeros((1, num_spanwise_panels + 1, 3))

                    for spanwise_position in range(num_spanwise_panels):
                        next_panel = next_wing.panels[chordwise_position, spanwise_position]
                        next_front_left_vertex = (
                            next_panel.back_left_vertex
                            + 0.25 * self.current_freestream_velocity * self.unsteady_problem.movement.delta_time
                        )

                        new_row_of_wake_ring_vertex_vertices[0, spanwise_position] = next_front_left_vertex

                        if spanwise_position == (num_spanwise_panels - 1):
                            next_front_right_vertex = (
                                next_panel.back_right_vertex
                                + 0.25 * self.current_freestream_velocity * self.unsteady_problem.movement.delta_time
                            )
                            new_row_of_wake_ring_vertex_vertices[0, spanwise_position + 1] = next_front_right_vertex

                    next_wing.wake_ring_vortex_vertices = np.vstack((next_wing.wake_ring_vortex_vertices,
                                                                     np.copy(new_row_of_wake_ring_vertex_vertices)))

                    num_chordwise_vertices = next_wing.wake_ring_vortex_vertices.shape[0]
                    num_spanwise_vertices = next_wing.wake_ring_vortex_vertices.shape[1]

                    new_row_of_wake_ring_vertex_vertices = np.zeros((1, num_spanwise_panels + 1, 3))

                    for chordwise_vertex_position in range(num_chordwise_vertices):
                        for spanwise_vertex_position in range(num_spanwise_vertices):
                            wing_wake_ring_vortex_vertex = next_wing.wake_ring_vortex_vertices[
                                chordwise_vertex_position,
                                spanwise_vertex_position
                            ]
                            velocity_at_wake_vortex_vertices = self.calculate_solution_velocity(
                                wing_wake_ring_vortex_vertex)

                            new_row_of_wake_ring_vertex_vertices[0, spanwise_vertex_position] = (
                                    wing_wake_ring_vortex_vertex
                                    + velocity_at_wake_vortex_vertices
                                    * delta_time
                            )

                    next_wing.wake_ring_vortex_vertices = np.vstack((next_wing.wake_ring_vortex_vertices,
                                                                     new_row_of_wake_ring_vertex_vertices))

                else:
                    next_wing.wake_ring_vortex_vertices = np.copy(this_wing.wake_ring_vortex_vertices)
                    wing_wake_ring_vortex_vertices = next_wing.wake_ring_vortex_vertices
                    num_chordwise_vertices = wing_wake_ring_vortex_vertices.shape[0]
                    num_spanwise_vertices = wing_wake_ring_vortex_vertices.shape[1]

                    for chordwise_vertex_position in range(num_chordwise_vertices):
                        for spanwise_vertex_position in range(num_spanwise_vertices):
                            wing_wake_ring_vortex_vertex = wing_wake_ring_vortex_vertices[chordwise_vertex_position,
                                                                                          spanwise_vertex_position]
                            velocity_at_wake_vortex_vertices = self.calculate_solution_velocity(
                                wing_wake_ring_vortex_vertex)
                            next_wing.wake_ring_vortex_vertices[chordwise_vertex_position,
                                                                spanwise_vertex_position] += (
                                    velocity_at_wake_vortex_vertices * delta_time
                            )

                    num_spanwise_panels = this_wing.num_spanwise_panels
                    num_chordwise_panels = this_wing.num_chordwise_panels
                    chordwise_position = num_chordwise_panels - 1
                    new_row_of_wake_ring_vertex_vertices = np.empty((1, num_spanwise_panels + 1, 3))

                    for spanwise_position in range(num_spanwise_panels):
                        next_panel = next_wing.panels[chordwise_position, spanwise_position]
                        next_front_left_vertex = (
                            next_panel.back_left_vertex
                            + 0.25 * self.current_freestream_velocity * self.unsteady_problem.movement.delta_time
                        )

                        new_row_of_wake_ring_vertex_vertices[0, spanwise_position] = next_front_left_vertex

                        if spanwise_position == (num_spanwise_panels - 1):
                            next_front_right_vertex = (
                                next_panel.back_right_vertex
                                + 0.25 * self.current_freestream_velocity * self.unsteady_problem.movement.delta_time
                            )
                            new_row_of_wake_ring_vertex_vertices[0, spanwise_position + 1] = next_front_right_vertex

                    next_wing.wake_ring_vortex_vertices = np.vstack(
                        (new_row_of_wake_ring_vertex_vertices, next_wing.wake_ring_vortex_vertices)
                    )
        else:
            return

    # ToDo: Properly document this method.
    def calculate_next_airplanes_wake_vortices(self):
        step = self.current_step
        num_steps = self.unsteady_problem.movement.num_steps
        this_airplane = self.current_airplane

        if step < num_steps - 1:
            next_airplane = self.unsteady_problem.movement.airplanes[step + 1]
            num_wings = len(this_airplane.wings)

            for wing_num in range(num_wings):
                this_wing = this_airplane.wings[wing_num]
                next_wing = next_airplane.wings[wing_num]

                next_wing_wake_ring_vortex_vertices = next_wing.wake_ring_vortex_vertices

                this_wing_wake_ring_vortices = this_wing.wake_ring_vortices

                num_chordwise_vertices = next_wing_wake_ring_vortex_vertices.shape[0]
                num_spanwise_vertices = next_wing_wake_ring_vortex_vertices.shape[1]

                new_row_of_wake_ring_vortices = np.empty((1, num_spanwise_vertices - 1), dtype=object)
                next_wing.wake_ring_vortices = np.vstack((new_row_of_wake_ring_vortices, np.copy(
                    this_wing_wake_ring_vortices)))

                for chordwise_vertex_position in range(num_chordwise_vertices):
                    for spanwise_vertex_position in range(num_spanwise_vertices):

                        has_right_vertex = (spanwise_vertex_position + 1) < num_spanwise_vertices
                        has_back_vertex = (chordwise_vertex_position + 1) < num_chordwise_vertices

                        if has_right_vertex and has_back_vertex:
                            front_left_vertex = next_wing_wake_ring_vortex_vertices[
                                chordwise_vertex_position, spanwise_vertex_position]
                            front_right_vertex = next_wing_wake_ring_vortex_vertices[
                                chordwise_vertex_position, spanwise_vertex_position + 1]
                            back_left_vertex = next_wing_wake_ring_vortex_vertices[
                                chordwise_vertex_position + 1, spanwise_vertex_position]
                            back_right_vertex = next_wing_wake_ring_vortex_vertices[
                                chordwise_vertex_position + 1, spanwise_vertex_position + 1]

                            if chordwise_vertex_position > 0:
                                next_wing.wake_ring_vortices[chordwise_vertex_position,
                                                             spanwise_vertex_position].update_position(
                                    front_left_vertex=front_left_vertex,
                                    front_right_vertex=front_right_vertex,
                                    back_left_vertex=back_left_vertex,
                                    back_right_vertex=back_right_vertex
                                )

                            if chordwise_vertex_position == 0:
                                this_strength = this_wing.panels[this_wing.num_chordwise_panels - 1,
                                                                 spanwise_vertex_position].ring_vortex.strength
                                next_wing.wake_ring_vortices[chordwise_vertex_position, spanwise_vertex_position] = (
                                    asmvp.aerodynamics.RingVortex(
                                        front_left_vertex=front_left_vertex,
                                        front_right_vertex=front_right_vertex,
                                        back_left_vertex=back_left_vertex,
                                        back_right_vertex=back_right_vertex,
                                        strength=this_strength
                                    )
                                )
        else:
            return

    # ToDo: Properly document this method.
    def debug_wake_vortices(self):
        for wing in self.current_airplane.wings:
            print("\nWing Vortex Strengths: ")
            print("\n\tPanel Vortex Strengths:\n")
            for chordwise_position in range(wing.num_chordwise_panels):
                for spanwise_position in range(wing.num_spanwise_panels):
                    panel = wing.panels[chordwise_position, spanwise_position]
                    ring_vortex = panel.ring_vortex
                    strength = ring_vortex.strength
                    print("\t\t" + str(round(strength, 2)), end="\t")
                print()
            print("\n\tWake Vortex Strengths:\n")
            for chordwise_position in range(wing.wake_ring_vortices.shape[0]):
                for spanwise_position in range(wing.wake_ring_vortices.shape[1]):
                    ring_vortex = wing.wake_ring_vortices[chordwise_position, spanwise_position]
                    strength = ring_vortex.strength
                    print("\t\t" + str(round(strength, 2)), end="\t")
                print()
        print()
