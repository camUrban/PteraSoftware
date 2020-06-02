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
        self.movement = self.unsteady_problem.movement

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.step = None
        self.airplane = None
        self.operating_point = None
        self.aerodynamic_influence_coefficients = None
        self.freestream_velocity = None
        self.normal_directions = None
        self.freestream_influences = None
        self.vortex_strengths = None
        self.total_near_field_force_wind_axes = None
        self.total_near_field_moment_wind_axes = None
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

        print("Initializing all airplanes' panel vortices...")
        for airplane in self.movement.airplanes:
            self.airplane = airplane
            # Initialize this problem's panels to have vortices congruent with this solver type.
            self.initialize_panel_vortices()
        print("All airplanes' panel vortices initialized.")

        for step in range(self.movement.num_steps):

            self.step = step

            print("\nBeginning time step " + str(self.step) + " out of " + str(self.movement.num_steps - 1) + "...")

            self.airplane = self.movement.airplanes[self.step]
            self.operating_point = self.movement.operating_points[self.step]

            # Initialize attributes to hold aerodynamic data that pertains to this problem.
            self.aerodynamic_influence_coefficients = None
            self.freestream_velocity = None
            self.normal_directions = None
            self.freestream_influences = None
            self.vortex_strengths = None
            self.total_near_field_force_wind_axes = None
            self.total_near_field_moment_wind_axes = None
            self.CL = None
            self.CDi = None
            self.CY = None
            self.Cl = None
            self.Cm = None
            self.Cn = None

            # Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
            print("Setting up geometry...")
            self.set_up_geometry()
            print("Geometry set up.")

            # Find the normal freestream speed at every collocation point without vortices.
            print("Setting up operating point...")
            self.set_up_operating_point()
            print("Operating point set up.")

            # Solve for each panel's vortex strength.
            print("Calculating vortex strengths...")
            self.calculate_vortex_strengths()
            print("Vortex strengths calculated.")

            if self.step < self.movement.num_steps:
                # Solve for the near field forces and moments on each panel.
                print("Shedding wake vortices...")
                self.calculate_wake_rollup()
                print("Wake vortices shed.")

            print("Finished time step " + str(self.step) + " out of " + str(self.movement.num_steps - 1) + ".")

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

        self.aerodynamic_influence_coefficients = np.zeros((self.airplane.num_panels, self.airplane.num_panels))

        for collocation_panel_wing_position in range(len(self.airplane.wings)):

            collocation_panel_wing = self.airplane.wings[collocation_panel_wing_position]

            for collocation_panel_chordwise_position in range(collocation_panel_wing.num_chordwise_panels):

                for collocation_panel_spanwise_position in range(collocation_panel_wing.num_spanwise_panels):

                    collocation_panel_position = (collocation_panel_chordwise_position
                                                  * collocation_panel_wing.num_chordwise_panels
                                                  + collocation_panel_spanwise_position)

                    collocation_panel = collocation_panel_wing.panels[collocation_panel_chordwise_position,
                                                                      collocation_panel_spanwise_position]

                    collocation_point = collocation_panel.collocation_point

                    collocation_panel_normal_direction = collocation_panel.normal_direction

                    for vortex_panel_wing_position in range(len(self.airplane.wings)):

                        vortex_panel_wing = self.airplane.wings[vortex_panel_wing_position]

                        for vortex_panel_chordwise_position in range(vortex_panel_wing.num_chordwise_panels):

                            for vortex_panel_spanwise_position in range(vortex_panel_wing.num_spanwise_panels):
                                vortex_panel_position = (vortex_panel_chordwise_position
                                                         * vortex_panel_wing.num_chordwise_panels
                                                         + vortex_panel_spanwise_position)

                                vortex_panel = vortex_panel_wing.panels[vortex_panel_chordwise_position,
                                                                        vortex_panel_spanwise_position]

                                normalized_induced_velocity_at_collocation_point = (
                                    vortex_panel.calculate_normalized_induced_velocity(collocation_point)
                                )

                                normal_normalized_velocity_at_collocation_point = np.dot(
                                    normalized_induced_velocity_at_collocation_point,
                                    collocation_panel_normal_direction
                                )

                                self.aerodynamic_influence_coefficients[collocation_panel_position,
                                                                        vortex_panel_position] = (
                                    normal_normalized_velocity_at_collocation_point
                                )

    def set_up_operating_point(self):

        self.normal_directions = np.zeros((self.airplane.num_panels, 3))

        self.freestream_influences = np.zeros(self.airplane.num_panels)

        self.freestream_velocity = np.expand_dims(self.operating_point.calculate_freestream_velocity_geometry_axes(), 0)

        for collocation_panel_wing_position in range(len(self.airplane.wings)):

            collocation_panel_wing = self.airplane.wings[collocation_panel_wing_position]

            for collocation_panel_chordwise_position in range(collocation_panel_wing.num_chordwise_panels):

                for collocation_panel_spanwise_position in range(collocation_panel_wing.num_spanwise_panels):

                    collocation_panel_position = (collocation_panel_chordwise_position
                                                  * collocation_panel_wing.num_chordwise_panels
                                                  + collocation_panel_spanwise_position)

                    collocation_panel = collocation_panel_wing.panels[collocation_panel_chordwise_position,
                                                                      collocation_panel_spanwise_position]

                    collocation_point = collocation_panel.collocation_point

                    collocation_panel_normal_direction = collocation_panel.normal_direction

                    self.normal_directions[collocation_panel_position] = collocation_panel.normal_direction

                    wake_velocity = np.zeros(3)

                    # for wake_ring_vortex_wing in self.airplane.wings:
                    #
                    #     wake_ring_vortices = np.ravel(wake_ring_vortex_wing.wake_ring_vortices)
                    #
                    #     for wake_ring_vortex in wake_ring_vortices:
                    #
                    #         if wake_ring_vortex is not None:
                    #             wake_velocity += (
                    #                 wake_ring_vortex.calculate_normalized_induced_velocity(collocation_point)
                    #             )

                    velocity_induced_by_flapping_at_collocation_point = (
                        self.movement.get_flapping_velocity_at_point_on_panel(
                            wing_position=collocation_panel_wing_position,
                            panel_chordwise_position=collocation_panel_chordwise_position,
                            panel_spanwise_position=collocation_panel_spanwise_position,
                            point_name='collocation',
                            current_step=self.step)
                    )

                    self.freestream_influences[collocation_panel_position] = (
                        np.dot(
                            (
                                    self.freestream_velocity
                                    + wake_velocity
                                    + velocity_induced_by_flapping_at_collocation_point
                            ),
                            collocation_panel_normal_direction
                        )
                    )

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

            # wake_ring_vortices = np.ravel(wing.wake_ring_vortices)
            # for wake_ring_vortex in wake_ring_vortices:
            #     if wake_ring_vortex is not None:
            #         velocity_induced_by_vortices += wake_ring_vortex.calculate_normalized_induced_velocity(point)

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
        if self.step < self.movement.num_steps:
            print('\n')
            self.calculate_next_airplanes_wake_vortex_vertices()
            print('\n')
        else:
            raise Exception('There is no need to modify the wake, because this is the last time step.')
    
    # ToDo: Properly document and cite this method.
    def convect_wake_ring_vortices(self):
        """
        
        :return: 
        """

        delta_time = self.movement.delta_time

        airplane = self.airplane

        if self.step < self.movement.num_steps:
            next_airplane = self.movement.airplanes[self.step + 1]
        else:
            raise Exception('There is no need to convect the wake, because this is the last time step.')

        for wing_num in range(len(airplane.wings)):

            wing = airplane.wings[wing_num]
            next_wing = next_airplane.wings[wing_num]

            num_wake_ring_vortex_rows = next_wing.wake_ring_vortices.shape[0]
            num_wake_ring_vortex_columns = next_wing.num_spanwise_panels

            for wake_ring_vortex_row in range(num_wake_ring_vortex_rows):
                for wake_ring_vortex_column in range(num_wake_ring_vortex_columns):

                    next_wake_ring_vortex = next_wing.wake_ring_vortices[wake_ring_vortex_row, wake_ring_vortex_column]

                    next_wake_ring_vortex_front_right_vertex = next_wake_ring_vortex.front_right_vertex
                    next_wake_ring_vortex_front_left_vertex = next_wake_ring_vortex.front_left_vertex
                    next_wake_ring_vortex_back_left_vertex = next_wake_ring_vortex.back_left_vertex
                    next_wake_ring_vortex_back_right_vertex = next_wake_ring_vortex.back_right_vertex

                    velocity_at_next_wake_ring_vortex_front_right_vertex = self.calculate_solution_velocity(
                        next_wake_ring_vortex_front_right_vertex
                    )
                    velocity_at_next_wake_ring_vortex_front_left_vertex = self.calculate_solution_velocity(
                        next_wake_ring_vortex_front_left_vertex
                    )
                    velocity_at_next_wake_ring_vortex_back_left_vertex = self.calculate_solution_velocity(
                        next_wake_ring_vortex_back_left_vertex
                    )
                    velocity_at_next_wake_ring_vortex_back_right_vertex = self.calculate_solution_velocity(
                        next_wake_ring_vortex_back_right_vertex
                    )

                    new_next_wake_ring_vortex_front_right_vertex = next_wake_ring_vortex_front_right_vertex + (
                            velocity_at_next_wake_ring_vortex_front_right_vertex * delta_time)
                    new_next_wake_ring_vortex_front_left_vertex = next_wake_ring_vortex_front_left_vertex + (
                            velocity_at_next_wake_ring_vortex_front_left_vertex * delta_time)
                    new_next_wake_ring_vortex_back_left_vertex = next_wake_ring_vortex_back_left_vertex + (
                            velocity_at_next_wake_ring_vortex_back_left_vertex * delta_time)
                    new_next_wake_ring_vortex_back_right_vertex = next_wake_ring_vortex_back_right_vertex + (
                            velocity_at_next_wake_ring_vortex_back_right_vertex * delta_time)

                    next_wake_ring_vortex.update_position(
                        front_left_vertex=new_next_wake_ring_vortex_front_left_vertex,
                        front_right_vertex=new_next_wake_ring_vortex_front_right_vertex,
                        back_right_vertex=new_next_wake_ring_vortex_back_right_vertex,
                        back_left_vertex=new_next_wake_ring_vortex_back_left_vertex
                    )

        print("Convecting wake vortices:")
        asmvp.output.draw(next_airplane, show_delta_pressures=False)

    # ToDo: Properly document and cite this method.
    def shed_new_wake_ring_vortices(self):
        """
        
        :return: 
        """

        delta_time = self.movement.delta_time

        airplane = self.airplane

        if self.step < self.movement.num_steps:
            next_airplane = self.movement.airplanes[self.step + 1]
        else:
            raise Exception('There is no need to shed new panels into the wake, because this is the last time step.')

        for wing_num in range(len(airplane.wings)):

            wing = airplane.wings[wing_num]
            next_wing = next_airplane.wings[wing_num]

            chordwise_position = wing.num_chordwise_panels - 1

            new_wake_ring_vortices = np.empty(next_wing.num_spanwise_panels, dtype=object)

            for spanwise_position in range(next_wing.num_spanwise_panels):

                panel = wing.panels[chordwise_position, spanwise_position]
                next_panel = next_wing.panels[chordwise_position, spanwise_position]

                ring_vortex = panel.ring_vortex
                next_ring_vortex = next_panel.ring_vortex

                new_wake_ring_vortex_front_left_vertex = next_ring_vortex.back_left_vertex
                new_wake_ring_vortex_front_right_vertex = next_ring_vortex.back_right_vertex

                velocity_at_new_wake_ring_vortex_front_right_vertex = self.calculate_solution_velocity(
                    new_wake_ring_vortex_front_right_vertex
                )
                velocity_at_new_wake_ring_vortex_front_left_vertex = self.calculate_solution_velocity(
                    new_wake_ring_vortex_front_left_vertex
                )

                new_wake_ring_vortex_back_left_vertex = (next_ring_vortex.front_left_vertex
                                                         + velocity_at_new_wake_ring_vortex_front_left_vertex
                                                         * delta_time)
                new_wake_ring_vortex_back_right_vertex = (next_ring_vortex.front_right_vertex
                                                          + velocity_at_new_wake_ring_vortex_front_right_vertex
                                                          * delta_time)

                new_wake_ring_vortex_strength = ring_vortex.strength

                new_wake_ring_vortex = asmvp.aerodynamics.RingVortex(
                    front_left_vertex=new_wake_ring_vortex_front_left_vertex,
                    front_right_vertex=new_wake_ring_vortex_front_right_vertex,
                    back_left_vertex=new_wake_ring_vortex_back_left_vertex,
                    back_right_vertex=new_wake_ring_vortex_back_right_vertex,
                    strength=new_wake_ring_vortex_strength
                )

                new_wake_ring_vortices[spanwise_position] = new_wake_ring_vortex

            next_wing.wake_ring_vortices = np.vstack((new_wake_ring_vortices, next_wing.wake_ring_vortices))
            print("Shedding wake vortices:")
            asmvp.output.draw(next_airplane, show_delta_pressures=False)

    def calculate_next_airplanes_wake_vortex_vertices(self):

        delta_time = self.movement.delta_time
        step = self.step
        num_steps = self.movement.num_steps
        this_airplane = self.airplane

        if step < num_steps - 1:
            next_airplane = self.movement.airplanes[step + 1]
            num_wings = len(this_airplane.wings)

            for wing_num in range(num_wings):
                this_wing = this_airplane.wings[wing_num]
                next_wing = next_airplane.wings[wing_num]
                next_wing.wake_ring_vortex_vertices = np.copy(this_wing.wake_ring_vortex_vertices)

                if step == 0:
                    # Shed wake vortex vertices:
                    num_spanwise_panels = this_wing.num_spanwise_panels
                    num_chordwise_panels = this_wing.num_chordwise_panels
                    chordwise_position = num_chordwise_panels - 1
                    new_row_of_wake_ring_vertex_vertices = np.empty((1, num_spanwise_panels + 1, 3))

                    for spanwise_position in range(num_spanwise_panels):
                        next_panel = next_wing.panels[chordwise_position, spanwise_position]
                        next_ring_vortex = next_panel.ring_vortex
                        next_front_left_vertex = next_ring_vortex.back_left_vertex

                        new_row_of_wake_ring_vertex_vertices[0, spanwise_position] = next_front_left_vertex

                        if spanwise_position == (num_spanwise_panels - 1):
                            next_front_right_vertex = next_ring_vortex.back_right_vertex
                            new_row_of_wake_ring_vertex_vertices[0, spanwise_position + 1] = next_front_right_vertex

                    next_wing.wake_ring_vortex_vertices = np.vstack(
                        (new_row_of_wake_ring_vertex_vertices, next_wing.wake_ring_vortex_vertices)
                    )

                else:

                    wing_wake_ring_vortex_vertices = np.reshape(np.copy(next_wing.wake_ring_vortex_vertices), (-1, 3))

                    for wing_wake_ring_vortex_vertex_position in range(wing_wake_ring_vortex_vertices.shape[0]):

                        wing_wake_ring_vortex_vertex = wing_wake_ring_vortex_vertices[wing_wake_ring_vortex_vertex_position]

                        velocity_at_wake_vortex_vertex = self.calculate_solution_velocity(wing_wake_ring_vortex_vertex)

                        next_wing.wake_ring_vortex_vertices[:, :] += velocities_at_wake_vortex_vertices * delta_time

                    num_spanwise_panels = this_wing.num_spanwise_panels
                    num_chordwise_panels = this_wing.num_chordwise_panels
                    chordwise_position = num_chordwise_panels - 1
                    new_row_of_wake_ring_vertex_vertices = np.empty((1, num_spanwise_panels + 1, 3))

                    for spanwise_position in range(num_spanwise_panels):
                        next_panel = next_wing.panels[chordwise_position, spanwise_position]
                        next_ring_vortex = next_panel.ring_vortex
                        next_front_left_vertex = next_ring_vortex.back_left_vertex

                        new_row_of_wake_ring_vertex_vertices[0, spanwise_position] = next_front_left_vertex

                        if spanwise_position == (num_spanwise_panels - 1):
                            next_front_right_vertex = next_ring_vortex.back_right_vertex
                            new_row_of_wake_ring_vertex_vertices[0, spanwise_position + 1] = next_front_right_vertex

                    next_wing.wake_ring_vortex_vertices = np.vstack(
                        (new_row_of_wake_ring_vertex_vertices, next_wing.wake_ring_vortex_vertices)
                    )

            print('Latest wake vortex vertices:')
            asmvp.output.draw(next_airplane, show_delta_pressures=False)

        else:
            return
