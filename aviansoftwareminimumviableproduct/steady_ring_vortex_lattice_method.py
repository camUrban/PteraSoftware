"""This module contains the class definition of this package's steady ring vortex lattice solver.

This module contains the following classes:
    SteadyRingVortexLatticeMethodSolver: This is an aerodynamics solver that uses a steady ring vortex lattice method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly cite document this method.
class SteadyRingVortexLatticeMethodSolver:
    """This is an aerodynamics solver that uses a steady ring vortex lattice method.

    Citation:
        Adapted from:         aerodynamics.vlm3.py in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    This class contains the following public methods:
        run: Run the solver on the steady problem.
        set_up_geometry: Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
        set_up_operating_point: Find the normal freestream speed at every collocation point without vortices.
        calculate_vortex_strengths: Solve for each panel's vortex strength.
        calculate_solution_velocity: Find the velocity at a given point due to the freestream and the vortices.
        calculate_velocity_influences: Find the velocity at a given point due to the vorticity of every vortex if their
                                       strengths were all set to 1.0 meters squared per second.
        calculate_delta_cp: Find the change in the pressure coefficient between the upper and lower surfaces of a panel.
        calculate_near_field_forces_and_moments: Find the the forces and moments calculated from the near field.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # ToDo: Properly cite document this method.
    def __init__(self, steady_problem):
        """This is the initialization method.

        :param steady_problem: SteadyProblem
            This is the steady problem's to be solved.
        """

        self.steady_problem = steady_problem
        self.airplane = self.steady_problem.airplane
        self.operating_point = self.steady_problem.operating_point

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.aerodynamic_influence_coefficients = np.zeros((self.airplane.num_panels, self.airplane.num_panels))
        self.downwash_influence_coefficients = np.zeros((self.airplane.num_panels, self.airplane.num_panels))
        self.freestream_velocity = np.zeros(3)
        self.normal_directions = np.zeros((self.airplane.num_panels, 3))
        self.freestream_influences = np.zeros(self.airplane.num_panels)
        self.vortex_strengths = np.zeros(self.airplane.num_panels)
        self.downwashes = np.zeros(self.airplane.num_panels)

    # ToDo: Properly document this method.
    def run(self):
        """Run the solver on the steady problem.

        :return: None
        """

        # Initialize this problem's panels to have vortices congruent with this solver type.
        self.initialize_panel_vortices()
        print("Panel vortices initialized.")
        # Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
        self.set_up_geometry()
        print("Geometry set up.")
        # Find the normal freestream speed at every collocation point without vortices.
        self.set_up_operating_point()
        print("Operating point set up.")
        # Solve for each panel's vortex strength.
        self.calculate_vortex_strengths()
        print("Vortex strength's calculated.")
        self.calculate_near_field_forces()
        print("Near field forces calculated.")
        self.calculate_streamlines()
        print("Streamlines calculated.")

    # ToDo: Properly document this method.
    def initialize_panel_vortices(self):
        """This method calculates the locations of the vortex vertices, and then initializes the panel's vortices.

        This function takes in the type of problem this panel will be used in, and initializes the appropriate vortices.

        For the "steady horseshoe vortex lattice method" problem type:
            Every panel has a horseshoe vortex. The vortex's finite leg runs along the panel's quarter chord from right
            to left. It's infinite legs point backwards in the positive x direction.

        For the "steady ring vortex lattice method" problem type:
            The panel's ring vortex is a quadrangle whose front vortex leg is at the panel's quarter chord. The left and
            right vortex legs run along the panel's left and right legs. They extend backwards and meet the back vortex
            leg at one quarter chord back from the panel's back leg.

            Panels that are at the trailing edge of a wing have a horseshoe vortex in addition to their ring vortex. The
            horseshoe vortex's finite leg runs along the ring vortex's back leg but in the opposite direction. It's
            infinite legs point backwards in the positive x direction. The ring vortex and horseshoe vortex have the
            same strength, so the back leg of the ring vortex is cancelled.

        :return: None
        """

        for wing in self.airplane.wings:
            # Increment through the wing's chordwise and spanwise positions.
            for chordwise_position in range(wing.num_chordwise_panels):
                for spanwise_position in range(wing.num_spanwise_panels):
                    # Pull the panel object out of the wing's list of panels.
                    panel = wing.panels[chordwise_position, spanwise_position]

                    front_left_vortex_vertex = panel.front_left_vortex_vertex
                    front_right_vortex_vertex = panel.front_right_vortex_vertex

                    if not panel.is_trailing_edge:
                        next_chordwise_panel = wing.panels[chordwise_position + 1, spanwise_position]
                        back_left_vortex_vertex = next_chordwise_panel.front_left_vortex_vertex
                        back_right_vortex_vertex = next_chordwise_panel.front_right_vortex_vertex
                    else:
                        back_left_vortex_vertex = front_left_vortex_vertex + (
                                panel.back_left_vertex - panel.front_left_vertex)
                        back_right_vortex_vertex = front_right_vortex_vertex + (
                                panel.back_right_vertex - panel.front_right_vertex)
                        panel.horseshoe_vortex = asmvp.aerodynamics.HorseshoeVortex(
                                finite_leg_origin=back_right_vortex_vertex,
                                finite_leg_termination=back_left_vortex_vertex,
                                strength=None
                            )

                    # If the panel has a ring vortex, initialize it.
                    panel.ring_vortex = asmvp.aerodynamics.RingVortex(
                        front_left_vertex=front_left_vortex_vertex,
                        front_right_vertex=front_right_vortex_vertex,
                        back_left_vertex=back_left_vortex_vertex,
                        back_right_vertex=back_right_vortex_vertex,
                        strength=None
                    )

    # ToDo: Properly document this method.
    def set_up_geometry(self):

        for collocation_panel_wing in self.airplane.wings:

            collocation_panel_wings_panels = np.ravel(collocation_panel_wing.panels)
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panel_wings_panels):

                for vortex_panel_wing in self.airplane.wings:

                    vortex_panel_wings_panels = np.ravel(vortex_panel_wing.panels)
                    for vortex_panel_index, vortex_panel in np.ndenumerate(vortex_panel_wings_panels):
                        normalized_induced_velocity_at_collocation_point = (
                            vortex_panel.calculate_normalized_induced_velocity(collocation_panel.collocation_point))

                        collocation_panel_normal_direction = collocation_panel.normal_direction

                        normal_normalized_induced_velocity_at_collocation_point = np.dot(
                            normalized_induced_velocity_at_collocation_point, collocation_panel_normal_direction)

                        self.aerodynamic_influence_coefficients[collocation_panel_index, vortex_panel_index] = (
                            normal_normalized_induced_velocity_at_collocation_point)

    # ToDo: Properly document this method.
    def set_up_operating_point(self):
        # This calculates and updates the direction the wind is going to, in geometry axes coordinates.
        self.freestream_velocity = np.expand_dims(
            self.operating_point.calculate_freestream_velocity_geometry_axes(),
            0)

        for collocation_panel_wing in self.airplane.wings:

            collocation_panel_wings_panels = np.ravel(collocation_panel_wing.panels)
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panel_wings_panels):
                self.normal_directions[collocation_panel_index] = (
                    collocation_panel.normal_direction)
                self.freestream_influences[collocation_panel_index] = (
                    np.dot(self.freestream_velocity, collocation_panel.normal_direction))

    # ToDo: Properly document this method.
    def calculate_vortex_strengths(self):
        """

        :return:
        """

        self.vortex_strengths = np.linalg.solve(self.aerodynamic_influence_coefficients, -self.freestream_influences)
        self.downwashes = self.downwash_influence_coefficients @ self.vortex_strengths
        for wing in self.airplane.wings:

            wing_panels = np.ravel(wing.panels)
            for panel_index, panel in np.ndenumerate(wing_panels):
                panel.ring_vortex.update_strength(self.vortex_strengths[panel_index])
                panel.downwash = self.downwashes[panel_index]
                if panel.horseshoe_vortex is not None:
                    panel.horseshoe_vortex.update_strength(self.vortex_strengths[panel_index])

    # ToDo: Properly document this method.
    def calculate_solution_velocity(self, point):
        velocity_induced_by_vortices = np.zeros(3)

        for wing in self.airplane.wings:

            wing_panels = np.ravel(wing.panels)
            for panel in wing_panels:
                velocity_induced_by_vortices += panel.calculate_induced_velocity(point)

        freestream = self.operating_point.calculate_freestream_velocity_geometry_axes()

        return velocity_induced_by_vortices + freestream

    # ToDo: Properly cite and document this method.
    def calculate_near_field_forces(self):

        airplane = self.airplane
        density = self.operating_point.density

        for wing in airplane.wings:
            panels = wing.panels
            num_chordwise_panels = wing.num_chordwise_panels
            num_spanwise_panels = wing.num_spanwise_panels
            for chordwise_location in range(num_chordwise_panels):
                for spanwise_location in range(num_spanwise_panels):
                    panel = panels[chordwise_location, spanwise_location]
                    ring_vortex = panel.ring_vortex
                    horseshoe_vortex = panel.horseshoe_vortex
                    bound_vortices = np.array([ring_vortex.front_leg, ring_vortex.left_leg, ring_vortex.back_leg,
                                               ring_vortex.right_leg])
                    if horseshoe_vortex is not None:
                        bound_vortices = np.append(bound_vortices, horseshoe_vortex.finite_leg)

                    for bound_vortex in bound_vortices:
                        velocity_at_center = self.calculate_solution_velocity(bound_vortex.center)
                        bound_vortex.near_field_force = (density * bound_vortex.strength
                                                         * np.cross(velocity_at_center, bound_vortex.vector))

                        bound_vortex.near_field_moment = np.cross(bound_vortex.near_field_force, bound_vortex.center)

                    panel.update_force_moment_and_pressure()

    def calculate_streamlines(self):

        airplane = self.airplane

        num_steps = 100
        delta_time = 0.01

        for wing in airplane.wings:
            wing.stream_line_points = np.zeros((num_steps + 1, wing.num_spanwise_panels, 3))
            chordwise_position = wing.num_chordwise_panels - 1
            # Increment through the wing's chordwise and spanwise positions.
            for spanwise_position in range(wing.num_spanwise_panels):
                # Pull the panel object out of the wing's list of panels.
                panel = wing.panels[chordwise_position, spanwise_position]
                seed_point = panel.back_left_vertex + 0.5 * (panel.back_right_vertex - panel.back_left_vertex)
                wing.stream_line_points[0, spanwise_position, :] = seed_point
                for step in range(num_steps):
                    last_point = wing.stream_line_points[step, spanwise_position, :]

                    wing.stream_line_points[step + 1, spanwise_position, :] = (
                            last_point
                            + delta_time
                            * self.calculate_solution_velocity(last_point)
                    )
