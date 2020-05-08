"""This module contains the class definition of this package's steady vortex lattice solver.

This module contains the following classes:
    SteadyRingVortexLatticeMethodSolver: This is an aerodynamics solver that uses a steady vortex lattice method.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly cite and document this method.
class SteadyRingVortexLatticeMethodSolver(asmvp.problems.SteadyProblem):
    """This is an aerodynamics solver that uses a steady vortex lattice method.

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

    # ToDo: Properly cite and document this method.
    def __init__(self, airplane, operating_point):
        """This is the initialization method.

        :param airplane: Airplane
            This is the problem's airplane object to be analyzed.
        :param operating_point: OperatingPoint
            This is the problem's operating point object to be analyzed.
        """

        # Call the parent class initialization method.
        super().__init__(airplane, operating_point)

        # Initialize attributes to hold aerodynamic data that pertains to this problem.
        self.aerodynamic_influence_coefficients = np.zeros((self.airplane.num_panels, self.airplane.num_panels))
        self.downwash_influence_coefficients = np.zeros((self.airplane.num_panels, self.airplane.num_panels))
        self.freestream_velocity = np.zeros(3)
        self.normal_directions = np.zeros((self.airplane.num_panels, 3))
        self.freestream_influences = np.zeros(self.airplane.num_panels)
        self.vortex_strengths = np.zeros(self.airplane.num_panels)
        self.downwashes = np.zeros(self.airplane.num_panels)

    # ToDo: Properly cite and document this method.
    def run(self):
        """Run the solver on the steady problem.

        :return: None
        """

        # Find the matrix of aerodynamic influence coefficients associated with this problem's geometry.
        self.set_up_geometry()

        # Find the normal freestream speed at every collocation point without vortices.
        self.set_up_operating_point()

        # Solve for each panel's vortex strength.
        self.calculate_vortex_strengths()

        origin = np.zeros(3)
        print(self.calculate_solution_velocity(origin))

        self.calculate_near_field_forces()

        # Find the the forces and moments calculated from the near field.
        # self.calculate_near_field_forces_and_moments()

        # Find the change in the pressure coefficient between the upper and lower surfaces of a panel.
        # self.calculate_delta_cp()

    # ToDo: Properly cite and document this method.
    def set_up_geometry(self):

        for collocation_panel_wing in self.airplane.wings:

            collocation_panel_wings_panels = np.ravel(collocation_panel_wing.panels)
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panel_wings_panels):

                for vortex_panel_wing in self.airplane.wings:

                    vortex_panel_wings_panels = np.ravel(vortex_panel_wing.panels)
                    for vortex_panel_index, vortex_panel in np.ndenumerate(vortex_panel_wings_panels):
                        normalized_induced_velocity_at_collocation_point = (
                            vortex_panel.calculate_normalized_induced_velocity(collocation_panel.collocation_point))
                        normalized_induced_downwash_at_collocation_point = (
                            vortex_panel.calculate_normalized_induced_downwash(collocation_panel.collocation_point))

                        normal_direction_at_collocation_point = collocation_panel.normal_direction_at_collocation_point

                        normal_normalized_induced_velocity_at_collocation_point = np.dot(
                            normalized_induced_velocity_at_collocation_point, normal_direction_at_collocation_point)
                        normal_normalized_induced_downwash_at_collocation_point = np.dot(
                            normalized_induced_downwash_at_collocation_point, normal_direction_at_collocation_point)

                        self.aerodynamic_influence_coefficients[collocation_panel_index, vortex_panel_index] = (
                            normal_normalized_induced_velocity_at_collocation_point)
                        self.downwash_influence_coefficients[collocation_panel_index, vortex_panel_index] = (
                            normal_normalized_induced_downwash_at_collocation_point)

    # ToDo: Properly cite and document this method.
    def set_up_operating_point(self):
        # This calculates and updates the direction the wind is going to, in geometry axes coordinates.
        self.freestream_velocity = np.expand_dims(
            self.operating_point.calculate_freestream_velocity_geometry_axes(),
            0)

        for collocation_panel_wing in self.airplane.wings:

            collocation_panel_wings_panels = np.ravel(collocation_panel_wing.panels)
            for collocation_panel_index, collocation_panel in np.ndenumerate(collocation_panel_wings_panels):
                self.normal_directions[collocation_panel_index] = (
                    collocation_panel.normal_direction_at_collocation_point)
                self.freestream_influences[collocation_panel_index] = (
                    np.dot(self.freestream_velocity, collocation_panel.normal_direction_at_collocation_point))

    # ToDo: Properly cite and document this method.
    def calculate_vortex_strengths(self):
        """

        :return:
        """
        # # Calculate Vortex Strengths
        # ----------------------------
        # Governing Equation: AIC @ Gamma + freestream_influence = 0

        self.vortex_strengths = np.linalg.solve(self.aerodynamic_influence_coefficients, -self.freestream_influences)
        self.downwashes = self.downwash_influence_coefficients @ self.vortex_strengths
        for wing in self.airplane.wings:

            wing_panels = np.ravel(wing.panels)
            for panel_index, panel in np.ndenumerate(wing_panels):
                panel.ring_vortex.update_strength(self.vortex_strengths[panel_index])
                panel.downwash = self.downwashes[panel_index]
                if panel.horseshoe_vortex is not None:
                    panel.horseshoe_vortex.update_strength(self.vortex_strengths[panel_index])

    # ToDo: Properly cite and document this method.
    def calculate_solution_velocity(self, point):
        velocity_induced_by_vortices = np.zeros(3)

        for wing in self.airplane.wings:

            wing_panels = np.ravel(wing.panels)
            for panel in wing_panels:
                velocity_induced_by_vortices += panel.calculate_induced_velocity(point)

        return velocity_induced_by_vortices

    # ToDo: Properly cite and document this method.
    def calculate_near_field_forces(self):

        airplane = self.airplane
        density = self.operating_point.density
        freestream_velocity = self.freestream_velocity

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

                        # if panel.is_leading_edge:
                        #     bound_vortex.far_field_induced_drag = (density
                        #                                            * panel.downwash
                        #                                            * bound_vortex.strength
                        #                                            * panel.width)
                        # else:
                        #     previous_chordwise_location = chordwise_location - 1
                        #     previous_panel = panels[previous_chordwise_location, spanwise_location]
                        #     previous_panel_ring_vortex = previous_panel.ring_vortex
                        #     previous_panel_vortex_strength = previous_panel_ring_vortex.strength
                        #     bound_vortex.induced_drag = (density
                        #                                  * panel.downwash
                        #                                  * (bound_vortex.strength - previous_panel_vortex_strength)
                        #                                  * panel.width)
                    panel.update_force_moment_and_pressure()
