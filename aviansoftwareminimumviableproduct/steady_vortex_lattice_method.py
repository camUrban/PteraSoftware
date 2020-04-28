# ToDo: Properly document this module.
"""This module contains useful aerodynamics functions.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np
import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly document this class.
class SteadyVortexLatticeMethod(asmvp.problems.SteadyProblem):
    """

    """

    # ToDo: Properly document this method.
    def __init__(self, airplane, operating_point):
        """

        :param airplane:
        :param operating_point:
        """

        super().__init__(airplane, operating_point)
        self.verbose = True

        self.front_left_vertices = None
        self.front_right_vertices = None
        self.back_left_vertices = None
        self.back_right_vertices = None

        self.front_left_vortex_vertices = None
        self.front_right_vortex_vertices = None
        self.back_left_vortex_vertices = None
        self.back_right_vortex_vertices = None

        self.front_vortex_legs = None
        self.right_vortex_legs = None
        self.back_vortex_legs = None
        self.left_vortex_legs = None

        self.front_vortex_leg_centers = None
        self.back_vortex_leg_centers = None
        self.left_vortex_leg_centers = None
        self.right_vortex_leg_centers = None

        self.areas = None
        self.is_trailing_edge = None
        self.is_leading_edge = None
        self.collocation_points = None
        self.normal_directions = None
        self.n_panels = None

        self.initial_front_left_vertices = None
        self.initial_front_right_vertices = None
        self.initial_back_left_vertices = None
        self.initial_back_right_vertices = None

        self.velocity_influences_at_collocations = None
        self.aerodynamic_influence_coefficients = None

        self.vortex_strengths = None

        self.steady_freestream_velocity = None
        self.rotation_freestream_velocities = None
        self.freestream_velocity = None
        self.freestream_influences = None

        self.panel_centers = None

    # ToDo: Properly document this method.
    def run(self):
        """

        :return:
        """
        asmvp.meshing.mesh_problem(self)
        self.setup_geometry()
        self.setup_operating_point()
        self.calculate_vortex_strengths()
        self.calculate_near_field_forces_and_moments()
        self.calculate_delta_cp()
        # asmvp.output_tools.draw(self, "collocation_points")

    # ToDo: Properly cite and document this method.
    def setup_geometry(self):
        """

        :return:
        """
        # Calculate AIC matrix
        if self.verbose:
            print("Calculating the collocation influence matrices...")

        airplane_num_panels = 0
        for wing in self.airplane.wings:
            airplane_num_panels = airplane_num_panels + wing.num_panels

        self.aerodynamic_influence_coefficients = np.empty(())

        this_wings_collocation_points_aerodynamics_influence_coefficients = np.empty((0, airplane_num_panels))

        for wing_with_collocation_points in self.airplane.wings:
            reshaped_collocation_points = np.reshape(wing_with_collocation_points.collocation_points, (-1, 3))
            reshaped_normal_directions = np.reshape(wing_with_collocation_points.normal_directions, (-1, 3))
            wings_collocation_points_aerodynamics_influence_coefficients_from_wings_vortices = np.empty(
                (self.airplane.wings[0].num_panels, 0))
            for wing_with_vortices in self.airplane.wings:
                reshaped_front_left_vortex_vertices = np.reshape(wing_with_vortices.front_left_vortex_vertices, (-1, 3))
                reshaped_front_right_vortex_vertices = np.reshape(wing_with_vortices.front_right_vortex_vertices,
                                                                  (-1, 3))
                reshaped_back_left_vortex_vertices = np.reshape(wing_with_vortices.back_left_vortex_vertices, (-1, 3))
                reshaped_back_right_vortex_vertices = np.reshape(wing_with_vortices.back_right_vortex_vertices, (-1, 3))
                reshaped_is_trailing_edge = np.reshape(wing_with_vortices.is_trailing_edge, (-1))
                wings_collocation_points_aerodynamics_influence_coefficients_from_this_wings_vortices = np.zeros(
                    (wing_with_collocation_points.num_panels, wing_with_vortices.num_panels))
                for collocation_point_num in range(wing_with_collocation_points.num_panels):
                    collocation_point = reshaped_collocation_points[collocation_point_num]
                    for vortex_num in range(wing_with_vortices.num_panels):
                        front_left_vortex_vertex = reshaped_front_left_vortex_vertices[vortex_num]
                        front_right_vortex_vertex = reshaped_front_right_vortex_vertices[vortex_num]
                        back_left_vortex_vertex = reshaped_back_left_vortex_vertices[vortex_num]
                        back_right_vortex_vertex = reshaped_back_right_vortex_vertices[vortex_num]
                        velocity_induced_by_ring_vortex = (asmvp.aerodynamics.induced_velocity_from_ring_vortex
                                                           (point=collocation_point,
                                                            front_left_vortex_point=front_left_vortex_vertex,
                                                            front_right_vortex_point=front_right_vortex_vertex,
                                                            back_left_vortex_point=back_left_vortex_vertex,
                                                            back_right_vortex_point=back_right_vortex_vertex,
                                                            vortex_strength=1))
                        if reshaped_is_trailing_edge[vortex_num]:
                            velocity_induced_by_horseshoe_vortex = (asmvp.aerodynamics.induced_velocity_from_horseshoe_vortex
                                                                    (point=collocation_point,
                                                                     finite_leg_origin=back_right_vortex_vertex,
                                                                     finite_leg_termination=back_left_vortex_vertex,
                                                                     vortex_strength=1))
                        else:
                            velocity_induced_by_horseshoe_vortex = [0, 0, 0]

                        velocity_induced = velocity_induced_by_ring_vortex + velocity_induced_by_horseshoe_vortex
                        normal_velocity_induced = np.dot(velocity_induced,
                                                         reshaped_normal_directions[collocation_point_num])
                        wings_collocation_points_aerodynamics_influence_coefficients_from_this_wings_vortices[
                            collocation_point_num, vortex_num] = normal_velocity_induced

                wings_collocation_points_aerodynamics_influence_coefficients_from_wings_vortices = np.hstack(
                    (wings_collocation_points_aerodynamics_influence_coefficients_from_wings_vortices,
                     wings_collocation_points_aerodynamics_influence_coefficients_from_this_wings_vortices))

            this_wings_collocation_points_aerodynamics_influence_coefficients = np.vstack(
                (this_wings_collocation_points_aerodynamics_influence_coefficients,
                 wings_collocation_points_aerodynamics_influence_coefficients_from_wings_vortices))

        self.aerodynamic_influence_coefficients = this_wings_collocation_points_aerodynamics_influence_coefficients

        if self.verbose:
            print("Collocation influence matrix calculated!")

    # ToDo: Properly cite and document this method.
    # ToDo: Check this method.
    # ToDo: Make this method actually compute the rotational velocity.
    def setup_operating_point(self):
        """

        :return:
        """

        if self.verbose:
            print("Calculating the freestream influence vector...")

        # This calculates and updates the direction the wind is going to, in geometry axes coordinates.
        self.steady_freestream_velocity = np.expand_dims(
            self.operating_point.compute_freestream_velocity_geometry_axes(),
            0)

        # This represents the freestream velocity at each panel collocation point. It is size N x 3 where N is the
        # number of collocation points.
        self.freestream_velocity = self.steady_freestream_velocity

        self.normal_directions = np.empty((0, 3))

        for wing in self.airplane.wings:
            reshaped_normal_directions = np.reshape(wing.normal_directions, (-1, 3))
            self.normal_directions = np.vstack((self.normal_directions, reshaped_normal_directions))

        self.freestream_influences = np.swapaxes(np.dot(self.freestream_velocity, np.swapaxes(self.normal_directions, 0, 1)), 0, 1)

        if self.verbose:
            print("Freestream influence vector calculated!")

    # ToDo: Properly cite and document this method.
    def calculate_vortex_strengths(self):
        """

        :return:
        """
        # # Calculate Vortex Strengths
        # ----------------------------
        # Governing Equation: AIC @ Gamma + freestream_influence = 0
        if self.verbose:
            print("Calculating vortex strengths...")
        self.vortex_strengths = np.linalg.solve(self.aerodynamic_influence_coefficients, -self.freestream_influences)

        vortex_strengths_to_be_sliced = self.vortex_strengths

        for wing in self.airplane.wings:
            wings_vortex_strengths_listed = vortex_strengths_to_be_sliced[:wing.num_panels, 0]
            wing.vortex_strengths = np.reshape(wings_vortex_strengths_listed, (wing.num_chordwise_panels,
                                                                               wing.num_spanwise_panels))
            vortex_strengths_to_be_sliced = vortex_strengths_to_be_sliced[wing.num_panels:, 0]

    # ToDo: Properly cite and document this method.
    def solution_velocity_from_vortices(self, point):
        """

        :param point:
        :return velocity_induced_by_vortices:
        """
        velocity_induced_by_vortices = [0, 0, 0]

        for wing in self.airplane.wings:
            for chordwise_panel_num in range(wing.num_chordwise_panels):
                for spanwise_panel_num in range(wing.num_spanwise_panels):

                    front_left_vortex_vertex = wing.front_left_vortex_vertices[chordwise_panel_num, spanwise_panel_num]
                    front_right_vortex_vertex = wing.front_right_vortex_vertices[chordwise_panel_num, spanwise_panel_num]
                    back_left_vortex_vertex = wing.back_left_vortex_vertices[chordwise_panel_num, spanwise_panel_num]
                    back_right_vortex_vertex = wing.back_right_vortex_vertices[chordwise_panel_num, spanwise_panel_num]
                    vortex_strength = wing.vortex_strengths[chordwise_panel_num, spanwise_panel_num]
                    velocity_induced_by_ring_vortex = (asmvp.aerodynamics.induced_velocity_from_ring_vortex
                                                       (point=point,
                                                        front_left_vortex_point=front_left_vortex_vertex,
                                                        front_right_vortex_point=front_right_vortex_vertex,
                                                        back_left_vortex_point=back_left_vortex_vertex,
                                                        back_right_vortex_point=back_right_vortex_vertex,
                                                        vortex_strength=vortex_strength))
                    if wing.is_trailing_edge[chordwise_panel_num, spanwise_panel_num]:
                        velocity_induced_by_horseshoe_vortex = (asmvp.aerodynamics.induced_velocity_from_horseshoe_vortex
                                                                (point=point,
                                                                 finite_leg_origin=back_right_vortex_vertex,
                                                                 finite_leg_termination=back_left_vortex_vertex,
                                                                 vortex_strength=vortex_strength))
                    else:
                        velocity_induced_by_horseshoe_vortex = [0, 0, 0]

                    velocity_induced_by_this_panel = velocity_induced_by_ring_vortex + velocity_induced_by_horseshoe_vortex
                    velocity_induced_by_vortices = velocity_induced_by_vortices + velocity_induced_by_this_panel

        return velocity_induced_by_vortices

    # ToDo: Properly cite and document this method.
    def calculate_velocity_influences(self, points):
        """

        :param points:
        :return:
        """
        # Calculates the  velocity influence matrix. The first index is collocation point number, second index is vortex
        # number.
        # Points is the ndarray of points, of shape N x 3, to calculate the velocity influence at.

        velocity_influences_ring_vortices = self.calculate_velocity_influences_without_wake(
            points)  # Calculates the section of the velocity influence matrix corresponding to the ring vortices.
        velocity_influences_horseshoe_vortices = self.calculate_velocity_influences_with_wake(
            points)  # Calculates the section of velocity_influences corresponding to the trailing horseshoe vortices
        # velocity_influences = np.hstack((velocity_influences_ring_vortices, velocity_influences_horseshoe_vortices))

        velocity_influences = np.zeros((points.shape[0], self.n_panels, 3))
        mask = np.tile(np.expand_dims(np.expand_dims(self.is_trailing_edge, 0), 2), (points.shape[0], 1, 3))

        np.place(
            velocity_influences,
            np.logical_not(mask),
            velocity_influences_ring_vortices
        )
        np.place(
            velocity_influences,
            mask,
            velocity_influences_horseshoe_vortices
        )

        return velocity_influences

    # ToDo: Properly cite and document this method.
    def calculate_delta_cp(self):
        """

        :return:
        """
        for wing in self.airplane.wings:
            for chordwise_panel_num in range(wing.num_chordwise_panels):
                for spanwise_panel_num in range(wing.num_spanwise_panels):
                    area = wing.areas[chordwise_panel_num, spanwise_panel_num]
                    total_force_on_panel = wing.total_force_on_panel_in_geometry_axes[chordwise_panel_num, spanwise_panel_num]
                    normal_direction = wing.normal_directions[chordwise_panel_num, spanwise_panel_num]
                    normal_force_on_panel = np.dot(total_force_on_panel, normal_direction)
                    pressure_on_panel = normal_force_on_panel / area
                    panel_delta_pressure_coefficient = pressure_on_panel / self.operating_point.dynamic_pressure()

                    wing.normal_force_on_panels[chordwise_panel_num, spanwise_panel_num] = normal_force_on_panel
                    wing.pressure_on_panels[chordwise_panel_num, spanwise_panel_num] = pressure_on_panel
                    wing.panels_delta_pressure_coefficient[chordwise_panel_num, spanwise_panel_num] = panel_delta_pressure_coefficient

    # ToDo: Properly cite and document this method.
    def calculate_near_field_forces_and_moments(self):
        """

        :return:
        """
        if self.verbose:
            print("Calculating forces on panels...")

        density = self.operating_point.density
        self.total_force_on_airplane_in_geometry_axes = np.array([0, 0, 0])
        self.total_moment_on_airplane_in_geometry_axes = np.array([0, 0, 0])
        for wing in self.airplane.wings:
            for chordwise_panel_num in range(wing.num_chordwise_panels):
                for spanwise_panel_num in range(wing.num_spanwise_panels):
                    front_vortex_leg_center = wing.front_vortex_leg_centers[chordwise_panel_num, spanwise_panel_num]
                    back_vortex_leg_center = wing.back_vortex_leg_centers[chordwise_panel_num, spanwise_panel_num]
                    left_vortex_leg_center = wing.left_vortex_leg_centers[chordwise_panel_num, spanwise_panel_num]
                    right_vortex_leg_center = wing.right_vortex_leg_centers[chordwise_panel_num, spanwise_panel_num]

                    front_vortex_leg = wing.front_vortex_legs[chordwise_panel_num, spanwise_panel_num]
                    back_vortex_leg = wing.back_vortex_legs[chordwise_panel_num, spanwise_panel_num]
                    left_vortex_leg = wing.left_vortex_legs[chordwise_panel_num, spanwise_panel_num]
                    right_vortex_leg = wing.right_vortex_legs[chordwise_panel_num, spanwise_panel_num]

                    panel_center = wing.panel_centers[chordwise_panel_num, spanwise_panel_num]

                    velocity_at_front_vortex_leg_center = self.solution_velocity_from_vortices(front_vortex_leg_center)
                    velocity_at_back_vortex_leg_center = self.solution_velocity_from_vortices(back_vortex_leg_center)
                    velocity_at_left_vortex_leg_center = self.solution_velocity_from_vortices(left_vortex_leg_center)
                    velocity_at_right_vortex_leg_center = self.solution_velocity_from_vortices(right_vortex_leg_center)

                    velocity_at_front_vortex_centers_cross_front_vortex_leg = np.cross(
                        velocity_at_front_vortex_leg_center, front_vortex_leg)
                    velocity_at_back_vortex_centers_cross_back_vortex_leg = np.cross(
                        velocity_at_back_vortex_leg_center, back_vortex_leg)
                    velocity_at_left_vortex_centers_cross_left_vortex_leg = np.cross(
                        velocity_at_left_vortex_leg_center, left_vortex_leg)
                    velocity_at_right_vortex_centers_cross_right_vortex_leg = np.cross(
                        velocity_at_right_vortex_leg_center, right_vortex_leg)

                    vortex_strength_expanded = wing.vortex_strengths[chordwise_panel_num, spanwise_panel_num]
                    # vortex_strength_expanded = np.expand_dims(self.vortex_strengths, axis=1)

                    wing.forces_on_front_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num] = (
                            density * velocity_at_front_vortex_centers_cross_front_vortex_leg
                            * vortex_strength_expanded)
                    wing.forces_on_back_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num] = (
                            density * velocity_at_back_vortex_centers_cross_back_vortex_leg
                            * vortex_strength_expanded)
                    wing.forces_on_left_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num] = (
                            density * velocity_at_left_vortex_centers_cross_left_vortex_leg
                            * vortex_strength_expanded)
                    wing.forces_on_right_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num] = (
                            density * velocity_at_right_vortex_centers_cross_right_vortex_leg
                            * vortex_strength_expanded)

                    wing.total_force_on_panel_in_geometry_axes[chordwise_panel_num, spanwise_panel_num] = (
                        wing.forces_on_front_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num]
                        + wing.forces_on_back_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num]
                        + wing.forces_on_left_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num]
                        + wing.forces_on_right_vortices_in_geometry_axes[chordwise_panel_num, spanwise_panel_num])

                    wing.total_moment_on_panel_in_geometry_axes[chordwise_panel_num, spanwise_panel_num] = (
                        np.cross(panel_center - self.airplane.xyz_ref,
                                 wing.total_force_on_panel_in_geometry_axes[chordwise_panel_num, spanwise_panel_num])
                    )

                    self.total_force_on_airplane_in_geometry_axes = (self.total_force_on_airplane_in_geometry_axes
                                                                     + wing.total_force_on_panel_in_geometry_axes
                                                                     [chordwise_panel_num, spanwise_panel_num])
                    self.total_moment_on_airplane_in_geometry_axes = (self.total_moment_on_airplane_in_geometry_axes
                                                                      + wing.total_moment_on_panel_in_geometry_axes
                                                                      [chordwise_panel_num, spanwise_panel_num])

        self.total_force_on_airplane_in_wind_axes = (np.transpose(
            self.operating_point.compute_rotation_matrix_wind_to_geometry())
                                                     @ self.total_force_on_airplane_in_geometry_axes)
        self.total_moment_on_airplane_in_wing_axes = (np.transpose(
            self.operating_point.compute_rotation_matrix_wind_to_geometry())
                                                     @ self.total_moment_on_airplane_in_geometry_axes)

        # Calculate nondimensional forces
        q = self.operating_point.dynamic_pressure()
        s_ref = self.airplane.s_ref
        b_ref = self.airplane.b_ref
        c_ref = self.airplane.c_ref

        self.airplane_coefficient_of_lift = -self.total_force_on_airplane_in_wind_axes[2] / q / s_ref
        self.airplane_coefficient_of_induced_drag = -self.total_force_on_airplane_in_wind_axes[0] / q / s_ref
        self.airplane_coefficient_of_side_force = self.total_force_on_airplane_in_wind_axes[1] / q / s_ref
        self.airplane_coefficient_of_rolling_moment = self.total_moment_on_airplane_in_wing_axes[0] / q / b_ref
        self.airplane_coefficient_of_pitching_moment = self.total_moment_on_airplane_in_wing_axes[1] / q / c_ref
        self.airplane_coefficient_of_yawing_moment = self.total_moment_on_airplane_in_wing_axes[2] / q / b_ref

        # Solves divide by zero error
        if self.airplane_coefficient_of_induced_drag == 0:
            self.airplane_coefficient_of_lift_over_coefficient_of_induced_drag = 0
        else:
            self.airplane_coefficient_of_lift_over_coefficient_of_induced_drag = (self.airplane_coefficient_of_lift / self.airplane_coefficient_of_induced_drag)

        if self.verbose:
            print("\nForces\n-----")
            print("airplane_coefficient_of_lift: ", self.airplane_coefficient_of_lift)
            print("airplane_coefficient_of_induced_drag: ", self.airplane_coefficient_of_induced_drag)
            print("airplane_coefficient_of_side_force: ", self.airplane_coefficient_of_side_force)
            print("CL/CDi: ", self.airplane_coefficient_of_lift_over_coefficient_of_induced_drag)
            print("\nMoments\n-----")
            print("Cl: ", self.airplane_coefficient_of_rolling_moment)
            print("Cm: ", self.airplane_coefficient_of_pitching_moment)
            print("Cn: ", self.airplane_coefficient_of_yawing_moment)

            print("Finished calculating forces!")

    # ToDo: Properly cite and document this method.
    def get_velocity_at_point(self, point):
        """

        :param point:
        :return:
        """
        # Input: a Nx3 numpy array of points that you would like to know the velocities at.
        # Output: a Nx3 numpy array of the velocities at those points.
        point = np.reshape(point, (-1, 3))

        Vi = self.get_induced_velocity_at_point(point)

        freestream = self.operating_point.compute_freestream_velocity_geometry_axes()

        V = Vi + freestream

        return V
