# ToDo: Properly document this module.
"""

"""
import numpy as np

import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly document this class.
class SteadyVortexLatticeMethod(asmvp.aerodynamics_problems.SteadyAerodynamicsProblem):
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
        asmvp.mesh_tools.mesh_problem(self)
        self.setup_geometry()
        self.setup_operating_point()
        self.calculate_vortex_strengths()
        # self.calculate_near_field_forces()
        # self.calculate_delta_cp()
        asmvp.output_tools.draw(self, "collocation_points")

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
                        velocity_induced_by_ring_vortex = (self.velocity_induced_by_ring_vortex
                                                           (point=collocation_point,
                                                            front_left_vortex_point=front_left_vortex_vertex,
                                                            front_right_vortex_point=front_right_vortex_vertex,
                                                            back_left_vortex_point=back_left_vortex_vertex,
                                                            back_right_vortex_point=back_right_vortex_vertex,
                                                            vortex_strength=1))
                        if reshaped_is_trailing_edge[vortex_num]:
                            velocity_induced_by_horseshoe_vortex = (self.velocity_induced_by_horseshoe_vortex
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

    def velocity_induced_by_vortices(self, point):

        velocity_induced_by_vortices = [0, 0, 0]

        for wing in self.airplane.wings:
            num_panels = wing.num_panels
            for panel_num in range(num_panels):
                front_left_vortex_vertex = wing.front_left_vortex_vertices[panel_num]
                front_right_vortex_vertex = wing.front_right_vortex_vertices[panel_num]
                back_left_vortex_vertex = wing.back_left_vortex_vertices[panel_num]
                back_right_vortex_vertex = wing.back_right_vortex_vertices[panel_num]
                velocity_induced_by_ring_vortex = (self.velocity_induced_by_ring_vortex
                                                   (point=point,
                                                    front_left_vortex_point=front_left_vortex_vertex,
                                                    front_right_vortex_point=front_right_vortex_vertex,
                                                    back_left_vortex_point=back_left_vortex_vertex,
                                                    back_right_vortex_point=back_right_vortex_vertex,
                                                    vortex_strength=1))
                if wing.is_trailing_edge[panel_num]:
                    velocity_induced_by_horseshoe_vortex = (self.velocity_induced_by_horseshoe_vortex
                                                            (point=point,
                                                             finite_leg_origin=back_right_vortex_vertex,
                                                             finite_leg_termination=back_left_vortex_vertex,
                                                             vortex_strength=1))
                else:
                    velocity_induced_by_horseshoe_vortex = [0, 0, 0]

                velocity_induced_by_this_panel = velocity_induced_by_ring_vortex + velocity_induced_by_horseshoe_vortex
                velocity_induced_by_vortices = velocity_induced_by_vortices + velocity_induced_by_this_panel

        return velocity_induced_by_vortices

    # ToDo: Properly cite and document this method.
    def velocity_induced_by_line_vortex(self, point, vortex_origin, vortex_termination,
                                        vortex_strength):
        """

        :param point:
        :param vortex_origin:
        :param vortex_termination:
        :param vortex_strength:
        :return: velocity_induced_by_vortex:
        """
        r_1 = point - vortex_origin
        r_2 = point - vortex_termination
        r_0 = r_1 - r_2

        r_1_cross_r_2 = np.cross(r_1, r_2)
        r_1_cross_r_2_absolute_value = r_1_cross_r_2[0] ** 2 + r_1_cross_r_2[1] ** 2 + r_1_cross_r_2[2] ** 2

        r_1_length = np.linalg.norm(r_1)
        r_2_length = np.linalg.norm(r_2)

        line_vortex_radius = 3.0e-16
        if (r_1_length < line_vortex_radius
                or r_2_length < line_vortex_radius
                or r_1_cross_r_2_absolute_value < line_vortex_radius):
            velocity_induced_by_vortex = [0, 0, 0]
        else:
            r_0_dot_r_1 = np.dot(r_0, r_1)
            r_0_dot_r_2 = np.dot(r_0, r_2)

            k = (vortex_strength / (4 * np.pi * r_1_cross_r_2_absolute_value)
                 * (r_0_dot_r_1 / r_1_length - r_0_dot_r_2 / r_2_length))

            u = k * r_1_cross_r_2[0]
            v = k * r_1_cross_r_2[1]
            w = k * r_1_cross_r_2[2]
            velocity_induced_by_vortex = np.array([u, v, w])

        return velocity_induced_by_vortex

    # ToDo: Properly cite and document this method.
    def velocity_induced_by_horseshoe_vortex(self, point, finite_leg_origin, finite_leg_termination, vortex_strength):
        """

        :param vortex_strength:
        :param point:
        :param finite_leg_origin:
        :param finite_leg_termination:
        :param vortex_strength:
        :return: velocity_induced_by_vortex:
        """
        a = point - finite_leg_origin
        b = point - finite_leg_termination

        a_length = np.linalg.norm(a)
        b_length = np.linalg.norm(b)

        point_on_first_infinite_leg = finite_leg_origin + [0, 0, 1]
        distance_along_first_infinite_leg = (np.linalg.norm(np.cross((point - finite_leg_origin),
                                                                     (point - point_on_first_infinite_leg)))
                                             / np.linalg.norm(point_on_first_infinite_leg - finite_leg_origin))

        point_on_second_infinite_leg = finite_leg_termination + [0, 0, 1]
        distance_along_second_infinite_leg = (np.linalg.norm(np.cross((point - finite_leg_termination),
                                                                      (point - point_on_second_infinite_leg)))
                                              / np.linalg.norm(point_on_second_infinite_leg
                                                               - finite_leg_termination))

        horseshoe_vortex_radius = 3.0e-16
        if (a_length < horseshoe_vortex_radius or b_length < horseshoe_vortex_radius
                or distance_along_first_infinite_leg < horseshoe_vortex_radius
                or distance_along_second_infinite_leg < horseshoe_vortex_radius):
            velocity_induced_by_vortex = [0, 0, 0]
        else:
            a_cross_b = np.cross(a, b)
            b_cross_a = np.cross(b, a)

            a_dot_b = np.dot(a, b)

            x_unit = [1, 0, 0]

            a_cross_x_unit = np.cross(a, x_unit)
            b_cross_x_unit = np.cross(b, x_unit)

            a_dot_x_unit = np.dot(a, x_unit)
            b_dot_x_unit = np.dot(b, x_unit)

            velocity_induced_by_vortex = (vortex_strength / (4 * np.pi)
                                          * ((a_cross_b / (a_length * b_length + a_dot_b))
                                             * (1 / a_length + 1 / b_length)
                                             + (a_cross_x_unit / (a_length - a_dot_x_unit))
                                             * (1 / a_length)
                                             - (b_cross_x_unit / (b_length - b_dot_x_unit))
                                             * (1 / b_length)))

        return velocity_induced_by_vortex

    # ToDo: Properly cite and document this method.
    def velocity_induced_by_ring_vortex(self, point, front_left_vortex_point, front_right_vortex_point,
                                        back_left_vortex_point, back_right_vortex_point, vortex_strength):
        """

        :param point:
        :param front_left_vortex_point:
        :param front_right_vortex_point:
        :param back_left_vortex_point:
        :param back_right_vortex_point:
        :param vortex_strength:
        :return: velocity_induced_by_vortex:
        """
        velocity_induced_by_upper_leg = self.velocity_induced_by_line_vortex(
            point=point, vortex_origin=front_right_vortex_point, vortex_termination=front_left_vortex_point,
            vortex_strength=vortex_strength)
        velocity_induced_by_left_leg = self.velocity_induced_by_line_vortex(
            point=point, vortex_origin=front_left_vortex_point, vortex_termination=back_left_vortex_point,
            vortex_strength=vortex_strength)
        velocity_induced_by_lower_leg = self.velocity_induced_by_line_vortex(
            point=point, vortex_origin=back_left_vortex_point, vortex_termination=back_right_vortex_point,
            vortex_strength=vortex_strength)
        velocity_induced_by_right_leg = self.velocity_induced_by_line_vortex(
            point=point, vortex_origin=back_right_vortex_point, vortex_termination=front_right_vortex_point,
            vortex_strength=vortex_strength)

        velocity_induced_by_vortex = (velocity_induced_by_upper_leg
                                      + velocity_induced_by_left_leg
                                      + velocity_induced_by_lower_leg
                                      + velocity_induced_by_right_leg)

        return velocity_induced_by_vortex

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
    def calculate_velocity_influences_without_wake(self, points):
        """

        :param points:
        :return:
        """
        # Calculates the doublet part of Vij, the velocity influence matrix (First index is collocation point number,
        # second index is vortex number).
        # points: the list of points (Nx3) to calculate the velocity influence at.

        # Data cleanup
        points = np.reshape(points, (-1, 3))

        # Make v1, v2, v3, and v4 vectors.
        # Each vector goes from all collocation points to one type of vertex (front left, front right, etc.). NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # v1: corresponds to front left vertices
        # v2: corresponds to front right vertices
        # v3: corresponds to back right vertices
        # v4: corresponds to back left vertices
        # Example: v1[i,j,:] = collocation_points[i,:] - front_left_vertices[j,:]

        points = np.expand_dims(points, 1)
        v1 = points - self.front_left_vortex_vertices
        v2 = points - self.front_right_vortex_vertices
        v3 = points - self.back_right_vortex_vertices
        v4 = points - self.back_left_vortex_vertices

        # Do some useful arithmetic
        v1_cross_v2 = np.cross(v1, v2, axis=2)
        v2_cross_v3 = np.cross(v2, v3, axis=2)
        v3_cross_v4 = np.cross(v3, v4, axis=2)
        v4_cross_v1 = np.cross(v4, v1, axis=2)
        v1_dot_v2 = np.einsum('ijk,ijk->ij', v1, v2)
        v2_dot_v3 = np.einsum('ijk,ijk->ij', v2, v3)
        v3_dot_v4 = np.einsum('ijk,ijk->ij', v3, v4)
        v4_dot_v1 = np.einsum('ijk,ijk->ij', v4, v1)
        norm_v1 = np.linalg.norm(v1, axis=2)
        norm_v2 = np.linalg.norm(v2, axis=2)
        norm_v3 = np.linalg.norm(v3, axis=2)
        norm_v4 = np.linalg.norm(v4, axis=2)
        norm_v1_inv = 1 / norm_v1
        norm_v2_inv = 1 / norm_v2
        norm_v3_inv = 1 / norm_v3
        norm_v4_inv = 1 / norm_v4

        # Check for the special case where the collocation point is along the bound vortex leg
        # Find where cross product is near zero, and set the dot product to infinity so that the value of the bound term
        # is zero.
        v1_v2_singularity_indices = (
                np.einsum('ijk,ijk->ij', v1_cross_v2, v1_cross_v2)  # norm(cross_product)^2
                < 3.0e-16)
        v1_dot_v2 = v1_dot_v2 + v1_v2_singularity_indices

        v2_v3_singularity_indices = (
                np.einsum('ijk,ijk->ij', v2_cross_v3, v2_cross_v3)  # norm(cross_product)^2
                < 3.0e-16)
        v2_dot_v3 = v2_dot_v3 + v2_v3_singularity_indices

        v3_v4_singularity_indices = (
                np.einsum('ijk,ijk->ij', v3_cross_v4, v3_cross_v4)  # norm(cross_product)^2
                < 3.0e-16)
        v3_dot_v4 = v3_dot_v4 + v3_v4_singularity_indices

        v4_v1_singularity_indices = (
                np.einsum('ijk,ijk->ij', v4_cross_v1, v4_cross_v1)  # norm(cross_product)^2
                < 3.0e-16)
        v4_dot_v1 = v4_dot_v1 + v4_v1_singularity_indices

        # Calculate Vij
        term1 = (norm_v1_inv + norm_v2_inv) / (norm_v1 * norm_v2 + v1_dot_v2)
        term1 = np.expand_dims(term1, 2)
        term2 = (norm_v2_inv + norm_v3_inv) / (norm_v2 * norm_v3 + v2_dot_v3)
        term2 = np.expand_dims(term2, 2)
        term3 = (norm_v3_inv + norm_v4_inv) / (norm_v3 * norm_v4 + v3_dot_v4)
        term3 = np.expand_dims(term3, 2)
        term4 = (norm_v4_inv + norm_v1_inv) / (norm_v4 * norm_v1 + v4_dot_v1)
        term4 = np.expand_dims(term4, 2)

        velocity_influences_ring_vortices = 1 / (4 * np.pi) * (
                v1_cross_v2 * term1 +
                v2_cross_v3 * term2 +
                v3_cross_v4 * term3 +
                v4_cross_v1 * term4
        )

        return velocity_influences_ring_vortices

    # ToDo: Properly cite and document this method.
    def calculate_velocity_influences_with_wake(self, points):
        """

        :param points:
        :return:
        """
        # Calculates the doublet part of Vij, the velocity influence matrix (First index is collocation point number,
        # second index is vortex number).
        # points: the list of points (Nx3) to calculate the velocity influence at.

        # Data cleanup
        points = np.reshape(points, (-1, 3))

        # Make v1, v2, v3, and v4 vectors.
        # Each vector goes from all collocation points to one type of vertex (front left, front right, etc.). NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # v1: corresponds to front left vertices
        # v2: corresponds to front right vertices
        # v3: corresponds to back right vertices
        # v4: corresponds to back left vertices
        # Example: v1[i,j,:] = collocation_points[i,:] - front_left_vertices[j,:]

        points = np.expand_dims(points, 1)
        v1 = points - self.front_left_vortex_vertices
        v2 = points - self.front_right_vortex_vertices
        v3 = points - self.back_right_vortex_vertices
        v4 = points - self.back_left_vortex_vertices

        # Do some useful arithmetic
        v1_cross_v2 = np.cross(v1, v2, axis=2)
        v2_cross_v3 = np.cross(v2, v3, axis=2)
        v3_cross_v4 = np.cross(v3, v4, axis=2)
        v4_cross_v1 = np.cross(v4, v1, axis=2)
        v1_dot_v2 = np.einsum('ijk,ijk->ij', v1, v2)
        v2_dot_v3 = np.einsum('ijk,ijk->ij', v2, v3)
        v3_dot_v4 = np.einsum('ijk,ijk->ij', v3, v4)
        v4_dot_v1 = np.einsum('ijk,ijk->ij', v4, v1)
        norm_v1 = np.linalg.norm(v1, axis=2)
        norm_v2 = np.linalg.norm(v2, axis=2)
        norm_v3 = np.linalg.norm(v3, axis=2)
        norm_v4 = np.linalg.norm(v4, axis=2)
        norm_v1_inv = 1 / norm_v1
        norm_v2_inv = 1 / norm_v2
        norm_v3_inv = 1 / norm_v3
        norm_v4_inv = 1 / norm_v4

        # Check for the special case where the collocation point is along the bound vortex leg
        # Find where cross product is near zero, and set the dot product to infinity so that the value of the bound term
        # is zero.
        v1_v2_singularity_indices = (
                np.einsum('ijk,ijk->ij', v1_cross_v2, v1_cross_v2)  # norm(cross_product)^2
                < 3.0e-16)
        v1_dot_v2 = v1_dot_v2 + v1_v2_singularity_indices

        v2_v3_singularity_indices = (
                np.einsum('ijk,ijk->ij', v2_cross_v3, v2_cross_v3)  # norm(cross_product)^2
                < 3.0e-16)
        v2_dot_v3 = v2_dot_v3 + v2_v3_singularity_indices

        v3_v4_singularity_indices = (
                np.einsum('ijk,ijk->ij', v3_cross_v4, v3_cross_v4)  # norm(cross_product)^2
                < 3.0e-16)
        v3_dot_v4 = v3_dot_v4 + v3_v4_singularity_indices

        v4_v1_singularity_indices = (
                np.einsum('ijk,ijk->ij', v4_cross_v1, v4_cross_v1)  # norm(cross_product)^2
                < 3.0e-16)
        v4_dot_v1 = v4_dot_v1 + v4_v1_singularity_indices

        # Calculate Vij
        term1 = (norm_v1_inv + norm_v2_inv) / (norm_v1 * norm_v2 + v1_dot_v2)
        term1 = np.expand_dims(term1, 2)
        term2 = (norm_v2_inv + norm_v3_inv) / (norm_v2 * norm_v3 + v2_dot_v3)
        term2 = np.expand_dims(term2, 2)
        term3 = (norm_v3_inv + norm_v4_inv) / (norm_v3 * norm_v4 + v3_dot_v4)
        term3 = np.expand_dims(term3, 2)
        term4 = (norm_v4_inv + norm_v1_inv) / (norm_v4 * norm_v1 + v4_dot_v1)
        term4 = np.expand_dims(term4, 2)

        velocity_influences_ring_vortices = 1 / (4 * np.pi) * (
                v1_cross_v2 * term1 +
                v2_cross_v3 * term2 +
                v3_cross_v4 * term3 +
                v4_cross_v1 * term4
        )

        # Calculates Vij, the velocity influence matrix (First index is collocation point number, second index is vortex
        # number).
        # points: the list of points (Nx3) to calculate the velocity influence at.

        # Make lv and rv
        left_horseshoe_vortex_vertices = self.back_left_vortex_vertices
        right_horseshoe_vortex_vertices = self.back_right_vortex_vertices

        points = np.reshape(points, (-1, 3))
        n_points = len(points)
        n_horseshoe_vortices = len(left_horseshoe_vortex_vertices)

        # Make a and b vectors.
        # a: Vector from all collocation points to all horseshoe vortex left  vertices, NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # b: Vector from all collocation points to all horseshoe vortex right vertices, NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # a[i,j,:] = c[i,:] - lv[j,:]
        # b[i,j,:] = c[i,:] - rv[j,:]
        points = np.expand_dims(points, 1)
        a = points - left_horseshoe_vortex_vertices
        b = points - right_horseshoe_vortex_vertices
        # x_hat = np.zeros([n_points, n_horseshoe_vortices, 3])
        # x_hat[:, :, 0] = 1

        # Do some useful arithmetic
        a_cross_b = np.cross(a, b, axis=2)
        a_dot_b = np.einsum('ijk,ijk->ij', a, b)

        a_cross_x = np.stack((
            np.zeros((n_points, n_horseshoe_vortices)),
            a[:, :, 2],
            -a[:, :, 1]
        ), axis=2)
        a_dot_x = a[:, :, 0]

        b_cross_x = np.stack((
            np.zeros((n_points, n_horseshoe_vortices)),
            b[:, :, 2],
            -b[:, :, 1]
        ), axis=2)
        b_dot_x = b[:, :, 0]  # np.sum(b * x_hat,axis=2)

        norm_a = np.linalg.norm(a, axis=2)
        norm_b = np.linalg.norm(b, axis=2)
        norm_a_inv = 1 / norm_a
        norm_b_inv = 1 / norm_b

        # Check for the special case where the collocation point is along the bound vortex leg
        # Find where cross product is near zero, and set the dot product to infinity so that the value of the bound term
        # is zero.
        bound_vortex_singularity_indices = (
                np.einsum('ijk,ijk->ij', a_cross_b, a_cross_b)  # norm(a_cross_b) ** 2
                < 3.0e-16)
        a_dot_b = a_dot_b + bound_vortex_singularity_indices
        left_vortex_singularity_indices = (
                np.einsum('ijk,ijk->ij', a_cross_x, a_cross_x)
                < 3.0e-16
        )
        a_dot_x = a_dot_x + left_vortex_singularity_indices
        right_vortex_singularity_indices = (
                np.einsum('ijk,ijk->ij', b_cross_x, b_cross_x)
                < 3.0e-16
        )
        b_dot_x = b_dot_x + right_vortex_singularity_indices

        # Calculate Vij
        term1 = (norm_a_inv + norm_b_inv) / (norm_a * norm_b + a_dot_b)
        term2 = norm_a_inv / (norm_a - a_dot_x)
        term3 = norm_b_inv / (norm_b - b_dot_x)
        term1 = np.expand_dims(term1, 2)
        term2 = np.expand_dims(term2, 2)
        term3 = np.expand_dims(term3, 2)

        velocity_influences_horseshoe_vortices = 1 / (4 * np.pi) * (
                a_cross_b * term1 +
                a_cross_x * term2 -
                b_cross_x * term3
        )

        return velocity_influences_ring_vortices + velocity_influences_horseshoe_vortices

    # ToDo: Properly cite and document this method.
    def calculate_velocity_influences_ring_vortices(self, points):
        """

        :param points:
        :return:
        """
        # Calculates the doublet part of Vij, the velocity influence matrix (First index is collocation point number,
        # second index is vortex number).
        # points: the list of points (Nx3) to calculate the velocity influence at.

        # Data cleanup
        points = np.reshape(points, (-1, 3))

        # Make v1, v2, v3, and v4 vectors.
        # Each vector goes from all collocation points to one type of vertex (front left, front right, etc.). NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # v1: corresponds to front left vertices
        # v2: corresponds to front right vertices
        # v3: corresponds to back right vertices
        # v4: corresponds to back left vertices
        # Example: v1[i,j,:] = collocation_points[i,:] - front_left_vertices[j,:]

        points = np.expand_dims(points, 1)
        v1 = points - self.front_left_vortex_vertices
        v2 = points - self.front_right_vortex_vertices
        v3 = points - self.back_right_vortex_vertices
        v4 = points - self.back_left_vortex_vertices

        # Do some useful arithmetic
        v1_cross_v2 = np.cross(v1, v2, axis=2)
        v2_cross_v3 = np.cross(v2, v3, axis=2)
        v3_cross_v4 = np.cross(v3, v4, axis=2)
        v4_cross_v1 = np.cross(v4, v1, axis=2)
        v1_dot_v2 = np.einsum('ijk,ijk->ij', v1, v2)
        v2_dot_v3 = np.einsum('ijk,ijk->ij', v2, v3)
        v3_dot_v4 = np.einsum('ijk,ijk->ij', v3, v4)
        v4_dot_v1 = np.einsum('ijk,ijk->ij', v4, v1)
        norm_v1 = np.linalg.norm(v1, axis=2)
        norm_v2 = np.linalg.norm(v2, axis=2)
        norm_v3 = np.linalg.norm(v3, axis=2)
        norm_v4 = np.linalg.norm(v4, axis=2)
        norm_v1_inv = 1 / norm_v1
        norm_v2_inv = 1 / norm_v2
        norm_v3_inv = 1 / norm_v3
        norm_v4_inv = 1 / norm_v4

        # Check for the special case where the collocation point is along the bound vortex leg
        # Find where cross product is near zero, and set the dot product to infinity so that the value of the bound term
        # is zero.
        v1_v2_singularity_indices = (
                np.einsum('ijk,ijk->ij', v1_cross_v2, v1_cross_v2)  # norm(cross_product)^2
                < 3.0e-16)
        v1_dot_v2 = v1_dot_v2 + v1_v2_singularity_indices

        v2_v3_singularity_indices = (
                np.einsum('ijk,ijk->ij', v2_cross_v3, v2_cross_v3)  # norm(cross_product)^2
                < 3.0e-16)
        v2_dot_v3 = v2_dot_v3 + v2_v3_singularity_indices

        v3_v4_singularity_indices = (
                np.einsum('ijk,ijk->ij', v3_cross_v4, v3_cross_v4)  # norm(cross_product)^2
                < 3.0e-16)
        v3_dot_v4 = v3_dot_v4 + v3_v4_singularity_indices

        v4_v1_singularity_indices = (
                np.einsum('ijk,ijk->ij', v4_cross_v1, v4_cross_v1)  # norm(cross_product)^2
                < 3.0e-16)
        v4_dot_v1 = v4_dot_v1 + v4_v1_singularity_indices

        # Calculate Vij
        term1 = (norm_v1_inv + norm_v2_inv) / (norm_v1 * norm_v2 + v1_dot_v2)
        term1 = np.expand_dims(term1, 2)
        term2 = (norm_v2_inv + norm_v3_inv) / (norm_v2 * norm_v3 + v2_dot_v3)
        term2 = np.expand_dims(term2, 2)
        term3 = (norm_v3_inv + norm_v4_inv) / (norm_v3 * norm_v4 + v3_dot_v4)
        term3 = np.expand_dims(term3, 2)
        term4 = (norm_v4_inv + norm_v1_inv) / (norm_v4 * norm_v1 + v4_dot_v1)
        term4 = np.expand_dims(term4, 2)

        velocity_influences_ring_vortices = 1 / (4 * np.pi) * (
                v1_cross_v2 * term1 +
                v2_cross_v3 * term2 +
                v3_cross_v4 * term3 +
                v4_cross_v1 * term4
        )

        return velocity_influences_ring_vortices

    # ToDo: Properly cite and document this method.
    def calculate_velocity_influences_horseshoe_vortices(self, points):
        """

        :param points:
        :return:
        """
        # Calculates Vij, the velocity influence matrix (First index is collocation point number, second index is vortex
        # number).
        # points: the list of points (Nx3) to calculate the velocity influence at.

        # Make lv and rv
        left_horseshoe_vortex_vertices = self.back_left_vortex_vertices[self.is_trailing_edge]
        right_horseshoe_vortex_vertices = self.back_right_vortex_vertices[self.is_trailing_edge]

        points = np.reshape(points, (-1, 3))
        n_points = len(points)
        n_horseshoe_vortices = len(left_horseshoe_vortex_vertices)

        # Make a and b vectors.
        # a: Vector from all collocation points to all horseshoe vortex left  vertices, NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # b: Vector from all collocation points to all horseshoe vortex right vertices, NxNx3.
        #   # First index is collocation point #, second is vortex #, and third is xyz. N=num_panels
        # a[i,j,:] = c[i,:] - lv[j,:]
        # b[i,j,:] = c[i,:] - rv[j,:]
        points = np.expand_dims(points, 1)
        a = points - left_horseshoe_vortex_vertices
        b = points - right_horseshoe_vortex_vertices
        # x_hat = np.zeros([n_points, n_horseshoe_vortices, 3])
        # x_hat[:, :, 0] = 1

        # Do some useful arithmetic
        a_cross_b = np.cross(a, b, axis=2)
        a_dot_b = np.einsum('ijk,ijk->ij', a, b)

        a_cross_x = np.stack((
            np.zeros((n_points, n_horseshoe_vortices)),
            a[:, :, 2],
            -a[:, :, 1]
        ), axis=2)
        a_dot_x = a[:, :, 0]

        b_cross_x = np.stack((
            np.zeros((n_points, n_horseshoe_vortices)),
            b[:, :, 2],
            -b[:, :, 1]
        ), axis=2)
        b_dot_x = b[:, :, 0]  # np.sum(b * x_hat,axis=2)

        norm_a = np.linalg.norm(a, axis=2)
        norm_b = np.linalg.norm(b, axis=2)
        norm_a_inv = 1 / norm_a
        norm_b_inv = 1 / norm_b

        # Check for the special case where the collocation point is along the bound vortex leg
        # Find where cross product is near zero, and set the dot product to infinity so that the value of the bound term
        # is zero.
        bound_vortex_singularity_indices = (
                np.einsum('ijk,ijk->ij', a_cross_b, a_cross_b)  # norm(a_cross_b) ** 2
                < 3.0e-16)
        a_dot_b = a_dot_b + bound_vortex_singularity_indices
        left_vortex_singularity_indices = (
                np.einsum('ijk,ijk->ij', a_cross_x, a_cross_x)
                < 3.0e-16
        )
        a_dot_x = a_dot_x + left_vortex_singularity_indices
        right_vortex_singularity_indices = (
                np.einsum('ijk,ijk->ij', b_cross_x, b_cross_x)
                < 3.0e-16
        )
        b_dot_x = b_dot_x + right_vortex_singularity_indices

        # Calculate Vij
        term1 = (norm_a_inv + norm_b_inv) / (norm_a * norm_b + a_dot_b)
        term2 = norm_a_inv / (norm_a - a_dot_x)
        term3 = norm_b_inv / (norm_b - b_dot_x)
        term1 = np.expand_dims(term1, 2)
        term2 = np.expand_dims(term2, 2)
        term3 = np.expand_dims(term3, 2)

        velocity_influences_horseshoe_vortices = 1 / (4 * np.pi) * (
                a_cross_b * term1 +
                a_cross_x * term2 -
                b_cross_x * term3
        )

        return velocity_influences_horseshoe_vortices

    # ToDo: Properly cite and document this method.
    def calculate_delta_cp(self):
        """

        :return:
        """
        # Find the area of each panel
        diag1 = self.front_left_vertices - self.back_right_vertices
        diag2 = self.front_right_vertices - self.back_left_vertices
        self.areas = np.linalg.norm(np.cross(diag1, diag2, axis=1), axis=1) / 2

        # Calculate panel data
        self.Fi_normal = np.einsum('ij,ij->i', self.forces_on_panels_in_geometry_axes, self.normal_directions)
        self.pressure_normal = self.Fi_normal / self.areas
        self.delta_cp = self.pressure_normal / self.operating_point.dynamic_pressure()

    # ToDo: Properly cite and document this method.
    def calculate_near_field_forces(self):
        """

        :return:
        """
        if self.verbose:
            print("Calculating forces on panels...")

        # Calculate local velocity at the ith horseshoe vortex center point
        velocity_at_front_vortex_leg_centers = self.get_velocity_at_point(self.front_vortex_leg_centers)
        velocity_at_back_vortex_leg_centers = self.get_velocity_at_point(self.back_vortex_leg_centers)
        velocity_at_left_vortex_leg_centers = self.get_velocity_at_point(self.left_vortex_leg_centers)
        velocity_at_right_vortex_leg_centers = self.get_velocity_at_point(self.right_vortex_leg_centers)

        # Calculate Fi_geometry, the force on the ith horseshoe vortex. Note that this is in GEOMETRY AXES,
        # not WIND AXES or BODY AXES.
        density = self.operating_point.density

        velocity_at_front_vortex_centers_cross_front_vortex_leg = np.cross(velocity_at_front_vortex_leg_centers,
                                                                           self.front_vortex_legs, axis=1)
        velocity_at_back_vortex_centers_cross_back_vortex_leg = np.cross(velocity_at_back_vortex_leg_centers,
                                                                         self.back_vortex_legs, axis=1)
        velocity_at_left_vortex_centers_cross_left_vortex_leg = np.cross(velocity_at_left_vortex_leg_centers,
                                                                         self.left_vortex_legs, axis=1)
        velocity_at_right_vortex_centers_cross_right_vortex_leg = np.cross(velocity_at_right_vortex_leg_centers,
                                                                           self.right_vortex_legs, axis=1)

        vortex_strengths_expanded = np.expand_dims(self.vortex_strengths, axis=1)

        self.forces_on_front_vortices_in_geometry_axes = (density
                                                          * velocity_at_front_vortex_centers_cross_front_vortex_leg
                                                          * vortex_strengths_expanded)
        self.forces_on_back_vortices_in_geometry_axes = (density
                                                         * velocity_at_back_vortex_centers_cross_back_vortex_leg
                                                         * vortex_strengths_expanded)
        self.forces_on_left_vortices_in_geometry_axes = (density
                                                         * velocity_at_left_vortex_centers_cross_left_vortex_leg
                                                         * vortex_strengths_expanded)
        self.forces_on_right_vortices_in_geometry_axes = (density
                                                          * velocity_at_right_vortex_centers_cross_right_vortex_leg
                                                          * vortex_strengths_expanded)

        if self.verbose:
            print("Finished calculating forces on each vortex leg!")

        # Calculate total forces
        if self.verbose:
            print("Calculating total forces...")

        self.forces_on_panels_in_geometry_axes = np.empty((0, 3))
        for panel in range(self.n_panels):
            forces_or_front_vortex_in_geometry_axes = self.forces_on_front_vortices_in_geometry_axes[panel]
            forces_or_back_vortex_in_geometry_axes = self.forces_on_back_vortices_in_geometry_axes[panel]
            forces_or_left_vortex_in_geometry_axes = self.forces_on_left_vortices_in_geometry_axes[panel]
            forces_or_right_vortex_in_geometry_axes = self.forces_on_right_vortices_in_geometry_axes[panel]

            forces_on_panel_legs = np.vstack((forces_or_front_vortex_in_geometry_axes,
                                              forces_or_back_vortex_in_geometry_axes,
                                              forces_or_left_vortex_in_geometry_axes,
                                              forces_or_right_vortex_in_geometry_axes))

            total_forces_on_panel = np.sum(forces_on_panel_legs, axis=0)

            self.forces_on_panels_in_geometry_axes = np.vstack((self.forces_on_panels_in_geometry_axes,
                                                                total_forces_on_panel))

        self.Ftotal_geometry = np.sum(self.forces_on_panels_in_geometry_axes,
                                      axis=0)  # Remember, this is in geometry axes, not wind axes or body axes

        self.Ftotal_wind = np.transpose(
            self.operating_point.compute_rotation_matrix_wind_to_geometry()) @ self.Ftotal_geometry

        # Calculate nondimensional forces
        q = self.operating_point.dynamic_pressure()
        s_ref = self.airplane.s_ref
        self.CL = -self.Ftotal_wind[2] / q / s_ref
        self.CDi = -self.Ftotal_wind[0] / q / s_ref
        self.CY = self.Ftotal_wind[1] / q / s_ref

        if self.verbose:
            print("\nForces\n-----")
            print("CL: ", self.CL)
            print("CDi: ", self.CDi)
            print("CY: ", self.CY)
            print("Finished calculating total forces!")

    # ToDo: Properly cite and document this method.
    def get_induced_velocity_at_point(self, point):
        """

        :param point:
        :return:
        """
        # Input: a Nx3 numpy array of points that you would like to know the induced velocities at.
        # Output: a Nx3 numpy array of the induced velocities at those points.
        point = np.reshape(point, (-1, 3))

        Vij = self.calculate_velocity_influences(point)

        vortex_strengths_expanded = np.expand_dims(self.vortex_strengths, 1)

        Vi_x = Vij[:, :, 0] @ vortex_strengths_expanded
        Vi_y = Vij[:, :, 1] @ vortex_strengths_expanded
        Vi_z = Vij[:, :, 2] @ vortex_strengths_expanded

        Vi = np.hstack((Vi_x, Vi_y, Vi_z))

        return Vi

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
