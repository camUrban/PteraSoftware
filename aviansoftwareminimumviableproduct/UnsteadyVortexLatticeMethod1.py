from aviansoftwareminimumviableproduct import *
import numpy as np
import pyvista as pv
import aerosandbox_legacy_v0 as asl


class UnsteadyVortexLatticeMethod1(UnsteadyAerodynamicsProblem):

    def __init__(self, airplane, operating_point, movement, simulation_duration, simulation_time_step):
        super().__init__(airplane, operating_point, movement, simulation_duration, simulation_time_step)
        self.verbose = True
        self.front_left_vertices = None
        self.front_right_vertices = None
        self.back_left_vertices = None
        self.back_right_vertices = None
        self.front_wing_vortex_centers = None
        self.back_wing_vortex_centers = None
        self.left_wing_vortex_centers = None
        self.right_wing_vortex_centers = None
        self.front_wing_vortex_legs = None
        self.right_wing_vortex_legs = None
        self.back_wing_vortex_legs = None
        self.left_wing_vortex_legs = None
        self.areas = None
        self.is_trailing_edge = None
        self.is_leading_edge = None
        self.collocation_points = None
        self.normal_directions = None
        self.n_panels = None

        self.time_steps = None
        self.sweep_angles = None
        self.current_time = 0
        self.current_sweep_angle = 0

        self.initial_front_left_vertices = None
        self.initial_front_right_vertices = None
        self.initial_back_left_vertices = None
        self.initial_back_right_vertices = None

    def run(self):

        self.time_steps = np.arange(0, self.simulation_duration + self.simulation_time_step, self.simulation_time_step)
        self.sweep_angles = np.arcsin(np.sin(self.movement.sweeping_amplitude / 2) * np.sin(
            2 * np.pi * self.time_steps / self.movement.movement_period))

        self.make_gif()

        # for i in range(len(self.time_steps)):
        #     self.current_time = self.time_steps[i]
        #     self.current_sweep_angle = self.sweep_angles[i]
        #     if self.current_time == 0:
        #         self.make_panels()
        #
        #         self.initial_front_left_vertices = self.front_left_vertices
        #         self.initial_front_right_vertices = self.front_right_vertices
        #         self.initial_back_left_vertices = self.back_left_vertices
        #         self.initial_back_right_vertices = self.back_right_vertices
        #         self.draw()
        #     else:
        #         self.move_panels()
        #         self.draw()

    # ToDo: Properly comment and cite the code in the following method
    def make_gif(self):

        plotter = pv.Plotter()

        for i in range(len(self.time_steps)):
            self.current_time = self.time_steps[i]
            self.current_sweep_angle = self.sweep_angles[i]
            if self.current_time == 0:
                self.make_panels()

                vertices = np.vstack((
                    self.front_left_vertices,
                    self.front_right_vertices,
                    self.back_right_vertices,
                    self.back_left_vertices
                ))
                faces = np.transpose(np.vstack((
                    4 * np.ones(self.n_panels),
                    np.arange(self.n_panels),
                    np.arange(self.n_panels) + self.n_panels,
                    np.arange(self.n_panels) + 2 * self.n_panels,
                    np.arange(self.n_panels) + 3 * self.n_panels,
                )))
                faces = np.reshape(faces, (-1), order='C')
                wing_surfaces = pv.PolyData(vertices, faces)

                plotter.add_mesh(wing_surfaces, color='white', show_edges=True,
                                 smooth_shading=False)

                print("Press \"q\" to close window and create an animation")

                plotter.show(auto_close=False)
                plotter.open_gif("flapping.gif")

                self.initial_front_left_vertices = self.front_left_vertices
                self.initial_front_right_vertices = self.front_right_vertices
                self.initial_back_left_vertices = self.back_left_vertices
                self.initial_back_right_vertices = self.back_right_vertices

                plotter.write_frame()
                plotter.clear()
            else:
                self.move_panels()

                vertices = np.vstack((
                    self.front_left_vertices,
                    self.front_right_vertices,
                    self.back_right_vertices,
                    self.back_left_vertices
                ))
                faces = np.transpose(np.vstack((
                    4 * np.ones(self.n_panels),
                    np.arange(self.n_panels),
                    np.arange(self.n_panels) + self.n_panels,
                    np.arange(self.n_panels) + 2 * self.n_panels,
                    np.arange(self.n_panels) + 3 * self.n_panels,
                )))
                faces = np.reshape(faces, (-1), order='C')
                wing_surfaces = pv.PolyData(vertices, faces)

                plotter.add_mesh(wing_surfaces, color='white', show_edges=True,
                                 smooth_shading=False)

                plotter.write_frame()
                plotter.clear()

        # Close movie and delete object
        plotter.close()

    # Adapted from:         vlm3.make_panels() in AeroSandbox
    # Author:               Peter Sharpe
    # Date of Retrieval:    03/28/2020
    def make_panels(self):
        # Creates self.panel_coordinates_structured_list and self.wing_mcl_normals.

        if self.verbose:
            print("Meshing...")

        collocation_points = np.empty((0, 3))
        normal_directions = np.empty((0, 3))
        left_vortex_vertices = np.empty((0, 3))
        right_vortex_vertices = np.empty((0, 3))
        front_left_vertices = np.empty((0, 3))
        front_right_vertices = np.empty((0, 3))
        back_left_vertices = np.empty((0, 3))
        back_right_vertices = np.empty((0, 3))
        areas = np.empty(0)
        is_trailing_edge = np.empty(0, dtype=bool)
        is_leading_edge = np.empty(0, dtype=bool)
        xsec = None

        for wing_num in range(len(self.airplane.wings)):
            # For each wing, we want (where M is the number of chordwise panels, N is the number of spanwise panels):
            #   panel_coordinates_structured_list: M+1 x N+1 x 3; corners of every panel.
            #   normals_structured_list: M x N x 3; normal direction of each panel

            # Get the wing
            wing = self.airplane.wings[wing_num]

            # Define number of chordwise points
            n_chordwise_coordinates = wing.chordwise_panels + 1

            # Get the chordwise coordinates
            if wing.chordwise_spacing == 'uniform':
                nondim_chordwise_coordinates = np.linspace(0, 1, n_chordwise_coordinates)
            elif wing.chordwise_spacing == 'cosine':
                nondim_chordwise_coordinates = asl.cosspace(0, 1, n_chordwise_coordinates)
            else:
                raise Exception("Bad value of wing.chordwise_spacing!")

            # Get corners of xsecs
            xsec_xyz_le = np.empty((0, 3))  # Nx3 array of leading edge points
            xsec_xyz_te = np.empty((0, 3))  # Nx3 array of trailing edge points
            for xsec in wing.xsecs:
                xsec_xyz_le = np.vstack((xsec_xyz_le, xsec.xyz_le + wing.xyz_le))
                xsec_xyz_te = np.vstack((xsec_xyz_te, xsec.xyz_te() + wing.xyz_le))

            # Get quarter-chord vector
            xsec_xyz_quarter_chords = 0.75 * xsec_xyz_le + 0.25 * xsec_xyz_te  # Nx3 array of quarter-chord points
            section_quarter_chords = (
                    xsec_xyz_quarter_chords[1:, :] -
                    xsec_xyz_quarter_chords[:-1, :]
            )  # Nx3 array of vectors connecting quarter-chords

            # Get directions for transforming 2D airfoil data to 3D

            # First, project quarter chords onto YZ plane and normalize.
            section_quarter_chords_proj = (section_quarter_chords[:, 1:] /
                                           np.expand_dims(np.linalg.norm(section_quarter_chords[:, 1:], axis=1),
                                                          axis=1))  # Nx2 array of quarter-chord vectors projected on YZ
            section_quarter_chords_proj = np.hstack(
                (np.zeros((section_quarter_chords_proj.shape[0], 1)),
                 section_quarter_chords_proj)
            )  # Convert back to a Nx3 array, since that's what we'll need later.

            # Then, construct the normal directions for each xsec.
            if len(wing.xsecs) > 2:  # Make normals for the inner xsecs, where we need to merge directions
                xsec_local_normal_inners = section_quarter_chords_proj[:-1, :] + section_quarter_chords_proj[1:, :]
                xsec_local_normal_inners = (xsec_local_normal_inners /
                                            np.expand_dims(np.linalg.norm(xsec_local_normal_inners, axis=1), axis=1)
                                            )
                xsec_local_normal = np.vstack((
                    section_quarter_chords_proj[0, :],
                    xsec_local_normal_inners,
                    section_quarter_chords_proj[-1, :]
                ))
            else:
                xsec_local_normal = np.vstack((
                    section_quarter_chords_proj[0, :],
                    section_quarter_chords_proj[-1, :]
                ))

            # xsec_local_normal is now a Nx3 array that represents the normal direction at each xsec.

            # Then, construct the back directions for each xsec.
            xsec_local_back = xsec_xyz_te - xsec_xyz_le  # aligned with chord
            xsec_chord = np.linalg.norm(xsec_local_back, axis=1)  # 1D vector, one per xsec
            xsec_local_back = (xsec_local_back /
                               np.expand_dims(xsec_chord, axis=1)
                               )

            # Then, construct the up direction for each xsec.
            xsec_local_up = np.cross(xsec_local_back, xsec_local_normal,
                                     axis=1)  # Nx3 array that represents the upwards direction at each xsec.

            # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to compensate)
            xsec_scaling_factor = 1 / np.sqrt(
                (1 + np.sum(section_quarter_chords_proj[1:, :] * section_quarter_chords_proj[:-1, :], axis=1)
                 ) / 2
            )
            xsec_scaling_factor = np.hstack((1, xsec_scaling_factor, 1))  # TODO Make sure this is always right.

            # Make the panels for each section.

            for section_num in range(len(wing.xsecs) - 1):
                # Define the relevant cross sections
                inner_xsec = wing.xsecs[section_num]  # type is WingXSec
                outer_xsec = wing.xsecs[section_num + 1]  # type is WingXSec

                # Define the airfoils at each cross section
                inner_airfoil = inner_xsec.airfoil.add_control_surface(
                    deflection=inner_xsec.control_surface_deflection,
                    hinge_point=inner_xsec.control_surface_hinge_point
                )
                outer_airfoil = outer_xsec.airfoil.add_control_surface(
                    deflection=inner_xsec.control_surface_deflection,
                    # inner xsec dictates control surface deflections.
                    hinge_point=inner_xsec.control_surface_hinge_point
                )

                # Make the mean camber lines for each.
                inner_xsec_mcl_nondim = inner_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)
                outer_xsec_mcl_nondim = outer_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)
                # inner_xsec_mcl: First index is point number, second index is xyz.
                inner_xsec_mcl = xsec_xyz_le[section_num, :] + (xsec_local_back[section_num, :] *
                                                                np.expand_dims(inner_xsec_mcl_nondim[:, 0], 1) *
                                                                xsec_chord[section_num] + xsec_local_up[section_num, :]
                                                                * np.expand_dims(inner_xsec_mcl_nondim[:, 1], 1) *
                                                                xsec_chord[section_num] *
                                                                xsec_scaling_factor[section_num]
                                                                )
                outer_xsec_mcl = xsec_xyz_le[section_num + 1, :] + (
                        xsec_local_back[section_num + 1, :] * np.expand_dims(outer_xsec_mcl_nondim[:, 0], 1) *
                        xsec_chord[section_num + 1] +
                        xsec_local_up[section_num + 1, :] * np.expand_dims(outer_xsec_mcl_nondim[:, 1], 1) * xsec_chord[
                            section_num + 1] * xsec_scaling_factor[
                            section_num + 1]
                )

                # Define number of spanwise points
                n_spanwise_coordinates = xsec.spanwise_panels + 1

                # Get the spanwise coordinates
                if xsec.spanwise_spacing == 'uniform':
                    nondim_spanwise_coordinates = np.linspace(0, 1, n_spanwise_coordinates)
                elif xsec.spanwise_spacing == 'cosine':
                    nondim_spanwise_coordinates = asl.cosspace(n_points=n_spanwise_coordinates)
                else:
                    raise Exception("Bad value of section.spanwise_spacing!")

                # Make section_mcl_coordinates: MxNx3 array of mean camberline coordinates.
                # First index is chordwise location, second index is spanwise location, third is xyz.
                section_mcl_coordinates = (
                        np.expand_dims(np.expand_dims((1 - nondim_spanwise_coordinates), 0),
                                       2) * np.expand_dims(inner_xsec_mcl, 1)
                        + np.expand_dims(np.expand_dims(nondim_spanwise_coordinates, 0), 2)
                        * np.expand_dims(outer_xsec_mcl, 1)
                )  # TODO Make this work for large twist angles.

                # Compute corners of each panel
                front_inner_coordinates = section_mcl_coordinates[:-1, :-1, :]
                front_outer_coordinates = section_mcl_coordinates[:-1, 1:, :]
                back_inner_coordinates = section_mcl_coordinates[1:, :-1, :]
                back_outer_coordinates = section_mcl_coordinates[1:, 1:, :]
                section_is_trailing_edge = np.vstack((
                    np.zeros((wing.chordwise_panels - 1, xsec.spanwise_panels), dtype=bool),
                    np.ones((1, xsec.spanwise_panels), dtype=bool)
                ))
                section_is_leading_edge = np.vstack((
                    np.ones((1, xsec.spanwise_panels), dtype=bool),
                    np.zeros((wing.chordwise_panels - 1, xsec.spanwise_panels), dtype=bool)
                ))

                # Reshape
                front_inner_coordinates = np.reshape(front_inner_coordinates, (-1, 3), order='F')
                front_outer_coordinates = np.reshape(front_outer_coordinates, (-1, 3), order='F')
                back_inner_coordinates = np.reshape(back_inner_coordinates, (-1, 3), order='F')
                back_outer_coordinates = np.reshape(back_outer_coordinates, (-1, 3), order='F')
                section_is_trailing_edge = np.reshape(section_is_trailing_edge, (-1), order='F')
                section_is_leading_edge = np.reshape(section_is_leading_edge, (-1), order='F')

                # Calculate panel normals and areas via diagonals
                diag1 = front_outer_coordinates - back_inner_coordinates
                diag2 = front_inner_coordinates - back_outer_coordinates
                diag_cross = np.cross(diag1, diag2, axis=1)
                diag_cross_norm = np.linalg.norm(diag_cross, axis=1)
                normals_to_add = diag_cross / np.expand_dims(diag_cross_norm, axis=1)
                areas_to_add = diag_cross_norm / 2

                # Make the panel data
                collocations_to_add = (
                        0.5 * (0.25 * front_inner_coordinates + 0.75 * back_inner_coordinates) +
                        0.5 * (0.25 * front_outer_coordinates + 0.75 * back_outer_coordinates)
                )
                inner_vortex_vertices_to_add = 0.75 * front_inner_coordinates + 0.25 * back_inner_coordinates
                outer_vortex_vertices_to_add = 0.75 * front_outer_coordinates + 0.25 * back_outer_coordinates

                # Append to the lists of panel data (c, n, lv, rv, etc.)
                front_left_vertices = np.vstack((
                    front_left_vertices,
                    front_inner_coordinates
                ))
                front_right_vertices = np.vstack((
                    front_right_vertices,
                    front_outer_coordinates
                ))
                back_left_vertices = np.vstack((
                    back_left_vertices,
                    back_inner_coordinates
                ))
                back_right_vertices = np.vstack((
                    back_right_vertices,
                    back_outer_coordinates
                ))
                areas = np.hstack((
                    areas,
                    areas_to_add
                ))
                is_trailing_edge = np.hstack((
                    is_trailing_edge,
                    section_is_trailing_edge
                ))
                is_leading_edge = np.hstack((
                    is_leading_edge,
                    section_is_leading_edge
                ))
                collocation_points = np.vstack((
                    collocation_points,
                    collocations_to_add
                ))
                normal_directions = np.vstack((
                    normal_directions,
                    normals_to_add
                ))
                left_vortex_vertices = np.vstack((
                    left_vortex_vertices,
                    inner_vortex_vertices_to_add
                ))
                right_vortex_vertices = np.vstack((
                    right_vortex_vertices,
                    outer_vortex_vertices_to_add
                ))

                # Handle symmetry
                if wing.symmetric:
                    if inner_xsec.control_surface_type == "asymmetric":
                        # Define the airfoils at each cross section
                        inner_airfoil = inner_xsec.airfoil.add_control_surface(
                            deflection=-inner_xsec.control_surface_deflection,
                            hinge_point=inner_xsec.control_surface_hinge_point
                        )
                        outer_airfoil = outer_xsec.airfoil.add_control_surface(
                            deflection=-inner_xsec.control_surface_deflection,
                            # inner xsec dictates control surface deflections.
                            hinge_point=inner_xsec.control_surface_hinge_point
                        )

                        # Make the mean camber lines for each.
                        inner_xsec_mcl_nondim = inner_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)
                        outer_xsec_mcl_nondim = outer_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)

                        # inner_xsec_mcl: First index is point number, second index is xyz.
                        inner_xsec_mcl = xsec_xyz_le[section_num, :] + (
                                xsec_local_back[section_num, :] * np.expand_dims(inner_xsec_mcl_nondim[:, 0], 1) *
                                xsec_chord[
                                    section_num] +
                                xsec_local_up[section_num, :] * np.expand_dims(inner_xsec_mcl_nondim[:, 1], 1) *
                                xsec_chord[
                                    section_num] * xsec_scaling_factor[section_num]
                        )

                        outer_xsec_mcl = xsec_xyz_le[section_num + 1, :] + (
                                xsec_local_back[section_num + 1, :] * np.expand_dims(outer_xsec_mcl_nondim[:, 0], 1) *
                                xsec_chord[section_num + 1] +
                                xsec_local_up[section_num + 1, :] * np.expand_dims(outer_xsec_mcl_nondim[:, 1], 1) *
                                xsec_chord[
                                    section_num + 1] * xsec_scaling_factor[
                                    section_num + 1]
                        )

                        # Make section_mcl_coordinates: MxNx3 array of mean camberline coordinates.
                        # First index is chordwise location, second index is spanwise location, third is xyz.
                        section_mcl_coordinates = (
                                np.expand_dims(np.expand_dims((1 - nondim_spanwise_coordinates), 0),
                                               2) * np.expand_dims(inner_xsec_mcl, 1)
                                + np.expand_dims(np.expand_dims(nondim_spanwise_coordinates, 0), 2)
                                * np.expand_dims(outer_xsec_mcl, 1)
                        )  # TODO Make this work for large twist angles.

                        # Compute corners of each panel
                        front_inner_coordinates = section_mcl_coordinates[:-1, :-1, :]
                        front_outer_coordinates = section_mcl_coordinates[:-1, 1:, :]
                        back_inner_coordinates = section_mcl_coordinates[1:, :-1, :]
                        back_outer_coordinates = section_mcl_coordinates[1:, 1:, :]

                        # Reshape
                        front_inner_coordinates = np.reshape(front_inner_coordinates, (-1, 3), order='F')
                        front_outer_coordinates = np.reshape(front_outer_coordinates, (-1, 3), order='F')
                        back_inner_coordinates = np.reshape(back_inner_coordinates, (-1, 3), order='F')
                        back_outer_coordinates = np.reshape(back_outer_coordinates, (-1, 3), order='F')

                        # Calculate panel normals and areas via diagonals
                        diag1 = front_outer_coordinates - back_inner_coordinates
                        diag2 = front_inner_coordinates - back_outer_coordinates
                        diag_cross = np.cross(diag1, diag2, axis=1)
                        diag_cross_norm = np.linalg.norm(diag_cross, axis=1)
                        normals_to_add = diag_cross / np.expand_dims(diag_cross_norm, axis=1)
                        areas_to_add = diag_cross_norm / 2

                        # Make the panels and append them to the lists of panel data (c, n, lv, rv, etc.)
                        collocations_to_add = (
                                0.5 * (0.25 * front_inner_coordinates + 0.75 * back_inner_coordinates) +
                                0.5 * (0.25 * front_outer_coordinates + 0.75 * back_outer_coordinates)
                        )
                        inner_vortex_vertices_to_add = 0.75 * front_inner_coordinates + 0.25 * back_inner_coordinates
                        outer_vortex_vertices_to_add = 0.75 * front_outer_coordinates + 0.25 * back_outer_coordinates

                    front_left_vertices = np.vstack((
                        front_left_vertices,
                        asl.reflect_over_XZ_plane(front_outer_coordinates)
                    ))
                    front_right_vertices = np.vstack((
                        front_right_vertices,
                        asl.reflect_over_XZ_plane(front_inner_coordinates)
                    ))
                    back_left_vertices = np.vstack((
                        back_left_vertices,
                        asl.reflect_over_XZ_plane(back_outer_coordinates)
                    ))
                    back_right_vertices = np.vstack((
                        back_right_vertices,
                        asl.reflect_over_XZ_plane(back_inner_coordinates)
                    ))
                    areas = np.hstack((
                        areas,
                        areas_to_add
                    ))
                    is_trailing_edge = np.hstack((
                        is_trailing_edge,
                        section_is_trailing_edge
                    ))
                    is_leading_edge = np.hstack((
                        is_leading_edge,
                        section_is_leading_edge
                    ))
                    collocation_points = np.vstack((
                        collocation_points,
                        asl.reflect_over_XZ_plane(collocations_to_add)
                    ))
                    normal_directions = np.vstack((
                        normal_directions,
                        asl.reflect_over_XZ_plane(normals_to_add)
                    ))
                    left_vortex_vertices = np.vstack((
                        left_vortex_vertices,
                        asl.reflect_over_XZ_plane(outer_vortex_vertices_to_add)
                    ))
                    right_vortex_vertices = np.vstack((
                        right_vortex_vertices,
                        asl.reflect_over_XZ_plane(inner_vortex_vertices_to_add)
                    ))

        # Write to self object
        self.front_left_vertices = front_left_vertices
        self.front_right_vertices = front_right_vertices
        self.back_left_vertices = back_left_vertices
        self.back_right_vertices = back_right_vertices

        self.front_wing_vortex_centers = (front_left_vertices + front_right_vertices) / 2
        self.back_wing_vortex_centers = (back_left_vertices + back_right_vertices) / 2
        self.left_wing_vortex_centers = (front_left_vertices + back_left_vertices) / 2
        self.right_wing_vortex_centers = (front_right_vertices + back_right_vertices) / 2

        self.front_wing_vortex_legs = (front_right_vertices - front_left_vertices)
        self.right_wing_vortex_legs = (back_right_vertices - front_right_vertices)
        self.back_wing_vortex_legs = (back_left_vertices - back_right_vertices)
        self.left_wing_vortex_legs = (front_left_vertices - back_left_vertices)

        self.areas = areas
        self.is_trailing_edge = is_trailing_edge
        self.is_leading_edge = is_leading_edge
        self.collocation_points = collocation_points
        self.normal_directions = normal_directions
        self.n_panels = len(self.collocation_points)

        if self.verbose:
            print("Meshing complete!")
        # -----------------------------------------------------
        # Review of the important things that have been done up to this point:
        #   We made panel_coordinates_structured_list, a MxNx3 array describing a structured quadrilateral mesh of the
        #   wing's mean camber surface.
        #       For reference: first index is chordwise coordinate, second index is spanwise coordinate, and third index
        #       is xyz.
        #   We made normals_structured_list, a MxNx3 array describing the normal direction of the mean camber surface at
        #   the collocation point.
        #       For reference: first index is chordwise coordinate, second index is spanwise coordinate, and third index
        #       is xyz.
        #       Takes into account control surface deflections
        #   Both panel_coordinates_structured_list and normals_structured_list have been appended to lists of ndarrays
        #   within the vlm2 class, accessible at self.panel_coordinates_structured_list and
        #   self.normals_structured_list, respectively.
        #   Control surface handling:
        #       Control surfaces are implemented into normal directions as intended.
        #   Symmetry handling:
        #       All symmetric wings have been split into separate halves.
        #       All wing halves have their spanwise coordinates labeled from the left side of the airplane to the right.
        #       Control surface deflection symmetry has been handled; this is encoded into the normal directions.

    # Adapted from:         vlm3.make_panels() in AeroSandbox
    # Author:               Peter Sharpe
    # Date of Retrieval:    03/28/2020
    def draw(self, points_type=None):

        if self.verbose:
            print("Drawing...")

        # Initialize Plotter
        plotter = pv.Plotter()

        # Make airplane geometry
        vertices = np.vstack((
            self.front_left_vertices,
            self.front_right_vertices,
            self.back_right_vertices,
            self.back_left_vertices
        ))
        faces = np.transpose(np.vstack((
            4 * np.ones(self.n_panels),
            np.arange(self.n_panels),
            np.arange(self.n_panels) + self.n_panels,
            np.arange(self.n_panels) + 2 * self.n_panels,
            np.arange(self.n_panels) + 3 * self.n_panels,
        )))
        faces = np.reshape(faces, (-1), order='C')
        wing_surfaces = pv.PolyData(vertices, faces)

        plotter.add_mesh(wing_surfaces, color='white', show_edges=True,
                         smooth_shading=True)

        # Points
        if points_type is not None:
            points = getattr(self, points_type)

            plotter.add_points(points)

        # Do the plotting
        plotter.set_background(color="black")
        plotter.show(cpos=(-1, -1, 1), full_screen=True)

        if self.verbose:
            print("Drawing finished!")

    # ToDo: Properly comment code in the following method
    def move_panels(self):

        new_front_left_vertices = self.initial_front_left_vertices
        new_front_right_vertices = self.initial_front_right_vertices
        new_back_left_vertices = self.initial_back_left_vertices
        new_back_right_vertices = self.initial_back_right_vertices

        front_left_vertex_sweep_angles = []
        front_right_vertex_sweep_angles = []
        back_left_vertex_sweep_angles = []
        back_right_vertex_sweep_angles = []

        for vertex in range(self.n_panels):
            if self.initial_front_left_vertices[vertex, 1] < 0:
                front_left_vertex_sweep_angles.append(np.pi - self.current_sweep_angle)
            else:
                front_left_vertex_sweep_angles.append(self.current_sweep_angle)

            if self.initial_front_right_vertices[vertex, 1] < 0:
                front_right_vertex_sweep_angles.append(np.pi - self.current_sweep_angle)
            else:
                front_right_vertex_sweep_angles.append(self.current_sweep_angle)

            if self.initial_back_left_vertices[vertex, 1] < 0:
                back_left_vertex_sweep_angles.append(np.pi - self.current_sweep_angle)
            else:
                back_left_vertex_sweep_angles.append(self.current_sweep_angle)

            if self.initial_front_right_vertices[vertex, 1] < 0:
                back_right_vertex_sweep_angles.append(np.pi - self.current_sweep_angle)
            else:
                back_right_vertex_sweep_angles.append(self.current_sweep_angle)

        front_left_radii = np.sqrt(
            self.initial_front_left_vertices[:, 1] ** 2 + self.initial_front_left_vertices[:, 2] ** 2)
        front_right_radii = np.sqrt(
            self.initial_front_right_vertices[:, 1] ** 2 + self.initial_front_right_vertices[:, 2] ** 2)
        back_left_radii = np.sqrt(
            self.initial_back_left_vertices[:, 1] ** 2 + self.initial_back_left_vertices[:, 2] ** 2)
        back_right_radii = np.sqrt(
            self.initial_back_right_vertices[:, 1] ** 2 + self.initial_back_right_vertices[:, 2] ** 2)

        new_front_left_vertices[:, 1] = front_left_radii * np.cos(front_left_vertex_sweep_angles)
        new_front_left_vertices[:, 2] = front_left_radii * np.sin(front_left_vertex_sweep_angles)

        new_front_right_vertices[:, 1] = front_right_radii * np.cos(front_right_vertex_sweep_angles)
        new_front_right_vertices[:, 2] = front_right_radii * np.sin(front_right_vertex_sweep_angles)

        new_back_left_vertices[:, 1] = back_left_radii * np.cos(back_left_vertex_sweep_angles)
        new_back_left_vertices[:, 2] = back_left_radii * np.sin(back_left_vertex_sweep_angles)

        new_back_right_vertices[:, 1] = back_right_radii * np.cos(back_right_vertex_sweep_angles)
        new_back_right_vertices[:, 2] = back_right_radii * np.sin(back_right_vertex_sweep_angles)

        self.front_left_vertices = new_front_left_vertices
        self.front_right_vertices = new_front_right_vertices
        self.back_left_vertices = new_back_left_vertices
        self.back_right_vertices = new_back_right_vertices

        self.front_wing_vortex_centers = (new_front_left_vertices + new_front_right_vertices) / 2
        self.back_wing_vortex_centers = (new_back_left_vertices + new_back_right_vertices) / 2
        self.left_wing_vortex_centers = (new_front_left_vertices + new_back_left_vertices) / 2
        self.right_wing_vortex_centers = (new_front_right_vertices + new_back_right_vertices) / 2

        self.front_wing_vortex_legs = (new_front_right_vertices - new_front_left_vertices)
        self.right_wing_vortex_legs = (new_back_right_vertices - new_front_right_vertices)
        self.back_wing_vortex_legs = (new_back_left_vertices - new_back_right_vertices)
        self.left_wing_vortex_legs = (new_front_left_vertices - new_back_left_vertices)
