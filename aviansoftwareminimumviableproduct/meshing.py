
"""This module contains useful functions for creating meshes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    mesh_wing: This function takes in an object of the Wing class and creates a quadrilateral mesh of its geometry,
               and then populates the object's panels with the mesh data.
    move_panels: This function takes in a problem of the UnsteadyAerodynamicsProblem class, and modifies it's panel
                 locations as it time steps through the simulation.
"""

import aviansoftwareminimumviableproduct as asmvp
import numpy as np


def mesh_wing(wing):
    """This function takes in an object of the Wing class and creates a quadrilateral mesh of its geometry,
    and then populates the object's panels with the mesh data.

    Citation:
        Adapted from:         vlm3.make_panels in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    05/01/2020

    :param wing: Wing
        This is the wing to be meshed.
    :return: None
    """

    # Define the number of chordwise panels and points.
    num_chordwise_panels = wing.num_chordwise_panels
    num_chordwise_coordinates = num_chordwise_panels + 1

    # Initialize the list of cross sections to None.
    wing_cross_section = None

    # Initialize an empty ndarray that will hold the panels of this wing. It currently has 0 columns and M rows,
    # where M is the number of the wing's chordwise panels.
    panels = np.empty((num_chordwise_panels, 0), dtype=object)

    # Get the chordwise coordinates.
    if wing.chordwise_spacing == 'uniform':
        nondim_chordwise_coordinates = np.linspace(0, 1, num_chordwise_coordinates)
    elif wing.chordwise_spacing == 'cosine':
        nondim_chordwise_coordinates = asmvp.geometry.cosspace(0, 1, num_chordwise_coordinates)
    else:
        raise Exception("Bad value of wing.chordwise_spacing!")

    # Initialize two empty 0 x 3 ndarrays to hold the corners of each cross section. They will eventually be L x 3
    # ndarrays, where L is number of cross sections.
    cross_section_xyz_le = np.empty((0, 3))
    cross_section_xyz_te = np.empty((0, 3))

    # Iterate through the meshed wing cross sections and vertically stack the global location each cross sections
    # leading and trailing edges. cross_section.xyz_te is a method that returns the cross section's trailing edge's
    # coordinates.
    for wing_cross_section in wing.wing_cross_sections:
        cross_section_xyz_le = np.vstack((cross_section_xyz_le, wing_cross_section.xyz_le + wing.xyz_le))
        cross_section_xyz_te = np.vstack((cross_section_xyz_te, wing_cross_section.xyz_te() + wing.xyz_le))

    # Get the quarter chord vectors, which are a L x 3 ndarray of points which are the quarter-chord points of cross
    # section, where L is the number of cross sections.
    cross_section_xyz_quarter_chords = cross_section_xyz_le + 0.25 * (cross_section_xyz_te - cross_section_xyz_le)

    # Get a (L - 1) x 3 ndarray of vectors connecting the cross section quarter chord points, where L is the number of
    # cross sections.
    section_quarter_chords = (cross_section_xyz_quarter_chords[1:, :] - cross_section_xyz_quarter_chords[:-1, :])

    # Get directions for transforming 2D airfoil data to 3D:
    #   Project quarter chords onto yz plane and normalize.
    #   Create a L x 2 ndarray with just the y and z components of the the section quarter chord vectors.
    section_quarter_chords_yz = section_quarter_chords[:, 1:]

    # Create a list of the magnitudes of each row of the section_quarter_chords_yz ndarray.
    section_quarter_chords_yz_magnitude_list = np.linalg.norm(section_quarter_chords_yz, axis=1)

    # Convert section_quarter_chords_yz_magnitude_list into a column vector.
    section_quarter_chords_yz_magnitude_column_vector = np.expand_dims(section_quarter_chords_yz_magnitude_list,
                                                                       axis=1)
    # Normalize the y and z components by the magnitudes
    section_quarter_chords_yz_norm_magnitudes = (section_quarter_chords_yz /
                                                 section_quarter_chords_yz_magnitude_column_vector)

    # Calculate the number of quarter chord vectors
    num_quarter_chords = section_quarter_chords_yz_magnitude_column_vector.shape[0]

    # Create a column vector of all zeros with height equal to the number of quarter chord vectors
    zero_column_vector_stand_in_for_quarter_chords_x_values = np.zeros((num_quarter_chords, 1))

    # Horizontally stack the zero column vector with the section_quarter_chords_yz_norm_magnitudes to
    # produce the normalized section quarter chords projected onto the yz plane.
    section_quarter_chords_proj_yz_norm = np.hstack((zero_column_vector_stand_in_for_quarter_chords_x_values,
                                                     section_quarter_chords_yz_norm_magnitudes))

    # Then, construct the normal directions for each cross_section. Make the normals for the inner
    # wing_cross_sections, where we need to merge directions.
    if len(wing.wing_cross_sections) > 2:
        # Add together the adjacent normalized section quarter chords projected onto the the yz plane.
        cross_sections_local_normal_inners_non_norm = (section_quarter_chords_proj_yz_norm[:-1, :] +
                                                       section_quarter_chords_proj_yz_norm[1:, :])

        # Create a list of the magnitudes of the summed adjacent normalized section quarter chords projected onto the yz
        # plane.
        cross_sections_local_normal_inners_mag_list = np.linalg.norm(cross_sections_local_normal_inners_non_norm,
                                                                     axis=1)

        # Convert the list to a column vector.
        cross_section_local_normal_inners_mag_column_vector = np.expand_dims(
            cross_sections_local_normal_inners_mag_list, axis=1)

        # Normalize the summed adjacent normalized section quarter chords projected onto the yz plane by their
        # magnitudes.
        cross_section_local_normal_inners_norm = (cross_sections_local_normal_inners_non_norm /
                                                  cross_section_local_normal_inners_mag_column_vector)

        # Vertically stack the first normalized section quarter chord, the inner normalized section quarter chords, and
        # the last normalized section quarter chord.
        cross_sections_local_normal = np.vstack((
            section_quarter_chords_proj_yz_norm[0, :],
            cross_section_local_normal_inners_norm,
            section_quarter_chords_proj_yz_norm[-1, :]
        ))
    else:
        # Vertically stack the first normalized section quarter chord, and the last normalized section quarter chord.
        cross_sections_local_normal = np.vstack((
            section_quarter_chords_proj_yz_norm[0, :],
            section_quarter_chords_proj_yz_norm[-1, :]
        ))
        # cross_sections_local_normal is now a L x 3 array that represents the normal direction at each cross section.

    # Then, construct the back directions for each cross section.
    cross_section_local_back_non_norm = cross_section_xyz_te - cross_section_xyz_le

    # Create a list of the cross section chord lengths.
    cross_section_chord_length_list = np.linalg.norm(cross_section_local_back_non_norm, axis=1)

    # Convert the list to a column vector.
    cross_section_chord_length_column_vector = np.expand_dims(cross_section_chord_length_list, axis=1)

    # Normalize the cross section back vectors by their magnitudes.
    cross_section_local_back_norm = (cross_section_local_back_non_norm / cross_section_chord_length_column_vector)

    # Then, construct the up direction for each cross sections.
    cross_section_local_up = np.cross(cross_section_local_back_norm, cross_sections_local_normal, axis=1)

    # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to compensate).
    cross_section_scaling_factor = 1 / np.sqrt(
        (1 + np.sum(section_quarter_chords_proj_yz_norm[1:, :] * section_quarter_chords_proj_yz_norm[:-1, :],
                    axis=1)
         ) / 2
    )
    cross_section_scaling_factor = np.hstack((1, cross_section_scaling_factor, 1))

    # Make the panels for each section.
    for section_num in range(len(wing.wing_cross_sections) - 1):
        # Define the relevant cross sections.
        inner_cross_section = wing.wing_cross_sections[section_num]
        outer_cross_section = wing.wing_cross_sections[section_num + 1]

        # Define the airfoils at each cross section.
        inner_airfoil = inner_cross_section.airfoil.add_control_surface(
            deflection=inner_cross_section.control_surface_deflection,
            hinge_point=inner_cross_section.control_surface_hinge_point
        )
        outer_airfoil = outer_cross_section.airfoil.add_control_surface(
            deflection=inner_cross_section.control_surface_deflection,

            # The inner cross section dictates control surface deflections.
            hinge_point=inner_cross_section.control_surface_hinge_point
        )

        # Make the mean camber lines for each cross section. First index is point number, second index is xyz.
        inner_cross_section_mcl_nondim = inner_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)
        outer_cross_section_mcl_nondim = outer_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)

        # Put the inner cross section's local up airfoil frame coordinates in a column vector.
        inner_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(inner_cross_section_mcl_nondim[:, 1], 1)

        # Put the inner cross section's local back airfoil frame coordinates in a column vector.
        inner_cross_section_mcl_nondim_local_back_column_vector = np.expand_dims(inner_cross_section_mcl_nondim[:, 0],
                                                                                 1)
        # Put the outer cross section's local up airfoil frame coordinates in a column vector.
        outer_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(outer_cross_section_mcl_nondim[:, 1], 1)

        # Put the outer cross section's local back airfoil frame coordinates in a column vector.
        outer_cross_section_mcl_nondim_local_back_column_vector = np.expand_dims(outer_cross_section_mcl_nondim[:, 0],
                                                                                 1)

        # Convert the inner cross section's non dimensional local back airfoil frame coordinates to meshed wing
        # coordinates.
        inner_cross_section_mcl_local_back = (cross_section_local_back_norm[section_num, :] *
                                              inner_cross_section_mcl_nondim_local_back_column_vector *
                                              cross_section_chord_length_list[section_num])

        # Convert the inner cross section's non dimensional local up airfoil frame coordinates to meshed wing
        # coordinates.
        inner_cross_section_mcl_local_up = (cross_section_local_up[section_num, :] *
                                            inner_cross_section_mcl_nondim_local_up_column_vector *
                                            cross_section_chord_length_list[section_num] *
                                            cross_section_scaling_factor[section_num])

        # Convert the outer cross section's non dimensional local back airfoil frame coordinates to meshed wing
        # coordinates.
        outer_cross_section_mcl_local_back = (cross_section_local_back_norm[section_num + 1, :] *
                                              outer_cross_section_mcl_nondim_local_back_column_vector *
                                              cross_section_chord_length_list[section_num + 1])

        # Convert the outer cross section's non dimensional local up airfoil frame coordinates to meshed wing
        # coordinates.
        outer_cross_section_mcl_local_up = (cross_section_local_up[section_num + 1, :] *
                                            outer_cross_section_mcl_nondim_local_up_column_vector *
                                            cross_section_chord_length_list[section_num + 1] *
                                            cross_section_scaling_factor[section_num + 1])

        # Convert the inner cross section's meshed wing coordinates to absolute coordinates. This is size M x 3, where M
        # is the number of chordwise points.
        inner_cross_section_mcl = (cross_section_xyz_le[section_num, :] + inner_cross_section_mcl_local_back
                                   + inner_cross_section_mcl_local_up)

        # Convert the outer cross section's meshed wing coordinates to absolute coordinates. This is size M x 3, where M
        # is the number of chordwise points.
        outer_cross_section_mcl = (cross_section_xyz_le[section_num + 1, :] + outer_cross_section_mcl_local_back
                                   + outer_cross_section_mcl_local_up)

        # Define number of spanwise points and panels.
        num_spanwise_panels = wing_cross_section.num_spanwise_panels
        num_spanwise_coordinates = num_spanwise_panels + 1

        # Get the spanwise coordinates.
        if wing_cross_section.spanwise_spacing == 'uniform':
            nondim_spanwise_coordinates = np.linspace(0, 1, num_spanwise_coordinates)
        elif wing_cross_section.spanwise_spacing == 'cosine':
            nondim_spanwise_coordinates = asmvp.geometry.cosspace(n_points=num_spanwise_coordinates)
        else:
            raise Exception("Bad value of section.spanwise_spacing!")

        # Make section_mcl_coordinates: M x N x 3 array of mean camberline coordinates. The first index is chordwise
        # point number, second index is spanwise point number, third is the x, y, and z coordinates. M is the number of
        # chordwise points. N is the number of spanwise points. Put a reversed version (from 1 to 0) of the non
        # dimensional spanwise coordinates in a row vector. This is size 1 x N, where N is the number of spanwise
        # points.
        reversed_nondim_spanwise_coordinates_row_vector = np.expand_dims((1 - nondim_spanwise_coordinates), 0)

        # Convert the reversed non dimensional spanwise coordinate row vector (from 1 to 0) to a matrix. This is size
        # 1 x N x 1, where N is the number of spanwise points.
        reversed_nondim_spanwise_coordinates_matrix = np.expand_dims(
            reversed_nondim_spanwise_coordinates_row_vector, 2)

        # Convert the inner and outer cross section's mean camberline coordinates column vectors to matrices. These are
        # size M x 1 x 3, where M is the number of chordwise points.
        inner_cross_section_mcl_matrix = np.expand_dims(inner_cross_section_mcl, 1)
        outer_cross_section_mcl_matrix = np.expand_dims(outer_cross_section_mcl, 1)

        # Put the non dimensional spanwise coordinates (from 0 to 1) in a row vector. This is size 1 x N, where N is the
        # number of spanwise points.
        nondim_spanwise_coordinates_row_vector = np.expand_dims(nondim_spanwise_coordinates, 0)

        # Convert the non dimensional spanwise coordinate row vector (from to 0 to 1) to a matrix. This is size
        # 1 x N x 1, where N is the number of spanwise points.
        nondim_spanwise_coordinates_matrix = np.expand_dims(nondim_spanwise_coordinates_row_vector, 2)

        # Linearly interpolate between inner and outer cross sections.
        #   This uses the following equation:
        #       f(a, b, i) = i * a + (1 - i) * b
        #       "a" is an N x 3 ndarray of the coordinates points along the outer cross section's mean
        #           camber line.
        #       "b" is an N x 3 ndarray of the coordinates of points along the inner cross section's mean
        #           camber line.
        #       "i" is a 1D array (or vector) of length M that holds the nondimensionalized spanwise panel
        #           spacing from 0 to 1.
        #       This produces a M x N x 3 ndarray where each slot holds the coordinates of a point on the
        #           surface between the inner and outer cross sections.
        section_mcl_vertices = (reversed_nondim_spanwise_coordinates_matrix * inner_cross_section_mcl_matrix
                                + nondim_spanwise_coordinates_matrix * outer_cross_section_mcl_matrix)

        # Compute the corners of each panel.
        front_inner_vertices = section_mcl_vertices[:-1, :-1, :]
        front_outer_vertices = section_mcl_vertices[:-1, 1:, :]
        back_inner_vertices = section_mcl_vertices[1:, :-1, :]
        back_outer_vertices = section_mcl_vertices[1:, 1:, :]

        # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels. The values are
        # either 1 if the panel at that location is a trailing edge, or 0 if not.
        section_is_trailing_edge = np.vstack((
            np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
            np.ones((1, num_spanwise_panels), dtype=bool)
        ))

        # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels. The values are
        # either 1 if the panel at that location is a leading edge, or 0 if not.
        section_is_leading_edge = np.vstack((
            np.ones((1, num_spanwise_panels), dtype=bool),
            np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool)
        ))

        # Initialize an empty ndarray to hold this sections. The matrix is size M x N, where M and N are the number
        # of chordwise and spanwise panels.
        section_panels = np.empty((num_chordwise_panels, num_spanwise_panels), dtype=object)

        # Loop through the empty section panels matrix and create a new panel object in each slot.
        for chordwise_position in range(num_chordwise_panels):
            for spanwise_position in range(num_spanwise_panels):
                section_panels[chordwise_position, spanwise_position] = asmvp.geometry.Panel(
                    front_left_vertex=front_inner_vertices[chordwise_position, spanwise_position],
                    front_right_vertex=front_outer_vertices[chordwise_position, spanwise_position],
                    back_left_vertex=back_inner_vertices[chordwise_position, spanwise_position],
                    back_right_vertex=back_outer_vertices[chordwise_position, spanwise_position],
                    is_trailing_edge=section_is_trailing_edge[chordwise_position, spanwise_position],
                    is_leading_edge=section_is_leading_edge[chordwise_position, spanwise_position])

        # This section's panel matrix is stack horizontally, to the right of the wing's panel matrix.
        panels = np.hstack((panels, section_panels))

        # Handle symmetry.
        if wing.symmetric:
            # Define the airfoils at each cross section.
            inner_airfoil = inner_cross_section.airfoil.add_control_surface(
                deflection=-inner_cross_section.control_surface_deflection,
                hinge_point=inner_cross_section.control_surface_hinge_point
            )
            outer_airfoil = outer_cross_section.airfoil.add_control_surface(
                deflection=-inner_cross_section.control_surface_deflection,
                # The inner cross section dictates control surface deflections.
                hinge_point=inner_cross_section.control_surface_hinge_point
            )

            # Make the mean camber lines for each cross section. First index is point number, second index is xyz.
            inner_cross_section_mcl_nondim = inner_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)
            outer_cross_section_mcl_nondim = outer_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)

            # Put the inner cross section's local up airfoil frame coordinates in a column vector.
            inner_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(inner_cross_section_mcl_nondim[:, 1],
                                                                                   1)

            # Put the inner cross section's local back airfoil frame coordinates in a column vector.
            inner_cross_section_mcl_nondim_local_back_column_vector = np.expand_dims(
                inner_cross_section_mcl_nondim[:, 0], 1)

            # Put the outer cross section's local up airfoil frame coordinates in a column vector.
            outer_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(outer_cross_section_mcl_nondim[:, 1],
                                                                                   1)

            # Put the outer cross section's local back airfoil frame coordinates in a column vector.
            outer_cross_section_mcl_nondim_local_back_column_vector = np.expand_dims(
                outer_cross_section_mcl_nondim[:, 0], 1)

            # Convert the inner cross section's non dimensional local back airfoil frame coordinates to meshed wing
            # coordinates.
            inner_cross_section_mcl_local_back = (cross_section_local_back_norm[section_num, :] *
                                                  inner_cross_section_mcl_nondim_local_back_column_vector *
                                                  cross_section_chord_length_list[section_num])

            # Convert the inner cross section's non dimensional local up airfoil frame coordinates to meshed wing
            # coordinates.
            inner_cross_section_mcl_local_up = (cross_section_local_up[section_num, :] *
                                                inner_cross_section_mcl_nondim_local_up_column_vector *
                                                cross_section_chord_length_list[section_num] *
                                                cross_section_scaling_factor[section_num])

            # Convert the outer cross section's non dimensional local back airfoil frame coordinates to meshed wing
            # coordinates.
            outer_cross_section_mcl_local_back = (cross_section_local_back_norm[section_num + 1, :] *
                                                  outer_cross_section_mcl_nondim_local_back_column_vector *
                                                  cross_section_chord_length_list[section_num + 1])

            # Convert the outer cross section's non dimensional local up airfoil frame coordinates to meshed wing
            # coordinates.
            outer_cross_section_mcl_local_up = (cross_section_local_up[section_num + 1, :] *
                                                outer_cross_section_mcl_nondim_local_up_column_vector *
                                                cross_section_chord_length_list[section_num + 1] *
                                                cross_section_scaling_factor[section_num + 1])

            # Convert the inner cross section's meshed wing coordinates to absolute coordinates. This is size M x 3,
            # where M is the number of chordwise points.
            inner_cross_section_mcl = (cross_section_xyz_le[section_num, :] + inner_cross_section_mcl_local_back
                                       + inner_cross_section_mcl_local_up)

            # Convert the outer cross section's meshed wing coordinates to absolute coordinates. This is size M x 3,
            # where M is the number of chordwise points.
            outer_cross_section_mcl = (cross_section_xyz_le[section_num + 1, :] + outer_cross_section_mcl_local_back
                                       + outer_cross_section_mcl_local_up)

            # Make section_mcl_coordinates: M x N x 3 array of mean camberline coordinates. First index is chordwise
            # point number, second index is spanwise point number, third are the x, y, and z coordinates. M is the
            # number of chordwise points. N is the number of spanwise points. Put a reversed version (from 1 to 0) of
            # the non dimensional spanwise coordinates in a row vector. This is size 1 x N, where N is the number of
            # spanwise points.
            reversed_nondim_spanwise_coordinates_row_vector = np.expand_dims((1 - nondim_spanwise_coordinates), 0)

            # Convert the reversed non dimensional spanwise coordinate row vector (from 1 to 0) to a matrix. This is
            # size 1 x N x 1, where N is the number of spanwise points.
            reversed_nondim_spanwise_coordinates_matrix = np.expand_dims(
                reversed_nondim_spanwise_coordinates_row_vector, 2)

            # Convert the inner and outer cross section's mean camberline coordinates column vectors to matrices. These
            # are size M x 1 x 3, where M is the number of chordwise points.
            inner_cross_section_mcl_matrix = np.expand_dims(inner_cross_section_mcl, 1)
            outer_cross_section_mcl_matrix = np.expand_dims(outer_cross_section_mcl, 1)

            # Put the non dimensional spanwise coordinates (from 0 to 1) in a row vector. This is size 1 x N, where N is
            # the number of spanwise points.
            nondim_spanwise_coordinates_row_vector = np.expand_dims(nondim_spanwise_coordinates, 0)

            # Convert the non dimensional spanwise coordinate row vector (from to 0 to 1) to a matrix. This is size
            # 1 x N x 1, where N is the number of spanwise points.
            nondim_spanwise_coordinates_matrix = np.expand_dims(nondim_spanwise_coordinates_row_vector, 2)

            # Linearly interpolate between inner and outer cross sections.
            #   This uses the following equation:
            #       f(a, b, i) = i * a + (1 - i) * b
            #       "a" is an N x 3 ndarray of the coordinates points along the outer cross section's mean
            #           camber line.
            #       "b" is an N x 3 ndarray of the coordinates of points along the inner cross section's mean
            #           camber line.
            #       "i" is a 1D array (or vector) of length M that holds the nondimensionalized spanwise panel
            #           spacing from 0 to 1.
            #       This produces a M x N x 3 ndarray where each slot holds the coordinates of a point on the
            #           surface between the inner and outer cross sections.
            section_mcl_vertices = (reversed_nondim_spanwise_coordinates_matrix * inner_cross_section_mcl_matrix
                                    + nondim_spanwise_coordinates_matrix * outer_cross_section_mcl_matrix)

            # Compute the corners of each panel.
            front_inner_vertices = section_mcl_vertices[:-1, :-1, :]
            front_outer_vertices = section_mcl_vertices[:-1, 1:, :]
            back_inner_vertices = section_mcl_vertices[1:, :-1, :]
            back_outer_vertices = section_mcl_vertices[1:, 1:, :]

            # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels. The values
            # are either 1 if the panel at that location is a trailing edge, or 0 if not.
            section_is_trailing_edge = np.vstack((
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
                np.ones((1, num_spanwise_panels), dtype=bool)
            ))

            # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels. The values
            # are either 1 if the panel at that location is a leading edge, or 0 if not.
            section_is_leading_edge = np.vstack((
                np.ones((1, num_spanwise_panels), dtype=bool),
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool)
            ))

            # Initialize an empty ndarray to hold this sections. The matrix is size M x N, where M and N are the
            # number of chordwise and spanwise panels.
            section_panels = np.empty((num_chordwise_panels, num_spanwise_panels), dtype=object)

            # Loop through the empty section panels matrix and create a new panel object in each slot.
            for chordwise_position in range(num_chordwise_panels):
                for spanwise_position in range(num_spanwise_panels):
                    # Reflect the vertices to create the reflected wing for the symmetric case.
                    front_inner_vertices_reflected = asmvp.geometry.reflect_over_xz_plane(front_inner_vertices[
                                                                                              chordwise_position,
                                                                                              spanwise_position])
                    front_outer_vertices_reflected = asmvp.geometry.reflect_over_xz_plane(front_outer_vertices[
                                                                                              chordwise_position,
                                                                                              spanwise_position])
                    back_inner_vertices_reflected = asmvp.geometry.reflect_over_xz_plane(back_inner_vertices[
                                                                                              chordwise_position,
                                                                                              spanwise_position])
                    back_outer_vertices_reflected = asmvp.geometry.reflect_over_xz_plane(back_outer_vertices[
                                                                                              chordwise_position,
                                                                                              spanwise_position])

                    section_panels[chordwise_position, spanwise_position] = asmvp.geometry.Panel(
                        front_left_vertex=front_outer_vertices_reflected,
                        front_right_vertex=front_inner_vertices_reflected,
                        back_left_vertex=back_outer_vertices_reflected,
                        back_right_vertex=back_inner_vertices_reflected,
                        is_trailing_edge=section_is_trailing_edge[chordwise_position, spanwise_position],
                        is_leading_edge=section_is_leading_edge[chordwise_position, spanwise_position])

            # This section's panel matrix is stack horizontally, to the left of the wing's panel matrix.
            panels = np.hstack((np.flip(section_panels, axis=1), panels))

    # Populate the wing's panels attribute.
    wing.panels = panels
