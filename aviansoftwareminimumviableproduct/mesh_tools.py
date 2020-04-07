import aerosandbox_legacy_v0 as asl
import numpy as np


# Adapted from:         vlm3.make_panels() in AeroSandbox
# Author:               Peter Sharpe
# Date of Retrieval:    03/28/2020
def make_panels(aerodynamics_problem):
    if aerodynamics_problem.verbose:
        print("Meshing...")

    # Initialize variables that will hold geometry information. These are empty ndarrays with 3 columns and 0 rows.
    collocation_points = np.empty((0, 3))
    normal_directions = np.empty((0, 3))
    left_vortex_vertices = np.empty((0, 3))
    right_vortex_vertices = np.empty((0, 3))
    front_left_vertices = np.empty((0, 3))
    front_right_vertices = np.empty((0, 3))
    back_left_vertices = np.empty((0, 3))
    back_right_vertices = np.empty((0, 3))

    # Initialize the areas variable which will hold the area of each panel.
    # It is an empty ndarray with 0 rows and columns.
    areas = np.empty(0)

    # Initialize the is_trailing_edge and is_leading_edge identifier variables.
    # They are empty ndarrays of booleans with 0 rows and columns.
    is_trailing_edge = np.empty(0, dtype=bool)
    is_leading_edge = np.empty(0, dtype=bool)

    # Initialize the list of cross sections, xsecs, to None.
    xsec = None

    # Iterate through the wings defined in this airplane.
    for wing_num in range(len(aerodynamics_problem.airplane.wings)):
        # For each wing, we want (where M is the number of chordwise panels, N is the number of spanwise panels)
        # the panel_coordinates_structured_list: (M + 1) x (N + 1) x 3; corners of every panel, and
        # the normals_structured_list: M x N x 3; normal direction of each panel.

        # Get the wing.
        wing = aerodynamics_problem.airplane.wings[wing_num]

        # Define number of chordwise points.
        n_chordwise_coordinates = wing.chordwise_panels + 1

        # Get the chordwise coordinates.
        if wing.chordwise_spacing == 'uniform':
            nondim_chordwise_coordinates = np.linspace(0, 1, n_chordwise_coordinates)
        elif wing.chordwise_spacing == 'cosine':
            nondim_chordwise_coordinates = asl.cosspace(0, 1, n_chordwise_coordinates)
        else:
            raise Exception("Bad value of wing.chordwise_spacing!")

        # Initialize two empty 0 x 3 ndarrays to hold the corners of each cross section. They will eventually be
        # L x 3 ndarrays, where L is number of cross sections.
        xsec_xyz_le = np.empty((0, 3))
        xsec_xyz_te = np.empty((0, 3))
        # Iterate through the wing cross sections and vertically stack the global location each cross sections
        # leading and trailing edges. xsec.xyz_te is a method that returns the cross section's trailing edge's
        # coordinates.
        for xsec in wing.xsecs:
            xsec_xyz_le = np.vstack((xsec_xyz_le, xsec.xyz_le + wing.xyz_le))
            xsec_xyz_te = np.vstack((xsec_xyz_te, xsec.xyz_te() + wing.xyz_le))

        # Get the quarter chord vectors.
        # Get a L x 3 ndarray of points which are the quarter-chord points of cross section, where L is the number
        # of cross sections.
        xsec_xyz_quarter_chords = xsec_xyz_le + 0.25 * (xsec_xyz_te - xsec_xyz_le)
        # Get a (L - 1) x 3 ndarray of vectors connecting the cross section quarter chord points, where L is the
        # number of cross sections.
        section_quarter_chords = (xsec_xyz_quarter_chords[1:, :] - xsec_xyz_quarter_chords[:-1, :])

        # Get directions for transforming 2D airfoil data to 3D

        # Project quarter chords onto YZ plane and normalize.
        # Create a L x 2 ndarray with just the y and z components of the the section quarter chord vectors.
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

        # Then, construct the normal directions for each xsec.
        # Make normals for the inner xsecs, where we need to merge directions.
        if len(wing.xsecs) > 2:
            # Add together the adjacent normalized section quarter chords projected onto the the yz plane.
            xsec_local_normal_inners_non_norm = (section_quarter_chords_proj_yz_norm[:-1, :] +
                                                 section_quarter_chords_proj_yz_norm[1:, :])
            # Create a list of the magnitudes of the summed adjacent normalized section quarter chords projected
            # onto the yz plane.
            xsec_local_normal_inners_mag_list = np.linalg.norm(xsec_local_normal_inners_non_norm, axis=1)
            # Convert the list to a column vector.
            xsec_local_normal_inners_mag_column_vector = np.expand_dims(xsec_local_normal_inners_mag_list, axis=1)
            # Normalize the summed adjacent normalized section quarter chords projected onto the yz plane by
            # their magnitudes.
            xsec_local_normal_inners_norm = (xsec_local_normal_inners_non_norm /
                                             xsec_local_normal_inners_mag_column_vector)
            # Vertically stack the first normalized section quarter chord, the inner normalized section quarter
            # chords, and the last normalized section quarter chord.
            xsec_local_normal = np.vstack((
                section_quarter_chords_proj_yz_norm[0, :],
                xsec_local_normal_inners_norm,
                section_quarter_chords_proj_yz_norm[-1, :]
            ))
        else:
            # Vertically stack the first normalized section quarter chord, and the last normalized section quarter
            # chord.
            xsec_local_normal = np.vstack((
                section_quarter_chords_proj_yz_norm[0, :],
                section_quarter_chords_proj_yz_norm[-1, :]
            ))
        # xsec_local_normal is now a L x 3 array that represents the normal direction at each cross section.

        # Then, construct the back directions for each cross section.
        xsec_local_back_non_norm = xsec_xyz_te - xsec_xyz_le
        # Create a list of the cross section chord lengths.
        xsec_chord_length_list = np.linalg.norm(xsec_local_back_non_norm, axis=1)
        # Convert the list to a column vector.
        xsec_chord_length_column_vector = np.expand_dims(xsec_chord_length_list, axis=1)
        # Normalize the cross section back vectors by their magnitudes.
        xsec_local_back_norm = (xsec_local_back_non_norm / xsec_chord_length_column_vector)

        # Then, construct the up direction for each cross sections.
        xsec_local_up = np.cross(xsec_local_back_norm, xsec_local_normal, axis=1)

        # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to compensate)
        xsec_scaling_factor = 1 / np.sqrt(
            (1 + np.sum(section_quarter_chords_proj_yz_norm[1:, :] * section_quarter_chords_proj_yz_norm[:-1, :],
                        axis=1)
             ) / 2
        )
        xsec_scaling_factor = np.hstack((1, xsec_scaling_factor, 1))

        # Make the panels for each section.
        for section_num in range(len(wing.xsecs) - 1):
            # Define the relevant cross sections.
            inner_xsec = wing.xsecs[section_num]
            outer_xsec = wing.xsecs[section_num + 1]

            # Define the airfoils at each cross section.
            inner_airfoil = inner_xsec.airfoil.add_control_surface(
                deflection=inner_xsec.control_surface_deflection,
                hinge_point=inner_xsec.control_surface_hinge_point
            )
            outer_airfoil = outer_xsec.airfoil.add_control_surface(
                deflection=inner_xsec.control_surface_deflection,
                # Inner cross section dictates control surface deflections.
                hinge_point=inner_xsec.control_surface_hinge_point
            )

            # Make the mean camber lines for each cross section. First index is point number, second index is xyz.
            inner_xsec_mcl_nondim = inner_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)
            outer_xsec_mcl_nondim = outer_airfoil.get_downsampled_mcl(nondim_chordwise_coordinates)

            # Put the inner cross section's local up airfoil frame coordinates in a column vector.
            inner_xsec_mcl_nondim_local_up_column_vector = np.expand_dims(inner_xsec_mcl_nondim[:, 1], 1)
            # Put the inner cross section's local back airfoil frame coordinates in a column vector.
            inner_xsec_mcl_nondim_local_back_column_vector = np.expand_dims(inner_xsec_mcl_nondim[:, 0], 1)
            # Put the outer cross section's local up airfoil frame coordinates in a column vector.
            outer_xsec_mcl_nondim_local_up_column_vector = np.expand_dims(outer_xsec_mcl_nondim[:, 1], 1)
            # Put the outer cross section's local back airfoil frame coordinates in a column vector.
            outer_xsec_mcl_nondim_local_back_column_vector = np.expand_dims(outer_xsec_mcl_nondim[:, 0], 1)

            # Convert the inner cross section's non dimensional local back airfoil frame coordinates to wing
            # coordinates.
            inner_xsec_mcl_local_back = (xsec_local_back_norm[section_num, :] *
                                         inner_xsec_mcl_nondim_local_back_column_vector *
                                         xsec_chord_length_list[section_num])
            # Convert the inner cross section's non dimensional local up airfoil frame coordinates to wing
            # coordinates.
            inner_xsec_mcl_local_up = (xsec_local_up[section_num, :] *
                                       inner_xsec_mcl_nondim_local_up_column_vector *
                                       xsec_chord_length_list[section_num] *
                                       xsec_scaling_factor[section_num])
            # Convert the outer cross section's non dimensional local back airfoil frame coordinates to wing
            # coordinates.
            outer_xsec_mcl_local_back = (xsec_local_back_norm[section_num + 1, :] *
                                         outer_xsec_mcl_nondim_local_back_column_vector *
                                         xsec_chord_length_list[section_num + 1])
            # Convert the outer cross section's non dimensional local up airfoil frame coordinates to wing
            # coordinates.
            outer_xsec_mcl_local_up = (xsec_local_up[section_num + 1, :] *
                                       outer_xsec_mcl_nondim_local_up_column_vector *
                                       xsec_chord_length_list[section_num + 1] *
                                       xsec_scaling_factor[section_num + 1])

            # Convert the inner cross section's wing coordinates to absolute coordinates.
            # This is size M x 3, where M is the number of chordwise points.
            inner_xsec_mcl = xsec_xyz_le[section_num, :] + inner_xsec_mcl_local_back + inner_xsec_mcl_local_up
            # Convert the outer cross section's wing coordinates to absolute coordinates.
            # This is size M x 3, where M is the number of chordwise points.
            outer_xsec_mcl = xsec_xyz_le[section_num + 1, :] + outer_xsec_mcl_local_back + outer_xsec_mcl_local_up

            # Define number of spanwise points.
            n_spanwise_coordinates = xsec.spanwise_panels + 1

            # Get the spanwise coordinates.
            if xsec.spanwise_spacing == 'uniform':
                nondim_spanwise_coordinates = np.linspace(0, 1, n_spanwise_coordinates)
            elif xsec.spanwise_spacing == 'cosine':
                nondim_spanwise_coordinates = asl.cosspace(n_points=n_spanwise_coordinates)
            else:
                raise Exception("Bad value of section.spanwise_spacing!")

            # Make section_mcl_coordinates: M x N x 3 array of mean camberline coordinates.
            # First index is chordwise point number, second index is spanwise point number, third is xyz.
            # M is the number of chordwise points.
            # N is the number of spanwise points

            # Put a reversed version (from 1 to 0) of the non dimensional spanwise coordinates in a row vector.
            # This is size 1 x N, where N is the number of spanwise points.
            reversed_nondim_spanwise_coordinates_row_vector = np.expand_dims((1 - nondim_spanwise_coordinates), 0)
            # Convert the reversed non dimensional spanwise coordinate row vector (from 1 to 0) to a matrix.
            # This is size 1 x N x 1, where N is the number of spanwise points.
            reversed_nondim_spanwise_coordinates_matrix = np.expand_dims(
                reversed_nondim_spanwise_coordinates_row_vector, 2)

            # Convert the inner and outer cross section's mean camberline coordinates column vectors to matrices.
            # These are size M x 1 x 3, where M is the number of chordwise points.
            inner_xsec_mcl_matrix = np.expand_dims(inner_xsec_mcl, 1)
            outer_xsec_mcl_matrix = np.expand_dims(outer_xsec_mcl, 1)

            # Put the non dimensional spanwise coordinates (from 0 to 1) in a row vector.
            # This is size 1 x N, where N is the number of spanwise points.
            nondim_spanwise_coordinates_row_vector = np.expand_dims(nondim_spanwise_coordinates, 0)
            # Convert the non dimensional spanwise coordinate row vector (from to 0 to 1) to a matrix.
            # This is size 1 x N x 1, where N is the number of spanwise points.
            nondim_spanwise_coordinates_matrix = np.expand_dims(nondim_spanwise_coordinates_row_vector, 2)

            # ToDo: Write comments to describe the following two variables.
            # This is size M x N x 3, where M and N are the number of chordwise points and spanwise points.
            section_mcl_coordinates_from_inner_xsec = (reversed_nondim_spanwise_coordinates_matrix *
                                                       inner_xsec_mcl_matrix)
            section_mcl_coordinates_from_outer_xsec = (nondim_spanwise_coordinates_matrix *
                                                       outer_xsec_mcl_matrix)

            # ToDo: Ask Peter how the following operation produces the section_mcl_coordinates.
            # TODO Make this work for large twist angles.
            section_mcl_coordinates = (section_mcl_coordinates_from_inner_xsec +
                                       section_mcl_coordinates_from_outer_xsec)

            # Compute corners of each panel.
            front_inner_coordinates = section_mcl_coordinates[:-1, :-1, :]
            front_outer_coordinates = section_mcl_coordinates[:-1, 1:, :]
            back_inner_coordinates = section_mcl_coordinates[1:, :-1, :]
            back_outer_coordinates = section_mcl_coordinates[1:, 1:, :]

            # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels.
            # The values are either 1 if the panel at that location is a trailing edge, or 0 if not.
            section_is_trailing_edge = np.vstack((
                np.zeros((wing.chordwise_panels - 1, xsec.spanwise_panels), dtype=bool),
                np.ones((1, xsec.spanwise_panels), dtype=bool)
            ))
            # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels.
            # The values are either 1 if the panel at that location is a leading edge, or 0 if not.
            section_is_leading_edge = np.vstack((
                np.ones((1, xsec.spanwise_panels), dtype=bool),
                np.zeros((wing.chordwise_panels - 1, xsec.spanwise_panels), dtype=bool)
            ))

            # Reshape the (M - 1) x (N - 1) x 3 corner coordinate matrices to be of shape ((M - 1) * (N-1)) x 3.
            # The panels are scanned moving down each chord first and then moving across the span.
            front_inner_coordinates = np.reshape(front_inner_coordinates, newshape=(-1, 3), order='F')
            front_outer_coordinates = np.reshape(front_outer_coordinates, newshape=(-1, 3), order='F')
            back_inner_coordinates = np.reshape(back_inner_coordinates, newshape=(-1, 3), order='F')
            back_outer_coordinates = np.reshape(back_outer_coordinates, newshape=(-1, 3), order='F')

            # Reshape the M x N x 1 leading edge identifying matrices to be 1D arrays of shape (M * N) x 1.
            section_is_trailing_edge = np.reshape(section_is_trailing_edge, newshape=(-1), order='F')
            section_is_leading_edge = np.reshape(section_is_leading_edge, newshape=(-1), order='F')

            # Calculate panel normals and areas via diagonals
            diag1 = front_outer_coordinates - back_inner_coordinates
            diag2 = front_inner_coordinates - back_outer_coordinates
            diag_cross = np.cross(diag1, diag2, axis=1)
            diag_cross_norm = np.linalg.norm(diag_cross, axis=1)
            normals_to_add = diag_cross / np.expand_dims(diag_cross_norm, axis=1)
            areas_to_add = diag_cross_norm / 2

            # ToDo: Continue 1.2.1. from here onward.
            # Make the panel data
            collocations_to_add = (
                    0.5 * (0.25 * front_inner_coordinates + 0.75 * back_inner_coordinates) +
                    0.5 * (0.25 * front_outer_coordinates + 0.75 * back_outer_coordinates)
            )
            inner_vortex_vertices_to_add = 0.75 * front_inner_coordinates + 0.25 * back_inner_coordinates
            outer_vortex_vertices_to_add = 0.75 * front_outer_coordinates + 0.25 * back_outer_coordinates

            # Append to the lists of panel data. The data for each new section is stacked vertically
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
                            xsec_local_back_norm[section_num, :] * np.expand_dims(inner_xsec_mcl_nondim[:, 0], 1) *
                            xsec_chord_length_list[
                                section_num] +
                            xsec_local_up[section_num, :] * np.expand_dims(inner_xsec_mcl_nondim[:, 1], 1) *
                            xsec_chord_length_list[
                                section_num] * xsec_scaling_factor[section_num]
                    )

                    outer_xsec_mcl = xsec_xyz_le[section_num + 1, :] + (
                            xsec_local_back_norm[section_num + 1, :] * np.expand_dims(outer_xsec_mcl_nondim[:, 0],
                                                                                      1) *
                            xsec_chord_length_list[section_num + 1] +
                            xsec_local_up[section_num + 1, :] * np.expand_dims(outer_xsec_mcl_nondim[:, 1], 1) *
                            xsec_chord_length_list[
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
                    )

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

    # Write to aerodynamics_problem object
    aerodynamics_problem.front_left_vertices = front_left_vertices
    aerodynamics_problem.front_right_vertices = front_right_vertices
    aerodynamics_problem.back_left_vertices = back_left_vertices
    aerodynamics_problem.back_right_vertices = back_right_vertices

    aerodynamics_problem.front_wing_vortex_centers = (front_left_vertices + front_right_vertices) / 2
    aerodynamics_problem.back_wing_vortex_centers = (back_left_vertices + back_right_vertices) / 2
    aerodynamics_problem.left_wing_vortex_centers = (front_left_vertices + back_left_vertices) / 2
    aerodynamics_problem.right_wing_vortex_centers = (front_right_vertices + back_right_vertices) / 2

    aerodynamics_problem.front_wing_vortex_legs = (front_right_vertices - front_left_vertices)
    aerodynamics_problem.right_wing_vortex_legs = (back_right_vertices - front_right_vertices)
    aerodynamics_problem.back_wing_vortex_legs = (back_left_vertices - back_right_vertices)
    aerodynamics_problem.left_wing_vortex_legs = (front_left_vertices - back_left_vertices)

    aerodynamics_problem.areas = areas
    aerodynamics_problem.is_trailing_edge = is_trailing_edge
    aerodynamics_problem.is_leading_edge = is_leading_edge
    aerodynamics_problem.collocation_points = collocation_points
    aerodynamics_problem.normal_directions = normal_directions
    aerodynamics_problem.n_panels = len(aerodynamics_problem.collocation_points)

    if aerodynamics_problem.verbose:
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
    #   within the vlm2 class, accessible at aerodynamics_problem.panel_coordinates_structured_list and
    #   aerodynamics_problem.normals_structured_list, respectively.
    #   Control surface handling:
    #       Control surfaces are implemented into normal directions as intended.
    #   Symmetry handling:
    #       All symmetric wings have been split into separate halves.
    #       All wing halves have their spanwise coordinates labeled from the left side of the airplane to the right.
    #       Control surface deflection symmetry has been handled; this is encoded into the normal directions.


# ToDo: Properly comment code in the following method
def move_panels(unsteady_aerodynamics_problem):
    new_front_left_vertices = unsteady_aerodynamics_problem.initial_front_left_vertices
    new_front_right_vertices = unsteady_aerodynamics_problem.initial_front_right_vertices
    new_back_left_vertices = unsteady_aerodynamics_problem.initial_back_left_vertices
    new_back_right_vertices = unsteady_aerodynamics_problem.initial_back_right_vertices

    front_left_vertex_sweep_angles = []
    front_right_vertex_sweep_angles = []
    back_left_vertex_sweep_angles = []
    back_right_vertex_sweep_angles = []

    for vertex in range(unsteady_aerodynamics_problem.n_panels):
        if unsteady_aerodynamics_problem.initial_front_left_vertices[vertex, 1] < 0:
            front_left_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            front_left_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

        if unsteady_aerodynamics_problem.initial_front_right_vertices[vertex, 1] < 0:
            front_right_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            front_right_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

        if unsteady_aerodynamics_problem.initial_back_left_vertices[vertex, 1] < 0:
            back_left_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            back_left_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

        if unsteady_aerodynamics_problem.initial_front_right_vertices[vertex, 1] < 0:
            back_right_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            back_right_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

    front_left_radii = np.sqrt(
        unsteady_aerodynamics_problem.initial_front_left_vertices[:, 1] ** 2 +
        unsteady_aerodynamics_problem.initial_front_left_vertices[:, 2] ** 2)
    front_right_radii = np.sqrt(
        unsteady_aerodynamics_problem.initial_front_right_vertices[:, 1] ** 2 +
        unsteady_aerodynamics_problem.initial_front_right_vertices[:, 2] ** 2)
    back_left_radii = np.sqrt(
        unsteady_aerodynamics_problem.initial_back_left_vertices[:, 1] ** 2 +
        unsteady_aerodynamics_problem.initial_back_left_vertices[:, 2] ** 2)
    back_right_radii = np.sqrt(
        unsteady_aerodynamics_problem.initial_back_right_vertices[:, 1] ** 2 +
        unsteady_aerodynamics_problem.initial_back_right_vertices[:, 2] ** 2)

    new_front_left_vertices[:, 1] = front_left_radii * np.cos(front_left_vertex_sweep_angles)
    new_front_left_vertices[:, 2] = front_left_radii * np.sin(front_left_vertex_sweep_angles)

    new_front_right_vertices[:, 1] = front_right_radii * np.cos(front_right_vertex_sweep_angles)
    new_front_right_vertices[:, 2] = front_right_radii * np.sin(front_right_vertex_sweep_angles)

    new_back_left_vertices[:, 1] = back_left_radii * np.cos(back_left_vertex_sweep_angles)
    new_back_left_vertices[:, 2] = back_left_radii * np.sin(back_left_vertex_sweep_angles)

    new_back_right_vertices[:, 1] = back_right_radii * np.cos(back_right_vertex_sweep_angles)
    new_back_right_vertices[:, 2] = back_right_radii * np.sin(back_right_vertex_sweep_angles)

    unsteady_aerodynamics_problem.front_left_vertices = new_front_left_vertices
    unsteady_aerodynamics_problem.front_right_vertices = new_front_right_vertices
    unsteady_aerodynamics_problem.back_left_vertices = new_back_left_vertices
    unsteady_aerodynamics_problem.back_right_vertices = new_back_right_vertices

    unsteady_aerodynamics_problem.front_wing_vortex_centers = (new_front_left_vertices + new_front_right_vertices) / 2
    unsteady_aerodynamics_problem.back_wing_vortex_centers = (new_back_left_vertices + new_back_right_vertices) / 2
    unsteady_aerodynamics_problem.left_wing_vortex_centers = (new_front_left_vertices + new_back_left_vertices) / 2
    unsteady_aerodynamics_problem.right_wing_vortex_centers = (new_front_right_vertices + new_back_right_vertices) / 2

    unsteady_aerodynamics_problem.front_wing_vortex_legs = (new_front_right_vertices - new_front_left_vertices)
    unsteady_aerodynamics_problem.right_wing_vortex_legs = (new_back_right_vertices - new_front_right_vertices)
    unsteady_aerodynamics_problem.back_wing_vortex_legs = (new_back_left_vertices - new_back_right_vertices)
    unsteady_aerodynamics_problem.left_wing_vortex_legs = (new_front_left_vertices - new_back_left_vertices)
