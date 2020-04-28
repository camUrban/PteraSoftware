# ToDo: Properly document this module.
"""This module contains useful aerodynamics functions.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import aerosandbox_legacy_v0 as asl
import numpy as np


# Adapted from:         vlm3.make_panels() in AeroSandbox
# Author:               Peter Sharpe
# Date of Retrieval:    03/28/2020
def make_panels(aerodynamics_problem):
    """This function takes in an object of the aerodynamics_problem class, creates a quadrilateral mesh of its geometry,
    and then populates the objects variables with the mesh data.

    Specifically, this function:
        makes structured lists of the panels' vertices, MxNx3 ndarrays describing a structured quadrilateral mesh of
            the wing's mean camber surface. For reference: The first index is the chordwise panel number, the second
            index is the spanwise panel number, and the third index contains the xyz coordinates,
        makes a structured list of the panels' normal directions, a MxNx3 ndarray containing a unit vector denoting the
            normal direction of the mean camber surface at the collocation point. For reference: The first index is the
            chordwise panel number, the second index is the spanwise panel number, and the third index contains the xyz
            coordinates,
        and takes into account control surface deflections.

    Notes:
        Control surface handling:
            Control surfaces are implemented into normal directions as intended.
        Symmetry handling:
            All symmetric wings have been split into separate halves.
            All wing halves have their spanwise vertices labeled from the left side of the airplane to the right.
            Control surface deflection symmetry has been handled; this is encoded into the normal directions.

    :param aerodynamics_problem: The problem to be analyzed.
    :return:
    """

    if aerodynamics_problem.verbose:
        print("Meshing...")

    # Initialize variables that will hold geometry information. These are empty ndarrays with 3 columns and 0 rows.
    collocation_points = np.empty((0, 3))
    normal_directions = np.empty((0, 3))
    front_left_vortex_vertices = np.empty((0, 3))
    front_right_vortex_vertices = np.empty((0, 3))
    back_left_vortex_vertices = np.empty((0, 3))
    back_right_vortex_vertices = np.empty((0, 3))
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

    # Initialize the list of cross sections to None.
    xsec = None

    # Iterate through the wings defined in this airplane.
    for wing_num in range(len(aerodynamics_problem.airplane.wings)):
        # For each wing, we want (where M is the number of chordwise panels, N is the number of spanwise panels)
        # the panel_vertices_structured_list: (M + 1) x (N + 1) x 3; corners of every panel, and
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
        #   Project quarter chords onto YZ plane and normalize.
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

        # Then, construct the normal directions for each xsec.
        # Make normals for the inner cross_sections, where we need to merge directions.
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
            #   chord.
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

        # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to compensate).
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
                # The inner cross section dictates control surface deflections.
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
            #   First index is chordwise point number, second index is spanwise point number, third is xyz.
            #   M is the number of chordwise points.
            #   N is the number of spanwise points.
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
            section_mcl_vertices = (reversed_nondim_spanwise_coordinates_matrix * inner_xsec_mcl_matrix
                                    + nondim_spanwise_coordinates_matrix * outer_xsec_mcl_matrix)

            # Compute the corners of each panel.
            front_inner_vertices = section_mcl_vertices[:-1, :-1, :]
            front_outer_vertices = section_mcl_vertices[:-1, 1:, :]
            back_inner_vertices = section_mcl_vertices[1:, :-1, :]
            back_outer_vertices = section_mcl_vertices[1:, 1:, :]

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
            front_inner_vertices = np.reshape(front_inner_vertices, newshape=(-1, 3), order='F')
            front_outer_vertices = np.reshape(front_outer_vertices, newshape=(-1, 3), order='F')
            back_inner_vertices = np.reshape(back_inner_vertices, newshape=(-1, 3), order='F')
            back_outer_vertices = np.reshape(back_outer_vertices, newshape=(-1, 3), order='F')

            # Reshape the M x N x 1 leading edge identifying matrices to be 1D arrays of shape (M * N) x 1.
            section_is_trailing_edge = np.reshape(section_is_trailing_edge, newshape=(-1), order='F')
            section_is_leading_edge = np.reshape(section_is_leading_edge, newshape=(-1), order='F')

            # Calculate the vortex vertices positions by performing vector math on the panel vertex locations.
            front_inner_vortex_vertices_to_add = front_inner_vertices + 0.25 * (back_inner_vertices -
                                                                                front_inner_vertices)
            front_outer_vortex_vertices_to_add = front_outer_vertices + 0.25 * (back_outer_vertices -
                                                                                front_outer_vertices)
            back_inner_vortex_vertices_to_add = front_inner_vortex_vertices_to_add + (back_inner_vertices -
                                                                                      front_inner_vertices)
            back_outer_vortex_vertices_to_add = front_outer_vortex_vertices_to_add + (back_outer_vertices -
                                                                                      front_outer_vertices)
            collocation_points_to_add = front_inner_vortex_vertices_to_add + 0.5 * (back_outer_vortex_vertices_to_add -
                                                                                    front_inner_vortex_vertices_to_add)

            # Calculate vortex ring normals and areas via diagonals.
            vortex_ring_first_diagonal = front_outer_vortex_vertices_to_add - back_inner_vortex_vertices_to_add
            vortex_ring_second_diagonal = front_inner_vortex_vertices_to_add - back_outer_vortex_vertices_to_add
            vortex_ring_cross_product = np.cross(vortex_ring_first_diagonal, vortex_ring_second_diagonal, axis=1)
            vortex_ring_cross_product_magnitude = np.linalg.norm(vortex_ring_cross_product, axis=1)
            normal_directions_to_add = (vortex_ring_cross_product
                                        / np.expand_dims(vortex_ring_cross_product_magnitude, axis=1))

            # Calculate panel normals and areas via diagonals.
            panel_first_diagonal = front_outer_vertices - back_inner_vertices
            panel_second_diagonal = front_inner_vertices - back_outer_vertices
            panel_cross_product = np.cross(panel_first_diagonal, panel_second_diagonal, axis=1)
            panel_cross_product_magnitude = np.linalg.norm(panel_cross_product, axis=1)
            areas_to_add = panel_cross_product_magnitude / 2

            # Append to the lists of mesh data. The data for each new section is stacked vertically.
            front_left_vertices = np.vstack((
                front_left_vertices,
                front_inner_vertices
            ))
            front_right_vertices = np.vstack((
                front_right_vertices,
                front_outer_vertices
            ))
            back_left_vertices = np.vstack((
                back_left_vertices,
                back_inner_vertices
            ))
            back_right_vertices = np.vstack((
                back_right_vertices,
                back_outer_vertices
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
            normal_directions = np.vstack((
                normal_directions,
                normal_directions_to_add
            ))
            front_left_vortex_vertices = np.vstack((
                front_left_vortex_vertices,
                front_inner_vortex_vertices_to_add
            ))
            front_right_vortex_vertices = np.vstack((
                front_right_vortex_vertices,
                front_outer_vortex_vertices_to_add
            ))
            back_left_vortex_vertices = np.vstack((
                back_left_vortex_vertices,
                back_inner_vortex_vertices_to_add
            ))
            back_right_vortex_vertices = np.vstack((
                back_right_vortex_vertices,
                back_outer_vortex_vertices_to_add
            ))
            collocation_points = np.vstack((
                collocation_points,
                collocation_points_to_add
            ))

            # Handle symmetry
            if wing.symmetric:
                # Define the airfoils at each cross section.
                inner_airfoil = inner_xsec.airfoil.add_control_surface(
                    deflection=-inner_xsec.control_surface_deflection,
                    hinge_point=inner_xsec.control_surface_hinge_point
                )
                outer_airfoil = outer_xsec.airfoil.add_control_surface(
                    deflection=-inner_xsec.control_surface_deflection,
                    # The inner cross section dictates control surface deflections.
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

                # Make section_mcl_coordinates: M x N x 3 array of mean camberline coordinates.
                #   First index is chordwise point number, second index is spanwise point number, third is xyz.
                #   M is the number of chordwise points.
                #   N is the number of spanwise points
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
                section_mcl_vertices = (reversed_nondim_spanwise_coordinates_matrix * inner_xsec_mcl_matrix
                                        + nondim_spanwise_coordinates_matrix * outer_xsec_mcl_matrix)

                # Compute the corners of each panel.
                front_inner_vertices = section_mcl_vertices[:-1, :-1, :]
                front_outer_vertices = section_mcl_vertices[:-1, 1:, :]
                back_inner_vertices = section_mcl_vertices[1:, :-1, :]
                back_outer_vertices = section_mcl_vertices[1:, 1:, :]

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
                front_inner_vertices = np.reshape(front_inner_vertices, newshape=(-1, 3), order='F')
                front_outer_vertices = np.reshape(front_outer_vertices, newshape=(-1, 3), order='F')
                back_inner_vertices = np.reshape(back_inner_vertices, newshape=(-1, 3), order='F')
                back_outer_vertices = np.reshape(back_outer_vertices, newshape=(-1, 3), order='F')

                # Reshape the M x N x 1 leading edge identifying matrices to be 1D arrays of shape (M * N) x 1.
                section_is_trailing_edge = np.reshape(section_is_trailing_edge, newshape=(-1), order='F')
                section_is_leading_edge = np.reshape(section_is_leading_edge, newshape=(-1), order='F')

                # Calculate the vortex vertices positions by performing vector math on the panel vertex locations.
                front_inner_vortex_vertices_to_add = front_inner_vertices + 0.25 * (back_inner_vertices -
                                                                                    front_inner_vertices)
                front_outer_vortex_vertices_to_add = front_outer_vertices + 0.25 * (back_outer_vertices -
                                                                                    front_outer_vertices)
                back_inner_vortex_vertices_to_add = front_inner_vortex_vertices_to_add + (back_inner_vertices -
                                                                                          front_inner_vertices)
                back_outer_vortex_vertices_to_add = front_outer_vortex_vertices_to_add + (back_outer_vertices -
                                                                                          front_outer_vertices)
                collocation_points_to_add = front_inner_vortex_vertices_to_add + 0.5 * (
                        back_outer_vortex_vertices_to_add -
                        front_inner_vortex_vertices_to_add)

                # Calculate vortex ring normals and areas via diagonals.
                vortex_ring_first_diagonal = front_outer_vortex_vertices_to_add - back_inner_vortex_vertices_to_add
                vortex_ring_second_diagonal = front_inner_vortex_vertices_to_add - back_outer_vortex_vertices_to_add
                vortex_ring_cross_product = np.cross(vortex_ring_first_diagonal, vortex_ring_second_diagonal, axis=1)
                vortex_ring_cross_product_magnitude = np.linalg.norm(vortex_ring_cross_product, axis=1)
                normal_directions_to_add = (vortex_ring_cross_product
                                            / np.expand_dims(vortex_ring_cross_product_magnitude, axis=1))

                # Calculate panel normals and areas via diagonals.
                panel_first_diagonal = front_outer_vertices - back_inner_vertices
                panel_second_diagonal = front_inner_vertices - back_outer_vertices
                panel_cross_product = np.cross(panel_first_diagonal, panel_second_diagonal, axis=1)
                panel_cross_product_magnitude = np.linalg.norm(panel_cross_product, axis=1)
                areas_to_add = panel_cross_product_magnitude / 2

            # Append to the lists of mesh data. The data for each new section is stacked vertically.
            front_left_vertices = np.vstack((
                front_left_vertices,
                asl.reflect_over_XZ_plane(front_outer_vertices)
            ))
            front_right_vertices = np.vstack((
                front_right_vertices,
                asl.reflect_over_XZ_plane(front_inner_vertices)
            ))
            back_left_vertices = np.vstack((
                back_left_vertices,
                asl.reflect_over_XZ_plane(back_outer_vertices)
            ))
            back_right_vertices = np.vstack((
                back_right_vertices,
                asl.reflect_over_XZ_plane(back_inner_vertices)
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
            normal_directions = np.vstack((
                normal_directions,
                asl.reflect_over_XZ_plane(normal_directions_to_add)
            ))
            front_left_vortex_vertices = np.vstack((
                front_left_vortex_vertices,
                asl.reflect_over_XZ_plane(front_outer_vortex_vertices_to_add)
            ))
            front_right_vortex_vertices = np.vstack((
                front_right_vortex_vertices,
                asl.reflect_over_XZ_plane(front_inner_vortex_vertices_to_add)
            ))
            back_left_vortex_vertices = np.vstack((
                back_left_vortex_vertices,
                asl.reflect_over_XZ_plane(back_outer_vortex_vertices_to_add)
            ))
            back_right_vortex_vertices = np.vstack((
                back_right_vortex_vertices,
                asl.reflect_over_XZ_plane(back_inner_vortex_vertices_to_add)
            ))
            collocation_points = np.vstack((
                collocation_points,
                asl.reflect_over_XZ_plane(collocation_points_to_add)
            ))

    # Update the panel and vortex vertex location variables in the aerodynamics_problem object.
    aerodynamics_problem.front_left_vertices = front_left_vertices
    aerodynamics_problem.front_right_vertices = front_right_vertices
    aerodynamics_problem.back_left_vertices = back_left_vertices
    aerodynamics_problem.back_right_vertices = back_right_vertices
    aerodynamics_problem.front_left_vortex_vertices = front_left_vortex_vertices
    aerodynamics_problem.front_right_vortex_vertices = front_right_vortex_vertices
    aerodynamics_problem.back_left_vortex_vertices = back_left_vortex_vertices
    aerodynamics_problem.back_right_vortex_vertices = back_right_vortex_vertices

    # Calculate the legs of each vortex ring. Update these variables in the aerodynamics_problem object.
    aerodynamics_problem.front_vortex_legs = (front_left_vortex_vertices - front_right_vortex_vertices)
    aerodynamics_problem.left_vortex_legs = (back_left_vortex_vertices - front_left_vortex_vertices)
    aerodynamics_problem.back_vortex_legs = (back_right_vortex_vertices - back_left_vortex_vertices)
    aerodynamics_problem.right_vortex_legs = (front_right_vortex_vertices - back_right_vortex_vertices)

    # Calculate the center of each leg in each vortex ring. Update these variables in the aerodynamics_problem object.
    aerodynamics_problem.front_vortex_leg_centers = (front_right_vortex_vertices + 0.5
                                                     * aerodynamics_problem.front_vortex_legs)
    aerodynamics_problem.left_vortex_leg_centers = (front_left_vortex_vertices + 0.5
                                                    * aerodynamics_problem.left_vortex_legs)
    aerodynamics_problem.back_vortex_leg_centers = (back_left_vortex_vertices + 0.5
                                                    * aerodynamics_problem.back_vortex_legs)
    aerodynamics_problem.right_vortex_leg_centers = (back_right_vortex_vertices + 0.5
                                                     * aerodynamics_problem.right_vortex_legs)

    # Update the following miscellaneous variables in the aerodynamics_problem object.
    aerodynamics_problem.areas = areas
    aerodynamics_problem.is_trailing_edge = is_trailing_edge
    aerodynamics_problem.is_leading_edge = is_leading_edge
    aerodynamics_problem.collocation_points = collocation_points
    aerodynamics_problem.normal_directions = normal_directions
    aerodynamics_problem.n_panels = len(aerodynamics_problem.collocation_points)
    aerodynamics_problem.panel_centers = front_left_vortex_vertices + 0.5 * (back_right_vertices - front_left_vertices)

    if aerodynamics_problem.verbose:
        print("Meshing complete!")


# ToDo: Properly document this function.
def mesh_problem(aerodynamics_problem):
    """

    :param aerodynamics_problem:
    :return:
    """
    airplane = aerodynamics_problem.airplane
    mesh_airplane(airplane)


# ToDo: Properly document this function.
def mesh_airplane(airplane):
    """

    :param airplane:
    :return:
    """

    for wing_num in range(len(airplane.wings)):
        wing = airplane.wings[wing_num]
        mesh_wing(wing)


# Adapted from:         vlm3.make_panels() in AeroSandbox
# Author:               Peter Sharpe
# Date of Retrieval:    03/28/2020
# ToDo: Properly document this function.
def mesh_wing(wing):
    """

    :param wing: The wing object to be meshed.
    :return:
    """

    # Define number of chordwise panels.
    num_chordwise_panels = wing.chordwise_panels
    # Define number of chordwise points.
    num_chordwise_coordinates = num_chordwise_panels + 1

    # Initialize variables that will hold geometry information. These are empty ndarrays of shape M x 0 x 3, where M is
    # the number of chordwise panels.
    collocation_points = np.empty((num_chordwise_panels, 0, 3))
    normal_directions = np.empty((num_chordwise_panels, 0, 3))
    front_left_vortex_vertices = np.empty((num_chordwise_panels, 0, 3))
    front_right_vortex_vertices = np.empty((num_chordwise_panels, 0, 3))
    back_left_vortex_vertices = np.empty((num_chordwise_panels, 0, 3))
    back_right_vortex_vertices = np.empty((num_chordwise_panels, 0, 3))
    front_left_vertices = np.empty((num_chordwise_panels, 0, 3))
    front_right_vertices = np.empty((num_chordwise_panels, 0, 3))
    back_left_vertices = np.empty((num_chordwise_panels, 0, 3))
    back_right_vertices = np.empty((num_chordwise_panels, 0, 3))

    # Initialize the areas variable which will hold the area of each panel.
    # It is an empty ndarray of shape M x 0, where M is the number of chordwise panels.
    areas = np.empty((num_chordwise_panels, 0))

    # Initialize the is_trailing_edge and is_leading_edge identifier variables.
    # They are empty ndarrays of shape M x 0, where M is the number of chordwise panels.
    is_trailing_edge = np.empty((num_chordwise_panels, 0), dtype=bool)
    is_leading_edge = np.empty((num_chordwise_panels, 0), dtype=bool)

    # Initialize the list of cross sections to None.
    xsec = None

    # Get the chordwise coordinates.
    if wing.chordwise_spacing == 'uniform':
        nondim_chordwise_coordinates = np.linspace(0, 1, num_chordwise_coordinates)
    elif wing.chordwise_spacing == 'cosine':
        nondim_chordwise_coordinates = asl.cosspace(0, 1, num_chordwise_coordinates)
    else:
        raise Exception("Bad value of wing.chordwise_spacing!")

    # Initialize two empty 0 x 3 ndarrays to hold the corners of each cross section. They will eventually be
    # L x 3 ndarrays, where L is number of cross sections.
    xsec_xyz_le = np.empty((0, 3))
    xsec_xyz_te = np.empty((0, 3))
    # Iterate through the meshed wing cross sections and vertically stack the global location each cross sections
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
    #   Project quarter chords onto YZ plane and normalize.
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

    # Then, construct the normal directions for each xsec.
    # Make normals for the inner cross_sections, where we need to merge directions.
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
        #   chord.
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

    # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to compensate).
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
            # The inner cross section dictates control surface deflections.
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

        # Convert the inner cross section's non dimensional local back airfoil frame coordinates to meshed wing
        # coordinates.
        inner_xsec_mcl_local_back = (xsec_local_back_norm[section_num, :] *
                                     inner_xsec_mcl_nondim_local_back_column_vector *
                                     xsec_chord_length_list[section_num])
        # Convert the inner cross section's non dimensional local up airfoil frame coordinates to meshed wing
        # coordinates.
        inner_xsec_mcl_local_up = (xsec_local_up[section_num, :] *
                                   inner_xsec_mcl_nondim_local_up_column_vector *
                                   xsec_chord_length_list[section_num] *
                                   xsec_scaling_factor[section_num])
        # Convert the outer cross section's non dimensional local back airfoil frame coordinates to meshed wing
        # coordinates.
        outer_xsec_mcl_local_back = (xsec_local_back_norm[section_num + 1, :] *
                                     outer_xsec_mcl_nondim_local_back_column_vector *
                                     xsec_chord_length_list[section_num + 1])
        # Convert the outer cross section's non dimensional local up airfoil frame coordinates to meshed wing
        # coordinates.
        outer_xsec_mcl_local_up = (xsec_local_up[section_num + 1, :] *
                                   outer_xsec_mcl_nondim_local_up_column_vector *
                                   xsec_chord_length_list[section_num + 1] *
                                   xsec_scaling_factor[section_num + 1])

        # Convert the inner cross section's meshed wing coordinates to absolute coordinates.
        # This is size M x 3, where M is the number of chordwise points.
        inner_xsec_mcl = xsec_xyz_le[section_num, :] + inner_xsec_mcl_local_back + inner_xsec_mcl_local_up
        # Convert the outer cross section's meshed wing coordinates to absolute coordinates.
        # This is size M x 3, where M is the number of chordwise points.
        outer_xsec_mcl = xsec_xyz_le[section_num + 1, :] + outer_xsec_mcl_local_back + outer_xsec_mcl_local_up

        # Define number of spanwise points.
        num_spanwise_coordinates = xsec.spanwise_panels + 1

        # Get the spanwise coordinates.
        if xsec.spanwise_spacing == 'uniform':
            nondim_spanwise_coordinates = np.linspace(0, 1, num_spanwise_coordinates)
        elif xsec.spanwise_spacing == 'cosine':
            nondim_spanwise_coordinates = asl.cosspace(n_points=num_spanwise_coordinates)
        else:
            raise Exception("Bad value of section.spanwise_spacing!")

        # Make section_mcl_coordinates: M x N x 3 array of mean camberline coordinates.
        #   First index is chordwise point number, second index is spanwise point number, third is xyz.
        #   M is the number of chordwise points.
        #   N is the number of spanwise points.
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
        section_mcl_vertices = (reversed_nondim_spanwise_coordinates_matrix * inner_xsec_mcl_matrix
                                + nondim_spanwise_coordinates_matrix * outer_xsec_mcl_matrix)

        # Compute the corners of each panel.
        front_inner_vertices = section_mcl_vertices[:-1, :-1, :]
        front_outer_vertices = section_mcl_vertices[:-1, 1:, :]
        back_inner_vertices = section_mcl_vertices[1:, :-1, :]
        back_outer_vertices = section_mcl_vertices[1:, 1:, :]

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

        # Calculate the vortex vertices positions by performing vector math on the panel vertex locations.
        front_inner_vortex_vertices_to_add = front_inner_vertices + 0.25 * (back_inner_vertices -
                                                                            front_inner_vertices)
        front_outer_vortex_vertices_to_add = front_outer_vertices + 0.25 * (back_outer_vertices -
                                                                            front_outer_vertices)
        back_inner_vortex_vertices_to_add = front_inner_vortex_vertices_to_add + (back_inner_vertices -
                                                                                  front_inner_vertices)
        back_outer_vortex_vertices_to_add = front_outer_vortex_vertices_to_add + (back_outer_vertices -
                                                                                  front_outer_vertices)
        collocation_points_to_add = front_inner_vortex_vertices_to_add + 0.5 * (back_outer_vortex_vertices_to_add -
                                                                                front_inner_vortex_vertices_to_add)

        # Calculate vortex ring normals and areas via diagonals.
        vortex_ring_first_diagonal = front_outer_vortex_vertices_to_add - back_inner_vortex_vertices_to_add
        vortex_ring_second_diagonal = front_inner_vortex_vertices_to_add - back_outer_vortex_vertices_to_add
        vortex_ring_cross_product = np.cross(vortex_ring_first_diagonal, vortex_ring_second_diagonal)
        vortex_ring_cross_product_magnitude = np.linalg.norm(vortex_ring_cross_product, axis=2)
        normal_directions_to_add = vortex_ring_cross_product / np.expand_dims(vortex_ring_cross_product_magnitude,
                                                                              axis=2)

        # Calculate panel normals and areas via diagonals.
        panel_first_diagonal = front_outer_vertices - back_inner_vertices
        panel_second_diagonal = front_inner_vertices - back_outer_vertices
        panel_cross_product = np.cross(panel_first_diagonal, panel_second_diagonal)
        panel_cross_product_magnitude = np.linalg.norm(panel_cross_product, axis=2)
        areas_to_add = panel_cross_product_magnitude / 2

        # Append to the lists of mesh data. The data for each new section is stacked vertically.
        front_left_vertices = np.hstack((
            front_left_vertices,
            front_inner_vertices
        ))
        front_right_vertices = np.hstack((
            front_right_vertices,
            front_outer_vertices
        ))
        back_left_vertices = np.hstack((
            back_left_vertices,
            back_inner_vertices
        ))
        back_right_vertices = np.hstack((
            back_right_vertices,
            back_outer_vertices
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
        normal_directions = np.hstack((
            normal_directions,
            normal_directions_to_add
        ))
        front_left_vortex_vertices = np.hstack((
            front_left_vortex_vertices,
            front_inner_vortex_vertices_to_add
        ))
        front_right_vortex_vertices = np.hstack((
            front_right_vortex_vertices,
            front_outer_vortex_vertices_to_add
        ))
        back_left_vortex_vertices = np.hstack((
            back_left_vortex_vertices,
            back_inner_vortex_vertices_to_add
        ))
        back_right_vortex_vertices = np.hstack((
            back_right_vortex_vertices,
            back_outer_vortex_vertices_to_add
        ))
        collocation_points = np.hstack((
            collocation_points,
            collocation_points_to_add
        ))

        # Handle symmetry
        if wing.symmetric:
            # Define the airfoils at each cross section.
            inner_airfoil = inner_xsec.airfoil.add_control_surface(
                deflection=-inner_xsec.control_surface_deflection,
                hinge_point=inner_xsec.control_surface_hinge_point
            )
            outer_airfoil = outer_xsec.airfoil.add_control_surface(
                deflection=-inner_xsec.control_surface_deflection,
                # The inner cross section dictates control surface deflections.
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

            # Convert the inner cross section's non dimensional local back airfoil frame coordinates to meshed wing
            # coordinates.
            inner_xsec_mcl_local_back = (xsec_local_back_norm[section_num, :] *
                                         inner_xsec_mcl_nondim_local_back_column_vector *
                                         xsec_chord_length_list[section_num])
            # Convert the inner cross section's non dimensional local up airfoil frame coordinates to meshed wing
            # coordinates.
            inner_xsec_mcl_local_up = (xsec_local_up[section_num, :] *
                                       inner_xsec_mcl_nondim_local_up_column_vector *
                                       xsec_chord_length_list[section_num] *
                                       xsec_scaling_factor[section_num])
            # Convert the outer cross section's non dimensional local back airfoil frame coordinates to meshed wing
            # coordinates.
            outer_xsec_mcl_local_back = (xsec_local_back_norm[section_num + 1, :] *
                                         outer_xsec_mcl_nondim_local_back_column_vector *
                                         xsec_chord_length_list[section_num + 1])
            # Convert the outer cross section's non dimensional local up airfoil frame coordinates to meshed wing
            # coordinates.
            outer_xsec_mcl_local_up = (xsec_local_up[section_num + 1, :] *
                                       outer_xsec_mcl_nondim_local_up_column_vector *
                                       xsec_chord_length_list[section_num + 1] *
                                       xsec_scaling_factor[section_num + 1])

            # Convert the inner cross section's meshed wing coordinates to absolute coordinates.
            # This is size M x 3, where M is the number of chordwise points.
            inner_xsec_mcl = xsec_xyz_le[section_num, :] + inner_xsec_mcl_local_back + inner_xsec_mcl_local_up
            # Convert the outer cross section's meshed wing coordinates to absolute coordinates.
            # This is size M x 3, where M is the number of chordwise points.
            outer_xsec_mcl = xsec_xyz_le[section_num + 1, :] + outer_xsec_mcl_local_back + outer_xsec_mcl_local_up

            # Make section_mcl_coordinates: M x N x 3 array of mean camberline coordinates.
            #   First index is chordwise point number, second index is spanwise point number, third is xyz.
            #   M is the number of chordwise points.
            #   N is the number of spanwise points
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
            section_mcl_vertices = (reversed_nondim_spanwise_coordinates_matrix * inner_xsec_mcl_matrix
                                    + nondim_spanwise_coordinates_matrix * outer_xsec_mcl_matrix)

            # Compute the corners of each panel.
            front_inner_vertices = section_mcl_vertices[:-1, :-1, :]
            front_outer_vertices = section_mcl_vertices[:-1, 1:, :]
            back_inner_vertices = section_mcl_vertices[1:, :-1, :]
            back_outer_vertices = section_mcl_vertices[1:, 1:, :]

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

            # Reshape the M x N x 3 corner coordinate matrices to be of shape (M * N) x 3, where M and N are the
            #   number of chordwise and spanwise panels.
            # The panels are scanned moving down each chord first and then moving across the span.
            # ToDo: Delete the following code if it is unnecessary.
            # front_inner_vertices = np.reshape(front_inner_vertices, newshape=(-1, 3), order='F')
            # front_outer_vertices = np.reshape(front_outer_vertices, newshape=(-1, 3), order='F')
            # back_inner_vertices = np.reshape(back_inner_vertices, newshape=(-1, 3), order='F')
            # back_outer_vertices = np.reshape(back_outer_vertices, newshape=(-1, 3), order='F')

            # Reshape the M x N x 1 leading edge identifying matrices to be 1D arrays of shape (M * N) x 1.
            # ToDo: Delete the following code if it is unnecessary.
            # section_is_trailing_edge = np.reshape(section_is_trailing_edge, newshape=(-1), order='F')
            # section_is_leading_edge = np.reshape(section_is_leading_edge, newshape=(-1), order='F')

            # Calculate the vortex vertices positions by performing vector math on the panel vertex locations.
            front_inner_vortex_vertices_to_add = front_inner_vertices + 0.25 * (back_inner_vertices -
                                                                                front_inner_vertices)
            front_outer_vortex_vertices_to_add = front_outer_vertices + 0.25 * (back_outer_vertices -
                                                                                front_outer_vertices)
            back_inner_vortex_vertices_to_add = front_inner_vortex_vertices_to_add + (back_inner_vertices -
                                                                                      front_inner_vertices)
            back_outer_vortex_vertices_to_add = front_outer_vortex_vertices_to_add + (back_outer_vertices -
                                                                                      front_outer_vertices)
            collocation_points_to_add = front_inner_vortex_vertices_to_add + 0.5 * (
                    back_outer_vortex_vertices_to_add -
                    front_inner_vortex_vertices_to_add)

            # Calculate vortex ring normals and areas via diagonals.
            vortex_ring_first_diagonal = front_outer_vortex_vertices_to_add - back_inner_vortex_vertices_to_add
            vortex_ring_second_diagonal = front_inner_vortex_vertices_to_add - back_outer_vortex_vertices_to_add
            vortex_ring_cross_product = np.cross(vortex_ring_first_diagonal, vortex_ring_second_diagonal)
            vortex_ring_cross_product_magnitude = np.linalg.norm(vortex_ring_cross_product, axis=2)
            normal_directions_to_add = vortex_ring_cross_product / np.expand_dims(vortex_ring_cross_product_magnitude,
                                                                                  axis=2)

            # Calculate panel normals and areas via diagonals.
            panel_first_diagonal = front_outer_vertices - back_inner_vertices
            panel_second_diagonal = front_inner_vertices - back_outer_vertices
            panel_cross_product = np.cross(panel_first_diagonal, panel_second_diagonal)
            panel_cross_product_magnitude = np.linalg.norm(panel_cross_product, axis=2)
            areas_to_add = panel_cross_product_magnitude / 2

            # Append to the lists of mesh data. The data for each new section is stacked horizontally.
            front_left_vertices = np.hstack((
                asl.reflect_over_XZ_plane(front_outer_vertices),
                front_left_vertices
            ))
            front_right_vertices = np.hstack((
                asl.reflect_over_XZ_plane(front_inner_vertices),
                front_right_vertices
            ))
            back_left_vertices = np.hstack((
                asl.reflect_over_XZ_plane(back_outer_vertices),
                back_left_vertices
            ))
            back_right_vertices = np.hstack((
                asl.reflect_over_XZ_plane(back_inner_vertices),
                back_right_vertices
            ))
            areas = np.hstack((
                areas_to_add,
                areas
            ))
            is_trailing_edge = np.hstack((
                section_is_trailing_edge,
                is_trailing_edge
            ))
            is_leading_edge = np.hstack((
                section_is_leading_edge,
                is_leading_edge
            ))
            normal_directions = np.hstack((
                asl.reflect_over_XZ_plane(normal_directions_to_add),
                normal_directions
            ))
            front_left_vortex_vertices = np.hstack((
                asl.reflect_over_XZ_plane(front_outer_vortex_vertices_to_add),
                front_left_vortex_vertices
            ))
            front_right_vortex_vertices = np.hstack((
                asl.reflect_over_XZ_plane(front_inner_vortex_vertices_to_add),
                front_right_vortex_vertices
            ))
            back_left_vortex_vertices = np.hstack((
                asl.reflect_over_XZ_plane(back_outer_vortex_vertices_to_add),
                back_left_vortex_vertices
            ))
            back_right_vortex_vertices = np.hstack((
                asl.reflect_over_XZ_plane(back_inner_vortex_vertices_to_add),
                back_right_vortex_vertices
            ))
            collocation_points = np.hstack((
                asl.reflect_over_XZ_plane(collocation_points_to_add),
                collocation_points
            ))

    # Update the panel and vortex vertex location variables in the aerodynamics_problem object.
    wing.front_left_vertices = front_left_vertices
    wing.front_right_vertices = front_right_vertices
    wing.back_left_vertices = back_left_vertices
    wing.back_right_vertices = back_right_vertices
    wing.front_left_vortex_vertices = front_left_vortex_vertices
    wing.front_right_vortex_vertices = front_right_vortex_vertices
    wing.back_left_vortex_vertices = back_left_vortex_vertices
    wing.back_right_vortex_vertices = back_right_vortex_vertices

    # Calculate the legs of each vortex ring. Update these variables in the aerodynamics_problem object.
    wing.front_vortex_legs = (front_left_vortex_vertices - front_right_vortex_vertices)
    wing.left_vortex_legs = (back_left_vortex_vertices - front_left_vortex_vertices)
    wing.back_vortex_legs = (back_right_vortex_vertices - back_left_vortex_vertices)
    wing.right_vortex_legs = (front_right_vortex_vertices - back_right_vortex_vertices)

    # Calculate the center of each leg in each vortex ring. Update these variables in the aerodynamics_problem object.
    wing.front_vortex_leg_centers = (front_right_vortex_vertices + 0.5
                                     * wing.front_vortex_legs)
    wing.left_vortex_leg_centers = (front_left_vortex_vertices + 0.5
                                    * wing.left_vortex_legs)
    wing.back_vortex_leg_centers = (back_left_vortex_vertices + 0.5
                                    * wing.back_vortex_legs)
    wing.right_vortex_leg_centers = (back_right_vortex_vertices + 0.5
                                     * wing.right_vortex_legs)

    # Update the following miscellaneous variables in the aerodynamics_problem object.
    wing.areas = areas
    wing.is_trailing_edge = is_trailing_edge
    wing.is_leading_edge = is_leading_edge
    wing.collocation_points = collocation_points
    wing.normal_directions = normal_directions
    wing.num_panels = len(np.reshape(wing.collocation_points, (-1, 3))[:, 0])
    wing.num_chordwise_panels = num_chordwise_panels
    wing.num_spanwise_panels = int(wing.num_panels / wing.num_chordwise_panels)
    wing.panel_centers = front_left_vortex_vertices + 0.5 * (back_right_vertices - front_left_vertices)

    # Initialize stuff.
    wing.forces_on_front_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.forces_on_back_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.forces_on_left_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.forces_on_right_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.total_force_on_panel_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.total_moment_on_panel_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))

    wing.normal_force_on_panels = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.pressure_on_panels = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.panels_delta_pressure_coefficient = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))


# ToDo: Properly document this method.
def move_panels(unsteady_aerodynamics_problem):
    """

    :param unsteady_aerodynamics_problem:
    :return:
    """
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
