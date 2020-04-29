"""This module contains useful functions for creating meshes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    mesh_problem: This function takes in an object of the AerodynamicsProblem class, and calls mesh_airplane on its
                  airplane object.
    mesh_airplane: This function takes in an object of the Airplane class and calls mesh_wing on each of its wings.
    mesh_wing: This function takes in an object of the Wing class and creates a quadrilateral mesh of its geometry, and
               then populates the objects variables with the mesh data.
    move_panels: This function takes in a problem of the UnsteadyAerodynamicsProblem class, and modifies it's panel
                 locations as it time steps through the simulation.
"""

import aviansoftwareminimumviableproduct as asmvp
import numpy as np


def mesh_problem(aerodynamics_problem):
    """This function takes in an object of the AerodynamicsProblem class, and calls mesh_airplane on its airplane
    object.

    :param aerodynamics_problem: AerodynamicsProblem
        This is the problem whose geometry is to be meshed.
    :return: None
    """

    # Initialize a variable to hold the problem's airplane object.
    airplane = aerodynamics_problem.airplane

    # Mesh the airplane.
    mesh_airplane(airplane)


def mesh_airplane(airplane):
    """This function takes in an object of the Airplane class and calls mesh_wing on each of its wings.

    :param airplane: Airplane.
        The airplane whose geometry is to be meshed.
    :return: None
    """

    # Iterate through the wings in the airplane object and mesh each wing.
    for wing_num in range(len(airplane.wings)):
        wing = airplane.wings[wing_num]
        mesh_wing(wing)


def mesh_wing(wing):
    """This function takes in an object of the Wing class and creates a quadrilateral mesh of its geometry, and then
    populates the objects variables with the mesh data.

    Citation:
        Adapted from:         vlm3.make_panels in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    03/28/2020

    :param wing: Wing
        The wing object to be meshed.
    :return: None
    """

    # Define the number of chordwise panels and points.
    num_chordwise_panels = wing.chordwise_panels
    num_chordwise_coordinates = num_chordwise_panels + 1

    # Initialize variables that will hold geometry information. These are empty numpy arrays of shape M x 0 x 3, where M
    # is the number of chordwise panels.
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

    # Initialize the areas variable which will hold the area of each panel. It is an empty ndarray of shape M x 0, where
    # M is the number of chordwise panels.
    areas = np.empty((num_chordwise_panels, 0))

    # Initialize the is_trailing_edge and is_leading_edge identifier variables. They are empty ndarrays of shape M x 0,
    # where M is the number of chordwise panels.
    is_trailing_edge = np.empty((num_chordwise_panels, 0), dtype=bool)
    is_leading_edge = np.empty((num_chordwise_panels, 0), dtype=bool)

    # Initialize the list of cross sections to None.
    cross_section = None

    # Get the chordwise coordinates.
    if wing.chordwise_spacing == 'uniform':
        nondim_chordwise_coordinates = np.linspace(0, 1, num_chordwise_coordinates)
    elif wing.chordwise_spacing == 'cosine':
        nondim_chordwise_coordinates = asmvp.geometry.cosspace(0, 1, num_chordwise_coordinates)
    else:
        raise Exception("Bad value of wing.chordwise_spacing!")

    # Initialize two empty 0 x 3 ndarrays to hold the corners of each cross section. They will eventually be L x 3
    # numpy arrays, where L is number of cross sections.
    cross_section_xyz_le = np.empty((0, 3))
    cross_section_xyz_te = np.empty((0, 3))

    # Iterate through the meshed wing cross sections and vertically stack the global location each cross sections
    # leading and trailing edges. cross_section.xyz_te is a method that returns the cross section's trailing edge's
    # coordinates.
    for cross_section in wing.cross_sections:
        cross_section_xyz_le = np.vstack((cross_section_xyz_le, cross_section.xyz_le + wing.xyz_le))
        cross_section_xyz_te = np.vstack((cross_section_xyz_te, cross_section.xyz_te() + wing.xyz_le))

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

    # Create a list of the magnitudes of each row of the section_quarter_chords_yz numpy array.
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

    # Then, construct the normal directions for each cross_section. Make the normals for the inner cross_sections, where
    # we need to merge directions.
    if len(wing.cross_sections) > 2:
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
    for section_num in range(len(wing.cross_sections) - 1):
        # Define the relevant cross sections.
        inner_cross_section = wing.cross_sections[section_num]
        outer_cross_section = wing.cross_sections[section_num + 1]

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

        # Define number of spanwise points.
        num_spanwise_coordinates = cross_section.spanwise_panels + 1

        # Get the spanwise coordinates.
        if cross_section.spanwise_spacing == 'uniform':
            nondim_spanwise_coordinates = np.linspace(0, 1, num_spanwise_coordinates)
        elif cross_section.spanwise_spacing == 'cosine':
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
            np.zeros((wing.chordwise_panels - 1, cross_section.spanwise_panels), dtype=bool),
            np.ones((1, cross_section.spanwise_panels), dtype=bool)
        ))

        # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels. The values are
        # either 1 if the panel at that location is a leading edge, or 0 if not.
        section_is_leading_edge = np.vstack((
            np.ones((1, cross_section.spanwise_panels), dtype=bool),
            np.zeros((wing.chordwise_panels - 1, cross_section.spanwise_panels), dtype=bool)
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
                np.zeros((wing.chordwise_panels - 1, cross_section.spanwise_panels), dtype=bool),
                np.ones((1, cross_section.spanwise_panels), dtype=bool)
            ))

            # Compute a matrix that is M x N, where M and N are the number of chordwise and spanwise panels. The values
            # are either 1 if the panel at that location is a leading edge, or 0 if not.
            section_is_leading_edge = np.vstack((
                np.ones((1, cross_section.spanwise_panels), dtype=bool),
                np.zeros((wing.chordwise_panels - 1, cross_section.spanwise_panels), dtype=bool)
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
                asmvp.geometry.reflect_over_xz_plane(front_outer_vertices),
                front_left_vertices
            ))
            front_right_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(front_inner_vertices),
                front_right_vertices
            ))
            back_left_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(back_outer_vertices),
                back_left_vertices
            ))
            back_right_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(back_inner_vertices),
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
                asmvp.geometry.reflect_over_xz_plane(normal_directions_to_add),
                normal_directions
            ))
            front_left_vortex_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(front_outer_vortex_vertices_to_add),
                front_left_vortex_vertices
            ))
            front_right_vortex_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(front_inner_vortex_vertices_to_add),
                front_right_vortex_vertices
            ))
            back_left_vortex_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(back_outer_vortex_vertices_to_add),
                back_left_vortex_vertices
            ))
            back_right_vortex_vertices = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(back_inner_vortex_vertices_to_add),
                back_right_vortex_vertices
            ))
            collocation_points = np.hstack((
                asmvp.geometry.reflect_over_xz_plane(collocation_points_to_add),
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

    # Initialize variables that will be filled in after the solution is found.
    wing.forces_on_front_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.forces_on_back_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.forces_on_left_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.forces_on_right_vortices_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.total_force_on_panel_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.total_moment_on_panel_in_geometry_axes = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.normal_force_on_panels = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.pressure_on_panels = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))
    wing.panels_delta_pressure_coefficient = np.zeros((wing.num_chordwise_panels, wing.num_spanwise_panels, 3))


def move_panels(unsteady_aerodynamics_problem):
    """This function takes in a problem of the UnsteadyAerodynamicsProblem class, and modifies it's panel locations as
    it time steps through the simulation.

    :param unsteady_aerodynamics_problem: UnsteadyAerodynamicsProblem
        This is the problem whose geometry is to be moved.
    :return: None
    """

    # Get the panel vertices from the unsteady aerodynamics problem.
    new_front_left_vertices = unsteady_aerodynamics_problem.initial_front_left_vertices
    new_front_right_vertices = unsteady_aerodynamics_problem.initial_front_right_vertices
    new_back_left_vertices = unsteady_aerodynamics_problem.initial_back_left_vertices
    new_back_right_vertices = unsteady_aerodynamics_problem.initial_back_right_vertices

    # Initialize empty lists to hold the sweep angles of each panel vertex.
    front_left_vertex_sweep_angles = []
    front_right_vertex_sweep_angles = []
    back_left_vertex_sweep_angles = []
    back_right_vertex_sweep_angles = []

    # Iterate through the unsteady aerodynamics problem's vertices.
    for vertex in range(unsteady_aerodynamics_problem.num_panels):
        # Check if the front left vertex's y coordinate is negative.
        if unsteady_aerodynamics_problem.initial_front_left_vertices[vertex, 1] < 0:
            # If so, modify the angle and append the angle to the list.
            front_left_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            # If not, append the angle to the list.
            front_left_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

        # Check if the front right vertex's y coordinate is negative.
        if unsteady_aerodynamics_problem.initial_front_right_vertices[vertex, 1] < 0:
            # If so, modify the angle and append the angle to the list.
            front_right_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            # If not, append the angle to the list.
            front_right_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

        # Check if the back left vertex's y coordinate is negative.
        if unsteady_aerodynamics_problem.initial_back_left_vertices[vertex, 1] < 0:
            # If so, modify the angle and append the angle to the list.
            back_left_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            # If not, append the angle to the list.
            back_left_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

        # Check if the back right vertex's y coordinate is negative.
        if unsteady_aerodynamics_problem.initial_front_right_vertices[vertex, 1] < 0:
            # If so, modify the angle and append the angle to the list.
            back_right_vertex_sweep_angles.append(np.pi - unsteady_aerodynamics_problem.current_sweep_angle)
        else:
            # If not, append the angle to the list.
            back_right_vertex_sweep_angles.append(unsteady_aerodynamics_problem.current_sweep_angle)

    # Find the distances from the origin to each panel vertex.
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

    # Rotate the coordinates of the front left vertex.
    new_front_left_vertices[:, 1] = front_left_radii * np.cos(front_left_vertex_sweep_angles)
    new_front_left_vertices[:, 2] = front_left_radii * np.sin(front_left_vertex_sweep_angles)

    # Rotate the coordinates of the front right vertex.
    new_front_right_vertices[:, 1] = front_right_radii * np.cos(front_right_vertex_sweep_angles)
    new_front_right_vertices[:, 2] = front_right_radii * np.sin(front_right_vertex_sweep_angles)

    # Rotate the coordinates of the back left vertex.
    new_back_left_vertices[:, 1] = back_left_radii * np.cos(back_left_vertex_sweep_angles)
    new_back_left_vertices[:, 2] = back_left_radii * np.sin(back_left_vertex_sweep_angles)

    # Rotate the coordinates of the back right vertex.
    new_back_right_vertices[:, 1] = back_right_radii * np.cos(back_right_vertex_sweep_angles)
    new_back_right_vertices[:, 2] = back_right_radii * np.sin(back_right_vertex_sweep_angles)

    # Update the class attributes with the new vertex locations.
    unsteady_aerodynamics_problem.front_left_vertices = new_front_left_vertices
    unsteady_aerodynamics_problem.front_right_vertices = new_front_right_vertices
    unsteady_aerodynamics_problem.back_left_vertices = new_back_left_vertices
    unsteady_aerodynamics_problem.back_right_vertices = new_back_right_vertices
