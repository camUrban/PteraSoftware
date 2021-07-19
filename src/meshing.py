# ToDo: Update this module's documentation.
"""This module contains useful functions for creating meshes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    mesh_wing: This function takes in an object of the Wing class and creates a
    quadrilateral mesh of its geometry, and then populates the object's panels with
    the mesh data.
"""
import numpy as np

from . import functions
from . import panel


def mesh_wing(wing):
    """This function takes in an object of the Wing class and creates a quadrilateral
    mesh of its geometry, and then populates the object's panels with the mesh data.

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

    # Get the chordwise coordinates.
    if wing.chordwise_spacing == "uniform":
        chordwise_coordinates = np.linspace(0, 1, num_chordwise_coordinates)
    else:
        chordwise_coordinates = functions.cosspace(0, 1, num_chordwise_coordinates)

    # Initialize two empty 0 x 3 arrays to hold the corners of each wing cross
    # section. They will eventually be L x 3 arrays, where L is number of wing cross
    # sections.
    wing_cross_sections_leading_edges = np.empty((0, 3))
    wing_cross_sections_trailing_edges = np.empty((0, 3))

    # Iterate through the meshed wing cross sections and vertically stack the global
    # location of each wing cross sections leading and trailing edges.
    # wing_cross_section.trailing_edge is a method that returns the wing cross section's
    # trailing edge's coordinates.
    for wing_cross_section in wing.wing_cross_sections:
        wing_cross_sections_leading_edges = np.vstack(
            (
                wing_cross_sections_leading_edges,
                wing_cross_section.leading_edge + wing.leading_edge,
            )
        )
        wing_cross_sections_trailing_edges = np.vstack(
            (
                wing_cross_sections_trailing_edges,
                wing_cross_section.trailing_edge() + wing.leading_edge,
            )
        )

    normalized_projected_quarter_chords = get_normalized_projected_quarter_chords(
        wing_cross_sections_leading_edges, wing_cross_sections_trailing_edges
    )

    # Get the number of wing cross sections.
    num_wing_cross_sections = len(wing.wing_cross_sections)
    num_wing_sections = num_wing_cross_sections - 1

    # Then, construct the normal directions for each wing cross section. Make the
    # normals for the inner wing cross sections, where we need to merge directions.
    if num_wing_cross_sections > 2:
        # Add together the adjacent normalized wing section quarter chords projected
        # onto the the YZ plane.
        wing_sections_local_normals = (
            normalized_projected_quarter_chords[:-1, :]
            + normalized_projected_quarter_chords[1:, :]
        )

        # Create a list of the magnitudes of the summed adjacent normalized wing
        # section quarter chords projected onto the YZ plane.
        wing_sections_local_normals_len = np.linalg.norm(
            wing_sections_local_normals, axis=1
        )

        # Convert the list to a column vector.
        transpose_wing_sections_local_normals_len = np.expand_dims(
            wing_sections_local_normals_len, axis=1
        )

        # Normalize the summed adjacent normalized wing section quarter chords projected
        # onto the YZ plane by their magnitudes.
        wing_sections_local_unit_normals = (
            wing_sections_local_normals / transpose_wing_sections_local_normals_len
        )

        # Vertically stack the first normalized wing section quarter chord, the inner
        # normalized wing section quarter chords, and the last normalized wing
        # section quarter chord.
        wing_sections_local_unit_normals = np.vstack(
            (
                normalized_projected_quarter_chords[0, :],
                wing_sections_local_unit_normals,
                normalized_projected_quarter_chords[-1, :],
            )
        )
    else:
        # Vertically stack the first normalized wing section quarter chord, and the
        # last normalized wing section quarter chord.
        wing_sections_local_unit_normals = np.vstack(
            (
                normalized_projected_quarter_chords[0, :],
                normalized_projected_quarter_chords[-1, :],
            )
        )

    # Then, construct the back directions for each wing cross section.
    wing_cross_sections_local_back_vectors = (
        wing_cross_sections_trailing_edges - wing_cross_sections_leading_edges
    )

    # Create a list of the wing cross section chord lengths.
    wing_cross_sections_chord_lengths = np.linalg.norm(
        wing_cross_sections_local_back_vectors, axis=1
    )

    # Convert the list to a column vector.
    transpose_wing_cross_sections_chord_lengths = np.expand_dims(
        wing_cross_sections_chord_lengths, axis=1
    )

    # Normalize the wing cross section back vectors by their magnitudes.
    wing_cross_sections_local_back_unit_vectors = (
        wing_cross_sections_local_back_vectors
        / transpose_wing_cross_sections_chord_lengths
    )

    # Then, construct the up direction for each wing cross section.
    wing_cross_sections_local_up_unit_vectors = np.cross(
        wing_cross_sections_local_back_unit_vectors,
        wing_sections_local_unit_normals,
        axis=1,
    )

    # If the wing is symmetric, set the local up position of the root cross section
    # to be the in local Z direction.
    if wing.symmetric:
        wing_cross_sections_local_up_unit_vectors[0] = np.array([0, 0, 1])

    # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to
    # compensate).
    wing_cross_sections_scaling_factors = get_wing_cross_section_scaling_factors(
        wing.symmetric, normalized_projected_quarter_chords
    )

    # Initialize an empty array that will hold the panels of this wing. It currently
    # has 0 columns and M rows, where M is the number of the wing's chordwise panels.
    wing_panels = np.empty((num_chordwise_panels, 0), dtype=object)

    # Make the panels for each wing section.
    for wing_section_num in range(num_wing_sections):
        # Define variables to hold the indices of this wing section's inner wing cross
        # section.
        inner_wing_cross_section_num = wing_section_num

        # Define the relevant wing cross sections.
        inner_wing_cross_section = wing.wing_cross_sections[
            inner_wing_cross_section_num
        ]
        outer_wing_cross_section = wing.wing_cross_sections[
            inner_wing_cross_section_num + 1
        ]

        # Define the airfoils at each wing cross section.
        inner_airfoil = inner_wing_cross_section.airfoil.add_control_surface(
            deflection=inner_wing_cross_section.control_surface_deflection,
            hinge_point=inner_wing_cross_section.control_surface_hinge_point,
        )
        outer_airfoil = outer_wing_cross_section.airfoil.add_control_surface(
            # The inner wing cross section dictates control surface deflections.
            deflection=inner_wing_cross_section.control_surface_deflection,
            hinge_point=inner_wing_cross_section.control_surface_hinge_point,
        )

        # Get the transposed mean camber line coordinates for the inner and outer
        # wing cross sections.
        transpose_mcl_vectors = get_transpose_mcl_vectors(
            inner_airfoil, outer_airfoil, chordwise_coordinates
        )

        # Define number of spanwise points and panels. This is based on the inner
        # wing cross section.
        num_spanwise_panels = inner_wing_cross_section.num_spanwise_panels
        num_spanwise_coordinates = num_spanwise_panels + 1

        # Get the spanwise coordinates. This is based on the inner wing cross section.
        if inner_wing_cross_section.spanwise_spacing == "uniform":
            spanwise_coordinates = np.linspace(0, 1, num_spanwise_coordinates)
        else:
            spanwise_coordinates = functions.cosspace(0, 1, num_spanwise_coordinates)

        # Get the panel vertices.
        [
            front_inner_vertices,
            front_outer_vertices,
            back_inner_vertices,
            back_outer_vertices,
        ] = get_panel_vertices(
            inner_wing_cross_section_num,
            wing_cross_sections_local_back_unit_vectors,
            wing_cross_sections_local_up_unit_vectors,
            wing_cross_sections_chord_lengths,
            wing_cross_sections_scaling_factors,
            wing_cross_sections_leading_edges,
            transpose_mcl_vectors,
            spanwise_coordinates,
        )

        # Compute a matrix that is M x N, where M and N are the number of chordwise
        # and spanwise panels. The values are either 1 if the panel at that location
        # is a trailing edge, or 0 if not.
        wing_section_is_trailing_edge = np.vstack(
            (
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
                np.ones((1, num_spanwise_panels), dtype=bool),
            )
        )

        # Compute a matrix that is M x N, where M and N are the number of chordwise
        # and spanwise panels. The values are either 1 if the panel at that location
        # is a leading edge, or 0 if not.
        wing_section_is_leading_edge = np.vstack(
            (
                np.ones((1, num_spanwise_panels), dtype=bool),
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
            )
        )

        # Get this wing section's panels.
        wing_section_panels = get_wing_section_panels(
            num_chordwise_panels=num_chordwise_panels,
            num_spanwise_panels=num_spanwise_panels,
            front_left_vertices=front_inner_vertices,
            front_right_vertices=front_outer_vertices,
            back_left_vertices=back_inner_vertices,
            back_right_vertices=back_outer_vertices,
            is_trailing_edge=wing_section_is_trailing_edge,
            is_leading_edge=wing_section_is_leading_edge,
        )

        # This wing section's panel matrix is stacked horizontally, to the right of the
        # wing's panel matrix.
        wing_panels = np.hstack((wing_panels, wing_section_panels))

        # Handle symmetry.
        if wing.symmetric:

            # The inner airfoil control surface type dictates this wing section's
            # control surface type.
            if inner_wing_cross_section.control_surface_type == "symmetric":

                # Define the airfoils at each wing cross section with flap-like control
                # surfaces.
                inner_airfoil = inner_wing_cross_section.airfoil.add_control_surface(
                    deflection=inner_wing_cross_section.control_surface_deflection,
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )
                outer_airfoil = outer_wing_cross_section.airfoil.add_control_surface(
                    deflection=inner_wing_cross_section.control_surface_deflection,
                    # The inner wing cross section dictates control surface deflections.
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )
            else:

                # Define the airfoils at each wing cross section with aileron-like
                # control surfaces.
                inner_airfoil = inner_wing_cross_section.airfoil.add_control_surface(
                    deflection=-inner_wing_cross_section.control_surface_deflection,
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )
                outer_airfoil = outer_wing_cross_section.airfoil.add_control_surface(
                    deflection=-inner_wing_cross_section.control_surface_deflection,
                    # The inner wing cross section dictates control surface deflections.
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )

            # Get the transposed mean camber line coordinates for the inner and outer
            # wing cross sections.
            transpose_mcl_vectors = get_transpose_mcl_vectors(
                inner_airfoil, outer_airfoil, chordwise_coordinates
            )

            # Get the panel vertices.
            [
                front_inner_vertices,
                front_outer_vertices,
                back_inner_vertices,
                back_outer_vertices,
            ] = get_panel_vertices(
                inner_wing_cross_section_num,
                wing_cross_sections_local_back_unit_vectors,
                wing_cross_sections_local_up_unit_vectors,
                wing_cross_sections_chord_lengths,
                wing_cross_sections_scaling_factors,
                wing_cross_sections_leading_edges,
                transpose_mcl_vectors,
                spanwise_coordinates,
            )

            # Compute a matrix that is M x N, where M and N are the number of
            # chordwise and spanwise panels. The values are either 1 if the panel at
            # that location is a trailing edge, or 0 if not.
            wing_section_is_trailing_edge = np.vstack(
                (
                    np.zeros(
                        (num_chordwise_panels - 1, num_spanwise_panels), dtype=bool
                    ),
                    np.ones((1, num_spanwise_panels), dtype=bool),
                )
            )

            # Compute a matrix that is M x N, where M and N are the number of
            # chordwise and spanwise panels. The values are either 1 if the panel at
            # that location is a leading edge, or 0 if not.
            wing_section_is_leading_edge = np.vstack(
                (
                    np.ones((1, num_spanwise_panels), dtype=bool),
                    np.zeros(
                        (num_chordwise_panels - 1, num_spanwise_panels), dtype=bool
                    ),
                )
            )

            # Reflect the vertices across the XZ plane.
            front_inner_vertices_reflected = functions.reflect_over_xz_plane(
                front_inner_vertices
            )
            front_outer_vertices_reflected = functions.reflect_over_xz_plane(
                front_outer_vertices
            )
            back_inner_vertices_reflected = functions.reflect_over_xz_plane(
                back_inner_vertices
            )
            back_outer_vertices_reflected = functions.reflect_over_xz_plane(
                back_outer_vertices
            )

            # Get the reflected wing section's panels.
            wing_section_panels = get_wing_section_panels(
                num_chordwise_panels=num_chordwise_panels,
                num_spanwise_panels=num_spanwise_panels,
                front_left_vertices=front_outer_vertices_reflected,
                front_right_vertices=front_inner_vertices_reflected,
                back_left_vertices=back_outer_vertices_reflected,
                back_right_vertices=back_inner_vertices_reflected,
                is_trailing_edge=wing_section_is_trailing_edge,
                is_leading_edge=wing_section_is_leading_edge,
            )

            # This wing section's panel matrix is stacked horizontally, to the left
            # of the wing's panel matrix.
            wing_panels = np.hstack((np.flip(wing_section_panels, axis=1), wing_panels))

    # Iterate through the panels and populate their left and right edge flags. Also
    # populate their local position attributes.
    for chordwise_position in range(wing.num_chordwise_panels):
        for spanwise_position in range(wing.num_spanwise_panels):
            this_panel = wing_panels[chordwise_position, spanwise_position]
            this_panel.local_chordwise_position = chordwise_position
            this_panel.local_spanwise_position = spanwise_position
            if spanwise_position == 0:
                this_panel.is_left_edge = True
            else:
                this_panel.is_left_edge = False
            if spanwise_position == wing.num_spanwise_panels - 1:
                this_panel.is_right_edge = True
            else:
                this_panel.is_right_edge = False

    # Populate the wing's panels attribute.
    wing.panels = wing_panels


# ToDo: Document the following function.
def get_wing_cross_section_scaling_factors(
    symmetric, wing_section_quarter_chords_proj_yz_norm
):
    """Get the scaling factors (airfoils at dihedral breaks need to be "taller" to
    compensate).

    :param symmetric:
    :param wing_section_quarter_chords_proj_yz_norm:
    :return:
    """
    num_wing_cross_sections = len(wing_section_quarter_chords_proj_yz_norm) + 1

    # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to
    # compensate).
    wing_cross_section_scaling_factors = np.ones(num_wing_cross_sections)

    for i in range(num_wing_cross_sections):
        if i == 0:
            if symmetric:
                first_chord_norm = wing_section_quarter_chords_proj_yz_norm[0]
                mirrored_first_chord_norm = first_chord_norm * np.array([1, 1, -1])

                product = first_chord_norm * mirrored_first_chord_norm
                collapsed_product = np.sum(product)
                this_scaling_factor = 1 / np.sqrt((1 + collapsed_product) / 2)
            else:
                this_scaling_factor = 1
        elif i == num_wing_cross_sections - 1:
            this_scaling_factor = 1
        else:
            this_chord_norm = wing_section_quarter_chords_proj_yz_norm[i - 1, :]
            next_chord_norm = wing_section_quarter_chords_proj_yz_norm[i, :]

            product = this_chord_norm * next_chord_norm
            collapsed_product = np.sum(product)
            this_scaling_factor = 1 / np.sqrt((1 + collapsed_product) / 2)

        wing_cross_section_scaling_factors[i] = this_scaling_factor

    return wing_cross_section_scaling_factors


# ToDo: Document the following function.
def get_panel_vertices(
    inner_wing_cross_section_num,
    wing_cross_sections_local_back_unit_vectors,
    wing_cross_sections_local_up_unit_vectors,
    wing_cross_sections_chord_lengths,
    wing_cross_sections_scaling_factors,
    wing_cross_sections_leading_edges,
    transpose_mcl_vectors,
    spanwise_coordinates,
):
    """

    :param inner_wing_cross_section_num:
    :param wing_cross_sections_local_back_unit_vectors:
    :param wing_cross_sections_local_up_unit_vectors:
    :param wing_cross_sections_chord_lengths:
    :param wing_cross_sections_scaling_factors:
    :param wing_cross_sections_leading_edges:
    :param transpose_mcl_vectors:
    :param spanwise_coordinates:
    :return:
    """
    [
        transpose_inner_mcl_up_vector,
        transpose_inner_mcl_back_vector,
        transpose_outer_mcl_up_vector,
        transpose_outer_mcl_back_vector,
    ] = transpose_mcl_vectors[:]

    # Convert the inner wing cross section's non dimensional local back airfoil frame
    # coordinates to meshed wing coordinates.
    inner_wing_cross_section_mcl_local_back = (
        wing_cross_sections_local_back_unit_vectors[inner_wing_cross_section_num, :]
        * transpose_inner_mcl_back_vector
        * wing_cross_sections_chord_lengths[inner_wing_cross_section_num]
    )

    # Convert the inner wing cross section's non dimensional local up airfoil frame
    # coordinates to meshed wing coordinates.
    inner_wing_cross_section_mcl_local_up = (
        wing_cross_sections_local_up_unit_vectors[inner_wing_cross_section_num, :]
        * transpose_inner_mcl_up_vector
        * wing_cross_sections_chord_lengths[inner_wing_cross_section_num]
        * wing_cross_sections_scaling_factors[inner_wing_cross_section_num]
    )

    # Define the index of this wing section's outer wing cross section.
    outer_wing_cross_section_num = inner_wing_cross_section_num + 1

    # Convert the outer wing cross section's non dimensional local back airfoil frame
    # coordinates to meshed wing coordinates.
    outer_wing_cross_section_mcl_local_back = (
        wing_cross_sections_local_back_unit_vectors[outer_wing_cross_section_num, :]
        * transpose_outer_mcl_back_vector
        * wing_cross_sections_chord_lengths[outer_wing_cross_section_num]
    )

    # Convert the outer wing cross section's non dimensional local up airfoil frame
    # coordinates to meshed wing coordinates.
    outer_wing_cross_section_mcl_local_up = (
        wing_cross_sections_local_up_unit_vectors[outer_wing_cross_section_num, :]
        * transpose_outer_mcl_up_vector
        * wing_cross_sections_chord_lengths[outer_wing_cross_section_num]
        * wing_cross_sections_scaling_factors[outer_wing_cross_section_num]
    )

    # Convert the inner wing cross section's meshed wing coordinates to absolute
    # coordinates. This is size M x 3, where M is the number of chordwise points.
    inner_wing_cross_section_mcl = (
        wing_cross_sections_leading_edges[inner_wing_cross_section_num, :]
        + inner_wing_cross_section_mcl_local_back
        + inner_wing_cross_section_mcl_local_up
    )

    # Convert the outer wing cross section's meshed wing coordinates to absolute
    # coordinates. This is size M x 3, where M is the number of chordwise points.
    outer_wing_cross_section_mcl = (
        wing_cross_sections_leading_edges[outer_wing_cross_section_num, :]
        + outer_wing_cross_section_mcl_local_back
        + outer_wing_cross_section_mcl_local_up
    )

    # Make section_mcl_coordinates: M x N x 3 array of mean camberline
    # coordinates. The first index is chordwise point number, second index is
    # spanwise point number, third is the x, y, and z coordinates. M is the
    # number of chordwise points. N is the number of spanwise points. Put a
    # reversed version (from 1 to 0) of the non dimensional spanwise coordinates
    # in a row vector. This is size 1 x N, where N is the number of spanwise
    # points.
    reversed_nondim_spanwise_coordinates_row_vector = np.expand_dims(
        (1 - spanwise_coordinates), 0
    )

    # Convert the reversed non dimensional spanwise coordinate row vector (from 1
    # to 0) to a matrix. This is size 1 x N x 1, where N is the number of
    # spanwise points.
    reversed_nondim_spanwise_coordinates_matrix = np.expand_dims(
        reversed_nondim_spanwise_coordinates_row_vector, 2
    )

    # Convert the inner and outer wing cross section's mean camberline coordinates
    # column vectors to matrices. These are size M x 1 x 3, where M is the number
    # of chordwise points.
    inner_wing_cross_section_mcl_matrix = np.expand_dims(
        inner_wing_cross_section_mcl, 1
    )
    outer_wing_cross_section_mcl_matrix = np.expand_dims(
        outer_wing_cross_section_mcl, 1
    )

    # Put the non dimensional spanwise coordinates (from 0 to 1) in a row vector.
    # This is size 1 x N, where N is the number of spanwise points.
    nondim_spanwise_coordinates_row_vector = np.expand_dims(spanwise_coordinates, 0)

    # Convert the non dimensional spanwise coordinate row vector (from to 0 to 1)
    # to a matrix. This is size 1 x N x 1, where N is the number of spanwise
    # points.
    nondim_spanwise_coordinates_matrix = np.expand_dims(
        nondim_spanwise_coordinates_row_vector, 2
    )

    # Linearly interpolate between inner and outer wing cross sections. This uses
    # the following equation:
    #
    # f(a, b, i) = i * a + (1 - i) * b
    #
    # "a" is an N x 3 array of the coordinates points along the outer wing cross
    # section's mean camber line.
    #
    # "b" is an N x 3 array of the coordinates of points along the inner wing
    # cross section's mean camber line.
    #
    # "i" is a 1D array (or vector) of length M that holds the nondimensionalized
    # spanwise panel spacing from 0 to 1.
    #
    # This produces a M x N x 3 array where each slot holds the coordinates of a
    # point on the surface between the inner and outer wing cross sections.
    wing_section_mcl_vertices = (
        reversed_nondim_spanwise_coordinates_matrix
        * inner_wing_cross_section_mcl_matrix
        + nondim_spanwise_coordinates_matrix * outer_wing_cross_section_mcl_matrix
    )

    # Compute the corners of each panel.
    front_inner_vertices = wing_section_mcl_vertices[:-1, :-1, :]
    front_outer_vertices = wing_section_mcl_vertices[:-1, 1:, :]
    back_inner_vertices = wing_section_mcl_vertices[1:, :-1, :]
    back_outer_vertices = wing_section_mcl_vertices[1:, 1:, :]

    return [
        front_inner_vertices,
        front_outer_vertices,
        back_inner_vertices,
        back_outer_vertices,
    ]


# ToDo: Document this function.
def get_normalized_projected_quarter_chords(
    wing_cross_sections_leading_edges, wing_cross_sections_trailing_edges
):
    """

    :param wing_cross_sections_leading_edges:
    :param wing_cross_sections_trailing_edges:
    :return:
    """
    # Get the location of each wing cross section's quarter chord point.
    wing_cross_sections_quarter_chord_points = (
        wing_cross_sections_leading_edges
        + 0.25
        * (wing_cross_sections_trailing_edges - wing_cross_sections_leading_edges)
    )

    # Get a (L - 1) x 3 array of vectors connecting the wing cross section quarter chord
    # points, where L is the number of wing cross sections.
    quarter_chords = (
        wing_cross_sections_quarter_chord_points[1:, :]
        - wing_cross_sections_quarter_chord_points[:-1, :]
    )

    # Get directions for transforming 2D airfoil data to 3D by the following steps.
    #
    # Project quarter chords onto YZ plane and normalize.
    #
    # Create a L x 2 array with just the y and z components of this wing section's
    # quarter chord vectors.
    projected_quarter_chords = quarter_chords[:, 1:]

    # Create a list of the lengths of each row of the projected_quarter_chords
    # array.
    projected_quarter_chords_len = np.linalg.norm(projected_quarter_chords, axis=1)

    # Convert projected_quarter_chords_len into a column vector.
    transpose_projected_quarter_chords_len = np.expand_dims(
        projected_quarter_chords_len, axis=1
    )
    # Normalize the coordinates by the magnitudes
    normalized_projected_quarter_chords = (
        projected_quarter_chords / transpose_projected_quarter_chords_len
    )

    # Create a column vector of all zeros with height equal to the number of quarter
    # chord vectors
    column_of_zeros = np.zeros((len(quarter_chords), 1))

    # Horizontally stack the zero column vector with the
    # normalized_projected_quarter_chords to give each normalized projected quarter
    # chord an X coordinate.
    normalized_projected_quarter_chords = np.hstack(
        (column_of_zeros, normalized_projected_quarter_chords)
    )

    return normalized_projected_quarter_chords


# ToDo: Document this function.
def get_transpose_mcl_vectors(inner_airfoil, outer_airfoil, chordwise_coordinates):
    """

    :param inner_airfoil:
    :param outer_airfoil:
    :param chordwise_coordinates:
    :return:
    """
    # Make the mean camber lines for each wing cross section. First index is
    # point number, second index is the coordinates in the airfoil frame.
    inner_mcl = inner_airfoil.get_downsampled_mcl(chordwise_coordinates)
    outer_mcl = outer_airfoil.get_downsampled_mcl(chordwise_coordinates)

    # Put the inner wing cross section's local up airfoil frame coordinates
    # in a column vector.
    transpose_inner_mcl_up_vector = np.expand_dims(inner_mcl[:, 1], 1)

    # Put the inner wing cross section's local back airfoil frame coordinates
    # in a column vector.
    transpose_inner_mcl_back_vector = np.expand_dims(inner_mcl[:, 0], 1)

    # Put the outer wing cross section's local up airfoil frame coordinates
    # in a column vector.
    transpose_outer_mcl_up_vector = np.expand_dims(outer_mcl[:, 1], 1)

    # Put the outer wing cross section's local back airfoil frame coordinates
    # in a column vector.
    transpose_outer_mcl_back_vector = np.expand_dims(outer_mcl[:, 0], 1)

    return [
        transpose_inner_mcl_up_vector,
        transpose_inner_mcl_back_vector,
        transpose_outer_mcl_up_vector,
        transpose_outer_mcl_back_vector,
    ]


# ToDo: Document this function.
def get_wing_section_panels(
    num_chordwise_panels,
    num_spanwise_panels,
    front_left_vertices,
    front_right_vertices,
    back_left_vertices,
    back_right_vertices,
    is_trailing_edge,
    is_leading_edge,
):
    """

    :param num_chordwise_panels:
    :param num_spanwise_panels:
    :param front_left_vertices:
    :param front_right_vertices:
    :param back_left_vertices:
    :param back_right_vertices:
    :param is_trailing_edge:
    :param is_leading_edge:
    :return:
    """
    # Initialize an empty array to hold the wing section's panels. The matrix is
    # size M x N, where M and N are the number of chordwise and spanwise panels.
    wing_section_panels = np.empty(
        (num_chordwise_panels, num_spanwise_panels), dtype=object
    )

    # Loop through the empty panels matrix and create a new panel object in each
    # slot.
    for chordwise_position in range(num_chordwise_panels):
        for spanwise_position in range(num_spanwise_panels):
            wing_section_panels[chordwise_position, spanwise_position] = panel.Panel(
                front_left_vertex=front_left_vertices[
                    chordwise_position, spanwise_position
                ],
                front_right_vertex=front_right_vertices[
                    chordwise_position, spanwise_position
                ],
                back_left_vertex=back_left_vertices[
                    chordwise_position, spanwise_position
                ],
                back_right_vertex=back_right_vertices[
                    chordwise_position, spanwise_position
                ],
                is_trailing_edge=is_trailing_edge[
                    chordwise_position, spanwise_position
                ],
                is_leading_edge=is_leading_edge[chordwise_position, spanwise_position],
            )

    return wing_section_panels
