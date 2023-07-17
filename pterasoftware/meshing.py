"""This module contains useful functions for creating meshes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    mesh_wing: This function takes in an object of the Wing class and creates a
    quadrilateral mesh of its geometry, and then populates the object's panels with
    the mesh data.

    get_panel_vertices: This function calculates the vertices of the panels on a wing.

    ToDo: Update this method's documentation.
    get_normalized_projected_quarter_chords: This method returns the quarter chords
    of a collection of wing cross sections based on the coordinates of their leading
    and trailing edges. These quarter chords are also projected on to the YZ plane
    and normalized by their magnitudes.

    get_transpose_mcl_vectors: This function takes in the inner and outer airfoils of
    a wing cross section and its chordwise coordinates. It returns a list of four
    column vectors. They are, in order, the inner airfoil's local up direction,
    the inner airfoil's local back direction, the outer airfoil's local up direction,
    and the outer airfoil's local back direction.

    get_wing_section_panels: This function takes in arrays panel attributes and
    returns a 2D array of panel objects.
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
                wing_cross_section.trailing_edge + wing.leading_edge,
            )
        )

    normalized_projected_quarter_chords = get_normalized_projected_quarter_chords(
        wing_cross_sections_leading_edges,
        wing_cross_sections_trailing_edges,
        wing.unit_chordwise_vector,
    )

    # Get the number of wing cross sections and wing sections.
    num_wing_cross_sections = len(wing.wing_cross_sections)
    num_wing_sections = num_wing_cross_sections - 1

    # ToDo: Change this part's language to remove YZ plane mentions.
    # Then, construct the adjusted normal directions for each wing cross section.
    # Make the normals for the inner wing cross sections, where we need to merge
    # directions.
    if num_wing_cross_sections > 2:
        # Add together the adjacent normalized wing section quarter chords projected
        # onto the YZ plane.
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

        # Normalize the summed adjacent normalized wing section quarter chords
        # projected onto the YZ plane by their magnitudes.
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
        # Vertically stack the first and last normalized wing section quarter chords.
        wing_sections_local_unit_normals = np.vstack(
            (
                normalized_projected_quarter_chords[0, :],
                normalized_projected_quarter_chords[-1, :],
            )
        )

    # Create a list of the wing cross section chord lengths.
    wing_cross_sections_chord_lengths = np.array(
        [wing_cross_section.chord for wing_cross_section in wing.wing_cross_sections]
    )

    # Normalize the wing cross section back vectors by their magnitudes.
    wing_cross_sections_unit_chordwise_vectors = np.vstack(
        [
            wing_cross_section.unit_chordwise_vector
            for wing_cross_section in wing.wing_cross_sections
        ]
    )

    # ToDo: Clarify these calculations' documentation.
    # Then, construct the up direction for each wing cross section.
    wing_cross_sections_unit_up_vectors = np.vstack(
        [
            wing_cross_section.unit_up_vector
            for wing_cross_section in wing.wing_cross_sections
        ]
    )
    adj_wing_cross_sections_unit_up_vectors = np.cross(
        wing_cross_sections_unit_chordwise_vectors,
        wing_sections_local_unit_normals,
        axis=1,
    )

    # ToDo: This isn't correct for all situations now. Fix it.
    # If the wing is symmetric, set the local up position of the root cross section
    # to be the in local Z direction.
    if wing.symmetric:
        adj_wing_cross_sections_unit_up_vectors[0] = np.array([0, 0, 1])

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
            wing_cross_sections_unit_chordwise_vectors,
            wing_cross_sections_unit_up_vectors,
            wing_cross_sections_chord_lengths,
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
                wing_cross_sections_unit_chordwise_vectors,
                adj_wing_cross_sections_unit_up_vectors,
                wing_cross_sections_chord_lengths,
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

            # Reflect the vertices across the symmetry plane.
            front_inner_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                front_inner_vertices,
                wing.symmetry_unit_normal_vector,
                wing.leading_edge,
            )
            front_outer_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                front_outer_vertices,
                wing.symmetry_unit_normal_vector,
                wing.leading_edge,
            )
            back_inner_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                back_inner_vertices,
                wing.symmetry_unit_normal_vector,
                wing.leading_edge,
            )
            back_outer_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                back_outer_vertices,
                wing.symmetry_unit_normal_vector,
                wing.leading_edge,
            )

            # Get the reflected wing section's panels.
            wing_section_panels = get_wing_section_panels(
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

    # Iterate through the panels and populate their local position attributes.
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


def get_panel_vertices(
    inner_wing_cross_section_num,
    wing_cross_sections_unit_chordwise_vectors,
    wing_cross_sections_unit_up_vectors,
    wing_cross_sections_chord_lengths,
    wing_cross_sections_leading_edges,
    transpose_mcl_vectors,
    spanwise_coordinates,
):
    # ToDo: Update this documentation.
    """This function calculates the vertices of the panels on a wing.

    :param inner_wing_cross_section_num: int
        This parameter is the integer index of this wing's section's inner wing cross
        section.
    :param wing_cross_sections_unit_chordwise_vectors: array
        This parameter is an array of floats with size (X, 3), where X is this wing's
        number of wing cross sections. It holds two unit vectors that correspond to
        the wing cross sections' local-back directions, written in the body frame.
    :param wing_cross_sections_unit_up_vectors: array
        This parameter is an array of floats with size (X, 3), where X is this wing's
        number of wing cross sections. It holds two unit vectors that correspond to
        the wing cross sections' local-up directions, written in the body frame.
    :param wing_cross_sections_chord_lengths: array
        This parameter is a 1D array of floats with length X, where X is this wing's
        number of wing cross sections. It holds the chord lengths of this wing's wing
        cross section in meters.
    :param wing_cross_sections_leading_edges: array
        This parameter is an array of floats with size (Xx3), where X is this wing's
        number of wing cross sections. It holds the coordinates of the leading edge
        points of this wing's wing cross sections. The units are in meters.
    :param transpose_mcl_vectors: list
        This parameter is a list of 4 (M x 1) arrays of floats, where M is the number
        of chordwise points. The first array contains the local-up component of the
        mean-camber-line slope at each of the chordwise points along the inner wing
        cross section. The second array contains the local-back component of the
        mean-camber-line slope at each of the chordwise points along the inner wing
        cross section. The third and fourth arrays are the same but for the outer
        wing cross section.
    :param spanwise_coordinates: array
        This parameter is a 1D array of floats with length N, where N is the number
        of spanwise points. It holds the distances of each spanwise point along the
        wing cross section and is normalized from 0 to 1.
    :return: list
        This function returns a list with four arrays. Each array is size (MxNx3),
        where M is the number of chordwise points and N is the number of spanwise
        points. The arrays are the body frame coordinates of this wing's panels'
        front-inner, front-outer, back-inner, and back-outer vertices. The units are
        in meters.
    """
    [
        transpose_inner_mcl_up_vector,
        transpose_inner_mcl_back_vector,
        transpose_outer_mcl_up_vector,
        transpose_outer_mcl_back_vector,
    ] = transpose_mcl_vectors[:]

    # Convert the inner wing cross section's non dimensional local back airfoil frame
    # coordinates to meshed wing coordinates.
    inner_wing_cross_section_mcl_back = (
        wing_cross_sections_unit_chordwise_vectors[inner_wing_cross_section_num, :]
        * transpose_inner_mcl_back_vector
        * wing_cross_sections_chord_lengths[inner_wing_cross_section_num]
    )

    # Convert the inner wing cross section's non dimensional local up airfoil frame
    # coordinates to meshed wing coordinates.
    inner_wing_cross_section_mcl_up = (
        wing_cross_sections_unit_up_vectors[inner_wing_cross_section_num, :]
        * transpose_inner_mcl_up_vector
        * wing_cross_sections_chord_lengths[inner_wing_cross_section_num]
        # * wing_cross_sections_scaling_factors[inner_wing_cross_section_num]
    )

    # Define the index of this wing section's outer wing cross section.
    outer_wing_cross_section_num = inner_wing_cross_section_num + 1

    # Convert the outer wing cross section's non dimensional local back airfoil frame
    # coordinates to meshed wing coordinates.
    outer_wing_cross_section_mcl_back = (
        wing_cross_sections_unit_chordwise_vectors[outer_wing_cross_section_num, :]
        * transpose_outer_mcl_back_vector
        * wing_cross_sections_chord_lengths[outer_wing_cross_section_num]
    )

    # Convert the outer wing cross section's non dimensional local up airfoil frame
    # coordinates to meshed wing coordinates.
    outer_wing_cross_section_mcl_up = (
        wing_cross_sections_unit_up_vectors[outer_wing_cross_section_num, :]
        * transpose_outer_mcl_up_vector
        * wing_cross_sections_chord_lengths[outer_wing_cross_section_num]
    )

    # Convert the inner wing cross section's meshed wing coordinates to absolute
    # coordinates. This is size M x 3, where M is the number of chordwise points.
    inner_wing_cross_section_mcl = (
        wing_cross_sections_leading_edges[inner_wing_cross_section_num, :]
        + inner_wing_cross_section_mcl_back
        + inner_wing_cross_section_mcl_up
    )

    # Convert the outer wing cross section's meshed wing coordinates to absolute
    # coordinates. This is size M x 3, where M is the number of chordwise points.
    outer_wing_cross_section_mcl = (
        wing_cross_sections_leading_edges[outer_wing_cross_section_num, :]
        + outer_wing_cross_section_mcl_back
        + outer_wing_cross_section_mcl_up
    )

    # Find the vertices of the points on this wing section with interpolation. This
    # returns an M x N x 3 matrix, where M and N are the number of chordwise points
    # and spanwise points.
    wing_section_mcl_vertices = functions.interp_between_points(
        inner_wing_cross_section_mcl,
        outer_wing_cross_section_mcl,
        spanwise_coordinates,
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


def get_normalized_projected_quarter_chords(
    wing_cross_sections_leading_edges,
    wing_cross_sections_trailing_edges,
    wing_unit_chordwise_vector,
):
    # ToDo: Update this docstring to swap mentions of YZ plane for the custom plane.
    """This method returns the quarter chords of a collection of wing cross sections
    based on the coordinates of their leading and trailing edges. These quarter
    chords are also projected on to the YZ plane and normalized by their magnitudes.

    :param wing_cross_sections_leading_edges: array
        This parameter is an array of floats with size (X, 3), where X is this wing's
        number of wing cross sections. For each cross section, this array holds the
        body-frame coordinates of its leading edge point in meters.
    :param wing_cross_sections_trailing_edges: array
        This parameter is an array of floats with size (X, 3), where X is this wing's
        number of wing cross sections. For each cross section, this array holds the
        body-frame coordinates of its trailing edge point in meters.
    :param wing_unit_chordwise_vector
    :return normalized_projected_quarter_chords: array
        This functions returns an array of floats with size (X - 1, 3), where X is
        this wing's number of wing cross sections. This array holds each wing
        section's quarter chords projected on to the YZ plane and normalized by their
        magnitudes.
    """
    # Get the location of each wing cross section's quarter chord point.
    wing_cross_sections_quarter_chord_points = (
        wing_cross_sections_leading_edges
        + 0.25
        * (wing_cross_sections_trailing_edges - wing_cross_sections_leading_edges)
    )

    # Get a (L - 1) x 3 array of vectors connecting the wing cross section quarter
    # chord points, where L is the number of wing cross sections.
    quarter_chords = (
        wing_cross_sections_quarter_chord_points[1:, :]
        - wing_cross_sections_quarter_chord_points[:-1, :]
    )

    # ToDo: Update this documentation.
    # Get directions for transforming 2D airfoil data to 3D by the following steps.
    # Project quarter chords onto YZ plane and normalize. Create an L x 2 array with
    # just the y and z components of this wing section's quarter chord vectors.
    quarter_chords_dot_wing_unit_chordwise_vector = np.einsum(
        "ij,j->i", quarter_chords, wing_unit_chordwise_vector
    )
    c = np.expand_dims(wing_unit_chordwise_vector, axis=0)
    b = np.einsum(
        "j,ji->ji",
        quarter_chords_dot_wing_unit_chordwise_vector,
        c,
    )
    projected_quarter_chords = quarter_chords - b

    # Create a list of the lengths of each row of the projected_quarter_chords array.
    projected_quarter_chords_len = np.linalg.norm(projected_quarter_chords, axis=1)

    # Convert projected_quarter_chords_len into a column vector.
    transpose_projected_quarter_chords_len = np.expand_dims(
        projected_quarter_chords_len, axis=1
    )
    # Normalize the coordinates by the magnitudes
    normalized_projected_quarter_chords = (
        projected_quarter_chords / transpose_projected_quarter_chords_len
    )

    return normalized_projected_quarter_chords


def get_transpose_mcl_vectors(inner_airfoil, outer_airfoil, chordwise_coordinates):
    """This function takes in the inner and outer airfoils of a wing cross section
    and its chordwise coordinates. It returns a list of four column vectors. They
    are, in order, the inner airfoil's local up direction, the inner airfoil's local
    back direction, the outer airfoil's local up direction, and the outer airfoil's
    local back direction.

    :param inner_airfoil: Airfoil
        This is the wing cross section's inner airfoil object.
    :param outer_airfoil:
        This is the wing cross section's inner airfoil object.
    :param chordwise_coordinates: 1D array of floats
        This is a 1D array of the normalized chordwise coordinates where we'd like to
        sample each airfoil's mean camber line.
    :return: list of 4 (2x1) arrays
        This is a list of four column vectors. They are, in order, the inner
        airfoil's local up direction, the inner airfoil's local back direction,
        the outer airfoil's local up direction, and the outer airfoil's local back
        direction.
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


def get_wing_section_panels(
    front_left_vertices,
    front_right_vertices,
    back_left_vertices,
    back_right_vertices,
    is_trailing_edge,
    is_leading_edge,
):
    """This function takes in arrays panel attributes and returns a 2D array of panel
    objects.

    :param front_left_vertices: array of floats
        This is 3D array of size (MxNx3), where M is the number of chordwise panels,
        N is the number of spanwise panels, and the last dimension contains the x, y,
        and z coordinates of each panel's front left vertex.
    :param front_right_vertices: array of floats
        This is 3D array of size (MxNx3), where M is the number of chordwise panels,
        N is the number of spanwise panels, and the last dimension contains the x, y,
        and z coordinates of each panel's front right vertex.
    :param back_left_vertices: array of floats
        This is 3D array of size (MxNx3), where M is the number of chordwise panels,
        N is the number of spanwise panels, and the last dimension contains the x, y,
        and z coordinates of each panel's back left vertex.
    :param back_right_vertices: array of floats
        This is 3D array of size (MxNx3), where M is the number of chordwise panels,
        N is the number of spanwise panels, and the last dimension contains the x, y,
        and z coordinates of each panel's back right vertex.
    :param is_trailing_edge: 2D array of Booleans
        This is 2D array of True or False values that correspond to if the panel in
        each location is on the trailing edge of the wing.
    :param is_leading_edge: 2D array of Booleans
        This is 2D array of True or False values that correspond to if the panel in
        each location is on the trailing edge of the wing.
    :return panel_array: 2D array of Panels
        This is a 2D array of panel objects with the requested attributes.
    """
    num_chordwise_panels = front_left_vertices.shape[0]
    num_spanwise_panels = front_left_vertices.shape[1]

    # Initialize an empty array to hold the wing section's panels. The matrix is
    # size M x N, where M and N are the number of chordwise and spanwise panels.
    panels = np.empty((num_chordwise_panels, num_spanwise_panels), dtype=object)

    # Loop through the empty panels matrix and create a new panel object in each
    # slot.
    for chordwise_position in range(num_chordwise_panels):
        for spanwise_position in range(num_spanwise_panels):
            panels[chordwise_position, spanwise_position] = panel.Panel(
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

    return panels
