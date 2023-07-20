"""This module contains useful functions for creating meshes.

This module contains the following classes:
    None

This module contains the following exceptions:
    None

This module contains the following functions:
    mesh_wing: This function takes in an object of the Wing class and creates a
    quadrilateral mesh of its geometry, and then populates the object's panels with
    the mesh data.

    get_panel_vertices: This function calculates the vertices of the panels on a wing
    section.

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

    # Get the number of wing cross sections and wing sections.
    num_wing_cross_sections = len(wing.wing_cross_sections)
    num_wing_sections = num_wing_cross_sections - 1

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
        ] = get_wing_section_panel_vertices(
            wing.leading_edge,
            inner_wing_cross_section,
            outer_wing_cross_section,
            transpose_mcl_vectors,
            spanwise_coordinates,
        )

        # Compute a matrix that is (M, N), where M and N are the number of chordwise
        # and spanwise panels. The values are either 1 if the panel at that location
        # is a trailing edge, or 0 if not.
        wing_section_is_trailing_edge = np.vstack(
            (
                np.zeros((num_chordwise_panels - 1, num_spanwise_panels), dtype=bool),
                np.ones((1, num_spanwise_panels), dtype=bool),
            )
        )

        # Compute a matrix that is (M, N), where M and N are the number of chordwise
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
            ] = get_wing_section_panel_vertices(
                wing.leading_edge,
                inner_wing_cross_section,
                outer_wing_cross_section,
                transpose_mcl_vectors,
                spanwise_coordinates,
            )

            # Compute a matrix that is (M, N), where M and N are the number of
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

            # Compute a matrix that is (M, N), where M and N are the number of
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
                wing.unit_normal_vector,
                wing.leading_edge,
            )
            front_outer_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                front_outer_vertices,
                wing.unit_normal_vector,
                wing.leading_edge,
            )
            back_inner_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                back_inner_vertices,
                wing.unit_normal_vector,
                wing.leading_edge,
            )
            back_outer_vertices_reflected = np.apply_along_axis(
                functions.reflect_point_across_plane,
                -1,
                back_outer_vertices,
                wing.unit_normal_vector,
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


def get_wing_section_panel_vertices(
    wing_leading_edge,
    inner_wing_cross_section,
    outer_wing_cross_section,
    transpose_mcl_vectors,
    spanwise_coordinates,
):
    """This function calculates the vertices of the panels on a wing section.

    :param wing_leading_edge: (3,) array of floats
        This is an array of the wing's leading edge coordinates. The units are meters.
    :param inner_wing_cross_section: WingCrossSection
        This is this wing section's inner Wing Cross Section object.
    :param outer_wing_cross_section: WingCrossSection
        This is this wing section's outer Wing Cross Section object.
    :param transpose_mcl_vectors: list of 4 (M, 1) arrays of floats
        This parameter is a list of 4 (M, 1) arrays where M is the number of
        chordwise points. The first array contains the local-up component of the
        mean-camber-line's slope at each of the chordwise points along the inner wing
        cross section. The second array contains the local-back component of the
        mean-camber-line's slope at each of the chordwise points along the inner wing
        cross section. The third and fourth arrays are the same but for the outer
        wing cross section instead of the inner wing cross section. The units are
        meters.
    :param spanwise_coordinates: (N, 1) array of floats
        This parameter is a (N, 1) array of floats, where N is the number of spanwise
        points. It holds the distances of each spanwise point along the wing section
        and is normalized from 0 to 1. These values are unitless.
    :return: list of 4 (M, N, 3) arrays of floats
        This function returns a list with four (M, N, 3) arrays, where M is the
        number of chordwise points and N is the number of spanwise points. The arrays
        are the coordinates of this wing's panels' front-inner, front-outer,
        back-inner, and back-outer vertices. The units are in meters.
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
        inner_wing_cross_section.unit_chordwise_vector
        * inner_wing_cross_section.chord
        * transpose_inner_mcl_back_vector
    )

    # Convert the inner wing cross section's non dimensional local up airfoil frame
    # coordinates to meshed wing coordinates.
    inner_wing_cross_section_mcl_up = (
        inner_wing_cross_section.unit_up_vector
        * inner_wing_cross_section.chord
        * transpose_inner_mcl_up_vector
    )

    # Convert the outer wing cross section's non dimensional local back airfoil frame
    # coordinates to meshed wing coordinates.
    outer_wing_cross_section_mcl_back = (
        outer_wing_cross_section.unit_chordwise_vector
        * outer_wing_cross_section.chord
        * transpose_outer_mcl_back_vector
    )

    # Convert the outer wing cross section's non dimensional local up airfoil frame
    # coordinates to meshed wing coordinates.
    outer_wing_cross_section_mcl_up = (
        outer_wing_cross_section.unit_up_vector
        * outer_wing_cross_section.chord
        * transpose_outer_mcl_up_vector
    )

    # Convert the inner wing cross section's meshed wing coordinates to absolute
    # coordinates. This is size (M, 3) where M is the number of chordwise points.
    inner_wing_cross_section_mcl = (
        wing_leading_edge
        + inner_wing_cross_section.leading_edge
        + inner_wing_cross_section_mcl_back
        + inner_wing_cross_section_mcl_up
    )

    # Convert the outer wing cross section's meshed wing coordinates to absolute
    # coordinates. This is size (M, 3) where M is the number of chordwise points.
    outer_wing_cross_section_mcl = (
        wing_leading_edge
        + outer_wing_cross_section.leading_edge
        + outer_wing_cross_section_mcl_back
        + outer_wing_cross_section_mcl_up
    )

    # Find the vertices of the points on this wing section with interpolation. This
    # returns an (M, N, 3) array, where M and N are the number of chordwise points
    # and spanwise points.
    wing_section_mcl_vertices = functions.interp_between_points(
        inner_wing_cross_section_mcl,
        outer_wing_cross_section_mcl,
        spanwise_coordinates,
    )

    # Extract the coordinates for corners of each panel.
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
