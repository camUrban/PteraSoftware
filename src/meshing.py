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

    # Initialize an empty array that will hold the panels of this wing. It currently
    # has 0 columns and M rows, where M is the number of the wing's chordwise panels.
    panels = np.empty((num_chordwise_panels, 0), dtype=object)

    # Get the chordwise coordinates.
    if wing.chordwise_spacing == "uniform":
        nondim_chordwise_coordinates = np.linspace(
            0,
            1,
            num_chordwise_coordinates,
            endpoint=True,
        )
    elif wing.chordwise_spacing == "cosine":
        nondim_chordwise_coordinates = functions.cosspace(
            0,
            1,
            num_chordwise_coordinates,
            endpoint=True,
        )
    else:
        raise Exception("Bad value of wing.chordwise_spacing!")

    # Initialize two empty 0 x 3 ndarrays to hold the corners of each wing cross
    # section. They will eventually be L x 3 ndarrays, where L is number of wing
    # cross sections.
    wing_cross_section_xyz_le = np.empty((0, 3))
    wing_cross_section_xyz_te = np.empty((0, 3))

    # Iterate through the meshed wing cross sections and vertically stack the global
    # location of each wing cross sections leading and trailing edges.
    # wing_cross_section.xyz_te is a method that returns the wing cross section's
    # trailing edge's coordinates.
    for wing_cross_section in wing.wing_cross_sections:
        wing_cross_section_xyz_le = np.vstack(
            (wing_cross_section_xyz_le, wing_cross_section.xyz_le + wing.xyz_le)
        )
        wing_cross_section_xyz_te = np.vstack(
            (wing_cross_section_xyz_te, wing_cross_section.xyz_te() + wing.xyz_le)
        )

    # Get the quarter chord vectors, which are a L x 3 array of points which are the
    # quarter-chord points of each wing cross section, where L is the number of wing
    # cross sections.
    wing_cross_section_xyz_quarter_chords = wing_cross_section_xyz_le + 0.25 * (
        wing_cross_section_xyz_te - wing_cross_section_xyz_le
    )

    # Get a (L - 1) x 3 array of vectors connecting the wing cross section quarter chord
    # points, where L is the number of wing cross sections.
    wing_section_quarter_chords = (
        wing_cross_section_xyz_quarter_chords[1:, :]
        - wing_cross_section_xyz_quarter_chords[:-1, :]
    )

    """Get directions for transforming 2D airfoil data to 3D by the following steps. 
    
    Project quarter chords onto YZ plane and normalize. 
    
    Create a L x 2 array with just the y and z components of this wing section's 
    quarter chord vectors."""
    wing_section_quarter_chords_yz = wing_section_quarter_chords[:, 1:]

    # Create a list of the magnitudes of each row of the wing_section_quarter_chords_yz
    # array.
    wing_section_quarter_chords_yz_magnitude_list = np.linalg.norm(
        wing_section_quarter_chords_yz, axis=1
    )

    # Convert wing_section_quarter_chords_yz_magnitude_list into a column vector.
    wing_section_quarter_chords_yz_magnitude_column_vector = np.expand_dims(
        wing_section_quarter_chords_yz_magnitude_list, axis=1
    )
    # Normalize the coordinates by the magnitudes
    wing_section_quarter_chords_yz_norm_magnitudes = (
        wing_section_quarter_chords_yz
        / wing_section_quarter_chords_yz_magnitude_column_vector
    )

    # Calculate the number of quarter chord vectors
    num_quarter_chords = wing_section_quarter_chords_yz_magnitude_column_vector.shape[0]

    # Create a column vector of all zeros with height equal to the number of quarter
    # chord vectors
    zero_column_vector_stand_in_for_quarter_chords_x_values = np.zeros(
        (num_quarter_chords, 1)
    )

    # Horizontally stack the zero column vector with the
    # wing_section_quarter_chords_yz_norm_magnitudes to produce the normalized wing
    # section quarter chords projected onto the YZ plane.
    wing_section_quarter_chords_proj_yz_norm = np.hstack(
        (
            zero_column_vector_stand_in_for_quarter_chords_x_values,
            wing_section_quarter_chords_yz_norm_magnitudes,
        )
    )

    # Get the number of wing cross sections.
    num_wing_cross_sections = len(wing.wing_cross_sections)
    num_wing_sections = num_wing_cross_sections - 1

    # Then, construct the normal directions for each wing cross section. Make the
    # normals for the inner wing cross sections, where we need to merge directions.
    if num_wing_cross_sections > 2:
        # Add together the adjacent normalized wing section quarter chords projected
        # onto the the YZ plane.
        wing_cross_sections_local_normal_inners_non_norm = (
            wing_section_quarter_chords_proj_yz_norm[:-1, :]
            + wing_section_quarter_chords_proj_yz_norm[1:, :]
        )

        # Create a list of the magnitudes of the summed adjacent normalized wing
        # section quarter chords projected onto the YZ plane.
        wing_cross_sections_local_normal_inners_mag_list = np.linalg.norm(
            wing_cross_sections_local_normal_inners_non_norm, axis=1
        )

        # Convert the list to a column vector.
        wing_cross_section_local_normal_inners_mag_column_vector = np.expand_dims(
            wing_cross_sections_local_normal_inners_mag_list, axis=1
        )

        # Normalize the summed adjacent normalized wing section quarter chords projected
        # onto the YZ plane by their magnitudes.
        wing_cross_section_local_normal_inners_norm = (
            wing_cross_sections_local_normal_inners_non_norm
            / wing_cross_section_local_normal_inners_mag_column_vector
        )

        # Vertically stack the first normalized wing section quarter chord, the inner
        # normalized wing section quarter chords, and the last normalized wing
        # section quarter chord.
        wing_cross_sections_local_normal = np.vstack(
            (
                wing_section_quarter_chords_proj_yz_norm[0, :],
                wing_cross_section_local_normal_inners_norm,
                wing_section_quarter_chords_proj_yz_norm[-1, :],
            )
        )
    else:
        # Vertically stack the first normalized wing section quarter chord, and the
        # last normalized wing section quarter chord.
        wing_cross_sections_local_normal = np.vstack(
            (
                wing_section_quarter_chords_proj_yz_norm[0, :],
                wing_section_quarter_chords_proj_yz_norm[-1, :],
            )
        )

    # Then, construct the back directions for each wing cross section.
    wing_cross_section_local_back_non_norm = (
        wing_cross_section_xyz_te - wing_cross_section_xyz_le
    )

    # Create a list of the wing cross section chord lengths.
    wing_cross_section_chord_lengths = np.linalg.norm(
        wing_cross_section_local_back_non_norm, axis=1
    )

    # Convert the list to a column vector.
    wing_cross_section_chord_length_column_vector = np.expand_dims(
        wing_cross_section_chord_lengths, axis=1
    )

    # Normalize the wing cross section back vectors by their magnitudes.
    wing_cross_section_local_back_norm = (
        wing_cross_section_local_back_non_norm
        / wing_cross_section_chord_length_column_vector
    )

    # Then, construct the up direction for each wing cross section.
    wing_cross_section_local_up = np.cross(
        wing_cross_section_local_back_norm, wing_cross_sections_local_normal, axis=1
    )

    # If the wing is symmetric, set the local up position of the root cross section
    # to be the in local Z direction.
    if wing.symmetric:
        wing_cross_section_local_up[0] = np.array([0, 0, 1])

    # Get the scaling factor (airfoils at dihedral breaks need to be "taller" to
    # compensate).
    wing_cross_section_scaling_factors = np.ones(num_wing_cross_sections)

    # ToDo: Document the following algorithm.
    for i in range(num_wing_cross_sections):
        if i == 0:
            if wing.symmetric:
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

    # Make the panels for each wing section.
    for wing_section_num in range(num_wing_sections):
        # Define variables to hold the indices of the wing cross sections associated
        # with each wing section.
        inner_wing_cross_section_num = wing_section_num
        outer_wing_cross_section_num = inner_wing_cross_section_num + 1

        # Define the relevant wing cross sections.
        inner_wing_cross_section = wing.wing_cross_sections[
            inner_wing_cross_section_num
        ]
        outer_wing_cross_section = wing.wing_cross_sections[
            outer_wing_cross_section_num
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

        # Make the mean camber lines for each wing cross section. First index is point
        # number, second index is xyz.
        inner_wing_cross_section_mcl_nondim = inner_airfoil.get_downsampled_mcl(
            nondim_chordwise_coordinates
        )
        outer_wing_cross_section_mcl_nondim = outer_airfoil.get_downsampled_mcl(
            nondim_chordwise_coordinates
        )

        # Put the inner wing cross section's local up airfoil frame coordinates in a
        # column vector.
        inner_wing_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(
            inner_wing_cross_section_mcl_nondim[:, 1], 1
        )

        # Put the inner wing cross section's local back airfoil frame coordinates in a
        # column vector.
        inner_wing_cross_section_mcl_nondim_local_back_column_vector = np.expand_dims(
            inner_wing_cross_section_mcl_nondim[:, 0], 1
        )
        # Put the outer wing cross section's local up airfoil frame coordinates in a
        # column vector.
        outer_wing_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(
            outer_wing_cross_section_mcl_nondim[:, 1], 1
        )

        # Put the outer wing cross section's local back airfoil frame coordinates in a
        # column vector.
        outer_wing_cross_section_mcl_nondim_local_back_column_vector = np.expand_dims(
            outer_wing_cross_section_mcl_nondim[:, 0], 1
        )

        # Convert the inner wing cross section's non dimensional local back airfoil frame
        # coordinates to meshed wing coordinates.
        inner_wing_cross_section_mcl_local_back = (
            wing_cross_section_local_back_norm[inner_wing_cross_section_num, :]
            * inner_wing_cross_section_mcl_nondim_local_back_column_vector
            * wing_cross_section_chord_lengths[inner_wing_cross_section_num]
        )

        # Convert the inner wing cross section's non dimensional local up airfoil frame
        # coordinates to meshed wing coordinates.
        inner_wing_cross_section_mcl_local_up = (
            wing_cross_section_local_up[inner_wing_cross_section_num, :]
            * inner_wing_cross_section_mcl_nondim_local_up_column_vector
            * wing_cross_section_chord_lengths[inner_wing_cross_section_num]
            * wing_cross_section_scaling_factors[inner_wing_cross_section_num]
        )

        # Convert the outer wing cross section's non dimensional local back airfoil frame
        # coordinates to meshed wing coordinates.
        outer_wing_cross_section_mcl_local_back = (
            wing_cross_section_local_back_norm[outer_wing_cross_section_num, :]
            * outer_wing_cross_section_mcl_nondim_local_back_column_vector
            * wing_cross_section_chord_lengths[outer_wing_cross_section_num]
        )

        # Convert the outer wing cross section's non dimensional local up airfoil frame
        # coordinates to meshed wing coordinates.
        outer_wing_cross_section_mcl_local_up = (
            wing_cross_section_local_up[outer_wing_cross_section_num, :]
            * outer_wing_cross_section_mcl_nondim_local_up_column_vector
            * wing_cross_section_chord_lengths[outer_wing_cross_section_num]
            * wing_cross_section_scaling_factors[outer_wing_cross_section_num]
        )

        # Convert the inner wing cross section's meshed wing coordinates to absolute
        # coordinates. This is size M x 3, where M is the number of chordwise points.
        inner_wing_cross_section_mcl = (
            wing_cross_section_xyz_le[inner_wing_cross_section_num, :]
            + inner_wing_cross_section_mcl_local_back
            + inner_wing_cross_section_mcl_local_up
        )

        # Convert the outer wing cross section's meshed wing coordinates to absolute
        # coordinates. This is size M x 3, where M is the number of chordwise points.
        outer_wing_cross_section_mcl = (
            wing_cross_section_xyz_le[outer_wing_cross_section_num, :]
            + outer_wing_cross_section_mcl_local_back
            + outer_wing_cross_section_mcl_local_up
        )

        # Define number of spanwise points and panels. This is based on the inner
        # wing cross section.
        num_spanwise_panels = inner_wing_cross_section.num_spanwise_panels
        num_spanwise_coordinates = num_spanwise_panels + 1

        # Get the spanwise coordinates. This is based on the inner wing cross section.
        if inner_wing_cross_section.spanwise_spacing == "uniform":
            nondim_spanwise_coordinates = np.linspace(
                0,
                1,
                num_spanwise_coordinates,
                endpoint=True,
            )
        elif inner_wing_cross_section.spanwise_spacing == "cosine":
            nondim_spanwise_coordinates = functions.cosspace(
                n_points=num_spanwise_coordinates,
                endpoint=True,
            )
        else:
            raise Exception(
                "This wing cross section's spanwise_spacing variable is invalid!"
            )

        # Make section_mcl_coordinates: M x N x 3 array of mean camberline
        # coordinates. The first index is chordwise point number, second index is
        # spanwise point number, third is the x, y, and z coordinates. M is the
        # number of chordwise points. N is the number of spanwise points. Put a
        # reversed version (from 1 to 0) of the non dimensional spanwise coordinates
        # in a row vector. This is size 1 x N, where N is the number of spanwise
        # points.
        reversed_nondim_spanwise_coordinates_row_vector = np.expand_dims(
            (1 - nondim_spanwise_coordinates), 0
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
        nondim_spanwise_coordinates_row_vector = np.expand_dims(
            nondim_spanwise_coordinates, 0
        )

        # Convert the non dimensional spanwise coordinate row vector (from to 0 to 1)
        # to a matrix. This is size 1 x N x 1, where N is the number of spanwise
        # points.
        nondim_spanwise_coordinates_matrix = np.expand_dims(
            nondim_spanwise_coordinates_row_vector, 2
        )

        """Linearly interpolate between inner and outer wing cross sections. This uses the 
        following equation. 
        
        f(a, b, i) = i * a + (1 - i) * b 
        
        "a" is an N x 3 array of the coordinates points along the outer wing cross 
        section's mean camber line. 
        
        "b" is an N x 3 array of the coordinates of points along the inner wing cross 
        section's mean camber line. 
        
        "i" is a 1D array (or vector) of length M that holds the nondimensionalized 
        spanwise panel spacing from 0 to 1. 
        
        This produces a M x N x 3 array where each slot holds the coordinates of a 
        point on the surface between the inner and outer wing cross sections."""
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

        # Initialize an empty array to hold this wing section's panels. The matrix is
        # size M x N, where M and N are the number of chordwise and spanwise panels.
        wing_section_panels = np.empty(
            (num_chordwise_panels, num_spanwise_panels), dtype=object
        )

        # Loop through the empty panels matrix and create a new panel object in each
        # slot.
        for chordwise_position in range(num_chordwise_panels):
            for spanwise_position in range(num_spanwise_panels):
                wing_section_panels[
                    chordwise_position, spanwise_position
                ] = panel.Panel(
                    front_left_vertex=front_inner_vertices[
                        chordwise_position, spanwise_position
                    ],
                    front_right_vertex=front_outer_vertices[
                        chordwise_position, spanwise_position
                    ],
                    back_left_vertex=back_inner_vertices[
                        chordwise_position, spanwise_position
                    ],
                    back_right_vertex=back_outer_vertices[
                        chordwise_position, spanwise_position
                    ],
                    is_trailing_edge=wing_section_is_trailing_edge[
                        chordwise_position, spanwise_position
                    ],
                    is_leading_edge=wing_section_is_leading_edge[
                        chordwise_position, spanwise_position
                    ],
                )

        # This wing section's panel matrix is stacked horizontally, to the right of the
        # wing's panel matrix.
        panels = np.hstack((panels, wing_section_panels))

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
                    # The inner wing cross section dictates control surface
                    # deflections.
                    hinge_point=inner_wing_cross_section.control_surface_hinge_point,
                )

            # Make the mean camber lines for each wing cross section. First index is
            # point number, second index is xyz.
            inner_wing_cross_section_mcl_nondim = inner_airfoil.get_downsampled_mcl(
                nondim_chordwise_coordinates
            )
            outer_wing_cross_section_mcl_nondim = outer_airfoil.get_downsampled_mcl(
                nondim_chordwise_coordinates
            )

            # Put the inner wing cross section's local up airfoil frame coordinates
            # in a column vector.
            inner_wing_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(
                inner_wing_cross_section_mcl_nondim[:, 1], 1
            )

            # Put the inner wing cross section's local back airfoil frame coordinates
            # in a column vector.
            inner_wing_cross_section_mcl_nondim_local_back_column_vector = (
                np.expand_dims(inner_wing_cross_section_mcl_nondim[:, 0], 1)
            )

            # Put the outer wing cross section's local up airfoil frame coordinates
            # in a column vector.
            outer_wing_cross_section_mcl_nondim_local_up_column_vector = np.expand_dims(
                outer_wing_cross_section_mcl_nondim[:, 1], 1
            )

            # Put the outer wing cross section's local back airfoil frame coordinates
            # in a column vector.
            outer_wing_cross_section_mcl_nondim_local_back_column_vector = (
                np.expand_dims(outer_wing_cross_section_mcl_nondim[:, 0], 1)
            )

            # Convert the inner wing cross section's non dimensional local back
            # airfoil frame coordinates to meshed wing coordinates.
            inner_wing_cross_section_mcl_local_back = (
                wing_cross_section_local_back_norm[inner_wing_cross_section_num, :]
                * inner_wing_cross_section_mcl_nondim_local_back_column_vector
                * wing_cross_section_chord_lengths[inner_wing_cross_section_num]
            )

            # Convert the inner wing cross section's non dimensional local up airfoil
            # frame coordinates to meshed wing coordinates.
            inner_wing_cross_section_mcl_local_up = (
                wing_cross_section_local_up[inner_wing_cross_section_num, :]
                * inner_wing_cross_section_mcl_nondim_local_up_column_vector
                * wing_cross_section_chord_lengths[inner_wing_cross_section_num]
                * wing_cross_section_scaling_factors[inner_wing_cross_section_num]
            )

            # Convert the outer wing cross section's non dimensional local back
            # airfoil frame coordinates to meshed wing coordinates.
            outer_wing_cross_section_mcl_local_back = (
                wing_cross_section_local_back_norm[outer_wing_cross_section_num, :]
                * outer_wing_cross_section_mcl_nondim_local_back_column_vector
                * wing_cross_section_chord_lengths[outer_wing_cross_section_num]
            )

            # Convert the outer wing cross section's non dimensional local up airfoil
            # frame coordinates to meshed wing coordinates.
            outer_wing_cross_section_mcl_local_up = (
                wing_cross_section_local_up[outer_wing_cross_section_num, :]
                * outer_wing_cross_section_mcl_nondim_local_up_column_vector
                * wing_cross_section_chord_lengths[outer_wing_cross_section_num]
                * wing_cross_section_scaling_factors[outer_wing_cross_section_num]
            )

            # Convert the inner wing cross section's meshed wing coordinates to
            # absolute coordinates. This is size M x 3, where M is the number of
            # chordwise points.
            inner_wing_cross_section_mcl = (
                wing_cross_section_xyz_le[inner_wing_cross_section_num, :]
                + inner_wing_cross_section_mcl_local_back
                + inner_wing_cross_section_mcl_local_up
            )

            # Convert the outer wing cross section's meshed wing coordinates to
            # absolute coordinates. This is size M x 3, where M is the number of
            # chordwise points.
            outer_wing_cross_section_mcl = (
                wing_cross_section_xyz_le[outer_wing_cross_section_num, :]
                + outer_wing_cross_section_mcl_local_back
                + outer_wing_cross_section_mcl_local_up
            )

            # Make section_mcl_coordinates: M x N x 3 array of mean camberline
            # coordinates. First index is chordwise point number, second index is
            # spanwise point number, third are the x, y, and z coordinates. M is the
            # number of chordwise points. N is the number of spanwise points. Put a
            # reversed version (from 1 to 0) of the non dimensional spanwise
            # coordinates in a row vector. This is size 1 x N, where N is the number
            # of spanwise points.
            reversed_nondim_spanwise_coordinates_row_vector = np.expand_dims(
                (1 - nondim_spanwise_coordinates), 0
            )

            # Convert the reversed non dimensional spanwise coordinate row vector (
            # from 1 to 0) to a matrix. This is size 1 x N x 1, where N is the number
            # of spanwise points.
            reversed_nondim_spanwise_coordinates_matrix = np.expand_dims(
                reversed_nondim_spanwise_coordinates_row_vector, 2
            )

            # Convert the inner and outer wing cross section's mean camberline coordinates
            # column vectors to matrices. These are size M x 1 x 3, where M is the
            # number of chordwise points.
            inner_wing_cross_section_mcl_matrix = np.expand_dims(
                inner_wing_cross_section_mcl, 1
            )
            outer_wing_cross_section_mcl_matrix = np.expand_dims(
                outer_wing_cross_section_mcl, 1
            )

            # Put the non dimensional spanwise coordinates (from 0 to 1) in a row
            # vector. This is size 1 x N, where N is the number of spanwise points.
            nondim_spanwise_coordinates_row_vector = np.expand_dims(
                nondim_spanwise_coordinates, 0
            )

            # Convert the non dimensional spanwise coordinate row vector (from to 0
            # to 1) to a matrix. This is size 1 x N x 1, where N is the number of
            # spanwise points.
            nondim_spanwise_coordinates_matrix = np.expand_dims(
                nondim_spanwise_coordinates_row_vector, 2
            )

            """Linearly interpolate between inner and outer wing cross sections. This 
            uses the following equation. 

            f(a, b, i) = i * a + (1 - i) * b 

            "a" is an N x 3 array of the coordinates points along the outer wing cross 
            section's mean camber line. 

            "b" is an N x 3 array of the coordinates of points along the inner wing 
            cross section's mean camber line. 

            "i" is a 1D array (or vector) of length M that holds the nondimensionalized 
            spanwise panel spacing from 0 to 1. 

            This produces a M x N x 3 array where each slot holds the coordinates of a 
            point on the surface between the inner and outer wing cross sections."""
            wing_section_mcl_vertices = (
                reversed_nondim_spanwise_coordinates_matrix
                * inner_wing_cross_section_mcl_matrix
                + nondim_spanwise_coordinates_matrix
                * outer_wing_cross_section_mcl_matrix
            )

            # Compute the corners of each panel.
            front_inner_vertices = wing_section_mcl_vertices[:-1, :-1, :]
            front_outer_vertices = wing_section_mcl_vertices[:-1, 1:, :]
            back_inner_vertices = wing_section_mcl_vertices[1:, :-1, :]
            back_outer_vertices = wing_section_mcl_vertices[1:, 1:, :]

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

            # Initialize an empty array to hold this wing section's panels. The
            # matrix is size M x N, where M and N are the number of chordwise and
            # spanwise panels.
            wing_section_panels = np.empty(
                (num_chordwise_panels, num_spanwise_panels), dtype=object
            )

            # Loop through the empty panels matrix and create a new panel object in
            # each slot.
            for chordwise_position in range(num_chordwise_panels):
                for spanwise_position in range(num_spanwise_panels):
                    # Reflect the vertices to create the reflected wing for the
                    # symmetric case.
                    front_inner_vertices_reflected = functions.reflect_over_xz_plane(
                        front_inner_vertices[chordwise_position, spanwise_position]
                    )
                    front_outer_vertices_reflected = functions.reflect_over_xz_plane(
                        front_outer_vertices[chordwise_position, spanwise_position]
                    )
                    back_inner_vertices_reflected = functions.reflect_over_xz_plane(
                        back_inner_vertices[chordwise_position, spanwise_position]
                    )
                    back_outer_vertices_reflected = functions.reflect_over_xz_plane(
                        back_outer_vertices[chordwise_position, spanwise_position]
                    )

                    wing_section_panels[
                        chordwise_position, spanwise_position
                    ] = panel.Panel(
                        front_left_vertex=front_outer_vertices_reflected,
                        front_right_vertex=front_inner_vertices_reflected,
                        back_left_vertex=back_outer_vertices_reflected,
                        back_right_vertex=back_inner_vertices_reflected,
                        is_trailing_edge=wing_section_is_trailing_edge[
                            chordwise_position, spanwise_position
                        ],
                        is_leading_edge=wing_section_is_leading_edge[
                            chordwise_position, spanwise_position
                        ],
                    )

            # This wing section's panel matrix is stacked horizontally, to the left
            # of the wing's panel matrix.
            panels = np.hstack((np.flip(wing_section_panels, axis=1), panels))

    # Iterate through the panels and populate their left and right edge flags. Also
    # populate their local position attributes.
    for chordwise_position in range(wing.num_chordwise_panels):
        for spanwise_position in range(wing.num_spanwise_panels):
            this_panel = panels[chordwise_position, spanwise_position]
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
    wing.panels = panels
