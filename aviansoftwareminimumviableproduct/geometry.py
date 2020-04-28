"""This module contains useful functions that relate to geometry, and the class definitions for different types of
geometries.

This module contains the following classes:
    Airplane: This is a class used to contain airplanes.
    Wing: This is a class used to contain the wings of an airplane.
    WingCrossSection: This class is used to contain the cross sections of the wings of an airplane.
    Airfoil: This class is used to contain the airfoil of a cross section of a wing of an airplane.
    Panel: This class is for the panels of a meshed wing.

This module contains the following exceptions:
    None

This module contains the following functions:
    cosspace: This function is used to create a numpy array containing a specified number of values between a specified
              minimum and maximum value that are spaced via a cosine function.
    reflect_over_xz_plane: This function is used to flip a the y coordinate of a coordinate vector.
    angle_axis_rotation_matrix: This function is used to find the rotation matrix for a given axis and angle.
"""

import numpy as np
import pyvista as pv
import scipy.interpolate as sp_interp


class Airplane:
    """This is a class used to contain airplanes.

    Citation:
        Adapted from:         geometry.Airplane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/23/2020

    This class contains the following public methods:
        set_reference_dimensions_from_wing: This method sets the reference dimensions of the airplane from measurements
                                            obtained from the main wing.
        set_paneling_everywhere: This method sets the chordwise and spanwise paneling everywhere to a specified value.
        draw: This method draws the airplane in a new window.
        is_symmetric: This method checks if the airplane is symmetric about the xz plane.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, name="Untitled", x_ref=0.0, y_ref=0.0, z_ref=0.0, wings=None, s_ref=None, c_ref=None,
                 b_ref=None):
        """This is the initialization method.

        :param name: str, optional
            A sensible name for your airplane. The default is "Untitled".
        :param x_ref: float, optional
            This is the x coordinate of the moment reference point. It should be the x coordinate of the center of
            gravity. The default is 0.0.
        :param y_ref: float, optional
            This is the y coordinate of the moment reference point. It should be the y coordinate of the center of
            gravity. The default is 0.0.
        :param z_ref: float, optional
            This is the z coordinate of the moment reference point. It should be the z coordinate of the center of
            gravity. The default is 0.0.
        :param wings: list of Wing objects, optional
            This is a list of the airplane's wings defined as Wing objects. The default is None, which this method
            converts to an empty list.
        :param s_ref: float, optional if more than one wing is in the wings list.
            This is the reference wetted area. If not set, it populates from first wing object.
        :param c_ref: float, optional if more than one wing is in the wings list.
            This is the reference chord length. If not set, it populates from first wing object.
        :param b_ref: float, optional if more than one wing is in the wings list.
            This is the reference span. If not set, it populates from first wing object.
        """

        # Initialize the name and the moment reference point.
        self.name = name
        self.xyz_ref = np.array([x_ref, y_ref, z_ref])

        # If wings was passed as None, set wings to an empty list.
        if wings is None:
            wings = []
        self.wings = wings

        # If the the wing list is not empty, set the wing reference dimensions to be the main wing's reference
        # dimensions.
        if len(self.wings) > 0:
            self.set_reference_dimensions_from_main_wing()

        # If any of the passed reference dimensions are not None, set that reference dimension to be what was passed.
        if s_ref is not None:
            self.s_ref = s_ref
        if c_ref is not None:
            self.c_ref = c_ref
        if b_ref is not None:
            self.b_ref = b_ref

        # Check that everything was set right.
        assert self.xyz_ref is not None
        assert self.s_ref is not None
        assert self.c_ref is not None
        assert self.b_ref is not None

    def set_reference_dimensions_from_main_wing(self):
        """This method sets the reference dimensions of the airplane from measurements obtained from the main wing.

        This method assumes the main wing to be the first wing in the wings list passed by the user.

        :return: None
        """

        # Define the main wing to be the first wing in the wings list.
        main_wing = self.wings[0]

        # Set the objects reference dimension attributes to be the reference dimension attributes of the main wing.
        # These attributes are calculated via methods in the Wing class.
        self.s_ref = main_wing.wetted_area()
        self.b_ref = main_wing.span()
        self.c_ref = main_wing.wetted_area() / main_wing.span()

    def set_paneling_everywhere(self, num_chordwise_panels, num_spanwise_panels):
        """This method sets the chordwise and spanwise paneling everywhere to a specified value.

        This is a quick way for the user to change the fidelity of the problem solution.

        :param num_chordwise_panels: int
            This parameter is the number of chordwise panels for all of the airplane's wings.
        :param num_spanwise_panels: int
            This parameter is the number of spanwise panels for all of the airplane's wings.
        :return: None
        """

        # Iterate through all the airplane's wings.
        for wing in self.wings:
            # Modify each wing's number of chordwise panels.
            wing.num_chordwise_panels = num_chordwise_panels
            # Iterate through each wing's cross sections.
            for cross_section in wing.cross_sections:
                cross_section.num_spanwise_panels = num_spanwise_panels

    def draw(self):
        """This method draws the airplane in a new window.

        This method uses the PyVista toolkit for visualization.

        :return: None
        """

        # Initialize variables in PyVista's PolyData format.
        vertices = np.empty((0, 3))
        faces = np.empty(0)

        # Iterate through the airplane's wings.
        for wing in self.wings:
            # Initialize variables for each wing in PyVista's PolyData format.
            wing_vertices = np.empty((0, 3))
            wing_quad_faces = np.empty((0, 5))
            # Iterate through the wing's cross sections, excluding the last cross section.
            for i in range(len(wing.cross_sections) - 1):
                is_second_to_last_section = i == len(wing.cross_sections) - 2

                # Define variables that are the current cross section's leading edge and trailing edge coordinates.
                le_start = wing.cross_sections[i].xyz_le + wing.xyz_le
                te_start = wing.cross_sections[i].xyz_te() + wing.xyz_le

                # Vertically stack this wing's vertices (that have been collected from previous cross sections) with
                # this cross section's leading edge and trailing edge coordinates.
                wing_vertices = np.vstack((wing_vertices, le_start, te_start))

                # Vertically stack this wing's face indices (that have been collected from previous cross sections) with
                # this cross section's face indices.
                wing_quad_faces = np.vstack((
                    wing_quad_faces,
                    np.expand_dims(np.array([4, 2 * i + 0, 2 * i + 1, 2 * i + 3, 2 * i + 2]), axis=0)
                ))

                # If the current cross section is the second to last cross section, then define the last cross section's
                # leading and trailing edge coordinates to be the last vertices, and stack them below the rest of the
                # wing's vertices.
                if is_second_to_last_section:
                    le_end = wing.xsecs[i + 1].xyz_le + wing.xyz_le
                    te_end = wing.xsecs[i + 1].xyz_te() + wing.xyz_le
                    wing_vertices = np.vstack((wing_vertices, le_end, te_end))

            # Initialize a variable to hold the number of rows in the vertices numpy array. This is equivalent to the
            # number of wings already iterated through.
            vertices_starting_index = len(vertices)
            # Make a copy of the wing_quad_faces numpy array.
            wing_quad_faces_reformatted = np.ndarray.copy(wing_quad_faces)

            # Modify our copy of the wing_quad_faces numpy array. For all the rows, and for all columns except the
            # first, add the number of wings already iterated through to the value in the wing_quad_faces_reformatted.
            wing_quad_faces_reformatted[:, 1:] = wing_quad_faces[:, 1:] + vertices_starting_index
            # In "C-like" order, read unravel the wing_quad_faces_reformatted numpy array into a 1D numpy array. This is
            # the format PyVista's PolyData object takes in faces.
            wing_quad_faces_reformatted = np.reshape(wing_quad_faces_reformatted, (-1), order='C')

            # Stack this wing's vertices with the vertices collected from other wings.
            vertices = np.vstack((vertices, wing_vertices))
            # Stack this wing's faces with the faces collected from other wings.
            faces = np.hstack((faces, wing_quad_faces_reformatted))

            # If the wing is symmetric, repeat this process to create vertices and faces of the wing reflected over the
            # xy plane.
            if wing.symmetric:
                vertices_starting_index = len(vertices)
                wing_vertices = reflect_over_xz_plane(wing_vertices)
                wing_quad_faces_reformatted = np.ndarray.copy(wing_quad_faces)
                wing_quad_faces_reformatted[:, 1:] = wing_quad_faces[:, 1:] + vertices_starting_index
                wing_quad_faces_reformatted = np.reshape(wing_quad_faces_reformatted, (-1), order='C')
                vertices = np.vstack((vertices, wing_vertices))
                faces = np.hstack((faces, wing_quad_faces_reformatted))

        # Initialize the PyVista plotter.
        plotter = pv.Plotter()

        # Add the wing surfaces to the plotter object.
        wing_surfaces = pv.PolyData(vertices, faces)
        plotter.add_mesh(wing_surfaces, color='#7EFC8F', show_edges=True, smooth_shading=True)

        # Add the moment reference point to the plotter object.
        xyz_ref = pv.PolyData(self.xyz_ref)
        plotter.add_points(xyz_ref, color='#50C7C7', point_size=10)

        # Format and show the visualization.
        plotter.show_grid(color='#444444')
        plotter.set_background(color="black")
        plotter.show(cpos=(-1, -1, 1), full_screen=False)


class Wing:
    """This is a class used to contain the wings of an airplane.

    If the wing is symmetric across the XZ plane, just define the right half and supply "symmetric=True" in
    the constructor. If the wing is not symmetric across the XZ plane, just define the wing.

    Citation:
        Adapted from:         geometry.Wing in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/24/2020

    This class contains the following public methods:
        wetted_area: This method calculates the wetted area of the wing.
        span: This method calculates the span of the wing.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, name="Untitled Wing", x_le=0.0, y_le=0.0, z_le=0.0, cross_sections=None, symmetric=False,
                 num_chordwise_panels=8, chordwise_spacing="cosine"):
        """This is the initialization method.

        :param name: str, optional
            This is a sensible name for the wing. The default is "Untitled Wing".
        :param x_le: float, optional
            This is the x coordinate of the leading edge of the wing, relative to the airplane's reference point. The
            default is 0.0.
        :param y_le: float, optional
            This is the y coordinate of the leading edge of the wing, relative to the airplane's reference point. The
            default is 0.0.
        :param z_le: float, optional
            This is the z coordinate of the leading edge of the wing, relative to the airplane's reference point. The
            default is 0.0.
        :param cross_sections: list of CrossSection objects, optional
            This is a list of WingCrossSection objects, that represent the wing's cross sections. The default is None.
        :param symmetric: bool, optional
            Set this to true if the wing is across the xz plane. Set it to false if not. The default is false.
        :param num_chordwise_panels: int, optional
            This is the number of chordwise panels to be used on this wing. The default is 8.
        :param chordwise_spacing: str, optional
            This is the type of spacing between the wing's chordwise panels. It can be set to "cosine" or "uniform".
            Cosine is highly recommended. The default is cosine.
        """

        # Initialize the name and the position of the wing's leading edge.
        self.name = name
        self.xyz_le = np.array([x_le, y_le, z_le])

        # If cross_sections is set to None, set it to an empty list.
        if cross_sections is None:
            cross_sections = []

        # Initialize the other attributes.
        self.cross_sections = cross_sections
        self.symmetric = symmetric
        self.chordwise_panels = num_chordwise_panels
        self.chordwise_spacing = chordwise_spacing

    def wetted_area(self):
        """This method calculates the wetted area of the wing.

        This method neglects the effects of camber and twist.

        :return wetted_area: float
            This method returns the wetted area of the wing.
        """

        # Initialize the wetted area of the wing.
        wetted_area = 0
        # Iterate through all except the last of the wing's cross sections.
        for i in range(len(self.cross_sections) - 1):
            # Calculate the average chord between the current and the next cross section.
            chord_eff = (self.cross_sections[i].chord + self.cross_sections[i + 1].chord) / 2
            # Calculate the location of the trailing edge for the current and next cross section.
            this_xyz_te = self.cross_sections[i].xyz_te()
            that_xyz_te = self.cross_sections[i + 1].xyz_te()
            # Calculate the effective span of the leading edge and trailing edge.
            span_le_eff = np.sqrt(
                (self.cross_sections[i].xyz_le[1] - self.cross_sections[i + 1].xyz_le[1]) ** 2 +
                (self.cross_sections[i].xyz_le[2] - self.cross_sections[i + 1].xyz_le[2]) ** 2
            )
            span_te_eff = np.sqrt(
                (this_xyz_te[1] - that_xyz_te[1]) ** 2 +
                (this_xyz_te[2] - that_xyz_te[2]) ** 2
            )
            # Average the effective leading and trailing edge spans to calculate the effective overall span.
            span_eff = (span_le_eff + span_te_eff) / 2
            # Calculate the wetted area of this section and add it to the total wetted area.
            wetted_area += chord_eff * span_eff
        # If the wing is symmetric, multiply the wetted area by two.
        if self.symmetric:
            wetted_area *= 2
        # Return the wetted area.
        return wetted_area

    def span(self):
        """This method calculates the span of the wing.

        :return span: float
            This method returns the span of the wing.
        """

        # Calculate the span (y-distance between the root and the tip) of the entire wing.
        span = self.cross_sections[-1].xyz_le[1] - self.cross_sections[0].xyz_le[1]

        # If the wing is symmetric, multiply the span by two.
        if self.symmetric:
            span *= 2

        # Return the span.
        return span


class WingCrossSection:
    """This class is used to contain the cross sections of the wings of an airplane.

    Citation:
        Adapted from:         geometry.WingXSec in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/26/2020

    This class contains the following public methods:
        xyz_te: This method calculates the coordinates of the trailing edge of the cross section.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, x_le=0.0, y_le=0.0, z_le=0.0, chord=0.0, twist=0.0, airfoil=None,
                 control_surface_type="symmetric", control_surface_hinge_point=0.75,
                 control_surface_deflection=0, spanwise_panels=8, spanwise_spacing="cosine"):
        """This is the initialization method.

        :param x_le: float, optional
            This is the x coordinate of the leading edge of the cross section relative to the wing's datum. The default
            value is 0.0.
        :param y_le: float, optional
            This is the y coordinate of the leading edge of the cross section relative to the wing's datum. The default
            value is 0.0.
        :param z_le: float, optional
            This is the z coordinate of the leading edge of the cross section relative to the wing's datum. The default
            value is 0.0.
        :param chord: float, optional
            This is the chord of the wing at this cross section. The default value is 0.0.
        :param twist: float, optional
            This is the twist of the cross section about the leading edge in degrees. The default value is 0.0.
        :param airfoil: Airfoil, optional
            This is the airfoil to be used at this cross section. The default value is None.
        :param control_surface_type: str, optional
            This is type of control surfaces for this cross section. It can be "symmetric" or "asymmetric". An example
            of symmetric control surfaces are flaps. An example of asymmetric control surfaces are ailerons. The default
            value is "symmetric".
        :param control_surface_hinge_point: float, optional
            This is the The location of the control surface hinge from the leading edge as a fraction of chord. The
            default value is 0.75.
        :param control_surface_deflection: float, optional
            This is the Control deflection in degrees. Deflection downwards is positive. The default value is 0.0.
        :param spanwise_panels: int, optional
            This is the number of spanwise panels to be used between this cross section and the next one. The default
            value is 8.
        :param spanwise_spacing: str, optional
            This is the Can be 'cosine' or 'uniform'. Highly recommended to be cosine. The default value is
        """

        # Initialize all the class attributes.
        self.x_le = x_le
        self.y_le = y_le
        self.z_le = z_le
        self.chord = chord
        self.twist = twist
        self.airfoil = airfoil
        self.control_surface_type = control_surface_type
        self.control_surface_hinge_point = control_surface_hinge_point
        self.control_surface_deflection = control_surface_deflection
        self.spanwise_panels = spanwise_panels
        self.spanwise_spacing = spanwise_spacing
        self.xyz_le = np.array([x_le, y_le, z_le])

    def xyz_te(self):
        """This method calculates the coordinates of the trailing edge of the cross section.

        :return xyz_te: numpy array
            This is a 1D numpy array that contains the coordinates of the cross section's trailing edge.
        """

        # Find the rotation matrix given the cross section's twist.
        rot = angle_axis_rotation_matrix(self.twist * np.pi / 180, np.array([0, 1, 0]))

        # Use the rotation matrix and the leading edge coordinates to calculate the trailing edge coordinates.
        xyz_te = self.xyz_le + rot @ np.array([self.chord, 0, 0])

        # Return the 1D numpy array that contains the trailing edge's coordinates.
        return xyz_te


class Airfoil:
    """This class is used to contain the airfoil of a cross section of a wing of an airplane.

    Citation:
        Adapted from:         geometry.Airfoil in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/27/2020

    This class contains the following public methods:
        populate_coordinates: This method populates a variable with the coordinates of the airfoil.
        populate_mcl_coordinates: This method creates a list of the airfoil's mean camber line coordinates. It also
                                  creates two lists of the vectors needed to go from the mcl coordinates to the upper
                                  and lower surfaces. It also creates list of the thicknesses at the x coordinates along
                                  the mean camber line.
        leading_edge_index: This method returns the index of the point along the leading edge.
        lower_coordinates: This method returns a matrix of x and y coordinates that describe the lower surface of the
                           airfoil.
        upper_coordinates: This method returns a matrix of x and y coordinates that describe the upper surface of the
                           airfoil.
        get_downsampled_mcl: This method returns the mean camber line in a downsampled form.
        get_camber_at_chord_fraction: This method returns the camber of the airfoil at a given fraction of the chord.
        repanel_current_airfoil: This method returns a repaneled version of the airfoil with cosine-spaced coordinates
                                 on the upper and lower surfaces.
        add_control_surface: This method returns a version of the airfoil with a control surface added at a given point.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, name="Untitled Airfoil", coordinates=None, repanel=True, n_points_per_side=400):
        """This is the initialization method.

        :param name: str, optional
            This is the name of the airfoil. It should correspond to the name in the airfoils directory unless you are
            passing in your own coordinates. The default is "Untitled Airfoil".
        :param coordinates: numpy array, optional
            This is a N x 2 numpy array of the airfoil's coordinates, where N is the number of coordinates. Treat this
            as an immutable, don't edit directly after initialization. If you wish to load coordinates from the airfoil
            directory, leave this as None. The default is None. Make sure that any airfoil coordinates used range in x
            from 0 to 1.
        :param repanel: bool, optional
            This is the variable that determines whether or not you would like to repanel the airfoil coordinates. This
            applies to coordinates passed in by the user or to the directory coordinates. It is highly recommended to
            set this to True. The default is True.
        :param n_points_per_side: int, optional
            This is number of points to use when repaneling the airfoil. It is ignored if the repanel is False. The
            default is 400.
        """

        # Initialize the airfoil name.
        self.name = name

        # Check if the user supplied coordinates.
        if coordinates is not None:
            self.coordinates = coordinates
        else:
            # If not, populate the coordinates from the directory.
            self.populate_coordinates()  # populates self.coordinates

        # Check that the coordinates have been set.
        assert hasattr(self, 'coordinates')

        # Initialize other attributes.
        self.repanel = repanel
        self.mcl_coordinates = None
        self.upper_minus_mcl = None
        self.thickness = None

        # If repanel is True, repanel the airfoil.
        if self.repanel:
            self.repanel_current_airfoil(n_points_per_side=n_points_per_side)

        # Populate the mean camber line attributes.
        self.populate_mcl_coordinates()

    def populate_coordinates(self):
        """This method populates a variable with the coordinates of the airfoil.

        The airfoil coordinates will either be generated, if the airfoil is a NACA 4-series airfoil, or loaded from the
        the airfoil database (a folder named "airfoils" in this directory, that contains a library of dat files for
        airfoil coordinates). NACA 4-series airfoil generation is an adaptation of:
        https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil.

        :return: None
        """

        # Sanitize the name input.
        name = self.name.lower().strip()

        # Check if the airfoil name is a NACA 4-series airfoil. If so, generate it.
        if "naca" in name:
            naca_number = name.split("naca")[1]
            if naca_number.isdigit():
                if len(naca_number) == 4:

                    # Parse the characteristics from the name.
                    max_camber = int(naca_number[0]) * 0.01
                    camber_loc = int(naca_number[1]) * 0.1
                    thickness = int(naca_number[2:]) * 0.01

                    # Set the number of points per side.
                    n_points_per_side = 100

                    # Make uncambered coordinates.
                    # Generate cosine-spaced points.
                    x_t = cosspace(n_points=n_points_per_side)
                    y_t = 5 * thickness * (
                            + 0.2969 * np.power(x_t, 0.5)
                            - 0.1260 * x_t
                            - 0.3516 * np.power(x_t, 2)
                            + 0.2843 * np.power(x_t, 3)
                            - 0.1015 * np.power(x_t, 4)
                    )

                    # Prevent divide by zero errors for airfoils like the NACA 0012.
                    if camber_loc == 0:
                        camber_loc = 0.5

                    # Get camber.
                    y_c_piece1 = max_camber / camber_loc ** 2 * (
                            2 * camber_loc * x_t[x_t <= camber_loc]
                            - x_t[x_t <= camber_loc] ** 2
                    )
                    y_c_piece2 = max_camber / (1 - camber_loc) ** 2 * (
                            (1 - 2 * camber_loc) +
                            2 * camber_loc * x_t[x_t > camber_loc]
                            - x_t[x_t > camber_loc] ** 2
                    )
                    y_c = np.hstack((y_c_piece1, y_c_piece2))

                    # Get camber slope.
                    first_piece_slope = 2 * max_camber / camber_loc ** 2 * (
                            camber_loc - x_t[x_t <= camber_loc]
                    )
                    second_piece_slope = 2 * max_camber / (1 - camber_loc) ** 2 * (
                            camber_loc - x_t[x_t > camber_loc]
                    )
                    slope = np.hstack((first_piece_slope, second_piece_slope))
                    theta = np.arctan(slope)

                    # Combine everything.
                    x_U = x_t - y_t * np.sin(theta)
                    x_L = x_t + y_t * np.sin(theta)
                    y_U = y_c + y_t * np.cos(theta)
                    y_L = y_c - y_t * np.cos(theta)

                    # Flip upper surface so it's back to front.
                    x_U, y_U = np.flipud(x_U), np.flipud(y_U)

                    # Trim 1 point from lower surface so there's no overlap.
                    x_L, y_L = x_L[1:], y_L[1:]

                    # Combine and format the coordinates.
                    x = np.hstack((x_U, x_L))
                    y = np.hstack((y_U, y_L))
                    coordinates = np.column_stack((x, y))

                    # Populate the coordinates attribute and return.
                    self.coordinates = coordinates
                    return
                else:
                    # If the airfoil is a NACA airfoil but not a NACA 4-series, print an error message.
                    print("Unfortunately, only 4-series NACA airfoils can be generated at this time.")

        # Try to read from the airfoil directory.
        try:
            import importlib.resources

            # Import the airfoil directory.
            from . import airfoils

            # Read the text from the airfoil file.
            raw_text = importlib.resources.read_text(str(airfoils), name + '.dat')

            # Trim the text at the return characters.
            trimmed_text = raw_text[raw_text.find('\n'):]

            # Input the coordinates into a 1D numpy array.
            coordinates1D = np.fromstring(trimmed_text, sep='\n')

            # Check to make sure the number of elements in the numpy array is even.
            assert len(
                coordinates1D) % 2 == 0, 'File was found in airfoil database, but it could not be read correctly!'

            # Reshape the 1D coordinates numpy array into a N x 2 numpy array, where N is the number of rows.
            coordinates = np.reshape(coordinates1D, (-1, 2))

            # Populate the coordinates attribute and return.
            self.coordinates = coordinates
            return

        except FileNotFoundError:
            # If the airfoil was not a NACA 4-series and was not found in the database, throw an error.
            print("File was not found in airfoil database!")

    def populate_mcl_coordinates(self):
        """This method creates a list of the airfoil's mean camber line coordinates. It also creates two lists of the
        vectors needed to go from the mcl coordinates to the upper and lower surfaces. It also creates list of the
        thicknesses at the x coordinates along the mean camber line.

        All vectors are listed from the leading edge to the trailing edge of the airfoil.

        :return: None
        """

        # Get the upper and lower coordinates. Flip the upper coordinates so that it is ordered from the leading edge to
        # the trailing edge.
        upper = np.flipud(self.upper_coordinates())
        lower = self.lower_coordinates()

        # Calculate the approximate mean camber line and populate the class attribute.
        mcl_coordinates = (upper + lower) / 2
        self.mcl_coordinates = mcl_coordinates

        # Find the vectors from each mean camber line coordinate to its upper coordinate.
        self.upper_minus_mcl = upper - self.mcl_coordinates

        # Create a list of values that are the thickness of the airfoil at each mean camber line.
        thickness = np.sqrt(
            np.sum(
                np.power(self.upper_minus_mcl, 2),
                axis=1
            )
        ) * 2

        # Populate the class attribute with the thicknesses at their associated x coordinates.
        self.thickness = np.column_stack((self.mcl_coordinates[:, 0], thickness))

    def leading_edge_index(self):
        """Returns the index of the leading edge point.

        :return leading_edge_index: int
            This is the index of the leading edge point.
        """

        # Find the index of the coordinate pair with the minimum value of the x coordinate. This is the leading edge
        # index.
        leading_edge_index = np.argmin(self.coordinates[:, 0])

        # Return the leading edge index.
        return leading_edge_index

    def lower_coordinates(self):
        """This method returns a matrix of x and y coordinates that describe the lower surface of the airfoil.

        The order of the returned matrix is from leading edge to trailing edge. This matrix includes the leading edge
        point so be careful about duplicates if using this method in conjunction with self.upper_coordinates.

        :return lower_coordinates: numpy array
            This is a N x 2 numpy array of x and y coordinates that describe the lower surface of the airfoil, where N
            is the number of points.
        """

        # Find the lower coordinates.
        lower_coordinates = self.coordinates[self.leading_edge_index():, :]

        # Return the lower coordinates.
        return lower_coordinates

    def upper_coordinates(self):
        """This method returns a matrix of x and y coordinates that describe the upper surface of the airfoil.

        The order of the returned matrix is from trailing edge to leading edge. This matrix includes the leading edge
        point so be careful about duplicates if using this method in conjunction with self.lower_coordinates.

        :return upper_coordinates: numpy array
            This is a N x 2 numpy array of x and y coordinates that describe the upper surface of the airfoil, where N
            is the number of points.
        """

        # Find the upper coordinates.
        upper_coordinates = self.coordinates[:self.leading_edge_index() + 1, :]

        # Return the upper coordinates.
        return upper_coordinates

    def get_downsampled_mcl(self, mcl_fractions):
        """This method returns the mean camber line in a downsampled form.

        :param mcl_fractions: 1D numpy array
            This is a 1D numpy array that lists the points along the mean camber line (normalized from 0 to 1) at which
            to return the mean camber line coordinates.
        :return mcl_downsampled: 2D numpy array
            This is a 2D numpy array that contains the coordinates of the downsampled mean camber line.
        """

        mcl = self.mcl_coordinates

        # Find the distances between points along the mean camber line, assuming linear interpolation.
        mcl_distances_between_points = np.sqrt(
            np.power(mcl[:-1, 0] - mcl[1:, 0], 2) +
            np.power(mcl[:-1, 1] - mcl[1:, 1], 2)
        )

        # Create a horizontal 1D numpy array that contains the distance along the mean camber line of each point.
        mcl_distances_cumulative = np.hstack((0, np.cumsum(mcl_distances_between_points)))

        # Normalize the 1D numpy array so that it ranges from 0 to 1.
        mcl_distances_cumulative_normalized = mcl_distances_cumulative / mcl_distances_cumulative[-1]

        # Linearly interpolate to find the x coordinates of the mean camber line at the given mean camber line
        # fractions.
        mcl_downsampled_x = np.interp(
            x=mcl_fractions,
            xp=mcl_distances_cumulative_normalized,
            fp=mcl[:, 0]
        )

        # Linearly interpolate to find the y coordinates of the mean camber line at the given mean camber line
        # fractions.
        mcl_downsampled_y = np.interp(
            x=mcl_fractions,
            xp=mcl_distances_cumulative_normalized,
            fp=mcl[:, 1]
        )

        # Combine the x and y coordinates of the downsampled mean camber line.
        mcl_downsampled = np.column_stack((mcl_downsampled_x, mcl_downsampled_y))

        # Return the coordinates of the downsampled mean camber line.
        return mcl_downsampled

    def get_camber_at_chord_fraction(self, chord_fraction):
        """This method returns the camber of the airfoil at a given fraction of the chord.

        :param chord_fraction: float
            This is a float of the fraction along the chord (normalized from 0 to 1) at which to return the camber.
        :return camber: float
            This is the camber of the airfoil at the requested fraction along the chord.
        """

        # Create a function that interpolates between the x and y coordinates of the mean camber line.
        camber_function = sp_interp.interp1d(
            x=self.mcl_coordinates[:, 0],
            y=self.mcl_coordinates[:, 1],
            copy=False,
            fill_value='extrapolate'
        )

        # Find the value of the camber (the y coordinate) of the airfoil at the requested chord fraction.
        camber = camber_function(chord_fraction)

        # Return the camber of the airfoil at the requested chord fraction.
        return camber

    def repanel_current_airfoil(self, n_points_per_side=100):
        """This method returns a repaneled version of the airfoil with cosine-spaced coordinates on the upper and lower
        surfaces.

        The number of points defining the final airfoil will be (n_points_per_side * 2 - 1), since the leading edge
        point is shared by both the upper and lower surfaces.

        :param n_points_per_side: int, optional
            This is the number of points on the upper and lower surfaces. The default value is 100.
        :return: None
        """

        # Get the upper and lower surface coordinates. These both contain the leading edge point.
        upper_original_coordinates = self.upper_coordinates()
        lower_original_coordinates = self.lower_coordinates()

        # Generate a cosine-spaced list of points from 0 to 1.
        cosine_spaced_x_values = cosspace(n_points=n_points_per_side)

        # Create interpolated functions for the x and y values of the upper and lower surfaces as a function of the
        # chord fractions
        x_upper_func = sp_interp.PchipInterpolator(x=cosine_spaced_x_values, y=upper_original_coordinates[:, 0])
        y_upper_func = sp_interp.PchipInterpolator(x=cosine_spaced_x_values, y=upper_original_coordinates[:, 1])
        x_lower_func = sp_interp.PchipInterpolator(x=cosine_spaced_x_values, y=lower_original_coordinates[:, 0])
        y_lower_func = sp_interp.PchipInterpolator(x=cosine_spaced_x_values, y=lower_original_coordinates[:, 1])

        # Find the x and y coordinates of the upper and lower surfaces at each of the cosine-spaced x values.
        x_coordinates = np.hstack((x_upper_func(cosine_spaced_x_values), x_lower_func(cosine_spaced_x_values)[1:]))
        y_coordinates = np.hstack((y_upper_func(cosine_spaced_x_values), y_lower_func(cosine_spaced_x_values)[1:]))

        # Stack the coordinates together and return them.
        coordinates = np.column_stack((x_coordinates, y_coordinates))
        self.coordinates = coordinates

    def add_control_surface(self, deflection=0.0, hinge_point=0.75):
        """This method returns a version of the airfoil with a control surface added at a given point.

        :param deflection: float, optional
            This is the deflection angle in degrees. Deflection downwards is positive. The default value is 0.0.
        :param hinge_point: float, optional
            This is the location of the hinge as a fraction of chord length. The default value is 0.75.
        :return flapped_airfoil: Airfoil
            This is the new airfoil with the control surface added.
        """

        # Make the rotation matrix for the given angle.
        sin_theta = np.sin(np.radians(-deflection))
        cos_theta = np.cos(np.radians(-deflection))
        rotation_matrix = np.array(
            [[cos_theta, -sin_theta],
             [sin_theta, cos_theta]]
        )

        # Find y coordinate at the hinge point x coordinate and make it a vector.
        hinge_point = np.array((hinge_point, self.get_camber_at_chord_fraction(hinge_point)))

        # Split the airfoil into the sections before and after the hinge.
        split_index = np.where(self.mcl_coordinates[:, 0] > hinge_point[0])[0][0]
        mcl_coordinates_before = self.mcl_coordinates[:split_index, :]
        mcl_coordinates_after = self.mcl_coordinates[split_index:, :]
        upper_minus_mcl_before = self.upper_minus_mcl[:split_index, :]
        upper_minus_mcl_after = self.upper_minus_mcl[split_index:, :]

        # Rotate the mean camber line coordinates and upper minus mean camber line vectors.
        new_mcl_coordinates_after = np.transpose(
            rotation_matrix @ np.transpose(mcl_coordinates_after - hinge_point)) + hinge_point
        new_upper_minus_mcl_after = np.transpose(rotation_matrix @ np.transpose(upper_minus_mcl_after))

        # Assemble the new, flapped airfoil.
        new_mcl_coordinates = np.vstack((mcl_coordinates_before, new_mcl_coordinates_after))
        new_upper_minus_mcl = np.vstack((upper_minus_mcl_before, new_upper_minus_mcl_after))
        upper_coordinates = np.flipud(new_mcl_coordinates + new_upper_minus_mcl)
        lower_coordinates = new_mcl_coordinates - new_upper_minus_mcl
        coordinates = np.vstack((upper_coordinates, lower_coordinates[1:, :]))

        # ToDo: Fix self-intersecting airfoils at high deflections.
        # Initialize the new, flapped airfoil and return it.
        flapped_airfoil = Airfoil(name=self.name + " flapped", coordinates=coordinates, repanel=False)
        return flapped_airfoil


# ToDo: Document and cite this class.
class Panel:
    """

    """
    # ToDo: Document this method.
    def __init__(self, front_left_panel_vertex, front_right_panel_vertex, back_left_panel_vertex,
                 back_right_panel_vertex, is_leading_edge, is_trailing_edge):
        """

        """
        self.front_left_panel_vertex = front_left_panel_vertex
        self.front_right_panel_vertex = front_right_panel_vertex
        self.back_left_panel_vertex = back_left_panel_vertex
        self.back_right_panel_vertex = back_right_panel_vertex
        self.is_leading_edge = is_leading_edge
        self.is_trailing_edge = is_trailing_edge

        self.front_left_vortex_vertex = None
        self.front_right_vortex_vertex = None
        self.back_left_vortex_vertex = None
        self.back_right_vortex_vertex = None
        self.calculate_vortex_vertex_locations()

        self.collocation_point = None
        self.calculate_collocation_point_location()

        self.front_vortex_leg = None
        self.left_vortex_leg = None
        self.back_vortex_leg = None
        self.right_vortex_leg = None
        self.calculate_vortex_leg_vectors()

        self.front_vortex_leg_center = None
        self.left_vortex_leg_center = None
        self.back_vortex_leg_center = None
        self.right_vortex_leg_center = None
        self.calculate_vortex_leg_center_locations()

        self.panel_area = None
        self.panel_normal_direction_at_collocation_point = None
        self.calculate_panel_area_and_normal()

        self.vortex_strength = None

        self.force_on_front_vortex_in_geometry_axes = None
        self.force_on_back_vortex_in_geometry_axes = None
        self.force_on_left_vortex_in_geometry_axes = None
        self.force_on_right_vortex_in_geometry_axes = None
        self.total_force_on_panel_in_geometry_axes = None

        self.total_moment_on_panel_in_geometry_axes = None

        self.normal_force_on_panel = None
        self.pressure_on_panel = None
        self.panel_delta_pressure_coefficient = None

    # ToDo: Document this method.
    def calculate_vortex_vertex_locations(self):
        """

        :return:
        """
        self.front_left_vortex_vertex = self.front_left_panel_vertex + 0.25 * (self.back_left_panel_vertex
                                                                               - self.front_right_panel_vertex)
        self.front_right_vortex_vertex = self.front_right_panel_vertex + 0.25 * (self.back_right_panel_vertex
                                                                                 - self.front_right_panel_vertex)
        self.back_left_vortex_vertex = self.front_left_vortex_vertex + (self.back_left_panel_vertex
                                                                        - self.front_left_panel_vertex)
        self.back_right_vortex_vertex = self.front_right_vortex_vertex + (self.back_right_panel_vertex
                                                                          - self.front_right_panel_vertex)

    # ToDo: Document this method.
    def calculate_collocation_point_location(self):
        """

        :return:
        """
        self.collocation_point = self.front_left_vortex_vertex + 0.5 * (self.back_right_vortex_vertex
                                                                        - self.front_left_vortex_vertex)

    # ToDo: Document this method.
    def calculate_vortex_leg_vectors(self):
        """

        :return:
        """
        self.front_vortex_leg = (self.front_left_vortex_vertex - self.front_right_vortex_vertex)
        self.left_vortex_leg = (self.back_left_vortex_vertex - self.front_left_vortex_vertex)
        self.back_vortex_leg = (self.back_right_vortex_vertex - self.back_left_vortex_vertex)
        self.right_vortex_leg = (self.front_right_vortex_vertex - self.back_right_vortex_vertex)

    # ToDo: Document this method.
    def calculate_vortex_leg_center_locations(self):
        """

        :return:
        """
        self.front_vortex_leg_center = (self.front_right_vortex_vertex + 0.5 * self.front_vortex_leg)
        self.left_vortex_leg_center = (self.front_left_vortex_vertex + 0.5 * self.left_vortex_leg)
        self.back_vortex_leg_center = (self.back_left_vortex_vertex + 0.5 * self.back_vortex_leg)
        self.right_vortex_leg_center = (self.back_right_vortex_vertex + 0.5 * self.right_vortex_leg)

    # ToDo: Document this method.
    def calculate_panel_area_and_normal(self):
        """

        :return:
        """
        # Calculate vortex ring normals and areas via diagonals.
        vortex_ring_first_diagonal = self.front_right_vortex_vertex - self.back_left_vortex_vertex
        vortex_ring_second_diagonal = self.front_left_vortex_vertex - self.back_right_vortex_vertex
        vortex_ring_cross_product = np.cross(vortex_ring_first_diagonal, vortex_ring_second_diagonal)
        vortex_ring_cross_product_magnitude = np.linalg.norm(vortex_ring_cross_product, axis=2)
        self.panel_normal_direction_at_collocation_point = (vortex_ring_cross_product
                                                            / np.expand_dims(vortex_ring_cross_product_magnitude,
                                                                             axis=2))

        # Calculate panel normals and areas via diagonals.
        panel_first_diagonal = self.front_right_panel_vertex - self.back_left_panel_vertex
        panel_second_diagonal = self.front_left_panel_vertex - self.back_right_panel_vertex
        panel_cross_product = np.cross(panel_first_diagonal, panel_second_diagonal)
        panel_cross_product_magnitude = np.linalg.norm(panel_cross_product, axis=2)
        self.panel_area = panel_cross_product_magnitude / 2


# ToDo: Document and cite this class.
def cosspace(minimum=0, maximum=1, n_points=50):
    """

    :param minimum:
    :param maximum:
    :param n_points:
    :return:
    """
    mean = (maximum + minimum) / 2
    amp = (maximum - minimum) / 2

    return mean + amp * np.cos(np.linspace(np.pi, 0, n_points))


# ToDo: Document and cite this function.
def reflect_over_xz_plane(input_vector):
    """
    Takes in a vector or an array and flips the y-coordinates.
    :param input_vector:
    :return:
    """
    output_vector = input_vector
    shape = np.shape(output_vector)
    if len(shape) == 1 and shape[0] == 3:  # Vector of 3 items
        output_vector = output_vector * np.array([1, -1, 1])
    elif len(shape) == 2 and shape[1] == 3:  # 2D Nx3 vector
        output_vector = output_vector * np.array([1, -1, 1])
    elif len(shape) == 3 and shape[2] == 3:  # 3D MxNx3 vector
        output_vector = output_vector * np.array([1, -1, 1])
    else:
        raise Exception("Invalid input for reflect_over_XZ_plane!")

    return output_vector


# ToDo: Document and cite this function.
def angle_axis_rotation_matrix(angle, axis, axis_already_normalized=False):
    """
    Gives the rotation matrix from an angle and an axis.
    An implementation of https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    Inputs:
      * angle: can be one angle or a vector (1d ndarray) of angles. Given in radians.
      * axis: a 1d numpy array of length 3 (p,y,z). Represents the angle.
      * axis_already_normalized: boolean, skips normalization for speed if you flag this true.
    Outputs:
      * If angle is a scalar, returns a 3x3 rotation matrix.
      * If angle is a vector, returns a 3x3xN rotation matrix.
    :param angle:
    :param axis:
    :param axis_already_normalized:
    :return:
    """
    if not axis_already_normalized:
        axis = axis / np.linalg.norm(axis)

    sin_theta = np.sin(angle)
    cos_theta = np.cos(angle)
    cpm = np.array(
        [[0, -axis[2], axis[1]],
         [axis[2], 0, -axis[0]],
         [-axis[1], axis[0], 0]]
    )
    # The cross product matrix of the rotation axis vector
    outer_axis = axis @ np.transpose(axis)

    angle = np.array(angle)  # make sure angle is a ndarray
    if len(angle.shape) == 0:  # is a scalar
        rot_matrix = cos_theta * np.eye(3) + sin_theta * cpm + (1 - cos_theta) * outer_axis
        return rot_matrix
    else:  # angle is assumed to be a 1d ndarray
        rot_matrix = cos_theta * np.expand_dims(np.eye(3), 2) + sin_theta * np.expand_dims(cpm, 2) + (
                1 - cos_theta) * np.expand_dims(outer_axis, 2)
        return rot_matrix
