"""This module contains useful functions that relate to geometry, and the class definitions for different types of
geometries.

This module contains the following classes:
    Airplane: This is a class used to contain airplanes.
    Wing: This is a class used to contain the wings of an airplane.
    WingCrossSection: This class is used to contain the cross sections of the wings of an airplane.
    Airfoil: This class is used to contain the airfoil of a cross section of a wing of an airplane.
    Panel: This class is used to contain the panels of a wing.

This module contains the following exceptions:
    None

This module contains the following functions:
    cosspace: This function is used to create a ndarray containing a specified number of values between a specified
              minimum and maximum value that are spaced via a cosine function.
    sinspace: This function is used to create a ndarray containing a specified number of values between a specified
              minimum and maximum value that are spaced via a sine function.
    reflect_over_xz_plane: This function is used to flip a the y coordinate of a coordinate vector.
    angle_axis_rotation_matrix: This function is used to find the rotation matrix for a given axis and angle.
    centroid_of_quadrilateral: This function is used to find the centroid of a quadrilateral.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sp_interp
import aviansoftwareminimumviableproduct as asmvp


class Airplane:
    """This is a class used to contain airplanes.

    Citation:
        Adapted from:         geometry.Airplane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/23/2020

    This class contains the following public methods:
        set_reference_dimensions_from_wing: This method sets the reference dimensions of the airplane from measurements
                                            obtained from the main wing.

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
            This is the reference calculate_span. If not set, it populates from first wing object.
        """

        # Initialize the name and the moment reference point.
        self.name = name
        self.x_ref = x_ref
        self.y_ref = y_ref
        self.z_ref = z_ref
        self.xyz_ref = np.array([float(self.x_ref), float(self.y_ref), float(self.z_ref)])

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
            self.s_ref = float(s_ref)
        if c_ref is not None:
            self.c_ref = float(c_ref)
        if b_ref is not None:
            self.b_ref = float(b_ref)

        # Calculate the number of panels in the entire airplane.
        self.num_panels = 0
        for wing in self.wings:
            self.num_panels += wing.num_panels

    def set_reference_dimensions_from_main_wing(self):
        """This method sets the reference dimensions of the airplane from measurements obtained from the main wing.

        This method assumes the main wing to be the first wing in the wings list passed by the user.

        :return: None
        """

        # Define the main wing to be the first wing in the wings list.
        main_wing = self.wings[0]

        # Set the objects reference dimension attributes to be the reference dimension attributes of the main wing.
        # These attributes are calculated via methods in the Wing class.
        self.s_ref = float(main_wing.wetted_area)
        self.b_ref = float(main_wing.span)
        self.c_ref = float(main_wing.wetted_area / main_wing.span)


class Wing:
    """This is a class used to contain the wings of an airplane.

    If the wing is symmetric across the XZ plane, just define the right half and supply "symmetric=True" in
    the constructor. If the wing is not symmetric across the XZ plane, just define the wing.

    Citation:
        Adapted from:         geometry.Wing in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/24/2020

    This class contains the following public methods:
        calculate_wetted_area: This method calculates the wetted area of the wing based on the areas of its panels.
        calculate_span: This method calculates the span of the wing.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, name="Untitled Wing", x_le=0.0, y_le=0.0, z_le=0.0, wing_cross_sections=None, symmetric=False,
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
        :param wing_cross_sections: list of WingCrossSection objects, optional
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
        self.x_le = x_le
        self.y_le = y_le
        self.z_le = z_le
        self.xyz_le = np.array([float(self.x_le), float(self.y_le), float(self.z_le)])

        # If wing_cross_sections is set to None, set it to an empty list.
        if wing_cross_sections is None:
            wing_cross_sections = []

        # Initialize the other attributes.
        self.wing_cross_sections = wing_cross_sections
        self.symmetric = symmetric
        self.num_chordwise_panels = num_chordwise_panels
        self.chordwise_spacing = chordwise_spacing

        # Initialize the the panels attribute. Then mesh the wing, which will populate this attribute.
        self.panels = None
        asmvp.meshing.mesh_wing(self)

        # Find the number of spanwise panels on the wing by adding each cross section's number of spanwise panels.
        # Exclude the last cross section's number of spanwise panels as this is irrelevant. If the wing is symmetric,
        # multiple the summation by two.
        self.num_spanwise_panels = 0
        for cross_section in self.wing_cross_sections[:-1]:
            self.num_spanwise_panels += cross_section.num_spanwise_panels
        if self.symmetric:
            self.num_spanwise_panels *= 2

        # Calculate the number of panels on this wing.
        self.num_panels = self.num_spanwise_panels * self.num_chordwise_panels

        # Initialize and calculate the wing's wetted area. If the wing is symmetrical, this includes the area of the
        # mirrored half.
        self.wetted_area = None
        self.calculate_wetted_area()

        # Initialize and calculate the wing's calculate_span. If the wing is symmetrical, this includes the length of
        # the mirrored half.
        self.span = None
        self.calculate_span()

        # Initialize an empty ndarray to hold this wing's wake ring vortices and its wake ring vortex vertices.
        self.wake_ring_vortex_vertices = np.empty((0, self.num_spanwise_panels + 1, 3))
        self.wake_ring_vortices = np.zeros((0, self.num_spanwise_panels), dtype=object)

    def calculate_wetted_area(self):
        """This method calculates the wetted area of the wing based on the areas of its panels.

        This method also updates the class's wetted area attribute. If the wing is symmetrical, it includes the area of
        the mirrored half.

        :return: None
        """

        wetted_area = 0

        # Iterate through the chordwise and spanwise indices of the panels.
        for chordwise_location in range(self.num_chordwise_panels):
            for spanwise_location in range(self.num_spanwise_panels):
                # Add each panel's area to the total wetted area of the wing.
                wetted_area += self.panels[chordwise_location, spanwise_location].area

        self.wetted_area = wetted_area

    def calculate_span(self):
        """This method calculates the calculate_span of the wing.

        This method also updates the class's span attribute. If the wing is symmetrical, it includes the length of the
        mirrored half.

        :return: None
        """

        # Calculate the span (y-distance between the root and the tip) of the entire wing.
        span = self.wing_cross_sections[-1].xyz_le[1] - self.wing_cross_sections[0].xyz_le[1]

        # If the wing is symmetric, multiply the span by two.
        if self.symmetric:
            span *= 2

        self.span = span


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
                 control_surface_type='symmetric', control_surface_hinge_point=0.75,
                 control_surface_deflection=0, num_spanwise_panels=8, spanwise_spacing='cosine'):
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
        :param num_spanwise_panels: int, optional
            This is the number of spanwise panels to be used between this cross section and the next one. The default
            value is 8.
        :param spanwise_spacing: str, optional
            This is the Can be 'cosine' or 'uniform'. Highly recommended to be cosine. The default value is
        """

        # Initialize all the class attributes.
        self.x_le = float(x_le)
        self.y_le = float(y_le)
        self.z_le = float(z_le)
        self.chord = float(chord)
        self.twist = float(twist)
        self.airfoil = airfoil
        self.control_surface_type = control_surface_type
        self.control_surface_hinge_point = float(control_surface_hinge_point)
        self.control_surface_deflection = float(control_surface_deflection)
        self.num_spanwise_panels = num_spanwise_panels
        self.spanwise_spacing = spanwise_spacing
        self.xyz_le = np.array([x_le, y_le, z_le])

    def xyz_te(self):
        """This method calculates the coordinates of the trailing edge of the cross section.

        :return xyz_te: ndarray
            This is a 1D ndarray that contains the coordinates of the cross section's trailing edge.
        """

        # Find the rotation matrix given the cross section's twist.
        rot = angle_axis_rotation_matrix(self.twist * np.pi / 180, np.array([0, 1, 0]))

        # Use the rotation matrix and the leading edge coordinates to calculate the trailing edge coordinates.
        xyz_te = self.xyz_le + rot @ np.array([self.chord, 0.0, 0.0])

        # Return the 1D ndarray that contains the trailing edge's coordinates.
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
        :param coordinates: ndarray, optional
            This is a N x 2 ndarray of the airfoil's coordinates, where N is the number of coordinates. Treat this
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
                # If the airfoil is a NACA airfoil but not a NACA 4-series, throw an error.
                raise Exception("Unfortunately, only 4-series NACA airfoils can be generated at this time.")

        # Try to read from the airfoil directory.
        try:
            import importlib.resources

            # Import the airfoil directory.
            from . import airfoils

            # Read the text from the airfoil file.
            raw_text = importlib.resources.read_text(str(airfoils), name + '.dat')

            # Trim the text at the return characters.
            trimmed_text = raw_text[raw_text.find('\n'):]

            # Input the coordinates into a 1D ndarray.
            coordinates1D = np.fromstring(trimmed_text, sep='\n')

            # Check to make sure the number of elements in the ndarray is even.
            assert len(
                coordinates1D) % 2 == 0, "File was found in airfoil database, but it could not be read correctly."

            # Reshape the 1D coordinates ndarray into a N x 2 ndarray, where N is the number of rows.
            coordinates = np.reshape(coordinates1D, (-1, 2))

            # Populate the coordinates attribute and return.
            self.coordinates = coordinates
            return

        except FileNotFoundError:
            # If the airfoil was not a NACA 4-series and was not found in the database, throw an error.
            raise Exception("File was not found in airfoil database.")

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

        :return lower_coordinates: ndarray
            This is a N x 2 ndarray of x and y coordinates that describe the lower surface of the airfoil, where N
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

        :return upper_coordinates: ndarray
            This is a N x 2 ndarray of x and y coordinates that describe the upper surface of the airfoil, where N
            is the number of points.
        """

        # Find the upper coordinates.
        upper_coordinates = self.coordinates[:self.leading_edge_index() + 1, :]

        # Return the upper coordinates.
        return upper_coordinates

    def get_downsampled_mcl(self, mcl_fractions):
        """This method returns the mean camber line in a downsampled form.

        :param mcl_fractions: 1D ndarray
            This is a 1D ndarray that lists the points along the mean camber line (normalized from 0 to 1) at which
            to return the mean camber line coordinates.
        :return mcl_downsampled: 2D ndarray
            This is a 2D ndarray that contains the coordinates of the downsampled mean camber line.
        """

        mcl = self.mcl_coordinates

        # Find the distances between points along the mean camber line, assuming linear interpolation.
        mcl_distances_between_points = np.sqrt(
            np.power(mcl[:-1, 0] - mcl[1:, 0], 2) +
            np.power(mcl[:-1, 1] - mcl[1:, 1], 2)
        )

        # Create a horizontal 1D ndarray that contains the distance along the mean camber line of each point.
        mcl_distances_cumulative = np.hstack((0, np.cumsum(mcl_distances_between_points)))

        # Normalize the 1D ndarray so that it ranges from 0 to 1.
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
        upper_func = sp_interp.PchipInterpolator(x=np.flip(upper_original_coordinates[:, 0]),
                                                 y=np.flip(upper_original_coordinates[:, 1]))
        lower_func = sp_interp.PchipInterpolator(x=lower_original_coordinates[:, 0],
                                                 y=lower_original_coordinates[:, 1])

        # Find the x and y coordinates of the upper and lower surfaces at each of the cosine-spaced x values.
        x_coordinates = np.hstack((np.flip(cosine_spaced_x_values), cosine_spaced_x_values[1:]))
        y_coordinates = np.hstack((upper_func(np.flip(cosine_spaced_x_values)), lower_func(cosine_spaced_x_values[1:])))

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

    def draw(self):
        x = self.coordinates[:, 0]
        y = self.coordinates[:, 1]
        plt.plot(x, y)
        plt.xlim(0, 1)
        plt.ylim(-0.5, 0.5)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()


class Panel:
    """This class is used to contain the panels of a wing.

    This class contains the following public methods:
        calculate_collocation_point_location: This method calculates the location of the collocation point.
        calculate_area_and_normal: This method calculates the panel's area and the panel's normal unit vector.
        calculate_normalized_induced_velocity: This method calculates the velocity induced at a point by this panel's
                                               vortices, assuming a unit vortex strength.
        calculate_induced_velocity: This method calculates the velocity induced at a point by this panel's vortices with
                                    their given vortex strengths.
        update_force_moment_and_pressure: This method updates the force, moment, and pressure on this panel.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, front_right_vertex, front_left_vertex, back_left_vertex, back_right_vertex, is_leading_edge,
                 is_trailing_edge):
        """ This is the initialization method.

        :param front_right_vertex: 1D ndarray with three elements
            This is an array containing the x, y, and z coordinates of the panel's front right vertex.
        :param front_left_vertex: 1D ndarray with three elements
            This is an array containing the x, y, and z coordinates of the panel's front left vertex.
        :param back_left_vertex: 1D ndarray with three elements
            This is an array containing the x, y, and z coordinates of the panel's back left vertex.
        :param back_right_vertex: 1D ndarray with three elements
            This is an array containing the x, y, and z coordinates of the panel's back right vertex.
        :param is_leading_edge: bool
            This is true if the panel is the leading edge panel on a wing, and false otherwise.
        :param is_trailing_edge: bool
            This is true if the panel is the trailing edge panel on a wing, and false otherwise.
        """

        # Initialize the attributes.
        self.front_right_vertex = front_right_vertex
        self.front_left_vertex = front_left_vertex
        self.back_left_vertex = back_left_vertex
        self.back_right_vertex = back_right_vertex
        self.is_leading_edge = is_leading_edge
        self.is_trailing_edge = is_trailing_edge

        # Initialize variables to hold the panel's ring and horseshoe vortices. These will be populated by the solver.
        self.ring_vortex = None
        self.horseshoe_vortex = None

        # Initialize a variable to hold the collocation point location and then populate it.
        self.collocation_point = None
        self.calculate_collocation_point_location()

        # Initialize variables to hold the panel area and the panel normal vector at the collocation point. Then
        # populate them.
        self.area = None
        self.normal_direction = None
        self.calculate_area_and_normal()

        # Calculate the center of the panel.
        self.center = asmvp.geometry.centroid_of_quadrilateral(front_right_vertex, front_left_vertex, back_left_vertex,
                                                               back_right_vertex)

        # Calculate the front and back leg lengths, then use them to find and populate the average panel width.
        front_leg_length = np.linalg.norm(front_left_vertex - front_right_vertex)
        back_leg_length = np.linalg.norm(back_right_vertex - back_left_vertex)
        self.width = (front_leg_length + back_leg_length) / 2

        # Initialize two variables that are along the panel's left and right legs at the quarter chord. These points
        # are used for all types of solvers, so we will define them here.
        self.front_right_vortex_vertex = self.front_right_vertex + 0.25 * (self.back_right_vertex
                                                                           - self.front_right_vertex)
        self.front_left_vortex_vertex = self.front_left_vertex + 0.25 * (self.back_left_vertex
                                                                         - self.front_left_vertex)

        # Initialize variables to hold attributes of the panel that will be defined after the solver finds a solution.
        self.near_field_force_geometry_axes = None
        self.near_field_moment_geometry_axes = None
        self.delta_pressure = None

    def calculate_collocation_point_location(self):
        """This method calculates the location of the collocation point.

        The collocation point is at the panel's three quarter chord point.

        :return: None
        """

        # Find the location of points three quarters of the way down the left and right legs of the panel.
        right_three_quarter_chord_mark = self.front_right_vertex + 0.75 * (self.back_right_vertex
                                                                           - self.front_right_vertex)
        left_three_quarter_chord_mark = self.front_left_vertex + 0.75 * (self.back_left_vertex
                                                                         - self.front_left_vertex)

        # Find the vector between the points three quarters of the way down the left and right legs of the panel.
        three_quarter_chord_vector = left_three_quarter_chord_mark - right_three_quarter_chord_mark

        # Find the collocation point, which is halfway between the points three quarters of the way down the left and
        # right legs of the panel. Then populate the class attribute.
        self.collocation_point = right_three_quarter_chord_mark + 0.5 * three_quarter_chord_vector

    def calculate_area_and_normal(self):
        """This method calculates the panel's area and the panel's normal unit vector.

        This method makes the assumption that the panel is planar. This is technically incorrect for wing's with twist
        but is a good approximation for small panels.

        :return: None
        """

        # Calculate panel's normal unit vector and its area via its diagonals.
        first_diagonal = self.front_right_vertex - self.back_left_vertex
        second_diagonal = self.front_left_vertex - self.back_right_vertex
        cross_product = np.cross(first_diagonal, second_diagonal)
        cross_product_magnitude = np.linalg.norm(cross_product)
        self.normal_direction = (cross_product / cross_product_magnitude)
        self.area = cross_product_magnitude / 2

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this panel's vortices, assuming a unit vortex
        strength.

        This method does not include the effect of the panel's wake vortices.

        :param point:  1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        normalized_induced_velocity = np.zeros(3)

        if self.ring_vortex is not None:
            normalized_induced_velocity += self.ring_vortex.calculate_normalized_induced_velocity(point=point)
        if self.horseshoe_vortex is not None:
            normalized_induced_velocity += self.horseshoe_vortex.calculate_normalized_induced_velocity(point=point)

        return normalized_induced_velocity

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this panel's vortices with their given vortex
        strengths.

        This method does not include the effect of the panel's wake vortices.

        :param point: 1D ndarray
            This is a vector containing the x, y, and z coordinates of the point to find the induced velocity at.
        :return: 1D ndarray
            This is a vector containing the x, y, and z components of the induced velocity.
        """

        induced_velocity = np.zeros(3)

        if self.ring_vortex is not None:
            induced_velocity += self.ring_vortex.calculate_induced_velocity(point=point)
        if self.horseshoe_vortex is not None:
            induced_velocity += self.horseshoe_vortex.calculate_induced_velocity(point=point)

        return induced_velocity

    def update_pressure(self):
        """This method updates the pressure across this panel.

        :return: None
        """

        self.delta_pressure = np.dot(self.near_field_force_geometry_axes, self.normal_direction) / self.area


def cosspace(minimum=0.0, maximum=1.0, n_points=50):
    """This function is used to create a ndarray containing a specified number of values between a specified minimum
    and maximum value that are spaced via a cosine function.

    Citation:
        Adapted from:         geometry.cosspace in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    :param minimum: float, optional
        This is the minimum value of the range of numbers you would like spaced. The default is 0.0.
    :param maximum: float, optional
        This is the maximum value of the range of numbers you would like spaced. The default is 1.0.
    :param n_points: int, optional
        This is the number of points to space. The default is 50.
    :return cosine_spaced_points: 1D ndarray
        This is a 1D ndarray of the points, ranging from the minimum to the maximum value (inclusive), spaced via a
        cosine function.
    """

    # Find the mean and the amplitude of the cosine function.
    mean = (maximum + minimum) / 2
    amp = (maximum - minimum) / 2

    # Space the points by applying cosine to the linspace function. Then return the points.
    cosine_spaced_points = mean + amp * np.cos(np.linspace(np.pi, 0, n_points))
    return cosine_spaced_points


def reflect_over_xz_plane(input_vector):
    """This function is used to flip a the y coordinate of a coordinate vector.

    Citation:
        Adapted from:         geometry.reflect_over_xz_plane in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    :param input_vector: ndarray
        This can either be a 1D ndarray of three items, a M x 3 2D ndarray, or a M x N x 3 3D ndarray. N and
        represent arbitrary numbers of rows or columns.
    :return output vector: ndarray
        This is a ndarray with each vertex's y variable flipped.
    """

    # Initialize the output vector.
    output_vector = input_vector

    # Find the shape of the input vector.
    shape = np.shape(output_vector)

    # Characterize the input vector.
    if len(shape) == 1 and shape[0] == 3:
        # The input vector is a 1D ndarray of 3 items. Flip the vertex's y variable.
        output_vector = output_vector * np.array([1, -1, 1])
    elif len(shape) == 2 and shape[1] == 3:
        # The input vector is a 2D ndarray of shape M x 3. Where M is some arbitrary number of rows. Flip each
        # vertex's y variable.
        output_vector = output_vector * np.array([1, -1, 1])
    elif len(shape) == 3 and shape[2] == 3:  # 3D MxNx3 vector
        # The input vector is a 3D ndarray of shape M x N x 3. Where M is some arbitrary number of rows, and N is
        # some arbitrary number of columns. Flip each vertex's y variable.
        output_vector = output_vector * np.array([1, -1, 1])
    else:
        # The input vector is an unacceptable shape. Throw an error.
        raise Exception("Invalid input for reflect_over_xz_plane.")

    # Return the output vector.
    return output_vector


def angle_axis_rotation_matrix(angle, axis, axis_already_normalized=False):
    """This function is used to find the rotation matrix for a given axis and angle.

    Citation:
        Adapted from:         geometry.angle_axis_rotation_matrix in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/28/2020

    For more information on the math behind this method, see:
    https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    :param angle: float or 1D ndarray of 3 angles
        This is the angle. Provide it in radians.
    :param axis: 1D ndarray of 3 elements.
        This is the axis.
    :param axis_already_normalized: bool, optional
        Set this as True if the axis given is already normalized, which will skip internal normalization and speed up
        the method. If not, set as False. The default is False.
    :return rotation matrix: ndarray
        This is the rotation matrix. If the given angle is a scalar, this will be a 3 x 3 ndarray. If the given
        angle is a vector, this will be a 3 x 3 x N ndarray.
    """

    # Normalize the axis is it is not already normalized.
    if not axis_already_normalized:
        axis = axis / np.linalg.norm(axis)

    sin_theta = np.sin(angle)
    cos_theta = np.cos(angle)
    cpm = np.array(
        [[0, -axis[2], axis[1]],
         [axis[2], 0, -axis[0]],
         [-axis[1], axis[0], 0]]
    )

    # Find the cross product matrix of the rotation axis vector.
    outer_axis = axis @ np.transpose(axis)

    # Make sure angle is a ndarray.
    angle = np.array(angle)

    # Check if the angle is a scalar.
    if len(angle.shape) == 0:
        rotation_matrix = cos_theta * np.eye(3) + sin_theta * cpm + (1 - cos_theta) * outer_axis
    else:
        # Otherwise, the angle is assumed to be a 1D ndarray.
        rotation_matrix = cos_theta * np.expand_dims(np.eye(3), 2) + sin_theta * np.expand_dims(cpm, 2) + (
                1 - cos_theta) * np.expand_dims(outer_axis, 2)

    # Return the rotation matrix.
    return rotation_matrix


def centroid_of_quadrilateral(front_left_vertex, front_right_vertex, back_left_vertex, back_right_vertex):
    """This function is used to find the centroid of a quadrilateral.

    :param front_left_vertex: 1D ndarray
        This is an array containing the x, y, and z components of the front left vertex of the quadrilateral.
    :param front_right_vertex: 1D ndarray
        This is an array containing the x, y, and z components of the front right vertex of the quadrilateral.
    :param back_left_vertex: 1D ndarray
        This is an array containing the x, y, and z components of the back left vertex of the quadrilateral.
    :param back_right_vertex: 1D ndarray
        This is an array containing the x, y, and z components of the back right vertex of the quadrilateral.
    :return: 1D ndarray
        This is an array containing the x, y, and z components of the centroid of the quadrilateral.
    """
    x_values = np.hstack((front_left_vertex[0], front_right_vertex[0], back_left_vertex[0], back_right_vertex[0]))
    x_average = np.average(x_values)

    y_values = np.hstack((front_left_vertex[1], front_right_vertex[1], back_left_vertex[1], back_right_vertex[1]))
    y_average = np.average(y_values)

    z_values = np.hstack((front_left_vertex[2], front_right_vertex[2], back_left_vertex[2], back_right_vertex[2]))
    z_average = np.average(z_values)

    centroid = np.array([x_average, y_average, z_average])

    return centroid
