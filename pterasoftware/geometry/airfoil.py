"""This module contains the Airfoil class.

This module contains the following classes:
    Airfoil: This is a class used to contain airfoils of a WingCrossSection.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import importlib.resources

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sp_interp

from .. import functions


# ToDo: Update this class's nomenclature to comply with the new standards discussed
#  in docs\ANGLE_VECTORS_AND_TRANSFORMATIONS.md and docs\AXES_POINTS_AND_FRAMES.md.
#  It works fine as is, so I'm leaving it alone for now.
class Airfoil:
    """This class is used to contain the airfoils of a WingCrossSection.

    Citation:
        Adapted from:         geometry.Airfoil in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/27/2020

    This class contains the following public methods:
        populate_coordinates: This method populates a variable with the coordinates
        of the airfoil.

        populate_mcl_coordinates: This method creates a list of the airfoil's mean
        camber line coordinates. It also creates two lists of the vectors needed to
        go from the mcl coordinates to the upper and lower surfaces. It also creates
        list of the thicknesses at the x coordinates along the mean camber line.

        leading_edge_index: This method returns the index of the point along the
        leading edge.

        lower_coordinates: This method returns a matrix of x and y coordinates that
        describe the lower surface of the airfoil.

        upper_coordinates: This method returns a matrix of x and y coordinates that
        describe the upper surface of the airfoil.

        get_downsampled_mcl: This method returns the mean camber line in a
        downsampled form.

        get_camber_at_chord_fraction: This method returns the camber of the airfoil
        at a given fraction of the chord.

        repanel_current_airfoil: This method returns a repaneled version of the
        airfoil with cosine-spaced coordinates on the upper and lower surfaces.

        add_control_surface: This method returns a version of the airfoil with a
        control surface added at a given point.

        draw: This method plots this Airfoil's coordinates using PyPlot.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        name="Untitled Airfoil",
        coordinates=None,
        repanel=True,
        n_points_per_side=400,
    ):
        """This is the initialization method.

        :param name: str, optional

            This is the name of the airfoil. It should correspond to the name in the
            airfoils directory unless you are passing in your own coordinates. The
            default is "Untitled Airfoil".

        :param coordinates: (N,2) array-like of numbers, optional

            This is a Nx2 array of the airfoil's coordinates, where N is the number
            of coordinates. Treat this as an immutable, don't edit directly after
            initialization. If you wish to load coordinates from the airfoil
            directory, leave this as None, which is the default. If not, it must be a
            Nx2 numpy array of numbers (int or float). Make sure that any airfoil
            coordinates used range in x from 0 to 1.

        :param repanel: bool, optional

            This is the variable that determines whether you would like to repanel
            the airfoil coordinates. This applies to coordinates passed in by the
            user or to the directory coordinates. I highly recommended setting this
            to True. The default is True.

        :param n_points_per_side: int, optional

            This is number of points to use when repaneling the airfoil. It must be a
            positive int. It is ignored if repanel=False. The default is 400.
        """

        # Initialize the airfoil name.
        self.name = name

        # Check if the user supplied coordinates.
        if coordinates is not None:
            self.coordinates = coordinates
        else:
            # If not, populate the coordinates from the directory.
            self.populate_coordinates()

        # Initialize other user-supplied attributes.
        self.repanel = repanel
        self.n_points_per_side = n_points_per_side

        # Check that the coordinates have been set.
        assert hasattr(self, "coordinates")

        # If repanel is True, repanel the airfoil.
        if self.repanel:
            self.repanel_current_airfoil(n_points_per_side=self.n_points_per_side)

        # Initialize other attributes that will be set by populate_mcl_coordinates.
        self.mcl_coordinates = None
        self.upper_minus_mcl = None
        self.thickness = None

        # Populate the mean camber line attributes.
        self.populate_mcl_coordinates()

    def populate_coordinates(self):
        """This method populates a variable with the coordinates of the airfoil.

        The airfoil coordinates will either be generated, if the airfoil is a NACA
        4-series airfoil, or loaded from the airfoil database (a folder named
        "airfoils" in this directory, that contains a library of dat files for
        airfoil coordinates). NACA 4-series airfoil generation is an adaptation of:
        https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4
        -digit_NACA_airfoil.

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

                    # Make uncambered coordinates and generate cosine-spaced points.
                    x_t = functions.cosspace(0, 1, n_points_per_side)
                    y_t = (
                        5
                        * thickness
                        * (
                            +0.2969 * np.power(x_t, 0.5)
                            - 0.1260 * x_t
                            - 0.3516 * np.power(x_t, 2)
                            + 0.2843 * np.power(x_t, 3)
                            - 0.1015 * np.power(x_t, 4)
                        )
                    )

                    # Prevent divide by zero errors for airfoils like the NACA 0012.
                    if camber_loc == 0:
                        camber_loc = 0.5

                    # Get the camber.
                    y_c_piece1 = (
                        max_camber
                        / camber_loc**2
                        * (
                            2 * camber_loc * x_t[x_t <= camber_loc]
                            - x_t[x_t <= camber_loc] ** 2
                        )
                    )
                    y_c_piece2 = (
                        max_camber
                        / (1 - camber_loc) ** 2
                        * (
                            (1 - 2 * camber_loc)
                            + 2 * camber_loc * x_t[x_t > camber_loc]
                            - x_t[x_t > camber_loc] ** 2
                        )
                    )
                    y_c = np.hstack((y_c_piece1, y_c_piece2))

                    # Get camber slope.
                    first_piece_slope = (
                        2
                        * max_camber
                        / camber_loc**2
                        * (camber_loc - x_t[x_t <= camber_loc])
                    )
                    second_piece_slope = (
                        2
                        * max_camber
                        / (1 - camber_loc) ** 2
                        * (camber_loc - x_t[x_t > camber_loc])
                    )
                    slope = np.hstack((first_piece_slope, second_piece_slope))
                    theta = np.arctan(slope)

                    # Combine everything.
                    x_u = x_t - y_t * np.sin(theta)
                    x_l = x_t + y_t * np.sin(theta)
                    y_u = y_c + y_t * np.cos(theta)
                    y_l = y_c - y_t * np.cos(theta)

                    # Flip upper surface so it's back to front.
                    x_u, y_u = np.flipud(x_u), np.flipud(y_u)

                    # Trim 1 point from lower surface so there's no overlap.
                    x_l, y_l = x_l[1:], y_l[1:]

                    # Combine and format the coordinates.
                    x = np.hstack((x_u, x_l))
                    y = np.hstack((y_u, y_l))
                    coordinates = np.column_stack((x, y))

                    # Populate the coordinates attribute and return.
                    self.coordinates = coordinates
                    return

        # Try to read from the airfoil directory.
        try:

            # Import the airfoils package as "airfoils".
            airfoils = importlib.import_module(
                name=".airfoils",
                package="pterasoftware",
            )

            # Read the text from the airfoil file.
            raw_text = importlib.resources.read_text(airfoils, name + ".dat")

            # Trim the text at the return characters.
            trimmed_text = raw_text[raw_text.find("\n") :]

            # Input the coordinates into a 1D array.
            coordinates_1d = np.fromstring(trimmed_text, sep="\n")

            # Check to make sure the number of elements in the array is even.
            assert len(coordinates_1d) % 2 == 0, (
                "File was found in airfoil database, "
                "but it could not be read correctly."
            )

            # Reshape the 1D coordinates array into an N x 2 array, where N is the
            # number of rows.
            coordinates = np.reshape(coordinates_1d, (-1, 2))

            # Populate the coordinates attribute and return.
            self.coordinates = coordinates
            return

        # If the airfoil was not a NACA 4-series and was not found in the database,
        # throw an error.
        except FileNotFoundError:
            raise Exception("Airfoil not in database!")

    def populate_mcl_coordinates(self):
        """This method creates a list of the airfoil's mean camber line coordinates.
        It also creates two lists of the vectors needed to go from the mcl
        coordinates to the upper and lower surfaces. It also creates list of the
        thicknesses at the x coordinates along the mean camber line.

        All vectors are listed from the leading edge to the trailing edge of the
        airfoil.

        :return: None
        """

        # Get the upper and lower coordinates. Flip the upper coordinates so that it
        # is ordered from the leading edge to the trailing edge.
        upper = np.flipud(self.upper_coordinates())
        lower = self.lower_coordinates()

        # Calculate the approximate mean camber line and populate the class attribute.
        mcl_coordinates = (upper + lower) / 2
        self.mcl_coordinates = mcl_coordinates

        # Find the vectors from each mean camber line coordinate to its upper
        # coordinate.
        self.upper_minus_mcl = upper - self.mcl_coordinates

        # Create a list of values that are the thickness of the airfoil at each mean
        # camber line.
        thickness = np.sqrt(np.sum(np.power(self.upper_minus_mcl, 2), axis=1)) * 2

        # Populate the class attribute with the thicknesses at their associated x
        # coordinates.
        self.thickness = np.column_stack((self.mcl_coordinates[:, 0], thickness))

    def leading_edge_index(self):
        """Returns the index of the leading edge point.

        :return leading_edge_index: int
            This is the index of the leading edge point.
        """

        # Find the index of the coordinate pair with the minimum value of the x
        # coordinate. This is the leading edge index.
        leading_edge_index = np.argmin(self.coordinates[:, 0])

        # Return the leading edge index.
        return leading_edge_index

    def lower_coordinates(self):
        """This method returns a matrix of x and y coordinates that describe the
        lower surface of the airfoil.

        The order of the returned matrix is from leading edge to trailing edge. This
        matrix includes the leading edge point so be careful about duplicates if
        using this method in conjunction with self.upper_coordinates.

        :return lower_coordinates: array
            This is an N x 2 array of x and y coordinates that describe the lower
            surface of the airfoil, where N is the number of points.
        """

        # Find the lower coordinates.
        lower_coordinates = self.coordinates[self.leading_edge_index() :, :]

        # Return the lower coordinates.
        return lower_coordinates

    def upper_coordinates(self):
        """This method returns a matrix of x and y coordinates that describe the
        upper surface of the airfoil.

        The order of the returned matrix is from trailing edge to leading edge. This
        matrix includes the leading edge point so be careful about duplicates if
        using this method in conjunction with self.lower_coordinates.

        :return upper_coordinates: array
            This is an N x 2 array of x and y coordinates that describe the upper
            surface of the airfoil, where N is the number of points.
        """

        # Find the upper coordinates.
        upper_coordinates = self.coordinates[: self.leading_edge_index() + 1, :]

        # Return the upper coordinates.
        return upper_coordinates

    def get_downsampled_mcl(self, mcl_fractions):
        """This method returns the mean camber line in a downsampled form.

        :param mcl_fractions: 1D array
            This is a 1D array that lists the points along the mean camber line (
            normalized from 0 to 1) at which to return the mean camber line
            coordinates.
        :return mcl_downsampled: 2D array
            This is a 2D array that contains the coordinates of the downsampled mean
            camber line.
        """

        mcl = self.mcl_coordinates

        # Find the distances between points along the mean camber line, assuming
        # linear interpolation.
        mcl_distances_between_points = np.sqrt(
            np.power(mcl[:-1, 0] - mcl[1:, 0], 2)
            + np.power(mcl[:-1, 1] - mcl[1:, 1], 2)
        )

        # Create a horizontal 1D array that contains the distance along the mean
        # camber line of each point.
        mcl_distances_cumulative = np.hstack(
            (0, np.cumsum(mcl_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0 to 1.
        mcl_distances_cumulative_normalized = (
            mcl_distances_cumulative / mcl_distances_cumulative[-1]
        )

        # Linearly interpolate to find the x coordinates of the mean camber line at
        # the given mean camber line fractions.
        mcl_downsampled_x = np.interp(
            x=mcl_fractions, xp=mcl_distances_cumulative_normalized, fp=mcl[:, 0]
        )

        # Linearly interpolate to find the y coordinates of the mean camber line at
        # the given mean camber line fractions.
        mcl_downsampled_y = np.interp(
            x=mcl_fractions, xp=mcl_distances_cumulative_normalized, fp=mcl[:, 1]
        )

        # Combine the x and y coordinates of the downsampled mean camber line.
        mcl_downsampled = np.column_stack((mcl_downsampled_x, mcl_downsampled_y))

        # Return the coordinates of the downsampled mean camber line.
        return mcl_downsampled

    def get_camber_at_chord_fraction(self, chord_fraction):
        """This method returns the camber of the airfoil at a given fraction of the
        chord.

        :param chord_fraction: float
            This is a float of the fraction along the chord (normalized from 0 to 1)
            at which to return the camber.
        :return camber: float
            This is the camber of the airfoil at the requested fraction along the
            chord.
        """

        # Create a function that interpolates between the x and y coordinates of the
        # mean camber line.
        camber_function = sp_interp.interp1d(
            x=self.mcl_coordinates[:, 0],
            y=self.mcl_coordinates[:, 1],
            copy=False,
            fill_value="extrapolate",
        )

        # Find the value of the camber (the y coordinate) of the airfoil at the
        # requested chord fraction.
        camber = camber_function(chord_fraction)

        # Return the camber of the airfoil at the requested chord fraction.
        return camber

    def repanel_current_airfoil(self, n_points_per_side=100):
        """This method returns a repaneled version of the airfoil with cosine-spaced
        coordinates on the upper and lower surfaces.

        The number of points defining the final airfoil will be (n_points_per_side *
        2 - 1), since the leading edge point is shared by both the upper and lower
        surfaces.

        :param n_points_per_side: int, optional
            This is the number of points on the upper and lower surfaces. The default
            value is 100.
        :return: None
        """

        # Get the upper and lower surface coordinates. These both contain the leading
        # edge point.
        upper_original_coordinates = self.upper_coordinates()
        lower_original_coordinates = self.lower_coordinates()

        # Generate a cosine-spaced list of points from 0 to 1.
        cosine_spaced_x_values = functions.cosspace(0, 1, n_points_per_side)

        # Create interpolated functions for the x and y values of the upper and lower
        # surfaces as a function of the chord fractions
        upper_func = sp_interp.PchipInterpolator(
            x=np.flip(upper_original_coordinates[:, 0]),
            y=np.flip(upper_original_coordinates[:, 1]),
        )
        lower_func = sp_interp.PchipInterpolator(
            x=lower_original_coordinates[:, 0], y=lower_original_coordinates[:, 1]
        )

        # Find the x and y coordinates of the upper and lower surfaces at each of the
        # cosine-spaced x values.
        x_coordinates = np.hstack(
            (np.flip(cosine_spaced_x_values), cosine_spaced_x_values[1:])
        )
        y_coordinates = np.hstack(
            (
                upper_func(np.flip(cosine_spaced_x_values)),
                lower_func(cosine_spaced_x_values[1:]),
            )
        )

        # Stack the coordinates together and return them.
        coordinates = np.column_stack((x_coordinates, y_coordinates))
        self.coordinates = coordinates

    def add_control_surface(self, deflection=0.0, hinge_point=0.75):
        """This method returns a version of the airfoil with a control surface added
        at a given point.

        :param deflection: float, optional
            This is the deflection angle in degrees. Deflection downwards is
            positive. The default value is 0.0.
        :param hinge_point: float, optional
            This is the location of the hinge as a fraction of chord length. The
            default value is 0.75.
        :return flapped_airfoil: Airfoil
            This is the new airfoil with the control surface added.
        """

        # Ensure that the airfoil's deflection is not too high, which increases the
        # risk of self intersection.
        if deflection > 90 or deflection < -90:
            raise Exception("Invalid value for deflection!")

        # Make the rotation matrix for the given angle.
        sin_theta = np.sin(np.radians(-deflection))
        cos_theta = np.cos(np.radians(-deflection))
        rotation_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])

        # Find y coordinate at the hinge point x coordinate and make it a vector.
        hinge_point = np.array(
            (hinge_point, self.get_camber_at_chord_fraction(hinge_point))
        )

        # Split the airfoil into the sections before and after the hinge.
        split_index = np.where(self.mcl_coordinates[:, 0] > hinge_point[0])[0][0]
        mcl_coordinates_before = self.mcl_coordinates[:split_index, :]
        mcl_coordinates_after = self.mcl_coordinates[split_index:, :]
        upper_minus_mcl_before = self.upper_minus_mcl[:split_index, :]
        upper_minus_mcl_after = self.upper_minus_mcl[split_index:, :]

        # Rotate the mean camber line coordinates and upper minus mean camber line
        # vectors.
        new_mcl_coordinates_after = (
            np.transpose(
                rotation_matrix @ np.transpose(mcl_coordinates_after - hinge_point)
            )
            + hinge_point
        )
        new_upper_minus_mcl_after = np.transpose(
            rotation_matrix @ np.transpose(upper_minus_mcl_after)
        )

        # Assemble the new, flapped airfoil.
        new_mcl_coordinates = np.vstack(
            (mcl_coordinates_before, new_mcl_coordinates_after)
        )
        new_upper_minus_mcl = np.vstack(
            (upper_minus_mcl_before, new_upper_minus_mcl_after)
        )
        upper_coordinates = np.flipud(new_mcl_coordinates + new_upper_minus_mcl)
        lower_coordinates = new_mcl_coordinates - new_upper_minus_mcl
        coordinates = np.vstack((upper_coordinates, lower_coordinates[1:, :]))

        # Initialize the new, flapped airfoil and return it.
        flapped_airfoil = Airfoil(
            name=self.name + " flapped", coordinates=coordinates, repanel=False
        )
        return flapped_airfoil

    def draw(self):
        """This method plots this Airfoil's coordinates using PyPlot.

        :return: None
        """
        x = self.coordinates[:, 0]
        y = self.coordinates[:, 1]
        plt.plot(x, y)
        plt.xlim(0, 1)
        plt.ylim(-0.5, 0.5)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.show()
