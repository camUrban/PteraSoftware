"""This module contains the Airfoil class.

This module contains the following classes:
    Airfoil: This is a class used to contain the Airfoil of a WingCrossSection.

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
from .. import parameter_validation


# NOTE: I'm in the process of refactoring this class.
class Airfoil:
    """This class is used to contain the Airfoil of a WingCrossSection.

    Citation:
        Adapted from:         geometry.Airfoil in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/27/2020

    This class contains the following public methods:

        add_control_surface: This method returns a version of the Airfoil with a
        control surface added at a given point.

        draw: This method plots this Airfoil using PyPlot.

        get_downsampled_mcl: This method returns the mean camber line (MCL) in a
        downsampled form.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I'm in the process of refactoring this method.
    def __init__(
        self,
        name="NACA0012",
        outline_A_lp=None,
        repanel=True,
        n_points_per_side=400,
    ):
        """This is the initialization method.

        :param name: str, optional

            This is the name of the Airfoil. It should correspond to the name of a
            file the airfoils directory (once converted to lower-case and stripped of
            leading and trailing whitespace) unless you are passing in your own
            array of points using outline_A_lp. The default is "NACA0012".

        :param outline_A_lp: array-like of shape (N,2), optional

            This is an array of the 2D points making up the Airfoil's outline (in
            airfoil axes, relative to the leading poit). If you wish to load
            coordinates from the airfoils directory, leave this as None, which is the
            default. If not, it must be an array-like object of numbers (int or
            float) with shape (N,2). It can be a tuple, list, or numpy array. Values
            are converted to floats internally.it must be a Nx2 numpy array of
            numbers (int or float). Make sure all x-component values are in the range
            [ 0.0, 1.0]. The default value is None.

        :param repanel: bool, optional

            This is the variable that determines whether you would like to repanel
            the airfoil coordinates. This applies to points passed in by the user or
            to those from the airfoils directory. I highly recommended setting this
            to True. The default is True.

        :param n_points_per_side: int or None, optional

            This is number of points to use when repaneling the Airfoil. If repanel
            is True, it must be a positive int. If repanel is False, it must be None.
            The default value is 400.
        """
        self.name = parameter_validation.validate_string(name, "name")

        if outline_A_lp is not None:
            self.outline_A_lp = parameter_validation.validate_2d_vector_array_float(
                outline_A_lp, "coordinates"
            )
            # ToDo: Also check if all x-component values are between 0.0 and 1.0.
        else:
            self._populate_points()

        self.repanel = parameter_validation.validate_boolean(repanel, "repanel")

        if self.repanel:
            if n_points_per_side is None:
                raise ValueError("n_points_per_side must be set if repanel is True")
        else:
            if n_points_per_side is not None:
                raise ValueError("n_points_per_side must be None if repanel is False")
        self.n_points_per_side = parameter_validation.validate_positive_scalar_int(
            n_points_per_side, "n_points_per_side"
        )

        # If repanel is True, repanel the Airfoil.
        if self.repanel:
            # NOTE: I've refactored this up to here. I'm now jumping to this method.
            self._repanel_current_airfoil(n_points_per_side=self.n_points_per_side)

        # Initialize other attributes that will be set by populate_mcl_coordinates.
        self.mcl_coordinates = None
        self.upper_minus_mcl = None
        self.thickness = None

        # Populate the mean camber line attributes.
        self._populate_mcl_coordinates()

    # NOTE: I've haven't yet started refactoring this method.
    def add_control_surface(self, deflection=0.0, hinge_point=0.75):
        """This method returns a version of the Airfoil with a control surface added
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
            (hinge_point, self._get_camber_at_chord_fraction(hinge_point))
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
        outline_A_lp = np.vstack((upper_coordinates, lower_coordinates[1:, :]))

        # Initialize the new, flapped airfoil and return it.
        flapped_airfoil = Airfoil(
            name=self.name + " flapped", outline_A_lp=outline_A_lp, repanel=False
        )
        return flapped_airfoil

    # NOTE: I've haven't yet started refactoring this method.
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

    # NOTE: I've haven't yet started refactoring this method.
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

    # NOTE: I've haven't yet started refactoring this method.
    def _populate_mcl_coordinates(self):
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
        upper = np.flipud(self._upper_coordinates())
        lower = self._lower_coordinates()

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

    # NOTE: I've finished refactoring this method.
    def _populate_points(self):
        """This method populates a variable with the points of the Airfoil's outline
        (in airfoil axes, relative to the leading point).

        The points will be generated if the Airfoil is a NACA 4-series airfoil,
        or loaded from the "airfoils" directory inside "pterasoftware", which is a
        database of dat files containing Airfoil points). NACA 4-series airfoil
        generation is an adaptation of:
        https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: None
        """
        # Sanitize the name input.
        sanitized_name = self.name.lower().strip()

        # Check if the sanitized Airfoil's name matches a name for a NACA 4-series
        # airfoil. If so, generate it.
        if "naca" in sanitized_name:
            naca_number = sanitized_name.split("naca")[1]
            if naca_number.isdigit():
                if len(naca_number) == 4:

                    # Parse the characteristics from the name.
                    max_camber = int(naca_number[0]) * 0.01
                    camber_loc = int(naca_number[1]) * 0.1
                    thickness = int(naca_number[2:]) * 0.01

                    # Set the number of points per side.
                    n_points_per_side = 100

                    # Get the x component of the MCL.
                    mclX_A_lp = functions.cosspace(0, 1, n_points_per_side)

                    # Find the half-thickness of the outline perpendicular to the MCL
                    # (in airfoil axes).
                    halfThickness_A = (
                        5
                        * thickness
                        * (
                            +0.2969 * np.power(mclX_A_lp, 0.5)
                            - 0.1260 * mclX_A_lp
                            - 0.3516 * np.power(mclX_A_lp, 2)
                            + 0.2843 * np.power(mclX_A_lp, 3)
                            - 0.1015 * np.power(mclX_A_lp, 4)
                        )
                    )

                    # Prevent divide by zero errors for airfoils like the NACA 0012.
                    if camber_loc == 0:
                        camber_loc = 0.5

                    # Get the y components of the MCL (in airfoil axes, relative to
                    # the leading point).
                    mclY1_A_lp = (
                        max_camber
                        / camber_loc**2
                        * (
                            2 * camber_loc * mclX_A_lp[mclX_A_lp <= camber_loc]
                            - mclX_A_lp[mclX_A_lp <= camber_loc] ** 2
                        )
                    )
                    mclY2_A_lp = (
                        max_camber
                        / (1 - camber_loc) ** 2
                        * (
                            (1 - 2 * camber_loc)
                            + 2 * camber_loc * mclX_A_lp[mclX_A_lp > camber_loc]
                            - mclX_A_lp[mclX_A_lp > camber_loc] ** 2
                        )
                    )
                    mclY_A_lp = np.hstack((mclY1_A_lp, mclY2_A_lp))

                    # Get the slope of the MCL (in airfoil axes).
                    mclSlope1_A = (
                        2
                        * max_camber
                        / camber_loc**2
                        * (camber_loc - mclX_A_lp[mclX_A_lp <= camber_loc])
                    )
                    mclSlope2_A = (
                        2
                        * max_camber
                        / (1 - camber_loc) ** 2
                        * (camber_loc - mclX_A_lp[mclX_A_lp > camber_loc])
                    )
                    mclSlope_A = np.hstack((mclSlope1_A, mclSlope2_A))

                    # Convert the slope of the MCL to the angle between the airfoil
                    # x-axis and the MCL tangent line.
                    thetaSlope_Ax_to_MCL = np.arctan(mclSlope_A)

                    # Find the upper and lower points of the Airfoil's outline (in
                    # airfoil axes, relative to the leading point) using the MCL
                    # points, the perpendicular half-thickness, and the angle between
                    # the x-axis and the MCL tangent line.
                    upperX_A_lp = mclX_A_lp - halfThickness_A * np.sin(
                        thetaSlope_Ax_to_MCL
                    )
                    lowerX_A_lp = mclX_A_lp + halfThickness_A * np.sin(
                        thetaSlope_Ax_to_MCL
                    )
                    upperY_A_lp = mclY_A_lp + halfThickness_A * np.cos(
                        thetaSlope_Ax_to_MCL
                    )
                    lowerY_A_lp = mclY_A_lp - halfThickness_A * np.cos(
                        thetaSlope_Ax_to_MCL
                    )

                    # Flip upper surface so it's back to front.
                    upperX_A_lp, upperY_A_lp = np.flipud(upperX_A_lp), np.flipud(
                        upperY_A_lp
                    )

                    # Trim one point from lower surface so there's no overlap.
                    lowerX_A_lp, lowerY_A_lp = lowerX_A_lp[1:], lowerY_A_lp[1:]

                    # Combine the points.
                    outlineX_A_lp = np.hstack((upperX_A_lp, lowerX_A_lp))
                    outlineY_A_lp = np.hstack((upperY_A_lp, lowerY_A_lp))

                    # Populate the outline_A_lp attribute and return.
                    self.outline_A_lp = np.column_stack((outlineX_A_lp, outlineY_A_lp))
                    return

        # Try to read from the airfoil directory.
        try:

            # Import the airfoils package as "airfoils".
            airfoils = importlib.import_module(
                name=".airfoils",
                package="pterasoftware",
            )

            # Read the text from the airfoil file.
            raw_text = importlib.resources.read_text(airfoils, sanitized_name + ".dat")

            # Trim the text at the return characters.
            trimmed_text = raw_text[raw_text.find("\n") :]

            # Input the coordinates into a 1D array. This represents the upper and
            # lower points of the Airfoil's outline (in airfoil axes, relative to the
            # leading point).
            outline1D_A_lp = np.fromstring(trimmed_text, sep="\n")

            # Check to make sure the number of elements in the array is even.
            if len(outline1D_A_lp) % 2 != 0:
                raise Exception(
                    "Airfoil file was in airfoil database, but it could not be read correctly."
                )

            # Populate the outline_A_lp attribute and return.
            self.outline_A_lp = np.reshape(outline1D_A_lp, (-1, 2))
            return

        # ToDo: Determine if I should call a particular error here instead of a
        #  general exception.
        # If the Airfoil was not a NACA 4-series and was not found in the database,
        # throw an error.
        except FileNotFoundError:
            raise Exception("Airfoil not in database!")

    # NOTE: I've haven't yet started refactoring this method.
    def _leading_edge_index(self):
        """Returns the index of the leading edge point.

        :return leading_edge_index: int
            This is the index of the leading edge point.
        """

        # Find the index of the coordinate pair with the minimum value of the x
        # coordinate. This is the leading edge index.
        leading_edge_index = np.argmin(self.coordinates[:, 0])

        # Return the leading edge index.
        return leading_edge_index

    # NOTE: I've haven't yet started refactoring this method.
    def _lower_coordinates(self):
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
        lower_coordinates = self.coordinates[self._leading_edge_index() :, :]

        # Return the lower coordinates.
        return lower_coordinates

    # NOTE: I've haven't yet started refactoring this method.
    def _upper_coordinates(self):
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
        upper_coordinates = self.coordinates[: self._leading_edge_index() + 1, :]

        # Return the upper coordinates.
        return upper_coordinates

    # NOTE: I've haven't yet started refactoring this method.
    def _get_camber_at_chord_fraction(self, chord_fraction):
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

    # NOTE: I've haven't yet started refactoring this method.
    def _repanel_current_airfoil(self, n_points_per_side=100):
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
        upper_original_coordinates = self._upper_coordinates()
        lower_original_coordinates = self._lower_coordinates()

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
