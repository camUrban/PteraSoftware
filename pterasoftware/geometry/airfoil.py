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
from .. import transformations


class Airfoil:
    """This class is used to contain the Airfoil of a WingCrossSection.

    Citation:
        Adapted from:         geometry.Airfoil in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/27/2020

    This class contains the following public methods:

        add_control_surface: This method returns a version of the Airfoil with a
        control surface added at a given point.

        draw: This method plots this Airfoil's outlines and mean camber line (MCL)
        using PyPlot.

        get_resampled_mcl: This method returns an array of points along the mean
        camber line (MCL), resampled from the mcl_A_outline attribute. It is used to
        discretize the MCL for meshing.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        name="NACA0012",
        outline_A_lp=None,
        resample=True,
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

        :param resample: bool, optional

            This is the variable that determines whether you would like to resample
            the points defining the Airfoil's outline. This applies to points passed
            in by the user or to those from the airfoils directory. I highly
            recommended setting this to True. The default is True.

        :param n_points_per_side: int or None, optional

            This is number of points to use when resampling the Airfoil's outline. If
            resample is True, it must be a positive int. If resample is False,
            it must be None. The default value is 400.
        """
        self.name = parameter_validation.validate_string(name, "name")

        if outline_A_lp is not None:
            outline_A_lp = parameter_validation.validate_2d_vector_array_float(
                outline_A_lp, "coordinates"
            )
            if not np.all((outline_A_lp[:, 0] >= 0.0) & (outline_A_lp[:, 0] <= 1.0)):
                raise ValueError(
                    "All x-component values in outline_A_lp must be between 0.0 and "
                    "1.0."
                )
            self.outline_A_lp = outline_A_lp
        else:
            self._populate_outline()

        self.resample = parameter_validation.validate_boolean(resample, "resample")

        if self.resample:
            if n_points_per_side is None:
                raise ValueError("n_points_per_side must be set if resample is True")
        else:
            if n_points_per_side is not None:
                raise ValueError("n_points_per_side must be None if resample is False")
        self.n_points_per_side = parameter_validation.validate_positive_scalar_int(
            n_points_per_side, "n_points_per_side"
        )

        # If resample is True, resample the Airfoil's outline points.
        if self.resample:
            self._resample_outline(self.n_points_per_side)

        # Initialize an attribute for an array of points along the MCL (in airfoil
        # axes, relative to the leading point). It will be set by _populate_mcl.
        self.mcl_A_lp = None
        self._populate_mcl()

    def add_control_surface(self, deflection, hinge_point):
        """This method returns a version of the Airfoil with a control surface added
        at a given point. It is called during meshing.

        :param deflection: number
            This is the control deflection in degrees. Deflection downwards is
            positive. It must be a number (int or float) in the range (-90.0,
            90.0) degrees. Values are converted to floats internally.
        :param hinge_point: float
            This is the location of the hinge as a fraction of chord length. It must
            be in the range (0.0, 1.0).
        :return flapped_airfoil: Airfoil
            This is the new airfoil with the control surface added.
        """
        # Validate the deflection and hinge_point inputs.
        deflection = parameter_validation.validate_scalar_in_range_float(
            deflection,
            "deflection",
            -90.0,
            False,
            90.0,
            False,
        )
        hinge_point = parameter_validation.validate_scalar_in_range_float(
            hinge_point,
            "hinge_point",
            0.0,
            False,
            1.0,
            False,
        )

        # Make the active rotational homogeneous transformation matrix for the given
        # angle.
        rot_T = transformations.generate_rot_T(
            (0, 0, -deflection), passive=False, intrinsic=False, order="zyx"
        )

        # ToDo: Delete these old lines after confirming the the refactoring was
        #  successful.
        # sin_theta = np.sin(np.radians(-deflection))
        # cos_theta = np.cos(np.radians(-deflection))
        # rotation_matrix = np.array([[cos_theta, -sin_theta], [sin_theta, cos_theta]])
        # # Rotate the mean camber line coordinates and upper minus mean camber line
        # # vectors.
        # new_mcl_coordinates_after = (
        #     np.transpose(
        #         rotation_matrix @ np.transpose(mcl_coordinates_after - hingePoint_A_lp)
        #     )
        #     + hingePoint_A_lp
        # )
        # new_upper_minus_mcl_after = np.transpose(
        #     rotation_matrix @ np.transpose(upper_minus_mcl_after)
        # )

        # Find the position of the hinge point on the MCL (in airfoil axes, relative
        # to the leading point)
        hingePoint_A_lp = np.array((hinge_point, self._get_mclY(hinge_point)))

        # Resample the upper outline to have points with the same x-components as the
        # MCL.
        upperFlippedOutline_A_lp = np.flipud(self._upper_outline())
        upper_func = sp_interp.PchipInterpolator(
            x=upperFlippedOutline_A_lp[:, 0],
            y=upperFlippedOutline_A_lp[:, 1],
            extrapolate=False,
        )
        upperFlippedOutline_A_lp = np.column_stack(
            [self.mcl_A_lp[:, 0], upper_func(self.mcl_A_lp[:, 0])]
        )

        # Find the vectors from each mean camber line coordinate to the upper portion
        # of the Airfoil's outline (in airfoil axes)
        mclToUpper_A = upperFlippedOutline_A_lp - self.mcl_A_lp

        # Split the airfoil into the sections before and after the hinge point.
        split_index = np.where(self.mcl_A_lp[:, 0] > hingePoint_A_lp[0])[0][0]
        preHingeMcl_A_lp = self.mcl_A_lp[:split_index, :]
        postHingeMcl_A_lp = self.mcl_A_lp[split_index:, :]
        preHingeMclToUpper_A = mclToUpper_A[:split_index, :]
        postHingeMclToUpper_A = mclToUpper_A[split_index:, :]

        postHingeMcl_Wcs_lp = np.hstack(
            [postHingeMcl_A_lp, np.zeros((postHingeMcl_A_lp.shape[0], 1))]
        )
        postHingeMclToUpper_Wcs = np.hstack(
            [postHingeMclToUpper_A, np.zeros((postHingeMclToUpper_A.shape[0], 1))]
        )

        flappedPostHingeMcl_A_lp = transformations.apply_T_to_vectors(
            rot_T, postHingeMcl_Wcs_lp, has_point=True
        )[:, :2]
        flappedPostHingeMclToUpper_A = transformations.apply_T_to_vectors(
            rot_T, postHingeMclToUpper_Wcs, has_point=False
        )[:, :2]

        # ToDo: Determine if this method is forcing Airfoil's to be symmetrical.
        # Assemble the new flapped Airfoil.
        flappedMcl_A_lp = np.vstack((preHingeMcl_A_lp, flappedPostHingeMcl_A_lp))
        flappedMclToUpper_A = np.vstack(
            (preHingeMclToUpper_A, flappedPostHingeMclToUpper_A)
        )
        flappedUpperOutline_A_lp = np.flipud(flappedMcl_A_lp + flappedMclToUpper_A)
        flappedLowerOutline_A_lp = flappedMcl_A_lp - flappedMclToUpper_A
        flappedOutline_A_lp = np.vstack(
            (flappedUpperOutline_A_lp, flappedLowerOutline_A_lp[1:, :])
        )

        # Initialize and return the new flapped Airfoil.
        flapped_airfoil = Airfoil(
            name=self.name + " flapped",
            outline_A_lp=flappedOutline_A_lp,
            resample=False,
            n_points_per_side=None,
        )
        return flapped_airfoil

    def draw(self):
        """This method plots this Airfoil's outlines and mean camber line (MCL) using
        PyPlot.

        :return: None
        """
        outlineX_A_lp = self.outline_A_lp[:, 0]
        outlineY_A_lp = self.outline_A_lp[:, 1]
        mclX_A_lp = self.mcl_A_lp[:, 0]
        mclY_A_lp = self.mcl_A_lp[:, 1]

        outlineYMin_A_lp = np.min(outlineY_A_lp)
        outlineYMax_A_lp = np.max(outlineY_A_lp)
        outlineYRange_A_lp = outlineYMax_A_lp - outlineYMin_A_lp
        y_padding = 0.1 * outlineYRange_A_lp

        plt.plot(outlineX_A_lp, outlineY_A_lp, "b-")
        plt.plot(mclX_A_lp, mclY_A_lp, "r-")

        plt.xlim(0.0, 1.0)
        plt.ylim(outlineYMin_A_lp - y_padding, outlineYMax_A_lp + y_padding)

        plt.xlabel("x (airfoil axes)")
        plt.ylabel("y (airfoil axes)")
        plt.title(f"Airfoil: {self.name}")
        plt.legend(["Outline", "Mean Camber Line (MCL)"])

        plt.gca().set_aspect("equal", adjustable="box")

        plt.show()

    def get_resampled_mcl(self, mcl_fractions):
        """This method returns an array of points along the mean camber line (MCL),
        resampled from the mcl_A_outline attribute. It is used to discretize the MCL
        for meshing.

        :param mcl_fractions: (N,) array-like of floats

            This is a (N,) array-like object of normalized distances along the MCL (
            from the leading to the trailing edge) at which to return the resampled
            MCL points. It can be a tuple, list, or ndarray. The first value must be
            0.0, the last must be 1.0, and the remaining must be in the range [0.0,
            1.0]. All values must be non-duplicated floats in ascending order.

        :return: (N,2) ndarray of floats

            This is a (N,2) ndarray of floats that contains the positions of the
            resampled MCL points (in airfoil axes, relative to the leading point).
        """
        # Validate the mcl_fractions input parameter.
        mcl_fractions = parameter_validation.validate_nd_vector_float(
            mcl_fractions, "mcl_fractions"
        )
        if len(mcl_fractions) < 2:
            raise ValueError("mcl_fractions must contain at least two values.")
        if not np.isclose(mcl_fractions[0], 0.0):
            raise ValueError("The first value in mcl_fractions must be 0.0.")
        if not np.isclose(mcl_fractions[-1], 1.0):
            raise ValueError("The last value in mcl_fractions must be 1.0.")
        if not np.all((mcl_fractions >= 0.0) & (mcl_fractions <= 1.0)):
            raise ValueError(
                "All values in mcl_fractions must be in the range[0.0, 1.0]."
            )
        if not np.all(np.diff(mcl_fractions) > 0):
            raise ValueError(
                "All values in mcl_fractions must be non-duplicated and in ascending "
                "order."
            )

        # Find the distance between points along the MCL.
        mcl_distances_between_points = np.sqrt(
            np.power(self.mcl_A_lp[:-1, 0] - self.mcl_A_lp[1:, 0], 2)
            + np.power(self.mcl_A_lp[:-1, 1] - self.mcl_A_lp[1:, 1], 2)
        )

        # Create a horizontal 1D array that contains the distance along the MCL of
        # each point on the MCL.
        mcl_distances_cumulative = np.hstack(
            (0, np.cumsum(mcl_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        mcl_distances_cumulative_normalized = (
            mcl_distances_cumulative / mcl_distances_cumulative[-1]
        )

        # Create interpolated functions for MCL's components as a function of
        # fractional distances along the MCL.
        mclX_func = sp_interp.PchipInterpolator(
            x=mcl_distances_cumulative_normalized,
            y=self.mcl_A_lp[:, 0],
            extrapolate=False,
        )
        mclY_func = sp_interp.PchipInterpolator(
            x=mcl_distances_cumulative_normalized,
            y=self.mcl_A_lp[:, 1],
            extrapolate=False,
        )

        return np.column_stack([mclX_func(mcl_fractions), mclY_func(mcl_fractions)])

    def _get_mclY(self, chord_fraction):
        """This method returns the y-component of the Airfoil's MCL (in airfoil axes,
        relative to the leading point) at a given fraction along the chord.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :param chord_fraction: number

            This number (int or float) is the fraction along the chord, from leading
            to trailing point, at which to return the y-component of the MCL (in
            airfoil axes, relative to the leading point). It should be in the range [
            0.0, 1.0].

        :return: float

            This is the y-component of the MCL (in airfoil axes, relative to the
            leading point) at the requested fraction along the chord.
        """
        mclY_func = sp_interp.PchipInterpolator(
            x=self.mcl_A_lp[:, 0],
            y=self.mcl_A_lp[:, 1],
            extrapolate=False,
        )

        return mclY_func(chord_fraction)

    def _leading_edge_index(self):
        """Returns the index of the leading point in the outline_A_lp attribute.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: int
            This is the index of the leading point.
        """
        return np.argmin(self.outline_A_lp[:, 0])

    def _lower_outline(self):
        """This method returns a 2D array of points on the lower portion of the
        Airfoil's outline (in airfoil axes, relative to the leading point).

        The order of the returned points is from leading edge to trailing edge.
        Included is the leading point, so be careful about duplicates if
        using this method in conjunction with _upper_outline.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: (N,2) ndarray of floats

            This is an (N,2) ndarray of floats that describe the position of N points
            on the Airfoil's lower outline (in airfoil axes, relative to the leading
            point).
        """
        return self.outline_A_lp[self._leading_edge_index() :, :]

    def _populate_mcl(self):
        """This method creates a 2D array of points along the Airfoil's MCL (in
        airfoil axes, relative to the leading point), which it uses to set the mcl_A_lp
        attribute. It is in order from the leading point to the trailing point.

        :return: None
        """
        # Split outline_A_lp into upper and lower sections. Flip the upper points so
        # that they are ordered from the leading point to the trailing point.
        upperFlippedOutline_A_lp = np.flipud(self._upper_outline())
        lowerOutline_A_lp = self._lower_outline()

        cosine_spaced_chord_fractions = functions.cosspace(
            0.0, 1.0, self.n_points_per_side
        )

        upper_func = sp_interp.PchipInterpolator(
            x=upperFlippedOutline_A_lp[:, 0],
            y=upperFlippedOutline_A_lp[:, 1],
            extrapolate=True,
        )
        lower_func = sp_interp.PchipInterpolator(
            x=lowerOutline_A_lp[:, 0],
            y=lowerOutline_A_lp[:, 1],
            extrapolate=True,
        )

        upperFlippedOutline_A_lp = upper_func(cosine_spaced_chord_fractions)
        lowerOutline_A_lp = lower_func(cosine_spaced_chord_fractions)

        # Calculate the approximate MCL points (in airfoil axes, relative to the
        # leading point) and set the class attribute.
        self.mcl_A_lp = np.column_stack(
            [
                cosine_spaced_chord_fractions,
                (upperFlippedOutline_A_lp + lowerOutline_A_lp) / 2,
            ]
        )

        # Resample the MCL points using cosine-spaced distances along the MCL.
        self.mcl_A_lp = self.get_resampled_mcl(
            mcl_fractions=cosine_spaced_chord_fractions
        )

    def _populate_outline(self):
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
                    n_points_per_side = 400

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
                    upperOutlineX_A_lp = mclX_A_lp - halfThickness_A * np.sin(
                        thetaSlope_Ax_to_MCL
                    )
                    lowerOutlineX_A_lp = mclX_A_lp + halfThickness_A * np.sin(
                        thetaSlope_Ax_to_MCL
                    )
                    upperOutlineY_A_lp = mclY_A_lp + halfThickness_A * np.cos(
                        thetaSlope_Ax_to_MCL
                    )
                    lowerOutlineY_A_lp = mclY_A_lp - halfThickness_A * np.cos(
                        thetaSlope_Ax_to_MCL
                    )

                    # Flip upper surface so it's back to front.
                    upperOutlineX_A_lp, upperOutlineY_A_lp = np.flipud(
                        upperOutlineX_A_lp
                    ), np.flipud(upperOutlineY_A_lp)

                    # Trim one point from lower surface so there's no overlap.
                    lowerOutlineX_A_lp, lowerOutlineY_A_lp = (
                        lowerOutlineX_A_lp[1:],
                        lowerOutlineY_A_lp[1:],
                    )

                    # Combine the points.
                    outlineX_A_lp = np.hstack((upperOutlineX_A_lp, lowerOutlineX_A_lp))
                    outlineY_A_lp = np.hstack((upperOutlineY_A_lp, lowerOutlineY_A_lp))

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

    def _resample_outline(self, n_points_per_side):
        """This method returns a resampled version of the points on the Airfoil's
        outline (in airfoil axes, relative to the leading point) with cosine-spaced
        points on the upper and lower surfaces.

        The number of points defining the final Airfoil's outline will be (
        n_points_per_side * 2 - 1), since the leading point is shared by both the
        upper and lower surfaces.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :param n_points_per_side: positive int
            This is the number of points on the upper and lower surfaces.
        :return: None
        """
        # Get the upper outline points. This contains the leading point.
        upperOutline_A_lp = self._upper_outline()

        upperFlippedOutlineX_A_lp = np.flip(upperOutline_A_lp[:, 0])
        upperFlippedOutlineY_A_lp = np.flip(upperOutline_A_lp[:, 1])

        # Find the distance between points along the upper flipped original outline.
        upperFlippedOutline_distances_between_points = np.sqrt(
            np.power(
                upperFlippedOutlineX_A_lp[:-1] - upperFlippedOutlineX_A_lp[1:],
                2,
            )
            + np.power(
                upperFlippedOutlineY_A_lp[:-1] - upperFlippedOutlineY_A_lp[1:],
                2,
            )
        )

        # Create a horizontal 1D array that contains the cumulative distance along
        # the upper flipped original outline of each point.
        upperFlippedOutline_distances_cumulative = np.hstack(
            (0, np.cumsum(upperFlippedOutline_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        upperFlippedOutline_distances_cumulative_normalized = (
            upperFlippedOutline_distances_cumulative
            / upperFlippedOutline_distances_cumulative[-1]
        )

        # Create interpolated functions for the x and y-components of points on the
        # upper outline as a function of distance along upper outline.
        upperX_func = sp_interp.PchipInterpolator(
            x=upperFlippedOutline_distances_cumulative_normalized,
            y=upperFlippedOutlineX_A_lp,
            extrapolate=False,
        )
        upperY_func = sp_interp.PchipInterpolator(
            x=upperFlippedOutline_distances_cumulative_normalized,
            y=upperFlippedOutlineY_A_lp,
            extrapolate=False,
        )

        # Get the lower outline points. This contains the leading point.
        lowerOutline_A_lp = self._lower_outline()

        lowerOutlineX_A_lp = lowerOutline_A_lp[:, 0]
        lowerOutlineY_A_lp = lowerOutline_A_lp[:, 1]

        # Find the distance between points along the lower original outline.
        lowerOutline_distances_between_points = np.sqrt(
            np.power(
                lowerOutlineX_A_lp[:-1] - lowerOutlineX_A_lp[1:],
                2,
            )
            + np.power(
                lowerOutlineY_A_lp[:-1] - lowerOutlineY_A_lp[1:],
                2,
            )
        )

        # Create a horizontal 1D array that contains the cumulative distance along
        # the lower original outline of each point.
        lowerOutline_distances_cumulative = np.hstack(
            (0, np.cumsum(lowerOutline_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        lowerOutline_distances_cumulative_normalized = (
            lowerOutline_distances_cumulative / lowerOutline_distances_cumulative[-1]
        )

        # Create interpolated functions for the x and y-components of points on the
        # lower outline as a function of distance along the lower outline.
        lowerX_func = sp_interp.PchipInterpolator(
            x=lowerOutline_distances_cumulative_normalized,
            y=lowerOutlineX_A_lp,
            extrapolate=False,
        )
        lowerY_func = sp_interp.PchipInterpolator(
            x=lowerOutline_distances_cumulative_normalized,
            y=lowerOutlineY_A_lp,
            extrapolate=False,
        )

        # Generate a cosine-spaced list of normalized distances from 0.0 to 1.0.
        cosine_spaced_normalized_distances = functions.cosspace(
            0.0, 1.0, n_points_per_side
        )

        # Find the x and y-components of the upper and lower outline points at each
        # of the resampled cosine-spaced normalized distances.
        upperResampledOutlineX_A_lp = np.flip(
            upperX_func(cosine_spaced_normalized_distances)
        )
        lowerResampledOutlineX_A_lp = lowerX_func(cosine_spaced_normalized_distances)[
            1:
        ]
        upperResampledOutlineY_A_lp = np.flip(
            upperY_func(cosine_spaced_normalized_distances)
        )
        lowerResampledOutlineY_A_lp = lowerY_func(cosine_spaced_normalized_distances)[
            1:
        ]

        resampledOutlineX_A_lp = np.hstack(
            (upperResampledOutlineX_A_lp, lowerResampledOutlineX_A_lp)
        )
        resampledOutlineY_A_lp = np.hstack(
            (upperResampledOutlineY_A_lp, lowerResampledOutlineY_A_lp)
        )

        self.outline_A_lp = np.column_stack(
            (resampledOutlineX_A_lp, resampledOutlineY_A_lp)
        )

    def _upper_outline(self):
        """This method returns a 2D array of points on the upper portion of the
        Airfoil's outline (in airfoil axes, relative to the leading point).

        The order of the returned points is from trailing edge to leading edge.
        Included is the leading point, so be careful about duplicates if
        using this method in conjunction with _lower_outline.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: (N,2) ndarray of floats

            This is an (N,2) ndarray of floats that describe the position of N points
            on the Airfoil's upper outline (in airfoil axes, relative to the leading
            point).
        """
        return self.outline_A_lp[: self._leading_edge_index() + 1, :]
