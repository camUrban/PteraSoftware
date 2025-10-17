"""This module contains the Airfoil class.

This module contains the following classes:
    Airfoil: This is a class used to contain the Airfoil of a WingCrossSection.

This module contains the following functions:
    None
"""

import importlib.resources

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sp_interp

from .. import _functions
from .. import _parameter_validation
from .. import _transformations

# Create a token object for bypassing outline_A_lp parameter validation in Airfoil's
# __init__ method.
_TRUST = object()


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

    # FIXME: Explicitly disallow NACA0000 Airfoils. They occasionally generate
    #  correctly but often causes issues due to their infinitesimal thickness.
    #  Mention the alternative of using any NACA00XX Airfoil, as these are the same
    #  thing for vortex lattice methods.

    def __init__(
        self,
        name="NACA0012",
        outline_A_lp=None,
        resample=True,
        n_points_per_side=400,
        _trust=None,
    ):
        """This is the initialization method.

        :param name: str, optional

            This is the name of the Airfoil. It should correspond to the name of a
            file the airfoils directory (once converted to lower-case and stripped of
            leading and trailing whitespace) unless you are passing in your own
            array of points using outline_A_lp. The default is "NACA0012".

        :param outline_A_lp: array-like of shape (N,2), optional

            This is an array of the 2D points making up the Airfoil's outline (in
            airfoil axes, relative to the leading point). If you wish to load
            coordinates from the airfoils directory, leave this as None, which is the
            default. If not, it must be an array-like object of numbers (int or
            float) with shape (N,2). It can be a tuple, list, or numpy array. Values
            are converted to floats internally.it must be a Nx2 numpy array of
            numbers (int or float). Make sure all x-component values are in the range
            [ 0.0, 1.0]. The default value is None.

        :param resample: boolLike, optional

            This is the variable that determines whether you would like to resample
            the points defining the Airfoil's outline. This applies to points passed
            in by the user or to those from the airfoils directory. I highly
            recommended setting this to True. It can be a boolean or a NumPy boolean
            and will be converted internally to a boolean. The default is True.

        :param n_points_per_side: int or None, optional

            This is number of points to use when creating the Airfoil's MCL and when
            resampling the upper and lower parts of the Airfoil's outline. It must be
            a positive int greater than or equal to 3. The resampled outline will
            have a total number of points equal to (2 * n_points_per_side) - 1. I
            highly recommend setting this to at least 100. The default value is 400.
        """
        self.name = _parameter_validation.string_return_string(name, "name")

        if outline_A_lp is not None:
            if _trust is not _TRUST:
                self.outline_A_lp = self._validate_outline(outline_A_lp)
            else:
                self.outline_A_lp = outline_A_lp
        else:
            self._populate_outline()

        self.resample = _parameter_validation.boolLike_return_bool(resample, "resample")

        self.n_points_per_side = _parameter_validation.int_in_range_return_int(
            n_points_per_side, "n_points_per_side", 3, True, None, None
        )

        # If resample is True, resample the Airfoil's outline points.
        if self.resample:
            self._resample_outline(self.n_points_per_side)

        # Initialize an attribute for an array of points along the MCL (in airfoil
        # axes, relative to the leading point). It will be set by _populate_mcl.
        self.mcl_A_lp = None
        self._populate_mcl()

    # TODO: In the future, if adding control surfaces becomes more important,
    #  we may want to rework this method. Using this method we need to artificially
    #  limit the maximum deflection to 5.0 degrees because higher values may cause
    #  the upper and lower outlines to intersect. This is because they each rotate
    #  about points on their respective outlines. Instead, it would be better to have
    #  everything rotate about the MCL's hinge point, however, this causes
    #  self-intersections for the upper and lower outlines, so we'd need that we'd
    #  need to write some logic to remove those.
    def add_control_surface(self, deflection, hinge_point):
        """This method returns a version of the Airfoil with a control surface added
        at a given point. It is called during meshing.

        :param deflection: number
            This is the control deflection in degrees. Deflection downwards is
            positive. It must be a number (int or float) in the range [-5.0,
            5.0] degrees. Values are converted to floats internally.
        :param hinge_point: float
            This is the location of the hinge as a fraction of chord length. It must
            be in the range (0.0, 1.0).
        :return flapped_airfoil: Airfoil
            This is the new airfoil with the control surface added.
        """
        # Validate the deflection and hinge_point inputs.
        deflection = _parameter_validation.number_in_range_return_float(
            deflection, "deflection", -5.0, True, 5.0, True
        )
        hinge_point = _parameter_validation.number_in_range_return_float(
            hinge_point, "hinge_point", 0.0, False, 1.0, False
        )

        flippedUpperOutline_A_lp = np.flipud(self._upper_outline())
        lowerOutline_A_lp = self._lower_outline()

        flippedUpperOutline_split_index = np.where(
            flippedUpperOutline_A_lp[:, 0] >= hinge_point
        )[0][0]
        lowerOutline_split_index = np.where(lowerOutline_A_lp[:, 0] >= hinge_point)[0][
            0
        ]

        preHingeFlippedUpperOutline_A_lp = flippedUpperOutline_A_lp[
            :flippedUpperOutline_split_index, :
        ]
        postHingeFlippedUpperOutline_A_lp = flippedUpperOutline_A_lp[
            flippedUpperOutline_split_index:, :
        ]
        preHingeLowerOutline_A_lp = lowerOutline_A_lp[:lowerOutline_split_index, :]
        postHingeLowerOutline_A_lp = lowerOutline_A_lp[lowerOutline_split_index:, :]

        flippedUpperOutlineHingePoint_A_lp = preHingeFlippedUpperOutline_A_lp[-1, :]
        lowerOutlineHingePoint_A_lp = preHingeLowerOutline_A_lp[-1, :]

        flippedUpperOutlineHingePoint_Wcs_lp = np.hstack(
            [
                flippedUpperOutlineHingePoint_A_lp[0],
                0.0,
                flippedUpperOutlineHingePoint_A_lp[1],
            ]
        )
        lowerOutlineHingePoint_Wcs_lp = np.hstack(
            [lowerOutlineHingePoint_A_lp[0], 0.0, lowerOutlineHingePoint_A_lp[1]]
        )

        flippedUpperOutlineToOrigin_T_act = _transformations.generate_trans_T(
            -flippedUpperOutlineHingePoint_Wcs_lp, passive=False
        )
        lowerOutlineToOrigin_T_act = _transformations.generate_trans_T(
            -lowerOutlineHingePoint_Wcs_lp, passive=False
        )

        # Make the active rotational homogeneous transformation matrix for the given
        # angle.
        rot_T = _transformations.generate_rot_T(
            (0, 0, -deflection), passive=False, intrinsic=False, order="zyx"
        )

        flippedUpperOutlineBack_T_act = _transformations.generate_trans_T(
            flippedUpperOutlineHingePoint_Wcs_lp, passive=False
        )
        lowerOutlineBack_T_act = _transformations.generate_trans_T(
            lowerOutlineHingePoint_Wcs_lp, passive=False
        )

        postHingeFlippedUpperOutline_T_act = _transformations.compose_T_act(
            flippedUpperOutlineToOrigin_T_act, rot_T, flippedUpperOutlineBack_T_act
        )
        postHingeLowerOutline_T_act = _transformations.compose_T_act(
            lowerOutlineToOrigin_T_act, rot_T, lowerOutlineBack_T_act
        )

        postHingeFlippedUpperOutline_Wcs_lp = np.column_stack(
            [
                postHingeFlippedUpperOutline_A_lp[:, 0],
                np.zeros_like(postHingeFlippedUpperOutline_A_lp[:, 0]),
                postHingeFlippedUpperOutline_A_lp[:, 1],
            ]
        )
        postHingeLowerOutline_Wcs_lp = np.column_stack(
            [
                postHingeLowerOutline_A_lp[:, 0],
                np.zeros_like(postHingeLowerOutline_A_lp[:, 0]),
                postHingeLowerOutline_A_lp[:, 1],
            ]
        )

        flappedPostHingeFlippedUpperOutline_A_lp = (
            _transformations.apply_T_to_vectors(
                postHingeFlippedUpperOutline_T_act,
                postHingeFlippedUpperOutline_Wcs_lp,
                has_point=True,
            )
        )[:, [0, 2]]
        flappedPostHingeLowerOutline_A_lp = _transformations.apply_T_to_vectors(
            postHingeLowerOutline_T_act,
            postHingeLowerOutline_Wcs_lp,
            has_point=True,
        )[:, [0, 2]]

        flappedFlippedUpperOutline_A_lp = np.vstack(
            [preHingeFlippedUpperOutline_A_lp, flappedPostHingeFlippedUpperOutline_A_lp]
        )
        flappedLowerOutline_A_lp = np.vstack(
            [preHingeLowerOutline_A_lp, flappedPostHingeLowerOutline_A_lp]
        )

        flappedOutline_A_lp = np.vstack(
            [
                np.flipud(flappedFlippedUpperOutline_A_lp),
                flappedLowerOutline_A_lp[1:, :],
            ]
        )

        # Return the new flapped Airfoil, with the _TRUST token so that we don't
        # re-validate the outline, which would fail because the validation requires
        # the trailing edge points be roughly at y=0.0 (in Airfoil axes, relative to
        # the leading point).
        return Airfoil(
            name=self.name + " flapped",
            outline_A_lp=flappedOutline_A_lp,
            resample=False,
            n_points_per_side=self.n_points_per_side,
            _trust=_TRUST,
        )

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

    # TEST: Consider adding unit tests for this method.
    # DOCUMENT: After testing it, document this method.
    def get_plottable_data(self, show=False):
        """

        :return:
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        outline_A_lp = self.outline_A_lp
        mcl_A_lp = self.mcl_A_lp

        if not show:
            return [outline_A_lp, mcl_A_lp]

        airfoil_figure, airfoil_axes = plt.subplots()

        outlineX_A_lp = self.outline_A_lp[:, 0]
        outlineY_A_lp = self.outline_A_lp[:, 1]
        mclX_A_lp = self.mcl_A_lp[:, 0]
        mclY_A_lp = self.mcl_A_lp[:, 1]

        outlineYMin_A_lp = np.min(outlineY_A_lp)
        outlineYMax_A_lp = np.max(outlineY_A_lp)
        outlineYRange_A_lp = outlineYMax_A_lp - outlineYMin_A_lp
        y_padding = 0.1 * outlineYRange_A_lp

        airfoil_axes.plot(outlineX_A_lp, outlineY_A_lp, "b-")
        airfoil_axes.plot(mclX_A_lp, mclY_A_lp, "r-")

        airfoil_axes.set_xlim(0.0, 1.0)
        airfoil_axes.set_ylim(
            outlineYMin_A_lp - y_padding, outlineYMax_A_lp + y_padding
        )

        airfoil_axes.set_xlabel("AX_lp")
        airfoil_axes.set_ylabel("AY_lp")
        airfoil_axes.set_title(f"{self.name} Airfoil")
        airfoil_axes.legend(
            ["Outline", "Mean Camber Line (MCL)"],
            loc="lower center",
            bbox_to_anchor=(0.5, -1.75),
        )

        airfoil_axes.set_aspect("equal", adjustable="box")

        airfoil_figure.show()

        return None

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
        mcl_fractions = _parameter_validation.nD_number_vectorLike_return_float(
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

    def _lp_index(self):
        """Returns the index of the leading point in the outline_A_lp attribute.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: int
            This is the index of the leading point.
        """
        return np.argmin(self.outline_A_lp[:, 0])

    def _lower_outline(self):
        """This method returns a 2D array of points on the lower portion of the
        Airfoil's outline (in airfoil axes, relative to the leading point).

        The order of the returned points is from leading point to trailing edge.
        Included is the leading point, so be careful about duplicates if
        using this method in conjunction with _upper_outline.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: (N,2) ndarray of floats

            This is an (N,2) ndarray of floats that describe the position of N points
            on the Airfoil's lower outline (in airfoil axes, relative to the leading
            point).
        """
        return self.outline_A_lp[self._lp_index() :, :]

    def _populate_mcl(self):
        """This method creates a 2D array of points along the Airfoil's MCL (in
        airfoil axes, relative to the leading point), which it uses to set the mcl_A_lp
        attribute. It is in order from the leading point to the trailing point.

        :return: None
        """
        # Split outline_A_lp into upper and lower sections. Flip the upper points so
        # that they are ordered from the leading point to the trailing point.
        flippedUpperOutline_A_lp = np.flipud(self._upper_outline())
        lowerOutline_A_lp = self._lower_outline()

        cosine_spaced_chord_fractions = _functions.cosspace(
            0.0, 1.0, self.n_points_per_side
        )

        upper_func = sp_interp.PchipInterpolator(
            x=flippedUpperOutline_A_lp[:, 0],
            y=flippedUpperOutline_A_lp[:, 1],
            extrapolate=True,
        )
        lower_func = sp_interp.PchipInterpolator(
            x=lowerOutline_A_lp[:, 0],
            y=lowerOutline_A_lp[:, 1],
            extrapolate=True,
        )

        flippedUpperOutlineY_A_lp = upper_func(cosine_spaced_chord_fractions)
        lowerOutlineY_A_lp = lower_func(cosine_spaced_chord_fractions)

        # Calculate the approximate MCL points (in airfoil axes, relative to the
        # leading point) and set the class attribute.
        self.mcl_A_lp = np.column_stack(
            [
                cosine_spaced_chord_fractions,
                (flippedUpperOutlineY_A_lp + lowerOutlineY_A_lp) / 2,
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
                    mclX_A_lp = _functions.cosspace(0, 1, n_points_per_side)

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

            # Get the path to the _airfoils data directory.
            airfoil_file = (
                importlib.resources.files("pterasoftware.geometry")
                .joinpath("_airfoils")
                .joinpath(sanitized_name + ".dat")
            )

            # Read the text from the airfoil file.
            raw_text = airfoil_file.read_text()

            # Trim the text at the return characters.
            trimmed_text = raw_text[raw_text.find("\n") :]

            # Input the coordinates into a 1D array. This represents the upper and
            # lower points of the Airfoil's outline (in airfoil axes, relative to the
            # leading point).
            outline1D_A_lp = np.fromstring(trimmed_text, sep="\n")

            # Check to make sure the number of elements in the array is even.
            if len(outline1D_A_lp) % 2 != 0:
                raise ValueError(
                    "name matched to an airfoil in the airfoils database, but it could not be read correctly."
                )

            # Populate the outline_A_lp attribute and return.
            self.outline_A_lp = np.reshape(outline1D_A_lp, (-1, 2))
            return

        # If the Airfoil was not a NACA 4-series and was not found in the database,
        # throw an error.
        except FileNotFoundError:
            raise ValueError(
                "name didn't match the NACA 4-series pattern nor was it found in the airfoils database."
            )

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

        flippedUpperOutlineX_A_lp = np.flipud(upperOutline_A_lp)[:, 0]
        flippedUpperOutlineY_A_lp = np.flipud(upperOutline_A_lp)[:, 1]

        # Find the distance between points along the upper flipped original outline.
        flippedUpperOutline_distances_between_points = np.sqrt(
            np.power(
                flippedUpperOutlineX_A_lp[:-1] - flippedUpperOutlineX_A_lp[1:],
                2,
            )
            + np.power(
                flippedUpperOutlineY_A_lp[:-1] - flippedUpperOutlineY_A_lp[1:],
                2,
            )
        )

        # Create a horizontal 1D array that contains the cumulative distance along
        # the upper flipped original outline of each point.
        flippedUpperOutline_distances_cumulative = np.hstack(
            (0, np.cumsum(flippedUpperOutline_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        flippedUpperOutline_distances_cumulative_normalized = (
            flippedUpperOutline_distances_cumulative
            / flippedUpperOutline_distances_cumulative[-1]
        )

        # Create interpolated functions for the x and y-components of points on the
        # upper outline as a function of distance along upper outline.
        upperX_func = sp_interp.PchipInterpolator(
            x=flippedUpperOutline_distances_cumulative_normalized,
            y=flippedUpperOutlineX_A_lp,
            extrapolate=False,
        )
        upperY_func = sp_interp.PchipInterpolator(
            x=flippedUpperOutline_distances_cumulative_normalized,
            y=flippedUpperOutlineY_A_lp,
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
        cosine_spaced_normalized_distances = _functions.cosspace(
            0.0, 1.0, n_points_per_side
        )

        # Find the x and y-components of the upper and lower outline points at each
        # of the resampled cosine-spaced normalized distances.
        upperResampledOutlineX_A_lp = np.flipud(
            upperX_func(cosine_spaced_normalized_distances)
        )
        lowerResampledOutlineX_A_lp = lowerX_func(cosine_spaced_normalized_distances)[
            1:
        ]
        upperResampledOutlineY_A_lp = np.flipud(
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

        The order of the returned points is from trailing edge to leading point.
        Included is the leading point, so be careful about duplicates if
        using this method in conjunction with _lower_outline.

        Note: This is a private method, so it doesn't perform any parameter validation.

        :return: (N,2) ndarray of floats

            This is an (N,2) ndarray of floats that describe the position of N points
            on the Airfoil's upper outline (in airfoil axes, relative to the leading
            point).
        """
        return self.outline_A_lp[: self._lp_index() + 1, :]

    @staticmethod
    def _validate_outline(outline_A_lp):
        """This method performs validates a user's provided outline_A_lp. However,
        it will fail for "flapped" airfoils.

        :param outline_A_lp: Any
            This is the input to validate.
        :return: (N,2) ndarray of floats
            This is the validated version of outline_A_lp
        """
        outline_A_lp = (
            _parameter_validation.arrayLike_of_twoD_number_vectorLikes_return_float(
                outline_A_lp, "outline_A_lp"
            )
        )

        n_outline_points = outline_A_lp.shape[0]

        # The outline must have at least 5 points.
        if n_outline_points < 5:
            raise ValueError("The Airfoil's outline must have at least five points")

        # The Airfoil must have roughly a unit chord length.
        allowance = 0.02
        outlineMaxX_A_lp = max(outline_A_lp[:, 0])
        outlineMinX_A_lp = min(outline_A_lp[:, 0])
        outlineChord = outlineMaxX_A_lp - outlineMinX_A_lp
        if outlineChord > 1 + allowance or outlineChord < 1 - allowance:
            raise ValueError(
                "The Airfoil's outline must have a chord length of roughly 1.0 m"
            )

        # The Airfoil's outline must have roughly a thickness of at least 0.1%.
        outlineMaxY_A_lp = max(outline_A_lp[:, 1])
        outlineMinY_A_lp = min(outline_A_lp[:, 1])
        outlineMaxThickness = outlineMaxY_A_lp - outlineMinY_A_lp
        if outlineMaxThickness < 0.001:
            raise ValueError(
                "The Airfoil's outline must have max thickness of at least 0.1%"
            )

        xAllowance = allowance * outlineChord
        yAllowance = allowance * outlineMaxThickness

        # Check that the upper outline's trailing point is at approximately [1.0,
        # 0.0] (in Airfoil axes, relative to the leading point).
        outlineUpperTp_A_lp = outline_A_lp[0, :]
        outlineUpperTpX_A_lp = outlineUpperTp_A_lp[0]
        outlineUpperTpY_A_lp = outlineUpperTp_A_lp[1]
        if (
            outlineUpperTpX_A_lp < 1.0 - xAllowance
            or outlineUpperTpX_A_lp > 1.0 + xAllowance
        ):
            raise ValueError(
                "The x-value of the Airfoil's upper outline's trailing point must be approximately 1.0 m (in Airfoil axes, relative to the leading point)"
            )
        if outlineUpperTpY_A_lp < -yAllowance or outlineUpperTpY_A_lp > yAllowance:
            raise ValueError(
                "The y-value of the Airfoil's upper outline's trailing point must be approximately 0.0 m (in Airfoil axes, relative to the leading point)"
            )

        # Check that the lower outline's trailing point is at approximately [1.0,
        # 0.0] (in Airfoil axes, relative to the leading point).
        outlineLowerTp_A_lp = outline_A_lp[-1, :]
        outlineLowerTpX_A_lp = outlineLowerTp_A_lp[0]
        outlineLowerTpY_A_lp = outlineLowerTp_A_lp[1]
        if (
            outlineLowerTpX_A_lp < 1.0 - xAllowance
            or outlineLowerTpX_A_lp > 1.0 + xAllowance
        ):
            raise ValueError(
                "The x-value of the Airfoil's lower outline's trailing point must be approximately 1.0 (in Airfoil axes, relative to the leading point)"
            )
        if outlineLowerTpY_A_lp < -yAllowance or outlineLowerTpY_A_lp > yAllowance:
            raise ValueError(
                "The y-value of the Airfoil's lower outline's trailing point must be approximately 0.0 (in Airfoil axes, relative to the leading point)"
            )

        # Check that the upper outline's trailing point has a y-value that's greater
        # than or equal to the lower outline's trailing point's y-value.
        if outlineLowerTpY_A_lp > outlineUpperTpY_A_lp:
            raise ValueError(
                "The upper outline's trailing point must have a y-value that's greater than or equal to the lower outline's trailing point's y-value (in Airfoil axes)"
            )

        # TODO: Consider moving this to another function for "normalizing" an
        #  Airfoil's outline such that its leading point is exactly at [0.0, 0.0],
        #  and its trailing point exactly at [1.0, 0.0]. However, we'd have to be
        #  careful we don't reject "flapped" Airfoils.
        # # If the upper outline's trailing point isn't exactly on the x-axis,
        # # check that it's slope is in the direction of the x-axis. If so, find the
        # # linearly extrapolated intersection with the x-axis (in Airfoil axes,
        # # relative to the leading point).
        # outlineUpperExtrapTpX_A_lp = outlineUpperTpX_A_lp
        # if outlineUpperTpY_A_lp != 0.0:
        #     outlineUpperPtp = outline_A_lp[1, :]
        #     outlineUpperSlopeTp_A = (outlineUpperTpY_A_lp - outlineUpperPtp[1]) / (
        #         outlineUpperTpX_A_lp - outlineUpperPtp[0]
        #     )
        #     if outlineUpperTpY_A_lp > 0.0:
        #         if outlineUpperSlopeTp_A >= 0.0:
        #             raise ValueError(
        #                 "If the y-value of the Airfoil's upper outline's trailing point is positive, the slope at the trailing point of the upper outline must be negative"
        #             )
        #     elif outlineUpperTpY_A_lp < 0.0:
        #         if outlineUpperSlopeTp_A <= 0.0:
        #             raise ValueError(
        #                 "If the y-value of the Airfoil's upper outline's trailing point is negative, the slope at the trailing point of the upper outline must be positive"
        #             )
        #     outlineUpperExtrapTpX_A_lp -= outlineUpperTpY_A_lp / outlineUpperSlopeTp_A
        #
        # # If the lower outline's trailing point isn't exactly on the x-axis,
        # # check that it's slope is in the direction of the x-axis. If so, find the
        # # linearly extrapolated intersection with the x-axis (in Airfoil axes,
        # # relative to the leading point).
        # outlineLowerExtrapTpX_A_lp = outlineLowerTpX_A_lp
        # if outlineLowerTpY_A_lp != 0.0:
        #     outlineLowerPtp = outline_A_lp[-2, :]
        #     outlineLowerSlopeTp_A = (outlineLowerTpY_A_lp - outlineLowerPtp[1]) / (
        #         outlineLowerTpX_A_lp - outlineLowerPtp[0]
        #     )
        #     if outlineLowerTpY_A_lp > 0.0:
        #         if outlineLowerSlopeTp_A >= 0.0:
        #             raise ValueError(
        #                 "If the y-value of the Airfoil's lower outline's trailing point is positive, the slope at the trailing point of the lower outline must be negative"
        #             )
        #     elif outlineLowerTpY_A_lp < 0.0:
        #         if outlineLowerSlopeTp_A <= 0.0:
        #             raise ValueError(
        #                 "If the y-value of the Airfoil's lower outline's trailing point is negative, the slope at the trailing point of the lower outline must be positive"
        #             )
        #     outlineLowerExtrapTpX_A_lp -= outlineLowerTpY_A_lp / outlineLowerSlopeTp_A
        #
        # # Find the average extrapolated trailing point.
        # outlineExtrapTpX = (outlineUpperExtrapTpX_A_lp + outlineLowerExtrapTpX_A_lp) / 2
        # outlineExtrapTp = np.array([outlineExtrapTpX, 0.0], dtype=float)
        #
        # # Add the extrapolated trailing point to outline_A_lp while avoiding duplicates.
        # if outlineUpperTpY_A_lp != 0.0:
        #     if outlineLowerTpY_A_lp != 0.0:
        #         outline_A_lp = np.hstack([outlineExtrapTp, outline_A_lp])
        #     else:
        #         outline_A_lp = np.hstack([outlineExtrapTp, outline_A_lp[:-1, :]])
        # else:
        #     if outlineLowerTpY_A_lp != 0.0:
        #         outline_A_lp = np.hstack([outlineExtrapTp, outline_A_lp[1:, :]])
        #     else:
        #         outline_A_lp = np.hstack([outlineExtrapTp, outline_A_lp[1:-1, :]])

        # Find the index of the outline's leading point.
        outlineLp_index = np.argmin(outline_A_lp[:, 0])

        # Check the outline's leading point is at approximately [0.0, 0.0] (in
        # Airfoil axes, relative to the leading point).
        outlineLp_A_lp = outline_A_lp[outlineLp_index, :]
        outlineLpX_A_lp = outlineLp_A_lp[0]
        outlineLpY_A_lp = outlineLp_A_lp[1]
        if outlineLpX_A_lp < -xAllowance or outlineLpX_A_lp > xAllowance:
            raise ValueError(
                "The x-value of the Airfoil's outline's leading point must be approximately 0.0 (in Airfoil axes, relative to the leading point)"
            )
        if outlineLpY_A_lp < -yAllowance or outlineLpY_A_lp > yAllowance:
            raise ValueError(
                "The y-value of the Airfoil's outline's leading point must be approximately 0.0 (in Airfoil axes, relative to the leading point)"
            )

        # Split the outline into its upper and lower sections.
        lowerOutline_A_lp = outline_A_lp[outlineLp_index:, :]
        upperOutline_A_lp = outline_A_lp[: outlineLp_index + 1, :]

        lowerOutlineDiffX_A = np.diff(lowerOutline_A_lp[:, 0])
        upperOutlineDiffX_A = np.diff(upperOutline_A_lp[:, 0])

        # Check that the upper outline is strictly decreasing in x and that the lower
        # outline is strictly increasing in x (in Airfoil axes).
        if not np.all(upperOutlineDiffX_A < 0.0):
            raise ValueError(
                "The every point in the Airfoil's outline's upper portion must have an x-value less than the point before it (in Airfoil axes)"
            )
        if not np.all(lowerOutlineDiffX_A > 0.0):
            raise ValueError(
                "The every point in the Airfoil's outline's lower portion must have an x-value greater than than the point before it (in Airfoil axes)"
            )

        # Check that the split portions both have at least three points.
        n_upper_points = upperOutline_A_lp.shape[0]
        n_lower_points = lowerOutline_A_lp.shape[0]
        if n_upper_points < 3:
            raise ValueError(
                "The upper portion of the Airfoil's outline must contain at least three points (including the outline's leading point)"
            )
        if n_lower_points < 3:
            raise ValueError(
                "The lower portion of the Airfoil's outline must contain at least three points (including the outline's leading point)"
            )

        outline_chord_fractions = np.linspace(
            outlineLpX_A_lp,
            max(outlineUpperTpX_A_lp, outlineUpperTpY_A_lp),
            2 * max(n_upper_points, n_lower_points),
        )

        flippedUpperOutline_A_lp = np.flipud(upperOutline_A_lp)

        upperY_func = sp_interp.PchipInterpolator(
            x=flippedUpperOutline_A_lp[:, 0],
            y=flippedUpperOutline_A_lp[:, 1],
            extrapolate=True,
        )
        lowerY_func = sp_interp.PchipInterpolator(
            x=lowerOutline_A_lp[:, 0],
            y=lowerOutline_A_lp[:, 1],
            extrapolate=True,
        )

        upperMinusLowerOutlineY = upperY_func(outline_chord_fractions) - lowerY_func(
            outline_chord_fractions
        )

        if not np.all(upperMinusLowerOutlineY[1:-1] >= 0.0):
            raise ValueError(
                "All points on the Airfoil's upper outline (excluding the outline's leading and trailing points) must have a y-value greater than the points on the lower outline at the same x-value (in Airfoil axes, relative to the leading point)"
            )

        return outline_A_lp
