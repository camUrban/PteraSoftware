"""Contains the Airfoil class.

**Contains the following classes:**

Airfoil: A class used to contain the Airfoil of a WingCrossSection.

**Contains the following functions:**

None
"""

from __future__ import annotations

import importlib.resources
import warnings
from collections.abc import Sequence
from typing import Any, cast

import matplotlib.pyplot as plt
import numpy as np

from .. import _functions, _parameter_validation, _transformations

# Create a token object for bypassing outline_A_Lp parameter validation in Airfoil's
# __init__ method.
_TRUST = object()

# Create a sentinel for detecting use of the deprecated outline_A_lp parameter. It is
# annotated as Any so the deprecated parameter can carry it as a default while keeping
# an outline-shaped type hint at the call sites mypy sees.
_UNSET: Any = object()


class Airfoil:
    """A class used to contain the Airfoil of a WingCrossSection.

    **Contains the following methods:**

    __deepcopy__: Returns an independent deep copy of this Airfoil.

    add_control_surface: Returns a version of the Airfoil with a control surface added
    at a given point.

    draw: Plots this Airfoil's outlines and mean camber line (MCL) using PyPlot.

    get_plottable_data: Returns plottable data for this Airfoil's outline and mean
    camber line.

    get_resampled_mcl: Returns a ndarray of points along the mean camber line (MCL),
    resampled from the mcl_A_Lp attribute. It is used to discretize the MCL for meshing.

    **Citation:**

    Adapted from: geometry.Airfoil in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/27/2020
    """

    __slots__ = (
        "_name",
        "_outline_A_Lp",
        "_resample",
        "_n_points_per_side",
        "_mcl_A_Lp",
    )

    def __init__(
        self,
        name: str = "NACA0012",
        outline_A_Lp: np.ndarray | Sequence[Sequence[float | int]] | None = None,
        resample: bool | np.bool_ = True,
        n_points_per_side: int = 400,
        _trust: object | None = None,
        outline_A_lp: Any = _UNSET,
    ) -> None:
        """The initialization method.

        :param name: The name of the Airfoil. It should correspond to the name of a file
            the airfoils directory, or to a valid NACA 4 series airfoil (once converted
            to lower-case and stripped of leading and trailing whitespace) unless you
            are passing in your own array of points using outline_A_Lp. Note that
            NACA0000 isn't a valid NACA 4 series airfoil, NACA 4 series airfoils with
            thickness above 30% are not supported, the first two digits must either both
            be zero (symmetric) or both be non zero (cambered), and for cambered
            airfoils the position of maximum camber must be greater than or equal to the
            maximum camber plus half the maximum thickness. The default is "NACA0012".
        :param outline_A_Lp: An array like object of numbers (int or float) with shape
            (N,2) representing the 2D points making up the Airfoil's outline. If you
            wish to load coordinates from the airfoils directory, leave this as None,
            which is the default. Can be a tuple, list, or ndarray. Values are converted
            to floats internally. The outline is automatically normalized to canonical
            form (leading point at origin, chord on x axis, unit chord length), so
            position and scale do not need to match this. Minor rotation offsets (less
            than 15 degrees, such as implicit angle of attack in airfoil data) are also
            corrected; larger rotations will raise a ValueError. The topology must be
            correct: points must be ordered from upper trailing point to leading point
            to lower trailing point. Also, after correcting for rotation, the upper
            portion's x values must be non increasing, the lower portion's x values must
            be non decreasing. Finally, both portions must have at least 3 unique
            points, and the outline must not self intersect (upper y must be strictly
            greater than lower y at all interior x positions). Open and blunt trailing
            edges are supported. The default value is None.
        :param resample: Determines whether to resample the points defining the
            Airfoil's outline. This applies to points passed in by the user or to those
            from the airfoils directory. I highly recommended setting this to True. Can
            be a bool or a numpy bool and will be converted internally to a bool. The
            default is True.
        :param n_points_per_side: The number of points to use when creating the
            Airfoil's MCL and when resampling the upper and lower parts of the Airfoil's
            outline. It must be a positive int greater than or equal to 3. The resampled
            outline will have a total number of points equal to (2 * n_points_per_side)
            - 1. I highly recommend setting this to at least 100. The default value is
            400.
        :param outline_A_lp: A deprecated alias for outline_A_Lp. Passing it emits a
            DeprecationWarning, and it will be removed in v6.0.0. It must not be passed
            together with outline_A_Lp. The default is a private sentinel meaning the
            parameter is unset.
        :return: None
        """
        if outline_A_lp is not _UNSET:
            warnings.warn(
                "The outline_A_lp parameter is deprecated and will be removed in "
                "v6.0.0. Use outline_A_Lp instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if outline_A_Lp is not None:
                raise ValueError(
                    "outline_A_Lp and the deprecated outline_A_lp cannot both be "
                    "passed."
                )
            outline_A_Lp = outline_A_lp

        self._name = _parameter_validation.str_return_str(name, "name")

        if outline_A_Lp is not None:
            if _trust is not _TRUST:
                # Validate, normalize, and final validate user provided outlines.
                self._outline_A_Lp = self._validate_outline_preliminary(outline_A_Lp)
                self._normalize_outline()
                self._validate_outline_final()
            else:
                # When _trust is _TRUST, we know outline_A_Lp is already validated.
                self._outline_A_Lp = cast(np.ndarray, outline_A_Lp)
        else:
            self._populate_outline()
            # Validate, normalize, and final validate database and generated NACA
            # outlines.
            self._outline_A_Lp = self._validate_outline_preliminary(self._outline_A_Lp)
            self._normalize_outline()
            self._validate_outline_final()

        self._resample = _parameter_validation.boolLike_return_bool(
            resample, "resample"
        )

        self._n_points_per_side = _parameter_validation.int_in_range_return_int(
            n_points_per_side, "n_points_per_side", 3, True, None, None
        )

        # If resample is True, resample the Airfoil's outline points.
        if self._resample:
            self._resample_outline(self._n_points_per_side)

        # Initialize an attribute for an array of points along the MCL (in airfoil axes,
        # relative to the leading point). It will be set by _populate_mcl.
        self._mcl_A_Lp: np.ndarray | None = None
        self._populate_mcl()

        # Set immutable numpy arrays to be read only.
        self._outline_A_Lp.flags.writeable = False
        if self._mcl_A_Lp is not None:
            self._mcl_A_Lp.flags.writeable = False

    # --- Deep copy method ---
    def __deepcopy__(self, memo: dict[int, Any]) -> Airfoil:
        """Creates a deep copy of this Airfoil.

        The copy preserves all attributes since they are all immutable.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new Airfoil with copied data.
        """
        # Create a new Airfoil instance without calling __init__.
        new_airfoil = object.__new__(Airfoil)

        # Store this Airfoil in memo to handle potential circular references.
        memo[id(self)] = new_airfoil

        # Copy immutable attributes. For those that are numpy arrays, make the copies
        # read only.
        new_airfoil._name = self._name
        new_airfoil._resample = self._resample
        new_airfoil._n_points_per_side = self._n_points_per_side

        new_airfoil._outline_A_Lp = self._outline_A_Lp.copy()
        new_airfoil._outline_A_Lp.flags.writeable = False

        if self._mcl_A_Lp is not None:
            new_airfoil._mcl_A_Lp = self._mcl_A_Lp.copy()
            new_airfoil._mcl_A_Lp.flags.writeable = False
        else:
            new_airfoil._mcl_A_Lp = None

        return new_airfoil

    # --- Immutable: read only properties ---
    @property
    def name(self) -> str:
        return self._name

    @property
    def outline_A_Lp(self) -> np.ndarray:
        return self._outline_A_Lp

    @property
    def resample(self) -> bool:
        return self._resample

    @property
    def n_points_per_side(self) -> int:
        return self._n_points_per_side

    @property
    def mcl_A_Lp(self) -> np.ndarray | None:
        return self._mcl_A_Lp

    # --- Deprecated aliases: remove in v6.0.0 ---
    @property
    def outline_A_lp(self) -> np.ndarray:
        """A deprecated alias for outline_A_Lp.

        Reading it emits a DeprecationWarning, and it will be removed in v6.0.0.

        :return: A (N,2) ndarray of floats representing the positions of the points
            making up the Airfoil's outline (in airfoil axes, relative to the leading
            point).
        """
        warnings.warn(
            "The outline_A_lp property is deprecated and will be removed in v6.0.0. "
            "Use outline_A_Lp instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._outline_A_Lp

    @property
    def mcl_A_lp(self) -> np.ndarray | None:
        """A deprecated alias for mcl_A_Lp.

        Reading it emits a DeprecationWarning, and it will be removed in v6.0.0.

        :return: A (N,2) ndarray of floats representing the positions of the points
            along the Airfoil's mean camber line (in airfoil axes, relative to the
            leading point), or None if the MCL has not yet been populated.
        """
        warnings.warn(
            "The mcl_A_lp property is deprecated and will be removed in v6.0.0. Use "
            "mcl_A_Lp instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._mcl_A_Lp

    # --- Other methods ---
    # TODO: In the future, if adding control surfaces becomes more important, we may
    #  want to rework this method. Using this method we need to artificially limit the
    #  maximum deflection to 5.0 degrees because higher values may cause the upper and
    #  lower outlines to intersect. This is because they each rotate about points on
    #  their respective outlines. Instead, it would be better to have everything rotate
    #  about the MCL's hinge point, however, this causes self intersections for the
    #  upper and lower outlines, so we'd need to write some logic to remove those.
    def add_control_surface(
        self, deflection: float | int, hinge_point: float | int
    ) -> Airfoil:
        """Returns a version of the Airfoil with a control surface added at a given
        point. It is called during meshing.

        :param deflection: The control deflection in degrees. Deflection downwards is
            positive. It must be a number (int or float) in the range [-5.0, 5.0]
            degrees. Values are converted to floats internally.
        :param hinge_point: The location of the hinge as a fraction of chord length. It
            must be a number (int or float) in the range (0.0, 1.0). Values are
            converted to floats internally.
        :return: The new Airfoil with the control surface added.
        """
        # Validate the deflection and hinge_point inputs.
        deflection = _parameter_validation.number_in_range_return_float(
            deflection, "deflection", -5.0, True, 5.0, True
        )

        # Early return optimization: no modification needed for zero deflection.
        if deflection == 0.0:
            return self

        hinge_point = _parameter_validation.number_in_range_return_float(
            hinge_point, "hinge_point", 0.0, False, 1.0, False
        )

        flippedUpperOutline_A_Lp = np.flipud(self._upper_outline())
        lowerOutline_A_Lp = self._lower_outline()

        flippedUpperOutline_split_index = np.where(
            flippedUpperOutline_A_Lp[:, 0] >= hinge_point
        )[0][0]
        lowerOutline_split_index = np.where(lowerOutline_A_Lp[:, 0] >= hinge_point)[0][
            0
        ]

        preHingeFlippedUpperOutline_A_Lp = flippedUpperOutline_A_Lp[
            :flippedUpperOutline_split_index, :
        ]
        postHingeFlippedUpperOutline_A_Lp = flippedUpperOutline_A_Lp[
            flippedUpperOutline_split_index:, :
        ]
        preHingeLowerOutline_A_Lp = lowerOutline_A_Lp[:lowerOutline_split_index, :]
        postHingeLowerOutline_A_Lp = lowerOutline_A_Lp[lowerOutline_split_index:, :]

        flippedUpperOutlineHingePoint_A_Lp = preHingeFlippedUpperOutline_A_Lp[-1, :]
        lowerOutlineHingePoint_A_Lp = preHingeLowerOutline_A_Lp[-1, :]

        flippedUpperOutlineHingePoint_Wcs_lp = np.hstack(
            [
                flippedUpperOutlineHingePoint_A_Lp[0],
                0.0,
                flippedUpperOutlineHingePoint_A_Lp[1],
            ]
        )
        lowerOutlineHingePoint_Wcs_lp = np.hstack(
            [lowerOutlineHingePoint_A_Lp[0], 0.0, lowerOutlineHingePoint_A_Lp[1]]
        )

        flippedUpperOutlineToOrigin_T_act = _transformations.generate_trans_T(
            -flippedUpperOutlineHingePoint_Wcs_lp, passive=False
        )
        lowerOutlineToOrigin_T_act = _transformations.generate_trans_T(
            -lowerOutlineHingePoint_Wcs_lp, passive=False
        )

        # Make the active rotational homogeneous transformation matrix that deflects the
        # aft portion about the spanwise hinge axis. The hinge points are in wing cross
        # section axes, where +x is chordwise, +y is spanwise (normal to the cross
        # section's plane), and +z is toward the wing's top surface. The hinge line
        # therefore runs along +y, so the deflection is a rotation about +y. A positive
        # deflection swings the trailing edge toward -z (downward, toward the lower
        # surface), matching the convention that downward deflection is positive.
        deflectionAngles_act_ezyx = np.array([0.0, deflection, 0.0])
        deflect_T_act = _transformations.generate_rot_T(
            deflectionAngles_act_ezyx,
            passive=False,
            intrinsic=False,
            order="zyx",
        )

        flippedUpperOutlineBack_T_act = _transformations.generate_trans_T(
            flippedUpperOutlineHingePoint_Wcs_lp, passive=False
        )
        lowerOutlineBack_T_act = _transformations.generate_trans_T(
            lowerOutlineHingePoint_Wcs_lp, passive=False
        )

        postHingeFlippedUpperOutline_T_act = _transformations.compose_T_act(
            flippedUpperOutlineToOrigin_T_act,
            deflect_T_act,
            flippedUpperOutlineBack_T_act,
        )
        postHingeLowerOutline_T_act = _transformations.compose_T_act(
            lowerOutlineToOrigin_T_act, deflect_T_act, lowerOutlineBack_T_act
        )

        postHingeFlippedUpperOutline_Wcs_lp = np.column_stack(
            [
                postHingeFlippedUpperOutline_A_Lp[:, 0],
                np.zeros_like(postHingeFlippedUpperOutline_A_Lp[:, 0]),
                postHingeFlippedUpperOutline_A_Lp[:, 1],
            ]
        )
        postHingeLowerOutline_Wcs_lp = np.column_stack(
            [
                postHingeLowerOutline_A_Lp[:, 0],
                np.zeros_like(postHingeLowerOutline_A_Lp[:, 0]),
                postHingeLowerOutline_A_Lp[:, 1],
            ]
        )

        flappedPostHingeFlippedUpperOutline_A_Lp = (
            _transformations.apply_T_to_vectors(
                postHingeFlippedUpperOutline_T_act,
                postHingeFlippedUpperOutline_Wcs_lp,
                is_position=True,
            )
        )[:, [0, 2]]
        flappedPostHingeLowerOutline_A_Lp = _transformations.apply_T_to_vectors(
            postHingeLowerOutline_T_act,
            postHingeLowerOutline_Wcs_lp,
            is_position=True,
        )[:, [0, 2]]

        flappedFlippedUpperOutline_A_Lp = np.vstack(
            [preHingeFlippedUpperOutline_A_Lp, flappedPostHingeFlippedUpperOutline_A_Lp]
        )
        flappedLowerOutline_A_Lp = np.vstack(
            [preHingeLowerOutline_A_Lp, flappedPostHingeLowerOutline_A_Lp]
        )

        flappedOutline_A_Lp = np.vstack(
            [
                np.flipud(flappedFlippedUpperOutline_A_Lp),
                flappedLowerOutline_A_Lp[1:, :],
            ]
        )

        # Return the new flapped Airfoil, with the _TRUST token so that we don't
        # re-validate the outline, which would fail because the validation requires the
        # trailing edge points be roughly at y = 0.0 (in airfoil axes, relative to the
        # leading point).
        return Airfoil(
            name=self.name + " flapped",
            outline_A_Lp=flappedOutline_A_Lp,
            resample=False,
            n_points_per_side=self.n_points_per_side,
            _trust=_TRUST,
        )

    def draw(self) -> None:
        """Plots this Airfoil's outlines and mean camber line (MCL) using PyPlot.

        :return: None
        """
        outlineX_A_Lp = self._outline_A_Lp[:, 0]
        outlineY_A_Lp = self._outline_A_Lp[:, 1]

        assert self.mcl_A_Lp is not None
        mclX_A_Lp = self.mcl_A_Lp[:, 0]
        mclY_A_Lp = self.mcl_A_Lp[:, 1]

        outlineYMin_A_Lp = float(np.min(outlineY_A_Lp))
        outlineYMax_A_Lp = float(np.max(outlineY_A_Lp))
        outlineYRange_A_Lp = outlineYMax_A_Lp - outlineYMin_A_Lp
        y_padding = 0.75 * outlineYRange_A_Lp

        outlineXMin_A_Lp = float(np.min(outlineX_A_Lp))
        outlineXMax_A_Lp = float(np.max(outlineX_A_Lp))
        outlineXRange_A_Lp = outlineXMax_A_Lp - outlineXMin_A_Lp
        x_padding = 0.05 * outlineXRange_A_Lp

        plt.plot(outlineX_A_Lp, outlineY_A_Lp, "b-")
        plt.plot(mclX_A_Lp, mclY_A_Lp, "r-")

        plt.xlim(outlineXMin_A_Lp - x_padding, outlineXMax_A_Lp + x_padding)
        plt.ylim(outlineYMin_A_Lp - y_padding, outlineYMax_A_Lp + y_padding)

        plt.xlabel("x (airfoil axes)")
        plt.ylabel("y (airfoil axes)")
        plt.title(f"Airfoil: {self.name}")
        plt.legend(["Outline", "Mean Camber Line (MCL)"])

        plt.gca().set_aspect("equal", adjustable="box")

        plt.show()

    def get_plottable_data(
        self, show: bool | np.bool_ = False
    ) -> list[np.ndarray] | None:
        """Returns plottable data for this Airfoil's outline and mean camber line.

        :param show: Determines whether to display the plot. Can be a bool or a numpy
            bool, and will be converted internally to a bool. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. The default is False.
        :return: A list of two ndarrays containing the outline and MCL data, or None if
            show is True.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        if not show:
            assert self.mcl_A_Lp is not None
            return [self.outline_A_Lp, self.mcl_A_Lp]

        airfoil_figure, airfoil_axes = plt.subplots()

        outlineX_A_Lp = self._outline_A_Lp[:, 0]
        outlineY_A_Lp = self._outline_A_Lp[:, 1]

        assert self.mcl_A_Lp is not None
        mclX_A_Lp = self.mcl_A_Lp[:, 0]
        mclY_A_Lp = self.mcl_A_Lp[:, 1]

        outlineYMin_A_Lp = float(np.min(outlineY_A_Lp))
        outlineYMax_A_Lp = float(np.max(outlineY_A_Lp))
        outlineYRange_A_Lp = outlineYMax_A_Lp - outlineYMin_A_Lp
        y_padding = 0.1 * outlineYRange_A_Lp

        airfoil_axes.plot(outlineX_A_Lp, outlineY_A_Lp, "b-")
        airfoil_axes.plot(mclX_A_Lp, mclY_A_Lp, "r-")

        airfoil_axes.set_xlim(0.0, 1.0)
        airfoil_axes.set_ylim(
            outlineYMin_A_Lp - y_padding, outlineYMax_A_Lp + y_padding
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

    def get_resampled_mcl(
        self, mcl_fractions: np.ndarray | Sequence[float]
    ) -> np.ndarray:
        """Returns a ndarray of points along the mean camber line (MCL), resampled from
        the mcl_A_Lp attribute. It is used to discretize the MCL for meshing.

        :param mcl_fractions: A (N,) array like object of floats representing normalized
            distances along the MCL (from the leading to the trailing edge) at which to
            return the resampled MCL points. Can be a tuple, list, or ndarray. The first
            value must be 0.0, the last must be 1.0, and the remaining must be in the
            range [0.0, 1.0]. All values must be non duplicated and in ascending order.
        :return: A (N,2) ndarray of floats that contains the positions of the resampled
            MCL points (in airfoil axes, relative to the leading point).
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
                "All values in mcl_fractions must be non duplicated and in ascending "
                "order."
            )

        # Find the distance between points along the MCL.
        assert self.mcl_A_Lp is not None
        mclX_A_Lp: np.ndarray = self.mcl_A_Lp[:, 0]
        mclY_A_Lp: np.ndarray = self.mcl_A_Lp[:, 1]
        mcl_distances_between_points = np.sqrt(
            np.power(mclX_A_Lp[:-1] - mclX_A_Lp[1:], 2)
            + np.power(mclY_A_Lp[:-1] - mclY_A_Lp[1:], 2)
        )

        # Create a horizontal 1D array that contains the distance along the MCL of each
        # point on the MCL.
        mcl_distances_cumulative = np.hstack(
            (0, np.cumsum(mcl_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        mcl_distances_cumulative_normalized = (
            mcl_distances_cumulative / mcl_distances_cumulative[-1]
        )

        # Use linear interpolation to resample the MCL.
        mclX_resampled = np.interp(
            mcl_fractions, mcl_distances_cumulative_normalized, self.mcl_A_Lp[:, 0]
        )
        mclY_resampled = np.interp(
            mcl_fractions, mcl_distances_cumulative_normalized, self.mcl_A_Lp[:, 1]
        )

        return np.column_stack([mclX_resampled, mclY_resampled])

    def _Lp_index(self) -> int:
        """Returns the index of the leading point in the _outline_A_Lp attribute.

        :return: The index of the leading point.
        """
        return int(np.argmin(self._outline_A_Lp[:, 0]))

    def _lower_outline(self) -> np.ndarray:
        """Returns a 2D ndarray of points on the lower portion of the Airfoil's outline
        (in airfoil axes, relative to the leading point).

        The order of the returned points is from leading point to trailing edge.
        Included is the leading point, so be careful about duplicates if using this
        method in conjunction with _upper_outline.

        :return: A (N,2) ndarray of floats that describe the position of N points on the
            Airfoil's lower outline (in airfoil axes, relative to the leading point).
        """
        return self._outline_A_Lp[self._Lp_index() :, :]

    def _upper_outline(self) -> np.ndarray:
        """Returns a 2D ndarray of points on the upper portion of the Airfoil's outline
        (in airfoil axes, relative to the leading point).

        The order of the returned points is from trailing edge to leading point.
        Included is the leading point, so be careful about duplicates if using this
        method in conjunction with _lower_outline.

        :return: A (N,2) ndarray of floats that describe the position of N points on the
            Airfoil's upper outline (in airfoil axes, relative to the leading point).
        """
        return self._outline_A_Lp[: self._Lp_index() + 1, :]

    def _populate_outline(self) -> None:
        """Populates a variable with the points of the Airfoil's outline (in airfoil
        axes, relative to the leading point).

        The points are generated if the Airfoil is a NACA 4 series airfoil, or loaded
        from the "airfoils" directory inside "pterasoftware", which is a database of dat
        files containing Airfoil points). NACA 4 series airfoil generation is an
        adaptation of:
        https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil.

        :return: None
        """
        # Sanitize the name input.
        sanitized_name = self._name.lower().strip()

        # Check if the sanitized Airfoil's name matches a name for a NACA 4 series
        # airfoil (NACA0000 is not valid, thickness must be at most 30%, first two
        # digits must both be zero or both be non zero, and for cambered airfoils the
        # position of maximum camber must be at least the maximum camber plus half the
        # maximum thickness). If so, generate it.
        if "naca" in sanitized_name:
            naca_number = sanitized_name.split("naca")[1]
            if naca_number.isdigit():
                if (len(naca_number) == 4) and (naca_number != "0000"):

                    # Parse the characteristics from the name.
                    max_camber = int(naca_number[0]) * 0.01
                    camber_loc = int(naca_number[1]) * 0.1
                    thickness = int(naca_number[2:]) * 0.01

                    # Validate that the thickness is at most 30%.
                    if thickness > 0.30:
                        raise ValueError(
                            "NACA 4 series airfoils with thickness above 30% are not "
                            "supported."
                        )

                    # Reject inconsistent camber parameters. The first two digits must
                    # either both be zero (symmetric) or both be non zero (cambered).
                    if (max_camber > 0) != (camber_loc > 0):
                        raise ValueError(
                            "NACA 4 series airfoils must have consistent camber "
                            "parameters: the first two digits must either both be "
                            "zero (symmetric) or both be non zero (cambered)."
                        )

                    # Reject cambered airfoils where the position of maximum camber is
                    # too close to the leading edge relative to the camber and
                    # thickness. This prevents geometric issues near the leading edge.
                    if max_camber > 0:
                        if camber_loc < max_camber + thickness / 2:
                            raise ValueError(
                                "NACA 4 series airfoils must have the position of "
                                "maximum camber (second digit x 10%) greater than "
                                "or equal to the maximum camber (first digit x 1%) "
                                "plus half the maximum thickness (last two digits "
                                "x 0.5%)."
                            )

                    # Set the number of points per side.
                    n_points_per_side = 400

                    # Get the x component of the MCL.
                    mclX_A_Lp = _functions.cosspace(0, 1, n_points_per_side)

                    # Find the half-thickness of the outline perpendicular to the MCL
                    # (in airfoil axes).
                    halfThickness_A = (
                        5
                        * thickness
                        * (
                            +0.2969 * np.power(mclX_A_Lp, 0.5)
                            - 0.1260 * mclX_A_Lp
                            - 0.3516 * np.power(mclX_A_Lp, 2)
                            + 0.2843 * np.power(mclX_A_Lp, 3)
                            - 0.1015 * np.power(mclX_A_Lp, 4)
                        )
                    )

                    # Prevent divide by zero errors for airfoils like the NACA 0012.
                    if camber_loc == 0:
                        camber_loc = 0.5

                    # Get the y components of the MCL (in airfoil axes, relative to the
                    # leading point).
                    mclY1_A_Lp = (
                        max_camber
                        / camber_loc**2
                        * (
                            2 * camber_loc * mclX_A_Lp[mclX_A_Lp <= camber_loc]
                            - mclX_A_Lp[mclX_A_Lp <= camber_loc] ** 2
                        )
                    )
                    mclY2_A_Lp = (
                        max_camber
                        / (1 - camber_loc) ** 2
                        * (
                            (1 - 2 * camber_loc)
                            + 2 * camber_loc * mclX_A_Lp[mclX_A_Lp > camber_loc]
                            - mclX_A_Lp[mclX_A_Lp > camber_loc] ** 2
                        )
                    )
                    mclY_A_Lp = np.hstack((mclY1_A_Lp, mclY2_A_Lp))

                    # Get the slope of the MCL (in airfoil axes).
                    mclSlope1_A = (
                        2
                        * max_camber
                        / camber_loc**2
                        * (camber_loc - mclX_A_Lp[mclX_A_Lp <= camber_loc])
                    )
                    mclSlope2_A = (
                        2
                        * max_camber
                        / (1 - camber_loc) ** 2
                        * (camber_loc - mclX_A_Lp[mclX_A_Lp > camber_loc])
                    )
                    mclSlope_A = np.hstack((mclSlope1_A, mclSlope2_A))

                    # Convert the slope of the MCL to the angle between the airfoil x
                    # axis and the MCL tangent line.
                    thetaSlopeRad_Ax_to_MCL = np.arctan(mclSlope_A)

                    # Find the upper and lower points of the Airfoil's outline (in
                    # airfoil axes, relative to the leading point) using the MCL points,
                    # the perpendicular half-thickness, and the angle between the x axis
                    # and the MCL tangent line.
                    upperOutlineX_A_Lp = mclX_A_Lp - halfThickness_A * np.sin(
                        thetaSlopeRad_Ax_to_MCL
                    )
                    lowerOutlineX_A_Lp = mclX_A_Lp + halfThickness_A * np.sin(
                        thetaSlopeRad_Ax_to_MCL
                    )
                    upperOutlineY_A_Lp = mclY_A_Lp + halfThickness_A * np.cos(
                        thetaSlopeRad_Ax_to_MCL
                    )
                    lowerOutlineY_A_Lp = mclY_A_Lp - halfThickness_A * np.cos(
                        thetaSlopeRad_Ax_to_MCL
                    )

                    # Flip upper surface so it's back to front.
                    upperOutlineX_A_Lp, upperOutlineY_A_Lp = np.flipud(
                        upperOutlineX_A_Lp
                    ), np.flipud(upperOutlineY_A_Lp)

                    # Trim one point from lower surface so there's no overlap.
                    lowerOutlineX_A_Lp, lowerOutlineY_A_Lp = (
                        lowerOutlineX_A_Lp[1:],
                        lowerOutlineY_A_Lp[1:],
                    )

                    # Combine the points.
                    outlineX_A_Lp = np.hstack((upperOutlineX_A_Lp, lowerOutlineX_A_Lp))
                    outlineY_A_Lp = np.hstack((upperOutlineY_A_Lp, lowerOutlineY_A_Lp))

                    # Populate the _outline_A_Lp attribute and return.
                    self._outline_A_Lp = np.column_stack((outlineX_A_Lp, outlineY_A_Lp))
                    return

        # Try to read from the airfoil directory.
        try:

            # Get the path to the _airfoils data directory.
            airfoils_dir = importlib.resources.files("pterasoftware.geometry").joinpath(
                "_airfoils"
            )

            # Find the file with a case-insensitive match. This is necessary because
            # some filesystems (e.g., Linux) are case-sensitive while others (e.g.,
            # Windows) are not.
            target_name = sanitized_name + ".dat"
            airfoil_file = None
            for item in airfoils_dir.iterdir():
                if item.name.lower() == target_name:
                    airfoil_file = item
                    break

            if airfoil_file is None:
                raise FileNotFoundError(f"Airfoil '{sanitized_name}' not found.")

            # Read the text from the airfoil file.
            raw_text = airfoil_file.read_text()

            # Split into lines, skip the header (first line), and filter out any comment
            # lines (starting with #).
            lines = raw_text.split("\n")[1:]
            data_lines = [
                line for line in lines if line.strip() and not line.startswith("#")
            ]
            trimmed_text = "\n".join(data_lines)

            # Input the coordinates into a 1D array. This represents the upper and lower
            # points of the Airfoil's outline (in airfoil axes, relative to the leading
            # point).
            outline1D_A_Lp = np.fromstring(trimmed_text, sep="\n")

            # Check to make sure the number of elements in the array is even.
            if len(outline1D_A_Lp) % 2 != 0:
                raise ValueError(
                    "name matched to an airfoil in the airfoils database, but it "
                    "could not be read correctly."
                )

            # Populate the _outline_A_Lp attribute and return.
            self._outline_A_Lp = np.reshape(outline1D_A_Lp, (-1, 2))
            return

        # If the Airfoil was not a NACA 4 series and was not found in the database,
        # throw an error.
        except FileNotFoundError:
            raise ValueError(
                "name didn't match a valid NACA 4 series pattern (4 digits, not 0000, "
                "thickness at most 30%, first two digits must both be zero or both be "
                "non zero, and for cambered airfoils the position of maximum camber "
                "must be at least the maximum camber plus half the maximum thickness) "
                "nor was it found in the airfoils database."
            )

    @staticmethod
    def _validate_outline_preliminary(outline_A_Lp: Any) -> np.ndarray:
        """Validates the basic structure of a user's provided outline_A_Lp.

        Only checks for issues that cannot be fixed by normalization. Orientation
        dependent checks (like x monotonicity) are deferred to _validate_outline_final()
        since they must be performed after rotation correction.

        :param outline_A_Lp: The input to validate (can be any type initially).
        :return: The validated version of outline_A_Lp as a (N,2) ndarray of floats.
            Always returns a copy to ensure the Airfoil owns its own data.
        """
        validated_outline_A_Lp = (
            _parameter_validation.arrayLike_of_twoD_number_vectorLikes_return_float(
                outline_A_Lp, "outline_A_Lp"
            )
        )
        # Make a copy to ensure the Airfoil owns its own data. This is important because
        # the input might be a read only array from another Airfoil.
        validated_outline_A_Lp = validated_outline_A_Lp.copy()

        n_outline_points = validated_outline_A_Lp.shape[0]

        # The outline must have at least 5 points.
        if n_outline_points < 5:
            raise ValueError("The Airfoil's outline must have at least five points.")

        # Find the index of the outline's leading point (minimum x).
        outlineLp_index = int(np.argmin(validated_outline_A_Lp[:, 0]))

        # Split the outline into its upper and lower sections.
        lowerOutline_A_Lp = validated_outline_A_Lp[outlineLp_index:, :]
        upperOutline_A_Lp = validated_outline_A_Lp[: outlineLp_index + 1, :]

        # Check that the split portions both have at least three unique points.
        upperOutlineUnique_A_Lp = np.unique(upperOutline_A_Lp, axis=0)
        lowerOutlineUnique_A_Lp = np.unique(lowerOutline_A_Lp, axis=0)
        n_upper_unique = upperOutlineUnique_A_Lp.shape[0]
        n_lower_unique = lowerOutlineUnique_A_Lp.shape[0]
        if n_upper_unique < 3:
            raise ValueError(
                "The upper portion of the Airfoil's outline must contain at least "
                "three unique points (including the outline's leading point)."
            )
        if n_lower_unique < 3:
            raise ValueError(
                "The lower portion of the Airfoil's outline must contain at least "
                "three unique points (including the outline's leading point)."
            )

        return validated_outline_A_Lp

    def _normalize_outline(self) -> None:
        """Transforms the Airfoil's outline to canonical form: leading point at origin,
        chord line on the x axis, and unit chord length.

        Uses an iterative approach to find the true leading point. If the airfoil data
        has an implicit angle of attack (up to 15 degrees), the initial minimum x point
        may not be the true aerodynamic leading edge. After rotating the chord onto the
        x axis, a different point may become the minimum x. The iteration continues
        until the leading point is stable (i.e., the same point remains at minimum x
        after rotation). Outlines with chord angles exceeding 15 degrees are rejected as
        they indicate unexpected data orientation rather than implicit angle of attack.

        :return: None
        """
        max_iterations = 10
        convergence_tol = 1e-9

        for iteration in range(max_iterations):
            # Find the current leading point (minimum x).
            Lp_index = self._Lp_index()
            Lp_A = self._outline_A_Lp[Lp_index, :]

            # Translate the leading point to the origin.
            self._outline_A_Lp[:, 0] -= Lp_A[0]
            self._outline_A_Lp[:, 1] -= Lp_A[1]

            # Find the trailing edge point (average of upper and lower TE points).
            upperTp_A_Lp = self._outline_A_Lp[0, :]
            lowerTp_A_Lp = self._outline_A_Lp[-1, :]
            te_A_Lp = (upperTp_A_Lp + lowerTp_A_Lp) / 2

            # Calculate the angle the chord makes with the x axis.
            chord_angle_rad = np.arctan2(te_A_Lp[1], te_A_Lp[0])

            # Check for excessive rotation on the first iteration. Minor rotation
            # offsets (implicit angle of attack) are acceptable, but large rotations
            # indicate the outline data is not in the expected orientation.
            max_rotation = 15.0
            if iteration == 0 and np.abs(chord_angle_rad) > np.deg2rad(max_rotation):
                raise ValueError(
                    f"The Airfoil's outline has excessive rotation "
                    f"({np.rad2deg(chord_angle_rad):#.3G} deg). The chord line must be "
                    f"within 15 deg of the x axis. Minor rotation offsets (such as "
                    f"implicit angle of attack) are corrected automatically, but the "
                    f"outline data appears to be in an unexpected orientation."
                )

            # Create an active rotation matrix to rotate the chord onto x axis. Convert
            # the angle to degrees to match the _transformations.py standard.
            rot_R_act = _transformations.generate_2D_rot_R(
                angle=np.rad2deg(-chord_angle_rad), passive=False
            )

            # Apply the active rotation to all points.
            self._outline_A_Lp = (rot_R_act @ self._outline_A_Lp.T).T

            # Check if the point at the origin is still the minimum x point.
            new_Lp_index = self._Lp_index()
            new_min_x = self._outline_A_Lp[new_Lp_index, 0]

            # If the minimum x is at or very close to 0.0, we have converged.
            if np.abs(new_min_x) < convergence_tol:
                break
        else:
            # If we exit the loop without breaking, normalization did not converge.
            raise ValueError(
                "Airfoil outline normalization did not converge. The outline may be "
                "malformed."
            )

        # Scale to unit chord.
        upperTp_A_Lp = self._outline_A_Lp[0, :]
        lowerTp_A_Lp = self._outline_A_Lp[-1, :]
        te_A_Lp = (upperTp_A_Lp + lowerTp_A_Lp) / 2
        chord_length = te_A_Lp[0]

        self._outline_A_Lp /= chord_length

        # Clamp x coordinates to [0.0, 1.0] and extrapolate trailing edge points to
        # exactly x = 1.0 if needed.
        self._outline_A_Lp[:, 0] = np.clip(self._outline_A_Lp[:, 0], 0.0, 1.0)

        # Extrapolate upper TE to x = 1.0 if it's not already there.
        if self._outline_A_Lp[0, 0] < 1.0:
            # Use linear extrapolation between the first two points to find y at x =
            # 1.0.
            x0, y0 = self._outline_A_Lp[1, :]
            x1, y1 = self._outline_A_Lp[0, :]
            if x1 > x0:
                y_at_1 = y0 + (y1 - y0) * (1.0 - x0) / (x1 - x0)
                self._outline_A_Lp[0, :] = [1.0, y_at_1]

        # Extrapolate lower TE to x = 1.0 if it's not already there.
        if self._outline_A_Lp[-1, 0] < 1.0:
            # Use linear extrapolation between the last two points to find y at x = 1.0.
            x0, y0 = self._outline_A_Lp[-2, :]
            x1, y1 = self._outline_A_Lp[-1, :]
            if x1 > x0:
                y_at_1 = y0 + (y1 - y0) * (1.0 - x0) / (x1 - x0)
                self._outline_A_Lp[-1, :] = [1.0, y_at_1]

        # Remove duplicate interior points while preserving first and last points. For
        # closed trailing edges, the upper and lower TE points may be identical, but
        # both must be kept to maintain proper outline topology.
        interior = self._outline_A_Lp[1:-1]
        _, unique_indices = np.unique(interior, axis=0, return_index=True)
        unique_indices = np.sort(unique_indices)
        unique_interior = interior[unique_indices]
        self._outline_A_Lp = np.vstack(
            [self._outline_A_Lp[0:1], unique_interior, self._outline_A_Lp[-1:]]
        )

        # Correct small trailing edge inversions. After normalization, floating-point
        # precision or minor data quality issues may cause the upper TE y to be slightly
        # below the lower TE y. If the inversion is small (<=0.1% chord), correct it by
        # averaging the y values to create a closed trailing edge. Larger inversions
        # indicate genuinely malformed data and are left for _validate_outline_final()
        # to catch.
        upperTpY_A_Lp = self._outline_A_Lp[0, 1]
        lowerTpY_A_Lp = self._outline_A_Lp[-1, 1]
        if lowerTpY_A_Lp > upperTpY_A_Lp:
            te_inversion = lowerTpY_A_Lp - upperTpY_A_Lp
            if te_inversion <= 1e-3:
                avg_te_y = (upperTpY_A_Lp + lowerTpY_A_Lp) / 2
                self._outline_A_Lp[0, 1] = avg_te_y
                self._outline_A_Lp[-1, 1] = avg_te_y

    def _validate_outline_final(self) -> None:
        """Verifies that outline normalization succeeded, checks x monotonicity, and
        checks for self intersection.

        :return: None
        """
        # Check for NaN values.
        if np.any(np.isnan(self._outline_A_Lp)):
            raise ValueError(
                "The Airfoil's outline contains NaN values after normalization."
            )

        # Check that the leading point is at approximately [0.0, 0.0].
        Lp_index = self._Lp_index()
        Lp_A_Lp = self._outline_A_Lp[Lp_index, :]
        tol = 1e-6
        if not np.isclose(Lp_A_Lp[0], 0.0, atol=tol):
            raise ValueError(
                "The Airfoil's outline's leading point x value must be approximately "
                "0.0 after normalization."
            )
        if not np.isclose(Lp_A_Lp[1], 0.0, atol=tol):
            raise ValueError(
                "The Airfoil's outline's leading point y value must be approximately "
                "0.0 after normalization."
            )

        # Check that trailing edge x values are approximately 1.0.
        upperTpX_A_Lp = self._outline_A_Lp[0, 0]
        lowerTpX_A_Lp = self._outline_A_Lp[-1, 0]
        te_tol = 1e-6
        if not np.isclose(upperTpX_A_Lp, 1.0, atol=te_tol):
            raise ValueError(
                "The Airfoil's upper outline's trailing point x value must be "
                "approximately 1.0 after normalization."
            )
        if not np.isclose(lowerTpX_A_Lp, 1.0, atol=te_tol):
            raise ValueError(
                "The Airfoil's lower outline's trailing point x value must be "
                "approximately 1.0 after normalization."
            )

        # Check that upper TE y >= lower TE y.
        upperTpY_A_Lp = self._outline_A_Lp[0, 1]
        lowerTpY_A_Lp = self._outline_A_Lp[-1, 1]
        if lowerTpY_A_Lp > upperTpY_A_Lp:
            raise ValueError(
                "The Airfoil's upper outline's trailing point y value must be "
                "greater than or equal to the lower outline's trailing point y value."
            )

        # Check x monotonicity. These checks are performed here (after normalization)
        # rather than in preliminary validation because airfoils with implicit angle of
        # attack may not have monotonic x values until after rotation correction.
        upperOutline_A_Lp = self._upper_outline()
        lowerOutline_A_Lp = self._lower_outline()

        upperOutlineDiffX_A = np.diff(upperOutline_A_Lp[:, 0])
        lowerOutlineDiffX_A = np.diff(lowerOutline_A_Lp[:, 0])

        # Check that the upper outline is non increasing in x and that the lower outline
        # is non decreasing in x (in airfoil axes). Adjacent points may have the same x
        # value.
        if not np.all(upperOutlineDiffX_A <= 0.0):
            raise ValueError(
                "Every point in the Airfoil's outline's upper portion must have "
                "an x value less than or equal to the point before it (in airfoil "
                "axes)."
            )
        if not np.all(lowerOutlineDiffX_A >= 0.0):
            raise ValueError(
                "Every point in the Airfoil's outline's lower portion must have "
                "an x value greater than or equal to the point before it (in airfoil "
                "axes)."
            )

        # Check for self intersection using linear interpolation.

        # Flip upper outline so it goes from LE to TE.
        flippedUpperOutline_A_Lp = np.flipud(upperOutline_A_Lp)

        # Remove duplicate x values from both outlines. Keep the first occurrence.
        _, upper_unique_indices = np.unique(
            flippedUpperOutline_A_Lp[:, 0], return_index=True
        )
        upper_unique_indices = np.sort(upper_unique_indices)
        flippedUpperOutlineUnique_A_Lp = flippedUpperOutline_A_Lp[upper_unique_indices]

        _, lower_unique_indices = np.unique(lowerOutline_A_Lp[:, 0], return_index=True)
        lower_unique_indices = np.sort(lower_unique_indices)
        lowerOutlineUnique_A_Lp = lowerOutline_A_Lp[lower_unique_indices]

        # Sample at all unique x values from both surfaces.
        all_x = np.union1d(
            flippedUpperOutlineUnique_A_Lp[:, 0], lowerOutlineUnique_A_Lp[:, 0]
        )

        # Use linear interpolation.
        upperY_interp = np.interp(
            all_x,
            flippedUpperOutlineUnique_A_Lp[:, 0],
            flippedUpperOutlineUnique_A_Lp[:, 1],
        )
        lowerY_interp = np.interp(
            all_x,
            lowerOutlineUnique_A_Lp[:, 0],
            lowerOutlineUnique_A_Lp[:, 1],
        )

        upperMinusLowerOutlineY = upperY_interp - lowerY_interp

        # Check for self intersection. Upper y must be strictly greater than lower y for
        # all interior points (excluding endpoints). This rejects zero thickness
        # airfoils.
        if not np.all(upperMinusLowerOutlineY[1:-1] > 0.0):
            raise ValueError(
                "All points on the Airfoil's upper outline (excluding the outline's "
                "leading and trailing points) must have a y value strictly greater "
                "than the points on the lower outline at the same x value "
                "(in airfoil axes, relative to the leading point)."
            )

    def _resample_outline(self, n_points_per_side: int) -> None:
        """Returns a resampled version of the points on the Airfoil's outline (in
        airfoil axes, relative to the leading point) with cosine spaced points on the
        upper and lower surfaces.

        The number of points defining the final Airfoil's outline is (n_points_per_side
        * 2 - 1), since the leading point is shared by both the upper and lower
        surfaces.

        :param n_points_per_side: The number of points (a positive int) on the upper and
            lower surfaces.
        :return: None
        """
        # Get the upper outline points. This contains the leading point.
        upperOutline_A_Lp = self._upper_outline()

        flippedUpperOutlineX_A_Lp = np.flipud(upperOutline_A_Lp)[:, 0]
        flippedUpperOutlineY_A_Lp = np.flipud(upperOutline_A_Lp)[:, 1]

        # Find the distance between points along the upper flipped original outline.
        flippedUpperOutline_distances_between_points = np.sqrt(
            np.power(
                flippedUpperOutlineX_A_Lp[:-1] - flippedUpperOutlineX_A_Lp[1:],
                2,
            )
            + np.power(
                flippedUpperOutlineY_A_Lp[:-1] - flippedUpperOutlineY_A_Lp[1:],
                2,
            )
        )

        # Create a horizontal 1D array that contains the cumulative distance along the
        # upper flipped original outline of each point.
        flippedUpperOutline_distances_cumulative = np.hstack(
            (0, np.cumsum(flippedUpperOutline_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        flippedUpperOutline_distances_cumulative_normalized = (
            flippedUpperOutline_distances_cumulative
            / flippedUpperOutline_distances_cumulative[-1]
        )

        # Get the lower outline points. This contains the leading point.
        lowerOutline_A_Lp = self._lower_outline()

        lowerOutlineX_A_Lp = lowerOutline_A_Lp[:, 0]
        lowerOutlineY_A_Lp = lowerOutline_A_Lp[:, 1]

        # Find the distance between points along the lower original outline.
        lowerOutline_distances_between_points = np.sqrt(
            np.power(
                lowerOutlineX_A_Lp[:-1] - lowerOutlineX_A_Lp[1:],
                2,
            )
            + np.power(
                lowerOutlineY_A_Lp[:-1] - lowerOutlineY_A_Lp[1:],
                2,
            )
        )

        # Create a horizontal 1D array that contains the cumulative distance along the
        # lower original outline of each point.
        lowerOutline_distances_cumulative = np.hstack(
            (0, np.cumsum(lowerOutline_distances_between_points))
        )

        # Normalize the 1D array so that it ranges from 0.0 to 1.0.
        lowerOutline_distances_cumulative_normalized = (
            lowerOutline_distances_cumulative / lowerOutline_distances_cumulative[-1]
        )

        # Generate a cosine spaced list of normalized distances from 0.0 to 1.0.
        cosine_spaced_normalized_distances = _functions.cosspace(
            0.0, 1.0, n_points_per_side
        )

        # Use linear interpolation to find the x and y components of the upper and lower
        # outline points at each of the resampled cosine spaced distances.
        upperResampledOutlineX_A_Lp = np.flipud(
            np.interp(
                cosine_spaced_normalized_distances,
                flippedUpperOutline_distances_cumulative_normalized,
                flippedUpperOutlineX_A_Lp,
            )
        )
        upperResampledOutlineY_A_Lp = np.flipud(
            np.interp(
                cosine_spaced_normalized_distances,
                flippedUpperOutline_distances_cumulative_normalized,
                flippedUpperOutlineY_A_Lp,
            )
        )
        lowerResampledOutlineX_A_Lp = np.interp(
            cosine_spaced_normalized_distances,
            lowerOutline_distances_cumulative_normalized,
            lowerOutlineX_A_Lp,
        )[1:]
        lowerResampledOutlineY_A_Lp = np.interp(
            cosine_spaced_normalized_distances,
            lowerOutline_distances_cumulative_normalized,
            lowerOutlineY_A_Lp,
        )[1:]

        resampledOutlineX_A_Lp = np.hstack(
            (upperResampledOutlineX_A_Lp, lowerResampledOutlineX_A_Lp)
        )
        resampledOutlineY_A_Lp = np.hstack(
            (upperResampledOutlineY_A_Lp, lowerResampledOutlineY_A_Lp)
        )

        self._outline_A_Lp = np.column_stack(
            (resampledOutlineX_A_Lp, resampledOutlineY_A_Lp)
        )

    def _populate_mcl(self) -> None:
        """Creates a 2D ndarray of points along the Airfoil's MCL (in airfoil axes,
        relative to the leading point), which it uses to set the mcl_A_Lp attribute. It
        is in order from the leading point to the trailing point.

        :return: None
        """
        # Split outline_A_Lp into upper and lower sections. Flip the upper points so
        # that they are ordered from the leading point to the trailing point.
        flippedUpperOutline_A_Lp = np.flipud(self._upper_outline())
        lowerOutline_A_Lp = self._lower_outline()

        # Generate cosine spaced x fractions from 0.0 to 1.0.
        cosine_spaced_chord_fractions = _functions.cosspace(
            0.0, 1.0, self._n_points_per_side
        )

        # Use linear interpolation to find the y values of the upper and lower outlines
        # at the cosine spaced x positions.
        flippedUpperOutlineY_A_Lp = np.interp(
            cosine_spaced_chord_fractions,
            flippedUpperOutline_A_Lp[:, 0],
            flippedUpperOutline_A_Lp[:, 1],
        )
        lowerOutlineY_A_Lp = np.interp(
            cosine_spaced_chord_fractions,
            lowerOutline_A_Lp[:, 0],
            lowerOutline_A_Lp[:, 1],
        )

        # Calculate the approximate MCL points (in airfoil axes, relative to the leading
        # point) and set the class attribute.
        self._mcl_A_Lp = np.column_stack(
            [
                cosine_spaced_chord_fractions,
                (flippedUpperOutlineY_A_Lp + lowerOutlineY_A_Lp) / 2,
            ]
        )

        # Resample the MCL points using cosine spaced distances along the MCL. The
        # fractions here represent normalized arc length (0.0 to 1.0).
        normalized_mcl_fractions = _functions.cosspace(
            0.0, 1.0, self._n_points_per_side
        )
        self._mcl_A_Lp = self.get_resampled_mcl(mcl_fractions=normalized_mcl_fractions)

        # Normalize the MCL so that x spans from 0.0 to 1.0. This corrects for slight
        # variations in the outline data where upper and lower surfaces may not extend
        # to exactly x = 0.0 and x = 1.0. We translate the leading edge to the origin
        # and scale the x values only (not y) to put the trailing edge at x = 1.0. We
        # intentionally do not rotate to put the trailing edge on the x axis, because
        # that would remove control surface deflection effects from the MCL.

        # Step 1: Translate so leading edge is at origin.
        self._mcl_A_Lp[:, 0] -= self._mcl_A_Lp[0, 0]
        self._mcl_A_Lp[:, 1] -= self._mcl_A_Lp[0, 1]

        # Step 2: Scale x only so trailing edge is at x = 1.0.
        x_scale_factor = 1.0 / self._mcl_A_Lp[-1, 0]
        self._mcl_A_Lp[:, 0] *= x_scale_factor
