"""Contains the FreeFlightWingCrossSectionMovement class.

**Contains the following classes:**

FreeFlightWingCrossSectionMovement: A class used to contain a WingCrossSection's
movement in a free flight simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from .. import _core, geometry


class FreeFlightWingCrossSectionMovement(_core.CoreWingCrossSectionMovement):
    """A class used to contain a WingCrossSection's movement in a free flight
    simulation.

    In free flight, wing cross section geometry is prescribed (the same oscillation
    based generation as WingCrossSectionMovement). This class exists so that a
    FreeFlightWingMovement always accepts FreeFlightWingCrossSectionMovements, keeping
    the free flight movement hierarchy consistently named.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this FreeFlightWingCrossSectionMovement.

    all_periods: All unique non zero periods from this
    FreeFlightWingCrossSectionMovement.

    max_period: FreeFlightWingCrossSectionMovement's longest period of motion.

    generate_wing_cross_section_at_time_step: Creates the WingCrossSection at a single
    time step.

    generate_wing_cross_sections: Creates the WingCrossSection at each time step, and
    returns them in a list.
    """

    __slots__ = ()

    def __init__(
        self,
        base_wing_cross_section: geometry.wing_cross_section.WingCrossSection,
        ampLp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseLp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        periodAngles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        spacingAngles_Wcsp_to_Wcs_ixyz: (
            np.ndarray | Sequence[str | Callable[[float], float]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseAngles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
    ) -> None:
        """The initialization method.

        :param base_wing_cross_section: The base WingCrossSection from which the
            WingCrossSection at each time step will be created.
        :param ampLp_Wcsp_Lpp: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            FreeFlightWingCrossSectionMovement's changes in its WingCrossSections'
            Lp_Wcsp_Lpp parameters. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Each amplitude must be low enough that it
            doesn't drive its base value out of the range of valid values. Otherwise,
            this FreeFlightWingCrossSectionMovement will try to create WingCrossSections
            with invalid parameter values. The units are in meters. The default is (0.0,
            0.0, 0.0).
        :param periodLp_Wcsp_Lpp: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            FreeFlightWingCrossSectionMovement's changes in its WingCrossSections'
            Lp_Wcsp_Lpp parameters. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampLp_Wcsp_Lpp is 0.0 and non zero if not. The
            units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingLp_Wcsp_Lpp: An array-like object of strs or callables with shape
            (3,) representing the spacing of the FreeFlightWingCrossSectionMovement's
            changes in its WingCrossSections' Lp_Wcsp_Lpp parameters. Can be a tuple,
            list, or ndarray. Each element can be the str "sine", the str "uniform", or
            a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0.0, return to 0.0 after one period of 2.0
            * pi radians, have amplitude of 1.0, be periodic, return finite values only,
            and accept a float as input and return a float. Custom functions are scaled
            by ampLp_Wcsp_Lpp, shifted horizontally and vertically by phaseLp_Wcsp_Lpp
            and the base value, and have a period set by periodLp_Wcsp_Lpp. The default
            is ("sine", "sine", "sine").
        :param phaseLp_Wcsp_Lpp: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's WingCrossSection's Lp_Wcsp_Lpp parameter relative to the base
            WingCrossSection's Lp_Wcsp_Lpp parameter. Can be a tuple, list, or ndarray.
            Elements must lie in the range (-180.0, 180.0]. Each element must be 0.0 if
            the corresponding element in ampLp_Wcsp_Lpp is 0.0 and non zero if not.
            Values are converted to floats internally. The units are in degrees. The
            default is (0.0, 0.0, 0.0).
        :param ampAngles_Wcsp_to_Wcs_ixyz: An array-like object of non negative numbers
            (int or float) with shape (3,) representing the amplitudes of the
            FreeFlightWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Values
            are converted to floats internally. Each amplitude must be low enough that
            it doesn't drive its base value out of the range of valid values. Otherwise,
            this FreeFlightWingCrossSectionMovement will try to create WingCrossSections
            with invalid parameter values. The units are in degrees. The default is
            (0.0, 0.0, 0.0).
        :param periodAngles_Wcsp_to_Wcs_ixyz: An array-like object of non negative
            numbers (int or float) with shape (3,) representing the periods of the
            FreeFlightWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Values
            are converted to floats internally. Each element must be 0.0 if the
            corresponding element in ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non zero if
            not. The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Wcsp_to_Wcs_ixyz: An array-like object of strs or callables
            with shape (3,) representing the spacing of the
            FreeFlightWingCrossSectionMovement's changes in its WingCrossSections'
            angles_Wcsp_to_Wcs_ixyz parameters. Can be a tuple, list, or ndarray. Each
            element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2.0 * pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a float
            as input and return a float. Custom functions are scaled by
            ampAngles_Wcsp_to_Wcs_ixyz, shifted horizontally and vertically by
            phaseAngles_Wcsp_to_Wcs_ixyz and the base value, and have a period set by
            periodAngles_Wcsp_to_Wcs_ixyz. The default is ("sine", "sine", "sine").
        :param phaseAngles_Wcsp_to_Wcs_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the phase offsets of the elements in the
            first time step's WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter
            relative to the base WingCrossSection's angles_Wcsp_to_Wcs_ixyz parameter.
            Can be a tuple, list, or ndarray. Elements must lie in the range (-180.0,
            180.0]. Each element must be 0.0 if the corresponding element in
            ampAngles_Wcsp_to_Wcs_ixyz is 0.0 and non zero if not. Values are converted
            to floats internally. The units are in degrees. The default is (0.0, 0.0,
            0.0).
        :return: None
        """
        super().__init__(
            base_wing_cross_section=base_wing_cross_section,
            ampLp_Wcsp_Lpp=ampLp_Wcsp_Lpp,
            periodLp_Wcsp_Lpp=periodLp_Wcsp_Lpp,
            spacingLp_Wcsp_Lpp=spacingLp_Wcsp_Lpp,
            phaseLp_Wcsp_Lpp=phaseLp_Wcsp_Lpp,
            ampAngles_Wcsp_to_Wcs_ixyz=ampAngles_Wcsp_to_Wcs_ixyz,
            periodAngles_Wcsp_to_Wcs_ixyz=periodAngles_Wcsp_to_Wcs_ixyz,
            spacingAngles_Wcsp_to_Wcs_ixyz=spacingAngles_Wcsp_to_Wcs_ixyz,
            phaseAngles_Wcsp_to_Wcs_ixyz=phaseAngles_Wcsp_to_Wcs_ixyz,
        )
