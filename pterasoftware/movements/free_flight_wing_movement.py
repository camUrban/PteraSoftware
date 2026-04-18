"""Contains the FreeFlightWingMovement class.

**Contains the following classes:**

FreeFlightWingMovement: A class used to contain a Wing's movement in a free flight
simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np

from .. import _core, geometry
from . import (
    free_flight_wing_cross_section_movement as free_flight_wing_cross_section_movement_mod,
)


class FreeFlightWingMovement(_core.CoreWingMovement):
    """A class used to contain a Wing's movement in a free flight simulation.

    In free flight, wing geometry is prescribed (the same oscillation based generation
    as WingMovement). This class exists so that a FreeFlightAirplaneMovement always
    accepts FreeFlightWingMovements, keeping the free flight movement hierarchy
    consistently named.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this FreeFlightWingMovement.

    all_periods: All unique non zero periods from this FreeFlightWingMovement and its
    FreeFlightWingCrossSectionMovements.

    max_period: The longest period of FreeFlightWingMovement's own motion and that of
    its sub movement objects.

    generate_wing_at_time_step: Creates the Wing at a single time step.

    generate_wings: Creates the Wing at each time step, and returns them in a list.

    **Notes:**

    Wings cannot undergo motion that causes them to switch symmetry types. A transition
    between types could change the number of Wings and the Panel structure, which is
    incompatible with the unsteady solver. This happens when a FreeFlightWingMovement
    defines motion that causes its base Wing's wing axes' yz plane and its symmetry
    plane to transition from coincident to non coincident, or vice versa. This is
    checked by this FreeFlightWingMovement's parent FreeFlightAirplaneMovement's parent
    FreeFlightMovement.
    """

    __slots__ = ()

    def __init__(
        self,
        base_wing: geometry.wing.Wing,
        wing_cross_section_movements: list[
            free_flight_wing_cross_section_movement_mod.FreeFlightWingCrossSectionMovement
        ],
        ampLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseLer_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        spacingAngles_Gs_to_Wn_ixyz: (
            np.ndarray | Sequence[str | Callable[[float], float]]
        ) = (
            "sine",
            "sine",
            "sine",
        ),
        phaseAngles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
        rotationPointOffset_Gs_Ler: np.ndarray | Sequence[float | int] = (
            0.0,
            0.0,
            0.0,
        ),
    ) -> None:
        """The initialization method.

        :param base_wing: The base Wing from which the Wing at each time step will be
            created.
        :param wing_cross_section_movements: A list of
            FreeFlightWingCrossSectionMovements associated with each of the base Wing's
            WingCrossSections. It must have the same length as the base Wing's list of
            WingCrossSections.
        :param ampLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            FreeFlightWingMovement's changes in its Wings' Ler_Gs_Cgs parameters. Can be
            a tuple, list, or ndarray. Values are converted to floats internally. Each
            amplitude must be low enough that it doesn't drive its base value out of the
            range of valid values. Otherwise, this FreeFlightWingMovement will try to
            create Wings with invalid parameters values. The units are in meters. The
            default is (0.0, 0.0, 0.0).
        :param periodLer_Gs_Cgs: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            FreeFlightWingMovement's changes in its Wings' Ler_Gs_Cgs parameters. Can be
            a tuple, list, or ndarray. Values are converted to floats internally. Each
            element must be 0.0 if the corresponding element in ampLer_Gs_Cgs is 0.0 and
            non zero if not. The units are in seconds. The default is (0.0, 0.0, 0.0).
        :param spacingLer_Gs_Cgs: An array-like object of strs or callables with shape
            (3,) representing the spacing of the FreeFlightWingMovement's change in its
            Wings' Ler_Gs_Cgs parameters. Can be a tuple, list, or ndarray. Each element
            can be the string "sine", the string "uniform", or a callable custom spacing
            function. Custom spacing functions are for advanced users and must start at
            0.0, return to 0.0 after one period of 2.0 * pi radians, have amplitude of
            1.0, be periodic, return finite values only, and accept a float as input and
            return a float. The custom function is scaled by ampLer_Gs_Cgs, shifted
            horizontally and vertically by phaseLer_Gs_Cgs and the base value, and have
            a period set by periodLer_Gs_Cgs. The default is ("sine", "sine", "sine").
        :param phaseLer_Gs_Cgs: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's Wing's Ler_Gs_Cgs parameter relative to the base Wing's Ler_Gs_Cgs
            parameter. Can be a tuple, list, or ndarray. Values must lie in the range
            (-180.0, 180.0] and will be converted to floats internally. Each element
            must be 0.0 if the corresponding element in ampLer_Gs_Cgs is 0.0 and non
            zero if not. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param ampAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or float)
            with shape (3,) representing the amplitudes of the FreeFlightWingMovement's
            changes in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Values must lie in the range [0.0, 180.0] and will be converted
            to floats internally. Each amplitude must be low enough that it doesn't
            drive its base value out of the range of valid values. Otherwise, this
            FreeFlightWingMovement will try to create Wings with invalid parameters
            values. The units are in degrees. The default is (0.0, 0.0, 0.0).
        :param periodAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or
            float) with shape (3,) representing the periods of the
            FreeFlightWingMovement's changes in its Wings' angles_Gs_to_Wn_ixyz
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampAngles_Gs_to_Wn_ixyz is 0.0 and non zero if not. The units are in
            seconds. The default is (0.0, 0.0, 0.0).
        :param spacingAngles_Gs_to_Wn_ixyz: An array-like object of strs or callables
            with shape (3,) representing the spacing of the FreeFlightWingMovement's
            change in its Wings' angles_Gs_to_Wn_ixyz parameters. Can be a tuple, list,
            or ndarray. Each element can be the string "sine", the string "uniform", or
            a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0.0, return to 0.0 after one period of 2.0
            * pi radians, have amplitude of 1.0, be periodic, return finite values only,
            and accept a float as input and return a float. The custom function is
            scaled by ampAngles_Gs_to_Wn_ixyz, shifted horizontally and vertically by
            phaseAngles_Gs_to_Wn_ixyz and the base value, with the period set by
            periodAngles_Gs_to_Wn_ixyz. The default is ("sine", "sine", "sine").
        :param phaseAngles_Gs_to_Wn_ixyz: An array-like object of numbers (int or float)
            with shape (3,) representing the phase offsets of the elements in the first
            time step's Wing's angles_Gs_to_Wn_ixyz parameter relative to the base
            Wing's angles_Gs_to_Wn_ixyz parameter. Can be a tuple, list, or ndarray.
            Values must lie in the range (-180.0, 180.0] and will be converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampAngles_Gs_to_Wn_ixyz is 0.0 and non zero if not. The units are in
            degrees. The default is (0.0, 0.0, 0.0).
        :param rotationPointOffset_Gs_Ler: An array-like object of 3 numbers (int or
            float) representing the position of the rotation point for the Wing's
            angular motion (in geometry axes after accounting for symmetry, relative to
            the leading edge root point). Can be a tuple, list, or ndarray. Values are
            converted to floats internally. This offset defines where the Wing rotates
            about when angles_Gs_to_Wn_ixyz oscillates. When set to (0, 0, 0), rotation
            occurs about the leading edge root point (default behavior). The units are
            in meters. The default is (0.0, 0.0, 0.0).
        :return: None
        """
        # Validate that every element is a FreeFlightWingCrossSectionMovement, not
        # just a CoreWingCrossSectionMovement. CoreWingMovement.__init__() validates
        # at the Core level, but FreeFlightWingMovement enforces the stricter type.
        for wing_cross_section_movement in wing_cross_section_movements:
            if not isinstance(
                wing_cross_section_movement,
                free_flight_wing_cross_section_movement_mod.FreeFlightWingCrossSectionMovement,
            ):
                raise TypeError(
                    "Every element in wing_cross_section_movements must "
                    "be a FreeFlightWingCrossSectionMovement."
                )

        super().__init__(
            base_wing=base_wing,
            wing_cross_section_movements=wing_cross_section_movements,
            ampLer_Gs_Cgs=ampLer_Gs_Cgs,
            periodLer_Gs_Cgs=periodLer_Gs_Cgs,
            spacingLer_Gs_Cgs=spacingLer_Gs_Cgs,
            phaseLer_Gs_Cgs=phaseLer_Gs_Cgs,
            ampAngles_Gs_to_Wn_ixyz=ampAngles_Gs_to_Wn_ixyz,
            periodAngles_Gs_to_Wn_ixyz=periodAngles_Gs_to_Wn_ixyz,
            spacingAngles_Gs_to_Wn_ixyz=spacingAngles_Gs_to_Wn_ixyz,
            phaseAngles_Gs_to_Wn_ixyz=phaseAngles_Gs_to_Wn_ixyz,
            rotationPointOffset_Gs_Ler=rotationPointOffset_Gs_Ler,
        )
