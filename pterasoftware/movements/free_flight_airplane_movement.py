"""Contains the FreeFlightAirplaneMovement class.

**Contains the following classes:**

FreeFlightAirplaneMovement: A class used to contain an Airplane's movement in a free
flight simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import cast

import numpy as np

from .. import _core, geometry
from . import free_flight_wing_movement as free_flight_wing_movement_mod


class FreeFlightAirplaneMovement(_core.CoreAirplaneMovement):
    """A class used to contain an Airplane's movement in a free flight simulation.

    In free flight, airplane geometry is prescribed (the same oscillation based
    generation as AirplaneMovement). This class exists so that a FreeFlightMovement
    always accepts FreeFlightAirplaneMovements, keeping the free flight movement
    hierarchy consistently named.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this FreeFlightAirplaneMovement.

    all_periods: All unique non zero periods from this FreeFlightAirplaneMovement, its
    FreeFlightWingMovement(s), and their FreeFlightWingCrossSectionMovements.

    max_period: The longest period of FreeFlightAirplaneMovement's own motion, the
    motion(s) of its sub movement object(s), and the motions of its sub sub movement
    objects.

    generate_airplane_at_time_step: Creates the Airplane at a single time step.

    generate_airplanes: Creates the Airplane at each time step, and returns them in a
    list.
    """

    __slots__ = ()

    def __init__(
        self,
        base_airplane: geometry.airplane.Airplane,
        wing_movements: list[free_flight_wing_movement_mod.FreeFlightWingMovement],
        ampCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        periodCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1: np.ndarray | Sequence[str | Callable[[float], float]] = (
            "sine",
            "sine",
            "sine",
        ),
        phaseCg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
    ) -> None:
        """The initialization method.

        :param base_airplane: The base Airplane from which the Airplane at each time
            step will be created.
        :param wing_movements: A list of the FreeFlightWingMovements associated with
            each of the base Airplane's Wings. It must have the same length as the base
            Airplane's list of Wings.
        :param ampCg_GP1_CgP1: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the amplitudes of the
            FreeFlightAirplaneMovement's changes in its Airplanes' Cg_GP1_CgP1
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each amplitude must be low enough that it doesn't drive its base
            value out of the range of valid values. Otherwise, this
            FreeFlightAirplaneMovement will try to create Airplanes with invalid
            parameter values. Because the first Airplane's Cg_GP1_CgP1 parameter must be
            all zeros, this means that the first Airplane's ampCg_GP1_CgP1 parameter
            must also be all zeros. The units are in meters. The default is (0.0, 0.0,
            0.0).
        :param periodCg_GP1_CgP1: An array-like object of non negative numbers (int or
            float) with shape (3,) representing the periods of the
            FreeFlightAirplaneMovement's changes in its Airplanes' Cg_GP1_CgP1
            parameters. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. Each element must be 0.0 if the corresponding element in
            ampCg_GP1_CgP1 is 0.0 and non zero if not. The units are in seconds. The
            default is (0.0, 0.0, 0.0).
        :param spacingCg_GP1_CgP1: An array-like object of strs or callables with shape
            (3,) representing the spacing of the FreeFlightAirplaneMovement's changes in
            its Airplanes' Cg_GP1_CgP1 parameters. Can be a tuple, list, or ndarray.
            Each element can be the str "sine", the str "uniform", or a callable custom
            spacing function. Custom spacing functions are for advanced users and must
            start at 0.0, return to 0.0 after one period of 2.0 * pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a float
            as input and return a float. Custom functions are scaled by ampCg_GP1_CgP1,
            shifted horizontally and vertically by phaseCg_GP1_CgP1 and the base value,
            and have a period set by periodCg_GP1_CgP1. The default is ("sine", "sine",
            "sine").
        :param phaseCg_GP1_CgP1: An array-like object of numbers (int or float) with
            shape (3,) representing the phase offsets of the elements in the first time
            step's Airplane's Cg_GP1_CgP1 parameter relative to the base Airplane's
            Cg_GP1_CgP1 parameter. Can be a tuple, list, or ndarray. Elements must lie
            in the range (-180.0, 180.0]. Each element must be 0.0 if the corresponding
            element in ampCg_GP1_CgP1 is 0.0 and non zero if not. Values are converted
            to floats internally. The units are in degrees. The default is (0.0, 0.0,
            0.0).
        :return: None
        """
        # Validate that every element is a FreeFlightWingMovement, not just a
        # CoreWingMovement. CoreAirplaneMovement.__init__() validates at the Core
        # level, but FreeFlightAirplaneMovement enforces the stricter type.
        for wing_movement in wing_movements:
            if not isinstance(
                wing_movement,
                free_flight_wing_movement_mod.FreeFlightWingMovement,
            ):
                raise TypeError(
                    "Every element in wing_movements must be a "
                    "FreeFlightWingMovement."
                )

        super().__init__(
            base_airplane=base_airplane,
            wing_movements=wing_movements,
            ampCg_GP1_CgP1=ampCg_GP1_CgP1,
            periodCg_GP1_CgP1=periodCg_GP1_CgP1,
            spacingCg_GP1_CgP1=spacingCg_GP1_CgP1,
            phaseCg_GP1_CgP1=phaseCg_GP1_CgP1,
        )

    # --- Immutable: read only properties ---
    @property
    def wing_movements(
        self,
    ) -> tuple[free_flight_wing_movement_mod.FreeFlightWingMovement, ...]:
        return cast(
            tuple[free_flight_wing_movement_mod.FreeFlightWingMovement, ...],
            self._wing_movements,
        )
