"""Contains the OperatingPointMovement class.

**Contains the following classes:**

OperatingPointMovement: A class used to contain an OperatingPoint's movements.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from .. import _parameter_validation
from .. import operating_point as operating_point_mod
from . import _functions


class OperatingPointMovement:
    """A class used to contain an OperatingPoint's movements.

    **Contains the following methods:**

    max_period: OperatingPointMovement's longest period of motion.

    generate_operating_points: Creates the OperatingPoint at each time step, and returns
    them in a list.
    """

    def __init__(
        self,
        base_operating_point: operating_point_mod.OperatingPoint,
        ampVCg__E: float | int = 0.0,
        periodVCg__E: float | int = 0.0,
        spacingVCg__E: str | Callable[[np.ndarray], np.ndarray] = "sine",
        phaseVCg__E: float | int = 0.0,
    ) -> None:
        """The initialization method.

        :param base_operating_point: The base OperatingPoint from which the
            OperatingPoint at each time step will be created.
        :param ampVCg__E: The amplitude of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameters. Must be a non negative number (int or
            float), and is converted to a float internally. The amplitude must be low
            enough that it doesn't drive its base value out of the range of valid
            values. Otherwise, this OperatingPointMovement will try to create
            OperatingPoints with invalid parameter values. The units are in meters per
            second. The default is 0.0.
        :param periodVCg__E: The period of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameter. Must be a non negative number (int or
            float), and is converted to a float internally. It must be 0.0 if ampVCg__E
            is 0.0 and non zero if not. The units are in seconds. The default is 0.0.
        :param spacingVCg__E: Determines the spacing of the OperatingPointMovement's
            change in its OperatingPoints' vCg__E parameters. Can be "sine", "uniform",
            or a callable custom spacing function. Custom spacing functions are for
            advanced users and must start at 0.0, return to 0.0 after one period of 2*pi
            radians, have amplitude of 1.0, be periodic, return finite values only, and
            accept a ndarray as input and return a ndarray of the same shape. The custom
            function is scaled by ampVCg__E, shifted horizontally and vertically by
            phaseVCg__E and the base value, and have a period set by periodVCg__E. The
            default is "sine".
        :param phaseVCg__E: The phase offset of the first time step's OperatingPoint's
            vCg__E parameter relative to the base OperatingPoint's vCg__E parameter.
            Must be a number (int or float) in the range (-180.0, 180.0], and will be
            converted to a float internally. It must be 0.0 if ampVCg__E is 0.0 and non
            zero if not. The units are in degrees. The default is 0.0.
        :return: None
        """
        # Validate and store immutable attributes.
        if not isinstance(base_operating_point, operating_point_mod.OperatingPoint):
            raise TypeError("base_operating_point must be an OperatingPoint")
        self._base_operating_point = base_operating_point

        self._ampVCg__E = _parameter_validation.number_in_range_return_float(
            ampVCg__E, "ampVCg__E", min_val=0.0, min_inclusive=True
        )

        periodVCg__E = _parameter_validation.number_in_range_return_float(
            periodVCg__E, "periodVCg__E", min_val=0.0, min_inclusive=True
        )
        if self._ampVCg__E == 0 and periodVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then periodVCg__E must also be 0.0.")
        self._periodVCg__E = periodVCg__E

        if isinstance(spacingVCg__E, str):
            if spacingVCg__E not in ["sine", "uniform"]:
                raise ValueError(
                    f"spacingVCg__E must be 'sine', 'uniform', or a callable, "
                    f"got string '{spacingVCg__E}'."
                )
        elif not callable(spacingVCg__E):
            raise TypeError(
                f"spacingVCg__E must be 'sine', 'uniform', or a callable, got "
                f"{type(spacingVCg__E).__name__}."
            )
        self._spacingVCg__E = spacingVCg__E

        phaseVCg__E = _parameter_validation.number_in_range_return_float(
            phaseVCg__E, "phaseVCg__E", -180.0, False, 180.0, True
        )
        if self._ampVCg__E == 0 and phaseVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then phaseVCg__E must also be 0.0.")
        self._phaseVCg__E = phaseVCg__E

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._max_period: float | None = None

    # --- Immutable: read only properties ---
    @property
    def base_operating_point(self) -> operating_point_mod.OperatingPoint:
        return self._base_operating_point

    @property
    def ampVCg__E(self) -> float:
        return self._ampVCg__E

    @property
    def periodVCg__E(self) -> float:
        return self._periodVCg__E

    @property
    def spacingVCg__E(self) -> str | Callable[[np.ndarray], np.ndarray]:
        return self._spacingVCg__E

    @property
    def phaseVCg__E(self) -> float:
        return self._phaseVCg__E

    # --- Immutable derived: manual lazy caching ---
    @property
    def max_period(self) -> float:
        """OperatingPointMovement's longest period of motion.

        :return: The longest period in seconds. If the motion is static, this will be
            0.0.
        """
        if self._max_period is None:
            self._max_period = self._periodVCg__E
        return self._max_period

    # --- Other methods ---
    def generate_operating_points(
        self, num_steps: int, delta_time: float | int
    ) -> list[operating_point_mod.OperatingPoint]:
        """Creates the OperatingPoint at each time step, and returns them in a list.

        :param num_steps: The number of time steps in this movement. It must be a
            positive int.
        :param delta_time: The time between each time step. It must be a positive number
            (int or float), and will be converted internally to a float. The units are
            in seconds.
        :return: The list of OperatingPoints associated with this
            OperatingPointMovement.
        """
        num_steps = _parameter_validation.int_in_range_return_int(
            num_steps,
            "num_steps",
            min_val=1,
            min_inclusive=True,
        )
        delta_time = _parameter_validation.number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )

        # Generate oscillating values for VCg__E.
        if self._spacingVCg__E == "sine":
            listVCg__E = _functions.oscillating_sinspaces(
                amps=self._ampVCg__E,
                periods=self._periodVCg__E,
                phases=self._phaseVCg__E,
                bases=self._base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self._spacingVCg__E == "uniform":
            listVCg__E = _functions.oscillating_linspaces(
                amps=self._ampVCg__E,
                periods=self._periodVCg__E,
                phases=self._phaseVCg__E,
                bases=self._base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif callable(self._spacingVCg__E):
            listVCg__E = _functions.oscillating_customspaces(
                amps=self._ampVCg__E,
                periods=self._periodVCg__E,
                phases=self._phaseVCg__E,
                bases=self._base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
                custom_function=self._spacingVCg__E,
            )
        else:
            raise ValueError(f"Invalid spacing value: {self._spacingVCg__E}")

        # Create an empty list to hold each time step's OperatingPoint.
        operating_points = []

        # Get the non changing OperatingPoint attributes.
        this_rho = self._base_operating_point.rho
        this_alpha = self._base_operating_point.alpha
        this_beta = self._base_operating_point.beta
        thisExternalFX_W = self._base_operating_point.externalFX_W
        this_nu = self._base_operating_point.nu
        thisAngles_E_to_BP1_izyx = self._base_operating_point.angles_E_to_BP1_izyx
        thisCgP1_E_Eo = self._base_operating_point.CgP1_E_Eo
        thisSurfaceNormal_E = self._base_operating_point.surfaceNormal_E
        thisSurfacePoint_E_Eo = self._base_operating_point.surfacePoint_E_Eo

        # Iterate through the time steps.
        for step in range(num_steps):
            thisVCg__E = listVCg__E[step]

            # Make a new operating point object for this time step.
            this_operating_point = operating_point_mod.OperatingPoint(
                rho=this_rho,
                vCg__E=thisVCg__E,
                alpha=this_alpha,
                beta=this_beta,
                externalFX_W=thisExternalFX_W,
                nu=this_nu,
                angles_E_to_BP1_izyx=thisAngles_E_to_BP1_izyx,
                CgP1_E_Eo=thisCgP1_E_Eo,
                surfaceNormal_E=thisSurfaceNormal_E,
                surfacePoint_E_Eo=thisSurfacePoint_E_Eo,
            )

            # Add this new OperatingPoint to the list of OperatingPoints.
            operating_points.append(this_operating_point)

        return operating_points
