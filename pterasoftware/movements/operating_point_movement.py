"""Contains the OperatingPointMovement class.

**Contains the following classes:**

OperatingPointMovement: A class used to contain an OperatingPoint's movements.

**Contains the following functions:**

None
"""

from collections.abc import Callable

from . import _functions

from .. import operating_point as operating_point_mod
from .. import _parameter_validation


class OperatingPointMovement:
    """A class used to contain an OperatingPoint's movements.

    **Contains the following methods:**

    generate_operating_points: Creates the OperatingPoint at each time step, and returns
    them in a list.

    max_period: Defines a property for OperatingPointMovement's longest period of
    motion.
    """

    def __init__(
        self,
        base_operating_point: operating_point_mod.OperatingPoint,
        ampVCg__E: float | int = 0.0,
        periodVCg__E: float | int = 0.0,
        spacingVCg__E: str | Callable = "sine",
        phaseVCg__E: float | int = 0.0,
    ) -> None:
        """The initialization method.

        :param base_operating_point: The base OperatingPoint from which the
            OperatingPoint at each time step will be created.
        :param ampVCg__E: The amplitude of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameters. Must be a non-negative number (int or
            float), and is converted to a float internally. The amplitude must be low
            enough that it doesn't drive its base value out of the range of valid
            values. Otherwise, this OperatingPointMovement will try to create
            OperatingPoints with invalid parameter values. The units are in meters per
            second. The default is 0.0.
        :param periodVCg__E: The period of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameter. Must be a non-negative number (int or
            float), and is converted to a float internally. It must be 0.0 if ampVCg__E
            is 0.0 and non-zero if not. The units are in seconds. The default is 0.0.
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
            converted to a float internally. It must be 0.0 if ampVCg__E is 0.0 and non-
            zero if not. The units are in degrees. The default is 0.0.
        :return: None
        """
        if not isinstance(base_operating_point, operating_point_mod.OperatingPoint):
            raise TypeError("base_operating_point must be an OperatingPoint")
        self.base_operating_point = base_operating_point

        self.ampVCg__E = _parameter_validation.non_negative_number_return_float(
            ampVCg__E, "ampVCg__E"
        )

        periodVCg__E = _parameter_validation.non_negative_number_return_float(
            periodVCg__E, "periodVCg__E"
        )
        if self.ampVCg__E == 0 and periodVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then periodVCg__E must also be 0.0.")
        self.periodVCg__E = periodVCg__E

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
        self.spacingVCg__E = spacingVCg__E

        phaseVCg__E = _parameter_validation.number_in_range_return_float(
            phaseVCg__E, "phaseVCg__E", -180.0, False, 180.0, True
        )
        if self.ampVCg__E == 0 and phaseVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then phaseVCg__E must also be 0.0.")
        self.phaseVCg__E = phaseVCg__E

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
        num_steps = _parameter_validation.positive_int_return_int(
            num_steps, "num_steps"
        )
        delta_time = _parameter_validation.positive_number_return_float(
            delta_time, "delta_time"
        )

        # Generate oscillating values for VCg__E.
        if self.spacingVCg__E == "sine":
            listVCg__E = _functions.oscillating_sinspaces(
                amps=self.ampVCg__E,
                periods=self.periodVCg__E,
                phases=self.phaseVCg__E,
                bases=self.base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.spacingVCg__E == "uniform":
            listVCg__E = _functions.oscillating_linspaces(
                amps=self.ampVCg__E,
                periods=self.periodVCg__E,
                phases=self.phaseVCg__E,
                bases=self.base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif callable(self.spacingVCg__E):
            listVCg__E = _functions.oscillating_customspaces(
                amps=self.ampVCg__E,
                periods=self.periodVCg__E,
                phases=self.phaseVCg__E,
                bases=self.base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
                custom_function=self.spacingVCg__E,
            )
        else:
            raise ValueError(f"Invalid spacing value: {self.spacingVCg__E}")

        # Create an empty list to hold each time step's OperatingPoint.
        operating_points = []

        # Get the non-changing OperatingPoint attributes.
        this_rho = self.base_operating_point.rho
        this_alpha = self.base_operating_point.alpha
        this_beta = self.base_operating_point.beta
        thisExternalFX_W = self.base_operating_point.externalFX_W
        this_nu = self.base_operating_point.nu

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
            )

            # Add this new OperatingPoint to the list of OperatingPoints.
            operating_points.append(this_operating_point)

        return operating_points

    @property
    def max_period(self) -> float:
        """Defines a property for OperatingPointMovement's longest period of motion.

        :return: The longest period in seconds. If the motion is static, this will be
            0.0.
        """
        return self.periodVCg__E
