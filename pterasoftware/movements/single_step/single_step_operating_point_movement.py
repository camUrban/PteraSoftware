"""This module contains the OperatingPointMovement class.

This module contains the following classes:
    OperatingPointMovement: This is a class used to contain the OperatingPoint
    movements.

This module contains the following functions:
    None
"""

from .._functions import (
    oscillating_sinspaces,
    oscillating_linspaces,
    oscillating_customspaces,
)

from ...operating_point import OperatingPoint
from ..._parameter_validation import (
    number_in_range_return_float,
    int_in_range_return_int,
)


class SingleStepOperatingPointMovement:
    """This is a class used to contain the OperatingPoint movements.

    This class contains the following public methods:

        generate_operating_points: Creates the OperatingPoint at each time step,
        and returns them in a list.

        max_period: Defines a property for the longest period of
        OperatingPointMovement's own motion.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    ):
        """This is the initialization method.

        :param base_operating_point: OperatingPoint

            This is the base OperatingPoint, from which the OperatingPoint at each
            time step will be created.

        :param ampVCg__E: number, optional

            The amplitude of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameters. Must be a non-negative number (int or
            float), and is converted to a float internally.  Also, the amplitude must
            be low enough that it doesn't drive its base value out of the range of
            valid values. Otherwise, this OperatingPointMovement will try to create
            OperatingPoints with invalid parameters values.The default value is 0.0.
            The units are in meters per second.

        :param periodVCg__E: number, optional

            The period of the OperatingPointMovement's changes in its
            OperatingPoints' vCg__E parameter. Must be a non-negative number (int or
            float), and is converted to a float internally. The default value is 0.0.
            It must be 0.0 if ampVCg__E 0.0 and non-zero if not. The units are in
            seconds.

        :param spacingVCg__E: string, optional

            The value determines the spacing of the OperatingPointMovement's change
            in its OperatingPoints' vCg__E parameters. Must be either "sine",
            "uniform", or a callable custom spacing function. Custom spacing
            functions are for advanced users and must start at 0, return to 0 after
            one period of 2*pi radians, have amplitude of 1, be periodic,
            return finite values only, and accept a ndarray as input and return a
            ndarray of the same shape. The custom function is scaled by ampVCg__E,
            shifted horizontally by phaseVCg__E, and vertically by the base value,
            with the period controlled by periodVCg__E. The default value is "sine".

        :param phaseVCg__E: number optional

            The phase offsets of the first time step's OperatingPoint's vCg__E
            parameter relative to the base OperatingPoint's vCg__E parameter. Must be
            a number (int or float) in the range (-180.0, 180.0], and is converted to a
            float internally. The default value is 0.0. It must be 0.0 if ampVCg__E
            is 0.0 and non-zero if not. The units are in degrees.
        """

        self.ampVCg__E = number_in_range_return_float(
            ampVCg__E, "ampVCg__E", min_val=0.0, min_inclusive=True
        )

        periodVCg__E = number_in_range_return_float(
            periodVCg__E, "periodVCg__E", min_val=0.0, min_inclusive=True
        )
        if self.ampVCg__E == 0 and periodVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then periodVCg__E must also be 0.0.")
        self.periodVCg__E = periodVCg__E

        if isinstance(spacingVCg__E, str):
            if spacingVCg__E not in ["sine", "uniform"]:
                raise ValueError(
                    f"spacingVCg__E must be 'sine', 'uniform', or a callable, got string '{spacingVCg__E}'."
                )
        elif not callable(spacingVCg__E):
            raise TypeError(
                f"spacingVCg__E must be 'sine', 'uniform', or a callable, got {type(spacingVCg__E).__name__}."
            )
        self.spacingVCg__E = spacingVCg__E

        phaseVCg__E = number_in_range_return_float(
            phaseVCg__E, "phaseVCg__E", -180.0, False, 180.0, True
        )
        if self.ampVCg__E == 0 and phaseVCg__E != 0:
            raise ValueError("If ampVCg__E is 0.0, then phaseVCg__E must also be 0.0.")
        self.phaseVCg__E = phaseVCg__E

        self.listVCg__E = None

    def generate_next_operating_point(self, delta_time, base_operating_point: OperatingPoint, num_steps, step):
        """Creates the OperatingPoint at each time step, and returns them in a list.

        :param num_steps: int

            This is the number of time steps in this movement. It must be a positive
            int.

        :param delta_time: number

            This is the time between each time step. It must be a positive number (
            int or float), and will be converted internally to a float. The units are
            in seconds.

        :return: list of OperatingPoints

            This is the list of OperatingPoints associated with this
            OperatingPointMovement.
        """
        num_steps = int_in_range_return_int(
            num_steps, "num_steps", min_val=1, min_inclusive=True
        )
        delta_time = number_in_range_return_float(
            delta_time, "delta_time", min_val=0.0, min_inclusive=False
        )

        if self.listVCg__E is None:
            self._initialize_oscillating_values(
                delta_time=delta_time,
                num_steps=num_steps,
                base_operating_point=base_operating_point,
            )

        # Get the non-changing OperatingPoint attributes.
        this_rho = base_operating_point.rho
        this_alpha = base_operating_point.alpha
        this_beta = base_operating_point.beta
        thisExternalFX_W = base_operating_point.externalFX_W
        this_nu = base_operating_point.nu

        # Make a new operating point object for this time step.
        this_operating_point = OperatingPoint(
            rho=this_rho,
            vCg__E=self.listVCg__E[step],
            alpha=this_alpha,
            beta=this_beta,
            externalFX_W=thisExternalFX_W,
            nu=this_nu,
        )

        return this_operating_point

    @property
    def max_period(self):
        """Defines a property for the longest period of OperatingPointMovement's own
        motion.

        :return: float

            The longest period in seconds. If the all the motion is static, this will
            be 0.0.
        """
        return self.periodVCg__E

    def _initialize_oscillating_values(self, delta_time, num_steps, base_operating_point):
        """Pre-computes the oscillating values for faster access later.

        :param delta_time: number  

            This is the time between each time step. It must be a positive number (
            int or float), and will be converted internally to a float. The units are
            in seconds.
        :param num_steps: int
            This is the number of time steps in this movement. It must be a positive
            int.
        :param base_operating_point: OperatingPoint
            This is the base OperatingPoint, from which the OperatingPoint at each
            time step will be created.
        """
        # Generate oscillating values for VCg__E.
        if self.spacingVCg__E == "sine":
            self.listVCg__E = oscillating_sinspaces(
                amps=self.ampVCg__E,
                periods=self.periodVCg__E,
                phases=self.phaseVCg__E,
                bases=base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.spacingVCg__E == "uniform":
            self.listVCg__E = oscillating_linspaces(
                amps=self.ampVCg__E,
                periods=self.periodVCg__E,
                phases=self.phaseVCg__E,
                bases=base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif callable(self.spacingVCg__E):
            self.listVCg__E = oscillating_customspaces(
                amps=self.ampVCg__E,
                periods=self.periodVCg__E,
                phases=self.phaseVCg__E,
                bases=base_operating_point.vCg__E,
                num_steps=num_steps,
                delta_time=delta_time,
                custom_function=self.spacingVCg__E,
            )
        else:
            raise ValueError(f"Invalid spacing value: {self.spacingVCg__E}")
