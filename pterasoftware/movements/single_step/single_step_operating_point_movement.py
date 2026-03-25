"""Contains the SingleStepOperatingPointMovement class.

**Contains the following classes:**

SingleStepOperatingPointMovement: A single step variant of OperatingPointMovement
that generates one OperatingPoint per time step instead of all at once. Uses
composition to wrap an OperatingPointMovement.

**Contains the following functions:**

None
"""
from ..operating_point_movement import OperatingPointMovement
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
    """A single step variant of OperatingPointMovement for coupled simulations.

    This class wraps an OperatingPointMovement via composition and generates one
    OperatingPoint per time step (via generate_next_operating_point) rather than
    generating all OperatingPoints at once. The composed OperatingPointMovement is
    accessible via corresponding_operating_point_movement.

    **Contains the following methods:**

    generate_next_operating_point: Creates the OperatingPoint at a single time step.

    max_period: The longest period of this SingleStepOperatingPointMovement's own
    motion.
    """

    __slots__ = (
        "ampVCg__E",
        "periodVCg__E",
        "spacingVCg__E",
        "phaseVCg__E",
        "listVCg__E",
        "corresponding_operating_point_movement",
    )

    def __init__(
        self,
        base_operating_point: OperatingPoint,
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    ):
        """The initialization method.

        :param base_operating_point: The base OperatingPoint from which the
            OperatingPoint at each time step will be created.
        :param ampVCg__E: The amplitude of the SingleStepOperatingPointMovement's
            changes in its OperatingPoints' vCg__E parameters. Must be a non negative
            number (int or float), and is converted to a float internally. The
            amplitude must be low enough that it doesn't drive its base value out of
            the range of valid values. Otherwise, this
            SingleStepOperatingPointMovement will try to create OperatingPoints with
            invalid parameter values. The units are in meters per second. The default
            is 0.0.
        :param periodVCg__E: The period of the SingleStepOperatingPointMovement's
            changes in its OperatingPoints' vCg__E parameter. Must be a non negative
            number (int or float), and is converted to a float internally. It must be
            0.0 if ampVCg__E is 0.0 and non zero if not. The units are in seconds. The
            default is 0.0.
        :param spacingVCg__E: Determines the spacing of the
            SingleStepOperatingPointMovement's change in its OperatingPoints' vCg__E
            parameters. Can be "sine", "uniform", or a callable custom spacing
            function. Custom spacing functions are for advanced users and must start at
            0.0, return to 0.0 after one period of 2*pi radians, have amplitude of
            1.0, be periodic, return finite values only, and accept a ndarray as input
            and return a ndarray of the same shape. The custom function is scaled by
            ampVCg__E, shifted horizontally and vertically by phaseVCg__E and the base
            value, and have a period set by periodVCg__E. The default is "sine".
        :param phaseVCg__E: The phase offset of the first time step's OperatingPoint's
            vCg__E parameter relative to the base OperatingPoint's vCg__E parameter.
            Must be a number (int or float) in the range (-180.0, 180.0], and will be
            converted to a float internally. It must be 0.0 if ampVCg__E is 0.0 and
            non zero if not. The units are in degrees. The default is 0.0.
        :return: None
        """

        # Create the corresponding OperatingPointMovement, which validates all
        # oscillation parameters and is also needed by coupled unsteady problems.
        self.corresponding_operating_point_movement = OperatingPointMovement(
            base_operating_point=base_operating_point,
            ampVCg__E=ampVCg__E,
            periodVCg__E=periodVCg__E,
            spacingVCg__E=spacingVCg__E,
            phaseVCg__E=phaseVCg__E,
        )

        # Copy validated attributes from the corresponding OperatingPointMovement.
        self.ampVCg__E = self.corresponding_operating_point_movement.ampVCg__E
        self.periodVCg__E = self.corresponding_operating_point_movement.periodVCg__E
        self.spacingVCg__E = self.corresponding_operating_point_movement.spacingVCg__E
        self.phaseVCg__E = self.corresponding_operating_point_movement.phaseVCg__E

        self.listVCg__E = None

    def generate_next_operating_point(self, delta_time, base_operating_point: OperatingPoint, num_steps, step):
        """Creates the OperatingPoint at a single time step.

        :param delta_time: The time between each time step in seconds.
        :param base_operating_point: The base OperatingPoint from which the new
            OperatingPoint will be created.
        :param num_steps: The total number of time steps in this movement.
        :param step: The index of the current time step.
        :return: The OperatingPoint at the specified time step.
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
        """The longest period of this SingleStepOperatingPointMovement's own motion.

        :return: The longest period in seconds. If all the motion is static, this
            will be 0.0.
        """
        return self.periodVCg__E

    def _initialize_oscillating_values(self, delta_time, num_steps, base_operating_point):
        """Pre computes the oscillating VCg__E values for all time steps.

        :param delta_time: The time between each time step in seconds.
        :param num_steps: The total number of time steps.
        :param base_operating_point: The base OperatingPoint providing the base
            VCg__E value.
        :return: None
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
