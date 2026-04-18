"""Contains the AeroelasticOperatingPointMovement class.

**Contains the following classes:**

AeroelasticOperatingPointMovement: A class used to contain an OperatingPoint's movements
in an aeroelastic simulation.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Callable

from .. import _core
from .. import operating_point as operating_point_mod


class AeroelasticOperatingPointMovement(_core.CoreOperatingPointMovement):
    """A class used to contain an OperatingPoint's movements in an aeroelastic
    simulation.

    In aeroelastic simulations, OperatingPoints are prescribed via oscillation
    parameters (the same oscillation based generation as OperatingPointMovement). This
    class exists so that an AeroelasticMovement always accepts
    AeroelasticOperatingPointMovements, keeping the aeroelastic movement hierarchy
    consistently named.

    **Contains the following methods:**

    max_period: AeroelasticOperatingPointMovement's longest period of motion.

    generate_operating_point_at_time_step: Creates the OperatingPoint at a single time
    step.

    generate_operating_points: Creates the OperatingPoint at each time step, and returns
    them in a list.
    """

    __slots__ = ()

    def __init__(
        self,
        base_operating_point: operating_point_mod.OperatingPoint,
        ampVCg__E: float | int = 0.0,
        periodVCg__E: float | int = 0.0,
        spacingVCg__E: str | Callable[[float], float] = "sine",
        phaseVCg__E: float | int = 0.0,
    ) -> None:
        """The initialization method.

        :param base_operating_point: The base OperatingPoint from which the
            OperatingPoint at each time step will be created.
        :param ampVCg__E: The amplitude of the AeroelasticOperatingPointMovement's
            change in its OperatingPoints' vCg__E parameters. It must be a non negative
            number (int or float). Values are converted to floats internally. The
            amplitude must be low enough that it doesn't drive its base value out of the
            range of valid values. Otherwise, this AeroelasticOperatingPointMovement
            will try to create OperatingPoints with invalid parameter values. The units
            are in meters per second. The default is 0.0.
        :param periodVCg__E: The period of the AeroelasticOperatingPointMovement's
            change in its OperatingPoints' vCg__E parameters. It must be a non negative
            number (int or float). Values are converted to floats internally. It must be
            0.0 if ampVCg__E is 0.0, and non zero if not. The units are in seconds. The
            default is 0.0.
        :param spacingVCg__E: The spacing type of the
            AeroelasticOperatingPointMovement's change in its OperatingPoints' vCg__E
            parameters. It can be the str "sine", the str "uniform", or a callable
            custom spacing function. Custom spacing functions are for advanced users and
            must start at 0.0, return to 0.0 after one period of 2.0 * pi radians, have
            amplitude of 1.0, be periodic, return finite values only, and accept a float
            as input and return a float. Custom functions are scaled by ampVCg__E,
            shifted horizontally and vertically by phaseVCg__E and the base value, and
            have a period set by periodVCg__E. The default is "sine".
        :param phaseVCg__E: The phase offset of the first time step's OperatingPoint's
            vCg__E parameter relative to the base OperatingPoint's vCg__E parameter. It
            must be a number (int or float) in the range (-180.0, 180.0]. It must be 0.0
            if ampVCg__E is 0.0 and non zero if not. Values are converted to floats
            internally. The units are in degrees. The default is 0.0.
        :return: None
        """
        super().__init__(
            base_operating_point=base_operating_point,
            ampVCg__E=ampVCg__E,
            periodVCg__E=periodVCg__E,
            spacingVCg__E=spacingVCg__E,
            phaseVCg__E=phaseVCg__E,
        )
