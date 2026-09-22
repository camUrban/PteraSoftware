"""Contains the FreeFlightOperatingPointMovement class."""

from __future__ import annotations

from .. import _core
from .. import operating_point as operating_point_mod


class FreeFlightOperatingPointMovement(_core.CoreOperatingPointMovement):
    """A class used to contain an OperatingPoint's movements in a free flight
    simulation.

    In free flight, OperatingPoints are not prescribed via oscillation parameters. They
    are dynamically determined by the solver as it integrates rigid body dynamics at
    each time step. FreeFlightOperatingPointMovement holds the initial OperatingPoint
    and provides a mutable list that the solver populates as dynamics integration
    produces new states.
    """

    __slots__ = ("operating_points",)

    def __init__(
        self,
        base_operating_point: operating_point_mod.OperatingPoint,
    ) -> None:
        """The initialization method.

        :param base_operating_point: The initial OperatingPoint representing the
            operating conditions at the start of the simulation.
        :return: None
        """
        super().__init__(base_operating_point=base_operating_point)

        # Mutable list of OperatingPoints. The solver appends new OperatingPoints from
        # dynamics integration at each time step. Starts with the base OperatingPoint at
        # step 0.
        self.operating_points: list[operating_point_mod.OperatingPoint] = [
            base_operating_point
        ]
