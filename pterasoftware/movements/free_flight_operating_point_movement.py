"""Contains the FreeFlightOperatingPointMovement class.

**Contains the following classes:**

FreeFlightOperatingPointMovement: A class used to contain an OperatingPoint's movements
in a free flight simulation.

**Contains the following functions:**

None
"""

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

    **Contains the following methods:**

    max_period: FreeFlightOperatingPointMovement's longest period of motion.

    generate_operating_point_at_time_step: Creates the OperatingPoint at a single time
    step.

    generate_operating_points: Creates the OperatingPoint at each time step, and returns
    them in a list.
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

        # Mutable list of OperatingPoints. The solver appends new
        # OperatingPoints from dynamics integration at each time step. Starts
        # with the base OperatingPoint at step 0.
        self.operating_points: list[operating_point_mod.OperatingPoint] = [
            base_operating_point
        ]
