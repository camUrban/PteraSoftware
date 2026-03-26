"""Contains the SingleStepMovement class.

**Contains the following classes:**

SingleStepMovement: A single step variant of Movement that generates Airplanes and
OperatingPoints one time step at a time instead of all at once. Uses composition to wrap
a Movement.

**Contains the following functions:**

None
"""

from ..._parameter_validation import (
    int_in_range_return_int,
    number_in_range_return_float,
)
from ..movement import Movement
from .single_step_airplane_movement import SingleStepAirplaneMovement
from .single_step_operating_point_movement import SingleStepOperatingPointMovement


class SingleStepMovement:
    """A single step variant of Movement for coupled simulations.

    This class wraps a Movement via composition and generates Airplanes and
    OperatingPoints one time step at a time (via generate_next_movement) rather than
    generating all of them at once. The composed Movement is accessible via
    corresponding_movement.

    **Contains the following methods:**

    generate_next_movement: Creates the Airplanes and OperatingPoint at a single time
    step.
    """

    __slots__ = (
        "airplane_movements",
        "operating_point_movement",
        "corresponding_movement",
        "delta_time",
        "num_steps",
    )

    def __init__(
        self,
        single_step_airplane_movements,
        single_step_operating_point_movement,
        delta_time: float | int | str | None = None,
        num_cycles: int | None = None,
        num_chords: int | None = None,
        num_steps: int | None = None,
    ):
        """The initialization method.

        :param single_step_airplane_movements: A list of SingleStepAirplaneMovements
            characterizing the movement of each Airplane.
        :param single_step_operating_point_movement: A SingleStepOperatingPointMovement
            characterizing any changes to the operating conditions.
        :param delta_time: The time between each time step. Accepts the following: None
            (default): SingleStepMovement analytically estimates the delta_time that
            produces wake RingVortices with roughly the same chord length as the bound
            trailing edge RingVortices, accounting for both freestream and geometry
            motion velocities. This provides good results across all Strouhal numbers.
            "optimize": SingleStepMovement first runs the analytical estimation, then
            uses that result as an initial guess for an iterative optimization that
            minimizes the area mismatch between wake RingVortices and their parent bound
            trailing edge RingVortices. This is slower but may produce slightly more
            accurate results. Positive number (int or float): Use the specified value
            directly. All values are converted internally to floats. The units are in
            seconds.
        :param num_cycles: The number of cycles of the maximum period motion used to
            calculate a num_steps parameter initialized as None if the
            SingleStepMovement isn't static. If num_steps is not None or if the
            SingleStepMovement is static, this must be None. If num_steps is initialized
            as None and the SingleStepMovement isn't static, num_cycles must be a
            positive int. In that case, I recommend setting num_cycles to 3. The default
            is None.
        :param num_chords: The number of chord lengths used to calculate a num_steps
            parameter initialized as None if the SingleStepMovement is static. If
            num_steps is not None or if the SingleStepMovement isn't static, this must
            be None. If num_steps is initialized as None and the SingleStepMovement is
            static, num_chords must be a positive int. In that case, I recommend setting
            num_chords to 10. For cases with multiple Airplanes, the num_chords will
            reference the largest reference chord length. The default is None.
        :param num_steps: The number of time steps of the unsteady simulation. If
            initialized as None, and the SingleStepMovement isn't static, it will
            calculate a value for num_steps such that the simulation will cover some
            number of cycles of the maximum period of all the motion described in the
            SingleStepMovement's sub movement objects, sub sub movement object(s), and
            sub sub sub movement objects. If num_steps is initialized as None, and the
            SingleStepMovement is static, it will calculate a value for num_steps such
            that the simulation will result in a wake extending back by some number of
            reference chord lengths.
        :return: None
        """
        if not isinstance(single_step_airplane_movements, list):
            raise TypeError("single_step_airplane_movements must be a list.")
        if len(single_step_airplane_movements) < 1:
            raise ValueError(
                "single_step_airplane_movements must have at least one element."
            )
        for airplane_movement in single_step_airplane_movements:
            if not isinstance(airplane_movement, SingleStepAirplaneMovement):
                raise TypeError(
                    "Every element in single_step_airplane_movements must be an SingleStepAirplaneMovement."
                )
        self.airplane_movements = single_step_airplane_movements

        if not isinstance(
            single_step_operating_point_movement, SingleStepOperatingPointMovement
        ):
            raise TypeError(
                "single_step_operating_point_movement must be an SingleStepOperatingPointMovement."
            )
        self.operating_point_movement = single_step_operating_point_movement

        corresponding_airplane_movements = [
            airplane_movement.corresponding_airplane_movement
            for airplane_movement in self.airplane_movements
        ]
        self.corresponding_movement = Movement(
            airplane_movements=corresponding_airplane_movements,
            operating_point_movement=self.operating_point_movement.corresponding_operating_point_movement,
            delta_time=delta_time,
            num_chords=num_chords,
            num_cycles=num_cycles,
            num_steps=num_steps,
        )

        self.delta_time = self.corresponding_movement.delta_time
        self.num_steps = self.corresponding_movement.num_steps

    def generate_next_movement(
        self, base_airplanes, base_operating_point, step, deformation_matrices=None
    ):
        """Creates the Airplanes and OperatingPoint at a single time step.

        :param base_airplanes: The list of Airplanes at the current time step.
        :param base_operating_point: The OperatingPoint at the current time step.
        :param step: The index of the time step to generate.
        :param deformation_matrices: Deformation matrices to apply to the Wings, or
            None. The default is None.
        :return: A tuple of (list of Airplanes, OperatingPoint) at the specified time
            step.
        """

        airplanes = []
        airplane_movement: SingleStepAirplaneMovement
        for airplane_id, airplane_movement in enumerate(self.airplane_movements):
            airplanes.append(
                airplane_movement.generate_next_airplane(
                    delta_time=self.delta_time,
                    base_airplane=base_airplanes[airplane_id],
                    num_steps=self.num_steps,
                    step=step,
                    deformation_matrices=deformation_matrices,
                )
            )

        operating_point = self.operating_point_movement.generate_next_operating_point(
            delta_time=self.delta_time,
            base_operating_point=base_operating_point,
            num_steps=self.num_steps,
            step=step,
        )
        return airplanes, operating_point
