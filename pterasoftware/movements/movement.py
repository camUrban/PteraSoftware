"""Contains the Movement class.

**Contains the following classes:**

Movement: A class used to contain an UnsteadyProblem's movement.

**Contains the following functions:**

None
"""

from __future__ import annotations

import math

from . import airplane_movement as airplane_movement_mod
from . import operating_point_movement as operating_point_movement_mod

from .. import _parameter_validation


class Movement:
    """A class used to contain an UnsteadyProblem's movement.

    **Contains the following methods:**

    max_period: The longest period of motion of Movement's sub movement objects, the
    motion(s) of its sub sub movement object(s), and the motions of its sub sub sub
    movement objects.

    static: Flags if Movement's sub movement objects, its sub sub movement object(s),
    and its sub sub sub movement objects all represent no motion.
    """

    def __init__(
        self,
        airplane_movements: list[airplane_movement_mod.AirplaneMovement],
        operating_point_movement: operating_point_movement_mod.OperatingPointMovement,
        delta_time: float | int | None = None,
        num_cycles: int | None = None,
        num_chords: int | None = None,
        num_steps: int | None = None,
    ) -> None:
        """The initialization method.

        This method checks that all Wings maintain their symmetry type across all time
        steps. See the WingMovement class documentation for more details on this
        requirement. See the Wing class documentation for more information on symmetry
        types.

        :param airplane_movements: A list of the AirplaneMovements associated with each
            of the UnsteadyProblem's Airplanes.
        :param operating_point_movement: An OperatingPointMovement characterizing any
            changes to the UnsteadyProblem's operating conditions.
        :param delta_time: The time between each time step. If set to None, Movement
            will calculate a value such that RingVortices shed from the first Wing will
            have roughly the same chord length as the RingVortices on the first Wing.
            This is based on first base Airplane's reference chord length, its first
            Wing's number of chordwise panels, and its base OperatingPoint's velocity.
            If not None, delta_time must be a positive number (int or float). It will be
            converted internally to a float. The units are in seconds. The default is
            None.
        :param num_cycles: The number of cycles of the maximum period motion used to
            calculate a num_steps parameter initialized as None if Movement isn't
            static. If num_steps is not None or if Movement is static, this must be
            None. If num_steps is initialized as None and the Movement isn't static,
            num_cycles must be a positive int. In that case, I recommend setting
            num_cycles to 3. The default is None.
        :param num_chords: The number of chord lengths used to calculate a num_steps
            parameter initialized as None if Movement is static. If num_steps is not
            None or if Movement isn't static, this must be None. If num_steps is
            initialized as None and Movement is static, num_chords must be a positive
            int. In that case, I recommend setting num_chords to 10. For cases with
            multiple Airplanes, the num_chords will reference the largest reference
            chord length. The default is None.
        :param num_steps: The number of time steps of the unsteady simulation. If
            initialized as None, and Movement isn't static, Movement will calculate a
            value for num_steps such that the simulation will cover some number of
            cycles of the maximum period of all the motion described in Movement's sub
            movement objects, sub sub movement object(s), and sub sub sub movement
            objects. If num_steps is initialized as None, and Movement is static,
            Movement will calculate a value for num_steps such that the simulation will
            result in a wake extending back by some number of reference chord lengths.
        :return: None
        """
        if not isinstance(airplane_movements, list):
            raise TypeError("airplane_movements must be a list.")
        if len(airplane_movements) < 1:
            raise ValueError("airplane_movements must have at least one element.")
        for airplane_movement in airplane_movements:
            if not isinstance(
                airplane_movement, airplane_movement_mod.AirplaneMovement
            ):
                raise TypeError(
                    "Every element in airplane_movements must be an AirplaneMovement."
                )
        self.airplane_movements = airplane_movements

        if not isinstance(
            operating_point_movement,
            operating_point_movement_mod.OperatingPointMovement,
        ):
            raise TypeError(
                "operating_point_movement must be an OperatingPointMovement."
            )
        self.operating_point_movement = operating_point_movement

        if delta_time is not None:
            delta_time = _parameter_validation.number_in_range_return_float(
                delta_time, "delta_time", min_val=0.0, min_inclusive=False
            )
        else:

            # FIXME: Automatic delta_time calculation gives very poor results if the
            #  motion has a high Strouhal number (i.e. a large ratio of
            #  flapping-motion to forward velocity). This is because the calculation
            #  assumes that the forward velocity is dominant. A better approach is
            #  needed.

            delta_times = []
            for airplane_movement in self.airplane_movements:
                # TODO: Consider making this also average across each Airplane's Wings.
                # For a given Airplane, the ideal time step length is that which
                # sheds RingVortices off the first Wing that have roughly the same
                # chord length as the RingVortices on the first Wing. This is based
                # on the base Airplane's reference chord length, its first Wing's
                # number of chordwise panels, and its base OperatingPoint's velocity.
                c_ref = airplane_movement.base_airplane.c_ref
                assert c_ref is not None
                delta_times.append(
                    c_ref
                    / airplane_movement.base_airplane.wings[0].num_chordwise_panels
                    / operating_point_movement.base_operating_point.vCg__E
                )

            # Set the delta_time to be the average of the Airplanes' ideal delta times.
            delta_time = sum(delta_times) / len(delta_times)
        self.delta_time: float = delta_time

        _static = self.static

        if num_steps is None:
            if _static:
                if num_cycles is not None:
                    raise ValueError(
                        "If num_steps is None and the Movement is static, num_cycles "
                        "must be left as None."
                    )
            else:
                if num_cycles is None:
                    raise ValueError(
                        "If num_steps is None and the Movement isn't static, "
                        "num_cycles must be set."
                    )
        if num_cycles is not None:
            num_cycles = _parameter_validation.int_in_range_return_int(
                num_cycles,
                "num_cycles",
                min_val=1,
                min_inclusive=True,
            )
        self.num_cycles = num_cycles

        if num_steps is None:
            if _static:
                if num_chords is None:
                    raise ValueError(
                        "If num_steps is None and the Movement is static, num_chords "
                        "must be set."
                    )
            else:
                if num_chords is not None:
                    raise ValueError(
                        "If num_steps is None and the Movement isn't static, "
                        "num_chords must be left as None."
                    )
        if num_chords is not None:
            num_chords = _parameter_validation.int_in_range_return_int(
                num_chords,
                "num_chords",
                min_val=1,
                min_inclusive=True,
            )
        self.num_chords = num_chords

        if self.num_cycles is not None or self.num_chords is not None:
            if num_steps is not None:
                raise ValueError(
                    "If either num_cycles or num_chords is not None, num_steps must "
                    "be None."
                )
        if num_steps is not None:
            num_steps = _parameter_validation.int_in_range_return_int(
                num_steps,
                "num_steps",
                min_val=1,
                min_inclusive=True,
            )
        else:
            if _static:
                # Find the value of the largest reference chord length of all the
                # base Airplanes.
                c_refs = []
                for airplane_movement in self.airplane_movements:
                    c_ref = airplane_movement.base_airplane.c_ref
                    assert c_ref is not None
                    c_refs.append(c_ref)
                max_c_ref = max(c_refs)

                # Set the number of time steps such that the wake extends back by
                # some number of reference chord lengths.
                assert self.num_chords is not None
                wake_length = self.num_chords * max_c_ref
                distance_per_time_step = (
                    delta_time
                    * self.operating_point_movement.base_operating_point.vCg__E
                )
                num_steps = math.ceil(wake_length / distance_per_time_step)
            else:
                # Set the number of time steps such that the simulation runs for some
                # number of cycles of the motion with the maximum period.
                assert self.num_cycles is not None
                num_steps = math.ceil(
                    self.num_cycles * self.max_period / self.delta_time
                )
        self.num_steps: int = num_steps

        # Generate a list of lists of Airplanes that are the steps through each
        # AirplaneMovement. The first index identifies the AirplaneMovement, and the
        # second index identifies the time step.
        self.airplanes = []
        for airplane_movement in self.airplane_movements:
            self.airplanes.append(
                airplane_movement.generate_airplanes(
                    num_steps=self.num_steps, delta_time=self.delta_time
                )
            )

        # Validate that all Wings maintain their symmetry type across all time steps.
        for airplane_movement_id, airplane_list in enumerate(self.airplanes):
            # Get the base Airplane (first time step).
            base_airplane = airplane_list[0]

            # Store the symmetry types of the base Wings.
            base_wing_symmetry_types = []
            for wing in base_airplane.wings:
                base_wing_symmetry_types.append(wing.symmetry_type)

            # Validate all subsequent time steps.
            for step_id, airplane in enumerate(airplane_list):
                # Check that Wings maintain their symmetry types.
                for wing_id, wing in enumerate(airplane.wings):
                    base_symmetry_type = base_wing_symmetry_types[wing_id]
                    if wing.symmetry_type != base_symmetry_type:
                        raise ValueError(
                            f"Wing {wing_id} in AirplaneMovement "
                            f"{airplane_movement_id} changed from type "
                            f"{base_symmetry_type} symmetry at time step 0 to type "
                            f"{wing.symmetry_type} symmetry at time step {step_id}. "
                            f"Wings cannot undergo motion that changes their symmetry "
                            f"type. This happens when a symmetric Wing moves such "
                            f"that its symmetry plane is no longer coincident with "
                            f"the wing axes' yz plane or vice versa."
                        )

        # Generate a lists of OperatingPoints that are the steps through the
        # OperatingPointMovement.
        self.operating_points = operating_point_movement.generate_operating_points(
            num_steps=self.num_steps, delta_time=self.delta_time
        )

    @property
    def max_period(self) -> float:
        """The longest period of motion of Movement's sub movement objects, the
        motion(s) of its sub sub movement object(s), and the motions of its sub sub sub
        movement objects.

        :return: The longest period in seconds. If all the motion is static, this will
            be 0.0.
        """
        # Iterate through the AirplaneMovements and find the one with the largest max
        # period.
        airplane_movement_max_periods = []
        for airplane_movement in self.airplane_movements:
            airplane_movement_max_periods.append(airplane_movement.max_period)
        max_airplane_period = max(airplane_movement_max_periods)

        # The global max period is the maximum of the max AirplaneMovement period and
        # the OperatingPointMovement max period.
        return max(
            max_airplane_period,
            self.operating_point_movement.max_period,
        )

    @property
    def static(self) -> bool:
        """Flags if the Movement's sub movement objects, its sub sub movement object(s),
        and its sub sub sub movement objects all represent no motion.

        :return: True if Movement's sub movement objects, its sub sub movement
            object(s), and its sub sub sub movement objects all represent no motion.
            False otherwise.
        """
        return self.max_period == 0
