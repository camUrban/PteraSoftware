"""This module contains the Movement class.

This module contains the following classes:
    Movement: This is a class used to contain an UnsteadyProblem's movement.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import math

from .airplane_movement import AirplaneMovement
from .operating_point_movement import OperatingPointMovement

from .. import parameter_validation


# TODO: Add unit tests for this Class.
class Movement:
    """This is a class used to contain an UnsteadyProblem's movement.

    This class contains the following public methods:

        max_period: Defines a property for the longest period of Movement's own
        motion and that of its sub-movement objects, sub-sub-movement objects, etc.

        static: Defines a property to flag if all the Movement itself, and all of its
        sub-movement objects, sub-sub-movement objects, etc. represent no motion.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        airplane_movements,
        operating_point_movement,
        delta_time=None,
        num_cycles=None,
        num_chords=None,
        num_steps=None,
    ):
        """This is the initialization method.

        :param airplane_movements: list of AirplaneMovements

            This is a list of objects which characterize the movement of each
            of the airplanes in the UnsteadyProblem.

        :param operating_point_movement: OperatingPointMovement

            This object characterizes changes to the UnsteadyProblem's the operating
            point.

        :param delta_time: number or None, optional

            delta_time is the time, in seconds, between each time step. If left as
            None, which is the default value, Movement will calculate a value such
            that RingVortices shed from the first Wing will have roughly the same
            chord length as the RingVortices on the first Wing. This is based on
            first base Airplane's reference chord length, its first Wing's number of
            chordwise panels, and its base OperatingPoint's velocity. If set,
            delta_time must be a positive number (int or float). It will be converted
            internally to a float.

        :param num_cycles: int or None, optional

            num_cycles is the number of cycles of the maximum period motion used to
            calculate a non-populated num_steps parameter if Movement isn't static.
            If num_steps is set or Movement is static, this must be left as None,
            which is the default value. If num_steps isn't set and Movement isn't
            static, num_cycles must be a positive int. In that case, I recommend
            setting num_cycles to 3.

        :param num_chords: int or None, optional

            num_chords is the number of chord lengths used to calculate a
            non-populated num_steps parameter if Movement is static. If num_steps is
            set or Movement isn't static, this must be left as None, which is the
            default value. If num_steps isn't set and Movement is static, num_chords
            must be a positive int. In that case, I recommend setting num_chords to
            10. For cases with multiple Airplanes, the num_chords will reference the
            largest reference chord length.

        :param num_steps: int or None, optional

            num_steps is the number of time steps of the unsteady simulation. It must
            be a positive int. The default value is None. If left as None,
            and Movement isn't static, Movement will calculate a value such that the
            simulation will cover some number of cycles of the maximum period of all
            the motion described in Movement's sub-movement objects, sub-sub-movement
            objects, etc. If num_steps is left as None, and Movement is static,
            it will default to the number of time steps such that the wake extends
            back by some number of reference chord lengths.
        """
        if not isinstance(airplane_movements, list):
            raise TypeError("airplane_movements must be a list.")
        if len(airplane_movements) < 1:
            raise ValueError("airplane_movements must have at least one element.")
        for airplane_movement in airplane_movements:
            if not isinstance(airplane_movement, AirplaneMovement):
                raise TypeError(
                    "Every element in airplane_movements must be an AirplaneMovement."
                )
        self.airplane_movements = airplane_movements

        if not isinstance(operating_point_movement, OperatingPointMovement):
            raise TypeError(
                "operating_point_movement must be an OperatingPointMovement."
            )
        self.operating_point_movement = operating_point_movement

        if delta_time is not None:
            delta_time = parameter_validation.positive_number_return_float(
                delta_time, "delta_time"
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
                delta_times.append(
                    airplane_movement.base_airplane.c_ref
                    / airplane_movement.base_airplane.wings[0].num_chordwise_panels
                    / operating_point_movement.base_operating_point.vCg__E
                )

            # Set the delta_time to be the average of the Airplanes' ideal delta times.
            delta_time = sum(delta_times) / len(delta_times)
        self.delta_time = delta_time

        _static = self.static

        if num_steps is None:
            if _static:
                if num_cycles is not None:
                    raise ValueError(
                        "If num_steps is None and the Movement is static, num_cycles must be left as None."
                    )
            else:
                if num_cycles is None:
                    raise ValueError(
                        "If num_steps is None and the Movement isn't static, num_cycles must be set."
                    )
        if num_cycles is not None:
            num_cycles = parameter_validation.positive_int_return_int(
                num_cycles, "num_cycles"
            )
        self.num_cycles = num_cycles

        if num_steps is None:
            if _static:
                if num_chords is None:
                    raise ValueError(
                        "If num_steps is None and the Movement is static, num_chords must be set."
                    )
            else:
                if num_chords is not None:
                    raise ValueError(
                        "If num_steps is None and the Movement isn't static, num_chords must be left as None."
                    )
        if num_chords is not None:
            num_chords = parameter_validation.positive_int_return_int(
                num_chords, "num_chords"
            )
        self.num_chords = num_chords

        if self.num_cycles is not None or self.num_chords is not None:
            if num_steps is not None:
                raise ValueError(
                    "If either num_cycles or num_chords is not None, num_steps must be None."
                )
        if num_steps is not None:
            num_steps = parameter_validation.positive_int_return_int(
                num_steps, "num_steps"
            )
        else:
            if _static:
                # Find the value of the largest reference chord length of all the
                # base Airplanes.
                c_refs = []
                for airplane_movement in self.airplane_movements:
                    c_refs.append(airplane_movement.base_airplane.c_ref)
                max_c_ref = max(c_refs)

                # Set the number of time steps such that the wake extends back by
                # some number of reference chord lengths.
                wake_length = self.num_chords * max_c_ref
                distance_per_time_step = (
                    delta_time
                    * self.operating_point_movement.base_operating_point.vCg__E
                )
                num_steps = math.ceil(wake_length / distance_per_time_step)
            else:
                # Set the number of time steps such that the simulation runs for some
                # number of cycles of the motion with the maximum period.
                num_steps = math.ceil(
                    self.num_cycles * self.max_period / self.delta_time
                )
        self.num_steps = num_steps

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

        # Generate a lists of OperatingPoints that are the steps through the
        # OperatingPointMovement.
        self.operating_points = operating_point_movement.generate_operating_points(
            num_steps=self.num_steps, delta_time=self.delta_time
        )

    # TODO: Add unit tests for this method.
    @property
    def max_period(self):
        """Defines a property for the longest period of Movement's own motion and
        that of its sub-movement objects, sub-sub-movement objects, etc.

        :return: float The longest period in seconds. If the all the motion is
        static, this will be 0.0.
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

    # TODO: Add unit tests for this method.
    @property
    def static(self):
        """Defines a property to flag if all the Movement itself, and all of its
        sub-movement objects, sub-sub-movement objects, etc. represent no motion.

        :return: bool

            True if all Movement and its sub-movement objects, sub-sub-movement
            objects, etc. represent no motion. False otherwise.
        """
        return self.max_period == 0
