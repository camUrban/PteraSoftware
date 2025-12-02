"""This module contains the Movement class.

This module contains the following classes:
    Movement: This is a class used to contain an UnsteadyProblem's movement.

This module contains the following functions:
    None
"""

import math

from .single_step_airplane_movement import SingleStepAirplaneMovement 
from .single_step_operating_point_movement import SingleStepOperatingPointMovement

from ..._parameter_validation import (
    positive_number_return_float,
    positive_int_return_int,
)

class SingleStepMovement:
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
        single_step_airplane_movements,
        single_step_operating_point_movement,
        delta_time=None,
        num_steps=40,
    ):
        """This is the initialization method.

        Note: This method checks that all Wings maintain their symmetry type across
        all time steps. See the WingMovement class documentation for more details on
        this requirement, and the Wing class documentation for more information on
        symmetry types.

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

        if delta_time is not None:
            delta_time = positive_number_return_float(
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
                    / single_step_operating_point_movement.base_operating_point.vCg__E
                )

            # Set the delta_time to be the average of the Airplanes' ideal delta times.
            delta_time = sum(delta_times) / len(delta_times)
        self.delta_time = delta_time

    
        num_steps = positive_int_return_int(
            num_steps, "num_steps"
        )
        
        self.num_steps = num_steps

        # Generate a list of lists of Airplanes that are the steps through each
        # AirplaneMovement. The first index identifies the AirplaneMovement, and the
        # second index identifies the time step.

    def generate_next_movement(self, base_airplanes, base_operating_point, step, deformation_matrices=None):
        """Creates the Airplanes and OperatingPoint at the next time step.
        :param base_airplanes: list of Airplanes

            This is the list of Airplanes at the base time step.
        :param base_operating_point: OperatingPoint

            This is the OperatingPoint at the base time step.
        :return: tuple (list of Airplanes, OperatingPoint)

            This is a tuple where the first element is the list of Airplanes at the
            next time step, and the second element is the OperatingPoint at the next
            time step.
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
