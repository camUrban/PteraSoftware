"""This module contains the class definitions for different types of problems.

This module contains the following classes:
    SteadyProblem: This is a class for steady aerodynamics problems.
    UnsteadyProblem: This is a class for unsteady aerodynamics problems.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import math

from . import geometry
from . import movement as mov
from . import operating_point as op
from . import parameter_validation


# TODO: Add unit tests for this class.
class SteadyProblem:
    """This is a class for steady aerodynamics problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, airplanes, operating_point):
        """This is the initialization method.

        :param airplanes: list of Airplanes

            This is the list of the Airplanes for this SteadyProblem.

        :param operating_point: OperatingPoint

            This is the OperatingPoint for this SteadyProblem.
        """
        if not isinstance(airplanes, list):
            raise TypeError("airplanes must be a list.")
        if len(airplanes) < 1:
            raise ValueError("airplanes must have at least one element.")
        for airplane in airplanes:
            if not isinstance(airplane, geometry.airplane.Airplane):
                raise TypeError("Every element in airplanes must be an Airplane.")
        self.airplanes = airplanes
        if not isinstance(operating_point, op.OperatingPoint):
            raise TypeError("operating_point must have be an OperatingPoint.")
        self.operating_point = operating_point


# TODO: Add unit tests for this class.
class UnsteadyProblem:
    """This is a class for unsteady aerodynamics problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, movement, only_final_results=False):
        """This is the initialization method.

        :param movement: Movement

            This is the Movement that contains this UnsteadyProblem's
            OperatingPointMovement and AirplaneMovements.

        :param only_final_results: boolLike, optional

            If set to True, the Solver will only calculate forces, moments,
            and pressures for the final complete cycle (of the Movement's
            sub-Movement with the longest period), which increases simulation speed.
            The default value is False.
        """
        if not isinstance(movement, mov.Movement):
            raise TypeError("movement must be a Movement.")
        self.movement = movement
        self.only_final_results = parameter_validation.boolLike_return_bool(
            only_final_results, "only_final_results"
        )

        self.num_steps = self.movement.num_steps
        self.delta_time = self.movement.delta_time

        # Find the maximum period of this UnsteadyProblem's Movement's sub-Movements.
        self.max_period = movement.get_max_period()

        # For UnsteadyProblems with a static Movement, users are typically interested
        # in the final time step's forces and moments, which, assuming convergence,
        # will be the most accurate. For UnsteadyProblems with cyclic movement,
        # (e.g. flapping wings) users are typically interested in the forces and
        # moments averaged over the last cycle simulated. Therefore, determine which
        # time step will be the first with relevant results based on if the Movement
        # is static or cyclic.
        if self.max_period == 0:
            self.first_averaging_step = self.num_steps - 1
        else:
            self.first_averaging_step = max(
                0, math.floor(self.num_steps - (self.max_period / self.delta_time))
            )

        # If the user only wants to calculate forces and moments for the final cycle
        # (for a cyclic Movement) or for the final time step (for a static Movement)
        # set the first step to calculate results to the first averaging step.
        # Otherwise, set it to the zero, which is the first time step.
        if self.only_final_results:
            self.first_results_step = self.first_averaging_step
        else:
            self.first_results_step = 0

        # Initialize empty lists to hold the final loads and load coefficients each
        # Airplane experiences. These will only be populated if this
        # UnsteadyProblem's Movement is static.
        self.final_near_field_forces_wind_axes = []
        self.final_near_field_force_coefficients_wind_axes = []
        self.final_near_field_moments_wind_axes = []
        self.final_near_field_moment_coefficients_wind_axes = []

        # Initialize empty lists to hold the final cycle-averaged loads and load
        # coefficients each Airplane experiences. These will only be populated if
        # this UnsteadyProblem's Movement is cyclic..
        self.final_mean_near_field_forces_wind_axes = []
        self.final_mean_near_field_force_coefficients_wind_axes = []
        self.final_mean_near_field_moments_wind_axes = []
        self.final_mean_near_field_moment_coefficients_wind_axes = []

        # Initialize empty lists to hold the final cycle-root-mean-squared loads and
        # load coefficients each airplane object experiences. These will only be
        # populated for variable geometry problems.
        self.final_rms_near_field_forces_wind_axes = []
        self.final_rms_near_field_force_coefficients_wind_axes = []
        self.final_rms_near_field_moments_wind_axes = []
        self.final_rms_near_field_moment_coefficients_wind_axes = []

        # Initialize an empty list to hold the SteadyProblems.
        self.steady_problems = []

        # Iterate through the UnsteadyProblem's time steps.
        for step_id in range(self.num_steps):

            # Get the Airplanes and the OperatingPoint associated with this time step.
            these_airplanes = []
            for this_base_airplane in movement.airplanes:
                these_airplanes.append(this_base_airplane[step_id])
            this_operating_point = movement.operating_points[step_id]

            # Initialize the SteadyProblem at this time step.
            this_steady_problem = SteadyProblem(
                airplanes=these_airplanes, operating_point=this_operating_point
            )

            # Append this SteadyProblem to the list of SteadyProblems.
            self.steady_problems.append(this_steady_problem)
