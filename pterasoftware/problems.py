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

        :param airplanes: list of Airplane objects
            This is a list of the airplane objects for this problem.
        :param operating_point: OperatingPoint
            This is the operating point object for this problem.
        """
        # Initialize the problem's attributes.
        self.airplanes = airplanes
        self.operating_point = operating_point


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
            This is the movement object that contains this problem's airplane and
            operating point objects.
        :param only_final_results: Bool
            If this is set true, the solver will only calculate forces, moments,
            and pressures for the final complete cycle (of the movement with the
            longest period), which increases speed. The default value is False.
        """
        # Initialize the class attributes.
        self.movement = movement
        self.num_steps = self.movement.num_steps
        self.delta_time = self.movement.delta_time
        self.only_final_results = only_final_results

        # Find the maximum period of this problem's movements.
        self.max_period = movement.get_max_period()

        # For unsteady problems with static movement, users are typically interested
        # in the final time step's forces and moments, which, assuming convergence,
        # will be the most accurate. For unsteady problems with cyclic movement,
        # users are typically interested in the forces and movements averaged over
        # the last cycle simulated. Based on if the movement is static or cyclic,
        # find the first time step with relevant results.
        if self.max_period == 0:
            self.first_averaging_step = self.num_steps - 1
        else:
            self.first_averaging_step = max(
                0, math.floor(self.num_steps - (self.max_period / self.delta_time))
            )

        # If the user only wants to calculate forces and moments for the final cycle
        # (for cyclic motion) or for the final time step (for static movement) set
        # the first step to calculate results to the first averaging step. Otherwise,
        # set it to the zero, which is the first time step.
        if self.only_final_results:
            self.first_results_step = self.first_averaging_step
        else:
            self.first_results_step = 0

        # Initialize empty lists to hold the final loads and load coefficients each
        # airplane object experiences. These will only be populated for static
        # geometry problems.
        self.final_near_field_forces_wind_axes = []
        self.final_near_field_force_coefficients_wind_axes = []
        self.final_near_field_moments_wind_axes = []
        self.final_near_field_moment_coefficients_wind_axes = []

        # Initialize empty lists to hold the final cycle-averaged loads and load
        # coefficients each airplane object experiences. These will only be populated
        # for variable geometry problems.
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

        # Initialize an empty list to hold the steady problems.
        self.steady_problems = []

        # Iterate through the problem's time steps.
        for step_id in range(self.num_steps):

            # Get the airplane objects and the operating point object associated with
            # this time step. This will be a list of the airplane snapshots associated
            # with each base airplane at this particular time step.
            these_airplanes = []
            for this_base_airplane in movement.airplanes:
                these_airplanes.append(this_base_airplane[step_id])

            this_operating_point = movement.operating_points[step_id]

            # Initialize the steady problem object at this time step.
            this_steady_problem = SteadyProblem(
                airplanes=these_airplanes, operating_point=this_operating_point
            )

            # Append this steady problem to the list of steady problems.
            self.steady_problems.append(this_steady_problem)
