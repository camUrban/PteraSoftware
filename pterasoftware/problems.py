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
        self.num_steps = movement.num_steps
        self.delta_time = movement.delta_time
        self.only_final_results = only_final_results

        # If the user only wants the results for the final cycle, find the first time
        # step index where the solver should start calculating results. Otherwise,
        # set the first time step index to 0.
        if self.only_final_results:
            self.max_period = movement.get_max_period()

            # If the movement is static, set the first time step index at which to
            # calculate results to 0. Otherwise, the solver will calculate over the
            # range which encompasses the last full cycle.
            if self.max_period == 0:
                first_results_step = 0
            else:
                first_results_step = max(
                    0, math.floor(self.num_steps - (self.max_period / self.delta_time))
                )
        else:
            first_results_step = 0

        self.first_results_step = first_results_step

        # Initialize an empty list to hold the steady problems.
        self.steady_problems = []

        # Iterate through the problem's time steps.
        for step_id in range(self.num_steps):

            # Get the airplane objects and the operating point object associated with
            # this time step. This will be a list of the airplane snapshots associated
            # with each base airplane at this particular time step.
            these_airplanes = []
            for base_airplane_id in range(len(movement.airplanes)):
                these_airplanes.append(movement.airplanes[base_airplane_id][step_id])

            this_operating_point = movement.operating_points[step_id]

            # Initialize the steady problem object at this time step.
            this_steady_problem = SteadyProblem(
                airplanes=these_airplanes, operating_point=this_operating_point
            )

            # Append this steady problem to the list of steady problems.
            self.steady_problems.append(this_steady_problem)
