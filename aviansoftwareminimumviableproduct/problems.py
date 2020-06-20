"""This module contains the class definitions for different types of problems.

This module contains the following classes:
    SteadyProblem: This is a class for steady aerodynamics problems.
    UnsteadyProblem: This is a class for unsteady aerodynamics problems.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""


class SteadyProblem:
    """This is a class for steady aerodynamics problems.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, airplane, operating_point):
        """This is the initialization method.

        :param airplane: Airplane
            This is the current_airplane object for this problem.
        :param operating_point: OperatingPoint
            This is the operating point object for this problem.
        """

        # Initialize the problem's attributes.
        self.airplane = airplane
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

    def __init__(self, movement):
        """This is the initialization method.

        :param movement: Movement
            This is the movement object that contains this problem's airplane and operating point objects.
        """

        # Initialize the class attributes for the number of time steps and the time in between these time steps.
        self.num_steps = movement.num_steps
        self.delta_time = movement.delta_time

        # Initialize an empty list to hold the steady problems.
        self.steady_problems = []

        # Iterate through the problem's time steps.
        for step in range(self.num_steps):

            # Get the airplane and operating point object at this time step.
            this_airplane = movement.airplanes[step]
            this_operating_point = movement.operating_points[step]

            # Initialize the steady problem object at this time step.
            this_steady_problem = SteadyProblem(
                airplane=this_airplane, operating_point=this_operating_point
            )

            # Append this steady problem to the list of steady problems.
            self.steady_problems.append(this_steady_problem)
