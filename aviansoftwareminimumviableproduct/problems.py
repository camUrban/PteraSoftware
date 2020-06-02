
"""This module contains the class definitions for different types of problems.

This module contains the following classes:
    Problem: This class is used to contain aerodynamics problems to be analyzed.
    SteadyProblem: This is a subclass of the Problems class and is for steady aerodynamics problems to be analyzed.
    UnsteadyProblem: This is a subclass of the Problems class and is for unsteady aerodynamics problems to be analyzed.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""


class Problem:
    """This class is used to contain aerodynamics problems to be analyzed.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is meant to be subclassed by SteadyProblem and UnsteadyProblem.
    """

    def __init__(self, airplane, operating_point):
        """This is the initialization method.

        :param airplane: Airplane
            This is the airplane object for this problem.
        :param operating_point: OperatingPoint
            This is the operating point object for this problem.
        """

        # Initialize the attributes.
        self.airplane = airplane
        self.operating_point = operating_point


class SteadyProblem(Problem):
    """This is a subclass of the Problems class and is for steady aerodynamics problems to be analyzed.

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
            This is the airplane object for this problem.
        :param operating_point: OperatingPoint
            This is the operating point object for this problem.
        """

        # Call the parent class's initialization method
        super().__init__(airplane=airplane, operating_point=operating_point)


class UnsteadyProblem(Problem):
    """This is a subclass of the Problems class and is for unsteady aerodynamics problems to be analyzed.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, airplane, operating_point, movement):
        """This is the initialization method.

        :param airplane: Airplane
            This is the airplane object for this problem.
        :param operating_point: OperatingPoint
            This is the operating point object for this problem.
        :param movement: Movement
            This is the movement object for this problem's airplane.
        """

        # Call the parent class's initialization method
        super().__init__(airplane=airplane, operating_point=operating_point)

        # Initialize the other attributes.
        self.movement = movement

