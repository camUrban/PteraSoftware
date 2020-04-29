"""This module contains the class definition for the geometry's movement.

This module contains the following classes:
    Movement: This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""


class Movement:
    """This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.

    """

    def __init__(self, movement_period, sweeping_amplitude):
        """This is the initialization method.

        :param movement_period: float
            This is the period of the problem's sweeping motion in seconds.
        :param sweeping_amplitude: float
            This is the amplitude of the problem's sweeping motion in radians.
        """

        # Initialize class attributes.
        self.movement_period = movement_period
        self.sweeping_amplitude = sweeping_amplitude
