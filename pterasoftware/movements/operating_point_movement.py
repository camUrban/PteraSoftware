# NOTE: I haven't yet started refactoring this module.
"""This module contains the OperatingPointMovement class.

This module contains the following classes:
    OperatingPointMovement: This is a class used to contain the OperatingPoint
    movements.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

from . import functions

from .. import operating_point


# NOTE: I haven't yet started refactoring this class.
class OperatingPointMovement:
    """This is a class used to contain the movement characteristics of an operating
    point.

    This class contains the following public methods:
        generate_operating_points: This method creates the operating point objects at
        each time current_step, and groups them into a list.

        max_period: This method returns the longest period of this movement
        object's cycles.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I haven't yet started refactoring this method.
    def __init__(
        self,
        base_operating_point,
        velocity_amplitude=0.0,
        velocity_period=0.0,
        velocity_spacing="sine",
    ):
        """This is the initialization method.

        :param base_operating_point: OperatingPoint
            This is the operating point object, from which the others will be created.
        :param velocity_amplitude: float, optional
            This is the amplitude of the operating point's change in velocity. Its
            units are meters per second and its default value is 0 meters per second.
        :param velocity_period: float, optional
            This is the period of the operating point's change in its velocity. Its
            units are seconds and its default value is 0 seconds.
        :param velocity_spacing: string, optional
            This value determines the spacing of the operating point's change in its
            velocity. The options are "sine", and "uniform". The default value is
            "sine".
        """

        # Initialize the class attributes.
        self.base_operating_point = base_operating_point
        self.velocity_base = self.base_operating_point.vCg__E
        self.velocity_amplitude = velocity_amplitude
        self.velocity_period = velocity_period
        self.velocity_spacing = velocity_spacing

    # NOTE: I haven't yet started refactoring this method.
    def generate_operating_points(self, num_steps=10, delta_time=0.1):
        """This method creates the operating point objects at each time current_step,
        and groups them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value
            is 0.1 seconds.
        :return operating_points: list of OperatingPoint objects
            This is the list of OperatingPoint objects that is associated with this
            OperatingPointMovement object.
        """

        # Check the velocity spacing value.
        if self.velocity_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            velocity_list = functions.oscillating_sinspace(
                amplitude=self.velocity_amplitude,
                period=self.velocity_period,
                base_value=self.velocity_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.velocity_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            velocity_list = functions.oscillating_linspace(
                amplitude=self.velocity_amplitude,
                period=self.velocity_period,
                base_value=self.velocity_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of velocity_spacing!")

        # Create an empty list of operating points.
        operating_points = []

        # Generate the non-changing operating point attributes.
        density = self.base_operating_point.density
        alpha = self.base_operating_point.alpha
        beta = self.base_operating_point.beta

        # Iterate through the time steps.
        for step in range(num_steps):
            # Get the velocity at this time step.
            velocity = velocity_list[step]

            # Make a new operating point object for this time step.
            this_operating_point = operating_point.OperatingPoint(
                density=density, vCg__E=velocity, alpha=alpha, beta=beta
            )

            # Add this new object to the list of operating points.
            operating_points.append(this_operating_point)

        # Return the list of operating points.
        return operating_points

    # NOTE: I haven't yet started refactoring this method.
    @property
    def max_period(self):
        """This method returns the longest period of this movement object's cycles.

        :return max_period: float
            The longest period in seconds.
        """

        max_period = self.velocity_period

        return max_period
