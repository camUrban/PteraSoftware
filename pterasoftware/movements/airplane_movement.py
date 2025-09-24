# NOTE: I haven't yet started refactoring this module.
"""This module contains the AirplaneMovement class.

This module contains the following classes:
    AirplaneMovement: This is a class used to contain the Airplane movements.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from . import functions

from .. import geometry


# NOTE: I haven't yet started refactoring this class.
class AirplaneMovement:
    """This is a class used to contain the movement characteristics of an airplane.

    This class contains the following public methods:
        generate_airplanes: This method creates the current_airplane object at each
        time current_step, and groups them into a list.

        max_period: This method returns the longest period of this movement
        object's sub-movement objects, sub-sub-movement objects, etc.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I haven't yet started refactoring this method.
    def __init__(
        self,
        base_airplane,
        wing_movements,
        x_ref_amplitude=0.0,
        x_ref_period=0.0,
        x_ref_spacing="sine",
        y_ref_amplitude=0.0,
        y_ref_period=0.0,
        y_ref_spacing="sine",
        z_ref_amplitude=0.0,
        z_ref_period=0.0,
        z_ref_spacing="sine",
    ):
        """This is the initialization method.

        :param base_airplane: Airplane
            This is the first airplane object, from which the others will be created.
        :param wing_movements: list of WingMovement objects
            This is a list of the WingMovement objects associated with each of the
            base airplane's wings.
        :param x_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its x reference point.
            Its units are meters and its default value is 0 meters.
        :param x_ref_period: float, optional
            This is the period of the airplane's change in its x reference point. Its
            units are seconds and its default value is 0 seconds.
        :param x_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its x
            reference point. The options are "sine", and "uniform". The default value
            is "sine".
        :param y_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its y reference point.
            Its units are meters and its default value is 0 meters.
        :param y_ref_period: float, optional
            This is the period of the airplane's change in its y reference point. Its
            units are seconds and its default value is 0 seconds.
        :param y_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its y
            reference point. The options are "sine", and "uniform". The default value
            is "sine".
        :param z_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its z reference point.
            Its units are meters and its default value is 0 meters.
        :param z_ref_period: float, optional
            This is the period of the airplane's change in its z reference point. Its
            units are seconds and its default value is 0 seconds.
        :param z_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its z
            reference point. The options are "sine", and "uniform". The default value
            is "sine".
        """

        # Initialize the class attributes.
        self.base_airplane = base_airplane
        self.wing_movements = wing_movements
        self.x_ref_base = self.base_airplane.x_ref
        self.x_ref_amplitude = x_ref_amplitude
        self.x_ref_period = x_ref_period
        self.x_ref_spacing = x_ref_spacing
        self.y_ref_base = self.base_airplane.y_ref
        self.y_ref_amplitude = y_ref_amplitude
        self.y_ref_period = y_ref_period
        self.y_ref_spacing = y_ref_spacing
        self.z_ref_base = self.base_airplane.z_ref
        self.z_ref_amplitude = z_ref_amplitude
        self.z_ref_period = z_ref_period
        self.z_ref_spacing = z_ref_spacing

    # NOTE: I haven't yet started refactoring this method.
    def generate_airplanes(self, num_steps=10, delta_time=0.1):
        """This method creates the current_airplane object at each time current_step,
        and groups them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value
            is 0.1 seconds.
        :return airplanes: list of Airplane objects
            This is the list of Airplane objects that is associated with this
            AirplaneMovement object.
        """

        # Check the x_ref spacing value.
        if self.x_ref_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            x_ref_list = functions.oscillating_sinspace(
                amplitude=self.x_ref_amplitude,
                period=self.x_ref_period,
                base_value=self.x_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.x_ref_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            x_ref_list = functions.oscillating_linspace(
                amplitude=self.x_ref_amplitude,
                period=self.x_ref_period,
                base_value=self.x_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of x_ref_spacing!")

        # Check the y_ref spacing value.
        if self.y_ref_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            y_ref_list = functions.oscillating_sinspace(
                amplitude=self.y_ref_amplitude,
                period=self.y_ref_period,
                base_value=self.y_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.y_ref_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            y_ref_list = functions.oscillating_linspace(
                amplitude=self.y_ref_amplitude,
                period=self.y_ref_period,
                base_value=self.y_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of y_ref_spacing!")

        # Check the z_ref spacing value.
        if self.z_ref_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            z_ref_list = functions.oscillating_sinspace(
                amplitude=self.z_ref_amplitude,
                period=self.z_ref_period,
                base_value=self.z_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.z_ref_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            z_ref_list = functions.oscillating_linspace(
                amplitude=self.z_ref_amplitude,
                period=self.z_ref_period,
                base_value=self.z_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of z_ref_spacing!")

        # Create an empty array that will hold each of the airplane's wing's vector
        # of other wing's based its movement.
        wings = np.empty((len(self.wing_movements), num_steps), dtype=object)

        # Iterate through the wing movement locations.
        for wing_movement_location, wing_movement in enumerate(self.wing_movements):

            # Generate this wing's vector of other wing's based on its movement.
            this_wings_list_of_wings = np.array(
                wing_movement.generate_wings(num_steps=num_steps, delta_time=delta_time)
            )

            # Add this vector the airplane's array of wing objects.
            wings[wing_movement_location, :] = this_wings_list_of_wings

        # Create an empty list of airplanes.
        airplanes = []

        # Generate the airplane name.
        name = self.base_airplane.name

        # Iterate through the time steps.
        for step in range(num_steps):
            # Get the reference position at this time step.
            x_ref = x_ref_list[step]
            y_ref = y_ref_list[step]
            z_ref = z_ref_list[step]
            these_wings = wings[:, step]

            # Make a new airplane object for this time step.
            this_airplane = geometry.airplane.Airplane(
                name=name, x_ref=x_ref, y_ref=y_ref, z_ref=z_ref, wings=these_wings
            )

            # Add this new object to the list of airplanes.
            airplanes.append(this_airplane)

        # Return the list of airplanes.
        return airplanes

    # NOTE: I haven't yet started refactoring this method.
    @property
    def max_period(self):
        """This method returns the longest period of this movement object's sub-
        movement objects, sub-sub-movement objects, etc.

        :return max_period: float
            The longest period in seconds.
        """

        wing_movement_max_periods = []
        for wing_movement in self.wing_movements:
            wing_movement_max_periods.append(wing_movement.max_period)
        max_wing_movement_period = max(wing_movement_max_periods)

        max_period = max(
            max_wing_movement_period,
            self.x_ref_period,
            self.y_ref_period,
            self.z_ref_period,
        )

        return max_period
