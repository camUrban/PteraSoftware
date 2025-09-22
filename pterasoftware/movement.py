# NOTE: I haven't yet started refactoring this module.
"""This module contains the class definitions for the problem's movement.

This module contains the following classes:
    Movement: This is a class used to contain the movement characteristics of an
    unsteady aerodynamics problem.

    AirplaneMovement: This is a class used to contain the movement characteristics of
    an airplane.

    WingMovement: This is a class used to contain the movement characteristics of a
    wing.

    WingCrossSectionMovement: This is a class used to contain the movement
    characteristics of a wing cross section.

    OperatingPointMovement: This is a class used to contain the movement
    characteristics of an operating point.

This module contains the following exceptions:
    None

This module contains the following functions:
    oscillating_sinspace: This function returns a 1D array of values that are
    calculated by inputting a vector of linearly spaced time steps into a sine
    function.

    oscillating_linspace: This function returns a 1D array of values that are
    calculated by inputting a vector of linearly spaced time steps into a triangle
    function.

    oscillating_customspace: This function returns a 1D array of values that are
    calculated by inputting a vector of linearly spaced time steps into a custom
    function.
"""

import math

import numpy as np
from scipy import signal

from . import geometry
from . import operating_point


# NOTE: I haven't yet started refactoring this class.
class Movement:
    """This is a class used to contain the movement characteristics of an unsteady
    aerodynamics problem.

    This class contains the following public methods:
        get_max_period: This method returns the longest period of this movement
        object's sub-movement objects, sub-sub-movement objects, etc.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I haven't yet started refactoring this method.
    def __init__(
        self,
        airplane_movements,
        operating_point_movement,
        num_steps=None,
        num_cycles=None,
        num_chords=None,
        delta_time=None,
    ):
        """This is the initialization method.

        :param airplane_movements: list of AirplaneMovement objects
            This is a list of objects which characterize the movement of
            each airplane in the problem.
        :param operating_point_movement: OperatingPointMovement
            This object characterizes the movement of the operating point.
        :param num_steps: int, optional
            This integer is the number of time steps of the unsteady simulation. If
            not given a value, and the movement is dynamic, this method will
            calculate one such that the simulation will cover some number of cycles
            of the maximum period movement. The number of cycles defaults to three,
            but can be changed with the num_cycles parameter. If not given a value,
            and the movement is static, the number of steps will default to the
            number of time steps such that the wake extends back by some number of
            reference chord lengths. The number of chord lengths defaults to ten,
            but can be changed with the num_chords parameter.
        :param num_cycles: int, optional
            This integer is the number of cycles of the maximum period movement used
            to calculate a non-populated num_steps parameter. This parameter is only
            used if the num_steps parameter is None, and the movement isn't static.
            The default value is None.
        :param num_chords: int, optional
            This integer is the number of reference chord lengths used to calculate a
            non-populated num_steps parameter. This parameter is only used if the
            num_steps parameter is None, and the movement is static. The default
            value is None.
        :param delta_time: float, optional
            This float is the time, in seconds, between each time current_step. If
            not given a value, this method will calculate one such the ring vortices
            shed off the main wing will have roughly the same chord length as the
            panels on the main wing. This is based on the base airplane's reference
            chord length, its main wing's number of chordwise panels, and its base
            operating point's velocity.
        """

        # Initialize the class attributes.
        self.airplane_movements = airplane_movements
        self.operating_point_movement = operating_point_movement

        # If the number of cycles were specified, make sure that the number of steps
        # isn't also specified and that the movement isn't static.
        if num_cycles is not None:
            if num_steps is not None or self.get_max_period() == 0:
                raise Exception(
                    "Only specify the number of cycles if you haven't specified the "
                    "number of steps and the movement isn't static!"
                )
            self.num_cycles = num_cycles
        else:
            self.num_cycles = None

        # If the number of chords were specified, make sure that the number of steps
        # isn't also specified and that the movement is static.
        if num_chords is not None:
            if num_steps is not None or self.get_max_period() != 0:
                raise Exception(
                    "Only specify the number of chords if you haven't specified the "
                    "number of steps and the movement is static!"
                )
            self.num_chords = num_chords
        else:
            self.num_chords = None

        # Calculate default num_steps and delta_time values if the user hasn't passed
        # one in.
        if delta_time is None:
            delta_times = []
            for airplane_movement in self.airplane_movements:

                # For a given airplane object, the ideal time step length is that
                # which sheds ring vortices off the main wing that have roughly the
                # same chord length as the panels on the main wing. This is based on
                # the base airplane's reference chord length, its main wing's number
                # of chordwise panels, and its base operating point's velocity.
                delta_times.append(
                    airplane_movement.base_airplane.c_ref
                    / airplane_movement.base_airplane.wings[0].num_chordwise_panels
                    / operating_point_movement.base_operating_point.vCg__E
                )

            # Set the delta time to be the average of the airplanes' ideal delta times.
            delta_time = sum(delta_times) / len(delta_times)

        if num_steps is None:

            # Get the maximum period of this movement's sub-movements.
            max_period = self.get_max_period()

            if max_period == 0:

                # Find the value of the largest reference chord length of all the
                # airplanes in this problem.
                c_refs = []
                for airplane_movement in self.airplane_movements:
                    c_refs.append(airplane_movement.base_airplane.c_ref)
                max_c_ref = max(c_refs)

                # If the user didn't specify the number of reference chord lengths
                # for which to convect the wake, set it to ten.
                if self.num_chords is None:
                    self.num_chords = 10

                # If the movement is static, then set the number of time steps such
                # that the wake extends back by some number of reference chord lengths.
                wake_length = self.num_chords * max_c_ref
                panel_length = (
                    delta_time
                    * self.operating_point_movement.base_operating_point.vCg__E
                )
                num_steps = math.ceil(wake_length / panel_length)
            else:

                # If the user didn't specify the number of cycles of the movement
                # with the maximum period over which to simulate, set it to three.
                if self.num_cycles is None:
                    self.num_cycles = 3

                # The default number of time steps is that simulates some number of
                # cycles of the movement with the maximum period.
                num_steps = math.ceil(self.num_cycles * max_period / delta_time)

        self.num_steps = num_steps
        self.delta_time = delta_time

        # Generate a list of lists of airplane objects that are the steps through the
        # movement of each of this problem's base airplanes. The first index
        # identifies the base airplane and the second index identifies the time step.
        self.airplanes = []
        for airplane_movement in self.airplane_movements:
            self.airplanes.append(
                airplane_movement.generate_airplanes(
                    num_steps=self.num_steps, delta_time=self.delta_time
                )
            )

        # Generate a lists of operating point objects that are the steps through the
        # movement of this problem's operating point.
        self.operating_points = operating_point_movement.generate_operating_points(
            num_steps=self.num_steps, delta_time=self.delta_time
        )

    # NOTE: I haven't yet started refactoring this method.
    def get_max_period(self):
        """This method returns the longest period of this movement object's sub-
        movement objects, sub-sub-movement objects, etc.

        :return: float
            The longest period in seconds.
        """
        # Iterate through the airplane movements and find the one with the largest
        # max period.
        max_airplane_periods = []
        for airplane_movement in self.airplane_movements:
            max_airplane_periods.append(airplane_movement.get_max_period())
        max_airplane_period = max(max_airplane_periods)

        # The global max period is the maximum of the max airplane period and the max
        # operating point period.
        return max(
            max_airplane_period,
            self.operating_point_movement.get_max_period(),
        )


# NOTE: I haven't yet started refactoring this class.
class AirplaneMovement:
    """This is a class used to contain the movement characteristics of an airplane.

    This class contains the following public methods:
        generate_airplanes: This method creates the current_airplane object at each
        time current_step, and groups them into a list.

        get_max_period: This method returns the longest period of this movement
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
            x_ref_list = oscillating_sinspace(
                amplitude=self.x_ref_amplitude,
                period=self.x_ref_period,
                base_value=self.x_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.x_ref_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            x_ref_list = oscillating_linspace(
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
            y_ref_list = oscillating_sinspace(
                amplitude=self.y_ref_amplitude,
                period=self.y_ref_period,
                base_value=self.y_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.y_ref_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            y_ref_list = oscillating_linspace(
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
            z_ref_list = oscillating_sinspace(
                amplitude=self.z_ref_amplitude,
                period=self.z_ref_period,
                base_value=self.z_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.z_ref_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            z_ref_list = oscillating_linspace(
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
    def get_max_period(self):
        """This method returns the longest period of this movement object's sub-
        movement objects, sub-sub-movement objects, etc.

        :return max_period: float
            The longest period in seconds.
        """

        wing_movement_max_periods = []
        for wing_movement in self.wing_movements:
            wing_movement_max_periods.append(wing_movement.get_max_period())
        max_wing_movement_period = max(wing_movement_max_periods)

        max_period = max(
            max_wing_movement_period,
            self.x_ref_period,
            self.y_ref_period,
            self.z_ref_period,
        )

        return max_period


# NOTE: I haven't yet started refactoring this class.
class WingMovement:
    """This is a class used to contain the movement characteristics of a wing.

    This class contains the following public methods:
        generate_wings: This method creates the wing object at each time step,
        and groups them into a list.

        get_max_period: This method returns the longest period of this movement
        object's sub-movement objects, sub-sub-movement objects, etc.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I haven't yet started refactoring this method.
    def __init__(
        self,
        base_wing,
        wing_cross_sections_movements,
        x_le_amplitude=0.0,
        x_le_period=0.0,
        x_le_spacing="sine",
        y_le_amplitude=0.0,
        y_le_period=0.0,
        y_le_spacing="sine",
        z_le_amplitude=0.0,
        z_le_period=0.0,
        z_le_spacing="sine",
    ):
        """This is the initialization method.

        :param base_wing: Wing
            This is the first wing object, from which the others will be created.
        :param wing_cross_sections_movements: list of WingCrossSectionMovement objects
            This is a list of the WingCrossSectionMovement objects associated with
            each of the base wing's cross sections.
        :param x_le_amplitude: float, optional
            This is the amplitude of the wing's change in its x reference point. Its
            units are meters and its default value is 0 meters.
        :param x_le_period: float, optional
            This is the period of the wing's change in its x reference point. Its
            units are seconds and its default value is 0 seconds.
        :param x_le_spacing: string, optional
            This value determines the spacing of the wing's change in its x reference
            point. The options are "sine", and "uniform". The default value is "sine".
        :param y_le_amplitude: float, optional
            This is the amplitude of the wing's change in its y reference point. Its
            units are meters and its default value is 0 meters.
        :param y_le_period: float, optional
            This is the period of the wing's change in its y reference point. Its
            units are seconds and its default value is 0 seconds.
        :param y_le_spacing: string, optional
            This value determines the spacing of the wing's change in its y reference
            point. The options are "sine", and "uniform". The default value is "sine".
        :param z_le_amplitude: float, optional
            This is the amplitude of the wing's change in its z reference point. Its
            units are meters and its default value is 0 meters.
        :param z_le_period: float, optional
            This is the period of the wing's change in its z reference point. Its
            units are seconds and its default value is 0 seconds.
        :param z_le_spacing: string, optional
            This value determines the spacing of the wing's change in its z reference
            point. The options are "sine", and "uniform". The default value is "sine".
        """
        # Initialize the class attributes.
        self.base_wing = base_wing
        self.wing_cross_section_movements = wing_cross_sections_movements
        self.x_le_base = self.base_wing.x_le
        self.x_le_amplitude = x_le_amplitude
        self.x_le_period = x_le_period
        self.x_le_spacing = x_le_spacing
        self.y_le_base = self.base_wing.y_le
        self.y_le_amplitude = y_le_amplitude
        self.y_le_period = y_le_period
        self.y_le_spacing = y_le_spacing
        self.z_le_base = self.base_wing.z_le
        self.z_le_amplitude = z_le_amplitude
        self.z_le_period = z_le_period
        self.z_le_spacing = z_le_spacing

    # NOTE: I haven't yet started refactoring this method.
    def generate_wings(self, num_steps=10, delta_time=0.1):
        """This method creates the wing object at each time current_step, and groups
        them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value
            is 0.1 seconds.
        :return wings: list of Wing objects
            This is the list of Wing objects that is associated with this
            WingMovement object.
        """
        # Check the x_le spacing value.
        if self.x_le_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            x_le_list = oscillating_sinspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.x_le_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            x_le_list = oscillating_linspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of x_le_spacing!")

        # Check the y_le spacing value.
        if self.y_le_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            y_le_list = oscillating_sinspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.y_le_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            y_le_list = oscillating_linspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of y_le_spacing!")

        # Check the z_le spacing value.
        if self.z_le_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            z_le_list = oscillating_sinspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.z_le_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            z_le_list = oscillating_linspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of z_le_spacing!")

        # Create an empty array that will hold each of the wing's wing cross
        # section's vector of other wing cross section's based its movement.
        wing_cross_sections = np.empty(
            (len(self.wing_cross_section_movements), num_steps), dtype=object
        )

        # Initialize a variable to hold the inner wing cross section's list of wing
        # cross sections for each time step.
        last_wing_cross_section_time_histories = None

        # Iterate through the wing cross section movement locations.
        for (
            wing_cross_section_movement_location,
            wing_cross_section_movement,
        ) in enumerate(self.wing_cross_section_movements):
            wing_is_vertical = False

            # Check if this is this wing's root cross section.
            if wing_cross_section_movement_location == 0:

                # Get the root cross section's sweeping and heaving attributes.
                first_wing_cross_section_movement_sweeping_amplitude = (
                    wing_cross_section_movement.sweeping_amplitude
                )
                first_wing_cross_section_movement_sweeping_period = (
                    wing_cross_section_movement.sweeping_period
                )
                first_wing_cross_section_movement_heaving_amplitude = (
                    wing_cross_section_movement.heaving_amplitude
                )
                first_wing_cross_section_movement_heaving_period = (
                    wing_cross_section_movement.heaving_period
                )

                # Check that the root cross section is not sweeping or heaving.
                assert first_wing_cross_section_movement_sweeping_amplitude == 0
                assert first_wing_cross_section_movement_sweeping_period == 0
                assert first_wing_cross_section_movement_heaving_amplitude == 0
                assert first_wing_cross_section_movement_heaving_period == 0

                # Set the variables relating this wing cross section to the inner
                # wing cross section to zero because this is the innermost wing cross
                # section
                wing_cross_section_span = 0.0
                base_wing_cross_section_sweep = 0.0
                base_wing_cross_section_heave = 0.0
                last_x_les = np.zeros(num_steps) * 0.0
                last_y_les = np.zeros(num_steps) * 0.0
                last_z_les = np.zeros(num_steps) * 0.0

            else:
                this_base_wing_cross_section = (
                    wing_cross_section_movement.base_wing_cross_section
                )

                this_x_le = this_base_wing_cross_section.x_le
                this_y_le = this_base_wing_cross_section.y_le
                this_z_le = this_base_wing_cross_section.z_le

                # Initialize variables to hold the inner wing cross section's time
                # histories of its leading edge coordinates.
                last_x_les = []
                last_y_les = []
                last_z_les = []

                # Iterate through the inner wing cross section's time history and
                # populate the leading edge coordinate variables.
                for last_wing_cross_section in last_wing_cross_section_time_histories:
                    last_x_les.append(last_wing_cross_section.x_le)
                    last_y_les.append(last_wing_cross_section.y_le)
                    last_z_les.append(last_wing_cross_section.z_le)

                # Find the span between this wing cross section and the inner wing
                # cross section.
                wing_cross_section_span = np.sqrt(
                    (this_x_le - last_x_les[0]) ** 2
                    + (this_y_le - last_y_les[0]) ** 2
                    + (this_z_le - last_z_les[0]) ** 2
                )

                if this_y_le != last_y_les[0]:
                    # Find the base sweep angle of this wing cross section compared
                    # to the inner wing cross section at the first time step.
                    base_wing_cross_section_sweep = (
                        np.arctan(
                            (this_z_le - last_z_les[0]) / (this_y_le - last_y_les[0])
                        )
                        * 180
                        / np.pi
                    )

                    # Find the base heave angle of this wing cross section compared
                    # to the inner wing cross section at the first time step.
                    base_wing_cross_section_heave = (
                        np.arctan(
                            (this_x_le - last_x_les[0]) / (this_y_le - last_y_les[0])
                        )
                        * 180
                        / np.pi
                    )
                else:
                    base_wing_cross_section_sweep = 0.0
                    base_wing_cross_section_heave = 0.0
                    wing_is_vertical = True

            # Generate this wing cross section's vector of wing cross sections at
            # each time step based on its movement.
            this_wing_cross_sections_list_of_wing_cross_sections = np.array(
                wing_cross_section_movement.generate_wing_cross_sections(
                    num_steps=num_steps,
                    delta_time=delta_time,
                    cross_section_span=wing_cross_section_span,
                    cross_section_sweep=base_wing_cross_section_sweep,
                    cross_section_heave=base_wing_cross_section_heave,
                    last_x_les=last_x_les,
                    last_y_les=last_y_les,
                    last_z_les=last_z_les,
                    wing_is_vertical=wing_is_vertical,
                )
            )

            # Add this vector the wing's array of wing cross section objects.
            wing_cross_sections[wing_cross_section_movement_location, :] = (
                this_wing_cross_sections_list_of_wing_cross_sections
            )

            # Update the inner wing cross section's list of wing cross sections for
            # each time step.
            last_wing_cross_section_time_histories = (
                this_wing_cross_sections_list_of_wing_cross_sections
            )

        # Create an empty list of wings.
        wings = []

        # Generate the non-changing wing attributes.
        name = self.base_wing.name
        symmetric = self.base_wing.symmetric
        num_chordwise_panels = self.base_wing.num_chordwise_panels
        chordwise_spacing = self.base_wing.chordwise_spacing

        # Iterate through the time steps.
        for step in range(num_steps):

            # Get the reference position at this time step.
            x_le = x_le_list[step]
            y_le = y_le_list[step]
            z_le = z_le_list[step]
            cross_sections = wing_cross_sections[:, step]

            # Make a new wing object for this time step.
            this_wing = geometry.wing.Wing(
                name=name,
                x_le=x_le,
                y_le=y_le,
                z_le=z_le,
                wing_cross_sections=cross_sections,
                symmetric=symmetric,
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing,
            )

            # Add this new object to the list of wings.
            wings.append(this_wing)

        # Return the list of wings.
        return wings

    # NOTE: I haven't yet started refactoring this method.
    def get_max_period(self):
        """This method returns the longest period of this movement object's
        sub-movement objects, sub-sub-movement objects, etc.

        :return max_period: float
            The longest period in seconds.
        """

        wing_cross_section_movement_max_periods = []
        for wing_cross_section_movement in self.wing_cross_section_movements:
            wing_cross_section_movement_max_periods.append(
                wing_cross_section_movement.get_max_period()
            )
        max_wing_cross_section_movement_period = max(
            wing_cross_section_movement_max_periods
        )

        max_period = max(
            max_wing_cross_section_movement_period,
            self.x_le_period,
            self.y_le_period,
            self.z_le_period,
        )

        return max_period


# NOTE: I haven't yet started refactoring this class.
class WingCrossSectionMovement:
    """This is a class used to contain the movement characteristics of a wing cross
    section.

    This class contains the following public methods:
        generate_wing_cross_sections: This method creates the wing cross section
        objects at each time current_step, and groups them into a list.

        get_max_period: This method returns the longest period of this movement
        object's cycles.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # NOTE: I haven't yet started refactoring this method.
    def __init__(
        self,
        base_wing_cross_section,
        sweeping_amplitude=0.0,
        sweeping_period=0.0,
        sweeping_spacing="sine",
        custom_sweep_function=None,
        pitching_amplitude=0.0,
        pitching_period=0.0,
        pitching_spacing="sine",
        custom_pitch_function=None,
        heaving_amplitude=0.0,
        heaving_period=0.0,
        heaving_spacing="sine",
        custom_heave_function=None,
    ):
        """This is the initialization method.

        :param base_wing_cross_section: WingCrossSection
            This is the first wing cross section object, from which the others will
            be created.
        :param sweeping_amplitude: float, optional
            This is the amplitude of the cross section's change in its sweep,
            relative to the vehicle's body axes. Its units are degrees and its
            default value is 0.0 degrees.
        :param sweeping_period: float, optional
            This is the period of the cross section's change in its sweep. Its units
            are seconds and its default value is 0.0 seconds.
        :param sweeping_spacing: string, optional
            This value determines the spacing of the cross section's change in its
            sweep. The options are "sine", "uniform", and "custom". The default value
            is "sine". If "custom", then the value of custom_sweep_function must not
            be none. If both sweeping_spacing and custom_sweep_function are not none,
            then the value of sweeping_spacing will take precedence.
        :param custom_sweep_function: function, optional
            This is a function that describes the motion of the sweeping. For
            example, it could be np.cos or np.sinh (assuming numpy had previously
            been imported as np). It will be horizontally scaled by the
            sweeping_period, vertically scaled by the sweeping_amplitude. For
            example, say the function has an amplitude of 2 units, a period of 3
            units, sweeping_amplitude is set to 4 units and sweeping_period is set to
            5 units. The sweeping motion will have a net amplitude of 8 units and a
            net period of 15 units.
        :param pitching_amplitude: float, optional
            This is the amplitude of the cross section's change in its pitch,
            relative to the vehicle's body axes. Its units are degrees and its
            default value is 0.0 degrees.
        :param pitching_period: float, optional
            This is the period of the cross section's change in its pitch. Its units
            are seconds and its default value is 0.0 seconds.
        :param pitching_spacing: string, optional
            This value determines the spacing of the cross section's change in its
            pitch. The options are "sine", "uniform", and "custom". The default value
            is "sine". If "custom", then the value of custom_pitch_function must not
            be none. If both pitching_spacing and custom_pitch_function are not none,
            then the value of pitching_spacing will take precedence.
        :param custom_pitch_function: function, optional
            This is a function that describes the motion of the pitching. For
            example, it could be np.cos or np.sinh (assuming numpy had previously
            been imported as np). It will be horizontally scaled by the
            pitching_period, vertically scaled by the pitching_amplitude. For
            example, say the function has an amplitude of 2 units, a period of 3
            units, pitching_amplitude is set to 4 units and pitching_period is set to
            5 units. The pitching motion will have a net amplitude of 8 units and a
            net period of 15 units.
        :param heaving_amplitude: float, optional
            This is the amplitude of the cross section's change in its heave,
            relative to the vehicle's body axes. Its units are degrees and its
            default value is 0.0 degrees.
        :param heaving_period: float, optional
            This is the period of the cross section's change in its heave. Its units
            are seconds and its default value is 0.0 seconds.
        :param heaving_spacing: string, optional
            This value determines the spacing of the cross section's change in its
            heave. The options are "sine", "uniform", and "custom". The default value
            is "sine". If "custom", then the value of custom_heave_function must not
            be none. If both heaving_spacing and custom_heave_function are not none,
            then the value of heaving_spacing will take precedence.
        :param custom_heave_function: function, optional
            This is a function that describes the motion of the heaving. For example,
            it could be np.cos or np.sinh (assuming numpy had previously been
            imported as np). It will be horizontally scaled by the heaving_period,
            vertically scaled by the heaving_amplitude. For example, say the function
            has an amplitude of 2 units, a period of 3 units, heaving_amplitude is
            set to 4 units and heaving_period is set to 5 units. The heaving motion
            will have a net amplitude of 8 units and a net period of 15 units.
        """

        # Initialize the class attributes.
        self.base_wing_cross_section = base_wing_cross_section
        self.sweeping_amplitude = sweeping_amplitude
        self.sweeping_period = sweeping_period
        self.sweeping_spacing = sweeping_spacing
        self.custom_sweep_function = custom_sweep_function
        self.sweeping_base = 0.0
        self.pitching_amplitude = pitching_amplitude
        self.pitching_period = pitching_period
        self.pitching_spacing = pitching_spacing
        self.custom_pitch_function = custom_pitch_function
        self.pitching_base = self.base_wing_cross_section.twist
        self.heaving_amplitude = heaving_amplitude
        self.heaving_period = heaving_period
        self.heaving_spacing = heaving_spacing
        self.custom_heave_function = custom_heave_function
        self.heaving_base = 0.0
        self.x_le_base = self.base_wing_cross_section.x_le
        self.y_le_base = self.base_wing_cross_section.y_le
        self.z_le_base = self.base_wing_cross_section.z_le
        self.twist_base = self.base_wing_cross_section.twist
        self.control_surface_deflection_base = (
            self.base_wing_cross_section.control_surface_deflection
        )

    # NOTE: I haven't yet started refactoring this method.
    def generate_wing_cross_sections(
        self,
        num_steps=10,
        delta_time=0.1,
        last_x_les=None,
        last_y_les=None,
        last_z_les=None,
        wing_is_vertical=False,
        cross_section_span=0.0,
        cross_section_sweep=0.0,
        cross_section_heave=0.0,
    ):
        """This method creates the wing cross section objects at each time
        current_step, and groups them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value
            is 0.1 seconds.
        :param last_x_les: float, optional
            This is an array of the x coordinates of the reference location of the
            previous cross section at each time step. Its units are in meters,
            and its default value is 0.0 meters.
        :param last_y_les: float, optional
            This is an array of the y coordinates of the reference location of the
            previous cross section at each time step. Its units are in meters,
            and its default value is 0.0 meters.
        :param last_z_les: float, optional
            This is an array of the z coordinates of the reference location of the
            previous cross section at each time step. Its units are in meters,
            and its default value is 0.0 meters.
        :param wing_is_vertical: bool, optional
            This flag is set to true if the wing containing this wing cross section
            is vertical. If true, the cross section's movement will automatically be
            eliminated. This is a temporary patch until vertical wing cross section
            movement is supported. The default value is false.
        :param cross_section_span: float, optional
            This is the length, in meters, of the leading edge stretching between
            this wing cross section at the previous wing cross section. If this is
            the first cross section, it should be 0.0 meters. The default value is
            0.0 meters.
        :param cross_section_sweep: float, optional
            This is the sweep, in degrees, of this wing cross section relative to the
            previous wing cross section. If this is the first cross section,
            it should be 0.0 degrees. The default value is 0.0 degrees.
        :param cross_section_heave: float, optional
            This is the heave, in degrees, of this wing cross section relative to the
            previous wing cross section. If this is the first cross section,
            it should be 0.0 degrees. The default value is 0.0 degrees.
        :return wing_cross_sections: list of WingCrossSection objects
            This is the list of WingCrossSection objects that is associated with this
            WingCrossSectionMovement object.
        """

        # Check the sweeping spacing value.
        if self.sweeping_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            sweeping_list = oscillating_sinspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=cross_section_sweep,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.sweeping_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            sweeping_list = oscillating_linspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=cross_section_sweep,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.sweeping_spacing == "custom":

            # Raise an exception if the user did not declare a custom sweep function.
            if self.custom_sweep_function is None:
                raise Exception(
                    "You can't declare custom sweep spacing without providing a "
                    "custom sweep function."
                )

            # Create an array of points with a uniform spacing.
            sweeping_list = oscillating_customspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=cross_section_sweep,
                num_steps=num_steps,
                delta_time=delta_time,
                custom_function=self.custom_sweep_function,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of sweeping_spacing!")

        # Check the pitching spacing value.
        if self.pitching_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            pitching_list = oscillating_sinspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.pitching_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            pitching_list = oscillating_linspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.pitching_spacing == "custom":

            # Raise an exception if the user did not declare a custom pitch function.
            if self.custom_pitch_function is None:
                raise Exception(
                    "You can't declare custom pitch spacing without providing a "
                    "custom pitch function."
                )

            # Create an array of points with a uniform spacing.
            pitching_list = oscillating_customspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time,
                custom_function=self.custom_pitch_function,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of pitching_spacing!")

        # Check the heaving spacing value.
        if self.heaving_spacing == "sine":

            # Create an array of points with a sinusoidal spacing.
            heaving_list = oscillating_sinspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=cross_section_heave,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.heaving_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            heaving_list = oscillating_linspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=cross_section_heave,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.heaving_spacing == "custom":

            # Raise an exception if the user did not declare a custom heave function.
            if self.custom_heave_function is None:
                raise Exception(
                    "You can't declare custom heave spacing without providing a "
                    "custom heave function."
                )

            # Create an array of points with custom spacing.
            heaving_list = oscillating_customspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=cross_section_heave,
                num_steps=num_steps,
                delta_time=delta_time,
                custom_function=self.custom_heave_function,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of heaving_spacing!")

        if wing_is_vertical:
            x_le_list = np.ones(num_steps) * self.x_le_base
            y_le_list = np.ones(num_steps) * self.y_le_base
            z_le_list = np.ones(num_steps) * self.z_le_base
            twist_list = np.ones(num_steps) * self.twist_base
        else:

            # Find the list of new leading edge points. This uses a spherical
            # coordinate transformation, referencing the previous wing cross
            # section's leading edge point (at each time step) as the origin. Also
            # convert the lists of sweep, pitch, and heave values to radians before
            # passing them into numpy's trigonometry functions.
            x_le_list = last_x_les + cross_section_span * np.cos(
                sweeping_list * np.pi / 180
            ) * np.sin(heaving_list * np.pi / 180)
            y_le_list = last_y_les + cross_section_span * np.cos(
                sweeping_list * np.pi / 180
            ) * np.cos(heaving_list * np.pi / 180)
            z_le_list = last_z_les + cross_section_span * np.sin(
                sweeping_list * np.pi / 180
            )
            twist_list = pitching_list

        # Create an empty list of wing cross sections.
        wing_cross_sections = []

        # Generate the non-changing wing cross section attributes.
        chord = self.base_wing_cross_section.chord
        airfoil = self.base_wing_cross_section.airfoil
        control_surface_deflection = self.control_surface_deflection_base
        control_surface_type = self.base_wing_cross_section.control_surface_type
        control_surface_hinge_point = (
            self.base_wing_cross_section.control_surface_hinge_point
        )
        num_spanwise_panels = self.base_wing_cross_section.num_spanwise_panels
        spanwise_spacing = self.base_wing_cross_section.spanwise_spacing

        # Iterate through the time steps.
        for step in range(num_steps):
            # Get the changing wing cross section attributes at this time step.
            x_le = x_le_list[step]
            y_le = y_le_list[step]
            z_le = z_le_list[step]
            twist = twist_list[step]

            # Make a new wing cross section object for this time step.
            this_wing_cross_section = geometry.wing_cross_section.WingCrossSection(
                x_le=x_le,
                y_le=y_le,
                z_le=z_le,
                chord=chord,
                twist=twist,
                airfoil=airfoil,
                control_surface_type=control_surface_type,
                control_surface_hinge_point=control_surface_hinge_point,
                control_surface_deflection=control_surface_deflection,
                num_spanwise_panels=num_spanwise_panels,
                spanwise_spacing=spanwise_spacing,
            )

            # Add this new object to the list of wing cross sections.
            wing_cross_sections.append(this_wing_cross_section)

        # Return the list of wing cross sections.
        return wing_cross_sections

    # NOTE: I haven't yet started refactoring this method.
    def get_max_period(self):
        """This method returns the longest period of this movement object's cycles.

        :return max_period: float
            The longest period in seconds.
        """

        max_period = max(
            self.sweeping_period, self.pitching_period, self.heaving_period
        )

        return max_period


# NOTE: I haven't yet started refactoring this class.
class OperatingPointMovement:
    """This is a class used to contain the movement characteristics of an operating
    point.

    This class contains the following public methods:
        generate_operating_points: This method creates the operating point objects at
        each time current_step, and groups them into a list.

        get_max_period: This method returns the longest period of this movement
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
            velocity_list = oscillating_sinspace(
                amplitude=self.velocity_amplitude,
                period=self.velocity_period,
                base_value=self.velocity_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.velocity_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            velocity_list = oscillating_linspace(
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
    def get_max_period(self):
        """This method returns the longest period of this movement object's cycles.

        :return max_period: float
            The longest period in seconds.
        """

        max_period = self.velocity_period

        return max_period


# NOTE: I haven't yet started refactoring this function.
def oscillating_sinspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D array of values that are calculated by inputting a
    vector of linearly spaced time steps into a sine function.

    :param amplitude: float
        This is the amplitude of the value fluctuation.
    :param period: float
        This is the period of the value fluctuation.
    :param base_value: float
        This is the starting value.
    :param num_steps: int
        This is the number of time steps to iterate through.
    :param delta_time: float
        This is the change in time between each time step.
    :return: 1D array of floats
        This is the resulting vector of sinusoidally spaced values
    """
    # If either the amplitude or the period are 0, return a vector with length equal
    # to the number of steps, and all the values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = 0
    k = base_value

    # Calculate and return the values.
    return a * np.sin(b * (times - h)) + k


# NOTE: I haven't yet started refactoring this function.
def oscillating_linspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D array of values that are calculated by inputting a
    vector of linearly spaced time steps into a triangle function.

    :param amplitude: float
        This is the amplitude of the value fluctuation.
    :param period: float
        This is the period of the value fluctuation.
    :param base_value: float
        This is the starting value.
    :param num_steps: int
        This is the number of time steps to iterate through.
    :param delta_time: float
        This is the change in time between each time step.
    :return: 1D array of floats
        This is the resulting vector of uniformly spaced values
    """
    # If either the amplitude or the period are 0, return a vector with length equal
    # to the number of steps, and all the values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = np.pi / 2
    k = base_value

    # Calculate and return the values.
    return a * signal.sawtooth((b * times + h), 0.5) + k


# NOTE: I haven't yet started refactoring this function.
def oscillating_customspace(
    amplitude, period, base_value, num_steps, delta_time, custom_function
):
    """This function returns a 1D array of values that are calculated by inputting a
    vector of linearly spaced time steps into a custom function.

    :param amplitude: float
        This is the amplitude of the value fluctuation.
    :param period: float
        This is the period of the value fluctuation.
    :param base_value: float
        This is the starting value.
    :param num_steps: int
        This is the number of time steps to iterate through.
    :param delta_time: float
        This is the change in time between each time step.
    :param custom_function: function
        This is a custom function used to return the values. For example, it could be
        np.cos or np.sinh (assuming numpy had previously been imported as np). It
        will be horizontally scaled by the period, vertically scaled by the
        amplitude. For example, say the function has an internal amplitude of 2
        units, an internal period of 3 units, amplitude is set to 4 units and period
        is set to 5 units. The result will have a net amplitude of 8 units and a net
        period of 15 units.
    :return: 1D array of floats
        This is the resulting vector of custom spaced values
    """
    # If either the amplitude or the period are 0, return a vector with length equal
    # to the number of steps, and all the values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps, endpoint=False)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = 0
    k = base_value

    # Calculate and return the values.
    return a * custom_function(b * (times - h)) + k
