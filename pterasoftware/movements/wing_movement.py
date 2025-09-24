# NOTE: I haven't yet started refactoring this module.
"""This module contains the WingMovement class.

This module contains the following classes:
    WingMovement: This is a class used to contain the Wing movements.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from . import functions

from .. import geometry


# NOTE: I haven't yet started refactoring this class.
class WingMovement:
    """This is a class used to contain the movement characteristics of a wing.

    This class contains the following public methods:
        generate_wings: This method creates the wing object at each time step,
        and groups them into a list.

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
            x_le_list = functions.oscillating_sinspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.x_le_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            x_le_list = functions.oscillating_linspace(
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
            y_le_list = functions.oscillating_sinspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.y_le_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            y_le_list = functions.oscillating_linspace(
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
            z_le_list = functions.oscillating_sinspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.z_le_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            z_le_list = functions.oscillating_linspace(
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
    @property
    def max_period(self):
        """This method returns the longest period of this movement object's
        sub-movement objects, sub-sub-movement objects, etc.

        :return max_period: float
            The longest period in seconds.
        """

        wing_cross_section_movement_max_periods = []
        for wing_cross_section_movement in self.wing_cross_section_movements:
            wing_cross_section_movement_max_periods.append(
                wing_cross_section_movement.max_period
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
