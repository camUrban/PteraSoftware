"""This module contains the class definitions for the problem's movement.

This module contains the following classes:
    Movement: This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.
    AirplaneMovement: This is a class used to contain the movement characteristics of an current_airplane.
    WingMovement: This is a class used to contain the movement characteristics of a wing.
    WingCrossSectionMovement: This is a class used to contain the movement characteristics of a wing cross section.
    OperatingPointMovement: This is a class used to contain the movement characteristics of an operating point.

This module contains the following exceptions:
    None

This module contains the following functions:
    oscillating_sinspace: This function returns a 1D ndarray of values that are calculated by inputting a vector of
                          linearly spaced time steps into a sine function.
    oscillating_linspace: This function returns a 1D ndarray of values that are calculated by inputting a vector of
                          linearly spaced time steps into a triangle function.
"""

import numpy as np
from scipy import signal

import pterasoftware as ps


class Movement:
    """This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.

    This class contains the following public methods:
        None

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self, airplane_movement, operating_point_movement, num_steps=10, delta_time=0.1
    ):
        """This is the initialization method.
        
        :param airplane_movement: AirplaneMovement
            This object characterizes the movement of the current_airplane.
        :param operating_point_movement: OperatingPointMovement
            This object characterizes the movement of the the operating point.
        :param num_steps: int, optional
            This integer is the number of time steps of the unsteady simulation. Its default value is 10.
        :param delta_time: float, optional
            This float is the time, in seconds, between each time current_step. Its default value is 0.1 seconds.
        """

        # Initialize the class attributes.
        self.num_steps = num_steps
        self.delta_time = delta_time

        # Generate a list of the airplanes and operating points that are the steps through this movement object.
        self.airplanes = airplane_movement.generate_airplanes(
            num_steps=self.num_steps, delta_time=self.delta_time
        )
        self.operating_points = operating_point_movement.generate_operating_points(
            num_steps=self.num_steps, delta_time=self.delta_time
        )


class AirplaneMovement:
    """This is a class used to contain the movement characteristics of an current_airplane.

    This class contains the following public methods:
        generate_airplanes: This method creates the current_airplane object at each time current_step, and groups them
                            into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

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
        """ This is the initialization method.
        
        :param base_airplane: Airplane
            This is the first airplane object, from which the others will be created.
        :param wing_movements: list of WingMovement objects
            This is a list of the WingMovement objects associated with each of the base airplane's wings.
        :param x_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its x reference point. Its units are meters and its
            default value is 0 meters.
        :param x_ref_period: float, optional
            This is the period of the airplane's change in its x reference point. Its units are seconds and its
            default value is 0 seconds.
        :param x_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its x reference point. The options are "sine",
            and "uniform". The default value is "sine".
        :param y_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its y reference point. Its units are meters and its
            default value is 0 meters.
        :param y_ref_period: float, optional
            This is the period of the airplane's change in its y reference point. Its units are seconds and its
            default value is 0 seconds.
        :param y_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its y reference point. The options are "sine",
            and "uniform". The default value is "sine".
        :param z_ref_amplitude: float, optional
            This is the amplitude of the airplane's change in its z reference point. Its units are meters and its
            default value is 0 meters.
        :param z_ref_period: float, optional
            This is the period of the airplane's change in its z reference point. Its units are seconds and its
            default value is 0 seconds.
        :param z_ref_spacing: string, optional
            This value determines the spacing of the airplane's change in its z reference point. The options are "sine",
            and "uniform". The default value is "sine".
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

    def generate_airplanes(self, num_steps=10, delta_time=0.1):
        """This method creates the current_airplane object at each time current_step, and groups them into a list.
        
        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value is 0.1 seconds.
        :return airplanes: list of Airplane objects
            This is the list of Airplane objects that is associated with this AirplaneMovement object.
        """

        # Check the x_ref spacing value.
        if self.x_ref_spacing == "sine":

            # Create an ndarray of points with a sinusoidal spacing.
            x_ref_list = oscillating_sinspace(
                amplitude=self.x_ref_amplitude,
                period=self.x_ref_period,
                base_value=self.x_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.x_ref_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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

            # Create an ndarray of points with a sinusoidal spacing.
            y_ref_list = oscillating_sinspace(
                amplitude=self.y_ref_amplitude,
                period=self.y_ref_period,
                base_value=self.y_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.y_ref_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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

            # Create an ndarray of points with a sinusoidal spacing.
            z_ref_list = oscillating_sinspace(
                amplitude=self.z_ref_amplitude,
                period=self.z_ref_period,
                base_value=self.z_ref_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.z_ref_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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

        # Create an empty ndarray that will hold each of the airplane's wing's vector of other wing's based its
        # movement.
        wings = np.empty((len(self.wing_movements), num_steps), dtype=object)

        # Iterate through the wing movement locations.
        for wing_movement_location in range(len(self.wing_movements)):

            # Get the wing movement.
            wing_movement = self.wing_movements[wing_movement_location]

            # Generate this wing's vector of other wing's based on its movement.
            this_wings_list_of_wings = np.array(
                wing_movement.generate_wings(num_steps=num_steps, delta_time=delta_time)
            )

            # Add this vector the airplane's ndarray of wing objects.
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
            this_airplane = ps.geometry.Airplane(
                name=name, x_ref=x_ref, y_ref=y_ref, z_ref=z_ref, wings=these_wings
            )

            # Add this new object to the list of airplanes.
            airplanes.append(this_airplane)

        # Return the list of airplanes.
        return airplanes


class WingMovement:
    """This is a class used to contain the movement characteristics of a wing.

    This class contains the following public methods:
        generate_wings: This method creates the wing object at each time current_step, and groups them into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

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
        """ This is the initialization method.

        :param base_wing: Wing
            This is the first wing object, from which the others will be created.
        :param wing_cross_sections_movements: list of WingCrossSectionMovement objects
            This is a list of the WingCrossSectionMovement objects associated with each of the base wing's cross
            sections.
        :param x_le_amplitude: float, optional
            This is the amplitude of the wing's change in its x reference point. Its units are meters and its
            default value is 0 meters.
        :param x_le_period: float, optional
            This is the period of the wing's change in its x reference point. Its units are seconds and its
            default value is 0 seconds.
        :param x_le_spacing: string, optional
            This value determines the spacing of the wing's change in its x reference point. The options are "sine",
            and "uniform". The default value is "sine".
        :param y_le_amplitude: float, optional
            This is the amplitude of the wing's change in its y reference point. Its units are meters and its
            default value is 0 meters.
        :param y_le_period: float, optional
            This is the period of the wing's change in its y reference point. Its units are seconds and its
            default value is 0 seconds.
        :param y_le_spacing: string, optional
            This value determines the spacing of the wing's change in its y reference point. The options are "sine",
            and "uniform". The default value is "sine".
        :param z_le_amplitude: float, optional
            This is the amplitude of the wing's change in its z reference point. Its units are meters and its
            default value is 0 meters.
        :param z_le_period: float, optional
            This is the period of the wing's change in its z reference point. Its units are seconds and its
            default value is 0 seconds.
        :param z_le_spacing: string, optional
            This value determines the spacing of the wing's change in its z reference point. The options are "sine",
            and "uniform". The default value is "sine".
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

    def generate_wings(self, num_steps=10, delta_time=0.1):
        """This method creates the wing object at each time current_step, and groups them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value is 0.1 seconds.
        :return wings: list of Wing objects
            This is the list of Wing objects that is associated with this WingMovement object.
        """

        # Check the x_le spacing value.
        if self.x_le_spacing == "sine":

            # Create an ndarray of points with a sinusoidal spacing.
            x_le_list = oscillating_sinspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.x_le_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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

            # Create an ndarray of points with a sinusoidal spacing.
            y_le_list = oscillating_sinspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.y_le_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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

            # Create an ndarray of points with a sinusoidal spacing.
            z_le_list = oscillating_sinspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.z_le_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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

        # Create an empty ndarray that will hold each of the wing's wing cross section's vector of other wing cross
        # section's based its movement.
        wing_cross_sections = np.empty(
            (len(self.wing_cross_section_movements), num_steps), dtype=object
        )

        # Iterate through the wing cross section movement locations.
        for wing_cross_section_movement_location in range(
            len(self.wing_cross_section_movements)
        ):

            # Get the wing cross section movement.
            wing_cross_section_movement = self.wing_cross_section_movements[
                wing_cross_section_movement_location
            ]

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

                cross_section_span = 0.0
                cross_section_sweep = 0.0
                cross_section_pitch = 0.0
                cross_section_heave = 0.0
                last_x_le = 0.0
                last_y_le = 0.0
                last_z_le = 0.0

            else:
                last_wing_cross_section_movement = self.wing_cross_section_movements[
                    wing_cross_section_movement_location - 1
                ]

                this_base_wing_cross_section = (
                    wing_cross_section_movement.base_wing_cross_section
                )
                last_base_wing_cross_section = (
                    last_wing_cross_section_movement.base_wing_cross_section
                )

                this_x_le = this_base_wing_cross_section.x_le
                this_y_le = this_base_wing_cross_section.y_le
                this_z_le = this_base_wing_cross_section.z_le
                this_twist = this_base_wing_cross_section.twist

                last_x_le = last_base_wing_cross_section.x_le
                last_y_le = last_base_wing_cross_section.y_le
                last_z_le = last_base_wing_cross_section.z_le
                last_twist = last_base_wing_cross_section.twist

                cross_section_span = np.sqrt(
                    (this_x_le - last_x_le) ** 2
                    + (this_y_le - last_y_le) ** 2
                    + (this_z_le - last_z_le) ** 2
                )

                try:
                    cross_section_sweep = (
                        np.arctan((this_z_le - last_z_le) / (this_y_le - last_y_le))
                        * 180
                        / np.pi
                    )
                except ZeroDivisionError:
                    cross_section_sweep = 0.0
                    wing_is_vertical = True

                cross_section_pitch = this_twist - last_twist

                try:
                    cross_section_heave = (
                        np.arctan((this_x_le - last_x_le) / (this_y_le - last_y_le))
                        * 180
                        / np.pi
                    )
                except ZeroDivisionError:
                    cross_section_heave = 0.0
                    wing_is_vertical = True

            # Generate this wing cross section's vector of other wing cross section's based on its movement.
            this_wing_cross_sections_list_of_wing_cross_sections = np.array(
                wing_cross_section_movement.generate_wing_cross_sections(
                    num_steps=num_steps,
                    delta_time=delta_time,
                    cross_section_span=cross_section_span,
                    cross_section_sweep=cross_section_sweep,
                    cross_section_pitch=cross_section_pitch,
                    cross_section_heave=cross_section_heave,
                    last_x_le=last_x_le,
                    last_y_le=last_y_le,
                    last_z_le=last_z_le,
                    wing_is_vertical=wing_is_vertical,
                )
            )

            # Add this vector the wing's ndarray of wing cross section objects.
            wing_cross_sections[
                wing_cross_section_movement_location, :
            ] = this_wing_cross_sections_list_of_wing_cross_sections

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
            this_wing = ps.geometry.Wing(
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


class WingCrossSectionMovement:
    """ This is a class used to contain the movement characteristics of a wing cross section.

    This class contains the following public methods:
        generate_wing_cross_sections: This method creates the wing cross section objects at each time current_step, and
                                      groups them into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        base_wing_cross_section,
        sweeping_amplitude=0.0,
        sweeping_period=0.0,
        sweeping_spacing="sine",
        pitching_amplitude=0.0,
        pitching_period=0.0,
        pitching_spacing="sine",
        heaving_amplitude=0.0,
        heaving_period=0.0,
        heaving_spacing="sine",
    ):
        """ This is the initialization method.

        :param base_wing_cross_section: WingCrossSection
            This is the first wing cross section object, from which the others will be created.
        :param sweeping_amplitude: float, optional
            This is the amplitude of the cross section's change in its sweep, relative to the previous cross section.
            Its units are degrees and its default value is 0.0 degrees.
        :param sweeping_period: float, optional
            This is the period of the cross section's change in its sweep. Its units are seconds and its default value
            is 0.0 seconds.
        :param sweeping_spacing: string, optional
            This value determines the spacing of the cross section's change in its sweep. The options are "sine", and
            "uniform". The default value is "sine".
        :param pitching_amplitude: float, optional
            This is the amplitude of the cross section's change in its pitch, relative to the previous cross section.
            Its units are degrees and its default value is 0.0 degrees.
        :param pitching_period: float, optional
            This is the period of the cross section's change in its pitch. Its units are seconds and its default value
            is 0.0 seconds.
        :param pitching_spacing: string, optional
            This value determines the spacing of the cross section's change in its pitch. The options are "sine", and
            "uniform". The default value is "sine".
        :param heaving_amplitude: float, optional
            This is the amplitude of the cross section's change in its heave, relative to the previous cross section.
            Its units are degrees and its default value is 0.0 degrees.
        :param heaving_period: float, optional
            This is the period of the cross section's change in its heave. Its units are seconds and its default value
            is 0.0 seconds.
        :param heaving_spacing: string, optional
            This value determines the spacing of the cross section's change in its heave. The options are "sine", and
            "uniform". The default value is "sine".
        """

        # Initialize the class attributes.
        self.base_wing_cross_section = base_wing_cross_section
        self.sweeping_amplitude = sweeping_amplitude
        self.sweeping_period = sweeping_period
        self.sweeping_spacing = sweeping_spacing
        self.sweeping_base = 0.0
        self.pitching_amplitude = pitching_amplitude
        self.pitching_period = pitching_period
        self.pitching_spacing = pitching_spacing
        self.pitching_base = self.base_wing_cross_section.twist
        self.heaving_amplitude = heaving_amplitude
        self.heaving_period = heaving_period
        self.heaving_spacing = heaving_spacing
        self.heaving_base = 0.0
        self.x_le_base = self.base_wing_cross_section.x_le
        self.y_le_base = self.base_wing_cross_section.y_le
        self.z_le_base = self.base_wing_cross_section.z_le
        self.twist_base = self.base_wing_cross_section.twist
        self.control_surface_deflection_base = (
            self.base_wing_cross_section.control_surface_deflection
        )

    def generate_wing_cross_sections(
        self,
        num_steps=10,
        delta_time=0.1,
        cross_section_span=0.0,
        cross_section_sweep=0.0,
        cross_section_pitch=0.0,
        cross_section_heave=0.0,
        last_x_le=0.0,
        last_y_le=0.0,
        last_z_le=0.0,
        wing_is_vertical=False,
    ):
        """This method creates the wing cross section objects at each time current_step, and groups them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value is 0.1 seconds.
        :param cross_section_span: float, optional
            This is the length, in meters, of the leading edge stretching between this cross section at the previous
            cross section. If this is the first cross section, it should be 0.0 meters. The default value is 0.0 meters.
        :param last_x_le: float, optional
            This is the x coordinate of the reference location of the previous cross section. Its units are in meters,
            and its default value is 0.0 meters.
        :param last_y_le: float, optional
            This is the y coordinate of the reference location of the previous cross section. Its units are in meters,
            and its default value is 0.0 meters.
        :param last_z_le: float, optional
            This is the z coordinate of the reference location of the previous cross section. Its units are in meters,
            and its default value is 0.0 meters.
        :return wing_cross_sections: list of WingCrossSection objects
            This is the list of WingCrossSection objects that is associated with this WingCrossSectionMovement object.
        """

        # Check the sweeping spacing value.
        if self.sweeping_spacing == "sine":

            # Create an ndarray of points with a sinusoidal spacing.
            sweeping_list = oscillating_sinspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=cross_section_sweep,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.sweeping_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
            sweeping_list = oscillating_linspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=cross_section_sweep,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of sweeping_spacing!")

        # Check the pitching spacing value.
        if self.pitching_spacing == "sine":

            # Create an ndarray of points with a sinusoidal spacing.
            pitching_list = oscillating_sinspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.pitching_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
            pitching_list = oscillating_linspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of pitching_spacing!")

        # Check the heaving spacing value.
        if self.heaving_spacing == "sine":

            # Create an ndarray of points with a sinusoidal spacing.
            heaving_list = oscillating_sinspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=cross_section_heave,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.heaving_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
            heaving_list = oscillating_linspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=cross_section_heave,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        else:

            # Throw an exception if the spacing value is not "sine" or "uniform".
            raise Exception("Bad value of heaving_spacing!")

        if wing_is_vertical:
            print("Wing is vertical!")

            x_le_list = np.ones(num_steps) * self.x_le_base
            y_le_list = np.ones(num_steps) * self.y_le_base
            z_le_list = np.ones(num_steps) * self.z_le_base
            twist_list = np.ones(num_steps) * self.twist_base
        else:

            # Find the list of new leading edge points. This uses a spherical coordinate transformation, referencing the
            # previous wing cross section's leading edge point as the origin. Also convert the lists of sweep, pitch, and
            # heave values to radians before passing them into numpy's trigonometry functions.
            x_le_list = last_x_le + cross_section_span * np.cos(
                sweeping_list * np.pi / 180
            ) * np.sin(heaving_list * np.pi / 180)
            y_le_list = last_y_le + cross_section_span * np.cos(
                sweeping_list * np.pi / 180
            ) * np.cos(heaving_list * np.pi / 180)
            z_le_list = last_z_le + cross_section_span * np.sin(
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
            this_wing_cross_section = ps.geometry.WingCrossSection(
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


class OperatingPointMovement:
    """This is a class used to contain the movement characteristics of an operating point.

    This class contains the following public methods:
        generate_operating_points: This method creates the operating point objects at each time current_step, and groups
                                   them into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        base_operating_point,
        velocity_amplitude=0.0,
        velocity_period=0.0,
        velocity_spacing="sine",
    ):
        """ This is the initialization method.

        :param base_operating_point: OperatingPoint
            This is the operating point object, from which the others will be created.
        :param velocity_amplitude: float, optional
            This is the amplitude of the operating point's change in velocity. Its units are meters per second and its
            default value is 0 meters per second.
        :param velocity_period: float, optional
            This is the period of the operating point's change in its velocity. Its units are seconds and its
            default value is 0 seconds.
        :param velocity_spacing: string, optional
            This value determines the spacing of the operating point's change in its velocity. The options are "sine",
            and "uniform". The default value is "sine".
        """

        # Initialize the class attributes.
        self.base_operating_point = base_operating_point
        self.velocity_base = self.base_operating_point.velocity
        self.velocity_amplitude = velocity_amplitude
        self.velocity_period = velocity_period
        self.velocity_spacing = velocity_spacing

    def generate_operating_points(self, num_steps=10, delta_time=0.1):
        """This method creates the operating point objects at each time current_step, and groups them into a list.

        :param num_steps: int, optional
            This is the number of time steps in this movement. The default value is 10.
        :param delta_time: float, optional
            This is the time, in seconds, between each time step. The default value is 0.1 seconds.
        :return operating_points: list of OperatingPoint objects
            This is the list of OperatingPoint objects that is associated with this OperatingPointMovement object.
        """

        # Check the velocity spacing value.
        if self.velocity_spacing == "sine":

            # Create an ndarray of points with a sinusoidal spacing.
            velocity_list = oscillating_sinspace(
                amplitude=self.velocity_amplitude,
                period=self.velocity_period,
                base_value=self.velocity_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.velocity_spacing == "uniform":

            # Create an ndarray of points with a uniform spacing.
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
            this_operating_point = ps.operating_point.OperatingPoint(
                density=density, velocity=velocity, alpha=alpha, beta=beta
            )

            # Add this new object to the list of operating points.
            operating_points.append(this_operating_point)

        # Return the list of operating points.
        return operating_points


def oscillating_sinspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D ndarray of values that are calculated by inputting a vector of linearly spaced time
    steps into a sine function.

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
    :return values: 1D ndarray of floats
        This is the resulting vector of sinusoidally spaced values
    """

    # If either the amplitude or the period are 0, return a vector with length equal to the number of steps, and all the
    # values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = 0
    k = base_value

    # Calculate and return the values.
    values = a * np.sin(b * (times - h)) + k
    return values


def oscillating_linspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D ndarray of values that are calculated by inputting a vector of linearly spaced time
    steps into a triangle function.

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
    :return values: 1D ndarray of floats
        This is the resulting vector of uniformly spaced values
    """

    # If either the amplitude or the period are 0, return a vector with length equal to the number of steps, and all the
    # values equal to the base value.
    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    # Calculate the total time.
    total_time = num_steps * delta_time

    # Get the time at each time step.
    times = np.linspace(0, total_time, num_steps)

    # Convert the function characteristics into classic wave function constants.
    a = amplitude
    b = 2 * np.pi / period
    h = np.pi / 2
    k = base_value

    # Calculate and return the values.
    values = a * signal.sawtooth((b * times + h), 0.5) + k
    return values
