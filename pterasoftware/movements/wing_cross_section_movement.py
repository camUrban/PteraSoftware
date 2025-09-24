# NOTE: I haven't yet started refactoring this module.
"""This module contains the WingCrossSectionMovement class.

This module contains the following classes:
    WingCrossSectionMovement: This is a class used to contain the WingCrossSection
    movements.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from . import functions

from .. import geometry


# NOTE: I haven't yet started refactoring this class.
class WingCrossSectionMovement:
    """This is a class used to contain the movement characteristics of a wing cross
    section.

    This class contains the following public methods:
        generate_wing_cross_sections: This method creates the wing cross section
        objects at each time current_step, and groups them into a list.

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
            sweeping_list = functions.oscillating_sinspace(
                amplitude=self.sweeping_amplitude,
                period=self.sweeping_period,
                base_value=cross_section_sweep,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.sweeping_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            sweeping_list = functions.oscillating_linspace(
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
            sweeping_list = functions.oscillating_customspace(
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
            pitching_list = functions.oscillating_sinspace(
                amplitude=self.pitching_amplitude,
                period=self.pitching_period,
                base_value=self.pitching_base,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.pitching_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            pitching_list = functions.oscillating_linspace(
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
            pitching_list = functions.oscillating_customspace(
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
            heaving_list = functions.oscillating_sinspace(
                amplitude=self.heaving_amplitude,
                period=self.heaving_period,
                base_value=cross_section_heave,
                num_steps=num_steps,
                delta_time=delta_time,
            )
        elif self.heaving_spacing == "uniform":

            # Create an array of points with a uniform spacing.
            heaving_list = functions.oscillating_linspace(
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
            heaving_list = functions.oscillating_customspace(
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
    @property
    def max_period(self):
        """This method returns the longest period of this movement object's cycles.

        :return max_period: float
            The longest period in seconds.
        """

        max_period = max(
            self.sweeping_period, self.pitching_period, self.heaving_period
        )

        return max_period
