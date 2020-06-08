
"""This module contains the class definitions for the problem's movement and the problem's operating point.

This module contains the following classes:
    Movement: This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.
    AirplaneMovement: This is a class used to contain the movement characteristics of an airplane.
    WingMovement: This is a class used to contain the movement characteristics of a wing.
    WingCrossSectionMovement: This is a class used to contain the movement characteristics of a wing cross section.
    OperatingPoint: This is a class used to contain the problem's operating point characteristics.
    OperatingPointMovement: This is a class used to contain the movement characteristics of an operating point.

This module contains the following exceptions:
    None

This module contains the following functions:
    oscillating_sinspace: This function returns a 1D ndarray of values that are calculated by inputting a vector of
                          linearly spaced time steps into a sine function.
    oscillating_linspace: This function returns a 1D ndarray of values that are calculated by inputting a vector of
                          linearly spaced time steps into a triangle function.
"""

from scipy import signal
import numpy as np
import aviansoftwareminimumviableproduct as asmvp


class Movement:
    """This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.

    This class contains the following public methods:
        get_flapping_velocity_at_point_on_panel: This method gets the velocity due to flapping at a point on the panel
                                                 of a given airplane based the time step, its current position, and its
                                                 last position.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, airplane_movement, operating_point_movement, num_steps=10, delta_time=0.1):
        """This is the initialization method.
        
        :param airplane_movement: AirplaneMovement
            This object characterizes the movement of the airplane.
        :param operating_point_movement: OperatingPointMovement
            This object characterizes the movement of the the operating point.
        :param num_steps: int, optional
            This integer is the number of time steps of the unsteady simulation. Its default value is 10.
        :param delta_time: float, optional
            This float is the time, in seconds, between each time step. Its default value is 0.1 seconds.
        """

        # Initialize the class attributes.
        self.airplane_movement = airplane_movement
        self.operating_point_movement = operating_point_movement
        self.num_steps = num_steps
        self.delta_time = delta_time

        # Generate a list of the airplanes and operating points that are the steps through this movement object.
        self.airplanes = self.airplane_movement.generate_airplanes(
            num_steps=self.num_steps,
            delta_time=self.delta_time
        )
        self.operating_points = self.operating_point_movement.generate_operating_points(
            num_steps=self.num_steps,
            delta_time=self.delta_time
        )

    def get_flapping_velocity_at_point_on_panel(self, wing_position, panel_chordwise_position, panel_spanwise_position,
                                                point_name, current_step):
        """This method gets the velocity due to flapping at a point on the panel of a given airplane based the time
        step, its current position, and its last position.

        :param wing_position: int
            This is the position of the panel's wing in the airplane's list of wings.
        :param panel_chordwise_position: int
            This is the chordwise position of the panel in the wing's ndarray of panels.
        :param panel_spanwise_position: int
            This is the spanwise position of the panel in the wing's ndarray of panels.
        :param point_name: string
            This is the name of the point at which to find the velocity due to flapping. It can be 'center',
            'collocation', 'front right vertex', 'front left vertex', 'back left vertex', 'back right vertex',
            'ring vortex front right vertex', 'ring vortex front left vertex', 'ring vortex back left vertex',
            'ring vortex back right vertex', 'ring vortex front leg center', 'ring vortex left leg center',
            'ring vortex back leg center', or 'ring vortex right leg center'.
        :param current_step: int
            This is the step number at which to evaluate the flapping velocity at the point.
        :return flapping_velocity: 1D ndarray
            This is the flapping velocity at the current time step at the given point. It is a (,3) ndarray with units
            of meters per second.
        """

        # If this is the first time step, there is no previous geometry, so the flapping velocity is (0, 0, 0)
        # meters per second.
        if current_step < 1:
            flapping_velocity = np.zeros(3)
            return flapping_velocity

        # Find the current, and the last airplane.
        airplane = self.airplanes[current_step]
        last_airplane = self.airplanes[current_step - 1]

        # Find the current, and the last wing.
        wing = airplane.wings[wing_position]
        last_wing = last_airplane.wings[wing_position]

        # Find the current, and the last panel.
        panel = wing.panels[panel_chordwise_position, panel_spanwise_position]
        last_panel = last_wing.panels[panel_chordwise_position, panel_spanwise_position]

        # Initialize variables to hold the current and the last position.
        position = None
        last_position = None

        # Populate the current and last position based on the point's name.
        if point_name == 'center':
            position = panel.center
            last_position = last_panel.center
        elif point_name == 'collocation':
            position = panel.collocation_point
            last_position = last_panel.collocation_point
        elif point_name == 'front right vertex':
            position = panel.front_right_vertex
            last_position = last_panel.front_right_vertex
        elif point_name == 'front left vertex':
            position = panel.front_left_vertex
            last_position = last_panel.front_left_vertex
        elif point_name == 'back left vertex':
            position = panel.back_left_vertex
            last_position = last_panel.back_left_vertex
        elif point_name == 'back right vertex':
            position = panel.back_right_vertex
            last_position = last_panel.back_right_vertex
        elif point_name == 'ring vortex front right vertex':
            position = panel.ring_vortex.front_right_vertex
            last_position = last_panel.ring_vortex.front_right_vertex
        elif point_name == 'ring vortex front left vertex':
            position = panel.ring_vortex.back_left_vertex
            last_position = last_panel.ring_vortex.back_left_vertex
        elif point_name == 'ring vortex back left vertex':
            position = panel.ring_vortex.back_left_vertex
            last_position = last_panel.ring_vortex.back_left_vertex
        elif point_name == 'ring vortex back right vertex':
            position = panel.ring_vortex.back_right_vertex
            last_position = last_panel.ring_vortex.back_right_vertex
        elif point_name == 'ring vortex front leg center':
            position = panel.ring_vortex.front_leg.center
            last_position = last_panel.ring_vortex.front_leg.center
        elif point_name == 'ring vortex left leg center':
            position = panel.ring_vortex.left_leg.center
            last_position = last_panel.ring_vortex.left_leg.center
        elif point_name == 'ring vortex back leg center':
            position = panel.ring_vortex.back_leg.center
            last_position = last_panel.ring_vortex.back_leg.center
        elif point_name == 'ring vortex right leg center':
            position = panel.ring_vortex.right_leg.center
            last_position = last_panel.ring_vortex.right_leg.center

        # Calculate and return the flapping velocity.
        flapping_velocity = (position - last_position) / self.delta_time
        return flapping_velocity


class AirplaneMovement:
    """This is a class used to contain the movement characteristics of an airplane.

    This class contains the following public methods:
        generate_airplanes: This method creates the airplane object at each time step, and groups them into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # ToDo: Properly document this method.
    def __init__(self, base_airplane, wing_movements, x_ref_amplitude=0.0, x_ref_period=0.0, x_ref_spacing='sine',
                 y_ref_amplitude=0.0, y_ref_period=0.0, y_ref_spacing='sine', z_ref_amplitude=0.0, z_ref_period=0.0,
                 z_ref_spacing='sine'):
        """
        
        :param base_airplane: 
        :param wing_movements: 
        :param x_ref_amplitude: 
        :param x_ref_period: 
        :param x_ref_spacing: 
        :param y_ref_amplitude: 
        :param y_ref_period: 
        :param y_ref_spacing: 
        :param z_ref_amplitude: 
        :param z_ref_period: 
        :param z_ref_spacing: 
        """

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

    # ToDo: Properly document this method.
    def generate_airplanes(self, num_steps=10, delta_time=0.1):
        """This method creates the airplane object at each time step, and groups them into a list.
        
        :param num_steps: 
        :param delta_time: 
        :return: 
        """

        # Create an ndarray of x_ref points.
        if self.x_ref_spacing == 'sine':
            x_ref_list = oscillating_sinspace(
                amplitude=self.x_ref_amplitude,
                period=self.x_ref_period,
                base_value=self.x_ref_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.x_ref_spacing == 'uniform':
            x_ref_list = oscillating_linspace(
                amplitude=self.x_ref_amplitude,
                period=self.x_ref_period,
                base_value=self.x_ref_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of x_ref_spacing!")

        # Create an ndarray of y_ref points.
        if self.y_ref_spacing == 'sine':
            y_ref_list = oscillating_sinspace(
                amplitude=self.y_ref_amplitude,
                period=self.y_ref_period,
                base_value=self.y_ref_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.y_ref_spacing == 'uniform':
            y_ref_list = oscillating_linspace(
                amplitude=self.y_ref_amplitude,
                period=self.y_ref_period,
                base_value=self.y_ref_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of y_ref_spacing!")

        # Create an ndarray of z_ref points.
        if self.z_ref_spacing == 'sine':
            z_ref_list = oscillating_sinspace(
                amplitude=self.z_ref_amplitude,
                period=self.z_ref_period,
                base_value=self.z_ref_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.z_ref_spacing == 'uniform':
            z_ref_list = oscillating_linspace(
                amplitude=self.z_ref_amplitude,
                period=self.z_ref_period,
                base_value=self.z_ref_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of z_ref_spacing!")

        wings = np.empty((len(self.wing_movements), num_steps), dtype=object)

        for wing_movement_location in range(len(self.wing_movements)):
            wing_movement = self.wing_movements[wing_movement_location]
            this_wings_list_of_wings = np.array(wing_movement.generate_wings(
                num_steps=num_steps, delta_time=delta_time))
            wings[wing_movement_location, :] = this_wings_list_of_wings

        # Create an empty list of airplanes.
        airplanes = []

        name = self.base_airplane.name

        for step in range(num_steps):
            x_ref = x_ref_list[step]
            y_ref = y_ref_list[step]
            z_ref = z_ref_list[step]
            these_wings = wings[:, step]

            this_airplane = asmvp.geometry.Airplane(
                name=name,
                x_ref=x_ref,
                y_ref=y_ref,
                z_ref=z_ref,
                wings=these_wings
            )

            airplanes.append(this_airplane)

        # Return the list of airplanes.
        return airplanes


class WingMovement:
    """This is a class used to contain the movement characteristics of a wing.

    This class contains the following public methods:
        generate_wings: This method creates the wing object at each time step, and groups them into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # ToDo: Properly document this method.
    def __init__(self, base_wing, wing_cross_sections_movements, x_le_amplitude=0.0, x_le_period=0.0,
                 x_le_spacing='sine', y_le_amplitude=0.0, y_le_period=0.0, y_le_spacing='sine', z_le_amplitude=0.0,
                 z_le_period=0.0, z_le_spacing='sine'):
        """

        :param base_wing:
        :param wing_cross_sections_movements:
        :param x_le_amplitude:
        :param x_le_period:
        :param x_le_spacing:
        :param y_le_amplitude:
        :param y_le_period:
        :param y_le_spacing:
        :param z_le_amplitude:
        :param z_le_period:
        :param z_le_spacing:
        """

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

    # ToDo: Properly document this method.
    def generate_wings(self, num_steps=10, delta_time=0.1):
        """This method creates the wing object at each time step, and groups them into a list.

        :param num_steps:
        :param delta_time:
        :return:
        """

        # Create an ndarray of x_le points.
        if self.x_le_spacing == 'sine':
            x_le_list = oscillating_sinspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.x_le_spacing == 'uniform':
            x_le_list = oscillating_linspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of x_le_spacing!")

        # Create an ndarray of y_le points.
        if self.y_le_spacing == 'sine':
            y_le_list = oscillating_sinspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.y_le_spacing == 'uniform':
            y_le_list = oscillating_linspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of y_le_spacing!")

        # Create an ndarray of z_le points.
        if self.z_le_spacing == 'sine':
            z_le_list = oscillating_sinspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.z_le_spacing == 'uniform':
            z_le_list = oscillating_linspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of z_le_spacing!")

        wing_cross_sections = np.empty((len(self.wing_cross_section_movements), num_steps), dtype=object)

        for wing_cross_section_movement_location in range(len(self.wing_cross_section_movements)):
            wing_cross_section_movement = self.wing_cross_section_movements[wing_cross_section_movement_location]
            x = np.array(wing_cross_section_movement.generate_wing_cross_sections(
                num_steps=num_steps, delta_time=delta_time))
            wing_cross_sections[wing_cross_section_movement_location, :] = x

        # Create an empty list of wings.
        wings = []

        name = self.base_wing.name
        symmetric = self.base_wing.symmetric
        num_chordwise_panels = self.base_wing.num_chordwise_panels
        chordwise_spacing = self.base_wing.chordwise_spacing

        for step in range(num_steps):
            x_le = x_le_list[step]
            y_le = y_le_list[step]
            z_le = z_le_list[step]
            cross_sections = wing_cross_sections[:, step]

            this_wing = asmvp.geometry.Wing(
                name=name,
                x_le=x_le,
                y_le=y_le,
                z_le=z_le,
                wing_cross_sections=cross_sections,
                symmetric=symmetric,
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing
            )

            wings.append(this_wing)

        # Return the list of wings.
        return wings


class WingCrossSectionMovement:
    """This is a class used to contain the movement characteristics of a wing cross section.

    This class contains the following public methods:
        generate_wing_cross_sections: This method creates the wing cross section objects at each time step, and groups
                                      them into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # ToDo: Properly document this method.
    def __init__(self, base_wing_cross_section, x_le_amplitude=0.0, x_le_period=0.0, x_le_spacing='sine',
                 y_le_amplitude=0.0, y_le_period=0.0, y_le_spacing='sine', z_le_amplitude=0.0, z_le_period=0.0,
                 z_le_spacing='sine', twist_amplitude=0.0, twist_period=0.0, twist_spacing='sine',
                 control_surface_deflection_amplitude=0.0, control_surface_deflection_period=0.0,
                 control_surface_deflection_spacing='sine'):
        """

        :param base_wing_cross_section:
        :param x_le_amplitude:
        :param x_le_period:
        :param x_le_spacing:
        :param y_le_amplitude:
        :param y_le_period:
        :param y_le_spacing:
        :param z_le_amplitude:
        :param z_le_period:
        :param z_le_spacing:
        :param twist_amplitude:
        :param twist_period:
        :param twist_spacing:
        :param control_surface_deflection_amplitude:
        :param control_surface_deflection_period:
        :param control_surface_deflection_spacing:
        """

        self.base_wing_cross_section = base_wing_cross_section

        self.x_le_base = self.base_wing_cross_section.x_le
        self.x_le_amplitude = x_le_amplitude
        self.x_le_period = x_le_period
        self.x_le_spacing = x_le_spacing

        self.y_le_base = self.base_wing_cross_section.y_le
        self.y_le_amplitude = y_le_amplitude
        self.y_le_period = y_le_period
        self.y_le_spacing = y_le_spacing

        self.z_le_base = self.base_wing_cross_section.z_le
        self.z_le_amplitude = z_le_amplitude
        self.z_le_period = z_le_period
        self.z_le_spacing = z_le_spacing

        self.twist_base = self.base_wing_cross_section.twist
        self.twist_amplitude = twist_amplitude
        self.twist_period = twist_period
        self.twist_spacing = twist_spacing

        self.control_surface_deflection_base = self.base_wing_cross_section.control_surface_deflection
        self.control_surface_deflection_amplitude = control_surface_deflection_amplitude
        self.control_surface_deflection_period = control_surface_deflection_period
        self.control_surface_deflection_spacing = control_surface_deflection_spacing

    # ToDo: Properly document this method.
    def generate_wing_cross_sections(self, num_steps=10, delta_time=0.1):
        """This method creates the wing cross section objects at each time step, and groups them into a list.

        :param num_steps:
        :param delta_time:
        :return:
        """

        # Create an ndarray of x_le points.
        if self.x_le_spacing == 'sine':
            x_le_list = oscillating_sinspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.x_le_spacing == 'uniform':
            x_le_list = oscillating_linspace(
                amplitude=self.x_le_amplitude,
                period=self.x_le_period,
                base_value=self.x_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of x_le_spacing!")

        # Create an ndarray of y_le points.
        if self.y_le_spacing == 'sine':
            y_le_list = oscillating_sinspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.y_le_spacing == 'uniform':
            y_le_list = oscillating_linspace(
                amplitude=self.y_le_amplitude,
                period=self.y_le_period,
                base_value=self.y_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of y_le_spacing!")

        # Create an ndarray of z_le points.
        if self.z_le_spacing == 'sine':
            z_le_list = oscillating_sinspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.z_le_spacing == 'uniform':
            z_le_list = oscillating_linspace(
                amplitude=self.z_le_amplitude,
                period=self.z_le_period,
                base_value=self.z_le_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of z_le_spacing!")

        # Create an ndarray of twist values.
        if self.twist_spacing == 'sine':
            twist_list = oscillating_sinspace(
                amplitude=self.twist_amplitude,
                period=self.twist_period,
                base_value=self.twist_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.twist_spacing == 'uniform':
            twist_list = oscillating_linspace(
                amplitude=self.twist_amplitude,
                period=self.twist_period,
                base_value=self.twist_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of twist_spacing!")

        # Create an ndarray of control surface deflection values.
        if self.control_surface_deflection_spacing == 'sine':
            control_surface_deflection_list = oscillating_sinspace(
                amplitude=self.control_surface_deflection_amplitude,
                period=self.control_surface_deflection_period,
                base_value=self.control_surface_deflection_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.control_surface_deflection_spacing == 'uniform':
            control_surface_deflection_list = oscillating_linspace(
                amplitude=self.control_surface_deflection_amplitude,
                period=self.control_surface_deflection_period,
                base_value=self.control_surface_deflection_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of control_surface_deflection_spacing!")

        # Create an empty list of wing cross sections.
        wing_cross_sections = []

        chord = self.base_wing_cross_section.chord
        airfoil = self.base_wing_cross_section.airfoil
        control_surface_type = self.base_wing_cross_section.control_surface_type
        control_surface_hinge_point = self.base_wing_cross_section.control_surface_hinge_point
        num_spanwise_panels = self.base_wing_cross_section.num_spanwise_panels
        spanwise_spacing = self.base_wing_cross_section.spanwise_spacing

        for step in range(num_steps):
            x_le = x_le_list[step]
            y_le = y_le_list[step]
            z_le = z_le_list[step]
            twist = twist_list[step]
            control_surface_deflection = control_surface_deflection_list[step]

            this_wing_cross_section = asmvp.geometry.WingCrossSection(
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
                spanwise_spacing=spanwise_spacing
            )

            wing_cross_sections.append(this_wing_cross_section)

        # Return the list of wing cross sections.
        return wing_cross_sections


class OperatingPoint:
    """This is a class used to contain the problem's operating point characteristics.

    Citation:
        Adapted from:         performance.OperatingPoint in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/29/2020

    This class contains the following public methods:
        calculate_dynamic_pressure: This method calculates the freestream dynamic pressure of the working fluid.
        calculate_rotation_matrix_wind_to_geometry: This method computes the 3 x 3 rotation matrix for converting from
                                                    wind axes to geometry axes.
        calculate_freestream_direction_geometry_axes: This method computes the freestream direction (the direction the
                                                      wind is going to) in geometry axes.
        calculate_freestream_velocity_geometry_axes: This method computes the freestream velocity vector (in the
                                                     direction the wind is going to) in geometry axes.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, density=1.225, velocity=10.0, alpha=5.0, beta=0.0):
        """This is the initialization method.

        :param density: float, optional
            This parameter is the density. The units are kilograms per meters cubed. The default value is 1.225.
        :param velocity: float, optional
            This parameter is the freestream speed in the positive x direction. The units are meters per second. The
            default value is 10.0.
        :param alpha: float, optional
            This parameter is the angle of attack. The units are degrees. The default value is 5.0.
        :param beta: float, optional
            This parameter is the side slip angle. The units are degrees. The default value is 0.0.
        """

        # Initialize the attributes.
        self.density = density
        self.velocity = velocity
        self.alpha = alpha
        self.beta = beta

    def calculate_dynamic_pressure(self):
        """This method calculates the freestream dynamic pressure of the working fluid.

        :return dynamic_pressure: float
            This is the freestream dynamic pressure. Its units are pascals.
        """

        # Calculate and return the freestream dynamic pressure
        dynamic_pressure = 0.5 * self.density * self.velocity ** 2
        return dynamic_pressure

    def calculate_rotation_matrix_wind_axes_to_geometry_axes(self):
        """This method computes the 3 x 3 rotation matrix for converting from wind axes to geometry axes.

        :return rotation_matrix_wind_axes_to_geometry_axes: 3 x 3 ndarray
            This is the rotation matrix to convert wind axes to geometry axes.
        """

        sin_alpha = np.sin(np.radians(self.alpha))
        cos_alpha = np.cos(np.radians(self.alpha))
        sin_beta = np.sin(np.radians(self.beta))
        cos_beta = np.cos(np.radians(self.beta))
        eye = np.eye(3)

        alpha_rotation = np.array([
            [cos_alpha, 0, -sin_alpha],
            [0, 1, 0],
            [sin_alpha, 0, cos_alpha]
        ])
        beta_rotation = np.array([
            [cos_beta, -sin_beta, 0],
            [sin_beta, cos_beta, 0],
            [0, 0, 1]
        ])

        # Flip the axes because in geometry axes x is downstream by convention, while in wind axes x is upstream by
        # convention. Same with z being up/down respectively.
        axes_flip = np.array([
            [-1, 0, 0],
            [0, 1, 0, ],
            [0, 0, -1]
        ])

        # Calculate and return the rotation matrix to convert wind axes to geometry axes.
        rotation_matrix_wind_axes_to_geometry_axes = axes_flip @ alpha_rotation @ beta_rotation @ eye
        return rotation_matrix_wind_axes_to_geometry_axes

    def calculate_freestream_direction_geometry_axes(self):
        """This method computes the freestream direction (the direction the wind is going to) in geometry axes.

        :return velocity_direction_geometry_axes: 1D ndarray
            This is the freestream velocity direction in geometry axes.
        """

        velocity_direction_wind_axes = np.array([-1, 0, 0])
        velocity_direction_geometry_axes = (self.calculate_rotation_matrix_wind_axes_to_geometry_axes()
                                            @ velocity_direction_wind_axes)
        return velocity_direction_geometry_axes

    def calculate_freestream_velocity_geometry_axes(self):
        """This method computes the freestream velocity vector (in the direction the wind is going to) in geometry axes.

        :return freestream_velocity_geometry_axes: 1D ndarray
            This is the freestream velocity vector in geometry axes.
        """

        freestream_velocity_geometry_axes = self.calculate_freestream_direction_geometry_axes() * self.velocity
        return freestream_velocity_geometry_axes


class OperatingPointMovement:
    """This is a class used to contain the movement characteristics of an operating point.

    This class contains the following public methods:
        generate_operating_points: This method creates the operating point objects at each time step, and groups them
                                   into a list.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    # ToDo: Properly document this method.
    def __init__(self, base_operating_point, velocity_amplitude=0.0, velocity_period=0.0, velocity_spacing='sine'):
        """
        
        :param base_operating_point: 
        :param velocity_amplitude: 
        :param velocity_period: 
        :param velocity_spacing:
        """

        self.base_operating_point = base_operating_point

        self.velocity_base = self.base_operating_point.velocity
        self.velocity_amplitude = velocity_amplitude
        self.velocity_period = velocity_period
        self.velocity_spacing = velocity_spacing

    # ToDo: Properly document this method.
    def generate_operating_points(self, num_steps=10, delta_time=0.1):
        """This method creates the operating point objects at each time step, and groups them into a list.
        
        :param num_steps: 
        :param delta_time: 
        :return: 
        """

        # Create an ndarray of velocities.
        if self.velocity_spacing == 'sine':
            velocity_list = oscillating_sinspace(
                amplitude=self.velocity_amplitude,
                period=self.velocity_period,
                base_value=self.velocity_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.velocity_spacing == 'uniform':
            velocity_list = oscillating_linspace(
                amplitude=self.velocity_amplitude,
                period=self.velocity_period,
                base_value=self.velocity_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of velocity_spacing!")

        # Create an empty list of operating points.
        operating_points = []

        density = self.base_operating_point.density
        alpha = self.base_operating_point.alpha
        beta = self.base_operating_point.beta

        for step in range(num_steps):
            velocity = velocity_list[step]

            this_operating_point = OperatingPoint(
                density=density,
                velocity=velocity,
                alpha=alpha,
                beta=beta
            )

            operating_points.append(this_operating_point)

        # Return the list of operating points.
        return operating_points


# ToDo: Properly document this function.
def oscillating_sinspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D ndarray of values that are calculated by inputting a vector of linearly spaced time
    steps into a sine function.

    :param amplitude:
    :param period:
    :param base_value:
    :param num_steps:
    :param delta_time:
    :return:
    """

    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    total_time = num_steps * delta_time

    times = np.linspace(0, total_time, num_steps)

    a = amplitude
    b = 2 * np.pi / period
    h = 0
    k = base_value

    values = a * np.sin(b * (times - h)) + k

    return values


# ToDo: Properly document this function.
def oscillating_linspace(amplitude, period, base_value, num_steps, delta_time):
    """This function returns a 1D ndarray of values that are calculated by inputting a vector of linearly spaced time
    steps into a triangle function.

    :param amplitude:
    :param period:
    :param base_value:
    :param num_steps:
    :param delta_time:
    :return:
    """

    if amplitude == 0 or period == 0:
        return np.ones(num_steps) * base_value

    total_time = num_steps * delta_time

    times = np.linspace(0, total_time, num_steps)

    a = amplitude
    b = 2 * np.pi / period
    h = np.pi / 2
    k = base_value

    values = a * signal.sawtooth((b * times + h), 0.5) + k

    return values
