# ToDo: Properly document this module.
"""This module contains the class definitions for the geometry's movement and the problem's operating point.

This module contains the following classes:
    Movement: This is a class used to contain the movement characteristics of an unsteady aerodynamics problem.
    OperatingPoint: This is a class used to contain the problem's operating point characteristics.

This module contains the following exceptions:
    None

This module contains the following functions:
    oscillating_sinspace:
    oscillating_linspace:
"""

from scipy import signal
import numpy as np
import aviansoftwareminimumviableproduct as asmvp


# ToDo: Properly document this class.
class Movement:
    """
    
    """
    
    # ToDo: Properly document this method.
    def __init__(self, airplane_movement, operating_point_movement, num_steps=10, delta_time=0.1):
        """
        
        :param airplane_movement: 
        :param operating_point_movement: 
        :param num_steps: 
        :param delta_time: 
        """
        
        self.airplane_movement = airplane_movement
        self.operating_point_movement = operating_point_movement
        self.num_steps = num_steps
        self.delta_time = delta_time

        self.airplanes = self.airplane_movement.generate_airplanes(
            num_steps=self.num_steps,
            delta_time=self.delta_time
        )
        self.operating_points = self.operating_point_movement.generate_operating_points(
            num_steps=self.num_steps,
            delta_time=self.delta_time
        )


# ToDo: Properly document this class.
class AirplaneMovement:
    """
    
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
        """
        
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

        wings = np.empty((len(self.wing_movements), num_steps))

        for wing_movement_location in range(len(self.wing_movements)):
            wing_movement = self.wing_movements[wing_movement_location]
            wings[wing_movement_location, :] = (
                wing_movement.generate_cross_sections(num_steps=num_steps, delta_time=delta_time))

        # Create an ndarray of airplanes.
        airplanes = np.empty(num_steps)

        name = self.base_airplane.name

        for step in range(num_steps):
            x_ref = x_ref_list[step]
            y_ref = y_ref_list[step]
            z_ref = z_ref_list[step]
            wings = wings[:, step]

            airplanes[step] = asmvp.geometry.Airplane(
                name=name,
                x_ref=x_ref,
                y_ref=y_ref,
                z_ref=z_ref,
                wings=wings
            )

        # Return the ndarray of airplanes.
        return airplanes


# ToDo: Properly document this class.
class WingMovement:
    """

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
        """

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

        wing_cross_sections = np.empty((len(self.wing_cross_section_movements), num_steps))

        for wing_cross_section_movement_location in range(len(self.wing_cross_section_movements)):
            wing_cross_section_movement = self.wing_cross_section_movements[wing_cross_section_movement_location]
            wing_cross_sections[wing_cross_section_movement_location, :] = (
                wing_cross_section_movement.generate_cross_sections(num_steps=num_steps, delta_time=delta_time))

        # Create an ndarray of wings.
        wings = np.empty(num_steps)

        name = self.base_wing.name
        symmetric = self.base_wing.symmetric
        num_chordwise_panels = self.base_wing.num_chordwise_panels
        chordwise_spacing = self.base_wing.chordwise_spacing

        for step in range(num_steps):
            x_le = x_le_list[step]
            y_le = y_le_list[step]
            z_le = z_le_list[step]
            cross_sections = wing_cross_sections[:, step]

            wings[step] = asmvp.geometry.Wing(
                name=name,
                x_le=x_le,
                y_le=y_le,
                z_le=z_le,
                cross_sections=cross_sections,
                symmetric=symmetric,
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing
            )

        # Return the ndarray of wing cross sections.
        return wings


# ToDo: Properly document this class.
class WingCrossSectionMovement:
    """

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
        """

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

        # Create an ndarray of wing cross sections.
        wing_cross_sections = np.empty(num_steps)

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

            wing_cross_sections[step] = asmvp.geometry.WingCrossSection(
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

        # Return the ndarray of wing cross sections.
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


# ToDo: Properly document this class.
class OperatingPointMovement:
    """
    
    """
    
    # ToDo: Properly document this method.
    def __init__(self, base_operating_point, velocity_amplitude=0.0, velocity_period=0.0, velocity_spacing='sine',
                 alpha_amplitude=0.0, alpha_period=0.0, alpha_spacing='sine', beta_amplitude=0.0, beta_period=0.0,
                 beta_spacing='sine'):
        """
        
        :param base_operating_point: 
        :param velocity_amplitude: 
        :param velocity_period: 
        :param velocity_spacing: 
        :param alpha_amplitude: 
        :param alpha_period: 
        :param alpha_spacing: 
        :param beta_amplitude: 
        :param beta_period: 
        :param beta_spacing: 
        """

        self.base_operating_point = base_operating_point

        self.velocity_base = self.base_operating_point.velocity
        self.velocity_amplitude = velocity_amplitude
        self.velocity_period = velocity_period
        self.velocity_spacing = velocity_spacing

        self.alpha_base = self.base_operating_point.alpha
        self.alpha_amplitude = alpha_amplitude
        self.alpha_period = alpha_period
        self.alpha_spacing = alpha_spacing

        self.beta_base = self.base_operating_point.beta
        self.beta_amplitude = beta_amplitude
        self.beta_period = beta_period
        self.beta_spacing = beta_spacing

    # ToDo: Properly document this method.
    def generate_operating_points(self, num_steps=10, delta_time=0.1):
        """
        
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

        # Create an ndarray of alpha values.
        if self.alpha_spacing == 'sine':
            alpha_list = oscillating_sinspace(
                amplitude=self.alpha_amplitude,
                period=self.alpha_period,
                base_value=self.alpha_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.alpha_spacing == 'uniform':
            alpha_list = oscillating_linspace(
                amplitude=self.alpha_amplitude,
                period=self.alpha_period,
                base_value=self.alpha_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of alpha_spacing!")

        # Create an ndarray of beta values.
        if self.beta_spacing == 'sine':
            beta_list = oscillating_sinspace(
                amplitude=self.beta_amplitude,
                period=self.beta_period,
                base_value=self.beta_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        elif self.beta_spacing == 'uniform':
            beta_list = oscillating_linspace(
                amplitude=self.beta_amplitude,
                period=self.beta_period,
                base_value=self.beta_base,
                num_steps=num_steps,
                delta_time=delta_time
            )
        else:
            raise Exception("Bad value of beta_spacing!")

        # Create an ndarray of operating points.
        operating_points = np.empty(num_steps)

        density = self.base_operating_point.density

        for step in range(num_steps):
            velocity = velocity_list[step]
            alpha = alpha_list[step]
            beta = beta_list[step]

            operating_points[step] = OperatingPoint(
                density=density,
                velocity=velocity,
                alpha=alpha,
                beta=beta
            )

        # Return the ndarray of operating points.
        return operating_points


# ToDo: Properly document this function.
def oscillating_sinspace(amplitude, period, base_value, num_steps, delta_time):
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
    total_time = num_steps * delta_time

    times = np.linspace(0, total_time, num_steps)

    a = amplitude
    b = 2 * np.pi / period
    h = np.pi / 2
    k = base_value

    values = a * signal.sawtooth((b * times + h), 0.5) + k

    return values
