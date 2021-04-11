"""This module contains the class definition for the problem's operating point.

This module contains the following classes:
    OperatingPoint: This is a class used to contain the problem's operating point
    characteristics.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np


class OperatingPoint:
    """This is a class used to contain the problem's operating point characteristics.

    Citation:
        Adapted from:         performance.OperatingPoint in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/29/2020

    This class contains the following public methods:
        calculate_dynamic_pressure: This method calculates the freestream dynamic
        pressure of the working fluid.
        calculate_rotation_matrix_wind_to_geometry: This method computes the 3 x 3
        rotation matrix for converting from
                                                    wind axes to geometry axes.
        calculate_freestream_direction_geometry_axes: This method computes the
        freestream direction (the direction the
                                                      wind is going to) in geometry
                                                      axes.
        calculate_freestream_velocity_geometry_axes: This method computes the
        freestream velocity vector (in the
                                                     direction the wind is going to)
                                                     in geometry axes.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, density=1.225, velocity=10.0, alpha=5.0, beta=0.0):
        """This is the initialization method.

        :param density: float, optional
            This parameter is the density. The units are kilograms per meters cubed.
            The default value is 1.225.
        :param velocity: float, optional
            This parameter is the freestream speed in the positive x direction. The
            units are meters per second. The
            default value is 10.0.
        :param alpha: float, optional
            This parameter is the angle of attack. The units are degrees. The default
            value is 5.0.
        :param beta: float, optional
            This parameter is the sideslip angle. The units are degrees. The default
            value is 0.0.
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
        """This method computes the 3 x 3 rotation matrix for converting from wind
        axes to geometry axes.

        :return rotation_matrix_wind_axes_to_geometry_axes: 3 x 3 array
            This is the rotation matrix to convert wind axes to geometry axes.
        """

        sin_alpha = np.sin(np.radians(self.alpha))
        cos_alpha = np.cos(np.radians(self.alpha))
        sin_beta = np.sin(np.radians(self.beta))
        cos_beta = np.cos(np.radians(self.beta))
        eye = np.eye(3)

        alpha_rotation = np.array(
            [[cos_alpha, 0, -sin_alpha], [0, 1, 0], [sin_alpha, 0, cos_alpha]]
        )
        beta_rotation = np.array(
            [[cos_beta, -sin_beta, 0], [sin_beta, cos_beta, 0], [0, 0, 1]]
        )

        # Flip the axes because in geometry axes x is downstream by convention,
        # while in wind axes x is upstream by
        # convention. Same with z being up/down respectively.
        axes_flip = np.array(
            [
                [-1, 0, 0],
                [
                    0,
                    1,
                    0,
                ],
                [0, 0, -1],
            ]
        )

        # Calculate and return the rotation matrix to convert wind axes to geometry
        # axes.
        rotation_matrix_wind_axes_to_geometry_axes = (
            axes_flip @ alpha_rotation @ beta_rotation @ eye
        )
        return rotation_matrix_wind_axes_to_geometry_axes

    def calculate_freestream_direction_geometry_axes(self):
        """This method computes the freestream direction (the direction the wind is
        going to) in geometry axes.

        :return velocity_direction_geometry_axes: 1D array
            This is the freestream velocity direction in geometry axes.
        """

        velocity_direction_wind_axes = np.array([-1, 0, 0])
        velocity_direction_geometry_axes = (
            self.calculate_rotation_matrix_wind_axes_to_geometry_axes()
            @ velocity_direction_wind_axes
        )
        return velocity_direction_geometry_axes

    def calculate_freestream_velocity_geometry_axes(self):
        """This method computes the freestream velocity vector (in the direction the
        wind is going to) in geometry axes.

        :return freestream_velocity_geometry_axes: 1D array
            This is the freestream velocity vector in geometry axes.
        """

        freestream_velocity_geometry_axes = (
            self.calculate_freestream_direction_geometry_axes() * self.velocity
        )
        return freestream_velocity_geometry_axes
