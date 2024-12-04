"""This module contains the Panel class.

This module contains the following classes:
    Panel: This class is used to contain the panels of a wing.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np


# ToDo: Update the list of methods for this class.
class Panel:
    """This class is used to contain the panels of a wing.

    This class contains the following public methods:
        calculate_collocation_point_location: This method calculates the location of
        the collocation point.

        calculate_area_and_normal: This method calculates the panel's area and the
        panel's normal unit vector.

        calculate_normalized_induced_velocity: This method calculates the velocity
        induced at a point by this panel's vortices, assuming a unit vortex strength.

        calculate_induced_velocity: This method calculates the velocity induced at a
        point by this panel's vortices with their given vortex strengths.

        update_coefficients: This method updates the panel's force coefficients.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(self, front_right_vertex, front_left_vertex, back_left_vertex,
            back_right_vertex, is_leading_edge, is_trailing_edge, ):
        """This is the initialization method.

        :param front_right_vertex: 1D array with three elements
            This is an array containing the x, y, and z coordinates of the panel's
            front right vertex.
        :param front_left_vertex: 1D array with three elements
            This is an array containing the x, y, and z coordinates of the panel's
            front left vertex.
        :param back_left_vertex: 1D array with three elements
            This is an array containing the x, y, and z coordinates of the panel's
            back left vertex.
        :param back_right_vertex: 1D array with three elements
            This is an array containing the x, y, and z coordinates of the panel's
            back right vertex.
        :param is_leading_edge: bool
            This is true if the panel is a leading edge panel on a wing, and false
            otherwise.
        :param is_trailing_edge: bool
            This is true if the panel is a trailing edge panel on a wing, and false
            otherwise.
        """

        # Initialize the attributes.
        self.front_right_vertex = front_right_vertex
        self.front_left_vertex = front_left_vertex
        self.back_left_vertex = back_left_vertex
        self.back_right_vertex = back_right_vertex
        self.is_leading_edge = is_leading_edge
        self.is_trailing_edge = is_trailing_edge

        # Initialize variables to hold attributes that describe the panel's position
        # in its wing's panel matrix. They will be populated by the meshing function.
        self.is_right_edge = None
        self.is_left_edge = None
        self.local_chordwise_position = None
        self.local_spanwise_position = None

        # Initialize variables to hold the panel's ring and horseshoe vortices. These
        # will be populated by the solver.
        self.ring_vortex = None
        self.horseshoe_vortex = None

        # Initialize variables to hold attributes of the panel that will be defined
        # after the solver finds a solution.
        self.near_field_force_geometry_axes = None
        self.near_field_moment_geometry_axes = None
        self.near_field_force_wind_axes = None
        self.near_field_moment_wind_axes = None
        self.induced_drag_coefficient = None
        self.side_force_coefficient = None
        self.lift_coefficient = None

    # ToDo: Update this method's documentation.
    @property
    def right_leg(self):
        return self.front_right_vertex - self.back_right_vertex

    # ToDo: Update this method's documentation.
    @property
    def front_leg(self):
        return self.front_left_vertex - self.front_right_vertex

    # ToDo: Update this method's documentation.
    @property
    def left_leg(self):
        return self.back_left_vertex - self.front_left_vertex

    # ToDo: Update this method's documentation.
    @property
    def back_leg(self):
        return self.back_right_vertex - self.back_left_vertex

    # ToDo: Update this method's documentation.
    @property
    def front_right_vortex_vertex(self):
        return self.back_right_vertex + 0.75 * self.right_leg

    # ToDo: Update this method's documentation.
    @property
    def front_left_vortex_vertex(self):
        return self.front_left_vertex + 0.25 * self.left_leg

    # ToDo: Update this method's documentation.
    @property
    def collocation_point(self):
        # Find the location of points three quarters of the way down the left and
        # right legs of the panel.
        right_three_quarter_chord_mark = self.back_right_vertex + 0.25 * self.right_leg
        left_three_quarter_chord_mark = self.front_left_vertex + 0.75 * self.left_leg

        # Find the vector between the points three quarters of the way down the left
        # and right legs of the panel.
        three_quarter_chord_vector = (
                left_three_quarter_chord_mark - right_three_quarter_chord_mark)

        # Find the collocation point, which is halfway between the points three
        # quarters of the way down the left and right legs of the panel. Then
        # populate the class attribute.
        return right_three_quarter_chord_mark + 0.5 * three_quarter_chord_vector

    # ToDo: Update this method's documentation.
    @property
    def area(self):
        return np.linalg.norm(self._cross) / 2

    # ToDo: Update this method's documentation.
    @property
    def unit_normal(self):
        return self._cross / np.linalg.norm(self._cross)

    # ToDo: Update this method's documentation.
    @property
    def unit_spanwise(self):
        front_spanwise = -self.front_leg
        back_spanwise = self.back_leg
        spanwise = (front_spanwise + back_spanwise) / 2
        return spanwise / np.linalg.norm(spanwise)

    # ToDo: Update this method's documentation.
    @property
    def average_span(self):
        front_leg_length = np.linalg.norm(self.front_leg)
        back_leg_length = np.linalg.norm(self.back_leg)
        return (front_leg_length + back_leg_length) / 2

    # ToDo: Update this method's documentation.
    @property
    def unit_chordwise(self):
        right_chordwise = -self.right_leg
        left_chordwise = self.left_leg
        chordwise = (right_chordwise + left_chordwise) / 2
        return chordwise / np.linalg.norm(chordwise)

    # ToDo: Update this method's documentation.
    @property
    def average_chord(self):
        right_leg_length = np.linalg.norm(self.right_leg)
        left_leg_length = np.linalg.norm(self.left_leg)
        return (right_leg_length + left_leg_length) / 2

    # ToDo: Update this method's documentation.
    @property
    def _first_diagonal(self):
        return self.front_right_vertex - self.back_left_vertex

    # ToDo: Update this method's documentation.
    @property
    def _second_diagonal(self):
        return self.front_left_vertex - self.back_right_vertex

    # ToDo: Update this method's documentation.
    @property
    def _cross(self):
        return np.cross(self._first_diagonal, self._second_diagonal)

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this panel's
        vortices, assuming a unit vortex strength.

        This method does not include the effect of the panel's wake vortices.

        :param point:  1D array
            This is a vector containing the x, y, and z coordinates of the point to
            find the induced velocity at.
        :return: 1D array
            This is a vector containing the x, y, and z components of the induced
            velocity.
        """

        normalized_induced_velocity = np.zeros(3)

        if self.ring_vortex is not None:
            normalized_induced_velocity += (
                self.ring_vortex.calculate_normalized_induced_velocity(point=point))
        if self.horseshoe_vortex is not None:
            normalized_induced_velocity += (
                self.horseshoe_vortex.calculate_normalized_induced_velocity(
                    point=point))

        return normalized_induced_velocity

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this panel's
        vortices with their given vortex strengths.

        This method does not include the effect of the panel's wake vortices.

        :param point: 1D array
            This is a vector containing the x, y, and z coordinates of the point to
            find the induced velocity at.
        :return: 1D array
            This is a vector containing the x, y, and z components of the induced
            velocity.
        """

        induced_velocity = np.zeros(3)

        if self.ring_vortex is not None:
            induced_velocity += self.ring_vortex.calculate_induced_velocity(point=point)
        if self.horseshoe_vortex is not None:
            induced_velocity += self.horseshoe_vortex.calculate_induced_velocity(
                point=point)

        return induced_velocity

    def update_coefficients(self, dynamic_pressure):
        """This method updates the panel's force coefficients.

        :return: None
        """

        induced_drag = -self.near_field_force_wind_axes[0]
        side_force = self.near_field_force_wind_axes[1]
        lift = -self.near_field_force_wind_axes[2]

        self.induced_drag_coefficient = induced_drag / self.area / dynamic_pressure
        self.side_force_coefficient = side_force / self.area / dynamic_pressure
        self.lift_coefficient = lift / self.area / dynamic_pressure
