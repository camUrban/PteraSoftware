# ToDo: Properly document this module.
import numpy as np

from . import functions


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

        update_force_moment_and_pressure: This method updates the force, moment,
        and pressure on this panel.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        front_right_vertex,
        front_left_vertex,
        back_left_vertex,
        back_right_vertex,
        is_leading_edge,
        is_trailing_edge,
    ):
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
        # in its wing's panel matrix. They
        # will be populated by the meshing function.
        self.is_right_edge = None
        self.is_left_edge = None
        self.local_chordwise_position = None
        self.local_spanwise_position = None

        # Initialize a variable to hold the position of the wing (who this panel
        # belongs to) in the airplane's wing
        # list. This will be populated by the airplane object.
        self.wing_position = None

        # Initialize variables to hold the panel's ring and horseshoe vortices. These
        # will be populated by the solver.
        self.ring_vortex = None
        self.horseshoe_vortex = None

        # Initialize a variable to hold the collocation point location and then
        # populate it.
        self.collocation_point = None
        self.calculate_collocation_point_location()

        # Initialize variables to hold the panel area and the panel normal vector at
        # the collocation point. Then
        # populate them.
        self.area = None
        self.normal_direction = None
        self.calculate_area_and_normal()

        # Calculate the center of the panel.
        self.center = functions.numba_centroid_of_quadrilateral(
            front_right_vertex, front_left_vertex, back_left_vertex, back_right_vertex
        )

        # Calculate the front and back leg lengths, then use them to find and
        # populate the average panel width.
        front_leg_length = np.linalg.norm(front_left_vertex - front_right_vertex)
        back_leg_length = np.linalg.norm(back_right_vertex - back_left_vertex)
        self.width = (front_leg_length + back_leg_length) / 2

        # Initialize two variables that are along the panel's left and right legs at
        # the quarter chord. These points
        # are used for all types of solvers, so we will define them here.
        self.front_right_vortex_vertex = self.front_right_vertex + 0.25 * (
            self.back_right_vertex - self.front_right_vertex
        )
        self.front_left_vortex_vertex = self.front_left_vertex + 0.25 * (
            self.back_left_vertex - self.front_left_vertex
        )

        # Initialize variables to hold attributes of the panel that will be defined
        # after the solver finds a solution.
        self.near_field_force_geometry_axes = None
        self.near_field_moment_geometry_axes = None
        self.delta_pressure = None

    def calculate_collocation_point_location(self):
        """This method calculates the location of the collocation point.

        The collocation point is at the panel's three quarter chord point.

        :return: None
        """

        # Find the location of points three quarters of the way down the left and
        # right legs of the panel.
        right_three_quarter_chord_mark = self.front_right_vertex + 0.75 * (
            self.back_right_vertex - self.front_right_vertex
        )
        left_three_quarter_chord_mark = self.front_left_vertex + 0.75 * (
            self.back_left_vertex - self.front_left_vertex
        )

        # Find the vector between the points three quarters of the way down the left
        # and right legs of the panel.
        three_quarter_chord_vector = (
            left_three_quarter_chord_mark - right_three_quarter_chord_mark
        )

        # Find the collocation point, which is halfway between the points three
        # quarters of the way down the left and
        # right legs of the panel. Then populate the class attribute.
        self.collocation_point = (
            right_three_quarter_chord_mark + 0.5 * three_quarter_chord_vector
        )

    def calculate_area_and_normal(self):
        """This method calculates the panel's area and the panel's normal unit vector.

        This method makes the assumption that the panel is planar. This is
        technically incorrect for wing's with twist
        but is a good approximation for small panels.

        :return: None
        """

        # Calculate panel's normal unit vector and its area via its diagonals.
        first_diagonal = self.front_right_vertex - self.back_left_vertex
        second_diagonal = self.front_left_vertex - self.back_right_vertex
        cross_product = np.cross(first_diagonal, second_diagonal)
        cross_product_magnitude = np.linalg.norm(cross_product)
        self.normal_direction = cross_product / cross_product_magnitude
        self.area = cross_product_magnitude / 2

    def calculate_normalized_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this panel's
        vortices, assuming a unit vortex
        strength.

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
                self.ring_vortex.calculate_normalized_induced_velocity(point=point)
            )
        if self.horseshoe_vortex is not None:
            normalized_induced_velocity += (
                self.horseshoe_vortex.calculate_normalized_induced_velocity(point=point)
            )

        return normalized_induced_velocity

    def calculate_induced_velocity(self, point):
        """This method calculates the velocity induced at a point by this panel's
        vortices with their given vortex
        strengths.

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
                point=point
            )

        return induced_velocity

    def update_pressure(self):
        """This method updates the pressure across this panel.

        :return: None
        """

        self.delta_pressure = (
            np.dot(self.near_field_force_geometry_axes, self.normal_direction)
            / self.area
        )
