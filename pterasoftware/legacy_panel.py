"""This module contains the Panel class.

This module contains the following classes:
    Panel: This class is used to contain the panels of a wing.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np


class Panel:
    """This class is used to contain the panels of a wing.

    This class contains the following public methods:
        rightLeg_G: This method defines a property for the panel's right leg vector as
        a (3,) array.

        frontLeg_G: This method defines a property for the panel's front leg vector as
        a (3,) array.

        leftLeg_G: This method defines a property for the panel's left leg vector as a
        (3,) array.

        backLeg_G: This method defines a property for the panel's back leg vector as a
        (3,) array.

        Frbvp_G_Cg: This method defines a property for the coordinates
        of the front-right vertex of the ring vortex as a (3,) array.

        Flbvp_G_Cg: This method defines a property for the coordinates
        of the front-left vertex of the ring vortex as a (3,) array.

        Cpp_G_Cg: This method defines a property for the coordinates of the
        panel's collocation point as a (3,) array.

        area: This method defines a property which is an estimate of the panel's area.

        unitNormal_G: This method defines a property for an estimate of the panel's
        unit normal vector as a (3,) array.

        unitSpanwise_G: This method defines a property for the panel's unit spanwise
        vector as a ( 3,) array.

        unitChordwise_G: This method defines a property for the panel's unit chordwise
        vector as a (3,) array.

        average_span: This method defines a property for the average span of the panel.

        average_chord: This method defines a property for the average chord of the
        panel.

        calculate_normalized_induced_velocity: This method calculates the velocity
        induced at a point by this panel's vortices, assuming a unit vortex strength.

        calculate_induced_velocity: This method calculates the velocity induced at a
        point by this panel's vortices with their given vortex strengths.

        calculate_projected_area: This method calculates the area of the panel
        projected on some plane defined by its unit normal vector.

        update_coefficients: This method updates the panel's force coefficients.

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

        :param front_right_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            front right vertex.
        :param front_left_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            front left vertex.
        :param back_left_vertex: (3,) array
            This is an array containing the x, y, and z coordinates of the panel's
            back left vertex.
        :param back_right_vertex: (3,) array
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

    @property
    def right_leg(self):
        """This method defines a property for the panel's right leg vector as a (3,
        ) array.

        :return: (3,) array of floats
            This is the panel's right leg vector, which is defined from back to
            front. The units are in meters.
        """
        return self.front_right_vertex - self.back_right_vertex

    @property
    def front_leg(self):
        """This method defines a property for the panel's front leg vector as a (3,
        ) array.

        :return: (3,) array of floats
            This is the panel's front leg vector, which is defined from right to
            left. The units are in meters.
        """
        return self.front_left_vertex - self.front_right_vertex

    @property
    def left_leg(self):
        """This method defines a property for the panel's left leg vector as a (3,
        ) array.

        :return: (3,) array of floats
            This is the panel's left leg vector, which is defined from front to
            back. The units are in meters.
        """
        return self.back_left_vertex - self.front_left_vertex

    @property
    def back_leg(self):
        """This method defines a property for the panel's back leg vector as a (3,
        ) array.

        :return: (3,) array of floats
            This is the panel's back leg vector, which is defined from left to
            right. The units are in meters.
        """
        return self.back_right_vertex - self.back_left_vertex

    @property
    def front_right_vortex_vertex(self):
        """This method defines a property for the coordinates of the front-right
        vertex of the ring vortex as a (3,) array.

        :return: (3,) array of floats
            This is the coordinates of the ring vortex's front-right vertex. The
            units are in meters.
        """
        return self.back_right_vertex + 0.75 * self.right_leg

    @property
    def front_left_vortex_vertex(self):
        """This method defines a property for the coordinates of the front-left
        vertex of the ring vortex as a (3,) array.

        :return: (3,) array of floats
            This is the coordinates of the ring vortex's front-left vertex. The units
            are in meters.
        """
        return self.front_left_vertex + 0.25 * self.left_leg

    @property
    def collocation_point(self):
        """This method defines a property for the coordinates of the panel's
        collocation point as a (3,) array.

        :return: (3,) array of floats
            This is the coordinates of the panel's collocation point. The units are
            in meters.
        """
        # Find the location of points three quarters of the way down the left and
        # right legs of the panel.
        right_three_quarter_chord_mark = self.back_right_vertex + 0.25 * self.right_leg
        left_three_quarter_chord_mark = self.front_left_vertex + 0.75 * self.left_leg

        # Find the vector between the points three quarters of the way down the left
        # and right legs of the panel.
        three_quarter_chord_vector = (
            left_three_quarter_chord_mark - right_three_quarter_chord_mark
        )

        # Find the collocation point, which is halfway between the points three
        # quarters of the way down the left and right legs of the panel. Then
        # populate the class attribute.
        return right_three_quarter_chord_mark + 0.5 * three_quarter_chord_vector

    @property
    def area(self):
        """This method defines a property which is an estimate of the panel's area.

        This is only an estimate because the surface defined by four line segments in
        3-space is a hyperboloid, and there doesn't seem to be a closed-form equation
        for the surface area of a hyperboloid between four points. Instead,
        we estimate the area using the cross product of panel's diagonal vectors,
        which should be relatively accurate if the panel can be approximated as a
        planar, convex quadrilateral.

        :return: float
            This is an estimate of the panel's area. The units are square meters.
        """
        return np.linalg.norm(self._cross) / 2

    @property
    def unit_normal(self):
        """This method defines a property for an estimate of the panel's unit
        normal vector as a (3,) array.

        :return: (3,) array of floats
            This is an estimate of the panel's unit normal vector as a (3,) array.
            The sign is determined via the right-hand rule given the orientation of
            panel's leg vectors (front-right to front-left to back-left to
            back-right). The units are in meters.
        """
        return self._cross / np.linalg.norm(self._cross)

    @property
    def unit_spanwise(self):
        """This method defines a property for the panel's unit spanwise vector as a (
        3,) array.

        :return: (3,) array of floats
            This is the panel's unit spanwise vector as a (3,) array. The positive
            direction is defined as left to right, which is opposite the direction of
            the front leg. The units are in meters.
        """
        front_spanwise = -self.front_leg
        back_spanwise = self.back_leg

        spanwise = (front_spanwise + back_spanwise) / 2

        return spanwise / np.linalg.norm(spanwise)

    @property
    def unit_chordwise(self):
        """This method defines a property for the panel's unit chordwise vector as a
        (3,) array.

        :return: (3,) array of floats
            This is the panel's unit chordwise vector as a (3,) array. The positive
            direction is defined as front to back. The units are in meters.
        """
        right_chordwise = -self.right_leg
        left_chordwise = self.left_leg

        chordwise = (right_chordwise + left_chordwise) / 2

        return chordwise / np.linalg.norm(chordwise)

    @property
    def average_span(self):
        """This method defines a property for the average span of the panel.

        :return: float
            This is the average span, which is defined as the average of the front
            and back leg lengths. The units are meters.
        """
        front_leg_length = np.linalg.norm(self.front_leg)
        back_leg_length = np.linalg.norm(self.back_leg)

        return (front_leg_length + back_leg_length) / 2

    @property
    def average_chord(self):
        """This method defines a property for the average chord of the panel.

        :return: float
            This is the average chord, which is defined as the average of the right
            and left leg lengths. The units are meters.
        """
        right_leg_length = np.linalg.norm(self.right_leg)
        left_leg_length = np.linalg.norm(self.left_leg)

        return (right_leg_length + left_leg_length) / 2

    @property
    def _first_diagonal(self):
        """This method defines a property for the panel's first diagonal vector.

        :return: (3,) array of floats
            This is the first diagonal vector, which is defined as the vector from
            the back-left vertex to the front-right vertex. The units are meters.
        """
        return self.front_right_vertex - self.back_left_vertex

    @property
    def _second_diagonal(self):
        """This method defines a property for the panel's second diagonal vector.

        :return: (3,) array of floats
            This is the second diagonal vector, which is defined as the vector from
            the back-right vertex to the front-left vertex. The units are meters.
        """
        return self.front_left_vertex - self.back_right_vertex

    @property
    def _cross(self):
        """This method defines a property for cross product of the panel's first and
        second diagonal vectors.

        :return: (3,) array of floats
            This is the cross product of the panel's first and second diagonal
            vectors. The units are meters.
        """
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
                self.ring_vortex.normalizedInducedVelocity_G()
            )
        if self.horseshoe_vortex is not None:
            normalized_induced_velocity += (
                self.horseshoe_vortex.normalizedInducedVelocity_G()
            )

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
            induced_velocity += self.ring_vortex.inducedVelocity_G(point=point)
        if self.horseshoe_vortex is not None:
            induced_velocity += self.horseshoe_vortex.inducedVelocity_G(point=point)

        return induced_velocity

    def calculate_projected_area(self, n_hat):
        """This method calculates the area of the panel projected on some plane
        defined by its unit normal vector.

        :param n_hat: (3,) array of floats
            This is a (3,) array of the components of the projection plane's unit
            normal vector. The vector must have a magnitude of one. The units are
            meters.
        :return: float
            This is the area of the panel projected onto the plane defined by the
            normal vector. The units are square meters.
        """
        # Find the projections of the first and second diagonal vectors onto the
        # plane's unit normal vector.
        proj_n_hat_first_diag = np.dot(self._first_diagonal, n_hat) * n_hat
        proj_n_hat_second_diag = np.dot(self._second_diagonal, n_hat) * n_hat

        # Find the projection of the first and second diagonal onto the plane.
        proj_plane_first_diag = self._first_diagonal - proj_n_hat_first_diag
        proj_plane_second_diag = self._second_diagonal - proj_n_hat_second_diag

        # The projected area is found by dividing the magnitude of cross product of
        # the diagonal vectors by two. Read the area method for a more detailed
        # explanation.
        proj_cross = np.cross(proj_plane_first_diag, proj_plane_second_diag)
        return np.linalg.norm(proj_cross) / 2

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
