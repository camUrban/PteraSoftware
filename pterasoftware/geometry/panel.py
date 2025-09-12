"""This module contains the Panel class.

This module contains the following classes:
    Panel: This class is used to contain the panels of a wing.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from .. import parameter_validation


# ToDo: Add unit tests for this class.
class Panel:
    """This class is used to contain the panels of a wing.

    This class contains the following public methods:
        rightLeg_G: This method sets a property for the Panel's right leg vector (in
        geometry axes).

        frontLeg_G: This method sets a property for the Panel's front leg vector (in
        geometry axes).

        leftLeg_G: This method sets a property for the Panel's left leg vector (in
        geometry axes).

        backLeg_G: This method sets a property for the Panel's back leg vector (in
        geometry axes).

        Frbvp_G_Cg: This method sets a property for the position of this Panel's
        front-right bound vortex point (in geometry axes, relative to the Cg).

        Flbvp_G_Cg: This method sets a property for the position of this Panel's
        front-right bound vortex point (in geometry axes, relative to the Cg).

        Cpp_G_Cg: This method sets a property for the position of this Panel's
        collocation point (in geometry axes, relative to the Cg).

        area: This method sets a property which is an estimate of the Panel's area.

        unitNormal_G: This method sets a property for an estimate of the Panel's unit
        normal vector (in geometry axes).

        unitSpanwise_G: This method sets a property for the Panel's unit spanwise
        vector (in geometry axes).

        unitChordwise_G: This method sets a property for the Panel's unit chordwise
        vector (in geometry axes).

        average_span: This method sets a property for the average span of the Panel.

        average_chord: This method sets a property for the average chord of the Panel.

        normalizedInducedVelocity_G: This method calculates the velocity (in geometry
        axes) induced at a point (in geometry axes, relative to the Cg) by this
        Panel's vortices, assuming a unit vortex strength.

        inducedVelocity_G: This method calculates the velocity (in geometry axes)
        induced at a point (in geometry axes, relative to the Cg) by this Panel's
        vortices with their given vortex strengths.

        calculate_projected_area: This method calculates the area of the Panel
        projected on a plane defined by a given normal vector (in geometry axes).

        update_coefficients: This method updates the Panel's force coefficients.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        Frpp_G_Cg,
        Flpp_G_Cg,
        Blpp_G_Cg,
        Brpp_G_Cg,
        is_leading_edge,
        is_trailing_edge,
    ):
        """This is the initialization method.

        :param Frpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's front right vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param Flpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's front left vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param Blpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's back left vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param Brpp_G_Cg: array-like of 3 numbers

            Position [x, y, z] of the Panel's back right vertex (in geometry axes,
            relative to the CG). Can be a list, tuple, or numpy array of numbers (int
            or float). Values are converted to floats internally. Front, back, left,
            and right are defined with respect to this Panel's position on its wing
            section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.

        :param is_leading_edge: boolLike

            This is True if the Panel is a leading edge Panel on a Wing, and False
            otherwise. It can be a boolean or a NumPy boolean and will be converted
            internally to a boolean.

        :param is_trailing_edge: boolLike

            This is True if the Panel is a trailing edge Panel on a Wing, and False
            otherwise. It can be a boolean or a NumPy boolean and will be converted
            internally to a boolean.
        """

        # Validate and initialize the attributes.
        self.Frpp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Frpp_G_Cg, "Frpp_G_Cg"
        )
        self.Flpp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Flpp_G_Cg, "Flpp_G_Cg"
        )
        self.Blpp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Blpp_G_Cg, "Blpp_G_Cg"
        )
        self.Brpp_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            Brpp_G_Cg, "Brpp_G_Cg"
        )
        self.is_leading_edge = parameter_validation.boolLike_return_bool(
            is_leading_edge, "is_leading_edge"
        )
        self.is_trailing_edge = parameter_validation.boolLike_return_bool(
            is_trailing_edge, "is_trailing_edge"
        )

        # Initialize variables to hold attributes that describe the Panel's position
        # in its Wing's Panel matrix. They will be populated by the meshing function.
        self.is_right_edge = None
        self.is_left_edge = None
        self.local_chordwise_position = None
        self.local_spanwise_position = None

        # Initialize variables to hold the Panel's RingVortex and HorseshoeVortex.
        # These will be populated by the solver.
        self.ring_vortex = None
        self.horseshoe_vortex = None

        # ToDo: Update and standardize these definitions for force coefficients.
        # Initialize variables to hold attributes of the Panel that will be defined
        # after the solver finds a solution.
        self.forces_G = None
        self.moments_G_Cg = None
        self.forces_W = None
        self.moments_W_Cg = None
        self.induced_drag_coefficient = None
        self.side_force_coefficient = None
        self.lift_coefficient = None

    @property
    def rightLeg_G(self):
        """This method sets a property for the Panel's right leg vector (in geometry
        axes).

        :return: (3,) ndarray of floats
            This is the Panel's right leg vector, which is defined from back to
            front. The units are in meters.
        """
        return self.Frpp_G_Cg - self.Brpp_G_Cg

    @property
    def frontLeg_G(self):
        """This method sets a property for the Panel's front leg vector (in geometry
        axes).

        :return: (3,) ndarray of floats
            This is the Panel's front leg vector, which is defined from right to
            left. The units are in meters.
        """
        return self.Flpp_G_Cg - self.Frpp_G_Cg

    @property
    def leftLeg_G(self):
        """This method sets a property for the Panel's left leg vector (in geometry
        axes).

        :return: (3,) ndarray of floats
            This is the Panel's left leg vector, which is defined from front to
            back. The units are in meters.
        """
        return self.Blpp_G_Cg - self.Flpp_G_Cg

    @property
    def backLeg_G(self):
        """This method sets a property for the Panel's back leg vector (in geometry
        axes).

        :return: (3,) ndarray of floats
            This is the Panel's back leg vector, which is defined from left to
            right. The units are in meters.
        """
        return self.Brpp_G_Cg - self.Blpp_G_Cg

    @property
    def Frbvp_G_Cg(self):
        """This method sets a property for the position of this Panel's front-right
        bound vortex point (in geometry axes, relative to the Cg).

        :return: (3,) ndarray of floats
            This is the position of this Panel's front-right bound vortex point. The
            units are in meters.
        """
        return self.Brpp_G_Cg + 0.75 * self.rightLeg_G

    @property
    def Flbvp_G_Cg(self):
        """This method sets a property for the position of this Panel's front-left
        bound vortex point (in geometry axes, relative to the Cg).

        :return: (3,) ndarray of floats
            This is the position of this Panel's front-right bound vortex point. The
            units are in meters.
        """
        return self.Flpp_G_Cg + 0.25 * self.leftLeg_G

    @property
    def Cpp_G_Cg(self):
        """This method sets a property for the position of this Panel's collocation
        point (in geometry axes, relative to the Cg).

        :return: (3,) ndarray of floats
            This is the position of this Panel's collocation point. The units are in
            meters.
        """
        # Find the positions of points three quarters of the way down the left and
        # right legs of the Panel (in geometry axes, relative to the CG).
        rightThreeQuarterChord_G_Cg = self.Brpp_G_Cg + 0.25 * self.rightLeg_G
        leftThreeQuarterChord_G_Cg = self.Flpp_G_Cg + 0.75 * self.leftLeg_G

        # Find the vector (in geometry axes) between the points three quarters of the
        # way down the left and right legs of the Panel.
        threeQuarterChord_G = leftThreeQuarterChord_G_Cg - rightThreeQuarterChord_G_Cg

        # Find the collocation point (in geometry axes, relative to the CG), which is
        # halfway between the points three quarters of the way down the left and
        # right legs of the Panel. Then populate the class attribute.
        return rightThreeQuarterChord_G_Cg + 0.5 * threeQuarterChord_G

    @property
    def area(self):
        """This method sets a property which is an estimate of the Panel's area.

        This is only an estimate because the surface defined by four line segments in
        3-space is a hyperboloid, and there doesn't seem to be a closed-form equation
        for the surface area of a hyperboloid between four points. Instead,
        we estimate the area using the cross product of Panel's diagonal vectors,
        which should be relatively accurate if the Panel can be approximated as a
        planar, convex quadrilateral.

        :return: float
            This is an estimate of the Panel's area. The units are square meters.
        """
        return np.linalg.norm(self._cross_G) / 2

    @property
    def unitNormal_G(self):
        """This method sets a property for an estimate of the Panel's unit normal
        vector (in geometry axes).

        :return: (3,) ndarray of floats
            This is an estimate of the Panel's unit normal vector as a (3,) ndarray
            of floats. The sign is determined via the right-hand rule given the
            orientation of Panel's leg vectors (front-right to front-left to
            back-left to back-right). The units are in meters.
        """
        return self._cross_G / np.linalg.norm(self._cross_G)

    @property
    def unitSpanwise_G(self):
        """This method sets a property for the Panel's unit spanwise vector (in
        geometry axes).

        :return: (3,) ndarray of floats
            This is the Panel's unit spanwise vector as a (3,) ndarray of floats. The
            positive direction is defined as left to right, which is opposite the
            direction of the front leg. The units are in meters.
        """
        frontSpanwise_G = -self.frontLeg_G
        backSpanwise_G = self.backLeg_G

        spanwise_G = (frontSpanwise_G + backSpanwise_G) / 2

        return spanwise_G / np.linalg.norm(spanwise_G)

    @property
    def unitChordwise_G(self):
        """This method sets a property for the Panel's unit chordwise vector (in
        geometry axes).

        :return: (3,) ndarray of floats
            This is the Panel's unit chordwise vector as a (3,) ndarray of floats.
            The positive direction is defined as front to back. The units are in meters.
        """
        rightChordwise_G = -self.rightLeg_G
        leftChordwise_G = self.leftLeg_G

        chordwise_G = (rightChordwise_G + leftChordwise_G) / 2

        return chordwise_G / np.linalg.norm(chordwise_G)

    @property
    def average_span(self):
        """This method sets a property for the average span of the Panel.

        :return: float
            This is the average span, which is defined as the average of the front
            and back leg lengths. The units are meters.
        """
        front_leg_length = np.linalg.norm(self.frontLeg_G)
        back_leg_length = np.linalg.norm(self.backLeg_G)

        return (front_leg_length + back_leg_length) / 2

    @property
    def average_chord(self):
        """This method sets a property for the average chord of the Panel.

        :return: float
            This is the average chord, which is defined as the average of the right
            and left leg lengths. The units are meters.
        """
        right_leg_length = np.linalg.norm(self.rightLeg_G)
        left_leg_length = np.linalg.norm(self.leftLeg_G)

        return (right_leg_length + left_leg_length) / 2

    @property
    def _firstDiagonal_G(self):
        """This method sets a property for the Panel's first diagonal vector (in
        geometry axes).

        :return: (3,) ndarray of floats
            This is the first diagonal vector, which is defined as the vector from
            the back-left panel point to the front-right panel point. The units are
            meters.
        """
        return self.Frpp_G_Cg - self.Blpp_G_Cg

    @property
    def _secondDiagonal_G(self):
        """This method sets a property for the Panel's second diagonal vector (in
        geometry axes).

        :return: (3,) ndarray of floats
            This is the second diagonal vector, which is defined as the vector from
            the back-right panel point to the front-left panel point. The units are
            meters.
        """
        return self.Flpp_G_Cg - self.Brpp_G_Cg

    @property
    def _cross_G(self):
        """This method sets a property for cross product (in geometry axes) of the
        Panel's first and second diagonal vectors.

        :return: (3,) ndarray of floats
            This is the cross product of the Panel's first and second diagonal
            vectors. The units are meters.
        """
        return np.cross(self._firstDiagonal_G, self._secondDiagonal_G)

    def normalizedInducedVelocity_G(self, point_G_Cg):
        """This method calculates the velocity (in geometry axes) induced at a point
        (in geometry axes, relative to the Cg) by this Panel's vortices, assuming a
        unit vortex strength.

        This method does not include the effect any wake vortices.

        :param point_G_Cg: array-like of 3 numbers
            Position [x, y, z] of the point (in geometry axes, relative to the CG) to
            find the induced velocity at. Can be a list, tuple, or numpy array of
            numbers (int or float). Values are converted to floats internally. The
            units are in meters.
        :return: (3,) ndarray of floats
            This is the normalized induced velocity (in geometry axes). The units are
            in meters per second.
        """
        point_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            point_G_Cg, "point_G_Cg"
        )

        normalizedInducedVelocity_G = np.zeros(3)

        if self.ring_vortex is not None:
            normalizedInducedVelocity_G += self.ring_vortex.normalizedInducedVelocity_G(
                point_G_Cg
            )
        if self.horseshoe_vortex is not None:
            normalizedInducedVelocity_G += (
                self.horseshoe_vortex.normalizedInducedVelocity_G(point_G_Cg)
            )

        return normalizedInducedVelocity_G

    def inducedVelocity_G(self, point_G_Cg):
        """This method calculates the velocity (in geometry axes) induced at a point
        (in geometry axes, relative to the Cg) by this Panel's vortices with their
        given vortex strengths.

        This method does not include the effect of any wake vortices.

        :param point_G_Cg: array-like of 3 numbers
            Position [x, y, z] of the point (in geometry axes, relative to the CG) to
            find the induced velocity at. Can be a list, tuple, or numpy array of
            numbers (int or float). Values are converted to floats internally. The
            units are in meters.
        :return: (3,) ndarray of floats
            This is the induced velocity (in geometry axes). The units are in meters
            per second.
        """
        point_G_Cg = parameter_validation.threeD_number_vectorLike_return_float(
            point_G_Cg, "point_G_Cg"
        )

        inducedVelocity_G = np.zeros(3)

        if self.ring_vortex is not None:
            inducedVelocity_G += self.ring_vortex.inducedVelocity_G(point_G_Cg)
        if self.horseshoe_vortex is not None:
            inducedVelocity_G += self.horseshoe_vortex.inducedVelocity_G(point_G_Cg)

        return inducedVelocity_G

    def calculate_projected_area(self, normal_G):
        """This method calculates the area of the Panel projected on a plane
        defined by a given normal vector (in geometry axes).

        :param normal_G: array-like of 3 numbers
            Normal vector (in geometry axes) defining the plane that will be used to
            calculate the projected area. Can be a list, tuple, or numpy array of
            numbers (int or float). Values are converted to floats internally. If not
            already a unit vector, it will also be converted to one. The units are in
            meters.
        :return: float
            This is the area of the Panel projected onto the plane defined by the given
            normal vector. The units are square meters.
        """
        # Validate and normalize the normal vector.
        unitNormal_G = (
            parameter_validation.threeD_number_vectorLike_return_float_unit_vector(
                normal_G, "normal_G"
            )
        )

        # Find the projections of the first and second diagonal vectors (in geometry
        # axes) onto the plane's unit normal vector.
        projFirstDiagonalOnNormal_G = (
            np.dot(self._firstDiagonal_G, unitNormal_G) * unitNormal_G
        )
        projSecondDiagonalOnNormal_G = (
            np.dot(self._secondDiagonal_G, unitNormal_G) * unitNormal_G
        )

        # Find the projection (in geometry axes) of the first and second diagonal
        # vectors onto the plane.
        projFirstDiagonalOnPlane_G = self._firstDiagonal_G - projFirstDiagonalOnNormal_G
        projSecondDiagonalOnPlane_G = (
            self._secondDiagonal_G - projSecondDiagonalOnNormal_G
        )

        # The projected area is found by dividing the magnitude of cross product of
        # the diagonal vectors (in geometry axes) by two. Read the area method for a
        # more detailed explanation.
        projDiagonalsOnPlaneCross_G = np.cross(
            projFirstDiagonalOnPlane_G, projSecondDiagonalOnPlane_G
        )
        return np.linalg.norm(projDiagonalsOnPlaneCross_G) / 2

    def update_coefficients(self, freestream_q):
        """This method updates the Panel's force coefficients.

        :param freestream_q: number
            This is the freestream dynamic pressure, as calculated for the current
            operating point. It must be a positive number, and will be converted to a
            float internally. The units are in Pascals.
        :return: None
        """
        freestream_q = parameter_validation.positive_number_return_float(
            freestream_q, "freestream_q"
        )

        # ToDo: Update and standardize these definitions for force coefficients.
        induced_drag = -self.forces_W[0]
        side_force = self.forces_W[1]
        lift = -self.forces_W[2]

        self.induced_drag_coefficient = induced_drag / self.area / freestream_q
        self.side_force_coefficient = side_force / self.area / freestream_q
        self.lift_coefficient = lift / self.area / freestream_q
