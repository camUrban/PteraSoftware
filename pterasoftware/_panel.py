"""Contains the Panel class."""

from __future__ import annotations

from typing import cast

import numpy as np

from . import _aerodynamics


class Panel:
    """A class used to contain the panels of a Wing.

    **Contains the following methods:**

    rightLeg_G: This Panel's right leg vector (in geometry axes).

    frontLeg_G: This Panel's front leg vector (in geometry axes).

    leftLeg_G: This Panel's left leg vector (in geometry axes).

    backLeg_G: This Panel's back leg vector (in geometry axes).

    Frbvp_G_Cg: The position of this Panel's front right bound vortex point (in geometry
    axes, relative to the CG).

    Flbvp_G_Cg: The position of this Panel's front left bound vortex point (in geometry
    axes, relative to the CG).

    Cpp_G_Cg: The position of this Panel's collocation point (in geometry axes, relative
    to the CG).

    rightLeg_GP1: This Panel's right leg vector ( in the first Airplane's geometry
    axes).

    frontLeg_GP1: This Panel's front leg vector ( in the first Airplane's geometry
    axes).

    leftLeg_GP1: This Panel's left leg vector (in the first Airplane's geometry axes).

    backLeg_GP1: This Panel's back leg vector (in the first Airplane's geometry axes).

    Frbvp_GP1_CgP1: The position of this Panel's front right bound vortex point (in the
    first Airplane's geometry axes, relative to the first Airplane's CG).

    Flbvp_GP1_CgP1: The position of this Panel's front left bound vortex point (in the
    first Airplane's geometry axes, relative to the first Airplane's CG).

    Cpp_GP1_CgP1: The position of this Panel's collocation point (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    area: An estimate of this Panel's area.

    unitNormal_G: An estimate of this Panel's unit normal vector (in geometry axes).

    unitNormal_GP1: An estimate of this Panel's unit normal vector (in the first
    Airplane's geometry axes).

    aspect_ratio: The aspect ratio of this Panel.

    calculate_projected_area: The area of this Panel projected on a plane defined by a
    given normal vector (in geometry axes).
    """

    def __init__(
        self,
        Frpp_G_Cg: np.ndarray,
        Flpp_G_Cg: np.ndarray,
        Blpp_G_Cg: np.ndarray,
        Brpp_G_Cg: np.ndarray,
        is_leading_edge: bool,
        is_trailing_edge: bool,
    ) -> None:
        """The initialization method.

        :param Frpp_G_Cg: A (3,) ndarray of floats representing the position of the
            Panel's front right vertex (in geometry axes, relative to the CG). Front,
            back, left, and right are defined with respect to this Panel's position on
            its wing section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.
        :param Flpp_G_Cg: A (3,) ndarray of floats representing the position of the
            Panel's front left vertex (in geometry axes, relative to the CG). Front,
            back, left, and right are defined with respect to this Panel's position on
            its wing section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.
        :param Blpp_G_Cg: A (3,) ndarray of floats representing the position of the
            Panel's back left vertex (in geometry axes, relative to the CG). Front,
            back, left, and right are defined with respect to this Panel's position on
            its wing section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.
        :param Brpp_G_Cg: A (3,) ndarray of floats representing the position of the
            Panel's back right vertex (in geometry axes, relative to the CG). Front,
            back, left, and right are defined with respect to this Panel's position on
            its wing section, with front meaning towards the leading edge, right meaning
            towards the next wing section or the Wing's tip, etc. The units are in
            meters.
        :param is_leading_edge: Flags if this Panel is at its parent Wing's leading
            edge.
        :param is_trailing_edge: Flags if this Panel is at its parent Wing's trailing
            edge.
        :return: None
        """
        # Initialize the attributes.
        self.Frpp_G_Cg = Frpp_G_Cg
        self.Flpp_G_Cg = Flpp_G_Cg
        self.Blpp_G_Cg = Blpp_G_Cg
        self.Brpp_G_Cg = Brpp_G_Cg
        self.is_leading_edge = is_leading_edge
        self.is_trailing_edge = is_trailing_edge

        # Initialize variables to hold attributes that describe the Panel's position
        # in its Wing's Panel matrix. They will be populated by the meshing function.
        self.is_right_edge: bool | None = None
        self.is_left_edge: bool | None = None
        self.local_chordwise_position: int | None = None
        self.local_spanwise_position: int | None = None

        # Initialize variables to hold the Panel's corner positions in formation
        # flight coordinates (in the first Airplane's geometry axes, relative to the
        # first Airplane's CG). These will be populated by the SteadyProblem or
        # UnsteadyProblem during initialization.
        self.Frpp_GP1_CgP1: np.ndarray | None = None
        self.Flpp_GP1_CgP1: np.ndarray | None = None
        self.Blpp_GP1_CgP1: np.ndarray | None = None
        self.Brpp_GP1_CgP1: np.ndarray | None = None

        # Initialize variables to hold the Panel's RingVortex and HorseshoeVortex.
        # These will be populated by the solver.
        self.ring_vortex: _aerodynamics.RingVortex | None = None
        self.horseshoe_vortex: _aerodynamics.HorseshoeVortex | None = None

        # Initialize variables to hold attributes of the Panel that will be defined
        # after the solver finds a solution.
        self.forces_GP1: np.ndarray | None = None
        self.moments_GP1_CgP1: np.ndarray | None = None
        self.forces_W: np.ndarray | None = None
        self.moments_W_CgP1: np.ndarray | None = None

    @property
    def rightLeg_G(self) -> np.ndarray:
        """This Panel's right leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's right leg vector,
            which is defined from back to front. The units are in meters.
        """
        return cast(np.ndarray, self.Frpp_G_Cg - self.Brpp_G_Cg)

    @property
    def frontLeg_G(self) -> np.ndarray:
        """This Panel's front leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's front leg vector,
            which is defined from right to left. The units are in meters.
        """
        return cast(np.ndarray, self.Flpp_G_Cg - self.Frpp_G_Cg)

    @property
    def leftLeg_G(self) -> np.ndarray:
        """This Panel's left leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's left leg vector,
            which is defined from front to back. The units are in meters.
        """
        return cast(np.ndarray, self.Blpp_G_Cg - self.Flpp_G_Cg)

    @property
    def backLeg_G(self) -> np.ndarray:
        """This Panel's back leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's back leg vector,
            which is defined from left to right. The units are in meters.
        """
        return cast(np.ndarray, self.Brpp_G_Cg - self.Blpp_G_Cg)

    @property
    def Frbvp_G_Cg(self) -> np.ndarray:
        """The position of this Panel's front right bound vortex point (in geometry
        axes, relative to the CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front right bound vortex point. The units are in meters.
        """
        return cast(np.ndarray, self.Brpp_G_Cg + 0.75 * self.rightLeg_G)

    @property
    def Flbvp_G_Cg(self) -> np.ndarray:
        """The position of this Panel's front left bound vortex point (in geometry axes,
        relative to the CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front left bound vortex point. The units are in meters.
        """
        return cast(np.ndarray, self.Flpp_G_Cg + 0.25 * self.leftLeg_G)

    @property
    def Cpp_G_Cg(self) -> np.ndarray:
        """The position of this Panel's collocation point (in geometry axes, relative to
        the CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            collocation point. The units are in meters.
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

    # TEST: Consider adding unit tests for this method.
    @property
    def rightLeg_GP1(self) -> np.ndarray | None:
        """This Panel's right leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's right leg vector,
            which is defined from back to front. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self.Frpp_GP1_CgP1 is None or self.Brpp_GP1_CgP1 is None:
            return None
        return cast(np.ndarray, self.Frpp_GP1_CgP1 - self.Brpp_GP1_CgP1)

    # TEST: Consider adding unit tests for this method.
    @property
    def frontLeg_GP1(self) -> np.ndarray | None:
        """This Panel's front leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's left leg vector,
            which is defined from right to left. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self.Flpp_GP1_CgP1 is None or self.Frpp_GP1_CgP1 is None:
            return None
        return cast(np.ndarray, self.Flpp_GP1_CgP1 - self.Frpp_GP1_CgP1)

    # TEST: Consider adding unit tests for this method.
    @property
    def leftLeg_GP1(self) -> np.ndarray | None:
        """This Panel's left leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's left leg vector,
            which is defined from front to back. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self.Blpp_GP1_CgP1 is None or self.Flpp_GP1_CgP1 is None:
            return None
        return cast(np.ndarray, self.Blpp_GP1_CgP1 - self.Flpp_GP1_CgP1)

    # TEST: Consider adding unit tests for this method.
    @property
    def backLeg_GP1(self) -> np.ndarray | None:
        """This Panel's back leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's back leg vector,
            which is defined from left to right. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self.Brpp_GP1_CgP1 is None or self.Blpp_GP1_CgP1 is None:
            return None
        return cast(np.ndarray, self.Brpp_GP1_CgP1 - self.Blpp_GP1_CgP1)

    # TEST: Consider adding unit tests for this method.
    @property
    def Frbvp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of this Panel's front right bound vortex point (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front right bound vortex point. The units are in meters. Returns None if
            this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self.Brpp_GP1_CgP1 is None or self.rightLeg_GP1 is None:
            return None
        return self.Brpp_GP1_CgP1 + 0.75 * self.rightLeg_GP1

    # TEST: Consider adding unit tests for this method.
    @property
    def Flbvp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of this Panel's front left bound vortex point (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front left bound vortex point. The units are in meters. Returns None if this
            Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self.Flpp_GP1_CgP1 is None or self.leftLeg_GP1 is None:
            return None
        return self.Flpp_GP1_CgP1 + 0.25 * self.leftLeg_GP1

    # TEST: Consider adding unit tests for this method.
    @property
    def Cpp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of this Panel's collocation point (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            collocation point. The units are in meters. Returns None if this Panel is
            not part of a SteadyProblem or UnsteadyProblem.
        """
        if (
            self.Brpp_GP1_CgP1 is None
            or self.rightLeg_GP1 is None
            or self.Flpp_GP1_CgP1 is None
            or self.leftLeg_GP1 is None
        ):
            return None

        # Find the positions of points three quarters of the way down the left and
        # right legs of the Panel (in the first Airplane's geometry axes, relative to
        # the first Airplane's CG).
        rightThreeQuarterChord_GP1_CgP1 = self.Brpp_GP1_CgP1 + 0.25 * self.rightLeg_GP1
        leftThreeQuarterChord_GP1_CgP1 = self.Flpp_GP1_CgP1 + 0.75 * self.leftLeg_GP1

        # Find the vector (in the first Airplane's geometry axes) between the points
        # three quarters of the way down the left and right legs of the Panel.
        threeQuarterChord_GP1 = (
            leftThreeQuarterChord_GP1_CgP1 - rightThreeQuarterChord_GP1_CgP1
        )

        # Find the collocation point (in the first Airplane's geometry axes, relative
        # to the first Airplane's CG), which is halfway between the points three
        # quarters of the way down the left and right legs of the Panel.
        return rightThreeQuarterChord_GP1_CgP1 + 0.5 * threeQuarterChord_GP1

    @property
    def area(self) -> float:
        """An estimate of this Panel's area.

        This is only an estimate because the surface defined by four line segments in
        3-space is a hyperboloid, and there doesn't seem to be a closed-form equation
        for the surface area of a hyperboloid between four points. Instead, we estimate
        the area using the cross product of Panel's diagonal vectors, which should be
        relatively accurate if the Panel can be approximated as a planar, convex
        quadrilateral.

        :return: An estimate of the Panel's area. The units are square meters.
        """
        return float(np.linalg.norm(self._cross_G) / 2)

    @property
    def unitNormal_G(self) -> np.ndarray:
        """An estimate of this Panel's unit normal vector (in geometry axes).

        :return: A (3,) ndarray of floats representing an estimate of this Panel's unit
            normal vector. The sign is determined via the right-hand rule given the
            orientation of Panel's leg vectors (front right to front left to back left
            to back right).
        """
        return cast(np.ndarray, self._cross_G / np.linalg.norm(self._cross_G))

    # TEST: Consider adding unit tests for this method.
    @property
    def unitNormal_GP1(self) -> np.ndarray | None:
        """An estimate of this Panel's unit normal vector (in the first Airplane's
        geometry axes).

        :return: A (3,) ndarray of floats representing an estimate of this Panel's unit
            normal vector. The sign is determined via the right-hand rule given the
            orientation of Panel's leg vectors (front right to front left to back left
            to back right). Returns None if this Panel is not part of a SteadyProblem or
            UnsteadyProblem.
        """
        if (
            self.Frpp_GP1_CgP1 is None
            or self.Flpp_GP1_CgP1 is None
            or self.Blpp_GP1_CgP1 is None
            or self.Brpp_GP1_CgP1 is None
        ):
            return None

        # Compute diagonal vectors (in the first Airplane's geometry axes).
        firstDiagonal_GP1 = self.Frpp_GP1_CgP1 - self.Blpp_GP1_CgP1
        secondDiagonal_GP1 = self.Flpp_GP1_CgP1 - self.Brpp_GP1_CgP1

        # Compute the cross product and normalize.
        cross_GP1 = np.cross(firstDiagonal_GP1, secondDiagonal_GP1)
        return cast(np.ndarray, cross_GP1 / np.linalg.norm(cross_GP1))

    # TEST: Consider adding unit tests for this method.
    @property
    def aspect_ratio(self) -> float:
        """The aspect ratio of this Panel.

        :return: The Panel's aspect ratio, which is defined as the distance between the
            right and left legs' center points divided by the distance between the front
            and back legs' center points.
        """
        frontCenterPoint_G_Cg = self.Frpp_G_Cg + self.frontLeg_G / 2
        leftCenterPoint_G_Cg = self.Flpp_G_Cg + self.leftLeg_G / 2
        backCenterPoint_G_Cg = self.Blpp_G_Cg + self.backLeg_G / 2
        rightCenterPoint_G_Cg = self.Brpp_G_Cg + self.rightLeg_G / 2

        right_left_distance = float(
            np.linalg.norm(rightCenterPoint_G_Cg - leftCenterPoint_G_Cg)
        )
        front_back_distance = float(
            np.linalg.norm(frontCenterPoint_G_Cg - backCenterPoint_G_Cg)
        )

        return right_left_distance / front_back_distance

    @property
    def _firstDiagonal_G(self) -> np.ndarray:
        """This Panel's first diagonal vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the Panel's first diagonal
            vector, which is defined as the vector from the back left panel point to the
            front right panel point. The units are in meters.
        """
        return cast(np.ndarray, self.Frpp_G_Cg - self.Blpp_G_Cg)

    @property
    def _secondDiagonal_G(self) -> np.ndarray:
        """This Panel's second diagonal vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the Panel's second diagonal
            vector, which is defined as the vector from the back right panel point to
            the front left panel point. The units are in meters.
        """
        return cast(np.ndarray, self.Flpp_G_Cg - self.Brpp_G_Cg)

    @property
    def _cross_G(self) -> np.ndarray:
        """The cross product (in geometry axes) of this Panel's first and second
        diagonal vectors.

        :return: A (3,) ndarray of floats representing the cross product of this Panel's
            first and second diagonal vectors.
        """
        return cast(np.ndarray, np.cross(self._firstDiagonal_G, self._secondDiagonal_G))

    def calculate_projected_area(self, normal_G: np.ndarray) -> float:
        """Calculates the area of this Panel projected on a plane defined by a given
        normal vector (in geometry axes).

        :param normal_G: A (3,) ndarray of floats representing the normal vector
            defining the plane that will be used to calculate the projected area.
        :return: The area of the Panel projected onto the plane defined by the given
            normal vector. The units are in square meters.
        """
        # Normalize the normal vector.
        unitNormal_G = normal_G / np.linalg.norm(normal_G)

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
        return float(np.linalg.norm(projDiagonalsOnPlaneCross_G) / 2)
