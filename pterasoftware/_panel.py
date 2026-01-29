"""Contains the Panel class."""

from __future__ import annotations

from typing import cast

import numpy as np

from . import _vortices


class Panel:
    """A class used to contain the panels of a Wing.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this Panel, preserving mesh geometry but
    resetting solver state.

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

    unitNormal_G: An estimate of this Panel's unit normal vector (in geometry axes).

    area: An estimate of this Panel's area.

    aspect_ratio: The aspect ratio of this Panel.

    Frpp_GP1_CgP1: The position of the Panel's front right vertex (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    Flpp_GP1_CgP1: The position of the Panel's front left vertex (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    Blpp_GP1_CgP1: The position of the Panel's back left vertex (in the first Airplane's
    geometry axes, relative to the first Airplane's CG).

    Brpp_GP1_CgP1: The position of the Panel's back right vertex (in the first
    Airplane's geometry axes, relative to the first Airplane's CG).

    right_edge: Flags if this Panel is at its parent Wing's right edge.

    left_edge: Flags if this Panel is at its parent Wing's left edge.

    local_chordwise_position: This Panel's local chordwise position in its parent Wing's
    mesh.

    local_spanwise_position: This Panel's local spanwise position in its parent Wing's
    mesh.

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

    unitNormal_GP1: An estimate of this Panel's unit normal vector (in the first
    Airplane's geometry axes).

    calculate_projected_area: The area of this Panel projected on a plane defined by a
    given normal vector (in geometry axes).

    **Notes:**

    Computed geometric properties (leg vectors, bound vortex points, collocation points,
    unit normals, area, and aspect ratio) are lazily evaluated and cached.
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
        # Initialize the immutable attributes. Set those that are numpy arrays to be
        # read only.
        self._Frpp_G_Cg = Frpp_G_Cg
        self._Frpp_G_Cg.flags.writeable = False
        self._Flpp_G_Cg = Flpp_G_Cg
        self._Flpp_G_Cg.flags.writeable = False
        self._Blpp_G_Cg = Blpp_G_Cg
        self._Blpp_G_Cg.flags.writeable = False
        self._Brpp_G_Cg = Brpp_G_Cg
        self._Brpp_G_Cg.flags.writeable = False
        self._is_leading_edge = is_leading_edge
        self._is_trailing_edge = is_trailing_edge

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._rightLeg_G: np.ndarray | None = None
        self._frontLeg_G: np.ndarray | None = None
        self._leftLeg_G: np.ndarray | None = None
        self._backLeg_G: np.ndarray | None = None
        self._Frbvp_G_Cg: np.ndarray | None = None
        self._Flbvp_G_Cg: np.ndarray | None = None
        self._Cpp_G_Cg: np.ndarray | None = None
        self._unitNormal_G: np.ndarray | None = None
        self._area: float | None = None
        self._aspect_ratio: float | None = None

        # Initialize the set once attributes. The attributes for the Panel's corner
        # positions in formation flight coordinates (in the first Airplane's geometry
        # axes, relative to the first Airplane's CG) will be populated by the
        # SteadyProblem or UnsteadyProblem during initialization. The attributes that
        # describe the Panel's position in its Wing's Panel matrix will be populated by
        # the meshing function.
        self._Frpp_GP1_CgP1: np.ndarray | None = None
        self._Flpp_GP1_CgP1: np.ndarray | None = None
        self._Blpp_GP1_CgP1: np.ndarray | None = None
        self._Brpp_GP1_CgP1: np.ndarray | None = None
        self._is_right_edge: bool | None = None
        self._is_left_edge: bool | None = None
        self._local_chordwise_position: int | None = None
        self._local_spanwise_position: int | None = None

        # Initialize the caches for the properties derived from the set once properties.
        self._rightLeg_GP1: np.ndarray | None = None
        self._frontLeg_GP1: np.ndarray | None = None
        self._leftLeg_GP1: np.ndarray | None = None
        self._backLeg_GP1: np.ndarray | None = None
        self._Frbvp_GP1_CgP1: np.ndarray | None = None
        self._Flbvp_GP1_CgP1: np.ndarray | None = None
        self._Cpp_GP1_CgP1: np.ndarray | None = None
        self._unitNormal_GP1: np.ndarray | None = None

        # Initialize mutable attributes to hold the Panel's RingVortex and
        # HorseshoeVortex and the loads on this Panel.
        self.ring_vortex: _vortices.ring_vortex.RingVortex | None = None
        self.horseshoe_vortex: _vortices.horseshoe_vortex.HorseshoeVortex | None = None
        self.forces_GP1: np.ndarray | None = None
        self.moments_GP1_CgP1: np.ndarray | None = None
        self.forces_W: np.ndarray | None = None
        self.moments_W_CgP1: np.ndarray | None = None

    def __deepcopy__(self, memo: dict) -> Panel:
        """Creates a deep copy of this Panel, preserving mesh geometry but resetting
        solver state.

        The copy preserves: (1) local corner positions (in geometry axes, relative to
        the CG), (2) mesh metadata (edge flags, grid positions), and (3) cached local
        geometric properties (leg vectors, bound vortex points, the collocation point,
        the unit normal vector, the area, and the aspect ratio).

        The copy resets to None: (1) global positions (in the first Airplane's geometry
        axes), (2) cached global geometric properties, (3) vortex objects (ring_vortex,
        horseshoe_vortex), and (4) loads.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new Panel with preserved mesh geometry and reset solver state.
        """
        # Create a new Panel instance without calling __init__ to avoid redundant
        # cache invalidation.
        new_panel = object.__new__(Panel)

        # Store this Panel in memo to handle potential circular references.
        memo[id(self)] = new_panel

        # Copy immutable attributes. For those that are numpy arrays, make the copies
        # read only.
        new_panel._Frpp_G_Cg = self._Frpp_G_Cg.copy()
        new_panel._Frpp_G_Cg.flags.writeable = False
        new_panel._Flpp_G_Cg = self._Flpp_G_Cg.copy()
        new_panel._Flpp_G_Cg.flags.writeable = False
        new_panel._Blpp_G_Cg = self._Blpp_G_Cg.copy()
        new_panel._Blpp_G_Cg.flags.writeable = False
        new_panel._Brpp_G_Cg = self._Brpp_G_Cg.copy()
        new_panel._Brpp_G_Cg.flags.writeable = False
        new_panel._is_leading_edge = self._is_leading_edge
        new_panel._is_trailing_edge = self._is_trailing_edge

        # Copy local set once attributes.
        new_panel._is_right_edge = self._is_right_edge
        new_panel._is_left_edge = self._is_left_edge
        new_panel._local_chordwise_position = self._local_chordwise_position
        new_panel._local_spanwise_position = self._local_spanwise_position

        # Copy cached local geometric properties. This preserves computation from
        # meshing. For those that are numpy arrays, make the copies read only.
        new_panel._rightLeg_G = (
            self._rightLeg_G.copy() if self._rightLeg_G is not None else None
        )
        if new_panel._rightLeg_G is not None:
            new_panel._rightLeg_G.flags.writeable = False

        new_panel._frontLeg_G = (
            self._frontLeg_G.copy() if self._frontLeg_G is not None else None
        )
        if new_panel._frontLeg_G is not None:
            new_panel._frontLeg_G.flags.writeable = False

        new_panel._leftLeg_G = (
            self._leftLeg_G.copy() if self._leftLeg_G is not None else None
        )
        if new_panel._leftLeg_G is not None:
            new_panel._leftLeg_G.flags.writeable = False

        new_panel._backLeg_G = (
            self._backLeg_G.copy() if self._backLeg_G is not None else None
        )
        if new_panel._backLeg_G is not None:
            new_panel._backLeg_G.flags.writeable = False

        new_panel._Frbvp_G_Cg = (
            self._Frbvp_G_Cg.copy() if self._Frbvp_G_Cg is not None else None
        )
        if new_panel._Frbvp_G_Cg is not None:
            new_panel._Frbvp_G_Cg.flags.writeable = False

        new_panel._Flbvp_G_Cg = (
            self._Flbvp_G_Cg.copy() if self._Flbvp_G_Cg is not None else None
        )
        if new_panel._Flbvp_G_Cg is not None:
            new_panel._Flbvp_G_Cg.flags.writeable = False

        new_panel._Cpp_G_Cg = (
            self._Cpp_G_Cg.copy() if self._Cpp_G_Cg is not None else None
        )
        if new_panel._Cpp_G_Cg is not None:
            new_panel._Cpp_G_Cg.flags.writeable = False

        new_panel._unitNormal_G = (
            self._unitNormal_G.copy() if self._unitNormal_G is not None else None
        )
        if new_panel._unitNormal_G is not None:
            new_panel._unitNormal_G.flags.writeable = False

        new_panel._area = self._area
        new_panel._aspect_ratio = self._aspect_ratio

        # Set global positions, cached global geometric attributes, vortex objects, and
        # loads to None (the solver will set/create/compute these).
        new_panel._Frpp_GP1_CgP1 = None
        new_panel._Flpp_GP1_CgP1 = None
        new_panel._Blpp_GP1_CgP1 = None
        new_panel._Brpp_GP1_CgP1 = None
        new_panel._rightLeg_GP1 = None
        new_panel._frontLeg_GP1 = None
        new_panel._leftLeg_GP1 = None
        new_panel._backLeg_GP1 = None
        new_panel._Frbvp_GP1_CgP1 = None
        new_panel._Flbvp_GP1_CgP1 = None
        new_panel._Cpp_GP1_CgP1 = None
        new_panel._unitNormal_GP1 = None
        new_panel.ring_vortex = None
        new_panel.horseshoe_vortex = None
        new_panel.forces_GP1 = None
        new_panel.moments_GP1_CgP1 = None
        new_panel.forces_W = None
        new_panel.moments_W_CgP1 = None

        return new_panel

    # --- Immutable: read only properties ---
    @property
    def Frpp_G_Cg(self) -> np.ndarray:
        return self._Frpp_G_Cg

    @property
    def Flpp_G_Cg(self) -> np.ndarray:
        return self._Flpp_G_Cg

    @property
    def Blpp_G_Cg(self) -> np.ndarray:
        return self._Blpp_G_Cg

    @property
    def Brpp_G_Cg(self) -> np.ndarray:
        return self._Brpp_G_Cg

    @property
    def is_leading_edge(self) -> bool:
        return self._is_leading_edge

    @property
    def is_trailing_edge(self) -> bool:
        return self._is_trailing_edge

    # --- Immutable derived: manual lazy caching ---
    @property
    def rightLeg_G(self) -> np.ndarray:
        """This Panel's right leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's right leg vector,
            which is defined from back to front. The units are in meters.
        """
        if self._rightLeg_G is None:
            self._rightLeg_G = cast(np.ndarray, self._Frpp_G_Cg - self._Brpp_G_Cg)
            self._rightLeg_G.flags.writeable = False
        return self._rightLeg_G

    @property
    def frontLeg_G(self) -> np.ndarray:
        """This Panel's front leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's front leg vector,
            which is defined from right to left. The units are in meters.
        """
        if self._frontLeg_G is None:
            self._frontLeg_G = cast(np.ndarray, self._Flpp_G_Cg - self._Frpp_G_Cg)
            self._frontLeg_G.flags.writeable = False
        return self._frontLeg_G

    @property
    def leftLeg_G(self) -> np.ndarray:
        """This Panel's left leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's left leg vector,
            which is defined from front to back. The units are in meters.
        """
        if self._leftLeg_G is None:
            self._leftLeg_G = cast(np.ndarray, self._Blpp_G_Cg - self._Flpp_G_Cg)
            self._leftLeg_G.flags.writeable = False
        return self._leftLeg_G

    @property
    def backLeg_G(self) -> np.ndarray:
        """This Panel's back leg vector (in geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's back leg vector,
            which is defined from left to right. The units are in meters.
        """
        if self._backLeg_G is None:
            self._backLeg_G = cast(np.ndarray, self._Brpp_G_Cg - self._Blpp_G_Cg)
            self._backLeg_G.flags.writeable = False
        return self._backLeg_G

    @property
    def Frbvp_G_Cg(self) -> np.ndarray:
        """The position of this Panel's front right bound vortex point (in geometry
        axes, relative to the CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front right bound vortex point. The units are in meters.
        """
        if self._Frbvp_G_Cg is None:
            self._Frbvp_G_Cg = self._Brpp_G_Cg + 0.75 * self.rightLeg_G
            self._Frbvp_G_Cg.flags.writeable = False
        return self._Frbvp_G_Cg

    @property
    def Flbvp_G_Cg(self) -> np.ndarray:
        """The position of this Panel's front left bound vortex point (in geometry axes,
        relative to the CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front left bound vortex point. The units are in meters.
        """
        if self._Flbvp_G_Cg is None:
            self._Flbvp_G_Cg = self._Flpp_G_Cg + 0.25 * self.leftLeg_G
            self._Flbvp_G_Cg.flags.writeable = False
        return self._Flbvp_G_Cg

    @property
    def Cpp_G_Cg(self) -> np.ndarray:
        """The position of this Panel's collocation point (in geometry axes, relative to
        the CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            collocation point. The units are in meters.
        """
        if self._Cpp_G_Cg is None:
            # Find the positions of points three quarters of the way down the left and
            # right legs of the Panel (in geometry axes, relative to the CG).
            rightThreeQuarterChord_G_Cg = self._Brpp_G_Cg + 0.25 * self.rightLeg_G
            leftThreeQuarterChord_G_Cg = self._Flpp_G_Cg + 0.75 * self.leftLeg_G

            # Find the vector (in geometry axes) between the points three quarters of
            # the way down the left and right legs of the Panel.
            threeQuarterChord_G = (
                leftThreeQuarterChord_G_Cg - rightThreeQuarterChord_G_Cg
            )

            # Find the collocation point (in geometry axes, relative to the CG), which
            # is halfway between the points three quarters of the way down the left and
            # right legs of the Panel. Then populate the class attribute.
            self._Cpp_G_Cg = rightThreeQuarterChord_G_Cg + 0.5 * threeQuarterChord_G
            self._Cpp_G_Cg.flags.writeable = False
        return self._Cpp_G_Cg

    @property
    def unitNormal_G(self) -> np.ndarray:
        """An estimate of this Panel's unit normal vector (in geometry axes).

        :return: A (3,) ndarray of floats representing an estimate of this Panel's unit
            normal vector. The sign is determined via the right-hand rule given the
            orientation of Panel's leg vectors (front right to front left to back left
            to back right).
        """
        if self._unitNormal_G is None:
            firstDiagonal_G = cast(np.ndarray, self._Frpp_G_Cg - self._Blpp_G_Cg)
            secondDiagonal_G = cast(np.ndarray, self._Flpp_G_Cg - self._Brpp_G_Cg)

            cross_G = cast(np.ndarray, np.cross(firstDiagonal_G, secondDiagonal_G))

            self._unitNormal_G = cast(np.ndarray, cross_G / np.linalg.norm(cross_G))
            self._unitNormal_G.flags.writeable = False
        return self._unitNormal_G

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
        if self._area is None:
            firstDiagonal_G = cast(np.ndarray, self._Frpp_G_Cg - self._Blpp_G_Cg)
            secondDiagonal_G = cast(np.ndarray, self._Flpp_G_Cg - self._Brpp_G_Cg)

            cross_G = cast(np.ndarray, np.cross(firstDiagonal_G, secondDiagonal_G))

            self._area = float(np.linalg.norm(cross_G) / 2)
        return self._area

    @property
    def aspect_ratio(self) -> float:
        """The aspect ratio of this Panel.

        :return: The Panel's aspect ratio, which is defined as the distance between the
            right and left legs' center points divided by the distance between the front
            and back legs' center points.
        """
        if self._aspect_ratio is None:
            frontCenterPoint_G_Cg = self._Frpp_G_Cg + self.frontLeg_G / 2
            leftCenterPoint_G_Cg = self._Flpp_G_Cg + self.leftLeg_G / 2
            backCenterPoint_G_Cg = self._Blpp_G_Cg + self.backLeg_G / 2
            rightCenterPoint_G_Cg = self._Brpp_G_Cg + self.rightLeg_G / 2

            right_left_distance = float(
                np.linalg.norm(rightCenterPoint_G_Cg - leftCenterPoint_G_Cg)
            )
            front_back_distance = float(
                np.linalg.norm(frontCenterPoint_G_Cg - backCenterPoint_G_Cg)
            )

            self._aspect_ratio = right_left_distance / front_back_distance
        return self._aspect_ratio

    # --- Set once: properties with single assignment enforcement ---
    @property
    def Frpp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of the Panel's front right vertex (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the Panel's front
            right vertex (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters. Returns None if not yet set, and
            can only be set once.
        """
        return self._Frpp_GP1_CgP1

    @Frpp_GP1_CgP1.setter
    def Frpp_GP1_CgP1(self, newFrpp_GP1_CgP1: np.ndarray) -> None:
        if self._Frpp_GP1_CgP1 is not None:
            raise AttributeError("Frpp_GP1_CgP1 can only be set once")
        self._Frpp_GP1_CgP1 = newFrpp_GP1_CgP1
        self._Frpp_GP1_CgP1.flags.writeable = False

    @property
    def Flpp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of the Panel's front left vertex (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the Panel's front
            left vertex (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters. Returns None if not yet set, and
            can only be set once.
        """
        return self._Flpp_GP1_CgP1

    @Flpp_GP1_CgP1.setter
    def Flpp_GP1_CgP1(self, newFlpp_GP1_CgP1: np.ndarray) -> None:
        if self._Flpp_GP1_CgP1 is not None:
            raise AttributeError("Flpp_GP1_CgP1 can only be set once")
        self._Flpp_GP1_CgP1 = newFlpp_GP1_CgP1
        self._Flpp_GP1_CgP1.flags.writeable = False

    @property
    def Blpp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of the Panel's back left vertex (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the Panel's back
            left vertex (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters. Returns None if not yet set, and
            can only be set once.
        """
        return self._Blpp_GP1_CgP1

    @Blpp_GP1_CgP1.setter
    def Blpp_GP1_CgP1(self, newBlpp_GP1_CgP1: np.ndarray) -> None:
        if self._Blpp_GP1_CgP1 is not None:
            raise AttributeError("Blpp_GP1_CgP1 can only be set once")
        self._Blpp_GP1_CgP1 = newBlpp_GP1_CgP1
        self._Blpp_GP1_CgP1.flags.writeable = False

    @property
    def Brpp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of the Panel's back right vertex (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of the Panel's back
            right vertex (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters. Returns None if not yet set, and
            can only be set once.
        """
        return self._Brpp_GP1_CgP1

    @Brpp_GP1_CgP1.setter
    def Brpp_GP1_CgP1(self, newBrpp_GP1_CgP1: np.ndarray) -> None:
        if self._Brpp_GP1_CgP1 is not None:
            raise AttributeError("Brpp_GP1_CgP1 can only be set once")
        self._Brpp_GP1_CgP1 = newBrpp_GP1_CgP1
        self._Brpp_GP1_CgP1.flags.writeable = False

    @property
    def is_right_edge(self) -> bool | None:
        """Flags if this Panel is at its parent Wing's right edge.

        :return: True if this Panel is at its parent Wing's right edge, False if not.
            Returns None if not yet set, and can only be set once.
        """
        return self._is_right_edge

    @is_right_edge.setter
    def is_right_edge(self, new_is_right_edge: bool) -> None:
        if self._is_right_edge is not None:
            raise AttributeError("is_right_edge can only be set once")
        self._is_right_edge = new_is_right_edge

    @property
    def is_left_edge(self) -> bool | None:
        """Flags if this Panel is at its parent Wing's left edge.

        :return: True if this Panel is at its parent Wing's left edge, False if not.
            Returns None if not yet set, and can only be set once.
        """
        return self._is_left_edge

    @is_left_edge.setter
    def is_left_edge(self, new_is_left_edge: bool) -> None:
        if self._is_left_edge is not None:
            raise AttributeError("is_left_edge can only be set once")
        self._is_left_edge = new_is_left_edge

    @property
    def local_chordwise_position(self) -> int | None:
        """This Panel's local chordwise position in its parent Wing's mesh.

        :return: The local chordwise position of this Panel. Returns None if not yet
            set, and can only be set once.
        """
        return self._local_chordwise_position

    @local_chordwise_position.setter
    def local_chordwise_position(self, new_local_chordwise_position: int) -> None:
        if self._local_chordwise_position is not None:
            raise AttributeError("local_chordwise_position can only be set once")
        self._local_chordwise_position = new_local_chordwise_position

    @property
    def local_spanwise_position(self) -> int | None:
        """This Panel's local spanwise position in its parent Wing's mesh.

        :return: The local spanwise position of this Panel. Returns None if not yet set,
            and can only be set once.
        """
        return self._local_spanwise_position

    @local_spanwise_position.setter
    def local_spanwise_position(self, new_local_spanwise_position: int) -> None:
        if self._local_spanwise_position is not None:
            raise AttributeError("local_spanwise_position can only be set once")
        self._local_spanwise_position = new_local_spanwise_position

    # --- Set once derived: manual lazy caching ---
    @property
    def rightLeg_GP1(self) -> np.ndarray | None:
        """This Panel's right leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's right leg vector,
            which is defined from back to front. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self._Frpp_GP1_CgP1 is None or self._Brpp_GP1_CgP1 is None:
            return None

        if self._rightLeg_GP1 is None:
            self._rightLeg_GP1 = cast(
                np.ndarray, self._Frpp_GP1_CgP1 - self._Brpp_GP1_CgP1
            )
            self._rightLeg_GP1.flags.writeable = False
        return self._rightLeg_GP1

    @property
    def frontLeg_GP1(self) -> np.ndarray | None:
        """This Panel's front leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's left leg vector,
            which is defined from right to left. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self._Flpp_GP1_CgP1 is None or self._Frpp_GP1_CgP1 is None:
            return None

        if self._frontLeg_GP1 is None:
            self._frontLeg_GP1 = cast(
                np.ndarray, self._Flpp_GP1_CgP1 - self._Frpp_GP1_CgP1
            )
            self._frontLeg_GP1.flags.writeable = False
        return self._frontLeg_GP1

    @property
    def leftLeg_GP1(self) -> np.ndarray | None:
        """This Panel's left leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's left leg vector,
            which is defined from front to back. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self._Blpp_GP1_CgP1 is None or self._Flpp_GP1_CgP1 is None:
            return None

        if self._leftLeg_GP1 is None:
            self._leftLeg_GP1 = cast(
                np.ndarray, self._Blpp_GP1_CgP1 - self._Flpp_GP1_CgP1
            )
            self._leftLeg_GP1.flags.writeable = False
        return self._leftLeg_GP1

    @property
    def backLeg_GP1(self) -> np.ndarray | None:
        """This Panel's back leg vector (in the first Airplane's geometry axes).

        :return: A (3,) ndarray of floats representing this Panel's back leg vector,
            which is defined from left to right. The units are in meters. Returns None
            if this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        if self._Brpp_GP1_CgP1 is None or self._Blpp_GP1_CgP1 is None:
            return None

        if self._backLeg_GP1 is None:
            self._backLeg_GP1 = cast(
                np.ndarray, self._Brpp_GP1_CgP1 - self._Blpp_GP1_CgP1
            )
            self._backLeg_GP1.flags.writeable = False
        return self._backLeg_GP1

    @property
    def Frbvp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of this Panel's front right bound vortex point (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front right bound vortex point. The units are in meters. Returns None if
            this Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        # No need to check if self._Brpp_GP1_CgP1 is None, as this is already done by
        # calling self.rightLeg_GP1.
        if self.rightLeg_GP1 is None:
            return None

        if self._Frbvp_GP1_CgP1 is None:
            self._Frbvp_GP1_CgP1 = cast(np.ndarray, self._Brpp_GP1_CgP1) + 0.75 * cast(
                np.ndarray, self._rightLeg_GP1
            )
            self._Frbvp_GP1_CgP1.flags.writeable = False
        return self._Frbvp_GP1_CgP1

    @property
    def Flbvp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of this Panel's front left bound vortex point (in the first
        Airplane's geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            front left bound vortex point. The units are in meters. Returns None if this
            Panel is not part of a SteadyProblem or UnsteadyProblem.
        """
        # No need to check if self._Flpp_GP1_CgP1 is None, as this is already done by
        # calling self.leftLeg_GP1.
        if self.leftLeg_GP1 is None:
            return None

        if self._Flbvp_GP1_CgP1 is None:
            self._Flbvp_GP1_CgP1 = cast(np.ndarray, self._Flpp_GP1_CgP1) + 0.25 * cast(
                np.ndarray, self._leftLeg_GP1
            )
            self._Flbvp_GP1_CgP1.flags.writeable = False
        return self._Flbvp_GP1_CgP1

    @property
    def Cpp_GP1_CgP1(self) -> np.ndarray | None:
        """The position of this Panel's collocation point (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of this Panel's
            collocation point. The units are in meters. Returns None if this Panel is
            not part of a SteadyProblem or UnsteadyProblem.
        """
        # No need to check if self._Brpp_GP1_CgP1 or self._Flpp_GP1_CgP1 is None, as
        # this is already done by calling self.rightLeg_GP1 and self.leftLeg_GP1.
        if self.rightLeg_GP1 is None or self.leftLeg_GP1 is None:
            return None

        if self._Cpp_GP1_CgP1 is None:
            # Find the positions of points three quarters of the way down the left and
            # right legs of the Panel (in the first Airplane's geometry axes, relative
            # to the first Airplane's CG).
            rightThreeQuarterChord_GP1_CgP1 = cast(
                np.ndarray, self._Brpp_GP1_CgP1
            ) + 0.25 * cast(np.ndarray, self._rightLeg_GP1)
            leftThreeQuarterChord_GP1_CgP1 = cast(
                np.ndarray, self._Flpp_GP1_CgP1
            ) + 0.75 * cast(np.ndarray, self._leftLeg_GP1)

            # Find the vector (in the first Airplane's geometry axes) between the points
            # three quarters of the way down the left and right legs of the Panel.
            threeQuarterChord_GP1 = (
                leftThreeQuarterChord_GP1_CgP1 - rightThreeQuarterChord_GP1_CgP1
            )

            # Find the collocation point (in the first Airplane's geometry axes,
            # relative to the first Airplane's CG), which is halfway between the points
            # three quarters of the way down the left and right legs of the Panel.
            self._Cpp_GP1_CgP1 = (
                rightThreeQuarterChord_GP1_CgP1 + 0.5 * threeQuarterChord_GP1
            )
            self._Cpp_GP1_CgP1.flags.writeable = False
        return self._Cpp_GP1_CgP1

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
            self._Frpp_GP1_CgP1 is None
            or self._Flpp_GP1_CgP1 is None
            or self._Blpp_GP1_CgP1 is None
            or self._Brpp_GP1_CgP1 is None
        ):
            return None

        if self._unitNormal_GP1 is None:
            # Compute diagonal vectors (in the first Airplane's geometry axes).
            firstDiagonal_GP1 = self._Frpp_GP1_CgP1 - self._Blpp_GP1_CgP1
            secondDiagonal_GP1 = self._Flpp_GP1_CgP1 - self._Brpp_GP1_CgP1

            # Compute the cross product and normalize.
            cross_GP1 = np.cross(firstDiagonal_GP1, secondDiagonal_GP1)

            self._unitNormal_GP1 = cast(
                np.ndarray, cross_GP1 / np.linalg.norm(cross_GP1)
            )
            self._unitNormal_GP1.flags.writeable = False
        return self._unitNormal_GP1

    # --- Non property methods ---
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

        firstDiagonal_G = cast(np.ndarray, self._Frpp_G_Cg - self._Blpp_G_Cg)
        secondDiagonal_G = cast(np.ndarray, self._Flpp_G_Cg - self._Brpp_G_Cg)

        # Find the projections of the first and second diagonal vectors (in geometry
        # axes) onto the plane's unit normal vector.
        projFirstDiagonalOnNormal_G = (
            np.dot(firstDiagonal_G, unitNormal_G) * unitNormal_G
        )
        projSecondDiagonalOnNormal_G = (
            np.dot(secondDiagonal_G, unitNormal_G) * unitNormal_G
        )

        # Find the projection (in geometry axes) of the first and second diagonal
        # vectors onto the plane.
        projFirstDiagonalOnPlane_G = firstDiagonal_G - projFirstDiagonalOnNormal_G
        projSecondDiagonalOnPlane_G = secondDiagonal_G - projSecondDiagonalOnNormal_G

        # The projected area is found by dividing the magnitude of cross product of
        # the diagonal vectors (in geometry axes) by two. Read the area method for a
        # more detailed explanation.
        projDiagonalsOnPlaneCross_G = np.cross(
            projFirstDiagonalOnPlane_G, projSecondDiagonalOnPlane_G
        )
        return float(np.linalg.norm(projDiagonalsOnPlaneCross_G) / 2)
