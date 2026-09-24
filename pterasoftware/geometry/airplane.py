"""Contains the Airplane class."""

from __future__ import annotations

import copy
import time
from collections.abc import Sequence
from typing import Any, cast

import numpy as np
import pyvista as pv
import scipy.interpolate as sp_interp
import webp

from .. import _parameter_validation, _transformations
from . import wing as wing_mod
from . import wing_cross_section as wing_cross_section_mod

# The relative tolerance used by the checks on the projected planform that sets the
# default reference dimensions. A strip is edge-on when its projected area is at most
# this fraction of the magnitude of its vector area, and two triangles overlap when
# their intersection area exceeds this fraction of the smaller triangle's area.
_PLANFORM_RELATIVE_TOLERANCE = 1e-9

# The number of evenly spaced points at which the planform samples an edge_defined
# Wing's edge curves between each neighboring pair of stored points.
_NUM_EDGE_CURVE_SAMPLES_BETWEEN_POINTS = 32


class Airplane:
    """A class used to contain airplanes.

    The Airplane class is responsible for: (1) Defining the local body axes and geometry
    axes, (2) managing Wings and their coordinate transformations, (3) processing
    symmetric Wings and converting them to separate wings when the symmetry plane is not
    coincident with the Wing's axes xz plane (type 5 symmetry), and (4) providing
    reference dimensions for aerodynamic calculations.

    Every Airplane has a body axis system, where +x points forward along fuselage, +y
    points to the right (starboard direction), and +z points downward (completing a
    right-handed system).

    Every Airplane also has a geometry axis system, where +x points aft along fuselage,
    +y points to the right (starboard direction), and +z points upward (completing a
    right-handed system).

    The reference planform supplies the default reference dimensions. It is the planform
    of the first element of the wings passed to the Airplane, including its mirrored
    half if it has type 4 symmetry, or the reflected Wing that symmetry processing
    creates for it if it has type 5 symmetry. It is built from the set geometry, in
    geometry axes, as a set of strips. For a Wing built from WingCrossSections, each
    strip joins two neighboring WingCrossSections' leading points and undeflected
    trailing points, so control surface deflections have no effect. An edge_defined
    Wing's whole half is one strip, bounded by its stored, untrimmed leading and
    trailing edge curves. For type 5 symmetry, one more strip joins the two halves' root
    chords, filling any gap between them. The strips are then projected onto the
    geometry axes' xy plane. The default s_ref is the sum of the strips' projected
    areas, the default b_ref is the extent of all the strips' points along the geometry
    axes' y axis, and the default c_ref is the mean aerodynamic chord: the integral of
    the square of the projected chord along the geometry axes' y axis, divided by the
    default s_ref. A strip that is edge-on in the projection, such as a vertical
    winglet, contributes nothing to s_ref or c_ref.

    Computing the defaults raises a ValueError in two cases. The first is when the
    reference planform is too steeply inclined: its area projected onto the geometry
    axes' xy plane is less than 1 / sqrt(2) times its area projected onto the xy plane
    of the wing axes of the first element of wings (for type 5 symmetry, both areas
    leave out the strip joining the halves). The second is when the reference planform's
    projection is ill-formed: a strip crosses over itself, two strips overlap, or, for
    type 5 symmetry, the two halves overlap along the geometry axes' y axis. In either
    case, s_ref, c_ref, and b_ref must all be passed explicitly.

    Immutable attributes (wings, name, Cg_GP1_CgP1, weight, s_ref, c_ref, and b_ref) are
    set during initialization and cannot be modified afterward. The numpy array
    Cg_GP1_CgP1 is made read only to prevent in place mutation. The wings attribute is
    stored as a tuple to prevent external mutation.

    Derived properties (num_panels and T_pas_G_Cg_to_GP1_CgP1) are lazily evaluated and
    cached since they depend only on immutable attributes.

    The load attributes remain mutable as they are set by the solver during simulation.
    The forces_W, forceCoefficients_W, moments_W_CgP1, and momentCoefficients_W_CgP1
    attributes hold the total force (in wind axes), the total moment (in wind axes,
    relative to the first Airplane's CG), and their coefficients. The forces_G,
    forceCoefficients_G, moments_G_Cg, momentCoefficients_G_Cg, moments_W_Cg, and
    momentCoefficients_W_Cg attributes hold the total force (in geometry axes), the
    total moment (in geometry axes, relative to the CG), the total moment (in wind axes,
    relative to the CG), and their coefficients. Because every Airplane's geometry axes
    are parallel to the first Airplane's geometry axes, the components of forces_G equal
    the components of the total force (in the first Airplane's geometry axes). For the
    first Airplane, moments_W_Cg equals moments_W_CgP1, and momentCoefficients_W_Cg
    equals momentCoefficients_W_CgP1.

    **Citations:**

    Adapted from: geometry.Airplane in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/23/2020

    Reference area, span, and mean aerodynamic chord conventions adapted from: Section
    6.1 of "Flight Vehicle Aerodynamics" (2014)

    Author: Mark Drela

    Date of retrieval: 09/23/2026

    Reference area and span conventions adapted from: Section 2.2 of "Aircraft Control
    and Simulation" (third edition, 2016)

    Authors: Brian L. Stevens, Frank L. Lewis, and Eric N. Johnson

    Date of retrieval: 09/23/2026
    """

    __slots__ = (
        "_wings",
        "_name",
        "_Cg_GP1_CgP1",
        "_weight",
        "_s_ref",
        "_c_ref",
        "_b_ref",
        "_num_panels",
        "_T_pas_G_Cg_to_GP1_CgP1",
        "forces_W",
        "forceCoefficients_W",
        "moments_W_CgP1",
        "momentCoefficients_W_CgP1",
        "forces_G",
        "forceCoefficients_G",
        "moments_G_Cg",
        "momentCoefficients_G_Cg",
        "moments_W_Cg",
        "momentCoefficients_W_Cg",
    )

    def __init__(
        self,
        wings: list[wing_mod.Wing],
        name: str = "Untitled Airplane",
        Cg_GP1_CgP1: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        weight: float | int = 0.0,
        s_ref: float | int | None = None,
        c_ref: float | int | None = None,
        b_ref: float | int | None = None,
    ) -> None:
        """The initialization method.

        :param wings: A list of the airplane's wings defined as Wings. It must contain
            at least one Wing. Wings with symmetric=True and non coincident symmetry
            planes will be automatically processed into separate Wings during
            initialization (type 5 symmetry).
        :param name: A sensible name for your airplane. The default is "Untitled
            Airplane".
        :param Cg_GP1_CgP1: An array-like object of 3 numbers representing the position
            of this Airplane's CG (in the first Airplane's geometry axes, relative to
            the first Airplane's CG). Can be a list, tuple, or ndarray. Values are
            converted to floats internally. For the first Airplane in a simulation, this
            must be equivalent to (0.0, 0.0, 0.0) by definition. The units are in
            meters. The default is (0.0, 0.0, 0.0).
        :param weight: A number (int or float) representing the weight of the aircraft
            in Newtons. This is used by the trim functions. It must be greater than or
            equal to zero. The default is 0.0. In free flight, it must also be
            consistent with the FreeFlightUnsteadyProblem's mass and the
            OperatingPoint's gravitational acceleration, satisfying weight == mass *
            np.linalg.norm(g_E) within floating point tolerance.
        :param s_ref: A number (int or float) representing the reference area. If not
            set or set to None (the default), it populates from the projected area of
            the reference planform, as described in the class docstring. If set, it must
            be greater than zero, and will be converted to a float internally. The units
            are square meters. Airplanes derived from this one (those at each time step
            of an unsteady simulation, those refined by a convergence analysis, and
            those tried by a trim analysis) inherit the resulting value rather than
            recomputing it from their own Wings, so all of their load coefficients share
            one normalization.
        :param c_ref: A number (int or float) representing the reference chord length.
            If not set or set to None (the default), it populates from the mean
            aerodynamic chord of the projected reference planform, as described in the
            class docstring. If set, it must be greater than zero, and will be converted
            to a float internally. The units are meters. Derived Airplanes inherit the
            resulting value in the same way as s_ref.
        :param b_ref: A number (int or float) representing the reference span. If not
            set or set to None (the default value), it populates from the span of the
            projected reference planform, as described in the class docstring. If set,
            it must be greater than zero, and will be converted to a float internally.
            The units are meters. Derived Airplanes inherit the resulting value in the
            same way as s_ref.
        """
        # Initialize the immutable attributes. Set those that are numpy arrays to be
        # read only. Store wings as a tuple to prevent external mutation.
        wings = _parameter_validation.non_empty_list_return_list(wings, "wings")
        processed_wings: list[wing_mod.Wing] = []
        first_wings: list[wing_mod.Wing] = []
        for wing in wings:
            if not isinstance(wing, wing_mod.Wing):
                raise TypeError("Every element in wings must be a Wing")
            wing_symmetry_wings = self.process_wing_symmetry(wing)
            if not first_wings:
                first_wings = wing_symmetry_wings
            processed_wings.extend(wing_symmetry_wings)
        self._wings = tuple(processed_wings)

        self._name = _parameter_validation.str_return_str(name, "name")

        self._Cg_GP1_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            Cg_GP1_CgP1, "Cg_GP1_CgP1"
        )
        self._Cg_GP1_CgP1.flags.writeable = False

        self._weight = _parameter_validation.number_in_range_return_float(
            weight,
            "weight",
            min_val=0.0,
            min_inclusive=True,
        )

        # If any of the passed reference dimensions are None, set them to the
        # corresponding dimension of the projected reference planform. Otherwise, set
        # them to the passed dimension after checking that it is valid.
        if s_ref is None or c_ref is None or b_ref is None:
            default_s_ref, default_c_ref, default_b_ref = (
                _get_planform_reference_dimensions(first_wings)
            )
        if s_ref is None:
            self._s_ref = default_s_ref
        else:
            self._s_ref = _parameter_validation.number_in_range_return_float(
                s_ref, "s_ref", min_val=0.0, min_inclusive=False
            )
        if c_ref is None:
            self._c_ref = default_c_ref
        else:
            self._c_ref = _parameter_validation.number_in_range_return_float(
                c_ref, "c_ref", min_val=0.0, min_inclusive=False
            )
        if b_ref is None:
            self._b_ref = default_b_ref
        else:
            self._b_ref = _parameter_validation.number_in_range_return_float(
                b_ref, "b_ref", min_val=0.0, min_inclusive=False
            )

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._num_panels: int | None = None
        self._T_pas_G_Cg_to_GP1_CgP1: np.ndarray | None = None

        # Initialize mutable attributes to hold the loads and load coefficients this
        # Airplane experiences.
        self.forces_W: np.ndarray | None = None
        self.forceCoefficients_W: np.ndarray | None = None
        self.moments_W_CgP1: np.ndarray | None = None
        self.momentCoefficients_W_CgP1: np.ndarray | None = None
        self.forces_G: np.ndarray | None = None
        self.forceCoefficients_G: np.ndarray | None = None
        self.moments_G_Cg: np.ndarray | None = None
        self.momentCoefficients_G_Cg: np.ndarray | None = None
        self.moments_W_Cg: np.ndarray | None = None
        self.momentCoefficients_W_Cg: np.ndarray | None = None

    # --- Deep copy methods ---
    def __deepcopy__(self, memo: dict) -> Airplane:
        """Creates a deep copy of this Airplane, preserving mesh geometry but resetting
        solver state.

        The copy preserves: (1) Wings tuple (each Wing is deep copied, preserving mesh
        and Panels) (2) Airplane parameters (name, Cg_GP1_CgP1, weight, reference
        dimensions), and (3) cached derived properties (num_panels,
        T_pas_G_Cg_to_GP1_CgP1).

        The copy resets to None: (1) loads and load coefficients.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new Airplane with preserved mesh geometry and reset solver state.
        """
        # Create a new Airplane instance without calling __init__ to avoid redundant
        # validation and Wing symmetry processing.
        new_airplane = object.__new__(Airplane)

        # Store this Airplane in memo to handle potential circular references.
        memo[id(self)] = new_airplane

        # Deep copy the Wings into a new tuple.
        new_airplane._wings = tuple(copy.deepcopy(wing, memo) for wing in self._wings)

        # Copy immutable attributes. For those that are numpy arrays, make the copies
        # read only.
        new_airplane._name = self._name
        new_airplane._Cg_GP1_CgP1 = self._Cg_GP1_CgP1.copy()
        new_airplane._Cg_GP1_CgP1.flags.writeable = False
        new_airplane._weight = self._weight
        new_airplane._s_ref = self._s_ref
        new_airplane._c_ref = self._c_ref
        new_airplane._b_ref = self._b_ref

        # Copy cached derived properties. For those that are numpy arrays, make the
        # copies read only.
        new_airplane._num_panels = self._num_panels
        if self._T_pas_G_Cg_to_GP1_CgP1 is not None:
            new_airplane._T_pas_G_Cg_to_GP1_CgP1 = self._T_pas_G_Cg_to_GP1_CgP1.copy()
            new_airplane._T_pas_G_Cg_to_GP1_CgP1.flags.writeable = False
        else:
            new_airplane._T_pas_G_Cg_to_GP1_CgP1 = None

        # Reset loads and load coefficients to None (solver will compute these).
        new_airplane.forces_W = None
        new_airplane.forceCoefficients_W = None
        new_airplane.moments_W_CgP1 = None
        new_airplane.momentCoefficients_W_CgP1 = None
        new_airplane.forces_G = None
        new_airplane.forceCoefficients_G = None
        new_airplane.moments_G_Cg = None
        new_airplane.momentCoefficients_G_Cg = None
        new_airplane.moments_W_Cg = None
        new_airplane.momentCoefficients_W_Cg = None

        return new_airplane

    def deep_copy_with_Cg_GP1_CgP1(
        self, new_Cg_GP1_CgP1: np.ndarray | Sequence[float | int]
    ) -> Airplane:
        """Creates a deep copy of this Airplane with a different Cg_GP1_CgP1 position.

        This method is used by AirplaneMovement to create Airplanes at different time
        steps that share the same geometry but have different positions in the
        formation. It maintains immutability by returning a new Airplane rather than
        modifying the existing one. The Wings and their child objects are deep copied
        unchanged, since their geometry is defined relative to this Airplane's own CG
        rather than its position in the formation.

        :param new_Cg_GP1_CgP1: An array-like object of 3 numbers representing the
            position of the new Airplane's CG (in the first Airplane's geometry axes,
            relative to the first Airplane's CG). Can be a list, tuple, or ndarray.
            Values are converted to floats internally. The units are in meters.
        :return: A new Airplane with the specified position and deep copied geometry.
        """
        # Only Cg_GP1_CgP1 and its derived cache (_T_pas_G_Cg_to_GP1_CgP1) need to
        # differ from a standard deep copy because (1) Wing geometry (Ler_Gs_Cgs,
        # panels, etc.) is defined relative to this Airplane's own CG, not the formation
        # position, so it remains valid, (2) Panel local coordinates (_G_Cg) are
        # independent of formation position (global coordinates (_GP1_CgP1) are reset to
        # None by Panel's __deepcopy__ and will be recomputed by the Problem using the
        # new transformation matrix), and (3) all other child objects
        # (WingCrossSections, Airfoils, vortices) have no dependency on Cg_GP1_CgP1.

        # Validate the new position.
        validated_Cg_GP1_CgP1 = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                new_Cg_GP1_CgP1, "new_Cg_GP1_CgP1"
            )
        )
        validated_Cg_GP1_CgP1.flags.writeable = False

        # Create a new Airplane instance without calling __init__ to avoid redundant
        # validation and Wing symmetry processing.
        new_airplane = object.__new__(Airplane)

        # Deep copy the Wings into a new tuple.
        memo: dict = {id(self): new_airplane}
        new_airplane._wings = tuple(copy.deepcopy(wing, memo) for wing in self._wings)

        # Copy immutable attributes, using the new Cg_GP1_CgP1.
        new_airplane._name = self._name
        new_airplane._Cg_GP1_CgP1 = validated_Cg_GP1_CgP1
        new_airplane._weight = self._weight
        new_airplane._s_ref = self._s_ref
        new_airplane._c_ref = self._c_ref
        new_airplane._b_ref = self._b_ref

        # Copy _num_panels cache (depends only on Wings, not position). Reset
        # _T_pas_G_Cg_to_GP1_CgP1 to None (depends on Cg_GP1_CgP1).
        new_airplane._num_panels = self._num_panels
        new_airplane._T_pas_G_Cg_to_GP1_CgP1 = None

        # Reset loads and load coefficients to None (solver will compute these).
        new_airplane.forces_W = None
        new_airplane.forceCoefficients_W = None
        new_airplane.moments_W_CgP1 = None
        new_airplane.momentCoefficients_W_CgP1 = None
        new_airplane.forces_G = None
        new_airplane.forceCoefficients_G = None
        new_airplane.moments_G_Cg = None
        new_airplane.momentCoefficients_G_Cg = None
        new_airplane.moments_W_Cg = None
        new_airplane.momentCoefficients_W_Cg = None

        return new_airplane

    # --- Immutable: read only properties ---
    @property
    def wings(self) -> tuple[wing_mod.Wing, ...]:
        return self._wings

    @property
    def name(self) -> str:
        return self._name

    @property
    def Cg_GP1_CgP1(self) -> np.ndarray:
        return self._Cg_GP1_CgP1

    @property
    def weight(self) -> float:
        return self._weight

    @property
    def s_ref(self) -> float:
        return self._s_ref

    @property
    def c_ref(self) -> float:
        return self._c_ref

    @property
    def b_ref(self) -> float:
        return self._b_ref

    # --- Immutable derived: manual lazy caching ---
    @property
    def num_panels(self) -> int:
        """The total number of Panels across all Wings.

        :return: The total number of Panels.
        """
        if self._num_panels is None:
            self._num_panels = sum(
                wing.num_panels if wing.num_panels is not None else 0
                for wing in self._wings
            )
        return self._num_panels

    @property
    def T_pas_G_Cg_to_GP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix from this Airplane's geometry axes,
        relative to this Airplane's CG to the first Airplane's geometry axes, relative
        to the first Airplane's CG.

        Computes the transformation chain: G_Cg > GP1_CgP1. This transformation matrix
        is used to position Airplanes relative to one another, in problems with more
        than one Airplane. If this Airplane is the first Airplane (where Cg_GP1_CgP1 =
        [0, 0, 0]), it returns an identity transformation.

        :return: A (4,4) ndarray of floats representing the passive transformation
            matrix from this Airplane's geometry axes, relative to its CG to the first
            Airplane's geometry axes, relative to its CG.
        """
        if self._T_pas_G_Cg_to_GP1_CgP1 is None:
            # generate_trans_T with passive=True expects the translations parameter to
            # be the position of the target reference point (CgP1) relative to the
            # source reference point (Cg). Using the notation from
            # AXES_POINTS_AND_FRAMES.md: translations = CgP1_G_Cg. However, we have
            # Cg_GP1_CgP1 (position of Cg, in GP1 axes, relative to CgP1). Since
            # geometry axes G and GP1 are parallel (pure translation, no rotation):
            # CgP1_G_Cg = -Cg_GP1_CgP1.
            self._T_pas_G_Cg_to_GP1_CgP1 = _transformations.generate_trans_T(
                translations=-self._Cg_GP1_CgP1, passive=True
            )
            self._T_pas_G_Cg_to_GP1_CgP1.flags.writeable = False
        return self._T_pas_G_Cg_to_GP1_CgP1

    # --- Mutable derived: read only properties, no backing slots ---
    @property
    def inducedDrag_W(self) -> float | None:
        """The total induced drag force experienced by this Airplane (in wind axes).

        Induced drag points along the wind axes' -x basis direction, so it is the
        negative of the wind axes' x force component.

        :return: The induced drag force in Newtons, or None if forces_W has not been
            set.
        """
        if self.forces_W is None:
            return None
        return float(-self.forces_W[0])

    @property
    def crosswindForce_W(self) -> float | None:
        """The total crosswind force experienced by this Airplane (in wind axes).

        Crosswind force points along the wind axes' -y basis direction, so it is the
        negative of the wind axes' y force component.

        :return: The crosswind force in Newtons, or None if forces_W has not been set.
        """
        if self.forces_W is None:
            return None
        return float(-self.forces_W[1])

    @property
    def lift_W(self) -> float | None:
        """The total lift force experienced by this Airplane (in wind axes).

        Lift points along the wind axes' -z basis direction, so it is the negative of
        the wind axes' z force component.

        :return: The lift force in Newtons, or None if forces_W has not been set.
        """
        if self.forces_W is None:
            return None
        return float(-self.forces_W[2])

    @property
    def inducedDragCoefficient_W(self) -> float | None:
        """The total induced drag force coefficient experienced by this Airplane (in
        wind axes).

        Induced drag coefficient corresponds to the wind axes' -x basis direction, so it
        is the negative of the wind axes' x force coefficient component.

        :return: The induced drag coefficient, or None if forceCoefficients_W has not
            been set.
        """
        if self.forceCoefficients_W is None:
            return None
        return float(-self.forceCoefficients_W[0])

    @property
    def crosswindForceCoefficient_W(self) -> float | None:
        """The total crosswind force coefficient experienced by this Airplane (in wind
        axes).

        Crosswind force coefficient corresponds to the wind axes' -y basis direction, so
        it is the negative of the wind axes' y force coefficient component.

        :return: The crosswind force coefficient, or None if forceCoefficients_W has not
            been set.
        """
        if self.forceCoefficients_W is None:
            return None
        return float(-self.forceCoefficients_W[1])

    @property
    def liftCoefficient_W(self) -> float | None:
        """The total lift force coefficient experienced by this Airplane (in wind axes).

        Lift coefficient corresponds to the wind axes' -z basis direction, so it is the
        negative of the wind axes' z force coefficient component.

        :return: The lift coefficient, or None if forceCoefficients_W has not been set.
        """
        if self.forceCoefficients_W is None:
            return None
        return float(-self.forceCoefficients_W[2])

    @property
    def rollingMoment_W_Cg(self) -> float | None:
        """The total rolling moment experienced by this Airplane (in wind axes, relative
        to the CG).

        Rolling moment acts about the wind axes' +x basis direction, so it equals the
        wind axes' x moment component.

        :return: The rolling moment in Newton-meters, or None if moments_W_Cg has not
            been set.
        """
        if self.moments_W_Cg is None:
            return None
        return float(self.moments_W_Cg[0])

    @property
    def pitchingMoment_W_Cg(self) -> float | None:
        """The total pitching moment experienced by this Airplane (in wind axes,
        relative to the CG).

        Pitching moment acts about the wind axes' +y basis direction, so it equals the
        wind axes' y moment component.

        :return: The pitching moment in Newton-meters, or None if moments_W_Cg has not
            been set.
        """
        if self.moments_W_Cg is None:
            return None
        return float(self.moments_W_Cg[1])

    @property
    def yawingMoment_W_Cg(self) -> float | None:
        """The total yawing moment experienced by this Airplane (in wind axes, relative
        to the CG).

        Yawing moment acts about the wind axes' +z basis direction, so it equals the
        wind axes' z moment component.

        :return: The yawing moment in Newton-meters, or None if moments_W_Cg has not
            been set.
        """
        if self.moments_W_Cg is None:
            return None
        return float(self.moments_W_Cg[2])

    @property
    def rollingMomentCoefficient_W_Cg(self) -> float | None:
        """The total rolling moment coefficient experienced by this Airplane (in wind
        axes, relative to the CG).

        Rolling moment coefficient corresponds to the wind axes' +x basis direction, so
        it equals the wind axes' x moment coefficient component.

        :return: The rolling moment coefficient, or None if momentCoefficients_W_Cg has
            not been set.
        """
        if self.momentCoefficients_W_Cg is None:
            return None
        return float(self.momentCoefficients_W_Cg[0])

    @property
    def pitchingMomentCoefficient_W_Cg(self) -> float | None:
        """The total pitching moment coefficient experienced by this Airplane (in wind
        axes, relative to the CG).

        Pitching moment coefficient corresponds to the wind axes' +y basis direction, so
        it equals the wind axes' y moment coefficient component.

        :return: The pitching moment coefficient, or None if momentCoefficients_W_Cg has
            not been set.
        """
        if self.momentCoefficients_W_Cg is None:
            return None
        return float(self.momentCoefficients_W_Cg[1])

    @property
    def yawingMomentCoefficient_W_Cg(self) -> float | None:
        """The total yawing moment coefficient experienced by this Airplane (in wind
        axes, relative to the CG).

        Yawing moment coefficient corresponds to the wind axes' +z basis direction, so
        it equals the wind axes' z moment coefficient component.

        :return: The yawing moment coefficient, or None if momentCoefficients_W_Cg has
            not been set.
        """
        if self.momentCoefficients_W_Cg is None:
            return None
        return float(self.momentCoefficients_W_Cg[2])

    # --- Other methods ---
    def draw(
        self, save: bool | np.bool = False, testing: bool | np.bool = False
    ) -> None:
        """Draws the 3D geometry of this Airplane.

        This method provides a convenient way to visualize the Airplane's Panels without
        needing to create a solver object first. It shows the Panel's surfaces in 3D
        using PyVista.

        :param save: Set to True to save the image as a WebP. Can be a bool or a numpy
            bool and will be converted internally to bool. The default value is False.
        :param testing: Set to True to close the image after 1 second, which is useful
            for running test suites. Can be a bool or a numpy bool and will be converted
            internally to bool. The default value is False.
        :return: None
        """
        save = _parameter_validation.boolLike_return_bool(save, "save")
        testing = _parameter_validation.boolLike_return_bool(testing, "testing")

        # Define visualization constants.
        panel_color = "chartreuse"
        plotter_background_color = "black"
        window_size = [1024, 768]
        quality = 75

        # Initialize the plotter and set it to use parallel projection.
        plotter = pv.Plotter(window_size=window_size, lighting=None)
        plotter.enable_parallel_projection()  # type: ignore[call-arg]

        # Initialize empty arrays to hold the Panels' vertices and faces.
        panel_vertices = np.empty((0, 3), dtype=float)
        panel_faces = np.empty(0, dtype=int)

        # Initialize a variable to keep track of how many Panels' data has been added to
        # the ndarrays.
        panel_num = 0

        # Iterate through this Airplane's Wings.
        for wing in self._wings:
            # Unravel the Wing's Panel matrix and iterate through it.
            if wing.panels is None:
                continue
            panels = np.ravel(wing.panels)
            for panel in panels:
                # Stack this Panel's vertices and faces.
                panel_vertices_to_add = np.vstack(
                    (
                        panel.Flpp_G_Cg,
                        panel.Frpp_G_Cg,
                        panel.Brpp_G_Cg,
                        panel.Blpp_G_Cg,
                    )
                )
                panel_face_to_add = np.array(
                    [
                        4,
                        (panel_num * 4),
                        (panel_num * 4) + 1,
                        (panel_num * 4) + 2,
                        (panel_num * 4) + 3,
                    ]
                )

                # Stack this Panel's vertices and faces with the array of all vertices
                # and faces.
                panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                panel_faces = np.hstack((panel_faces, panel_face_to_add))

                # Update the number of previous Panels.
                panel_num += 1

        # Convert the Panel vertices and faces to PolyData.
        panel_surfaces = pv.PolyData(panel_vertices, panel_faces)

        # Add the Panels to the plotter.
        plotter.add_mesh(
            panel_surfaces,
            show_edges=True,
            color=panel_color,
            smooth_shading=False,
        )

        # Set the plotter's background color.
        plotter.set_background(color=plotter_background_color)  # type: ignore[call-arg]

        if not testing:
            # Show the plotter so the user can adjust the camera position and window
            plotter.show(
                title=f"Airplane: {self._name}",
                cpos=(-1, -1, 1),
                full_screen=False,
                auto_close=False,
            )
        else:
            # Show the plotter for 1 second, then proceed automatically (for testing).
            plotter.show(
                title=f"Airplane: {self._name}",
                cpos=(-1, -1, 1),
                full_screen=False,
                interactive=False,
                auto_close=False,
            )
            time.sleep(1)

        # If the user wants to save the image, take a screenshot and save as WebP.
        if save:
            screenshot = plotter.screenshot(
                filename=None,
                transparent_background=True,
                return_img=True,
            )
            image = webp.Image.fromarray(
                cast(np.ndarray[Any, Any], screenshot),
            )
            # The compression method matches _output_rendering.WEBP_METHOD, which this
            # module cannot import without a circular import. It trades file size for
            # encode time and leaves the image quality alone.
            webp.save_image(
                img=image,
                file_path=f"{self._name.lower().replace(' ', '_')}_geometry.webp",
                lossless=False,
                quality=quality,
                method=0,
            )

        # Close all the plotters.
        pv.close_all()

    def get_plottable_data(
        self, show: bool | np.bool = False
    ) -> list[list[list[np.ndarray]]] | None:
        """Returns plottable data for this Airplane's Airfoils' outlines and mean camber
        lines.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: If show is True, returns None. If show is False, returns a list of sub
            lists (one sub list for each of this Airplane's Wings). Each sub list
            contains sub sub lists (one for each of this Wing's WingCrossSections). Each
            sub sub list contains two ndarrays. The first ndarray contains points on
            that WingCrossSection's Airfoil's outline and the second contains points on
            its mean camber line. The points are in geometry axes, relative to the CG.
            The units are in meters.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        airfoilOutlines_G_Cg = []
        airfoilMcls_G_Cg = []
        for wing_id, wing in enumerate(self._wings):
            plottable_data = wing.get_plottable_data(show=False)

            assert plottable_data is not None
            [airfoilOutlines_Wn_Ler, airfoilMcls_Wn_Ler] = plottable_data

            these_airfoilOutlines_G_Cg = []
            these_airfoilMcls_G_Cg = []
            for airfoil_id in range(len(airfoilOutlines_Wn_Ler)):
                airfoilOutline_Wn_Ler = airfoilOutlines_Wn_Ler[airfoil_id]
                airfoilMcl_Wn_Ler = airfoilMcls_Wn_Ler[airfoil_id]

                assert wing.T_pas_Wn_Ler_to_G_Cg is not None
                airfoilOutline_G_Cg = _transformations.apply_T_to_vectors(
                    wing.T_pas_Wn_Ler_to_G_Cg, airfoilOutline_Wn_Ler, is_position=True
                )
                airfoilMcl_G_Cg = _transformations.apply_T_to_vectors(
                    wing.T_pas_Wn_Ler_to_G_Cg, airfoilMcl_Wn_Ler, is_position=True
                )

                these_airfoilOutlines_G_Cg.append(airfoilOutline_G_Cg)
                these_airfoilMcls_G_Cg.append(airfoilMcl_G_Cg)

            airfoilOutlines_G_Cg.append(these_airfoilOutlines_G_Cg)
            airfoilMcls_G_Cg.append(these_airfoilMcls_G_Cg)

        if not show:
            return [airfoilOutlines_G_Cg, airfoilMcls_G_Cg]

        plotter = pv.Plotter()

        AxesGCg = pv.AxesAssembly(
            x_label="GX@Cg",
            y_label="GY@Cg",
            z_label="GZ@Cg",
            # labels=None,
            label_color="black",
            show_labels=True,
            # label_position=(1, 1, 1),
            label_size=15,
            x_color="red",
            y_color="green",
            z_color="blue",
            # position=(0.0, 0.0, 0.0),
            # orientation=(0.0, 0.0, 0.0),
            # origin=(0.0, 0.0, 0.0),
            scale=(0.25, 0.25, 0.25),
            user_matrix=np.eye(4, dtype=float),
            name="G",
            shaft_type="cylinder",
            shaft_radius=0.025,
            shaft_length=(0.8, 0.8, 0.8),
            tip_type="cone",
            tip_radius=0.1,
            tip_length=(0.2, 0.2, 0.2),
            symmetric_bounds=False,
        )

        plotter.add_actor(AxesGCg)

        for wing_id, wing in enumerate(self._wings):
            wing_num = wing_id + 1

            assert wing.T_pas_G_Cg_to_Wn_Ler is not None
            AxesWLerWcs1Lp1_G_Cg = pv.AxesAssembly(
                x_label=f"W{wing_num}X@Ler/Wcs1XLp1",
                y_label=f"W{wing_num}Y@Ler/Wcs1YLp1",
                z_label=f"W{wing_num}Z@Ler/Wcs1ZLp1",
                # labels=None,
                label_color="black",
                show_labels=True,
                # label_position=(1, 1, 1),
                label_size=15,
                x_color="red",
                y_color="green",
                z_color="blue",
                # position=(0.0, 0.0, 0.0),
                # orientation=(0.0, 0.0, 0.0),
                # origin=(0.0, 0.0, 0.0),
                scale=(0.25, 0.25, 0.25),
                user_matrix=np.linalg.inv(wing.T_pas_G_Cg_to_Wn_Ler),
                # user_matrix=wingAxes_T_act,
                name=f"W{wing_num}/Wcs1",
                shaft_type="cylinder",
                shaft_radius=0.025,
                shaft_length=(0.8, 0.8, 0.8),
                tip_type="cone",
                tip_radius=0.1,
                tip_length=(0.2, 0.2, 0.2),
                symmetric_bounds=False,
            )

            plotter.add_actor(AxesWLerWcs1Lp1_G_Cg)

            these_airfoilOutlines_G_Cg = airfoilOutlines_G_Cg[wing_id]
            these_airfoilMcls_G_Cg = airfoilMcls_G_Cg[wing_id]

            for wing_cross_section_id, wing_cross_section in enumerate(
                wing.wing_cross_sections
            ):
                airfoilOutline_G_Cg = these_airfoilOutlines_G_Cg[wing_cross_section_id]
                airfoilMcl_G_Cg = these_airfoilMcls_G_Cg[wing_cross_section_id]

                airfoilOutline_faces = np.hstack(
                    [
                        airfoilOutline_G_Cg.shape[0],
                        np.arange(airfoilOutline_G_Cg.shape[0]),
                    ]
                )
                airfoilOutline_mesh = pv.PolyData(
                    airfoilOutline_G_Cg, faces=airfoilOutline_faces
                )
                plotter.add_mesh(airfoilOutline_mesh)
                plotter.add_lines(airfoilMcl_G_Cg)

                if wing_cross_section_id != 0:
                    wing_cross_section_num = wing_cross_section_id + 1

                    AxesWcsLp_G_Cg = pv.AxesAssembly(
                        x_label=f"Wcs{wing_cross_section_num}Wn{wing_num}X@Lp"
                        f"{wing_cross_section_num}Wn{wing_num}",
                        y_label=f"Wcs{wing_cross_section_num}Wn{wing_num}Y@Lp"
                        f"{wing_cross_section_num}Wn{wing_num}",
                        z_label=f"Wcs{wing_cross_section_num}Wn{wing_num}Z@Lp"
                        f"{wing_cross_section_num}Wn{wing_num}",
                        # labels=None,
                        label_color="black",
                        show_labels=True,
                        # label_position=(1, 1, 1),
                        label_size=15,
                        x_color="red",
                        y_color="green",
                        z_color="blue",
                        # position=(0.0, 0.0, 0.0),
                        # orientation=(0.0, 0.0, 0.0),
                        # origin=(0.0, 0.0, 0.0),
                        scale=(0.25, 0.25, 0.25),
                        user_matrix=np.linalg.inv(
                            wing.children_T_pas_G_Cg_to_Wcs_Lp[wing_cross_section_id]
                        ),
                        name=f"Wcs{wing_cross_section_id}Wn{wing_num}",
                        shaft_type="cylinder",
                        shaft_radius=0.025,
                        shaft_length=(0.8, 0.8, 0.8),
                        tip_type="cone",
                        tip_radius=0.1,
                        tip_length=(0.2, 0.2, 0.2),
                        symmetric_bounds=False,
                    )

                    plotter.add_actor(AxesWcsLp_G_Cg)

            if wing.panels is not None:
                # Initialize empty arrays to hold the Panels' vertices and faces.
                panel_vertices = np.empty((0, 3), dtype=float)
                panel_faces = np.empty(0, dtype=int)

                # Initialize a variable to keep track of how many Panels' data has been
                # added to the arrays.
                panel_num = 0

                # Unravel the Wing's Panel matrix and iterate through it.
                panels = np.ravel(wing.panels)
                for panel in panels:
                    # Stack this Panel's vertices and faces.
                    panel_vertices_to_add = np.vstack(
                        (
                            panel.Flpp_G_Cg,
                            panel.Frpp_G_Cg,
                            panel.Brpp_G_Cg,
                            panel.Blpp_G_Cg,
                        )
                    )
                    panel_face_to_add = np.array(
                        [
                            4,
                            (panel_num * 4),
                            (panel_num * 4) + 1,
                            (panel_num * 4) + 2,
                            (panel_num * 4) + 3,
                        ]
                    )

                    # Stack this Panel's vertices and faces with the array of all
                    # vertices and faces.
                    panel_vertices = np.vstack((panel_vertices, panel_vertices_to_add))
                    panel_faces = np.hstack((panel_faces, panel_face_to_add))

                    # Update the number of previous Panels.
                    panel_num += 1

                    # Convert the Panel vertices and faces to PolyData.
                    panel_surfaces = pv.PolyData(panel_vertices, panel_faces)

                    # Add the Panels to the plotter.
                    plotter.add_mesh(
                        panel_surfaces,
                        show_edges=True,
                        color="chartreuse",
                        smooth_shading=False,
                    )

        plotter.enable_parallel_projection()  # type: ignore[call-arg]

        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )

        return None

    def validate_first_airplane_constraints(self) -> None:
        """Validates that the first Airplane in a simulation has Cg_GP1_CgP1 set to
        zeros.

        This method should be called by SteadyProblem or UnsteadyProblem.

        :return: None
        """
        if not np.allclose(self._Cg_GP1_CgP1, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The first Airplane in a simulation must have Cg_GP1_CgP1 set to ("
                "0.0, 0.0, 0.0) by definition."
            )

    @staticmethod
    def process_wing_symmetry(wing: wing_mod.Wing) -> list[wing_mod.Wing]:
        """Processes a Wing to determine what type of symmetry it has.

        If necessary, it then modifies the Wing. If type 5 symmetry is detected, it also
        creates a second reflected Wing. Finally, it returns a list of Wings.

        :param wing: The Wing to process for symmetry analysis and potential
            modification.
        :return: The list of processed Wings. For types 1-4 symmetry it contains only
            the one modified Wing, but for type 5 symmetry it contains the modified Wing
            followed by the new reflected Wing. Before returning them, it also calls
            each Wing's generate_mesh method, preparing them for use simulation.
        """
        # Determine if the symmetry plane is coincident with the wing axes' xz plane. If
        # symmetryNormal_G or symmetryPoint_G_Cg is None, then there is no symmetry and
        # the symmetry plane doesn't exist. Otherwise, the symmetry plane is coincident
        # to the wing axes' xz plane if Ler_Gs_Cgs lies on the symmetry plane, and if
        # symmetryNormal_G is parallel with WnY_G. We don't need to check types, values,
        # or normalize because this is done in Wing's init method.
        coincident_symmetry_plane = True
        if wing.symmetryPoint_G_Cg is None or wing.symmetryNormal_G is None:
            coincident_symmetry_plane = False
        else:
            # If the symmetry plane exists, we first need to check if its normal vector
            # is parallel with the wing axes' y axis vector.

            # Actively transform geometry axes' second basis vector (in geometry axes)
            # to this Wing's axes' second basis vector (in geometry axes). We can skip
            # the translation step (step 2) as we are only transforming a direction
            # vector, not a position vector.
            GY_G = np.array([0.0, 1.0, 0.0], dtype=float)
            GsY_G = _transformations.apply_T_to_vectors(
                _transformations.generate_reflect_T(
                    plane_point_A_a=wing.symmetryPoint_G_Cg,
                    plane_normal_A=wing.symmetryNormal_G,
                    passive=False,
                ),
                GY_G,
                is_position=False,
            )
            WnY_G = _transformations.apply_T_to_vectors(
                _transformations.generate_rot_T(
                    wing.angles_Gs_to_Wn_ixyz,
                    passive=False,
                    intrinsic=True,
                    order="xyz",
                ),
                GsY_G,
                is_position=False,
            )

            # If symmetryNormal_G is parallel with WnY_G, their cross product will be
            # the zero vector.
            is_parallel = np.allclose(
                np.cross(wing.symmetryNormal_G, WnY_G),
                np.array([0.0, 0.0, 0.0], dtype=float),
            )

            if not is_parallel:
                coincident_symmetry_plane = False
            else:
                # If the symmetry plane's normal vector and the wing axes' y axis vector
                # are parallel, then the last check for a coincident symmetry plane is
                # to check if the Ler is on the symmetry plane.

                # To do this, we first find the symmetry plane's normal vector (in
                # geometry axes after accounting for symmetry) and the symmetry plane's
                # point (in geometry axes after accounting for symmetry, relative to the
                # CG after accounting for symmetry). As the symmetry plane is defined
                # using these quantities, they don't change after reflection.
                symmetryPoint_Gs_Cgs = wing.symmetryPoint_G_Cg
                symmetryNormal_Gs_Cgs = wing.symmetryNormal_G

                # The leading edge root point is on the symmetry plane if the distance
                # between it and the symmetry plane is zero.
                Ler_on_plane = np.allclose(
                    np.dot(
                        symmetryNormal_Gs_Cgs, (wing.Ler_Gs_Cgs - symmetryPoint_Gs_Cgs)
                    ),
                    0.0,
                )

                if not Ler_on_plane:
                    coincident_symmetry_plane = False

        # See the Wing class docstring for the interpretation of the different symmetry
        # types.
        if not wing.symmetric:
            if not wing.mirror_only:
                # Type 1 Symmetry:
                # symmetric = False, mirror_only = False
                symmetry_type = 1
            else:
                if coincident_symmetry_plane:
                    # Type 2 Symmetry:
                    # symmetric = False, mirror_only = True, coincident_symmetry_plane = True
                    symmetry_type = 2
                else:
                    # Type 3 Symmetry:
                    # symmetric = False, mirror_only = True, coincident_symmetry_plane = False
                    symmetry_type = 3
        else:
            if coincident_symmetry_plane:
                # Type 4 Symmetry:
                # symmetric = True, coincident_symmetry_plane = True
                symmetry_type = 4
            else:
                # Type 5 Symmetry:
                # symmetric = True, coincident_symmetry_plane = False
                symmetry_type = 5

        # Based on the determined symmetry type, validate the Wing's WingCrossSections'
        # control_surface_symmetry types. From the validation done during each
        # WingCrossSection's initialization method, we already know that
        # control_surface_symmetry type is None or a valid string.
        for wing_cross_section in wing.wing_cross_sections:
            control_surface_symmetry_type = (
                wing_cross_section.control_surface_symmetry_type
            )
            if symmetry_type in [1, 2, 3]:
                if control_surface_symmetry_type is not None:
                    raise ValueError(
                        f"control_surface_symmetry_type must be None for symmetry type "
                        f"{symmetry_type}"
                    )
            else:
                if wing_cross_section.control_surface_symmetry_type is None:
                    raise ValueError(
                        f"control_surface_symmetry_type must be specified for symmetry "
                        f"type {symmetry_type}"
                    )

        # For a type 4 Wing, the root WingCrossSection lies on the coincident symmetry
        # plane and is shared between the original and mirrored halves, which are meshed
        # into a single Panel grid joined at that centerline seam. An asymmetric control
        # surface there would deflect the two halves in opposite directions, so their
        # Panels would meet at the seam with non-coincident edges. The solvers assume
        # every internal edge within a Wing is shared by two coincident Panels (the
        # Kutta-Joukowski effective-strength subtraction and the wake shedding both
        # identify neighbors by grid index, not geometry), so such a torn seam would
        # silently corrupt the loads and the wake. A zero deflection is harmless because
        # it produces no geometric offset.
        if symmetry_type == 4:
            root_wing_cross_section = wing.wing_cross_sections[0]
            if (
                root_wing_cross_section.control_surface_symmetry_type == "asymmetric"
                and root_wing_cross_section.control_surface_deflection != 0.0
            ):
                raise ValueError(
                    'control_surface_symmetry_type cannot be "asymmetric" with a '
                    "nonzero control_surface_deflection on the root WingCrossSection "
                    "of a Wing with a coincident symmetry plane"
                )

        # Based on symmetry type, generate the mesh and return the Wing(s).
        if symmetry_type in [1, 2, 3, 4]:
            wing.generate_mesh(symmetry_type)
            return [wing]
        else:
            assert wing.symmetryNormal_G is not None
            assert wing.symmetryPoint_G_Cg is not None

            # Build the reflected Wing the same way the original was built so it carries
            # the same spanwise mesh provenance, which the convergence tools dispatch
            # on.
            if wing.spanwise_mesh == "edge_defined":
                # The stored edge curves are in wing axes, which for a type 3 Wing
                # already include the reflection, so they are reused verbatim.
                leadingEdgePoints_Wn_Ler = wing.leadingEdgePoints_Wn_Ler
                trailingEdgePoints_Wn_Ler = wing.trailingEdgePoints_Wn_Ler
                tip_trim_fraction = wing.tip_trim_fraction
                assert leadingEdgePoints_Wn_Ler is not None
                assert trailingEdgePoints_Wn_Ler is not None
                assert tip_trim_fraction is not None

                reflected_wing = wing_mod.Wing.from_edge_points(
                    leadingEdgePoints_Wn_Ler=np.copy(leadingEdgePoints_Wn_Ler),
                    trailingEdgePoints_Wn_Ler=np.copy(trailingEdgePoints_Wn_Ler),
                    num_wing_cross_sections=len(wing.wing_cross_sections),
                    airfoil=copy.deepcopy(wing.wing_cross_sections[0].airfoil),
                    name=f"Reflected {wing.name}",
                    Ler_Gs_Cgs=np.copy(wing.Ler_Gs_Cgs),
                    angles_Gs_to_Wn_ixyz=np.copy(wing.angles_Gs_to_Wn_ixyz),
                    symmetric=False,
                    mirror_only=True,
                    symmetryNormal_G=np.copy(wing.symmetryNormal_G),
                    symmetryPoint_G_Cg=np.copy(wing.symmetryPoint_G_Cg),
                    num_chordwise_panels=wing.num_chordwise_panels,
                    chordwise_spacing=wing.chordwise_spacing,
                    tip_trim_fraction=tip_trim_fraction,
                )
            else:
                reflected_wing_cross_sections = []
                for wing_cross_section in wing.wing_cross_sections:
                    airfoil = wing_cross_section.airfoil

                    reflected_airfoil = copy.deepcopy(airfoil)

                    if wing_cross_section.control_surface_symmetry_type == "asymmetric":
                        reflected_control_surface_deflection = (
                            -1 * wing_cross_section.control_surface_deflection
                        )
                    else:
                        reflected_control_surface_deflection = (
                            wing_cross_section.control_surface_deflection
                        )

                    reflected_wing_cross_sections.append(
                        wing_cross_section_mod.WingCrossSection(
                            airfoil=reflected_airfoil,
                            num_spanwise_panels=wing_cross_section.num_spanwise_panels,
                            chord=wing_cross_section.chord,
                            Lp_Wcsp_Lpp=np.copy(wing_cross_section.Lp_Wcsp_Lpp),
                            angles_Wcsp_to_Wcs_ixyz=np.copy(
                                wing_cross_section.angles_Wcsp_to_Wcs_ixyz
                            ),
                            control_surface_symmetry_type=None,
                            control_surface_hinge_point=wing_cross_section.control_surface_hinge_point,
                            control_surface_deflection=reflected_control_surface_deflection,
                            spanwise_spacing=wing_cross_section.spanwise_spacing,
                        )
                    )

                # Re-exploding an already exploded Wing reproduces it exactly.
                reflected_wing = wing_mod.Wing(
                    wing_cross_sections=reflected_wing_cross_sections,
                    name=f"Reflected {wing.name}",
                    Ler_Gs_Cgs=np.copy(wing.Ler_Gs_Cgs),
                    angles_Gs_to_Wn_ixyz=np.copy(wing.angles_Gs_to_Wn_ixyz),
                    symmetric=False,
                    mirror_only=True,
                    symmetryNormal_G=np.copy(wing.symmetryNormal_G),
                    symmetryPoint_G_Cg=np.copy(wing.symmetryPoint_G_Cg),
                    explode_into_strips=wing.spanwise_mesh == "exploded",
                    num_chordwise_panels=wing.num_chordwise_panels,
                    chordwise_spacing=wing.chordwise_spacing,
                )

            wing.symmetric = False
            wing.mirror_only = False
            wing.symmetryNormal_G = None
            wing.symmetryPoint_G_Cg = None

            # Reset control_surface_symmetry_type to None for Type 1 symmetry.
            for wing_cross_section in wing.wing_cross_sections:
                wing_cross_section.control_surface_symmetry_type = None

            wing.generate_mesh(symmetry_type=1)
            reflected_wing.generate_mesh(symmetry_type=3)
            return [wing, reflected_wing]


def _get_planform_reference_dimensions(
    first_wings: list[wing_mod.Wing],
) -> tuple[float, float, float]:
    """Calculates the default reference area, reference chord, and reference span from
    an Airplane's projected reference planform.

    The Airplane class docstring defines the reference planform, the three reference
    dimensions, and the conditions under which this function raises. A type 4 Wing's
    mirrored half is found by reflecting the original half, and for type 5 symmetry, the
    strip joining the two halves is called the bridge strip.

    :param first_wings: The list of meshed Wings that process_wing_symmetry returns for
        the first element of the wings passed to an Airplane. It holds one Wing, or two
        for type 5 symmetry.
    :return: A tuple of three floats, which are the reference area, the reference chord,
        and the reference span. Their units are square meters, meters, and meters.
    :raises ValueError: If the planform is too steeply inclined relative to the geometry
        axes' xy plane, or if its projection is ill-formed.
    """
    first_wing = first_wings[0]
    first_WnZ_G = first_wing.WnZ_G
    assert first_WnZ_G is not None

    # Build each half's strips (in geometry axes, relative to the CG). Each strip is a
    # (N, 2, 3) ndarray of N stations, ordered from root to tip, where each station
    # holds a leading point and then a trailing point. The list is indexed first by
    # half, and then by strip.
    listGridStripPoints_G_Cg = [_get_wing_strip_points(wing) for wing in first_wings]
    if first_wing.symmetry_type == 4:
        assert first_wing.symmetryPoint_G_Cg is not None
        assert first_wing.symmetryNormal_G is not None
        T_act_reflect = _transformations.generate_reflect_T(
            plane_point_A_a=first_wing.symmetryPoint_G_Cg,
            plane_normal_A=first_wing.symmetryNormal_G,
            passive=False,
        )
        listGridStripPoints_G_Cg.append(
            [
                _transformations.apply_T_to_vectors(
                    T_act_reflect, gridStripPoints_G_Cg.reshape(-1, 3), is_position=True
                ).reshape(-1, 2, 3)
                for gridStripPoints_G_Cg in listGridStripPoints_G_Cg[0]
            ]
        )

    # Collect every strip, with each half's strips first. For type 5 symmetry, the
    # bridge strip comes last. It joins the halves' root chords, which are the first
    # stations of each half's first strip.
    listAllGridStripPoints_G_Cg = [
        gridStripPoints_G_Cg
        for half in listGridStripPoints_G_Cg
        for gridStripPoints_G_Cg in half
    ]
    if len(first_wings) == 2:
        listAllGridStripPoints_G_Cg.append(
            np.array(
                [listGridStripPoints_G_Cg[0][0][0], listGridStripPoints_G_Cg[1][0][0]],
                dtype=float,
            )
        )
    half_ends = np.cumsum([len(half) for half in listGridStripPoints_G_Cg])

    # Find each strip's vector area (in geometry axes) as half the sum of the cross
    # products of its outline's consecutive points. The outline runs from root to tip
    # along the leading points, and then back from tip to root along the trailing
    # points. Taking the points relative to the outline's first point limits round off.
    # Also find each strip's range along the geometry axes' y axis.
    listVectorAreas_G = []
    for gridStripPoints_G_Cg in listAllGridStripPoints_G_Cg:
        outlinePoints_G_Cg = np.concatenate(
            (gridStripPoints_G_Cg[:, 0], gridStripPoints_G_Cg[::-1, 1])
        )
        outlinePoints_G_Cg = outlinePoints_G_Cg - outlinePoints_G_Cg[0]
        listVectorAreas_G.append(
            0.5
            * np.sum(
                np.cross(outlinePoints_G_Cg, np.roll(outlinePoints_G_Cg, -1, axis=0)),
                axis=0,
            )
        )
    stackVectorAreas_G = np.array(listVectorAreas_G, dtype=float)
    stripMinimumsY_G_Cg = np.array(
        [np.min(strip[:, :, 1]) for strip in listAllGridStripPoints_G_Cg], dtype=float
    )
    stripMaximumsY_G_Cg = np.array(
        [np.max(strip[:, :, 1]) for strip in listAllGridStripPoints_G_Cg], dtype=float
    )

    # Check that the halves aren't too steeply inclined relative to the geometry axes'
    # xy plane. The bridge strip is left out because it only fills the gap between the
    # halves by convention.
    stackHalfVectorAreas_G = stackVectorAreas_G[: half_ends[-1]]
    geometry_projected_area = float(np.sum(np.abs(stackHalfVectorAreas_G[:, 2])))
    wing_projected_area = float(np.sum(np.abs(stackHalfVectorAreas_G @ first_WnZ_G)))
    vector_area_magnitude = float(
        np.sum(np.linalg.norm(stackHalfVectorAreas_G, axis=1))
    )
    if (
        geometry_projected_area < wing_projected_area / np.sqrt(2.0)
        or geometry_projected_area
        <= _PLANFORM_RELATIVE_TOLERANCE * vector_area_magnitude
    ):
        raise ValueError(
            "The default reference dimensions come from the reference planform (defined "
            "in the Airplane class docstring) projected onto the geometry axes' xy "
            "plane, but that planform is too "
            "steeply inclined relative to that plane to serve as a reference. Its "
            "projected area is less than 1 / sqrt(2) times its area projected onto the "
            "xy plane of the wing axes of the first element of wings. Pass s_ref, "
            "c_ref, and b_ref explicitly."
        )

    # A strip is edge-on when its projected area is negligible compared to the magnitude
    # of its vector area. Edge-on strips play no further part, except in setting the
    # reference span.
    signed_areas = stackVectorAreas_G[:, 2]
    is_not_edge_on = np.abs(signed_areas) > _PLANFORM_RELATIVE_TOLERANCE * (
        np.linalg.norm(stackVectorAreas_G, axis=1)
    )

    # Check that the projected planform is well-formed. A strip whose projected area has
    # the opposite sign to its neighbors' isn't ill-formed by itself, since a Wing can
    # legitimately turn back toward its root (a winglet canted inward, for example), so
    # the checks below look for actual crossings and overlaps instead of sign changes.
    ill_formed_message = (
        "The default reference dimensions come from the reference planform (defined in "
        "the Airplane class docstring) projected onto the geometry axes' xy plane, but "
        "that projection is "
        "ill-formed because {}. Pass s_ref, c_ref, and b_ref explicitly. This geometry "
        "passed parameter validation, but it is likely a degenerate edge case, such as "
        "a Wing that folds or spirals back over itself. Other parts of the simulation "
        "may break or produce inaccurate results for it, so its definition is worth "
        "rechecking."
    )

    # For type 5 symmetry, check that the halves' ranges along the geometry axes' y axis
    # don't overlap.
    if len(first_wings) == 2:
        firstHalfMinimumY_G_Cg = np.min(stripMinimumsY_G_Cg[: half_ends[0]])
        firstHalfMaximumY_G_Cg = np.max(stripMaximumsY_G_Cg[: half_ends[0]])
        secondHalfMinimumY_G_Cg = np.min(
            stripMinimumsY_G_Cg[half_ends[0] : half_ends[1]]
        )
        secondHalfMaximumY_G_Cg = np.max(
            stripMaximumsY_G_Cg[half_ends[0] : half_ends[1]]
        )
        half_overlap = min(firstHalfMaximumY_G_Cg, secondHalfMaximumY_G_Cg) - max(
            firstHalfMinimumY_G_Cg, secondHalfMinimumY_G_Cg
        )
        half_extent = min(
            firstHalfMaximumY_G_Cg - firstHalfMinimumY_G_Cg,
            secondHalfMaximumY_G_Cg - secondHalfMinimumY_G_Cg,
        )
        if half_overlap > _PLANFORM_RELATIVE_TOLERANCE * half_extent:
            raise ValueError(
                ill_formed_message.format(
                    "its two halves overlap along the geometry axes' y axis"
                )
            )

    # Split each strip that isn't edge-on into triangles whose projected areas share its
    # sign. This fails for a strip that crosses over itself.
    listStackTrianglePointsXY_G_Cg = []
    triangle_strip_ids: list[int] = []
    for strip_id in np.flatnonzero(is_not_edge_on):
        stackStripTrianglePointsXY_G_Cg = _triangulate_strip(
            listAllGridStripPoints_G_Cg[strip_id][:, :, :2],
            float(signed_areas[strip_id]),
        )
        if stackStripTrianglePointsXY_G_Cg is None:
            raise ValueError(
                ill_formed_message.format("part of it crosses over itself")
            )
        listStackTrianglePointsXY_G_Cg.append(stackStripTrianglePointsXY_G_Cg)
        triangle_strip_ids.extend(
            [int(strip_id)] * len(stackStripTrianglePointsXY_G_Cg)
        )
    stackTrianglePointsXY_G_Cg = np.concatenate(listStackTrianglePointsXY_G_Cg)

    # Check that no two strips' interiors overlap.
    if _triangles_overlap(
        stackTrianglePointsXY_G_Cg, np.array(triangle_strip_ids, dtype=int)
    ):
        raise ValueError(ill_formed_message.format("two of its parts overlap"))

    # Find the reference dimensions.
    b_ref = float(np.max(stripMaximumsY_G_Cg) - np.min(stripMinimumsY_G_Cg))
    s_ref = float(np.sum(np.abs(signed_areas[is_not_edge_on])))
    c_ref = _get_chord_squared_integral(stackTrianglePointsXY_G_Cg) / s_ref

    return s_ref, c_ref, b_ref


def _get_wing_strip_points(wing: wing_mod.Wing) -> list[np.ndarray]:
    """Returns the strips that make up one meshed Wing's planform (in geometry axes,
    relative to the CG).

    For a Wing built from WingCrossSections, each strip joins two neighboring
    WingCrossSections' leading points and undeflected trailing points. A
    WingCrossSection's undeflected trailing point is (chord, 0.0, 0.0) in its own axes,
    relative to its leading point, so control surface deflections don't affect it.

    An edge_defined Wing is planar, so its whole half is one strip, bounded by the
    stored, untrimmed edge curves. Each curve is the monotone PCHIP interpolation of its
    x component as a function of its y component, which is the same curve used when
    resampling the edges for meshing. Both curves are sampled at every stored point of
    either curve, plus a fixed number of evenly spaced points between each neighboring
    pair.

    :param wing: The meshed Wing whose strips to return.
    :return: A list of (N, 2, 3) ndarrays of floats, one per strip, ordered from root to
        tip. Each holds N stations, ordered from root to tip, where each station holds a
        leading point and then a trailing point (in geometry axes, relative to the CG).
        The units are in meters.
    """
    if wing.spanwise_mesh == "edge_defined":
        leadingEdgePoints_Wn_Ler = wing.leadingEdgePoints_Wn_Ler
        trailingEdgePoints_Wn_Ler = wing.trailingEdgePoints_Wn_Ler
        T_pas_Wn_Ler_to_G_Cg = wing.T_pas_Wn_Ler_to_G_Cg
        assert leadingEdgePoints_Wn_Ler is not None
        assert trailingEdgePoints_Wn_Ler is not None
        assert T_pas_Wn_Ler_to_G_Cg is not None

        # Sample both curves at the same y values, so that the strip's stations each
        # hold a leading point and a trailing point at the same y component.
        listSampleYs_Wn_Ler = [
            np.linspace(start_y, end_y, _NUM_EDGE_CURVE_SAMPLES_BETWEEN_POINTS + 2)
            for edgePoints_Wn_Ler in (
                leadingEdgePoints_Wn_Ler,
                trailingEdgePoints_Wn_Ler,
            )
            for start_y, end_y in zip(
                edgePoints_Wn_Ler[:-1, 1], edgePoints_Wn_Ler[1:, 1]
            )
        ]
        sampleYs_Wn_Ler = np.unique(np.concatenate(listSampleYs_Wn_Ler))

        gridStripPoints_Wn_Ler = np.zeros((len(sampleYs_Wn_Ler), 2, 3), dtype=float)
        for edge_id, edgePoints_Wn_Ler in enumerate(
            (leadingEdgePoints_Wn_Ler, trailingEdgePoints_Wn_Ler)
        ):
            gridStripPoints_Wn_Ler[:, edge_id, 0] = sp_interp.PchipInterpolator(
                edgePoints_Wn_Ler[:, 1], edgePoints_Wn_Ler[:, 0]
            )(sampleYs_Wn_Ler)
            gridStripPoints_Wn_Ler[:, edge_id, 1] = sampleYs_Wn_Ler

        gridStripPoints_G_Cg = _transformations.apply_T_to_vectors(
            T_pas_Wn_Ler_to_G_Cg,
            gridStripPoints_Wn_Ler.reshape(-1, 3),
            is_position=True,
        ).reshape(-1, 2, 3)
        return [gridStripPoints_G_Cg]

    gridSectionPoints_G_Cg = np.array(
        [
            _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_G_Cg,
                np.array(
                    [[0.0, 0.0, 0.0], [wing_cross_section.chord, 0.0, 0.0]],
                    dtype=float,
                ),
                is_position=True,
            )
            for T_pas_Wcs_Lp_to_G_Cg, wing_cross_section in zip(
                wing.children_T_pas_Wcs_Lp_to_G_Cg, wing.wing_cross_sections
            )
        ],
        dtype=float,
    )
    return [
        gridSectionPoints_G_Cg[section_id : section_id + 2]
        for section_id in range(len(gridSectionPoints_G_Cg) - 1)
    ]


def _triangulate_strip(
    gridStripPointsXY_G_Cg: np.ndarray, signed_area: float
) -> np.ndarray | None:
    """Splits a projected strip into triangles whose projected areas all share the
    strip's sign.

    Each quadrilateral between neighboring stations is split along whichever of its
    diagonals gives two triangles that share the strip's sign. Triangles with negligible
    area are dropped.

    :param gridStripPointsXY_G_Cg: A (N, 2, 2) ndarray of floats holding the strip's N
        stations, where each station holds a leading point and then a trailing point (in
        geometry axes, projected onto its xy plane, relative to the CG). The units are
        in meters.
    :param signed_area: The strip's signed projected area. It must not be zero. The
        units are square meters.
    :return: A (K, 3, 2) ndarray of floats holding the K triangles' points, each in
        counterclockwise order when viewed from above the geometry axes' xy plane (in
        geometry axes, projected onto its xy plane, relative to the CG), or None if a
        quadrilateral can't be split into triangles that share the strip's sign, which
        means the strip crosses over itself. The units are in meters.
    """
    # Each quadrilateral's corners run in the same direction as the strip's outline,
    # from its inner leading point to its outer leading point, then to its outer
    # trailing point, and then to its inner trailing point.
    gridQuadrilateralPointsXY_G_Cg = np.stack(
        (
            gridStripPointsXY_G_Cg[:-1, 0],
            gridStripPointsXY_G_Cg[1:, 0],
            gridStripPointsXY_G_Cg[1:, 1],
            gridStripPointsXY_G_Cg[:-1, 1],
        ),
        axis=1,
    )

    # Each of a quadrilateral's two diagonals splits it into two triangles. The first
    # diagonal joins its first and third corners, and the second joins its second and
    # fourth corners. The candidate triangles are indexed by quadrilateral, then by
    # diagonal, and then by triangle.
    triangle_corner_ids = np.array(
        [[[0, 1, 2], [0, 2, 3]], [[0, 1, 3], [1, 2, 3]]], dtype=int
    )
    gridCandidateTrianglePointsXY_G_Cg = gridQuadrilateralPointsXY_G_Cg[
        :, triangle_corner_ids
    ]
    candidate_areas = _get_signed_triangle_areas(
        gridCandidateTrianglePointsXY_G_Cg.reshape(-1, 3, 2)
    ).reshape(-1, 2, 2)

    # A diagonal is valid for a quadrilateral if neither of its triangles has a non
    # negligible area with the opposite sign to the strip's.
    tolerance = _PLANFORM_RELATIVE_TOLERANCE * abs(signed_area)
    sign = np.sign(signed_area)
    is_valid = np.all(sign * candidate_areas >= -tolerance, axis=2)
    if not np.all(np.any(is_valid, axis=1)):
        return None

    quadrilateral_ids = np.arange(len(gridQuadrilateralPointsXY_G_Cg))
    diagonal_ids = np.where(is_valid[:, 0], 0, 1)
    stackTrianglePointsXY_G_Cg: np.ndarray = gridCandidateTrianglePointsXY_G_Cg[
        quadrilateral_ids, diagonal_ids
    ].reshape(-1, 3, 2)
    areas = candidate_areas[quadrilateral_ids, diagonal_ids].ravel()
    stackTrianglePointsXY_G_Cg = stackTrianglePointsXY_G_Cg[np.abs(areas) > tolerance]

    # Reverse the triangles' points if needed to put them in counterclockwise order.
    if sign < 0.0:
        stackTrianglePointsXY_G_Cg = stackTrianglePointsXY_G_Cg[:, ::-1]
    return stackTrianglePointsXY_G_Cg


def _get_signed_triangle_areas(stackTrianglePointsXY_G_Cg: np.ndarray) -> np.ndarray:
    """Returns the signed areas of projected triangles.

    :param stackTrianglePointsXY_G_Cg: A (K, 3, 2) ndarray of floats holding the K
        triangles' points (in geometry axes, projected onto its xy plane, relative to
        the CG). The units are in meters.
    :return: A (K,) ndarray of floats holding the triangles' signed areas, which are
        positive for triangles whose points are in counterclockwise order when viewed
        from above the geometry axes' xy plane. The units are square meters.
    """
    firstLegsXY_G = stackTrianglePointsXY_G_Cg[:, 1] - stackTrianglePointsXY_G_Cg[:, 0]
    secondLegsXY_G = stackTrianglePointsXY_G_Cg[:, 2] - stackTrianglePointsXY_G_Cg[:, 0]
    signed_areas: np.ndarray = 0.5 * (
        firstLegsXY_G[:, 0] * secondLegsXY_G[:, 1]
        - firstLegsXY_G[:, 1] * secondLegsXY_G[:, 0]
    )
    return signed_areas


def _triangles_overlap(
    stackTrianglePointsXY_G_Cg: np.ndarray, strip_ids: np.ndarray
) -> bool:
    """Checks whether the interiors of any two triangles from different strips overlap.

    Each pair of triangles from different strips whose bounding boxes overlap is clipped
    against each other, and the pair overlaps if the intersection's area exceeds a small
    fraction of the smaller triangle's area. The threshold can't be zero, because
    clipping two triangles that only touch along an edge can leave a sliver whose area
    is on the order of round off.

    :param stackTrianglePointsXY_G_Cg: A (K, 3, 2) ndarray of floats holding the K
        triangles' points, each in counterclockwise order when viewed from above the
        geometry axes' xy plane (in geometry axes, projected onto its xy plane, relative
        to the CG). The units are in meters.
    :param strip_ids: A (K,) ndarray of ints holding the ID of the strip each triangle
        came from.
    :return: True if the interiors of two triangles from different strips overlap, and
        False otherwise.
    """
    if len(np.unique(strip_ids)) < 2:
        return False

    # Sort the triangles by the minimum y components of their bounding boxes. Each
    # triangle then only needs to be compared with the later triangles whose bounding
    # boxes start below the end of its own.
    stackMinimumsXY_G_Cg = np.min(stackTrianglePointsXY_G_Cg, axis=1)
    stackMaximumsXY_G_Cg = np.max(stackTrianglePointsXY_G_Cg, axis=1)
    order = np.argsort(stackMinimumsXY_G_Cg[:, 1], kind="stable")
    stackTrianglePointsXY_G_Cg = stackTrianglePointsXY_G_Cg[order]
    stackMinimumsXY_G_Cg = stackMinimumsXY_G_Cg[order]
    stackMaximumsXY_G_Cg = stackMaximumsXY_G_Cg[order]
    strip_ids = strip_ids[order]
    areas = _get_signed_triangle_areas(stackTrianglePointsXY_G_Cg)
    end_ids = np.searchsorted(
        stackMinimumsXY_G_Cg[:, 1], stackMaximumsXY_G_Cg[:, 1], side="left"
    )

    for triangle_id, end_id in enumerate(end_ids):
        other_ids = np.arange(triangle_id + 1, end_id)
        candidate_ids = other_ids[
            (strip_ids[other_ids] != strip_ids[triangle_id])
            & (
                stackMinimumsXY_G_Cg[other_ids, 0]
                < stackMaximumsXY_G_Cg[triangle_id, 0]
            )
            & (
                stackMaximumsXY_G_Cg[other_ids, 0]
                > stackMinimumsXY_G_Cg[triangle_id, 0]
            )
        ]
        for other_id in candidate_ids:
            intersection_area = _get_triangle_intersection_area(
                stackTrianglePointsXY_G_Cg[triangle_id],
                stackTrianglePointsXY_G_Cg[other_id],
            )
            if intersection_area > _PLANFORM_RELATIVE_TOLERANCE * min(
                areas[triangle_id], areas[other_id]
            ):
                return True
    return False


def _get_triangle_intersection_area(
    trianglePointsXY_G_Cg: np.ndarray, clipTrianglePointsXY_G_Cg: np.ndarray
) -> float:
    """Returns the area of the intersection of two projected triangles whose points are
    in counterclockwise order.

    The first triangle is clipped by each edge of the second triangle in turn
    (Sutherland-Hodgman clipping), and the area of the remaining polygon is returned.

    **Citation:**

    Algorithm adapted from: "Reentrant Polygon Clipping"

    Authors: Ivan E. Sutherland and Gary W. Hodgman

    Date of retrieval: 09/23/2026

    :param trianglePointsXY_G_Cg: A (3, 2) ndarray of floats holding the points of the
        triangle to clip, in counterclockwise order (in geometry axes, projected onto
        its xy plane, relative to the CG). The units are in meters.
    :param clipTrianglePointsXY_G_Cg: A (3, 2) ndarray of floats holding the points of
        the triangle to clip by, in counterclockwise order (in geometry axes, projected
        onto its xy plane, relative to the CG). The units are in meters.
    :return: The intersection's area. The units are square meters.
    """
    polygonPointsXY_G_Cg = trianglePointsXY_G_Cg
    for edge_id in range(3):
        edgeStartXY_G_Cg = clipTrianglePointsXY_G_Cg[edge_id]
        edgeXY_G = clipTrianglePointsXY_G_Cg[(edge_id + 1) % 3] - edgeStartXY_G_Cg

        # A point is inside the clipping edge when it's on the edge's left side, where a
        # counterclockwise triangle's interior is.
        distances = edgeXY_G[0] * (
            polygonPointsXY_G_Cg[:, 1] - edgeStartXY_G_Cg[1]
        ) - edgeXY_G[1] * (polygonPointsXY_G_Cg[:, 0] - edgeStartXY_G_Cg[0])
        is_inside = distances >= 0.0

        listClippedPointsXY_G_Cg = []
        num_points = len(polygonPointsXY_G_Cg)
        for point_id in range(num_points):
            next_point_id = (point_id + 1) % num_points
            if is_inside[point_id]:
                listClippedPointsXY_G_Cg.append(polygonPointsXY_G_Cg[point_id])
            if is_inside[point_id] != is_inside[next_point_id]:
                fraction = distances[point_id] / (
                    distances[point_id] - distances[next_point_id]
                )
                listClippedPointsXY_G_Cg.append(
                    polygonPointsXY_G_Cg[point_id]
                    + fraction
                    * (
                        polygonPointsXY_G_Cg[next_point_id]
                        - polygonPointsXY_G_Cg[point_id]
                    )
                )
        if len(listClippedPointsXY_G_Cg) < 3:
            return 0.0
        polygonPointsXY_G_Cg = np.array(listClippedPointsXY_G_Cg, dtype=float)

    nextPolygonPointsXY_G_Cg = np.roll(polygonPointsXY_G_Cg, -1, axis=0)
    return max(
        0.5
        * float(
            np.sum(
                polygonPointsXY_G_Cg[:, 0] * nextPolygonPointsXY_G_Cg[:, 1]
                - polygonPointsXY_G_Cg[:, 1] * nextPolygonPointsXY_G_Cg[:, 0]
            )
        ),
        0.0,
    )


def _get_chord_squared_integral(stackTrianglePointsXY_G_Cg: np.ndarray) -> float:
    """Returns the integral of the square of the projected chord along the geometry
    axes' y axis, for a set of non overlapping projected triangles.

    The projected chord at each position along the geometry axes' y axis is the total
    length, along the geometry axes' x axis, of the parts of that slice inside the
    triangles. Each triangle's slice length varies linearly between the y components of
    its points, so the chord varies linearly between neighboring y components of all the
    triangles' points. The integral is therefore evaluated exactly, one piece at a time.

    :param stackTrianglePointsXY_G_Cg: A (K, 3, 2) ndarray of floats holding the K
        triangles' points (in geometry axes, projected onto its xy plane, relative to
        the CG). Every triangle must have a non zero area. The units are in meters.
    :return: The integral. The units are cubic meters.
    """
    # Each triangle's slice length is zero at its lowest and highest points and peaks at
    # the y component of its middle point. A triangle's area is half its height times
    # that peak length.
    sortedPointsY_G_Cg = np.sort(stackTrianglePointsXY_G_Cg[:, :, 1], axis=1)
    lowPointsY_G_Cg = sortedPointsY_G_Cg[:, 0]
    middlePointsY_G_Cg = sortedPointsY_G_Cg[:, 1]
    highPointsY_G_Cg = sortedPointsY_G_Cg[:, 2]
    middle_lengths = (
        2.0
        * np.abs(_get_signed_triangle_areas(stackTrianglePointsXY_G_Cg))
        / (highPointsY_G_Cg - lowPointsY_G_Cg)
    )

    # Split each triangle's slice length into two linear pieces, and keep those with a
    # non zero extent along the geometry axes' y axis.
    zeros = np.zeros_like(middle_lengths)
    start_ys = np.concatenate((lowPointsY_G_Cg, middlePointsY_G_Cg))
    end_ys = np.concatenate((middlePointsY_G_Cg, highPointsY_G_Cg))
    start_lengths = np.concatenate((zeros, middle_lengths))
    end_lengths = np.concatenate((middle_lengths, zeros))
    has_extent = end_ys > start_ys
    start_ys = start_ys[has_extent]
    end_ys = end_ys[has_extent]
    start_lengths = start_lengths[has_extent]
    end_lengths = end_lengths[has_extent]

    # Add each piece's contribution to the chord at both ends of every interval, between
    # neighboring breakpoints, that it covers.
    breakpoint_ys = np.unique(np.concatenate((start_ys, end_ys)))
    first_interval_ids = np.searchsorted(breakpoint_ys, start_ys)
    num_intervals = np.searchsorted(breakpoint_ys, end_ys) - first_interval_ids
    piece_ids = np.repeat(np.arange(len(start_ys)), num_intervals)
    interval_ids = np.repeat(first_interval_ids, num_intervals) + (
        np.arange(len(piece_ids))
        - np.repeat(np.cumsum(num_intervals) - num_intervals, num_intervals)
    )
    piece_heights = end_ys[piece_ids] - start_ys[piece_ids]
    piece_length_changes = end_lengths[piece_ids] - start_lengths[piece_ids]
    start_chords = np.zeros(len(breakpoint_ys) - 1, dtype=float)
    end_chords = np.zeros(len(breakpoint_ys) - 1, dtype=float)
    for end_offset, chords in ((0, start_chords), (1, end_chords)):
        np.add.at(
            chords,
            interval_ids,
            start_lengths[piece_ids]
            + piece_length_changes
            * (breakpoint_ys[interval_ids + end_offset] - start_ys[piece_ids])
            / piece_heights,
        )

    # The chord is linear on each interval, so the integral of its square over an
    # interval is the interval's width times the mean of the squares and the product of
    # the chords at its ends.
    return float(
        np.sum(
            np.diff(breakpoint_ys)
            * (start_chords**2 + start_chords * end_chords + end_chords**2)
            / 3.0
        )
    )
