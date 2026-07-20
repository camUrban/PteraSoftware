"""Contains the Wing class.

**Contains the following classes:**

Wing: A class used to contain wings of an Airplane.

**Contains the following functions:**

None
"""

from __future__ import annotations

import copy
from collections.abc import Sequence

import numpy as np
import pyvista as pv
import scipy.interpolate as sp_interp

from .. import _panel, _parameter_validation, _transformations
from . import _meshing
from . import airfoil as airfoil_mod
from . import wing_cross_section as wing_cross_section_mod


class Wing:
    """A class used to contain the wings of an Airplane.

    **Contains the following methods:**

    from_edge_points: Builds a planar, untwisted Wing directly from leading edge and
    trailing edge curves, approximating an arbitrary planform with a sequence of
    WingCrossSections joined by single-panel strips.

    __deepcopy__: Creates a deep copy of this Wing, preserving mesh geometry but
    resetting wake state.

    T_pas_G_Cg_to_Wn_Ler: The passive transformation matrix which maps in homogeneous
    coordinates from geometry axes relative to the CG to wing axes relative to the
    leading edge root point. It is None if the Wing's symmetry type hasn't been defined
    yet.

    T_pas_Wn_Ler_to_G_Cg: The passive transformation matrix which maps in homogeneous
    coordinates from wing axes relative to the leading edge root point to geometry axes
    relative to the CG point. It is None if the Wing's symmetry type hasn't been defined
    yet.

    WnX_G: The wing axes' first basis vector (in geometry axes).

    WnY_G: The wing axes' second basis vector (in geometry axes).

    WnZ_G: The wing axes' third basis vector (in geometry axes).

    children_T_pas_Wn_Ler_to_Wcs_Lp: A list of passive transformation matrices which map
    in homogeneous coordinates from wing axes, relative to the leading edge root point,
    to each of this Wing's WingCrossSection's axes, relative to their respective leading
    points.

    children_T_pas_Wcs_Lp_to_Wn_Ler: A list of passive transformation matrices which map
    in homogeneous coordinates from each of this Wing's WingCrossSection's axes,
    relative to their respective leading points, to wing axes, relative to the leading
    edge root point.

    children_T_pas_G_Cg_to_Wcs_Lp: A list of passive transformation matrices which map
    in homogeneous coordinates from geometry axes, relative to the CG, to each of this
    Wing's WingCrossSection's axes, relative to their respective leading points.

    children_T_pas_Wcs_Lp_to_G_Cg: A list of passive transformation matrices which map
    in homogeneous coordinates from each of this Wing's WingCrossSection's axes,
    relative to their respective leading points, to geometry axes, relative to the CG.

    symmetry_type: The symmetry type of this Wing.

    num_spanwise_panels: The number of spanwise Panels on this Wing.

    num_panels: The total number of Panels on this Wing.

    panels: The 2D array of Panels on this Wing.

    projected_area: The area of the Wing projected onto the plane defined by the wing
    axes' xy plane.

    wetted_area: The Wing's wetted area.

    average_panel_aspect_ratio: The average aspect ratio of the Wing's Panels.

    span: The Wing's span.

    standard_mean_chord: The Wing's standard mean chord.

    mean_aerodynamic_chord: The Wing's mean aerodynamic chord.

    spanwise_mesh: How this Wing's spanwise mesh was defined.

    leadingEdgePoints_Wn_Ler: The original leading edge curve this Wing was built from,
    or None.

    trailingEdgePoints_Wn_Ler: The original trailing edge curve this Wing was built
    from, or None.

    tip_trim_fraction: The fraction of the span dropped off the tip when this Wing was
    built, or None.

    generate_mesh: Generates this Wing's mesh, which finishes the process of preparing
    the Wing to be used in a simulation. It is called by the Wing's parent Airplane,
    after it's determined its symmetry type.

    get_plottable_data: Returns plottable data for this Wing's Airfoils' outlines and
    mean camber lines.

    **Notes:**

    Immutable attributes (wing_cross_sections, name, Ler_Gs_Cgs, angles_Gs_to_Wn_ixyz,
    num_chordwise_panels, chordwise_spacing, spanwise_mesh, leadingEdgePoints_Wn_Ler,
    and trailingEdgePoints_Wn_Ler) are set during initialization and cannot be modified
    afterward. The numpy arrays Ler_Gs_Cgs and angles_Gs_to_Wn_ixyz are made read only
    to prevent in place mutation. The wing_cross_sections attribute is stored as a tuple
    to prevent external mutation. The spanwise_mesh attribute is not a constructor
    parameter; the standard constructor derives it from explode_into_strips, and
    from_edge_points sets it to "edge_defined". The leadingEdgePoints_Wn_Ler,
    trailingEdgePoints_Wn_Ler, and tip_trim_fraction attributes are None unless the Wing
    was built by from_edge_points, in which case they hold its original edge curves (as
    read only arrays) and the tip trim fraction applied while resampling them.

    Derived transformation matrices and basis vectors (T_pas_G_Cg_to_Wn_Ler,
    T_pas_Wn_Ler_to_G_Cg, WnX_G, WnY_G, WnZ_G, and the children_T_pas_* properties) are
    lazily evaluated and cached. Derived geometric properties (projected_area,
    wetted_area, average_panel_aspect_ratio, span, standard_mean_chord, and
    mean_aerodynamic_chord) are also lazily evaluated and cached.

    The symmetry_type, num_spanwise_panels, num_panels, and panels attributes are set
    once by generate_mesh and cannot be modified after being set.

    The symmetric, mirror_only, symmetryNormal_G, and symmetryPoint_G_Cg attributes
    remain mutable as they may be modified by Airplane.process_wing_symmetry() for type
    5 symmetry handling. The gridWrvp_GP1_CgP1 attribute is mutable as it is modified
    during simulation.

    Every Wing has its own axis system, known as wing axes. The user sets the
    relationship between these axes and geometry axes with the Ler_Gs_Cgs and
    angles_Gs_to_Wn_ixyz parameters. However, the steps for transforming a vector from
    geometry axes to wing axes, and the interpretation of the wing axes orientation and
    position relative to an Airplane's geometry axes, also depend on the parameters
    symmetric, mirror_only, symmetryNormal_G, and symmetryPoint_G_Cg. In all cases, the
    order of transformations from geometry axes to wing axes is reflection (if
    applicable), translation, and then rotation.

    There are five symmetry types. Type 1: symmetric=False, mirror_only=False, and the
    symmetry plane must be undefined. Type 2: symmetric=False, mirror_only=True, and the
    symmetry plane is coincident with the wing axes' xz plane. Type 3: symmetric=False,
    mirror_only=True, and the symmetry plane is not coincident with the wing axes' xz
    plane. Type 4: symmetric=True, mirror_only=False, and the symmetry plane is
    coincident with the wing axes' xz plane. Type 5: symmetric=True, mirror_only=False,
    and the symmetry plane is not coincident with the wing axes' xz plane.

    **Citation:**

    Adapted from: geometry.Wing in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/24/2020
    """

    __slots__ = (
        # Immutable
        "_wing_cross_sections",
        "_name",
        "_Ler_Gs_Cgs",
        "_angles_Gs_to_Wn_ixyz",
        "_num_chordwise_panels",
        "_chordwise_spacing",
        "_spanwise_mesh",
        "_leadingEdgePoints_Wn_Ler",
        "_trailingEdgePoints_Wn_Ler",
        "_tip_trim_fraction",
        # Mutable (type 5 symmetry)
        "symmetric",
        "mirror_only",
        "symmetryNormal_G",
        "symmetryPoint_G_Cg",
        # Set once
        "_symmetry_type",
        "_num_spanwise_panels",
        "_num_panels",
        "_panels",
        # Mutable (wake)
        "gridWrvp_GP1_CgP1",
        # Caches from immutable
        "_T_pas_G_Cg_to_Wn_Ler",
        "_T_pas_Wn_Ler_to_G_Cg",
        "_WnX_G",
        "_WnY_G",
        "_WnZ_G",
        "_children_T_pas_Wn_Ler_to_Wcs_Lp",
        "_children_T_pas_Wcs_Lp_to_Wn_Ler",
        "_children_T_pas_G_Cg_to_Wcs_Lp",
        "_children_T_pas_Wcs_Lp_to_G_Cg",
        # Caches from set once
        "_projected_area",
        "_wetted_area",
        "_average_panel_aspect_ratio",
        "_span",
        "_standard_mean_chord",
        "_mean_aerodynamic_chord",
    )

    def __init__(
        self,
        wing_cross_sections: list[wing_cross_section_mod.WingCrossSection],
        name: str = "Untitled Wing",
        Ler_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        angles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        symmetric: bool | np.bool_ = False,
        mirror_only: bool | np.bool_ = False,
        symmetryNormal_G: None | np.ndarray | Sequence[float | int] = None,
        symmetryPoint_G_Cg: None | np.ndarray | Sequence[float | int] = None,
        explode_into_strips: bool | np.bool_ = False,
        num_chordwise_panels: int = 8,
        chordwise_spacing: str = "cosine",
    ) -> None:
        """The initialization method.

        :param wing_cross_sections: A list of WingCrossSections representing the wing's
            cross sections in order from root to tip. It must contain at least two
            WingCrossSections.
        :param name: A sensible name for the Wing. The default is "Untitled Wing".
        :param Ler_Gs_Cgs: An array-like object of 3 numbers (int or float) representing
            the position of the origin of this Wing's axes (in geometry axes after
            accounting for symmetry, relative to the CG after accounting for symmetry).
            Can be a tuple, list, or ndarray. Values are converted to floats internally.
            The units are meters. The default is (0.0, 0.0, 0.0).
        :param angles_Gs_to_Wn_ixyz: An array-like object of 3 numbers (int or float)
            representing the angle vector that defines the orientation of this Wing's
            axes relative to the geometry axes (after accounting for symmetry). Can be a
            tuple, list, or ndarray. Values are converted to floats internally. All
            angles must be in the range [-90, 90] degrees. Rotations are intrinsic, and
            proceed in the xy'z" order. The units are degrees. The default is (0.0, 0.0,
            0.0).
        :param symmetric: Set this to True if the Wing's geometry should be mirrored
            across the symmetry plane while retaining the non mirrored side. If
            mirror_only is True, symmetric must be False. If symmetric is True, then
            neither symmetryNormal_G nor symmetryPoint_G_Cg can be None. If the symmetry
            plane is coincident with this Wing's axes' xz plane, the mirrored and non
            mirrored geometry will be meshed as a single wing. If not, this Wing's
            Airplane will automatically create another Wing with the mirrored geometry,
            modify both Wings' parameters, and add the reflected Wing to its list of
            Wings immediately following this one. For more details on that process, and
            how this parameter interacts with symmetryNormal_G, symmetryPoint_G_Cg, and
            mirror_only, see the class docstring. Can be a bool or a numpy bool and will
            be converted internally to a bool. The default is False.
        :param mirror_only: Set this to True if the Wing's geometry should be reflected
            about the symmetry plane without retaining the non reflected geometry. If
            symmetric is True, mirror_only must be False. If mirror_only is True, then
            neither symmetryNormal_G nor symmetryPoint_G_Cg can be None. For more
            details on how this parameter interacts with symmetryNormal_G,
            symmetryPoint_G_Cg, and symmetric, see the class docstring. Can be a bool or
            a numpy bool and will be converted internally to a bool. The default is
            False.
        :param symmetryNormal_G: None, or an array-like of 3 numbers (int or float)
            representing the unit normal vector (in geometry axes) that, together with
            symmetryPoint_G_Cg, defines the plane used for symmetry or mirroring. Can be
            None, or a tuple, list, or ndarray. If not None, values are converted to
            floats and normalized internally. Note that reversing the normal direction
            (using the antiparallel vector) defines the same plane and produces the same
            result. This value must be None if both symmetric and mirror_only are False,
            and cannot be None if either are True. For more details on how this
            parameter interacts with symmetryPoint_G_Cg, symmetric, and mirror_only, see
            the class docstring. The default is None.
        :param symmetryPoint_G_Cg: None or an array-like object of 3 numbers (int or
            float) representing a point (in geometry axes, relative to the CG) that,
            along with symmetryNormal_G, defines the location of the plane about which
            symmetry or mirroring is applied. Can be None, or a list, tuple, or ndarray.
            If not None, values are converted to floats internally. This value must be
            None if both symmetric and mirror_only are False, and cannot be None if
            either are True. For more details on how this parameter interacts with
            symmetryNormal_G, symmetric, and mirror_only, see the class docstring. The
            units are meters. The default is None.
        :param explode_into_strips: Set this to True to replace wing_cross_sections
            during initialization with a new list in which every panel is broken into
            single spanwise strips for deformation. When True, every non tip
            WingCrossSection in wing_cross_sections must have
            spanwise_spacing="uniform"; the explosion assumes uniformly spaced
            intermediates and rejects other spacings rather than silently overriding
            them. This parameter is consumed during initialization and is not stored as
            an attribute.
        :param num_chordwise_panels: The number of chordwise panels to be used on this
            Wing, which must be set to a positive integer. The default is 8.
        :param chordwise_spacing: The type of spacing between the Wing's chordwise
            panels. Can be "cosine" or "uniform". Using cosine spacing is highly
            recommended for steady simulations and uniform spacing is highly recommended
            for unsteady simulations. The default is "cosine".
        :return: None
        """
        # Validate wing_cross_sections.
        wing_cross_sections = _parameter_validation.non_empty_list_return_list(
            wing_cross_sections, "wing_cross_sections"
        )
        explode_into_strips = _parameter_validation.boolLike_return_bool(
            explode_into_strips, "explode_into_strips"
        )
        if explode_into_strips:
            wing_cross_sections = self._explode_wing(wing_cross_sections)
        num_wing_cross_sections = len(wing_cross_sections)
        if num_wing_cross_sections < 2:
            raise ValueError("wing_cross_sections must contain at least two elements.")
        for wing_cross_section_id, wing_cross_section in enumerate(wing_cross_sections):
            if not isinstance(
                wing_cross_section, wing_cross_section_mod.WingCrossSection
            ):
                raise TypeError(
                    "Every element in wing_cross_sections must be a WingCrossSection."
                )
            if wing_cross_section_id == 0:
                # Validate root WingCrossSection constraints.
                wing_cross_section.validate_root_constraints()
            elif wing_cross_section_id == num_wing_cross_sections - 1:
                # Validate tip WingCrossSection constraints.
                wing_cross_section.validate_tip_constraints()
            else:
                wing_cross_section.validate_mid_constraints()
            # Set the validated flag for this WingCrossSection.
            wing_cross_section.validated = True
        # Store as tuple to prevent external mutation.
        self._wing_cross_sections: tuple[
            wing_cross_section_mod.WingCrossSection, ...
        ] = tuple(wing_cross_sections)

        # Record how this Wing's spanwise mesh was defined (immutable). A Wing built
        # with explode_into_strips has every non tip WingCrossSection as a single
        # spanwise panel strip, so it cannot be refined by the standard convergence
        # tools, which sweep the number of spanwise Panels while holding the
        # WingCrossSections fixed. The standard convergence tools dispatch on this
        # marker to reject such a Wing. The marker is set by provenance rather than
        # detected from the geometry, since an exploded Wing is structurally identical
        # to one a user defines with single strip WingCrossSections directly, yet the
        # latter is a valid coarse trapezoidal Wing that the standard tools should
        # refine.
        self._spanwise_mesh: str = "exploded" if explode_into_strips else "trapezoidal"

        # The original, full resolution leading edge and trailing edge curves a Wing was
        # built from. These are populated only by from_edge_points, which sets the
        # spanwise mesh marker to "edge_defined" and stores the curves (together with
        # the tip trim fraction applied while resampling them) so a future non
        # trapezoidal convergence tool can resample them to a different number of
        # WingCrossSections. A Wing built any other way carries no such curves and a
        # None tip trim fraction.
        self._leadingEdgePoints_Wn_Ler: np.ndarray | None = None
        self._trailingEdgePoints_Wn_Ler: np.ndarray | None = None
        self._tip_trim_fraction: float | None = None

        # Validate name and store as immutable.
        self._name = _parameter_validation.str_return_str(name, "name")

        # Validate Ler_Gs_Cgs, store as immutable, and make read-only.
        self._Ler_Gs_Cgs = _parameter_validation.threeD_number_vectorLike_return_float(
            Ler_Gs_Cgs, "Ler_Gs_Cgs"
        )
        self._Ler_Gs_Cgs.flags.writeable = False

        # Validate angles_Gs_to_Wn_ixyz, store as immutable, and make read-only.
        self._angles_Gs_to_Wn_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                angles_Gs_to_Wn_ixyz, "angles_Gs_to_Wn_ixyz"
            )
        )
        if not np.all(
            (-90.0 <= self._angles_Gs_to_Wn_ixyz) & (self._angles_Gs_to_Wn_ixyz <= 90.0)
        ):
            raise ValueError(
                "All elements of angles_Gs_to_Wn_ixyz must lie in the range [-90, "
                "90] degrees."
            )
        self._angles_Gs_to_Wn_ixyz.flags.writeable = False

        # Validate symmetric and mirror_only. These are mutable because
        # Airplane.process_wing_symmetry modifies them for type 5 symmetry.
        symmetric = _parameter_validation.boolLike_return_bool(symmetric, "symmetric")
        mirror_only = _parameter_validation.boolLike_return_bool(
            mirror_only, "mirror_only"
        )
        if symmetric and mirror_only:
            raise ValueError("symmetric and mirror_only cannot both be True.")
        self.symmetric = symmetric
        self.mirror_only = mirror_only

        # Validate symmetryNormal_G and symmetryPoint_G_Cg. These are mutable because
        # Airplane.process_wing_symmetry modifies them for type 5 symmetry.
        if self.symmetric or self.mirror_only:
            if symmetryNormal_G is None:
                raise ValueError(
                    "symmetryNormal_G cannot be None when symmetric or mirror_only is "
                    "True."
                )
            symmetryNormal_G = (
                _parameter_validation.threeD_number_vectorLike_return_float_unit_vector(
                    symmetryNormal_G, "symmetryNormal_G"
                )
            )
            if symmetryPoint_G_Cg is None:
                raise ValueError(
                    "symmetryPoint_G_Cg cannot be None when symmetric or mirror_only "
                    "is True."
                )
            symmetryPoint_G_Cg = (
                _parameter_validation.threeD_number_vectorLike_return_float(
                    symmetryPoint_G_Cg, "symmetryPoint_G_Cg"
                )
            )
        else:
            if symmetryNormal_G is not None:
                raise ValueError(
                    "symmetryNormal_G must be None when both symmetric and "
                    "mirror_only are False."
                )
            if symmetryPoint_G_Cg is not None:
                raise ValueError(
                    "symmetryPoint_G_Cg must be None when both symmetric and "
                    "mirror_only are False."
                )
        self.symmetryNormal_G = symmetryNormal_G
        self.symmetryPoint_G_Cg = symmetryPoint_G_Cg

        # Validate num_chordwise_panels and chordwise_spacing. Store as immutable.
        self._num_chordwise_panels = _parameter_validation.int_in_range_return_int(
            num_chordwise_panels,
            "num_chordwise_panels",
            min_val=1,
            min_inclusive=True,
        )
        if chordwise_spacing not in ["cosine", "uniform"]:
            raise ValueError('chordwise_spacing must be "cosine" or "uniform".')
        self._chordwise_spacing = chordwise_spacing

        # Set once attributes: will be initialized or populated once this Wing's parent
        # Airplane calls generate_mesh.
        self._symmetry_type: int | None = None
        self._num_spanwise_panels: int | None = None
        self._num_panels: int | None = None
        self._panels: np.ndarray | None = None

        # Mutable wake state.
        self.gridWrvp_GP1_CgP1: np.ndarray | None = None

        # Caches for properties derived from immutable attributes. These are populated
        # on first access and preserved in deepcopy.
        self._T_pas_G_Cg_to_Wn_Ler: np.ndarray | None = None
        self._T_pas_Wn_Ler_to_G_Cg: np.ndarray | None = None
        self._WnX_G: np.ndarray | None = None
        self._WnY_G: np.ndarray | None = None
        self._WnZ_G: np.ndarray | None = None
        self._children_T_pas_Wn_Ler_to_Wcs_Lp: list[np.ndarray] | None = None
        self._children_T_pas_Wcs_Lp_to_Wn_Ler: list[np.ndarray] | None = None
        self._children_T_pas_G_Cg_to_Wcs_Lp: list[np.ndarray] | None = None
        self._children_T_pas_Wcs_Lp_to_G_Cg: list[np.ndarray] | None = None

        # Caches for properties derived from set once attributes. These are populated on
        # first access and reset to None in deepcopy.
        self._projected_area: float | None = None
        self._wetted_area: float | None = None
        self._average_panel_aspect_ratio: float | None = None
        self._span: float | None = None
        self._standard_mean_chord: float | None = None
        self._mean_aerodynamic_chord: float | None = None

    # --- Alternate constructor ---
    @classmethod
    def from_edge_points(
        cls,
        leadingEdgePoints_Wn_Ler: np.ndarray | Sequence[Sequence[float | int]],
        trailingEdgePoints_Wn_Ler: np.ndarray | Sequence[Sequence[float | int]],
        num_wing_cross_sections: int,
        airfoil: airfoil_mod.Airfoil,
        name: str = "Untitled Wing",
        Ler_Gs_Cgs: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        angles_Gs_to_Wn_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        symmetric: bool | np.bool_ = False,
        mirror_only: bool | np.bool_ = False,
        symmetryNormal_G: None | np.ndarray | Sequence[float | int] = None,
        symmetryPoint_G_Cg: None | np.ndarray | Sequence[float | int] = None,
        num_chordwise_panels: int = 8,
        chordwise_spacing: str = "cosine",
        tip_trim_fraction: float | int = 0.0,
    ) -> Wing:
        """Builds a planar, untwisted Wing directly from leading edge and trailing edge
        curves, approximating an arbitrary planform with a sequence of WingCrossSections
        joined by single-panel strips.

        The user supplies samples of the leading edge and trailing edge curves (in wing
        axes, relative to the leading edge root point) at whatever resolution they have.
        This method resamples both curves at num_wing_cross_sections points spaced
        uniformly in the spanwise direction, then builds one WingCrossSection from each
        pair of leading and trailing edge points. The standard Wing parameters are
        forwarded to the normal constructor unchanged.

        The resulting Wing records that its spanwise mesh was defined from edge curves
        (its spanwise_mesh is "edge_defined") and stores the original curves so they can
        be resampled to a different number of WingCrossSections later. There is
        intentionally no explode_into_strips parameter: a Wing built from edge points is
        already in single-panel-strip form.

        **Notes:**

        This builds a planar, untwisted Wing. Every leading edge and trailing edge point
        must have a zero z component, every WingCrossSection keeps a zero angle vector,
        and the chord line stays aligned with the wing axes' x direction. Dihedral,
        twist, a three dimensional chord line, per-span Airfoil variation, and control
        surfaces are not supported.

        A planform that tapers to a point at the tip cannot be built directly, because
        the outermost WingCrossSection would have a zero chord, which is invalid. Use
        tip_trim_fraction to drop a fraction of the span off the tip so the outermost
        WingCrossSection has a finite chord. Root trimming is not supported, so the
        leading edge root point must already have a finite root chord.

        :param leadingEdgePoints_Wn_Ler: An array-like of leading edge points (in wing
            axes, relative to the leading edge root point) with shape (M, 3), where M is
            at least 2. These are samples of the leading edge curve and need not match
            num_wing_cross_sections. The first point must be the origin (the leading
            edge root point), the y components must strictly increase from root to tip,
            and every z component must be zero. The units are meters.
        :param trailingEdgePoints_Wn_Ler: An array-like of trailing edge points (in wing
            axes, relative to the leading edge root point) with shape (M, 3), where M is
            at least 2. The leading edge and trailing edge curves are independent and
            may have different point counts. The first point must have a zero y
            component, a zero z component, and a positive x component (the root chord),
            the y components must strictly increase from root to tip, and every z
            component must be zero. The curve must span the same maximum y as the
            leading edge curve. The units are meters.
        :param num_wing_cross_sections: The number of WingCrossSections to generate by
            resampling the curves. It must be an integer of at least 2.
        :param airfoil: The single Airfoil used at every WingCrossSection.
        :param name: A sensible name for the Wing. The default is "Untitled Wing".
        :param Ler_Gs_Cgs: An array-like object of 3 numbers (int or float) representing
            the position of the origin of this Wing's axes (in geometry axes after
            accounting for symmetry, relative to the CG after accounting for symmetry).
            The units are meters. The default is (0.0, 0.0, 0.0).
        :param angles_Gs_to_Wn_ixyz: An array-like object of 3 numbers (int or float)
            representing the angle vector that defines the orientation of this Wing's
            axes relative to the geometry axes (after accounting for symmetry). All
            angles must be in the range [-90, 90] degrees. The units are degrees. The
            default is (0.0, 0.0, 0.0).
        :param symmetric: Set this to True if the Wing's geometry should be mirrored
            across the symmetry plane while retaining the non mirrored side. See the
            class docstring for details. The default is False.
        :param mirror_only: Set this to True if the Wing's geometry should be reflected
            about the symmetry plane without retaining the non reflected geometry. See
            the class docstring for details. The default is False.
        :param symmetryNormal_G: None, or an array-like of 3 numbers (int or float)
            representing the unit normal vector (in geometry axes) of the symmetry
            plane. See the class docstring for details. The default is None.
        :param symmetryPoint_G_Cg: None, or an array-like of 3 numbers (int or float)
            representing a point (in geometry axes, relative to the CG) on the symmetry
            plane. See the class docstring for details. The units are meters. The
            default is None.
        :param num_chordwise_panels: The number of chordwise panels to be used on this
            Wing, which must be a positive integer. The default is 8.
        :param chordwise_spacing: The type of spacing between the Wing's chordwise
            panels. Can be "cosine" or "uniform". The default is "cosine".
        :param tip_trim_fraction: The fraction of the span to drop off the tip before
            resampling, which lets a planform that tapers to a point at the tip be built
            with a finite tip chord. The curves are resampled over the spanwise range
            from the root to (1 - tip_trim_fraction) of the maximum y, so the outermost
            WingCrossSection sits inboard of the geometric tip. It must be a number in
            the range [0, 1). The original, untrimmed curves are still stored. The
            default is 0.0 (no trimming).
        :return: A fully built, validated Wing whose spanwise_mesh is "edge_defined".
        """
        # Stage 1: validate the edge-point inputs structurally and against the planar,
        # untwisted conventions. The remaining invariants (positive chord, the correct
        # single-panel-strip pattern) are enforced downstream by the WingCrossSection
        # and Wing constructors.
        leadingEdgePoints_Wn_Ler, trailingEdgePoints_Wn_Ler = cls._validate_edge_points(
            leadingEdgePoints_Wn_Ler, trailingEdgePoints_Wn_Ler
        )
        num_wing_cross_sections = _parameter_validation.int_in_range_return_int(
            num_wing_cross_sections,
            "num_wing_cross_sections",
            min_val=2,
            min_inclusive=True,
        )
        if not isinstance(airfoil, airfoil_mod.Airfoil):
            raise TypeError("airfoil must be an Airfoil.")
        tip_trim_fraction = _parameter_validation.number_in_range_return_float(
            tip_trim_fraction,
            "tip_trim_fraction",
            min_val=0.0,
            min_inclusive=True,
            max_val=1.0,
            max_inclusive=False,
        )

        # Stage 2: resample both curves at a common set of points spaced uniformly in
        # the spanwise (y) direction, from the root (y of 0) to the trimmed tip. Without
        # trimming the trimmed tip is the shared maximum y. The tip_trim_fraction pulls
        # it inboard so a planform that tapers to a point at the tip ends at a finite
        # chord. The z component is zero at every point in this planar version.
        yTip_Wn_Ler = float(leadingEdgePoints_Wn_Ler[-1, 1])
        yTrimmedTip_Wn_Ler = (1.0 - tip_trim_fraction) * yTip_Wn_Ler
        ys_Wn_Ler = np.linspace(0.0, yTrimmedTip_Wn_Ler, num_wing_cross_sections)
        resampledLeadingEdgePoints_Wn_Ler = cls._resample_edge_points(
            leadingEdgePoints_Wn_Ler, ys_Wn_Ler
        )
        resampledTrailingEdgePoints_Wn_Ler = cls._resample_edge_points(
            trailingEdgePoints_Wn_Ler, ys_Wn_Ler
        )

        # The chord at each WingCrossSection is the trailing edge x value minus the
        # leading edge x value. Reject a non positive chord here with a
        # WingCrossSection-specific message that the downstream WingCrossSection
        # constructor could not phrase as clearly.
        chords = (
            resampledTrailingEdgePoints_Wn_Ler[:, 0]
            - resampledLeadingEdgePoints_Wn_Ler[:, 0]
        )
        for wing_cross_section_id, chord in enumerate(chords):
            if chord <= 0.0:
                raise ValueError(
                    f"The chord at WingCrossSection {wing_cross_section_id} is not "
                    "positive. The trailing edge x value must exceed the leading edge x "
                    "value at every WingCrossSection."
                )

        # A symmetric Wing always resolves to symmetry type 4 or 5, both of which
        # require every WingCrossSection to carry a non None
        # control_surface_symmetry_type, while a non symmetric Wing (types 1 through 3)
        # requires None. These Wings build no control surfaces, so "symmetric" with the
        # default zero deflection adds no geometry. It is just the value that lets a
        # symmetric Wing mesh.
        control_surface_symmetry_type = "symmetric" if symmetric else None

        # Stage 3: build one single-panel WingCrossSection from each resampled point.
        # The parent relative leading point offset is the difference of consecutive
        # resampled leading edge points taken directly in wing axes, which is valid only
        # because every angle vector is zero, so the wing cross section parent axes
        # coincide with the wing axes. The root offset is the zero vector.
        wing_cross_sections = []
        for wing_cross_section_id in range(num_wing_cross_sections):
            is_tip = wing_cross_section_id == num_wing_cross_sections - 1
            if wing_cross_section_id == 0:
                Lp_Wcsp_Lpp = np.array([0.0, 0.0, 0.0])
            else:
                Lp_Wcsp_Lpp = (
                    resampledLeadingEdgePoints_Wn_Ler[wing_cross_section_id]
                    - resampledLeadingEdgePoints_Wn_Ler[wing_cross_section_id - 1]
                )
            wing_cross_sections.append(
                wing_cross_section_mod.WingCrossSection(
                    airfoil=airfoil,
                    num_spanwise_panels=None if is_tip else 1,
                    chord=float(chords[wing_cross_section_id]),
                    Lp_Wcsp_Lpp=Lp_Wcsp_Lpp,
                    angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                    spanwise_spacing=None if is_tip else "uniform",
                    control_surface_symmetry_type=control_surface_symmetry_type,
                )
            )

        # Construct a normal Wing, which marks it "trapezoidal", then record its
        # edge-defined provenance and store the original full resolution curves. These
        # private slots are assigned in-class, the same way __deepcopy__ and
        # deserialization set slots outside __init__. The Wing is never exposed to the
        # user in the intermediate "trapezoidal" state.
        wing = cls(
            wing_cross_sections=wing_cross_sections,
            name=name,
            Ler_Gs_Cgs=Ler_Gs_Cgs,
            angles_Gs_to_Wn_ixyz=angles_Gs_to_Wn_ixyz,
            symmetric=symmetric,
            mirror_only=mirror_only,
            symmetryNormal_G=symmetryNormal_G,
            symmetryPoint_G_Cg=symmetryPoint_G_Cg,
            num_chordwise_panels=num_chordwise_panels,
            chordwise_spacing=chordwise_spacing,
        )
        wing._spanwise_mesh = "edge_defined"
        leadingEdgePoints_Wn_Ler.flags.writeable = False
        trailingEdgePoints_Wn_Ler.flags.writeable = False
        wing._leadingEdgePoints_Wn_Ler = leadingEdgePoints_Wn_Ler
        wing._trailingEdgePoints_Wn_Ler = trailingEdgePoints_Wn_Ler
        wing._tip_trim_fraction = tip_trim_fraction
        return wing

    # --- Deep copy method ---
    def __deepcopy__(self, memo: dict) -> Wing:
        """Creates a deep copy of this Wing, preserving mesh geometry but resetting wake
        state.

        The copy preserves: (1) wing parameters (name, position, angles, symmetry
        settings, panel counts), (2) WingCrossSections (deep copied), (3) mesh metadata
        (symmetry_type, num_spanwise_panels, num_panels), (4) Panels array (each Panel
        is deep copied), and (5) caches for properties derived from immutable attributes
        (transformation matrices, basis vectors, children transformation matrices).

        The copy resets: (1) wake state (gridWrvp_GP1_CgP1 is reset to an empty array
        with correct shape if meshed, or None if not meshed), and (2) caches for
        properties derived from set once attributes (projected_area, wetted_area,
        average_panel_aspect_ratio, span, standard_mean_chord, mean_aerodynamic_chord).

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new Wing with preserved mesh geometry and reset wake state.
        """
        # Create a new Wing instance without calling __init__ to avoid redundant
        # validation and meshing.
        new_wing = object.__new__(Wing)

        # Store this Wing in memo to handle potential circular references.
        memo[id(self)] = new_wing

        # Deep copy the WingCrossSections into a new tuple.
        new_wing._wing_cross_sections = tuple(
            copy.deepcopy(wing_cross_section, memo)
            for wing_cross_section in self._wing_cross_sections
        )

        # Copy immutable Wing parameters (primitive types).
        new_wing._name = self._name
        new_wing._num_chordwise_panels = self._num_chordwise_panels
        new_wing._chordwise_spacing = self._chordwise_spacing
        new_wing._spanwise_mesh = self._spanwise_mesh

        # Copy the stored edge curves (present only for a from_edge_points Wing) and
        # keep the copies read-only, or carry the None.
        if self._leadingEdgePoints_Wn_Ler is not None:
            new_wing._leadingEdgePoints_Wn_Ler = np.copy(self._leadingEdgePoints_Wn_Ler)
            new_wing._leadingEdgePoints_Wn_Ler.flags.writeable = False
        else:
            new_wing._leadingEdgePoints_Wn_Ler = None
        if self._trailingEdgePoints_Wn_Ler is not None:
            new_wing._trailingEdgePoints_Wn_Ler = np.copy(
                self._trailingEdgePoints_Wn_Ler
            )
            new_wing._trailingEdgePoints_Wn_Ler.flags.writeable = False
        else:
            new_wing._trailingEdgePoints_Wn_Ler = None
        new_wing._tip_trim_fraction = self._tip_trim_fraction

        # Copy mutable symmetry attributes (these may be modified by
        # process_wing_symmetry for type 5 symmetry).
        new_wing.symmetric = self.symmetric
        new_wing.mirror_only = self.mirror_only

        # Copy immutable numpy arrays and make them read-only.
        new_wing._Ler_Gs_Cgs = np.copy(self._Ler_Gs_Cgs)
        new_wing._Ler_Gs_Cgs.flags.writeable = False
        new_wing._angles_Gs_to_Wn_ixyz = np.copy(self._angles_Gs_to_Wn_ixyz)
        new_wing._angles_Gs_to_Wn_ixyz.flags.writeable = False

        # Copy mutable symmetry arrays (may be None, may be modified by
        # process_wing_symmetry).
        new_wing.symmetryNormal_G = (
            self.symmetryNormal_G.copy() if self.symmetryNormal_G is not None else None
        )
        new_wing.symmetryPoint_G_Cg = (
            self.symmetryPoint_G_Cg.copy()
            if self.symmetryPoint_G_Cg is not None
            else None
        )

        # Copy set once mesh metadata directly to private attributes (bypassing setters)
        # since we're copying, not setting for the first time.
        new_wing._symmetry_type = self._symmetry_type
        new_wing._num_spanwise_panels = self._num_spanwise_panels
        new_wing._num_panels = self._num_panels

        # Deep copy the Panels array if it exists (directly to private attribute).
        if self._panels is not None:
            new_wing._panels = np.empty_like(self._panels, dtype=object)
            for i in range(self._panels.shape[0]):
                for j in range(self._panels.shape[1]):
                    new_wing._panels[i, j] = copy.deepcopy(self._panels[i, j], memo)
        else:
            new_wing._panels = None

        # Reset wake state to empty arrays with correct shape (if meshed).
        if self._num_spanwise_panels is not None:
            new_wing.gridWrvp_GP1_CgP1 = np.empty(
                (0, self._num_spanwise_panels + 1, 3), dtype=float
            )
        else:
            new_wing.gridWrvp_GP1_CgP1 = None

        # Preserve caches for properties derived from immutable attributes. Copy numpy
        # arrays and make them read-only.
        if self._T_pas_G_Cg_to_Wn_Ler is not None:
            new_wing._T_pas_G_Cg_to_Wn_Ler = self._T_pas_G_Cg_to_Wn_Ler.copy()
            new_wing._T_pas_G_Cg_to_Wn_Ler.flags.writeable = False
        else:
            new_wing._T_pas_G_Cg_to_Wn_Ler = None

        if self._T_pas_Wn_Ler_to_G_Cg is not None:
            new_wing._T_pas_Wn_Ler_to_G_Cg = self._T_pas_Wn_Ler_to_G_Cg.copy()
            new_wing._T_pas_Wn_Ler_to_G_Cg.flags.writeable = False
        else:
            new_wing._T_pas_Wn_Ler_to_G_Cg = None

        if self._WnX_G is not None:
            new_wing._WnX_G = self._WnX_G.copy()
            new_wing._WnX_G.flags.writeable = False
        else:
            new_wing._WnX_G = None

        if self._WnY_G is not None:
            new_wing._WnY_G = self._WnY_G.copy()
            new_wing._WnY_G.flags.writeable = False
        else:
            new_wing._WnY_G = None

        if self._WnZ_G is not None:
            new_wing._WnZ_G = self._WnZ_G.copy()
            new_wing._WnZ_G.flags.writeable = False
        else:
            new_wing._WnZ_G = None

        # Copy list caches (lists of numpy arrays).
        if self._children_T_pas_Wn_Ler_to_Wcs_Lp is not None:
            new_wing._children_T_pas_Wn_Ler_to_Wcs_Lp = []
            for T in self._children_T_pas_Wn_Ler_to_Wcs_Lp:
                T_copy = T.copy()
                T_copy.flags.writeable = False
                new_wing._children_T_pas_Wn_Ler_to_Wcs_Lp.append(T_copy)
        else:
            new_wing._children_T_pas_Wn_Ler_to_Wcs_Lp = None

        if self._children_T_pas_Wcs_Lp_to_Wn_Ler is not None:
            new_wing._children_T_pas_Wcs_Lp_to_Wn_Ler = []
            for T in self._children_T_pas_Wcs_Lp_to_Wn_Ler:
                T_copy = T.copy()
                T_copy.flags.writeable = False
                new_wing._children_T_pas_Wcs_Lp_to_Wn_Ler.append(T_copy)
        else:
            new_wing._children_T_pas_Wcs_Lp_to_Wn_Ler = None

        if self._children_T_pas_G_Cg_to_Wcs_Lp is not None:
            new_wing._children_T_pas_G_Cg_to_Wcs_Lp = []
            for T in self._children_T_pas_G_Cg_to_Wcs_Lp:
                T_copy = T.copy()
                T_copy.flags.writeable = False
                new_wing._children_T_pas_G_Cg_to_Wcs_Lp.append(T_copy)
        else:
            new_wing._children_T_pas_G_Cg_to_Wcs_Lp = None

        if self._children_T_pas_Wcs_Lp_to_G_Cg is not None:
            new_wing._children_T_pas_Wcs_Lp_to_G_Cg = []
            for T in self._children_T_pas_Wcs_Lp_to_G_Cg:
                T_copy = T.copy()
                T_copy.flags.writeable = False
                new_wing._children_T_pas_Wcs_Lp_to_G_Cg.append(T_copy)
        else:
            new_wing._children_T_pas_Wcs_Lp_to_G_Cg = None

        # Reset caches for properties derived from set once attributes.
        new_wing._projected_area = None
        new_wing._wetted_area = None
        new_wing._average_panel_aspect_ratio = None
        new_wing._span = None
        new_wing._standard_mean_chord = None
        new_wing._mean_aerodynamic_chord = None

        return new_wing

    # --- Immutable: read only properties ---
    @property
    def wing_cross_sections(
        self,
    ) -> tuple[wing_cross_section_mod.WingCrossSection, ...]:
        return self._wing_cross_sections

    @property
    def name(self) -> str:
        return self._name

    @property
    def Ler_Gs_Cgs(self) -> np.ndarray:
        return self._Ler_Gs_Cgs

    @property
    def angles_Gs_to_Wn_ixyz(self) -> np.ndarray:
        return self._angles_Gs_to_Wn_ixyz

    @property
    def num_chordwise_panels(self) -> int:
        return self._num_chordwise_panels

    @property
    def chordwise_spacing(self) -> str:
        return self._chordwise_spacing

    @property
    def spanwise_mesh(self) -> str:
        """How this Wing's spanwise mesh was defined.

        It is "trapezoidal" for a Wing whose adjacent WingCrossSections form trapezoidal
        panel strips, "exploded" for a Wing built with explode_into_strips=True, or
        "edge_defined" for a Wing built by from_edge_points. An "exploded" or
        "edge_defined" Wing has every non tip WingCrossSection as a single spanwise
        panel strip. Only an "edge_defined" Wing additionally stores its leading edge
        and trailing edge curves. The standard convergence tools support only
        "trapezoidal" Wings.

        :return: A string, one of "trapezoidal", "exploded", or "edge_defined",
            recording how this Wing's spanwise mesh was defined.
        """
        return self._spanwise_mesh

    @property
    def leadingEdgePoints_Wn_Ler(self) -> np.ndarray | None:
        """The original leading edge curve this Wing was built from, or None.

        It is the full resolution leading edge curve (in wing axes, relative to the
        leading edge root point) passed to from_edge_points, held as a read-only (M, 3)
        ndarray. It is None for any Wing not built by from_edge_points.

        :return: A read-only (M, 3) ndarray of floats representing the leading edge
            points (in wing axes, relative to the leading edge root point), or None. The
            units are meters.
        """
        return self._leadingEdgePoints_Wn_Ler

    @property
    def trailingEdgePoints_Wn_Ler(self) -> np.ndarray | None:
        """The original trailing edge curve this Wing was built from, or None.

        It is the full resolution trailing edge curve (in wing axes, relative to the
        leading edge root point) passed to from_edge_points, held as a read-only (M, 3)
        ndarray. It is None for any Wing not built by from_edge_points.

        :return: A read-only (M, 3) ndarray of floats representing the trailing edge
            points (in wing axes, relative to the leading edge root point), or None. The
            units are meters.
        """
        return self._trailingEdgePoints_Wn_Ler

    @property
    def tip_trim_fraction(self) -> float | None:
        """The fraction of the span dropped off the tip when this Wing was built, or
        None.

        It is the tip_trim_fraction passed to from_edge_points, in the range [0, 1). It
        is None for any Wing not built by from_edge_points.

        :return: A float in the range [0, 1) representing the fraction of the span
            dropped off the tip during resampling, or None.
        """
        return self._tip_trim_fraction

    # --- Immutable derived: manual lazy caching ---
    @property
    def T_pas_G_Cg_to_Wn_Ler(self) -> None | np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        geometry axes relative to the CG to wing axes relative to the leading edge root
        point. Is None if the Wing's symmetry type hasn't been defined yet.

        :return: A (4,4) ndarray of floats representing the transformation matrix or
            None if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._T_pas_G_Cg_to_Wn_Ler is not None:
            return self._T_pas_G_Cg_to_Wn_Ler

        # Step 1: Create T_reflect_pas_G_Cg_to_Gs_Cgs, which maps in homogeneous
        # coordinates from geometry axes relative to the CG to reflected geometry axes
        # (after accounting for symmetry) relative to the CG (after accounting for
        # symmetry). This is the reflection step. Only apply reflection for mirror-only
        # Wings (types 2 and 3), not for symmetric Wings (type 4).
        if self.symmetry_type in (2, 3):
            assert self.symmetryPoint_G_Cg is not None
            assert self.symmetryNormal_G is not None
            T_reflect_pas_G_Cg_to_Gs_Cgs = _transformations.generate_reflect_T(
                plane_point_A_a=self.symmetryPoint_G_Cg,
                plane_normal_A=self.symmetryNormal_G,
                passive=True,
            )
        else:
            T_reflect_pas_G_Cg_to_Gs_Cgs = np.eye(4, dtype=float)

        # Step 2: Create T_trans_pas_Gs_Cgs_to_Gs_Ler, which maps in homogeneous
        # coordinates from geometry axes (after accounting for symmetry) relative to the
        # CG (after accounting for symmetry) to geometry axes (after accounting for
        # symmetry) relative to the leading edge root point. This is the translation
        # step.
        T_trans_pas_Gs_Cgs_to_Gs_Ler = _transformations.generate_trans_T(
            self.Ler_Gs_Cgs, passive=True
        )

        # Step 3: Create T_rot_pas_Gs_to_Wn, which maps in homogeneous coordinates from
        # geometry axes (after accounting for symmetry) to wing axes. This is the
        # rotation step.
        T_rot_pas_Gs_to_Wn = _transformations.generate_rot_T(
            self.angles_Gs_to_Wn_ixyz, passive=True, intrinsic=True, order="xyz"
        )

        # Cache and return the result.
        self._T_pas_G_Cg_to_Wn_Ler = _transformations.compose_T_pas(
            T_reflect_pas_G_Cg_to_Gs_Cgs,
            T_trans_pas_Gs_Cgs_to_Gs_Ler,
            T_rot_pas_Gs_to_Wn,
        )
        self._T_pas_G_Cg_to_Wn_Ler.flags.writeable = False
        return self._T_pas_G_Cg_to_Wn_Ler

    @property
    def T_pas_Wn_Ler_to_G_Cg(self) -> None | np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        wing axes relative to the leading edge root point to geometry axes relative to
        the CG. Is None if the Wing's symmetry type hasn't been defined yet.

        :return: A (4,4) ndarray of floats representing the transformation matrix or
            None if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._T_pas_Wn_Ler_to_G_Cg is not None:
            return self._T_pas_Wn_Ler_to_G_Cg

        _T_pas_G_Cg_to_Wn_Ler = self.T_pas_G_Cg_to_Wn_Ler
        assert _T_pas_G_Cg_to_Wn_Ler is not None

        # Cache and return the result.
        self._T_pas_Wn_Ler_to_G_Cg = _transformations.invert_T_pas(
            _T_pas_G_Cg_to_Wn_Ler
        )
        self._T_pas_Wn_Ler_to_G_Cg.flags.writeable = False
        return self._T_pas_Wn_Ler_to_G_Cg

    @property
    def WnX_G(self) -> None | np.ndarray:
        """The wing axes' first basis vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the wing axes' first basis vector
            (in geometry axes) or None if the Wing's symmetry type hasn't been defined
            yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._WnX_G is not None:
            return self._WnX_G

        WnX_Wn = np.array([1.0, 0.0, 0.0])

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        # Cache and return the result.
        self._WnX_G = _transformations.apply_T_to_vectors(
            _T_pas_Wn_Ler_to_G_Cg, WnX_Wn, is_position=False
        )
        self._WnX_G.flags.writeable = False
        return self._WnX_G

    @property
    def WnY_G(self) -> None | np.ndarray:
        """The wing axes' second basis vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the wing axes' second basis
            vector (in geometry axes) or None if the Wing's symmetry type hasn't been
            defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._WnY_G is not None:
            return self._WnY_G

        WnY_Wn = np.array([0.0, 1.0, 0.0])

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        # Cache and return the result.
        self._WnY_G = _transformations.apply_T_to_vectors(
            _T_pas_Wn_Ler_to_G_Cg, WnY_Wn, is_position=False
        )
        self._WnY_G.flags.writeable = False
        return self._WnY_G

    @property
    def WnZ_G(self) -> None | np.ndarray:
        """The wing axes' third basis vector (in geometry axes).

        :return: A (3,) ndarray of floats representing the wing axes' third basis vector
            (in geometry axes) or None if the Wing's symmetry type hasn't been defined
            yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._WnZ_G is not None:
            return self._WnZ_G

        WnZ_Wn = np.array([0.0, 0.0, 1.0])

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        # Cache and return the result.
        self._WnZ_G = _transformations.apply_T_to_vectors(
            _T_pas_Wn_Ler_to_G_Cg, WnZ_Wn, is_position=False
        )
        self._WnZ_G.flags.writeable = False
        return self._WnZ_G

    @property
    def children_T_pas_Wn_Ler_to_Wcs_Lp(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from wing axes, relative to the leading edge root point, to each of
        this Wing's WingCrossSection's axes, relative to their respective leading
        points.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """
        # Return cached value if available.
        if self._children_T_pas_Wn_Ler_to_Wcs_Lp is not None:
            return self._children_T_pas_Wn_Ler_to_Wcs_Lp

        # Compute, cache, and return the result.
        result = []
        for i in range(len(self.wing_cross_sections)):
            T = _transformations.compose_T_pas(
                *(
                    _assert_T_not_none(wing_cross_section.T_pas_Wcsp_Lpp_to_Wcs_Lp)
                    for wing_cross_section in self.wing_cross_sections[: i + 1]
                )
            )
            T.flags.writeable = False
            result.append(T)
        self._children_T_pas_Wn_Ler_to_Wcs_Lp = result
        return self._children_T_pas_Wn_Ler_to_Wcs_Lp

    @property
    def children_T_pas_Wcs_Lp_to_Wn_Ler(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from each of this Wing's WingCrossSection's axes, relative to their
        respective leading points, to wing axes, relative to the leading edge root
        point.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """
        # Return cached value if available.
        if self._children_T_pas_Wcs_Lp_to_Wn_Ler is not None:
            return self._children_T_pas_Wcs_Lp_to_Wn_Ler

        # Compute, cache, and return the result.
        result = []
        for i in range(len(self.wing_cross_sections)):
            T = _transformations.invert_T_pas(self.children_T_pas_Wn_Ler_to_Wcs_Lp[i])
            T.flags.writeable = False
            result.append(T)
        self._children_T_pas_Wcs_Lp_to_Wn_Ler = result
        return self._children_T_pas_Wcs_Lp_to_Wn_Ler

    @property
    def children_T_pas_G_Cg_to_Wcs_Lp(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from geometry axes, relative to the CG, to each of this Wing's
        WingCrossSection's axes, relative to their respective leading points.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """
        # Return cached value if available.
        if self._children_T_pas_G_Cg_to_Wcs_Lp is not None:
            return self._children_T_pas_G_Cg_to_Wcs_Lp

        _T_pas_G_Cg_to_Wn_Ler = self.T_pas_G_Cg_to_Wn_Ler
        assert _T_pas_G_Cg_to_Wn_Ler is not None

        # Compute, cache, and return the result.
        result = []
        for i in range(len(self.wing_cross_sections)):
            T = _transformations.compose_T_pas(
                _T_pas_G_Cg_to_Wn_Ler, self.children_T_pas_Wn_Ler_to_Wcs_Lp[i]
            )
            T.flags.writeable = False
            result.append(T)
        self._children_T_pas_G_Cg_to_Wcs_Lp = result
        return self._children_T_pas_G_Cg_to_Wcs_Lp

    @property
    def children_T_pas_Wcs_Lp_to_G_Cg(self) -> list[np.ndarray]:
        """A list of passive transformation matrices which map in homogeneous
        coordinates from each of this Wing's WingCrossSection's axes, relative to their
        respective leading points, to geometry axes, relative to the CG.

        :return: A list of (4,4) ndarrays of floats representing the homogeneous
            transformation matrices.
        """
        # Return cached value if available.
        if self._children_T_pas_Wcs_Lp_to_G_Cg is not None:
            return self._children_T_pas_Wcs_Lp_to_G_Cg

        # Compute, cache, and return the result.
        result = []
        for i in range(len(self.wing_cross_sections)):
            T = _transformations.invert_T_pas(self.children_T_pas_G_Cg_to_Wcs_Lp[i])
            T.flags.writeable = False
            result.append(T)
        self._children_T_pas_Wcs_Lp_to_G_Cg = result
        return self._children_T_pas_Wcs_Lp_to_G_Cg

    # --- Set once: properties with single assignment enforcement ---
    @property
    def symmetry_type(self) -> int | None:
        """The symmetry type of this Wing.

        :return: An integer from 1-4 representing the symmetry type, or None if the
            Wing's symmetry type hasn't been determined yet.
        """
        return self._symmetry_type

    @symmetry_type.setter
    def symmetry_type(self, value: int) -> None:
        if self._symmetry_type is not None:
            raise AttributeError("symmetry_type can only be set once")
        self._symmetry_type = value

    @property
    def num_spanwise_panels(self) -> int | None:
        """The total number of spanwise Panels on this Wing.

        :return: A positive integer representing the number of spanwise Panels, or None
            if the Wing hasn't been meshed yet.
        """
        return self._num_spanwise_panels

    @num_spanwise_panels.setter
    def num_spanwise_panels(self, value: int) -> None:
        if self._num_spanwise_panels is not None:
            raise AttributeError("num_spanwise_panels can only be set once")
        self._num_spanwise_panels = value

    @property
    def num_panels(self) -> int | None:
        """The total number of Panels on this Wing.

        :return: A positive integer representing the total number of Panels, or None if
            the Wing hasn't been meshed yet.
        """
        return self._num_panels

    @num_panels.setter
    def num_panels(self, value: int) -> None:
        if self._num_panels is not None:
            raise AttributeError("num_panels can only be set once")
        self._num_panels = value

    @property
    def panels(self) -> np.ndarray | None:
        """The 2D array of Panels on this Wing.

        :return: A (num_chordwise_panels, num_spanwise_panels) ndarray of Panel objects,
            or None if the Wing hasn't been meshed yet.
        """
        return self._panels

    @panels.setter
    def panels(self, value: np.ndarray) -> None:
        if self._panels is not None:
            raise AttributeError("panels can only be set once")
        self._panels = value

    # --- Set once derived: manual lazy caching ---
    @property
    def projected_area(self) -> None | float:
        """The area of the Wing projected onto the plane defined by the wing axes' xy
        plane.

        **Notes:**

        If the Wing is symmetric and continuous, the area of the mirrored half is
        included.

        :return: The projected area of the Wing. It has units of square meters. If the
            Wing hasn't been meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self._panels is None:
            return None

        # Return cached value if available.
        if self._projected_area is not None:
            return self._projected_area

        projected_area = 0.0

        # Get the wing Z-axis once before iterating.
        WnZ_G = self.WnZ_G
        assert WnZ_G is not None

        # Iterate through the chordwise and spanwise indices of the Panels and add their
        # area to the total projected area.
        assert self._num_spanwise_panels is not None
        for chordwise_location in range(self._num_chordwise_panels):
            for spanwise_location in range(self._num_spanwise_panels):
                this_panel: _panel.Panel = self._panels[
                    chordwise_location, spanwise_location
                ]
                projected_area += this_panel.calculate_projected_area(WnZ_G)

        # Cache the computed value.
        self._projected_area = projected_area
        return projected_area

    @property
    def wetted_area(self) -> None | float:
        """The Wing's wetted area.

        **Notes:**

        If the Wing is symmetric and continuous, the area of the mirrored half is
        included.

        :return: The wetted area of the Wing. It has units of square meters. If the Wing
            hasn't been meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self._panels is None:
            return None

        # Return cached value if available.
        if self._wetted_area is not None:
            return self._wetted_area

        wetted_area = 0.0

        # Iterate through the chordwise and spanwise indices of the panels and add their
        # area to the total wetted area.
        assert self._num_spanwise_panels is not None
        for chordwise_location in range(self._num_chordwise_panels):
            for spanwise_location in range(self._num_spanwise_panels):
                this_panel: _panel.Panel = self._panels[
                    chordwise_location, spanwise_location
                ]
                wetted_area += this_panel.area

        # Cache the computed value.
        self._wetted_area = wetted_area
        return wetted_area

    @property
    def average_panel_aspect_ratio(self) -> None | float:
        """The average aspect ratio of the Wing's Panels.

        :return: The average aspect ratio of the Wing's Panels. If the Wing hasn't been
            meshed yet, None is returned instead.
        """
        # Return None if the Wing hasn't been meshed yet.
        if self._panels is None:
            return None

        # Return cached value if available.
        if self._average_panel_aspect_ratio is not None:
            return self._average_panel_aspect_ratio

        aspect_ratio_sum = 0.0

        # Iterate through the chordwise and spanwise indices of the Panels and sum all
        # the Panels' aspect ratios.
        assert self._num_spanwise_panels is not None
        for chordwise_location in range(self._num_chordwise_panels):
            for spanwise_location in range(self._num_spanwise_panels):
                this_panel: _panel.Panel = self._panels[
                    chordwise_location, spanwise_location
                ]
                aspect_ratio_sum += this_panel.aspect_ratio

        assert self._num_panels is not None

        # Cache and return the result.
        self._average_panel_aspect_ratio = aspect_ratio_sum / self._num_panels
        return self._average_panel_aspect_ratio

    @property
    def span(self) -> None | float:
        """The Wing's span.

        **Notes:**

        The span is derived by first finding the vector connecting the leading edges of
        the root and tip WingCrossSections. Then, this vector is projected onto the wing
        axes' second basis vector. The span is defined as the magnitude of this
        projection.

        If the Wing is symmetric and continuous, this method includes the span of the
        mirrored half.

        :return: The Wing's span. It has units of meters. None is returned if the Wing's
            symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self._symmetry_type is None:
            return None

        # Return cached value if available.
        if self._span is not None:
            return self._span

        tipLp_Wcsp_Lpp = self._wing_cross_sections[-1].Lp_Wcsp_Lpp

        tip_T_pas_Wcsp_Lpp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[-2]

        tipLp_Wn_Ler = _transformations.apply_T_to_vectors(
            tip_T_pas_Wcsp_Lpp_to_Wn_Ler, tipLp_Wcsp_Lpp, is_position=True
        )

        # Project the tip position onto the wing axes' y direction (spanwise direction).
        projectedTipLp_Wn_Ler = np.dot(
            tipLp_Wn_Ler, np.array([0.0, 1.0, 0.0])
        ) * np.array([0.0, 1.0, 0.0])

        span = float(np.linalg.norm(projectedTipLp_Wn_Ler))

        # If the wing is symmetric and continuous, multiply the span by two.
        if self._symmetry_type == 4:
            span *= 2

        # Cache the computed value.
        self._span = span
        return span

    @property
    def standard_mean_chord(self) -> None | float:
        """The Wing's standard mean chord.

        **Notes:**

        The standard mean chord is defined as the projected area divided by the span.
        See their respective methods for the definitions of span and projected area.

        :return: The standard mean chord of the Wing. It has units of meters. None is
            returned if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._standard_mean_chord is not None:
            return self._standard_mean_chord

        _projected_area = self.projected_area
        assert _projected_area is not None

        _span = self.span
        assert _span is not None

        # Cache and return the result.
        self._standard_mean_chord = _projected_area / _span
        return self._standard_mean_chord

    @property
    def mean_aerodynamic_chord(self) -> None | float:
        """The Wing's mean aerodynamic chord.

        :return: The mean aerodynamic chord of the Wing. It has units of meters. None is
            returned if the Wing's symmetry type hasn't been defined yet.
        """
        # If the Wing's symmetry type hasn't been set yet, return None to avoid
        # incorrect symmetry handling.
        if self.symmetry_type is None:
            return None

        # Return cached value if available.
        if self._mean_aerodynamic_chord is not None:
            return self._mean_aerodynamic_chord

        # This method is based on the equation for the mean aerodynamic chord of a wing,
        # which can be found here:
        # https://en.wikipedia.org/wiki/Chord_(aeronautics)#Mean_aerodynamic_chord. This
        # equation integrates the squared chord from the Wing's center to the Wing's
        # tip. We will perform this integral piecewise for each section of the Wing.
        integral = 0.0

        # Iterate through the WingCrossSections to add the contribution of their
        # corresponding Wing section to the piecewise integral.
        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections[:-1]
        ):
            next_wing_cross_section = self.wing_cross_sections[
                wing_cross_section_id + 1
            ]

            chord = wing_cross_section.chord
            next_chord = next_wing_cross_section.chord

            # Find this section's span by calculating the positions of both
            # WingCrossSections in wing axes, then finding the distance between them.

            # Calculate current WingCrossSection's position in wing axes.
            Lp_Wcs_Lp = np.array([0.0, 0.0, 0.0])

            T_pas_Wcs_Lp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[
                wing_cross_section_id
            ]

            Lp_Wn_Ler = _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_Wn_Ler, Lp_Wcs_Lp, is_position=True
            )

            # Calculate next WingCrossSection's position in wing axes.
            nextLp_nextWcs_nextLp = np.array([0.0, 0.0, 0.0])

            T_pas_nextWcs_nextLp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[
                wing_cross_section_id + 1
            ]

            nextLp_Wn_Ler = _transformations.apply_T_to_vectors(
                T_pas_nextWcs_nextLp_to_Wn_Ler, nextLp_nextWcs_nextLp, is_position=True
            )

            # Find the section vector and project it onto spanwise direction (wing axes
            # y direction).
            nextLp_Wn_Lp = nextLp_Wn_Ler - Lp_Wn_Ler

            nextLpProj_Wn_Lp = np.dot(
                nextLp_Wn_Lp, np.array([0.0, 1.0, 0.0])
            ) * np.array([0.0, 1.0, 0.0])

            section_span = float(np.linalg.norm(nextLpProj_Wn_Lp))

            # Each Wing section is, by definition, trapezoidal (at least when projected
            # on to the wing axes' xy plane). For a trapezoid, the integral from the
            # cited equation can be shown to evaluate to the following.
            integral += (
                section_span * (chord**2 + chord * next_chord + next_chord**2) / 3
            )

        _projected_area = self.projected_area
        assert _projected_area is not None

        # Multiply the integral's value by the coefficients from the cited equation.
        # Double if the wing is symmetric and continuous.
        if self.symmetry_type == 4:
            self._mean_aerodynamic_chord = 2 * integral / _projected_area
        else:
            self._mean_aerodynamic_chord = integral / _projected_area
        return self._mean_aerodynamic_chord

    # --- Other methods ---
    def generate_mesh(self, symmetry_type: int) -> None:
        """Generates this Wing's mesh, which finishes the process of preparing the Wing
        to be used in a simulation. It is called by the Wing's parent Airplane, after
        it's determined its symmetry type.

        :param symmetry_type: The symmetry type of this Wing as an int from 1-4. See the
            class docstring for details on how to interpret the symmetry types.
        :return: None
        """
        # Validate and apply symmetry_type. 5 isn't a valid symmetry type, because the
        # parent Airplane should have modified a Wing that initially had type 5 symmetry
        # to have type 1 symmetry, and then made a new reflected Wing with type 3
        # symmetry.
        validated_symmetry_type = _parameter_validation.int_in_range_return_int(
            symmetry_type,
            "symmetry_type",
            min_val=1,
            min_inclusive=True,
            max_val=4,
            max_inclusive=True,
        )
        self.symmetry_type = validated_symmetry_type

        # Set this Wing's children WingCrossSections' symmetry type parameters.
        for wing_cross_section in self.wing_cross_sections:
            wing_cross_section.symmetry_type = validated_symmetry_type

        # Find the number of spanwise Panels on this Wing by adding each
        # WingCrossSection's number of spanwise Panels. Exclude the last
        # WingCrossSection's number of spanwise Panels as this is irrelevant. If this
        # Wing has type 4 symmetry multiply the summation by two.
        computed_num_spanwise_panels = 0
        for wing_cross_section in self.wing_cross_sections[:-1]:
            assert wing_cross_section.num_spanwise_panels is not None
            computed_num_spanwise_panels += wing_cross_section.num_spanwise_panels
        if validated_symmetry_type == 4:
            computed_num_spanwise_panels *= 2
        self.num_spanwise_panels = computed_num_spanwise_panels

        # Calculate the number of panels on this wing.
        self.num_panels = computed_num_spanwise_panels * self.num_chordwise_panels

        # Initialize an empty array to hold this Wing's wake ring vortex points.
        self.gridWrvp_GP1_CgP1 = np.empty(
            (0, computed_num_spanwise_panels + 1, 3), dtype=float
        )

        # Generate the wing's mesh, which populates the Panels attribute.
        _meshing.mesh_wing(self)

    def get_plottable_data(
        self, show: bool | np.bool_ = False
    ) -> list[list[np.ndarray]] | None:
        """Returns plottable data for this Wing's Airfoils' outlines and mean camber
        lines.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: None if the Wing's symmetry type hasn't been set yet, or if show is
            True. Otherwise, returns a list of two lists, each containing one ndarray
            for every one of this Wing's Airfoils. These ndarrays represent points on
            each Airfoil's outline and mean camber lines, respectively. The points are
            in wing axes, relative to the leading edge root point. The units are in
            meters.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        # If this Wing hasn't had its symmetry type set, return None.
        if self.symmetry_type is None:
            return None

        airfoilOutlines_Wn_ler = []
        airfoilMcls_Wn_ler = []
        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections
        ):
            plottable_data = wing_cross_section.get_plottable_data(show=False)
            assert plottable_data is not None

            [airfoilOutline_Wcs_lp, airfoilMcl_Wcs_lp] = plottable_data

            T_pas_Wcs_Lp_to_Wn_Ler = self.children_T_pas_Wcs_Lp_to_Wn_Ler[
                wing_cross_section_id
            ]

            airfoilOutline_Wn_ler = _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_Wn_Ler, airfoilOutline_Wcs_lp, is_position=True
            )
            airfoilMcl_Wn_ler = _transformations.apply_T_to_vectors(
                T_pas_Wcs_Lp_to_Wn_Ler, airfoilMcl_Wcs_lp, is_position=True
            )

            airfoilOutlines_Wn_ler.append(airfoilOutline_Wn_ler)
            airfoilMcls_Wn_ler.append(airfoilMcl_Wn_ler)

        if not show:
            return [airfoilOutlines_Wn_ler, airfoilMcls_Wn_ler]

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

        _T_pas_G_Cg_to_Wn_Ler = self.T_pas_G_Cg_to_Wn_Ler
        assert _T_pas_G_Cg_to_Wn_Ler is not None

        AxesWLerWcs1Lp1_G_Cg = pv.AxesAssembly(
            x_label="WX@Ler/Wcs1XLp1",
            y_label="WY@Ler/Wcs1YLp1",
            z_label="WZ@Ler/Wcs1ZLp1",
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
            user_matrix=np.linalg.inv(_T_pas_G_Cg_to_Wn_Ler),
            # user_matrix=wingAxes_T_act,
            name="W/Wcs1",
            shaft_type="cylinder",
            shaft_radius=0.025,
            shaft_length=(0.8, 0.8, 0.8),
            tip_type="cone",
            tip_radius=0.1,
            tip_length=(0.2, 0.2, 0.2),
            symmetric_bounds=False,
        )

        plotter.add_actor(AxesWLerWcs1Lp1_G_Cg)

        _T_pas_Wn_Ler_to_G_Cg = self.T_pas_Wn_Ler_to_G_Cg
        assert _T_pas_Wn_Ler_to_G_Cg is not None

        for wing_cross_section_id, wing_cross_section in enumerate(
            self.wing_cross_sections
        ):
            airfoilOutline_Wn_ler = airfoilOutlines_Wn_ler[wing_cross_section_id]
            airfoilMcl_Wn_ler = airfoilMcls_Wn_ler[wing_cross_section_id]

            airfoilOutline_G_Cg = _transformations.apply_T_to_vectors(
                _T_pas_Wn_Ler_to_G_Cg, airfoilOutline_Wn_ler, is_position=True
            )
            airfoilMcl_G_Cg = _transformations.apply_T_to_vectors(
                _T_pas_Wn_Ler_to_G_Cg, airfoilMcl_Wn_ler, is_position=True
            )

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
                    x_label=f"Wcs{wing_cross_section_num}X@Lp{wing_cross_section_num}",
                    y_label=f"Wcs{wing_cross_section_num}Y@Lp{wing_cross_section_num}",
                    z_label=f"Wcs{wing_cross_section_num}Z@Lp{wing_cross_section_num}",
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
                        self.children_T_pas_G_Cg_to_Wcs_Lp[wing_cross_section_id]
                    ),
                    name=f"Wcs{wing_cross_section_id}",
                    shaft_type="cylinder",
                    shaft_radius=0.025,
                    shaft_length=(0.8, 0.8, 0.8),
                    tip_type="cone",
                    tip_radius=0.1,
                    tip_length=(0.2, 0.2, 0.2),
                    symmetric_bounds=False,
                )

                plotter.add_actor(AxesWcsLp_G_Cg)

        if self.panels is not None:
            # Initialize empty arrays to hold the Panels' vertices and faces.
            panel_vertices = np.empty((0, 3), dtype=float)
            panel_faces = np.empty(0, dtype=int)

            # Initialize a variable to keep track of how many Panels' data has been
            # added to the arrays.
            panel_num = 0

            # Unravel the Wing's Panel matrix and iterate through it.
            panels = np.ravel(self.panels)
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

    @staticmethod
    def _validate_edge_points(
        leadingEdgePoints_Wn_Ler: np.ndarray | Sequence[Sequence[float | int]],
        trailingEdgePoints_Wn_Ler: np.ndarray | Sequence[Sequence[float | int]],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Validates the leading edge and trailing edge curves passed to
        from_edge_points, returning them as (M, 3) float arrays.

        This is the preliminary, structural layer of validation plus the planar,
        untwisted convention checks that the downstream WingCrossSection and Wing
        constructors cannot phrase in terms the user would recognize.

        :param leadingEdgePoints_Wn_Ler: The leading edge curve, an array-like of points
            (in wing axes, relative to the leading edge root point) of shape (M, 3).
        :param trailingEdgePoints_Wn_Ler: The trailing edge curve, an array-like of
            points (in wing axes, relative to the leading edge root point) of shape (M,
            3).
        :return: A tuple of the validated leading edge and trailing edge curves, each a
            writeable (M, 3) ndarray of floats.
        :raises ValueError: If either curve is not a valid (M, 3) array with M at least
            two, has non strictly increasing y components, is not anchored at the root,
            has a non zero z component anywhere, or the two curves do not span the same
            maximum y.
        """
        leadingEdgePoints_Wn_Ler = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                leadingEdgePoints_Wn_Ler, "leadingEdgePoints_Wn_Ler"
            )
        )
        trailingEdgePoints_Wn_Ler = (
            _parameter_validation.arrayLike_of_threeD_number_vectorLikes_return_float(
                trailingEdgePoints_Wn_Ler, "trailingEdgePoints_Wn_Ler"
            )
        )

        for points, points_name in (
            (leadingEdgePoints_Wn_Ler, "leadingEdgePoints_Wn_Ler"),
            (trailingEdgePoints_Wn_Ler, "trailingEdgePoints_Wn_Ler"),
        ):
            if points.ndim != 2:
                raise ValueError(f"{points_name} must have shape (M, 3).")
            if points.shape[0] < 2:
                raise ValueError(f"{points_name} must contain at least two points.")
            if not np.all(np.diff(points[:, 1]) > 0.0):
                raise ValueError(
                    f"{points_name}'s y components must strictly increase from root to "
                    "tip."
                )
            if not np.all(points[:, 2] == 0.0):
                raise ValueError(
                    f"Every point in {points_name} must have a zero z component, "
                    "because from_edge_points builds only planar Wings."
                )

        # The leading edge curve starts at the leading edge root point, which is the
        # origin by definition.
        if not np.all(leadingEdgePoints_Wn_Ler[0] == 0.0):
            raise ValueError(
                "leadingEdgePoints_Wn_Ler's first point must be the origin (the leading "
                "edge root point)."
            )

        # The trailing edge curve starts at the root chord: a zero y component, a zero z
        # component, and a positive x component.
        rootTp_Wn_Ler = trailingEdgePoints_Wn_Ler[0]
        if not (
            rootTp_Wn_Ler[1] == 0.0
            and rootTp_Wn_Ler[2] == 0.0
            and rootTp_Wn_Ler[0] > 0.0
        ):
            raise ValueError(
                "trailingEdgePoints_Wn_Ler's first point must have a zero y component, a "
                "zero z component, and a positive x component (the root chord)."
            )

        # Both curves must reach the same tip so they can be sampled at a common set of
        # points.
        if leadingEdgePoints_Wn_Ler[-1, 1] != trailingEdgePoints_Wn_Ler[-1, 1]:
            raise ValueError(
                "leadingEdgePoints_Wn_Ler and trailingEdgePoints_Wn_Ler must span the "
                "same maximum y."
            )

        return leadingEdgePoints_Wn_Ler, trailingEdgePoints_Wn_Ler

    @staticmethod
    def _resample_edge_points(
        edgePoints_Wn_Ler: np.ndarray, ys_Wn_Ler: np.ndarray
    ) -> np.ndarray:
        """Resamples an edge curve at the requested spanwise y values using monotone
        PCHIP interpolation of the x component as a function of the y component.

        :param edgePoints_Wn_Ler: A validated (M, 3) ndarray of edge points (in wing
            axes, relative to the leading edge root point), with strictly increasing y
            components and zero z components.
        :param ys_Wn_Ler: A (N,) ndarray of the spanwise (wing axes y) values to sample
            at, relative to the leading edge root point.
        :return: A (N, 3) ndarray of resampled edge points (in wing axes, relative to
            the leading edge root point). The z component is zero at every point.
        """
        x_interpolator = sp_interp.PchipInterpolator(
            edgePoints_Wn_Ler[:, 1], edgePoints_Wn_Ler[:, 0]
        )
        xs_Wn_Ler = x_interpolator(ys_Wn_Ler)
        zs_Wn_Ler = np.zeros_like(ys_Wn_Ler)
        return np.column_stack((xs_Wn_Ler, ys_Wn_Ler, zs_Wn_Ler))

    @staticmethod
    def _interpolate_between_wing_cross_sections(
        first_wing_cross_section: wing_cross_section_mod.WingCrossSection,
        second_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    ) -> list[wing_cross_section_mod.WingCrossSection]:
        """Build the single-panel WingCrossSections spanning from just past
        first_wing_cross_section through second_wing_cross_section.

        When exploding a wing to 1 spanwise panel per cross section, we need to
        interpolate the intermediate cross sections. This returns the N sections
        downstream of first_wing_cross_section, where N is first_wing_cross_section's
        spanwise panel count. The section at first_wing_cross_section itself is seeded
        once by _explode_wing, so it is not repeated here.

        :param first_wing_cross_section: The first WingCrossSection.
        :param second_wing_cross_section: The second WingCrossSection.
        :return: A list of WingCrossSections representing the interpolated cross
            sections
        """
        interpolated = []

        N = first_wing_cross_section.num_spanwise_panels
        assert (
            N is not None
        ), "first_wing_cross_section.num_spanwise_panels must not be None"

        for i in range(N):
            t = (i + 1) / N  # interpolation parameter between 0 and 1

            chord = (
                1 - t
            ) * first_wing_cross_section.chord + t * second_wing_cross_section.chord
            # Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_ixyz are parent-relative deltas (see
            # WingCrossSection's constructor docstring), so each of the N intermediates
            # carries 1.0 / N of the second WingCrossSection's delta and the chain
            # composes back to the second WingCrossSection's offset and twist. This is
            # exact only when the intermediates are uniformly spaced, which
            # _explode_wing enforces by validating the input.
            Lp_Wcsp_Lpp = tuple(np.array(second_wing_cross_section.Lp_Wcsp_Lpp) / N)
            angles_Wcsp_to_Wcs_ixyz = (
                second_wing_cross_section.angles_Wcsp_to_Wcs_ixyz / N
            )
            is_final_section = (
                second_wing_cross_section.num_spanwise_panels is None and i == N - 1
            )

            interpolated.append(
                wing_cross_section_mod.WingCrossSection(
                    num_spanwise_panels=None if is_final_section else 1,
                    chord=chord,
                    Lp_Wcsp_Lpp=Lp_Wcsp_Lpp,
                    angles_Wcsp_to_Wcs_ixyz=angles_Wcsp_to_Wcs_ixyz,
                    control_surface_symmetry_type=first_wing_cross_section.control_surface_symmetry_type,
                    control_surface_hinge_point=first_wing_cross_section.control_surface_hinge_point,
                    control_surface_deflection=first_wing_cross_section.control_surface_deflection,
                    spanwise_spacing=None if is_final_section else "uniform",
                    airfoil=first_wing_cross_section.airfoil,
                )
            )
        return interpolated

    def _explode_wing(
        self, wing_cross_sections: list[wing_cross_section_mod.WingCrossSection]
    ) -> list[wing_cross_section_mod.WingCrossSection]:
        """Takes a list of WingCrossSections and returns a new list where all cross
        sections have num_spanwise_panels = 1.

        The interpolation distributes each non-tip WingCrossSection's parent-relative
        offset and twist uniformly across its N intermediates, which is exact only
        when the source WingCrossSection specifies uniform spanwise spacing. Non
        uniform spacings (for example "cosine") would silently be replaced with
        uniform spacing in the output mesh, so this method rejects them up front.

        :param wing_cross_sections: The list of wing cross sections to explode. Every
            non tip WingCrossSection must have spanwise_spacing="uniform".
        :return: A new list of exploded wing cross sections.
        :raises ValueError: If any non tip WingCrossSection has a spanwise_spacing
            other than "uniform".
        """
        for i, wing_cross_section in enumerate(wing_cross_sections[:-1]):
            if wing_cross_section.spanwise_spacing != "uniform":
                raise ValueError(
                    f"wing_cross_sections[{i}].spanwise_spacing is "
                    f'"{wing_cross_section.spanwise_spacing}", but exploding a Wing '
                    f'(explode_into_strips=True) requires "uniform" on every non tip '
                    f"WingCrossSection."
                )

        # Seed the result with the root strip, then fill in the single-panel sections
        # spanning each adjacent pair. Seeding the root here, rather than inside the
        # per-pair interpolation, keeps that helper independent of its position in the
        # wing.
        root = wing_cross_sections[0]
        new_cross_sections = [
            wing_cross_section_mod.WingCrossSection(
                num_spanwise_panels=1,
                chord=root.chord,
                Lp_Wcsp_Lpp=root.Lp_Wcsp_Lpp,
                angles_Wcsp_to_Wcs_ixyz=root.angles_Wcsp_to_Wcs_ixyz,
                control_surface_symmetry_type=root.control_surface_symmetry_type,
                control_surface_hinge_point=root.control_surface_hinge_point,
                control_surface_deflection=root.control_surface_deflection,
                spanwise_spacing="uniform",
                airfoil=root.airfoil,
            )
        ]

        for i in range(len(wing_cross_sections) - 1):
            new_cross_sections.extend(
                self._interpolate_between_wing_cross_sections(
                    wing_cross_sections[i], wing_cross_sections[i + 1]
                )
            )

        return new_cross_sections


def _assert_T_not_none(T: np.ndarray | None) -> np.ndarray:
    """Assert that a transformation matrix is not None and return it.

    :param T: None, or a (4,4) ndarray of floats representing the transformation matrix.
    :return: A (4,4) ndarray of floats representing the transformation matrix
        (guaranteed not to be None).
    """
    assert T is not None
    return T
