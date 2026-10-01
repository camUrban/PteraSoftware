"""Contains the WingCrossSection class."""

from __future__ import annotations

import copy
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pyvista as pv

from .. import _output_rendering, _parameter_validation, _transformations
from . import airfoil as airfoil_mod


class WingCrossSection:
    """A class used to contain the wing cross sections of a Wing.

    Immutable attributes (airfoil, num_spanwise_panels, chord, Lp_Wcsp_Lpp,
    angles_Wcsp_to_Wcs_ixyz, control_surface_hinge_point, control_surface_deflection,
    and spanwise_spacing) are set during initialization and cannot be modified
    afterward. The numpy arrays Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_ixyz are made read
    only to prevent in place mutation.

    Derived transformation matrices (T_pas_Wcsp_Lpp_to_Wcs_Lp and
    T_pas_Wcs_Lp_to_Wcsp_Lpp) are lazily evaluated and cached.

    The validated and symmetry_type attributes are set once by the parent Wing and
    cannot be modified after being set.

    The control_surface_symmetry_type attribute remains mutable as it may be modified by
    Airplane.process_wing_symmetry() for type 5 symmetry handling.

    The first WingCrossSection in a Wing's wing_cross_section list is known as the root
    WingCrossSection. The last is known as the tip WingCrossSection.

    Every WingCrossSection has its own wing cross section axes. For root
    WingCrossSections, their wing cross section axes are identical in position,
    orientation, and handedness to their Wing's wing axes. For all other
    WingCrossSections, their wing cross section axes are defined relative to the axes of
    the previous WingCrossSection. Locally, the x axis points from a cross section's
    leading point to its trailing point, the y axis points spanwise in the general
    direction of the next WingCrossSection, and the z axis points upwards.

    Things can get a little confusing with respect to WingCrossSections for Wings with
    symmetric or mirror_only set to True. For more details, look in the Wing class's
    docstring. Also, remember that WingCrossSections themselves aren't used by the
    solvers, they are merely one of the Wing attributes that tell the meshing function
    how we'd like to generate the Wing's Panels.

    **Citation:**

    Adapted from: geometry.WingXSec in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/26/2020
    """

    __slots__ = (
        "_airfoil",
        "_num_spanwise_panels",
        "_chord",
        "_Lp_Wcsp_Lpp",
        "_angles_Wcsp_to_Wcs_ixyz",
        "control_surface_symmetry_type",
        "_control_surface_hinge_point",
        "_control_surface_deflection",
        "_spanwise_spacing",
        "_validated",
        "_symmetry_type",
        "_T_pas_Wcsp_Lpp_to_Wcs_Lp",
        "_T_pas_Wcs_Lp_to_Wcsp_Lpp",
    )

    def __init__(
        self,
        airfoil: airfoil_mod.Airfoil,
        num_spanwise_panels: int | None,
        chord: float | int = 1.0,
        Lp_Wcsp_Lpp: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_ixyz: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        control_surface_symmetry_type: str | None = None,
        control_surface_hinge_point: float = 0.75,
        control_surface_deflection: float | int = 0.0,
        spanwise_spacing: str | None = None,
    ) -> None:
        """The initialization method.

        :param airfoil: The Airfoil to be used at this WingCrossSection.
        :param num_spanwise_panels: The number of spanwise Panels to be used between
            this WingCrossSection and the next one. For tip WingCrossSections, it must
            be None. For all other WingCrossSections, it must be a positive integer.
        :param chord: The Wing's chord at this WingCrossSection. It must be greater than
            0.0 and a number (int or float), and will be converted internally to a
            float. The units are in meters. The default value is 1.0.
        :param Lp_Wcsp_Lpp: An array-like object of 3 numbers (int or float)
            representing the position in meters of this WingCrossSection's leading edge
            in parent wing cross section axes, relative to the parent leading edge
            point. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. If this is the root WingCrossSection, the parent wing cross
            section axes are the wing axes and the parent leading point is the Wing's
            leading edge root point. If not, the parent axes and point are those of the
            previous WingCrossSection. If this is the root WingCrossSection, it must be
            a zero vector. The second component must be non negative. For non root
            WingCrossSections, the second component must be strictly positive. The units
            are in meters. The default is (0.0, 0.0, 0.0).
        :param angles_Wcsp_to_Wcs_ixyz: An array-like object of 3 numbers (int or float)
            representing the angle vector of rotation angles that define the orientation
            of this WingCrossSection's axes relative to the parent wing cross section
            axes. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. If this is a root WingCrossSection, these are the wing axes. If
            not, the parent axes are the previous WingCrossSection's axes. For the root
            WingCrossSection, this must be a zero vector. For other WingCrossSections,
            all angles must be in the range [-90, 90] degrees. Rotations are intrinsic,
            and proceed in the xy'z" order. The units are in degrees. The default is
            (0.0, 0.0, 0.0).
        :param control_surface_symmetry_type: Determines how control surfaces behave
            when the Wing has symmetry. Can be "symmetric", "asymmetric", or None. With
            "symmetric", mirrored control surfaces have the same deflection (like
            flaps). With "asymmetric", mirrored control surfaces have opposite
            deflections (like ailerons). The default is None. For Wings with type 4 or 5
            symmetry, this parameter must be specified. For Wings with type 1, 2, or 3
            symmetry, this parameter must be None. The root WingCrossSection of a type 4
            Wing cannot be "asymmetric" with a nonzero deflection, since the mirrored
            halves would tear the mesh at the shared centerline. This validation is
            performed by the parent Airplane during Wing processing.
        :param control_surface_hinge_point: The location of the control surface hinge
            from the leading edge as a fraction of chord. It must be a float in the
            range (0.0, 1.0). The default is 0.75.
        :param control_surface_deflection: The control deflection in degrees. Deflection
            downwards is positive. It must be a number (int or float) in the range
            [-5.0, 5.0] degrees. It will be converted to a float internally. The default
            is 0.0 degrees.
        :param spanwise_spacing: For non tip WingCrossSections, this can be "cosine" or
            "uniform". I highly recommend using cosine spacing. For tip
            WingCrossSections it must be None. If the parent Wing has
            explode_into_strips=True, every non tip WingCrossSection must use "uniform"
            because the explosion distributes parent-relative offsets and twists
            uniformly across each WingCrossSection's intermediates.
        :return: None
        """
        # Validate airfoil (immutable).
        if not isinstance(airfoil, airfoil_mod.Airfoil):
            raise TypeError("airfoil must be an Airfoil.")
        self._airfoil = airfoil

        # Perform a preliminary validation for num_spanwise_panels (immutable). The
        # parent Wing will later check that this is None if this WingCrossSection is a
        # tip WingCrossSection.
        if num_spanwise_panels is not None:
            num_spanwise_panels = _parameter_validation.int_in_range_return_int(
                num_spanwise_panels,
                "Non-None num_spanwise",
                min_val=1,
                min_inclusive=True,
            )
        self._num_spanwise_panels = num_spanwise_panels

        # Validate chord (immutable).
        self._chord = _parameter_validation.number_in_range_return_float(
            chord, "chord", min_val=0.0, min_inclusive=False
        )

        # Perform a preliminary validation for Lp_Wcsp_Lpp (immutable). The parent Wing
        # will later check that this is a zero vector if this WingCrossSection is a root
        # WingCrossSection.
        Lp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            Lp_Wcsp_Lpp, "Lp_Wcsp_Lpp"
        )
        Lp_Wcsp_Lpp[1] = _parameter_validation.number_in_range_return_float(
            Lp_Wcsp_Lpp[1], "Lp_Wcsp_Lpp[1]", min_val=0.0, min_inclusive=True
        )
        self._Lp_Wcsp_Lpp = Lp_Wcsp_Lpp
        self._Lp_Wcsp_Lpp.flags.writeable = False

        # Perform a preliminary validation for angles_Wcsp_to_Wcs_ixyz (immutable). The
        # parent Wing will later check that this is a zero vector if this
        # WingCrossSection is a root WingCrossSection.
        angles_Wcsp_to_Wcs_ixyz = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                angles_Wcsp_to_Wcs_ixyz, "angles_Wcsp_to_Wcs_ixyz"
            )
        )
        for angle_id, angle in enumerate(angles_Wcsp_to_Wcs_ixyz):
            angles_Wcsp_to_Wcs_ixyz[angle_id] = (
                _parameter_validation.number_in_range_return_float(
                    angle,
                    f"angles_Wcsp_to_Wcs_ixyz[{angle_id}]",
                    -90.0,
                    True,
                    90.0,
                    True,
                )
            )
        self._angles_Wcsp_to_Wcs_ixyz = angles_Wcsp_to_Wcs_ixyz
        self._angles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        # Validate control surface symmetry type (mutable, may be modified by
        # Airplane.process_wing_symmetry for type 5 symmetry).
        if control_surface_symmetry_type is not None:
            control_surface_symmetry_type = _parameter_validation.str_return_str(
                control_surface_symmetry_type, "control_surface_symmetry_type"
            )
            valid_control_surface_symmetry_types = ["symmetric", "asymmetric"]
            if (
                control_surface_symmetry_type
                not in valid_control_surface_symmetry_types
            ):
                raise ValueError(
                    f"control_surface_symmetry_type must be one of "
                    f"{valid_control_surface_symmetry_types} or None."
                )
        self.control_surface_symmetry_type = control_surface_symmetry_type

        # Validate control_surface_hinge_point and control_surface_deflection
        # (immutable).
        self._control_surface_hinge_point = (
            _parameter_validation.number_in_range_return_float(
                control_surface_hinge_point,
                "control_surface_hinge_point",
                0.0,
                False,
                1.0,
                False,
            )
        )
        self._control_surface_deflection = (
            _parameter_validation.number_in_range_return_float(
                control_surface_deflection,
                "control_surface_deflection",
                -5.0,
                True,
                5.0,
                True,
            )
        )

        # Perform a preliminary validation for spanwise_spacing (immutable). The parent
        # Wing will later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if spanwise_spacing is not None:
            spanwise_spacing = _parameter_validation.str_return_str(
                spanwise_spacing, "spanwise_spacing"
            )
            valid_non_none_spanwise_spacings = ["cosine", "uniform"]
            if spanwise_spacing not in valid_non_none_spanwise_spacings:
                raise ValueError(
                    f"Values for non None spanwise_spacing must be one of "
                    f"{valid_non_none_spanwise_spacings}."
                )
        self._spanwise_spacing = spanwise_spacing

        # Define a flag for if this WingCrossSection has been fully validated (set
        # once). This will be set by the parent Wing after calling its additional
        # validation methods.
        self._validated: bool = False

        # Define a flag for this WingCrossSection's parent Wing's symmetry type (set
        # once). This will be set by its parent Wing immediately after it has its own
        # symmetry_type parameter set by its parent Airplane.
        self._symmetry_type: int | None = None

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._T_pas_Wcsp_Lpp_to_Wcs_Lp: np.ndarray | None = None
        self._T_pas_Wcs_Lp_to_Wcsp_Lpp: np.ndarray | None = None

    # --- Deep copy method ---
    def __deepcopy__(self, memo: dict) -> WingCrossSection:
        """Creates a deep copy of this WingCrossSection.

        All attributes are copied. The Airfoil is deep copied to ensure independence.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new WingCrossSection with copied attributes.
        """
        # Create a new WingCrossSection instance without calling __init__ to avoid
        # redundant validation.
        new_wing_cross_section = object.__new__(WingCrossSection)

        # Store this WingCrossSection in memo to handle potential circular references.
        memo[id(self)] = new_wing_cross_section

        # Deep copy the Airfoil to ensure independence (immutable).
        new_wing_cross_section._airfoil = copy.deepcopy(self._airfoil, memo)

        # Copy simple immutable attributes (primitive types).
        new_wing_cross_section._num_spanwise_panels = self._num_spanwise_panels
        new_wing_cross_section._chord = self._chord
        new_wing_cross_section._control_surface_hinge_point = (
            self._control_surface_hinge_point
        )
        new_wing_cross_section._control_surface_deflection = (
            self._control_surface_deflection
        )
        new_wing_cross_section._spanwise_spacing = self._spanwise_spacing

        # Copy mutable attribute.
        new_wing_cross_section.control_surface_symmetry_type = (
            self.control_surface_symmetry_type
        )

        # Copy set once attributes directly to private fields to preserve their state.
        new_wing_cross_section._validated = self._validated
        new_wing_cross_section._symmetry_type = self._symmetry_type

        # Copy numpy arrays and make them read only.
        new_wing_cross_section._Lp_Wcsp_Lpp = self._Lp_Wcsp_Lpp.copy()
        new_wing_cross_section._Lp_Wcsp_Lpp.flags.writeable = False
        new_wing_cross_section._angles_Wcsp_to_Wcs_ixyz = (
            self._angles_Wcsp_to_Wcs_ixyz.copy()
        )
        new_wing_cross_section._angles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        # Copy cached derived properties. This preserves computation from validation.
        # For those that are numpy arrays, make the copies read only.
        new_wing_cross_section._T_pas_Wcsp_Lpp_to_Wcs_Lp = (
            self._T_pas_Wcsp_Lpp_to_Wcs_Lp.copy()
            if self._T_pas_Wcsp_Lpp_to_Wcs_Lp is not None
            else None
        )
        if new_wing_cross_section._T_pas_Wcsp_Lpp_to_Wcs_Lp is not None:
            new_wing_cross_section._T_pas_Wcsp_Lpp_to_Wcs_Lp.flags.writeable = False

        new_wing_cross_section._T_pas_Wcs_Lp_to_Wcsp_Lpp = (
            self._T_pas_Wcs_Lp_to_Wcsp_Lpp.copy()
            if self._T_pas_Wcs_Lp_to_Wcsp_Lpp is not None
            else None
        )
        if new_wing_cross_section._T_pas_Wcs_Lp_to_Wcsp_Lpp is not None:
            new_wing_cross_section._T_pas_Wcs_Lp_to_Wcsp_Lpp.flags.writeable = False

        return new_wing_cross_section

    # --- Immutable: read only properties ---
    @property
    def airfoil(self) -> airfoil_mod.Airfoil:
        return self._airfoil

    @property
    def num_spanwise_panels(self) -> int | None:
        return self._num_spanwise_panels

    @property
    def chord(self) -> float:
        return self._chord

    @property
    def Lp_Wcsp_Lpp(self) -> np.ndarray:
        return self._Lp_Wcsp_Lpp

    @property
    def angles_Wcsp_to_Wcs_ixyz(self) -> np.ndarray:
        return self._angles_Wcsp_to_Wcs_ixyz

    @property
    def control_surface_hinge_point(self) -> float:
        return self._control_surface_hinge_point

    @property
    def control_surface_deflection(self) -> float:
        return self._control_surface_deflection

    @property
    def spanwise_spacing(self) -> str | None:
        return self._spanwise_spacing

    # --- Immutable derived: manual lazy caching ---
    @property
    def T_pas_Wcsp_Lpp_to_Wcs_Lp(self) -> np.ndarray | None:
        """Defines a property for the passive transformation matrix which maps in
        homogeneous coordinates from parent wing cross section axes, relative to the
        parent leading point, to wing cross section axes, relative to the leading point.
        Is None if the WingCrossSection hasn't been fully validated yet.

        :return: A (4,4) ndarray of floats representing the transformation matrix or
            None if self.validated=False.
        """
        if not self.validated:
            return None

        if self._T_pas_Wcsp_Lpp_to_Wcs_Lp is None:
            # Step 1: Create T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp, which maps in homogeneous
            # coordinates from parent wing cross section axes relative to the parent
            # leading point to parent wing cross section axes relative to the leading
            # point. This is the translation step.
            T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp = _transformations.generate_trans_T(
                self._Lp_Wcsp_Lpp, passive=True
            )

            # Step 2: Create T_rot_pas_Wcsp_to_Wcs, which maps in homogeneous
            # coordinates from parent wing cross section axes to wing cross section
            # axes. This is the rotation step.
            T_rot_pas_Wcsp_to_Wcs = _transformations.generate_rot_T(
                self._angles_Wcsp_to_Wcs_ixyz, passive=True, intrinsic=True, order="xyz"
            )

            self._T_pas_Wcsp_Lpp_to_Wcs_Lp = _transformations.compose_T_pas(
                T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp, T_rot_pas_Wcsp_to_Wcs
            )
            self._T_pas_Wcsp_Lpp_to_Wcs_Lp.flags.writeable = False
        return self._T_pas_Wcsp_Lpp_to_Wcs_Lp

    @property
    def T_pas_Wcs_Lp_to_Wcsp_Lpp(self) -> np.ndarray | None:
        """Defines a property for the passive transformation matrix which maps in
        homogeneous coordinates from wing cross section axes, relative to the leading
        point, to parent wing cross section axes, relative to the parent leading point.
        Is None if the WingCrossSection hasn't been fully validated yet.

        :return: A (4,4) ndarray of floats representing the transformation matrix or
            None if self.validated=False.
        """
        if not self.validated:
            return None

        if self._T_pas_Wcs_Lp_to_Wcsp_Lpp is None:
            _T_pas_Wcsp_Lpp_to_Wcs_Lp = self.T_pas_Wcsp_Lpp_to_Wcs_Lp
            assert _T_pas_Wcsp_Lpp_to_Wcs_Lp is not None

            self._T_pas_Wcs_Lp_to_Wcsp_Lpp = _transformations.invert_T_pas(
                _T_pas_Wcsp_Lpp_to_Wcs_Lp
            )
            self._T_pas_Wcs_Lp_to_Wcsp_Lpp.flags.writeable = False
        return self._T_pas_Wcs_Lp_to_Wcsp_Lpp

    # --- Set once: properties with single assignment enforcement ---
    @property
    def validated(self) -> bool:
        """A flag indicating if this WingCrossSection has been fully validated by its
        parent Wing.

        :return: True if validated, False otherwise.
        """
        return self._validated

    @validated.setter
    def validated(self, value: bool) -> None:
        if self._validated:
            raise AttributeError("validated can only be set once")
        self._validated = value

    @property
    def symmetry_type(self) -> int | None:
        """The symmetry type inherited from the parent Wing.

        :return: The symmetry type (1, 2, 3, or 4) or None if not yet set.
        """
        return self._symmetry_type

    @symmetry_type.setter
    def symmetry_type(self, value: int) -> None:
        if self._symmetry_type is not None:
            raise AttributeError("symmetry_type can only be set once")
        self._symmetry_type = value

    # --- Other methods ---
    def diagram(
        self,
        *,
        show_airfoil_axes_and_points: bool | np.bool = True,
        show_airfoil: bool | np.bool = True,
        show_mcl: bool | np.bool = True,
        math_labels: bool | np.bool = False,
        save: bool | np.bool = False,
        path: str | Path = "diagram.webp",
        quality: int | float = 75.0,
    ) -> None:
        """Displays a diagram of this WingCrossSection's Airfoil, placed and scaled by
        this WingCrossSection, along with its axes and its parent axes.

        The diagram is drawn in parent wing cross section axes, relative to the parent
        leading point. It shows this WingCrossSection's axes at its leading point and
        its parent axes at its parent leading point. For a WingCrossSection whose
        symmetry type is 2 or 3, the diagram is reflected across the parent wing cross
        section axes' xz plane, so the WingCrossSection appears as it does on its
        reflected Wing. The units are in meters.

        :param show_airfoil_axes_and_points: Determines whether to draw the Airfoil's
            axes at its leading point, which coincides with this WingCrossSection's
            leading point. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is True.
        :param show_airfoil: Determines whether to draw the Airfoil's outline and mean
            camber line. Can be a bool or a numpy bool and will be converted internally
            to a bool. The default is True.
        :param show_mcl: Determines whether to draw the Airfoil's mean camber line. It
            has no effect if show_airfoil is False. Can be a bool or a numpy bool and
            will be converted internally to a bool. The default is True.
        :param math_labels: Determines whether to write the axes and point labels as
            math, set in the STIX font. Each basis direction arrow is then labeled with
            a unit vector whose superscript lists its axes' abbreviations, and each
            point with its name in capitals, whose subscript lists what it belongs to,
            if anything. If False, each label is the plain ID, set in a monospaced font.
            Can be a bool or a numpy bool and will be converted internally to a bool.
            The default is False.
        :param save: Determines whether to save the diagram as a WebP with a white
            background once its window is closed, which keeps the view's orientation and
            any labels dragged by hand. Can be a bool or a numpy bool and will be
            converted internally to a bool. The default is False.
        :param path: The file path to save the diagram to. It can be a str or a Path,
            must end with ".webp", and its directory must already exist. It has no
            effect if save is False. The default is "diagram.webp".
        :param quality: The quality of the saved WebP, where 0.0 is the smallest file
            with the most compression artifacts and 100.0 is the largest file with the
            fewest. It can be an int or a float and will be converted internally to a
            float. It has no effect if save is False. The default is 75.0.
        :return: None
        """
        show_airfoil_axes_and_points = _parameter_validation.boolLike_return_bool(
            show_airfoil_axes_and_points, "show_airfoil_axes_and_points"
        )
        show_airfoil = _parameter_validation.boolLike_return_bool(
            show_airfoil, "show_airfoil"
        )
        show_mcl = _parameter_validation.boolLike_return_bool(show_mcl, "show_mcl")
        math_labels = _parameter_validation.boolLike_return_bool(
            math_labels, "math_labels"
        )
        save = _parameter_validation.boolLike_return_bool(save, "save")
        path = _parameter_validation.pathLike_return_path(path, "path", (".webp",))
        quality = _parameter_validation.number_in_range_return_float(
            quality, "quality", 0.0, True, 100.0, True
        )

        if self.symmetry_type is None or not self.validated:
            raise ValueError(
                "A WingCrossSection can only be diagrammed after its parent Wing has "
                "validated it and set its symmetry type."
            )

        _T_pas_Wcs_Lp_to_Wcsp_Lpp = self.T_pas_Wcs_Lp_to_Wcsp_Lpp
        assert _T_pas_Wcs_Lp_to_Wcsp_Lpp is not None

        # Map from parent wing cross section axes to diagram axes, which are reflected
        # across the parent axes' xz plane for symmetry types 2 and 3. Both are relative
        # to the parent leading point.
        if self.symmetry_type in (2, 3):
            T_pas_Wcsp_Lpp_to_D_Do = _transformations.generate_reflect_T(
                np.array([0.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), passive=True
            )
        else:
            T_pas_Wcsp_Lpp_to_D_Do = np.eye(4, dtype=float)
        T_pas_Wcs_Lp_to_D_Do = _transformations.compose_T_pas(
            _T_pas_Wcs_Lp_to_Wcsp_Lpp, T_pas_Wcsp_Lpp_to_D_Do
        )

        plotter = pv.Plotter()

        if show_airfoil:
            _output_rendering.add_airfoil(
                plotter, self, T_pas_Wcs_Lp_to_D_Do, show_mcl=show_mcl
            )

        # Draw all the axes with one call, so that the airfoil axes' arrows merge with
        # the wing cross section axes' arrows they coincide with. Size the axes relative
        # to this WingCrossSection's chord, so they stay legible regardless of the
        # geometry's absolute size.
        axes_ids = ["Wcs", "Wcsp"]
        point_ids = ["Lp", "Lpp"]
        transformations = [T_pas_Wcs_Lp_to_D_Do, T_pas_Wcsp_Lpp_to_D_Do]
        if show_airfoil_axes_and_points:
            axes_ids.append("A")
            point_ids.append("Lp")
            transformations.append(
                _output_rendering.get_airfoil_axes_transformation(T_pas_Wcs_Lp_to_D_Do)
            )
        _output_rendering.add_axes_and_points(
            plotter,
            axes_ids=axes_ids,
            point_ids=point_ids,
            transformations=transformations,
            axes_scale=0.5 * self.chord,
            two_dimensional_axes_ids=["A"],
            math_labels=math_labels,
        )

        _output_rendering.show_diagram(
            plotter, cpos=(-1, -1, 1), save=save, path=path, quality=quality
        )

    # Coverage ignores this method because it is deprecated in favor of diagram.
    def get_plottable_data(  # pragma: no cover
        self,
        show: bool | np.bool = False,
    ) -> list[np.ndarray] | None:
        """A deprecated method that returns plottable data for this WingCrossSection's
        Airfoil's outline and mean camber line.

        Calling it emits a DeprecationWarning, and it will be removed in v6.0.0. Use
        diagram instead.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: None if the WingCrossSection hasn't been validated, if its symmetry
            type hasn't been set, or if show is True. Otherwise, returns a list of two
            ndarrays. These ndarrays represent points on this WingCrossSection's
            Airfoil's outline and mean camber lines, respectively. The points are in
            wing cross section axes, relative to the leading point. The units are in
            meters.
        """
        warnings.warn(
            "The get_plottable_data method is deprecated and will be removed in "
            "v6.0.0. Use diagram instead.",
            DeprecationWarning,
            stacklevel=2,
        )

        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        # If this WingCrossSection hasn't been fully validated, or its symmetry type
        # hasn't been set, return None.
        if self.symmetry_type is None or not self.validated:
            return None

        plottable_data = self.airfoil.get_plottable_data(show=False)
        assert plottable_data is not None
        [airfoilOutline_A_Lp, airfoilMcl_A_Lp] = plottable_data

        airfoilNonScaledOutline_Wcs_Lp = np.column_stack(
            [
                airfoilOutline_A_Lp[:, 0],
                np.zeros_like(airfoilOutline_A_Lp[:, 0]),
                airfoilOutline_A_Lp[:, 1],
            ]
        )
        airfoilNonScaledMcl_Wcs_Lp = np.column_stack(
            [
                airfoilMcl_A_Lp[:, 0],
                np.zeros_like(airfoilMcl_A_Lp[:, 0]),
                airfoilMcl_A_Lp[:, 1],
            ]
        )

        airfoilScalingMatrix = np.array(
            [
                [self.chord, 0.0, 0.0, 0.0],
                [0.0, self.chord, 0.0, 0.0],
                [0.0, 0.0, self.chord, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )

        airfoilOutline_Wcs_Lp = _transformations.apply_T_to_vectors(
            airfoilScalingMatrix, airfoilNonScaledOutline_Wcs_Lp, is_position=True
        )
        airfoilMcl_Wcs_Lp = _transformations.apply_T_to_vectors(
            airfoilScalingMatrix, airfoilNonScaledMcl_Wcs_Lp, is_position=True
        )

        if not show:
            return [airfoilOutline_Wcs_Lp, airfoilMcl_Wcs_Lp]

        plotter = pv.Plotter()

        _T_pas_Wcs_Lp_to_Wcsp_Lpp = self.T_pas_Wcs_Lp_to_Wcsp_Lpp
        assert _T_pas_Wcs_Lp_to_Wcsp_Lpp is not None

        airfoilOutline_Wcsp_Lpp = _transformations.apply_T_to_vectors(
            _T_pas_Wcs_Lp_to_Wcsp_Lpp, airfoilOutline_Wcs_Lp, is_position=True
        )
        airfoilMcl_Wcsp_Lpp = _transformations.apply_T_to_vectors(
            _T_pas_Wcs_Lp_to_Wcsp_Lpp, airfoilMcl_Wcs_Lp, is_position=True
        )

        if self.symmetry_type in (2, 3):
            UserMatrixAxesWcspLpp = _transformations.generate_reflect_T(
                np.array([0.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), passive=False
            )

            airfoilOutline_WcspReflectY_Lpp = _transformations.apply_T_to_vectors(
                UserMatrixAxesWcspLpp, airfoilOutline_Wcsp_Lpp, is_position=True
            )
            airfoilMcl_WcspReflectY_Lpp = _transformations.apply_T_to_vectors(
                UserMatrixAxesWcspLpp, airfoilMcl_Wcsp_Lpp, is_position=True
            )

        else:
            UserMatrixAxesWcspLpp = np.eye(4, dtype=float)
            airfoilOutline_WcspReflectY_Lpp = airfoilOutline_Wcsp_Lpp
            airfoilMcl_WcspReflectY_Lpp = airfoilMcl_Wcsp_Lpp

        rot_T_act = _transformations.generate_rot_T(
            angles=self.angles_Wcsp_to_Wcs_ixyz,
            passive=False,
            intrinsic=True,
            order="xyz",
        )
        trans_T_act = _transformations.generate_trans_T(
            translations=self.Lp_Wcsp_Lpp,
            passive=False,
        )

        UserMatrixAxesWcsLp_WcspLpp = _transformations.compose_T_act(
            rot_T_act,
            trans_T_act,
            UserMatrixAxesWcspLpp,
        )

        airfoilOutline_faces = np.hstack(
            [
                airfoilOutline_WcspReflectY_Lpp.shape[0],
                np.arange(airfoilOutline_WcspReflectY_Lpp.shape[0]),
            ]
        )
        airfoilOutline_mesh = pv.PolyData(
            airfoilOutline_WcspReflectY_Lpp, faces=airfoilOutline_faces
        )
        plotter.add_mesh(airfoilOutline_mesh)
        plotter.add_lines(airfoilMcl_WcspReflectY_Lpp)

        if np.allclose(UserMatrixAxesWcsLp_WcspLpp, UserMatrixAxesWcspLpp):
            AxesWcsLpWcspLpp_Wcsp_Lpp = pv.AxesAssembly(
                x_label="WcsX@Lp/WcspX@Lpp",
                y_label="WcsY@Lp/WcspY@Lpp",
                z_label="WcsZ@Lp/WcspZ@Lpp",
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
                user_matrix=UserMatrixAxesWcsLp_WcspLpp,
                name="Wcs/Wcsp",
                shaft_type="cylinder",
                shaft_radius=0.025,
                shaft_length=(0.8, 0.8, 0.8),
                tip_type="cone",
                tip_radius=0.1,
                tip_length=(0.2, 0.2, 0.2),
                symmetric_bounds=False,
            )
            plotter.add_actor(AxesWcsLpWcspLpp_Wcsp_Lpp)
        else:
            AxesWcsLp_Wcsp_Lpp = pv.AxesAssembly(
                x_label="WcsX@Lp",
                y_label="WcsY@Lp",
                z_label="WcsZ@Lp",
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
                user_matrix=UserMatrixAxesWcsLp_WcspLpp,
                # user_matrix=self.T_pas_Wcs_Lp_to_Wcsp_Lpp,
                name="Wcs",
                shaft_type="cylinder",
                shaft_radius=0.025,
                shaft_length=(0.8, 0.8, 0.8),
                tip_type="cone",
                tip_radius=0.1,
                tip_length=(0.2, 0.2, 0.2),
                symmetric_bounds=False,
            )
            AxesWcspLpp = pv.AxesAssembly(
                x_label="WcspX@Lpp",
                y_label="WcspY@Lpp",
                z_label="WcspZ@Lpp",
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
                user_matrix=UserMatrixAxesWcspLpp,
                name="Wcsp",
                shaft_type="cylinder",
                shaft_radius=0.025,
                shaft_length=(0.8, 0.8, 0.8),
                tip_type="cone",
                tip_radius=0.1,
                tip_length=(0.2, 0.2, 0.2),
                symmetric_bounds=False,
            )
            plotter.add_actor(AxesWcsLp_Wcsp_Lpp)
            plotter.add_actor(AxesWcspLpp)

        plotter.enable_parallel_projection()  # type: ignore[call-arg]

        plotter.show(
            cpos=(-1, -1, 1),
            full_screen=False,
            auto_close=False,
        )

        return None

    def validate_root_constraints(self) -> None:
        """Called by the parent Wing to validate constraints specific to root
        WingCrossSections.

        Root WingCrossSections must have Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_ixyz set to
        zero vectors. They also must have num_spanwise_panels not None (it's previously
        been checked to be None or a positive int).

        :return: None
        """
        # These checks are sufficient because the types were already validated by the
        # initialization method.
        if not np.allclose(self.Lp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The root WingCrossSection's Lp_Wcsp_Lpp must be np.array([0.0, 0.0, "
                "0.0])."
            )
        if not np.allclose(self.angles_Wcsp_to_Wcs_ixyz, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The root WingCrossSection's angles_Wcsp_to_Wcs_ixyz must be "
                "np.array([0.0, 0.0, 0.0])."
            )
        if self.num_spanwise_panels is None:
            raise ValueError(
                "The root WingCrossSection cannot have num_spanwise_panels set to None."
            )

    def validate_mid_constraints(self) -> None:
        """Called by the parent Wing to validate constraints specific to middle
        WingCrossSections.

        Middle WingCrossSections must have num_spanwise_panels not None (it's previously
        been checked to be None or a positive int). They must also have a strictly
        positive second component of Lp_Wcsp_Lpp.

        :return: None
        """
        if self.Lp_Wcsp_Lpp[1] <= 0.0:
            raise ValueError(
                "Non root WingCrossSections must have a strictly positive second "
                "component of Lp_Wcsp_Lpp. To create a vertical surface, define "
                "the WingCrossSections with spanwise (y) offsets and rotate the "
                "Wing using angles_Gs_to_Wn_ixyz (e.g., (90.0, 0.0, 0.0) for a "
                "vertical tail). To create a winglet, rotate the previous "
                "WingCrossSection's axes using angles_Wcsp_to_Wcs_ixyz."
            )
        if self.num_spanwise_panels is None:
            raise ValueError(
                "Middle WingCrossSections cannot have num_spanwise_panels set to None."
            )

    def validate_tip_constraints(self) -> None:
        """Called by the parent Wing to validate constraints specific to tip
        WingCrossSections.

        Tip WingCrossSections must have num_spanwise_panels and spanwise_spacing set to
        None. They must also have a strictly positive second component of Lp_Wcsp_Lpp.

        :return: None
        """
        if self.Lp_Wcsp_Lpp[1] <= 0.0:
            raise ValueError(
                "Non root WingCrossSections must have a strictly positive second "
                "component of Lp_Wcsp_Lpp. To create a vertical surface, define "
                "the WingCrossSections with spanwise (y) offsets and rotate the "
                "Wing using angles_Gs_to_Wn_ixyz (e.g., (90.0, 0.0, 0.0) for a "
                "vertical tail). To create a winglet, rotate the previous "
                "WingCrossSection's axes using angles_Wcsp_to_Wcs_ixyz."
            )
        if self.num_spanwise_panels is not None:
            raise ValueError(
                "The tip WingCrossSection must have num_spanwise_panels=None."
            )
        if self.spanwise_spacing is not None:
            raise ValueError(
                "The tip WingCrossSection must have spanwise_spacing=None."
            )
