"""Contains the WingCrossSection class.

**Contains the following classes:**

WingCrossSection: A class used to contain wing cross sections of a Wing.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pyvista as pv

from . import airfoil as airfoil_mod
from .. import _parameter_validation
from .. import _transformations
from .._transformations import apply_T_to_vectors


class WingCrossSection:
    """A class used to contain the wing cross sections of a Wing.

    **Contains the following methods:**

    get_plottable_data: Returns plottable data for this WingCrossSection's Airfoil's
    outline and mean camber line.

    validate_root_constraints: Called by the parent Wing to validate constraints
    specific to root WingCrossSections.

    validate_mid_constraints: Called by the parent Wing to validate constraints specific
    to middle WingCrossSections.

    validate_tip_constraints: Called by the parent Wing to validate constraints specific
    to tip WingCrossSections.

    T_pas_Wcsp_Lpp_to_Wcs_Lp: Defines a property for the passive transformation matrix
    which maps in homogeneous coordinates from parent wing cross section axes, relative
    to the parent leading point, to wing cross section axes, relative to the leading
    point. Is None if the WingCrossSection hasn't been fully validated yet.

    T_pas_Wcs_Lp_to_Wcsp_Lpp: Defines a property for the passive transformation matrix
    which maps in homogeneous coordinates from wing cross section axes, relative to the
    leading point, to parent wing cross section axes, relative to the parent leading
    point. Is None if the WingCrossSection hasn't been fully validated yet.

    **Notes:**

    The first WingCrossSection in a Wing's wing_cross_section list is known as the root
    WingCrossSection. The last is known as the tip WingCrossSection.

    Every WingCrossSection has its own wing cross section axes. For root
    WingCrossSections, their wing cross section axes are identical in position,
    orientation, and handedness to their Wing's wing axes. For all other
    WingCrossSections, their wing cross section axes are defined relative to the axes of
    the previous WingCrossSection. Locally, the x-axis points from a cross section's
    leading point to its trailing point, the y-axis points spanwise in the general
    direction of the next WingCrossSection, and the z-axis points upwards.

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
            a zero vector. The second component must be non-negative. The units are in
            meters. The default is (0.0, 0.0, 0.0).
        :param angles_Wcsp_to_Wcs_ixyz: An array-like object of 3 numbers (int or float)
            representing the angle vector of rotation angles that define the orientation
            of this WingCrossSection's axes relative to the parent wing cross section
            axes. Can be a tuple, list, or ndarray. Values are converted to floats
            internally. If this is a root WingCrossSection, these are the wing axes. If
            not, the parent axes are the previous WingCrossSection's axes. For the root
            WingCrossSection, this must be a zero vector. For other WingCrossSections,
            all angles must be in the range [-90, 90] degrees. Rotations are intrinsic,
            and proceed in the x-y'-z" order. The units are in degrees. The default is
            (0.0, 0.0, 0.0).
        :param control_surface_symmetry_type: Determines how control surfaces behave
            when the Wing has symmetry. Can be "symmetric", "asymmetric", or None. With
            "symmetric", mirrored control surfaces have the same deflection (like
            flaps). With "asymmetric", mirrored control surfaces have opposite
            deflections (like ailerons). The default is None. For Wings with type 4 or 5
            symmetry, this parameter must be specified. For Wings with type 1, 2, or 3
            symmetry, this parameter must be None. This validation is performed by the
            parent Airplane during Wing processing.
        :param control_surface_hinge_point: The location of the control surface hinge
            from the leading edge as a fraction of chord. It must be a float in the
            range (0.0, 1.0). The default is 0.75.
        :param control_surface_deflection: The control deflection in degrees. Deflection
            downwards is positive. It must be a number (int or float) in the range
            [-5.0, 5.0] degrees. It will be converted to a float internally. The default
            is 0.0 degrees.
        :param spanwise_spacing: For non-tip WingCrossSections, this can be "cosine" or
            "uniform". I highly recommend using cosine spacing. For tip
            WingCrossSections it must be None.
        :return: None
        """
        # Validate airfoil.
        if not isinstance(airfoil, airfoil_mod.Airfoil):
            raise TypeError("airfoil must be an Airfoil.")
        self.airfoil = airfoil

        # Perform a preliminary validation for num_spanwise_panels. The parent Wing
        # will later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if num_spanwise_panels is not None:
            num_spanwise_panels = _parameter_validation.positive_int_return_int(
                num_spanwise_panels, "Non-None num_spanwise"
            )
        self.num_spanwise_panels = num_spanwise_panels

        # Validate chord.
        self.chord = _parameter_validation.positive_number_return_float(chord, "chord")

        # Perform a preliminary validation for Lp_Wcsp_Lpp. The parent Wing will
        # later check that this is a zero vector if this WingCrossSection is a root
        # WingCrossSection.
        Lp_Wcsp_Lpp = _parameter_validation.threeD_number_vectorLike_return_float(
            Lp_Wcsp_Lpp, "Lp_Wcsp_Lpp"
        )
        Lp_Wcsp_Lpp[1] = _parameter_validation.non_negative_number_return_float(
            Lp_Wcsp_Lpp[1], "Lp_Wcsp_Lpp[1]"
        )
        self.Lp_Wcsp_Lpp = Lp_Wcsp_Lpp

        # Perform a preliminary validation for angles_Wcsp_to_Wcs_ixyz. The parent
        # Wing will later check that this is a zero vector if this WingCrossSection
        # is a root WingCrossSection.
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
        self.angles_Wcsp_to_Wcs_ixyz = angles_Wcsp_to_Wcs_ixyz

        # Validate control surface symmetry type.
        if control_surface_symmetry_type is not None:
            control_surface_symmetry_type = _parameter_validation.string_return_string(
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

        # Validate control_surface_hinge_point and control_surface_deflection.
        self.control_surface_hinge_point = (
            _parameter_validation.number_in_range_return_float(
                control_surface_hinge_point,
                "control_surface_hinge_point",
                0.0,
                False,
                1.0,
                False,
            )
        )
        self.control_surface_deflection = (
            _parameter_validation.number_in_range_return_float(
                control_surface_deflection,
                "control_surface_deflection",
                -5.0,
                True,
                5.0,
                True,
            )
        )

        # Perform a preliminary validation for spanwise_spacing. The parent Wing will
        # later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if spanwise_spacing is not None:
            spanwise_spacing = _parameter_validation.string_return_string(
                spanwise_spacing, "spanwise_spacing"
            )
            valid_non_none_spanwise_spacings = ["cosine", "uniform"]
            if spanwise_spacing not in valid_non_none_spanwise_spacings:
                raise ValueError(
                    f"Values for non-None spanwise_spacing must be one of "
                    f"{valid_non_none_spanwise_spacings}."
                )
        self.spanwise_spacing = spanwise_spacing

        # Define a flag for if this WingCrossSection has been fully validated. This
        # will be set by the parent Wing after calling its additional validation
        # methods.
        self.validated: bool = False

        # Define a flag for this WingCrossSection's parent Wing's symmetry type. This
        # will be set by its parent Wing immediately after it has its own
        # symmetry_type parameter set by its parent Airplane.
        self.symmetry_type: int | None = None

    # TEST: Consider adding unit tests for this method.
    def get_plottable_data(
        self,
        show: bool | np.bool_ = False,
    ) -> list[np.ndarray] | None:
        """Returns plottable data for this WingCrossSection's Airfoil's outline and mean
        camber line.

        :param show: Determines whether to display the plot. If True, the method
            displays the plot and returns None. If False, the method returns the data
            without displaying. Can be a bool or a numpy bool and will be converted
            internally to a bool. The default is False.
        :return: If show is True, returns None. If show is False, returns a list of two
            ndarrays. These ndarrays represent points on this WingCrossSection's
            Airfoil's outline and mean camber lines, respectively. The points are in
            wing cross section axes, relative to the leading point. The units are in
            meters.
        """
        # Validate the input flag.
        show = _parameter_validation.boolLike_return_bool(show, "show")

        # If this WingCrossSection hasn't been fully validated, or its symmetry type
        # hasn't been set, return None.
        if self.symmetry_type is None or self.validated is None:
            return None

        plottable_data = self.airfoil.get_plottable_data(show=False)
        assert (
            plottable_data is not None
        ), "get_plottable_data with show=False should not return None"
        [airfoilOutline_A_lp, airfoilMcl_A_lp] = plottable_data

        airfoilNonScaledOutline_Wcs_lp = np.column_stack(
            [
                airfoilOutline_A_lp[:, 0],
                np.zeros_like(airfoilOutline_A_lp[:, 0]),
                airfoilOutline_A_lp[:, 1],
            ]
        )
        airfoilNonScaledMcl_Wcs_lp = np.column_stack(
            [
                airfoilMcl_A_lp[:, 0],
                np.zeros_like(airfoilMcl_A_lp[:, 0]),
                airfoilMcl_A_lp[:, 1],
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

        airfoilOutline_Wcs_lp = apply_T_to_vectors(
            airfoilScalingMatrix, airfoilNonScaledOutline_Wcs_lp, has_point=True
        )
        airfoilMcl_Wcs_lp = apply_T_to_vectors(
            airfoilScalingMatrix, airfoilNonScaledMcl_Wcs_lp, has_point=True
        )

        if not show:
            return [airfoilOutline_Wcs_lp, airfoilMcl_Wcs_lp]

        plotter = pv.Plotter()

        airfoilOutline_Wcsp_lpp = _transformations.apply_T_to_vectors(
            self.T_pas_Wcs_Lp_to_Wcsp_Lpp, airfoilOutline_Wcs_lp, has_point=True
        )
        airfoilMcl_Wcsp_lpp = _transformations.apply_T_to_vectors(
            self.T_pas_Wcs_Lp_to_Wcsp_Lpp, airfoilMcl_Wcs_lp, has_point=True
        )

        if self.symmetry_type in (2, 3):
            UserMatrixAxesWcspLpp = _transformations.generate_reflect_T(
                (0, 0, 0), (0, 1, 0), passive=False
            )

            airfoilOutline_WcspReflectY_lpp = _transformations.apply_T_to_vectors(
                UserMatrixAxesWcspLpp, airfoilOutline_Wcsp_lpp, has_point=True
            )
            airfoilMcl_WcspReflectY_lpp = _transformations.apply_T_to_vectors(
                UserMatrixAxesWcspLpp, airfoilMcl_Wcsp_lpp, has_point=True
            )

        else:
            UserMatrixAxesWcspLpp = np.eye(4, dtype=float)
            airfoilOutline_WcspReflectY_lpp = airfoilOutline_Wcsp_lpp
            airfoilMcl_WcspReflectY_lpp = airfoilMcl_Wcsp_lpp

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
                airfoilOutline_WcspReflectY_lpp.shape[0],
                np.arange(airfoilOutline_WcspReflectY_lpp.shape[0]),
            ]
        )
        airfoilOutline_mesh = pv.PolyData(
            airfoilOutline_WcspReflectY_lpp, faces=airfoilOutline_faces
        )
        plotter.add_mesh(airfoilOutline_mesh)
        plotter.add_lines(airfoilMcl_WcspReflectY_lpp)

        if np.allclose(UserMatrixAxesWcsLp_WcspLpp, UserMatrixAxesWcspLpp):
            AxesWcsLpWcspLpp_Wcsp_lpp = pv.AxesAssembly(
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
            plotter.add_actor(AxesWcsLpWcspLpp_Wcsp_lpp)  # type: ignore[arg-type]
        else:
            AxesWcsLp_Wcsp_lpp = pv.AxesAssembly(
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
            plotter.add_actor(AxesWcsLp_Wcsp_lpp)  # type: ignore[arg-type]
            plotter.add_actor(AxesWcspLpp)  # type: ignore[arg-type]

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

    # TODO: Check that tip WingCrossSections have self.Lp_Wcsp_Lpp[0] != 0.
    def validate_mid_constraints(self) -> None:
        """Called by the parent Wing to validate constraints specific to middle
        WingCrossSections.

        Middle WingCrossSections must have num_spanwise_panels not None (it's previously
        been checked to be None or a positive int).

        :return: None
        """
        if self.num_spanwise_panels is None:
            raise ValueError(
                "Middle WingCrossSections cannot have num_spanwise_panels set to None."
            )

    # TODO: Check that tip WingCrossSections have self.Lp_Wcsp_Lpp[0] != 0.
    def validate_tip_constraints(self) -> None:
        """Called by the parent Wing to validate constraints specific to tip
        WingCrossSections.

        Tip WingCrossSections must have num_spanwise_panels and spanwise_spacing set to
        None.

        :return: None
        """
        if self.num_spanwise_panels is not None:
            raise ValueError(
                "The tip WingCrossSection must have num_spanwise_panels=None."
            )
        if self.spanwise_spacing is not None:
            raise ValueError(
                "The tip WingCrossSection must have spanwise_spacing=None."
            )

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

        # Step 1: Create T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp, which maps in homogeneous
        # coordinates from parent wing cross section axes relative to the parent
        # leading point to parent wing cross section axes relative to the leading
        # point. This is the translation step.
        T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp = _transformations.generate_trans_T(
            self.Lp_Wcsp_Lpp, passive=True
        )

        # Step 2: Create T_rot_pas_Wcsp_to_Wcs, which maps in homogeneous coordinates
        # from parent wing cross section axes to wing cross section axes This is the
        # rotation step.
        T_rot_pas_Wcsp_to_Wcs = _transformations.generate_rot_T(
            self.angles_Wcsp_to_Wcs_ixyz, passive=True, intrinsic=True, order="xyz"
        )

        return _transformations.compose_T_pas(
            T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp, T_rot_pas_Wcsp_to_Wcs
        )

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

        return _transformations.invert_T_pas(self.T_pas_Wcsp_Lpp_to_Wcs_Lp)
