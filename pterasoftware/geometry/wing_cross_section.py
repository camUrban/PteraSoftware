"""This module contains the WingCrossSection class.

This module contains the following classes:
    WingCrossSection: This is a class used to contain wing cross sections of a Wing.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from .airfoil import Airfoil

from .. import parameter_validation
from .. import transformations


class WingCrossSection:
    """This class is used to contain the wing cross sections of a Wing.

    Citation:
        Adapted from:         geometry.WingXSec in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/26/2020

    This class contains the following public methods:

        validate_root_constraints: This method is called by the parent Wing to
        validate constraints specific to root WingCrossSections.

        validate_tip_constraints: This method is called by the parent Wing to
        validate constraints specific to tip WingCrossSections.

        T_pas_Wcsp_Lpp_to_Wcs_Lp: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from parent wing
        cross section axes relative to the parent leading point to wing cross section
        axes relative to the leading point. This is set to None if the
        WingCrossSection hasn't been fully validated yet.

        T_pas_Wcs_Lp_to_Wcsp_Lpp: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from wing cross
        section axes relative to the leading point to parent wing cross section axes
        relative to the parent leading point. This is set to None if the
        WingCrossSection hasn't been fully validated yet.

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.

    The first WingCrossSection in a Wing's wing_cross_section list is known as the
    root WingCrossSection. The last is known as the tip WingCrossSection.

    Every WingCrossSection has its own wing cross section axes. For root
    WingCrossSections, their wing cross section axes are identical in position,
    orientation, and handedness to their Wing's wing axes. For all other
    WingCrossSections, their wing cross section axes are defined relative to the axes
    of the previous WingCrossSection. Locally, the x-axis points from a cross
    section's leading point to its trailing point, the y-axis points spanwise in the
    general direction of the next WingCrossSection, and the z-axis points upwards.

    Things can get a little confusing with respect to WingCrossSections for Wings
    with symmetric or mirror_only set to True. For more details, look in the Wing
    class's docstring. Also remember that WingCrossSections themselves aren't used
    for any simulations, they are merely one of the Wing attributes that tell the
    meshing function how we'd like to generate the Wing's Panels.
    """

    # ToDo: Make control_surface_type have a default value of None. Also add
    #  validation that checks that it is only not equal to None if (1) this
    #  WingCrossSection is not a tip WingCrossSection, and (2) this
    #  WingCrossSection's Wing has type 4 or type 5 symmetry.
    def __init__(
        self,
        airfoil,
        num_spanwise_panels,
        chord=1.0,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_izyx=(0.0, 0.0, 0.0),
        control_surface_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing=None,
    ):
        """This is the initialization method.

        :param airfoil: Airfoil

            This is the Airfoil to be used at this WingCrossSection.

        :param num_spanwise_panels: int or None

            This is the number of spanwise panels to be used between this
            WingCrossSection and the next one. For tip WingCrossSections,
            this must be None. For all other WingCrossSections, this must be a
            positive integer.

        :param chord: number, optional

            This is the chord of the wing at this WingCrossSection. The units are
            meters. It must be greater than 0.0 and a number (int or float). The
            default value is 1.0.

        :param Lp_Wcsp_Lpp: array-like of 3 numbers, optional

            This is the position [x, y, z] in meters of this WingCrossSection's
            leading edge in parent wing cross section axes, relative to the parent
            leading edge point. Can be a tuple, list, or numpy array of numbers (int
            or float). Values are converted to floats internally. If this is the root
            WingCrossSection, the parent wing cross section axes are the wing axes
            and the parent leading point is the Wing's leading edge root point. If
            not, the parent axes and point are those of the previous
            WingCrossSection. If this is the root WingCrossSection, it must be a zero
            vector. The second component must be non-negative. The default is (0.0,
            0.0, 0.0).

        :param angles_Wcsp_to_Wcs_izyx: array-like of 3 numbers, optional

            This is the angle vector of rotation angles [roll, pitch, yaw] in degrees
            that define the orientation of this WingCrossSection's axes relative to
            the parent wing cross section axes. Can be a tuple, list, or numpy array
            of numbers (int or float). Values are converted to floats internally. If
            this is a root WingCrossSection, these are the wing axes. If not,
            the parent axes are the previous WingCrossSection's axes. For the root
            WingCrossSection, this must be a zero vector. For other
            WingCrossSections, all angles must be in the range (-90, 90) degrees.
            Roll is rotation about the x-axis, pitch is rotation about the y-axis,
            and yaw is rotation about the z-axis. Rotations are intrinsic,
            and proceed in the z-y'-x'' order conventional for Euler angles. The
            units are degrees. The default is (0.0, 0.0, 0.0).

        :param control_surface_type: str, optional

            This is type of control surfaces for this WingCrossSection. It can
            be "symmetric" or "asymmetric". An example of symmetric control surfaces
            are flaps. An example of asymmetric control surfaces are ailerons. The
            default value is "symmetric". This value only affects a Wing's geometry
            if it has type 4 or 5 symmetry.

        :param control_surface_hinge_point: number, optional

            This is the location of the control surface hinge from the leading edge
            as a fraction of chord. It must be a number (int or float) the range (
            0.0, 1.0). The default value is 0.75.

        :param control_surface_deflection: number, optional

            This is the control deflection in degrees. Deflection downwards is
            positive. It must be a number (int or float) in the range (-90.0,
            90.0) degrees. The default value is 0.0 degrees.

        :param spanwise_spacing: str or None, optional

            For non-tip WingCrossSections, this can be "cosine" or "uniform". Using
            cosine spacing is highly recommended. For tip WingCrossSections it must
            be None.
        """
        # Validate airfoil.
        if not isinstance(airfoil, Airfoil):
            raise TypeError("airfoil must be an Airfoil.")
        self.airfoil = airfoil

        # Perform a preliminary validation for num_spanwise_panels. The parent Wing
        # will later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if num_spanwise_panels is not None:
            num_spanwise_panels = parameter_validation.validate_positive_scalar_int(
                num_spanwise_panels, "Non-None num_spanwise"
            )
        self.num_spanwise_panels = num_spanwise_panels

        # Validate chord.
        self.chord = parameter_validation.validate_positive_scalar_float(chord, "chord")

        # Perform a preliminary validation for Lp_Wcsp_Lpp. The parent Wing will
        # later check that this is a zero vector if this WingCrossSection is a root
        # WingCrossSection.
        Lp_Wcsp_Lpp = parameter_validation.validate_3d_vector_float(
            Lp_Wcsp_Lpp, "Lp_Wcsp_Lpp"
        )
        Lp_Wcsp_Lpp[1] = parameter_validation.validate_non_negative_scalar_float(
            Lp_Wcsp_Lpp[1], "Lp_Wcsp_Lpp[1]"
        )
        self.Lp_Wcsp_Lpp = Lp_Wcsp_Lpp

        # Perform a preliminary validation for angles_Wcsp_to_Wcs_izyx. The parent
        # Wing will later check that this is a zero vector if this WingCrossSection
        # is a root WingCrossSection.
        angles_Wcsp_to_Wcs_izyx = parameter_validation.validate_3d_vector_float(
            angles_Wcsp_to_Wcs_izyx, "angles_Wcsp_to_Wcs_izyx"
        )
        for angle_id, angle in enumerate(angles_Wcsp_to_Wcs_izyx):
            angles_Wcsp_to_Wcs_izyx[angle_id] = (
                parameter_validation.validate_scalar_in_range_float(
                    angle,
                    f"angles_Wcsp_to_Wcs_izyx[{angle_id}]",
                    -90.0,
                    False,
                    90.0,
                    False,
                )
            )
        self.angles_Wcsp_to_Wcs_izyx = angles_Wcsp_to_Wcs_izyx

        # Validate control surface type.
        control_surface_type = parameter_validation.validate_string(
            control_surface_type, "control_surface_type"
        )
        valid_control_surface_types = ["symmetric", "asymmetric"]
        if control_surface_type not in valid_control_surface_types:
            raise ValueError(
                f"control_surface_type must be one of {valid_control_surface_types}."
            )
        self.control_surface_type = control_surface_type

        # Validate control_surface_hinge_point and control_surface_deflection.
        self.control_surface_hinge_point = (
            parameter_validation.validate_scalar_in_range_float(
                control_surface_hinge_point,
                "control_surface_hinge_point",
                0.0,
                False,
                1.0,
                False,
            )
        )
        self.control_surface_deflection = (
            parameter_validation.validate_scalar_in_range_float(
                control_surface_deflection,
                "control_surface_deflection",
                -90.0,
                False,
                90.0,
                False,
            )
        )

        # Perform a preliminary validation for spanwise_spacing. The parent Wing will
        # later check that this is None if this WingCrossSection is a tip
        # WingCrossSection.
        if spanwise_spacing is not None:
            spanwise_spacing = parameter_validation.validate_string(
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
        self.validated = False

    def validate_root_constraints(self):
        """This method is called by the parent Wing to validate constraints specific
        to root WingCrossSections.

        Root WingCrossSections must have Lp_Wcsp_Lpp and angles_Wcsp_to_Wcs_izyx
        set to zero vectors.

        :raises ValueError: If root constraints are violated.
        """
        # These checks are sufficient because the types were already validated by the
        # initialization method.
        if not np.allclose(self.Lp_Wcsp_Lpp, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The root WingCrossSection's Lp_Wcsp_Lpp must be np.array([0.0, 0.0, 0.0])."
            )
        if not np.allclose(self.angles_Wcsp_to_Wcs_izyx, np.array([0.0, 0.0, 0.0])):
            raise ValueError(
                "The root WingCrossSection's angles_Wcsp_to_Wcs_izyx must be np.array([0.0, 0.0, 0.0])."
            )

    def validate_tip_constraints(self):
        """This method is called by the parent Wing to validate constraints specific
        to tip WingCrossSections.

        Tip WingCrossSections must have num_spanwise_panels and spanwise_spacing set
        to None.

        :raises ValueError: If tip constraints are violated.
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
    def T_pas_Wcsp_Lpp_to_Wcs_Lp(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from parent wing cross section axes relative
        to the parent leading point to wing cross section axes relative to the
        leading point. This is set to None if the WingCrossSection hasn't been
        fully validated yet.

        :return: (4,4) ndarray of floats or None
            4x4 transformation matrix or None if self.validated=False.
        """
        if not self.validated:
            return None

        # Step 1: Create T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp, which maps in homogenous
        # coordinates from parent wing cross section axes relative to the parent
        # leading point to parent wing cross section axes relative to the leading
        # point. This is the translation step.
        T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp = transformations.generate_trans_T(
            self.Lp_Wcsp_Lpp, passive=True
        )

        # Step 2: Create T_rot_pas_Wcsp_to_Wcs, which maps in homogeneous coordinates
        # from parent wing cross section axes to wing cross section axes This is the
        # rotation step.
        T_rot_pas_Wcsp_to_Wcs = transformations.generate_rot_T(
            self.angles_Wcsp_to_Wcs_izyx, passive=True, intrinsic=True, order="zyx"
        )

        return T_rot_pas_Wcsp_to_Wcs @ T_trans_pas_Wcsp_Lpp_to_Wcsp_Lp

    @property
    def T_pas_Wcs_Lp_to_Wcsp_Lpp(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from wing cross section axes relative to the
        leading point to parent wing cross section axes relative to the parent
        leading point. This is set to None if the WingCrossSection hasn't been fully
        validated yet.

        :return: (4,4) ndarray of floats or None
            4x4 transformation matrix or None if self.validated=False.
        """
        if not self.validated:
            return None

        return np.linalg.inv(self.T_pas_Wcsp_Lpp_to_Wcs_Lp)
