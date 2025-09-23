"""This module contains the class definition for a Problem's operating point.

This module contains the following classes:
    OperatingPoint: This is a class used to contain a Problem's operating point
    characteristics.

This module contains the following exceptions:
    None

This module contains the following functions:
    None
"""

import numpy as np

from . import transformations
from . import parameter_validation


# TODO: Add unit tests for this class.
class OperatingPoint:
    """This is a class used to contain a Problem's operating point characteristics.

    Citation:
        Adapted from:         performance.OperatingPoint in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/29/2020

    This class contains the following public methods:
        qInf__E: This method calculates the freestream dynamic pressure experienced
        by the Airplane (observed in the Earth frame).

        T_pas_G_Cg_to_W_Cg: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from geometry
        axes relative to the CG to wind axes relative to the CG.

        T_pas_W_Cg_to_G_Cg: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from wind axes
        relative to the CG to geometry axes relative to the CG.

        vInfHat_G__E: This method computes the freestream direction (in geometry
        axes, observed from the Earth frame).

        vInf_G__E: This method computes the freestream velocity (in geometry axes,
        observed from the Earth frame).

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        density=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    ):
        """This is the initialization method.

        # TODO: Update the parameters' descriptions and types.
        :param density: float, optional
            This parameter is the density. The units are kilograms per meters cubed.
            The default value is 1.225.
        :param vCg__E: float, optional
            This parameter is the freestream speed in the positive x direction. The
            units are meters per second. The
            default value is 10.0.
        :param alpha: float, optional
            This parameter is the angle of attack. The units are degrees. The default
            value is 5.0.
        :param beta: float, optional
            This parameter is the sideslip angle. The units are degrees. The default
            value is 0.0.
        :param externalFX_W: float, optional
            This parameter is for any thrust that's due to the airplanes' wings. For
            example, this may hold thrust due to a non-modeled propeller or engine.
            The default value is 0.0.
        :param nu: float, optional
            This parameter is the air's kinematic viscosity. The units are meters
            squared per second. This parameter is only used in the unsteady ring
            vortex lattice method's vortex core growth model. The default value is
            15.06e-6 meters squared per second, which corresponds to air's kinematic
            viscosity at 20 degrees Celsius [source:
            https://www.engineeringtoolbox.com].
        """
        self.density = parameter_validation.positive_number_return_float(
            density, "density"
        )
        # TODO: In the future, test what happens with vCg__E = 0.
        self.vCg__E = parameter_validation.positive_number_return_float(
            vCg__E, "vCg__E"
        )
        # TODO: Restrict alpha and beta's range if testing reveals that high absolute
        #  magnitude values break things.
        self.alpha = parameter_validation.number_in_range_return_float(
            alpha, "alpha", -180.0, False, 180.0, True
        )
        self.beta = parameter_validation.number_in_range_return_float(
            beta, "beta", -180.0, False, 180.0, True
        )
        self.externalFX_W = parameter_validation.number_return_float(
            externalFX_W, "externalFX_W"
        )
        self.nu = parameter_validation.positive_number_return_float(nu, "nu")

    @property
    def qInf__E(self):
        """This method calculates the freestream dynamic pressure experienced by the
        Airplane (observed in the Earth frame).

        :return: float

            This is the freestream dynamic pressure (observed in the Earth frame).
            Its units are Pascals.
        """
        return 0.5 * self.density * self.vCg__E**2

    @property
    def T_pas_G_Cg_to_W_Cg(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from geometry axes relative to the CG to wind
        axes relative to the CG.

        :return: (4,4) ndarray of floats

            This is the passive transformation matrix which maps in homogeneous
            coordinates from geometry axes relative to the CG to wind axes relative
            to the CG.
        """
        alpha_rad = np.radians(self.alpha)
        beta_rad = np.radians(self.beta)

        T_pas_G_Cg_to_B_Cg = transformations.generate_reflect_T(
            plane_point_A_a=np.array([0.0, 0.0, 0.0]),
            plane_normal_A=np.array([1.0, 0.0, 0.0]),
            passive=True,
        )

        angles_B_to_W_exyz = np.array([0.0, alpha_rad, beta_rad])

        T_pas_B_Cg_to_W_Cg = transformations.generate_rot_T(
            angles=angles_B_to_W_exyz, passive=True, intrinsic=False, order="xyz"
        )

        return transformations.compose_T_pas(T_pas_G_Cg_to_B_Cg, T_pas_B_Cg_to_W_Cg)

    @property
    def T_pas_W_Cg_to_G_Cg(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from wind axes relative to the CG to geometry
        axes relative to the CG.

        :return: (4,4) ndarray of floats

            This is the passive transformation matrix which maps in homogeneous
            coordinates from wind axes relative to the CG to geometry axes relative
            to the CG.
        """
        return transformations.invert_T_pas(self.T_pas_G_Cg_to_W_Cg)

        # TODO: Delete the following older method of finding the rotation matrix once
        #  we've confirmed that the new method delivers equivalent results.
        # sin_alpha = np.sin(np.radians(self.alpha))
        # cos_alpha = np.cos(np.radians(self.alpha))
        # sin_beta = np.sin(np.radians(self.beta))
        # cos_beta = np.cos(np.radians(self.beta))
        # eye = np.eye(3)
        #
        # alpha_rotation = np.array(
        #     [[cos_alpha, 0, -sin_alpha], [0, 1, 0], [sin_alpha, 0, cos_alpha]]
        # )
        # beta_rotation = np.array(
        #     [[cos_beta, -sin_beta, 0], [sin_beta, cos_beta, 0], [0, 0, 1]]
        # )
        #
        # # Flip the axes because in geometry axes x is downstream by convention,
        # # while in wind axes x is upstream by convention. Same with z being up/down
        # # respectively.
        # axes_flip = np.array(
        #     [
        #         [-1, 0, 0],
        #         [0, 1, 0],
        #         [0, 0, -1],
        #     ]
        # )
        #
        # # Calculate and return the rotation matrix to convert wind axes to geometry
        # # axes.
        # rotation_matrix_wind_axes_to_geometry_axes = (
        #     axes_flip @ alpha_rotation @ beta_rotation @ eye
        # )
        # return rotation_matrix_wind_axes_to_geometry_axes

    @property
    def vInfHat_G__E(self):
        """This method computes the freestream direction (in geometry axes, observed
        from the Earth frame).

        Note: See the docstring for vInf_G__E for details on how to interpret this
        property.

        :return: (3,) ndarray of floats

            This is a unit vector along the freestream velocity vector (in geometry
            axes, observed from the Earth frame).
        """
        vInfHat_W__E = np.array([-1.0, 0.0, 0.0])

        return transformations.apply_T_to_vectors(
            self.T_pas_W_Cg_to_G_Cg, vInfHat_W__E, has_point=False
        )

    @property
    def vInf_G__E(self):
        """This method computes the freestream velocity (in geometry axes, observed
        from the Earth frame).

        Note: I'm defining vInf_G__E to be -1 * vCgX_G__E. This may seem obvious,
        but the important takeaways are that the freestream velocity is (1) entirely
        due to the body's motion (a still airmass), and (2) the freestream velocity
        is observed from the Earth frame, which is inertial. Given point 1,
        a possible interpretation is that vInf_G__E must be zero, which is why I'm
        being specific with the definition.

        :return: (3,) ndarray of floats

            This is the freestream velocity vector (in geometry axes, observed from
            the Earth frame).
        """
        freestream_velocity_geometry_axes = self.vInfHat_G__E * self.vCg__E

        return freestream_velocity_geometry_axes
