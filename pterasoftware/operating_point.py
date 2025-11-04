"""Contains the class definition for a Problem's operating point.

This module contains the following classes:
    OperatingPoint: This is a class used to contain a Problem's operating point
    characteristics.

This module contains the following functions:
    None
"""

import numpy as np

from . import _parameter_validation
from . import _transformations


class OperatingPoint:
    """This is a class used to contain a Problem's operating point characteristics.

    Citation:
        Adapted from:         performance.OperatingPoint in AeroSandbox
        Author:               Peter Sharpe
        Date of Retrieval:    04/29/2020

    This class contains the following public methods:
        qInf__E: This method calculates the freestream dynamic pressure experienced
        by the Airplane (observed in the Earth frame).

        T_pas_GP1_CgP1_to_W_CgP1: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from the first
        Airplane's geometry axes relative to the first Airplane's CG to wind axes
        relative to the first Airplane's CG.

        T_pas_W_CgP1_to_GP1_CgP1: This method defines a property for the passive
        transformation matrix which maps in homogeneous coordinates from wind axes
        relative to the first Airplane's CG to the first Airplane's geometry axes
        relative to the first Airplane's CG.

        vInfHat_GP1__E: This method computes the freestream direction (in the first
        Airplane's geometry axes, observed from the Earth frame).

        vInf_GP1__E: This method computes the freestream velocity (in the first
        Airplane's geometry axes, observed from the Earth frame).

    This class contains the following class attributes:
        None

    Subclassing:
        This class is not meant to be subclassed.
    """

    def __init__(
        self,
        rho: float | int = 1.225,
        vCg__E: float | int = 10.0,
        alpha: float | int = 5.0,
        beta: float | int = 0.0,
        externalFX_W: float | int = 0.0,
        nu: float | int = 15.06e-6,
    ) -> None:
        """This is the initialization method.

        :param rho: number, optional

            This parameter is the fluid's density. It must be a positive number and
            will be converted internally to a float. The units are kilograms per
            meters cubed. The default value is 1.225.

        :param vCg__E: number, optional

            This parameter is the speed of the Airplane's or Airplanes' CG(s)
            (observed from the Earth frame). In formation flight with multiple
            Airplanes, all Airplanes share the same velocity magnitude. Given that (1)
            this is the magnitude of a vector, and (2) we always assume a still fluid
            in our simulations, this value is equivalent to the freestream speed (the
            speed of the apparent wind, infinitely far away from the Airplane or
            Airplanes, observed while moving at the same speed as the non-accelerating
            CG or CGs). It must be a positive number and will be converted internally
            to a float. Its units are in meters per second. The default value is 10.0.

        :param alpha: number, optional

            This parameter is the angle of attack for the problem's Airplane(s). For
            more details on the exact interpretation of this value, see the
            description of wind axes in docs/AXES_POINTS_AND_FRAMES.md. It must be a
            number in the range (-180.0, 180.0] and will be converted internally to a
            float. The units are degrees. The default value is 5.0.

        :param beta: number, optional

            This parameter is the sideslip angle for the problem's Airplane(s). For
            more details on the exact interpretation of this value, see the
            description of wind axes in docs/AXES_POINTS_AND_FRAMES.md. It must be a
            number in the range (-180.0, 180.0] and will be converted internally to a
            float. The units are degrees. The default value is 0.0.

        :param externalFX_W: number, optional

            This parameter is for any additional thrust or drag on a problem's
            Airplane(s) (in wind axes) not due to the Airplanes' Wings. It is useful
            for trim analyses. It must be a number and will be converted internally
            to a float. The units are Newtons. The default value is 0.0.

        :param nu: number, optional

            This parameter is the fluid's kinematic viscosity. The units are meters
            squared per second. This parameter is only used in the unsteady ring
            vortex lattice method's vortex core growth model. It must be a positive
            number and will be converted internally to a float. Its units are in
            meters squared per second. The default value is 15.06e-6,
            which corresponds to air's kinematic viscosity at 20 degrees Celsius
            [source: https://www.engineeringtoolbox.com].
        """
        self.rho = _parameter_validation.positive_number_return_float(rho, "rho")
        # TODO: In the future, test what happens with vCg__E = 0.
        self.vCg__E = _parameter_validation.positive_number_return_float(
            vCg__E, "vCg__E"
        )
        # TODO: Restrict alpha and beta's range if testing reveals that high absolute
        #  magnitude values break things.
        self.alpha = _parameter_validation.number_in_range_return_float(
            alpha, "alpha", -180.0, False, 180.0, True
        )
        self.beta = _parameter_validation.number_in_range_return_float(
            beta, "beta", -180.0, False, 180.0, True
        )
        self.externalFX_W = _parameter_validation.number_return_float(
            externalFX_W, "externalFX_W"
        )
        self.nu = _parameter_validation.positive_number_return_float(nu, "nu")

    @property
    def qInf__E(self):
        """This method calculates the freestream dynamic pressure experienced by the
        Airplane (observed in the Earth frame).

        :return: float

            This is the freestream dynamic pressure (observed in the Earth frame).
            Its units are Pascals.
        """
        return 0.5 * self.rho * self.vCg__E**2

    @property
    def T_pas_GP1_CgP1_to_W_CgP1(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from the first Airplane's geometry axes
        relative to the first Airplane's CG to wind axes relative to the first
        Airplane's CG.

        :return: (4,4) ndarray of floats

            This is the passive transformation matrix which maps in homogeneous
            coordinates from the first Airplane's geometry axes relative to the first
            Airplane's CG to wind axes relative to the first Airplane's CG.
        """
        # Geometry axes to body axes transformation: flip x (aft to forward) and z (up
        # to down). This is equivalent to a 180-degree rotation about y.
        T_pas_GP1_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
            angles=np.array([0.0, 180.0, 0.0]),
            passive=True,
            intrinsic=False,
            order="xyz",
        )

        angles_B_to_W_exyz = np.array([0.0, -self.alpha, self.beta])

        T_pas_BP1_CgP1_to_W_CgP1 = _transformations.generate_rot_T(
            angles=angles_B_to_W_exyz, passive=True, intrinsic=False, order="xyz"
        )

        return _transformations.compose_T_pas(
            T_pas_GP1_CgP1_to_BP1_CgP1, T_pas_BP1_CgP1_to_W_CgP1
        )

    @property
    def T_pas_W_CgP1_to_GP1_CgP1(self):
        """This method defines a property for the passive transformation matrix which
        maps in homogeneous coordinates from wind axes relative to the first
        Airplane's CG to the first Airplane's geometry axes relative to the first
        Airplane's CG.

        :return: (4,4) ndarray of floats

            This is the passive transformation matrix which maps in homogeneous
            coordinates from wind axes relative to the first Airplane's CG to the
            first Airplane's geometry axes relative to the first Airplane's CG.
        """
        return _transformations.invert_T_pas(self.T_pas_GP1_CgP1_to_W_CgP1)

    @property
    def vInfHat_GP1__E(self):
        """This method computes the freestream direction (in the first Airplane's
        geometry axes, observed from the Earth frame).

        Note: See the docstring for vInf_GP1__E for details on how to interpret this
        property.

        :return: (3,) ndarray of floats

            This is a unit vector along the freestream velocity vector (in the first
            Airplane's geometry axes, observed from the Earth frame).
        """
        vInfHat_W__E = np.array([-1.0, 0.0, 0.0])

        return _transformations.apply_T_to_vectors(
            self.T_pas_W_CgP1_to_GP1_CgP1, vInfHat_W__E, has_point=False
        )

    @property
    def vInf_GP1__E(self):
        """This method computes the freestream velocity (in the first Airplane's
        geometry axes, observed from the Earth frame).

        Note: I'm defining vInf_GP1__E to be -1 * vCgX_GP1__E. This may seem obvious,
        but the important takeaways are that the freestream velocity is (1) entirely
        due to the Airplane's (or Airplanes') body's motion (a still airmass),
        and (2) the freestream velocity is observed from the Earth frame, which is
        inertial. Given point 1, a possible interpretation is that vInf_GP1__E must
        be zero, which is why I'm being specific with the definition.

        :return: (3,) ndarray of floats

            This is the freestream velocity vector (in the first Airplane's geometry
            axes, observed from the Earth frame).
        """
        return self.vInfHat_GP1__E * self.vCg__E
