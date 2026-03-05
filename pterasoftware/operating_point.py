"""Contains the OperatingPoint class.

**Contains the following classes:**

OperatingPoint: A class used to contain the operating conditions of an aerodynamic
problem.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from . import _parameter_validation, _transformations


class OperatingPoint:
    """A class used to contain the operating conditions of an aerodynamic problem.

    **Contains the following methods:**

    qInf__E: The freestream dynamic pressure experienced by the Airplane (observed in
    the Earth frame).

    T_pas_GP1_CgP1_to_W_CgP1: The passive transformation matrix which maps in
    homogeneous coordinates from the first Airplane's geometry axes relative to the
    first Airplane's CG to wind axes relative to the first Airplane's CG.

    T_pas_W_CgP1_to_GP1_CgP1: The passive transformation matrix which maps in
    homogeneous coordinates from wind axes relative to the first Airplane's CG to the
    first Airplane's geometry axes relative to the first Airplane's CG.

    vInfHat_GP1__E: The freestream direction (in the first Airplane's geometry axes,
    observed from the Earth frame).

    vInf_GP1__E: The freestream velocity (in the first Airplane's geometry axes,
    observed from the Earth frame).

    **Citation:**

    Adapted from: performance.OperatingPoint in AeroSandbox

    Author: Peter Sharpe

    Date of retrieval: 04/29/2020
    """

    def __init__(
        self,
        rho: float | int = 1.225,
        vCg__E: float | int = 10.0,
        alpha: float | int = 5.0,
        beta: float | int = 0.0,
        angles_E_to_BP1_izyx: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        Cg_E_Eo: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        surfaceNormal_E: None | np.ndarray | Sequence[float | int] = None,
        surfacePoint_E_Eo: None | np.ndarray | Sequence[float | int] = None,
        externalFX_W: float | int = 0.0,
        nu: float | int = 15.06e-6,
    ) -> None:
        """The initialization method.

        :param rho: The fluid's density. It must be a positive number (int or float) and
            will be converted internally to a float. The units are in kilograms per
            meters cubed. The default is 1.225.
        :param vCg__E: The speed of the Airplane's or Airplanes' CG(s) (observed from
            the Earth frame). In formation flight with multiple Airplanes, all Airplanes
            share the same velocity magnitude. Given that (1) this is the magnitude of a
            vector, and (2) we always assume a still fluid in our simulations, this
            value is equivalent to the freestream speed (the speed of the apparent wind,
            infinitely far away from the Airplane or Airplanes, observed while moving at
            the same speed as the non accelerating CG or CGs). It must be a positive
            number (int or float) and will be converted internally to a float. Its units
            are in meters per second. The default is 10.0.
        :param alpha: The angle of attack for the problem's Airplane(s). For more
            details on the exact interpretation of this value, see the description of
            wind axes in docs/AXES_POINTS_AND_FRAMES.md. It must be a number (int or
            float) in the range (-180.0, 180.0] and will be converted internally to a
            float. The units are in degrees. The default is 5.0.
        :param beta: The sideslip angle for the problem's Airplane(s). For more details
            on the exact interpretation of this value, see the description of wind axes
            in docs/AXES_POINTS_AND_FRAMES.md. It must be a number (int or float) in the
            range (-180.0, 180.0] and will be converted internally to a float. The units
            are in degrees. The default is 0.0.
        :param angles_E_to_BP1_izyx: An array-like object of 3 numbers representing the
            angles from Earth axes to the first Airplane's body axes using an intrinsic
            zy'x" sequence. Can be a tuple, list, or ndarray. Values are converted to
            floats internally. Note that body axes differ from geometry axes: body axes
            point forward/right/down while geometry axes point aft/right/up. The units
            are in degrees. All angles must lie in the range (-180.0, 180.0] degrees.
            The default is (0.0, 0.0, 0.0).
        :param Cg_E_Eo: An array-like object of 3 numbers representing the position of
            the first Airplane's CG (in Earth axes, relative to the Earth origin). Can
            be a tuple, list, or ndarray. Values are converted to floats internally. The
            units are in meters. The default is (0.0, 0.0, 0.0).
        :param surfaceNormal_E: None, or an array-like of 3 numbers (int or float)
            representing the unit normal vector (in Earth axes) that, together with
            surfacePoint_E_Eo, defines the image surface used for surface effect
            modeling via the method of images. Can be None, or a tuple, list, or
            ndarray. If not None, values are converted to floats and normalized
            internally. Note that reversing the normal direction (using the antiparallel
            vector) defines the same plane and produces the same result. This value must
            be None if surfacePoint_E_Eo is None, and cannot be None if
            surfacePoint_E_Eo is not None. The default is None.
        :param surfacePoint_E_Eo: None, or an array-like of 3 numbers (int or float)
            representing a point (in Earth axes, relative to the Earth origin) that,
            along with surfaceNormal_E, defines the location of the image surface used
            for surface effect modeling via the method of images. Can be None, or a
            tuple, list, or ndarray. If not None, values are converted to floats
            internally. This value must be None if surfaceNormal_E is None, and cannot
            be None if surfaceNormal_E is not None. The units are in meters. The default
            is None.
        :param externalFX_W: The additional thrust or drag on a problem's Airplane(s)
            (in wind axes) not due to the Airplanes' Wings. It is useful for trim
            analyses. It must be a number (int or float) and will be converted
            internally to a float. The units are in Newtons. The default is 0.0.
        :param nu: The fluid's kinematic viscosity. The units are in meters squared per
            second. This parameter is only used in the unsteady ring vortex lattice
            method's vortex core growth model. It must be a positive number and will be
            converted internally to a float. Its units are in meters squared per second.
            The default is 15.06e-6, which corresponds to air's kinematic viscosity at
            20 degrees Celsius [source: https://www.engineeringtoolbox.com].
        :return: None
        """
        # Initialize the immutable attributes.
        self._rho = _parameter_validation.number_in_range_return_float(
            rho, "rho", min_val=0.0, min_inclusive=False
        )
        # TODO: In the future, test what happens with vCg__E = 0.
        self._vCg__E = _parameter_validation.number_in_range_return_float(
            vCg__E, "vCg__E", min_val=0.0, min_inclusive=False
        )
        # TODO: Restrict alpha and beta's range if testing reveals that high absolute
        #  magnitude values break things.
        self._alpha = _parameter_validation.number_in_range_return_float(
            alpha, "alpha", -180.0, False, 180.0, True
        )
        self._beta = _parameter_validation.number_in_range_return_float(
            beta, "beta", -180.0, False, 180.0, True
        )
        angles_E_to_BP1_izyx = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                angles_E_to_BP1_izyx, "angles_E_to_BP1_izyx"
            )
        )
        angles_E_to_BP1_izyx[0] = _parameter_validation.number_in_range_return_float(
            angles_E_to_BP1_izyx[0],
            "angles_E_to_BP1_izyx[0]",
            -180.0,
            False,
            180.0,
            True,
        )
        angles_E_to_BP1_izyx[1] = _parameter_validation.number_in_range_return_float(
            angles_E_to_BP1_izyx[1],
            "angles_E_to_BP1_izyx[1]",
            -180.0,
            False,
            180.0,
            True,
        )
        angles_E_to_BP1_izyx[2] = _parameter_validation.number_in_range_return_float(
            angles_E_to_BP1_izyx[2],
            "angles_E_to_BP1_izyx[2]",
            -180.0,
            False,
            180.0,
            True,
        )
        self._angles_E_to_BP1_izyx = angles_E_to_BP1_izyx
        self._angles_E_to_BP1_izyx.flags.writeable = False
        self._Cg_E_Eo = _parameter_validation.threeD_number_vectorLike_return_float(
            Cg_E_Eo, "Cg_E_Eo"
        )
        self._Cg_E_Eo.flags.writeable = False
        if surfaceNormal_E is not None and surfacePoint_E_Eo is not None:
            surfaceNormal_E = (
                _parameter_validation.threeD_number_vectorLike_return_float_unit_vector(
                    surfaceNormal_E, "surfaceNormal_E"
                )
            )
            surfaceNormal_E.flags.writeable = False
            surfacePoint_E_Eo = (
                _parameter_validation.threeD_number_vectorLike_return_float(
                    surfacePoint_E_Eo, "surfacePoint_E_Eo"
                )
            )
            surfacePoint_E_Eo.flags.writeable = False
        elif surfaceNormal_E is None and surfacePoint_E_Eo is None:
            pass
        elif surfaceNormal_E is None:
            raise ValueError(
                "surfaceNormal_E cannot be None when surfacePoint_E_Eo is not None."
            )
        else:
            raise ValueError(
                "surfacePoint_E_Eo cannot be None when surfaceNormal_E is not None."
            )
        self._surfaceNormal_E = surfaceNormal_E
        self._surfacePoint_E_Eo = surfacePoint_E_Eo
        self._externalFX_W = _parameter_validation.number_in_range_return_float(
            externalFX_W, "externalFX_W"
        )
        self._nu = _parameter_validation.number_in_range_return_float(
            nu, "nu", min_val=0.0, min_inclusive=False
        )

        # Initialize the caches for the properties derived from the immutable
        # attributes.
        self._qInf__E: float | None = None
        self._T_pas_GP1_CgP1_to_BP1_CgP1: np.ndarray | None = None
        self._T_pas_BP1_CgP1_to_GP1_CgP1: np.ndarray | None = None
        self._T_pas_BP1_CgP1_to_W_CgP1: np.ndarray | None = None
        self._T_pas_W_CgP1_to_BP1_CgP1: np.ndarray | None = None
        self._T_pas_GP1_CgP1_to_W_CgP1: np.ndarray | None = None
        self._T_pas_W_CgP1_to_GP1_CgP1: np.ndarray | None = None
        self._T_pas_E_CgP1_to_BP1_CgP1: np.ndarray | None = None
        self._T_pas_BP1_CgP1_to_E_CgP1: np.ndarray | None = None
        self._T_pas_E_CgP1_to_GP1_CgP1: np.ndarray | None = None
        self._T_pas_GP1_CgP1_to_E_CgP1: np.ndarray | None = None
        self._vInfHat_GP1__E: np.ndarray | None = None
        self._vInf_GP1__E: np.ndarray | None = None

    # --- Immutable: read only properties ---
    @property
    def rho(self) -> float:
        return self._rho

    @property
    def vCg__E(self) -> float:
        return self._vCg__E

    @property
    def alpha(self) -> float:
        return self._alpha

    @property
    def beta(self) -> float:
        return self._beta

    @property
    def angles_E_to_BP1_izyx(self) -> np.ndarray:
        return self._angles_E_to_BP1_izyx

    @property
    def Cg_E_Eo(self) -> np.ndarray:
        return self._Cg_E_Eo

    @property
    def surfaceNormal_E(self) -> np.ndarray | None:
        return self._surfaceNormal_E

    @property
    def surfacePoint_E_Eo(self) -> np.ndarray | None:
        return self._surfacePoint_E_Eo

    @property
    def externalFX_W(self) -> float:
        return self._externalFX_W

    @property
    def nu(self) -> float:
        return self._nu

    # --- Immutable derived: manual lazy caching ---
    @property
    def qInf__E(self) -> float:
        """The freestream dynamic pressure experienced by the Airplane (observed in the
        Earth frame).

        :return: The freestream dynamic pressure (observed in the Earth frame). Its
            units are in Pascals.
        """
        if self._qInf__E is None:
            self._qInf__E = 0.5 * self._rho * self._vCg__E**2
        return self._qInf__E

    @property
    def T_pas_GP1_CgP1_to_BP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's geometry axes relative to the first Airplane's CG to the
        first Airplane's body axes relative to the first Airplane's CG.

        Geometry axes to body axes transformation: flip x (aft to forward) and z (up to
        down). This is equivalent to a 180-degree rotation about y.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's geometry axes relative to the first Airplane's CG
            to the first Airplane's body axes relative to the first Airplane's CG.
        """
        if self._T_pas_GP1_CgP1_to_BP1_CgP1 is None:
            self._T_pas_GP1_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
                angles=np.array([0.0, 180.0, 0.0]),
                passive=True,
                intrinsic=False,
                order="xyz",
            )
            self._T_pas_GP1_CgP1_to_BP1_CgP1.flags.writeable = False
        return self._T_pas_GP1_CgP1_to_BP1_CgP1

    @property
    def T_pas_BP1_CgP1_to_GP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's body axes relative to the first Airplane's CG to the first
        Airplane's geometry axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's body axes relative to the first Airplane's CG to
            the first Airplane's geometry axes relative to the first Airplane's CG.
        """
        if self._T_pas_BP1_CgP1_to_GP1_CgP1 is None:
            self._T_pas_BP1_CgP1_to_GP1_CgP1 = _transformations.invert_T_pas(
                self.T_pas_GP1_CgP1_to_BP1_CgP1
            )
            self._T_pas_BP1_CgP1_to_GP1_CgP1.flags.writeable = False
        return self._T_pas_BP1_CgP1_to_GP1_CgP1

    @property
    def T_pas_BP1_CgP1_to_W_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's body axes relative to the first Airplane's CG to wind axes
        relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's body axes relative to the first Airplane's CG to
            wind axes relative to the first Airplane's CG.
        """
        if self._T_pas_BP1_CgP1_to_W_CgP1 is None:
            angles_BP1_to_W_exyz = np.array([0.0, -self._alpha, self._beta])
            self._T_pas_BP1_CgP1_to_W_CgP1 = _transformations.generate_rot_T(
                angles=angles_BP1_to_W_exyz,
                passive=True,
                intrinsic=False,
                order="xyz",
            )
            self._T_pas_BP1_CgP1_to_W_CgP1.flags.writeable = False
        return self._T_pas_BP1_CgP1_to_W_CgP1

    @property
    def T_pas_W_CgP1_to_BP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        wind axes relative to the first Airplane's CG to the first Airplane's body axes
        relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from wind axes relative to the first Airplane's CG to the first Airplane's
            body axes relative to the first Airplane's CG.
        """
        if self._T_pas_W_CgP1_to_BP1_CgP1 is None:
            self._T_pas_W_CgP1_to_BP1_CgP1 = _transformations.invert_T_pas(
                self.T_pas_BP1_CgP1_to_W_CgP1
            )
            self._T_pas_W_CgP1_to_BP1_CgP1.flags.writeable = False
        return self._T_pas_W_CgP1_to_BP1_CgP1

    @property
    def T_pas_GP1_CgP1_to_W_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's geometry axes relative to the first Airplane's CG to wind
        axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's geometry axes relative to the first Airplane's CG
            to wind axes relative to the first Airplane's CG.
        """
        if self._T_pas_GP1_CgP1_to_W_CgP1 is None:
            self._T_pas_GP1_CgP1_to_W_CgP1 = _transformations.compose_T_pas(
                self.T_pas_GP1_CgP1_to_BP1_CgP1, self.T_pas_BP1_CgP1_to_W_CgP1
            )
            self._T_pas_GP1_CgP1_to_W_CgP1.flags.writeable = False
        return self._T_pas_GP1_CgP1_to_W_CgP1

    @property
    def T_pas_W_CgP1_to_GP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        wind axes relative to the first Airplane's CG to the first Airplane's geometry
        axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from wind axes relative to the first Airplane's CG to the first Airplane's
            geometry axes relative to the first Airplane's CG.
        """
        if self._T_pas_W_CgP1_to_GP1_CgP1 is None:
            self._T_pas_W_CgP1_to_GP1_CgP1 = _transformations.invert_T_pas(
                self.T_pas_GP1_CgP1_to_W_CgP1
            )
            self._T_pas_W_CgP1_to_GP1_CgP1.flags.writeable = False
        return self._T_pas_W_CgP1_to_GP1_CgP1

    @property
    def T_pas_E_CgP1_to_BP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        Earth axes relative to the first Airplane's CG to the first Airplane's body axes
        relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from Earth axes relative to the first Airplane's CG to the first Airplane's
            body axes relative to the first Airplane's CG.
        """
        if self._T_pas_E_CgP1_to_BP1_CgP1 is None:
            self._T_pas_E_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
                angles=self._angles_E_to_BP1_izyx,
                passive=True,
                intrinsic=True,
                order="zyx",
            )
            self._T_pas_E_CgP1_to_BP1_CgP1.flags.writeable = False
        return self._T_pas_E_CgP1_to_BP1_CgP1

    @property
    def T_pas_BP1_CgP1_to_E_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's body axes relative to the first Airplane's CG to Earth axes
        relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's body axes relative to the first Airplane's CG to
            Earth axes relative to the first Airplane's CG.
        """
        if self._T_pas_BP1_CgP1_to_E_CgP1 is None:
            self._T_pas_BP1_CgP1_to_E_CgP1 = _transformations.invert_T_pas(
                self.T_pas_E_CgP1_to_BP1_CgP1
            )
            self._T_pas_BP1_CgP1_to_E_CgP1.flags.writeable = False
        return self._T_pas_BP1_CgP1_to_E_CgP1

    @property
    def T_pas_E_CgP1_to_GP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        Earth axes relative to the first Airplane's CG to the first Airplane's geometry
        axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from Earth axes relative to the first Airplane's CG to the first Airplane's
            geometry axes relative to the first Airplane's CG.
        """
        if self._T_pas_E_CgP1_to_GP1_CgP1 is None:
            self._T_pas_E_CgP1_to_GP1_CgP1 = _transformations.compose_T_pas(
                self.T_pas_E_CgP1_to_BP1_CgP1, self.T_pas_BP1_CgP1_to_GP1_CgP1
            )
            self._T_pas_E_CgP1_to_GP1_CgP1.flags.writeable = False
        return self._T_pas_E_CgP1_to_GP1_CgP1

    @property
    def T_pas_GP1_CgP1_to_E_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's geometry axes relative to the first Airplane's CG to Earth
        axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's geometry axes relative to the first Airplane's CG
            to Earth axes relative to the first Airplane's CG.
        """
        if self._T_pas_GP1_CgP1_to_E_CgP1 is None:
            self._T_pas_GP1_CgP1_to_E_CgP1 = _transformations.invert_T_pas(
                self.T_pas_E_CgP1_to_GP1_CgP1
            )
            self._T_pas_GP1_CgP1_to_E_CgP1.flags.writeable = False
        return self._T_pas_GP1_CgP1_to_E_CgP1

    @property
    def vInfHat_GP1__E(self) -> np.ndarray:
        """The freestream direction (in the first Airplane's geometry axes, observed
        from the Earth frame).

        **Notes:**

        See the docstring for vInf_GP1__E for details on how to interpret this property.

        :return: The unit vector along the freestream velocity vector (in the first
            Airplane's geometry axes, observed from the Earth frame).
        """
        if self._vInfHat_GP1__E is None:
            vInfHat_W__E = np.array([-1.0, 0.0, 0.0])

            self._vInfHat_GP1__E = _transformations.apply_T_to_vectors(
                self.T_pas_W_CgP1_to_GP1_CgP1, vInfHat_W__E, has_point=False
            )
            self._vInfHat_GP1__E.flags.writeable = False
        return self._vInfHat_GP1__E

    @property
    def vInf_GP1__E(self) -> np.ndarray:
        """The freestream velocity (in the first Airplane's geometry axes, observed from
        the Earth frame).

        **Notes:**

        I'm defining vInf_GP1__E to be -1 * vCgX_GP1__E. This may seem obvious, but the
        important takeaways are that the freestream velocity is (1) entirely due to the
        Airplane's (or Airplanes') body's motion (a still airmass), and (2) the
        freestream velocity is observed from the Earth frame, which is inertial. Given
        point 1, a possible interpretation is that vInf_GP1__E must be zero, which is
        why I'm being specific with the definition.

        :return: The freestream velocity vector (in the first Airplane's geometry axes,
            observed from the Earth frame).
        """
        if self._vInf_GP1__E is None:
            self._vInf_GP1__E = self.vInfHat_GP1__E * self._vCg__E
            self._vInf_GP1__E.flags.writeable = False
        return self._vInf_GP1__E
