"""Contains the OperatingPoint class.

**Contains the following classes:**

OperatingPoint: A class used to contain the operating conditions of an aerodynamic
problem.

**Contains the following functions:**

None
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from . import _parameter_validation, _transformations


class _Unset:
    """A sentinel class used to mark that an argument was not passed to an
    OperatingPoint."""


# The module-level sentinel instance that marks an omitted argument.
_UNSET = _Unset()


class OperatingPoint:
    """A class used to contain the operating conditions of an aerodynamic problem.

    **Contains the following methods:**

    __deepcopy__: Creates a deep copy of this OperatingPoint.

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

    __slots__ = (
        "_rho",
        "_vCg__E",
        "_alpha",
        "_beta",
        "_angles_E_to_BP1_izyx",
        "_CgP1_E_Eo",
        "_surfaceNormal_E",
        "_surfacePoint_E_Eo",
        "_externalFX_W",
        "_nu",
        "_g_E",
        "_omegas_BP1__E",
        "_qInf__E",
        "_T_pas_GP1_CgP1_to_BP1_CgP1",
        "_T_pas_BP1_CgP1_to_GP1_CgP1",
        "_T_pas_BP1_CgP1_to_W_CgP1",
        "_T_pas_W_CgP1_to_BP1_CgP1",
        "_T_pas_GP1_CgP1_to_W_CgP1",
        "_T_pas_W_CgP1_to_GP1_CgP1",
        "_T_pas_E_CgP1_to_BP1_CgP1",
        "_T_pas_BP1_CgP1_to_E_CgP1",
        "_T_pas_E_CgP1_to_GP1_CgP1",
        "_T_pas_GP1_CgP1_to_E_CgP1",
        "_T_pas_W_CgP1_to_E_CgP1",
        "_T_pas_E_CgP1_to_W_CgP1",
        "_surfaceNormal_GP1",
        "_surfacePoint_GP1_CgP1",
        "_surfaceReflect_T_act_GP1_CgP1",
        "_vInfHat_GP1__E",
        "_vInf_GP1__E",
    )

    def __init__(
        self,
        rho: float | int = 1.225,
        vCg__E: float | int = 10.0,
        alpha: float | int | _Unset = _UNSET,
        beta: float | int = 0.0,
        angles_E_to_BP1_izyx: None | np.ndarray | Sequence[float | int] = None,
        CgP1_E_Eo: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        surfaceNormal_E: None | np.ndarray | Sequence[float | int] = None,
        surfacePoint_E_Eo: None | np.ndarray | Sequence[float | int] = None,
        externalFX_W: float | int | _Unset = _UNSET,
        nu: float | int = 15.06e-6,
        g_E: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        omegas_BP1__E: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
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
            float. The units are in degrees. If alpha is not passed, it defaults to 5.0
            and a FutureWarning is issued, because the default will change to None in
            v6.0.0. Pass alpha explicitly to silence the warning.
        :param beta: The sideslip angle for the problem's Airplane(s). For more details
            on the exact interpretation of this value, see the description of wind axes
            in docs/AXES_POINTS_AND_FRAMES.md. It must be a number (int or float) in the
            range (-180.0, 180.0] and will be converted internally to a float. The units
            are in degrees. The default is 0.0.
        :param angles_E_to_BP1_izyx: None, or an array-like object of 3 numbers
            representing the angles from Earth axes to the first Airplane's body axes
            using an intrinsic zy'x" sequence. Can be None, a tuple, list, or ndarray.
            If not None, values are converted to floats internally and all angles must
            lie in the range (-180.0, 180.0] degrees. Note that body axes differ from
            geometry axes: body axes point forward/right/down while geometry axes point
            aft/right/up. The units are in degrees. The default is None, which resolves
            at construction to the attitude that places the first Airplane in level
            flight along Earth +x at the given alpha and beta (equivalently, the
            attitude that makes wind axes coincide with Earth axes). This is the
            resolution most users expect: alpha and beta tilt the body relative to a
            level flight path fixed in the Earth frame rather than tilting the flight
            path relative to Earth. The distinction only affects results when an effect
            anchored to the Earth frame is active (an image surface, a gravitational
            field, or free flight). Pass (0.0, 0.0, 0.0) explicitly to instead align the
            body axes with the Earth axes regardless of alpha and beta.
        :param CgP1_E_Eo: An array-like object of 3 numbers representing the position of
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
            internally to a float. The units are in Newtons. The default is 0.0. The
            free-flight solver never applies externalFX_W and raises if it is non-zero;
            model thrust there with FreeFlightUnsteadyProblem's external_loads_fn
            instead. externalFX_W is deprecated: passing it issues a DeprecationWarning,
            and it will be removed in v6.0.0 in favor of a more general external_loads
            parameter.
        :param nu: The fluid's kinematic viscosity. The units are in meters squared per
            second. This parameter is only used in the unsteady ring vortex lattice
            method's vortex core growth model. It must be a positive number and will be
            converted internally to a float. Its units are in meters squared per second.
            The default is 15.06e-6, which corresponds to air's kinematic viscosity at
            20 degrees Celsius [source: https://www.engineeringtoolbox.com].
        :param g_E: An array-like of 3 numbers (int or float) representing the
            gravitational acceleration vector (in Earth axes). Can be a tuple, list, or
            ndarray. Values are converted to floats internally. The units are in meters
            per second squared. The default is (0.0, 0.0, 0.0), which corresponds to no
            gravitational field; set it explicitly (for example (0.0, 0.0, 9.80665) for
            standard gravity pointing along +z in Earth axes) to model a body in a
            gravitational field. This parameter is only used by the free-flight solver;
            other solvers ignore it.
        :param omegas_BP1__E: An array-like of 3 numbers (int or float) representing the
            angular velocity of the first Airplane's body axes (observed from the Earth
            frame, expressed in the first Airplane's body axes). Can be a tuple, list,
            or ndarray. Values are converted to floats internally. The units are in
            degrees per second. The default is (0.0, 0.0, 0.0). Only the free-flight
            solver accepts non-zero values; other solvers raise if any component is non-
            zero because they do not model body rotation.
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
        # Resolve the alpha sentinel. Alpha's implicit default is scheduled to change in
        # v6.0.0, so warn users who rely on it.
        if isinstance(alpha, _Unset):
            warnings.warn(
                "OperatingPoint was constructed without an explicit alpha, which "
                "currently defaults to 5.0. In v6.0.0, the default will change to "
                "None, which will resolve to 0.0 when vCg__E is positive and to NaN "
                "when vCg__E is 0.0. Pass alpha explicitly to keep the current "
                "behavior and silence this warning.",
                FutureWarning,
                stacklevel=2,
            )
            alpha = 5.0
        # TODO: Restrict alpha and beta's range if testing reveals that high absolute
        #  magnitude values break things.
        self._alpha = _parameter_validation.number_in_range_return_float(
            alpha, "alpha", -180.0, False, 180.0, True
        )
        self._beta = _parameter_validation.number_in_range_return_float(
            beta, "beta", -180.0, False, 180.0, True
        )
        if angles_E_to_BP1_izyx is None:
            # Resolve the default attitude to the one that makes wind axes coincide with
            # Earth axes, which places the first Airplane in level flight along Earth +x
            # at the given alpha and beta. The body-to-wind rotation depends only on
            # alpha and beta, so setting the Earth-to-body rotation equal to the
            # wind-to-body rotation makes Earth and wind coincide.
            T_pas_BP1_CgP1_to_W_CgP1 = _transformations.generate_rot_T(
                angles=np.array([0.0, -self._alpha, self._beta]),
                passive=True,
                intrinsic=False,
                order="xyz",
            )
            R_pas_W_to_BP1 = T_pas_BP1_CgP1_to_W_CgP1[:3, :3].T
            angles_E_to_BP1_izyx = _transformations.R_to_angles_izyx(R_pas_W_to_BP1)
        else:
            angles_E_to_BP1_izyx = (
                _parameter_validation.threeD_number_vectorLike_return_float(
                    angles_E_to_BP1_izyx, "angles_E_to_BP1_izyx"
                )
            )
            angles_E_to_BP1_izyx[0] = (
                _parameter_validation.number_in_range_return_float(
                    angles_E_to_BP1_izyx[0],
                    "angles_E_to_BP1_izyx[0]",
                    -180.0,
                    False,
                    180.0,
                    True,
                )
            )
            angles_E_to_BP1_izyx[1] = (
                _parameter_validation.number_in_range_return_float(
                    angles_E_to_BP1_izyx[1],
                    "angles_E_to_BP1_izyx[1]",
                    -180.0,
                    False,
                    180.0,
                    True,
                )
            )
            angles_E_to_BP1_izyx[2] = (
                _parameter_validation.number_in_range_return_float(
                    angles_E_to_BP1_izyx[2],
                    "angles_E_to_BP1_izyx[2]",
                    -180.0,
                    False,
                    180.0,
                    True,
                )
            )
        self._angles_E_to_BP1_izyx = angles_E_to_BP1_izyx
        self._angles_E_to_BP1_izyx.flags.writeable = False
        self._CgP1_E_Eo = _parameter_validation.threeD_number_vectorLike_return_float(
            CgP1_E_Eo, "CgP1_E_Eo"
        )
        self._CgP1_E_Eo.flags.writeable = False
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
        # Resolve the externalFX_W sentinel. The parameter is scheduled for removal in
        # v6.0.0, so warn users who pass it.
        if isinstance(externalFX_W, _Unset):
            externalFX_W = 0.0
        else:
            warnings.warn(
                "externalFX_W is deprecated and will be removed in v6.0.0, in favor "
                "of a more general external_loads parameter. Until then, stop "
                "passing externalFX_W to silence this warning.",
                DeprecationWarning,
                stacklevel=2,
            )
        self._externalFX_W = _parameter_validation.number_in_range_return_float(
            externalFX_W, "externalFX_W"
        )
        self._nu = _parameter_validation.number_in_range_return_float(
            nu, "nu", min_val=0.0, min_inclusive=False
        )
        self._g_E = _parameter_validation.threeD_number_vectorLike_return_float(
            g_E, "g_E"
        )
        self._g_E.flags.writeable = False
        self._omegas_BP1__E = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                omegas_BP1__E, "omegas_BP1__E"
            )
        )
        self._omegas_BP1__E.flags.writeable = False

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
        self._T_pas_W_CgP1_to_E_CgP1: np.ndarray | None = None
        self._T_pas_E_CgP1_to_W_CgP1: np.ndarray | None = None
        self._surfaceNormal_GP1: np.ndarray | None = None
        self._surfacePoint_GP1_CgP1: np.ndarray | None = None
        self._surfaceReflect_T_act_GP1_CgP1: np.ndarray | None = None
        self._vInfHat_GP1__E: np.ndarray | None = None
        self._vInf_GP1__E: np.ndarray | None = None

    def __deepcopy__(self, memo: dict[int, Any]) -> OperatingPoint:
        """Creates a deep copy of this OperatingPoint.

        Every attribute is preserved, since all of this OperatingPoint's attributes are
        either immutable or derived from immutable attributes and so remain valid in the
        copy. Each numpy array is copied into an independent read only array, and each
        cached derived property is copied when it has been computed or left as None when
        it has not. Each slot is handled explicitly rather than in a loop so that adding
        a slot of a different category (for example one the solver mutates, which would
        need resetting instead of copying) forces a deliberate choice here.

        :param memo: A dict used by the copy module to track already copied objects and
            avoid infinite recursion.
        :return: A new OperatingPoint with copied data.
        """
        # Create a new OperatingPoint instance without calling __init__ to avoid
        # redundant validation.
        new_operating_point = object.__new__(OperatingPoint)

        # Store this OperatingPoint in memo to handle potential circular references.
        memo[id(self)] = new_operating_point

        # Copy the immutable scalar attributes.
        new_operating_point._rho = self._rho
        new_operating_point._vCg__E = self._vCg__E
        new_operating_point._alpha = self._alpha
        new_operating_point._beta = self._beta
        new_operating_point._externalFX_W = self._externalFX_W
        new_operating_point._nu = self._nu

        # Copy the immutable array attributes into new read only arrays.
        new_operating_point._angles_E_to_BP1_izyx = self._angles_E_to_BP1_izyx.copy()
        new_operating_point._angles_E_to_BP1_izyx.flags.writeable = False
        new_operating_point._CgP1_E_Eo = self._CgP1_E_Eo.copy()
        new_operating_point._CgP1_E_Eo.flags.writeable = False
        new_operating_point._g_E = self._g_E.copy()
        new_operating_point._g_E.flags.writeable = False
        new_operating_point._omegas_BP1__E = self._omegas_BP1__E.copy()
        new_operating_point._omegas_BP1__E.flags.writeable = False

        # Copy the immutable nullable surface arrays into new read only arrays, or keep
        # them None.
        if self._surfaceNormal_E is not None:
            new_operating_point._surfaceNormal_E = self._surfaceNormal_E.copy()
            new_operating_point._surfaceNormal_E.flags.writeable = False
        else:
            new_operating_point._surfaceNormal_E = None
        if self._surfacePoint_E_Eo is not None:
            new_operating_point._surfacePoint_E_Eo = self._surfacePoint_E_Eo.copy()
            new_operating_point._surfacePoint_E_Eo.flags.writeable = False
        else:
            new_operating_point._surfacePoint_E_Eo = None

        # Copy the scalar dynamic pressure cache directly, preserving None when it has
        # not been computed.
        new_operating_point._qInf__E = self._qInf__E

        # Copy each cached transformation matrix into a new read only array, or keep it
        # None. These are all derived from immutable attributes, so the copies stay
        # valid.
        if self._T_pas_GP1_CgP1_to_BP1_CgP1 is not None:
            new_operating_point._T_pas_GP1_CgP1_to_BP1_CgP1 = (
                self._T_pas_GP1_CgP1_to_BP1_CgP1.copy()
            )
            new_operating_point._T_pas_GP1_CgP1_to_BP1_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_GP1_CgP1_to_BP1_CgP1 = None
        if self._T_pas_BP1_CgP1_to_GP1_CgP1 is not None:
            new_operating_point._T_pas_BP1_CgP1_to_GP1_CgP1 = (
                self._T_pas_BP1_CgP1_to_GP1_CgP1.copy()
            )
            new_operating_point._T_pas_BP1_CgP1_to_GP1_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_BP1_CgP1_to_GP1_CgP1 = None
        if self._T_pas_BP1_CgP1_to_W_CgP1 is not None:
            new_operating_point._T_pas_BP1_CgP1_to_W_CgP1 = (
                self._T_pas_BP1_CgP1_to_W_CgP1.copy()
            )
            new_operating_point._T_pas_BP1_CgP1_to_W_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_BP1_CgP1_to_W_CgP1 = None
        if self._T_pas_W_CgP1_to_BP1_CgP1 is not None:
            new_operating_point._T_pas_W_CgP1_to_BP1_CgP1 = (
                self._T_pas_W_CgP1_to_BP1_CgP1.copy()
            )
            new_operating_point._T_pas_W_CgP1_to_BP1_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_W_CgP1_to_BP1_CgP1 = None
        if self._T_pas_GP1_CgP1_to_W_CgP1 is not None:
            new_operating_point._T_pas_GP1_CgP1_to_W_CgP1 = (
                self._T_pas_GP1_CgP1_to_W_CgP1.copy()
            )
            new_operating_point._T_pas_GP1_CgP1_to_W_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_GP1_CgP1_to_W_CgP1 = None
        if self._T_pas_W_CgP1_to_GP1_CgP1 is not None:
            new_operating_point._T_pas_W_CgP1_to_GP1_CgP1 = (
                self._T_pas_W_CgP1_to_GP1_CgP1.copy()
            )
            new_operating_point._T_pas_W_CgP1_to_GP1_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_W_CgP1_to_GP1_CgP1 = None
        if self._T_pas_E_CgP1_to_BP1_CgP1 is not None:
            new_operating_point._T_pas_E_CgP1_to_BP1_CgP1 = (
                self._T_pas_E_CgP1_to_BP1_CgP1.copy()
            )
            new_operating_point._T_pas_E_CgP1_to_BP1_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_E_CgP1_to_BP1_CgP1 = None
        if self._T_pas_BP1_CgP1_to_E_CgP1 is not None:
            new_operating_point._T_pas_BP1_CgP1_to_E_CgP1 = (
                self._T_pas_BP1_CgP1_to_E_CgP1.copy()
            )
            new_operating_point._T_pas_BP1_CgP1_to_E_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_BP1_CgP1_to_E_CgP1 = None
        if self._T_pas_E_CgP1_to_GP1_CgP1 is not None:
            new_operating_point._T_pas_E_CgP1_to_GP1_CgP1 = (
                self._T_pas_E_CgP1_to_GP1_CgP1.copy()
            )
            new_operating_point._T_pas_E_CgP1_to_GP1_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_E_CgP1_to_GP1_CgP1 = None
        if self._T_pas_GP1_CgP1_to_E_CgP1 is not None:
            new_operating_point._T_pas_GP1_CgP1_to_E_CgP1 = (
                self._T_pas_GP1_CgP1_to_E_CgP1.copy()
            )
            new_operating_point._T_pas_GP1_CgP1_to_E_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_GP1_CgP1_to_E_CgP1 = None
        if self._T_pas_W_CgP1_to_E_CgP1 is not None:
            new_operating_point._T_pas_W_CgP1_to_E_CgP1 = (
                self._T_pas_W_CgP1_to_E_CgP1.copy()
            )
            new_operating_point._T_pas_W_CgP1_to_E_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_W_CgP1_to_E_CgP1 = None
        if self._T_pas_E_CgP1_to_W_CgP1 is not None:
            new_operating_point._T_pas_E_CgP1_to_W_CgP1 = (
                self._T_pas_E_CgP1_to_W_CgP1.copy()
            )
            new_operating_point._T_pas_E_CgP1_to_W_CgP1.flags.writeable = False
        else:
            new_operating_point._T_pas_E_CgP1_to_W_CgP1 = None

        # Copy each remaining cached derived array into a new read only array, or keep
        # it None. These are also derived from immutable attributes.
        if self._surfaceNormal_GP1 is not None:
            new_operating_point._surfaceNormal_GP1 = self._surfaceNormal_GP1.copy()
            new_operating_point._surfaceNormal_GP1.flags.writeable = False
        else:
            new_operating_point._surfaceNormal_GP1 = None
        if self._surfacePoint_GP1_CgP1 is not None:
            new_operating_point._surfacePoint_GP1_CgP1 = (
                self._surfacePoint_GP1_CgP1.copy()
            )
            new_operating_point._surfacePoint_GP1_CgP1.flags.writeable = False
        else:
            new_operating_point._surfacePoint_GP1_CgP1 = None
        if self._surfaceReflect_T_act_GP1_CgP1 is not None:
            new_operating_point._surfaceReflect_T_act_GP1_CgP1 = (
                self._surfaceReflect_T_act_GP1_CgP1.copy()
            )
            new_operating_point._surfaceReflect_T_act_GP1_CgP1.flags.writeable = False
        else:
            new_operating_point._surfaceReflect_T_act_GP1_CgP1 = None
        if self._vInfHat_GP1__E is not None:
            new_operating_point._vInfHat_GP1__E = self._vInfHat_GP1__E.copy()
            new_operating_point._vInfHat_GP1__E.flags.writeable = False
        else:
            new_operating_point._vInfHat_GP1__E = None
        if self._vInf_GP1__E is not None:
            new_operating_point._vInf_GP1__E = self._vInf_GP1__E.copy()
            new_operating_point._vInf_GP1__E.flags.writeable = False
        else:
            new_operating_point._vInf_GP1__E = None

        return new_operating_point

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
    def CgP1_E_Eo(self) -> np.ndarray:
        return self._CgP1_E_Eo

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

    @property
    def g_E(self) -> np.ndarray:
        return self._g_E

    @property
    def omegas_BP1__E(self) -> np.ndarray:
        return self._omegas_BP1__E

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
    def T_pas_W_CgP1_to_E_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        wind axes relative to the first Airplane's CG to Earth axes relative to the
        first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from wind axes relative to the first Airplane's CG to Earth axes relative to
            the first Airplane's CG.
        """
        if self._T_pas_W_CgP1_to_E_CgP1 is None:
            self._T_pas_W_CgP1_to_E_CgP1 = _transformations.compose_T_pas(
                self.T_pas_W_CgP1_to_BP1_CgP1, self.T_pas_BP1_CgP1_to_E_CgP1
            )
            self._T_pas_W_CgP1_to_E_CgP1.flags.writeable = False
        return self._T_pas_W_CgP1_to_E_CgP1

    @property
    def T_pas_E_CgP1_to_W_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        Earth axes relative to the first Airplane's CG to wind axes relative to the
        first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from Earth axes relative to the first Airplane's CG to wind axes relative to
            the first Airplane's CG.
        """
        if self._T_pas_E_CgP1_to_W_CgP1 is None:
            self._T_pas_E_CgP1_to_W_CgP1 = _transformations.invert_T_pas(
                self.T_pas_W_CgP1_to_E_CgP1
            )
            self._T_pas_E_CgP1_to_W_CgP1.flags.writeable = False
        return self._T_pas_E_CgP1_to_W_CgP1

    @property
    def surfaceNormal_GP1(self) -> np.ndarray | None:
        """The image surface's unit normal vector (in the first Airplane's geometry
        axes).

        :return: A (3,) ndarray of floats representing the image surface's unit normal
            vector (in the first Airplane's geometry axes), or None if no image surface
            is defined.
        """
        if self._surfaceNormal_E is None:
            return None
        if self._surfaceNormal_GP1 is None:
            self._surfaceNormal_GP1 = _transformations.apply_T_to_vectors(
                self.T_pas_E_CgP1_to_GP1_CgP1, self._surfaceNormal_E, is_position=False
            )
            self._surfaceNormal_GP1.flags.writeable = False
        return self._surfaceNormal_GP1

    @property
    def surfacePoint_GP1_CgP1(self) -> np.ndarray | None:
        """The position of a point on the image surface (in the first Airplane's
        geometry axes, relative to the first Airplane's CG).

        :return: A (3,) ndarray of floats representing the position of a point on the
            image surface (in the first Airplane's geometry axes, relative to the first
            Airplane's CG). The units are in meters. Returns None if no image surface is
            defined.
        """
        if self._surfacePoint_E_Eo is None:
            return None
        if self._surfacePoint_GP1_CgP1 is None:
            surfacePoint_E_CgP1 = self._surfacePoint_E_Eo - self._CgP1_E_Eo
            self._surfacePoint_GP1_CgP1 = _transformations.apply_T_to_vectors(
                self.T_pas_E_CgP1_to_GP1_CgP1, surfacePoint_E_CgP1, is_position=True
            )
            self._surfacePoint_GP1_CgP1.flags.writeable = False
        return self._surfacePoint_GP1_CgP1

    @property
    def surfaceReflect_T_act_GP1_CgP1(self) -> np.ndarray | None:
        """The active reflection transformation matrix for the image surface (in the
        first Airplane's geometry axes, relative to the first Airplane's CG).

        When applied with is_position=True, this matrix reflects a point across the
        image surface. When applied with is_position=False, it reflects a non-position
        vector (such as a velocity) across the image surface's normal direction, without
        any translational component.

        :return: A (4,4) ndarray of floats representing the active reflection
            transformation matrix (in the first Airplane's geometry axes, relative to
            the first Airplane's CG), or None if no image surface is defined.
        """
        if self._surfaceNormal_E is None or self._surfacePoint_E_Eo is None:
            return None
        if self._surfaceReflect_T_act_GP1_CgP1 is None:
            assert self.surfacePoint_GP1_CgP1 is not None
            assert self.surfaceNormal_GP1 is not None
            self._surfaceReflect_T_act_GP1_CgP1 = _transformations.generate_reflect_T(
                plane_point_A_a=self.surfacePoint_GP1_CgP1,
                plane_normal_A=self.surfaceNormal_GP1,
                passive=False,
            )
            self._surfaceReflect_T_act_GP1_CgP1.flags.writeable = False
        return self._surfaceReflect_T_act_GP1_CgP1

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
                self.T_pas_W_CgP1_to_GP1_CgP1, vInfHat_W__E, is_position=False
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
