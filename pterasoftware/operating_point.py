"""Contains classes for the operating conditions of aerodynamic problems.

**Contains the following classes:**

OperatingPoint: A class used to contain the operating conditions of an aerodynamic
problem.

CoupledOperatingPoint: A subclass of OperatingPoint used to contain the operating
conditions at one of the time steps in a coupled aerodynamic problem.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from . import _parameter_validation
from . import _transformations


# DOCUMENT: Add the new methods to this class's docstring.
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
        self.rho = _parameter_validation.number_in_range_return_float(
            rho, "rho", min_val=0.0, min_inclusive=False
        )
        # TODO: In the future, test what happens with vCg__E = 0.
        self.vCg__E = _parameter_validation.number_in_range_return_float(
            vCg__E, "vCg__E", min_val=0.0, min_inclusive=False
        )
        # TODO: Restrict alpha and beta's range if testing reveals that high absolute
        #  magnitude values break things.
        self.alpha = _parameter_validation.number_in_range_return_float(
            alpha, "alpha", -180.0, False, 180.0, True
        )
        self.beta = _parameter_validation.number_in_range_return_float(
            beta, "beta", -180.0, False, 180.0, True
        )
        self.externalFX_W = _parameter_validation.number_in_range_return_float(
            externalFX_W, "externalFX_W"
        )
        self.nu = _parameter_validation.number_in_range_return_float(
            nu, "nu", min_val=0.0, min_inclusive=False
        )

        self.T_pas_BP1_CgP1_to_GP1_CgP1 = _transformations.generate_rot_T(
            angles=np.array([0.0, -180.0, 0.0]),
            passive=True,
            intrinsic=True,
            order="xyz",
        )

        self.T_pas_GP1_CgP1_to_BP1_CgP1 = _transformations.invert_T_pas(
            self.T_pas_BP1_CgP1_to_GP1_CgP1
        )

    @property
    def qInf__E(self) -> float:
        """The freestream dynamic pressure experienced by the Airplane (observed in the
        Earth frame).

        :return: The freestream dynamic pressure (observed in the Earth frame). Its
            units are in Pascals.
        """
        return 0.5 * self.rho * self.vCg__E**2

    @property
    def T_pas_GP1_CgP1_to_W_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        the first Airplane's geometry axes relative to the first Airplane's CG to wind
        axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from the first Airplane's geometry axes relative to the first Airplane's CG
            to wind axes relative to the first Airplane's CG.
        """

        return _transformations.compose_T_pas(
            self.T_pas_GP1_CgP1_to_BP1_CgP1, self.T_pas_BP1_CgP1_to_W_CgP1
        )

    @property
    def T_pas_W_CgP1_to_GP1_CgP1(self) -> np.ndarray:
        """The passive transformation matrix which maps in homogeneous coordinates from
        wind axes relative to the first Airplane's CG to the first Airplane's geometry
        axes relative to the first Airplane's CG.

        :return: The passive transformation matrix which maps in homogeneous coordinates
            from wind axes relative to the first Airplane's CG to the first Airplane's
            geometry axes relative to the first Airplane's CG.
        """
        return _transformations.invert_T_pas(self.T_pas_GP1_CgP1_to_W_CgP1)

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_BP1_CgP1_to_W_CgP1(self) -> np.ndarray:
        angles_BP1_to_W_exyz = np.array([0.0, -self.alpha, self.beta])

        return _transformations.generate_rot_T(
            angles=angles_BP1_to_W_exyz, passive=True, intrinsic=False, order="xyz"
        )

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_W_CgP1_to_BP1_CgP1(self) -> np.ndarray:
        return _transformations.invert_T_pas(self.T_pas_BP1_CgP1_to_W_CgP1)

    @property
    def vInfHat_GP1__E(self) -> np.ndarray:
        """The freestream direction (in the first Airplane's geometry axes, observed
        from the Earth frame).

        **Notes:**

        See the docstring for vInf_GP1__E for details on how to interpret this property.

        :return: The unit vector along the freestream velocity vector (in the first
            Airplane's geometry axes, observed from the Earth frame).
        """
        vInfHat_W__E = np.array([-1.0, 0.0, 0.0])

        return _transformations.apply_T_to_vectors(
            self.T_pas_W_CgP1_to_GP1_CgP1, vInfHat_W__E, has_point=False
        )

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
        return self.vInfHat_GP1__E * self.vCg__E


# REFACTOR: Add a section to ANGLE_VECTORS_AND_TRANSFORMATIONS.md about angular speeds.
# DOCUMENT: Add this class's new methods to its docstring.
class CoupledOperatingPoint(OperatingPoint):
    """A subclass of OperatingPoint used to contain the operating conditions at one of
    the time steps in a coupled aerodynamic problem.

    **Notes:**

    Inherits all parameters and methods from OperatingPoint without modification.

    CoupledUnsteadyProblems problems are initialized with a CoupledMovement object. This
    CoupledMovement contains an AirplaneMovement and a CoupledOperatingPoint. During
    initialization, the CoupledUnsteadyProblem will create a CoupledSteadyProblem, and
    pass it an Airplane and a CoupledOperatingPoint. This CoupledSteadyProblem will thus
    represent the Airplane and its operating conditions at the first time step. During
    the simulation, a new CoupledSteadyProblem is created for each time step, and it
    will contain the Airplane associated with that time step along with a new
    CoupledOperatingPoint with the current operating conditions.

    Although the CoupledUnsteadyRingVortexLatticeMethodSolver can only simulate a single
    airplane, CoupledOperatingPoint's documentation still include phrases like "in the
    first Airplane's geometry axes" and variable suffices like "_CgP1_GP1" for
    consistency with the rest of the codebase.

    **Contains the following methods:**

    None
    """

    def __init__(
        self,
        rho: float | int = 1.225,
        vCg__E: float | int = 10.0,
        omegas_BP1__E: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        angles_E_to_BP1_izyx: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
        alpha: float | int = 5.0,
        beta: float | int = 0.0,
        externalFX_W: float | int = 0.0,
        nu: float | int = 15.06e-6,
        g_E: np.ndarray | Sequence[float | int] = (0.0, 0.0, 9.80665),
    ) -> None:
        """The initialization method.

        See OperatingPoint's initialization method for descriptions of inherited
        parameters.

        :param omegas_BP1__E: An array-like object of numbers (int or float), with shape
            (3,), representing the current angular speed of the first Airplane's body
            axes (observed from the Earth frame). Can be a tuple, list, or ndarray.
            Values are converted to floats internally. The units are in degrees per
            second. The default is (0.0, 0.0, 0.0).
        :param angles_E_to_BP1_izyx: An array-like object of 3 numbers representing the
            current angles from Earth axes to the first Airplane's body axes using an
            intrinsic zy'x" sequence. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. Note that body axes differ from geometry
            axes: body axes point forward/right/down while geometry axes point
            aft/right/up. The units are degrees. All angles must lie in the range
            (-180.0, 180.0] degrees. The default is (0.0, 0.0, 0.0).
        :param g_E: An array-like object of numbers (ints or floats), with shape (3,),
            representing the direction of gravitational acceleration (in Earth axes).
            Can be a tuple, list, or ndarray. Values are converted to floats internally.
            The units are in meters per second per second. The default is (0.0, 0.0,
            9.80665).
        :return: None
        """
        super().__init__(rho, vCg__E, alpha, beta, externalFX_W, nu)

        self.omegas_BP1__E = (
            _parameter_validation.threeD_number_vectorLike_return_float(
                omegas_BP1__E, "omegas_BP1__E"
            )
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
        self.angles_E_to_BP1_izyx = angles_E_to_BP1_izyx

        self.g_E = _parameter_validation.threeD_number_vectorLike_return_float(
            g_E, "g_E"
        )

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_E_CgP1_to_BP1_CgP1(self) -> np.ndarray:
        return _transformations.generate_rot_T(
            angles=self.angles_E_to_BP1_izyx, passive=True, intrinsic=True, order="zyx"
        )

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_BP1_CgP1_to_E_CgP1(self) -> np.ndarray:
        return _transformations.invert_T_pas(self.T_pas_E_CgP1_to_BP1_CgP1)

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_E_CgP1_to_GP1_CgP1(self) -> np.ndarray:
        return _transformations.compose_T_pas(
            self.T_pas_E_CgP1_to_BP1_CgP1, self.T_pas_BP1_CgP1_to_GP1_CgP1
        )

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_GP1_CgP1_to_E_CgP1(self) -> np.ndarray:
        return _transformations.invert_T_pas(self.T_pas_E_CgP1_to_GP1_CgP1)

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_E_CgP1_to_W_CgP1(self) -> np.ndarray:
        return _transformations.compose_T_pas(
            self.T_pas_E_CgP1_to_BP1_CgP1, self.T_pas_BP1_CgP1_to_W_CgP1
        )

    # TEST: Add unit tests for this method.
    # DOCUMENT: Fill out this method's docstring.
    @property
    def T_pas_W_CgP1_to_E_CgP1(self) -> np.ndarray:
        return _transformations.invert_T_pas(self.T_pas_E_CgP1_to_W_CgP1)
