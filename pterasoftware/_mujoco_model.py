"""Contains the class definition for interfacing with MuJoCo for free flight
simulations."""

from __future__ import annotations

from collections.abc import Sequence

import mujoco
import numpy as np

from . import movements
from . import _parameter_validation
from . import _transformations


# TEST: Add unit tests for this class's initialization.
class MuJoCoModel:
    """A class that wraps MuJoCo models and data objects to provide a clean interface
    for free flight simulations.

    Provides methods for applying aerodynamic loads to the first Airplane, advancing the
    MuJoCo simulation, and extracting the current state of the first Airplane.

    **Contains the following methods:**

    apply_loads: Applies loads to the model.

    step: Advances the MuJoCo simulation by one time step.

    get_state: Extracts the current position, orientation, velocity, and angular
    velocity from the model.

    reset: Resets the model's state to the initial conditions, time to zero seconds, and
    removes any applied loads.
    """

    def __init__(
        self,
        coupled_movement: movements.movement.CoupledMovement,
        I_BP1_CgP1: np.ndarray,
    ) -> None:
        """The initialization method.

        :param coupled_movement: The CoupledMovement this model is associated with.
        :param I_BP1_CgP1: A (3,3) ndarray of floats representing the inertia matrix of
            the airplane represented by coupled_movement's AirplaneMovement. It is in
            the first Airplane's body axes, relative to the first Airplane's CG.
        :return: None
        """
        start_key_frame_name: str = "start"

        initial_airplane = coupled_movement.airplanes[0]
        initial_coupled_operating_point = coupled_movement.coupled_operating_points[0]
        delta_time = coupled_movement.delta_time

        name = initial_airplane.name
        weight = initial_airplane.weight
        vCg__E = initial_coupled_operating_point.vCg__E
        omegasRad_BP1__E = np.deg2rad(initial_coupled_operating_point.omegas_BP1__E)
        g_E = initial_coupled_operating_point.g_E
        T_pas_E_CgP1_to_BP1_CgP1 = (
            initial_coupled_operating_point.T_pas_E_CgP1_to_BP1_CgP1
        )
        T_pas_W_CgP1_to_E_CgP1 = initial_coupled_operating_point.T_pas_W_CgP1_to_E_CgP1

        # REFACTOR: Determine if we are double-counting the gravitational force on
        #  the object by both setting gravity and mass in this class, and separately
        #  applying weight in the coupled solver.
        mass = weight / np.linalg.norm(g_E)
        vCg_W__E = np.array([vCg__E, 0.0, 0.0], dtype=float)
        omegaXRad_BP1__E, omegaYRad_BP1__E, omegaZRad_BP1__E = omegasRad_BP1__E[:]
        gX_E, gY_E, gZ_E = g_E[:]
        R_pas_E_to_BP1 = T_pas_E_CgP1_to_BP1_CgP1[:3, :3]

        IXX_BP1_CgP1, IXY_BP1_CgP1, IXZ_BP1_CgP1 = I_BP1_CgP1[0]
        IYX_BP1_CgP1, IYY_BP1_CgP1, IYZ_BP1_CgP1 = I_BP1_CgP1[1]
        IZX_BP1_CgP1, IZY_BP1_CgP1, IZZ_BP1_CgP1 = I_BP1_CgP1[2]

        IXY_BP1_CgP1 = (IXY_BP1_CgP1 + IYX_BP1_CgP1) / 2
        IXZ_BP1_CgP1 = (IXZ_BP1_CgP1 + IZX_BP1_CgP1) / 2
        IYZ_BP1_CgP1 = (IYZ_BP1_CgP1 + IZY_BP1_CgP1) / 2

        # REFACTOR: Add section on quaternions to ANGLES_VECTORS_AND_TRANSFORMATIONS.md.
        quat_E_to_BP1_wxyz = _transformations.R_to_quat_wxyz(R_pas_E_to_BP1)

        quatW_E_to_BP1, quatX_E_to_BP1, quatY_E_to_BP1, quatZ_E_to_BP1 = (
            quat_E_to_BP1_wxyz[:]
        )

        vCg_E__E = _transformations.apply_T_to_vectors(
            T_pas_W_CgP1_to_E_CgP1, vCg_W__E, has_point=False
        )

        vCgX_E__E, vCgY_E__E, vCgZ_E__E = vCg_E__E[:]

        gravity_str = f"{gX_E} {gY_E} {gZ_E}"
        inertia_str = (
            f"{IXX_BP1_CgP1} {IYY_BP1_CgP1} {IZZ_BP1_CgP1} {IXY_BP1_CgP1} "
            f"{IXZ_BP1_CgP1} {IYZ_BP1_CgP1}"
        )
        qpos_str = (
            f"0.0 0.0 0.0 {quatW_E_to_BP1} {quatX_E_to_BP1} {quatY_E_to_BP1} "
            f"{quatZ_E_to_BP1}"
        )
        qvel_str = (
            f"{vCgX_E__E} {vCgY_E__E} {vCgZ_E__E} {omegaXRad_BP1__E} "
            f"{omegaYRad_BP1__E} {omegaZRad_BP1__E}"
        )

        self.xml_str = f"""
        <mujoco model="{name}">
          <option timestep="{delta_time}" integrator="RK4" gravity="{gravity_str}"/>

          <worldbody>
            <body name="{name}" pos="0.0 0.0 0.0" >
              <freejoint/>
              <inertial pos="0.0 0.0 0.0" mass="{mass}" fullinertia="{inertia_str}"/>
            </body>
          </worldbody>

          <keyframe>
            <key name="{start_key_frame_name}" qpos="{qpos_str}" qvel="{qvel_str}"/>
          </keyframe>
        </mujoco>
        """

        # Create the internal MuJoCo model object from the XML str.
        # noinspection PyArgumentList
        self.model = mujoco.MjModel.from_xml_string(self.xml_str)

        # Set the internal model's time step to be the same as CoupledMovement's.
        self.model.opt.timestep = delta_time

        # Create the MuJoCo data structure.
        self.data: mujoco.MjData = mujoco.MjData(self.model)

        # Get and store the body ID and the initial conditions key frame ID.
        self.body_id: int = mujoco.mj_name2id(
            self.model, mujoco.mjtObj.mjOBJ_BODY, name
        )
        self.initial_key_frame_id: int = mujoco.mj_name2id(
            self.model, mujoco.mjtObj.mjOBJ_KEY, start_key_frame_name
        )

        # Set the internal model's state to the initial conditions.
        mujoco.mj_resetDataKeyframe(self.model, self.data, self.initial_key_frame_id)

        # Store initial state for reset functionality.
        self._initial_qpos: np.ndarray = np.copy(self.data.qpos)
        self._initial_qvel: np.ndarray = np.copy(self.data.qvel)

    # TEST: Add unit tests for this method.
    def apply_loads(
        self,
        forces_E: np.ndarray | Sequence[float | int],
        moments_E_CgP1: np.ndarray | Sequence[float | int],
    ) -> None:
        """Applies loads to the model.

        **Notes:**

        xfrc_applied[0:3] = forces_E: The current force applied to the first Airplane's
        CG (in Earth axes) in Newtons.

        xfrc_applied[3:6] = moments_E_CgP1: The current moment applied to the first
        Airplane's CG (in Earth axes, relative to the first Airplane's CG) in Newton
        meters.

        The loads will persist until the next call to apply_loads or until they are
        explicitly cleared.

        :param forces_E: A (3,) array-like object of numbers (int or float) representing
            the forces (in Earth axes) to apply to the first Airplane at the first
            Airplane's CG. Can be a tuple, list, or ndarray. Values are converted to
            floats internally. The units are in Newtons.
        :param moments_E_CgP1: A (3,) array-like object of numbers (int or float)
            representing the moments (in Earth axes, relative to the first Airplane's
            CG) to apply to the first Airplane at the first Airplane's CG. Can be a
            tuple, list, or ndarray. Values are converted to floats internally. The
            units are in Newton meters.
        :return: None
        """
        forces_E = _parameter_validation.threeD_number_vectorLike_return_float(
            forces_E, "forces_E"
        )
        moments_E_CgP1 = _parameter_validation.threeD_number_vectorLike_return_float(
            moments_E_CgP1, "moments_E_CgP1"
        )

        # Pack the force and moment into the model's 6-element xfrc_applied array.
        self.data.xfrc_applied[self.body_id][:] = np.hstack([forces_E, moments_E_CgP1])

    # TEST: Add unit tests for this method.
    def step(self) -> None:
        """Advances the MuJoCo simulation by one time step.

        Steps the equations of motion forward in time by one time step, taking into
        account all forces, moments, contacts, and constraints in the model.

        :return: None
        """
        mujoco.mj_step(self.model, self.data)

    # TEST: Add unit tests for this method.
    def get_state(self) -> dict[str, np.ndarray | float]:
        """Extracts the current position, orientation, velocity, and angular velocity of
        the model.

        **Notes:**

        qpos[0:3] = position_E_E: The current position of the first Airplane's CG (in
        Earth axes, relative to the Earth origin) in meters.

        qvel[0:3] = velocity_E__E: The current velocity of the first Airplane's CG (in
        Earth axes, observed from the Earth frame) in meters per second.

        np.rad2deg(qvel[3:6]) = omegas_BP1__E: The current angular velocity of the first
        Airplane's body axes (in the first Airplane's body axes, observed from the Earth
        frame) in degrees per second.

        xmat = R_pas_BP1_to_E: The current orientation of the first Airplane as a
        passive rotation matrix from the first Airplane's body axes to Earth axes.

        We define MuJoCo world coordinates to be identical to Ptera Software Earth axes.

        :return: A dictionary containing the following keys: ``position_E_E``, a (3,)
            ndarray of floats representing the current position of the first Airplane's
            CG (in Earth axes, relative to the Earth origin) in meters;
            ``R_pas_E_to_BP1``, a (3,3) ndarray of floats representing the current
            orientation of the first Airplane as a passive rotation matrix from Earth
            axes to first Airplane's body axes; ``velocity_E__E``, a (3,) ndarray of
            floats representing the current velocity of the first Airplane's CG (in
            Earth axes, observed from the Earth frame) in meters per second;
            ``omegas_BP1__E``, a (3,) ndarray of floats representing the current angular
            velocity of the first Airplane's body axes (in the first Airplane's body
            axes, observed from the Earth frame) in degrees per second; ``time``, a
            float representing the current simulation time in seconds.
        """
        # MuJoCo's xmat is R_pas_BP1_to_E: it transforms vectors from the first
        # Airplane's body axes to Earth axes. To get R_pas_E_to_BP1, we take the
        # transpose.
        R_pas_BP1_to_E = self.data.xmat[self.body_id].reshape(3, 3)
        # REFACTOR: Consider creating an invert_R_pas function in _transformations.py
        #  and calling it here.
        R_pas_E_to_BP1 = R_pas_BP1_to_E.T

        return {
            "position_E_E": np.copy(self.data.qpos[0:3]),
            "R_pas_E_to_BP1": np.copy(R_pas_E_to_BP1),
            "velocity_E__E": np.copy(self.data.qvel[0:3]),
            "omegas_BP1__E": np.rad2deg(np.copy(self.data.qvel[3:6])),
            "time": float(self.data.time),
        }

    # TEST: Add unit tests for this method.
    def reset(self) -> None:
        """Resets the model's state to the initial conditions, time to zero seconds, and
        removes any applied loads.

        :return: None
        """
        # Reset the model's state to the initial conditions.
        self.data.qpos[:] = self._initial_qpos
        self.data.qvel[:] = self._initial_qvel

        # Reset time to zero seconds.
        self.data.time = 0.0

        # Remove any applied loads.
        self.data.xfrc_applied[:] = 0.0

        # Run forward kinematics to update dependent quantities.
        mujoco.mj_forward(self.model, self.data)
