"""Contains the MuJoCoModel class."""

from __future__ import annotations

from typing import TypedDict

import mujoco
import numpy as np

from . import _transformations


class MuJoCoState(TypedDict):
    """The state returned by MuJoCoModel.get_state."""

    position_E_Eo: np.ndarray
    R_pas_E_to_BP1: np.ndarray
    velocity_E__E: np.ndarray
    omegas_BP1__E: np.ndarray
    time: float


class MuJoCoModel:
    """A class used to interface with MuJoCo for free flight simulations.

    **Contains the following methods:**

    apply_loads: Applies loads to the model.

    step: Advances the MuJoCo simulation by one time step.

    get_state: Extracts the current position, orientation, velocity, and angular
    velocity from the model.

    reset: Resets the model's state to the initial conditions, time to zero seconds, and
    removes any applied loads.

    save_state: Saves and returns a snapshot of the model's current integration state.

    restore_state: Restores the model to a previously saved integration state.

    **Notes:**

    Wraps MuJoCo models and data objects to provide a clean interface for applying
    aerodynamic loads to the first Airplane, advancing the MuJoCo simulation, and
    extracting the current state of the first Airplane.

    MuJoCoModel performs no input validation. It is private and constructed only by
    FreeFlightUnsteadyProblem, which supplies already-validated inputs; see that class
    for the validation of its constructor arguments.

    The model integrates the rigid-body dynamics freely and has no awareness of any
    aerodynamic image surface (the method-of-images ground or water plane the
    aerodynamics reflect across). Any surface effect reaches the body only through the
    aerodynamic loads passed to apply_loads; the surface is not modeled as a collision
    boundary, so nothing in the dynamics prevents the body from passing through it.
    """

    __slots__ = (
        "_xml_str",
        "_model",
        "_data",
        "_body_id",
        "_initial_key_frame_id",
        "_initial_qpos",
        "_initial_qvel",
        "_mujoco_assets",
    )

    def __init__(
        self,
        name: str,
        mass: float | int,
        omegas_BP1__E: np.ndarray,
        T_pas_BP1_CgP1_to_E_CgP1: np.ndarray,
        vCg_E__E: np.ndarray,
        I_BP1_CgP1: np.ndarray,
        delta_time: float | int,
        integrator: str = "RK4",
        extra_xml: dict[str, str] | None = None,
        mujoco_assets: dict[str, bytes] | None = None,
    ) -> None:
        """The initialization method.

        :param name: The name of the Airplane, used as the MuJoCo body name. Supplied by
            FreeFlightUnsteadyProblem from the Airplane.
        :param mass: The mass of the first Airplane in kilograms. Supplied by
            FreeFlightUnsteadyProblem.
        :param omegas_BP1__E: A (3,) ndarray of floats representing the initial angular
            velocity of the first Airplane's body axes (in the first Airplane's body
            axes, observed from the Earth frame), in degrees per second. Supplied by
            FreeFlightUnsteadyProblem from the OperatingPoint.
        :param T_pas_BP1_CgP1_to_E_CgP1: A (4,4) ndarray of floats representing the
            passive transformation matrix from the first Airplane's body axes, relative
            to the first Airplane's CG, to Earth axes, relative to the first Airplane's
            CG. Supplied by FreeFlightUnsteadyProblem from the OperatingPoint.
        :param vCg_E__E: A (3,) ndarray of floats representing the initial velocity of
            the first Airplane's CG (in Earth axes, observed from the Earth frame), in
            meters per second. Computed by FreeFlightUnsteadyProblem from the
            OperatingPoint.
        :param I_BP1_CgP1: A (3,3) ndarray of floats representing the inertia matrix of
            the first Airplane (in the first Airplane's body axes, relative to the first
            Airplane's CG), in kilogram square meters. Supplied by
            FreeFlightUnsteadyProblem, which validates that it is symmetric.
        :param delta_time: The time, in seconds, between each time step. Supplied by
            FreeFlightUnsteadyProblem from the Movement.
        :param integrator: A str naming the MuJoCo integrator to set in the generated
            model XML's option element. Validated by FreeFlightUnsteadyProblem before
            being passed here. The default is "RK4".
        :param extra_xml: A dict (or None) mapping injection point names to XML fragment
            strings to inject into the generated MuJoCo XML. Supported keys are
            "default", "asset", and "visual" (inserted as top level elements),
            "worldbody" (inserted inside the worldbody element, before the body), and
            "body" (inserted inside the body element, after the inertial element).
            Validated by FreeFlightUnsteadyProblem before being passed here. The default
            is None, which injects no extra XML.
        :param mujoco_assets: A dict (or None) mapping virtual filenames to their binary
            contents. These are passed to MuJoCo's from_xml_string as the assets
            parameter, allowing meshes and other binary files to be loaded without
            writing to disk. The dict is retained, which lets the serialization layer
            refuse to save an asset-based model. Validated by FreeFlightUnsteadyProblem
            before being passed here. The default is None, which provides no extra
            assets.
        :return: None
        """
        start_key_frame_name: str = "start"

        omegasRad_BP1__E = np.deg2rad(omegas_BP1__E)

        omegaXRad_BP1__E, omegaYRad_BP1__E, omegaZRad_BP1__E = omegasRad_BP1__E[:]

        R_pas_BP1_to_E = T_pas_BP1_CgP1_to_E_CgP1[:3, :3]

        R_act_BP1_to_E = np.linalg.inv(R_pas_BP1_to_E)

        R_act_E_to_BP1 = R_act_BP1_to_E.T

        quat_act_E_to_BP1_wxyz = _transformations.R_to_quat_wxyz(R_act_E_to_BP1)

        IXX_BP1_CgP1, IXY_BP1_CgP1, IXZ_BP1_CgP1 = I_BP1_CgP1[0]
        IYX_BP1_CgP1, IYY_BP1_CgP1, IYZ_BP1_CgP1 = I_BP1_CgP1[1]
        IZX_BP1_CgP1, IZY_BP1_CgP1, IZZ_BP1_CgP1 = I_BP1_CgP1[2]

        IXY_BP1_CgP1 = (IXY_BP1_CgP1 + IYX_BP1_CgP1) / 2
        IXZ_BP1_CgP1 = (IXZ_BP1_CgP1 + IZX_BP1_CgP1) / 2
        IYZ_BP1_CgP1 = (IYZ_BP1_CgP1 + IZY_BP1_CgP1) / 2

        (
            quatW_act_E_to_BP1,
            quatX_act_E_to_BP1,
            quatY_act_E_to_BP1,
            quatZ_act_E_to_BP1,
        ) = quat_act_E_to_BP1_wxyz[:]

        vCgX_E__E, vCgY_E__E, vCgZ_E__E = vCg_E__E[:]

        # Gravity in the MuJoCo model is turned off as it is applied by
        # FreeFlightUnsteadyProblem.initialize_next_problem.
        gravity_str = "0.0 0.0 0.0"
        inertia_str = (
            f"{IXX_BP1_CgP1} {IYY_BP1_CgP1} {IZZ_BP1_CgP1} {IXY_BP1_CgP1} "
            f"{IXZ_BP1_CgP1} {IYZ_BP1_CgP1}"
        )
        qpos_str = (
            f"0.0 0.0 0.0 {quatW_act_E_to_BP1} {quatX_act_E_to_BP1} {quatY_act_E_to_BP1} "
            f"{quatZ_act_E_to_BP1}"
        )
        qvel_str = (
            f"{vCgX_E__E} {vCgY_E__E} {vCgZ_E__E} {omegaXRad_BP1__E} "
            f"{omegaYRad_BP1__E} {omegaZRad_BP1__E}"
        )

        # Build the extra XML fragments to inject. If extra_xml is None, use an empty
        # dict so the .get calls below return empty strings.
        _extra = extra_xml if extra_xml is not None else {}
        extra_default = _extra.get("default", "")
        extra_asset = _extra.get("asset", "")
        extra_visual = _extra.get("visual", "")
        extra_worldbody = _extra.get("worldbody", "")
        extra_body = _extra.get("body", "")

        # Initialize the immutable attributes.
        self._xml_str: str = f"""
        <mujoco model="{name}">
          <option timestep="{delta_time}" integrator="{integrator}" gravity="{gravity_str}"/>

          {extra_default}
          {extra_asset}
          {extra_visual}

          <worldbody>
            {extra_worldbody}
            <body name="{name}" pos="0.0 0.0 0.0" >
              <freejoint/>
              <inertial pos="0.0 0.0 0.0" mass="{mass}" fullinertia="{inertia_str}"/>
              {extra_body}
            </body>
          </worldbody>

          <keyframe>
            <key name="{start_key_frame_name}" qpos="{qpos_str}" qvel="{qvel_str}"/>
          </keyframe>
        </mujoco>
        """

        # Create the internal MuJoCo model object from the XML str. If mujoco_assets
        # is provided, pass it so MuJoCo can resolve virtual filenames (e.g., STL
        # meshes) without writing them to disk.
        # noinspection PyArgumentList
        if mujoco_assets is not None:
            self._model: mujoco.MjModel = mujoco.MjModel.from_xml_string(
                self._xml_str, assets=mujoco_assets
            )
        else:
            # noinspection PyArgumentList
            self._model = mujoco.MjModel.from_xml_string(self._xml_str)

        # Retain the assets dict so the serialization layer can refuse to save an
        # asset-based model: the engine is rebuilt on load from the stored XML alone,
        # whose asset references would be unresolvable.
        self._mujoco_assets: dict[str, bytes] | None = mujoco_assets

        # Set the internal model's time step to be the same as the simulation's.
        self._model.opt.timestep = delta_time

        # Initialize the mutable attributes.
        self._data: mujoco.MjData = mujoco.MjData(self._model)

        # Get and store the body ID and the initial conditions key frame ID.
        self._body_id: int = mujoco.mj_name2id(
            self._model, mujoco.mjtObj.mjOBJ_BODY, name
        )
        self._initial_key_frame_id: int = mujoco.mj_name2id(
            self._model, mujoco.mjtObj.mjOBJ_KEY, start_key_frame_name
        )

        # Set the internal model's state to the initial conditions.
        mujoco.mj_resetDataKeyframe(self._model, self._data, self._initial_key_frame_id)

        # Run forward kinematics to compute derived quantities (xmat, xpos, etc.)
        # from the initial qpos/qvel. Without this, xmat would be zeros until the
        # first call to mj_step.
        mujoco.mj_forward(self._model, self._data)

        # Store initial state for reset functionality.
        self._initial_qpos: np.ndarray = np.copy(self._data.qpos)
        self._initial_qpos.flags.writeable = False
        self._initial_qvel: np.ndarray = np.copy(self._data.qvel)
        self._initial_qvel.flags.writeable = False

    # --- Immutable: read only properties ---
    @property
    def xml_str(self) -> str:
        return self._xml_str

    @property
    def body_id(self) -> int:
        return self._body_id

    @property
    def initial_key_frame_id(self) -> int:
        return self._initial_key_frame_id

    @property
    def initial_qpos(self) -> np.ndarray:
        return self._initial_qpos

    @property
    def initial_qvel(self) -> np.ndarray:
        return self._initial_qvel

    # --- Other methods ---
    def apply_loads(
        self,
        forces_E: np.ndarray,
        moments_E_CgP1: np.ndarray,
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

        :param forces_E: A (3,) ndarray of floats representing the forces (in Earth
            axes) to apply to the first Airplane at the first Airplane's CG, in Newtons.
            Supplied by FreeFlightUnsteadyProblem.initialize_next_problem.
        :param moments_E_CgP1: A (3,) ndarray of floats representing the moments (in
            Earth axes, relative to the first Airplane's CG) to apply to the first
            Airplane at the first Airplane's CG, in Newton meters. Supplied by
            FreeFlightUnsteadyProblem.initialize_next_problem.
        :return: None
        """
        # Pack the force and moment into the model's 6-element xfrc_applied array.
        self._data.xfrc_applied[self._body_id][:] = np.hstack(
            [forces_E, moments_E_CgP1]
        )

    def step(self) -> None:
        """Advances the MuJoCo simulation by one time step.

        Steps the equations of motion forward in time by one time step, taking into
        account all forces, moments, contacts, and constraints in the model.

        :return: None
        """
        mujoco.mj_step(self._model, self._data)

    def get_state(self) -> MuJoCoState:
        """Extracts the current position, orientation, velocity, and angular velocity of
        the model.

        **Notes:**

        qpos[0:3] = position_E_Eo: The current position of the first Airplane's CG (in
        Earth axes, relative to the Earth origin) in meters.

        qvel[0:3] = velocity_E__E: The current velocity of the first Airplane's CG (in
        Earth axes, observed from the Earth frame) in meters per second.

        np.rad2deg(qvel[3:6]) = omegas_BP1__E: The current angular velocity of the first
        Airplane's body axes (in the first Airplane's body axes, observed from the Earth
        frame) in degrees per second.

        xmat = R_pas_BP1_to_E: The current orientation of the first Airplane as a
        passive rotation matrix from the first Airplane's body axes to Earth axes.

        We define MuJoCo world coordinates to be identical to Ptera Software Earth axes.

        :return: A dictionary containing the following keys: ``position_E_Eo``, a (3,)
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
        R_pas_BP1_to_E = self._data.xmat[self._body_id].reshape(3, 3)
        # TODO: Consider creating an invert_R_pas function in _transformations.py
        #  and calling it here.
        R_pas_E_to_BP1 = R_pas_BP1_to_E.T

        return {
            "position_E_Eo": np.copy(self._data.qpos[0:3]),
            "R_pas_E_to_BP1": np.copy(R_pas_E_to_BP1),
            "velocity_E__E": np.copy(self._data.qvel[0:3]),
            "omegas_BP1__E": np.rad2deg(np.copy(self._data.qvel[3:6])),
            "time": float(self._data.time),
        }

    def reset(self) -> None:
        """Resets the model's state to the initial conditions, time to zero seconds, and
        removes any applied loads.

        :return: None
        """
        # Reset the model's state to the initial conditions.
        self._data.qpos[:] = self._initial_qpos
        self._data.qvel[:] = self._initial_qvel

        # Reset time to zero seconds.
        self._data.time = 0.0

        # Remove any applied loads.
        self._data.xfrc_applied[:] = 0.0

        # Run forward kinematics to update dependent quantities.
        mujoco.mj_forward(self._model, self._data)

    def save_state(self) -> np.ndarray:
        """Saves and returns a snapshot of the model's current integration state.

        **Notes:**

        The snapshot is MuJoCo's integration state (the generalized positions and
        velocities, the time, the applied loads in xfrc_applied, and the constraint
        solver warm start), captured through MuJoCo's native state API. Its size and
        buffer layout vary across MuJoCo versions, so the size is queried from
        mj_stateSize at runtime and the snapshot is treated as opaque. It is transient
        within a time step and must never be serialized. Pair it with restore_state to
        return the model to this state.

        :return: A (N,) ndarray of floats representing the model's integration state
            snapshot, where N is the integration state size reported by mj_stateSize.
        """
        state_size = mujoco.mj_stateSize(
            self._model, mujoco.mjtState.mjSTATE_INTEGRATION
        )
        state = np.empty(state_size, dtype=float)
        mujoco.mj_getState(
            self._model, self._data, state, mujoco.mjtState.mjSTATE_INTEGRATION
        )

        return state

    def restore_state(self, state: np.ndarray) -> None:
        """Restores the model to a previously saved integration state.

        **Notes:**

        Takes a snapshot produced by save_state and writes it back through MuJoCo's
        native state API, then runs forward kinematics so the derived quantities (such
        as the orientation matrix get_state reads) match the restored configuration.
        Restoring the constraint solver warm start alongside the generalized positions
        and velocities, the time, and the applied loads makes re-stepping from the
        restored state reproducible.

        :param state: A (N,) ndarray of floats representing an integration state
            snapshot produced by this model's save_state, where N is the integration
            state size reported by mj_stateSize.
        :return: None
        """
        mujoco.mj_setState(
            self._model, self._data, state, mujoco.mjtState.mjSTATE_INTEGRATION
        )

        # Run forward kinematics to update dependent quantities.
        mujoco.mj_forward(self._model, self._data)

    def _rebuild_engine(self) -> None:
        """Rebuilds the native MuJoCo model and data objects from the stored XML string.

        The serialization layer restores the picklable slots (the XML string, the body
        and keyframe ids, and the initial qpos and qvel) but cannot restore the live
        MuJoCo model and data objects, which wrap native MuJoCo state. This recreates
        them from the XML string and resets the data to the initial keyframe, matching
        the state established at construction. The live, post-run data state is not
        restored: the canonical per-step state lives in the FreeFlightUnsteadyProblem's
        steady problems. Models built with mujoco_assets (such as meshes) never reach
        this method: the serialization layer refuses to save them, since the rebuilt
        engine could not resolve the XML string's asset references.

        :return: None
        """
        # noinspection PyArgumentList
        model = mujoco.MjModel.from_xml_string(self._xml_str)
        data = mujoco.MjData(model)
        mujoco.mj_resetDataKeyframe(model, data, self._initial_key_frame_id)
        mujoco.mj_forward(model, data)
        self._model = model
        self._data = data
