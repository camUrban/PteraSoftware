"""Contains the class definition for interfacing with MuJoCo for free flight
simulations.

**Contains the following classes:**

MuJoCoModel: A class that wraps MuJoCo models and data objects to provide a clean
interface for free flight simulations.

**Contains the following functions:**

None
"""

from __future__ import annotations

from collections.abc import Sequence

import mujoco
import numpy as np

from . import _parameter_validation


# TEST: Add unit tests for this class's initialization.
class MuJoCoModel:
    """A class that wraps MuJoCo models and data objects to provide a clean interface
    for free flight simulations.

    Provides methods for applying aerodynamic loads to the first Airplane, advancing the
    MuJoCo simulation, and extracting the current state of the first Airplane.

    **Contains the following methods:**

    apply_loads: Applies loads to the first Airplane.

    step: Advances the MuJoCo simulation by one time step.

    get_state: Extracts the current position, orientation, velocity, and angular
    velocity of the first Airplane.

    reset: Resets the simulation to the initial conditions.
    """

    def __init__(
        self,
        xml: str,
        body_name: str,
        initial_key_frame_name: str,
        delta_time: float | int | None = None,
    ) -> None:
        """The initialization method.

        :param xml: A str that contains either the complete MJCF (MuJoCo XML) model
            definition, or a file path to an XML file containing the model. If the str
            ends with '.xml', it is treated as a file path. Otherwise, it is treated as
            an XML.
        :param body_name: The name of the body (as defined in the MJCF model) to which
            aerodynamic loads will be applied. This body should typically have a
            freejoint to allow 6-DOF motion.
        :param initial_key_frame_name: The name of the initial condition key frame (as
            defined in the MJCF model).
        :param delta_time: A positive number (int or float) representing the time step
            for the MuJoCo simulation or None. If None, the time step defined in the XML
            model will be used. If not None, it will be converted internally to a float.
            The units are seconds.
        :return: None
        """
        xml = _parameter_validation.str_return_str(xml, "xml")
        self.body_name = _parameter_validation.str_return_str(body_name, "body_name")
        self.initial_key_frame_name = _parameter_validation.str_return_str(
            initial_key_frame_name, "initial_key_frame_name"
        )

        # REFACTOR: Validate xml contents or file existence?
        # Load the MuJoCo model from XML string or file path.
        self.model: mujoco.MjModel
        if xml.endswith(".xml"):
            # noinspection PyArgumentList
            self.model = mujoco.MjModel.from_xml_path(xml)
        else:
            # noinspection PyArgumentList
            self.model = mujoco.MjModel.from_xml_string(xml)

        # Set the time step if provided.
        if delta_time is not None:
            delta_time = _parameter_validation.number_in_range_return_float(
                delta_time, "delta_time", min_val=0.0, min_inclusive=False
            )
            self.model.opt.timestep = delta_time

        # Create the MuJoCo data structure.
        self.data: mujoco.MjData = mujoco.MjData(self.model)

        # Get and store the body ID and the initial conditions key frame ID.
        self.body_id: int = mujoco.mj_name2id(
            self.model, mujoco.mjtObj.mjOBJ_BODY, self.body_name
        )
        self.initial_key_frame_id: int = mujoco.mj_name2id(
            self.model, mujoco.mjtObj.mjOBJ_KEY, self.initial_key_frame_name
        )

        if self.body_id == -1:
            raise ValueError(
                f"Body with name '{self.body_name}' not found in the MuJoCo model."
            )
        if self.initial_key_frame_id == -1:
            raise ValueError(
                f"Key frame with name '{self.initial_key_frame_name}' not found in "
                f"the MuJoCo model."
            )

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
        """Applies loads to the first Airplane.

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

        # Pack the force and moment into the 6-element xfrc_applied array.
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
        the first Airplane.

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
        """Resets the simulation to the initial conditions.

        :return: None
        """
        # Reset the positions and velocities to their initial values.
        self.data.qpos[:] = self._initial_qpos
        self.data.qvel[:] = self._initial_qvel

        # Reset time to zero.
        self.data.time = 0.0

        # Reset applied forces.
        self.data.xfrc_applied[:] = 0.0

        # Run forward kinematics to update dependent quantities.
        mujoco.mj_forward(self.model, self.data)
