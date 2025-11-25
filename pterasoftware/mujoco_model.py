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
from . import _transformations


# TEST: Add unit tests for this class's initialization.
class MuJoCoModel:
    """A class that wraps MuJoCo models and data objects to provide a clean interface
    for free flight simulations.

    Provides methods for applying aerodynamic loads to a rigid body, advancing the
    MuJoCo simulation, and extracting the current state of the body.

    **Contains the following methods:**

    apply_loads: Applies loads to the body.

    step: Advances the MuJoCo simulation by one time step.

    get_state: Extracts the current position, orientation, velocity, and angular
    velocity of the body.

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
        moments_E_Cg: np.ndarray | Sequence[float | int],
    ) -> None:
        """Applies loads to the body.

        The forces are in Earth axes and moments are in Earth axes relative to the
        body's CG. Both are applied at the body's CG. They will persist until the next
        call to apply_loads or until they are explicitly cleared.

        :param forces_E: A (3,) array-like object of numbers (int or float) representing
            the forces (in Earth axes) to apply to the body at its CG. Can be a tuple,
            list, or ndarray. Values are converted to floats internally. The units are
            in Newtons.
        :param moments_E_Cg: A (3,) array-like object of numbers (int or float)
            representing the moments (in Earth axes, relative to the body's CG) to apply
            to the body at its CG. Can be a tuple, list, or ndarray. Values are
            converted to floats internally. The units are in Newton-meters.
        :return: None
        """
        forces_E = _parameter_validation.threeD_number_vectorLike_return_float(
            forces_E, "forces_E"
        )
        moments_E_Cg = _parameter_validation.threeD_number_vectorLike_return_float(
            moments_E_Cg, "moments_E_Cg"
        )

        # Pack the force and moment into the 6-element xfrc_applied array.
        self.data.xfrc_applied[self.body_id][:] = np.hstack([forces_E, moments_E_Cg])

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
        the body.

        :return: A dictionary containing the following keys: ``position_E_E``, a (3,)
            ndarray of floats representing the current position of the body's CG in
            Earth axes relative to the Earth origin in meters; ``R_pas_E_to_B``, a (3,3)
            ndarray of floats representing the current orientation of the body as a
            passive rotation matrix from Earth axes to body axes; ``velocity_E__E``, a
            (3,) ndarray of floats representing the current velocity of the body's CG in
            Earth axes as observed from the Earth frame in meters per second;
            ``omegas_B__E``, a (3,) ndarray of floats representing the current angular
            velocity of the body axes in Earth axes as observed from the Earth frame in
            degrees per second; ``time``, a float representing the current simulation
            time in seconds.
        """
        # REFACTOR: This transformation logic is EXTREMELY janky. We should clean it
        #  up and get a better understanding of what is actually happening. I
        #  verified that it works with alpha > 0, beta = 0, and angles_E_to_B_ixyz =
        #  (0, 0, 0). However, I have no idea if it will work with other parameters.
        #  What does xmat[self.body_id].reshape(3, 3) actually represent? How do we
        #  incorporate both alpha and beta, and the orientation angles.
        R_mujoco = np.copy(self.data.xmat[self.body_id].reshape(3, 3))
        R_pas_W_Cg_to_G_Cg = (
            np.asarray(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], dtype=float
            )
            @ R_mujoco.T
        )
        T_pas_W_Cg_to_G_Cg = np.eye(4, dtype=float)
        T_pas_W_Cg_to_G_Cg[:3, :3] = R_pas_W_Cg_to_G_Cg
        T_pas_G_Cg_to_B_Cg = _transformations.generate_rot_T(
            angles=(0, 180, 0), passive=True, intrinsic=True, order="zyx"
        )
        T_pas_W_Cg_to_B_Cg = _transformations.compose_T_pas(
            T_pas_W_Cg_to_G_Cg, T_pas_G_Cg_to_B_Cg
        )

        return {
            "position_E_E": np.copy(self.data.qpos[0:3]),
            "R_pas_E_to_B": T_pas_W_Cg_to_B_Cg[0:3, 0:3],
            "velocity_E__E": np.copy(self.data.qvel[0:3]),
            "omegas_B__E": np.rad2deg(np.copy(self.data.qvel[3:6])),
            "time": self.data.time,
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
