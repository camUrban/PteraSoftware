"""This module contains functions to create MuJoCoModels for use in tests."""

from collections.abc import Sequence

import numpy as np

# noinspection PyProtectedMember
from pterasoftware import _mujoco_model, _transformations


def make_basic_mujoco_model_fixture() -> _mujoco_model.MuJoCoModel:
    """This method makes a fixture that is a MuJoCoModel with basic parameters
    representing a simple rigid body at rest with identity orientation.

    :return basic_mujoco_model_fixture: MuJoCoModel This is the MuJoCoModel with basic
        parameters.
    """
    basic_mujoco_model_fixture = _mujoco_model.MuJoCoModel(
        name="test_airplane",
        mass=1.0,
        omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
        T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
        vCg_E__E=np.array((10.0, 0.0, 0.0)),
        I_BP1_CgP1=np.eye(3, dtype=float),
        delta_time=0.01,
    )

    return basic_mujoco_model_fixture


def make_rotated_mujoco_model_fixture() -> _mujoco_model.MuJoCoModel:
    """This method makes a fixture that is a MuJoCoModel with a 90 degree rotation about
    the z axis and non zero initial angular velocity.

    :return rotated_mujoco_model_fixture: MuJoCoModel This is the MuJoCoModel with
        rotated initial orientation.
    """
    T_pas = np.eye(4, dtype=float)
    T_pas[:3, :3] = np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
    )

    rotated_mujoco_model_fixture = _mujoco_model.MuJoCoModel(
        name="rotated_airplane",
        mass=2.0,
        omegas_BP1__E=np.array((0.0, 0.0, 10.0)),
        T_pas_BP1_CgP1_to_E_CgP1=T_pas,
        vCg_E__E=np.array((0.0, 5.0, -1.0)),
        I_BP1_CgP1=np.diag([2.0, 3.0, 4.0]),
        delta_time=0.005,
    )

    return rotated_mujoco_model_fixture


def make_pitched_mujoco_model_fixture(
    omegas_BP1__E: np.ndarray | Sequence[float | int] = (0.0, 0.0, 0.0),
) -> _mujoco_model.MuJoCoModel:
    """This method makes a fixture that is a MuJoCoModel pitched 90 degrees about the y
    axis, with isotropic inertia and unit mass.

    A 90 degree pitch places the body's +x direction along Earth -z, which makes this
    fixture useful for verifying the MuJoCo axis conventions documented in
    MUJOCO_CONVENTIONS.md. The isotropic inertia and unit mass keep the rotational and
    translational responses independent of orientation. The orientation is built from an
    Euler angle vector through the transformation helpers rather than a hand-written
    matrix.

    :param omegas_BP1__E: A (3,) array-like of floats representing the initial angular
        velocity of the body axes (in body axes, observed from the Earth frame), in
        degrees per second. The default is no rotation, which leaves the body at rest.
    :return pitched_mujoco_model_fixture: MuJoCoModel This is the MuJoCoModel pitched 90
        degrees about the y axis.
    """
    angles_E_to_BP1_izyx = np.array([0.0, 90.0, 0.0], dtype=float)
    T_pas_E_CgP1_to_BP1_CgP1 = _transformations.generate_rot_T(
        angles=angles_E_to_BP1_izyx, passive=True, intrinsic=True, order="zyx"
    )
    T_pas_BP1_CgP1_to_E_CgP1 = _transformations.invert_T_pas(T_pas_E_CgP1_to_BP1_CgP1)

    pitched_mujoco_model_fixture = _mujoco_model.MuJoCoModel(
        name="pitched_airplane",
        mass=1.0,
        omegas_BP1__E=np.array(omegas_BP1__E, dtype=float),
        T_pas_BP1_CgP1_to_E_CgP1=T_pas_BP1_CgP1_to_E_CgP1,
        vCg_E__E=np.array((0.0, 0.0, 0.0)),
        I_BP1_CgP1=np.eye(3, dtype=float),
        delta_time=0.01,
    )

    return pitched_mujoco_model_fixture


def make_basic_mujoco_model_name_fixture() -> str:
    """This method makes a fixture that is the name used by the basic MuJoCoModel
    fixture.

    :return: str The name of the basic MuJoCoModel fixture.
    """
    return "test_airplane"


def make_basic_mujoco_model_mass_fixture() -> float:
    """This method makes a fixture that is the mass used by the basic MuJoCoModel
    fixture.

    :return: float The mass in kilograms.
    """
    return 1.0


def make_basic_mujoco_model_delta_time_fixture() -> float:
    """This method makes a fixture that is the delta_time used by the basic MuJoCoModel
    fixture.

    :return: float The time step in seconds.
    """
    return 0.01
