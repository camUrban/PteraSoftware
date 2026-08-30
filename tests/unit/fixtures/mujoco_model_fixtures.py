"""This module contains functions to create MuJoCoModels for use in tests."""

import struct
import tempfile
from collections.abc import Sequence
from pathlib import Path

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


def make_tetrahedron_vertices_fixture() -> np.ndarray:
    """This method makes a fixture that is the authored vertices of the tetrahedron
    whose STL bytes make_tetrahedron_stl_bytes_fixture returns.

    :return tetrahedron_vertices_fixture: ndarray This is a (4,3) ndarray of floats
        representing the tetrahedron's vertices (in the mesh file's axes, relative to
        the mesh file's origin).
    """
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def make_tetrahedron_stl_bytes_fixture() -> bytes:
    """This method makes a fixture that is the binary STL contents of a tetrahedron.

    The tetrahedron's vertices are those returned by make_tetrahedron_vertices_fixture.
    MuJoCo requires a mesh to have at least four vertices, so a tetrahedron is the
    smallest mesh a model can carry.

    :return tetrahedron_stl_bytes_fixture: bytes This is the binary STL contents of the
        tetrahedron.
    """
    vertices = make_tetrahedron_vertices_fixture()
    faces = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]

    header = b"\0" * 80
    triangles = b""
    for i0, i1, i2 in faces:
        normal = struct.pack("<fff", 0.0, 0.0, 0.0)
        triangle = b"".join(struct.pack("<fff", *vertices[i]) for i in (i0, i1, i2))
        triangles += normal + triangle + struct.pack("<H", 0)

    return header + struct.pack("<I", len(faces)) + triangles


def make_render_geometry_mujoco_model_fixture() -> _mujoco_model.MuJoCoModel:
    """This method makes a fixture that is a MuJoCoModel whose extra_xml injects one
    geom of every drawable type.

    The worldbody carries a finite plane, and the body carries a box, a sphere with a
    red material, a capsule, a cylinder, an ellipsoid, and a tetrahedron mesh loaded
    from mujoco_assets. The compiled model's geom order follows the document order, so
    get_render_geometry returns the geoms in the order just listed. The box declares its
    own rgba, the sphere takes its color from its material, and the mesh geom is posed
    with a translation and a 90 degree rotation about its parent's y axis.

    :return render_geometry_mujoco_model_fixture: MuJoCoModel This is the MuJoCoModel
        with one geom of every drawable type.
    """
    extra_xml = {
        "asset": (
            "<asset>"
            '<mesh name="tetrahedron" file="tetrahedron.stl"/>'
            '<material name="red_material" rgba="1 0 0 1"/>'
            "</asset>"
        ),
        "worldbody": '<geom name="ground" type="plane" size="5 4 0.1" pos="0 0 2"/>',
        "body": (
            '<geom name="box_geom" type="box" size="0.1 0.2 0.3" pos="1 2 3" '
            'rgba="0.1 0.2 0.3 0.4"/>'
            '<geom name="sphere_geom" type="sphere" size="0.5" '
            'material="red_material"/>'
            '<geom name="capsule_geom" type="capsule" size="0.1 0.4"/>'
            '<geom name="cylinder_geom" type="cylinder" size="0.2 0.5"/>'
            '<geom name="ellipsoid_geom" type="ellipsoid" size="0.1 0.2 0.3"/>'
            '<geom name="mesh_geom" type="mesh" mesh="tetrahedron" pos="1 2 3" '
            'euler="0 90 0"/>'
        ),
    }
    mujoco_assets = {"tetrahedron.stl": make_tetrahedron_stl_bytes_fixture()}

    render_geometry_mujoco_model_fixture = _mujoco_model.MuJoCoModel(
        name="render_geometry_airplane",
        mass=1.0,
        omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
        T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
        vCg_E__E=np.array((10.0, 0.0, 0.0)),
        I_BP1_CgP1=np.eye(3, dtype=float),
        delta_time=0.01,
        extra_xml=extra_xml,
        mujoco_assets=mujoco_assets,
    )

    return render_geometry_mujoco_model_fixture


def make_disk_mesh_mujoco_model_fixture() -> _mujoco_model.MuJoCoModel:
    """This method makes a fixture that is a MuJoCoModel whose XML references a mesh by
    an absolute on-disk path instead of through mujoco_assets.

    The tetrahedron STL bytes are written to a temporary file that MuJoCo reads at
    compile time, and the file is deleted before returning, which leaves the model
    holding a file reference that resolved only while the file existed on this machine.

    :return disk_mesh_mujoco_model_fixture: MuJoCoModel This is the MuJoCoModel whose
        XML references a mesh by an absolute on-disk path.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        stl_path = Path(temp_dir) / "tetrahedron.stl"
        stl_path.write_bytes(make_tetrahedron_stl_bytes_fixture())

        extra_xml = {
            "asset": f'<asset><mesh name="tetrahedron" file="{stl_path}"/></asset>',
            "body": '<geom name="mesh_geom" type="mesh" mesh="tetrahedron"/>',
        }
        disk_mesh_mujoco_model_fixture = _mujoco_model.MuJoCoModel(
            name="disk_mesh_airplane",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=np.eye(3, dtype=float),
            delta_time=0.01,
            extra_xml=extra_xml,
        )

    return disk_mesh_mujoco_model_fixture


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
