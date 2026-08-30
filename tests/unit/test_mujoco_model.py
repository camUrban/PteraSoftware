"""This module contains classes to test MuJoCoModels."""

import struct
import unittest

import mujoco
import numpy as np
import numpy.testing as npt
import pyvista as pv

# noinspection PyProtectedMember
from pterasoftware import _mujoco_model, _transformations
from tests.unit.fixtures import mujoco_model_fixtures


class TestMuJoCoModelInit(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel initialization."""

    def setUp(self) -> None:
        """Set up test fixtures for MuJoCoModel initialization tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_instantiation_returns_correct_type(self) -> None:
        """Test that MuJoCoModel instantiation returns a MuJoCoModel."""
        self.assertIsInstance(self.model, _mujoco_model.MuJoCoModel)

    def test_xml_str_contains_model_name(self) -> None:
        """Test that the generated XML contains the Airplane name."""
        name = mujoco_model_fixtures.make_basic_mujoco_model_name_fixture()
        self.assertIn(name, self.model.xml_str)

    def test_xml_str_contains_timestep(self) -> None:
        """Test that the generated XML contains the delta_time."""
        delta_time = mujoco_model_fixtures.make_basic_mujoco_model_delta_time_fixture()
        self.assertIn(str(delta_time), self.model.xml_str)

    def test_xml_str_contains_default_integrator(self) -> None:
        """Test that the generated XML defaults to the RK4 integrator."""
        self.assertIn('integrator="RK4"', self.model.xml_str)

    def test_accepts_integrator(self) -> None:
        """Test that MuJoCoModel accepts an integrator, injects it into the XML, and the
        compiled model uses it."""
        model = _mujoco_model.MuJoCoModel(
            name="integrator_test",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=np.eye(3, dtype=float),
            delta_time=0.01,
            integrator="implicitfast",
        )
        self.assertIn('integrator="implicitfast"', model.xml_str)
        self.assertEqual(
            model._model.opt.integrator, mujoco.mjtIntegrator.mjINT_IMPLICITFAST
        )

    def test_model_is_mj_model(self) -> None:
        """Test that the internal model object is an MjModel."""
        self.assertIsInstance(self.model._model, mujoco.MjModel)

    def test_data_is_mj_data(self) -> None:
        """Test that the internal data object is an MjData."""
        self.assertIsInstance(self.model._data, mujoco.MjData)

    def test_body_id_is_int(self) -> None:
        """Test that body_id is an int."""
        self.assertIsInstance(self.model.body_id, int)

    def test_initial_key_frame_id_is_int(self) -> None:
        """Test that initial_key_frame_id is an int."""
        self.assertIsInstance(self.model.initial_key_frame_id, int)

    def test_initial_qpos_shape(self) -> None:
        """Test that initial_qpos has the correct shape for a free joint."""
        self.assertEqual(self.model.initial_qpos.shape, (7,))

    def test_initial_qvel_shape(self) -> None:
        """Test that initial_qvel has the correct shape for a free joint."""
        self.assertEqual(self.model.initial_qvel.shape, (6,))

    def test_initial_qvel_contains_velocity(self) -> None:
        """Test that initial_qvel contains the initial velocity."""
        npt.assert_allclose(self.model.initial_qvel[0:3], [10.0, 0.0, 0.0], atol=1e-14)

    def test_initial_qvel_contains_angular_velocity(self) -> None:
        """Test that initial_qvel contains the initial angular velocity in radians per
        second."""
        npt.assert_allclose(self.model.initial_qvel[3:6], [0.0, 0.0, 0.0], atol=1e-14)

    def test_mass_set_on_model(self) -> None:
        """Test that the MuJoCo body mass equals the supplied mass."""
        expected_mass = mujoco_model_fixtures.make_basic_mujoco_model_mass_fixture()

        body_id = self.model.body_id
        actual_mass = self.model._model.body_mass[body_id]
        self.assertAlmostEqual(actual_mass, expected_mass, places=10)

    def test_timestep_set_on_model(self) -> None:
        """Test that the MuJoCo model timestep matches delta_time."""
        delta_time = mujoco_model_fixtures.make_basic_mujoco_model_delta_time_fixture()
        self.assertAlmostEqual(self.model._model.opt.timestep, delta_time, places=14)

    def test_gravity_disabled_in_mujoco(self) -> None:
        """Test that gravity is set to zero in the MuJoCo model."""
        npt.assert_array_equal(self.model._model.opt.gravity, [0.0, 0.0, 0.0])

    def test_accepts_extra_xml(self) -> None:
        """Test that MuJoCoModel accepts extra_xml and injects it into the XML."""
        extra_xml = {"visual": '<visual><quality shadowsize="2048"/></visual>'}
        model = _mujoco_model.MuJoCoModel(
            name="extra_xml_test",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=np.eye(3, dtype=float),
            delta_time=0.01,
            extra_xml=extra_xml,
        )
        self.assertIn("shadowsize", model.xml_str)

    def test_accepts_mujoco_assets(self) -> None:
        """Test that MuJoCoModel accepts mujoco_assets without error."""
        header = b"\x00" * 80
        verts = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ]
        faces = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
        triangles = b""
        for i0, i1, i2 in faces:
            normal = struct.pack("<fff", 0.0, 0.0, 0.0)
            tri = b"".join(struct.pack("<fff", *verts[i]) for i in (i0, i1, i2))
            triangles += normal + tri + struct.pack("<H", 0)
        stl_bytes = header + struct.pack("<I", len(faces)) + triangles

        mujoco_assets = {"dummy.stl": stl_bytes}
        extra_xml = {
            "asset": '<asset><mesh name="dummy" file="dummy.stl"/></asset>',
            "body": '<geom type="mesh" mesh="dummy"/>',
        }
        model = _mujoco_model.MuJoCoModel(
            name="assets_test",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=np.eye(3, dtype=float),
            delta_time=0.01,
            extra_xml=extra_xml,
            mujoco_assets=mujoco_assets,
        )
        self.assertIsInstance(model, _mujoco_model.MuJoCoModel)
        self.assertIs(model._mujoco_assets, mujoco_assets)

    def test_mujoco_assets_default_none(self) -> None:
        """Test that the retained mujoco_assets defaults to None."""
        self.assertIsNone(self.model._mujoco_assets)

    def test_rotated_initial_orientation(self) -> None:
        """Test that a rotated initial orientation produces correct initial state."""
        model = mujoco_model_fixtures.make_rotated_mujoco_model_fixture()
        state = model.get_state()

        R = state["R_pas_E_to_BP1"]
        expected_R = np.array(
            [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
        )
        npt.assert_allclose(R, expected_R, atol=1e-10)

    def test_non_zero_initial_angular_velocity(self) -> None:
        """Test that non zero initial angular velocity is stored correctly."""
        model = mujoco_model_fixtures.make_rotated_mujoco_model_fixture()
        state = model.get_state()
        npt.assert_allclose(state["omegas_BP1__E"], [0.0, 0.0, 10.0], atol=1e-10)

    def test_symmetrizes_inertia_matrix(self) -> None:
        """Test that an asymmetric inertia matrix is symmetrized."""
        I_asymmetric = np.array(
            [[1.0, 0.2, 0.0], [0.4, 2.0, 0.0], [0.0, 0.0, 3.0]], dtype=float
        )
        model = _mujoco_model.MuJoCoModel(
            name="sym_test",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=I_asymmetric,
            delta_time=0.01,
        )
        expected_Ixy = (0.2 + 0.4) / 2
        self.assertIn(str(expected_Ixy), model.xml_str)


class TestMuJoCoModelImmutability(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel immutability patterns."""

    def setUp(self) -> None:
        """Set up test fixtures for MuJoCoModel immutability tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_immutable_xml_str_raises_attribute_error(self) -> None:
        """Test that setting xml_str raises AttributeError."""
        with self.assertRaises(AttributeError):
            setattr(self.model, "xml_str", "new xml")

    def test_immutable_body_id_raises_attribute_error(self) -> None:
        """Test that setting body_id raises AttributeError."""
        with self.assertRaises(AttributeError):
            setattr(self.model, "body_id", 99)

    def test_immutable_initial_key_frame_id_raises_attribute_error(self) -> None:
        """Test that setting initial_key_frame_id raises AttributeError."""
        with self.assertRaises(AttributeError):
            setattr(self.model, "initial_key_frame_id", 99)

    def test_immutable_initial_qpos_raises_attribute_error(self) -> None:
        """Test that setting initial_qpos raises AttributeError."""
        with self.assertRaises(AttributeError):
            setattr(self.model, "initial_qpos", np.zeros(7))

    def test_immutable_initial_qvel_raises_attribute_error(self) -> None:
        """Test that setting initial_qvel raises AttributeError."""
        with self.assertRaises(AttributeError):
            setattr(self.model, "initial_qvel", np.zeros(6))

    def test_initial_qpos_is_read_only_array(self) -> None:
        """Test that initial_qpos array is not writeable."""
        with self.assertRaises(ValueError):
            self.model.initial_qpos[0] = 999.0

    def test_initial_qvel_is_read_only_array(self) -> None:
        """Test that initial_qvel array is not writeable."""
        with self.assertRaises(ValueError):
            self.model.initial_qvel[0] = 999.0


class TestMuJoCoModelApplyLoads(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.apply_loads."""

    def setUp(self) -> None:
        """Set up test fixtures for apply_loads tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_apply_loads_sets_forces(self) -> None:
        """Test that apply_loads sets the force on the body."""
        forces_E = np.array([1.0, 2.0, 3.0])
        moments_E_CgP1 = np.array([0.0, 0.0, 0.0])
        self.model.apply_loads(forces_E, moments_E_CgP1)

        applied = self.model._data.xfrc_applied[self.model.body_id]
        npt.assert_allclose(applied[0:3], forces_E, atol=1e-14)

    def test_apply_loads_sets_moments(self) -> None:
        """Test that apply_loads sets the moment on the body."""
        forces_E = np.array([0.0, 0.0, 0.0])
        moments_E_CgP1 = np.array([4.0, 5.0, 6.0])
        self.model.apply_loads(forces_E, moments_E_CgP1)

        applied = self.model._data.xfrc_applied[self.model.body_id]
        npt.assert_allclose(applied[3:6], moments_E_CgP1, atol=1e-14)

    def test_apply_loads_overwrites_previous(self) -> None:
        """Test that apply_loads overwrites previously applied loads."""
        self.model.apply_loads(np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0]))
        self.model.apply_loads(np.array([10.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))

        applied = self.model._data.xfrc_applied[self.model.body_id]
        npt.assert_allclose(applied[0:3], [10.0, 0.0, 0.0], atol=1e-14)
        npt.assert_allclose(applied[3:6], [0.0, 0.0, 0.0], atol=1e-14)


class TestMuJoCoModelStep(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.step."""

    def setUp(self) -> None:
        """Set up test fixtures for step tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_step_advances_time(self) -> None:
        """Test that step advances the simulation time by delta_time."""
        delta_time = mujoco_model_fixtures.make_basic_mujoco_model_delta_time_fixture()
        self.model.step()
        self.assertAlmostEqual(self.model._data.time, delta_time, places=14)

    def test_multiple_steps_advance_time(self) -> None:
        """Test that multiple steps advance time correctly."""
        delta_time = mujoco_model_fixtures.make_basic_mujoco_model_delta_time_fixture()
        num_steps = 10
        for _ in range(num_steps):
            self.model.step()
        self.assertAlmostEqual(self.model._data.time, num_steps * delta_time, places=10)

    def test_step_with_force_changes_velocity(self) -> None:
        """Test that stepping with an applied force changes the velocity."""
        initial_state = self.model.get_state()
        initial_velocity = initial_state["velocity_E__E"].copy()

        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        self.model.step()

        final_state = self.model.get_state()
        final_velocity = final_state["velocity_E__E"]

        self.assertFalse(np.allclose(initial_velocity, final_velocity))


class TestMuJoCoModelGetState(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.get_state."""

    def setUp(self) -> None:
        """Set up test fixtures for get_state tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_get_state_returns_dict(self) -> None:
        """Test that get_state returns a dict."""
        state = self.model.get_state()
        self.assertIsInstance(state, dict)

    def test_get_state_contains_expected_keys(self) -> None:
        """Test that get_state returns all expected keys."""
        state = self.model.get_state()
        expected_keys = {
            "position_E_Eo",
            "R_pas_E_to_BP1",
            "velocity_E__E",
            "omegas_BP1__E",
            "time",
        }
        self.assertEqual(set(state.keys()), expected_keys)

    def test_position_shape(self) -> None:
        """Test that position_E_Eo has shape (3,)."""
        state = self.model.get_state()
        self.assertEqual(state["position_E_Eo"].shape, (3,))

    def test_rotation_matrix_shape(self) -> None:
        """Test that R_pas_E_to_BP1 has shape (3,3)."""
        state = self.model.get_state()
        self.assertEqual(state["R_pas_E_to_BP1"].shape, (3, 3))

    def test_velocity_shape(self) -> None:
        """Test that velocity_E__E has shape (3,)."""
        state = self.model.get_state()
        self.assertEqual(state["velocity_E__E"].shape, (3,))

    def test_omegas_shape(self) -> None:
        """Test that omegas_BP1__E has shape (3,)."""
        state = self.model.get_state()
        self.assertEqual(state["omegas_BP1__E"].shape, (3,))

    def test_time_is_float(self) -> None:
        """Test that time is a float."""
        state = self.model.get_state()
        self.assertIsInstance(state["time"], float)

    def test_initial_position_at_origin(self) -> None:
        """Test that the initial position is at the origin."""
        state = self.model.get_state()
        npt.assert_allclose(state["position_E_Eo"], [0.0, 0.0, 0.0], atol=1e-14)

    def test_initial_time_is_zero(self) -> None:
        """Test that the initial time is zero."""
        state = self.model.get_state()
        self.assertAlmostEqual(state["time"], 0.0, places=14)

    def test_identity_rotation_at_init(self) -> None:
        """Test that identity T_pas produces identity R_pas_E_to_BP1."""
        state = self.model.get_state()
        npt.assert_allclose(state["R_pas_E_to_BP1"], np.eye(3), atol=1e-10)

    def test_initial_velocity_matches_input(self) -> None:
        """Test that the initial velocity matches the input vCg_E__E."""
        state = self.model.get_state()
        npt.assert_allclose(state["velocity_E__E"], [10.0, 0.0, 0.0], atol=1e-14)

    def test_get_state_returns_copies(self) -> None:
        """Test that get_state returns copies that do not alias internal data."""
        state1 = self.model.get_state()
        state1["position_E_Eo"][0] = 999.0
        state2 = self.model.get_state()
        self.assertNotAlmostEqual(state2["position_E_Eo"][0], 999.0)


class TestMuJoCoModelReset(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.reset."""

    def setUp(self) -> None:
        """Set up test fixtures for reset tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_reset_restores_time_to_zero(self) -> None:
        """Test that reset restores time to zero."""
        self.model.step()
        self.model.step()
        self.model.reset()
        self.assertAlmostEqual(self.model._data.time, 0.0, places=14)

    def test_reset_restores_initial_qpos(self) -> None:
        """Test that reset restores initial generalized positions."""
        initial_qpos = self.model.initial_qpos.copy()
        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        for _ in range(10):
            self.model.step()
        self.model.reset()
        npt.assert_allclose(self.model._data.qpos, initial_qpos, atol=1e-14)

    def test_reset_restores_initial_qvel(self) -> None:
        """Test that reset restores initial generalized velocities."""
        initial_qvel = self.model.initial_qvel.copy()
        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        for _ in range(10):
            self.model.step()
        self.model.reset()
        npt.assert_allclose(self.model._data.qvel, initial_qvel, atol=1e-14)

    def test_reset_clears_applied_loads(self) -> None:
        """Test that reset clears any applied loads."""
        self.model.apply_loads(
            np.array([100.0, 200.0, 300.0]), np.array([10.0, 20.0, 30.0])
        )
        self.model.reset()
        applied = self.model._data.xfrc_applied[self.model.body_id]
        npt.assert_array_equal(applied, np.zeros(6))

    def test_reset_produces_same_state_as_init(self) -> None:
        """Test that reset produces the same get_state output as after init."""
        initial_state = self.model.get_state()

        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 10.0]))
        for _ in range(20):
            self.model.step()

        self.model.reset()
        reset_state = self.model.get_state()

        npt.assert_allclose(
            reset_state["position_E_Eo"], initial_state["position_E_Eo"], atol=1e-14
        )
        npt.assert_allclose(
            reset_state["velocity_E__E"], initial_state["velocity_E__E"], atol=1e-14
        )
        npt.assert_allclose(
            reset_state["omegas_BP1__E"], initial_state["omegas_BP1__E"], atol=1e-14
        )
        npt.assert_allclose(
            reset_state["R_pas_E_to_BP1"],
            initial_state["R_pas_E_to_BP1"],
            atol=1e-10,
        )
        self.assertAlmostEqual(reset_state["time"], initial_state["time"], places=14)


class TestMuJoCoModelSaveState(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.save_state."""

    def setUp(self) -> None:
        """Set up test fixtures for save_state tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_save_state_returns_ndarray(self) -> None:
        """Test that save_state returns an ndarray."""
        snapshot = self.model.save_state()
        self.assertIsInstance(snapshot, np.ndarray)

    def test_save_state_is_one_dimensional(self) -> None:
        """Test that the saved snapshot is a one dimensional array."""
        snapshot = self.model.save_state()
        self.assertEqual(snapshot.ndim, 1)

    def test_save_state_is_float_array(self) -> None:
        """Test that the saved snapshot is a float array, as MuJoCo's state API
        requires."""
        snapshot = self.model.save_state()
        self.assertEqual(snapshot.dtype, np.float64)

    def test_save_state_length_matches_integration_state_size(self) -> None:
        """Test that the snapshot length matches MuJoCo's integration state size.

        The size is queried from mj_stateSize at runtime rather than hard coded, since
        the buffer layout varies across MuJoCo versions within the pin.
        """
        snapshot = self.model.save_state()
        expected_size = mujoco.mj_stateSize(
            self.model._model, mujoco.mjtState.mjSTATE_INTEGRATION
        )
        self.assertEqual(snapshot.shape, (expected_size,))

    def test_save_state_does_not_alias_internal_data(self) -> None:
        """Test that the snapshot does not alias the model's live data.

        Advancing the simulation after saving must leave a previously saved snapshot
        unchanged.
        """
        snapshot = self.model.save_state()
        snapshot_before = snapshot.copy()

        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        for _ in range(10):
            self.model.step()

        npt.assert_array_equal(snapshot, snapshot_before)

    def test_save_state_captures_distinct_states(self) -> None:
        """Test that snapshots taken at different states differ."""
        snapshot_initial = self.model.save_state()

        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        for _ in range(10):
            self.model.step()
        snapshot_advanced = self.model.save_state()

        self.assertFalse(np.array_equal(snapshot_initial, snapshot_advanced))


class TestMuJoCoModelRestoreState(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.restore_state."""

    def setUp(self) -> None:
        """Set up test fixtures for restore_state tests."""
        self.model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()

    def test_restore_state_restores_position(self) -> None:
        """Test that restore_state restores the position."""
        saved_state = self.model.get_state()
        snapshot = self.model.save_state()

        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        for _ in range(10):
            self.model.step()
        self.model.restore_state(snapshot)

        restored_state = self.model.get_state()
        npt.assert_allclose(
            restored_state["position_E_Eo"], saved_state["position_E_Eo"], atol=1e-14
        )

    def test_restore_state_restores_velocity(self) -> None:
        """Test that restore_state restores the velocity."""
        saved_state = self.model.get_state()
        snapshot = self.model.save_state()

        self.model.apply_loads(np.array([100.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        for _ in range(10):
            self.model.step()
        self.model.restore_state(snapshot)

        restored_state = self.model.get_state()
        npt.assert_allclose(
            restored_state["velocity_E__E"], saved_state["velocity_E__E"], atol=1e-14
        )

    def test_restore_state_restores_time(self) -> None:
        """Test that restore_state restores the simulation time."""
        for _ in range(5):
            self.model.step()
        saved_state = self.model.get_state()
        snapshot = self.model.save_state()

        for _ in range(10):
            self.model.step()
        self.model.restore_state(snapshot)

        restored_state = self.model.get_state()
        self.assertAlmostEqual(restored_state["time"], saved_state["time"], places=14)

    def test_restore_state_restores_orientation_and_angular_velocity(self) -> None:
        """Test that restore_state restores the orientation and angular velocity.

        Restoring must recompute the derived orientation matrix, so get_state reports
        the restored attitude even without an intervening step. A rotated model is
        stepped so both quantities evolve to a non trivial state before the snapshot.
        """
        model = mujoco_model_fixtures.make_rotated_mujoco_model_fixture()
        for _ in range(5):
            model.step()
        saved_state = model.get_state()
        snapshot = model.save_state()

        model.apply_loads(np.array([0.0, 0.0, 0.0]), np.array([10.0, 20.0, 30.0]))
        for _ in range(10):
            model.step()
        model.restore_state(snapshot)

        restored_state = model.get_state()
        npt.assert_allclose(
            restored_state["R_pas_E_to_BP1"], saved_state["R_pas_E_to_BP1"], atol=1e-12
        )
        npt.assert_allclose(
            restored_state["omegas_BP1__E"], saved_state["omegas_BP1__E"], atol=1e-12
        )

    def test_restore_state_restores_applied_loads(self) -> None:
        """Test that restore_state restores the applied loads.

        The integration state snapshot includes xfrc_applied, where apply_loads persists
        the applied loads, so restoring returns them.
        """
        self.model.apply_loads(np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0]))
        snapshot = self.model.save_state()

        self.model.apply_loads(np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))
        self.model.restore_state(snapshot)

        applied = self.model._data.xfrc_applied[self.model.body_id]
        npt.assert_array_equal(applied[0:3], [1.0, 2.0, 3.0])
        npt.assert_array_equal(applied[3:6], [4.0, 5.0, 6.0])

    def test_restore_state_enables_reproducible_restepping(self) -> None:
        """Test that re-stepping from a restored state reproduces the trajectory.

        Saving the state, stepping under a load, then restoring and stepping under the
        same load reproduces the first outcome exactly. This reproducible re-stepping is
        what the strongly coupled sub-iteration relies on.
        """
        self.model.apply_loads(np.array([5.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]))
        for _ in range(5):
            self.model.step()
        snapshot = self.model.save_state()

        self.model.apply_loads(np.array([20.0, 10.0, 0.0]), np.array([0.0, 1.0, 0.0]))
        self.model.step()
        first_state = self.model.get_state()

        self.model.restore_state(snapshot)

        self.model.apply_loads(np.array([20.0, 10.0, 0.0]), np.array([0.0, 1.0, 0.0]))
        self.model.step()
        second_state = self.model.get_state()

        npt.assert_array_equal(
            second_state["position_E_Eo"], first_state["position_E_Eo"]
        )
        npt.assert_array_equal(
            second_state["velocity_E__E"], first_state["velocity_E__E"]
        )
        npt.assert_array_equal(
            second_state["omegas_BP1__E"], first_state["omegas_BP1__E"]
        )
        npt.assert_array_equal(
            second_state["R_pas_E_to_BP1"], first_state["R_pas_E_to_BP1"]
        )
        self.assertEqual(second_state["time"], first_state["time"])


class TestMuJoCoModelConventions(unittest.TestCase):
    """This class contains methods verifying the MuJoCo axis conventions documented in
    MUJOCO_CONVENTIONS.md.

    Each test runs through MuJoCoModel's own translation between MuJoCo's freejoint
    state and Ptera Software's conventions (the qpos quaternion construction,
    get_state's xmat handling, and apply_loads's xfrc_applied packing), so it pins that
    contract rather than just MuJoCo. The model is pitched 90 degrees about the y axis,
    which places the body's +x direction along Earth -z, so a body-axes quantity and its
    Earth-axes counterpart point in clearly different directions. Assertions are
    directional rather than exact to stay robust to the integrator's small numerical
    drift.
    """

    def setUp(self) -> None:
        """Set up a pitched MuJoCoModel at rest for the convention tests."""
        self.model = mujoco_model_fixtures.make_pitched_mujoco_model_fixture()

    def test_xmat_columns_are_body_axes_in_earth_axes(self) -> None:
        """Test that xmat's columns are the body axes expressed in Earth axes.

        With a 90 degree pitch about the y axis, the body's +x direction points along
        Earth -z, its +y direction along Earth +y, and its +z direction along Earth +x.
        """
        xmat = self.model._data.xmat[self.model.body_id].reshape(3, 3)
        npt.assert_allclose(xmat[:, 0], [0.0, 0.0, -1.0], atol=1e-10)
        npt.assert_allclose(xmat[:, 1], [0.0, 1.0, 0.0], atol=1e-10)
        npt.assert_allclose(xmat[:, 2], [1.0, 0.0, 0.0], atol=1e-10)

    def test_angular_velocity_qvel_is_in_body_axes(self) -> None:
        """Test that qvel[3:6] is interpreted in body axes, not Earth axes.

        An angular velocity along the body's +x direction is applied to the pitched
        model, where the body's +x direction points along Earth -z. If qvel[3:6] were in
        Earth axes the body would spin about Earth +x; because it is in body axes, the
        body spins about Earth -z.
        """
        model = mujoco_model_fixtures.make_pitched_mujoco_model_fixture(
            omegas_BP1__E=(90.0, 0.0, 0.0)
        )

        for _ in range(20):
            model.step()

        state = model.get_state()
        R_pas_BP1_to_E = _transformations.invert_R_pas(state["R_pas_E_to_BP1"])
        omegas_E__E = _transformations.apply_R_to_vectors(
            R_pas_BP1_to_E, state["omegas_BP1__E"]
        )

        # The spin is about Earth -z, so the Earth-axes angular velocity is dominated by
        # a negative z component with negligible x and y components.
        self.assertLess(omegas_E__E[2], 0.0)
        self.assertGreater(abs(omegas_E__E[2]), 10.0 * abs(omegas_E__E[0]))
        self.assertGreater(abs(omegas_E__E[2]), 10.0 * abs(omegas_E__E[1]))

    def test_applied_torque_is_in_earth_axes(self) -> None:
        """Test that an applied torque is in Earth axes, independent of orientation.

        A torque about Earth +x is applied to the pitched model. If the torque were in
        body axes it would act about the body's +x direction (Earth -z); because it is
        in Earth axes, the body accelerates about Earth +x regardless of its
        orientation.
        """
        self.model.apply_loads(np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))

        for _ in range(20):
            self.model.step()

        state = self.model.get_state()
        R_pas_BP1_to_E = state["R_pas_E_to_BP1"].T
        omegas_E__E = R_pas_BP1_to_E @ state["omegas_BP1__E"]

        # The body spins up about Earth +x, so the Earth-axes angular velocity is
        # dominated by a positive x component with negligible y and z components.
        self.assertGreater(omegas_E__E[0], 0.0)
        self.assertGreater(abs(omegas_E__E[0]), 10.0 * abs(omegas_E__E[1]))
        self.assertGreater(abs(omegas_E__E[0]), 10.0 * abs(omegas_E__E[2]))

    def test_applied_force_is_in_earth_axes(self) -> None:
        """Test that an applied force is in Earth axes, independent of orientation.

        A force along Earth +x is applied to the pitched model. If the force were in
        body axes it would act along the body's +x direction (Earth -z); because it is
        in Earth axes, the body accelerates along Earth +x regardless of its
        orientation.
        """
        self.model.apply_loads(np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))

        for _ in range(20):
            self.model.step()

        velocity_E__E = self.model.get_state()["velocity_E__E"]

        # The body accelerates along Earth +x, so the velocity is dominated by a
        # positive x component with negligible y and z components.
        self.assertGreater(velocity_E__E[0], 0.0)
        self.assertGreater(abs(velocity_E__E[0]), 10.0 * abs(velocity_E__E[1]))
        self.assertGreater(abs(velocity_E__E[0]), 10.0 * abs(velocity_E__E[2]))


class TestMuJoCoModelGetRenderGeometry(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.get_render_geometry.

    The shared fixture's geoms arrive in document order: the worldbody plane, then the
    body's box, sphere, capsule, cylinder, ellipsoid, and tetrahedron mesh. The
    tessellated primitives (the sphere, capsule, cylinder, and ellipsoid) inscribe their
    ideal surfaces, so their bounds are compared with a loose tolerance, while the box,
    plane, and mesh vertices are exact.
    """

    render_geoms: list[_mujoco_model.RenderGeom]

    @classmethod
    def setUpClass(cls) -> None:
        """Set up the shared test fixtures."""
        cls.render_geoms = (
            mujoco_model_fixtures.make_render_geometry_mujoco_model_fixture()
        ).get_render_geometry()

    def test_returns_one_render_geom_per_drawable_geom(self) -> None:
        """Test that every drawable geom produces one RenderGeom."""
        self.assertEqual(len(self.render_geoms), 7)
        for render_geom in self.render_geoms:
            self.assertIsInstance(render_geom, _mujoco_model.RenderGeom)

    def test_meshes_are_polydata(self) -> None:
        """Test that every RenderGeom's mesh is a PolyData."""
        for render_geom in self.render_geoms:
            self.assertIsInstance(render_geom.mesh, pv.PolyData)

    def test_meshes_carry_no_point_normals(self) -> None:
        """Test that no RenderGeom's mesh carries stored point normals.

        The rendering layer recomputes shading normals from the posed triangles when it
        adds each geom actor, so stored normals are never read, and ones left behind by
        a PyVista primitive source would sit stale against the posed points.
        """
        for render_geom in self.render_geoms:
            self.assertNotIn("Normals", render_geom.mesh.point_data)

    def test_model_without_extra_geometry_returns_empty_list(self) -> None:
        """Test that a model with no extra_xml geoms returns an empty list."""
        model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()
        self.assertEqual(model.get_render_geometry(), [])

    def test_worldbody_geoms_are_not_body_attached(self) -> None:
        """Test that a worldbody geom's body_attached flag is False."""
        self.assertFalse(self.render_geoms[0].body_attached)

    def test_body_geoms_are_body_attached(self) -> None:
        """Test that every body geom's body_attached flag is True."""
        for render_geom in self.render_geoms[1:]:
            self.assertTrue(render_geom.body_attached)

    def test_rgba_defaults_to_the_geoms_rgba(self) -> None:
        """Test that a geom without a material takes its own rgba."""
        npt.assert_allclose(
            self.render_geoms[1].rgba, np.array([0.1, 0.2, 0.3, 0.4]), atol=1e-6
        )

    def test_rgba_uses_the_material_when_assigned(self) -> None:
        """Test that a geom with a material takes the material's rgba."""
        npt.assert_allclose(
            self.render_geoms[2].rgba, np.array([1.0, 0.0, 0.0, 1.0]), atol=1e-6
        )

    def test_plane_quad_matches_half_extents_and_position(self) -> None:
        """Test that a finite plane becomes a quad spanning its half-extents at its
        authored position."""
        npt.assert_allclose(
            self.render_geoms[0].mesh.bounds,
            (-5.0, 5.0, -4.0, 4.0, 2.0, 2.0),
            atol=1e-6,
        )

    def test_box_bounds_match_half_extents_and_position(self) -> None:
        """Test that a box's bounds span its half-extents about its authored
        position."""
        npt.assert_allclose(
            self.render_geoms[1].mesh.bounds,
            (0.9, 1.1, 1.8, 2.2, 2.7, 3.3),
            atol=1e-6,
        )

    def test_sphere_bounds_match_radius(self) -> None:
        """Test that a sphere's bounds span its radius."""
        npt.assert_allclose(
            self.render_geoms[2].mesh.bounds,
            (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
            atol=1e-2,
        )

    def test_capsule_bounds_match_radius_and_half_length(self) -> None:
        """Test that a capsule's axis is the local z axis, with bounds spanning its
        radius laterally and its half-length plus its radius axially."""
        npt.assert_allclose(
            self.render_geoms[3].mesh.bounds,
            (-0.1, 0.1, -0.1, 0.1, -0.5, 0.5),
            atol=1e-2,
        )

    def test_cylinder_bounds_match_radius_and_half_length(self) -> None:
        """Test that a cylinder's axis is the local z axis, with bounds spanning its
        radius laterally and its half-length axially."""
        npt.assert_allclose(
            self.render_geoms[4].mesh.bounds,
            (-0.2, 0.2, -0.2, 0.2, -0.5, 0.5),
            atol=1e-2,
        )

    def test_ellipsoid_bounds_match_radii(self) -> None:
        """Test that an ellipsoid's bounds span its three radii."""
        npt.assert_allclose(
            self.render_geoms[5].mesh.bounds,
            (-0.1, 0.1, -0.2, 0.2, -0.3, 0.3),
            atol=1e-2,
        )

    def test_mesh_geom_reconstructs_the_authored_placement(self) -> None:
        """Test that a mesh geom's points match the authored vertices under the authored
        pose.

        The compiler re-centers and re-orients a mesh's stored vertices and folds the
        offset into the geom's position and orientation, so recovering the authored
        placement proves that get_render_geometry composes the stored vertices with the
        compiled pose correctly. The fixture poses the tetrahedron with a 90 degree
        rotation about its parent's y axis followed by a translation of (1.0, 2.0, 3.0),
        and the mesh geom is attached to the body, so the expected points are in the
        first Airplane's body axes, relative to the first Airplane's CG.
        """
        rotate_T_act = _transformations.generate_rot_T(
            angles=np.array([0.0, 90.0, 0.0]),
            passive=False,
            intrinsic=False,
            order="xyz",
        )
        translate_T_act = _transformations.generate_trans_T(
            translations=np.array([1.0, 2.0, 3.0]), passive=False
        )
        pose_T_act = _transformations.compose_T_act(rotate_T_act, translate_T_act)
        stackExpectedPoints_BP1_CgP1 = _transformations.apply_T_to_vectors(
            pose_T_act,
            mujoco_model_fixtures.make_tetrahedron_vertices_fixture(),
            is_position=True,
        )

        stackRebuiltPoints_BP1_CgP1 = np.asarray(
            self.render_geoms[6].mesh.points, dtype=float
        )

        # The compiler may reorder the vertices, so match each authored vertex to its
        # nearest rebuilt point rather than comparing row by row. The authored vertices
        # are pairwise distinct, so requiring every one of them to sit within tolerance
        # of a rebuilt point, with the counts equal, proves the two sets match.
        self.assertEqual(
            len(stackRebuiltPoints_BP1_CgP1), len(stackExpectedPoints_BP1_CgP1)
        )
        for expectedPoint_BP1_CgP1 in stackExpectedPoints_BP1_CgP1:
            distances = np.linalg.norm(
                stackRebuiltPoints_BP1_CgP1 - expectedPoint_BP1_CgP1, axis=1
            )
            self.assertLess(float(distances.min()), 1e-6)

    def test_skips_infinite_plane_with_warning(self) -> None:
        """Test that an infinite plane is skipped with a warning naming the geom and
        stating that it still participates in the dynamics."""
        model = _mujoco_model.MuJoCoModel(
            name="infinite_plane_airplane",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=np.eye(3, dtype=float),
            delta_time=0.01,
            extra_xml={"worldbody": '<geom name="endless" type="plane" size="0 0 1"/>'},
        )

        with self.assertLogs("pterasoftware._mujoco_model", level="WARNING") as logs:
            render_geoms = model.get_render_geometry()

        self.assertEqual(render_geoms, [])
        self.assertEqual(len(logs.output), 1)
        self.assertIn("'endless'", logs.output[0])
        self.assertIn("an infinite plane", logs.output[0])
        self.assertIn("still participates in the dynamics", logs.output[0])

    def test_skips_heightfield_with_warning(self) -> None:
        """Test that a heightfield is skipped with a warning naming the geom and stating
        that it still participates in the dynamics."""
        model = _mujoco_model.MuJoCoModel(
            name="heightfield_airplane",
            mass=1.0,
            omegas_BP1__E=np.array((0.0, 0.0, 0.0)),
            T_pas_BP1_CgP1_to_E_CgP1=np.eye(4, dtype=float),
            vCg_E__E=np.array((10.0, 0.0, 0.0)),
            I_BP1_CgP1=np.eye(3, dtype=float),
            delta_time=0.01,
            extra_xml={
                "asset": (
                    '<asset><hfield name="terrain" nrow="2" ncol="2" '
                    'size="1 1 0.1 0.01"/></asset>'
                ),
                "worldbody": '<geom name="rough" type="hfield" hfield="terrain"/>',
            },
        )

        with self.assertLogs("pterasoftware._mujoco_model", level="WARNING") as logs:
            render_geoms = model.get_render_geometry()

        self.assertEqual(render_geoms, [])
        self.assertEqual(len(logs.output), 1)
        self.assertIn("'rough'", logs.output[0])
        self.assertIn("a heightfield", logs.output[0])
        self.assertIn("still participates in the dynamics", logs.output[0])


class TestMuJoCoModelUncoveredFileReferences(unittest.TestCase):
    """This class contains methods for testing MuJoCoModel.uncovered_file_references."""

    def test_no_file_references_returns_empty_list(self) -> None:
        """Test that a model whose XML has no file references returns an empty list."""
        model = mujoco_model_fixtures.make_basic_mujoco_model_fixture()
        self.assertEqual(model.uncovered_file_references(), [])

    def test_covered_reference_returns_empty_list(self) -> None:
        """Test that a mesh reference covered by the assets dict returns an empty
        list."""
        model = mujoco_model_fixtures.make_render_geometry_mujoco_model_fixture()
        self.assertEqual(model.uncovered_file_references(), [])

    def test_disk_reference_is_returned(self) -> None:
        """Test that a mesh referenced by an absolute on-disk path is returned."""
        model = mujoco_model_fixtures.make_disk_mesh_mujoco_model_fixture()
        uncovered = model.uncovered_file_references()
        self.assertEqual(len(uncovered), 1)
        self.assertTrue(uncovered[0].endswith("tetrahedron.stl"))
