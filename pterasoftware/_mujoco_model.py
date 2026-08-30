"""Contains the MuJoCoModel class."""

from __future__ import annotations

from typing import NamedTuple, TypedDict
from xml.etree import ElementTree

import mujoco
import numpy as np
import pyvista as pv

from . import _logging, _transformations

_logger = _logging.get_logger("_mujoco_model")


class MuJoCoState(TypedDict):
    """The state returned by MuJoCoModel.get_state."""

    position_E_Eo: np.ndarray
    R_pas_E_to_BP1: np.ndarray
    velocity_E__E: np.ndarray
    omegas_BP1__E: np.ndarray
    time: float


class RenderGeom(NamedTuple):
    """One geom's renderable form, extracted from the compiled MuJoCo model by
    MuJoCoModel.get_render_geometry.

    mesh holds the geom's surface with its points posed in the axes of the geom's
    parent: the first Airplane's body axes (relative to the first Airplane's CG) when
    body_attached is True, and Earth axes (relative to the Earth origin) when it is
    False. rgba holds the geom's display color as a (4,) ndarray of floats in the range
    [0.0, 1.0], taken from the geom's material when one is assigned and from the geom's
    own rgba otherwise.
    """

    mesh: pv.PolyData
    rgba: np.ndarray
    body_attached: bool


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

    get_render_geometry: Extracts the renderable geometry of every geom in the compiled
    model.

    uncovered_file_references: Returns the file references in the model XML that the
    assets dict does not cover.

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
            save it alongside the XML string and lets _rebuild_engine resolve the XML
            string's asset references after a load. Validated by
            FreeFlightUnsteadyProblem before being passed here. The default is None,
            which provides no extra assets.
        :return: None
        """
        start_key_frame_name: str = "start"

        omegasRad_BP1__E = np.deg2rad(omegas_BP1__E)

        omegaXRad_BP1__E, omegaYRad_BP1__E, omegaZRad_BP1__E = omegasRad_BP1__E[:]

        R_pas_BP1_to_E = _transformations.extract_R_from_T(T_pas_BP1_CgP1_to_E_CgP1)

        R_act_E_to_BP1 = _transformations.invert_R_act(R_pas_BP1_to_E)

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

        # Create the internal MuJoCo model object from the XML str. If mujoco_assets is
        # provided, pass it so MuJoCo can resolve virtual filenames (e.g., STL meshes)
        # without writing them to disk.
        # noinspection PyArgumentList
        if mujoco_assets is not None:
            self._model: mujoco.MjModel = mujoco.MjModel.from_xml_string(
                self._xml_str, assets=mujoco_assets
            )
        else:
            # noinspection PyArgumentList
            self._model = mujoco.MjModel.from_xml_string(self._xml_str)

        # Retain the assets dict so the serialization layer can save it alongside the
        # XML string, which lets _rebuild_engine resolve the XML string's asset
        # references after a load without touching the filesystem.
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

        # Run forward kinematics to compute derived quantities (xmat, xpos, etc.) from
        # the initial qpos/qvel. Without this, xmat would be zeros until the first call
        # to mj_step.
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

        R_pas_E_to_BP1 = _transformations.invert_R_pas(R_pas_BP1_to_E)

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

    def get_render_geometry(self) -> list[RenderGeom]:
        """Extracts the renderable geometry of every geom in the compiled model.

        **Notes:**

        The generated model XML contains no geoms of its own, so every geom in the
        compiled model comes from the extra_xml fragments, and any mesh assets have
        already been decoded by MuJoCo's compiler into flat vertex and face arrays. The
        returned meshes are therefore rebuilt from the compiled model alone, without
        parsing XML or re-reading asset bytes.

        Each geom's points are posed in the axes of its parent. For geoms attached to
        the body, that is the first Airplane's body axes, relative to the first
        Airplane's CG. For geoms attached to the worldbody, that is Earth axes, relative
        to the Earth origin, since MuJoCo's world coordinates are defined to be
        identical to Earth axes.

        Infinite planes, heightfields, and signed distance field geoms have no finite
        surface to triangulate, so they are skipped with a logged warning. A skipped
        geom still participates in the dynamics: it is just not drawn.

        :return: A list of RenderGeoms, one per drawable geom, in the compiled model's
            geom order.
        """
        render_geoms: list[RenderGeom] = []
        for geom_id in range(self._model.ngeom):
            geom_mesh = self._get_local_geom_mesh(geom_id)

            if geom_mesh is None:
                geom_type = int(self._model.geom_type[geom_id])
                if geom_type == mujoco.mjtGeom.mjGEOM_PLANE:
                    shape_description = "an infinite plane"
                elif geom_type == mujoco.mjtGeom.mjGEOM_HFIELD:
                    shape_description = "a heightfield"
                elif geom_type == mujoco.mjtGeom.mjGEOM_SDF:
                    shape_description = "a signed distance field"
                else:
                    shape_description = "an unsupported geom type"

                # Describe the skipped geom by its name when it has one and by its ID
                # otherwise.
                geom_name = mujoco.mj_id2name(
                    self._model, mujoco.mjtObj.mjOBJ_GEOM, geom_id
                )
                geom_label = f"'{geom_name}'" if geom_name else f"with ID {geom_id}"
                _logger.warning(
                    _logging.indent()
                    + "Not drawing the MuJoCo geom %s because it is %s, which this "
                    "renderer does not support. The skipped shape still participates "
                    "in the dynamics: it is just not drawn.",
                    geom_label,
                    shape_description,
                )
                continue

            # The geom's stored orientation quaternion (scalar first) encodes
            # R_pas_geom_to_parent, the passive rotation from the geom's local axes to
            # the axes of its parent. MUJOCO_CONVENTIONS.md declares the informal geom,
            # geomOrigin, parent, and parentOrigin IDs used here. For a mesh geom, the
            # compiler has already folded the mesh's re-centering and re-orientation
            # into the stored pose, so applying it to the stored vertices reconstructs
            # the authored placement.
            R_pas_geom_to_parent: np.ndarray = np.zeros(9, dtype=float)
            mujoco.mju_quat2Mat(R_pas_geom_to_parent, self._model.geom_quat[geom_id])
            R_pas_geom_to_parent = R_pas_geom_to_parent.reshape(3, 3)

            # Lift the rotation into the transformation that maps from geom axes to
            # parent axes while holding the reference point at the geom origin.
            T_pas_geom_geomOrigin_to_parent_geomOrigin = np.eye(4, dtype=float)
            T_pas_geom_geomOrigin_to_parent_geomOrigin[:3, :3] = R_pas_geom_to_parent

            # The geom's stored position is geomOrigin_parent_parentOrigin, the position
            # of the geom origin (in parent axes, relative to the parent origin). For a
            # passive translation, the parameter is the position of the final reference
            # point (the parent origin) relative to the initial one (the geom origin),
            # which is its negative.
            geomOrigin_parent_parentOrigin = self._model.geom_pos[geom_id]
            T_pas_parent_geomOrigin_to_parent_parentOrigin = (
                _transformations.generate_trans_T(
                    translations=-geomOrigin_parent_parentOrigin, passive=True
                )
            )

            T_pas_geom_geomOrigin_to_parent_parentOrigin = (
                _transformations.compose_T_pas(
                    T_pas_geom_geomOrigin_to_parent_geomOrigin,
                    T_pas_parent_geomOrigin_to_parent_parentOrigin,
                )
            )
            geom_mesh.points = _transformations.apply_T_to_vectors(
                T_pas_geom_geomOrigin_to_parent_parentOrigin,
                geom_mesh.points,
                is_position=True,
            )

            # Drop any point normals the mesh's PyVista source attached. The rendering
            # layer recomputes shading normals from the posed triangles when it adds
            # each geom actor, so stored normals are never read, and ones left behind
            # here would sit stale against the posed points.
            if "Normals" in geom_mesh.point_data:
                del geom_mesh.point_data["Normals"]

            # Take the display color from the geom's material when one is assigned and
            # from the geom's own rgba otherwise.
            mat_id = int(self._model.geom_matid[geom_id])
            if mat_id >= 0:
                rgba = self._model.mat_rgba[mat_id].astype(float)
            else:
                rgba = self._model.geom_rgba[geom_id].astype(float)

            render_geoms.append(
                RenderGeom(
                    mesh=geom_mesh,
                    rgba=rgba,
                    body_attached=int(self._model.geom_bodyid[geom_id])
                    == self._body_id,
                )
            )

        return render_geoms

    def _get_local_geom_mesh(self, geom_id: int) -> pv.PolyData | None:
        """Builds one geom's surface mesh (in geom axes, relative to the geom origin),
        or returns None for a geom with no finite surface to triangulate.

        MuJoCo's geom_size semantics differ per geom type, so each supported type reads
        its own slots: a plane's first two sizes are the half-extents of its rendered
        rectangle (zero meaning infinite), a sphere's first size is its radius, a
        capsule's and a cylinder's are their radius and half-length (with their axis
        along the local z axis), an ellipsoid's are its three radii, and a box's are its
        three half-extents. A mesh geom ignores geom_size and reads the compiled vertex
        and face arrays instead.

        :param geom_id: The index of the geom in the compiled model's geom arrays.
        :return: The geom's surface as a PolyData with its points (in geom axes,
            relative to the geom origin), or None for an infinite plane, a heightfield,
            a signed distance field, or an unrecognized geom type.
        """
        geom_type = int(self._model.geom_type[geom_id])
        geom_size = self._model.geom_size[geom_id]

        if geom_type == mujoco.mjtGeom.mjGEOM_PLANE:
            # A zero half-extent marks the plane as infinite in that direction.
            if geom_size[0] <= 0.0 or geom_size[1] <= 0.0:
                return None
            return pv.Plane(
                direction=(0.0, 0.0, 1.0),
                i_size=2.0 * float(geom_size[0]),
                j_size=2.0 * float(geom_size[1]),
            )

        if geom_type == mujoco.mjtGeom.mjGEOM_SPHERE:
            return pv.Sphere(radius=float(geom_size[0]))

        if geom_type == mujoco.mjtGeom.mjGEOM_CAPSULE:
            return pv.Capsule(
                direction=(0.0, 0.0, 1.0),
                radius=float(geom_size[0]),
                cylinder_length=2.0 * float(geom_size[1]),
            )

        if geom_type == mujoco.mjtGeom.mjGEOM_ELLIPSOID:
            return pv.ParametricEllipsoid(
                xradius=float(geom_size[0]),
                yradius=float(geom_size[1]),
                zradius=float(geom_size[2]),
            )

        if geom_type == mujoco.mjtGeom.mjGEOM_CYLINDER:
            return pv.Cylinder(
                direction=(0.0, 0.0, 1.0),
                radius=float(geom_size[0]),
                height=2.0 * float(geom_size[1]),
            )

        if geom_type == mujoco.mjtGeom.mjGEOM_BOX:
            return pv.Box(
                bounds=(
                    -float(geom_size[0]),
                    float(geom_size[0]),
                    -float(geom_size[1]),
                    float(geom_size[1]),
                    -float(geom_size[2]),
                    float(geom_size[2]),
                )
            )

        if geom_type == mujoco.mjtGeom.mjGEOM_MESH:
            mesh_id = int(self._model.geom_dataid[geom_id])
            vert_start = int(self._model.mesh_vertadr[mesh_id])
            vert_count = int(self._model.mesh_vertnum[mesh_id])
            face_start = int(self._model.mesh_faceadr[mesh_id])
            face_count = int(self._model.mesh_facenum[mesh_id])
            vertices = self._model.mesh_vert[
                vert_start : vert_start + vert_count
            ].astype(float)

            # The compiled face array holds each triangle's three vertex indices, local
            # to the mesh's own vertex block. VTK's flat cell format prefixes each face
            # with its point count.
            triangles = self._model.mesh_face[
                face_start : face_start + face_count
            ].astype(int)
            faces = np.hstack((np.full((face_count, 1), 3, dtype=int), triangles))

            return pv.PolyData(vertices, faces.ravel())

        return None

    def uncovered_file_references(self) -> list[str]:
        """Returns the file references in the model XML that the assets dict does not
        cover.

        MuJoCo resolves a file attribute from the assets dict when the attribute's value
        matches a key, and falls back to the local filesystem otherwise. A reference
        resolved from the filesystem ties the model to a machine-specific path and could
        not survive a save and load round trip, so FreeFlightUnsteadyProblem rejects a
        model for which this method returns a nonempty list at construction time. Every
        file attribute in the XML is collected regardless of its element, which covers
        meshes, heightfields, textures, skins, and includes.

        :return: A list of strs representing the file references that the assets dict
            does not cover, in document order and without duplicates.
        """
        covered = self._mujoco_assets if self._mujoco_assets is not None else {}
        uncovered: list[str] = []
        for element in ElementTree.fromstring(self._xml_str).iter():
            file_reference = element.get("file")
            if file_reference is None:
                continue
            if file_reference in covered or file_reference in uncovered:
                continue
            uncovered.append(file_reference)
        return uncovered

    def _rebuild_engine(self) -> None:
        """Rebuilds the native MuJoCo model and data objects from the stored XML string
        and the stored assets dict.

        The serialization layer restores the serializable slots (the XML string, the
        assets dict, the body and keyframe ids, and the initial qpos and qvel) but
        cannot restore the live MuJoCo model and data objects, which wrap native MuJoCo
        state. This recreates them from the XML string and resets the data to the
        initial keyframe, matching the state established at construction. The restored
        assets dict is passed to MuJoCo so the rebuilt engine can resolve the XML
        string's asset references (such as meshes) without touching the filesystem,
        mirroring the construction path. The live, post-run data state is not restored:
        the canonical per-step state lives in the FreeFlightUnsteadyProblem's steady
        problems.

        :return: None
        """
        # noinspection PyArgumentList
        if self._mujoco_assets is not None:
            model = mujoco.MjModel.from_xml_string(
                self._xml_str, assets=self._mujoco_assets
            )
        else:
            # noinspection PyArgumentList
            model = mujoco.MjModel.from_xml_string(self._xml_str)
        data = mujoco.MjData(model)
        mujoco.mj_resetDataKeyframe(model, data, self._initial_key_frame_id)
        mujoco.mj_forward(model, data)
        self._model = model
        self._data = data
