# MuJoCo Conventions

This document gives the definitive interpretation of MuJoCo's state variables and compiled geometry constants, and how they map onto Ptera Software's conventions for axes, points, frames, angle vectors, and transformations. Every interpretation below was verified empirically through systematic testing during the development of free flight. The mapping is implemented by MuJoCoModel in `pterasoftware/_mujoco_model.py`, which is the one place the package interfaces with MuJoCo.

Free flight support is single-airplane for now, so this document describes the state of the one free-jointed body in terms of the first Airplane's body axes (BP1). We define MuJoCo's world coordinates to be identical to Ptera Software's Earth axes (E), so no transformation is needed between MuJoCo's world and our Earth axes.

## MuJoCo Freejoint State Variables

A body with a freejoint is a six degree of freedom floating body. MuJoCo describes its configuration with a seven element generalized position array (`qpos`: three position components followed by a four component orientation quaternion) and its motion with a six element generalized velocity array (`qvel`: three linear components followed by three angular components). External loads are applied through the body's six element `xfrc_applied` entry (three force components followed by three torque components). The sections below interpret each block.

### Position (qpos\[0:3])

`qpos[0:3]` is the position of the first Airplane's CG (in Earth axes, relative to the Earth origin), in meters. No transformation is needed.

```python
position_E_Eo = qpos[0:3]
```

### Orientation Quaternion (qpos\[3:7])

`qpos[3:7]` is the body's orientation quaternion in scalar-first format, [w, x, y, z]. It encodes the active rotation from Earth axes to the first Airplane's body axes (the active rotation that carries a frame from the Earth orientation to the body orientation). MuJoCoModel sets it from the initial passive orientation matrix R_pas_BP1_to_E (the rotation block of T_pas_BP1_CgP1_to_E_CgP1):

```python
R_pas_BP1_to_E = T_pas_BP1_CgP1_to_E_CgP1[:3, :3]
R_act_BP1_to_E = np.linalg.inv(R_pas_BP1_to_E)
R_act_E_to_BP1 = R_act_BP1_to_E.T
quat_act_E_to_BP1_wxyz = _transformations.R_to_quat_wxyz(R_act_E_to_BP1)
qpos[3:7] = quat_act_E_to_BP1_wxyz
```

### Orientation Matrix (xmat)

MuJoCo also exposes the body orientation as a 3x3 rotation matrix, `xmat`. This is R_pas_BP1_to_E, the passive rotation matrix that maps from the first Airplane's body axes to Earth axes. Equivalently, the columns of `xmat` are the body axes' basis directions expressed in Earth axes: column 0 is the body's +x direction in Earth axes, column 1 is the body's +y direction, and column 2 is the body's +z direction. To obtain R_pas_E_to_BP1, take the transpose.

```python
R_pas_BP1_to_E = data.xmat[body_id].reshape(3, 3)
R_pas_E_to_BP1 = R_pas_BP1_to_E.T
```

Verification: with the body rotated 90 degrees about the Earth y axis, the body's +x direction points along Earth -z. The test showed `xmat[:, 0] = [0, 0, -1]`, confirming that the columns of `xmat` are the body axes expressed in Earth axes.

### Linear Velocity (qvel\[0:3])

`qvel[0:3]` is the linear velocity of the first Airplane's CG (in Earth axes, observed from the Earth frame), in meters per second. No transformation is needed.

```python
velocity_E__E = qvel[0:3]
```

Verification: when a force along Earth +x is applied to a rotated body, the velocity grows along Earth +x regardless of the body's orientation.

### Angular Velocity (qvel\[3:6])

This is the most error-prone block. `qvel[3:6]` is in body axes, not Earth axes. It is the angular velocity of the first Airplane's body axes (in the first Airplane's body axes, observed from the Earth frame), in radians per second. Ptera Software works in degrees per second, so MuJoCoModel converts:

```python
omegas_BP1__E = np.rad2deg(qvel[3:6])
```

To express the same angular velocity in Earth axes instead, left-multiply by `xmat` (which is R_pas_BP1_to_E):

```python
omegas_E__E = np.rad2deg(xmat @ qvel[3:6])
```

Verification: with the body rotated 90 degrees about the y axis (so the body's +x direction equals Earth -z), setting `qvel[3:6] = [1, 0, 0]` produced rotation about Earth -z, not Earth +x, proving that `qvel[3:6]` is expressed in body axes.

### Applied Forces (xfrc_applied\[body_id]\[0:3])

`xfrc_applied[body_id][0:3]` is the external force applied to the body at its CG (in Earth axes), in Newtons. No transformation is needed when the force is already in Earth axes.

```python
xfrc_applied[body_id][0:3] = forces_E
```

### Applied Torques (xfrc_applied\[body_id]\[3:6])

`xfrc_applied[body_id][3:6]` is the external torque applied to the body about its CG (in Earth axes, relative to the first Airplane's CG), in Newton meters. No transformation is needed when the moment is already in Earth axes.

```python
xfrc_applied[body_id][3:6] = moments_E_CgP1
```

Verification: a torque of [10, 0, 0] applied to a rotated body always produced angular acceleration about Earth +x, regardless of the body's orientation.

## MuJoCo Geom and Mesh Arrays

The generated model XML contains no geoms of its own, so every geom in the compiled model comes from `extra_xml`. MuJoCoModel.get_render_geometry rebuilds their renderable surfaces from the compiled model's constant arrays alone, using the interpretations below, rather than parsing XML fragments or asset bytes.

### Informal Axes and Point IDs (geom, geomOrigin, parent, and parentOrigin)

The geom arrays involve two systems that are MuJoCo-internal, so their variables take informal lowercase IDs rather than IDs registered in AXES_POINTS_AND_FRAMES.md. The axes ID geom denotes a geom's own local axes, and the point ID geomOrigin denotes the geom origin, the point where MuJoCo anchors the geom's local coordinates. The axes ID parent denotes the parent body's axes, and the point ID parentOrigin denotes the parent origin, the point where MuJoCo anchors the parent body's coordinates. For a geom attached to the free-jointed body, parent axes are the first Airplane's body axes and the parent origin is the first Airplane's CG, since the generated body is anchored at the CG. For a worldbody geom, parent axes are Earth axes and the parent origin is the Earth origin. A single variable over the geom arrays spans both cases, which is why no registered ID can be pinned on it.

### Geom Ownership (geom_bodyid)

`geom_bodyid` gives the index of the body a geom is attached to. The worldbody is body 0, and the generated model contains exactly one other body, so a geom is attached to the free-jointed body exactly when its `geom_bodyid` equals that body's ID. Worldbody geoms are static in Earth axes (MuJoCo's world coordinates are Earth axes), while body geoms ride with the first Airplane's body axes.

### Geom Pose (geom_pos and geom_quat)

`geom_pos` and `geom_quat` give the geom's position and orientation relative to its parent body. `geom_pos` is geomOrigin_parent_parentOrigin, the position of the geom origin (in parent axes, relative to the parent origin). `geom_quat` is scalar first, [w, x, y, z], and `mju_quat2Mat` converts it to R_pas_geom_to_parent, the passive rotation matrix which maps from geom axes to parent axes. Together they build the passive transformation which maps in homogeneous coordinates from geom axes relative to the geom origin to parent axes relative to the parent origin:

```python
R_pas_geom_to_parent = np.zeros(9, dtype=float)
mujoco.mju_quat2Mat(R_pas_geom_to_parent, model.geom_quat[geom_id])
R_pas_geom_to_parent = R_pas_geom_to_parent.reshape(3, 3)

T_pas_geom_geomOrigin_to_parent_geomOrigin = np.eye(4, dtype=float)
T_pas_geom_geomOrigin_to_parent_geomOrigin[:3, :3] = R_pas_geom_to_parent

geomOrigin_parent_parentOrigin = model.geom_pos[geom_id]
T_pas_parent_geomOrigin_to_parent_parentOrigin = _transformations.generate_trans_T(
    translations=-geomOrigin_parent_parentOrigin, passive=True
)

T_pas_geom_geomOrigin_to_parent_parentOrigin = _transformations.compose_T_pas(
    T_pas_geom_geomOrigin_to_parent_geomOrigin,
    T_pas_parent_geomOrigin_to_parent_parentOrigin,
)
stackPoints_parent_parentOrigin = _transformations.apply_T_to_vectors(
    T_pas_geom_geomOrigin_to_parent_parentOrigin,
    stackPoints_geom_geomOrigin,
    is_position=True,
)
```

Verification: a geom authored with euler="0 90 0" produced a matrix mapping local +z to parent +x, which is the 90 degree rotation about y applied to geom-local coordinates.

### Geom Sizes (geom_size)

`geom_size` holds three slots whose meaning depends on `geom_type`:

- plane: the half-extents of the rendered rectangle in the local x and y directions, then the spacing of its rendered grid lines. A zero half-extent marks the plane as infinite in that direction.
- sphere: the radius, then two unused slots.
- capsule: the radius, then the half-length of the cylindrical section (the axis is the local z axis), then an unused slot.
- ellipsoid: the three radii.
- cylinder: the radius, then the half-length (the axis is the local z axis), then an unused slot.
- box: the three half-extents.
- mesh: derived bounding information that the rendering path does not use.

### Mesh Data (geom_dataid, mesh_vert, mesh_face, and their address arrays)

`geom_dataid` gives a mesh geom's index into the mesh arrays, and is -1 for every other geom type. `mesh_vertadr` and `mesh_vertnum` locate a mesh's block in `mesh_vert`, and `mesh_faceadr` and `mesh_facenum` locate its block in `mesh_face`. Each face's three vertex indices are local to the mesh's own vertex block, not global into `mesh_vert`. The compiler re-centers and re-orients each mesh's stored vertices and folds the offset into the owning geom's `geom_pos` and `geom_quat`, so `mesh_vert` does not hold the authored vertices, and composing the stored vertices with the compiled geom pose is what reconstructs the authored placement.

Verification: a tetrahedron authored with known vertices and posed with pos="1 2 3" and euler="0 90 0" compiled to re-centered vertices and a different geom pose, and applying that pose to the stored vertices reproduced the authored vertices under the authored pose. In a model with two meshes, the second mesh's face indices began at 0, confirming that face indices are local.

### Geom Colors (geom_rgba, geom_matid, and mat_rgba)

`geom_rgba` holds the geom's own display color and defaults to (0.5, 0.5, 0.5, 1.0). `geom_matid` is the index of the geom's material in `mat_rgba`, or -1 when no material is assigned. When a material is assigned, its rgba is the geom's display color.

## Common Pitfalls

1. Angular velocity axes: the most common mistake is assuming `qvel[3:6]` is in Earth axes. It is in the body axes.
2. Orientation matrix direction: `xmat` maps from body axes to Earth axes (R_pas_BP1_to_E), not the other way around. Use its transpose for R_pas_E_to_BP1.
3. Torque axes versus angular velocity axes: even though the applied torque (`xfrc_applied[3:6]`) is in Earth axes, the resulting angular velocity (`qvel[3:6]`) is reported in body axes.
4. Quaternion ordering: MuJoCo uses scalar-first quaternions, [w, x, y, z]. Libraries that use scalar-last quaternions, [x, y, z, w], will silently misinterpret the orientation.
5. Angular velocity units: `qvel[3:6]` is in radians per second, while Ptera Software uses degrees per second for angular quantities.
