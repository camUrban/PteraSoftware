# MuJoCo Conventions

This document provides the definitive interpretation of MuJoCo state variables and their mapping to Ptera Software's conventions for axes, points, frames, angle vectors, and transformations. All interpretations have been empirically verified through systematic testing.

## MuJoCo Freejoint State Variables

For a body with a freejoint (6-DOF floating body), MuJoCo provides the following state variables:

### Position (`qpos[0:3]`)

**Interpretation:** Position of the body's origin in Earth axes.

```python
position_E_E = qpos[0:3]
```
No transformation needed.

### Quaternion (`qpos[3:7]`)

**Interpretation:** Orientation quaternion in scalar-first format `[w, x, y, z]`.

**Terminology note:** MuJoCo documentation may describe this as "the rotation from identity orientation to current body orientation." This phrasing describes the body's physical orientation (where its axes point in Earth axes) as if an active rotation were applied to move the body from identity to its current state. To create this quaternion from `xmat`, we do the following:

```python
R_pas_B_to_E = xmat.reshape(3, 3)

R_act_B_to_E = np.linalg.inv(R_pas_B_to_E)

R_act_E_to_B = R_act_B_to_E.T

quat_act_E_to_B_wxyz = _transformations.R_to_quat_wxyz(R_act_E_to_B)

qpos[3:7] = quat_act_E_to_B_wxyz
```

### Rotation Matrix (`xmat`)

**Interpretation:** `xmat` is `R_pas_B_to_E` (the passive rotation matrix that maps from body axes to Earth axes).

This means:
- `xmat` maps from body axes to Earth axes
- The columns of `xmat` are the body frame basis vectors expressed in Earth axes
  - Column 0: body +X direction in Earth axes
  - Column 1: body +Y direction in Earth axes
  - Column 2: body +Z direction in Earth axes

**Verification:** With a 90 degree rotation about the Earth Y axis, body +X points along Earth -Z. The test showed `xmat[:, 0] = [0, 0, -1]`, confirming that `xmat` columns represent body axes in Earth axes.

### Linear Velocity (`qvel[0:3]`)

**Interpretation:** Linear velocity of the body's center of mass in Earth axes.

```python
velocity_E__E = qvel[0:3]
```
No transformation needed. This is the velocity in Earth axes, observed from the Earth frame.

**Verification:** When force is applied along Earth +X to a rotated body, the velocity increases along Earth +X regardless of body orientation.

### Angular Velocity (`qvel[3:6]`)

**CRITICAL:** `qvel[3:6]` is in body axes, NOT Earth axes!

**Interpretation:** Angular velocity of the body expressed in body axes, in radians per second.

```python
omegas_B__E = np.rad2deg(qvel[3:6])
```
This is the angular velocity of body axes (expressed in body axes) observed from the Earth frame.

**To convert to Earth axes:**
```python
omegas_E__E = np.rad2deg(xmat @ qvel[3:6])
```
where `xmat = R_pas_B_to_E`.

**Verification:** With body rotated 90 deg about Y (so body +X = Earth -Z), setting `qvel[3:6] = [1, 0, 0]` caused rotation about Earth -Z (not Earth +X), proving that `qvel[3:6]` is in body axes.

### Applied Forces (`xfrc_applied[body_id][0:3]`)

**Interpretation:** External forces applied to the body at its center of mass, in Earth axes.

```python
xfrc_applied[body_id][0:3] = forces_E
```
No transformation needed if forces are already in Earth axes.

### Applied Torques (`xfrc_applied[body_id][3:6]`)

**Interpretation:** External torques applied to the body about its center of mass, in Earth axes.

```python
xfrc_applied[body_id][3:6] = moments_E_Cg
```
No transformation needed if moments are already in Earth axes.

**Verification:** Torque applied as `[10, 0, 0]` to a rotated body always produced angular acceleration about Earth +X, regardless of body orientation.

## Common Pitfalls

1. **Angular velocity axes:** The most common mistake is assuming `qvel[3:6]` is in Earth axes. It is in body axes.

2. **Rotation matrix direction:** `xmat` maps from body axes to Earth axes, not from Earth axes to body axes. Use `xmat.T` for the matrix that maps from Earth axes to body axes.

3. **Torque axes vs angular velocity axes:** Even though torques (`xfrc_applied[3:6]`) are in Earth axes, the resulting angular velocity (`qvel[3:6]`) is reported in body axes.

4. **Units:** `qvel[3:6]` is in radians per second, while Ptera Software typically uses degrees per second for angular quantities.