# Vortex Singularities Branch Log

This file tracks progress and plans for the `feature/vortex_singularities` branch.

## Biot-Savart Kernels in `_aerodynamics_functions.py`

### Module Level Constants

The module defines three constants used by the Biot-Savart kernels:

- `_squire` (10^-4): Squire's parameter, which relates to the size and growth rate of vortex cores. Cited for use in flapping wing vehicles from Ananthan and Leishman (2004). It is unitless.
- `_lamb` (1.25643): Lamb's constant, which also relates to core size and growth. From Nguyen et al. (2016). It is unitless.
- `_eps` (machine epsilon): Used to guard against removable discontinuities when a point lies on or very near a LineVortex.

### Core Kernel Functions (Private)

These two functions contain the actual Biot-Savart computation. Both are `@njit` compiled with `cache=True` and `fastmath=False`.

#### `_collapsed_velocities_from_line_vortices`

Takes N points and M LineVortices. Returns a (N, 3) ndarray: the cumulative induced velocity at each point summed over all M LineVortices.

#### `_expanded_velocities_from_line_vortices`

Takes N points and M LineVortices. Returns a (N, M, 3) ndarray: the induced velocity at each point due to each individual LineVortex (not summed).

#### Shared Kernel Logic

Both functions implement the same modified Biot-Savart law (adapted from Nguyen et al., 2016). For each LineVortex with a given point:

1. Compute the vortex core radius: `r_c = 2 * sqrt(_lamb * (nu + _squire * |strength|) * age)`. For bound vortices (age = 0), this evaluates to 0.
2. Define geometric vectors:
   - `r0`: from the LineVortex start to end.
   - `r1`: from the evaluation point to the LineVortex start.
   - `r2`: from the evaluation point to the LineVortex end.
   - `r3`: cross product of `r1` and `r2`.
3. Check for removable discontinuities: if `|r1| < eps`, `|r2| < eps`, or `|r3|^2 < eps`, the point is on or essentially touching the LineVortex, so the induced velocity is zero (skip).
4. Otherwise, compute the induced velocity using the regularized kernel: `c_4 = (strength / (4 * pi)) * (|r1| + |r2|) * (|r1| * |r2| - r1 . r2) / (|r1| * |r2| * (|r3|^2 + |r0|^2 * r_c^2))`

The `|r0|^2 * r_c^2` term in the denominator is the regularization that smoothly desingularizes the velocity field near the vortex core.

### Public Wrapper Functions

These functions decompose RingVortices and HorseshoeVortices into their constituent LineVortex legs and call the appropriate core kernel.

#### Ring Vortex Wrappers

- `collapsed_velocities_from_ring_vortices`: Calls the collapsed kernel for all 4 legs (back right to front right, front right to front left, front left to back left, back left to back right). Returns (N, 3).
- `collapsed_velocities_from_ring_vortices_chordwise_segments`: Calls the collapsed kernel for the 2 chordwise legs only (back right to front right, front left to back left). Returns (N, 3). **Not currently used by any solver.**
- `expanded_velocities_from_ring_vortices`: Calls the expanded kernel for all 4 legs. Returns (N, M, 3).

#### Horseshoe Vortex Wrappers

- `collapsed_velocities_from_horseshoe_vortices`: Calls the collapsed kernel for 3 legs (back right to front right, front right to front left, front left to back left). Returns (N, 3).
- `expanded_velocities_from_horseshoe_vortices`: Calls the expanded kernel for 3 legs. Returns (N, M, 3).

### Parameters Common to All Functions

- `stackP_GP1_CgP1`: (N, 3) ndarray of evaluation point positions (in the first Airplane's geometry axes, relative to the first Airplane's CG).
- `strengths`: (M,) ndarray of vortex strengths in meters squared per second.
- `ages`: None for bound vortices, or (M,) ndarray of ages in seconds for wake vortices. Defaults to None.
- `nu`: Kinematic viscosity of the fluid in meters squared per second. Defaults to 0.0.

## Where the Solvers Call These Functions

### SteadyHorseshoeVortexLatticeMethodSolver

This solver uses only the horseshoe vortex wrapper functions.

| Aerodynamics function                          | Caller method                     | Purpose                                                                             |
|------------------------------------------------|-----------------------------------|-------------------------------------------------------------------------------------|
| `expanded_velocities_from_horseshoe_vortices`  | `_calculate_wing_wing_influences` | Build the (N, N) influence coefficient matrix with unit strength HorseshoeVortices. |
| `collapsed_velocities_from_horseshoe_vortices` | `calculate_solution_velocity`     | Compute the total induced velocity at arbitrary evaluation points.                  |

`calculate_solution_velocity` is itself called from:

- `_calculate_loads`: evaluates velocity at bound vortex centers for Kutta-Joukowski force computation.
- `calculate_streamlines` in `_functions.py` (via `run`): evaluates velocity at streamline seed points during streamline tracing.

### SteadyRingVortexLatticeMethodSolver

This solver uses both ring vortex and horseshoe vortex wrapper functions. The HorseshoeVortices model the semi infinite wake behind trailing edge Panels.

| Aerodynamics function                          | Caller method                     | Purpose                                                                      |
|------------------------------------------------|-----------------------------------|------------------------------------------------------------------------------|
| `expanded_velocities_from_ring_vortices`       | `_calculate_wing_wing_influences` | Build the ring vortex contribution to the influence coefficient matrix.      |
| `expanded_velocities_from_horseshoe_vortices`  | `_calculate_wing_wing_influences` | Build the horseshoe vortex contribution to the influence coefficient matrix. |
| `collapsed_velocities_from_ring_vortices`      | `calculate_solution_velocity`     | Compute the ring vortex induced velocity at evaluation points.               |
| `collapsed_velocities_from_horseshoe_vortices` | `calculate_solution_velocity`     | Compute the horseshoe vortex induced velocity at evaluation points.          |

`calculate_solution_velocity` is itself called from:

- `_calculate_loads`: evaluates velocity at each of the four LineVortex leg centers for all Panels' RingVortices (4 separate calls).
- `calculate_streamlines` in `_functions.py` (via `run`): streamline tracing.

### UnsteadyRingVortexLatticeMethodSolver

This solver uses only ring vortex wrapper functions. There are no HorseshoeVortices because the wake is modeled explicitly with shed wake RingVortices.

| Aerodynamics function                     | Caller method                     | Purpose                                                                               |
|-------------------------------------------|-----------------------------------|---------------------------------------------------------------------------------------|
| `expanded_velocities_from_ring_vortices`  | `_calculate_wing_wing_influences` | Build the influence coefficient matrix from bound RingVortices.                       |
| `collapsed_velocities_from_ring_vortices` | `_calculate_wake_wing_influences` | Compute the cumulative induced velocity at collocation points from wake RingVortices. |
| `collapsed_velocities_from_ring_vortices` | `calculate_solution_velocity`     | Compute bound RingVortex induced velocity at evaluation points.                       |
| `collapsed_velocities_from_ring_vortices` | `calculate_solution_velocity`     | Compute wake RingVortex induced velocity at evaluation points.                        |

`calculate_solution_velocity` is itself called from:

- `_calculate_loads`: evaluates velocity at each of the four LineVortex leg centers for all Panels' RingVortices (4 separate calls).
- `_populate_next_airplanes_wake_vortex_points`: evaluates velocity at wake RingVortex grid points during free wake convection (only when `prescribed_wake` is False).
- `calculate_streamlines` in `_functions.py` (via `run`): streamline tracing.

### Summary of Kernel Usage Across Solvers

| Kernel                                     | Steady Horseshoe VLM | Steady Ring VLM        | Unsteady Ring UVLM |
|--------------------------------------------|----------------------|------------------------|--------------------|
| `_expanded_velocities_from_line_vortices`  | Yes (horseshoe)      | Yes (ring + horseshoe) | Yes (ring)         |
| `_collapsed_velocities_from_line_vortices` | Yes (horseshoe)      | Yes (ring + horseshoe) | Yes (ring)         |

The expanded kernel is only used for building the Wing Wing influence coefficient matrices. The collapsed kernel is used everywhere else: solution velocity, load computation, wake influence, streamlines, and free wake convection.

`collapsed_velocities_from_ring_vortices_chordwise_segments` is not called by any solver.
