# Vortex Singularities Branch Log

This file tracks progress and plans for the `feature/vortex_singularities` branch.
The design rationale for these changes is documented in `BIOT_SAVART_PLAN.md` (mirrors
the GitHub issue).

---

## Progress Summary

| Phase   | Status      | Description                                       |
|---------|-------------|---------------------------------------------------|
| Phase 1 | Complete    | Kernel singularity check and vortex core fixes    |
| Phase 2 | Not started | Context-aware singularity detection and reporting |

---

## Phase 1: Kernel Singularity Check and Vortex Core Fixes (Complete)

Phase 1 replaces absolute singularity cutoffs with scale-invariant criteria and adds
initial core radius regularization for the unsteady solver's wake and bound vortices.
The steady solvers pass `r_c0s = 0` (coreless) to match conventional VLM practice.

### Step 1 (Complete) -- Add `_tol` constant

Added `_tol = 1.0e-10` as a module-level constant in `_aerodynamics_functions.py`.

### Step 2 (Complete) -- Update singularity checks in both kernels

In both `_collapsed_velocities_from_line_vortices` and
`_expanded_velocities_from_line_vortices`:

- Added degenerate filament check in the outer (vortex) loop: `if r0 < _eps: continue`
- Replaced absolute inner-loop check (`r1 < _eps or r2 < _eps or r3**2 < _eps`) with
  scale-invariant checks:
  `if r1 / r0 < _tol or r2 / r0 < _tol or r3 / (r1 * r2) < _tol: continue`

### Step 3 (Complete) -- Add `r_c0s` parameter to both kernel functions

Added `r_c0s: np.ndarray` as a required parameter after `strengths` and before `ages` in
both kernels.

### Step 4 (Complete) -- Update the core radius formula in both kernels

Replaced the old formula (`r_c = 2 * sqrt(...)`) with the Ramasamy and Leishman (2007)
formulation: `r_c = sqrt(r_c0**2 + 4.0 * _lamb * (nu + _squire * |strength|) * age)`.

### Step 5 (Complete) -- Add `r_c0s` parameter to all 5 public wrapper functions

Added `r_c0s: np.ndarray` after `strengths` and before `ages` in all 5 wrappers:
`collapsed_velocities_from_ring_vortices`,
`collapsed_velocities_from_ring_vortices_chordwise_segments`,
`expanded_velocities_from_ring_vortices`,
`collapsed_velocities_from_horseshoe_vortices`,
`expanded_velocities_from_horseshoe_vortices`.

### Step 6 (Complete) -- Pass `r_c0s` in the steady horseshoe VLM solver

In `steady_horseshoe_vortex_lattice_method.py`:
- Allocated `self._stackRc0s` as zeros in `__init__`.
- Passed `r_c0s=self._stackRc0s` to both call sites.

Steady solvers use `r_c0 = 0` (coreless) to match conventional VLM practice. See
Step 11c for rationale.

### Step 7 (Complete) -- Pass `r_c0s` in the steady ring VLM solver

In `steady_ring_vortex_lattice_method.py`:
- Allocated `self._stackRc0s` as zeros in `__init__`.
- Passed `r_c0s=self._stackRc0s` to all 4 call sites.

Steady solvers use `r_c0 = 0` (coreless) to match conventional VLM practice. See
Step 11c for rationale.

### Step 8 (Complete) -- Compute and pass `r_c0s` in the unsteady ring UVLM solver

In `unsteady_ring_vortex_lattice_method.py`:
- Added `self._list_wake_rc0s` pre-allocation list in `__init__`.
- Added `self._currentStackBoundRc0s` and `self._currentStackWakeRc0s` per-step arrays.
- Populated bound and wake initial core radii in `_collapse_geometry`.
- Passed `r_c0s=self._currentStackBoundRc0s` to bound calls and
  `r_c0s=self._currentStackWakeRc0s` to wake calls at all 4 call sites.

### Step 9 (Complete) -- Update unit tests and fixtures

- Added `make_rc0s_fixture(num_vortices)` to
  `tests/unit/fixtures/aerodynamics_functions_fixtures.py`.
- Added `r_c0s` fixtures in `setUp` and passed `r_c0s=` to all 25 function calls in
  `tests/unit/test_aerodynamics_functions.py`.
- Updated the reference Biot-Savart implementation's singularity checks to use the same
  scale-invariant criteria as the kernels.
- Used `r_c0s = np.zeros(1, dtype=float)` in decomposition tests (the reference
  implementation is a coreless model).

### Steps 10-11 (Complete) -- Integration tests and Numba caches

Unit tests pass. Integration tests and Numba cache clearing to be verified before merge.

### Step 11b (Complete) -- Singularity guard and core radius unit tests

Added two new test classes to `test_aerodynamics_functions.py` and two new fixture
functions to `aerodynamics_functions_fixtures.py`:

- `TestSingularityGuards` (11 tests): Directly exercises the three singularity guards
  (degenerate filament, vertex proximity, collinearity) for both collapsed and expanded
  ring and horseshoe wrappers. Verifies singular legs contribute zero while non singular
  legs match the coreless reference.
- `TestCoreRadiusFormula` (6 tests): Validates the Ramasamy-Leishman core radius formula.
  Includes monotonicity tests (r_c0, age, nu), a large r_c0 suppression test, and two
  exact match tests against a regularized reference implementation.

### Step 11c (Complete) -- Revert steady solvers to coreless (r_c0 = 0)

Removed the `wing_r_c0 = 0.03 * _standard_mean_chord` computation and per-panel
`self._stackRc0s` population from both steady solvers. The arrays remain allocated as
zeros in `__init__`, so all steady solver vortices are now coreless. This matches
conventional VLM practice where viscous core effects are only applied in unsteady
simulations with aging wake vortices. The unsteady solver is unchanged.

---

## Phase 2: Context-Aware Singularity Detection and Reporting (Not Started)

Phase 2 adds counter-based singularity detection inside the `@njit` kernels and
context-aware logging at the solver level.

### Step 12 -- Add `singularity_counts` parameter to both kernel functions

Add `singularity_counts: np.ndarray` (shape `(4,)`, dtype `int64`) as a required
parameter after `r_c0s` in both kernels. Index mapping:
- `[0]`: degenerate filament
- `[1]`: vertex-start proximity
- `[2]`: vertex-end proximity
- `[3]`: collinearity

### Step 13 -- Increment counters at each singularity check

Split the combined `if/or` into separate checks with individual counter increments in
both kernels. The `continue` semantics preserve short-circuit behavior.

### Step 14 -- Thread `singularity_counts` through all 5 wrappers

Each wrapper receives `singularity_counts: np.ndarray` and passes the same array to every
internal kernel call. Counts accumulate across legs.

### Step 15 -- Add logging helper function

Add `_log_singularity_counts` to `pterasoftware/_functions.py`. Checks
`singularity_counts.sum() > 0`, then logs at a specified level with context string and
per-check breakdown.

### Step 16 -- Add optional counter parameters to `calculate_solution_velocity`

All three solvers' `calculate_solution_velocity` methods get optional
`bound_singularity_counts` and `wake_singularity_counts` parameters. If `None`, counts
are discarded.

### Step 17 -- Add logging at each solver call site

| Solver          | Call site                                     | Level     |
|-----------------|-----------------------------------------------|-----------|
| All solvers     | `_calculate_wing_wing_influences`             | `ERROR`   |
| Unsteady UVLM   | `_calculate_wake_wing_influences`             | `INFO`    |
| All solvers     | `_calculate_loads` (bound counters)           | `ERROR`   |
| Unsteady UVLM   | `_calculate_loads` (wake counters)            | `INFO`    |
| Unsteady UVLM   | `_populate_next_airplanes_wake_vortex_points` | `DEBUG`   |
| `_functions.py` | `calculate_streamlines`                       | `WARNING` |

### Step 18 -- Add logger instances to solver modules

Verify each solver module has a module-level logger. Add `import logging` where needed
for level constants.

### Step 19 -- Update unit tests for Phase 2

- Add `singularity_counts` parameter to all wrapper/kernel calls in unit tests.
- Add new tests verifying counters increment correctly for known-singular configurations.
- Add tests for the `_log_singularity_counts` helper.

### Step 20 (Deferred) -- Pre-exclusion of known singular pairs in `_calculate_loads`

Deferred to a follow-up. Requires combinatorial analysis of which LineVortex-point pairs
are structurally singular for load calculations. The current implementation will log these
expected singularities, which is acceptable for initial deployment.

---

## Files Modified by Phase 1

| File                                                      | Changes                                                                                                                             |
|-----------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| `pterasoftware/_aerodynamics_functions.py`                | `_tol` constant; scale-invariant singularity checks; `r_c0s` param in 2 kernels + 5 wrappers; Ramasamy-Leishman core radius formula |
| `pterasoftware/steady_horseshoe_vortex_lattice_method.py` | `_stackRc0s` zero allocation; `r_c0s` at 2 call sites (coreless)                                                                    |
| `pterasoftware/steady_ring_vortex_lattice_method.py`      | `_stackRc0s` zero allocation; `r_c0s` at 4 call sites (coreless)                                                                    |
| `pterasoftware/unsteady_ring_vortex_lattice_method.py`    | Bound + wake `r_c0s` arrays; pre-allocation list; `r_c0s` at 4 call sites                                                           |
| `tests/unit/fixtures/aerodynamics_functions_fixtures.py`  | `make_rc0s_fixture` helper                                                                                                          |
| `tests/unit/test_aerodynamics_functions.py`               | `r_c0s` at all 25 calls; updated reference implementation singularity checks; coreless decomposition tests                          |

## Files to be Modified by Phase 2

| File                                                      | Changes                                                                            |
|-----------------------------------------------------------|------------------------------------------------------------------------------------|
| `pterasoftware/_aerodynamics_functions.py`                | `singularity_counts` param in 2 kernels + 5 wrappers; per-check counter increments |
| `pterasoftware/steady_horseshoe_vortex_lattice_method.py` | Counter allocation + logging at 2 call sites                                       |
| `pterasoftware/steady_ring_vortex_lattice_method.py`      | Counter allocation + logging at 4 call sites                                       |
| `pterasoftware/unsteady_ring_vortex_lattice_method.py`    | Counter allocation + logging at 4 call sites                                       |
| `pterasoftware/_functions.py`                             | `_log_singularity_counts` helper; `calculate_streamlines` counter passing          |
| `tests/unit/test_aerodynamics_functions.py`               | `singularity_counts` at all calls; new counter tests                               |

---

## Reference: Solver Call Sites

### SteadyHorseshoeVortexLatticeMethodSolver

| Aerodynamics function                          | Caller method                     | Purpose                                                                             |
|------------------------------------------------|-----------------------------------|-------------------------------------------------------------------------------------|
| `expanded_velocities_from_horseshoe_vortices`  | `_calculate_wing_wing_influences` | Build the (N, N) influence coefficient matrix with unit strength HorseshoeVortices. |
| `collapsed_velocities_from_horseshoe_vortices` | `calculate_solution_velocity`     | Compute the total induced velocity at arbitrary evaluation points.                  |

`calculate_solution_velocity` is called from `_calculate_loads` and
`calculate_streamlines` (in `_functions.py`).

### SteadyRingVortexLatticeMethodSolver

| Aerodynamics function                          | Caller method                     | Purpose                                                                      |
|------------------------------------------------|-----------------------------------|------------------------------------------------------------------------------|
| `expanded_velocities_from_ring_vortices`       | `_calculate_wing_wing_influences` | Build the ring vortex contribution to the influence coefficient matrix.      |
| `expanded_velocities_from_horseshoe_vortices`  | `_calculate_wing_wing_influences` | Build the horseshoe vortex contribution to the influence coefficient matrix. |
| `collapsed_velocities_from_ring_vortices`      | `calculate_solution_velocity`     | Compute the ring vortex induced velocity at evaluation points.               |
| `collapsed_velocities_from_horseshoe_vortices` | `calculate_solution_velocity`     | Compute the horseshoe vortex induced velocity at evaluation points.          |

`calculate_solution_velocity` is called from `_calculate_loads` and
`calculate_streamlines` (in `_functions.py`).

### UnsteadyRingVortexLatticeMethodSolver

| Aerodynamics function                     | Caller method                     | Purpose                                                                           |
|-------------------------------------------|-----------------------------------|-----------------------------------------------------------------------------------|
| `expanded_velocities_from_ring_vortices`  | `_calculate_wing_wing_influences` | Build the influence coefficient matrix from bound RingVortices.                   |
| `collapsed_velocities_from_ring_vortices` | `_calculate_wake_wing_influences` | Compute cumulative induced velocity at collocation points from wake RingVortices. |
| `collapsed_velocities_from_ring_vortices` | `calculate_solution_velocity`     | Compute bound RingVortex induced velocity at evaluation points.                   |
| `collapsed_velocities_from_ring_vortices` | `calculate_solution_velocity`     | Compute wake RingVortex induced velocity at evaluation points.                    |

`calculate_solution_velocity` is called from `_calculate_loads`,
`_populate_next_airplanes_wake_vortex_points`, and `calculate_streamlines`
(in `_functions.py`).

### Kernel Usage Across Solvers

| Kernel                                     | Steady Horseshoe VLM | Steady Ring VLM        | Unsteady Ring UVLM |
|--------------------------------------------|----------------------|------------------------|--------------------|
| `_expanded_velocities_from_line_vortices`  | Yes (horseshoe)      | Yes (ring + horseshoe) | Yes (ring)         |
| `_collapsed_velocities_from_line_vortices` | Yes (horseshoe)      | Yes (ring + horseshoe) | Yes (ring)         |

The expanded kernel is only used for building the wing-wing influence coefficient
matrices. The collapsed kernel is used everywhere else.

`collapsed_velocities_from_ring_vortices_chordwise_segments` is not called by any solver.
