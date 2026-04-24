# Classes and Immutability

This document describes the consistent pattern of immutability and lazy caching across the following core data and geometry classes in the Ptera Software codebase:

- `CoreUnsteadyProblem` / `UnsteadyProblem`
- `_CoupledUnsteadyProblem`
- `CoreMovement` / `Movement`
- `CoreAirplaneMovement` / `AirplaneMovement`
- `CoreWingMovement` / `WingMovement`
- `CoreWingCrossSectionMovement` / `WingCrossSectionMovement`
- `CoreOperatingPointMovement` / `OperatingPointMovement`
- `SteadyProblem`
- `OperatingPoint`
- `Airplane`
- `Wing`
- `WingCrossSection`
- `Airfoil`
- `Panel`
- `RingVortex`
- `HorseshoeVortex`
- `LineVortex`

The `Core*` classes live in `pterasoftware/_core.py` and own the shared slots and properties. The public classes extend their core parents, sometimes adding additional slots and sometimes inheriting everything with an empty `__slots__`. See each class section below for details on which attributes are defined at which level.

## Design Principles

### Class Attribute Categories

Most attribute falls into one of these categories:

| Category                | Pattern                                                          |
|-------------------------|------------------------------------------------------------------|
| **Immutable**           | Read-only property (no setter), set once in `__init__`           |
| **Derived (Immutable)** | Manual lazy caching (check `None`, compute, cache, return)       |
| **Set Once**            | Property with setter that raises `AttributeError` if already set |
| **Mutable**             | Property with setter, or plain attribute                         |
| **Derived (Set Once)**  | Manual lazy caching, depends on set once attributes              |

### Key Decisions

1. **No cache invalidation for immutable/set once attributes**: Since these are only set once, we don't need invalidation logic in setters.

2. **Enforce set once semantics at runtime**: Set once properties raise `AttributeError` if assigned a second time. This catches bugs early where code incorrectly attempts to modify values that should be immutable after initial assignment.

3. **Use manual lazy caching for all derived properties**: This approach:
    - Works consistently for properties derived from both immutable and set once attributes
    - Is compatible with `__slots__` (no dependency on `__dict__`)
    - Simplifies `__deepcopy__` (cache variables can be copied directly)
    - Maintains the existing pattern already used throughout the codebase

4. **Solver mutable attributes remain mutable**: Properties that the solver needs to set keep their setters.

5. **`__slots__` on every class**: All classes in the package define `__slots__`, which eliminates per-instance `__dict__` overhead and prevents accidental dynamic attribute assignment. This catches typos like `self.num_panles = 5` at runtime with an `AttributeError` instead of silently creating a new attribute.

### NumPy Array Mutability

Even with read-only properties, numpy arrays can still be mutated in place via the getter (e.g., `panel.Frpp_G_Cg[0] = 999.0`). To prevent this, all numpy arrays that should be immutable are set to read-only using `arr.flags.writeable = False`:

1. **Immutable arrays**: Set in `__init__` immediately after assignment
2. **Set once arrays**: Set in the setter after assignment
3. **Derived cached arrays**: Set in the lazy property after computation (since numpy operations like subtraction create new writable arrays regardless of input writability)
4. **Deepcopy**: Use `.copy()` then set `flags.writeable = False` on the copy

### Deepcopy Cache Handling

When implementing `__deepcopy__`, handle cached derived properties based on their source:

1. **Derived from Immutable -> Preserve**: Copy the cached values (they remain valid since the source immutable attributes are also copied). For numpy arrays, use `.copy()` then set `flags.writeable = False`.
2. **Derived from Set Once -> Reset to None**: These depend on values that will be set fresh by the solver or meshing, so reset them.

### Serialization Attribute Handling

The `_serialization` module uses the same `object.__new__()` + `object.__setattr__()` pattern as the `__deepcopy__` methods but with a simpler strategy: the logical state of each object is preserved exactly as it was at save time, including all cached values on primary slots. For these primary attributes, no values are reset to `None` during deserialization.

Some redundant or alias slots (for example, solver alias slots and `Movement._airplanes` / `Movement._operating_points`) are treated specially: they are serialized as `null` and then deterministically rebuilt from their canonical sources during deserialization. This preserves object graph identity and avoids duplicating equivalent objects while still yielding the same effective state as the original instance.

This works because:

1. Deserialized objects are fully formed snapshots. There is no subsequent `__init__`, solver run, or meshing step that would populate set once attributes, so there is no conflict with set once guards.
2. Cached derived values remain valid because the immutable and set once attributes they depend on are also restored with their original values, and any alias slots are rebuilt to match.
3. NumPy array writeable flags are preserved through serialization, maintaining the same mutability guarantees as the original objects.

When adding or renaming `__slots__` on any class, both `__deepcopy__` and `_serialization` are affected. The serialization module discovers attributes generically via `__slots__`, so new slots are automatically serialized (subject to the intentional omission of redundant/alias slots described above). However, adding or removing slots requires incrementing `_FORMAT_VERSION` in `_serialization.py` to ensure old files are not loaded with incompatible code.

### List Collection Immutability

Store collections as tuples internally to prevent external mutation via `.append()`, `.pop()`, etc.

---

## CoreUnsteadyProblem / UnsteadyProblem Class (`_core.py`, `problems.py`)

`UnsteadyProblem` extends `CoreUnsteadyProblem`. `CoreUnsteadyProblem` owns all attributes except `movement` and `steady_problems`, which are defined on `UnsteadyProblem`.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type                        | Defined On              | Notes                 |
|------------------------|-----------------------------|-------------------------|-----------------------|
| `movement`             | `Movement`                  | `UnsteadyProblem`       | Movement definition   |
| `only_final_results`   | `bool`                      | `CoreUnsteadyProblem`   | Results flag          |
| `num_steps`            | `int`                       | `CoreUnsteadyProblem`   | Copied from Movement  |
| `delta_time`           | `float`                     | `CoreUnsteadyProblem`   | Copied from Movement  |
| `first_averaging_step` | `int`                       | `CoreUnsteadyProblem`   | Computed during init  |
| `first_results_step`   | `int`                       | `CoreUnsteadyProblem`   | Computed during init  |
| `max_wake_rows`        | `int \| None`               | `CoreUnsteadyProblem`   | Copied from Movement  |
| `steady_problems`      | `tuple[SteadyProblem, ...]` | `UnsteadyProblem`       | Generated during init |

#### Mutable (populated by solver)

| Attribute                            | Type               | Notes                              |
|--------------------------------------|--------------------|------------------------------------|
| `finalForces_W`                      | `list[np.ndarray]` | Final forces                       |
| `finalForceCoefficients_W`           | `list[np.ndarray]` | Final force coefficients           |
| `finalMoments_W_CgP1`                | `list[np.ndarray]` | Final moments                      |
| `finalMomentCoefficients_W_CgP1`     | `list[np.ndarray]` | Final moment coefficients          |
| `finalMeanForces_W`                  | `list[np.ndarray]` | Cycle averaged forces              |
| `finalMeanForceCoefficients_W`       | `list[np.ndarray]` | Cycle averaged coefficients        |
| `finalMeanMoments_W_CgP1`            | `list[np.ndarray]` | Cycle averaged moments             |
| `finalMeanMomentCoefficients_W_CgP1` | `list[np.ndarray]` | Cycle averaged moment coefficients |
| `finalRmsForces_W`                   | `list[np.ndarray]` | RMS forces                         |
| `finalRmsForceCoefficients_W`        | `list[np.ndarray]` | RMS force coefficients             |
| `finalRmsMoments_W_CgP1`             | `list[np.ndarray]` | RMS moments                        |
| `finalRmsMomentCoefficients_W_CgP1`  | `list[np.ndarray]` | RMS moment coefficients            |

**Note**: The mutable solver result lists are defined on `CoreUnsteadyProblem` and must remain mutable as they are populated after initialization by the solver. These are initialized as empty lists and appended to during the solve.

## _CoupledUnsteadyProblem Class (`problems.py`)

`_CoupledUnsteadyProblem` is a private middle-layer class that extends `CoreUnsteadyProblem`. It is the base for concrete subclasses (forthcoming `AeroelasticUnsteadyProblem` and `FreeFlightUnsteadyProblem`) whose geometry at each time step depends on the solver's results from the previous step. Unlike `UnsteadyProblem`, which builds all `SteadyProblem`s up front from a pre-generated `Movement`, the coupled subclasses grow their `SteadyProblem` collection one step at a time during the solve.

All `CoreUnsteadyProblem` attributes (documented in the section above) are inherited unchanged. The additions are:

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute         | Type                        | Notes                                                                                                                                                |
|-------------------|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| `movement`        | `CoreMovement`              | Source of `delta_time`, `num_steps`, `max_wake_rows`, and `lcm_period`                                                                               |
| `steady_problems` | `tuple[SteadyProblem, ...]` | Read-only view of the `_steady_problems` backing list; returned tuple is frozen, but successive calls may return different-length tuples (see below) |

**Note on `steady_problems`**: The parent class's `steady_problems` property is doubly immutable. The returned tuple is read-only and its value never changes over the lifetime of the `UnsteadyProblem`. On `_CoupledUnsteadyProblem`, the first guarantee still holds (callers cannot mutate the tuple), but the second does not. The backing slot `_steady_problems` is a `list[SteadyProblem]` seeded at init with a single entry built from `initial_airplanes` and `initial_operating_point`. Subclass `initialize_next_problem` overrides append to this list as each step is initialized during the solve, so calling `steady_problems` at different points can yield different-length tuples. External code that needs a consistent snapshot should read `steady_problems` once after the solver has completed.

## CoreMovement / Movement Class (`_core.py`, `movements/movement.py`)

`Movement` extends `CoreMovement`. `CoreMovement` owns the shared slots (`airplane_movements`, `operating_point_movement`, `delta_time`, `num_steps`, `max_wake_rows`) and derived properties (`lcm_period`, `max_period`, `min_period`, `static`). `Movement` adds cycle/chord counting, wake sizing parameters, and batch pre-generation of `Airplane`s and `OperatingPoint`s.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                  | Type                               | Defined On     | Notes                     |
|----------------------------|------------------------------------|----------------|---------------------------|
| `airplane_movements`       | `tuple[AirplaneMovement, ...]`     | `CoreMovement` | Tuple prevents mutation   |
| `operating_point_movement` | `OperatingPointMovement`           | `CoreMovement` | Operating point changes   |
| `delta_time`               | `float`                            | `CoreMovement` | Time step                 |
| `num_steps`                | `int`                              | `CoreMovement` | Total time steps          |
| `max_wake_rows`            | `int \| None`                      | `CoreMovement` | Max wake rows per Wing    |
| `num_cycles`               | `int \| None`                      | `Movement`     | Number of cycles          |
| `num_chords`               | `int \| None`                      | `Movement`     | Number of chord lengths   |
| `max_wake_chords`          | `int \| None`                      | `Movement`     | Max wake in chord lengths |
| `max_wake_cycles`          | `int \| None`                      | `Movement`     | Max wake in motion cycles |
| `airplanes`                | `tuple[tuple[Airplane, ...], ...]` | `Movement`     | Generated Airplanes       |
| `operating_points`         | `tuple[OperatingPoint, ...]`       | `Movement`     | Generated OperatingPoints |

#### Derived from Immutable (use manual lazy caching)

| Property     | Depends On                                       | Defined On     | Notes  |
|--------------|--------------------------------------------------|----------------|--------|
| `lcm_period` | `airplane_movements`, `operating_point_movement` | `CoreMovement` | Cached |
| `max_period` | `airplane_movements`, `operating_point_movement` | `CoreMovement` | Cached |
| `min_period` | `airplane_movements`, `operating_point_movement` | `CoreMovement` | Cached |
| `static`     | `max_period`                                     | `CoreMovement` | Cached |

**Note on `airplanes` and `operating_points`**: These are defined on `Movement` and generated during `__init__` by calling the child movements' `generate_*` methods. They are stored as nested tuples to prevent modification after generation.

## CoreAirplaneMovement / AirplaneMovement Class (`_core.py`, `movements/airplane_movement.py`)

`AirplaneMovement` extends `CoreAirplaneMovement`. All slots are defined on `CoreAirplaneMovement`; `AirplaneMovement` has empty `__slots__` and only narrows the `wing_movements` type to require `WingMovement` children.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute            | Type                       | Notes                   |
|----------------------|----------------------------|-------------------------|
| `base_airplane`      | `Airplane`                 | Base geometry           |
| `wing_movements`     | `tuple[WingMovement, ...]` | Tuple prevents mutation |
| `ampCg_GP1_CgP1`     | `np.ndarray`               | CG position amplitudes  |
| `periodCg_GP1_CgP1`  | `np.ndarray`               | CG position periods     |
| `spacingCg_GP1_CgP1` | `tuple`                    | CG position spacing     |
| `phaseCg_GP1_CgP1`   | `np.ndarray`               | CG position phases      |

#### Derived from Immutable (use manual lazy caching)

| Property      | Depends On                        | Notes                                     |
|---------------|-----------------------------------|-------------------------------------------|
| `all_periods` | Own periods + child `all_periods` | Tuple of unique non zero periods (cached) |
| `max_period`  | Own periods + child `max_period`  | Scalar float, longest period (cached)     |

## CoreWingMovement / WingMovement Class (`_core.py`, `movements/wing_movement.py`)

`WingMovement` extends `CoreWingMovement`. All slots are defined on `CoreWingMovement`; `WingMovement` has empty `__slots__` and only narrows the `wing_cross_section_movements` type to require `WingCrossSectionMovement` children.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                      | Type                                   | Notes                   |
|--------------------------------|----------------------------------------|-------------------------|
| `base_wing`                    | `Wing`                                 | Base geometry           |
| `wing_cross_section_movements` | `tuple[WingCrossSectionMovement, ...]` | Tuple prevents mutation |
| `ampLer_Gs_Cgs`                | `np.ndarray`                           | Position amplitudes     |
| `periodLer_Gs_Cgs`             | `np.ndarray`                           | Position periods        |
| `spacingLer_Gs_Cgs`            | `tuple`                                | Position spacing        |
| `phaseLer_Gs_Cgs`              | `np.ndarray`                           | Position phases         |
| `ampAngles_Gs_to_Wn_ixyz`      | `np.ndarray`                           | Angle amplitudes        |
| `periodAngles_Gs_to_Wn_ixyz`   | `np.ndarray`                           | Angle periods           |
| `spacingAngles_Gs_to_Wn_ixyz`  | `tuple`                                | Angle spacing           |
| `phaseAngles_Gs_to_Wn_ixyz`    | `np.ndarray`                           | Angle phases            |
| `rotationPointOffset_Gs_Ler`   | `np.ndarray`                           | Rotation point offset   |

#### Derived from Immutable (use manual lazy caching)

| Property      | Depends On                        | Notes                                     |
|---------------|-----------------------------------|-------------------------------------------|
| `all_periods` | Own periods + child `all_periods` | Tuple of unique non zero periods (cached) |
| `max_period`  | Own periods + child `max_period`  | Scalar float, longest period (cached)     |

## CoreWingCrossSectionMovement / WingCrossSectionMovement Class (`_core.py`, `movements/wing_cross_section_movement.py`)

`WingCrossSectionMovement` extends `CoreWingCrossSectionMovement`. All slots are defined on `CoreWingCrossSectionMovement`; `WingCrossSectionMovement` has empty `__slots__`.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                        | Type               | Notes               |
|----------------------------------|--------------------|---------------------|
| `base_wing_cross_section`        | `WingCrossSection` | Base geometry       |
| `ampLp_Wcsp_Lpp`                 | `np.ndarray`       | Position amplitudes |
| `periodLp_Wcsp_Lpp`              | `np.ndarray`       | Position periods    |
| `spacingLp_Wcsp_Lpp`             | `tuple`            | Position spacing    |
| `phaseLp_Wcsp_Lpp`               | `np.ndarray`       | Position phases     |
| `ampAngles_Wcsp_to_Wcs_ixyz`     | `np.ndarray`       | Angle amplitudes    |
| `periodAngles_Wcsp_to_Wcs_ixyz`  | `np.ndarray`       | Angle periods       |
| `spacingAngles_Wcsp_to_Wcs_ixyz` | `tuple`            | Angle spacing       |
| `phaseAngles_Wcsp_to_Wcs_ixyz`   | `np.ndarray`       | Angle phases        |

#### Derived from Immutable (use manual lazy caching)

| Property      | Depends On    | Notes                                     |
|---------------|---------------|-------------------------------------------|
| `all_periods` | Period arrays | Tuple of unique non zero periods (cached) |
| `max_period`  | Period arrays | Scalar float, longest period (cached)     |

## CoreOperatingPointMovement / OperatingPointMovement Class (`_core.py`, `movements/operating_point_movement.py`)

`OperatingPointMovement` extends `CoreOperatingPointMovement`. All slots are defined on `CoreOperatingPointMovement`; `OperatingPointMovement` has empty `__slots__`.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type              | Notes                     |
|------------------------|-------------------|---------------------------|
| `base_operating_point` | `OperatingPoint`  | Base operating conditions |
| `ampVCg__E`            | `float`           | Amplitude                 |
| `periodVCg__E`         | `float`           | Period                    |
| `spacingVCg__E`        | `str \| Callable` | Spacing function          |
| `phaseVCg__E`          | `float`           | Phase offset              |

#### Derived from Immutable (use manual lazy caching)

| Property     | Depends On     | Notes  |
|--------------|----------------|--------|
| `max_period` | `periodVCg__E` | Cached |

## SteadyProblem Class (`problems.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute         | Type                   | Notes                            |
|-------------------|------------------------|----------------------------------|
| `airplanes`       | `tuple[Airplane, ...]` | Tuple prevents external mutation |
| `operating_point` | `OperatingPoint`       | Operating conditions             |

#### Derived from Immutable (use manual lazy caching)

| Property           | Depends On                     | Notes                                  |
|--------------------|--------------------------------|----------------------------------------|
| `reynolds_numbers` | `airplanes`, `operating_point` | Tuple of Re for each Airplane (cached) |

## OperatingPoint Class (`operating_point.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type                 | Notes                     |
|------------------------|----------------------|---------------------------|
| `rho`                  | `float`              | Fluid density             |
| `vCg__E`               | `float`              | CG speed                  |
| `alpha`                | `float`              | Angle of attack           |
| `beta`                 | `float`              | Sideslip angle            |
| `angles_E_to_BP1_izyx` | `np.ndarray`         | Earth-to-body orientation |
| `CgP1_E_Eo`            | `np.ndarray`         | CG position in Earth axes |
| `surfaceNormal_E`      | `np.ndarray \| None` | Image surface normal      |
| `surfacePoint_E_Eo`    | `np.ndarray \| None` | Image surface point       |
| `externalFX_W`         | `float`              | External force            |
| `nu`                   | `float`              | Kinematic viscosity       |

#### Derived from Immutable (use manual lazy caching)

| Property                        | Depends On                                                   | Notes                             |
|---------------------------------|--------------------------------------------------------------|-----------------------------------|
| `qInf__E`                       | `rho`, `vCg__E`                                              | Dynamic pressure (cached)         |
| `T_pas_GP1_CgP1_to_BP1_CgP1`    | (constant)                                                   | Geometry-to-body matrix (cached)  |
| `T_pas_BP1_CgP1_to_GP1_CgP1`    | `T_pas_GP1_CgP1_to_BP1_CgP1`                                 | Inverse of above (cached)         |
| `T_pas_BP1_CgP1_to_W_CgP1`      | `alpha`, `beta`                                              | Body-to-wind matrix (cached)      |
| `T_pas_W_CgP1_to_BP1_CgP1`      | `T_pas_BP1_CgP1_to_W_CgP1`                                   | Inverse of above (cached)         |
| `T_pas_GP1_CgP1_to_W_CgP1`      | `T_pas_GP1_CgP1_to_BP1_CgP1`, `T_pas_BP1_CgP1_to_W_CgP1`     | Geometry-to-wind matrix (cached)  |
| `T_pas_W_CgP1_to_GP1_CgP1`      | `T_pas_GP1_CgP1_to_W_CgP1`                                   | Inverse of above (cached)         |
| `T_pas_E_CgP1_to_BP1_CgP1`      | `angles_E_to_BP1_izyx`                                       | Earth-to-body matrix (cached)     |
| `T_pas_BP1_CgP1_to_E_CgP1`      | `T_pas_E_CgP1_to_BP1_CgP1`                                   | Inverse of above (cached)         |
| `T_pas_E_CgP1_to_GP1_CgP1`      | `T_pas_E_CgP1_to_BP1_CgP1`, `T_pas_BP1_CgP1_to_GP1_CgP1`     | Earth-to-geometry matrix (cached) |
| `T_pas_GP1_CgP1_to_E_CgP1`      | `T_pas_E_CgP1_to_GP1_CgP1`                                   | Inverse of above (cached)         |
| `surfaceNormal_GP1`             | `surfaceNormal_E`, `T_pas_E_CgP1_to_GP1_CgP1`                | Surface normal in GP1 (cached)    |
| `surfacePoint_GP1_CgP1`         | `surfacePoint_E_Eo`, `CgP1_E_Eo`, `T_pas_E_CgP1_to_GP1_CgP1` | Surface point in GP1 (cached)     |
| `surfaceReflect_T_act_GP1_CgP1` | `surfacePoint_GP1_CgP1`, `surfaceNormal_GP1`                 | Active reflection matrix (cached) |
| `vInfHat_GP1__E`                | `T_pas_W_CgP1_to_GP1_CgP1`                                   | Freestream direction (cached)     |
| `vInf_GP1__E`                   | `vInfHat_GP1__E`, `vCg__E`                                   | Freestream velocity (cached)      |

**Note on caching**: While `vInfHat_GP1__E` and `vInf_GP1__E` are simple computations once the transformation matrix is cached, they are cached for consistency with the overall pattern and because they are called repeatedly during solver operations. The same is true for `qInf__E`.

**Note on `T_pas_GP1_CgP1_to_BP1_CgP1`**: This matrix is a constant (180-degree rotation about y) and does not depend on any `__init__` parameters. It is still lazily cached for consistency with the overall pattern and to avoid recomputing it on every access.

**Note on transformation decomposition**: The geometry-to-wind transformation (`T_pas_GP1_CgP1_to_W_CgP1`) is composed from `T_pas_GP1_CgP1_to_BP1_CgP1` and `T_pas_BP1_CgP1_to_W_CgP1` via `compose_T_pas`. Similarly, `T_pas_E_CgP1_to_GP1_CgP1` is composed from `T_pas_E_CgP1_to_BP1_CgP1` and `T_pas_BP1_CgP1_to_GP1_CgP1`. Decomposing through body axes enables the Earth transformation chain to reuse the geometry-to-body and body-to-geometry matrices without duplication.

## Airplane Class (`geometry/airplane.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute     | Type               | Notes                                                        |
|---------------|--------------------|--------------------------------------------------------------|
| `wings`       | `tuple[Wing, ...]` | Processed for symmetry during init (tuple prevents mutation) |
| `name`        | `str`              | Airplane identifier                                          |
| `Cg_GP1_CgP1` | `np.ndarray`       | CG position in formation coordinates                         |
| `weight`      | `float`            | Aircraft weight in Newtons                                   |
| `s_ref`       | `float`            | Reference wetted area                                        |
| `c_ref`       | `float`            | Reference chord length                                       |
| `b_ref`       | `float`            | Reference span                                               |

#### Derived from Immutable (use manual lazy caching)

| Property                 | Depends On    | Notes                    |
|--------------------------|---------------|--------------------------|
| `num_panels`             | `wings`       | Sum of wing panel counts |
| `T_pas_G_Cg_to_GP1_CgP1` | `Cg_GP1_CgP1` | Transformation matrix    |

#### Mutable (set by solver)

| Attribute                   | Type                 | Notes                |
|-----------------------------|----------------------|----------------------|
| `forces_W`                  | `np.ndarray \| None` | Forces in wind axes  |
| `forceCoefficients_W`       | `np.ndarray \| None` | Force coefficients   |
| `moments_W_CgP1`            | `np.ndarray \| None` | Moments in wind axes |
| `momentCoefficients_W_CgP1` | `np.ndarray \| None` | Moment coefficients  |

---

## Wing Class (`geometry/wing.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type                           | Notes                                         |
|------------------------|--------------------------------|-----------------------------------------------|
| `wing_cross_sections`  | `tuple[WingCrossSection, ...]` | Wing cross sections (tuple prevents mutation) |
| `name`                 | `str`                          | Wing identifier                               |
| `Ler_Gs_Cgs`           | `np.ndarray`                   | Leading edge root position                    |
| `angles_Gs_to_Wn_ixyz` | `np.ndarray`                   | Rotation angles                               |
| `num_chordwise_panels` | `int`                          | Chordwise panel count                         |
| `chordwise_spacing`    | `str`                          | "cosine" or "uniform"                         |

#### Derived from Immutable (use manual lazy caching)

| Property                  | Depends On      | Notes                  |
|---------------------------|-----------------|------------------------|
| `T_pas_G_Cg_to_Wn_Ler`    | Immutable attrs | Transformation matrix  |
| `T_pas_Wn_Ler_to_G_Cg`    | Above           | Inverse transformation |
| `WnX_G`, `WnY_G`, `WnZ_G` | Above           | Basis vectors          |
| `children_T_pas_*`        | Cross sections  | Child transformations  |

#### Set Once (set by `generate_mesh`, never modified after)

| Attribute             | Type                 | Set By          | Notes                |
|-----------------------|----------------------|-----------------|----------------------|
| `symmetry_type`       | `int \| None`        | `generate_mesh` | 1, 2, 3, or 4        |
| `num_spanwise_panels` | `int \| None`        | `generate_mesh` | Total spanwise count |
| `num_panels`          | `int \| None`        | `generate_mesh` | Total panel count    |
| `panels`              | `np.ndarray \| None` | `generate_mesh` | Panel matrix         |

#### Derived from Set Once (use manual lazy caching)

| Property                     | Depends On                           | Notes                |
|------------------------------|--------------------------------------|----------------------|
| `projected_area`             | `panels`                             | Projected area       |
| `wetted_area`                | `panels`                             | Wetted area          |
| `average_panel_aspect_ratio` | `panels`                             | Average aspect ratio |
| `span`                       | Wing cross sections, `symmetry_type` | Wing span            |
| `standard_mean_chord`        | `projected_area`, `span`             | Standard mean chord  |
| `mean_aerodynamic_chord`     | `projected_area`, `symmetry_type`    | MAC                  |

**Note on caching**: Most derived properties iterate over panels or wing cross sections. For large meshes, caching `projected_area`, `wetted_area`, and `span` provides meaningful performance gains if they're accessed multiple times. Since their source attributes are immutable or set once, these are cached after first computation without invalidation logic.

#### Mutable (modified by `process_wing_symmetry` for type 5 symmetry)

| Attribute            | Type                 | Notes                        |
|----------------------|----------------------|------------------------------|
| `symmetric`          | `bool`               | Modified to False for type 5 |
| `mirror_only`        | `bool`               | Modified to False for type 5 |
| `symmetryNormal_G`   | `np.ndarray \| None` | Modified to None for type 5  |
| `symmetryPoint_G_Cg` | `np.ndarray \| None` | Modified to None for type 5  |

**Note**: These are modified by `Airplane.process_wing_symmetry()` when type 5 symmetry is detected. The original Wing becomes a type 1 wing and a new reflected Wing is created with type 3 symmetry.

#### Mutable (modified during simulation for wake)

| Attribute            | Type                 | Notes                        |
|----------------------|----------------------|------------------------------|
| `wake_ring_vortices` | `np.ndarray \| None` | Wake vortex array, grows     |
| `gridWrvp_GP1_CgP1`  | `np.ndarray \| None` | Wake vortex positions, grows |

## WingCrossSection Class (`geometry/wing_cross_section.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                     | Type          | Notes                      |
|-------------------------------|---------------|----------------------------|
| `airfoil`                     | `Airfoil`     | Wing cross section airfoil |
| `num_spanwise_panels`         | `int \| None` | Spanwise panel count       |
| `chord`                       | `float`       | Chord length               |
| `Lp_Wcsp_Lpp`                 | `np.ndarray`  | Position in parent axes    |
| `angles_Wcsp_to_Wcs_ixyz`     | `np.ndarray`  | Rotation angles            |
| `control_surface_hinge_point` | `float`       | Hinge location (0-1)       |
| `control_surface_deflection`  | `float`       | Deflection in degrees      |
| `spanwise_spacing`            | `str \| None` | "cosine" or "uniform"      |

#### Derived from Immutable (use manual lazy caching)

| Property                   | Depends On                  | Notes                  |
|----------------------------|-----------------------------|------------------------|
| `T_pas_Wcsp_Lpp_to_Wcs_Lp` | `Lp_Wcsp_Lpp`, `angles_...` | Transformation matrix  |
| `T_pas_Wcs_Lp_to_Wcsp_Lpp` | Above                       | Inverse transformation |

#### Set Once (set by parent Wing)

| Attribute       | Type          | Set By               | Notes                   |
|-----------------|---------------|----------------------|-------------------------|
| `validated`     | `bool`        | `Wing.__init__`      | Validation flag         |
| `symmetry_type` | `int \| None` | `Wing.generate_mesh` | Inherited symmetry type |

#### Mutable (modified by `process_wing_symmetry`)

| Attribute                       | Type          | Notes                           |
|---------------------------------|---------------|---------------------------------|
| `control_surface_symmetry_type` | `str \| None` | Set to None for type 5 symmetry |

**Note**: Modified at when type 5 symmetry is split into two wings.

## Airfoil Class (`geometry/airfoil.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute           | Type         | Notes                          |
|---------------------|--------------|--------------------------------|
| `name`              | `str`        | Airfoil identifier             |
| `outline_A_lp`      | `np.ndarray` | Outline coordinates            |
| `resample`          | `bool`       | Resampling flag                |
| `n_points_per_side` | `int`        | Points per side for resampling |
| `mcl_A_lp`          | `np.ndarray` | Mean camber line coordinates   |

**Note**: The `add_control_surface` method creates and returns a new Airfoil instance rather than modifying the existing one. This is the correct immutable pattern.

## Panel Class (`_panel.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute          | Type         | Notes                              |
|--------------------|--------------|------------------------------------|
| `Frpp_G_Cg`        | `np.ndarray` | Front right corner (geometry axes) |
| `Flpp_G_Cg`        | `np.ndarray` | Front left corner (geometry axes)  |
| `Blpp_G_Cg`        | `np.ndarray` | Back left corner (geometry axes)   |
| `Brpp_G_Cg`        | `np.ndarray` | Back right corner (geometry axes)  |
| `is_leading_edge`  | `bool`       | Edge flag                          |
| `is_trailing_edge` | `bool`       | Edge flag                          |

#### Derived from Immutable (use manual lazy caching)

| Property       | Depends On                | Notes                          |
|----------------|---------------------------|--------------------------------|
| `rightLeg_G`   | `Frpp_G_Cg`, `Brpp_G_Cg`  | Right leg vector               |
| `frontLeg_G`   | `Flpp_G_Cg`, `Frpp_G_Cg`  | Front leg vector               |
| `leftLeg_G`    | `Blpp_G_Cg`, `Flpp_G_Cg`  | Left leg vector                |
| `backLeg_G`    | `Brpp_G_Cg`, `Blpp_G_Cg`  | Back leg vector                |
| `Frbvp_G_Cg`   | `Brpp_G_Cg`, `rightLeg_G` | Front right bound vortex point |
| `Flbvp_G_Cg`   | `Flpp_G_Cg`, `leftLeg_G`  | Front left bound vortex point  |
| `Cpp_G_Cg`     | Multiple corners          | Collocation point              |
| `unitNormal_G` | All corners               | Unit normal vector             |
| `area`         | All corners               | Panel area                     |
| `aspect_ratio` | All legs                  | Aspect ratio                   |

#### Set Once (set after construction by meshing or Problem)

| Attribute                  | Type                 | Set By  | Notes                        |
|----------------------------|----------------------|---------|------------------------------|
| `is_right_edge`            | `bool \| None`       | Meshing | Edge flag                    |
| `is_left_edge`             | `bool \| None`       | Meshing | Edge flag                    |
| `local_chordwise_position` | `int \| None`        | Meshing | Grid position                |
| `local_spanwise_position`  | `int \| None`        | Meshing | Grid position                |
| `Frpp_GP1_CgP1`            | `np.ndarray \| None` | Problem | Front right (formation axes) |
| `Flpp_GP1_CgP1`            | `np.ndarray \| None` | Problem | Front left (formation axes)  |
| `Blpp_GP1_CgP1`            | `np.ndarray \| None` | Problem | Back left (formation axes)   |
| `Brpp_GP1_CgP1`            | `np.ndarray \| None` | Problem | Back right (formation axes)  |

#### Derived from Set Once (use manual lazy caching)

| Property         | Depends On                       | Notes                         |
|------------------|----------------------------------|-------------------------------|
| `rightLeg_GP1`   | `Frpp_GP1_CgP1`, `Brpp_GP1_CgP1` | Right leg (formation axes)    |
| `frontLeg_GP1`   | `Flpp_GP1_CgP1`, `Frpp_GP1_CgP1` | Front leg (formation axes)    |
| `leftLeg_GP1`    | `Blpp_GP1_CgP1`, `Flpp_GP1_CgP1` | Left leg (formation axes)     |
| `backLeg_GP1`    | `Brpp_GP1_CgP1`, `Blpp_GP1_CgP1` | Back leg (formation axes)     |
| `Frbvp_GP1_CgP1` | `Brpp_GP1_CgP1`, `rightLeg_GP1`  | Front right BVP (formation)   |
| `Flbvp_GP1_CgP1` | `Flpp_GP1_CgP1`, `leftLeg_GP1`   | Front left BVP (formation)    |
| `Cpp_GP1_CgP1`   | Multiple GP1 corners             | Collocation point (formation) |
| `unitNormal_GP1` | All GP1 corners                  | Unit normal (formation)       |

#### Mutable (set by solver)

| Attribute          | Type                      | Notes                |
|--------------------|---------------------------|----------------------|
| `ring_vortex`      | `RingVortex \| None`      | Attached vortex      |
| `horseshoe_vortex` | `HorseshoeVortex \| None` | Attached vortex      |
| `forces_GP1`       | `np.ndarray \| None`      | Computed forces      |
| `moments_GP1_CgP1` | `np.ndarray \| None`      | Computed moments     |
| `forces_W`         | `np.ndarray \| None`      | Forces in wind axes  |
| `moments_W_CgP1`   | `np.ndarray \| None`      | Moments in wind axes |

---

## RingVortex Class (`_vortices/ring_vortex.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute        | Type         | Notes              |
|------------------|--------------|--------------------|
| `Frrvp_GP1_CgP1` | `np.ndarray` | Front right corner |
| `Flrvp_GP1_CgP1` | `np.ndarray` | Front left corner  |
| `Blrvp_GP1_CgP1` | `np.ndarray` | Back left corner   |
| `Brrvp_GP1_CgP1` | `np.ndarray` | Back right corner  |

#### Derived from Immutable (use manual lazy caching)

| Property        | Depends On  | Notes             |
|-----------------|-------------|-------------------|
| `Crvp_GP1_CgP1` | All corners | Centroid position |
| `area`          | All corners | Vortex area       |

#### Mutable (solver sets these)

| Attribute  | Type    | Notes                                         |
|------------|---------|-----------------------------------------------|
| `strength` | `float` | Vortex strength (solver finds this)           |
| `age`      | `float` | Age in simulation time (incremented for wake) |

#### Derived (special: child objects)

| Property    | Depends On                                     | Notes              |
|-------------|------------------------------------------------|--------------------|
| `front_leg` | `Frrvp_GP1_CgP1`, `Flrvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `left_leg`  | `Flrvp_GP1_CgP1`, `Blrvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `back_leg`  | `Blrvp_GP1_CgP1`, `Brrvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `right_leg` | `Brrvp_GP1_CgP1`, `Frrvp_GP1_CgP1`, `strength` | Child `LineVortex` |

**Note on legs**: The leg `LineVortex` objects depend on both geometry (immutable) and `strength` (mutable). Since `strength` is set by the solver AFTER the vortex is created, we propagate strength updates to existing legs. Since legs are accessed repeatedly during induced velocity calculations, option 2 (keeping the current propagation) is more efficient.

---

## HorseshoeVortex Class (`_vortices/horseshoe_vortex.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                | Type         | Notes                              |
|--------------------------|--------------|------------------------------------|
| `Frhvp_GP1_CgP1`         | `np.ndarray` | Front right point                  |
| `Flhvp_GP1_CgP1`         | `np.ndarray` | Front left point                   |
| `leftLegVector_GP1`      | `np.ndarray` | Direction of left leg (normalized) |
| `left_right_leg_lengths` | `float`      | Length of semi-infinite legs       |

#### Derived from Immutable (use manual lazy caching)

| Property         | Depends On                                                      | Notes            |
|------------------|-----------------------------------------------------------------|------------------|
| `Brhvp_GP1_CgP1` | `Frhvp_GP1_CgP1`, `leftLegVector_GP1`, `left_right_leg_lengths` | Back right point |
| `Blhvp_GP1_CgP1` | `Flhvp_GP1_CgP1`, `leftLegVector_GP1`, `left_right_leg_lengths` | Back left point  |

#### Mutable (solver sets these)

| Attribute  | Type    | Notes           |
|------------|---------|-----------------|
| `strength` | `float` | Vortex strength |

#### Derived (child objects)

| Property     | Depends On                                     | Notes              |
|--------------|------------------------------------------------|--------------------|
| `right_leg`  | `Brhvp_GP1_CgP1`, `Frhvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `finite_leg` | `Frhvp_GP1_CgP1`, `Flhvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `left_leg`   | `Flhvp_GP1_CgP1`, `Blhvp_GP1_CgP1`, `strength` | Child `LineVortex` |

---

## LineVortex Class (`_vortices/_line_vortex.py`)

### Attribute Classification

Since `LineVortex` is an internal class whose endpoints ARE updated by parent vortex classes when their corners change, we need to consider whether this mutation actually happens.

#### Immutable (set in `__init__`)

| Attribute       | Type         | Notes       |
|-----------------|--------------|-------------|
| `Slvp_GP1_CgP1` | `np.ndarray` | Start point |
| `Elvp_GP1_CgP1` | `np.ndarray` | End point   |

#### Derived from Immutable (use manual lazy caching)

| Property        | Depends On                       | Notes        |
|-----------------|----------------------------------|--------------|
| `vector_GP1`    | `Elvp_GP1_CgP1`, `Slvp_GP1_CgP1` | Line vector  |
| `Clvp_GP1_CgP1` | `Slvp_GP1_CgP1`, `vector_GP1`    | Center point |

#### Mutable

| Attribute  | Type    | Notes                    |
|------------|---------|--------------------------|
| `strength` | `float` | Updated by parent vortex |

---

## Solver Classes (Not Covered Above)

The four solver classes (`SteadyHorseshoeVortexLatticeMethodSolver`, `SteadyRingVortexLatticeMethodSolver`, `UnsteadyRingVortexLatticeMethodSolver`, and `CoupledUnsteadyRingVortexLatticeMethodSolver`) are intentionally omitted from the immutability and lazy caching patterns described in this document. Unlike the data and geometry classes above, the solver classes are algorithmic classes whose attributes are internal mutable working state in a procedural computation pipeline. They are not shared data that external code accesses or modifies, so immutable properties, set once enforcement, and lazy caching would add significant boilerplate with no meaningful safety benefit. The solver classes do still use `__slots__`, like all other classes in the package, to protect against dynamic attribute assignment typos.