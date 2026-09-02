# Classes and Immutability

This document describes the consistent pattern of immutability and lazy caching across the following core data and geometry classes in the Ptera Software codebase:

- `CoreUnsteadyProblem` / `UnsteadyProblem`
- `_CoupledUnsteadyProblem`
- `AeroelasticUnsteadyProblem`
- `FreeFlightUnsteadyProblem`
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
- `MuJoCoModel`

The `Core*` classes live in `pterasoftware/_core.py` and own the shared slots and properties. The public classes extend their core parents, sometimes adding additional slots and sometimes inheriting everything with an empty `__slots__`. See each class section below for details on which attributes are defined at which level.

## Design Principles

### Class Attribute Categories

Most attributes fall into one of these categories:

| Category                | Pattern                                                          |
|-------------------------|------------------------------------------------------------------|
| **Immutable**           | Read-only property (no setter), set once in `__init__`           |
| **Derived (Immutable)** | Manual lazy caching (check `None`, compute, cache, return)       |
| **Set Once**            | Property with setter that raises `AttributeError` if already set |
| **Mutable**             | Property with setter, or plain attribute                         |
| **Derived (Set Once)**  | Manual lazy caching, depends on set once attributes              |
| **Derived (Mutable)**   | Read-only property, no backing slot, recomputed on every access  |

### Construction-Only Parameters

A constructor parameter is not always retained as an attribute. Some parameters shape how an object is built but are deliberately discarded once construction finishes: they have no `__slots__` entry, appear in no attribute-category table, and are not serialized or deep copied, because nothing reads them after `__init__`. They are still validated exactly like their stored counterparts.

`Wing`'s `explode_into_strips` is one example. When True, it triggers `Wing._explode_wing` to replace the supplied wing cross sections with single-strip cross sections. The bool itself is never assigned to `self`: its validated value lives in a local variable inside `__init__`, passing through `boolLike_return_bool`, the same check applied to the stored bool-likes `symmetric` and `mirror_only`. Construction does, however, derive the immutable `spanwise_mesh` marker from it (see the Wing Class section), recording whether the spanwise mesh was defined as trapezoidal corners or as single-strip cross sections. That provenance is set explicitly here rather than detected later, precisely because an exploded Wing is structurally indistinguishable from one a user defines with single-strip cross sections directly: the standard convergence tools must refine the latter but reject the former, so the distinction cannot be recovered from the geometry alone.

### Key Decisions

1. **No cache invalidation for immutable/set once attributes**: Since these are only set once, we don't need invalidation logic in setters.

2. **Enforce set once semantics at runtime**: Set once properties raise `AttributeError` if assigned a second time. This catches bugs early where code incorrectly attempts to modify values that should be immutable after initial assignment.

3. **Use manual lazy caching for all cached derived properties**: This approach:
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

A saved file is a `.psz` file: a zip archive of compressed JSON members. The first member is a header holding the format version, provenance metadata, the saved object's class name, and a manifest of the members that follow. The last member holds the object itself. For an `UnsteadyProblem` or an unsteady solver, the per time step sequences (the problem's `_steady_problems` and the solver's per step bound and wake vortex lists) are split out of the last member into one member per time step, and the corresponding slots in the last member hold placeholders that `load` splices the step members back into. This keeps the serialized text of at most one time step in memory at once during a save or a load, on top of the object itself.

Shared references are preserved by an identity memo that spans every member of the file: the first encounter of an object serializes it in full, and every later encounter (for example, a solver's `airplanes` alias into its `SteadyProblem`, or `Movement._airplanes` inside an `UnsteadyProblem`) serializes as a reference to it. Deserialization resolves those references to the same rebuilt instance, so object graph identity survives the round trip without duplicating equivalent objects. The one exception is `MuJoCoModel`'s native engine objects, which cannot be serialized at all, so they are serialized as `null` and rebuilt from the model's XML string on deserialization.

Custom callables are the other thing a saved file cannot hold. The built in `"sine"` and `"uniform"` spacing functions serialize by name and are restored exactly. Every other callable (a custom spacing function, an `AeroelasticWingMovement`'s `spacingAnglesSecondDerivative_Gs_to_Wn_ixyz` entry, or a `FreeFlightUnsteadyProblem`'s `external_loads_fn`) serializes as an inert marker holding its qualified name, its dedented source text when `inspect.getsource` can retrieve it, and the SHA-256 hash of that text. Nothing stored in the file is ever executed. Deserialization restores the marker as an `UnboundCallable`, a private slotted placeholder in `_serialization.py` holding the marker's three fields as immutable attributes, which occupies the same slot, passes `callable(...)` checks, and re-emits the identical marker when saved again, so a loaded object's `hash_object` digest matches the original's. Because the source text is part of the serialized content, editing a custom function changes the digest, which is what lets the convergence cache self-invalidate. Loading never invokes these callables, so a solved simulation that holds them loads and visualizes normally, and an `UnboundCallable` raises only if something later tries to call it. To re-solve or regenerate motion, the user passes the original functions to `load`'s `callables` argument, keyed by the recorded qualified names, and `load` restores them in place of the placeholders during deserialization rather than mutating the loaded object afterwards, which would need an immutability bypass. A key that matches no marker raises, and a function whose source does not hash to a marker's recorded hash logs a warning, since the file was saved with a different definition.

This works because:

1. Deserialized objects are fully formed snapshots. There is no subsequent `__init__`, solver run, or meshing step that would populate set once attributes, so there is no conflict with set once guards.
2. Cached derived values remain valid because the immutable and set once attributes they depend on are also restored with their original values, and alias slots resolve to the same restored objects as their canonical sources.
3. NumPy array writeable flags are preserved through serialization, maintaining the same mutability guarantees as the original objects.

When adding or renaming `__slots__` on any class, both `__deepcopy__` and `_serialization` are affected. The serialization module discovers attributes generically via `__slots__`, so new slots are automatically serialized (subject to the `MuJoCoModel` exception described above). However, adding or removing slots requires incrementing `_FORMAT_VERSION` in `_serialization.py` to ensure old files are not loaded with incompatible code.

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

| Attribute                            | Type               | Notes                                            |
|--------------------------------------|--------------------|--------------------------------------------------|
| `finalForces_W`                      | `list[np.ndarray]` | Final forces                                     |
| `finalForceCoefficients_W`           | `list[np.ndarray]` | Final force coefficients                         |
| `finalMoments_W_CgP1`                | `list[np.ndarray]` | Final moments                                    |
| `finalMomentCoefficients_W_CgP1`     | `list[np.ndarray]` | Final moment coefficients                        |
| `finalForces_G`                      | `list[np.ndarray]` | Final geometry-axes forces                       |
| `finalForceCoefficients_G`           | `list[np.ndarray]` | Final geometry-axes force coefficients           |
| `finalMoments_G_Cg`                  | `list[np.ndarray]` | Final own-CG moments (geometry axes)             |
| `finalMomentCoefficients_G_Cg`       | `list[np.ndarray]` | Final own-CG moment coefficients (geometry axes) |
| `finalMoments_W_Cg`                  | `list[np.ndarray]` | Final own-CG moments (wind axes)                 |
| `finalMomentCoefficients_W_Cg`       | `list[np.ndarray]` | Final own-CG moment coefficients (wind axes)     |
| `finalMeanForces_W`                  | `list[np.ndarray]` | Cycle averaged counterparts of the above         |
| `finalMeanForceCoefficients_W`       | `list[np.ndarray]` |                                                  |
| `finalMeanMoments_W_CgP1`            | `list[np.ndarray]` |                                                  |
| `finalMeanMomentCoefficients_W_CgP1` | `list[np.ndarray]` |                                                  |
| `finalMeanForces_G`                  | `list[np.ndarray]` |                                                  |
| `finalMeanForceCoefficients_G`       | `list[np.ndarray]` |                                                  |
| `finalMeanMoments_G_Cg`              | `list[np.ndarray]` |                                                  |
| `finalMeanMomentCoefficients_G_Cg`   | `list[np.ndarray]` |                                                  |
| `finalMeanMoments_W_Cg`              | `list[np.ndarray]` |                                                  |
| `finalMeanMomentCoefficients_W_Cg`   | `list[np.ndarray]` |                                                  |
| `finalRmsForces_W`                   | `list[np.ndarray]` | RMS counterparts of the above                    |
| `finalRmsForceCoefficients_W`        | `list[np.ndarray]` |                                                  |
| `finalRmsMoments_W_CgP1`             | `list[np.ndarray]` |                                                  |
| `finalRmsMomentCoefficients_W_CgP1`  | `list[np.ndarray]` |                                                  |
| `finalRmsForces_G`                   | `list[np.ndarray]` |                                                  |
| `finalRmsForceCoefficients_G`        | `list[np.ndarray]` |                                                  |
| `finalRmsMoments_G_Cg`               | `list[np.ndarray]` |                                                  |
| `finalRmsMomentCoefficients_G_Cg`    | `list[np.ndarray]` |                                                  |
| `finalRmsMoments_W_Cg`               | `list[np.ndarray]` |                                                  |
| `finalRmsMomentCoefficients_W_Cg`    | `list[np.ndarray]` |                                                  |

**Note**: The mutable solver result lists are defined on `CoreUnsteadyProblem` and must remain mutable as they are populated after initialization by the solver. These are initialized as empty lists and appended to during the solve. In these per-Airplane lists, the G and Cg IDs are used without an Airplane index to implicitly mean "in the entry's own Airplane's geometry axes" and "relative to the entry's own Airplane's CG".

#### Derived from Mutable (read-only property, no backing slot)

These read-only properties expose the named load components and coefficients defined in `AXES_POINTS_AND_FRAMES.md` for the wind axes lists. Each returns one signed component per Airplane, or an empty list while its source list is empty, and stores nothing of its own. They are recomputed on every access, since a cached value would go stale when the solver populates the source lists.

| Property                                   | Depends On                         | Notes                                    |
|--------------------------------------------|------------------------------------|------------------------------------------|
| `finalInducedDrags_W`                      | `finalForces_W`                    | Negatives of the x components            |
| `finalSideForces_W`                        | `finalForces_W`                    | The y components                         |
| `finalLifts_W`                             | `finalForces_W`                    | Negatives of the z components            |
| `finalInducedDragCoefficients_W`           | `finalForceCoefficients_W`         | Negatives of the x components            |
| `finalSideForceCoefficients_W`             | `finalForceCoefficients_W`         | The y components                         |
| `finalLiftCoefficients_W`                  | `finalForceCoefficients_W`         | Negatives of the z components            |
| `finalRollingMoments_W_Cg`                 | `finalMoments_W_Cg`                | The x components                         |
| `finalPitchingMoments_W_Cg`                | `finalMoments_W_Cg`                | The y components                         |
| `finalYawingMoments_W_Cg`                  | `finalMoments_W_Cg`                | The z components                         |
| `finalRollingMomentCoefficients_W_Cg`      | `finalMomentCoefficients_W_Cg`     | The x components                         |
| `finalPitchingMomentCoefficients_W_Cg`     | `finalMomentCoefficients_W_Cg`     | The y components                         |
| `finalYawingMomentCoefficients_W_Cg`       | `finalMomentCoefficients_W_Cg`     | The z components                         |
| `finalMeanInducedDrags_W`                  | `finalMeanForces_W`                | Cycle averaged counterparts of the above |
| `finalMeanSideForces_W`                    | `finalMeanForces_W`                |                                          |
| `finalMeanLifts_W`                         | `finalMeanForces_W`                |                                          |
| `finalMeanInducedDragCoefficients_W`       | `finalMeanForceCoefficients_W`     |                                          |
| `finalMeanSideForceCoefficients_W`         | `finalMeanForceCoefficients_W`     |                                          |
| `finalMeanLiftCoefficients_W`              | `finalMeanForceCoefficients_W`     |                                          |
| `finalMeanRollingMoments_W_Cg`             | `finalMeanMoments_W_Cg`            |                                          |
| `finalMeanPitchingMoments_W_Cg`            | `finalMeanMoments_W_Cg`            |                                          |
| `finalMeanYawingMoments_W_Cg`              | `finalMeanMoments_W_Cg`            |                                          |
| `finalMeanRollingMomentCoefficients_W_Cg`  | `finalMeanMomentCoefficients_W_Cg` |                                          |
| `finalMeanPitchingMomentCoefficients_W_Cg` | `finalMeanMomentCoefficients_W_Cg` |                                          |
| `finalMeanYawingMomentCoefficients_W_Cg`   | `finalMeanMomentCoefficients_W_Cg` |                                          |

## _CoupledUnsteadyProblem Class (`problems.py`)

`_CoupledUnsteadyProblem` is a private middle-layer class that extends `CoreUnsteadyProblem`. It is the base for concrete subclasses (`AeroelasticUnsteadyProblem` and `FreeFlightUnsteadyProblem`, both documented below) whose per-step `SteadyProblem` depends on the solver's results from the previous step: deformed wing geometry for aeroelasticity, updated rigid body state for free flight. Unlike `UnsteadyProblem`, which builds all `SteadyProblem`s up front from a pre-generated `Movement`, the coupled subclasses grow their `SteadyProblem` collection one step at a time during the solve.

All `CoreUnsteadyProblem` attributes (documented in the section above) are inherited unchanged. The additions are:

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute         | Type                        | Notes                                                                                                                                                |
|-------------------|-----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| `movement`        | `CoreMovement`              | Source of `delta_time`, `num_steps`, `max_wake_rows`, and `lcm_period`                                                                               |
| `steady_problems` | `tuple[SteadyProblem, ...]` | Read-only view of the `_steady_problems` backing list; returned tuple is frozen, but successive calls may return different-length tuples (see below) |

**Note on `steady_problems`**: The parent class's `steady_problems` property is doubly immutable. The returned tuple is read-only and its value never changes over the lifetime of the `UnsteadyProblem`. On `_CoupledUnsteadyProblem`, the first guarantee still holds (callers cannot mutate the tuple), but the second does not. The backing slot `_steady_problems` is a `list[SteadyProblem]` seeded at init with a single entry built from `initial_airplanes` and `initial_operating_point`. Subclass `initialize_next_problem` overrides append to this list as each step is initialized during the solve, so calling `steady_problems` at different points can yield different-length tuples. External code that needs a consistent snapshot should read `steady_problems` once after the solver has completed.

## AeroelasticUnsteadyProblem Class (`problems.py`)

`AeroelasticUnsteadyProblem` extends `_CoupledUnsteadyProblem`. It couples aerodynamic loads with a torsional spring-mass-damper structural model so that each wing's deformation at a given time step is driven by the previous step's aerodynamic, inertial, and spring-restoring moments. All `_CoupledUnsteadyProblem` and `CoreUnsteadyProblem` attributes (documented in the sections above) are inherited unchanged. The additions are the structural configuration (set once at construction) and the per-wing structural state (populated as the solve advances).

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

Each is stored in a `_`-prefixed backing slot and exposed through a read-only property of the unprefixed name.

| Attribute              | Type    | Notes                                                  |
|------------------------|---------|--------------------------------------------------------|
| `wing_density`         | `float` | Mass per unit span area (kg/m^2)                       |
| `spring_constant_rad`  | `float` | Torsional spring stiffness (N*m/rad)                   |
| `damping_constant_rad` | `float` | Torsional damping coefficient (N*m*s/rad)              |
| `step_discards`        | `int`   | Number of initial steps discarded for stability        |

#### Derived from Immutable (read-only property, no backing slot)

| Property                | Depends On              | Notes                                                                                               |
|-------------------------|-------------------------|-----------------------------------------------------------------------------------------------------|
| `_aeroelastic_movement` | `_movement`             | Typed-narrow cast of the inherited `_movement` slot to `AeroelasticMovement` (recomputed, uncached) |

#### Mutable (populated by solver)

These lists are allocated in `__init__` with one entry per wing in the initial airplane, then appended to by the solver during the run. They are plain slots rather than read-only properties because the solver must update them after construction, mirroring the mutable-result-list treatment on `CoreUnsteadyProblem`.

| Attribute                                                | Type                     | Notes                                                                                 |
|----------------------------------------------------------|--------------------------|---------------------------------------------------------------------------------------|
| `listDeformationAnglesYRad_Wcsp_to_Wcs_ixyz`             | `list[list[np.ndarray]]` | Cumulative torsional angle time series (radians), indexed `[wing][entry]`, seeded     |
| `_listDeformationAnglesDerivativeYRad_Wcsp_to_Wcs_ixyz`  | `list[list[np.ndarray]]` | Torsional angle time derivative time series (rad/s), indexed `[wing][entry]`, seeded  |

## FreeFlightUnsteadyProblem Class (`problems.py`)

`FreeFlightUnsteadyProblem` extends `_CoupledUnsteadyProblem`. It couples aerodynamic loads with six-degree-of-freedom rigid body dynamics, integrated by a `MuJoCoModel`, so that the Airplane's position, orientation, and velocity at a given time step are driven by the previous step's aerodynamic loads, gravity, and any external loads. The wing geometry stays prescribed; it is the per-step `OperatingPoint` (body pose and rates) that the dynamics update. All `_CoupledUnsteadyProblem` and `CoreUnsteadyProblem` attributes (documented in the sections above) are inherited unchanged. The additions are the rigid body configuration (set once at construction) and the per-step aerodynamic load history (populated as the solve advances).

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

Each is stored in a `_`-prefixed backing slot and exposed through a read-only property of the unprefixed name.

| Attribute           | Type               | Notes                                                                                                        |
|---------------------|--------------------|--------------------------------------------------------------------------------------------------------------|
| `I_BP1_CgP1`        | `np.ndarray`       | Inertia matrix in the first Airplane's body axes about its CG (kg*m^2); the array itself is set read-only    |
| `mass`              | `float`            | Mass of the first Airplane (kg)                                                                              |
| `k_max`             | `int`              | Maximum strongly coupled sub-iterations per free-flight time step (a capped step is accepted with a warning) |
| `external_loads_fn` | `Callable \| None` | Optional callback returning additional (force, moment) loads to apply each step, or None                     |

The rigid body dynamics engine lives in the private `_mujoco_model` slot (a `MuJoCoModel`). Like the attributes above, it is set in `__init__` and never reassigned, and the reference stays fixed while the engine's own state advances during the solve. It has no public property: users never interact with the engine directly (see the MuJoCoModel section).

#### Derived from Immutable (read-only property, no backing slot)

| Property                | Depends On  | Notes                                                                                              |
|-------------------------|-------------|----------------------------------------------------------------------------------------------------|
| `_free_flight_movement` | `_movement` | Typed-narrow cast of the inherited `_movement` slot to `FreeFlightMovement` (recomputed, uncached) |

#### Mutable (populated by solver)

This is allocated in `__init__` (the validation guard as `False`) and updated by the problem's `initialize_next_problem` on the `external_loads_fn`'s first invocation. It is a plain slot rather than a read-only property because it must change after construction.

| Attribute                   | Type   | Notes                                                                                                      |
|-----------------------------|--------|------------------------------------------------------------------------------------------------------------|
| `_external_loads_validated` | `bool` | Once-only guard; flips to `True` after the `external_loads_fn` return is validated on its first invocation |

#### Construction-only parameters

`integrator`, `extra_xml`, and `mujoco_assets` are constructor parameters, not attributes: all three are validated here (`integrator` a str naming a supported MuJoCo integrator; `extra_xml` a dict or None with keys restricted to the permitted injection points and str values; `mujoco_assets` a dict or None mapping str filenames to bytes, where each filename must be a bare basename with a nonempty extension), then forwarded to the `MuJoCoModel` constructed in `__init__` and not stored on the problem, so none has a slot or an attribute-category entry above. They are the only raw user input reaching the `MuJoCoModel`, which performs no validation of its own; deeper XML and asset-reference correctness is left to MuJoCo, except that construction raises if the generated model XML references files that `mujoco_assets` does not cover, which keeps machine-specific paths out of saved files and every constructed problem saveable. See Construction-Only Parameters under Design Principles.

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

#### `FreeFlightMovement` and `AeroelasticMovement` variant additions

`FreeFlightMovement` and `AeroelasticMovement` extend `CoreMovement` directly, as feature-specific siblings of `Movement`. They inherit every attribute above unchanged and add the immutable state below, each stored in a `_`-prefixed backing slot and exposed through a read-only property of the unprefixed name. The pre-generation is asymmetric: `FreeFlightMovement` pre-generates only `Airplane`s (the solver produces `OperatingPoint`s from the rigid body dynamics), while `AeroelasticMovement` pre-generates only `OperatingPoint`s (deformed `Airplane` geometry depends on the solver's structural response).

| Attribute              | Type                               | Defined On            | Notes                                                                           |
|------------------------|------------------------------------|-----------------------|---------------------------------------------------------------------------------|
| `prescribed_num_steps` | `int`                              | `FreeFlightMovement`  | Prescribed-flight time steps before the free-flight phase                       |
| `free_num_steps`       | `int`                              | `FreeFlightMovement`  | Free-flight time steps after the prescribed phase                               |
| `airplanes`            | `tuple[tuple[Airplane, ...], ...]` | `FreeFlightMovement`  | Pre-generated prescribed Airplane geometry, indexed `[airplane_movement][step]` |
| `operating_points`     | `tuple[OperatingPoint, ...]`       | `AeroelasticMovement` | Pre-generated prescribed OperatingPoints                                        |

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

#### `AeroelasticWingMovement` variant additions

`AeroelasticWingMovement` extends `CoreWingMovement` directly, as a feature-specific sibling of `WingMovement`. It inherits every attribute above unchanged and adds the one immutable slot below, stored in a `_`-prefixed backing slot and exposed through a read-only property of the unprefixed name.

| Attribute                                     | Type                           | Notes                                                                                                                                                                                                           |
|-----------------------------------------------|--------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `spacingAnglesSecondDerivative_Gs_to_Wn_ixyz` | `tuple[Callable \| None, ...]` | Per-basis-direction (x, y, z) analytical second time derivative of the matching custom angular spacing; an entry is a callable when its `spacingAngles_Gs_to_Wn_ixyz` component is a custom callable, else None |

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

#### `FreeFlightOperatingPointMovement` variant additions

`FreeFlightOperatingPointMovement` extends `CoreOperatingPointMovement` directly, as a feature-specific sibling of `OperatingPointMovement`. It inherits every attribute above unchanged and adds the one mutable slot below, since its `OperatingPoint`s are produced by the solver's rigid body dynamics integration rather than prescribed.

| Attribute          | Type                   | Notes                                                                                                                                                                                                                                            |
|--------------------|------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `operating_points` | `list[OperatingPoint]` | Mutable OperatingPoint history; seeded with the base OperatingPoint at step 0, then the solver appends one per step. A plain mutable slot rather than a read-only property, mirroring the mutable-result-list treatment on `CoreUnsteadyProblem` |

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
| `g_E`                  | `np.ndarray`         | Gravity in Earth axes     |
| `omegas_BP1__E`        | `np.ndarray`         | Body angular velocity     |

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
| `T_pas_W_CgP1_to_E_CgP1`        | `T_pas_W_CgP1_to_BP1_CgP1`, `T_pas_BP1_CgP1_to_E_CgP1`       | Wind-to-Earth matrix (cached)     |
| `T_pas_E_CgP1_to_W_CgP1`        | `T_pas_W_CgP1_to_E_CgP1`                                     | Inverse of above (cached)         |
| `surfaceNormal_GP1`             | `surfaceNormal_E`, `T_pas_E_CgP1_to_GP1_CgP1`                | Surface normal in GP1 (cached)    |
| `surfacePoint_GP1_CgP1`         | `surfacePoint_E_Eo`, `CgP1_E_Eo`, `T_pas_E_CgP1_to_GP1_CgP1` | Surface point in GP1 (cached)     |
| `surfaceReflect_T_act_GP1_CgP1` | `surfacePoint_GP1_CgP1`, `surfaceNormal_GP1`                 | Active reflection matrix (cached) |
| `vInfHat_GP1__E`                | `T_pas_W_CgP1_to_GP1_CgP1`                                   | Freestream direction (cached)     |
| `vInf_GP1__E`                   | `vInfHat_GP1__E`, `vCg__E`                                   | Freestream velocity (cached)      |

**Note on caching**: While `vInfHat_GP1__E` and `vInf_GP1__E` are simple computations once the transformation matrix is cached, they are cached for consistency with the overall pattern and because they are called repeatedly during solver operations. The same is true for `qInf__E`.

**Note on `T_pas_GP1_CgP1_to_BP1_CgP1`**: This matrix is a constant (180-degree rotation about y) and does not depend on any `__init__` parameters. It is still lazily cached for consistency with the overall pattern and to avoid recomputing it on every access.

**Note on transformation decomposition**: The geometry-to-wind transformation (`T_pas_GP1_CgP1_to_W_CgP1`) is composed from `T_pas_GP1_CgP1_to_BP1_CgP1` and `T_pas_BP1_CgP1_to_W_CgP1` via `compose_T_pas`. Similarly, `T_pas_E_CgP1_to_GP1_CgP1` is composed from `T_pas_E_CgP1_to_BP1_CgP1` and `T_pas_BP1_CgP1_to_GP1_CgP1`, and `T_pas_W_CgP1_to_E_CgP1` is composed from `T_pas_W_CgP1_to_BP1_CgP1` and `T_pas_BP1_CgP1_to_E_CgP1`. Decomposing through body axes lets each non-body chain reuse the body-relative matrices rather than computing fresh rotations.

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

| Attribute                   | Type                 | Notes                                       |
|-----------------------------|----------------------|---------------------------------------------|
| `forces_W`                  | `np.ndarray \| None` | Forces in wind axes                         |
| `forceCoefficients_W`       | `np.ndarray \| None` | Force coefficients                          |
| `moments_W_CgP1`            | `np.ndarray \| None` | Moments relative to the first Airplane's CG |
| `momentCoefficients_W_CgP1` | `np.ndarray \| None` | Moment coefficients                         |
| `forces_G`                  | `np.ndarray \| None` | Forces in geometry axes                     |
| `forceCoefficients_G`       | `np.ndarray \| None` | Force coefficients                          |
| `moments_G_Cg`              | `np.ndarray \| None` | Moments relative to this Airplane's CG      |
| `momentCoefficients_G_Cg`   | `np.ndarray \| None` | Moment coefficients                         |
| `moments_W_Cg`              | `np.ndarray \| None` | Wind-axes rotation of `moments_G_Cg`        |
| `momentCoefficients_W_Cg`   | `np.ndarray \| None` | Moment coefficients                         |

#### Derived from Mutable (read-only property, no backing slot)

These read-only properties expose the named load components and coefficients defined in `AXES_POINTS_AND_FRAMES.md`. Each returns one signed component of its source array, or None while the source is None, and stores nothing of its own. They are recomputed on every access, since a cached value would go stale when the solver populates the source arrays.

| Property                         | Depends On                | Notes                       |
|----------------------------------|---------------------------|-----------------------------|
| `inducedDrag_W`                  | `forces_W`                | Negative of the x component |
| `sideForce_W`                    | `forces_W`                | The y component             |
| `lift_W`                         | `forces_W`                | Negative of the z component |
| `inducedDragCoefficient_W`       | `forceCoefficients_W`     | Negative of the x component |
| `sideForceCoefficient_W`         | `forceCoefficients_W`     | The y component             |
| `liftCoefficient_W`              | `forceCoefficients_W`     | Negative of the z component |
| `rollingMoment_W_Cg`             | `moments_W_Cg`            | The x component             |
| `pitchingMoment_W_Cg`            | `moments_W_Cg`            | The y component             |
| `yawingMoment_W_Cg`              | `moments_W_Cg`            | The z component             |
| `rollingMomentCoefficient_W_Cg`  | `momentCoefficients_W_Cg` | The x component             |
| `pitchingMomentCoefficient_W_Cg` | `momentCoefficients_W_Cg` | The y component             |
| `yawingMomentCoefficient_W_Cg`   | `momentCoefficients_W_Cg` | The z component             |

---

## Wing Class (`geometry/wing.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                   | Type                           | Notes                                                                                                                      |
|-----------------------------|--------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| `wing_cross_sections`       | `tuple[WingCrossSection, ...]` | Wing cross sections (tuple prevents mutation)                                                                              |
| `name`                      | `str`                          | Wing identifier                                                                                                            |
| `Ler_Gs_Cgs`                | `np.ndarray`                   | Leading edge root position                                                                                                 |
| `angles_Gs_to_Wn_ixyz`      | `np.ndarray`                   | Rotation angles                                                                                                            |
| `num_chordwise_panels`      | `int`                          | Chordwise panel count                                                                                                      |
| `chordwise_spacing`         | `str`                          | "cosine" or "uniform"                                                                                                      |
| `spanwise_mesh`             | `str`                          | "trapezoidal", "exploded", or "edge_defined", set by provenance during construction rather than as a constructor parameter |
| `leadingEdgePoints_Wn_Ler`  | `np.ndarray \| None`           | Original leading edge curve, a read-only array set only by `from_edge_points` (else None)                                  |
| `trailingEdgePoints_Wn_Ler` | `np.ndarray \| None`           | Original trailing edge curve, a read-only array set only by `from_edge_points` (else None)                                 |
| `tip_trim_fraction`         | `float \| None`                | Span fraction dropped off the tip while resampling, set only by `from_edge_points` (else None)                             |

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

#### Construction-only parameters

`explode_into_strips` is a constructor parameter, not an attribute: when True it triggers `_explode_wing` during initialization and the bool is then discarded, so it has no slot of its own. It does, however, set the derived immutable `spanwise_mesh` marker listed in the Immutable table above ("exploded" when True, "trapezoidal" otherwise). See Construction-Only Parameters under Design Principles.

The `from_edge_points` classmethod is the third source of the `spanwise_mesh` marker. It builds the single-strip cross sections from leading edge and trailing edge curves, constructs a normal Wing through `__init__` (which marks it "trapezoidal"), then sets the marker to "edge_defined" and fills the `leadingEdgePoints_Wn_Ler`, `trailingEdgePoints_Wn_Ler`, and `tip_trim_fraction` slots with the original curves and the tip trim fraction applied while resampling them before returning. These four immutable slots are assigned in-class outside `__init__`, the same way `__deepcopy__` and deserialization set slots. The Wing is never exposed in the intermediate "trapezoidal" state. Storing the original curves and the trim fraction (rather than the resampled cross sections) is what lets a future non-trapezoidal convergence tool resample them at a different number of WingCrossSections, which is why the marker is three-valued: only an "edge_defined" Wing carries curves, while an "exploded" Wing does not.

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
| `outline_A_Lp`      | `np.ndarray` | Outline coordinates            |
| `resample`          | `bool`       | Resampling flag                |
| `n_points_per_side` | `int`        | Points per side for resampling |
| `mcl_A_Lp`          | `np.ndarray` | Mean camber line coordinates   |

**Note**: The `add_control_surface` method creates and returns a new Airfoil instance rather than modifying the existing one. This is the correct immutable pattern.

**Note**: `outline_A_lp` and `mcl_A_lp` are deprecated read-only aliases for `outline_A_Lp` and `mcl_A_Lp`, and the constructor accepts a deprecated `outline_A_lp` alias for its `outline_A_Lp` parameter. Using any of the aliases emits a DeprecationWarning, and all three will be removed in v6.0.0.

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
| `forces_GP1`       | `np.ndarray \| None`      | Computed forces      |
| `moments_GP1_CgP1` | `np.ndarray \| None`      | Computed moments     |
| `forces_W`         | `np.ndarray \| None`      | Forces in wind axes  |
| `moments_W_CgP1`   | `np.ndarray \| None`      | Moments in wind axes |

#### Derived from Mutable (read-only property, no backing slot)

These read-only properties expose the named force components defined in `AXES_POINTS_AND_FRAMES.md`. Each returns one signed component of `forces_W`, or None while it is None, and stores nothing of its own. They are recomputed on every access, since a cached value would go stale when the solver populates `forces_W`.

| Property        | Depends On | Notes                       |
|-----------------|------------|-----------------------------|
| `inducedDrag_W` | `forces_W` | Negative of the x component |
| `sideForce_W`   | `forces_W` | The y component             |
| `lift_W`        | `forces_W` | Negative of the z component |

---

## MuJoCoModel Class (`_mujoco_model.py`)

`MuJoCoModel` is a private class that wraps MuJoCo's `MjModel` and `MjData` objects. It is constructed by `FreeFlightUnsteadyProblem` and provides methods for applying aerodynamic loads, stepping the rigid body dynamics, and extracting the updated state. Users pass raw scalars and arrays to `FreeFlightUnsteadyProblem` and never interact with `MuJoCoModel` directly.

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type                       | Notes                                                                                                                                        |
|------------------------|----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| `xml_str`              | `str`                      | Generated MuJoCo XML                                                                                                                         |
| `_model`               | `mujoco.MjModel`           | Compiled MuJoCo model; private slot with no property                                                                                         |
| `body_id`              | `int`                      | MuJoCo body ID for the Airplane                                                                                                              |
| `initial_key_frame_id` | `int`                      | MuJoCo key frame ID for initial conditions                                                                                                   |
| `initial_qpos`         | `np.ndarray`               | Initial generalized positions (computed during init)                                                                                         |
| `initial_qvel`         | `np.ndarray`               | Initial generalized velocities (computed during init)                                                                                        |
| `_mujoco_assets`       | `dict[str, bytes] \| None` | Retained assets dict (or None); private slot with no property; serialized into saved files so `_rebuild_engine` can resolve asset references |

#### Mutable

| Attribute | Type            | Notes                                                                                                   |
|-----------|-----------------|---------------------------------------------------------------------------------------------------------|
| `_data`   | `mujoco.MjData` | Mutated by `apply_loads`, `step`, `reset`, `restore_state`, `mj_forward`; private slot with no property |

#### Construction-only parameters

`integrator` and `extra_xml` are constructor parameters, not attributes: both shape the generated model during initialization, are folded into `xml_str` (so their content survives indirectly through the stored XML), and are then discarded, so neither has a slot or an attribute-category entry above. `mujoco_assets`, by contrast, is retained in the `_mujoco_assets` slot after being passed to MuJoCo's `from_xml_string`: an asset-based model cannot be rebuilt from `xml_str` alone, so the slot is serialized alongside the XML string (with the asset bytes encoded as base64) and `_rebuild_engine` passes the restored dict back to MuJoCo, which keeps saved files self-contained and load off the filesystem. `MuJoCoModel` does not validate any of the three: it is private and validates nothing, so they arrive already validated from `FreeFlightUnsteadyProblem` (the only constructor), with deeper XML and asset-reference correctness left to MuJoCo, except that `FreeFlightUnsteadyProblem` calls the model's `uncovered_file_references` method after construction and raises if the XML references files that the assets dict does not cover. See Construction-Only Parameters under Design Principles.

---


## Solver Classes (Not Covered Above)

The six solver classes (`SteadyHorseshoeVortexLatticeMethodSolver`, `SteadyRingVortexLatticeMethodSolver`, `UnsteadyRingVortexLatticeMethodSolver`, `CoupledUnsteadyRingVortexLatticeMethodSolver`, `AeroelasticUnsteadyRingVortexLatticeMethodSolver`, and `FreeFlightUnsteadyRingVortexLatticeMethodSolver`) are intentionally omitted from the immutability and lazy caching patterns described in this document. Unlike the data and geometry classes above, the solver classes are algorithmic classes whose attributes are internal mutable working state in a procedural computation pipeline. They are not shared data that external code accesses or modifies, so immutable properties, set once enforcement, and lazy caching would add significant boilerplate with no meaningful safety benefit. The solver classes do still use `__slots__`, like all other classes in the package, to protect against dynamic attribute assignment typos.

One narrow exception is `UnsteadyRingVortexLatticeMethodSolver.steady_problems`. It is not internal working state but results data that external code reads after a run, so rather than holding a separate mutable copy it is a read-only property that returns the underlying `UnsteadyProblem.steady_problems` directly. This keeps the problem as the single source of truth and removes any chance of the solver's view going stale, which matters for coupled problems whose `steady_problems` grows step by step during the solve.

The other exceptions are the solver attributes that describe a run rather than serve it, which are exposed through read-only properties over `_`-prefixed backing slots. Every solver has `ran`, a bool that starts False and is set True at the end of a completed `run()` call. The unsteady solver, and its coupled subclasses by inheritance, also has `prescribed_wake` and `force_method`, which report the configuration of the most recent `run()` call. Their backing slots are None until the first call to `run()` assigns them validated values, and reading either property while its backing slot is None raises a RuntimeError, so there is no pre-run default that could be mistaken for a real run's configuration. Unlike the immutable properties on the data classes, these values change over the solver's lifetime: `ran` flips once, and the run configuration is reassigned by each call to `run()`. The read-only property guarantees only that external code cannot write them.
