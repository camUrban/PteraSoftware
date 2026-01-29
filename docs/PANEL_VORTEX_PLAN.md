# Panel and Vortex Class Plan

This document describes the planned refactoring of class attributes to use a consistent
pattern of immutability and lazy caching across Panel and the vortex classes.

## Design Principles

### Attribute Categories

Each attribute falls into one of these categories:

| Category      | Pattern                                                          | Example               |
|---------------|------------------------------------------------------------------|-----------------------|
| **Immutable** | Read-only property (no setter), set once in `__init__`           | `Panel.Frpp_G_Cg`     |
| **Set-Once**  | Property with setter that raises `AttributeError` if already set | `Panel.Frpp_GP1_CgP1` |
| **Mutable**   | Property with setter, or plain attribute                         | `RingVortex.strength` |
| **Derived**   | Manual lazy caching (check `None`, compute, cache, return)       | `Panel.area`          |

### Key Decisions

1. **No cache invalidation for immutable/set-once attributes**: Since these are only
   set once, we don't need invalidation logic in setters.

2. **Enforce set-once semantics at runtime**: Set-once properties raise `AttributeError`
   if assigned a second time. This catches bugs early where code incorrectly attempts
   to modify values that should be immutable after initial assignment.

3. **Use manual lazy caching for all derived properties**: This approach:
   - Works consistently for properties derived from both immutable and set-once attributes
   - Enables future `__slots__` adoption (no dependency on `__dict__`)
   - Simplifies `__deepcopy__` (cache variables can be copied directly)
   - Maintains the existing pattern already used throughout the codebase

4. **Solver-mutable attributes remain mutable**: Properties like `strength` that the
   solver needs to set keep their setters.

### Numpy Array Mutability

Even with read-only properties, numpy arrays can still be mutated in place via the
getter (e.g., `panel.Frpp_G_Cg[0] = 999.0`). To prevent this, all numpy arrays that
should be immutable are set to read-only using `arr.flags.writeable = False`:

1. **Immutable arrays**: Set in `__init__` immediately after assignment
2. **Set-once arrays**: Set in the setter after assignment
3. **Derived cached arrays**: Set in the lazy property after computation (since numpy
   operations like subtraction create new writable arrays regardless of input writability)
4. **Deepcopy**: Use `.copy()` then set `flags.writeable = False` on the copy

**Status**: Implemented for `Panel` and `LineVortex`. Still needed for `RingVortex`
and `HorseshoeVortex`.

---

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

#### Set-Once (set after construction by meshing or Problem)

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

#### Mutable (set by solver)

| Attribute          | Type                      | Notes                |
|--------------------|---------------------------|----------------------|
| `ring_vortex`      | `RingVortex \| None`      | Attached vortex      |
| `horseshoe_vortex` | `HorseshoeVortex \| None` | Attached vortex      |
| `forces_GP1`       | `np.ndarray \| None`      | Computed forces      |
| `moments_GP1_CgP1` | `np.ndarray \| None`      | Computed moments     |
| `forces_W`         | `np.ndarray \| None`      | Forces in wind axes  |
| `moments_W_CgP1`   | `np.ndarray \| None`      | Moments in wind axes |

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

#### Derived from Set-Once (use manual lazy caching)

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

#### Mutable (solver sets these)

| Attribute  | Type    | Notes                                         |
|------------|---------|-----------------------------------------------|
| `strength` | `float` | Vortex strength (solver finds this)           |
| `age`      | `float` | Age in simulation time (incremented for wake) |

#### Derived from Immutable (use manual lazy caching)

| Property        | Depends On  | Notes             |
|-----------------|-------------|-------------------|
| `Crvp_GP1_CgP1` | All corners | Centroid position |
| `area`          | All corners | Vortex area       |

#### Derived (special: child objects)

| Property    | Depends On                                     | Notes              |
|-------------|------------------------------------------------|--------------------|
| `front_leg` | `Frrvp_GP1_CgP1`, `Flrvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `left_leg`  | `Flrvp_GP1_CgP1`, `Blrvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `back_leg`  | `Blrvp_GP1_CgP1`, `Brrvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `right_leg` | `Brrvp_GP1_CgP1`, `Frrvp_GP1_CgP1`, `strength` | Child `LineVortex` |

**Note on legs**: The leg `LineVortex` objects depend on both geometry (immutable)
and `strength` (mutable). Since `strength` is set by the solver AFTER the vortex is
created, we have two options:

1. **Recreate legs when strength is accessed** (simpler, slight overhead)
2. **Propagate strength updates to existing legs** (current approach)

Since legs are accessed repeatedly during induced velocity calculations, option 2
(keeping the current propagation) is more efficient.

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

#### Mutable (solver sets these)

| Attribute  | Type    | Notes           |
|------------|---------|-----------------|
| `strength` | `float` | Vortex strength |

#### Derived from Immutable (use manual lazy caching)

| Property         | Depends On                                                      | Notes            |
|------------------|-----------------------------------------------------------------|------------------|
| `Brhvp_GP1_CgP1` | `Frhvp_GP1_CgP1`, `leftLegVector_GP1`, `left_right_leg_lengths` | Back right point |
| `Blhvp_GP1_CgP1` | `Flhvp_GP1_CgP1`, `leftLegVector_GP1`, `left_right_leg_lengths` | Back left point  |

#### Derived (child objects)

| Property     | Depends On                                     | Notes              |
|--------------|------------------------------------------------|--------------------|
| `right_leg`  | `Brhvp_GP1_CgP1`, `Frhvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `finite_leg` | `Frhvp_GP1_CgP1`, `Flhvp_GP1_CgP1`, `strength` | Child `LineVortex` |
| `left_leg`   | `Flhvp_GP1_CgP1`, `Blhvp_GP1_CgP1`, `strength` | Child `LineVortex` |

**Note**: The `leftLegVector_GP1` setter currently normalizes the input vector. Since this
property becomes immutable, normalization must happen in `__init__`.

---

## LineVortex Class (`_vortices/_line_vortex.py`)

### Attribute Classification

Since `LineVortex` is an internal class whose endpoints ARE updated by parent vortex
classes when their corners change, we need to consider whether this mutation actually
happens.

**Analysis**: Looking at the parent classes, endpoint updates only happen in setters
that we're now removing (since corners are immutable). Therefore, `LineVortex`
endpoints are also effectively immutable after creation.

#### Immutable (set in `__init__`)

| Attribute       | Type         | Notes       |
|-----------------|--------------|-------------|
| `Slvp_GP1_CgP1` | `np.ndarray` | Start point |
| `Elvp_GP1_CgP1` | `np.ndarray` | End point   |

#### Mutable

| Attribute  | Type    | Notes                    |
|------------|---------|--------------------------|
| `strength` | `float` | Updated by parent vortex |

#### Derived from Immutable (use manual lazy caching)

| Property        | Depends On                       | Notes        |
|-----------------|----------------------------------|--------------|
| `vector_GP1`    | `Elvp_GP1_CgP1`, `Slvp_GP1_CgP1` | Line vector  |
| `Clvp_GP1_CgP1` | `Slvp_GP1_CgP1`, `vector_GP1`    | Center point |

---

## Summary of Changes

| Class               | Remove                                                          | Add                                                                       | Keep                                      |
|---------------------|-----------------------------------------------------------------|---------------------------------------------------------------------------|-------------------------------------------|
| **Panel**           | Setters on `_G_Cg` corners, cache invalidation in those setters | Set-once enforcement on meshing/`_GP1` attrs, `writeable=False` on caches | Manual lazy cache for all derived props   |
| **RingVortex**      | Setters on corners, child update logic in corner setters        | `writeable=False` on cached arrays                                        | Manual lazy cache, `strength` propagation |
| **HorseshoeVortex** | Setters on front points/vector/lengths, child update logic      | Normalization of `leftLegVector_GP1` in `__init__`, `writeable=False`     | Manual lazy cache, `strength` propagation |
| **LineVortex**      | Setters on endpoints                                            | `writeable=False` on cached arrays                                        | Manual lazy cache, plain `strength` attr  |

---

## Migration Steps

1. **Panel class**:
   - Remove setters from `_G_Cg` corner properties (make read-only)
   - Add `writeable = False` to all derived cached numpy arrays after computation
   - Convert meshing attributes (`is_right_edge`, `is_left_edge`, `local_chordwise_position`,
     `local_spanwise_position`) to properties with set-once enforcement
   - Convert `_GP1_CgP1` setters to enforce single assignment (raise `AttributeError` if already set)
   - Remove cache invalidation logic from the (now removed) `_G_Cg` setters

2. **RingVortex class**:
   - Remove setters from corner properties (make read-only)
   - Remove child `LineVortex` update logic from (now removed) corner setters
   - Add `writeable = False` to derived cached numpy arrays (`Crvp_GP1_CgP1`) after computation
   - Keep `strength` setter with child propagation

3. **HorseshoeVortex class**:
   - Remove setters from `Frhvp_GP1_CgP1`, `Flhvp_GP1_CgP1`, `leftLegVector_GP1`, `left_right_leg_lengths`
   - Move `leftLegVector_GP1` normalization from setter to `__init__`
   - Remove child update logic from (now removed) setters
   - Add `writeable = False` to derived cached numpy arrays (`Brhvp_GP1_CgP1`, `Blhvp_GP1_CgP1`)
     after computation
   - Keep `strength` setter with child propagation

4. **LineVortex class**:
   - Remove setters from `Slvp_GP1_CgP1` and `Elvp_GP1_CgP1` (make read-only)
   - Add `writeable = False` to derived cached numpy arrays (`vector_GP1`, `Clvp_GP1_CgP1`)
     after computation
   - Keep `strength` as plain attribute (no propagation needed)

5. **Update `__deepcopy__` (Panel only)**:
   - No structural changes needed since we're keeping manual lazy caching
   - Cache variables continue to be copied as before

6. **Update tests**:
   - Remove tests that verify cache invalidation on corner modification
   - Add tests verifying immutability (setting immutable property raises `AttributeError`)
   - Add tests verifying set-once enforcement (second assignment raises `AttributeError`)