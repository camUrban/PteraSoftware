# Geometry Class Plan

This document extends the panel and vortex class plan to cover the geometry classes:
Airplane, Wing, WingCrossSection, and Airfoil. It follows the same taxonomic categories
and logic established in PANEL_VORTEX_PLAN.md, covers the Panel, RingVortex,
HorseshoeVortex, and LineVortex classes, and has been fully implemented and tested.
Refer back to that plan and the corresponding code as the "gold standard" if any
questions arise here.

## Design Principles (Recap)

### Attribute Categories

| Category                | Pattern                                                          | Example                   |
|-------------------------|------------------------------------------------------------------|---------------------------|
| **Immutable**           | Read-only property (no setter), set once in `__init__`           | `Airplane.name`           |
| **Set-Once**            | Property with setter that raises `AttributeError` if already set | `Wing.symmetry_type`      |
| **Mutable**             | Property with setter, or plain attribute                         | `Wing.wake_ring_vortices` |
| **Derived (Immutable)** | Manual lazy caching, depends on immutable attributes             | `Wing.span`               |
| **Derived (Set-Once)**  | Manual lazy caching, depends on set-once attributes              | `Wing.projected_area`     |

### Numpy Array Mutability

Even with read-only properties, numpy arrays can still be mutated in place via the
getter (e.g., `airplane.Cg_GP1_CgP1[0] = 999.0`). To prevent this, all numpy arrays
that should be immutable are set to read-only using `arr.flags.writeable = False`:

1. **Immutable arrays**: Set in `__init__` immediately after assignment
2. **Set-once arrays**: Set in the setter after assignment
3. **Derived cached arrays**: Set in the lazy property after computation (since numpy
   operations like subtraction create new writable arrays regardless of input writability)
4. **Deepcopy**: Use `.copy()` then set `flags.writeable = False` on the copy

### Deepcopy Cache Handling

When implementing `__deepcopy__`, handle cached derived properties based on their source:

1. **Derived from Immutable → Preserve**: Copy the cached values (they remain valid since
   the source immutable attributes are also copied). For numpy arrays, use `.copy()` then
   set `flags.writeable = False`.

2. **Derived from Set-Once → Reset to None**: These depend on values that will be set
   fresh by the solver or meshing, so reset them.

| Class                | Preserve (Derived from Immutable)                        | Reset (Derived from Set-Once)                                     |
|----------------------|----------------------------------------------------------|-------------------------------------------------------------------|
| **Airplane**         | `_num_panels`, `_T_pas_G_Cg_to_GP1_CgP1`                 | (none)                                                            |
| **Wing**             | `_span`, `_T_pas_*` matrices, basis vectors              | `_projected_area`, `_wetted_area`, `_average_panel_aspect_ratio`  |
| **WingCrossSection** | `_T_pas_Wcsp_Lpp_to_Wcs_Lp`, `_T_pas_Wcs_Lp_to_Wcsp_Lpp` | (none)                                                            |
| **Airfoil**          | (no cached derived properties)                           | (none)                                                            |

---

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

#### Mutable (set by solver)

| Attribute                   | Type                 | Notes                |
|-----------------------------|----------------------|----------------------|
| `forces_W`                  | `np.ndarray \| None` | Forces in wind axes  |
| `forceCoefficients_W`       | `np.ndarray \| None` | Force coefficients   |
| `moments_W_CgP1`            | `np.ndarray \| None` | Moments in wind axes |
| `momentCoefficients_W_CgP1` | `np.ndarray \| None` | Moment coefficients  |

#### Derived from Immutable (use manual lazy caching)

| Property                  | Depends On      | Notes                    |
|---------------------------|-----------------|--------------------------|
| `num_panels`              | `wings`         | Sum of wing panel counts |
| `T_pas_G_Cg_to_GP1_CgP1`  | `Cg_GP1_CgP1`   | Transformation matrix    |

---

## Wing Class (`geometry/wing.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type                           | Notes                                    |
|------------------------|--------------------------------|------------------------------------------|
| `wing_cross_sections`  | `tuple[WingCrossSection, ...]` | Cross sections (tuple prevents mutation) |
| `name`                 | `str`                          | Wing identifier                          |
| `Ler_Gs_Cgs`           | `np.ndarray`                   | Leading edge root position               |
| `angles_Gs_to_Wn_ixyz` | `np.ndarray`                   | Rotation angles                          |
| `num_chordwise_panels` | `int`                          | Chordwise panel count                    |
| `chordwise_spacing`    | `str`                          | "cosine" or "uniform"                    |

#### Mutable (modified by `process_wing_symmetry` for type 5 symmetry)

| Attribute            | Type                 | Notes                           |
|----------------------|----------------------|---------------------------------|
| `symmetric`          | `bool`               | Modified to False for type 5    |
| `mirror_only`        | `bool`               | Modified to False for type 5    |
| `symmetryNormal_G`   | `np.ndarray \| None` | Modified to None for type 5     |
| `symmetryPoint_G_Cg` | `np.ndarray \| None` | Modified to None for type 5     |

**Note**: These are modified by `Airplane.process_wing_symmetry()` hen type 5 symmetry
is detected. The original Wing becomes a type 1 wing and a new reflected Wing is created
with type 3 symmetry.

#### Set-Once (set by `generate_mesh`, never modified after)

| Attribute             | Type                 | Set By          | Notes                |
|-----------------------|----------------------|-----------------|----------------------|
| `symmetry_type`       | `int \| None`        | `generate_mesh` | 1, 2, 3, or 4        |
| `num_spanwise_panels` | `int \| None`        | `generate_mesh` | Total spanwise count |
| `num_panels`          | `int \| None`        | `generate_mesh` | Total panel count    |
| `panels`              | `np.ndarray \| None` | `generate_mesh` | Panel matrix         |

#### Mutable (modified during simulation for wake)

| Attribute             | Type                 | Notes                           |
|-----------------------|----------------------|---------------------------------|
| `wake_ring_vortices`  | `np.ndarray \| None` | Wake vortex array, grows        |
| `gridWrvp_GP1_CgP1`   | `np.ndarray \| None` | Wake vortex positions, grows    |

#### Derived from Immutable (use manual lazy caching)

| Property                   | Depends On       | Notes                  |
|----------------------------|------------------|------------------------|
| `T_pas_G_Cg_to_Wn_Ler`     | Immutable attrs  | Transformation matrix  |
| `T_pas_Wn_Ler_to_G_Cg`     | Above            | Inverse transformation |
| `WnX_G`, `WnY_G`, `WnZ_G`  | Above            | Basis vectors          |
| `children_T_pas_*`         | Cross sections   | Child transformations  |
| `span`                     | Cross sections   | Wing span              |
| `standard_mean_chord`      | Cross sections   | Standard mean chord    |
| `mean_aerodynamic_chord`   | Cross sections   | MAC                    |

#### Derived from Set-Once (use manual lazy caching)

| Property                     | Depends On | Notes                |
|------------------------------|------------|----------------------|
| `projected_area`             | `panels`   | Projected area       |
| `wetted_area`                | `panels`   | Wetted area          |
| `average_panel_aspect_ratio` | `panels`   | Average aspect ratio |

**Note on caching**: Most derived properties iterate over panels or cross sections.
For large meshes, caching `projected_area`, `wetted_area`, and `span` could provide
meaningful performance gains if they're accessed multiple times. Since their source
attributes are immutable or set-once, these can be cached after first computation
without invalidation logic.

## WingCrossSection Class (`geometry/wing_cross_section.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                     | Type          | Notes                   |
|-------------------------------|---------------|-------------------------|
| `airfoil`                     | `Airfoil`     | Cross section airfoil   |
| `num_spanwise_panels`         | `int \| None` | Spanwise panel count    |
| `chord`                       | `float`       | Chord length            |
| `Lp_Wcsp_Lpp`                 | `np.ndarray`  | Position in parent axes |
| `angles_Wcsp_to_Wcs_ixyz`     | `np.ndarray`  | Rotation angles         |
| `control_surface_hinge_point` | `float`       | Hinge location (0-1)    |
| `control_surface_deflection`  | `float`       | Deflection in degrees   |
| `spanwise_spacing`            | `str \| None` | "cosine" or "uniform"   |

#### Mutable (modified by `process_wing_symmetry`)

| Attribute                       | Type          | Notes                           |
|---------------------------------|---------------|---------------------------------|
| `control_surface_symmetry_type` | `str \| None` | Set to None for type 5 symmetry |

**Note**: Modified at `airplane.py:810` when type 5 symmetry is split into two wings.

#### Set-Once (set by parent Wing)

| Attribute       | Type          | Set By               | Notes                   |
|-----------------|---------------|----------------------|-------------------------|
| `validated`     | `bool`        | `Wing.__init__`      | Validation flag         |
| `symmetry_type` | `int \| None` | `Wing.generate_mesh` | Inherited symmetry type |

#### Derived from Immutable (use manual lazy caching)

| Property                   | Depends On                  | Notes                   |
|----------------------------|-----------------------------|-------------------------|
| `T_pas_Wcsp_Lpp_to_Wcs_Lp` | `Lp_Wcsp_Lpp`, `angles_...` | Transformation matrix   |
| `T_pas_Wcs_Lp_to_Wcsp_Lpp` | Above                       | Inverse transformation  |

## Airfoil Class (`geometry/airfoil.py`)

### Attribute Classification

#### Immutable (all attributes)

| Attribute           | Type         | Notes                          |
|---------------------|--------------|--------------------------------|
| `name`              | `str`        | Airfoil identifier             |
| `outline_A_lp`      | `np.ndarray` | Outline coordinates            |
| `resample`          | `bool`       | Resampling flag                |
| `n_points_per_side` | `int`        | Points per side for resampling |
| `mcl_A_lp`          | `np.ndarray` | Mean camber line coordinates   |

**Note**: The `add_control_surface` method creates and returns a new Airfoil instance
rather than modifying the existing one. This is the correct immutable pattern.

## Summary of Changes

| Class                | Remove                      | Add                                         | Keep                      |
|----------------------|-----------------------------|---------------------------------------------|---------------------------|
| **Airplane**         | Direct attribute assignment | Read-only properties for immutable attrs    | Mutable solver state      |
| **Wing**             | Direct attribute assignment | Read-only for immutable, set-once for mesh  | Mutable symmetry & wake   |
| **WingCrossSection** | Direct attribute assignment | Read-only for immutable, set-once for flags | Mutable ctrl surface symm |
| **Airfoil**          | Direct attribute assignment | Read-only properties for all attrs          | (none mutable)            |

---

## Migration Steps

1. **Airfoil class** (simplest, no dependencies):
   - Convert all attributes to read-only properties
   - Internal methods use `_outline_A_lp` directly during init
   - Set `flags.writeable = False` on `_outline_A_lp` and `_mcl_A_lp` at end of init
   - Update `__deepcopy__` to use private attributes and set `flags.writeable = False`

2. **WingCrossSection class**:
   - Convert immutable attributes to read-only properties
   - Set `flags.writeable = False` on `_Lp_Wcsp_Lpp` and `_angles_Wcsp_to_Wcs_ixyz`
   - Add set-once enforcement for `validated` and `symmetry_type`
   - Keep `control_surface_symmetry_type` mutable (plain attribute)
   - Update `__deepcopy__` to use private attributes, set `flags.writeable = False`,
     and preserve derived-from-immutable caches

3. **Wing class**:
   - Convert immutable attributes to read-only properties
   - **Convert `wing_cross_sections` from list to tuple** to prevent external mutation
   - Set `flags.writeable = False` on `_Ler_Gs_Cgs` and `_angles_Gs_to_Wn_ixyz`
   - Add set-once enforcement for mesh-related attributes
   - Keep symmetry attributes and wake state mutable
   - Add caches for derived properties
   - Update `generate_mesh` to use setters
   - Update `__deepcopy__` to deepcopy into a new tuple, set `flags.writeable = False`,
     preserve derived-from-immutable caches, and reset derived-from-set-once caches

4. **Airplane class**:
   - Convert all init-time attributes to read-only properties
   - **Convert `wings` from list to tuple** to prevent external mutation
   - Set `flags.writeable = False` on `_Cg_GP1_CgP1`
   - Keep solver outputs mutable
   - Update `process_wing_symmetry` to work with Wing's mutable symmetry attrs
   - Update `__deepcopy__` to deepcopy into a new tuple, set `flags.writeable = False`,
     and preserve derived-from-immutable caches

5. **Update tests**:
   - Remove tests that verify cache invalidation on attribute modification (if any exist)
   - Remove tests that verify mutable behavior for now-immutable attributes
   - Add tests verifying immutability (setting read-only property raises `AttributeError`)
   - Add tests verifying set-once enforcement (second assignment raises `AttributeError`)
   - Add tests verifying numpy array immutability (in-place mutation raises `ValueError`)
   - **Add tests verifying tuple collection immutability** (calling `.append()` on `wings`
     or `wing_cross_sections` raises `AttributeError`)
   - Verify `process_wing_symmetry` still works with Wing's mutable symmetry attributes
   - Verify iteration over `wings` and `wing_cross_sections` still works correctly

---

## Special Considerations

### List Collection Immutability

Even with read-only properties, Python lists can still be mutated externally:

```python
# Without protection, this would modify the internal list:
airplane.wings.append(rogue_wing)  # Mutates the list!
airplane.wings.pop()               # Also mutates!
```

**Solution**: Store collections as tuples internally. Tuples are immutable, so any
attempt to modify them raises an `AttributeError`:

```python
# With tuple storage:
airplane.wings.append(rogue_wing)  # AttributeError: 'tuple' has no attribute 'append'
```

**Benefits:**
- No runtime overhead per access (unlike returning a copy)
- Clear error message when mutation is attempted
- Type checkers can catch incorrect usage (e.g., calling `.append()`)

**Note on `__deepcopy__`:** When deepcopying, create a new tuple from deepcopied
elements:

```python
def __deepcopy__(self, memo: dict) -> Airplane:
    ...
    new_airplane._wings = tuple(copy.deepcopy(wing, memo) for wing in self._wings)
    ...
```

**Note on iteration:** Tuples support iteration identically to lists, so existing
code that iterates over `wings` or `wing_cross_sections` will continue to work
without modification.

---

### Type 5 Symmetry Handling

The `Airplane.process_wing_symmetry` method modifies Wing attributes for type 5
symmetry:

```python
wing.symmetric = False
wing.mirror_only = False
wing.symmetryNormal_G = None
wing.symmetryPoint_G_Cg = None

for wing_cross_section in wing.wing_cross_sections:
    wing_cross_section.control_surface_symmetry_type = None
```

These attributes MUST remain mutable to support this use case. This is why they are
classified as **Mutable** rather than **Immutable** despite being set in `__init__`.

### Wing Caching Strategy

For derived properties that iterate over panels (`projected_area`, `wetted_area`,
`average_panel_aspect_ratio`), caching provides value since:
1. Panels are set-once, so no invalidation needed
2. These properties may be accessed multiple times during output/analysis
3. Iteration over panels has O(n) cost

The caches can be populated on first access and never invalidated.

### Numpy Array Mutability

Even with read-only properties, numpy arrays can still be mutated in place:
```python
airplane.Cg_GP1_CgP1[0] = 999.0  # Still works without protection!
```

To fully prevent mutation, we set immutable arrays as read-only using
`arr.flags.writeable = False`. This is done once in `__init__` after the array is
fully initialized. Any attempt to mutate the array will raise a `ValueError`.

```python
# In __init__, after validation and any processing:
self._Cg_GP1_CgP1 = validated_array
self._Cg_GP1_CgP1.flags.writeable = False
```

**Benefits:**
- No allocation overhead (unlike returning copies)
- Clear error message when mutation is attempted
- Code that copies first still works: `arr = airplane.Cg_GP1_CgP1.copy()` returns a
  writable copy

**Note on `__deepcopy__`:** When deepcopying, use `np.copy()` which returns a writable
array, then set `flags.writeable = False` on the copy:

```python
def __deepcopy__(self, memo: dict) -> Airplane:
    ...
    new_airplane._Cg_GP1_CgP1 = np.copy(self._Cg_GP1_CgP1)
    new_airplane._Cg_GP1_CgP1.flags.writeable = False
    ...
```