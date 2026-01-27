# Geometry Class Plan

This document extends the panel and vortex class plan to cover the geometry classes:
Airplane, Wing, WingCrossSection, and Airfoil. It follows the same taxonomic categories
established in PANEL_VORTEX_PLAN.md.

## Design Principles (Recap)

### Attribute Categories

| Category      | Pattern                                                          | Example                        |
|---------------|------------------------------------------------------------------|--------------------------------|
| **Immutable** | Read-only property (no setter), set once in `__init__`           | `Airplane.name`                |
| **Set-Once**  | Property with setter that raises `AttributeError` if already set | `Wing.symmetry_type`           |
| **Mutable**   | Property with setter, or plain attribute                         | `Wing.wake_ring_vortices`      |
| **Derived**   | Manual lazy caching (check `None`, compute, cache, return)       | `Wing.projected_area`          |

---

## Airplane Class (`geometry/airplane.py`)

### Current State

The Airplane class processes Wings for symmetry during initialization, stores reference
dimensions, and holds solver outputs. The `__deepcopy__` method preserves Wings and
parameters while resetting solver state.

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

#### Derived (computed properties, no caching needed)

| Property                  | Depends On      | Notes                              |
|---------------------------|-----------------|------------------------------------|
| `num_panels`              | `wings`         | Sum of wing panel counts           |
| `T_pas_G_Cg_to_GP1_CgP1`  | `Cg_GP1_CgP1`   | Transformation matrix              |

**Note on derived properties**: These are simple computations that don't benefit from
caching. `num_panels` iterates over wings (small list), and `T_pas_G_Cg_to_GP1_CgP1`
calls a transformation function. Both are fast enough to compute on access.

### Implementation Changes

```python
class Airplane:
    def __init__(self, ...):
        # Immutable attributes (store directly, no setter needed)
        # Convert list to tuple to prevent external mutation
        self._wings: tuple[Wing, ...] = tuple(processed_wings)
        self._name = name
        self._weight = weight
        self._s_ref = s_ref
        self._c_ref = c_ref
        self._b_ref = b_ref

        # Immutable numpy array: make read-only to prevent in-place mutation
        self._Cg_GP1_CgP1 = Cg_GP1_CgP1
        self._Cg_GP1_CgP1.flags.writeable = False

        # Mutable solver state
        self.forces_W: np.ndarray | None = None
        self.forceCoefficients_W: np.ndarray | None = None
        self.moments_W_CgP1: np.ndarray | None = None
        self.momentCoefficients_W_CgP1: np.ndarray | None = None

    # --- Immutable: read-only properties ---
    @property
    def wings(self) -> tuple[Wing, ...]:
        return self._wings

    @property
    def name(self) -> str:
        return self._name

    @property
    def Cg_GP1_CgP1(self) -> np.ndarray:
        return self._Cg_GP1_CgP1

    @property
    def weight(self) -> float:
        return self._weight

    @property
    def s_ref(self) -> float:
        return self._s_ref

    @property
    def c_ref(self) -> float:
        return self._c_ref

    @property
    def b_ref(self) -> float:
        return self._b_ref

    # --- Derived: computed on access (no caching) ---
    @property
    def num_panels(self) -> int:
        return sum(wing.num_panels or 0 for wing in self._wings)

    def __deepcopy__(self, memo: dict) -> Airplane:
        ...
        # Deepcopy wings into a new tuple
        new_airplane._wings = tuple(
            copy.deepcopy(wing, memo) for wing in self._wings
        )
        # Copy numpy array and make it read-only
        new_airplane._Cg_GP1_CgP1 = np.copy(self._Cg_GP1_CgP1)
        new_airplane._Cg_GP1_CgP1.flags.writeable = False
        ...
```

---

## Wing Class (`geometry/wing.py`)

### Current State

The Wing class has complex lifecycle: attributes set in `__init__`, some modified by
`Airplane.process_wing_symmetry()` for type 5 symmetry, and others set by
`generate_mesh()`. The `__deepcopy__` method preserves mesh geometry while resetting
wake state.

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

**Note**: These are modified by `Airplane.process_wing_symmetry()` at lines 803-806
when type 5 symmetry is detected. The original Wing becomes a type 1 wing and a new
reflected Wing is created with type 3 symmetry.

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

#### Derived (computed properties)

| Property                     | Depends On            | Notes                  |
|------------------------------|-----------------------|------------------------|
| `T_pas_G_Cg_to_Wn_Ler`       | `symmetry_type`, etc. | Transformation matrix  |
| `T_pas_Wn_Ler_to_G_Cg`       | Above                 | Inverse transformation |
| `WnX_G`, `WnY_G`, `WnZ_G`    | Above                 | Basis vectors          |
| `children_T_pas_*`           | Cross sections        | Child transformations  |
| `projected_area`             | `panels`              | Projected area         |
| `wetted_area`                | `panels`              | Wetted area            |
| `average_panel_aspect_ratio` | `panels`              | Average aspect ratio   |
| `span`                       | Cross sections        | Wing span              |
| `standard_mean_chord`        | Above                 | Standard mean chord    |
| `mean_aerodynamic_chord`     | Above                 | MAC                    |

**Note on caching**: Most derived properties iterate over panels or cross sections.
For large meshes, caching `projected_area`, `wetted_area`, and `span` could provide
meaningful performance gains if they're accessed multiple times. However, since panels
are set-once, these could be cached after first computation without invalidation logic.

### Implementation Changes

```python
class Wing:
    def __init__(self, ...):
        # Immutable attributes
        # Convert list to tuple to prevent external mutation
        self._wing_cross_sections: tuple[WingCrossSection, ...] = tuple(
            wing_cross_sections
        )
        self._name = name
        self._num_chordwise_panels = num_chordwise_panels
        self._chordwise_spacing = chordwise_spacing

        # Immutable numpy arrays: make read-only to prevent in-place mutation
        self._Ler_Gs_Cgs = Ler_Gs_Cgs
        self._Ler_Gs_Cgs.flags.writeable = False
        self._angles_Gs_to_Wn_ixyz = angles_Gs_to_Wn_ixyz
        self._angles_Gs_to_Wn_ixyz.flags.writeable = False

        # Mutable (modified by process_wing_symmetry for type 5)
        self.symmetric = symmetric
        self.mirror_only = mirror_only
        self.symmetryNormal_G = symmetryNormal_G
        self.symmetryPoint_G_Cg = symmetryPoint_G_Cg

        # Set-once by generate_mesh
        self._symmetry_type: int | None = None
        self._num_spanwise_panels: int | None = None
        self._num_panels: int | None = None
        self._panels: np.ndarray | None = None

        # Mutable wake state
        self.wake_ring_vortices: np.ndarray | None = None
        self.gridWrvp_GP1_CgP1: np.ndarray | None = None

        # Optional: caches for derived properties
        self._projected_area: float | None = None
        self._wetted_area: float | None = None
        self._span: float | None = None

    # --- Immutable: read-only properties ---
    @property
    def wing_cross_sections(self) -> tuple[WingCrossSection, ...]:
        return self._wing_cross_sections

    @property
    def name(self) -> str:
        return self._name

    @property
    def Ler_Gs_Cgs(self) -> np.ndarray:
        return self._Ler_Gs_Cgs

    @property
    def angles_Gs_to_Wn_ixyz(self) -> np.ndarray:
        return self._angles_Gs_to_Wn_ixyz

    @property
    def num_chordwise_panels(self) -> int:
        return self._num_chordwise_panels

    @property
    def chordwise_spacing(self) -> str:
        return self._chordwise_spacing

    # --- Set-once: properties with single-assignment enforcement ---
    @property
    def symmetry_type(self) -> int | None:
        return self._symmetry_type

    @symmetry_type.setter
    def symmetry_type(self, value: int) -> None:
        if self._symmetry_type is not None:
            raise AttributeError("symmetry_type can only be set once")
        self._symmetry_type = value

    @property
    def num_spanwise_panels(self) -> int | None:
        return self._num_spanwise_panels

    @num_spanwise_panels.setter
    def num_spanwise_panels(self, value: int) -> None:
        if self._num_spanwise_panels is not None:
            raise AttributeError("num_spanwise_panels can only be set once")
        self._num_spanwise_panels = value

    @property
    def num_panels(self) -> int | None:
        return self._num_panels

    @num_panels.setter
    def num_panels(self, value: int) -> None:
        if self._num_panels is not None:
            raise AttributeError("num_panels can only be set once")
        self._num_panels = value

    @property
    def panels(self) -> np.ndarray | None:
        return self._panels

    @panels.setter
    def panels(self, value: np.ndarray) -> None:
        if self._panels is not None:
            raise AttributeError("panels can only be set once")
        self._panels = value

    # --- Derived with optional caching ---
    @property
    def projected_area(self) -> float | None:
        if self._panels is None:
            return None
        if self._projected_area is None:
            # ... computation ...
            self._projected_area = computed_area
        return self._projected_area

    def __deepcopy__(self, memo: dict) -> Wing:
        ...
        # Deepcopy wing cross sections into a new tuple
        new_wing._wing_cross_sections = tuple(
            copy.deepcopy(wcs, memo) for wcs in self._wing_cross_sections
        )
        # Copy numpy arrays and make them read-only
        new_wing._Ler_Gs_Cgs = np.copy(self._Ler_Gs_Cgs)
        new_wing._Ler_Gs_Cgs.flags.writeable = False
        new_wing._angles_Gs_to_Wn_ixyz = np.copy(self._angles_Gs_to_Wn_ixyz)
        new_wing._angles_Gs_to_Wn_ixyz.flags.writeable = False
        ...
```

---

## WingCrossSection Class (`geometry/wing_cross_section.py`)

### Current State

WingCrossSection holds cross section geometry and control surface parameters. Some
attributes are modified by the parent Wing/Airplane during symmetry processing.

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

#### Derived (computed properties)

| Property                   | Depends On                  | Notes                   |
|----------------------------|-----------------------------|-------------------------|
| `T_pas_Wcsp_Lpp_to_Wcs_Lp` | `Lp_Wcsp_Lpp`, `angles_...` | Transformation matrix   |
| `T_pas_Wcs_Lp_to_Wcsp_Lpp` | Above                       | Inverse transformation  |

### Implementation Changes

```python
class WingCrossSection:
    def __init__(self, ...):
        # Immutable attributes
        self._airfoil = airfoil
        self._num_spanwise_panels = num_spanwise_panels
        self._chord = chord
        self._control_surface_hinge_point = control_surface_hinge_point
        self._control_surface_deflection = control_surface_deflection
        self._spanwise_spacing = spanwise_spacing

        # Immutable numpy arrays: make read-only to prevent in-place mutation
        self._Lp_Wcsp_Lpp = Lp_Wcsp_Lpp
        self._Lp_Wcsp_Lpp.flags.writeable = False
        self._angles_Wcsp_to_Wcs_ixyz = angles_Wcsp_to_Wcs_ixyz
        self._angles_Wcsp_to_Wcs_ixyz.flags.writeable = False

        # Mutable (modified by process_wing_symmetry)
        self.control_surface_symmetry_type = control_surface_symmetry_type

        # Set-once by parent Wing
        self._validated: bool = False
        self._symmetry_type: int | None = None

    # --- Immutable: read-only properties ---
    @property
    def airfoil(self) -> Airfoil:
        return self._airfoil

    @property
    def num_spanwise_panels(self) -> int | None:
        return self._num_spanwise_panels

    @property
    def chord(self) -> float:
        return self._chord

    @property
    def Lp_Wcsp_Lpp(self) -> np.ndarray:
        return self._Lp_Wcsp_Lpp

    @property
    def angles_Wcsp_to_Wcs_ixyz(self) -> np.ndarray:
        return self._angles_Wcsp_to_Wcs_ixyz

    @property
    def control_surface_hinge_point(self) -> float:
        return self._control_surface_hinge_point

    @property
    def control_surface_deflection(self) -> float:
        return self._control_surface_deflection

    @property
    def spanwise_spacing(self) -> str | None:
        return self._spanwise_spacing

    # --- Set-once: properties with single-assignment enforcement ---
    @property
    def validated(self) -> bool:
        return self._validated

    @validated.setter
    def validated(self, value: bool) -> None:
        if self._validated:
            raise AttributeError("validated can only be set once")
        self._validated = value

    @property
    def symmetry_type(self) -> int | None:
        return self._symmetry_type

    @symmetry_type.setter
    def symmetry_type(self, value: int) -> None:
        if self._symmetry_type is not None:
            raise AttributeError("symmetry_type can only be set once")
        self._symmetry_type = value

    def __deepcopy__(self, memo: dict) -> WingCrossSection:
        ...
        # Copy numpy arrays and make them read-only
        new_wcs._Lp_Wcsp_Lpp = np.copy(self._Lp_Wcsp_Lpp)
        new_wcs._Lp_Wcsp_Lpp.flags.writeable = False
        new_wcs._angles_Wcsp_to_Wcs_ixyz = np.copy(self._angles_Wcsp_to_Wcs_ixyz)
        new_wcs._angles_Wcsp_to_Wcs_ixyz.flags.writeable = False
        ...
```

---

## Airfoil Class (`geometry/airfoil.py`)

### Current State

Airfoil is the simplest geometry class. All attributes are set during `__init__` and
never modified. The `add_control_surface` method returns a NEW Airfoil rather than
modifying the existing one.

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

### Implementation Changes

```python
class Airfoil:
    def __init__(self, ...):
        # All attributes are immutable
        self._name = name
        self._resample = resample
        self._n_points_per_side = n_points_per_side

        # These are computed during init but then never modified.
        # Use private attributes during init methods, then make read-only at the end.
        self._outline_A_lp: np.ndarray  # Set by _populate_outline or parameter
        self._mcl_A_lp: np.ndarray | None = None  # Set by _populate_mcl

        # ... validation, normalization, resampling using self._outline_A_lp ...
        # ... _populate_mcl sets self._mcl_A_lp ...

        # After all processing is complete, make arrays read-only
        self._outline_A_lp.flags.writeable = False
        if self._mcl_A_lp is not None:
            self._mcl_A_lp.flags.writeable = False

    # --- Immutable: read-only properties ---
    @property
    def name(self) -> str:
        return self._name

    @property
    def outline_A_lp(self) -> np.ndarray:
        return self._outline_A_lp

    @property
    def resample(self) -> bool:
        return self._resample

    @property
    def n_points_per_side(self) -> int:
        return self._n_points_per_side

    @property
    def mcl_A_lp(self) -> np.ndarray | None:
        return self._mcl_A_lp

    def __deepcopy__(self, memo: dict) -> Airfoil:
        ...
        # Copy numpy arrays and make them read-only
        new_airfoil._outline_A_lp = np.copy(self._outline_A_lp)
        new_airfoil._outline_A_lp.flags.writeable = False
        if self._mcl_A_lp is not None:
            new_airfoil._mcl_A_lp = np.copy(self._mcl_A_lp)
            new_airfoil._mcl_A_lp.flags.writeable = False
        else:
            new_airfoil._mcl_A_lp = None
        ...
```

**Note on internal outline processing**: During `__init__`, several methods modify
`_outline_A_lp` (normalization, resampling). These methods use the private attribute
directly. The `flags.writeable = False` is set only at the end of `__init__`, after
all processing is complete. This allows the init methods to work normally while still
protecting the final result from external mutation.

---

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
   - Update `__deepcopy__` to use private attributes and set `flags.writeable = False`

3. **Wing class**:
   - Convert immutable attributes to read-only properties
   - **Convert `wing_cross_sections` from list to tuple** to prevent external mutation
   - Set `flags.writeable = False` on `_Ler_Gs_Cgs` and `_angles_Gs_to_Wn_ixyz`
   - Add set-once enforcement for mesh-related attributes
   - Keep symmetry attributes and wake state mutable
   - Consider adding caches for `projected_area`, `wetted_area`, `span`
   - Update `generate_mesh` to use setters
   - Update `__deepcopy__` to deepcopy into a new tuple and set `flags.writeable = False`

4. **Airplane class**:
   - Convert all init-time attributes to read-only properties
   - **Convert `wings` from list to tuple** to prevent external mutation
   - Set `flags.writeable = False` on `_Cg_GP1_CgP1`
   - Keep solver outputs mutable
   - Update `process_wing_symmetry` to work with Wing's mutable symmetry attrs
   - Update `__deepcopy__` to deepcopy into a new tuple and set `flags.writeable = False`

5. **Update tests**:
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

**Affected attributes:**

| Class        | Attribute              | Store As                        | Return As                       |
|--------------|------------------------|---------------------------------|---------------------------------|
| `Airplane`   | `wings`                | `tuple[Wing, ...]`              | `tuple[Wing, ...]`              |
| `Wing`       | `wing_cross_sections`  | `tuple[WingCrossSection, ...]`  | `tuple[WingCrossSection, ...]`  |

**Implementation pattern:**

```python
# In __init__, convert the validated list to a tuple:
self._wings: tuple[Wing, ...] = tuple(processed_wings)

# The property returns the tuple directly (no copy needed):
@property
def wings(self) -> tuple[Wing, ...]:
    return self._wings
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
symmetry (lines 803-810):

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

**Immutable numpy arrays requiring this protection:**

| Class              | Attribute                  |
|--------------------|----------------------------|
| `Airplane`         | `Cg_GP1_CgP1`              |
| `Wing`             | `Ler_Gs_Cgs`               |
| `Wing`             | `angles_Gs_to_Wn_ixyz`     |
| `WingCrossSection` | `Lp_Wcsp_Lpp`              |
| `WingCrossSection` | `angles_Wcsp_to_Wcs_ixyz`  |
| `Airfoil`          | `outline_A_lp`             |
| `Airfoil`          | `mcl_A_lp`                 |

**Note on `__deepcopy__`:** When deepcopying, use `np.copy()` which returns a writable
array, then set `flags.writeable = False` on the copy:

```python
def __deepcopy__(self, memo: dict) -> Airplane:
    ...
    new_airplane._Cg_GP1_CgP1 = np.copy(self._Cg_GP1_CgP1)
    new_airplane._Cg_GP1_CgP1.flags.writeable = False
    ...
```
