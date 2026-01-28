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
getter (e.g., `panel.Frpp_G_Cg[0] = 999.0`). To prevent this, immutable numpy arrays
are set to read-only using `arr.flags.writeable = False` in `__init__` and after
set-once assignment. When deepcopying, use `.copy()` then set `flags.writeable = False`
on the copy.

---

## Panel Class (`_panel.py`)

### Current State

The Panel class currently has setters on all corner properties that invalidate caches.
This invalidation is unnecessary since corners are never modified after
construction/meshing.

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

### Implementation Changes

```python
class Panel:
    def __init__(self, ...):
        # Immutable geometry (store directly, no setter needed)
        self._Frpp_G_Cg = Frpp_G_Cg
        self._Flpp_G_Cg = Flpp_G_Cg
        self._Blpp_G_Cg = Blpp_G_Cg
        self._Brpp_G_Cg = Brpp_G_Cg
        self._is_leading_edge = is_leading_edge
        self._is_trailing_edge = is_trailing_edge

        # Caches for G-derived properties
        self._rightLeg_G: np.ndarray | None = None
        self._area: float | None = None
        # ... etc

        # Set-once by meshing (enforced via properties)
        self._is_right_edge: bool | None = None
        self._is_left_edge: bool | None = None
        self._local_chordwise_position: int | None = None
        self._local_spanwise_position: int | None = None

        # Set-once by Problem
        self._Frpp_GP1_CgP1: np.ndarray | None = None
        self._Flpp_GP1_CgP1: np.ndarray | None = None
        self._Blpp_GP1_CgP1: np.ndarray | None = None
        self._Brpp_GP1_CgP1: np.ndarray | None = None

        # Caches for GP1-derived properties
        self._rightLeg_GP1: np.ndarray | None = None
        self._Cpp_GP1_CgP1: np.ndarray | None = None
        # ... etc

        # Mutable solver state
        self.ring_vortex: RingVortex | None = None
        self.horseshoe_vortex: HorseshoeVortex | None = None
        self.forces_GP1: np.ndarray | None = None
        # ... etc

    # --- Immutable: read-only properties ---
    @property
    def Frpp_G_Cg(self) -> np.ndarray:
        return self._Frpp_G_Cg

    @property
    def is_leading_edge(self) -> bool:
        return self._is_leading_edge

    # --- Immutable derived: manual lazy caching ---
    @property
    def rightLeg_G(self) -> np.ndarray:
        if self._rightLeg_G is None:
            self._rightLeg_G = cast(np.ndarray, self._Frpp_G_Cg - self._Brpp_G_Cg)
        return self._rightLeg_G

    @property
    def area(self) -> float:
        if self._area is None:
            # ... computation ...
        return self._area

    # --- Set-once by meshing: properties with single-assignment enforcement ---
    @property
    def is_right_edge(self) -> bool | None:
        return self._is_right_edge

    @is_right_edge.setter
    def is_right_edge(self, value: bool) -> None:
        if self._is_right_edge is not None:
            raise AttributeError("is_right_edge can only be set once")
        self._is_right_edge = value

    @property
    def is_left_edge(self) -> bool | None:
        return self._is_left_edge

    @is_left_edge.setter
    def is_left_edge(self, value: bool) -> None:
        if self._is_left_edge is not None:
            raise AttributeError("is_left_edge can only be set once")
        self._is_left_edge = value

    @property
    def local_chordwise_position(self) -> int | None:
        return self._local_chordwise_position

    @local_chordwise_position.setter
    def local_chordwise_position(self, value: int) -> None:
        if self._local_chordwise_position is not None:
            raise AttributeError("local_chordwise_position can only be set once")
        self._local_chordwise_position = value

    @property
    def local_spanwise_position(self) -> int | None:
        return self._local_spanwise_position

    @local_spanwise_position.setter
    def local_spanwise_position(self, value: int) -> None:
        if self._local_spanwise_position is not None:
            raise AttributeError("local_spanwise_position can only be set once")
        self._local_spanwise_position = value

    # --- Set-once by Problem: properties with single-assignment enforcement ---
    @property
    def Frpp_GP1_CgP1(self) -> np.ndarray | None:
        return self._Frpp_GP1_CgP1

    @Frpp_GP1_CgP1.setter
    def Frpp_GP1_CgP1(self, value: np.ndarray) -> None:
        if self._Frpp_GP1_CgP1 is not None:
            raise AttributeError("Frpp_GP1_CgP1 can only be set once")
        self._Frpp_GP1_CgP1 = value

    # (similar pattern for Flpp_GP1_CgP1, Blpp_GP1_CgP1, Brpp_GP1_CgP1)

    # --- Set-once derived: manual lazy caching ---
    @property
    def rightLeg_GP1(self) -> np.ndarray | None:
        if self._Frpp_GP1_CgP1 is None or self._Brpp_GP1_CgP1 is None:
            return None
        if self._rightLeg_GP1 is None:
            self._rightLeg_GP1 = cast(np.ndarray, self._Frpp_GP1_CgP1 - self._Brpp_GP1_CgP1)
        return self._rightLeg_GP1
```

---

## RingVortex Class (`_aerodynamics.py`)

### Current State

Has setters on corner positions that update child `_LineVortex` objects and invalidate
caches. However, corners are never modified after construction in practice.

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

| Property    | Depends On                                     | Notes               |
|-------------|------------------------------------------------|---------------------|
| `front_leg` | `Frrvp_GP1_CgP1`, `Flrvp_GP1_CgP1`, `strength` | Child `_LineVortex` |
| `left_leg`  | `Flrvp_GP1_CgP1`, `Blrvp_GP1_CgP1`, `strength` | Child `_LineVortex` |
| `back_leg`  | `Blrvp_GP1_CgP1`, `Brrvp_GP1_CgP1`, `strength` | Child `_LineVortex` |
| `right_leg` | `Brrvp_GP1_CgP1`, `Frrvp_GP1_CgP1`, `strength` | Child `_LineVortex` |

**Note on legs**: The leg `_LineVortex` objects depend on both geometry (immutable)
and `strength` (mutable). Since `strength` is set by the solver AFTER the vortex is
created, we have two options:

1. **Recreate legs when strength is accessed** (simpler, slight overhead)
2. **Propagate strength updates to existing legs** (current approach)

Since legs are accessed repeatedly during induced velocity calculations, option 2
(keeping the current propagation) is more efficient.

### Implementation Changes

```python
class RingVortex:
    def __init__(self, ...):
        # Immutable geometry
        self._Frrvp_GP1_CgP1 = Frrvp_GP1_CgP1
        self._Flrvp_GP1_CgP1 = Flrvp_GP1_CgP1
        self._Blrvp_GP1_CgP1 = Blrvp_GP1_CgP1
        self._Brrvp_GP1_CgP1 = Brrvp_GP1_CgP1

        # Caches for derived properties
        self._Crvp_GP1_CgP1: np.ndarray | None = None
        self._area: float | None = None

        # Mutable
        self._strength = strength
        self.age: float = 0.0

        # Child objects (lazily created)
        self._front_leg: _LineVortex | None = None
        self._left_leg: _LineVortex | None = None
        self._back_leg: _LineVortex | None = None
        self._right_leg: _LineVortex | None = None

    # --- Immutable: read-only ---
    @property
    def Frrvp_GP1_CgP1(self) -> np.ndarray:
        return self._Frrvp_GP1_CgP1

    # (similar for other corners)

    # --- Mutable: strength with propagation to children ---
    @property
    def strength(self) -> float:
        return self._strength

    @strength.setter
    def strength(self, value: float) -> None:
        self._strength = value
        # Propagate to existing children
        if self._front_leg is not None:
            self._front_leg.strength = value
        if self._left_leg is not None:
            self._left_leg.strength = value
        if self._back_leg is not None:
            self._back_leg.strength = value
        if self._right_leg is not None:
            self._right_leg.strength = value

    # --- Derived from immutable: manual lazy caching ---
    @property
    def Crvp_GP1_CgP1(self) -> np.ndarray:
        if self._Crvp_GP1_CgP1 is None:
            self._Crvp_GP1_CgP1 = _functions.numba_centroid_of_quadrilateral(...)
        return self._Crvp_GP1_CgP1

    @property
    def area(self) -> float:
        if self._area is None:
            # ... computation ...
        return self._area

    # --- Child objects: lazy creation ---
    @property
    def front_leg(self) -> _LineVortex:
        if self._front_leg is None:
            self._front_leg = _LineVortex(
                Slvp_GP1_CgP1=self._Frrvp_GP1_CgP1,
                Elvp_GP1_CgP1=self._Flrvp_GP1_CgP1,
                strength=self._strength,
            )
        return self._front_leg
```

---

## HorseshoeVortex Class (`_aerodynamics.py`)

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

| Property     | Depends On                                     | Notes               |
|--------------|------------------------------------------------|---------------------|
| `right_leg`  | `Brhvp_GP1_CgP1`, `Frhvp_GP1_CgP1`, `strength` | Child `_LineVortex` |
| `finite_leg` | `Frhvp_GP1_CgP1`, `Flhvp_GP1_CgP1`, `strength` | Child `_LineVortex` |
| `left_leg`   | `Flhvp_GP1_CgP1`, `Blhvp_GP1_CgP1`, `strength` | Child `_LineVortex` |

### Implementation Changes

Same pattern as `RingVortex`:
- Immutable geometry: read-only properties
- Derived geometry (`Brhvp_GP1_CgP1`, `Blhvp_GP1_CgP1`): manual lazy caching
- Strength: mutable with propagation to children
- Child legs: lazy creation

**Note**: The `leftLegVector_GP1` setter currently normalizes the input vector. Since this
property becomes immutable, normalization must happen in `__init__`:

```python
def __init__(self, ..., leftLegVector_GP1: np.ndarray, ...):
    # Normalize during initialization
    self._leftLegVector_GP1 = leftLegVector_GP1 / np.linalg.norm(leftLegVector_GP1)
```

---

## _LineVortex Class (`_aerodynamics.py`)

### Current State

Has setters on endpoints that invalidate caches. Endpoints are set by parent
`RingVortex`/`HorseshoeVortex` during propagation.

### Attribute Classification

Since `_LineVortex` is an internal class whose endpoints ARE updated by parent vortex
classes when their corners change, we need to consider whether this mutation actually
happens.

**Analysis**: Looking at the parent classes, endpoint updates only happen in setters
that we're now removing (since corners are immutable). Therefore, `_LineVortex`
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

### Implementation Changes

```python
class _LineVortex:
    def __init__(self, ...):
        # Immutable geometry
        self._Slvp_GP1_CgP1 = Slvp_GP1_CgP1
        self._Elvp_GP1_CgP1 = Elvp_GP1_CgP1

        # Caches for derived properties
        self._vector_GP1: np.ndarray | None = None
        self._Clvp_GP1_CgP1: np.ndarray | None = None

        # Mutable
        self.strength = strength

    # --- Immutable: read-only ---
    @property
    def Slvp_GP1_CgP1(self) -> np.ndarray:
        return self._Slvp_GP1_CgP1

    @property
    def Elvp_GP1_CgP1(self) -> np.ndarray:
        return self._Elvp_GP1_CgP1

    # --- Derived: manual lazy caching ---
    @property
    def vector_GP1(self) -> np.ndarray:
        if self._vector_GP1 is None:
            self._vector_GP1 = self._Elvp_GP1_CgP1 - self._Slvp_GP1_CgP1
        return self._vector_GP1

    @property
    def Clvp_GP1_CgP1(self) -> np.ndarray:
        if self._Clvp_GP1_CgP1 is None:
            self._Clvp_GP1_CgP1 = self._Slvp_GP1_CgP1 + 0.5 * self.vector_GP1
        return self._Clvp_GP1_CgP1
```

---

## Summary of Changes

| Class               | Remove                                                          | Add                                                    | Keep                                       |
|---------------------|-----------------------------------------------------------------|--------------------------------------------------------|--------------------------------------------|
| **Panel**           | Setters on `_G_Cg` corners, cache invalidation in those setters | Set-once enforcement on meshing and `_GP1` attrs       | Manual lazy cache for all derived props    |
| **RingVortex**      | Setters on corners, child update logic in corner setters        | (none)                                                 | Manual lazy cache, `strength` propagation  |
| **HorseshoeVortex** | Setters on front points/vector/lengths, child update logic      | Normalization of `leftLegVector_GP1` in `__init__`     | Manual lazy cache, `strength` propagation  |
| **_LineVortex**     | Setters on endpoints                                            | (none)                                                 | Manual lazy cache, plain `strength` attr   |

---

## Migration Steps

1. **Panel class**:
   - Remove setters from `_G_Cg` corner properties (make read-only)
   - Keep existing manual lazy caching for `_G` derived properties (no change needed)
   - Convert meshing attributes (`is_right_edge`, `is_left_edge`, `local_chordwise_position`,
     `local_spanwise_position`) to properties with set-once enforcement
   - Convert `_GP1_CgP1` setters to enforce single assignment (raise `AttributeError` if already set)
   - Remove cache invalidation logic from the (now removed) `_G_Cg` setters

2. **RingVortex class**:
   - Remove setters from corner properties (make read-only)
   - Remove child `_LineVortex` update logic from (now removed) corner setters
   - Keep existing manual lazy caching for `Crvp_GP1_CgP1` and `area` (no change needed)
   - Keep `strength` setter with child propagation

3. **HorseshoeVortex class**:
   - Remove setters from `Frhvp_GP1_CgP1`, `Flhvp_GP1_CgP1`, `leftLegVector_GP1`, `left_right_leg_lengths`
   - Move `leftLegVector_GP1` normalization from setter to `__init__`
   - Remove child update logic from (now removed) setters
   - Keep existing manual lazy caching for `Brhvp_GP1_CgP1` and `Blhvp_GP1_CgP1` (no change needed)
   - Keep `strength` setter with child propagation

4. **_LineVortex class**:
   - Remove setters from `Slvp_GP1_CgP1` and `Elvp_GP1_CgP1` (make read-only)
   - Keep existing manual lazy caching for `vector_GP1` and `Clvp_GP1_CgP1` (no change needed)
   - Keep `strength` as plain attribute (no propagation needed)

5. **Update `__deepcopy__` (Panel only)**:
   - No structural changes needed since we're keeping manual lazy caching
   - Cache variables continue to be copied as before

6. **Update tests**:
   - Remove tests that verify cache invalidation on corner modification
   - Add tests verifying immutability (setting immutable property raises `AttributeError`):
     - `Panel.Frpp_G_Cg`, `Panel.Flpp_G_Cg`, `Panel.Blpp_G_Cg`, `Panel.Brpp_G_Cg`
     - `Panel.is_leading_edge`, `Panel.is_trailing_edge`
     - `RingVortex.Frrvp_GP1_CgP1`, etc.
     - `HorseshoeVortex.Frhvp_GP1_CgP1`, `leftLegVector_GP1`, etc.
     - `_LineVortex.Slvp_GP1_CgP1`, `_LineVortex.Elvp_GP1_CgP1`
   - Add tests verifying set-once enforcement (second assignment raises `AttributeError`):
     - `Panel.is_right_edge`, `Panel.is_left_edge`
     - `Panel.local_chordwise_position`, `Panel.local_spanwise_position`
     - `Panel.Frpp_GP1_CgP1`, `Panel.Flpp_GP1_CgP1`, `Panel.Blpp_GP1_CgP1`, `Panel.Brpp_GP1_CgP1`