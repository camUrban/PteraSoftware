# OperatingPoint, Movement, and Problem Class Plan

This document describes the planned refactoring of class attributes to use a consistent
pattern of immutability and lazy caching across OperatingPoint, the movement classes,
SteadyProblem, and UnsteadyProblem. It follows the same taxonomic categories and design
principles established in PANEL_VORTEX_PLAN.md and GEOMETRY_PLAN.md.

## Design Principles (Recap)

### Attribute Categories

| Category                | Pattern                                                          | Example                            |
|-------------------------|------------------------------------------------------------------|------------------------------------|
| **Immutable**           | Read only property (no setter), set once in `__init__`           | `OperatingPoint.alpha`             |
| **Set Once**            | Property with setter that raises `AttributeError` if already set | (none in these classes)            |
| **Mutable**             | Property with setter, or plain attribute                         | `UnsteadyProblem.finalForces_W`    |
| **Derived (Immutable)** | Manual lazy caching, depends on immutable attributes             | `OperatingPoint.T_pas_GP1_to_W`    |

### Numpy Array Mutability

All numpy arrays that should be immutable are set to read only using
`arr.flags.writeable = False` in `__init__` immediately after assignment.

### List Collection Immutability

Store collections as tuples internally to prevent external mutation via `.append()`,
`.pop()`, etc.

### Stylistic Conventions

The following conventions from the previously refactored files (Panel, Wing,
WingCrossSection, Airfoil, etc.) should be carried forward. Read `_panel.py`, 

1. **Section header comments**: Use comments like
   `# --- Immutable: read only properties ---` to organize property groups within
   classes.

2. **Property organization order**:
   - Immutable read only properties (simple getters) (no docstrings as these are
     described in `__init__`'s docstring)
   - Immutable derived properties (lazy cached) (docstrings)
   - Set once properties (if any) (docstrings if they aren't parameters described in
     `__init__`)
   - Set once derived properties (if any) (docstrings)
   - Other methods (docstrings)
   - Class docstrings contains a list of all methods in the class. It only includes
     methods that are public, not parameter properties, and not setters. This list
     contains the name followed by the exact summary line from that method's docstring.
     The list must be in order of appearance in the class.

3. **Cache initialization**: All cache variables are initialized to `None` in `__init__`
   with a comment explaining they are caches.

4. **Read only enforcement**: For numpy arrays, set `flags.writeable = False`
   immediately after assignment in `__init__` or in the setter for set once properties.

5. **Tuple conversion**: When converting a list parameter to a tuple for immutability,
   validate the list first, then store as tuple with a comment explaining this prevents
   external mutation.

---

## OperatingPoint Class (`operating_point.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute      | Type    | Notes                    |
|----------------|---------|--------------------------|
| `rho`          | `float` | Fluid density            |
| `vCg__E`       | `float` | CG speed                 |
| `alpha`        | `float` | Angle of attack          |
| `beta`         | `float` | Sideslip angle           |
| `externalFX_W` | `float` | External force           |
| `nu`           | `float` | Kinematic viscosity      |

#### Derived from Immutable (use manual lazy caching)

| Property                   | Depends On                      | Notes                           |
|----------------------------|---------------------------------|---------------------------------|
| `qInf__E`                  | `rho`, `vCg__E`                 | Dynamic pressure (cached)       |
| `T_pas_GP1_CgP1_to_W_CgP1` | `alpha`, `beta`                 | Transformation matrix (cached)  |
| `T_pas_W_CgP1_to_GP1_CgP1` | `T_pas_GP1_CgP1_to_W_CgP1`      | Inverse transformation (cached) |
| `vInfHat_GP1__E`           | `T_pas_W_CgP1_to_GP1_CgP1`      | Freestream direction (cached)   |
| `vInf_GP1__E`              | `vInfHat_GP1__E`, `vCg__E`      | Freestream velocity (cached)    |

**Note on caching**: While `vInfHat_GP1__E` and `vInf_GP1__E` are simple computations
once the transformation matrix is cached, they are cached for consistency with the
overall pattern and because they are called repeatedly during solver operations. The
same is true for `qInf__E`.

---

## SteadyProblem Class (`problems.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute         | Type                     | Notes                                    |
|-------------------|--------------------------|------------------------------------------|
| `airplanes`       | `tuple[Airplane, ...]`   | Tuple prevents external mutation         |
| `operating_point` | `OperatingPoint`         | Operating conditions                     |

#### Derived from Immutable (use manual lazy caching)

| Property           | Depends On                     | Notes                                  |
|--------------------|--------------------------------|----------------------------------------|
| `reynolds_numbers` | `airplanes`, `operating_point` | Tuple of Re for each Airplane (cached) |

---

## OperatingPointMovement Class (`movements/operating_point_movement.py`)

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

| Property     | Depends On     | Notes    |
|--------------|----------------|----------|
| `max_period` | `periodVCg__E` | (cached) |

---

## WingCrossSectionMovement Class (`movements/wing_cross_section_movement.py`)

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

| Property      | Depends On    | Notes                                        |
|---------------|---------------|----------------------------------------------|
| `all_periods` | Period arrays | Tuple of unique non zero periods (cached)    |
| `max_period`  | Period arrays | Scalar float, longest period (cached)        |

---

## WingMovement Class (`movements/wing_movement.py`)

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

---

## AirplaneMovement Class (`movements/airplane_movement.py`)

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

---

## Movement Class (`movements/movement.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute                  | Type                               | Notes                     |
|----------------------------|------------------------------------|---------------------------|
| `airplane_movements`       | `tuple[AirplaneMovement, ...]`     | Tuple prevents mutation   |
| `operating_point_movement` | `OperatingPointMovement`           | Operating point changes   |
| `delta_time`               | `float`                            | Time step                 |
| `num_cycles`               | `int \| None`                      | Number of cycles          |
| `num_chords`               | `int \| None`                      | Number of chord lengths   |
| `num_steps`                | `int`                              | Total time steps          |
| `airplanes`                | `tuple[tuple[Airplane, ...], ...]` | Generated Airplanes       |
| `operating_points`         | `tuple[OperatingPoint, ...]`       | Generated OperatingPoints |

#### Derived from Immutable (use manual lazy caching)

| Property     | Depends On                                       | Notes    |
|--------------|--------------------------------------------------|----------|
| `lcm_period` | `airplane_movements`, `operating_point_movement` | (cached) |
| `max_period` | `airplane_movements`, `operating_point_movement` | (cached) |
| `static`     | `max_period`                                     | (cached) |

**Note on `airplanes` and `operating_points`**: These are generated during `__init__` by
calling the child movements' `generate_*` methods. They should be stored as nested
tuples to prevent modification after generation.

---

## UnsteadyProblem Class (`problems.py`)

### Attribute Classification

#### Immutable (set in `__init__`, never modified)

| Attribute              | Type                        | Notes                 |
|------------------------|-----------------------------|-----------------------|
| `movement`             | `Movement`                  | Movement definition   |
| `only_final_results`   | `bool`                      | Results flag          |
| `num_steps`            | `int`                       | Copied from Movement  |
| `delta_time`           | `float`                     | Copied from Movement  |
| `first_averaging_step` | `int`                       | Computed during init  |
| `first_results_step`   | `int`                       | Computed during init  |
| `steady_problems`      | `tuple[SteadyProblem, ...]` | Generated during init |

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

**Note**: The solver result lists must remain mutable as they are populated after
initialization by the solver. These are initialized as empty lists and appended to
during the solve.

---

## Summary of Changes

| Class                        | Changes                                                                                                                            |
|------------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| **OperatingPoint**           | Convert 6 scalars to read only properties; add lazy caching                                                                        |
| **SteadyProblem**            | Convert to read only properties; convert `airplanes` list to tuple; add lazy caching                                               |
| **OperatingPointMovement**   | Convert all attrs to read only properties; add lazy caching                                                                        |
| **WingCrossSectionMovement** | Convert to read only properties; set `flags.writeable = False` on numpy arrays; add lazy caching; add `__deepcopy__`               |
| **WingMovement**             | Convert to read only; tuple for `wing_cross_section_movements`; `flags.writeable = False`; add lazy caching; add `__deepcopy__`    |
| **AirplaneMovement**         | Convert to read only; tuple for `wing_movements`; `flags.writeable = False`; add lazy caching; add `__deepcopy__`                  |
| **Movement**                 | Convert to read only; tuples for all collections including generated `airplanes`/`operating_points`; add lazy caching              |
| **UnsteadyProblem**          | Convert config attrs to read only; convert `steady_problems` to tuple; keep result lists mutable                                   |

---

## Migration Steps

1. **OperatingPoint class**:
   - Convert all scalar attributes to read only properties with private storage
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments

2. **SteadyProblem class**:
   - Convert `airplanes` list to tuple after validation
   - Convert `operating_point` to read only property
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments

3. **OperatingPointMovement class**:
   - Convert all attributes to read only properties
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments

4. **WingCrossSectionMovement class**:
   - Convert all attributes to read only properties
   - Set `flags.writeable = False` on numpy arrays in `__init__`
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments
   - Add `__deepcopy__` method that copies numpy arrays and sets them read only

5. **WingMovement class**:
   - Convert `wing_cross_section_movements` list to tuple
   - Convert all other attributes to read only properties
   - Set `flags.writeable = False` on numpy arrays in `__init__`
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments
   - Add `__deepcopy__` method that copies numpy arrays and sets them read only

6. **AirplaneMovement class**:
   - Convert `wing_movements` list to tuple
   - Convert all other attributes to read only properties
   - Set `flags.writeable = False` on numpy arrays in `__init__`
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments
   - Add `__deepcopy__` method that copies numpy arrays and sets them read only

7. **Movement class**:
   - Convert `airplane_movements` list to tuple
   - Convert generated `airplanes` to tuple of tuples using nested comprehension:
     ```python
     # After generating the list of lists:
     self._airplanes = tuple(tuple(airplane_list) for airplane_list in airplanes_temp)
     ```
   - Convert generated `operating_points` to tuple:
     ```python
     self._operating_points = tuple(operating_points_temp)
     ```
   - Convert all other attributes to read only properties
   - Add lazy caching
   - Initialize cache variables to `None` in `__init__`
   - Add section header comments

8. **UnsteadyProblem class**:
   - Convert configuration attributes to read only properties
   - Convert `steady_problems` list to tuple after generation
   - Keep solver result lists as mutable plain attributes
   - Add section header comments

9. **Update tests**:
   - Add tests verifying immutability (setting read only property raises
     `AttributeError`)
   - Add tests verifying numpy array immutability (in place mutation raises
     `ValueError`)
   - Add tests verifying tuple collection immutability (calling `.append()` raises
     `AttributeError`)
   - Verify `generate_*` methods continue to work correctly
   - Verify solver can still populate UnsteadyProblem result lists
   - Verify `copy.deepcopy()` works correctly on Movement classes:
     - Deepcopied numpy arrays remain read only (`flags.writeable == False`)
     - Deepcopied objects are independent (modifying copy doesn't affect original)
     - Cached derived properties are reset to `None` in the copy
     - Used by `_compute_wake_area_mismatch` in `movement.py`

---

## Special Considerations

### Deepcopy Methods for Movement Classes with Numpy Arrays

The default `copy.deepcopy()` behavior does **not** preserve the `writeable = False` flag
on numpy arrays. When deepcopying an array, the resulting copy has `writeable = True` by
default. To maintain immutability invariants through deepcopy operations, classes with
immutable numpy array attributes require custom `__deepcopy__` methods.

**Classes that need custom `__deepcopy__` methods:**

- **WingCrossSectionMovement**: Has 6 numpy arrays (`ampLp_Wcsp_Lpp`, `periodLp_Wcsp_Lpp`,
  `phaseLp_Wcsp_Lpp`, `ampAngles_Wcsp_to_Wcs_ixyz`, `periodAngles_Wcsp_to_Wcs_ixyz`,
  `phaseAngles_Wcsp_to_Wcs_ixyz`)
- **WingMovement**: Has 7 numpy arrays (`ampLer_Gs_Cgs`, `periodLer_Gs_Cgs`,
  `phaseLer_Gs_Cgs`, `ampAngles_Gs_to_Wn_ixyz`, `periodAngles_Gs_to_Wn_ixyz`,
  `phaseAngles_Gs_to_Wn_ixyz`, `rotationPointOffset_Gs_Ler`)
- **AirplaneMovement**: Has 3 numpy arrays (`ampCg_GP1_CgP1`, `periodCg_GP1_CgP1`,
  `phaseCg_GP1_CgP1`)

**Classes that do NOT need custom `__deepcopy__` methods:**

- **OperatingPoint**: Only has `float` scalar immutable attributes. Derived cached
  properties (transformation matrices) are numpy arrays but can be recomputed.
- **OperatingPointMovement**: Only has `float`, `str | Callable`, and object reference
  attributes.
- **Movement**: Only has `tuple`, `float`, `int`, and object reference attributes.
- **SteadyProblem**: Only has `tuple` and object reference attributes.
- **UnsteadyProblem**: Only has scalars, `tuple`, and mutable `list` attributes.

The `__deepcopy__` implementation pattern follows the geometry classes (see
`WingCrossSection.__deepcopy__` for reference):

1. Create new instance via `object.__new__()` to skip `__init__` validation
2. Store in `memo` dict to handle circular references
3. Copy simple immutable attributes directly
4. Deepcopy object references (e.g., `base_wing_cross_section`) with memo
5. Copy tuples directly (they are immutable)
6. Copy numpy arrays via `.copy()` and set `flags.writeable = False`
7. Initialize cache variables to `None` (caches will be recomputed on access)

### Solver Result Lists

The `UnsteadyProblem` class has 12 result list attributes that are populated by the
solver after initialization. These must remain mutable (plain list attributes) since
the solver appends to them during the solve process. They are initialized as empty lists
in `__init__`.

### Generated Collections

The `Movement` class generates `airplanes` (a list of lists of Airplanes) and
`operating_points` (a list of OperatingPoints) during `__init__`. After generation,
these should be converted to tuples to prevent modification:

- `airplanes` becomes `tuple[tuple[Airplane, ...], ...]`
- `operating_points` becomes `tuple[OperatingPoint, ...]`

### Panel GP1_CgP1 Coordinates and SteadyProblem

During `SteadyProblem.__init__`, the Panel objects within the Airplanes have their
`*_GP1_CgP1` coordinates (e.g., `Frpp_GP1_CgP1`, `Flpp_GP1_CgP1`, etc.) populated via
coordinate transformations. These Panel attributes are **set once properties**
(implemented in PANEL_VORTEX_PLAN.md refactoring) that:

- Raise `AttributeError` if set more than once
- Set `flags.writeable = False` after assignment

This behavior is intentional and compatible with storing `airplanes` as a tuple. The
tuple immutability prevents replacing/adding/removing Airplane references, but does not
prevent the one time initialization of set once properties on objects referenced within
the tuple. This is the expected workflow: geometry is created, then global coordinates
are computed once when the problem is initialized.

### Deepcopy Compatibility

The `_compute_wake_area_mismatch` function in `movement.py` uses `copy.deepcopy()` on
`AirplaneMovement` and `OperatingPointMovement` objects. With the custom `__deepcopy__`
methods implemented for movement classes containing numpy arrays, deepcopy will correctly
preserve immutability:

- `AirplaneMovement.__deepcopy__` copies its numpy arrays and sets them read only
- `WingMovement.__deepcopy__` (called recursively) does the same for its arrays
- `WingCrossSectionMovement.__deepcopy__` (called recursively) does the same for its
  arrays
- `OperatingPointMovement` does not need a custom `__deepcopy__` as it has no numpy arrays
- The geometry objects (Airplane, Wing, etc.) already have proper `__deepcopy__`
  implementations from previous refactoring

Test cases should verify:
1. Deepcopied movement objects maintain read only numpy arrays
2. Deepcopied objects are independent (modifying one doesn't affect the other)
3. Cached derived properties are reset to `None` in the copy