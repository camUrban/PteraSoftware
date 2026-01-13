# bug/as6095_airfoil Branch Plan

## Goal

All NACA 4-series airfoils (within parameter constraints) and all database airfoils must pass validation checks. Use the validation scripts to verify:
```bash
cd experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe validate_airfoil_database.py
cd experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe validate_naca4_airfoils.py
```

### Current Validation Results

**Database Airfoils: 1256/1609 passed (353 failed)**
| Error                                 | Count |
|---------------------------------------|-------|
| Upper TE y not ≈ 0                    | 226   |
| Leading point y not ≈ 0               | 65    |
| Lower x not strictly increasing       | 29    |
| Lower TE y not ≈ 0                    | 19    |
| Upper TE y < Lower TE y               | 10    |
| Self-intersection (upper y < lower y) | 3     |
| Upper x not strictly decreasing       | 1     |

**NACA 4-Series: 1096/2260 passed (1164 failed)**
| Error                   | Count |
|-------------------------|-------|
| Leading point y not ≈ 0 | 1164  |

All NACA failures are the same error—cambered airfoils have their leading point on the camber line, not at y=0.

---

## Current Implementation

### Flow: Base Airfoil Creation
```
Airfoil(name, outline_A_lp, resample, n_points_per_side)
                    │
    ┌───────────────┴───────────────┐
    │                               │
outline_A_lp provided          outline_A_lp is None
    │                               │
    │                      _populate_outline()
    │                      (NACA gen or DB load)
    │                               │
    └───────────────┬───────────────┘
                    │
           _validate_outline()
                    │
            resample=True? ───yes───► _resample_outline()
                    │
             _populate_mcl()
```

### Flow: Control Surface Deflection
```
base_airfoil.add_control_surface(deflection, hinge_point)
                    │
                    ▼
         Split outline at hinge
                    │
                    ▼
      Rotate post-hinge portions
                    │
                    ▼
    Airfoil(..., _trust=_TRUST)  ◄── Bypasses _validate_outline()
                    │                (deflected TE would fail y≈0 check)
                    ▼
         _resample_outline() if needed
                    │
                    ▼
          _populate_mcl()
```

### Current `_validate_outline()` Checks
- At least 5 points
- Chord ≈ 1.0 (±2%)
- Thickness ≥ 0.1%
- Upper/lower trailing points at ≈ [1.0, 0.0]
- Leading point at ≈ [0.0, 0.0]
- Upper x strictly decreasing, lower x strictly increasing
- Upper y ≥ lower y at all x (no self-intersection)

### Current `_populate_mcl()` Normalization
- Computes MCL as average of upper/lower y values
- Resamples by arc-length
- Translates MCL LE to origin
- Scales MCL x to [0, 1] (y unchanged to preserve control surface effects)

---

## Planned Implementation: Three-Phase Validation

### Flow: Base Airfoil Creation (NEW)
```
Airfoil(name, outline_A_lp, resample, n_points_per_side)
                    │
    ┌───────────────┴───────────────┐
    │                               │
outline_A_lp provided          outline_A_lp is None
    │                               │
    │                      _populate_outline()
    │                      (NACA gen or DB load)
    │                               │
    └───────────────┬───────────────┘
                    │
      _validate_outline_preliminary()  ◄── Topology checks only
                    │
           _normalize_outline()        ◄── Translate, rotate, scale
                    │
        _validate_outline_final()      ◄── Position checks
                    │
            resample=True? ───yes───► _resample_outline()
                    │
             _populate_mcl()           ◄── MCL normalization still needed
```

### Flow: Control Surface Deflection (NEW)
```
base_airfoil.add_control_surface(deflection, hinge_point)
                    │
                    ▼
         Split outline at hinge
                    │
                    ▼
      Rotate post-hinge portions
                    │
                    ▼
    Airfoil(..., _trust=_TRUST)  ◄── Bypasses preliminary, normalize, final
                    │                (outline already normalized from base;
                    │                 deflected TE is intentional)
                    ▼
         _resample_outline() if needed
                    │
                    ▼
          _populate_mcl()        ◄── MCL normalization preserves deflection
                                     by not rotating (scales x only)
```

### Phase 1: `_validate_outline_preliminary()`
Topology checks for unfixable issues:
- At least 5 points
- Can identify a leading point (minimum x)
- Upper x strictly decreasing
- Lower x strictly increasing
- No self-intersection (upper y ≥ lower y)
- Both portions have ≥3 points

**Checks removed** (normalization handles these):
- ~~Chord ≈ 1.0~~ → normalization scales
- ~~Thickness ≥ 0.1%~~ → unnecessary for database/user airfoils
- ~~LE at [0,0]~~ → normalization translates
- ~~TE at [1,0]~~ → normalization rotates and scales

### Phase 2: `_normalize_outline()`
Transform to canonical form:
1. **Translate** leading point to origin
2. **Rotate** chord line onto x-axis
3. **Scale** to unit chord

### Phase 3: `_validate_outline_final()`
Verify normalization succeeded:
- LE at [0, 0] (within float tolerance)
- TE x ≈ 1.0 (both upper and lower)
- Upper TE y ≥ Lower TE y (catches inverted airfoils)
- No NaN values

**Note:** TE y ≈ 0 check relaxed to allow open/blunt trailing edges.

### Why MCL Normalization Is Still Needed
- Arc-length resampling shifts x values away from [0, 1]
- Must preserve y-values for control surface deflection effects
- Provides floating-point precision cleanup

### Expected Results After Implementation
| Source        | Before          | After             |
|---------------|-----------------|-------------------|
| NACA 4-series | 1096/2260 (48%) | ~2260/2260 (100%) |
| Database      | 1256/1609 (78%) | ~1566/1609 (97%)  |

Remaining ~43 database failures are truly invalid (non-monotonic, self-intersecting).

---

## Completed Changes (This Branch)

### Bug Fixes
1. **Reflected Airfoil resampling** - Added `__deepcopy__` to `Airfoil` class to prevent double-resampling when creating reflected wings
2. **Outline interpolation** - Fixed extrapolation errors by:
   - Interpolating only within overlapping x-range of upper/lower outlines
   - Disabling extrapolation in `PchipInterpolator`
   - Clamping interpolation inputs with `np.clip` to avoid NaN at boundaries
3. **MCL normalization** - Normalized mean camber line x-values to span [0.0, 1.0]

### Cleanup
- Removed unused `_get_mclY` method from `Airfoil`
- Increased padding in `Airfoil.draw()`
- Fixed linter complaints

### Database Changes
Removed 12 invalid airfoils with non-monotonic or discontinuous outlines:
- e664ex, n0009sm, n0012, n2414, n2415, n6409
- r1145msf, r1145msm, s1221-4deg-flap, ua79sff, ua79sfm

### Validation Tooling
Added experimental scripts:
- `validate_airfoil_database.py` - Validates all database airfoils
- `validate_naca4_airfoils.py` - Validates generated NACA 4-series airfoils