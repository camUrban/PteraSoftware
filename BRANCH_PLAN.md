# bug/as6095_airfoil Branch Plan

## Goal

All NACA 4-series airfoils (within parameter constraints) and all database airfoils must pass validation checks. Use the validation scripts to verify:
```bash
cd experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe validate_airfoil_database.py
cd experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe validate_naca4_airfoils.py
```

### Current Validation Results

**NACA 4-Series: 2260/2260 passed (100%)**

**Database Airfoils: 1570/1609 passed (39 failed)**
| Error                                      | Count |
|--------------------------------------------|-------|
| Upper TE y < Lower TE y (inverted TE)      | 31    |
| Self-intersection (upper y < lower y)      | 8     |

The 39 database failures are truly invalid airfoils:
- 31 have inverted trailing edges (as6093, as6096, as6098, dsma523a, etc.)
- 8 have self-intersecting outlines (e340, e378, fx38153, fx62k131, fx63147, etc.)

---

## Current Implementation: Three-Phase Validation

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
      _validate_outline_preliminary()  ◄── Topology checks only
                    │
           _normalize_outline()        ◄── Translate, rotate, scale
                    │
        _validate_outline_final()      ◄── Position checks + self-intersection
                    │
            resample=True? ───yes───► _resample_outline()
                    │
             _populate_mcl()           ◄── MCL normalization still needed
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

### Phase 1: `_validate_outline_preliminary()` (Implemented)
Topology checks for unfixable issues:
- At least 5 points
- Can identify a leading point (minimum x)
- Upper x non-increasing (allows adjacent points with same x)
- Lower x non-decreasing (allows adjacent points with same x)
- Both portions have at least 3 unique points

**Checks removed** (normalization handles these):
- ~~Chord approx 1.0~~ -> normalization scales
- ~~Thickness >= 0.1%~~ -> unnecessary for database/user airfoils
- ~~LE at [0,0]~~ -> normalization translates
- ~~TE at [1,0]~~ -> normalization rotates and scales

**Checks moved to final validation:**
- ~~No self-intersection~~ -> moved to after rotation for proper orientation

### Phase 2: `_normalize_outline()` (Implemented)
Transform to canonical form using iterative approach:
1. **Translate** leading point (minimum x) to origin
2. **Rotate** chord line (LE to average TE) onto x axis
3. **Check convergence**: If minimum x point changed, repeat from step 1
4. **Scale** to unit chord

The iterative approach handles airfoils with implicit angle of attack, where the initial
minimum x point is not the true aerodynamic leading edge. After rotation, a different
point may become minimum x; iteration continues until the leading point is stable.

### Phase 3: `_validate_outline_final()` (Implemented)
Verify normalization succeeded:
- No NaN values
- LE at [0, 0] (within 1e-6 tolerance)
- TE x = 1.0 (both upper and lower, within 1e-6 tolerance after extrapolation)
- Upper TE y >= Lower TE y (catches inverted airfoils)
- No self-intersection (upper y **strictly greater than** lower y using linear interpolation)

**Note:** TE y approx 0 check relaxed to allow open/blunt trailing edges.
**Note:** Self-intersection check uses strict inequality to reject zero-thickness airfoils.

---

## `_populate_mcl()` Normalization (Unchanged)
- Computes MCL as average of upper/lower y values
- Resamples by arc-length
- Translates MCL LE to origin
- Scales MCL x to [0, 1] (y unchanged to preserve control surface effects)

### Why MCL Normalization Is Still Needed
- Arc-length resampling shifts x values away from [0, 1]
- Must preserve y-values for control surface deflection effects
- Provides floating-point precision cleanup

---

## Completed Changes (This Branch)

### Three-Phase Validation (Complete)
1. **`_validate_outline_preliminary()`** - Topology checks only
2. **`_normalize_outline()`** - Iterative translate, rotate, scale to canonical form
3. **`_validate_outline_final()`** - Position checks and self-intersection
4. **Updated `__init__`** - Now calls three-phase validation instead of `_validate_outline()`

### Bug Fixes (Previous Work)
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
- Simplified self-intersection check to sample at union of all unique x values from both surfaces (`np.union1d`) instead of linearly-spaced points
- Fixed duplicate removal to only remove interior duplicates, preserving first/last points for closed trailing edges

### Database Changes
Removed 12 invalid airfoils with non-monotonic or discontinuous outlines:
- e664ex, n0009sm, n0012, n2414, n2415, n6409
- r1145msf, r1145msm, s1221-4deg-flap, ua79sff, ua79sfm

### Validation Tooling
Added experimental scripts:
- `validate_airfoil_database.py` - Validates all database airfoils
- `validate_naca4_airfoils.py` - Validates generated NACA 4-series airfoils

---

## Next Steps

1. Get rid of all PChip interpolation in airfoil code and replace with piecewise linear interpolation.
2. Consider if it would be useful to add debug-level logging messages. Note: It may be that these would clog up output, and aren't a good idea. However, if you determine that they would be useful, add them using the patterns established in other parts of the codebase.
3. Document the any changes to acceptable airfoil outlines in the docstrings within `airfoil.py`.
4. Create a script in experimental to more closely analyze the failing database airfoils to manually check if they are malformed.
5. Update unit tests to match new validation behavior
