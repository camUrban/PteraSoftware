# Ptera Software Development Guidelines for Claude

## Writing Style Guidelines

### Terminology
- **"Ptera Software"**: When writing as text, always write as two words without a 
hyphen, each being capitalized (never "ptera", "ptera software", or "PteraSoftware"). 
When writing as a package to be imported, use ```import pterasoftware as ps```  
- **"cross section"**: Always write as two words, never hyphenated 
(not "cross-section")  
- **Object references**: When referring to code objects, use proper class naming 
convention. This will illustrate that we are talking about a code object, not an 
abstraction. Therefore, don't include the "object" suffix (e.g. "this WingCrossSection 
object") unless it is needed for clarity. In summary, when talking about code objects:
  - ✅ "the previous WingCrossSection"
  - ❌ "the previous cross section"
  - ✅ "this Wing"
  - ❌ "this wing"
- **Abstract references**: When referring to abstractions, use lowercase and separate 
individual words with a space (e.g. "an airplane's wings are used to generate lift" and 
"the cross section of a wing typically has a streamlined shape known as an airfoil"). 
This is to distinguish them from code objects.

### Docstring Style
- Follow existing PteraSoftware docstring conventions  
- Use reStructuredText (rST) formatting guidelines  
- Maintain consistent parameter descriptions  
- Preserve existing documentation structure and completeness unless we are explicitly 
updating or improving it  
- Always include units in parameter descriptions where applicable  
- Keep parameter descriptions complete.  
- If a parameter is a numpy array, specify the expected shape and data type. For  
example, say "(3,) ndarray of floats" for a 1D array of 3 floats.

## Code Style Guidelines

### Code Formatting
- Follow existing code style (black) and conventions  
- Maintain consistent indentation and spacing  
- Preserve existing comment structure and detail level  

### Variable Naming
- Use descriptive variables names that clearly indicate their purpose  
- Use lowercase with underscores for variable names
- Variables that are coordinates should be named 1D ndarrays, and have a suffix 
abbreviation indicating their reference frame. The reference frames used are:  
  - The wind frame (wind_frame)
  - The geometry frame (geometry_frame)
  - The wing frame (wing_frame)
  - The wing cross section frame (wing_cross_section_frame)
  - The airfoil frame (airfoil_frame)
  - The for a 1D ndarray of x-coordinates in the body frame, or "x_wing" for a 1D 
  ndarray of x-coordinates in the wing frame.

## Miscellaneous Guidelines
- Use clear, descriptive variable names  
- Avoid abbreviations unless they are well-known in the context  
- In docstrings and comments, never use em-dashes (—) or en-dashes (–); always use 
hyphens (-) for clarity  
- Never use emojis in code, comments, or docstrings  

---

## Appendix: Axes, Points, and Frames

This section standardizes how to name, reference, and reason about **axes**, **points**, and **frames** in Ptera Software. It follows (and slightly extends) Mark Drela’s "Flight Vehicle Aerodynamics" conventions while matching Ptera Software’s object hierarchy.

### 1) Quick definitions (use these precisely)

- **Axis system (“axes”)**: three basis directions (Cartesian, polar, or spherical) used to express vector components.
- **Reference point (“point”)**: a specific location used as an origin (for positions) or a moment center (for moments).
- **Reference frame (“frame”)**: the observer’s frame of motion (needed for time derivatives like velocity/acceleration).

**Which ingredients are required?**
- Generic vector (e.g., **force**): **axes** only.
- **Position** or **moment**: **axes + point**. (For moments, the point is the moment center, not “the origin.”)
- **Velocity/acceleration** (time derivatives of position): **axes + frame** (no point).

Ptera Software’s simulations can involve many axes/points/frames simultaneously (e.g. a simulation with multiple Airplanes, each with several Wings, with each of those having several WingCrossSections). Therefore, naming must encode *what* (axes/point/frame) and *whose* (local vs. non-local) they are.

### 2) Variable naming patterns

Append identifiers to base variable names using underscores as separators.

1. **Axes only**  
   `[name]_[axesID]`  
   Example: `f_W` → "f" (local wind axes).

2. **Axes + frame (no point)**  
   `[name]_[axesID]__[frameID]`  (double underscore between axes and frame)  
   Example: `v_B__E` → "v" (local body axes) observed from the Earth frame.

3. **Axes + point (no frame)**  
   `[name]_[axesID]_[pointID]`  
   Example: `m_B_Cg` → "m" (local body axes) relative to the local CG point.

> Do **not** combine point and frame on the same variable (choose the physically correct formulation).  

### 3) Text reference style (for docstrings/comments)

Use short, tags immediately after the quantity name. Keep the order **([axes]) relative to [point], observed from [frame]**.

- **Axes only:**  
  Pattern: `"... ([axes])"`  
  Example: `Force F (local wind axes).`

- **Axes + point (position/moment):**  
  Pattern: `"... ([axes]) relative to [point]"`  
  Examples: `Moment M (local body axes) relative to the second Airplane's CG point.`  
           `Position r (first Airplane's geometry axes) relative to relative to the first Airplane’s second Wing’s leading edge root point).`

- **Axes + frame (time-derivatives):**  
  Pattern: `"... ([axes]), observed from [frame]"`  
  Example: `Velocity v (second Airplane's body axes), observed from the Earth frame.`

- **Multiple ownership indices:**  
  Include explicit owners: `"... (local WingCrossSection axes), observed in the first Airplane's second Wing's frame"`

### 4) ID building blocks (most-specific → least-specific)

Use these abbreviations to build IDs; prepend local indices when needed.
- `E` : Earth
- `B` : body
- `P` : Airplane (ownership index follows, e.g., `P1`)
- `W` : wind
- `G` : geometry
- `Wn` : Wing (ownership index follows in non-local contexts, e.g., `Wn2`)
- `Wcs` : WingCrossSection (ownership index follows in non-local contexts, e.g., `Wcs3`)
- `A` : Airfoil (used with a Wcs context)

**Points**
- `I` : simulation starting point (first Airplane’s CG at t0)
- `Cgi` : Airplane CG at t0 (use ownership index if non-local, e.g., `CgiP1`)
- `Cg` : Airplane CG (use `CgP2`, etc., if non-local)
- `Ler` : Wing leading edge *root* point (e.g., `Ler2P1`)
- `Lp` : WingCrossSection *leading* point (e.g., `Lp1Wn2P1`)

**Panel / Vortex points (VLM)**
- Panel points: `...pp...` with suffixes for quadrant tag `Fr/Fl/Bl/Br`, collocation `C`, and indices `r[m]c[n]`.
  - Example: `Frppr3c2Wn2P1` = the first Airplane's second Wing's (3, 2) Panel's front right point
- Bound horseshoe: `...bhvp...`, bound ring: `...brvp...`, wake ring: `...wrvp...`, wake horseshoe: `...whvp...`
  (same quadrant/indexing patterns).
- Line vortex points: `lvp` with `S`(start), `E`(end), `C`(center), and leg tags `f/l/b/r`.

### 5) Axis systems (bases are right-handed unless noted)

**Earth axes (`E`)**  
- Basis: North, East, Down.  
- Text reference: “(in Earth axes)”.  
- Variable suffix: `_E`.

**Body axes (`B`)**  
- Basis: +x forward, +y right, +z down.  
- Local text: “(local body axes)”; non-local: “(the first Airplane’s body axes)”.  
- Variables: `_B` (local), `_BP1` (non-local for Airplane 1).

**Wind axes (`W`)**  
- Assumption: still airmass; freestream in the body frame equals airplane’s Earth-frame motion.  
- Basis:  
  1) +x parallel (not anti-parallel) to freestream as seen from the local body frame,  
  2) +y, +z chosen to be orthogonal; they are set via a thought experiment: find α, β such that a 2–3 **extrinsic** rotation aligns geometry-axes +x with wind-axes +x; the resulting +y and +z define wind +y and +z.  
  - This convention makes “lift” the component along wind +z and independent of sideslip by construction.  
- Text reference: “(local wind axes)”.  
- Variables: `_W` (local), `_WP1` (non-local for Airplane 1).

**Geometry axes (`G`)**  
- Basis: +x aft along fuselage, +y right, +z up.  
- Text reference: “(local geometry axes)”.  
- Variables: `_G` (local), `_GP1` (non-local for Airplane 1).

**Wing axes (`Wn`)**  
- Basis:  
  1) +x aft in the first wing-section’s plane,  
  2) +y normal to the first section’s plane toward the next section,  
  3) +z toward wing top surface.  
- Handedness: right-handed for non-symmetric and symmetric-continuous wings; **left-handed** for mirror-only wings.  
- Text reference: “(local Wing axes)”.  
- Variables: `_Wn` (local), `_Wn1` (Airplane-local first Wing), `_Wn2P1` (non-local: Wing 2 of Airplane 1).

**WingCrossSection axes (`Wcs`)**  
- Basis:  
  1) +x toward trailing edge in the section plane,  
  2) +y normal to the section plane toward the next section,  
  3) +z toward wing top surface.  
- Handedness mirrors the parent wing rule above.  
- Text reference: “(local WingCrossSection axes)”.  
- Variables: `_Wcs` (local), `_Wcs3Wn2` (Airplane-local), `_Wcs1Wn2P1` (non-local).

**Airfoil axes (`A`)** (2-D, tied to a Wcs)  
- Basis: +x chordwise toward trailing point; +y normal to chord toward upper line.  
- Text reference: “(local Airfoil axes)”.  
- Variables: `_A` (local), `_AWcs3Wn2`, `_AWcs1Wn2P1`.

### 6) Reference points (how to refer)

- **Simulation starting point (`I`)**: first Airplane’s CG at t0.  
  - Text: “relative to the simulation starting point”.  
  - Var: `_I`.

- **Starting points (`Cgi`)**: this Airplane’s CG at t0.  
  - Text: “relative to the local starting point”; non-local: “…the first Airplane’s starting point”.  
  - Var: `_Cgi`, `_CgiP1`.

- **CG points (`Cg`)**  
  - Text: “relative to the local CG point”; non-local: “…the second Airplane’s CG point”.  
  - Var: `_Cg`, `_CgP2`.

- **Leading edge root points (`Ler`)**: one per Wing  
  - Local/airplane-local/non-local as needed.  
  - Text example: “relative to first Wing's leading edge root point" (in Airplane-context)  
  - Var examples: `_Ler`, `_Ler1`, `_Ler2P1`.

- **Leading point (`Lp`)**: one per WingCrossSection  
  - Text example: “relative to the local leading point”.  
  - Var examples: `_Lp`, `_Lp1`, `_Lp1Wn2`, `_Lp1Wn2P1`.

- **Panel / vortex points**  
  - Panels (`…pp…`): e.g., `_Frppr3c2`, `_Frppr3c2Wn2`, `_Frppr3c2Wn2P1`.  
  - Bound HorseshoeVortex (`…bhvp…`), bound RingVortex (`…brvp…`), wake HorseshoeVortex (`…whvp…`), wake RingVortex (`…wrvp…`): same patterns.  
  - LineVortex (`lvp`): start/end/center prefix and front/right/back/left suffix (e.g., `_Slvp`, `_ClvpfBhvr3c2Wn2P2`).

### 7) Reference frames (who is observing)

- **Earth frame (`__E`)**: inertial, rigidly attached to the Earth.  
  - Text: “, observed from the Earth frame”.  
  - Var: `__E`.

- **Body frame (`__B`)**: non-inertial, rigidly attached to an Airplane's body.  
  - Text Example: “, observed from the local body frame”.  
  - Var Examples: `__B`, `__BP2`.

- **Wing frame (`__Wn`)**: non-inertial, attached to a Wing’s leading edge root point.  
  - Text Example: “, observed from the second Airplane's first Wing's frame”.  
  - Var Examples: `__Wn`, `__Wn2`, `__Wn2P4`.

- **WingCrossSection frame (`__Wcs`)**: non-inertial, attached to a WingCrossSection's leading point.  
  - Text Example: “, observed in the local WingCrossSection frame”.  
  - Var Examples: `__Wcs`, `__Wcs3Wn2`, `__Wcs3Wn2P4`.

### 8) Context rules and examples

**Local vs non-local**  
- If a variable is *inside* `Wing`, you may use `Wn` or `Wcs` without an Airplane index when the ownership is unambiguous at that scope.  
- If ambiguity exists (e.g., referencing another Airplane’s Wing from within `Wing`), include Airplane/Wing indices: `...Wn2P1`.

**More examples**
- `F_W` → Text: “F (wind axes)”.
- `V_B__E` → Text: “V (local body axes), observed from the Earth frame”.
- `a_G__B` → Text: “a (local geometry axes), observed from the local body frame”.
- `r_BP2_Cg` → Text: “r (the second airplane's body axes) relative to the local CG point”.
- `M_W_Ler2P1` → Text: “M (local wind axes) relative to the first Airplane's second Wing's leading edge root point”.

> Prefer brevity **only** when scope is unambiguous in code; otherwise include indices.

### 9) Common pitfalls (avoid these)

- **Anti-parallel wind +x**: by convention, wind +x is **parallel** to the freestream as seen from the body frame.  
- **Right-handed assumptions for mirror-only Wings**: their Wing axes and WingCrossSection axes are **left-handed**, so sign assumptions and cross products require extra care.  
- **Ambiguous ownership**: always include `P#`, `Wn#`, `Wcs#` when referencing outside the local object.

---
