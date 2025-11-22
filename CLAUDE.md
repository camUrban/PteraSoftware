# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Ptera Software Development Guidelines for Claude

## Project Overview

Ptera Software is a fast, easy-to-use, and open-source package for analyzing flapping-wing flight using unsteady and steady vortex lattice methods. It supports steady and unsteady aerodynamic simulations with formation flight capabilities.

## Example Commands

### Running Tests

```
Bash(.venv/Scripts/python.exe:-m:unittest:discover:-s:${WORKSPACE}/tests:*)
```

If that fails, try:

```
Bash(.venv/bin/python:-m:unittest:discover:-s:${WORKSPACE}/tests:*)
```

To run a specific test module in tests\unit\:

```
".venv/Scripts/python.exe" -m unittest tests.unit.test_module_name -v
```

For example:

```
".venv/Scripts/python.exe" -m unittest tests.unit.test_wing_cross_section_movement -v
```

### Running Scripts That Import Ptera Software

When running scripts outside the main pterasoftware directory that import the package (e.g., scripts in `experimental/`), you need to set `PYTHONPATH` to the project root:

```bash
cd ${WORKSPACE}/experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe script_name.py
```

On Unix-like systems:

```bash
cd ${WORKSPACE}/experimental && PYTHONPATH="$PWD/.." ../.venv/bin/python script_name.py
```

This pattern:
1. Changes to the script's directory
2. Sets `PYTHONPATH` to the parent directory (project root)
3. Runs the script using the virtual environment's Python interpreter

## Architecture Overview

### Core Package Structure
- **pterasoftware/**: Main package with modular solver architecture
  - **Geometry**: `geometry/` package
    - `geometry/_airfoils/` Directory containing data files with airfoil coordinates
    - `geometry/_meshing.py`: Wing mesh generation
    - `geometry/airfoil.py`: Airfoil class with coordinate generation
    - `geometry/airplane.py`: Airplane class with coordinate transformations
    - `geometry/wing.py`: Wing class with symmetry processing
    - `geometry/wing_cross_section.py`: WingCrossSection class with validation
  - **Movement System**: `movements/` package (flapping dynamics and motion)
    - `movements/_functions.py`: Movement utility functions
    - `movements/airplane_movement.py`: Airplane motion definitions
    - `movements/movement.py`: Core Movement class
    - `movements/operating_point_movement.py`: Operating condition changes
    - `movements/wing_cross_section_movement.py`: Wing cross section motion
    - `movements/wing_movement.py`: Wing flapping motion
  - `_aerodynamics.py`: Vortex elements and velocity calculations
  - `_functions.py`: Shared utility functions
  - `_panel.py`: Panel class for discretized mesh elements
  - `_parameter_validation.py`: Input validation functions
  - `_transformations.py` (coordinate transformations and rotations)
  - `convergence.py`: Convergence analysis tools
  - `operating_point.py` (OperatingPoint class)
  - `output.py`: Visualization and results processing
  - `problems.py` (SteadyProblem, UnsteadyProblem classes)
  - `steady_horseshoe_vortex_lattice_method.py`: Steady horseshoe VLM solver
  - `steady_ring_vortex_lattice_method.py`: Steady ring VLM solver
  - `trim.py`: Trim analysis functionality
  - `unsteady_ring_vortex_lattice_method.py`: Unsteady ring UVLM solver

### Key Features
- **Multiple Simulation Methods**: Steady horseshoe VLM, steady ring VLM, unsteady ring UVLM
- **Customizable Aircraft Geometry**: Multi-wing aircraft with arbitrary wing cross sections and airfoils
- **Time-Dependent Motion**: Flapping motion defined by sweep, pitch, and heave functions
- **Formation Flight**: Multi-airplane simulations supported since v2.0.0
- **High-Speed Computing**: JIT compilation via Numba for fast simulations
- **Trim Analysis**: Automatic search for trim operating points
- **GUI Interface**: Basic GUI available (beta stage)

### Solver Architecture Pattern
1. **Problem Definition** - Combine geometry with operating conditions
2. **Automatic Meshing** - Panel discretization from geometry
3. **VLM Computation** - Matrix-based influence coefficient methods
4. **Post-processing** - Force/moment calculation and visualization

**Formation Flight Coordinate Handling:**
- For multi-Airplane simulations, all solver computations use the first Airplane's geometry axes and CG (GP1_CgP1)
- Panels store both local (G_Cg) and formation (GP1_CgP1) coordinates; vortices only use GP1_CgP1
- Coordinate transformations from local to first-Airplane-axes occur in Problem.__init__()
- Forces and moments are calculated in GP1_CgP1, then transformed to wind axes for output

### Key Dependencies
- **NumPy/SciPy**: Core numerical computations
- **Numba**: JIT compilation for performance-critical loops
- **PyVista**: 3D mesh processing and visualization
- **PySide6**: GUI framework
- **Matplotlib**: 2D plotting and analysis output

### Python Version Constraint
Requires Python 3.11, but active development is done in 3.13  

### Package Dependencies
**Core Runtime Dependencies:**
- matplotlib >= 3.10.7, < 4.0.0
- numpy >= 2.3.5, < 3.0.0
- pyvista >= 0.46.4, < 1.0.0
- scipy >= 1.16.3, < 2.0.0
- numba >= 0.62.1, < 1.0.0
- cmocean >= 4.0.3, < 5.0.0
- tqdm >= 4.67.1, < 5.0.0
- webp >= 0.4.0, < 1.0.0
- PySide6 >= 6.10.1, < 7.0.0

**Development Dependencies:**
- black >= 25.11.0, < 26.0.0
- build >= 1.3.0, < 2.0.0
- codecov >= 2.1.13, < 3.0.0
- codespell >= 2.4.1, < 3.0.0
- docformatter >= 1.7.7, < 2.0.0
- mypy >= 1.18.2, < 2.0.0
- pre-commit >= 4.4.0, < 5.0.0
- PyInstaller >= 6.16.0, < 7.0.0
- setuptools >= 80.9.0, < 81.0.0
- twine >= 6.2.0, < 7.0.0
- wheel >= 0.45.1, < 0.46.0

## Writing Style Guidelines

### Terminology
- **"Ptera Software"**: When writing as text, always write as two words without a hyphen, each being capitalized (never "ptera", "ptera software", or "PteraSoftware"). When writing as a package to be imported, use ```import pterasoftware as ps```  
- **"cross section"**: Always write as two words, never hyphenated (not "cross-section")  
- **Object references**: When referring to code objects, use proper class naming convention. The capitalization indicates that we are talking about a code object, not an abstraction. You don't need to add "object" or "objects" after the class name since the capitalization already makes this clear (e.g. "update the Wings" instead of "update the Wing objects"). In summary, when talking about code objects:
  - ✅ "the previous WingCrossSection"
  - ❌ "the previous cross section"
  - ✅ "this Wing"
  - ❌ "this wing"
  - ✅ "update the Wings"
  - ❌ "update the Wing objects" (unnecessary)
- **Abstract references**: When referring to abstractions, use lowercase and separate individual words with a space (e.g. "an airplane's wings are used to generate lift" and "the cross section of a wing typically has a streamlined shape known as an airfoil"). This is to distinguish them from code objects.

### Docstring Style
- See pterasoftware/geometry/_meshing.py and pterasoftware/geometry/airfoil.py for gold standard docstring style. 
- Use reStructuredText (rST) formatting guidelines  
- Maintain consistent parameter descriptions  
- Preserve existing documentation structure and completeness unless we are explicitly updating or improving it  
- Always include units in parameter descriptions where applicable  
- Keep parameter descriptions complete.  
- If a parameter is a numpy array, specify the expected shape and data type. For example, say "(3,) ndarray of floats" for a 1D array of 3 floats.

## Code Style Guidelines

### Code Formatting
- Follow existing code style (black) and conventions
- Maintain consistent indentation and spacing
- Preserve existing comment structure and detail level
- Write comments as complete sentences ending with a period
- When initializing numpy arrays (e.g., using np.zeros, np.ones, np.empty, etc.), always specify the dtype parameter (e.g., dtype=float, dtype=int, dtype=object, etc.)

### Variable Naming
- Use descriptive variable names that clearly indicate their purpose  
- Use lowercase with underscores for variable names
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the documentation

#### Coordinate System Variable Naming Patterns

Variables must follow one of these four patterns based on their requirements:

1. **Axes without a point and without a frame**
   `[variable name]_[axes ID]`
   - Example: `force_W` force (in wind axes)

2. **Axes without a point and with a frame**
   `[variable name]_[axes ID]__[frame ID]` (note double underscore)
   - Example: `velocity_W__B` velocity (in wind axes, observed from body frame)

3. **Axes with a point and without a frame**
   `[variable name]_[axes ID]_[point ID]`
   - Example: `position_G_Cg` position (in geometry axes, relative to CG point)

4. **Only a frame (for scalar values like speed)**
   `[variable name]__[frame ID]`
   - Example: `speed__E` speed (observed from Earth frame)

#### Standardized IDs

**Axis Systems:**
- E: Earth, B: body, W: wind, G: geometry, Gs: geometry axes (after accounting for symmetry), Wn: wing, Wcs: wing cross section, Wcsp: wing cross section parent, A: airfoil

**Reference Points:**
- Cg: CG point, Cgs: CG (after accounting for symmetry), Ler: leading edge root point, Lp: leading point, Lpp: leading point parent

**Reference Frames:**
- E: Earth frame, B: body frame, Wn: wing frame, Wcs: wing cross section frame, Wcsp: wing cross section parent frame

#### Angle Vector Naming

**Passive angle vectors (relating two axis systems):**
- `[variable name]_[source axes ID]_to_[target axes ID]_[sequence ID]`
- Example: `angles_E_to_B_izyx` (angles from Earth to body axes using intrinsic zy'x" sequence)

**Active angle vectors (for rotation within current axis system):**
- `[variable name]_act_[sequence ID]`
- Example: `rotation_angles_act_izyx` (rotation angles using intrinsic zy'x" sequence)

#### Transformation Matrix Naming

**Passive transformation matrices:**
- Rotation: `R_pas_[source]_to_[target]`
- Full transformation: `T_pas_[source axes]_[source point]_to_[target axes]_[target point]`
- Examples: `R_pas_G_to_Wn`, `T_pas_G_Cg_to_Wn_Ler`

**Active transformation matrices:**
- Rotation: `[name]_R_act`
- Full transformation: `[name]_T_act`
- Examples: `mirror_R_act`, `translate_T_act`

#### Text Referencing Conventions

When referencing coordinate variables in comments and docstrings, use the following patterns:

1. **Axes without a point and without a frame**
   "[variable name] (in [axes name])"
   - Example: "force (in wind axes)"

2. **Axes without a point and with a frame**
   "[variable name] (in [axes name], observed from the [frame name])"
   - Example: "velocity (in wind axes, observed from the body frame)"

3. **Axes with a point and without a frame**
   "[variable name] (in [axes name], relative to the [point name])"
   - Example: "position (in geometry axes, relative to the CG point)"

4. **Only a frame (for scalar values like speed)**
   "[variable name] (observed from the [frame name])"
   - Example: "speed (observed from the Earth frame)"

**Context-Dependent References:**
- **Local context**: "...in wing axes..." (when already within Wing class)
- **Airplane-local**: "...in the first Wing's axes..." (when within Airplane class)
- **Non-local**: "...in the first Airplane's second Wing's axes..." (when outside Airplane)

#### Vector and Angle Vector Format Conventions

**Vector Format:**
- All (non-homogeneous) vectors must take the form `[x, y, z]`
- This applies to position vectors, force vectors, velocity vectors, etc.
- Example: `position_G_Cg = np.array([1.5, 2.0, -0.3])`
- Example: `force_W = np.array([100.2, -5.1, 85.7])`

**Angle Vector Format:**
- All angle vectors must take the form `[angleX, angleY, angleZ]`
- This applies to both passive and active angle vectors
- Angles are always in degrees unless explicitly noted otherwise
- Example: `angles_E_to_B_izyx = np.array([15.0, -5.2, 30.0])`
- Example: `rotationAngles_act_izyx = np.array([0.0, 45.0, 0.0])`

**Homogeneous Coordinates:**
- For 4x4 transformation matrices, vectors are converted to homogeneous form
- Position vectors: `[x, y, z, 1]` (has_point=True)
- Direction vectors: `[x, y, z, 0]` (has_point=False)

## Coordinate System Conventions

### Axes, Points, and Frames Overview

Ptera Software simulates flapping-wing dynamics and aerodynamics using several different axis systems, reference points, and reference frames.

**Key Concepts:**
- **Axis system** ("axes"): Contains information about three cartesian directions. They do not have an origin, therefore it's incorrect to write something like "the wing axes' origin."
- **Reference points** ("points"): Contains information about the location of a particular point in space.
- **Reference frame** ("frame"): Contains information about the location of an "observer" and their motion relative to what is observed.

**Vector Requirements:**
- Force vectors: require only an axis system
- Position vectors: require both axes and a reference point (origin)
- Moment vectors: require axes and reference point (point about which moment acts)
- Velocity/acceleration vectors: require both axis system and reference frame
- Some scalars (e.g., speed): require a reference frame

### Variable Naming Patterns

There are four useful combinations of axes, points, and frames. Variables are denoted by appending information using **IDs**, and described in text using **names**:

1. **Axes without a point and without a frame**
   `[variable name]_[axes ID]`
   "[variable name] (in [axes name])"

2. **Axes without a point and with a frame**
   `[variable name]_[axes ID]__[frame ID]`
   "[variable name] (in [axes name], observed from the [frame name])"

3. **Axes with a point and without a frame**
   `[variable name]_[axes ID]_[point ID]`
   "[variable name] (in [axes name], relative to the [point name])"

4. **Only a frame (for scalar values like speed)**
   `[variable name]__[frame ID]`
   "[variable name] (observed from the [frame name])"

**Naming Convention:** IDs move in scope from most to least specific. Names move from least to most specific.

### ID Abbreviations and Names

- E: Earth
- B: body
- P: Airplane
- P1: the first Airplane
- W: wind
- Pr: problem
- G: geometry axes
- GP1: the first Airplane's geometry axes
- Gs: geometry axes (after accounting for symmetry)
- Wn: Wing
- Wcs: WingCrossSection
- Wcsp: WingCrossSection parent
- A: Airfoil
- Cg: CG
- CgP1: the first Airplane's CG
- Cgs: CG (after accounting for symmetry)
- Ler: leading edge root point
- Lp: leading point
- Lpp: leading point parent
- ...pp...: Panel point (Fr=front right, Fl=front left, Bl=back left, Br=back right, C=collocation)
- ...bhvp...: bound HorseshoeVortex point
- ...brvp...: bound RingVortex point
- ...wrvp...: wake RingVortex point
- ...whvp...: wake HorseshoeVortex point
- ...lvp...: LineVortex point (S=start, E=end, C=center)

**Reference Frames:**
- E: Earth reference frame (inertial)
- B: Body reference frame (non-inertial, attached to Airplane's body)
- Wn: Wing reference frame (non-inertial, attached to Wing's leading edge root)
- Wcs: Wing cross section reference frame (non-inertial, attached to WingCrossSection's leading point)
- Wcsp: Wing cross section parent reference frame (non-inertial, attached to WingCrossSection's leading point parent)

### Axis System Definitions

**1. Earth Axes**
- Basis: North, East, Down (right-handed)
- Variables: `..._E...`
- Text: "...in Earth axes..."

**2. Body Axes**
- Basis: Front of Airplane, Right of Airplane, Bottom of Airplane (right-handed)
- Variables: `..._B...` (local) or `..._BP1...` (non-local)
- Text: "...in body axes..." or "...in the first Airplane's body axes..."

**3. Wind Axes**
- Basis: Parallel to freestream velocity, perpendicular directions defined via angle of attack/sideslip
- Variables: `..._W...` (local) or `..._WPr1...` (non-local)
- Text: "...in wind axes..." or "...in the first Problem's wind axes..."

**4. Geometry Axes**
- Basis: Back of Airplane, Right of Airplane, Top of Airplane (right-handed)
- Variables: `..._G...` (local) or `..._GP1...` (non-local)
- Text: "...in geometry axes..." or "...in the first Airplane's geometry axes..."

**5. Geometry Axes (After Accounting for Symmetry)**
- Basis: For a given Wing, the basis directions are identical to that Wing's Airplane's geometry axes if the Wing is non-symmetric or symmetric-continuous. For mirror-only Wings, the basis directions are that Wing's Airplane's geometry axes reflected about that Wing's symmetry plane.
- Right-handed for non-symmetric and symmetric-continuous Wings. Left-handed for mirror-only Wings.
- Variables: `..._Gs...` (local), `..._Gs1...` (Airplane-local), `..._Gs2P1...` (non-local)
- Text: "...in geometry axes (after accounting for symmetry)...", "...in geometry axes (after accounting for the first Wing's symmetry)...", "...in the first Airplane's geometry axes (after accounting for its second Wing's symmetry)..."

**6. Wing Axes**
- Basis: Back of Wing in first cross section plane, normal to plane, top surface
- Right-handed for non-symmetric/symmetric-continuous Wings, left-handed for mirror-only Wings
- Variables: `..._Wn...` (local), `..._Wn1...` (Airplane-local), `..._Wn2P1...` (non-local)
- Text: "...in wing axes...", "...in the first Wing's axes...", "...in the first Airplane's second Wing's axes..."

**7. Wing Cross Section Axes**
- Basis: Trailing edge in plane, normal to plane, top surface
- Handedness same as wing axes
- Variables: `..._Wcs...` (local), `..._Wcs1...` (Wing-local), `..._Wcs3Wn2...` (Airplane-local), `..._Wcs1Wn2P1...` (non-local)
- Text: "...in wing cross section axes...", "...in the first WingCrossSection's axes...", etc.

**8. Wing Cross Section Parent Axes**
- Basis: Identical to wing axes for first cross section, identical to previous cross section axes for subsequent ones
- Variables: `..._Wcsp...` with similar local/non-local patterns
- Text: "...in wing cross section parent axes..."

**9. Airfoil Axes**
- Basis: Chordwise to trailing point, Normal to chord toward upper line (2D)
- Variables: `..._A...` (local), `..._AWcs2...` (Wing-local), etc.
- Text: "...in airfoil axes...", "...in the second WingCrossSection's Airfoil's axes...", etc.

### Angle Vectors and Transformations

**Rotation Sequences (Tait-Bryan only):**
- ixyz, ixzy, iyxz, iyzx, izxy, izyx: intrinsic sequences (xy'z", xz'y", etc.)
- exyz, exzy, eyxz, eyzx, ezxy, ezyx: extrinsic sequences (xyz, xzy, etc.)

**Passive Angle Vectors (relating axis systems):**
- Variable pattern: `[name]_[source axes ID]_to_[target axes ID]_[sequence ID]`
- Text pattern: "[name] describing the orientation of [target axes name] relative to [source axes name] using an [sequence name] sequence"
- Example: `angles_E_to_B_izyx` = "angles describing the orientation of body axes relative to Earth axes using an intrinsic zy'x" sequence"

**Active Angle Vectors (rotating within axis system):**
- Variable pattern: `[name]_act_[sequence ID]`
- Text pattern: "[name] for rotation using an [sequence name] sequence"
- Example: `angles_act_izyx` = "angles for rotation using an intrinsic zy'x" sequence"

**Rotation and Transformation Matrices:**
- **Passive rotation matrices:** `R_pas_[source]_to_[target]` (3x3)
- **Passive transformation matrices:** `T_pas_[source axes]_[source point]_to_[target axes]_[target point]` (4x4 homogeneous)
  - It only changes how the same physical quantity is expressed (axes and/or point) and never introduces scaling or shear.
- **Active rotation matrices:** `[name]_R_act` (3x3)  
- **Active transformation matrices:** `[name]_T_act` (4x4 homogeneous)
  - Active transforms operate within a single axis system; they never change which axes a vector is expressed in and never introduce scaling or shear. For free vectors (e.g., forces), translation has no effect.

**Implementation Notes:**
- All angles in degrees unless noted otherwise
- Angle wrapping to (-180, 180] range
- Vectors treated as column vectors with left-multiplication by matrices
- Different sequences have different gimbal lock singularities

## Miscellaneous Guidelines
- Use clear, descriptive variable names
- Avoid abbreviations unless they are well-known in the context
- Never hyphenate words in docstrings or comments. This is because they are often incorrectly wrapped by docformatter and incorrectly rendered in PyCharm's quick documentation. For example, even though not standard grammar, it's okay to write "non symmetric" instead of "non-symmetric".
- In docstrings and comments, represent subtraction using hyphens surrounded by spaces ( - ); never use em-dashes (—) or en-dashes (–).
- In docstrings and comments, never use a multiplication sign (×); always use a lowercase x or an asterisk surrounded by spaces (" x " or " * ")
- In docstrings and comments, never use the pi symbol (π); always write "pi" instead (e.g., "2 * pi" not "2 π")
- In docstrings and comments, never use the approximately-equal sign (≈); always write " ~ " instead (e.g., "a ~ b" not "a≈b")
- **Coordinate and axis references**: When referring to axes, coordinates, or planes, use  lowercase letters without hyphens between coordinate letters and descriptors (e.g., "x axis", "y component", "xz plane", "z direction"). Never use uppercase letters for axis references in text.
- Never use emojis in code, comments, or docstrings
- Always use straight single and double quotes, not curly ones
- In .md files, use "…" instead of "...". In all other files, use "...".
- Always end *.py file with an empty line.