# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Ptera Software Development Guidelines for Claude

## Project Overview

Ptera Software is a fast, easy-to-use, and open-source package for analyzing flapping-wing flight using unsteady and steady vortex lattice methods. It supports steady and unsteady aerodynamic simulations with formation flight capabilities.

## Development Commands

### Testing
```bash
# Run all tests with coverage
coverage run --source=pterasoftware -m unittest discover -s tests

# Run unit tests only
python -m unittest discover -s tests/unit

# Run integration tests only  
python -m unittest discover -s tests/integration

# Run a specific test file
python -m unittest tests.unit.test_vortex
```

### Code Quality
```bash
# Format code (handled by pre-commit hooks)
black .

# Build package for distribution
python -m build

# Install development dependencies
pip install -r requirements_dev.txt

# Install package in development mode
pip install -e .
```

### GUI Development
```bash
# Run GUI application
python main.py

# Build GUI executable
python -O -m PyInstaller --noconfirm "pterasoftware.spec"
```

## Architecture Overview

### Core Package Structure
- **pterasoftware/**: Main package with modular solver architecture
  - **Geometry**: `geometry/` package (Airplane → Wing → WingCrossSection → Airfoil hierarchy)
    - `geometry/airplane.py`: Airplane class with coordinate transformations
    - `geometry/wing.py`: Wing class with symmetry processing
    - `geometry/wing_cross_section.py`: WingCrossSection class with validation
    - `geometry/airfoil.py`: Airfoil class with coordinate generation
  - **Transformations**: `transformations.py` (coordinate system transformations and rotations)
  - **Solvers**: Three VLM implementations (steady horseshoe, steady ring, unsteady ring)
  - **Problems**: `problems.py` (SteadyProblem, UnsteadyProblem classes)
  - **Support**: meshing, movement, output, trim, convergence modules

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

### Key Dependencies
- **NumPy/SciPy**: Core numerical computations
- **Numba**: JIT compilation for performance-critical loops
- **PyVista**: 3D mesh processing and visualization
- **PySide6**: GUI framework
- **Matplotlib**: 2D plotting and analysis output

### Python Version Constraint
Requires Python 3.10.0 to < 3.11.0 (strict constraint for dependency compatibility)  

### Package Dependencies
**Core Runtime Dependencies:**
- matplotlib >= 3.10.3, < 4.0.0
- numpy >= 2.2.6, < 2.2.7  
- pyvista >= 0.45.2, < 1.0.0
- scipy >= 1.15.3, < 1.15.4
- numba >= 0.61.2, < 1.0.0
- cmocean >= 4.0.3, < 5.0.0
- tqdm >= 4.67.1, < 5.0.0
- webp >= 0.4.0, < 1.0.0
- PySide6 >= 6.9.1, < 7.0.0

**Development Dependencies:**
- codecov, black, codespell, pre-commit, build, twine, PyInstaller, setuptools, wheel

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
- **Abstract references**: When referring to abstractions, use lowercase and separate individual words with a space (e.g. "an airplane's wings are used to generate lift" and"the cross section of a wing typically has a streamlined shape known as an airfoil"). This is to distinguish them from code objects.

### Docstring Style
- Follow existing PteraSoftware docstring conventions  
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

### Variable Naming
- Use descriptive variable names that clearly indicate their purpose  
- Use lowercase with underscores for variable names
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the documentation

#### Coordinate System Variable Naming Patterns

Variables must follow one of these three patterns based on their requirements:

1. **Axes without a point and without a frame**  
   `[variable name]_[axes ID]`
   - Example: `force_W` (force in wind axes)

2. **Axes without a point and with a frame**  
   `[variable name]_[axes ID]__[frame ID]` (note double underscore)
   - Example: `velocity_W__B` (velocity in wind axes, observed from body frame)

3. **Axes with a point and without a frame**  
   `[variable name]_[axes ID]_[point ID]`
   - Example: `position_G_Cg` (position in geometry axes, relative to CG point)

#### Standardized IDs

**Axis Systems:**
- E: Earth, B: body, W: wind, G: geometry, Wn: wing, Wcs: wing cross section, Wcsp: wing cross section parent, A: airfoil

**Reference Points:**  
- I: simulation starting point, Cgi: starting point, Cg: CG point, Ler: leading edge root point, Lp: leading point, Lpp: leading point parent

**Reference Frames:**
- E: Earth frame, B: body frame, Wn: wing frame, Wcs: wing cross section frame, Wcsp: wing cross section parent frame

#### Angle Vector Naming

**Passive angle vectors (relating two axis systems):**
- `[variable name]_[source axes ID]_to_[target axes ID]_[sequence ID]`
- Example: `angles_E_to_B_izyx` (angles from Earth to body axes using intrinsic z-y'-x" sequence)

**Active angle vectors (for rotation within current axis system):**
- `[variable name]_act_[sequence ID]`
- Example: `rotation_angles_act_izyx` (rotation angles using intrinsic z-y'-x" sequence)

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
- Example: `rotation_angles_act_izyx = np.array([0.0, 45.0, 0.0])`

**Homogeneous Coordinates:**
- For 4x4 transformation matrices, vectors are converted to homogeneous form
- Position vectors: `[x, y, z, 1]` (has_point=True)
- Direction vectors: `[x, y, z, 0]` (has_point=False)

## Coordinate System Conventions

### Axes, Points, and Frames Overview

Ptera Software simulates flapping-wing dynamics and aerodynamics using several different axis systems, reference points, and reference frames. The notation and terminology used is an extended version of that introduced in "Flight Vehicle Aerodynamics" by Mark Drela.

**Key Concepts:**
- **Axis system** ("axes"): Contains information about three cartesian directions. They do not have an origin, therefore it's incorrect to write something like "the wing axes' origin."
- **Reference points** ("points"): Contains information about the location of a particular point in space.
- **Reference frame** ("frame"): Contains information about the location of an "observer" and their motion relative to what is observed.

**Vector Requirements:**
- Force vectors: require only an axis system
- Position vectors: require both axes and a reference point (origin)
- Moment vectors: require axes and reference point (point about which moment acts)
- Velocity/acceleration vectors: require both axis system and reference frame

### Variable Naming Patterns

There are three useful combinations of axes, points, and frames. Variables are denoted by appending information using **IDs**, and described in text using **names**:

1. **Axes without a point and without a frame**  
   `[variable name]_[axes ID]`  
   "[variable name] (in [axes name])"

2. **Axes without a point and with a frame**  
   `[variable name]_[axes ID]__[frame ID]`  
   "[variable name] (in [axes name], observed from the [frame name])"  

3. **Axes with a point and without a frame**  
   `[variable name]_[axes ID]_[point ID]`  
   "[variable name] (in [axes name], relative to the [point name])"

**Naming Convention:** IDs move in scope from most to least specific. Names move from least to most specific.

### ID Abbreviations and Names

- E: Earth  
- B: body
- P: airplane  
- W: wind
- Pr: problem
- G: geometry  
- Wn: wing  
- Wcs: wing cross section  
- Wcsp: wing cross section parent  
- A: airfoil
- I: simulation starting point  
- Cgi: starting point  
- Cg: CG point  
- Ler: leading edge root point  
- Lp: leading point  
- Lpp: leading point parent  
- ...pp...: panel point (Fr=front right, Fl=front left, Bl=back left, Br=back right, C=collocation)
- ...bhvp...: bound horseshoe vortex point  
- ...brvp...: bound ring vortex point  
- ...wrvp...: wake ring vortex point  
- ...whvp...: wake horseshoe vortex point  
- ...lvp...: line vortex point (S=start, E=end, C=center)

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

**5. Wing Axes**
- Basis: Back of Wing in first cross section plane, normal to plane, top surface
- Right-handed for non-symmetric/symmetric-continuous Wings, left-handed for mirror-only Wings
- Variables: `..._Wn...` (local), `..._Wn1...` (Airplane-local), `..._Wn2P1...` (non-local)
- Text: "...in wing axes...", "...in the first Wing's axes...", "...in the first Airplane's second Wing's axes..."

**6. Wing Cross Section Axes**
- Basis: Trailing edge in plane, normal to plane, top surface
- Handedness same as wing axes
- Variables: `..._Wcs...` (local), `..._Wcs1...` (Wing-local), `..._Wcs3Wn2...` (Airplane-local), `..._Wcs1Wn2P1...` (non-local)
- Text: "...in wing cross section axes...", "...in the first WingCrossSection's axes...", etc.

**7. Wing Cross Section Parent Axes**
- Basis: Identical to wing axes for first cross section, identical to previous cross section axes for subsequent ones
- Variables: `..._Wcsp...` with similar local/non-local patterns
- Text: "...in wing cross section parent axes..."

**8. Airfoil Axes**
- Basis: Chordwise to trailing point, Normal to chord toward upper line (2D)
- Variables: `..._A...` (local), `..._AWcs2...` (Wing-local), etc.
- Text: "...in airfoil axes...", "...in the second WingCrossSection's Airfoil's axes...", etc.

### Angle Vectors and Transformations

**Rotation Sequences (Tait-Bryan only):**
- ixyz, ixzy, iyxz, iyzx, izxy, izyx: intrinsic sequences (x-y'-z", x-z'-y", etc.)
- exyz, exzy, eyxz, eyzx, ezxy, ezyx: extrinsic sequences (x-y-z, x-z-y, etc.)

**Passive Angle Vectors (relating axis systems):**
- Variable pattern: `[name]_[source axes ID]_to_[target axes ID]_[sequence ID]`
- Text pattern: "[name] describing the orientation of [target axes name] relative to [source axes name] using an [sequence name] sequence"
- Example: `angles_E_to_B_izyx` = "angles describing the orientation of body axes relative to Earth axes using an intrinsic z-y'-x" sequence"

**Active Angle Vectors (rotating within axis system):**
- Variable pattern: `[name]_[sequence ID]`
- Text pattern: "[name] for rotation using an [sequence name] sequence"
- Example: `angles_izyx` = "angles for rotation using an intrinsic z-y'-x" sequence"

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

## Example Usage Pattern

Basic usage follows this pattern:
```python
import pterasoftware as ps

# Define geometry using the new package structure
airplane = ps.geometry.airplane.Airplane(
    wings=[
        ps.geometry.wing.Wing(
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.wing_cross_section.WingCrossSection(
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                    num_spanwise_panels=8,
                ),
                ps.geometry.wing_cross_section.WingCrossSection(
                    Lp_Wcsp_Lpp=[0.0, 5.0, 0.0],
                    airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
                    num_spanwise_panels=None,  # Tip cross section
                ),
            ],
        ),
    ],
)

# Define operating conditions
operating_point = ps.operating_point.OperatingPoint()

# Create problem
problem = ps.problems.SteadyProblem(
    airplanes=[airplane], operating_point=operating_point
)

# Choose and run solver
solver = ps.steady_horseshoe_vortex_lattice_method.SteadyHorseshoeVortexLatticeMethodSolver(
    steady_problem=problem
)
solver.run()

# Visualize results
ps.output.draw(solver=solver, scalar_type="lift", show_streamlines=True)
```

## Miscellaneous Guidelines
- Use clear, descriptive variable names  
- Avoid abbreviations unless they are well-known in the context  
- In docstrings and comments, never use em-dashes (—) or en-dashes (–); always use hyphens (-) for clarity  
- In docstrings and comments, never use a multiplication sign (×); always use a lowercase x (x)
- **Coordinate and axis references**: When referring to axes, coordinates, or planes, use  lowercase letters with hyphens between coordinate letters and descriptors (e.g., "x-axis", "y-coordinate", "xz-plane", "z-direction"). However, never add hyphens between a word and "axis" or "axes" (e.g., "wing axis" not "wing-axis", "body axes" not "body-axes").  Never use uppercase letters for axis references in text.
- Never use emojis in code, comments, or docstrings
- Always use straight single and double quotes, not curly ones
- In .md files, use "…" instead of "...". In all other files, use "...".