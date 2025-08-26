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
  - **Geometry**: `geometry.py` (Airplane → Wing → WingCrossSection → Airfoil hierarchy)
  - **Solvers**: Three VLM implementations (steady horseshoe, steady ring, unsteady ring)
  - **Problems**: `problems.py` (SteadyProblem, UnsteadyProblem classes)
  - **Support**: meshing, movement, output, trim, convergence modules

### Key Features
- **Multiple Simulation Methods**: Steady horseshoe VLM, steady ring VLM, unsteady ring UVLM
- **Customizable Aircraft Geometry**: Multi-wing aircraft with arbitrary cross sections and airfoils
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
- **"Ptera Software"**: When writing as text, always write as two words without a 
hyphen, each being capitalized (never "ptera", "ptera software", or "PteraSoftware"). 
When writing as a package to be imported, use ```import pterasoftware as ps```  
- **"cross section"**: Always write as two words, never hyphenated 
(not "cross-section")  
- **Object references**: When referring to code objects, use proper class naming 
convention. The capitalization indicates that we are talking about a code object, not an 
abstraction. You don't need to add "object" or "objects" after the class name since the 
capitalization already makes this clear (e.g. "update the Wings" instead of "update the 
Wing objects"). In summary, when talking about code objects:
  - ✅ "the previous WingCrossSection"
  - ❌ "the previous cross section"
  - ✅ "this Wing"
  - ❌ "this wing"
  - ✅ "update the Wings"
  - ❌ "update the Wing objects" (unnecessary)
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

## Coordinate System Conventions

### Axes, Points, and Frames Overview

Ptera Software simulates flapping-wing dynamics and aerodynamics using several different axis systems, reference points, and reference frames. The notation and terminology used is an extended version of that introduced in "Flight Vehicle Aerodynamics" by Mark Drela.

**Key Concepts:**
- **Axis system** ("axes"): Contains information about three directions (cartesian, polar, or spherical)
- **Reference points** ("points"): Contains information about the location of a particular point in space
- **Reference frame** ("frame"): Contains information about the location of an "observer" and their motion relative to what is observed

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

**Axis Systems:**
- E: Earth  
- B: body  
- P: airplane  
- W: wind  
- G: geometry  
- Wn: wing  
- Wcs: wing cross section  
- Wcsp: wing cross section's parent  
- A: airfoil  

**Reference Points:**
- I: Simulation's starting point  
- Cgi: starting point  
- Cg: CG point  
- Ler: leading edge's root point  
- Lp: leading point  
- Lpp: leading point's parent  
- ...pp...: panel point (Fr=front right, Fl=front left, Bl=back left, Br=back right, C=collocation)
- ...bhvp...: bound horseshoe vortex's point  
- ...brvp...: bound ring vortex point  
- ...wrvp...: wake ring vortex point  
- ...whvp...: wake horseshoe vortex point  
- ...lvp...: line vortex's point (S=start, E=end, C=center)

**Reference Frames:**
- E: Earth's reference frame (inertial)
- B: Body's reference frame (non-inertial, attached to airplane body)
- Wn: Wing's reference frame (non-inertial, attached to wing leading edge root)
- Wcs: Wing cross section's reference frame (non-inertial, attached to leading point)
- Wcsp: Wing cross section's parent reference frame (non-inertial, attached to leading parent point)

### Axis System Definitions

**1. Earth Axes**
- Basis: North, East, Down (right-handed)
- Variables: `..._E...`
- Text: "...in Earth axes..."

**2. Body Axes**
- Basis: Front of airplane, Right of airplane, Bottom of airplane (right-handed)
- Variables: `..._B...` (local) or `..._BP1...` (non-local)
- Text: "...in body axes..." or "...in the first airplane's body axes..."

**3. Wind Axes**
- Basis: Parallel to freestream velocity, perpendicular directions defined via angle of attack/sideslip
- Variables: `..._W...` (local) or `..._WP1...` (non-local)
- Text: "...in wind axes..." or "...in the first airplane's wind axes..."

**4. Geometry Axes**
- Basis: Back of airplane, Right of airplane, Top of airplane (right-handed)
- Variables: `..._G...` (local) or `..._GP1...` (non-local)
- Text: "...in geometry axes..." or "...in the first airplane's geometry axes..."

**5. Wing Axes**
- Basis: Back of wing in first cross section plane, normal to plane, top surface
- Right-handed for non-symmetric/symmetric-continuous wings, left-handed for mirror-only wings
- Variables: `..._Wn...` (local), `..._Wn1...` (airplane-local), `..._Wn2P1...` (non-local)
- Text: "...in wing axes...", "...in the first wing's axes...", "...in the first airplane's second wing's axes..."

**6. Wing Cross Section Axes**
- Basis: Trailing edge in plane, normal to plane, top surface
- Handedness same as wing axes
- Variables: `..._Wcs...` (local), `..._Wcs1...` (wing-local), `..._Wcs3Wn2...` (airplane-local), `..._Wcs1Wn2P1...` (non-local)
- Text: "...in wing cross section axes...", "...in the first wing cross section's axes...", etc.

**7. Wing Cross Section Parent Axes**
- Basis: Identical to wing axes for first cross section, identical to previous cross section axes for subsequent ones
- Variables: `..._Wcsp...` with similar local/non-local patterns
- Text: "...in the wing cross section's parent axes..."

**8. Airfoil Axes**
- Basis: Chordwise to trailing point, Normal to chord toward upper line (2D)
- Variables: `..._A...` (local), `..._AWcs2...` (wing-local), etc.
- Text: "...in airfoil axes...", "...in the second wing cross section's airfoil's axes...", etc.

### Angle Vectors and Transformations

**Rotation Sequences (Tait-Bryan only):**
- i123, i132, i213, i231, i312, i321: intrinsic sequences
- e123, e132, e213, e231, e312, e321: extrinsic sequences

**Passive Angle Vectors (relating axis systems):**
- Variable pattern: `[name]_pas_[source axes ID]_to_[target axes ID]_[sequence ID]`
- Text pattern: "[name] describing the orientation of [target axes name] relative to [source axes name] using an [sequence name] sequence"
- Example: `angles_pas_E_to_B_i321` = "angles describing the orientation of the Earth axes relative to the body axes using an intrinsic 3-2'-1" sequence"

**Active Angle Vectors (rotating within axis system):**
- Variable pattern: `[name]_act_[sequence ID]`
- Text pattern: "[name] for rotation using an [sequence name] sequence"
- Example: `angles_act_i321` = "angles for rotation using an intrinsic 3-2'-1" sequence"

**Rotation and Transformation Matrices:**
- **Passive rotation matrices:** `R_pas_[source]_to_[target]` (3×3)
- **Passive transformation matrices:** `T_pas_[source axes]_[source point]_to_[target axes]_[target point]` (4×4 homogeneous)
- **Active rotation matrices:** `[name]_R_act` (3×3)  
- **Active transformation matrices:** `[name]_T_act` (4×4 homogeneous)

**Implementation Notes:**
- All angles in degrees unless noted otherwise
- Angle wrapping to (-180, 180] range
- Vectors treated as column vectors with left-multiplication by matrices
- Different sequences have different gimbal lock singularities

## Miscellaneous Guidelines
- Use clear, descriptive variable names  
- Avoid abbreviations unless they are well-known in the context  
- In docstrings and comments, never use em-dashes (—) or en-dashes (–); always use 
hyphens (-) for clarity  
- **Coordinate and axis references**: When referring to axes, coordinates, or planes, use 
lowercase letters with hyphens between coordinate letters and descriptors (e.g., "x-axis", 
"y-coordinate", "xz-plane", "z-direction"). However, never add hyphens between a word and 
"axis" or "axes" (e.g., "wing axis" not "wing-axis", "body axes" not "body-axes"). 
Never use uppercase letters for axis references in text.
- Never use emojis in code, comments, or docstrings

## Example Usage Pattern

Basic usage follows this pattern:
```python
import pterasoftware as ps

# Define geometry
airplane = ps.geometry.Airplane(
    wings=[
        ps.geometry.Wing(
            symmetric=True,
            wing_cross_sections=[
                ps.geometry.WingCrossSection(
                    airfoil=ps.geometry.Airfoil(name="naca2412"),
                ),
                ps.geometry.WingCrossSection(
                    y_le=5.0,
                    airfoil=ps.geometry.Airfoil(name="naca2412"),
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
