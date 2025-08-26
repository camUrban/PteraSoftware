# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Ptera Software Development Guidelines for Claude

## Before Starting, Please Read
README.md  
docs\AXES_COORDINATES_AND_FRAMES.md  
docs\ANGLE_VECTORS_AND_TRANSFORMATIONS.md  

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
