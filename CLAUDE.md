# Ptera Software Development Guidelines for Claude

## Project Overview
Ptera Software is a fast, easy-to-use, and open-source package for analyzing flapping-wing flight using unsteady and steady vortex lattice methods.

### Key Features
- **Multiple Simulation Methods**: Steady horseshoe VLM, steady ring VLM, and unsteady ring UVLM
- **Customizable Aircraft Geometry**: Multi-wing aircraft with arbitrary wing cross sections and airfoils
- **Time-Dependent Motion**: Custom prescribed flapping motions
- **Formation Flight**: Multi-airplane simulations supported since v2.0.0
- **High-Speed Computing**: JIT compilation via Numba for fast simulations
- **Trim Analysis**: Automatic search for trim operating points
- **Convergence Analysis**: Automatic search for converged parameters
- **GUI Interface**: Basic GUI available (beta stage)
- **Visualization Tools**: 3D mesh visualization and 2D plotting of results
- **Extensive Testing**: Comprehensive unit and integration tests for reliability

### Python Version Constraint
Requires Python 3.11, but active development is done in 3.13  

### Key Runtime Dependencies
- **NumPy/SciPy**: Core numerical computations
- **Numba**: JIT compilation for performance-critical loops
- **PyVista**: 3D mesh processing and visualization
- **PySide6**: GUI framework
- **Matplotlib**: 2D plotting and analysis output

## Architecture Overview

### Relevant Directories, Packages, and Files
- `.github/`: Directory with GitHub configuration files
    - `ISSUE_TEMPLATE/`: Directory with issue templates
      - `bug_report.md`
      - `feature_request.md`
    - `workflows/`: Directory with GitHub Actions workflows
      - `black.yml`
      - `codespell.yml`
      - `mypy.yml`
      - `tests.yml`
    - `pull_request_template.md`
- `.venv/`: Directory for the Python virtual environment (not included in version control)
- `experimental/`: Directory with experimental scripts and prototypes (not included in version control)
- `docs/`: Directory with documentation files
  - `examples expected output/`: Example output files for verification
  - `ANGLES_VECTORS_AND_TRANSFORMATIONS.md`: Conventions and definitions for angle vectors and transformations **READ BEFORE CONTRIBUTING ANY CODE, PARTICIPATING IN DISCUSSIONS REGARDING, OR PLANNING RELATED TO VECTOR-VALUED VARIABLES**
  - `AXES_POINTS_AND_FRAMES.md`: Conventions and definitions for axis systems, points, and reference points: **READ BEFORE CONTRIBUTING ANY CODE, PARTICIPATING IN DISCUSSIONS REGARDING, OR PLANNING RELATED TO VECTOR-VALUED VARIABLES**
  - `CODE_STYLE.md`: Code style guidelines: **READ BEFORE CONTRIBUTING ANY CODE**
  - `RUNNING_TESTS_AND_TYPE_CHECKS.md`: Instructions for running tests and type checks **READ BEFORE RUNNING TESTS OR TYPE CHECKS LOCALLY**
  - `TYPE_HINT_AND_DOCSTRING_STYLE.md`: Guidelines for type hinting and docstring formatting: **READ BEFORE CONTRIBUTING ANY CODE OR WRITING ANY DOCSTRINGS**
  - `WRITING_STYLE.md`: Guidelines for writing style in comments and documentation: **READ BEFORE WRITING ANY DOCUMENTATION, DOCSTRINGS, OR COMMENTS**
- `examples/`: Directory with example scripts for users
- `gui/`: Directory with GUI source code
- `pterasoftware/`: Main package with modular solver architecture
  - `geometry/`: Package with aircraft geometry classes
    - `_airfoils/` Directory containing data files with airfoil coordinates
    - `_meshing.py`: Wing mesh generation
    - `airfoil.py`: Airfoil class with coordinate generation
    - `airplane.py`: Airplane class with coordinate transformations
    - `wing.py`: Wing class with symmetry processing
    - `wing_cross_section.py`: WingCrossSection class with validation
  - `movements/` Package with movement classes (definitions for time-dependent motion)
    - `_functions.py`: Movement utility functions
    - `airplane_movement.py`: Airplane motion definitions
    - `movement.py`: Core Movement class
    - `operating_point_movement.py`: Operating condition changes
    - `wing_cross_section_movement.py`: Wing cross section motion
    - `wing_movement.py`: Wing flapping motion
  - `_aerodynamics.py`: Vortex elements and velocity calculations
  - `_functions.py`: Shared utility functions
  - `_logging.py`: Contains function for setting up logging
  - `_panel.py`: Panel class for discretized mesh elements
  - `_parameter_validation.py`: Input validation functions
  - `_transformations.py` Coordinate transformations and rotations
  - `convergence.py`: Convergence analysis tools
  - `operating_point.py` OperatingPoint and CoupledOperatingPoint classes
  - `output.py`: Visualization and results processing
  - `problems.py` SteadyProblem and UnsteadyProblem classes
  - `steady_horseshoe_vortex_lattice_method.py`: Steady horseshoe VLM solver
  - `steady_ring_vortex_lattice_method.py`: Steady ring VLM solver
  - `trim.py`: Trim analysis functionality
  - `unsteady_ring_vortex_lattice_method.py`: Unsteady ring UVLM solver
- `tests/`: Directory with unit and integration tests
  - `integration/`: Integration tests for combined functionality
    - `fixtures/`: Fixtures for integration tests
      - `airplane_fixtures.py`
      - `movement_fixtures.py`
      - `operating_point_fixtures.py`
      - `problem_fixtures.py`
      - `solver_fixtures.py`
    - `test_output.py`
    - `test_steady_convergence.py`
    - `test_steady_horseshoe_vortex_lattice_method.py`
    - `test_steady_ring_vortex_lattice_method.py`
    - `test_steady_trim.py`
    - `test_unsteady_convergence.py`
    - `test_unsteady_ring_vortex_lattice_method_multiple_wing_static_geometry.py`
    - `test_unsteady_ring_vortex_lattice_method_multiple_wing_variable_geometry.py`
    - `test_unsteady_ring_vortex_lattice_method_static_geometry.py`
    - `test_unsteady_ring_vortex_lattice_method_variable_geometry.py`
  - `unit/`: Unit tests for individual classes and functions
    - `fixtures/`: Fixtures for unit tests
      - `aerodynamics_function_fixtures.py`
      - `airplane_movement_fixtures.py`
      - `geometry_fixtures.py`
      - `horseshoe_fixtures.py`
      - `movement_fixtures.py`
      - `movement_function_fixtures.py`
      - `operating_point_fixtures.py`
      - `operating_point_movement_fixtures.py`
      - `problem_fixtures.py`
      - `ring_vortex_fixtures.py`
      - `wing_cross_section_movement_fixtures.py`
      - `wing_movement_fixtures.py`
    - `test_aerodynamics_functions.py`
    - `test_airfoil.py`
    - `test_airplane.py`
    - `test_airplane_movement.py`
    - `test_horseshoe_vortex.py`
    - `test_movement.py`
    - `test_movement_functions.py`
    - `test_operating_point.py`
    - `test_operating_point_movement.py`
    - `test_panel.py`
    - `test_problems.py`
    - `test_ring_vortex.py`
    - `test_transformations.py`
    - `test_wing.py`
    - `test_wing_cross_section.py`
    - `test_wing_cross_section_movement.py`
    - `test_wing_movement.py`
- `.codespell-ignore`: File listing words to ignore in spell checking
- `.gitignore`: Git ignore file
- `.pre-commit-config.yaml`: Pre-commit configuration file
- `BUILD.md`: Instructions for building the GUI
- `CONTRIBUTING.md`: Contribution guidelines for developers
- `make_installer.iss`: Inno Setup script for building Windows installer
- `MANIFEST.in`: Manifest file for packaging
- `mypy.ini`: MyPy configuration file
- `pterasoftware.spec`: PyInstaller specification file for building executables
- `pyproject.toml`: Project configuration file
- `README.md`: Project overview and installation instructions for developers
- `requirements.txt`: Full list of runtime dependencies with version constraints
- `requirements_dev.txt`: Full list of development dependencies with version constraints
- `setup.cfg`: Setup configuration file

## Running Scripts That Import Ptera Software
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