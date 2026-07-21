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
- **Aeroelasticity**: First-order structural wing deformation coupled to the UVLM via a torsional spring-mass-damper model (beta)
- **Free Flight**: Six-degree-of-freedom flight dynamics from UVLM aerodynamics coupled to MuJoCo rigid-body dynamics (beta)
- **Visualization Tools**: 3D mesh visualization and 2D plotting of results
- **Save and Load**: JSON serialization of solved simulations without pickle security risks
- **Extensive Testing**: Comprehensive unit and integration tests for reliability

### Python Version Constraint

Requires Python 3.11, but active development is done in 3.13

### Key Runtime Dependencies

- **NumPy/SciPy**: Core numerical computations
- **Numba**: JIT compilation for performance-critical loops
- **PyVista**: 3D mesh processing and visualization
- **Matplotlib**: 2D plotting and analysis output

## Architecture Overview

### Relevant Directories, Packages, and Files

- `.github/`: Directory with GitHub configuration files: the issue and pull request templates, the label definitions (`labels.yml`), the Dependabot, funding, and code-owner configuration, and the GitHub Actions workflows (one per pre-commit hook or hook group, plus `label-sync.yml`, `publish.yml`, and `tests.yml`)
- `.venv/`: Directory for the Python virtual environment, configured for the host machine's OS (not included in version control)
- `.venv-wsl/`: Directory for the Python virtual environment configured for a WSL OS (not included in version control, may be missing if host machine doesn't use WSL for development)
- `experimental/`: Directory with experimental scripts and prototypes (not included in version control)
- `docs/`: Directory with documentation files
    - `examples_expected_output/`: Example output files for verification
    - `hero_graphics/`: Assets for the README hero graphic
    - `private/`: Directory with documentation not included in this repository's version control (may be missing if the private repo hasn't been cloned and linked to this local repo)
        - `katz_plotkin_12_2/`: A recreation of Chapter 12.2, which describes efficiently including the effects of symmetry and ground effect for vortex lattice methods, from the textbook "Low-Speed Aerodynamics" by Katz and Plotkin
        - `katz_plotkin_13_12/`: A recreation of Chapter 13.12, which describes the UVLM, from the textbook "Low-Speed Aerodynamics" by Katz and Plotkin
        - `katz_plotkin_d/`: A recreation of Appendix D, which includes example Fortran programs, from the textbook "Low-Speed Aerodynamics" by Katz and Plotkin
        - `lambert_2015_2_3__2_4/`: A recreation of Sections 2.3 and 2.4 from Thomas Lambert's thesis "Modeling of aerodynamic forces in flapping flight with the unsteady vortex lattice method"
    - `website/`: Directory with the source files for generating the documentation website
    - `ANGLE_VECTORS_AND_TRANSFORMATIONS.md`: Conventions and definitions for angle vectors and transformations **READ BEFORE CONTRIBUTING ANY CODE, PARTICIPATING IN DISCUSSIONS REGARDING, OR PLANNING RELATED TO VECTOR-VALUED VARIABLES**
    - `AXES_POINTS_AND_FRAMES.md`: Conventions and definitions for axis systems, points, and reference points: **READ BEFORE CONTRIBUTING ANY CODE, PARTICIPATING IN DISCUSSIONS REGARDING, OR PLANNING RELATED TO VECTOR-VALUED VARIABLES**
    - `CLASSES_AND_IMMUTABILITY.md`: Description of class structure and attribute immutability.
    - `CODE_STYLE.md`: Code style guidelines: **READ BEFORE CONTRIBUTING ANY CODE**
    - `MUJOCO_CONVENTIONS.md`: Definitive interpretation of MuJoCo state variables and their mapping to Ptera Software's axes, points, frames, and transformations
    - `RUNNING_TESTS_AND_TYPE_CHECKS.md`: Instructions for running tests and type checks **READ BEFORE RUNNING TESTS OR TYPE CHECKS LOCALLY**
    - `STRONG_COUPLING.md`: Mathematical framework for the strongly coupled free-flight UVLM-MuJoCo solver: the fixed-point sub-iteration, Aitken relaxation, the weighting matrix, and the convergence tolerances
    - `TYPE_HINT_AND_DOCSTRING_STYLE.md`: Guidelines for type hinting and docstring formatting: **READ BEFORE CONTRIBUTING ANY CODE OR WRITING ANY DOCSTRINGS**
    - `WRITING_STYLE.md`: Guidelines for writing style in comments and documentation: **READ BEFORE WRITING ANY DOCUMENTATION, DOCSTRINGS, OR COMMENTS**
- `examples/`: Directory with example scripts for users
- `pterasoftware/`: Main package with modular solver architecture
    - `geometry/`: Package with aircraft geometry classes
        - `_airfoils/`: Directory containing data files with airfoil coordinates
        - `_meshing.py`: Wing mesh generation
        - `airfoil.py`: Airfoil class with coordinate generation
        - `airplane.py`: Airplane class with coordinate transformations
        - `wing.py`: Wing class with symmetry processing
        - `wing_cross_section.py`: WingCrossSection class with validation
    - `movements/`: Package with movement classes (definitions for time-dependent motion)
        - `aeroelastic_airplane_movement.py`: AeroelasticAirplaneMovement class
        - `aeroelastic_movement.py`: AeroelasticMovement class
        - `aeroelastic_wing_cross_section_movement.py`: AeroelasticWingCrossSectionMovement class
        - `aeroelastic_wing_movement.py`: AeroelasticWingMovement class
        - `airplane_movement.py`: AirplaneMovement class
        - `free_flight_movement.py`: FreeFlightMovement class
        - `free_flight_operating_point_movement.py`: FreeFlightOperatingPointMovement class
        - `movement.py`: Movement class
        - `operating_point_movement.py`: OperatingPointMovement class
        - `wing_cross_section_movement.py`: WingCrossSectionMovement class
        - `wing_movement.py`: WingMovement class
    - `_colormap_data/`: Directory containing data files with the vendored color map and color palette colors, along with their licenses
    - `_aerodynamics_functions.py`: Induced velocity functions
    - `_colormaps.py`: Color maps and color palettes used by the visualization functions
    - `_convergence_cache.py`: JSON solve and memo cache for convergence analysis
    - `_convergence_meshing.py`: Mesh building and refinement for convergence iterations
    - `_core.py`: Core classes for the movement and problem hierarchies
    - `_coupled_unsteady_ring_vortex_lattice_method.py`: Coupled unsteady UVLM solver subclass with step-by-step geometry
    - `_fixed_point_relaxation.py`: Pure fixed-point relaxation helpers (weighted norm, convergence test, Aitken relaxation factor) for the strong-coupling sub-iteration
    - `_functions.py`: Shared utility functions
    - `_logging.py`: Contains function for setting up logging
    - `_mujoco_model.py`: Contains the MuJoCoModel class.
    - `_oscillation.py`: Oscillation functions for movement classes
    - `_panel.py`: Panel class for discretized mesh elements
    - `_parameter_validation.py`: Input validation functions
    - `_serialization.py`: JSON serialization and deserialization (save/load)
    - `_transformations.py`: Coordinate transformations and rotations
    - `aeroelastic_unsteady_ring_vortex_lattice_method.py`: Aeroelastic UVLM solver subclass with first-order structural deformation
    - `convergence.py`: Convergence analysis tools
    - `free_flight_unsteady_ring_vortex_lattice_method.py`: Free flight UVLM solver subclass with six-DOF MuJoCo coupling
    - `operating_point.py`: OperatingPoint class
    - `output.py`: Visualization and results processing
    - `problems.py`: SteadyProblem and UnsteadyProblem classes
    - `steady_horseshoe_vortex_lattice_method.py`: Steady horseshoe VLM solver
    - `steady_ring_vortex_lattice_method.py`: Steady ring VLM solver
    - `trim.py`: Trim analysis functionality
    - `unsteady_ring_vortex_lattice_method.py`: Unsteady ring UVLM solver
- `scripts/`: Directory with maintenance and tooling scripts
    - `hero_generation/`: Scripts for creating and finalizing the README hero graphic
        - `create_solve_and_save_hero.py`: Creates, solves, and saves the hero simulation
        - `finalize_and_save_hero.py`: Renames preview hero graphics to their permanent names
        - `load_and_visualize_hero.py`: Loads the saved hero simulation and generates preview graphics
    - `analyze_webp.py`: Renders WebP frames to PNG files for inspection (backs the `analyze-webp` slash command)
    - `check_ascii_only.py`: Pre-commit hook script that flags non-ASCII characters in text files
    - `find_unused_fixtures.py`: Finds and optionally deletes unused fixtures and dead `setUp` attributes across the test suite (backs the `delete-unused-fixtures` slash command)
    - `regenerate_example_outputs.py`: Runs all example scripts (or a single named example) and collects their outputs into `docs/examples_expected_output/`, re-rendering oversized WebP files at lower quality
- `tests/`: Directory with unit and integration tests
    - `integration/`: Integration tests for combined functionality, with shared fixtures in a `fixtures/` subpackage. There is one test module per solver configuration worth exercising end-to-end (each solver, plus its surface-effect, wake-truncation, multiple-wing, and variable-geometry variants), along with modules for convergence, trim, output, and serialized output
    - `unit/`: Unit tests for individual classes and functions, with shared fixtures in a `fixtures/` subpackage. The test modules mirror the package's modules one-to-one (for example, `test_wing.py` tests `geometry/wing.py`), plus a few suite-level modules such as `test_package_init.py`, `test_slots.py`, and `test_test_environment.py`. Fixture modules are named after the test modules they serve
    - `_test_environment.py`: Configures the test process to quiet known sources of test run noise (the serialization dirty-tree warnings, tqdm progress bars, and the headless-Linux VTK warning). It is imported as the first line of `tests/__init__.py` so the suppressions are in place before any pterasoftware, pyvista, or tqdm module loads
- `validation/`: Directory with the experimental validation study: `validation_study.py` simulates a flapping test stand from Yeo et al., 2011 and compares the UVLM results against the paper's published pressure measurements, which are stored alongside it as CSV files extracted from the paper, along with the accompanying report (`validation_paper.pdf`)
- `.codespell-ignore.txt`: File listing words to ignore in spell checking
- `.gitignore`: Git ignore file
- `.pre-commit-config.yaml`: Pre-commit configuration file
- `.readthedocs.yaml`: Read the Docs build configuration
- `CITATION.cff`: Citation metadata for the project (powers GitHub's cite-this-repository feature and the Zenodo release record)
- `codecov.yml`: Codecov configuration for test coverage reporting
- `CONTRIBUTING.md`: Contribution guidelines for developers
- `MANIFEST.in`: Manifest file for packaging
- `mypy.ini`: MyPy configuration file
- `pyproject.toml`: Project configuration file
- `README.md`: Project overview and installation instructions for developers
- `requirements.txt`: Full list of runtime dependencies with version constraints
- `requirements_dev.txt`: Full list of development dependencies with version constraints
- `requirements_min.txt`: Minimum-version runtime dependencies
- `setup.cfg`: Setup configuration file

## Git Command Permissions

The project's `.claude/settings.json` restricts which shell commands Claude Code may run. For git specifically, the subcommands fall into three tiers:

- Denied (cannot run at all): `git am`, `git apply`, `git bisect`, `git branch`, `git cherry-pick`, `git clean`, `git clone`, `git config`, `git fetch`, `git init`, `git lfs`, `git merge`, `git mv`, `git pull`, `git push`, `git rebase`, `git remote`, `git reset`, `git revert`, `git rm`, `git sparse-checkout`, `git stash`, `git submodule`, `git switch`, `git tag`, and `git worktree`.
- Ask (prompt for approval each time): `git add`, `git checkout`, `git commit`, and `git restore`.
- Allowed (run without a prompt): `git diff`, `git grep`, `git log`, `git ls-files`, `git show`, and `git status`.

The practical consequences are worth internalizing before planning any git workflow:

- `git push` is denied, so Claude cannot push. When a workflow needs the branch on the remote (for example, before opening a pull request), the user must push it manually.
- `git branch` is denied even for read-only listing, so use `git status` or `git log` to determine the current branch instead.
- A compound command fails if any single part is denied. For example, `git branch --show-current && git status` is rejected because `git branch` is denied, even though `git status` on its own is allowed. Split such commands into separately allowed invocations.

## Common Mistakes

- Forgetting to read RUNNING_TESTS_AND_TYPE_CHECKS.md before running tests and trying to use pytest (Ptera Software uses unittest)
- Forgetting to read CODE_STYLE.md before contributing code
- Forgetting to read TYPE_HINT_AND_DOCSTRING_STYLE.md before writing docstrings
- Forgetting to read ANGLE_VECTORS_AND_TRANSFORMATIONS.md and AXES_POINTS_AND_FRAMES.md before working with vector-valued variables. If in doubt, before writing code, read both of these documents
- Forgetting to read WRITING_STYLE.md before writing documentation, docstrings, or comments
