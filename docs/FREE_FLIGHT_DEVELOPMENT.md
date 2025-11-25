# Free Flight Development Tracking

This document tracks the development progress and goals for the free flight simulation feature in Ptera Software. Free flight simulations couple aerodynamic analysis (via Ptera Software's UVLM solver) with rigid body dynamics (via MuJoCo) to enable fully dynamic flapping-wing flight simulations.

## Feature Branch
- **Branch name**: `feature/free_flight`
- **Base branch**: `main`

## Current Status Summary

**🟢 DEVELOPMENT IN PROGRESS - Core 6-DOF Dynamics Now Functional**

**What Works:**
- ✅ Complete coupled solver infrastructure implemented and functional
- ✅ MuJoCo integration with proper load passing and state extraction
- ✅ Prescribed motion steps work correctly and match baseline simulations
- ✅ Solver runs through full coupling loop without crashing (300+ steps tested)
- ✅ Comprehensive debugging infrastructure with systematic test progression
- ✅ Free flight visualization tools (animate_free_flight, draw support)
- ✅ Converged, trimmed, stable test case (simple glider validated with XFLR5)
- ✅ MuJoCo rotation matrix convention understood and properly transformed
- ✅ Euler angle extraction from rotation matrix (with gimbal lock handling)
- ✅ Rotational velocity terms properly implemented with body-to-geometry transformation

**What Needs Testing/Validation:**
- ⚠️ Full 6-DOF dynamics need validation against known solutions
- ⚠️ Non-zero sideslip (beta) cases need testing
- ⚠️ Animation coordinate transformations may still have issues

**Immediate Focus:**
Validating the complete 6-DOF free flight simulation against known solutions and debugging any remaining coordinate transformation issues in visualization.

## Overview

Free flight mode allows Ptera Software to simulate flapping-wing vehicles where the aerodynamic forces feed back into the vehicle's motion, enabling simulation of natural flight dynamics. This requires bidirectional coupling:
1. Ptera Software computes aerodynamic forces and moments
2. Forces/moments are passed to MuJoCo for rigid body dynamics integration
3. MuJoCo returns updated vehicle state (position, orientation, velocities)
4. Ptera Software updates geometry and operating conditions for the next time step

## Architecture

### Coupled Class Hierarchy

A parallel class hierarchy has been created with "Coupled" prefixes to distinguish free flight classes from standard prescribed-motion classes:

```
CoupledOperatingPoint (extends OperatingPoint)
    └─ Adds: `omegas_E_to_B__E` (angular speeds) and `angles_E_to_B_izyx` (orientation)

CoupledMovement (parallel to Movement)
    └─ Contains: AirplaneMovements + initial CoupledOperatingPoint
    └─ Generates: time-varying CoupledOperatingPoints (populated by solver)

CoupledSteadyProblem (parallel to SteadyProblem)
    └─ Pairs: Airplanes + CoupledOperatingPoint for a single time step

CoupledUnsteadyProblem (parallel to UnsteadyProblem)
    └─ Contains: CoupledMovement
    └─ Generates: list of CoupledSteadyProblems

CoupledUnsteadyRingVortexLatticeMethodSolver (parallel to UnsteadyRingVortexLatticeMethodSolver)
    └─ Takes: CoupledUnsteadyProblem + mujoco_model
    └─ Implements: bidirectional coupling loop
```

### Key Architectural Differences from Standard UVLM

**Standard UVLM:**
- Geometry initialized once before time loop
- Operating point varies according to prescribed OperatingPointMovement
- No feedback from aerodynamic loads to motion

**Coupled UVLM:**
- Geometry reinitialized at EACH time step (based on updated state from MuJoCo)
- Operating point updated at each time step based on actual vehicle motion
- Aerodynamic loads directly influence vehicle dynamics

### Solver Run Loop Structure

The `CoupledUnsteadyRingVortexLatticeMethodSolver.run()` method follows this structure for each time step:

1. **Initialize geometric arrays** - Allocate arrays for current time step
2. **`_initialize_panel_vortices()`** - Set up bound RingVortices on panels
3. **`_collapse_geometry()`** - Flatten geometry matrices to 1D arrays
4. **`_calculate_wing_wing_influences()`** - Compute wing-wing influence coefficient matrix
5. **`_calculate_freestream_wing_influences()`** - Compute freestream-wing influences
6. **`_calculate_wake_wing_influences()`** - Compute wake-wing influences
7. **`_calculate_vortex_strengths()`** - Solve for bound vortex circulations
8. **`_calculate_loads()`** - Compute aerodynamic forces and moments
9. **`_pass_loads_to_mujoco()`** - Transform loads to Earth axes and pass to MuJoCo
10. **`mujoco_model.step()`** - Advance MuJoCo simulation by one time step
11. **`_process_new_states_from_mujoco()`** - Receive updated state from MuJoCo
12. **`_create_next_coupled_steady_problem()`** - Create next time step's problem
13. **`_populate_next_wake()`** - Shed wake vortices for next time step

## Implementation Status

### ✅ Completed

1. **Class Structure** (pterasoftware/operating_point.py, pterasoftware/problems.py, pterasoftware/movements/movement.py)
   - CoupledOperatingPoint class with angular velocity support
   - CoupledMovement class for free flight motion definition
   - CoupledSteadyProblem class for time step characterization
   - CoupledUnsteadyProblem class for overall problem definition

2. **Solver Framework** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - CoupledUnsteadyRingVortexLatticeMethodSolver class created
   - Initialization method with mujoco_model parameter and validation
   - Main run() method structure established
   - Per-time-step geometric reinitialization implemented
   - All standard UVLM calculation methods copied and functional

3. **MuJoCo Interface** (pterasoftware/mujoco_model.py)
   - MuJoCoModel class created to wrap MuJoCo model and data structures
   - Support for loading models from XML strings or file paths
   - Methods for applying loads: `apply_loads(forces_E, moments_E_Cg)`
   - Methods for advancing simulation: `step()`
   - Methods for extracting state: `get_state()` returns `position_E_E`, `R_pas_E_to_B`, `velocity_E__E`, `omegas_B__E`, time
   - Methods for resetting: `reset()`
   - Proper parameter validation and error handling
   - Follows Ptera Software coordinate system naming conventions

4. **Package Integration** (`pterasoftware/__init__.py`)
   - Coupled solver module imported and exposed
   - MuJoCo model module imported and exposed

5. **Reference Materials** (gammabot_simulations/, mujoco_examples/)
   - Non-free-flight GammaBot simulation example
   - MuJoCo tutorial notebooks and examples
   - Previous MuJoCo integration reference files

6. **Bug Fixes**
   - Fixed uninitialized `self.current_airplanes` in `CoupledUnsteadyRingVortexLatticeMethodSolver.run()`
   - Fixed missing `self.` assignment for `initial_coupled_operating_point` in `CoupledMovement.__init__()`

7. **Load Passing to MuJoCo** (`pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - `_pass_loads_to_mujoco()` method implemented
   - Transforms aerodynamic loads from wind axes to Earth axes
   - Transformation chain: W_CgP1 > GP1_CgP1 > BP1_CgP1 > E_CgP1
   - Uses `_transformations` module for coordinate transformations
   - Applies loads to MuJoCo via `mujoco_model.apply_loads(forces_E, moments_E_Cg)`
   - Comprehensive docstring documenting transformation approach

8. **Processing MuJoCo State Updates** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - `_process_new_states_from_mujoco()` method implemented
   - Receives updated state from MuJoCo (position, orientation, velocity, angular velocity)
   - Computes speed and freestream velocity from vehicle motion (assuming still air)
   - Transforms freestream velocity: E > BP1 > body axes for alpha/beta computation
   - Calculates angle of attack (alpha) and sideslip angle (beta) from body-axes freestream
   - Creates new CoupledOperatingPoint with updated aerodynamic state
   - Stores intermediate state (_next_position_E_E, _next_R_pas_E_to_B, etc.) for next method
   - Comprehensive docstring with transformation chain documentation

9. **Creating Next Time Step's Problem** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - `_create_next_coupled_steady_problem()` method implemented
   - Retrieves prescribed Airplane geometry (correct flapping Wing positions)
   - Extracts Euler angles from MuJoCo rotation matrix (z-y'-x" intrinsic sequence)
   - Creates new Airplane with updated position/orientation from MuJoCo
   - Preserves Wing geometry from prescribed motion (flapping kinematics)
   - Retrieves CoupledOperatingPoint created by _process_new_states_from_mujoco()
   - Creates and appends new CoupledSteadyProblem for next time step
   - Enforces `Cg_GP1_CgP1 = (0.0, 0.0, 0.0)` constraint (first Airplane's CG defines origin)
   - TODO: Verify wake handling (wake RingVortices have absolute GP1_CgP1 coordinates)

10. **Initial Free Flight Test Case** (examples/free_flight_simple_glider.py)
   - Simple gliding rectangular wing test case created
   - 4m wingspan (2m per side), 0.2m chord, NACA 0012 airfoil
   - 0.5 kg mass, 5° initial angle of attack, 10 m/s initial speed
   - 100 time steps at 0.01s each (1 second total simulation)
   - **Solver completes successfully** - coupling loop verified to work end-to-end
   - MuJoCo warnings indicate numerical instability (loads/initial conditions need tuning)
   - Demonstrates complete workflow: aerodynamics → MuJoCo → state extraction → next problem

11. **Bug Fixes During Initial Testing**
   - Fixed einsum signature in `_calculate_freestream_wing_influences()`
     - Changed from `"ij,j->i"` to `"ij,ij->i"` to handle rotational velocity component
     - `np.cross` produces `(num_panels, 3)` array, not `(3,)` vector
   - Fixed first Airplane CG position constraint in `_create_next_coupled_steady_problem()`
     - Enforced `Cg_GP1_CgP1 = (0.0, 0.0, 0.0)` by design
     - Motion captured in orientation and velocities, not absolute position
   - Fixed ring vortex initialization timing
     - First step: Initialize in main loop with conditional check
     - Subsequent steps: Initialize after creating next problem, before populating wake
     - Prevents `_populate_next_wake()` from accessing uninitialized ring_vortices

12. **Bug Fixes During Debugging Phase**
   - Fixed incorrect initial orientation in MuJoCo model by computing proper quaternion from rotation matrix (`coupled_unsteady_ring_vortex_lattice_method.py`, `4_simple_glider_free_flight.py`)
   - Fixed wake RingVortex coordinate transformations in `animate_free_flight()` (`output.py`)
   - Added external forces (weight + external wind axes x-force) to MuJoCo loads (`coupled_unsteady_ring_vortex_lattice_method.py`)
   - Modified solver to iterate several prescribed motion steps before free flight to prevent numerical issues from initial wake shedding spikes (`coupled_unsteady_ring_vortex_lattice_method.py`)
   - Updated `draw()` function to work with coupled solver (`output.py`)
   - Added state history arrays (`stackPosition_E_E`, `stackR_pas_E_to_BP1`) for visualization (`coupled_unsteady_ring_vortex_lattice_method.py`)
   - Added `initial_key_frame_name` parameter to MuJoCoModel for proper initialization (`mujoco_model.py`)
   - Improved MuJoCo rotation matrix extraction (semi-fixed, still brittle) (`mujoco_model.py`)

13. **Free Flight Visualization** (`pterasoftware/output.py`)
   - `animate_free_flight()` function implemented
   - Shows Airplane trajectory through 3D space with proper Earth axes transformations
   - Displays wake RingVortices in correct reference frame
   - Supports scalar coloring (lift, induced drag, side force)
   - Optional trajectory path visualization
   - Helper functions: `_get_panel_surfaces_free_flight()`, `_get_wake_ring_vortex_surfaces_free_flight()`
   - Updated `draw()` to support CoupledUnsteadyRingVortexLatticeMethodSolver

14. **Debugging Infrastructure** (`debugging scripts/`)
   - Systematic debugging case progression created:
     - `1_simple_glider_convergence.py`: Convergence study
     - `2_simple_glider_trim.py`: Trim analysis
     - `3_simple_glider_prescribed.py`: Prescribed motion validation
     - `4_simple_glider_free_flight.py`: Free flight simulation
   - Single-panel flat plate debugging cases:
     - `flat_plate_prescribed.py`: 1-Panel prescribed motion baseline
     - `flat_plate_free_flight.py`: 1-Panel free flight for isolation testing
   - DEBUGGING_LOG.md tracks active issues, likely culprits, and resolved problems
   - XFLR5 model (`simple_glider.xfl`) for validation with known stable configuration
   - Simple glider design verified: converged, trimmed, statically stable in pitch/yaw, realistic inertia

15. **MuJoCo Rotation Matrix Convention Resolved** (pterasoftware/mujoco_model.py)
   - Clarified that MuJoCo's `xmat` is `R_pas_B_to_E` (transforms vectors from body to Earth axes)
   - Fixed `get_state()` to properly compute `R_pas_E_to_B = R_pas_B_to_E.T`
   - Removed janky empirical transformation logic
   - Added clear documentation in docstring explaining `xmat` convention

16. **Euler Angle Extraction Implemented** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - Implemented proper extraction for intrinsic z-y'-x" (izyx) sequence from R_pas_E_to_B
   - Added gimbal lock detection and handling (when pitch near ±90°)
   - Documented matrix element formulas in code comments
   - REFACTOR note added for future extraction into `_transformations.py`

17. **Rotational Velocity Terms Re-enabled** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - Fixed `_calculate_freestream_wing_influences()` to properly include omega cross r terms
   - Added coordinate transformation from body axes to geometry axes for angular velocity
   - Body to geometry transformation: `omegas_GP1__E = omegas_BP1__E * [-1, 1, -1]` (180° rotation about Y)
   - Properly converts angular velocity from degrees/s to radians/s for cross product
   - Added detailed comments explaining the transformation and physics

18. **Code Quality Improvements**
   - Fixed sign convention for beta in initial orientation calculation (debugging scripts)
   - Added REFACTOR notes for functions that should be moved to `_transformations.py`
   - Improved docstrings in `MuJoCoModel` with explicit `xfrc_applied` and state variable documentation
   - Standardized use of `np.rad2deg()` and `np.deg2rad()` instead of manual `* 180 / pi`

### ⚠️ In Progress / Known Issues

**Current Status:**

The core 6-DOF dynamics infrastructure is now complete. All previously critical issues (rotation matrix convention, Euler angle extraction, rotational velocity terms) have been resolved. The focus shifts to validation and testing.

**⚠️ Needs Validation/Testing:**

1. **Full 6-DOF Dynamics Validation** (HIGH PRIORITY)
   - Rotational velocity terms are now enabled but need validation against known solutions
   - Euler angle extraction needs verification with various initial orientations
   - Energy conservation should be checked for gliding simulations

2. **Non-Zero Sideslip Cases** (MEDIUM PRIORITY)
   - Current implementation verified for alpha > 0, beta = 0
   - Need to test with non-zero sideslip angles
   - May reveal remaining coordinate transformation issues

3. **Wake Coordinate Handling in Free Flight** (MEDIUM PRIORITY)
   - Wake RingVortices stored in Wing have absolute coordinates (GP1_CgP1)
   - When Airplane position/orientation change, wake may need coordinate updates
   - Alternatively, consider storing wake in solver instead of Wing for coupled simulations
   - Lower priority since visualization fixes suggest current approach may work

**🔍 Lower Priority Issues:**

4. **Animation Coordinate Issues**
   - Airplane appears rotated 180° about body y-axis in `animate_free_flight()`
   - Airplane appears to move in opposite direction from expected
   - Both likely visualization transformation issues since wake convects correctly in `draw()`

5. **Airplane Mesh Disappears in Animation**
   - After several tens of time steps in `animate_free_flight()`, Airplane mesh vanishes
   - May be rendering issue or coordinate transformation problem

**📋 Future Issues (Post-Validation):**

6. **No MuJoCo Model Visualization**
   - Cannot easily verify MuJoCo XML setup is correct
   - Consider integration with MuJoCo's built-in visualizer

7. **Airplane Direction Reversal**
   - If Airplane switches direction, wake sheds from wrong edge
   - Detection exists but no recovery mechanism
   - Low priority for gliding simulations

8. **MuJoCo Model Not Programmatically Generated**
    - MuJoCo XML created separately from CoupledUnsteadyProblem
    - Changes to problem don't automatically update MuJoCo model
    - Consider programmatic model generation from problem definition

9. **Utility Functions Need Extraction** (REFACTOR)
   - `R_to_quat_wxyz` function duplicated in debugging scripts; should be in `_transformations.py`
   - Euler angle extraction from rotation matrix should be in `_transformations.py`
   - Alpha/beta extraction from velocity vector should be in `_transformations.py`
   - Body-to-geometry axes transformation should have dedicated function

### ❌ Not Yet Started

1. **Coordinate Transformation Utilities**
   - ✅ Transformation of forces/moments to MuJoCo's expected frames (completed in `_pass_loads_to_mujoco()`)
   - ✅ Transformation of MuJoCo state to compute alpha/beta (completed in `_process_new_states_from_mujoco()`)
   - ✅ Extraction of Euler angles from rotation matrix (completed in `_create_next_coupled_steady_problem()`)
   - ⚠️ Wake coordinate transformation (might need, requires investigation)

2. **Testing Infrastructure**
   - Unit tests for CoupledOperatingPoint
   - Unit tests for CoupledMovement
   - Unit tests for CoupledSteadyProblem
   - Unit tests for CoupledUnsteadyProblem
   - Integration tests for coupled solver
   - Validation against known free flight data

3. **Example Scripts**
   - Free flight GammaBot simulation example
   - Documentation on how to set up free flight simulations
   - Tutorial notebook demonstrating workflow

4. **Convergence and Stability Analysis**
   - Investigate time step size requirements for stable coupling
   - Study convergence behavior with different panel counts
   - Analyze effect of prescribed vs. free wake on free flight

## Key Design Decisions

### Coordinate System Handling

- All solver computations use the first Airplane's geometry axes and CG (GP1_CgP1)
- Forces and moments calculated in GP1_CgP1, then transformed to wind axes for output
- MuJoCo uses Earth axes (E) and Earth frame (E), corresponding to MuJoCo's "world coordinates"
- Earth origin point (E) defined as the initial CG position (documented in AXES_POINTS_AND_FRAMES.md)
- MuJoCo interface requires transformation from GP1_CgP1 (Ptera Software) to E (MuJoCo)

### Time Stepping

- Aerodynamic time step (delta_time in CoupledMovement) must match or subdivide MuJoCo time step
- Current approach: Ptera Software drives the time loop, calling MuJoCo for each step
- Alternative: Let MuJoCo drive, calling Ptera Software as needed (not yet explored)

### Wake Model

- Prescribed wake assumed for initial implementation (prescribed_wake=True)
- Free wake may be explored later for improved accuracy
- Wake vortices shed based on current geometry, propagate via freestream + induced velocity

## Current Goals

### Immediate (Next Steps)

**Priority 1: Rotational Motion and Coordinate Transformations (COMPLETED ✅)**

All previously critical blockers have been resolved:

1. ✅ **MuJoCo Rotation Matrix Convention** (COMPLETED)
   - Clarified that MuJoCo's `xmat` is `R_pas_B_to_E`
   - Implemented proper transformation: `R_pas_E_to_B = R_pas_B_to_E.T`
   - Removed empirical transformation logic

2. ✅ **Euler Angle Extraction** (COMPLETED)
   - Implemented proper extraction for izyx sequence with gimbal lock handling
   - Verified formulas against rotation matrix structure

3. ✅ **Rotational Velocity Terms** (COMPLETED)
   - Re-enabled omega cross r terms in `_calculate_freestream_wing_influences()`
   - Added proper body-to-geometry coordinate transformation for angular velocity
   - Documented physics and sign conventions in code comments

**Priority 2: Validation and Verification (CURRENT FOCUS)**

4. ✅ **Create Simple Free Flight Test Case** (COMPLETED)
   - ✅ Implemented simple gliding wing test case (`debugging scripts/4_simple_glider_free_flight.py`)
   - ✅ Set up MuJoCo model with converged, trimmed, stable configuration
   - ✅ Verified solver runs through complete coupling loop (300 steps)
   - ✅ Created systematic debugging progression (convergence → trim → prescribed → free)
   - ✅ Created single-panel flat plate isolation tests

5. ✅ **Develop Visualization for Free Flight Trajectories** (PARTIALLY COMPLETED)
   - ✅ Implemented `animate_free_flight()` function with trajectory support
   - ✅ Updated `draw()` to work with coupled solver
   - ✅ Basic plotting for debugging (position over time)
   - ⚠️ Animation has coordinate transformation issues (lower priority)
   - 🔲 Still need: comprehensive time-series plots (forces, moments, alpha, beta, omegas)

6. **Validate Full 6-DOF Dynamics** (IN PROGRESS)
   - Verify force/moment calculations produce physically sensible results
   - Check energy conservation (kinetic + potential should be approximately constant for gliding)
   - Compare simple glider free flight trajectory with XFLR5 dynamic stability predictions
   - Tune time step size for stability/accuracy trade-off
   - Validate wake convection is physically reasonable
   - Test with non-zero sideslip (beta) cases

### Medium-term

1. Implement free flight GammaBot simulation
2. Validate against experimental GammaBot data
3. Create comprehensive documentation and examples
4. Add unit and integration tests

### Long-term

1. Support multiple Airplanes in free flight (formation flight dynamics)

## Related Files

### Core Implementation
- `pterasoftware/operating_point.py` - CoupledOperatingPoint class
- `pterasoftware/problems.py` - CoupledSteadyProblem, CoupledUnsteadyProblem
- `pterasoftware/movements/movement.py` - CoupledMovement class
- `pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py` - Main solver
- `pterasoftware/mujoco_model.py` - MuJoCoModel class for MuJoCo interface

### Reference Materials
- `gammabot_simulations/GammaBot.py` - Example of standard (non-free-flight) unsteady simulation
- `mujoco_examples/` - MuJoCo tutorials and examples
- `docs/AXES_POINTS_AND_FRAMES.md` - Coordinate system documentation
- `docs/ANGLE_VECTORS_AND_TRANSFORMATIONS.md` - Transformation documentation

### Testing
- Tests to be added in `tests/unit/` and `tests/integration/`

## Notes and TODOs

### From Code Comments

- **`CoupledOperatingPoint.__init__()`**:
  - REFACTOR: Add section to ANGLE_VECTORS_AND_TRANSFORMATIONS.md about angular speeds

- **`CoupledMovement.__init__()`**:
  - FIXME: Automatic delta_time calculation gives poor results at high Strouhal numbers
  - Consider improving time step calculation for flapping-dominated motion

- **`CoupledUnsteadyRingVortexLatticeMethodSolver.run()`**:
  - TODO: Initialization steps at start of time loop may be redundant, consider optimizing

### Open Questions

1. How to handle control surface actuation in free flight (coupled with controller)?
2. Should gravity be included in CoupledOperatingPoint or handled entirely by MuJoCo?
3. How to best visualize free flight trajectories?

## Debugging Insights and Lessons Learned

### Key Insights from Debugging Process

**1. Importance of Validated Test Cases**
- Creating a simple glider with known converged, trimmed, and stable parameters (validated via XFLR5) was crucial
- Having an external reference (XFLR5 model with confirmed static stability) provided confidence that the test case itself was physically reasonable
- This prevented wasting time debugging a potentially unstable or poorly configured test case

**2. Value of Systematic Debugging Progression**
- Single-panel flat plate tests provided the simplest possible test case for comparing prescribed vs. free flight
- This methodical approach revealed that prescribed motion works correctly, narrowing the problem to free flight-specific logic

**3. Prescribed Motion Steps as Numerical Stabilizer**
- Initial wake shedding causes large force spikes that can destabilize MuJoCo integration
- Running several prescribed motion steps first allows wake to develop before enabling coupling
- This workaround is documented in `coupled_unsteady_ring_vortex_lattice_method.py`

**4. Coordinate Transformation Complexity**
- The interface between MuJoCo's coordinate conventions and Ptera Software's extensive coordinate system formalism is surprisingly difficult to get right
- Understanding the theoretical basis of MuJoCo's rotation matrices is essential before proceeding

**5. Iterative Debugging with Temporary Workarounds**
- Temporarily disabling rotational velocity terms (multiplying by zero) was a valid debugging strategy
- This isolated the problem and confirmed that rotational motion coupling is the remaining issue
- However, these workarounds must be documented clearly and eventually removed

**6. Importance of State History for Visualization**
- Adding `stackPosition_E_E` and `stackR_pas_E_to_BP1` arrays enabled simple debugging plots
- Being able to plot position and orientation over time was invaluable for understanding behavior
- Need to expand this to include forces, moments, and aerodynamic angles for complete debugging picture

### Debugging Methodology That Worked Well

1. **Create minimal reproducible case** (single-panel flat plate)
2. **Establish baseline** (prescribed motion that works correctly)
3. **Compare baseline to coupled case** (identify where they diverge)
4. **Isolate problematic terms** (disable rotational terms to identify culprit)
5. **Document findings** (DEBUGGING_LOG.md tracks progress and hypotheses)
6. **Validate assumptions** (XFLR5 model confirms test case is reasonable)

### Next Steps Based on Insights

The coordinate transformation issues identified during debugging have now been resolved:

**✅ Resolved Issues:**
- MuJoCo's `xmat` is `R_pas_B_to_E`; we take the transpose to get `R_pas_E_to_B`
- Euler angle extraction uses the izyx sequence with proper gimbal lock handling
- Rotational velocity cross products now use the correct body-to-geometry transformation

**Current Focus:**
The infrastructure is now complete for full 6-DOF free flight simulation. The next steps are:
1. Validate the simulation produces physically sensible results
2. Test with various initial conditions including non-zero sideslip
3. Compare against XFLR5 dynamic stability predictions
4. Extract utility functions into `_transformations.py` for reuse and testing
