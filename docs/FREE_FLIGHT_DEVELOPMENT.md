# Free Flight Development Tracking

This document tracks the development progress and goals for the free flight simulation feature in Ptera Software. Free flight simulations couple aerodynamic analysis (via Ptera Software's UVLM solver) with rigid body dynamics (via MuJoCo) to enable fully dynamic flapping-wing flight simulations.

## Feature Branch
- **Branch name**: `feature/free_flight`
- **Base branch**: `main`

## Current Status Summary (as of commit 331ec61)

**🟡 DEVELOPMENT IN PROGRESS - Core Infrastructure Complete, Debugging Active Issues**

**What Works:**
- ✅ Complete coupled solver infrastructure implemented and functional
- ✅ MuJoCo integration with proper load passing and state extraction
- ✅ Prescribed motion steps work correctly and match baseline simulations
- ✅ Solver runs through full coupling loop without crashing (300+ steps tested)
- ✅ Comprehensive debugging infrastructure with systematic test progression
- ✅ Free flight visualization tools (animate_free_flight, draw support)
- ✅ Converged, trimmed, stable test case (simple glider validated with XFLR5)

**What Doesn't Work Yet:**
- ❌ Rotational velocity terms disabled (set to zero to maintain stability)
- ❌ Euler angle extraction disabled (hardcoded to zero)
- ❌ Rotation matrix transformation from MuJoCo is brittle and empirical
- ❌ Free flight portion exhibits non-physical behavior when rotational terms enabled

**Immediate Focus:**
Understanding and fixing MuJoCo's rotation matrix convention and coordinate transformations to enable proper 6-DOF rotational dynamics. This is the primary blocker for realistic free flight simulation.

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
    └─ Adds: `omegas_E_to_B__E` (angular speeds)

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
   - Methods for extracting state: `get_state()` returns `position_E_E`, `R_pas_E_to_B`, `velocity_E__E`, `omegas_E_to_B__E`, time
   - Methods for resetting: `reset(qpos, qvel)`
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
   - Fixed uninitialized `self.current_airplanes` in `CoupledUnsteadyRingVortexLatticeMethodSolver.run()` (line 222)
   - Fixed missing `self.` assignment for `initial_coupled_operating_point` in `CoupledMovement.__init__()` (`movement.py`:413)

7. **Load Passing to MuJoCo** (`pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - `_pass_loads_to_mujoco()` method implemented (line 1158-1233)
   - Transforms aerodynamic loads from wind axes to Earth axes
   - Transformation chain: W_CgP1 > GP1_CgP1 > BP1_CgP1 > E_CgP1
   - Uses `_transformations` module for coordinate transformations
   - Applies loads to MuJoCo via `mujoco_model.apply_loads(forces_E, moments_E_Cg)`
   - Comprehensive docstring documenting transformation approach

8. **Processing MuJoCo State Updates** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - `_process_new_states_from_mujoco()` method implemented (line 1238-1313)
   - Receives updated state from MuJoCo (position, orientation, velocity, angular velocity)
   - Computes speed and freestream velocity from vehicle motion (assuming still air)
   - Transforms freestream velocity: E > BP1 > body axes for alpha/beta computation
   - Calculates angle of attack (alpha) and sideslip angle (beta) from body-axes freestream
   - Creates new CoupledOperatingPoint with updated aerodynamic state
   - Stores intermediate state (_next_position_E_E, _next_R_pas_E_to_B, etc.) for next method
   - Comprehensive docstring with transformation chain documentation

9. **Creating Next Time Step's Problem** (pterasoftware/coupled_unsteady_ring_vortex_lattice_method.py)
   - `_create_next_coupled_steady_problem()` method implemented (line 1341-1429)
   - Retrieves prescribed Airplane geometry (correct flapping Wing positions)
   - Extracts Euler angles from MuJoCo rotation matrix (z-y'-x" intrinsic sequence)
   - Creates new Airplane with updated position/orientation from MuJoCo
   - Preserves Wing geometry from prescribed motion (flapping kinematics)
   - Retrieves CoupledOperatingPoint created by _process_new_states_from_mujoco()
   - Creates and appends new CoupledSteadyProblem for next time step
   - Enforces `Cg_E_CgP1 = (0.0, 0.0, 0.0)` constraint (first Airplane's CG defines origin)
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
   - Fixed einsum signature in `_calculate_freestream_wing_influences()` (line 757)
     - Changed from `"ij,j->i"` to `"ij,ij->i"` to handle rotational velocity component
     - `np.cross` produces `(num_panels, 3)` array, not `(3,)` vector
   - Fixed first Airplane CG position constraint in `_create_next_coupled_steady_problem()`
     - Enforced `Cg_E_CgP1 = (0.0, 0.0, 0.0)` by design (line 1392)
     - Motion captured in orientation and velocities, not absolute position
   - Fixed ring vortex initialization timing (lines 439-502)
     - First step: Initialize in main loop with conditional check
     - Subsequent steps: Initialize after creating next problem, before populating wake
     - Prevents `_populate_next_wake()` from accessing uninitialized ring_vortices

12. **Bug Fixes During Debugging Phase** (multiple commits: 99d147d..331ec61)
   - Fixed incorrect initial orientation in MuJoCo model by computing proper quaternion from rotation matrix (`coupled_unsteady_ring_vortex_lattice_method.py`:75-78, `4_simple_glider_free_flight.py`:16-48)
   - Fixed wake RingVortex coordinate transformations in `animate_free_flight()` (`output.py`:753-761)
   - Added external forces (weight + external wind axes x-force) to MuJoCo loads (`coupled_unsteady_ring_vortex_lattice_method.py`:1216-1218, 1260-1263)
   - Modified solver to iterate several prescribed motion steps before free flight to prevent numerical issues from initial wake shedding spikes (`coupled_unsteady_ring_vortex_lattice_method.py`:1272-1277)
   - Updated `draw()` function to work with coupled solver (`output.py`:217-229)
   - Added state history arrays (`stackCgP1_E_E`, `stackR_pas_E_to_BP1`) for visualization (`coupled_unsteady_ring_vortex_lattice_method.py`:184-185, 1310-1311)
   - Added `initial_key_frame_name` parameter to MuJoCoModel for proper initialization (`mujoco_model.py`:43, 65-69, 98-114)
   - Improved MuJoCo rotation matrix extraction (semi-fixed, still brittle) (`mujoco_model.py`:192-205)

13. **Free Flight Visualization** (`pterasoftware/output.py`)
   - `animate_free_flight()` function implemented (line 667-995)
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

### ⚠️ In Progress / Known Issues

**Current Debugging Status (as of commit 331ec61):**

The coupled solver is now stable for prescribed motion steps and has been extensively debugged. The free flight portion completes without crashing but still exhibits non-physical behavior. Key achievements and remaining issues:

**✅ Recent Progress:**
- Prescribed and free flight simulations now match closely during prescribed time steps
- Eliminated wild instability by temporarily disabling omega-related velocity terms
- Successfully created converged, trimmed, stable simple glider based on XFLR5 validation
- Fixed initial orientation quaternion generation
- Improved alpha/beta calculation from MuJoCo state (works for alpha > 0, beta = 0, angles_E_to_B_izyx = (0, 0, 0))

**⚠️ Active Debugging Issues:**

1. **Rotational Motion Velocity Terms Disabled** (CRITICAL)
   - All velocity components due to angular velocities (`omegas_E_to_B__E`) are currently multiplied by zero
   - Lines affected: `coupled_unsteady_ring_vortex_lattice_method.py`:784, 1068, 1071, 1074
   - This is non-physical but necessary to maintain stability
   - Hints that remaining bugs exist in how rotational motion affects local velocities
   - Need to investigate: cross product sign conventions, reference frame transformations

2. **Euler Angle Extraction Disabled** (CRITICAL)
   - Euler angle extraction in `_create_next_coupled_steady_problem()` is hardcoded to zero
   - Line 1408-1410: `angleY = 0.0`, `angleX = 0.0`, `angleZ = 0.0`
   - Original extraction formulas commented out (lines 1405-1407)
   - This prevents proper orientation updates from MuJoCo
   - Related to issue with rotation matrix transformations

3. **Rotation Matrix Transformation is Brittle** (HIGH PRIORITY)
   - `MuJoCoModel.get_state()` uses janky transformation logic (`mujoco_model.py`:192-205)
   - Only verified to work for: alpha > 0, beta = 0, angles_E_to_B_izyx = (0, 0, 0)
   - Unclear what `xmat[self.body_id].reshape(3, 3)` actually represents in MuJoCo
   - Need to understand MuJoCo's rotation matrix convention and map it properly to Ptera Software's coordinate systems
   - Current approach uses empirical matrix transformations without theoretical justification

4. **Wake Coordinate Handling in Free Flight** (MEDIUM PRIORITY)
   - Wake RingVortices stored in Wing have absolute coordinates (GP1_CgP1)
   - When Airplane position/orientation change, wake may need coordinate updates
   - Alternatively, consider storing wake in solver instead of Wing for coupled simulations
   - Lower priority since visualization fixes suggest current approach may work

**🔍 Less-Likely Issues (Lower Priority):**

5. **Animation Coordinate Issues**
   - Airplane appears rotated 180° about body y-axis in `animate_free_flight()`
   - Airplane appears to move in opposite direction from expected
   - Both likely visualization transformation issues since wake convects correctly in `draw()`

6. **Initial Orientation Parameter Confusion**
   - `airplane.angles_E_to_B_izyx` parameter is confusing and may not affect anything
   - Currently forced to zero; unclear if it should exist at all
   - Consider removing and handling all orientation offsets via OperatingPoint

7. **Airplane Mesh Disappears in Animation**
   - After several tens of time steps in `animate_free_flight()`, Airplane mesh vanishes
   - May be rendering issue or coordinate transformation problem

**📋 Future Issues (Post-Stability):**

8. **No MuJoCo Model Visualization**
   - Cannot easily verify MuJoCo XML setup is correct
   - Consider integration with MuJoCo's built-in visualizer

9. **Airplane Direction Reversal**
   - If Airplane switches direction, wake sheds from wrong edge
   - Detection exists but no recovery mechanism
   - Low priority for gliding simulations

10. **MuJoCo Model Not Programmatically Generated**
    - MuJoCo XML created separately from CoupledUnsteadyProblem
    - Changes to problem don't automatically update MuJoCo model
    - Consider programmatic model generation from problem definition

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

**Priority 1: Fix Rotational Motion and Coordinate Transformations (CRITICAL)**

The current blocker for realistic free flight is the disabled rotational velocity terms and Euler angle extraction. These must be fixed to enable proper 6-DOF dynamics:

1. **Understand and Fix MuJoCo Rotation Matrix Convention**
   - Research MuJoCo's `xmat` documentation to understand what rotation it represents
   - Determine proper transformation from MuJoCo's rotation matrix to Ptera Software's R_pas_E_to_B
   - Test with non-zero beta and non-zero angles_E_to_B_izyx cases
   - Replace empirical transformation in `mujoco_model.py`:192-205 with theoretically justified approach

2. **Re-enable and Fix Euler Angle Extraction**
   - Once rotation matrix issue is resolved, re-enable angle extraction in `_create_next_coupled_steady_problem()`
   - Verify extracted angles produce correct Airplane orientations

3. **Re-enable and Fix Rotational Velocity Terms**
   - Investigate cross product sign conventions for omega-induced velocities
   - Verify reference frame transformations for rotational motion effects
   - Gradually re-enable terms (remove 0.0 multipliers) at lines 784, 1068, 1071, 1074
   - Test stability with single-panel flat plate case first, then simple glider

**Priority 2: Validation and Verification**

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

6. **Debug and Validate Once Core Issues Resolved**
   - Verify force/moment calculations produce physically sensible results
   - Check energy conservation (kinetic + potential should be approximately constant for gliding)
   - Compare simple glider free flight trajectory with XFLR5 dynamic stability predictions
   - Tune time step size for stability/accuracy trade-off
   - Validate wake convection is physically reasonable

### Medium-term

1. Implement free flight GammaBot simulation
2. Validate against experimental GammaBot data
3. Create comprehensive documentation and examples
4. Add unit and integration tests

### Long-term

1. Support multiple Airplanes in free flight (formation flight dynamics)

## Related Files

### Core Implementation
- `pterasoftware/operating_point.py` - CoupledOperatingPoint class (line 225+)
- `pterasoftware/problems.py` - CoupledSteadyProblem (line 207+), CoupledUnsteadyProblem (line 286+)
- `pterasoftware/movements/movement.py` - CoupledMovement class (line 306+)
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

- **`CoupledOperatingPoint.__init__()`** (operating_point.py:86):
  - REFACTOR: Add section to ANGLE_VECTORS_AND_TRANSFORMATIONS.md about angular speeds

- **`CoupledSteadyProblem.__init__()`** (problems.py:244-250):
  - REFACTOR: Consider how to handle reference points with multiple Airplanes in free flight
  - Suggestion: Disallow free flight with multiple Airplanes initially

- **`CoupledMovement.__init__()`** (movements/movement.py:419-425):
  - FIXME: Automatic delta_time calculation gives poor results at high Strouhal numbers
  - Consider improving time step calculation for flapping-dominated motion

- **`CoupledUnsteadyRingVortexLatticeMethodSolver.run()`** (line 325-326):
  - TODO: Initialization steps at start of time loop may be redundant, consider optimizing

### Open Questions

1. Should MuJoCo handle collision detection between Airplane components?
2. How to handle control surface actuation in free flight (coupled with controller)?
3. Should gravity be included in CoupledOperatingPoint or handled entirely by MuJoCo?
4. How to best visualize free flight trajectories?

## Debugging Insights and Lessons Learned

### Key Insights from Debugging Process (commits 99d147d..331ec61)

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
- This workaround is documented in `coupled_unsteady_ring_vortex_lattice_method.py`:1272-1277

**4. Coordinate Transformation Complexity**
- The interface between MuJoCo's coordinate conventions and Ptera Software's extensive coordinate system formalism is surprisingly difficult to get right
- Understanding the theoretical basis of MuJoCo's rotation matrices is essential before proceeding

**5. Iterative Debugging with Temporary Workarounds**
- Temporarily disabling rotational velocity terms (multiplying by zero) was a valid debugging strategy
- This isolated the problem and confirmed that rotational motion coupling is the remaining issue
- However, these workarounds must be documented clearly and eventually removed

**6. Importance of State History for Visualization**
- Adding `stackCgP1_E_E` and `stackR_pas_E_to_BP1` arrays enabled simple debugging plots
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

The debugging process has clarified that the core issue is **coordinate transformation between MuJoCo and Ptera Software**, specifically:
- Understanding what rotation matrix MuJoCo's `xmat` represents
- Properly mapping that to Ptera Software's `R_pas_E_to_B` convention
- Ensuring rotational velocity cross products use correct sign conventions and reference frames

Once this is resolved, the rest of the infrastructure is in place for full 6-DOF free flight simulation.
