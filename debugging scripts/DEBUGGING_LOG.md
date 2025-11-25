# Debugging Log for Free Flight Simple Glider Simulation

## Description of Problem
The free flight simple glider simulation is unstable and the Airplane accelerates wildly. 


## Current Work

I've made new debugging scripts to isolate the problem. They are two unsteady simulations of the same flat plate, one prescribed and one free flight. The flat plates have exactly one Panel, which should make it easier to see differences.

I finally got the prescribed and free flight simulations to match closely for the prescribed time steps. I also managed to eliminate the instability in the free flight time steps by setting any velocity components due to omegas to zero. This is non-physical, but hints at where in the code there may be remaining bugs.

## Active Issues

### Likely Culprits for Current Problem

None

### Less-Likely Culprits for Current Problem
1. When calling `ps.output.animate_free_flight`, the Airplane appears to be rotated 180 degrees about the body y-axis. However, I think this is an issue in the transformations in animate, as the wake appears to have convected in the correct direction (behind the Airplane), when I call the old draw method, which plots everything in the first Airplane's geometry axes.
2. When calling `ps.output.animate_free_flight`, the Airplane appears to be moving in the direction opposite what we expect. However, I think this is an issue in the transformations in animate, as the wake appears to have convected in the correct direction (behind the Airplane), when I call the old draw method, which plots everything in the first Airplane's geometry axes.
3. I need to correctly specify the initial orientation of the Airplane in the MuJoCo model. Right now, it's just aligned with the Earth axes. This is likely a small problem now, as we are using small angles.

### Future Issues
1. We have no way of visualizing the MuJoCo model itself, so it's hard to know if we've set up the XML correctly.
2. If the Airplane switches direction, then the wake will shed from the wrong edge. There isn't an easy way to fix this, but we check for this condition and terminate the simulation if it happens.
3. The MuJoCo model isn't intrinsically connected to the CoupledUnsteadyProblem. Therefore, if we change something in the model it doesn't automatically update the MuJoCo model. Perhaps, instead of using creating the XML separately, we could build the MuJoCo model programmatically from the CoupledUnsteadyProblem.'
4. When calling `ps.output.animate_free_flight`, the Airplane mesh disappears after several tens of time steps. It's almost like it's gone offscreen even though I've zoomed way out.

## Fixed Issues (Chronological Order)
1. I modified the simple glider Airplane, and the MuJoCo model, to:
    * Be based on an XFLR5 airplane that I know is statically stable in pitch and yaw
    * Use converged parameters
    * Operate at a trimmed condition
    * Use a known inertia matrix that is reasonable given its geometry (based on the XFLR project and verified via a simple Onshape model)
2. I fixed an issue where the free flight animation method was plotting the wake RingVortices using the wrong axes.
3. I incorporated the external wind axes x-axis force in the MuJoCo model.
4. I modified the coupled solver to iterate several time steps with prescribed motion. This is to prevent numerical issues relating to the large force spikes present when wake is first shed.
5. I modified the draw method to work with the coupled solver.
6. We were setting the initial orientation of the body axes in the MuJoCo model to be aligned with the Earth axes. This was incorrect because alpha wasn't 0.0. I fixed this by finding the correct initial rotation matrix, and then wrote a function to convert it into a quaternion.
7. We weren't correctly using the rotation matrix from MuJoCo. This resulted in incorrect angles of attack, which caused incorrect vInf_GP1__E vectors. I've fixed this for cases with alpha > 0 and beta = 0, and airplane.angles_E_to_B_izyx = (0, 0, 0), however, I think my solution is likely brittle.
8. I forced the solver to always have airplane.angles_E_to_B_izyx = (0, 0, 0) to see if there were any issues related to that parameter. I don't think it matters, but it is non-physical. We should revisit this once we've fixed the current problem, and come up with a better way of handling orientation offsets.
9. I modified PteraSoftware to remove the angles_E_to_B_izyx parameter from Airplane and instead have it be a parameter of CoupledOperatingPoint.
10. I added type-hints to new code and checked them with mypy.
11. I clarified that the coupled solver can only simulate a single airplane throughout the new code, and added validation checks for this condition.
12. I fixed the issue with Wing's span property (bug originating in main branch)
13. I fixed the issue with wake RingVortices being shed too far back from their Wing (bug originating in main branch)
