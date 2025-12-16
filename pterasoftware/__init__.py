"""Contains the source code for the Ptera Software.

**Contains the following subpackages:**

geometry: Contains the geometry classes.

movements: Contains the movement classes.

**Contains the following directories:**

None

**Contains the following modules:**

convergence.py: Contains functions for analyzing the convergence of SteadyProblems and
UnsteadyProblems.

coupled_unsteady_ring_vortex_lattice_method.py Contains the
CoupledUnsteadyRingVortexLatticeMethodSolver class.

operating_point.py: Contains classes for the operating conditions of aerodynamic
problems.

output.py: Contains functions for visualizing geometry and results.

problems.py: Contains classes for aerodynamic problems.

steady_horseshoe_vortex_lattice_method.py: Contains the
SteadyHorseshoeVortexLatticeMethodSolver class.

steady_ring_vortex_lattice_method.py: Contains the SteadyRingVortexLatticeMethodSolver
class.

trim.py: Contains functions to analyze the trim conditions of SteadyProblems and
UnsteadyProblems.

unsteady_ring_vortex_lattice_method.py: Contains the
UnsteadyRingVortexLatticeMethodSolver class.

**Contains the following functions:**

set_up_logging: Configures logging for the pterasoftware package that is compatible with
TQDM progress bars.
"""

import pterasoftware.convergence
import pterasoftware.coupled_unsteady_ring_vortex_lattice_method
import pterasoftware.geometry
import pterasoftware.movements
import pterasoftware.operating_point
import pterasoftware.output
import pterasoftware.problems
import pterasoftware.steady_horseshoe_vortex_lattice_method
import pterasoftware.steady_ring_vortex_lattice_method
import pterasoftware.trim
import pterasoftware.unsteady_ring_vortex_lattice_method
from pterasoftware._logging import set_up_logging
