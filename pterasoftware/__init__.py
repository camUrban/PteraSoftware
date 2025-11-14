"""Contains the source code for the Ptera Software.

**Contains the following subpackages:**

geometry: Contains the geometry classes.

movements: Contains the movement classes.

**Contains the following directories:**

None

**Contains the following modules:**

convergence.py: Contains functions for analyzing the convergence of SteadyProblems and
UnsteadyProblems.

operating_point.py: Contains the OperatingPoint class.

output.py: Contains useful functions for visualizing geometry and results.

problems.py: Contains the SteadyProblem and UnsteadyProblem classes.

steady_horseshoe_vortex_lattice_method.py: Contains the class definition of this
package's steady horseshoe vortex lattice solver.

steady_ring_vortex_lattice_method.py: Contains the class definition of this package's
steady ring vortex lattice solver.

trim.py: Contains functions to analyze the trim conditions of SteadyProblems and
UnsteadyProblems.

unsteady_ring_vortex_lattice_method.py: Contains the class definition of this package's
unsteady ring vortex lattice solver.
"""

import pterasoftware.geometry
import pterasoftware.movements
import pterasoftware.convergence
import pterasoftware.operating_point
import pterasoftware.output
import pterasoftware.problems
import pterasoftware.steady_horseshoe_vortex_lattice_method
import pterasoftware.steady_ring_vortex_lattice_method
import pterasoftware.trim
import pterasoftware.unsteady_ring_vortex_lattice_method
