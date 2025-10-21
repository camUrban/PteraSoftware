"""This package contains all the source code for the Ptera Software.

This package contains the following subpackages:
    geometry: This package contains the geometry classes.

    movements: This package contains the movement classes.

This package contains the following directories:
    None

This package contains the following modules:
    __init__.py: This module is this package's initialization script.

    convergence.py: This module contains functions for analyzing the convergence of
    SteadyProblems and UnsteadyProblems.

    operating_point.py: This module contains the class definition for a Problem's
    operating point.

    output.py: This module contains useful functions for visualizing geometry and
    results.

    problems.py: This module contains the class definitions for different types of
    problems.

    steady_horseshoe_vortex_lattice_method.py: This module contains the class
    definition of this package's steady horseshoe vortex lattice solver.

    steady_ring_vortex_lattice_method.py: This module contains the class definition
    of this package's steady ring vortex lattice solver.

    trim.py: This module contains functions to analyze the trim conditions of
    SteadyProblems and UnsteadyProblems.

    unsteady_ring_vortex_lattice_method.py: This module contains the class definition
    of this package's unsteady ring vortex lattice solver.
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
