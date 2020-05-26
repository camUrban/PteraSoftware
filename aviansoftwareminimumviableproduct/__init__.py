
"""This package contains all the source code for the Avian Software Minimum Viable Product project.

This package contains the following subpackages:
    None

This package contains the following directories:
    airfoils: This folder contains a collection of airfoils whose coordinates are stored in DAT files.

This package contains the following modules:
    __init__.py: This module is this package's initialization script. It imports all the modules from this package.
    aerodynamics.py: This module contains vortex class definitions.
    geometry.py: This module contains useful functions that relate to geometry, and the class definitions for different
                 types of geometries.
    meshing.py: This module contains useful functions for creating meshes.
    output.py: This module contains useful functions for visualizing solutions to problems.
    performance.py: This module contains the class definitions for the problem's movement and the problem's operating
                    point.
    problems.py: This module contains the class definitions for different types of problems.
    steady_horseshoe_vortex_lattice_method.py: This module contains the class definition of this package's steady
                                               horseshoe vortex lattice solver.
    steady_ring_vortex_lattice_method.py: This module contains the class definition of this package's steady ring vortex
                                          lattice solver.
    unsteady_ring_vortex_lattice_method.py: This module contains the class definition of this package's unsteady ring
                                            vortex lattice solver.

"""

from aviansoftwareminimumviableproduct import aerodynamics
from aviansoftwareminimumviableproduct import geometry
from aviansoftwareminimumviableproduct import meshing
from aviansoftwareminimumviableproduct import output
from aviansoftwareminimumviableproduct import performance
from aviansoftwareminimumviableproduct import problems
from aviansoftwareminimumviableproduct import steady_horseshoe_vortex_lattice_method
from aviansoftwareminimumviableproduct import steady_ring_vortex_lattice_method
from aviansoftwareminimumviableproduct import unsteady_ring_vortex_lattice_method
