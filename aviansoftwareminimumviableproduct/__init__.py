"""This package contains all the source code for the Avian Software Minimum Viable Product project.

This package contains the following subpackages:
    None

This package contains the following directories:
    airfoils: This folder contains a collection of airfoils whose coordinates are stored in DAT files.

This package contains the following modules:
    __init__.py: This module is this packages initialization script. It imports all the modules from this package.
    aerodynamics.py: This module contains useful aerodynamics functions.
    geometry.py: This module contains useful functions that relate to geometry, and the class definitions for different
                 types of geometries.
    meshing.py: This module contains useful functions for creating meshes.
    movement.py: This module contains the class definition for the geometry's movement.
    output.py: This module contains useful functions for visualizing solutions to problems.
    problems.py: This module contains the class definitions for different types of problems.
    steady_vortex_lattice_method.py: This module contains the class definition of this package's steady vortex lattice
                                     solver.
    unsteady_vortex_lattice_method.py: This module contains the class definition for this package's unsteady vortex
                                       lattice solver.
"""

from aviansoftwareminimumviableproduct import aerodynamics
from aviansoftwareminimumviableproduct import problems
from aviansoftwareminimumviableproduct import geometry
from aviansoftwareminimumviableproduct import meshing
from aviansoftwareminimumviableproduct import movement
from aviansoftwareminimumviableproduct import output
from aviansoftwareminimumviableproduct import steady_vortex_lattice_method
from aviansoftwareminimumviableproduct import unsteady_vortex_lattice_method
