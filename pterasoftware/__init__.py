""" This package contains all the source code for the Ptera Software.

This package contains the following subpackages:
    None

This package contains the following directories:
    airfoils: This folder contains a collection of airfoils whose coordinates are
    stored in DAT files.

This package contains the following modules:
    __init__.py: This module is this package's initialization script.
    aerodynamics.py: This module contains vortex class definitions.
    geometry.py: This module contains useful functions that relate to geometry,
    and the class definitions for different
                 types of geometries.
    meshing.py: This module contains useful functions for creating meshes.
    output.py: This module contains useful functions for visualizing solutions to
    problems.
    movement.py: This module contains the class definitions for the problem's movement.
    current_operating_point.py: This module contains the class definition for the
    problem's operating point.
    problems.py: This module contains the class definitions for different types of
    problems.
    steady_horseshoe_vortex_lattice_method.py: This module contains the class
    definition of this package's steady
                                               horseshoe vortex lattice solver.
    steady_ring_vortex_lattice_method.py: This module contains the class definition
    of this package's steady ring vortex
                                          lattice solver.
    unsteady_ring_vortex_lattice_method.py: This module contains the class definition
    of this package's unsteady ring
                                            vortex lattice solver.

"""

from pterasoftware import aerodynamics
from pterasoftware import airfoils
from pterasoftware import geometry
from pterasoftware import meshing
from pterasoftware import movement
from pterasoftware import operating_point
from pterasoftware import output
from pterasoftware import problems
from pterasoftware import steady_horseshoe_vortex_lattice_method
from pterasoftware import steady_ring_vortex_lattice_method
from pterasoftware import unsteady_ring_vortex_lattice_method
