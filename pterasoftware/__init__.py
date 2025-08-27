"""This package contains all the source code for the Ptera Software.

This package contains the following subpackages:
    None

This package contains the following directories:
    airfoils: This folder contains a collection of airfoils whose coordinates are
    stored in DAT files.

This package contains the following modules:
    __init__.py: This module is this package's initialization script.

    aerodynamics.py: This module contains vortex class definitions.

    convergence.py: This module contains functions for analyzing the convergence of
    steady and unsteady problems.

    functions.py: This module contains functions shared by other modules in the
    pterasoftware package.

    geometry.py: This module contains useful class definitions for different types of
    geometries.

    meshing.py: This module contains useful functions for creating meshes.

    movement.py: This module contains the class definitions for the problem's movement.

    operating_point.py: This module contains the class definition for the problem's
    operating point.

    output.py: This module contains useful functions for visualizing solutions to
    problems.

    panel.py: This module contains the Panel class.

    parameter_validation.py: This module contains common parameter validation functions.

    problems.py: This module contains the class definitions for different types of
    problems.

    steady_horseshoe_vortex_lattice_method.py: This module contains the class
    definition of this package's steady horseshoe vortex lattice solver.

    steady_ring_vortex_lattice_method.py: This module contains the class definition
    of this package's steady ring vortex lattice solver.

    transformations.py This module contains functions use for geometric transformations.

    trim.py: This module contains functions to analyze the trim conditions of steady
    and unsteady solvers.

    unsteady_ring_vortex_lattice_method.py: This module contains the class definition
    of this package's unsteady ring vortex lattice solver.
"""

import pterasoftware.aerodynamics
import pterasoftware.airfoils
import pterasoftware.convergence
import pterasoftware.functions
import pterasoftware.geometry
import pterasoftware.meshing
import pterasoftware.movement
import pterasoftware.operating_point
import pterasoftware.output
import pterasoftware.panel
import pterasoftware.parameter_validation
import pterasoftware.problems
import pterasoftware.steady_horseshoe_vortex_lattice_method
import pterasoftware.steady_ring_vortex_lattice_method
import pterasoftware.transformations
import pterasoftware.trim
import pterasoftware.unsteady_ring_vortex_lattice_method
