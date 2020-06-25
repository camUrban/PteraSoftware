""" This package contains integration tests.

This package contains the following subpackages:
    None

This package contains the following directories:
    fixtures: This package contains modules which create fixtures for use in integration tests.

This package contains the following modules:
    __init__: This module is this package's initialization script.
    test_output: This module is a testing case for the output module.
    test_steady_horseshoe_vortex_lattice_method: This module is a testing case for the steady horseshoe vortex lattice
                                                 method solver.
    test_steady_ring_vortex_lattice_method: This module is a testing case for the steady ring vortex lattice method
                                            solver.
    test_unsteady_ring_vortex_lattice_method_static_geometry: This module is a testing case for the unsteady ring vortex
                                                              lattice method solver with static geometry.
    test_unsteady_ring_vortex_lattice_method_variable_geometry: This module is a testing case for the unsteady ring
                                                                vortex lattice method solver with variable geometry.
"""

from tests.integration import fixtures
