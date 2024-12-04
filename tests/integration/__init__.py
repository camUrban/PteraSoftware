"""This package contains the integration tests.

This package contains the following subpackages:
    fixtures: This package contains the fixtures for integration testing.

This package contains the following directories:
    None

This package contains the following modules:
    __init__.py: This module is this package's initialization script.

    test_output.py: This module is a testing case for the output module.

    test_steady_convergence.py: This module contains testing cases for the steady
    convergence function.

    test_steady_horseshoe_vortex_lattice_method.py: This module is a testing case for
    the steady horseshoe vortex lattice method solver.

    test_steady_ring_vortex_lattice_method.py: This module is a testing case for the
    steady ring vortex lattice method solver.

    test_steady_trim.py: This module contains a testing case for the steady trim
    function.

    test_unsteady_convergence.py This module contains a testing case for the unsteady
    convergence function.

    test_unsteady_ring_vortex_lattice_method_multiple_wing_static_geometry.py: This
    is a testing case for the unsteady ring vortex lattice method solver with static,
    multi-wing geometry.

    test_unsteady_ring_vortex_lattice_method_multiple_wing_variable_geometry.py: This
    is a testing case for the unsteady ring vortex lattice method solver with
    variable, multi-wing geometry.

    test_unsteady_ring_vortex_lattice_method_static_geometry.py: This is a testing
    case for the unsteady ring vortex lattice method solver with static geometry.

    test_unsteady_ring_vortex_lattice_method_variable_geometry.py: This is a testing
    case for the unsteady ring vortex lattice method solver with variable geometry.
"""

import tests.integration.fixtures
import tests.integration.test_output
import tests.integration.test_steady_convergence
import tests.integration.test_steady_horseshoe_vortex_lattice_method
import tests.integration.test_steady_ring_vortex_lattice_method
import tests.integration.test_steady_trim
import tests.integration.test_unsteady_convergence
import tests.integration.test_unsteady_ring_vortex_lattice_method_multiple_wing_static_geometry
import tests.integration.test_unsteady_ring_vortex_lattice_method_multiple_wing_variable_geometry
import tests.integration.test_unsteady_ring_vortex_lattice_method_static_geometry
import tests.integration.test_unsteady_ring_vortex_lattice_method_variable_geometry
