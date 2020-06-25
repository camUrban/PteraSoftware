""" This package contains modules which create fixtures for use in integration tests.

This package contains the following subpackages:
    None

This package contains the following directories:
    None

This package contains the following modules:
    __init__.py: This module is this package's initialization script. It imports all the modules from this package.
    airplane_fixtures: This module creates airplane objects to be used as fixtures.
    movement_fixtures: This module creates movement objects to be used as fixtures.
    operating_point_fixtures: This module creates operating point objects to be used as fixtures.
    problem_fixtures: This module creates problem objects to be used as fixtures.
    solver_fixtures: This module creates solver objects to be used as fixtures.
"""

from tests.integration.fixtures import airplane_fixtures
from tests.integration.fixtures import movement_fixtures
from tests.integration.fixtures import operating_point_fixtures
from tests.integration.fixtures import problem_fixtures
from tests.integration.fixtures import solver_fixtures
