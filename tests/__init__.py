"""This package contains the tests for Ptera Software.

This package contains the following subpackages:
    integration: This package contains integration tests.
    unit: This package contains the unit tests.

This package contains the following directories:
    benchmarks: This directory contains benchmark scripts and results.
    references: This directory contains reference files for tests.

This package contains the following modules:
    __init__.py: This module is this package's initialization script.

    _test_environment.py: This module quiets known sources of test run noise
    (serialization dirty tree warnings, tqdm progress bars, and the VTK
    headless display warning). It is imported first so its side effects take
    effect before any pterasoftware, pyvista, or tqdm module loads.
"""

# Must be the first import so environment variables are set before any other
# module reads them.
import tests._test_environment  # noqa: F401
import tests.integration
import tests.unit
