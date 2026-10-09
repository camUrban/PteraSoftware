"""Contains the deprecated module path for the analyze_steady_convergence and
analyze_unsteady_convergence functions.

They are now available at the package top level, as
pterasoftware.analyze_steady_convergence and pterasoftware.analyze_unsteady_convergence.
Accessing one through this module emits a DeprecationWarning, and this module will be
removed in v6.0.0.
"""

from typing import Any

from ._deprecation import get_deprecated_attribute

NAMES = ("analyze_steady_convergence", "analyze_unsteady_convergence")


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(__name__, "pterasoftware._convergence", NAMES, name)


def __dir__() -> list[str]:
    return list(NAMES)
