"""Contains the deprecated module path for the analyze_steady_trim and
analyze_unsteady_trim functions.

They are now available at the package top level, as pterasoftware.analyze_steady_trim
and pterasoftware.analyze_unsteady_trim. Accessing one through this module emits a
DeprecationWarning, and this module will be removed in v6.0.0.
"""

from typing import Any

from ._deprecation import get_deprecated_attribute

NAMES = ("analyze_steady_trim", "analyze_unsteady_trim")


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(__name__, "pterasoftware._trim", NAMES, name)


def __dir__() -> list[str]:
    return list(NAMES)
