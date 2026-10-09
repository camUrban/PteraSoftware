"""Contains the deprecated module path for the Airfoil class.

It is now available at the package top level, as pterasoftware.Airfoil. Accessing it
through this module emits a DeprecationWarning, and this module will be removed in
v6.0.0.
"""

from typing import Any

from .._deprecation import get_deprecated_attribute

NAMES = ("Airfoil",)


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(
        __name__, "pterasoftware._geometry.airfoil", NAMES, name
    )


def __dir__() -> list[str]:
    return list(NAMES)
