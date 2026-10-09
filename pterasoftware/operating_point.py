"""Contains the deprecated module path for the OperatingPoint class.

It is now available at the package top level, as pterasoftware.OperatingPoint. Accessing
it through this module emits a DeprecationWarning, and this module will be removed in
v6.0.0.
"""

from typing import Any

from ._deprecation import get_deprecated_attribute

NAMES = ("OperatingPoint",)


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(
        __name__, "pterasoftware._operating_point", NAMES, name
    )


def __dir__() -> list[str]:
    return list(NAMES)
