"""Contains the deprecated module path for the Movement class.

It is now available at the package top level, as pterasoftware.Movement. Accessing it
through this module emits a DeprecationWarning, and this module will be removed in
v6.0.0.
"""

from typing import Any

from .._deprecation import get_deprecated_attribute

NAMES = ("Movement",)


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(
        __name__, "pterasoftware._movements.movement", NAMES, name
    )


def __dir__() -> list[str]:
    return list(NAMES)
