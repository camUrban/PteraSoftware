"""Contains the deprecated module path for the OperatingPointMovement class.

It is now available at the package top level, as pterasoftware.OperatingPointMovement.
Accessing it through this module emits a DeprecationWarning, and this module will be
removed in v6.0.0.
"""

from .._deprecation import make_deprecated_module

__getattr__, __dir__, __all__ = make_deprecated_module(
    __name__,
    "pterasoftware._movements.operating_point_movement",
    ("OperatingPointMovement",),
)
