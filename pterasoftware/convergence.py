"""Contains the deprecated module path for the analyze_steady_convergence and
analyze_unsteady_convergence functions.

They are now available at the package top level, as
pterasoftware.analyze_steady_convergence and pterasoftware.analyze_unsteady_convergence.
Accessing one through this module emits a DeprecationWarning, and this module will be
removed in v6.0.0.
"""

from ._deprecation import make_deprecated_module

__getattr__, __dir__, __all__ = make_deprecated_module(
    __name__,
    "pterasoftware._convergence",
    ("analyze_steady_convergence", "analyze_unsteady_convergence"),
)
