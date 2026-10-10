"""Contains the deprecated module path for the SteadyProblem, UnsteadyProblem,
FreeFlightUnsteadyProblem, and AeroelasticUnsteadyProblem classes.

They are now available at the package top level, as pterasoftware.SteadyProblem,
pterasoftware.UnsteadyProblem, pterasoftware.FreeFlightUnsteadyProblem, and
pterasoftware.AeroelasticUnsteadyProblem. Accessing one through this module emits a
DeprecationWarning, and this module will be removed in v6.0.0.
"""

from ._deprecation import make_deprecated_module

__getattr__, __dir__, __all__ = make_deprecated_module(
    __name__,
    "pterasoftware._problems",
    (
        "SteadyProblem",
        "UnsteadyProblem",
        "FreeFlightUnsteadyProblem",
        "AeroelasticUnsteadyProblem",
    ),
)
