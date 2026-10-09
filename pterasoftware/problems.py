"""Contains the deprecated module path for the SteadyProblem, UnsteadyProblem,
FreeFlightUnsteadyProblem, and AeroelasticUnsteadyProblem classes.

They are now available at the package top level, as pterasoftware.SteadyProblem,
pterasoftware.UnsteadyProblem, pterasoftware.FreeFlightUnsteadyProblem, and
pterasoftware.AeroelasticUnsteadyProblem. Accessing one through this module emits a
DeprecationWarning, and this module will be removed in v6.0.0.
"""

from typing import Any

from ._deprecation import get_deprecated_attribute

NAMES = (
    "SteadyProblem",
    "UnsteadyProblem",
    "FreeFlightUnsteadyProblem",
    "AeroelasticUnsteadyProblem",
)


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(__name__, "pterasoftware._problems", NAMES, name)


def __dir__() -> list[str]:
    return list(NAMES)
