"""Contains the deprecated module path for the draw, animate, plot_results_versus_time,
and log_results functions.

They are now available at the package top level, as pterasoftware.draw,
pterasoftware.animate, pterasoftware.plot_results_versus_time, and
pterasoftware.log_results. Accessing one through this module emits a DeprecationWarning,
and this module will be removed in v6.0.0.
"""

from typing import Any

from ._deprecation import get_deprecated_attribute

NAMES = ("draw", "animate", "plot_results_versus_time", "log_results")


def __getattr__(name: str) -> Any:
    return get_deprecated_attribute(__name__, "pterasoftware._output", NAMES, name)


def __dir__() -> list[str]:
    return list(NAMES)
