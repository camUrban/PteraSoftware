"""Contains the deprecated module path for the draw, animate, plot_results_versus_time,
and log_results functions.

They are now available at the package top level, as pterasoftware.draw,
pterasoftware.animate, pterasoftware.plot_results_versus_time, and
pterasoftware.log_results. Accessing one through this module emits a DeprecationWarning,
and this module will be removed in v6.0.0.
"""

from ._deprecation import make_deprecated_module

__getattr__, __dir__, __all__ = make_deprecated_module(
    __name__,
    "pterasoftware._output",
    ("draw", "animate", "plot_results_versus_time", "log_results"),
)
