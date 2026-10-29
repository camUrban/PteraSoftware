"""Contains deprecated module paths for camera state and output functions.

They are now available at the package top level, as pterasoftware.CameraState,
pterasoftware.draw, pterasoftware.animate, pterasoftware.plot_results_versus_time, and
pterasoftware.log_results. Accessing a name through this module emits a
DeprecationWarning, and this module will be removed in v6.0.0.
"""

from ._deprecation import make_deprecated_module

__getattr__, __dir__, __all__ = make_deprecated_module(
    __name__,
    "pterasoftware._output",
    (
        "CameraState",
        "draw",
        "animate",
        "plot_results_versus_time",
        "log_results",
    ),
)
