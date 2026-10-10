"""Contains the deprecated module path for the UnsteadyRingVortexLatticeMethodSolver
class.

It is now available at the package top level, as
pterasoftware.UnsteadyRingVortexLatticeMethodSolver. Accessing it through this module
emits a DeprecationWarning, and this module will be removed in v6.0.0.
"""

from ._deprecation import make_deprecated_module

__getattr__, __dir__, __all__ = make_deprecated_module(
    __name__,
    "pterasoftware._unsteady_ring_vortex_lattice_method",
    ("UnsteadyRingVortexLatticeMethodSolver",),
)
