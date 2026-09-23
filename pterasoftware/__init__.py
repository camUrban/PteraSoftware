"""Contains the source code for Ptera Software."""

from typing import TYPE_CHECKING, Any

# Eager imports: core modules always needed to define simulations.
import pterasoftware.geometry
import pterasoftware.movements
import pterasoftware.operating_point
import pterasoftware.problems

# Static imports of every lazily loaded name, visible only to type checkers. Without
# these, any name resolved through __getattr__ is typed as Any, so every use of the lazy
# modules and callables through the package namespace goes unchecked.
if TYPE_CHECKING:
    from pterasoftware import (
        aeroelastic_unsteady_ring_vortex_lattice_method,
        convergence,
        free_flight_unsteady_ring_vortex_lattice_method,
        output,
        steady_horseshoe_vortex_lattice_method,
        steady_ring_vortex_lattice_method,
        trim,
        unsteady_ring_vortex_lattice_method,
    )
    from pterasoftware._logging import set_up_logging
    from pterasoftware._serialization import load, save

# Lazy imports configuration: modules loaded on first access.
_LAZY_MODULES = {
    "aeroelastic_unsteady_ring_vortex_lattice_method": "pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method",
    "convergence": "pterasoftware.convergence",
    "free_flight_unsteady_ring_vortex_lattice_method": "pterasoftware.free_flight_unsteady_ring_vortex_lattice_method",
    "output": "pterasoftware.output",
    "steady_horseshoe_vortex_lattice_method": "pterasoftware.steady_horseshoe_vortex_lattice_method",
    "steady_ring_vortex_lattice_method": "pterasoftware.steady_ring_vortex_lattice_method",
    "trim": "pterasoftware.trim",
    "unsteady_ring_vortex_lattice_method": "pterasoftware.unsteady_ring_vortex_lattice_method",
}

# Lazy callable imports: functions that need special handling.
_LAZY_CALLABLES = {
    "load": ("pterasoftware._serialization", "load"),
    "save": ("pterasoftware._serialization", "save"),
    "set_up_logging": ("pterasoftware._logging", "set_up_logging"),
}


def __getattr__(name: str) -> Any:
    if name in _LAZY_MODULES:
        import importlib

        module = importlib.import_module(_LAZY_MODULES[name])
        globals()[name] = module
        return module
    if name in _LAZY_CALLABLES:
        import importlib

        module_path, attr_name = _LAZY_CALLABLES[name]
        module = importlib.import_module(module_path)
        attr = getattr(module, attr_name)
        globals()[name] = attr
        return attr
    raise AttributeError(f'module "pterasoftware" has no attribute "{name}"')


def __dir__() -> list[str]:
    # Include lazy modules in dir() for discoverability.
    return (
        list(globals().keys())
        + list(_LAZY_MODULES.keys())
        + list(_LAZY_CALLABLES.keys())
    )
