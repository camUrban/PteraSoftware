"""Contains the source code for Ptera Software."""

from typing import TYPE_CHECKING, Any

# Eager imports: core modules always needed to define simulations.
import pterasoftware._geometry
import pterasoftware._movements
import pterasoftware._operating_point
import pterasoftware._problems

# Static imports of every lazily loaded name, visible only to type checkers. Without
# these, any name resolved through __getattr__ is typed as Any, so every use of the lazy
# modules and callables through the package namespace goes unchecked. Each deprecated
# module name is typed as the internal module that now defines its public names.
if TYPE_CHECKING:
    from pterasoftware import (
        _aeroelastic_unsteady_ring_vortex_lattice_method as aeroelastic_unsteady_ring_vortex_lattice_method,
    )
    from pterasoftware import _convergence as convergence
    from pterasoftware import (
        _free_flight_unsteady_ring_vortex_lattice_method as free_flight_unsteady_ring_vortex_lattice_method,
    )
    from pterasoftware import _geometry as geometry
    from pterasoftware import _movements as movements
    from pterasoftware import _operating_point as operating_point
    from pterasoftware import _output as output
    from pterasoftware import _problems as problems
    from pterasoftware import (
        _steady_horseshoe_vortex_lattice_method as steady_horseshoe_vortex_lattice_method,
    )
    from pterasoftware import (
        _steady_ring_vortex_lattice_method as steady_ring_vortex_lattice_method,
    )
    from pterasoftware import _trim as trim
    from pterasoftware import (
        _unsteady_ring_vortex_lattice_method as unsteady_ring_vortex_lattice_method,
    )
    from pterasoftware._logging import set_up_logging
    from pterasoftware._serialization import load, save

# Lazy imports configuration: modules loaded on first access. Each is a deprecated
# module at an old public module path, which forwards that module's public names to
# their new locations.
_LAZY_MODULES = {
    "aeroelastic_unsteady_ring_vortex_lattice_method": "pterasoftware.aeroelastic_unsteady_ring_vortex_lattice_method",
    "convergence": "pterasoftware.convergence",
    "free_flight_unsteady_ring_vortex_lattice_method": "pterasoftware.free_flight_unsteady_ring_vortex_lattice_method",
    "geometry": "pterasoftware.geometry",
    "movements": "pterasoftware.movements",
    "operating_point": "pterasoftware.operating_point",
    "output": "pterasoftware.output",
    "problems": "pterasoftware.problems",
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
