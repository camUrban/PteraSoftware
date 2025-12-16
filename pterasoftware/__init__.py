"""Contains the source code for Ptera Software.

**Contains the following subpackages:**

geometry: Contains the geometry classes.

movements: Contains the movement classes.

**Contains the following directories:**

None

**Contains the following modules:**

convergence.py: Contains functions for analyzing the convergence of SteadyProblems and
UnsteadyProblems.

operating_point.py: Contains the OperatingPoint class.

output.py: Contains functions for visualizing geometry and results.

problems.py: Contains the SteadyProblem and UnsteadyProblem classes.

steady_horseshoe_vortex_lattice_method.py: Contains the
SteadyHorseshoeVortexLatticeMethodSolver class.

steady_ring_vortex_lattice_method.py: Contains the SteadyRingVortexLatticeMethodSolver
class.

trim.py: Contains functions to analyze the trim conditions of SteadyProblems and
UnsteadyProblems.

unsteady_ring_vortex_lattice_method.py: Contains the
UnsteadyRingVortexLatticeMethodSolver class.

**Contains the following functions:**

set_up_logging: Configures logging for the pterasoftware package that is compatible with
TQDM progress bars.
"""

# Eager imports: core modules always needed to define simulations.
import pterasoftware.geometry
import pterasoftware.movements
import pterasoftware.operating_point
import pterasoftware.problems

# Lazy imports configuration: modules loaded on first access.
_LAZY_MODULES = {
    "convergence": "pterasoftware.convergence",
    "output": "pterasoftware.output",
    "steady_horseshoe_vortex_lattice_method": "pterasoftware.steady_horseshoe_vortex_lattice_method",
    "steady_ring_vortex_lattice_method": "pterasoftware.steady_ring_vortex_lattice_method",
    "trim": "pterasoftware.trim",
    "unsteady_ring_vortex_lattice_method": "pterasoftware.unsteady_ring_vortex_lattice_method",
}

# Lazy callable imports: functions that need special handling.
_LAZY_CALLABLES = {
    "set_up_logging": ("pterasoftware._logging", "set_up_logging"),
}


def __getattr__(name: str):
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
    raise AttributeError(f"module 'pterasoftware' has no attribute {name!r}")


def __dir__():
    # Include lazy modules in dir() for discoverability.
    return (
        list(globals().keys())
        + list(_LAZY_MODULES.keys())
        + list(_LAZY_CALLABLES.keys())
    )
