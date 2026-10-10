"""Contains the source code for Ptera Software."""

import importlib
from typing import TYPE_CHECKING, Any

# The public API. Every public class and function is available here, at the package top
# level, and nowhere else.
__all__ = [
    "AeroelasticAirplaneMovement",
    "AeroelasticMovement",
    "AeroelasticUnsteadyProblem",
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver",
    "AeroelasticWingCrossSectionMovement",
    "AeroelasticWingMovement",
    "Airfoil",
    "Airplane",
    "AirplaneMovement",
    "FreeFlightMovement",
    "FreeFlightOperatingPointMovement",
    "FreeFlightUnsteadyProblem",
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver",
    "Movement",
    "OperatingPoint",
    "OperatingPointMovement",
    "SteadyHorseshoeVortexLatticeMethodSolver",
    "SteadyProblem",
    "SteadyRingVortexLatticeMethodSolver",
    "UnsteadyProblem",
    "UnsteadyRingVortexLatticeMethodSolver",
    "Wing",
    "WingCrossSection",
    "WingCrossSectionMovement",
    "WingMovement",
    "analyze_steady_convergence",
    "analyze_steady_trim",
    "analyze_unsteady_convergence",
    "analyze_unsteady_trim",
    "animate",
    "draw",
    "load",
    "log_results",
    "plot_results_versus_time",
    "save",
    "set_up_logging",
]

# Static imports of every public name, visible only to type checkers. Without these, any
# name resolved through __getattr__ is typed as Any, so every use of the public names
# through the package namespace goes unchecked.
if TYPE_CHECKING:
    from pterasoftware._aeroelastic_unsteady_ring_vortex_lattice_method import (
        AeroelasticUnsteadyRingVortexLatticeMethodSolver,
    )
    from pterasoftware._convergence import (
        analyze_steady_convergence,
        analyze_unsteady_convergence,
    )
    from pterasoftware._free_flight_unsteady_ring_vortex_lattice_method import (
        FreeFlightUnsteadyRingVortexLatticeMethodSolver,
    )
    from pterasoftware._geometry.airfoil import Airfoil
    from pterasoftware._geometry.airplane import Airplane
    from pterasoftware._geometry.wing import Wing
    from pterasoftware._geometry.wing_cross_section import WingCrossSection
    from pterasoftware._logging import set_up_logging
    from pterasoftware._movements.aeroelastic_airplane_movement import (
        AeroelasticAirplaneMovement,
    )
    from pterasoftware._movements.aeroelastic_movement import AeroelasticMovement
    from pterasoftware._movements.aeroelastic_wing_cross_section_movement import (
        AeroelasticWingCrossSectionMovement,
    )
    from pterasoftware._movements.aeroelastic_wing_movement import (
        AeroelasticWingMovement,
    )
    from pterasoftware._movements.airplane_movement import AirplaneMovement
    from pterasoftware._movements.free_flight_movement import FreeFlightMovement
    from pterasoftware._movements.free_flight_operating_point_movement import (
        FreeFlightOperatingPointMovement,
    )
    from pterasoftware._movements.movement import Movement
    from pterasoftware._movements.operating_point_movement import (
        OperatingPointMovement,
    )
    from pterasoftware._movements.wing_cross_section_movement import (
        WingCrossSectionMovement,
    )
    from pterasoftware._movements.wing_movement import WingMovement
    from pterasoftware._operating_point import OperatingPoint
    from pterasoftware._output import (
        animate,
        draw,
        log_results,
        plot_results_versus_time,
    )
    from pterasoftware._problems import (
        AeroelasticUnsteadyProblem,
        FreeFlightUnsteadyProblem,
        SteadyProblem,
        UnsteadyProblem,
    )
    from pterasoftware._serialization import load, save
    from pterasoftware._steady_horseshoe_vortex_lattice_method import (
        SteadyHorseshoeVortexLatticeMethodSolver,
    )
    from pterasoftware._steady_ring_vortex_lattice_method import (
        SteadyRingVortexLatticeMethodSolver,
    )
    from pterasoftware._trim import analyze_steady_trim, analyze_unsteady_trim
    from pterasoftware._unsteady_ring_vortex_lattice_method import (
        UnsteadyRingVortexLatticeMethodSolver,
    )

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

# Lazy callable imports: every name in __all__, mapped to the internal module that
# defines it and its name there. Each is loaded on first access.
_LAZY_CALLABLES = {
    "AeroelasticAirplaneMovement": (
        "pterasoftware._movements.aeroelastic_airplane_movement",
        "AeroelasticAirplaneMovement",
    ),
    "AeroelasticMovement": (
        "pterasoftware._movements.aeroelastic_movement",
        "AeroelasticMovement",
    ),
    "AeroelasticUnsteadyProblem": (
        "pterasoftware._problems",
        "AeroelasticUnsteadyProblem",
    ),
    "AeroelasticUnsteadyRingVortexLatticeMethodSolver": (
        "pterasoftware._aeroelastic_unsteady_ring_vortex_lattice_method",
        "AeroelasticUnsteadyRingVortexLatticeMethodSolver",
    ),
    "AeroelasticWingCrossSectionMovement": (
        "pterasoftware._movements.aeroelastic_wing_cross_section_movement",
        "AeroelasticWingCrossSectionMovement",
    ),
    "AeroelasticWingMovement": (
        "pterasoftware._movements.aeroelastic_wing_movement",
        "AeroelasticWingMovement",
    ),
    "Airfoil": ("pterasoftware._geometry.airfoil", "Airfoil"),
    "Airplane": ("pterasoftware._geometry.airplane", "Airplane"),
    "AirplaneMovement": (
        "pterasoftware._movements.airplane_movement",
        "AirplaneMovement",
    ),
    "FreeFlightMovement": (
        "pterasoftware._movements.free_flight_movement",
        "FreeFlightMovement",
    ),
    "FreeFlightOperatingPointMovement": (
        "pterasoftware._movements.free_flight_operating_point_movement",
        "FreeFlightOperatingPointMovement",
    ),
    "FreeFlightUnsteadyProblem": (
        "pterasoftware._problems",
        "FreeFlightUnsteadyProblem",
    ),
    "FreeFlightUnsteadyRingVortexLatticeMethodSolver": (
        "pterasoftware._free_flight_unsteady_ring_vortex_lattice_method",
        "FreeFlightUnsteadyRingVortexLatticeMethodSolver",
    ),
    "Movement": ("pterasoftware._movements.movement", "Movement"),
    "OperatingPoint": ("pterasoftware._operating_point", "OperatingPoint"),
    "OperatingPointMovement": (
        "pterasoftware._movements.operating_point_movement",
        "OperatingPointMovement",
    ),
    "SteadyHorseshoeVortexLatticeMethodSolver": (
        "pterasoftware._steady_horseshoe_vortex_lattice_method",
        "SteadyHorseshoeVortexLatticeMethodSolver",
    ),
    "SteadyProblem": ("pterasoftware._problems", "SteadyProblem"),
    "SteadyRingVortexLatticeMethodSolver": (
        "pterasoftware._steady_ring_vortex_lattice_method",
        "SteadyRingVortexLatticeMethodSolver",
    ),
    "UnsteadyProblem": ("pterasoftware._problems", "UnsteadyProblem"),
    "UnsteadyRingVortexLatticeMethodSolver": (
        "pterasoftware._unsteady_ring_vortex_lattice_method",
        "UnsteadyRingVortexLatticeMethodSolver",
    ),
    "Wing": ("pterasoftware._geometry.wing", "Wing"),
    "WingCrossSection": (
        "pterasoftware._geometry.wing_cross_section",
        "WingCrossSection",
    ),
    "WingCrossSectionMovement": (
        "pterasoftware._movements.wing_cross_section_movement",
        "WingCrossSectionMovement",
    ),
    "WingMovement": ("pterasoftware._movements.wing_movement", "WingMovement"),
    "analyze_steady_convergence": (
        "pterasoftware._convergence",
        "analyze_steady_convergence",
    ),
    "analyze_steady_trim": ("pterasoftware._trim", "analyze_steady_trim"),
    "analyze_unsteady_convergence": (
        "pterasoftware._convergence",
        "analyze_unsteady_convergence",
    ),
    "analyze_unsteady_trim": ("pterasoftware._trim", "analyze_unsteady_trim"),
    "animate": ("pterasoftware._output", "animate"),
    "draw": ("pterasoftware._output", "draw"),
    "load": ("pterasoftware._serialization", "load"),
    "log_results": ("pterasoftware._output", "log_results"),
    "plot_results_versus_time": (
        "pterasoftware._output",
        "plot_results_versus_time",
    ),
    "save": ("pterasoftware._serialization", "save"),
    "set_up_logging": ("pterasoftware._logging", "set_up_logging"),
}


def __getattr__(name: str) -> Any:
    if name in _LAZY_CALLABLES:
        module_path, attr_name = _LAZY_CALLABLES[name]
        module = importlib.import_module(module_path)
        attr = getattr(module, attr_name)
        globals()[name] = attr
        return attr
    if name in _LAZY_MODULES:
        module = importlib.import_module(_LAZY_MODULES[name])
        globals()[name] = module
        return module
    raise AttributeError(f'module "pterasoftware" has no attribute "{name}"')


def __dir__() -> list[str]:
    return list(__all__)
