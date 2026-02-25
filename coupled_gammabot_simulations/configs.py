"""Configuration parameters for coupled GammaBot simulations.

Re-exports wing motion configurations from gammabot_simulations/configs.py and adds
coupled simulation specific parameters (mass, inertia, drag, etc.).
"""

import importlib.util
from pathlib import Path

# Load gammabot_simulations/configs.py by file path to avoid circular import (both
# modules are named "configs").
_gammabot_configs_path = (
    Path(__file__).parent.parent / "gammabot_simulations" / "configs.py"
)
_spec = importlib.util.spec_from_file_location(
    "gammabot_configs", _gammabot_configs_path
)
assert _spec is not None and _spec.loader is not None
_gammabot_configs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_gammabot_configs)

CONFIGURATIONS = _gammabot_configs.CONFIGURATIONS
MESH_PARAMS = _gammabot_configs.MESH_PARAMS
SHARED_PARAMS = _gammabot_configs.SHARED_PARAMS
get_config = _gammabot_configs.get_config
list_configs = _gammabot_configs.list_configs

# Coupled simulation specific parameters.
COUPLED_PARAMS = {
    # TODO: Replace placeholder mass/inertia with values from CAD.
    "weight": 0.01 * 9.80665,  # Placeholder weight (N).
    "I_BP1_CgP1": [  # Placeholder inertia tensor (kg*m^2).
        [1.0e-8, 0.0, 0.0],
        [0.0, 1.0e-8, 0.0],
        [0.0, 0.0, 1.0e-8],
    ],
    "body_drag_c1": 2.5526,  # From water_drag_constant_estimation.py (1/m).
    "rho": 1.225,  # Air density (kg/m^3). Wings are in air.
    "nu": 15.06e-6,  # Air kinematic viscosity (m^2/s).
    "g_E": (0.0, 0.0, 9.80665),  # Gravity (+z down).
}

__all__ = [
    "CONFIGURATIONS",
    "COUPLED_PARAMS",
    "MESH_PARAMS",
    "SHARED_PARAMS",
    "get_config",
    "list_configs",
]
