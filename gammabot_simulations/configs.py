"""Configuration parameters for GammaBot simulations.

Each configuration represents a different voltage setting for left and right wings
at 170Hz flapping frequency. The naming convention is L{voltage}V_R{voltage}V_170Hz.
"""

# Shared parameters across all configurations
SHARED_PARAMS = {
    "velocity": 0.9,
    "alpha": 0.0,
    "flapping_frequency": 170.0,
    "wing_spacing": 0.02746,
    "x_offset": 1.215e-6,
    "y_offset": -0.00215,
    "chordwise_spacing": "uniform",
}

# Mesh settings
MESH_PARAMS = {
    "coarse": {
        "num_flaps": 2,
        "num_chordwise_panels": 12,
        "num_spanwise_sections": 13,
        "delta_time": 7.261728395061729e-05,
        "prescribed_wake": False,
    },
    "fine": {
        "num_flaps": 3,
        "num_chordwise_panels": 13,
        "num_spanwise_sections": 14,
        "delta_time": 6.68409090909091e-05,
        "prescribed_wake": False,
    },
}

# Wing motion configurations
# Keys: phi_max, phi_v_shift, psi_max, psi_v_shift, delta
CONFIGURATIONS = {
    "L0V_R150V_170Hz": {
        "left": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
        "right": {
            "phi_max": 15.63,
            "phi_v_shift": 3.58,
            "psi_max": 37.07,
            "psi_v_shift": 2.72,
            "delta": 15.40,
        },
    },
    "L0V_R160V_170Hz": {
        "left": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
        "right": {
            "phi_max": 17.49,
            "phi_v_shift": 3.09,
            "psi_max": 40.44,
            "psi_v_shift": 1.71,
            "delta": 17.56,
        },
    },
    "L0V_R170V_170Hz": {
        "left": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
        "right": {
            "phi_max": 18.90,
            "phi_v_shift": 3.06,
            "psi_max": 42.79,
            "psi_v_shift": 1.1,
            "delta": 16.25,
        },
    },
    "L0V_R180V_170Hz": {
        "left": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
        "right": {
            "phi_max": 20.22,
            "phi_v_shift": 3.17,
            "psi_max": 44.02,
            "psi_v_shift": 1.91,
            "delta": 15.65,
        },
    },
    "L150V_R0V_170Hz": {
        "left": {
            "phi_max": 15.58,
            "phi_v_shift": 1.12,
            "psi_max": 30.60,
            "psi_v_shift": 0.30,
            "delta": 36.06,
        },
        "right": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
    },
    "L160V_R0V_170Hz": {
        "left": {
            "phi_max": 16.53,
            "phi_v_shift": 0.98,
            "psi_max": 32.37,
            "psi_v_shift": 0.60,
            "delta": 34.64,
        },
        "right": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
    },
    "L170V_R0V_170Hz": {
        "left": {
            "phi_max": 17.54,
            "phi_v_shift": 1.04,
            "psi_max": 34.18,
            "psi_v_shift": 0.70,
            "delta": 32.00,
        },
        "right": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
    },
    "L170V_R180V_170Hz": {
        "left": {
            "phi_max": 17.81,
            "phi_v_shift": 1.26,
            "psi_max": 33.83,
            "psi_v_shift": 0.51,
            "delta": 35.22,
        },
        "right": {
            "phi_max": 20.20,
            "phi_v_shift": 3.06,
            "psi_max": 43.68,
            "psi_v_shift": 1.78,
            "delta": 16.91,
        },
    },
    "L180V_R0V_170Hz": {
        "left": {
            "phi_max": 18.59,
            "phi_v_shift": 1.23,
            "psi_max": 35.94,
            "psi_v_shift": 0.86,
            "delta": 29.28,
        },
        "right": {
            "phi_max": 0.0,
            "phi_v_shift": 0.0,
            "psi_max": 0.0,
            "psi_v_shift": 0.0,
            "delta": 0.0,
        },
    },
}


def get_config(name: str) -> dict:
    """Get a configuration by name.

    :param name: Configuration name (e.g., "L0V_R150V_170Hz").

    :return: Dictionary with left and right wing parameters.

    :raises KeyError: If configuration name is not found.
    """
    if name not in CONFIGURATIONS:
        available = ", ".join(sorted(CONFIGURATIONS.keys()))
        raise KeyError(f"Unknown configuration: {name}. Available: {available}")
    return CONFIGURATIONS[name]


def list_configs() -> list[str]:
    """List all available configuration names.

    :return: Sorted list of configuration names.
    """
    return sorted(CONFIGURATIONS.keys())
