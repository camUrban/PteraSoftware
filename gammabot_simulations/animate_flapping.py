"""Animate wing flapping for NC=12,13,14,15 at AR=1 (no simulation needed).

Usage:
    python animate_flapping.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent))

import json
import os

import configs
import dxf_to_csv
import numpy as np

import pterasoftware as ps
from pterasoftware.convergence.unsteady_non_trapezoidal import (
    _get_num_cross_sections_for_panel_ar,
)

# -- Parameters ----------------------------------------------------------------
TARGET_AR = 1
CHORDWISE_VALUES = [12, 13, 14, 15]
CONFIG_NAME = "L170V_R180V_170Hz"
NUM_CYCLES = 1
_DT_OVERRIDES = {"1,15": 5.823762e-05}
# ------------------------------------------------------------------------------


def build_problem(nc: int):
    """Build an UnsteadyProblem at the given chordwise panel count with AR=1."""
    wing_config = configs.get_config(CONFIG_NAME)
    shared = configs.SHARED_PARAMS

    velocity = float(shared["velocity"])
    alpha = float(shared["alpha"])
    flapping_frequency = float(shared["flapping_frequency"])
    wing_spacing = float(shared["wing_spacing"])
    x_offset = shared["x_offset"]
    y_offset = shared["y_offset"]
    chordwise_spacing = str(shared["chordwise_spacing"])

    flapping_period = 1.0 / flapping_frequency
    dxf_filepath = Path(__file__).parent / "gammabot_approximate_wing.dxf"

    ref_data = dxf_to_csv.process_dxf_to_wing_section_data(str(dxf_filepath), 8)
    span = float(np.sum(np.sqrt(np.sum(ref_data[1:, :3] ** 2, axis=1))))
    avg_chord = float(np.mean(ref_data[:, 3]))
    num_sections = _get_num_cross_sections_for_panel_ar(span, avg_chord, TARGET_AR, nc)

    cache_path = Path(__file__).parent / "convergence_cache" / f"{CONFIG_NAME}.json"
    with open(cache_path, "r") as f:
        cache_data = json.load(f)
    dt_cache = cache_data.get("_optimized_dt", {})
    dt_key = f"{NUM_CYCLES},{nc}"
    delta_time = dt_cache.get(dt_key, _DT_OVERRIDES.get(dt_key, flapping_period / 30))

    wing_section_data = dxf_to_csv.process_dxf_to_wing_section_data(
        str(dxf_filepath), num_sections
    )
    num_wcs = num_sections + 1

    def make_cross_sections():
        cross_sections = []
        for i in range(num_wcs):
            num_spanwise_panels = 1 if i < num_sections else None
            wcs = ps.geometry.wing_cross_section.WingCrossSection(
                Lp_Wcsp_Lpp=tuple(wing_section_data[i, :3]),
                angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                chord=float(wing_section_data[i, 3]),
                airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
                num_spanwise_panels=num_spanwise_panels,
            )
            cross_sections.append(wcs)
        return cross_sections

    left_params = wing_config["left"]
    right_params = wing_config["right"]

    def compute_flapping(params):
        phi_max = params["phi_max"]
        psi_max = params["psi_max"]
        delta = params["delta"]
        return {
            "amp_x": phi_max,
            "amp_y": psi_max,
            "period_x": flapping_period if phi_max != 0.0 else 0.0,
            "period_y": flapping_period if psi_max != 0.0 else 0.0,
            "phase_y": (90.0 + delta) if psi_max != 0.0 else 0.0,
        }

    left_flap = compute_flapping(left_params)
    right_flap = compute_flapping(right_params)

    airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=make_cross_sections(),
                Ler_Gs_Cgs=(0.0, wing_spacing / 2, 0.0),
                angles_Gs_to_Wn_ixyz=(
                    left_params["phi_v_shift"],
                    left_params["psi_v_shift"],
                    0.0,
                ),
                symmetric=False,
                mirror_only=False,
                num_chordwise_panels=nc,
                chordwise_spacing=chordwise_spacing,
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=make_cross_sections(),
                Ler_Gs_Cgs=(0.0, wing_spacing / 2, 0.0),
                angles_Gs_to_Wn_ixyz=(
                    right_params["phi_v_shift"],
                    right_params["psi_v_shift"],
                    0.0,
                ),
                symmetric=False,
                mirror_only=True,
                symmetryNormal_G=(0, 1, 0),
                symmetryPoint_G_Cg=(0, 0, 0),
                num_chordwise_panels=nc,
                chordwise_spacing=chordwise_spacing,
            ),
        ],
        name="GammaBot",
    )

    def make_wcs_movements(wing):
        return [
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
            )
            for wcs in wing.wing_cross_sections
        ]

    left_wcsm = make_wcs_movements(airplane.wings[0])
    right_wcsm = make_wcs_movements(airplane.wings[1])

    left_wm = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=left_wcsm,
        rotationPointOffset_Gs_Ler=(x_offset, y_offset, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(left_flap["amp_x"], left_flap["amp_y"], 0.0),
        periodAngles_Gs_to_Wn_ixyz=(left_flap["period_x"], left_flap["period_y"], 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, left_flap["phase_y"], 0.0),
    )
    right_wm = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[1],
        wing_cross_section_movements=right_wcsm,
        rotationPointOffset_Gs_Ler=(x_offset, y_offset, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(right_flap["amp_x"], right_flap["amp_y"], 0.0),
        periodAngles_Gs_to_Wn_ixyz=(
            right_flap["period_x"],
            right_flap["period_y"],
            0.0,
        ),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, right_flap["phase_y"], 0.0),
    )

    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[left_wm, right_wm],
    )

    operating_point = ps.operating_point.OperatingPoint(vCg__E=velocity, alpha=alpha)
    op_movement = ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=operating_point,
    )

    movement = ps.movements.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=op_movement,
        num_cycles=NUM_CYCLES,
        delta_time=delta_time,
    )

    return ps.problems.UnsteadyProblem(movement=movement)


def main():
    for nc in CHORDWISE_VALUES:
        print(f"Animating NC={nc}...")
        problem = build_problem(nc)

        solver = ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_problem=problem
        )

        print(f"  {solver.num_panels} panels, {solver.num_steps} steps")

        ps.output.animate(
            unsteady_solver=solver,
            scalar_type=None,
            show_wake_vortices=False,
            save=True,
        )

        # Rename the hardcoded output file to include the NC value.
        src = "Animate.webp"
        dst = str(Path(__file__).parent / f"animate_NC{nc}_AR{TARGET_AR}.webp")
        if os.path.exists(src):
            os.replace(src, dst)
            print(f"  Saved to {dst}")


if __name__ == "__main__":
    main()
