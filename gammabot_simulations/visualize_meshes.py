"""Visualize wing meshes at specific convergence parameters.

Creates mesh visualizations for comparing panel geometries across different
chordwise panel counts at a fixed panel aspect ratio.

Usage:
    python visualize_meshes.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
sys.path.insert(0, str(Path(__file__).parent))

import dxf_to_csv
import numpy as np

import pterasoftware as ps
from pterasoftware.convergence._functions import (
    _verify_panel_aspect_ratio,
    _visualize_wing_mesh,
)
from pterasoftware.convergence.unsteady_non_trapezoidal import (
    _get_num_cross_sections_for_panel_ar,
)

import configs

# ── Parameters to visualize ──────────────────────────────────────────────────
TARGET_AR = 1
CHORDWISE_VALUES = [5, 10, 15]
CONFIG_NAME = "L170V_R180V_170Hz"
# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    # Load configuration.
    wing_config = configs.get_config(CONFIG_NAME)
    shared = configs.SHARED_PARAMS

    wing_spacing = float(shared["wing_spacing"])
    chordwise_spacing = str(shared["chordwise_spacing"])

    # Load DXF geometry to get span and avg_chord.
    dxf_filepath = Path(__file__).parent / "gammabot_approximate_wing.dxf"

    # Get wing geometry info at a reference resolution.
    ref_sections = 8
    ref_data = dxf_to_csv.process_dxf_to_wing_section_data(
        str(dxf_filepath), ref_sections
    )
    span = float(np.sum(np.sqrt(np.sum(ref_data[1:, :3] ** 2, axis=1))))
    avg_chord = float(np.mean(ref_data[:, 3]))

    print(f"Wing span: {span:.6f} m")
    print(f"Avg chord: {avg_chord:.6f} m")
    print()

    output_dir = Path(__file__).parent / "convergence_meshes" / CONFIG_NAME
    output_dir.mkdir(parents=True, exist_ok=True)

    for nc in CHORDWISE_VALUES:
        num_sections = _get_num_cross_sections_for_panel_ar(
            span, avg_chord, TARGET_AR, nc
        )
        print(f"NC={nc}: num_spanwise_sections={num_sections}")

        # Get resampled geometry.
        wing_section_data = dxf_to_csv.process_dxf_to_wing_section_data(
            str(dxf_filepath), num_sections
        )

        num_wcs = num_sections + 1

        # Create wing cross sections for left and right wings.
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

        # Verify AR.
        ar_ok, actual_ar = _verify_panel_aspect_ratio(airplane, TARGET_AR)
        print(f"  Target AR: {TARGET_AR}, Actual AR: {actual_ar:.2f}, OK: {ar_ok}")

        # Count panels.
        total_panels = 0
        for wing in airplane.wings:
            if wing.panels is not None:
                total_panels += np.ravel(wing.panels).size
        print(f"  Total panels: {total_panels}")

        # Visualize and save.
        _visualize_wing_mesh(
            airplane,
            title=f"AR={TARGET_AR}, Chordwise={nc}",
            show=True,
            save_path=str(output_dir / f"mesh_ar{TARGET_AR}_chord{nc}.png"),
        )
        print()


if __name__ == "__main__":
    main()
