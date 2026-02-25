"""Run coupled GammaBot UVLM/MuJoCo simulations with configurable parameters.

Couples Ptera Software's unsteady ring vortex lattice method aerodynamic solver with
MuJoCo rigid body dynamics to simulate the GammaBot's free flight on the water surface.

Usage:
    python run_coupled_simulation.py <config_name> [options]
    python run_coupled_simulation.py --list

Examples:
    python run_coupled_simulation.py L170V_R180V_170Hz
    python run_coupled_simulation.py L170V_R180V_170Hz --fine --tag baseline
    python run_coupled_simulation.py L170V_R180V_170Hz --prescribed-num-flaps 3
    python run_coupled_simulation.py --list

Results are saved to:
    simulation_results/<run_id>/<config_name>/<mesh_type>/
"""

import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import cast

# Add project root to path for pterasoftware import.
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import the local configs module first (which re-exports from gammabot_simulations).
import configs

# Add gammabot_simulations to path for helper imports (dxf_to_csv, run_simulation).
# Use append so the local configs.py keeps priority over gammabot_simulations/configs.py.
sys.path.append(str(Path(__file__).parent.parent / "gammabot_simulations"))

import dxf_to_csv

import matplotlib.pyplot as plt
import numpy as np

import pterasoftware as ps

from run_simulation import (
    compute_flapping_parameters,
    create_wing_cross_section_movements,
    create_wing_cross_sections,
)

# Configure logging.
logging.root.setLevel(logging.INFO)
for handler in logging.root.handlers:
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("[%(name)s] %(levelname)s: %(message)s"))

logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("numba").setLevel(logging.WARNING)

script_logger = logging.getLogger("coupled_gammabot")
script_logger.setLevel(logging.DEBUG)


def gammabot_body_drag(
    coupled_operating_point: ps.operating_point.CoupledOperatingPoint,
    airplane: ps.geometry.airplane.Airplane,
) -> tuple[np.ndarray, np.ndarray]:
    """Computes the speed dependent body drag force on the GammaBot.

    The drag model is F_drag = -c1 * mass * vCg__E^2, acting in the negative wind
    x axis direction. This represents hydrodynamic drag from the water surface.

    :param coupled_operating_point: The current CoupledOperatingPoint.
    :param airplane: The current Airplane.
    :return: A tuple of two (3,) ndarrays of floats. The first is the additional
        force (in wind axes, in Newtons). The second is the additional moment (in
        wind axes, relative to the first Airplane's CG, in Newton meters).
    """
    c1 = configs.COUPLED_PARAMS["body_drag_c1"]
    mass = airplane.weight / np.linalg.norm(coupled_operating_point.g_E)
    drag_FX_W = -c1 * mass * coupled_operating_point.vCg__E**2
    return np.array([drag_FX_W, 0.0, 0.0], dtype=float), np.zeros(3, dtype=float)


def build_extra_xml() -> dict[str, str]:
    """Builds the extra_xml dict for the GammaBot MuJoCo model.

    Includes default contact parameters, the gammabot.stl mesh asset, a ground plane,
    and a mesh geom for contact.

    :return: A dict mapping injection point names to XML fragment strings.
    """
    extra_xml = {
        "default": (
            "<default>\n"
            '  <geom condim="1"\n'
            '        solref="0.01 1"\n'
            '        solimp="0.99 0.999 0.001 0.5 2"/>\n'
            "</default>"
        ),
        "asset": (
            "<asset>\n"
            '  <mesh name="gammabot" file="gammabot.stl"'
            ' scale="0.001 0.001 0.001"/>\n'
            "</asset>"
        ),
        "worldbody": (
            '<geom size="1 1 .01" type="plane" margin="0.001"/>\n'
            '<light pos="0 0 .6"/>'
        ),
        "body": (
            '<geom type="mesh" mesh="gammabot" rgba="0.7 0.7 0.9 1"'
            ' contype="1" conaffinity="1" mass="0"/>'
        ),
        "visual": '<visual><global offwidth="1920" offheight="1088"/></visual>',
    }
    return extra_xml


def load_mujoco_assets() -> dict[str, bytes]:
    """Loads binary assets for the MuJoCo model.

    :return: A dict mapping virtual filenames to their binary contents.
    """
    stl_path = (
        Path(__file__).parent.parent / "mujoco_examples" / "assets" / "gammabot.stl"
    )
    return {"gammabot.stl": stl_path.read_bytes()}


def run_coupled_simulation(
    config_name: str,
    run_id: str,
    mesh_type: str = "coarse",
    prescribed_num_flaps: int = 2,
    free_num_flaps: int = 2,
    show_results: bool = True,
) -> None:
    """Run a coupled GammaBot simulation with the specified configuration.

    :param config_name: Name of the configuration to use.
    :param run_id: Identifier for this run (tag or timestamp). Used in output path.
    :param mesh_type: Mesh type to use ("coarse" or "fine").
    :param prescribed_num_flaps: Number of prescribed wake development flaps before
        free flight.
    :param free_num_flaps: Number of free flight flaps after the prescribed portion.
    :param show_results: If True, display results interactively.
    :return: None
    """
    script_logger.info(
        f"Initializing coupled simulation: {config_name} ({mesh_type} mesh)"
    )

    # Load configuration.
    wing_config = configs.get_config(config_name)
    shared = configs.SHARED_PARAMS
    mesh = configs.MESH_PARAMS[mesh_type]
    coupled = configs.COUPLED_PARAMS

    # Extract parameters.
    velocity = cast(float, shared["velocity"])
    alpha = cast(float, shared["alpha"])
    flapping_frequency = cast(float, shared["flapping_frequency"])
    wing_spacing = cast(float, shared["wing_spacing"])
    x_offset = shared["x_offset"]
    y_offset = shared["y_offset"]
    chordwise_spacing = cast(str, shared["chordwise_spacing"])

    num_chordwise_panels = cast(int, mesh["num_chordwise_panels"])
    num_spanwise_sections = cast(int, mesh["num_spanwise_sections"])
    delta_time: float = mesh["delta_time"]

    flapping_period: float = 1.0 / flapping_frequency
    num_wing_cross_sections: int = num_spanwise_sections + 1

    # Compute the number of time steps from flap counts.
    steps_per_flap = round(flapping_period / delta_time)
    prescribed_num_steps = prescribed_num_flaps * steps_per_flap
    free_num_steps = free_num_flaps * steps_per_flap

    script_logger.info(f"  Steps per flap: {steps_per_flap}")
    script_logger.info(f"  Prescribed steps: {prescribed_num_steps}")
    script_logger.info(f"  Free steps: {free_num_steps}")
    script_logger.info(f"  Total steps: {prescribed_num_steps + free_num_steps}")

    # Compute flapping parameters for each wing.
    left_params = wing_config["left"]
    right_params = wing_config["right"]

    left_flapping = compute_flapping_parameters(left_params, flapping_period)
    right_flapping = compute_flapping_parameters(right_params, flapping_period)

    # Load wing geometry from DXF.
    dxf_filepath = (
        Path(__file__).parent.parent
        / "gammabot_simulations"
        / "gammabot_approximate_wing.dxf"
    )
    wing_section_data = dxf_to_csv.process_dxf_to_wing_section_data(
        str(dxf_filepath), num_spanwise_sections
    )

    # Create wing cross sections.
    left_wing_cross_sections = create_wing_cross_sections(
        wing_section_data, num_wing_cross_sections
    )
    right_wing_cross_sections = create_wing_cross_sections(
        wing_section_data, num_wing_cross_sections
    )

    # Define the GammaBot Airplane with weight from coupled params.
    airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=left_wing_cross_sections,
                Ler_Gs_Cgs=(0.0, wing_spacing / 2, 0.0),
                angles_Gs_to_Wn_ixyz=(
                    left_params["phi_v_shift"],
                    left_params["psi_v_shift"],
                    0.0,
                ),
                symmetric=False,
                mirror_only=False,
                symmetryNormal_G=None,
                symmetryPoint_G_Cg=None,
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing,
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=right_wing_cross_sections,
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
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing=chordwise_spacing,
            ),
        ],
        name="GammaBot",
        weight=coupled["weight"],
    )

    # Create wing cross section movements.
    left_wing_cross_section_movements = create_wing_cross_section_movements(
        airplane.wings[0], num_wing_cross_sections
    )
    right_wing_cross_section_movements = create_wing_cross_section_movements(
        airplane.wings[1], num_wing_cross_sections
    )

    # Define wing movements.
    left_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=left_wing_cross_section_movements,
        rotationPointOffset_Gs_Ler=(x_offset, y_offset, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(
            left_flapping["amplitude_x"],
            left_flapping["amplitude_y"],
            0.0,
        ),
        periodAngles_Gs_to_Wn_ixyz=(
            left_flapping["period_x"],
            left_flapping["period_y"],
            0.0,
        ),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, left_flapping["phase_y"], 0.0),
    )

    right_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[1],
        wing_cross_section_movements=right_wing_cross_section_movements,
        rotationPointOffset_Gs_Ler=(x_offset, y_offset, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(
            right_flapping["amplitude_x"],
            right_flapping["amplitude_y"],
            0.0,
        ),
        periodAngles_Gs_to_Wn_ixyz=(
            right_flapping["period_x"],
            right_flapping["period_y"],
            0.0,
        ),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, right_flapping["phase_y"], 0.0),
    )

    # Define airplane movement.
    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[left_wing_movement, right_wing_movement],
    )

    # Define the initial coupled operating point.
    initial_coupled_operating_point = ps.operating_point.CoupledOperatingPoint(
        rho=coupled["rho"],
        vCg__E=velocity,
        alpha=alpha,
        beta=0.0,
        angles_E_to_BP1_izyx=(0.0, alpha, 0.0),
        externalFX_W=0.0,
        nu=coupled["nu"],
        g_E=coupled["g_E"],
    )

    # Define the coupled movement.
    coupled_movement = ps.movements.movement.CoupledMovement(
        airplane_movement=airplane_movement,
        initial_coupled_operating_point=initial_coupled_operating_point,
        delta_time=delta_time,
        prescribed_num_steps=prescribed_num_steps,
        free_num_steps=free_num_steps,
    )

    # Build the MuJoCo XML extensions and assets.
    extra_xml = build_extra_xml()
    mujoco_assets = load_mujoco_assets()

    # Define the coupled unsteady problem.
    coupled_problem = ps.problems.CoupledUnsteadyProblem(
        coupled_movement=coupled_movement,
        I_BP1_CgP1=coupled["I_BP1_CgP1"],
        external_forces_fn=gammabot_body_drag,
        extra_xml=extra_xml,
        mujoco_assets=mujoco_assets,
    )

    # Create the solver.
    solver = (
        ps.coupled_unsteady_ring_vortex_lattice_method.CoupledUnsteadyRingVortexLatticeMethodSolver(
            coupled_unsteady_problem=coupled_problem,
        )
    )

    script_logger.info("Starting coupled simulation.")
    script_logger.info(f"  Total time steps: {coupled_problem.num_steps}")
    script_logger.info(f"  Delta time: {coupled_problem.delta_time} s")
    script_logger.info(
        f"  Total duration: "
        f"{coupled_problem.num_steps * coupled_problem.delta_time:.4f} s"
    )

    # Verify the generated XML contains our injected elements.
    xml_str = solver.mujoco_model.xml_str
    assert "gammabot.stl" in xml_str, "Mesh asset not found in generated XML."
    assert 'type="plane"' in xml_str, "Ground plane not found in generated XML."
    script_logger.info("MuJoCo XML verified: ground plane and mesh geom present.")

    # Run the coupled simulation.
    prescribed_wake = cast(bool, mesh.get("prescribed_wake", False))
    solver.run(prescribed_wake=prescribed_wake, show_progress=True)

    script_logger.info("Coupled simulation completed successfully.")

    # Create output directory.
    mesh_dir = mesh_type.title()
    output_dir = (
        Path(__file__).parent / "simulation_results" / run_id / config_name / mesh_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    original_dir = Path.cwd()
    os.chdir(output_dir)
    script_logger.info(f"Saving results to: {output_dir}")

    try:
        # Plot position time history.
        total_steps = prescribed_num_steps + free_num_steps
        times = np.linspace(0, (total_steps - 1) * delta_time, total_steps)

        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        labels = ["x (m)", "y (m)", "z (m)"]
        for i, (ax, label) in enumerate(zip(axes, labels)):
            ax.plot(times, solver.stackPosition_E_E[:, i])
            ax.set_ylabel(label)
            ax.grid(True)
        axes[-1].set_xlabel("Time (s)")
        fig.suptitle(f"CG Position - {config_name}")
        fig.tight_layout()
        fig.savefig("position_history.png", dpi=150)
        if show_results:
            plt.show()
        else:
            plt.close(fig)

        # Print final position summary.
        final_pos = solver.stackPosition_E_E[-1]
        print(f"Final CG position (E): x={final_pos[0]:.6f}, "
              f"y={final_pos[1]:.6f}, z={final_pos[2]:.6f} m")

        # Animate free flight.
        ps.output.animate_free_flight(
            coupled_solver=solver,
            scalar_type="lift",
            show_wake_vortices=True,
            save=True,
            testing=not show_results,
        )

        # Draw the final state.
        ps.output.draw(
            solver=solver,
            scalar_type="lift",
            show_wake_vortices=True,
            save=True,
            testing=not show_results,
        )

    finally:
        os.chdir(original_dir)


def main() -> None:
    """Parse arguments and run the coupled simulation."""
    parser = argparse.ArgumentParser(
        description="Run coupled GammaBot UVLM/MuJoCo simulations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "config",
        nargs="?",
        help="Configuration name (e.g., L170V_R180V_170Hz)",
    )
    parser.add_argument(
        "--fine",
        action="store_true",
        help="Use fine mesh settings instead of coarse",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display plots interactively",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available configurations and exit",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Tag for this run (used in output path). Defaults to timestamp.",
    )
    parser.add_argument(
        "--prescribed-num-flaps",
        type=int,
        default=2,
        help="Number of prescribed wake development flaps (default: 2)",
    )
    parser.add_argument(
        "--free-num-flaps",
        type=int,
        default=2,
        help="Number of free flight flaps (default: 2)",
    )

    args = parser.parse_args()

    if args.list:
        print("Available configurations:")
        for name in configs.list_configs():
            print(f"  {name}")
        return

    if not args.config:
        parser.error("Configuration name required (or use --list)")

    run_id = args.tag if args.tag else datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    mesh_type = "fine" if args.fine else "coarse"

    run_coupled_simulation(
        config_name=args.config,
        run_id=run_id,
        mesh_type=mesh_type,
        prescribed_num_flaps=args.prescribed_num_flaps,
        free_num_flaps=args.free_num_flaps,
        show_results=not args.no_show,
    )


if __name__ == "__main__":
    main()
