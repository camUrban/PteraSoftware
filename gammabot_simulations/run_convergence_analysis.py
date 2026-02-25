"""Run convergence analysis for GammaBot simulations.

This script demonstrates how to use the analyze_unsteady_convergence_non_trapezoidal
function to find converged simulation parameters for GammaBot wings.

Usage:
    python run_convergence_analysis.py [config_name] [options]

Examples:
    python run_convergence_analysis.py L170V_R180V_170Hz --visualize
    python run_convergence_analysis.py L170V_R180V_170Hz --num-cycles 2 5
    python run_convergence_analysis.py L170V_R180V_170Hz --num-cycles 6 10

    # Run two ranges in parallel (both write to the same cache file):
    python run_convergence_analysis.py L170V_R180V_170Hz --num-cycles 2 5 &
    python run_convergence_analysis.py L170V_R180V_170Hz --num-cycles 6 10 &
"""

import argparse
import sys
from pathlib import Path
from typing import cast

# Add project root to path for pterasoftware import
sys.path.insert(0, str(Path(__file__).parent.parent))

import configs
from convergence_helpers import create_gammabot_geometry_resampler

import pterasoftware as ps


def compute_flapping_parameters(wing_params: dict, flapping_period: float) -> dict:
    """Compute derived flapping parameters from wing configuration.

    :param wing_params: Dictionary with phi_max, psi_max, delta, etc.
    :param flapping_period: Flapping period in seconds.
    :return: Dictionary with amplitude, period, and phase values.
    """
    phi_max = wing_params["phi_max"]
    psi_max = wing_params["psi_max"]
    delta = wing_params["delta"]

    amplitude_x = phi_max
    amplitude_y = psi_max

    period_x = flapping_period if amplitude_x != 0.0 else 0.0
    period_y = flapping_period if amplitude_y != 0.0 else 0.0

    phase_y = (90.0 + delta) if amplitude_y != 0.0 else 0.0

    return {
        "amplitude_x": amplitude_x,
        "amplitude_y": amplitude_y,
        "period_x": period_x,
        "period_y": period_y,
        "phase_y": phase_y,
    }


def create_wing_cross_sections(
    wing_section_data, num_wing_cross_sections: int
) -> list[ps.geometry.wing_cross_section.WingCrossSection]:
    """Create WingCrossSection objects from wing section data.

    :param wing_section_data: Array with position and chord data.
    :param num_wing_cross_sections: Number of cross sections to create.
    :return: List of WingCrossSection objects.
    """
    cross_sections = []

    for i in range(num_wing_cross_sections):
        num_spanwise_panels = 1 if i < (num_wing_cross_sections - 1) else None

        cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            Lp_Wcsp_Lpp=tuple(wing_section_data[i, :3]),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            chord=float(wing_section_data[i, 3]),
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            num_spanwise_panels=num_spanwise_panels,
            control_surface_symmetry_type=None,
            control_surface_hinge_point=0.75,
            control_surface_deflection=0.0,
        )
        cross_sections.append(cross_section)

    return cross_sections


def create_wing_cross_section_movements(
    wing: ps.geometry.wing.Wing, num_cross_sections: int
) -> list[ps.movements.wing_cross_section_movement.WingCrossSectionMovement]:
    """Create WingCrossSectionMovement objects for each cross section.

    :param wing: The Wing containing the cross sections.
    :param num_cross_sections: Number of cross sections in the wing.
    :return: List of WingCrossSectionMovement objects.
    """
    movements = []
    for i in range(num_cross_sections):
        movement = ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=wing.wing_cross_sections[i],
        )
        movements.append(movement)
    return movements


def create_ref_problem(
    config_name: str,
    num_spanwise_sections: int = 8,
    num_chordwise_panels: int = 4,
    num_cycles: int = 2,
) -> ps.problems.UnsteadyProblem:
    """Create a reference UnsteadyProblem for convergence analysis.

    :param config_name: Name of the GammaBot configuration to use.
    :param num_spanwise_sections: Initial number of spanwise sections.
    :param num_chordwise_panels: Initial number of chordwise panels.
    :param num_cycles: Number of flapping cycles.
    :return: UnsteadyProblem configured for the specified parameters.
    """
    import dxf_to_csv

    # Load configuration
    wing_config = configs.get_config(config_name)
    shared = configs.SHARED_PARAMS

    # Extract parameters
    velocity = cast(float, shared["velocity"])
    alpha = cast(float, shared["alpha"])
    flapping_frequency = cast(float, shared["flapping_frequency"])
    wing_spacing = cast(float, shared["wing_spacing"])
    x_offset = shared["x_offset"]
    y_offset = shared["y_offset"]
    chordwise_spacing = cast(str, shared["chordwise_spacing"])

    flapping_period: float = 1.0 / flapping_frequency
    num_wing_cross_sections: int = num_spanwise_sections + 1

    # Compute flapping parameters for each wing
    left_params = wing_config["left"]
    right_params = wing_config["right"]

    left_flapping = compute_flapping_parameters(left_params, flapping_period)
    right_flapping = compute_flapping_parameters(right_params, flapping_period)

    # Load wing geometry from DXF
    dxf_filepath = Path(__file__).parent / "gammabot_approximate_wing.dxf"
    wing_section_data = dxf_to_csv.process_dxf_to_wing_section_data(
        str(dxf_filepath), num_spanwise_sections
    )

    # Create wing cross sections
    left_wing_cross_sections = create_wing_cross_sections(
        wing_section_data, num_wing_cross_sections
    )
    right_wing_cross_sections = create_wing_cross_sections(
        wing_section_data, num_wing_cross_sections
    )

    # Define the GammaBot Airplane
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
    )

    # Create wing cross section movements
    left_wcs_movements = create_wing_cross_section_movements(
        airplane.wings[0], num_wing_cross_sections
    )
    right_wcs_movements = create_wing_cross_section_movements(
        airplane.wings[1], num_wing_cross_sections
    )

    # Define wing movements
    left_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=left_wcs_movements,
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
        wing_cross_section_movements=right_wcs_movements,
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

    # Define airplane movement
    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[left_wing_movement, right_wing_movement],
    )

    # Define operating point
    operating_point = ps.operating_point.OperatingPoint(vCg__E=velocity, alpha=alpha)
    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point,
        )
    )

    # Define overall movement
    # Use 50 steps per cycle for the reference problem
    delta_time = flapping_period / 30  # 30 steps per cycle
    movement = ps.movements.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_cycles=num_cycles,
        delta_time=delta_time,
    )

    # Create and return the problem
    return ps.problems.UnsteadyProblem(movement=movement)


def run_convergence_analysis(
    config_name: str,
    visualize: bool = False,
    clear_cache: bool = False,
    prescribed_wake: bool = True,
    free_wake: bool = False,
    num_cycles_bounds: tuple[int, int] = (1, 1),
    panel_ar_bounds: tuple[int, int] = (4, 1),
    chordwise_bounds: tuple[int, int] = (5, 10),
) -> tuple[bool | None, int | None, int | None, int | None]:
    """Run convergence analysis for a GammaBot configuration.

    :param config_name: Name of the GammaBot configuration to analyze.
    :param visualize: If True, save mesh visualizations.
    :param clear_cache: If True, delete the cache file before running.
    :param prescribed_wake: If True, analyze prescribed wake. Default is True.
    :param free_wake: If True, analyze free wake. Default is False.
    :param num_cycles_bounds: Range of wake lengths in cycles. Default is (1, 1).
    :param panel_ar_bounds: Range of panel aspect ratios (coarsest, finest). Default is
        (4, 1).
    :param chordwise_bounds: Range of chordwise panel counts. Default is (5, 10).
    :return: Tuple of (converged_wake, converged_cycles, converged_ar, converged_chordwise).
    """
    print(f"Running convergence analysis for configuration: {config_name}")
    print(f"  prescribed_wake={prescribed_wake}, free_wake={free_wake}")
    print(f"  num_cycles_bounds={num_cycles_bounds}")
    print(f"  panel_ar_bounds={panel_ar_bounds}")
    print(f"  chordwise_bounds={chordwise_bounds}")
    print("-" * 60)

    # Set up logging
    ps._logging.set_up_logging(level="INFO")

    # Set up cache file
    cache_dir = Path(__file__).parent / "convergence_cache"
    cache_file = cache_dir / f"{config_name}.json"
    if clear_cache and cache_file.exists():
        cache_file.unlink()
        print(f"Cleared cache file: {cache_file}")

    # Create reference problem
    print("Creating reference problem...")
    ref_problem = create_ref_problem(config_name)

    # Create geometry resampler
    dxf_filepath = Path(__file__).parent / "gammabot_approximate_wing.dxf"
    resampler = create_gammabot_geometry_resampler(dxf_filepath)

    # Set up visualization directory if needed
    vis_dir = None
    if visualize:
        vis_dir = str(Path(__file__).parent / "convergence_meshes" / config_name)

    # Run convergence analysis
    print("Starting convergence analysis...")
    print("This may take a while depending on the parameter ranges.")
    print()

    result = ps.convergence.analyze_unsteady_convergence_non_trapezoidal_optimized_dt(
        ref_problem=ref_problem,
        wing_geometry_resampler=resampler,
        prescribed_wake=prescribed_wake,
        free_wake=free_wake,
        num_cycles_bounds=num_cycles_bounds,
        panel_aspect_ratio_bounds=panel_ar_bounds,
        num_chordwise_panels_bounds=chordwise_bounds,
        visualize_meshes=visualize,
        visualization_dir=vis_dir,
        coefficient_mask=(True, False, True, False, True, False),
        rtol=0.05,
        atol=0.1,
        cache_file=cache_file,
    )

    # Print results
    print()
    print("=" * 60)
    print("CONVERGENCE ANALYSIS RESULTS")
    print("=" * 60)

    if result[0] is not None:
        (
            converged_wake,
            converged_cycles,
            converged_ar,
            converged_chordwise,
        ) = result
        print(f"Wake type:       {'Prescribed' if converged_wake else 'Free'}")
        print(f"Wake length:     {converged_cycles} cycles")
        print(f"Panel AR:        {converged_ar}")
        print(f"Chordwise panels: {converged_chordwise}")
    else:
        print("No converged parameters found within the specified bounds.")
        print(
            "Consider expanding the parameter ranges or relaxing the convergence criteria."
        )

    return result


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Run convergence analysis for GammaBot simulations."
    )
    parser.add_argument(
        "config_name",
        nargs="?",
        default="L170V_R180V_170Hz",
        help="Configuration name (default: L170V_R180V_170Hz)",
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Save mesh visualizations for each parameter combination",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available configurations",
    )
    parser.add_argument(
        "--clear-cache",
        action="store_true",
        help="Delete the cache file before running",
    )
    parser.add_argument(
        "--prescribed-wake",
        action="store_true",
        default=None,
        help="Analyze prescribed wake (default if neither wake flag is set)",
    )
    parser.add_argument(
        "--free-wake",
        action="store_true",
        default=None,
        help="Analyze free wake",
    )
    parser.add_argument(
        "--num-cycles",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=(1, 1),
        help="Range of wake lengths in cycles (default: 1 1)",
    )
    parser.add_argument(
        "--panel-ar",
        type=int,
        nargs=2,
        metavar=("COARSEST", "FINEST"),
        default=(4, 1),
        help="Range of panel aspect ratios, coarsest to finest (default: 4 1)",
    )
    parser.add_argument(
        "--chordwise",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=(5, 10),
        help="Range of chordwise panel counts (default: 5 10)",
    )

    args = parser.parse_args()

    if args.list:
        print("Available configurations:")
        for name in configs.CONFIGURATIONS.keys():
            print(f"  {name}")
        return

    # Default to prescribed_wake=True, free_wake=False if neither flag is set.
    if args.prescribed_wake is None and args.free_wake is None:
        prescribed_wake = True
        free_wake = False
    else:
        prescribed_wake = bool(args.prescribed_wake)
        free_wake = bool(args.free_wake)

    run_convergence_analysis(
        config_name=args.config_name,
        visualize=args.visualize,
        clear_cache=args.clear_cache,
        prescribed_wake=prescribed_wake,
        free_wake=free_wake,
        num_cycles_bounds=tuple(args.num_cycles),
        panel_ar_bounds=tuple(args.panel_ar),
        chordwise_bounds=tuple(args.chordwise),
    )


if __name__ == "__main__":
    main()
