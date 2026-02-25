"""Run GammaBot UVLM simulations with configurable parameters.

Usage:     python run_simulation.py <config_name> [--fine | --coarse-and-fine]
[--no-show] [--tag TAG]     python run_simulation.py --all --no-show
[--fine | --coarse-and-fine] [--tag TAG] [--sequential] [--workers N]
python run_simulation.py --list

Examples:     python run_simulation.py L0V_R150V_170Hz     python run_simulation.py
L170V_R180V_170Hz --fine --tag baseline     python run_simulation.py --all --no-show
--tag baseline     python run_simulation.py --all --no-show --sequential
python run_simulation.py --all --no-show --workers 2     python run_simulation.py --all
--coarse-and-fine --no-show --tag full_study     python run_simulation.py --list

Results are saved to: simulation_results/<run_id>/<config_name>/<mesh_type>/ where
run_id is either the --tag value or a timestamp (YYYY-MM-DD_HH-MM-SS), and mesh_type is
one of Coarse or Fine.

Parallel Execution:     --all runs simulations in parallel by default (5 workers max).
Use --sequential to disable parallel execution.     Use --workers N to override the
default worker count.
"""

import argparse
import io
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TextIO, cast

# Add project root to path for pterasoftware import
sys.path.insert(0, str(Path(__file__).parent.parent))

import configs
import dxf_to_csv

import pterasoftware as ps


@dataclass
class SimulationTask:
    """Represents a single simulation task to be executed."""

    config_name: str
    mesh_type: str  # "coarse" or "fine"
    run_id: str
    show_results: bool


@dataclass
class SimulationResult:
    """Result of executing a simulation task."""

    task: SimulationTask
    success: bool
    elapsed_time: float
    error_message: str | None = None


class TeeWriter:
    """Write to both a file and the original stdout."""

    def __init__(self, file: TextIO, original_stdout: TextIO) -> None:
        """Initialize the TeeWriter.

        :param file: File to write to.
        :param original_stdout: Original stdout to also write to.
        """
        self.file = file
        self.original_stdout = original_stdout

    def write(self, data: str) -> int:
        """Write data to both file and stdout.

        :param data: Data to write.
        :return: Number of characters written.
        """
        self.file.write(data)
        self.original_stdout.write(data)
        return len(data)

    def flush(self) -> None:
        """Flush both file and stdout."""
        self.file.flush()
        self.original_stdout.flush()


def create_codebase_snapshot(output_dir: Path) -> None:
    """Create a snapshot of the codebase in the output directory.

    Copies: - pterasoftware/ (excluding __pycache__ and _airfoils directories) - Files
    directly in gammabot_simulations/ (not subdirectories)

    :param output_dir: Directory to create the snapshot in.
    """
    snapshot_dir = output_dir / "codebase_snapshot"
    snapshot_dir.mkdir(parents=True, exist_ok=True)

    project_root = Path(__file__).parent.parent
    gammabot_dir = Path(__file__).parent

    # Copy pterasoftware package (excluding __pycache__ and _airfoils)
    pterasoftware_src = project_root / "pterasoftware"
    pterasoftware_dst = snapshot_dir / "pterasoftware"

    def ignore_patterns(_directory: str, files: list[str]) -> list[str]:
        """Return files/directories to ignore during copy."""
        ignored = []
        for f in files:
            if f == "__pycache__":
                ignored.append(f)
            elif f == "_airfoils":
                ignored.append(f)
        return ignored

    if pterasoftware_src.exists():
        shutil.copytree(pterasoftware_src, pterasoftware_dst, ignore=ignore_patterns)

    # Copy files directly in gammabot_simulations (not subdirectories)
    gammabot_snapshot_dir = snapshot_dir / "gammabot_simulations"
    gammabot_snapshot_dir.mkdir(parents=True, exist_ok=True)

    # Exclude plotting script and the script with its experimental thrust data. These
    # are saved once plotting is run.
    excluded_files = {"plot_single_wing_thrust.py", "experimental_thrust.py"}
    for item in gammabot_dir.iterdir():
        if item.is_file() and item.name not in excluded_files:
            shutil.copy2(item, gammabot_snapshot_dir / item.name)


def format_elapsed_time(seconds: float) -> str:
    """Format elapsed time in a human-readable format.

    :param seconds: Elapsed time in seconds.
    :return: Formatted string (e.g., "2h 15m 30s" or "45.2s").
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}h {minutes}m {secs:.1f}s"


def is_single_wing_config(config_name: str) -> bool:
    """Check if a configuration has only one wing flapping (the other at 0V).

    :param config_name: Configuration name (e.g., "L150V_R0V_170Hz").
    :return: True if only one wing is flapping, False if both are flapping.
    """
    import re

    pattern = r"L(\d+)V_R(\d+)V_170Hz"
    match = re.match(pattern, config_name)
    if not match:
        return False
    left_v, right_v = int(match.group(1)), int(match.group(2))
    return left_v == 0 or right_v == 0


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


def create_wing_cross_sections(wing_section_data, num_wing_cross_sections: int) -> list:
    """Create WingCrossSection objects from wing section data.

    :param wing_section_data: Array with position and chord data.
    :param num_wing_cross_sections: Number of cross sections to create.
    :return: List of WingCrossSection objects.
    """
    cross_sections = []

    for i in range(num_wing_cross_sections):
        num_spanwise_panels = 1 if i < (num_wing_cross_sections - 1) else None

        cross_section = ps.geometry.wing_cross_section.WingCrossSection(
            Lp_Wcsp_Lpp=wing_section_data[i, :3],
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            chord=wing_section_data[i, 3],
            airfoil=ps.geometry.airfoil.Airfoil(name="naca0012"),
            num_spanwise_panels=num_spanwise_panels,
            control_surface_symmetry_type=None,
            control_surface_hinge_point=0.75,
            control_surface_deflection=0.0,
        )
        cross_sections.append(cross_section)

    return cross_sections


def create_wing_cross_section_movements(wing, num_cross_sections: int) -> list:
    """Create WingCrossSectionMovement objects for a wing.

    :param wing: Wing object.
    :param num_cross_sections: Number of cross sections.
    :return: List of WingCrossSectionMovement objects.
    """
    movements = []
    for j in range(num_cross_sections):
        movement = ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=wing.wing_cross_sections[j]
        )
        movements.append(movement)
    return movements


def run_simulation(
    config_name: str,
    run_id: str,
    mesh_type: str = "coarse",
    show_results: bool = True,
) -> None:
    """Run a GammaBot simulation with the specified configuration.

    :param config_name: Name of the configuration to use.
    :param run_id: Identifier for this run (tag or timestamp). Used in output path.
    :param mesh_type: Mesh type to use ("coarse" or "fine").
    :param show_results: If True, display results interactively. If False, doesn't show
        plots and automatically closes renders after 1 second.
    :return: None
    """
    print(f"Initializing simulation: {config_name} ({mesh_type} mesh)")

    # Load configuration
    wing_config = configs.get_config(config_name)
    shared = configs.SHARED_PARAMS
    mesh = configs.MESH_PARAMS[mesh_type]

    # Extract parameters with explicit type casts for mypy
    velocity = cast(float, shared["velocity"])
    alpha = cast(float, shared["alpha"])
    flapping_frequency = cast(float, shared["flapping_frequency"])
    wing_spacing = cast(float, shared["wing_spacing"])
    x_offset = shared["x_offset"]
    y_offset = shared["y_offset"]
    chordwise_spacing = cast(str, shared["chordwise_spacing"])

    num_flaps = cast(int, mesh["num_flaps"])
    num_chordwise_panels = cast(int, mesh["num_chordwise_panels"])
    num_spanwise_sections = cast(int, mesh["num_spanwise_sections"])
    delta_time: float = mesh["delta_time"]
    prescribed_wake = cast(bool, mesh["prescribed_wake"])

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
    movement = ps.movements.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        num_cycles=num_flaps,
        delta_time=delta_time,
    )

    # Define and run the solver
    problem = ps.problems.UnsteadyProblem(movement=movement)
    solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_problem=problem,
        )
    )

    solver.run(prescribed_wake=prescribed_wake)

    # Create output directory: simulation_results/<run_id>/<config_name>/<mesh_dir>/
    # Convert mesh_type to directory name (e.g., "coarse" -> "Coarse")
    mesh_dir = mesh_type.title().replace(" ", "-")
    output_dir = (
        Path(__file__).parent / "simulation_results" / run_id / config_name / mesh_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    # Change to output directory for saving results
    original_dir = Path.cwd()
    os.chdir(output_dir)
    print(f"Saving results to: {output_dir}")

    try:
        # Generate outputs
        ps.output.plot_results_versus_time(
            unsteady_solver=solver,
            show=show_results,
            save=True,
        )

        ps.output.draw(
            solver=solver,
            show_wake_vortices=True,
            scalar_type="induced drag",
            save=True,
            testing=not show_results,
        )

        ps.output.animate(
            unsteady_solver=solver,
            show_wake_vortices=True,
            scalar_type="induced drag",
            save=True,
            testing=not show_results,
        )

        # Capture print_results output and save to file
        results_buffer = io.StringIO()
        with redirect_stdout(results_buffer):
            ps.output.print_results(solver=solver)
        results_text = results_buffer.getvalue()

        # Print to console
        print(results_text)

        # Save to file
        results_file = output_dir / "print_results.txt"
        results_file.write_text(results_text)
        print(f"Results saved to: {results_file}")
    finally:
        # Restore original working directory
        os.chdir(original_dir)


def _calculate_max_workers() -> int:
    """Calculate the maximum number of workers for parallel execution.

    Uses a conservative default based on available CPUs, capped at 5 to avoid memory
    issues (each fine simulation uses ~1.5 GB RAM).

    :return: Number of workers to use.
    """
    cpu_count = os.cpu_count() or 5
    return min(cpu_count // 2, 5)


def _run_simulation_worker(task: SimulationTask) -> SimulationResult:
    """Worker function for parallel simulation execution.

    Suppresses stdout from the simulation to avoid cluttering the console.

    :param task: SimulationTask to execute.
    :return: SimulationResult with success status and timing.
    """
    # Print start message before redirecting stdout
    print(f"[STARTED] {task.config_name} ({task.mesh_type})", flush=True)

    start_time = time.perf_counter()

    try:
        # Suppress stdout and stderr from the simulation to avoid console clutter
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            run_simulation(
                config_name=task.config_name,
                run_id=task.run_id,
                mesh_type=task.mesh_type,
                show_results=task.show_results,
            )

        elapsed = time.perf_counter() - start_time
        return SimulationResult(task=task, success=True, elapsed_time=elapsed)

    except Exception as e:
        elapsed = time.perf_counter() - start_time
        return SimulationResult(
            task=task, success=False, elapsed_time=elapsed, error_message=str(e)
        )


def run_batch_parallel(
    tasks: list[SimulationTask], max_workers: int
) -> list[SimulationResult]:
    """Run a batch of simulations in parallel using ProcessPoolExecutor.

    :param tasks: List of SimulationTask objects to execute.
    :param max_workers: Maximum number of parallel workers.
    :return: List of SimulationResult objects.
    :raises KeyboardInterrupt: Re-raised after cleanup if user interrupts.
    """
    results: list[SimulationResult] = []
    total = len(tasks)

    print(f"Starting parallel execution with {max_workers} workers")
    print(f"{'=' * 60}")

    executor = ProcessPoolExecutor(max_workers=max_workers)
    try:
        # Submit all tasks
        future_to_task = {
            executor.submit(_run_simulation_worker, task): task for task in tasks
        }

        # Process completions as they happen
        for future in as_completed(future_to_task):
            result = future.result()
            results.append(result)

            completed = len(results)

            if result.success:
                print(
                    f"[OK] {result.task.config_name} ({result.task.mesh_type}) - "
                    f"{format_elapsed_time(result.elapsed_time)} "
                    f"[{completed}/{total}]"
                )
            else:
                print(
                    f"[FAILED] {result.task.config_name} ({result.task.mesh_type}) - "
                    f"{format_elapsed_time(result.elapsed_time)} "
                    f"[{completed}/{total}]"
                )
                print(f"         Error: {result.error_message}")

    except KeyboardInterrupt:
        print("\n\nInterrupted! Cancelling pending tasks...")
        executor.shutdown(wait=False, cancel_futures=True)
        raise

    finally:
        executor.shutdown(wait=True)

    return results


def run_batch_sequential(tasks: list[SimulationTask]) -> list[SimulationResult]:
    """Run a batch of simulations sequentially.

    :param tasks: List of SimulationTask objects to execute.
    :return: List of SimulationResult objects.
    """
    results: list[SimulationResult] = []
    total = len(tasks)

    for i, task in enumerate(tasks, 1):
        print(f"\n{'=' * 60}")
        print(f"Running simulation {i}/{total}: {task.config_name} ({task.mesh_type})")
        print(f"{'=' * 60}\n")

        start_time = time.perf_counter()
        try:
            run_simulation(
                config_name=task.config_name,
                run_id=task.run_id,
                mesh_type=task.mesh_type,
                show_results=task.show_results,
            )
            elapsed = time.perf_counter() - start_time
            print(f"Simulation completed in {format_elapsed_time(elapsed)}")
            results.append(
                SimulationResult(task=task, success=True, elapsed_time=elapsed)
            )
        except Exception as e:
            elapsed = time.perf_counter() - start_time
            print(
                f"ERROR: {task.config_name} ({task.mesh_type}) failed after "
                f"{format_elapsed_time(elapsed)}: {e}"
            )
            results.append(
                SimulationResult(
                    task=task, success=False, elapsed_time=elapsed, error_message=str(e)
                )
            )

    return results


def print_batch_summary(
    results: list[SimulationResult], batch_elapsed: float, parallel: bool
) -> None:
    """Print a summary of batch execution results.

    :param results: List of SimulationResult objects.
    :param batch_elapsed: Total batch execution time in seconds.
    :param parallel: Whether parallel execution was used.
    """
    total = len(results)
    succeeded = sum(1 for r in results if r.success)
    failed = [r for r in results if not r.success]

    total_sim_time = sum(r.elapsed_time for r in results)
    avg_sim_time = total_sim_time / total if total > 0 else 0

    print(f"\n{'=' * 60}")
    print("BATCH EXECUTION SUMMARY")
    print(f"{'=' * 60}")
    print(f"Completed: {succeeded}/{total}")
    print(f"Total time: {format_elapsed_time(batch_elapsed)}")
    print(f"Average simulation time: {format_elapsed_time(avg_sim_time)}")

    if parallel and total > 0:
        speedup = total_sim_time / batch_elapsed if batch_elapsed > 0 else 0
        print(f"Parallel speedup: {speedup:.2f}x")

    if failed:
        print(f"\nFAILED SIMULATIONS:")
        for result in failed:
            print(
                f"  - {result.task.config_name} ({result.task.mesh_type}): "
                f"{result.error_message}"
            )


def main():
    """Parse arguments and run simulation."""
    parser = argparse.ArgumentParser(
        description="Run GammaBot UVLM simulations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "config",
        nargs="?",
        help="Configuration name (e.g., L0V_R150V_170Hz)",
    )

    # Mesh selection (mutually exclusive)
    mesh_group = parser.add_mutually_exclusive_group()
    mesh_group.add_argument(
        "--fine",
        action="store_true",
        help="Use fine mesh settings instead of coarse",
    )
    mesh_group.add_argument(
        "--coarse-and-fine",
        action="store_true",
        help="Run both coarse and fine mesh for each configuration",
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
        "--all",
        action="store_true",
        help="Run all configurations (parallel by default, use --sequential for serial)",
    )
    parser.add_argument(
        "--all-single-wing",
        action="store_true",
        help="Run all single-wing configurations (exclude combined left+right flapping)",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Tag for this run (used in output path). Defaults to timestamp.",
    )
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="Disable parallel execution, run simulations one-by-one",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Override auto-detected worker count (default: 5)",
    )

    args = parser.parse_args()

    # Validation
    run_all = args.all or args.all_single_wing
    if run_all and not args.no_show:
        parser.error(
            "--all/--all-single-wing requires --no-show "
            "(parallel execution cannot show plots)"
        )

    if args.all and args.all_single_wing:
        parser.error("--all and --all-single-wing are mutually exclusive")

    if args.sequential and args.workers is not None:
        parser.error("--workers cannot be combined with --sequential")

    if args.list:
        print("Available configurations:")
        for name in configs.list_configs():
            print(f"  {name}")
        return

    if not run_all and not args.config:
        parser.error("Configuration name required (or use --all or --list)")

    # Generate run_id once for the entire batch (tag or timestamp)
    run_id = args.tag if args.tag else datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Determine which mesh types to run
    coarse_and_fine = getattr(args, "coarse_and_fine", False)
    if coarse_and_fine:
        mesh_settings = ["coarse", "fine"]
    elif args.fine:
        mesh_settings = ["fine"]
    else:
        mesh_settings = ["coarse"]

    # Create run output directory and set up logging/snapshot
    run_output_dir = Path(__file__).parent / "simulation_results" / run_id
    run_output_dir.mkdir(parents=True, exist_ok=True)

    # Set up console output capture
    console_log_path = run_output_dir / "console_output.txt"
    console_log_file = open(console_log_path, "w", encoding="utf-8")
    original_stdout = sys.stdout
    sys.stdout = TeeWriter(console_log_file, original_stdout)

    try:
        # Create codebase snapshot
        print(
            f"Creating codebase snapshot in {run_output_dir / 'codebase_snapshot'}..."
        )
        create_codebase_snapshot(run_output_dir)
        print("Codebase snapshot created.")

        if args.all:
            all_configs = configs.list_configs()
            num_configs = len(all_configs)
            num_meshes = len(mesh_settings)
            total = num_configs * num_meshes

            # Build task list for all mesh types
            all_tasks: list[SimulationTask] = []
            for name in all_configs:
                for mesh_type in mesh_settings:
                    task = SimulationTask(
                        config_name=name,
                        mesh_type=mesh_type,
                        run_id=run_id,
                        show_results=False,  # Always False for batch runs
                    )
                    all_tasks.append(task)

            print(f"Run ID: {run_id}")
            if num_meshes > 1:
                mesh_label = " + ".join(mesh_settings)
                print(
                    f"Running {num_configs} configurations x {num_meshes} mesh types "
                    f"= {total} simulations ({mesh_label})"
                )
            else:
                print(f"Running {num_configs} configurations ({mesh_settings[0]} mesh)")

            batch_start_time = time.perf_counter()
            results: list[SimulationResult] = []

            # Run simulations with warmup strategy
            use_parallel = not args.sequential
            if use_parallel:
                max_workers = args.workers if args.workers else _calculate_max_workers()

                # Run first task sequentially to warm up Numba cache
                print(f"\n{'=' * 60}")
                print("Running first simulation to warm up Numba cache...")
                print(f"{'=' * 60}")
                warmup_results = run_batch_sequential([all_tasks[0]])
                results.extend(warmup_results)

                # Run remaining tasks in parallel
                remaining_tasks = all_tasks[1:]
                if remaining_tasks:
                    print(f"\n{'=' * 60}")
                    print(
                        f"Running remaining {len(remaining_tasks)} simulations in parallel..."
                    )
                    print(f"{'=' * 60}")
                    parallel_results = run_batch_parallel(remaining_tasks, max_workers)
                    results.extend(parallel_results)
            else:
                # Sequential mode - just run all tasks in order
                results = run_batch_sequential(all_tasks)

            batch_elapsed = time.perf_counter() - batch_start_time

            # Print summary
            print_batch_summary(results, batch_elapsed, parallel=use_parallel)
            return

        if args.all_single_wing:
            # Filter to only single-wing configurations
            all_configs = [
                name for name in configs.list_configs() if is_single_wing_config(name)
            ]
            num_configs = len(all_configs)
            num_meshes = len(mesh_settings)
            total = num_configs * num_meshes

            # Build task list for all mesh types
            all_tasks = []
            for name in all_configs:
                for mesh_type in mesh_settings:
                    task = SimulationTask(
                        config_name=name,
                        mesh_type=mesh_type,
                        run_id=run_id,
                        show_results=False,  # Always False for batch runs
                    )
                    all_tasks.append(task)

            print(f"Run ID: {run_id}")
            if num_meshes > 1:
                mesh_label = " + ".join(mesh_settings)
                print(
                    f"Running {num_configs} single-wing configurations x {num_meshes} "
                    f"mesh types = {total} simulations ({mesh_label})"
                )
            else:
                print(
                    f"Running {num_configs} single-wing configurations "
                    f"({mesh_settings[0]} mesh)"
                )

            batch_start_time = time.perf_counter()
            results = []

            # Run simulations with warmup strategy
            use_parallel = not args.sequential
            if use_parallel:
                max_workers = args.workers if args.workers else _calculate_max_workers()

                # Run first task sequentially to warm up Numba cache
                print(f"\n{'=' * 60}")
                print("Running first simulation to warm up Numba cache...")
                print(f"{'=' * 60}")
                warmup_results = run_batch_sequential([all_tasks[0]])
                results.extend(warmup_results)

                # Run remaining tasks in parallel
                remaining_tasks = all_tasks[1:]
                if remaining_tasks:
                    print(f"\n{'=' * 60}")
                    print(
                        f"Running remaining {len(remaining_tasks)} simulations in parallel..."
                    )
                    print(f"{'=' * 60}")
                    parallel_results = run_batch_parallel(remaining_tasks, max_workers)
                    results.extend(parallel_results)
            else:
                # Sequential mode - just run all tasks in order
                results = run_batch_sequential(all_tasks)

            batch_elapsed = time.perf_counter() - batch_start_time

            # Print summary
            print_batch_summary(results, batch_elapsed, parallel=use_parallel)
            return

        for mesh_type in mesh_settings:
            sim_start_time = time.perf_counter()
            run_simulation(
                config_name=args.config,
                run_id=run_id,
                mesh_type=mesh_type,
                show_results=not args.no_show,
            )
            sim_elapsed = time.perf_counter() - sim_start_time
            print(f"Simulation completed in {format_elapsed_time(sim_elapsed)}")

    finally:
        # Restore stdout and close log file
        sys.stdout = original_stdout
        console_log_file.close()
        print(f"Console output saved to: {console_log_path}")


if __name__ == "__main__":
    main()
