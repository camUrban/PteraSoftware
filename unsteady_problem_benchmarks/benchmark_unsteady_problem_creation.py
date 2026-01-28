"""Benchmark script for measuring UnsteadyProblem creation time.

This script measures how long it takes to create UnsteadyProblem objects with:
- Different numbers of panels (controlled by spanwise and chordwise panel counts)
- Static geometry (no flapping motion)
- Variable geometry (with flapping motion)
- Different numbers of cycles

Results are saved to JSON files in the benchmark_results/ directory for comparison
after optimizations.

Usage:
    cd experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe benchmark_unsteady_problem_creation.py

    # With a custom label for the run:
    cd experimental && PYTHONPATH="$PWD/.." ../.venv/Scripts/python.exe benchmark_unsteady_problem_creation.py --label "before_optimization"
"""

import argparse
import gc
import json
import platform
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import pterasoftware as ps


@dataclass
class BenchmarkConfig:
    """Configuration for a single benchmark run."""

    name: str
    num_spanwise_panels: int
    num_chordwise_panels: int
    num_cycles: int
    is_variable_geometry: bool
    num_cross_sections: int = 2

    @property
    def expected_panels_per_wing(self) -> int:
        """Expected number of panels for a single wing (before symmetry)."""
        num_sections = max(2, self.num_cross_sections) - 1
        return self.num_spanwise_panels * self.num_chordwise_panels * num_sections


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    config: BenchmarkConfig
    airplane_creation_time: float
    movement_creation_time: float
    problem_creation_time: float
    total_time: float
    actual_num_panels: int
    num_steps: int

    def __str__(self) -> str:
        return (
            f"{self.config.name:40s} | "
            f"panels={self.actual_num_panels:5d} | "
            f"steps={self.num_steps:4d} | "
            f"airplane={self.airplane_creation_time:6.3f}s | "
            f"movement={self.movement_creation_time:6.3f}s | "
            f"problem={self.problem_creation_time:6.3f}s | "
            f"total={self.total_time:7.3f}s"
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "config": asdict(self.config),
            "airplane_creation_time": self.airplane_creation_time,
            "movement_creation_time": self.movement_creation_time,
            "problem_creation_time": self.problem_creation_time,
            "total_time": self.total_time,
            "actual_num_panels": self.actual_num_panels,
            "num_steps": self.num_steps,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "BenchmarkResult":
        """Create from dictionary."""
        return cls(
            config=BenchmarkConfig(**data["config"]),
            airplane_creation_time=data["airplane_creation_time"],
            movement_creation_time=data["movement_creation_time"],
            problem_creation_time=data["problem_creation_time"],
            total_time=data["total_time"],
            actual_num_panels=data["actual_num_panels"],
            num_steps=data["num_steps"],
        )


@dataclass
class BenchmarkRun:
    """Complete benchmark run with metadata."""

    timestamp: str
    label: str
    python_version: str
    platform_info: str
    git_commit: str
    results: list[BenchmarkResult]

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "timestamp": self.timestamp,
            "label": self.label,
            "python_version": self.python_version,
            "platform_info": self.platform_info,
            "git_commit": self.git_commit,
            "results": [r.to_dict() for r in self.results],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "BenchmarkRun":
        """Create from dictionary."""
        return cls(
            timestamp=data["timestamp"],
            label=data["label"],
            python_version=data["python_version"],
            platform_info=data["platform_info"],
            git_commit=data["git_commit"],
            results=[BenchmarkResult.from_dict(r) for r in data["results"]],
        )

    def save(self, output_dir: Path) -> Path:
        """Save benchmark run to JSON file.

        :param output_dir: Directory to save the file.
        :return: Path to the saved file.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create filename from timestamp and label
        safe_label = self.label.replace(" ", "_").replace("/", "-")
        filename = f"benchmark_{self.timestamp}_{safe_label}.json"
        filepath = output_dir / filename

        with open(filepath, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

        return filepath

    @classmethod
    def load(cls, filepath: Path) -> "BenchmarkRun":
        """Load benchmark run from JSON file.

        :param filepath: Path to the JSON file.
        :return: BenchmarkRun object.
        """
        with open(filepath) as f:
            data = json.load(f)
        return cls.from_dict(data)


def get_git_commit() -> str:
    """Get the current git commit hash."""
    import subprocess

    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return "unknown"


def get_system_info() -> str:
    """Get system information string."""
    return f"{platform.system()} {platform.release()} ({platform.machine()})"


def create_simple_airplane(
    num_spanwise_panels: int,
    num_chordwise_panels: int,
    use_symmetric_wing: bool = True,
    num_cross_sections: int = 2,
) -> ps.geometry.airplane.Airplane:
    """Create a simple airplane with a single wing for benchmarking.

    :param num_spanwise_panels: Number of spanwise panels per wing section.
    :param num_chordwise_panels: Number of chordwise panels.
    :param use_symmetric_wing: If True, creates a symmetric wing (type 4 symmetry).
    :param num_cross_sections: Number of wing cross sections (minimum 2).
    :return: An Airplane object.
    """
    if use_symmetric_wing:
        # Type 4 symmetry: symmetric wing with coincident symmetry plane
        symmetric = True
        symmetry_normal = (0.0, 1.0, 0.0)
        symmetry_point = (0.0, 0.0, 0.0)
        control_surface_symmetry_type: str | None = "symmetric"
    else:
        symmetric = False
        symmetry_normal = None
        symmetry_point = None
        control_surface_symmetry_type = None

    # Create wing cross sections
    num_cross_sections = max(2, num_cross_sections)
    wing_cross_sections = []
    span_per_section = 5.0 / (num_cross_sections - 1)  # Total span of 5.0

    for i in range(num_cross_sections):
        is_last = i == num_cross_sections - 1
        # Taper from chord 1.0 at root to 0.6 at tip
        taper_ratio = i / (num_cross_sections - 1)
        chord = 1.0 - 0.4 * taper_ratio
        # Sweep: x position increases toward tip
        x_pos = 0.2 * taper_ratio
        y_pos = span_per_section * i

        wcs = ps.geometry.wing_cross_section.WingCrossSection(
            num_spanwise_panels=None if is_last else num_spanwise_panels,
            chord=chord,
            Lp_Wcsp_Lpp=(x_pos, y_pos, 0.0),
            angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            control_surface_symmetry_type=control_surface_symmetry_type,
            control_surface_hinge_point=0.75,
            control_surface_deflection=0.0,
            spanwise_spacing=None if is_last else "uniform",
            airfoil=ps.geometry.airfoil.Airfoil(name="naca2412"),
        )
        wing_cross_sections.append(wcs)

    airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=wing_cross_sections,
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=symmetric,
                mirror_only=False,
                symmetryNormal_G=symmetry_normal,
                symmetryPoint_G_Cg=symmetry_point,
                num_chordwise_panels=num_chordwise_panels,
                chordwise_spacing="uniform",
            ),
        ],
        name="Benchmark Airplane",
    )

    return airplane


def create_movement(
    airplane: ps.geometry.airplane.Airplane,
    num_cycles: int,
    is_variable_geometry: bool,
) -> ps.movements.movement.Movement:
    """Create a Movement for the airplane.

    :param airplane: The base Airplane object.
    :param num_cycles: Number of flapping cycles.
    :param is_variable_geometry: If True, adds flapping motion to the wing.
    :return: A Movement object.
    """
    # Create WingCrossSectionMovements for each cross section
    wing_cross_section_movements = []
    for wcs in airplane.wings[0].wing_cross_sections:
        wcs_movement = (
            ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
                base_wing_cross_section=wcs,
                ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
                phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
                phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
            )
        )
        wing_cross_section_movements.append(wcs_movement)

    # Set flapping parameters based on geometry type
    if is_variable_geometry:
        # Variable geometry: wing flaps with 15 degree amplitude, 1 second period
        amp_angles = (15.0, 0.0, 0.0)
        period_angles = (1.0, 0.0, 0.0)
    else:
        # Static geometry: no flapping
        amp_angles = (0.0, 0.0, 0.0)
        period_angles = (0.0, 0.0, 0.0)

    # Create WingMovement
    wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=wing_cross_section_movements,
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=amp_angles,
        periodAngles_Gs_to_Wn_ixyz=period_angles,
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    # Create AirplaneMovement
    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[wing_movement],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    # Create OperatingPoint and OperatingPointMovement
    operating_point = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=operating_point,
            ampVCg__E=0.0,
            periodVCg__E=0.0,
            spacingVCg__E="sine",
            phaseVCg__E=0.0,
        )
    )

    # Create Movement
    if is_variable_geometry:
        # Variable geometry: use num_cycles
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=None,
            num_cycles=num_cycles,
            num_chords=None,
            num_steps=None,
        )
    else:
        # Static geometry: use num_chords (simulates wake development)
        # For static geometry, num_steps = num_chords * num_chordwise_panels
        num_chords = num_cycles * 10  # Scale to get similar step counts
        movement = ps.movements.movement.Movement(
            airplane_movements=[airplane_movement],
            operating_point_movement=operating_point_movement,
            delta_time=None,
            num_cycles=None,
            num_chords=num_chords,
            num_steps=None,
        )

    return movement


def run_benchmark(config: BenchmarkConfig) -> BenchmarkResult:
    """Run a single benchmark with the given configuration.

    :param config: The benchmark configuration.
    :return: The benchmark result.
    """
    # Force garbage collection before starting
    gc.collect()

    # Time airplane creation
    # For variable geometry, we must use non-symmetric wings because flapping
    # motion changes the symmetry type (plane becomes non-coincident with wing axes).
    use_symmetric = not config.is_variable_geometry

    start_time = time.perf_counter()
    airplane = create_simple_airplane(
        num_spanwise_panels=config.num_spanwise_panels,
        num_chordwise_panels=config.num_chordwise_panels,
        use_symmetric_wing=use_symmetric,
        num_cross_sections=config.num_cross_sections,
    )
    airplane_time = time.perf_counter() - start_time

    # Time movement creation
    start_time = time.perf_counter()
    movement = create_movement(
        airplane=airplane,
        num_cycles=config.num_cycles,
        is_variable_geometry=config.is_variable_geometry,
    )
    movement_time = time.perf_counter() - start_time

    # Time problem creation
    start_time = time.perf_counter()
    problem = ps.problems.UnsteadyProblem(movement=movement)
    problem_time = time.perf_counter() - start_time

    total_time = airplane_time + movement_time + problem_time

    result = BenchmarkResult(
        config=config,
        airplane_creation_time=airplane_time,
        movement_creation_time=movement_time,
        problem_creation_time=problem_time,
        total_time=total_time,
        actual_num_panels=airplane.num_panels,
        num_steps=problem.num_steps,
    )

    # Clean up to free memory
    del problem
    del movement
    del airplane
    gc.collect()

    return result


def compare_runs(run1: BenchmarkRun, run2: BenchmarkRun) -> None:
    """Compare two benchmark runs and print the differences.

    :param run1: The baseline (older) run.
    :param run2: The comparison (newer) run.
    """
    print()
    print("=" * 120)
    print("Benchmark Comparison")
    print("=" * 120)
    print(f"Baseline: {run1.label} ({run1.timestamp}, commit {run1.git_commit})")
    print(f"Current:  {run2.label} ({run2.timestamp}, commit {run2.git_commit})")
    print()

    # Create lookup by config name
    run1_by_name = {r.config.name: r for r in run1.results}
    run2_by_name = {r.config.name: r for r in run2.results}

    # Find common configs
    common_names = set(run1_by_name.keys()) & set(run2_by_name.keys())

    if not common_names:
        print("No common configurations found between the two runs.")
        return

    print(
        f"{'Configuration':40s} | {'Baseline':>9s} | {'Current':>9s} | "
        f"{'Diff':>9s} | {'Change':>8s}"
    )
    print("-" * 120)

    total_baseline = 0.0
    total_current = 0.0

    for name in sorted(common_names):
        r1 = run1_by_name[name]
        r2 = run2_by_name[name]

        baseline_time = r1.total_time
        current_time = r2.total_time
        diff = current_time - baseline_time

        if baseline_time > 0:
            pct_change = (diff / baseline_time) * 100
            pct_str = f"{pct_change:+.1f}%"
        else:
            pct_str = "N/A"

        # Color indicators (using text since we're in terminal)
        if diff < -0.01:
            indicator = "faster"
        elif diff > 0.01:
            indicator = "SLOWER"
        else:
            indicator = "same"

        print(
            f"{name:40s} | {baseline_time:8.3f}s | {current_time:8.3f}s | "
            f"{diff:+8.3f}s | {pct_str:>7s} {indicator}"
        )

        total_baseline += baseline_time
        total_current += current_time

    print("-" * 120)
    total_diff = total_current - total_baseline
    if total_baseline > 0:
        total_pct = (total_diff / total_baseline) * 100
        total_pct_str = f"{total_pct:+.1f}%"
    else:
        total_pct_str = "N/A"

    print(
        f"{'TOTAL':40s} | {total_baseline:8.3f}s | {total_current:8.3f}s | "
        f"{total_diff:+8.3f}s | {total_pct_str:>7s}"
    )


def list_saved_runs(output_dir: Path) -> list[Path]:
    """List all saved benchmark runs.

    :param output_dir: Directory containing benchmark files.
    :return: List of paths to benchmark files, sorted by timestamp.
    """
    if not output_dir.exists():
        return []

    files = list(output_dir.glob("benchmark_*.json"))
    return sorted(files)


def main() -> None:
    """Run the benchmark suite."""
    parser = argparse.ArgumentParser(
        description="Benchmark UnsteadyProblem creation time."
    )
    parser.add_argument(
        "--label",
        type=str,
        default="default",
        help="Label for this benchmark run (e.g., 'before_optimization')",
    )
    parser.add_argument(
        "--compare",
        type=str,
        metavar="FILE",
        help="Compare results against a previous benchmark file",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List all saved benchmark runs",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "benchmark_results",
        help="Directory to save/load benchmark results",
    )

    args = parser.parse_args()

    # Handle --list command
    if args.list:
        runs = list_saved_runs(args.output_dir)
        if not runs:
            print(f"No benchmark runs found in {args.output_dir}")
        else:
            print(f"Saved benchmark runs in {args.output_dir}:")
            for filepath in runs:
                run = BenchmarkRun.load(filepath)
                print(f"  {filepath.name}")
                print(f"    Label: {run.label}")
                print(f"    Timestamp: {run.timestamp}")
                print(f"    Commit: {run.git_commit}")
                print(f"    Results: {len(run.results)} configurations")
                print()
        return

    print("=" * 120)
    print("UnsteadyProblem Creation Benchmark")
    print("=" * 120)
    print()
    print(f"Label: {args.label}")
    print(f"Python: {sys.version}")
    print(f"Platform: {get_system_info()}")
    print(f"Git commit: {get_git_commit()}")
    print()

    # Define benchmark configurations
    configs: list[BenchmarkConfig] = []

    # Panel count scaling (static geometry, 2 cycles equivalent)
    panel_counts = [
        (4, 4),  # 16 panels per half-wing
        (8, 6),  # 48 panels per half-wing
        (12, 8),  # 96 panels per half-wing
        (16, 10),  # 160 panels per half-wing
        (20, 12),  # 240 panels per half-wing
        (24, 14),  # 336 panels per half-wing
    ]

    for spanwise, chordwise in panel_counts:
        configs.append(
            BenchmarkConfig(
                name=f"Static {spanwise}x{chordwise} panels, 2 cycles eq",
                num_spanwise_panels=spanwise,
                num_chordwise_panels=chordwise,
                num_cycles=2,
                is_variable_geometry=False,
            )
        )

    # Cycle count scaling (static geometry, medium panel count)
    for num_cycles in [1, 2, 3, 5]:
        configs.append(
            BenchmarkConfig(
                name=f"Static 12x8 panels, {num_cycles} cycles eq",
                num_spanwise_panels=12,
                num_chordwise_panels=8,
                num_cycles=num_cycles,
                is_variable_geometry=False,
            )
        )

    # Variable geometry: panel count scaling
    for spanwise, chordwise in panel_counts:
        configs.append(
            BenchmarkConfig(
                name=f"Variable {spanwise}x{chordwise} panels, 2 cycles",
                num_spanwise_panels=spanwise,
                num_chordwise_panels=chordwise,
                num_cycles=2,
                is_variable_geometry=True,
            )
        )

    # Variable geometry: cycle count scaling
    for num_cycles in [1, 2, 3, 5]:
        configs.append(
            BenchmarkConfig(
                name=f"Variable 12x8 panels, {num_cycles} cycles",
                num_spanwise_panels=12,
                num_chordwise_panels=8,
                num_cycles=num_cycles,
                is_variable_geometry=True,
            )
        )

    # Many cross sections: test scaling with 10 wing cross sections
    # This creates 9 spanwise sections (between 10 cross sections)
    configs.append(
        BenchmarkConfig(
            name="Static 4x6 panels, 10 WCS, 5 cycles eq",
            num_spanwise_panels=4,
            num_chordwise_panels=6,
            num_cycles=5,
            is_variable_geometry=False,
            num_cross_sections=10,
        )
    )
    configs.append(
        BenchmarkConfig(
            name="Static 4x6 panels, 10 WCS, 10 cycles eq",
            num_spanwise_panels=4,
            num_chordwise_panels=6,
            num_cycles=10,
            is_variable_geometry=False,
            num_cross_sections=10,
        )
    )
    configs.append(
        BenchmarkConfig(
            name="Variable 4x6 panels, 10 WCS, 5 cycles",
            num_spanwise_panels=4,
            num_chordwise_panels=6,
            num_cycles=5,
            is_variable_geometry=True,
            num_cross_sections=10,
        )
    )
    configs.append(
        BenchmarkConfig(
            name="Variable 4x6 panels, 10 WCS, 10 cycles",
            num_spanwise_panels=4,
            num_chordwise_panels=6,
            num_cycles=10,
            is_variable_geometry=True,
            num_cross_sections=10,
        )
    )

    # Run benchmarks
    results: list[BenchmarkResult] = []

    print(
        f"{'Configuration':40s} | {'panels':>6s} | {'steps':>5s} | "
        f"{'airplane':>9s} | {'movement':>9s} | {'problem':>9s} | {'total':>8s}"
    )
    print("-" * 120)

    for config in configs:
        result = run_benchmark(config)
        results.append(result)
        print(result)

    print()
    print("=" * 120)
    print("Summary")
    print("=" * 120)

    # Group results by geometry type
    static_results = [r for r in results if not r.config.is_variable_geometry]
    variable_results = [r for r in results if r.config.is_variable_geometry]

    print()
    print("Static Geometry Results:")
    print("-" * 60)
    if static_results:
        avg_airplane = sum(r.airplane_creation_time for r in static_results) / len(
            static_results
        )
        avg_movement = sum(r.movement_creation_time for r in static_results) / len(
            static_results
        )
        avg_problem = sum(r.problem_creation_time for r in static_results) / len(
            static_results
        )
        print(f"  Average airplane creation time: {avg_airplane:.3f}s")
        print(f"  Average movement creation time: {avg_movement:.3f}s")
        print(f"  Average problem creation time:  {avg_problem:.3f}s")

    print()
    print("Variable Geometry Results:")
    print("-" * 60)
    if variable_results:
        avg_airplane = sum(r.airplane_creation_time for r in variable_results) / len(
            variable_results
        )
        avg_movement = sum(r.movement_creation_time for r in variable_results) / len(
            variable_results
        )
        avg_problem = sum(r.problem_creation_time for r in variable_results) / len(
            variable_results
        )
        print(f"  Average airplane creation time: {avg_airplane:.3f}s")
        print(f"  Average movement creation time: {avg_movement:.3f}s")
        print(f"  Average problem creation time:  {avg_problem:.3f}s")

    # Find slowest operations
    print()
    print("Slowest Problem Creations:")
    print("-" * 60)
    sorted_results = sorted(
        results, key=lambda r: r.problem_creation_time, reverse=True
    )
    for result in sorted_results[:5]:
        print(
            f"  {result.config.name}: {result.problem_creation_time:.3f}s "
            f"({result.actual_num_panels} panels, {result.num_steps} steps)"
        )

    # Create and save the benchmark run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    benchmark_run = BenchmarkRun(
        timestamp=timestamp,
        label=args.label,
        python_version=sys.version,
        platform_info=get_system_info(),
        git_commit=get_git_commit(),
        results=results,
    )

    saved_path = benchmark_run.save(args.output_dir)
    print()
    print("=" * 120)
    print(f"Results saved to: {saved_path}")
    print("=" * 120)

    # Handle --compare option
    if args.compare:
        compare_path = Path(args.compare)
        if not compare_path.exists():
            # Try relative to output_dir
            compare_path = args.output_dir / args.compare

        if compare_path.exists():
            baseline_run = BenchmarkRun.load(compare_path)
            compare_runs(baseline_run, benchmark_run)
        else:
            print(f"Warning: Comparison file not found: {args.compare}")


if __name__ == "__main__":
    main()
