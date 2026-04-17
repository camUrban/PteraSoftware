"""Benchmarks save/load times and file sizes for each serializable class.

Uses private serialization functions directly because the public save() and load()
restrict top level objects to public classes. The benchmark covers all classes in the
registry, including internal ones like LineVortex and Panel.
"""

import gzip
import json
import tempfile
import time
from pathlib import Path

import numpy as np

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware._panel import Panel

# noinspection PyProtectedMember
from pterasoftware._serialization import _object_from_dict, _object_to_dict

# noinspection PyProtectedMember
from pterasoftware._vortices._line_vortex import LineVortex

# noinspection PyProtectedMember
from pterasoftware._vortices.horseshoe_vortex import HorseshoeVortex

# noinspection PyProtectedMember
from pterasoftware._vortices.ring_vortex import RingVortex
from pterasoftware.geometry.airfoil import Airfoil
from pterasoftware.geometry.airplane import Airplane
from pterasoftware.geometry.wing import Wing
from pterasoftware.geometry.wing_cross_section import WingCrossSection
from pterasoftware.movements.airplane_movement import AirplaneMovement
from pterasoftware.movements.movement import Movement
from pterasoftware.movements.operating_point_movement import OperatingPointMovement
from pterasoftware.movements.wing_cross_section_movement import (
    WingCrossSectionMovement,
)
from pterasoftware.movements.wing_movement import WingMovement
from pterasoftware.operating_point import OperatingPoint
from pterasoftware.problems import SteadyProblem, UnsteadyProblem
from pterasoftware.steady_horseshoe_vortex_lattice_method import (
    SteadyHorseshoeVortexLatticeMethodSolver,
)
from pterasoftware.steady_ring_vortex_lattice_method import (
    SteadyRingVortexLatticeMethodSolver,
)
from pterasoftware.unsteady_ring_vortex_lattice_method import (
    UnsteadyRingVortexLatticeMethodSolver,
)


def _benchmark(name: str, obj: object) -> None:
    """Benchmarks save/load for a single object and prints the results.

    Uses private serialization functions directly to bypass the public API restriction
    that limits save/load to public classes.
    """
    serialized_dict = _object_to_dict(obj)

    with tempfile.TemporaryDirectory() as tmp:
        compact_path = Path(tmp) / "compact.json"
        indented_path = Path(tmp) / "indented.json"
        gz_path = Path(tmp) / "test.json.gz"

        # Benchmark compact JSON save.
        start = time.perf_counter()
        compact_bytes = json.dumps(serialized_dict).encode("utf-8")
        with open(compact_path, "wb") as f:
            f.write(compact_bytes)
        compact_save_time = time.perf_counter() - start

        compact_size = compact_path.stat().st_size

        # Benchmark compact JSON load.
        start = time.perf_counter()
        with open(compact_path, "rb") as f:
            _object_from_dict(json.loads(f.read()))
        compact_load_time = time.perf_counter() - start

        # Measure indented JSON size (used by save() for .json files).
        indented_bytes = json.dumps(serialized_dict, indent=2).encode("utf-8")
        with open(indented_path, "wb") as f:
            f.write(indented_bytes)
        indented_size = indented_path.stat().st_size

        # Benchmark gzip save (compact format).
        start = time.perf_counter()
        with gzip.open(gz_path, "wb") as f:
            f.write(compact_bytes)
        gz_save_time = time.perf_counter() - start

        gz_size = gz_path.stat().st_size

        # Benchmark gzip load.
        start = time.perf_counter()
        with gzip.open(gz_path, "rb") as f:
            _object_from_dict(json.loads(f.read()))
        gz_load_time = time.perf_counter() - start

    indent_overhead = indented_size / compact_size
    print(f"  {name}")
    print(
        f"    JSON (compact):  {compact_size:>12,} B | save {compact_save_time:.4f}s"
        f" | load {compact_load_time:.4f}s"
    )
    print(
        f"    JSON (indented): {indented_size:>12,} B | {indent_overhead:.2f}x overhead"
    )
    print(
        f"    GZIP:            {gz_size:>12,} B | save {gz_save_time:.4f}s"
        f" | load {gz_load_time:.4f}s"
    )
    print(f"    Compression ratio (compact / gzip): {compact_size / gz_size:.1f}x")
    print()


def _make_line_vortex() -> LineVortex:
    """Creates a typical LineVortex."""
    return LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0]),
        Elvp_GP1_CgP1=np.array([1.0, 0.5, 0.0]),
        strength=1.0,
    )


def _make_ring_vortex() -> RingVortex:
    """Creates a typical RingVortex."""
    ring_vortex = RingVortex(
        Frrvp_GP1_CgP1=np.array([1.0, 0.0, 0.0]),
        Flrvp_GP1_CgP1=np.array([1.0, 1.0, 0.0]),
        Blrvp_GP1_CgP1=np.array([0.0, 1.0, 0.0]),
        Brrvp_GP1_CgP1=np.array([0.0, 0.0, 0.0]),
        strength=1.0,
    )
    ring_vortex.age = 0.05
    return ring_vortex


def _make_horseshoe_vortex() -> HorseshoeVortex:
    """Creates a typical HorseshoeVortex."""
    return HorseshoeVortex(
        Frhvp_GP1_CgP1=np.array([1.0, 0.0, 0.0]),
        Flhvp_GP1_CgP1=np.array([1.0, 1.0, 0.0]),
        leftLegVector_GP1=np.array([-1.0, 0.0, 0.0]),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )


def _make_airfoil() -> Airfoil:
    """Creates a typical Airfoil."""
    return Airfoil(name="NACA2412")


def _make_operating_point() -> OperatingPoint:
    """Creates a typical OperatingPoint."""
    return OperatingPoint(rho=1.225, vCg__E=10.0, alpha=5.0, beta=0.0)


def _make_wing_cross_section() -> WingCrossSection:
    """Creates a typical WingCrossSection."""
    return WingCrossSection(
        airfoil=Airfoil(name="NACA0012"),
        num_spanwise_panels=8,
        chord=1.0,
    )


def _make_panel() -> Panel:
    """Creates a typical Panel."""
    return Panel(
        Frpp_G_Cg=np.array([1.0, 0.0, 0.0]),
        Flpp_G_Cg=np.array([1.0, 1.0, 0.0]),
        Blpp_G_Cg=np.array([0.0, 1.0, 0.0]),
        Brpp_G_Cg=np.array([0.0, 0.0, 0.0]),
        is_leading_edge=True,
        is_trailing_edge=False,
    )


def _make_meshed_wing() -> Wing:
    """Creates a typical meshed Wing with 8 x 16 = 128 Panels."""
    wing = Wing(
        wing_cross_sections=[
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=16,
                chord=1.0,
            ),
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=0.6,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
            ),
        ],
        num_chordwise_panels=8,
        chordwise_spacing="uniform",
    )
    airplane = Airplane(wings=[wing], name="Temp")
    operating_point = OperatingPoint()
    SteadyProblem(airplanes=[airplane], operating_point=operating_point)
    return airplane.wings[0]


def _make_meshed_airplane() -> Airplane:
    """Creates a typical meshed Airplane with two Wings, each 8 x 16 = 128 Panels (256
    total).
    """
    wing = Wing(
        wing_cross_sections=[
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=16,
                chord=1.0,
                control_surface_symmetry_type="symmetric",
            ),
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=0.6,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                control_surface_symmetry_type="symmetric",
            ),
        ],
        num_chordwise_panels=8,
        chordwise_spacing="uniform",
        symmetric=True,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
    )
    airplane = Airplane(wings=[wing], name="Benchmark")
    operating_point = OperatingPoint()
    SteadyProblem(airplanes=[airplane], operating_point=operating_point)
    return airplane


def _make_steady_problem() -> SteadyProblem:
    """Creates a typical SteadyProblem with a 256 panel Airplane."""
    wing = Wing(
        wing_cross_sections=[
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=16,
                chord=1.0,
                control_surface_symmetry_type="symmetric",
            ),
            WingCrossSection(
                airfoil=Airfoil(name="NACA0012"),
                num_spanwise_panels=None,
                chord=0.6,
                Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                control_surface_symmetry_type="symmetric",
            ),
        ],
        num_chordwise_panels=8,
        chordwise_spacing="uniform",
        symmetric=True,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
    )
    airplane = Airplane(wings=[wing], name="Benchmark")
    operating_point = OperatingPoint()
    return SteadyProblem(airplanes=[airplane], operating_point=operating_point)


def _make_operating_point_movement() -> OperatingPointMovement:
    """Creates a typical OperatingPointMovement."""
    return OperatingPointMovement(base_operating_point=OperatingPoint())


def _make_wing_cross_section_movement() -> WingCrossSectionMovement:
    """Creates a typical WingCrossSectionMovement."""
    return WingCrossSectionMovement(
        base_wing_cross_section=WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=8,
            chord=1.0,
        ),
    )


def _make_wing_movement() -> WingMovement:
    """Creates a typical WingMovement with a base Wing of 8 x 16 = 128 Panels."""
    wing_cross_section_movement_root = WingCrossSectionMovement(
        base_wing_cross_section=WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=16,
            chord=1.0,
        ),
    )
    wing_cross_section_movement_tip = WingCrossSectionMovement(
        base_wing_cross_section=WingCrossSection(
            airfoil=Airfoil(name="NACA0012"),
            num_spanwise_panels=None,
            chord=0.6,
            Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
        ),
    )
    return WingMovement(
        base_wing=Wing(
            wing_cross_sections=[
                wing_cross_section_movement_root.base_wing_cross_section,
                wing_cross_section_movement_tip.base_wing_cross_section,
            ],
            num_chordwise_panels=8,
            chordwise_spacing="uniform",
        ),
        wing_cross_section_movements=[
            wing_cross_section_movement_root,
            wing_cross_section_movement_tip,
        ],
    )


def _make_airplane_movement() -> AirplaneMovement:
    """Creates a typical AirplaneMovement with two Wings (256 panels total)."""
    wing_movement1 = _make_wing_movement()
    wing_movement2 = _make_wing_movement()
    return AirplaneMovement(
        base_airplane=Airplane(
            wings=[wing_movement1.base_wing, wing_movement2.base_wing]
        ),
        wing_movements=[wing_movement1, wing_movement2],
    )


def _make_movement() -> Movement:
    """Creates a typical Movement with 150 time steps and 256 panels."""
    airplane_movement = _make_airplane_movement()
    return Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=_make_operating_point_movement(),
        num_steps=150,
    )


def _make_unsteady_problem() -> UnsteadyProblem:
    """Creates a typical UnsteadyProblem with 150 time steps."""
    return UnsteadyProblem(movement=_make_movement())


def _make_unsteady_solver_problem() -> UnsteadyProblem:
    """Creates an UnsteadyProblem with two symmetric wings (256 panels total) for the
    unsteady solver benchmark. Follows the same pattern as the examples.
    """
    # Build the Airplane with two symmetric wings (main wing + v-tail).
    airplane = Airplane(
        wings=[
            Wing(
                wing_cross_sections=[
                    WingCrossSection(
                        airfoil=Airfoil(name="NACA0012"),
                        num_spanwise_panels=8,
                        chord=1.0,
                        control_surface_symmetry_type="symmetric",
                    ),
                    WingCrossSection(
                        airfoil=Airfoil(name="NACA0012"),
                        num_spanwise_panels=None,
                        chord=0.6,
                        Lp_Wcsp_Lpp=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                    ),
                ],
                name="Main Wing",
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
                symmetric=True,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            ),
            Wing(
                wing_cross_sections=[
                    WingCrossSection(
                        airfoil=Airfoil(name="NACA0012"),
                        num_spanwise_panels=8,
                        chord=0.8,
                        control_surface_symmetry_type="symmetric",
                    ),
                    WingCrossSection(
                        airfoil=Airfoil(name="NACA0012"),
                        num_spanwise_panels=None,
                        chord=0.5,
                        Lp_Wcsp_Lpp=(0.0, 2.0, 0.5),
                        control_surface_symmetry_type="symmetric",
                    ),
                ],
                name="V-Tail",
                Ler_Gs_Cgs=(4.0, 0.0, 0.0),
                num_chordwise_panels=8,
                chordwise_spacing="uniform",
                symmetric=True,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
            ),
        ],
    )

    # Build WingCrossSectionMovements referencing the post symmetry Airplane wings.
    main_wing_cross_section_movements = [
        WingCrossSectionMovement(base_wing_cross_section=wing_cross_section)
        for wing_cross_section in airplane.wings[0].wing_cross_sections
    ]
    tail_wing_cross_section_movements = [
        WingCrossSectionMovement(base_wing_cross_section=wing_cross_section)
        for wing_cross_section in airplane.wings[1].wing_cross_sections
    ]

    main_wing_movement = WingMovement(
        base_wing=airplane.wings[0],
        wing_cross_section_movements=main_wing_cross_section_movements,
    )
    tail_wing_movement = WingMovement(
        base_wing=airplane.wings[1],
        wing_cross_section_movements=tail_wing_cross_section_movements,
    )

    airplane_movement = AirplaneMovement(
        base_airplane=airplane,
        wing_movements=[main_wing_movement, tail_wing_movement],
    )
    movement = Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=OperatingPointMovement(
            base_operating_point=OperatingPoint(),
        ),
        num_steps=150,
    )
    return UnsteadyProblem(movement=movement)


if __name__ == "__main__":
    ps.set_up_logging()

    print("Serialization Benchmark")
    print("=" * 60)
    print()

    _benchmark("LineVortex", _make_line_vortex())
    _benchmark("RingVortex", _make_ring_vortex())
    _benchmark("HorseshoeVortex", _make_horseshoe_vortex())
    _benchmark("Airfoil", _make_airfoil())
    _benchmark("OperatingPoint", _make_operating_point())
    _benchmark("WingCrossSection", _make_wing_cross_section())
    _benchmark("Panel", _make_panel())
    _benchmark("Wing (meshed, 128 panels)", _make_meshed_wing())
    _benchmark("Airplane (meshed, 256 panels)", _make_meshed_airplane())
    _benchmark("SteadyProblem (256 panels)", _make_steady_problem())

    problem = _make_steady_problem()
    horseshoe_solver_unsolved = SteadyHorseshoeVortexLatticeMethodSolver(problem)
    _benchmark("Horseshoe solver (unsolved, 256 panels)", horseshoe_solver_unsolved)

    problem = _make_steady_problem()
    horseshoe_solver = SteadyHorseshoeVortexLatticeMethodSolver(problem)
    horseshoe_solver.run()
    _benchmark("Horseshoe solver (solved, 256 panels)", horseshoe_solver)

    problem = _make_steady_problem()
    ring_solver_unsolved = SteadyRingVortexLatticeMethodSolver(problem)
    _benchmark("Ring solver (unsolved, 256 panels)", ring_solver_unsolved)

    problem = _make_steady_problem()
    ring_solver = SteadyRingVortexLatticeMethodSolver(problem)
    ring_solver.run()
    _benchmark("Ring solver (solved, 256 panels)", ring_solver)

    _benchmark("OperatingPointMovement", _make_operating_point_movement())
    _benchmark("WingCrossSectionMovement", _make_wing_cross_section_movement())
    _benchmark("WingMovement (base: 128 panels)", _make_wing_movement())
    _benchmark("AirplaneMovement (base: 256 panels)", _make_airplane_movement())
    _benchmark("Movement (150 steps, 256 panels)", _make_movement())
    _benchmark("UnsteadyProblem (150 steps, 256 panels)", _make_unsteady_problem())

    # The unsteady solver benchmarks use a single wing (128 panels) because the
    # two identical wing geometry produces a singular influence matrix.
    unsolved_unsteady_solver = UnsteadyRingVortexLatticeMethodSolver(
        _make_unsteady_solver_problem()
    )
    _benchmark(
        "Unsteady solver (unsolved, 150 steps, 256 panels)", unsolved_unsteady_solver
    )

    solved_unsteady_solver = UnsteadyRingVortexLatticeMethodSolver(
        _make_unsteady_solver_problem()
    )
    solved_unsteady_solver.run()
    _benchmark(
        "Unsteady solver (solved, 150 steps, 256 panels)", solved_unsteady_solver
    )
