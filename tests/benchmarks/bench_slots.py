"""Contains benchmarks for measuring the impact of adding __slots__ to all classes.

Run this script before adding __slots__ (to record baseline results) and again after all
__slots__ are added (to measure the improvement). Results are printed to stdout and can
be redirected to a file for comparison.

Usage:
    PYTHONPATH="$PWD" .venv/bin/python tests/benchmarks/bench_slots.py
"""

from __future__ import annotations

import sys
import timeit
import tracemalloc

import numpy as np

import pterasoftware as ps

# noinspection PyProtectedMember
from pterasoftware import _panel, _vortices

# noinspection PyProtectedMember
from pterasoftware._vortices import _line_vortex

# ---------------------------------------------------------------------------
# Helper functions for creating instances.
# ---------------------------------------------------------------------------


def _make_line_vortex() -> _line_vortex.LineVortex:
    """Creates a basic LineVortex for benchmarking.

    :return: A LineVortex with unit strength.
    """
    return _line_vortex.LineVortex(
        Slvp_GP1_CgP1=np.array([0.0, 0.0, 0.0], dtype=float),
        Elvp_GP1_CgP1=np.array([1.0, 0.0, 0.0], dtype=float),
        strength=1.0,
    )


def _make_ring_vortex() -> _vortices.ring_vortex.RingVortex:
    """Creates a basic RingVortex for benchmarking.

    :return: A RingVortex with unit strength.
    """
    return _vortices.ring_vortex.RingVortex(
        Frrvp_GP1_CgP1=np.array([0.0, 0.5, 0.0], dtype=float),
        Flrvp_GP1_CgP1=np.array([0.0, -0.5, 0.0], dtype=float),
        Blrvp_GP1_CgP1=np.array([1.0, -0.5, 0.0], dtype=float),
        Brrvp_GP1_CgP1=np.array([1.0, 0.5, 0.0], dtype=float),
        strength=1.0,
    )


def _make_horseshoe_vortex() -> _vortices.horseshoe_vortex.HorseshoeVortex:
    """Creates a basic HorseshoeVortex for benchmarking.

    :return: A HorseshoeVortex with unit strength.
    """
    return _vortices.horseshoe_vortex.HorseshoeVortex(
        Frhvp_GP1_CgP1=np.array([0.0, 0.5, 0.0], dtype=float),
        Flhvp_GP1_CgP1=np.array([0.0, -0.5, 0.0], dtype=float),
        leftLegVector_GP1=np.array([1.0, 0.0, 0.0], dtype=float),
        left_right_leg_lengths=20.0,
        strength=1.0,
    )


def _make_panel() -> _panel.Panel:
    """Creates a basic Panel for benchmarking.

    :return: A Panel with 1.0 m chord and 0.5 m span.
    """
    return _panel.Panel(
        Frpp_G_Cg=np.array([0.0, 0.5, 0.0], dtype=float),
        Flpp_G_Cg=np.array([0.0, 0.0, 0.0], dtype=float),
        Blpp_G_Cg=np.array([1.0, 0.0, 0.0], dtype=float),
        Brpp_G_Cg=np.array([1.0, 0.5, 0.0], dtype=float),
        is_leading_edge=False,
        is_trailing_edge=False,
    )


def _make_airfoil() -> ps.geometry.airfoil.Airfoil:
    """Creates a basic Airfoil for benchmarking.

    :return: An Airfoil using the NACA 2412 profile.
    """
    return ps.geometry.airfoil.Airfoil(
        name="naca2412",
        outline_A_lp=None,
        resample=True,
        n_points_per_side=50,
    )


def _make_wing_cross_section() -> ps.geometry.wing_cross_section.WingCrossSection:
    """Creates a basic root WingCrossSection for benchmarking.

    :return: A root WingCrossSection with 10 spanwise panels.
    """
    return ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=_make_airfoil(),
        num_spanwise_panels=10,
        chord=2.0,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )


def _make_wing() -> ps.geometry.wing.Wing:
    """Creates a basic symmetric Wing for benchmarking.

    :return: A symmetric Wing with 8 chordwise panels.
    """
    return ps.geometry.wing.Wing(
        wing_cross_sections=[
            ps.geometry.wing_cross_section.WingCrossSection(
                airfoil=_make_airfoil(),
                num_spanwise_panels=8,
                chord=1.5,
                Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                control_surface_symmetry_type="symmetric",
                control_surface_hinge_point=0.75,
                control_surface_deflection=0.0,
                spanwise_spacing="cosine",
            ),
            ps.geometry.wing_cross_section.WingCrossSection(
                airfoil=_make_airfoil(),
                num_spanwise_panels=None,
                chord=1.0,
                Lp_Wcsp_Lpp=(0.5, 5.0, 0.0),
                angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                control_surface_symmetry_type="symmetric",
                control_surface_hinge_point=0.75,
                control_surface_deflection=0.0,
                spanwise_spacing=None,
            ),
        ],
        name="Benchmark Wing",
        Ler_Gs_Cgs=(0.0, 0.0, 0.0),
        angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
        symmetric=True,
        mirror_only=False,
        symmetryNormal_G=(0.0, 1.0, 0.0),
        symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
        num_chordwise_panels=8,
        chordwise_spacing="cosine",
    )


def _make_airplane() -> ps.geometry.airplane.Airplane:
    """Creates a basic single Wing Airplane for benchmarking.

    :return: An Airplane with one symmetric Wing.
    """
    return ps.geometry.airplane.Airplane(
        wings=[_make_wing()],
        name="Benchmark Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
    )


def _make_operating_point() -> ps.operating_point.OperatingPoint:
    """Creates a basic OperatingPoint for benchmarking.

    :return: An OperatingPoint with standard sea level conditions.
    """
    return ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=5.0,
        beta=0.0,
        externalFX_W=0.0,
        nu=15.06e-6,
    )


def _make_steady_problem() -> ps.problems.SteadyProblem:
    """Creates a basic SteadyProblem for benchmarking.

    :return: A SteadyProblem with one Airplane.
    """
    return ps.problems.SteadyProblem(
        airplanes=[_make_airplane()],
        operating_point=_make_operating_point(),
    )


def _make_operating_point_movement() -> (
    ps.movements.operating_point_movement.OperatingPointMovement
):
    """Creates a basic OperatingPointMovement for benchmarking.

    :return: An OperatingPointMovement with no movement.
    """
    return ps.movements.operating_point_movement.OperatingPointMovement(
        base_operating_point=_make_operating_point(),
        ampVCg__E=0.0,
        periodVCg__E=0.0,
        spacingVCg__E="sine",
        phaseVCg__E=0.0,
    )


def _make_wing_cross_section_movement() -> (
    ps.movements.wing_cross_section_movement.WingCrossSectionMovement
):
    """Creates a basic WingCrossSectionMovement for benchmarking.

    :return: A WingCrossSectionMovement with no movement.
    """
    return ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
        base_wing_cross_section=_make_wing_cross_section(),
        ampLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        periodLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        spacingLp_Wcsp_Lpp=("sine", "sine", "sine"),
        phaseLp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        ampAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        periodAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        spacingAngles_Wcsp_to_Wcs_ixyz=("sine", "sine", "sine"),
        phaseAngles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
    )


def _make_wing_movement() -> ps.movements.wing_movement.WingMovement:
    """Creates a basic WingMovement for benchmarking.

    :return: A WingMovement with no movement.
    """
    root_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=_make_airfoil(),
        num_spanwise_panels=8,
        chord=1.5,
        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing="cosine",
    )
    tip_wing_cross_section = ps.geometry.wing_cross_section.WingCrossSection(
        airfoil=_make_airfoil(),
        num_spanwise_panels=None,
        chord=1.0,
        Lp_Wcsp_Lpp=(0.5, 5.0, 0.0),
        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
        control_surface_symmetry_type="symmetric",
        control_surface_hinge_point=0.75,
        control_surface_deflection=0.0,
        spanwise_spacing=None,
    )
    root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=root_wing_cross_section,
        )
    )
    tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=tip_wing_cross_section,
        )
    )
    return ps.movements.wing_movement.WingMovement(
        base_wing=_make_wing(),
        wing_cross_section_movements=[
            root_wing_cross_section_movement,
            tip_wing_cross_section_movement,
        ],
    )


def _make_airplane_movement() -> ps.movements.airplane_movement.AirplaneMovement:
    """Creates a basic AirplaneMovement for benchmarking.

    :return: An AirplaneMovement with no movement.
    """
    return ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=_make_airplane(),
        wing_movements=[_make_wing_movement()],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )


def _make_movement() -> ps.movements.movement.Movement:
    """Creates a basic Movement for benchmarking.

    :return: A Movement with 3 chords and no movement.
    """
    return ps.movements.movement.Movement(
        airplane_movements=[_make_airplane_movement()],
        operating_point_movement=_make_operating_point_movement(),
        num_chords=3,
    )


def _make_unsteady_problem() -> ps.problems.UnsteadyProblem:
    """Creates a basic UnsteadyProblem for benchmarking.

    :return: An UnsteadyProblem with all results retained.
    """
    return ps.problems.UnsteadyProblem(
        movement=_make_movement(),
        only_final_results=False,
    )


def _make_multi_wing_airplane() -> ps.geometry.airplane.Airplane:
    """Creates a multi Wing Airplane for aggregate memory benchmarking.

    Uses a two Wing configuration (main Wing + horizontal stabilizer) with symmetric
    Wings for a realistic panel count.

    :return: An Airplane with two symmetric Wings and ~750 panels.
    """
    return ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=20,
                        chord=1.0,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(1.0, 5.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=14,
                chordwise_spacing="cosine",
            ),
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=12,
                        chord=0.75,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca0010",
                            outline_A_lp=None,
                            resample=True,
                            n_points_per_side=50,
                        ),
                        num_spanwise_panels=None,
                        chord=0.50,
                        Lp_Wcsp_Lpp=(0.25, 2.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                    ),
                ],
                name="Horizontal Stabilizer",
                Ler_Gs_Cgs=(5.0, 0.0, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, -5.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=8,
                chordwise_spacing="cosine",
            ),
        ],
        name="Benchmark Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )


# ---------------------------------------------------------------------------
# 1. Per instance memory (sys.getsizeof).
# ---------------------------------------------------------------------------


def bench_per_instance_memory() -> None:
    """Measures per instance memory using sys.getsizeof for all 16 classes.

    :return: None
    """
    print("=" * 70)
    print("PER INSTANCE MEMORY (sys.getsizeof)")
    print("=" * 70)

    # Build instances in dependency order. The higher level classes are more
    # expensive to construct, but we only need one instance of each.
    instances = {
        "LineVortex": _make_line_vortex(),
        "RingVortex": _make_ring_vortex(),
        "HorseshoeVortex": _make_horseshoe_vortex(),
        "Panel": _make_panel(),
        "Airfoil": _make_airfoil(),
        "WingCrossSection": _make_wing_cross_section(),
        "Wing": _make_wing(),
        "Airplane": _make_airplane(),
        "OperatingPoint": _make_operating_point(),
        "SteadyProblem": _make_steady_problem(),
        "UnsteadyProblem": _make_unsteady_problem(),
        "OperatingPointMovement": _make_operating_point_movement(),
        "WingCrossSectionMovement": _make_wing_cross_section_movement(),
        "WingMovement": _make_wing_movement(),
        "AirplaneMovement": _make_airplane_movement(),
        "Movement": _make_movement(),
    }

    print(f"{'Class':<28} {'Size (bytes)':>15} {'Dict (bytes)':>15}")
    print("-" * 58)
    for name, instance in instances.items():
        instance_size = sys.getsizeof(instance)
        has_dict = hasattr(instance, "__dict__")

        if has_dict:
            dict_size_str = f"{sys.getsizeof(instance.__dict__):>15,}"
        else:
            dict_size_str = f"{'N/A':>15}"

        print(f"{name:<28} {instance_size:>15,} {dict_size_str}")

    print()


# ---------------------------------------------------------------------------
# 2. Aggregate memory (tracemalloc) for a realistic workload.
# ---------------------------------------------------------------------------


def bench_aggregate_memory() -> None:
    """Measures total memory for constructing a multi Wing Airplane using tracemalloc.

    A warm up Airplane is constructed and discarded before measuring so that one time
    costs (Numba JIT compilation, Airfoil file I/O, scipy initialization) are not
    included in the measurement.

    :return: None
    """
    print("=" * 70)
    print("AGGREGATE MEMORY (tracemalloc)")
    print("=" * 70)

    # Warm up: trigger all one time initialization (Numba JIT, file I/O, etc.).
    _make_multi_wing_airplane()

    tracemalloc.start()
    snapshot_before = tracemalloc.take_snapshot()

    airplane = _make_multi_wing_airplane()

    snapshot_after = tracemalloc.take_snapshot()
    tracemalloc.stop()

    stats = snapshot_after.compare_to(snapshot_before, "lineno")

    total_allocated = sum(stat.size_diff for stat in stats if stat.size_diff > 0)

    total_panels = airplane.num_panels
    print(f"Airplane panels:          {total_panels}")
    print(f"Total memory allocated:   {total_allocated:,} bytes")
    print(f"                          {total_allocated / 1024:.1f} KB")
    if total_panels > 0:
        print(f"Bytes per panel:          {total_allocated / total_panels:.0f}")

    print()


# ---------------------------------------------------------------------------
# 3. Attribute access speed (timeit) for hot path classes.
# ---------------------------------------------------------------------------


def bench_attribute_access() -> None:
    """Measures attribute read and write speed for hot path classes.

    :return: None
    """
    print("=" * 70)
    print("ATTRIBUTE ACCESS SPEED (timeit)")
    print("=" * 70)

    num_iterations = 100_000

    # Panel: read Frpp_G_Cg (immutable cached property), read area (computed cache).
    panel = _make_panel()
    panel_read_time = timeit.timeit(
        "panel.Frpp_G_Cg",
        globals={"panel": panel},
        number=num_iterations,
    )

    # Access area once to populate cache, then benchmark cached reads.
    _ = panel.area
    panel_cache_read_time = timeit.timeit(
        "panel.area",
        globals={"panel": panel},
        number=num_iterations,
    )

    # RingVortex: read strength (property with getter), write age.
    ring_vortex = _make_ring_vortex()
    rv_read_time = timeit.timeit(
        "ring_vortex.strength",
        globals={"ring_vortex": ring_vortex},
        number=num_iterations,
    )
    rv_write_time = timeit.timeit(
        "ring_vortex.age = 5",
        globals={"ring_vortex": ring_vortex},
        number=num_iterations,
    )

    # HorseshoeVortex: read strength.
    horseshoe = _make_horseshoe_vortex()
    hv_read_time = timeit.timeit(
        "horseshoe.strength",
        globals={"horseshoe": horseshoe},
        number=num_iterations,
    )

    print(f"{'Benchmark':<45} {'Time (ms)':>12} {'Per op (ns)':>12}")
    print("-" * 69)
    print(
        f"{'Panel.Frpp_G_Cg read':<45} "
        f"{panel_read_time * 1000:>12.2f} "
        f"{panel_read_time / num_iterations * 1e9:>12.1f}"
    )
    print(
        f"{'Panel.area cached read':<45} "
        f"{panel_cache_read_time * 1000:>12.2f} "
        f"{panel_cache_read_time / num_iterations * 1e9:>12.1f}"
    )
    print(
        f"{'RingVortex.strength read':<45} "
        f"{rv_read_time * 1000:>12.2f} "
        f"{rv_read_time / num_iterations * 1e9:>12.1f}"
    )
    print(
        f"{'RingVortex.age write':<45} "
        f"{rv_write_time * 1000:>12.2f} "
        f"{rv_write_time / num_iterations * 1e9:>12.1f}"
    )
    print(
        f"{'HorseshoeVortex.strength read':<45} "
        f"{hv_read_time * 1000:>12.2f} "
        f"{hv_read_time / num_iterations * 1e9:>12.1f}"
    )

    print()


# ---------------------------------------------------------------------------
# 4. Object creation speed (timeit) for high volume classes.
# ---------------------------------------------------------------------------


def bench_object_creation() -> None:
    """Measures object construction speed for high volume classes.

    :return: None
    """
    print("=" * 70)
    print("OBJECT CREATION SPEED (timeit)")
    print("=" * 70)

    num_iterations = 10_000

    lv_time = timeit.timeit(
        "_make_line_vortex()",
        globals={"_make_line_vortex": _make_line_vortex},
        number=num_iterations,
    )

    rv_time = timeit.timeit(
        "_make_ring_vortex()",
        globals={"_make_ring_vortex": _make_ring_vortex},
        number=num_iterations,
    )

    panel_time = timeit.timeit(
        "_make_panel()",
        globals={"_make_panel": _make_panel},
        number=num_iterations,
    )

    print(f"{'Class':<25} {'Total (ms)':>12} {'Per object (us)':>16}")
    print("-" * 53)
    print(
        f"{'LineVortex':<25} "
        f"{lv_time * 1000:>12.2f} "
        f"{lv_time / num_iterations * 1e6:>16.2f}"
    )
    print(
        f"{'RingVortex':<25} "
        f"{rv_time * 1000:>12.2f} "
        f"{rv_time / num_iterations * 1e6:>16.2f}"
    )
    print(
        f"{'Panel':<25} "
        f"{panel_time * 1000:>12.2f} "
        f"{panel_time / num_iterations * 1e6:>16.2f}"
    )

    print()


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    print()
    print("Ptera Software __slots__ Benchmark")
    print(f"Python {sys.version}")
    print(f"Has __slots__: {hasattr(_line_vortex.LineVortex, '__slots__')}")
    print()

    bench_per_instance_memory()
    bench_attribute_access()
    bench_object_creation()
    bench_aggregate_memory()

    print("Done.")
