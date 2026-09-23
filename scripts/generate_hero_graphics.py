"""Generates the README hero graphics.

Creates and solves the hero simulation, a flapping wing with surface effects (ground
effect), then renders the static and animated hero graphics directly over the committed
files in docs/hero_graphics/. Review the new graphics with git diff, and discard them
with git restore if they are not an improvement.

The simulation is solved once per run, and every render and re-render uses that same
solver.

There is not yet a way to set the camera position programmatically, so draw and animate
each open a window in which to orient the view, and pressing any key starts the render.

Any WebP file larger than the size ceiling is re-rendered at progressively lower quality
until it fits, though never below the quality floor where the visualizations' text stops
being readable. Each re-render opens its window again.

Because the view has to be oriented by hand, the script stops before solving if it has
no display to open windows on, such as when run over SSH. Nothing is overwritten in that
case.
"""

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pyvista as pv

import pterasoftware as ps

_HERO_GRAPHICS_DIR = Path(__file__).resolve().parent.parent / "docs" / "hero_graphics"
_STATIC_PATH = _HERO_GRAPHICS_DIR / "hero_static.webp"
_ANIMATED_PATH = _HERO_GRAPHICS_DIR / "hero_animated.webp"

_MAX_WEBP_BYTES = 5 * 1024 * 1024
_INITIAL_QUALITY = 75.0
_QUALITY_STEP = 25.0
# The lowest quality a re-render will try. Below this, WebP compression makes the
# visualizations' overlay text hard to read. The floor was chosen by inspecting the
# aeroelastic example's animation rendered at qualities from 5 to 95.
_MIN_QUALITY = 25.0
_MAX_RERENDER_ATTEMPTS = 2

_DRAW_KWARGS: dict[str, Any] = {
    "scalar_type": "induced drag",
    "show_wake_vortices": True,
}
_ANIMATE_KWARGS: dict[str, Any] = {
    "scalar_type": "induced drag",
    "show_wake_vortices": True,
}

# The render window classes VTK falls back to when it cannot open an on-screen window,
# after printing a warning that it cannot connect to an X server.
_OFF_SCREEN_WINDOW_CLASSES = ("vtkEGLRenderWindow", "vtkOSOpenGLRenderWindow")


def _create_and_solve_hero() -> (
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver
):
    """Creates and solves the hero simulation.

    :return: The solved hero solver.
    """
    hero_airplane = ps.geometry.airplane.Airplane(
        wings=[
            ps.geometry.wing.Wing(
                wing_cross_sections=[
                    ps.geometry.wing_cross_section.WingCrossSection(
                        num_spanwise_panels=12,
                        chord=1.75,
                        Lp_Wcsp_Lpp=(0.0, 0.0, 0.0),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 0.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing="cosine",
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_Lp=None,
                            resample=True,
                            n_points_per_side=400,
                        ),
                    ),
                    ps.geometry.wing_cross_section.WingCrossSection(
                        num_spanwise_panels=None,
                        chord=1.25,
                        Lp_Wcsp_Lpp=(0.75, 7.0, 0.5),
                        angles_Wcsp_to_Wcs_ixyz=(0.0, 5.0, 0.0),
                        control_surface_symmetry_type="symmetric",
                        control_surface_hinge_point=0.75,
                        control_surface_deflection=0.0,
                        spanwise_spacing=None,
                        airfoil=ps.geometry.airfoil.Airfoil(
                            name="naca2412",
                            outline_A_Lp=None,
                            resample=True,
                            n_points_per_side=400,
                        ),
                    ),
                ],
                name="Main Wing",
                Ler_Gs_Cgs=(0.0, 0.5, 0.0),
                angles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
                symmetric=True,
                mirror_only=False,
                symmetryNormal_G=(0.0, 1.0, 0.0),
                symmetryPoint_G_Cg=(0.0, 0.0, 0.0),
                num_chordwise_panels=6,
                chordwise_spacing="uniform",
            ),
        ],
        name="Hero Airplane",
        Cg_GP1_CgP1=(0.0, 0.0, 0.0),
        weight=0.0,
        s_ref=None,
        c_ref=None,
        b_ref=None,
    )

    main_wing_root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=hero_airplane.wings[0].wing_cross_sections[0],
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
    main_wing_tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=hero_airplane.wings[0].wing_cross_sections[1],
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

    reflected_main_wing_root_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=hero_airplane.wings[1].wing_cross_sections[0],
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
    reflected_main_wing_tip_wing_cross_section_movement = (
        ps.movements.wing_cross_section_movement.WingCrossSectionMovement(
            base_wing_cross_section=hero_airplane.wings[1].wing_cross_sections[1],
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

    main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=hero_airplane.wings[0],
        wing_cross_section_movements=[
            main_wing_root_wing_cross_section_movement,
            main_wing_tip_wing_cross_section_movement,
        ],
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )
    reflected_main_wing_movement = ps.movements.wing_movement.WingMovement(
        base_wing=hero_airplane.wings[1],
        wing_cross_section_movements=[
            reflected_main_wing_root_wing_cross_section_movement,
            reflected_main_wing_tip_wing_cross_section_movement,
        ],
        ampLer_Gs_Cgs=(0.0, 0.0, 0.0),
        periodLer_Gs_Cgs=(0.0, 0.0, 0.0),
        spacingLer_Gs_Cgs=("sine", "sine", "sine"),
        phaseLer_Gs_Cgs=(0.0, 0.0, 0.0),
        ampAngles_Gs_to_Wn_ixyz=(15.0, 0.0, 0.0),
        periodAngles_Gs_to_Wn_ixyz=(1.0, 0.0, 0.0),
        spacingAngles_Gs_to_Wn_ixyz=("sine", "sine", "sine"),
        phaseAngles_Gs_to_Wn_ixyz=(0.0, 0.0, 0.0),
    )

    airplane_movement = ps.movements.airplane_movement.AirplaneMovement(
        base_airplane=hero_airplane,
        wing_movements=[main_wing_movement, reflected_main_wing_movement],
        ampCg_GP1_CgP1=(0.0, 0.0, 0.0),
        periodCg_GP1_CgP1=(0.0, 0.0, 0.0),
        spacingCg_GP1_CgP1=("sine", "sine", "sine"),
        phaseCg_GP1_CgP1=(0.0, 0.0, 0.0),
    )

    hero_operating_point = ps.operating_point.OperatingPoint(
        rho=1.225,
        vCg__E=10.0,
        alpha=30.0,
        beta=0.0,
        CgP1_E_Eo=(0.0, 0.0, -6.0),
        surfaceNormal_E=(0.0, 0.0, 1.0),
        surfacePoint_E_Eo=(0.0, 0.0, 0.0),
        externalFX_W=0.0,
        nu=15.06e-6,
    )

    operating_point_movement = (
        ps.movements.operating_point_movement.OperatingPointMovement(
            base_operating_point=hero_operating_point,
            periodVCg__E=0.0,
            spacingVCg__E="sine",
        )
    )

    movement = ps.movements.movement.Movement(
        airplane_movements=[airplane_movement],
        operating_point_movement=operating_point_movement,
        delta_time=None,
        num_cycles=3,
        num_chords=None,
        num_steps=None,
    )

    hero_problem = ps.problems.UnsteadyProblem(
        movement=movement,
    )

    hero_solver = (
        ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver(
            unsteady_problem=hero_problem,
        )
    )

    hero_solver.run(
        prescribed_wake=False,
        show_progress=True,
    )

    return hero_solver


def _render_within_size_ceiling(
    webp_path: Path,
    render_func: Callable[..., None],
    render_kwargs: dict[str, Any],
) -> None:
    """Renders a WebP file, then re-renders it at progressively lower quality while it
    is larger than the size ceiling, never going below the quality floor.

    :param webp_path: The path to write the WebP file to.
    :param render_func: The output function that renders the file, either ps.output.draw
        or ps.output.animate.
    :param render_kwargs: The keyword arguments to pass to render_func, other than save,
        path, and quality.
    :return: None
    """
    render_func(**render_kwargs, save=True, path=webp_path)

    original_bytes = webp_path.stat().st_size
    if original_bytes <= _MAX_WEBP_BYTES:
        return

    name = webp_path.name
    quality = _INITIAL_QUALITY
    for _ in range(_MAX_RERENDER_ATTEMPTS):
        quality = max(quality - _QUALITY_STEP, _MIN_QUALITY)
        print(
            f"  {name} is {original_bytes / 1024:.0f} KB, "
            f"re-rendering at quality={quality:.0f}..."
        )
        render_func(**render_kwargs, save=True, path=webp_path, quality=quality)

        new_bytes = webp_path.stat().st_size
        if new_bytes <= _MAX_WEBP_BYTES or quality == _MIN_QUALITY:
            break

    new_bytes = webp_path.stat().st_size
    if new_bytes <= _MAX_WEBP_BYTES:
        print(
            f"Re-rendered {name}: "
            f"{original_bytes / 1024:.0f} KB -> {new_bytes / 1024:.0f} KB "
            f"(quality={quality:.0f})"
        )
    else:
        print(
            f"Warning: {name} is still {new_bytes / 1024:.0f} KB after "
            f"re-rendering down to quality={quality:.0f}. Quality is never "
            f"reduced below {_MIN_QUALITY:.0f}, where the text stops being "
            f"readable."
        )


def main() -> int:
    """Solves the hero simulation and renders the hero graphics.

    :return: An int representing the exit code, which is always 0.
    """
    # Creating a Plotter picks its render window class without opening a window, so this
    # detects a missing display before any time is spent solving.
    probe = pv.Plotter()
    render_window = probe.ren_win
    assert render_window is not None
    window_class = render_window.GetClassName()
    probe.close()
    if window_class in _OFF_SCREEN_WINDOW_CLASSES:
        raise RuntimeError(
            f"VTK could not open an on-screen window and fell back to {window_class}, "
            f"so there is no way to orient the hero graphics' view. Run this script "
            f"from a session with a display, not over SSH."
        )

    ps.set_up_logging()

    hero_solver = _create_and_solve_hero()

    _render_within_size_ceiling(
        _STATIC_PATH, ps.output.draw, {"solver": hero_solver, **_DRAW_KWARGS}
    )
    _render_within_size_ceiling(
        _ANIMATED_PATH,
        ps.output.animate,
        {"unsteady_solver": hero_solver, **_ANIMATE_KWARGS},
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
