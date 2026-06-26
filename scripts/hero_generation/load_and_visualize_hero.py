"""Loads the saved hero simulation and generates preview graphics.

This produces Animate.webp and Draw.webp in the current directory. Run this script
repeatedly while iterating on the simulation or visualization parameters. Once
satisfied, run finalize_and_save_hero.py to promote the previews to the permanent hero
graphics.

Any WebP file larger than the size ceiling is re-rendered at progressively lower quality
(dropping by 35 each attempt) until it fits.
"""

from pathlib import Path

import pterasoftware as ps

_hero_graphics_dir = (
    Path(__file__).resolve().parent.parent.parent / "docs" / "hero_graphics"
)

_MAX_WEBP_BYTES = 5 * 1024 * 1024
_INITIAL_QUALITY = 75.0
_QUALITY_STEP = 35.0
_MAX_RERENDER_ATTEMPTS = 2

_DRAW_KWARGS = {
    "scalar_type": "induced drag",
    "show_wake_vortices": True,
}
_ANIMATE_KWARGS = {
    "scalar_type": "induced drag",
    "show_wake_vortices": True,
}

ps.set_up_logging()

loaded_hero_solver = ps.load(_hero_graphics_dir / "hero_solver.json.gz")

ps.output.draw(
    solver=loaded_hero_solver,
    **_DRAW_KWARGS,
    save=True,
)

ps.output.animate(
    unsteady_solver=loaded_hero_solver,
    **_ANIMATE_KWARGS,
    save=True,
)

for name, render_func, render_kwargs in [
    ("Draw.webp", ps.output.draw, {"solver": loaded_hero_solver, **_DRAW_KWARGS}),
    (
        "Animate.webp",
        ps.output.animate,
        {"unsteady_solver": loaded_hero_solver, **_ANIMATE_KWARGS},
    ),
]:
    webp_path = Path(name)
    if not webp_path.exists():
        continue
    original_bytes = webp_path.stat().st_size
    if original_bytes <= _MAX_WEBP_BYTES:
        continue

    quality = _INITIAL_QUALITY
    for _ in range(_MAX_RERENDER_ATTEMPTS):
        quality -= _QUALITY_STEP
        print(
            f"  {name} is {original_bytes / 1024:.0f} KB, "
            f"re-rendering at quality={quality:.0f}..."
        )
        ps.output._quality = quality
        render_func(**render_kwargs, save=True)
        ps.output._quality = _INITIAL_QUALITY

        new_bytes = webp_path.stat().st_size
        if new_bytes <= _MAX_WEBP_BYTES:
            print(
                f"Re-rendered {name}: "
                f"{original_bytes / 1024:.0f} KB -> {new_bytes / 1024:.0f} KB "
                f"(quality={quality:.0f})"
            )
            break
    else:
        final_bytes = webp_path.stat().st_size
        print(
            f"Warning: {name} is still {final_bytes / 1024:.0f} KB "
            f"after {_MAX_RERENDER_ATTEMPTS} re-render attempts "
            f"(final quality={quality:.0f})."
        )
