"""Loads the saved hero simulation and generates preview graphics.

This produces animate.webp and draw.webp in the current directory. Run this script
repeatedly while iterating on the simulation or visualization parameters. Once
satisfied, run finalize_and_save_hero.py to promote the previews to the permanent hero
graphics.

Any WebP file larger than the size ceiling is re-rendered at progressively lower quality
until it fits, though never below the quality floor where the visualizations' text stops
being readable.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pterasoftware as ps

_hero_graphics_dir = (
    Path(__file__).resolve().parent.parent.parent / "docs" / "hero_graphics"
)

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

ps.set_up_logging()

loaded_hero_solver = ps.load(_hero_graphics_dir / "hero_solver.psz")
assert isinstance(
    loaded_hero_solver,
    ps.unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver,
)

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

_RENDER_TARGETS: list[tuple[str, Callable[..., None], dict[str, Any]]] = [
    ("draw.webp", ps.output.draw, {"solver": loaded_hero_solver, **_DRAW_KWARGS}),
    (
        "animate.webp",
        ps.output.animate,
        {"unsteady_solver": loaded_hero_solver, **_ANIMATE_KWARGS},
    ),
]

for name, render_func, render_kwargs in _RENDER_TARGETS:
    webp_path = Path(name)
    if not webp_path.exists():
        continue
    original_bytes = webp_path.stat().st_size
    if original_bytes <= _MAX_WEBP_BYTES:
        continue

    quality = _INITIAL_QUALITY
    for _ in range(_MAX_RERENDER_ATTEMPTS):
        quality = max(quality - _QUALITY_STEP, _MIN_QUALITY)
        print(
            f"  {name} is {original_bytes / 1024:.0f} KB, "
            f"re-rendering at quality={quality:.0f}..."
        )
        render_func(**render_kwargs, save=True, quality=quality)

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
