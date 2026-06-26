"""Loads the saved hero simulation and generates preview graphics.

This produces Animate.webp and Draw.webp in the current directory. Run this script
repeatedly while iterating on the simulation or visualization parameters. Once
satisfied, run finalize_and_save_hero.py to promote the previews to the permanent hero
graphics.
"""

from pathlib import Path

import pterasoftware as ps

_hero_graphics_dir = (
    Path(__file__).resolve().parent.parent.parent / "docs" / "hero_graphics"
)

ps.set_up_logging()

loaded_hero_solver = ps.load(_hero_graphics_dir / "hero_solver.json.gz")

ps.output.draw(
    solver=loaded_hero_solver,
    scalar_type="induced drag",
    show_wake_vortices=True,
    save=True,
)

ps.output.animate(
    unsteady_solver=loaded_hero_solver,
    scalar_type="induced drag",
    show_wake_vortices=True,
    save=True,
)
