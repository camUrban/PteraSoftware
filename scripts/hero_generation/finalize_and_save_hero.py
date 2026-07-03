"""Renames the preview hero graphics to their permanent names.

Run this script after you are satisfied with the current animate.webp and draw.webp
produced by load_and_visualize_hero.py. It overwrites the saved hero graphics with the
new previews.
"""

from pathlib import Path

_hero_graphics_dir = (
    Path(__file__).resolve().parent.parent.parent / "docs" / "hero_graphics"
)

preview_to_final = {
    "animate.webp": "hero_animated.webp",
    "draw.webp": "hero_static.webp",
}

for preview_name, final_name in preview_to_final.items():
    preview = Path(preview_name)
    final = _hero_graphics_dir / final_name
    if not preview.exists():
        print(f"Skipping {preview_name}: file not found.")
        continue
    preview.replace(final)
    print(f"Renamed {preview_name} to {final_name}.")
