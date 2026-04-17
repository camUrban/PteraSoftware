"""Renames the preview hero graphics to their permanent names.

Run this script after you are satisfied with the current Animate.webp and Draw.webp
produced by load_and_visualize_hero.py. It overwrites the saved hero graphics with the
new previews.
"""

from pathlib import Path

here = Path(__file__).parent

preview_to_final = {
    "Animate.webp": "hero_animated.webp",
    "Draw.webp": "hero_static.webp",
}

for preview_name, final_name in preview_to_final.items():
    preview = here / preview_name
    final = here / final_name
    if not preview.exists():
        print(f"Skipping {preview_name}: file not found.")
        continue
    preview.replace(final)
    print(f"Renamed {preview_name} to {final_name}.")
