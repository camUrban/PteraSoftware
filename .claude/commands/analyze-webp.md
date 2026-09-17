---
description: Analyze a static or animated WebP file by viewing its metadata and rendered frames
---

# Analyze a WebP

Inspect the WebP file given in the arguments and answer any question that accompanies it. If the arguments contain no question, summarize what the file shows and call out anything that looks wrong.

The Read tool only ever shows the first frame of a WebP file, so never judge an animated WebP by Reading it directly. Use `scripts/analyze_webp.py`, which renders frames to PNG files for viewing.

Arguments: $ARGUMENTS

## Steps

1. **Probe the file** to learn whether it is static or animated, its frame count, its dimensions, its duration, and whether it has transparency:
   ```bash
   python -u scripts/analyze_webp.py <file> info
   ```
   If multiple files were given, analyze each one in turn with these same steps.
2. **Render and view the frames.** The script prints the path of each PNG it saves. View them with the Read tool.
    - For a static file (1 frame), render it on its own:
      ```bash
      python -u scripts/analyze_webp.py <file> frames 0
      ```
    - For an animated file, render a contact sheet of frames sampled evenly across the animation:
      ```bash
      python -u scripts/analyze_webp.py <file> sheet
      ```
      The default samples 9 frames into a 3 x 3 grid, which is usually sufficient on its own to answer the question. To inspect more frames, render multiple sheets over sub-ranges with `--first` and `--last` rather than raising `--count` on one sheet. Tiles much smaller than the default grid's are illegible.
3. **Drill down only where needed.** The sheet's tiles render at about two-thirds of these models' native viewing resolution, so they carry nearly as much visible detail as full-resolution renders. Do not rerender frames individually just because they appeared as tiles. Reserve full-resolution renders for when a specific tile shows something that genuinely needs closer inspection (an apparent artifact, fine text, a subtle geometry question) and name that reason before rendering:
   ```bash
   python -u scripts/analyze_webp.py <file> frames <index> [<index> ...]
   ```
4. **Check contrast against a specific background if asked.** Both `sheet` and `frames` accept `--background <color>` (any matplotlib color specification, such as `black` or `#1a2b3c`), which composites the pixels inside each frame over that color instead of keeping their alpha. Use this when the question involves how the content will look on a particular background, for example whether a transparent animation stays legible on a dark page.
5. **Answer.** Address the question from the arguments (or summarize the file), citing specific frame indices and timestamps as evidence.

## Interpreting the rendered PNGs

- Light red strictly means "not part of any frame": it fills the margins, the labels' background, and any unused grid cells. It never appears inside a frame's rectangle unless the frame itself contains that color.
- Without `--background`, pixels inside each frame keep the original file's alpha values. Transparent regions inside a frame therefore display as whatever background the viewer composites them over (typically white in the Read tool), and are NOT visually distinguishable from opaque pixels of that same color. To determine whether the file has transparency, rely on the `info` subcommand's alpha statistics, not on the rendered PNGs.
- With `--background`, pixels inside each frame are opaque and show exactly how the frame composites over that color.
- Each tile's label gives its frame index and, for animated files, its timestamp in seconds.
