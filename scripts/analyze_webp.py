"""Analyze a static or animated WebP file by rendering its frames to PNG.

This helper backs the analyze-webp slash command, but it can also be run by hand. It
decodes a WebP file with the webp package's animation decoder (which treats a static
image as a one-frame animation), reports metadata including transparency statistics, and
renders frames to PNG files for inspection.

Rendered PNGs composite an opaque light red backdrop everywhere outside the frames
themselves, so light red strictly means "not part of any frame." By default, pixels
inside each frame keep their original alpha values, which means transparent regions
inside a frame display as whatever background the PNG viewer uses. Check the info
subcommand's alpha statistics to determine whether a file has transparency. The sheet
and frames subcommands also accept a background color over which to composite the pixels
inside each frame, which shows how the frames would look against that background.
"""

import argparse
import io
import math
import sys
import tempfile
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np
import webp
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure

# The backdrop marks every pixel that is not part of a frame. Light red is used
# because it is unlikely to appear in real rendered output, so it is visually
# unambiguous.
BACKDROP_RGB = np.array([1.0, 0.815, 0.815], dtype=float)

DPI = 110.0

# At this module's DPI constant, this width makes a default 9-tile sheet of 4:3
# frames render at about 2000 x 1600 pixels, which fits within the resolution
# at which current Claude models view images natively (2576 pixels on the long
# edge and roughly 3.6 megapixels) with no downscaling.
SHEET_TILE_WIDTH_INCHES = 6.0


def _decode_frames(webp_bytes: bytes) -> list[tuple[np.ndarray, float]]:
    """Decodes a WebP file's bytes into a list of RGBA frames with timestamps.

    :param webp_bytes: The raw bytes of the WebP file.
    :return: A list of tuples, one per frame in order. Each tuple holds a (height,
        width, 4) ndarray of uint8s representing the frame's RGBA pixel values and a
        float representing the frame's timestamp. The units of the timestamps are
        seconds. A static WebP file decodes to a single frame with a timestamp of 0.0.
    """
    decoder = webp.WebPAnimDecoder.new(webp.WebPData.from_buffer(webp_bytes))
    return [(frame, timestamp_ms / 1000.0) for frame, timestamp_ms in decoder.frames()]


def _sample_indices(first: int, last: int, count: int) -> list[int]:
    """Returns evenly spaced frame indices spanning a range, without duplicates.

    The first and last indices of the range are always included.

    :param first: The first frame index in the range, inclusive.
    :param last: The last frame index in the range, inclusive. It must be greater than
        or equal to first.
    :param count: The requested number of indices. It is clamped to the size of the
        range.
    :return: A sorted list of unique ints representing the sampled frame indices.
    """
    span = last - first
    count = min(count, span + 1)
    if count <= 1:
        return [first]
    positions = [first + round(index * span / (count - 1)) for index in range(count)]
    return sorted(set(positions))


def _render_tiles(
    tiles: list[tuple[int, np.ndarray, float]],
    n_frames: int,
    columns: int,
    tile_width_inches: float,
    out_path: Path,
    background_rgb: np.ndarray | None = None,
) -> None:
    """Renders frames as a labeled grid of tiles and saves the grid as a PNG.

    Every pixel outside the frames themselves (margins, labels, and unused grid cells)
    is composited over an opaque light red backdrop. Pixels inside each frame keep their
    original alpha values unless a background color is given, in which case they are
    composited over that color.

    :param tiles: A list of tuples, one per tile to render, in order. Each tuple holds
        an int representing the frame's index, a (height, width, 4) ndarray of uint8s
        representing the frame's RGBA pixel values, and a float representing the frame's
        timestamp. The units of the timestamps are seconds.
    :param n_frames: The total number of frames in the source file, used for the tile
        labels.
    :param columns: The number of tiles per row of the grid.
    :param tile_width_inches: The rendered width of each tile. The units are inches at
        this module's DPI constant.
    :param out_path: The path at which to save the rendered PNG.
    :param background_rgb: A (3,) ndarray of floats representing the RGB components of
        the color over which to composite the pixels inside each frame. The values are
        normalized from 0.0 to 1.0 and are unitless. If None, which is the default, the
        pixels inside each frame keep their original alpha values.
    :return: None
    """
    columns = min(columns, len(tiles))
    rows = math.ceil(len(tiles) / columns)
    frame_height, frame_width = tiles[0][1].shape[:2]
    aspect_ratio = frame_height / frame_width
    fig_width = columns * tile_width_inches
    fig_height = rows * (tile_width_inches * aspect_ratio + 0.45)
    fig = Figure(figsize=(fig_width, fig_height), dpi=DPI)
    canvas = FigureCanvasAgg(fig)
    axs = fig.subplots(rows, columns, squeeze=False)
    artists = []
    for ax in axs.flat:
        ax.set_visible(False)
    for ax, (index, frame, timestamp) in zip(axs.flat, tiles):
        ax.set_visible(True)
        artists.append(ax.imshow(frame))
        if n_frames > 1:
            label = f"frame {index}/{n_frames - 1} (t = {timestamp:.2f} s)"
        else:
            label = "static image"
        ax.set_title(label, fontsize=10)
        ax.axis("off")
    fig.tight_layout()

    # Render the figure with a fully transparent background so the pixels
    # inside each frame keep their original alpha values.
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=DPI, transparent=True)
    buffer.seek(0)
    rendered = mpimg.imread(buffer)
    height, width = rendered.shape[:2]

    # Composite the whole render over the opaque backdrop, then refill each
    # frame's exact extent: with the raw pixels (keeping their original alpha
    # values) by default, or composited over the given background color. Labels
    # land in the margins, so they survive the backdrop compositing.
    alpha = rendered[:, :, 3:4]
    out = np.empty((height, width, 4), dtype=float)
    out[:, :, :3] = rendered[:, :, :3] * alpha + BACKDROP_RGB * (1.0 - alpha)
    out[:, :, 3] = 1.0
    if background_rgb is None:
        refill = rendered
    else:
        refill = np.empty((height, width, 4), dtype=float)
        refill[:, :, :3] = rendered[:, :, :3] * alpha + background_rgb * (1.0 - alpha)
        refill[:, :, 3] = 1.0
    renderer = canvas.get_renderer()
    for artist in artists:
        bbox = artist.get_window_extent(renderer)
        col_0 = int(np.floor(bbox.x0))
        col_1 = int(np.ceil(bbox.x1))
        row_0 = height - int(np.ceil(bbox.y1))
        row_1 = height - int(np.floor(bbox.y0))
        out[row_0:row_1, col_0:col_1] = refill[row_0:row_1, col_0:col_1]
    mpimg.imsave(out_path, np.clip(out, 0.0, 1.0))


def _run_info(frames: list[tuple[np.ndarray, float]], file_path: Path) -> None:
    """Prints metadata about a decoded WebP file.

    :param frames: The decoded frames, as returned by _decode_frames.
    :param file_path: The path of the source WebP file, used for display only.
    :return: None
    """
    n_frames = len(frames)
    frame_height, frame_width = frames[0][0].shape[:2]
    print(f"file: {file_path}")
    if n_frames > 1:
        print(f"type: animated ({n_frames} frames)")
    else:
        print("type: static (1 frame)")
    print(f"size: {frame_width} x {frame_height} pixels")
    if n_frames > 1:
        duration = frames[-1][1]
        if duration > 0.0:
            print(f"duration: {duration:.2f} s ({n_frames / duration:.1f} fps)")
        else:
            print("duration: unavailable (no timestamps)")
    min_alpha = min(int(frame[:, :, 3].min()) for frame, _ in frames)
    non_opaque = float(np.mean([np.mean(frame[:, :, 3] < 255) for frame, _ in frames]))
    if non_opaque > 0.0:
        print(
            f"alpha: has transparency (min alpha is {min_alpha}; "
            f"{100.0 * non_opaque:.1f} percent of pixels are non-opaque)"
        )
    else:
        print("alpha: fully opaque")


def main(argv: list[str]) -> int:
    """Parses the command-line arguments and runs the requested subcommand.

    :param argv: The command-line arguments, excluding the program name.
    :return: An int representing the exit code, which is 0 on success and 1 on a usage
        error.
    """
    parser = argparse.ArgumentParser(
        description="Analyze a static or animated WebP file."
    )
    parser.add_argument("file", type=Path, help="The path to the WebP file.")
    subparsers = parser.add_subparsers(dest="mode", required=True)
    subparsers.add_parser("info", help="Print metadata about the file.")
    sheet_parser = subparsers.add_parser(
        "sheet", help="Render sampled frames as a labeled contact sheet PNG."
    )
    sheet_parser.add_argument(
        "--count",
        type=int,
        default=9,
        help="The number of frames to sample. The default is 9.",
    )
    sheet_parser.add_argument(
        "--columns",
        type=int,
        default=3,
        help="The number of tiles per row. The default is 3.",
    )
    sheet_parser.add_argument(
        "--first",
        type=int,
        default=0,
        help="The first frame index of the range to sample. The default is 0.",
    )
    sheet_parser.add_argument(
        "--last",
        type=int,
        default=None,
        help="The last frame index of the range to sample. The default is the "
        "last frame.",
    )
    sheet_parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="The path at which to save the PNG. The default is a file named "
        "after the input in the system's temporary directory.",
    )
    sheet_parser.add_argument(
        "--background",
        type=str,
        default=None,
        help="A color over which to composite the pixels inside each frame, "
        "given as any matplotlib color specification (for example, 'black' or "
        "'#1a2b3c'). The default is to keep the frames' original alpha values.",
    )
    frames_parser = subparsers.add_parser(
        "frames", help="Render specific frames as full-resolution PNGs."
    )
    frames_parser.add_argument(
        "indices",
        type=int,
        nargs="+",
        help="The frame indices to render.",
    )
    frames_parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="The directory in which to save the PNGs. The default is a "
        "directory named after the input in the system's temporary directory.",
    )
    frames_parser.add_argument(
        "--background",
        type=str,
        default=None,
        help="A color over which to composite the pixels inside each frame, "
        "given as any matplotlib color specification (for example, 'black' or "
        "'#1a2b3c'). The default is to keep the frames' original alpha values.",
    )
    args = parser.parse_args(argv)

    background_rgb = None
    if getattr(args, "background", None) is not None:
        try:
            background_rgb = np.array(to_rgb(args.background), dtype=float)
        except ValueError:
            print(
                f"error: {args.background!r} is not a valid matplotlib color "
                "specification",
                file=sys.stderr,
            )
            return 1

    try:
        webp_bytes = args.file.read_bytes()
    except OSError as exc:
        print(f"error: could not read {args.file} ({exc})", file=sys.stderr)
        return 1
    frames = _decode_frames(webp_bytes)
    n_frames = len(frames)

    if args.mode == "info":
        _run_info(frames, args.file)
        return 0

    if args.mode == "sheet":
        first = args.first
        last = args.last if args.last is not None else n_frames - 1
        if not 0 <= first <= last <= n_frames - 1:
            print(
                f"error: the frame range [{first}, {last}] is invalid for a "
                f"file with {n_frames} frame(s)",
                file=sys.stderr,
            )
            return 1
        if args.count < 1 or args.columns < 1:
            print("error: --count and --columns must be positive", file=sys.stderr)
            return 1
        indices = _sample_indices(first, last, args.count)
        out_path = args.out
        if out_path is None:
            out_path = Path(tempfile.gettempdir()) / f"{args.file.stem}_sheet.png"
        tiles = [(index, frames[index][0], frames[index][1]) for index in indices]
        _render_tiles(
            tiles,
            n_frames,
            args.columns,
            SHEET_TILE_WIDTH_INCHES,
            out_path,
            background_rgb,
        )
        print(f"sampled frame indices: {indices}")
        print(f"saved sheet to: {out_path}")
        return 0

    # The mode is "frames".
    invalid = [index for index in args.indices if not 0 <= index <= n_frames - 1]
    if invalid:
        print(
            f"error: the frame indices {invalid} are invalid for a file with "
            f"{n_frames} frame(s)",
            file=sys.stderr,
        )
        return 1
    out_dir = args.out_dir
    if out_dir is None:
        out_dir = Path(tempfile.gettempdir()) / f"{args.file.stem}_frames"
    out_dir.mkdir(parents=True, exist_ok=True)
    for index in args.indices:
        frame, timestamp = frames[index]
        frame_width = frame.shape[1]
        out_path = out_dir / f"frame_{index:04d}.png"
        _render_tiles(
            [(index, frame, timestamp)],
            n_frames,
            1,
            frame_width / DPI,
            out_path,
            background_rgb,
        )
        print(f"saved frame {index} to: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
