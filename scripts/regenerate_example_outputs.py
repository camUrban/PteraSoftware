"""Regenerates the expected outputs for the example scripts.

Runs each example script in examples/ and collects its outputs (WebPs, PNGs, CSVs, and
JSON files) into the corresponding subdirectory under docs/examples_expected_output/.
When run without arguments, wipes and rebuilds the entire output tree. When given an
example file name, regenerates only that example's subdirectory.

An example chooses its own output destinations, so a file may land in a subdirectory of
the one the example runs in. Everything here that walks an example's output therefore
walks it recursively.

After each example runs, any WebP file larger than the size ceiling is re-rendered from
the saved solver at progressively lower quality until it fits. This avoids generation
loss from recompression and keeps the expected output tree suitable for display on
GitHub and ReadTheDocs without manual size management.
"""

import argparse
import ast
import os
import shutil
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "examples"
OUTPUT_DIR = PROJECT_ROOT / "docs" / "examples_expected_output"

_MAX_WEBP_BYTES = 5 * 1024 * 1024
_INITIAL_QUALITY = 75.0
_QUALITY_STEP = 35.0
_MAX_RERENDER_ATTEMPTS = 2

# The keyword arguments this script supplies itself, which override whatever an example
# passes. Every other keyword is forwarded to the re-render so that it reproduces the
# example's output. path is deliberately absent: a re-render has to write back over the
# file it is shrinking, so it needs the destination the example chose.
_MANAGED_KWARGS = {"solver", "unsteady_solver", "save", "quality", "testing"}

# The destinations draw and animate write to when an example names none.
_DEFAULT_DRAW_PATH = "draw.webp"
_DEFAULT_ANIMATE_PATH = "animate.webp"


def _extract_output_kwargs(
    script_path: Path,
) -> tuple[dict[str, object] | None, dict[str, object] | None]:
    """Extracts the keyword arguments from ps.output.draw and ps.output.animate calls.

    Parses the example script's AST and returns the keyword arguments (excluding the
    managed kwargs solver, unsteady_solver, save, quality, and testing) for each call.
    Only literal values are extracted, so a call that builds an argument at runtime
    yields no entry for it and is treated as though it had used the default.

    :param script_path: The path to the example script to parse.
    :return: A tuple of (draw_kwargs, animate_kwargs). Either value is None if the
        corresponding call is not found in the script.
    """
    tree = ast.parse(script_path.read_text(), filename=str(script_path))

    draw_kwargs: dict[str, object] | None = None
    animate_kwargs: dict[str, object] | None = None

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        func = node.func
        if not (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Attribute)
            and func.value.attr == "output"
            and isinstance(func.value.value, ast.Name)
            and func.value.value.id == "ps"
        ):
            continue

        if func.attr not in ("draw", "animate"):
            continue

        kwargs: dict[str, object] = {}
        for kw in node.keywords:
            if kw.arg is None or kw.arg in _MANAGED_KWARGS:
                continue
            try:
                kwargs[kw.arg] = ast.literal_eval(kw.value)
            except (ValueError, TypeError):
                continue

        if func.attr == "draw":
            draw_kwargs = kwargs
        elif func.attr == "animate":
            animate_kwargs = kwargs

    return draw_kwargs, animate_kwargs


def _find_solver_file(directory: Path) -> Path | None:
    """Finds the saved solver JSON file in a directory or any of its subdirectories.

    :param directory: The directory to search.
    :return: The path to the solver file, or None if the directory tree does not hold
        exactly one.
    """
    candidates = sorted(directory.rglob("*.json.gz"))
    if len(candidates) == 1:
        return candidates[0]
    return None


def _output_destination(
    output_subdir: Path, kwargs: dict[str, object], default_path: str
) -> Path:
    """Returns the absolute path an output call writes its file to.

    :param output_subdir: The directory the example runs in, which a relative
        destination is resolved against.
    :param kwargs: The keyword arguments extracted from the output call.
    :param default_path: The destination the output function uses when the call names
        none.
    :return: The resolved absolute path.
    """
    return (output_subdir / Path(str(kwargs.get("path", default_path)))).resolve()


def _rerender_oversized_webps(output_subdir: Path, script_path: Path) -> None:
    """Re-renders any oversized WebP files by loading the saved solver and calling draw
    or animate at progressively lower quality.

    :param output_subdir: The directory containing the example's output files.
    :param script_path: The path to the example script, used to extract the original
        draw/animate keyword arguments.
    :return: None
    """
    oversized = [
        p
        for p in sorted(output_subdir.rglob("*.webp"))
        if p.stat().st_size > _MAX_WEBP_BYTES
    ]
    if not oversized:
        return

    solver_file = _find_solver_file(output_subdir)
    if solver_file is None:
        print("    Cannot re-render: no saved solver found.")
        return

    draw_kwargs, animate_kwargs = _extract_output_kwargs(script_path)

    import pterasoftware as ps

    loaded_solver = ps.load(solver_file)

    # Match each oversized file to the call that produced it by destination rather than
    # by name, since an example can name its own. A call whose destination is not a
    # literal has none to match against, so its file is reported rather than re-rendered
    # under a name it never used.
    renderers: dict[Path, tuple[Callable[..., None], dict[str, Any]]] = {}
    if draw_kwargs is not None:
        draw_destination = _output_destination(
            output_subdir, draw_kwargs, _DEFAULT_DRAW_PATH
        )
        renderers[draw_destination] = (
            ps.output.draw,
            {**draw_kwargs, "solver": loaded_solver},
        )
    if animate_kwargs is not None:
        animate_destination = _output_destination(
            output_subdir, animate_kwargs, _DEFAULT_ANIMATE_PATH
        )
        renderers[animate_destination] = (
            ps.output.animate,
            {**animate_kwargs, "unsteady_solver": loaded_solver},
        )

    for webp_path in oversized:
        original_bytes = webp_path.stat().st_size
        name = str(webp_path.relative_to(output_subdir))

        renderer = renderers.get(webp_path.resolve())
        if renderer is None:
            print(f"    Cannot re-render {name}: no matching output call found.")
            continue
        render_func, render_kwargs = renderer

        original_cwd = os.getcwd()
        quality = _INITIAL_QUALITY

        try:
            os.chdir(output_subdir)

            for _ in range(_MAX_RERENDER_ATTEMPTS):
                quality -= _QUALITY_STEP
                render_func(**render_kwargs, save=True, quality=quality, testing=True)

                new_bytes = webp_path.stat().st_size
                if new_bytes <= _MAX_WEBP_BYTES:
                    print(
                        f"    Re-rendered {name}: "
                        f"{original_bytes / 1024:.0f} KB -> "
                        f"{new_bytes / 1024:.0f} KB "
                        f"(quality={quality:.1f})"
                    )
                    break
            else:
                final_bytes = webp_path.stat().st_size
                print(
                    f"    Warning: {name} is still "
                    f"{final_bytes / 1024:.0f} KB after "
                    f"{_MAX_RERENDER_ATTEMPTS} re-render attempts "
                    f"(final quality={quality:.1f})."
                )
        finally:
            os.chdir(original_cwd)


def _discover_examples() -> list[Path]:
    """Discovers all example scripts in the examples directory.

    :return: A sorted list of Paths to example scripts, excluding __init__.py.
    """
    return sorted(p for p in EXAMPLES_DIR.glob("*.py") if p.name != "__init__.py")


_SUBPROCESS_WRAPPER = """\
import runpy
import sys
from unittest.mock import patch

import matplotlib.pyplot as plt
import pterasoftware as ps

original_draw = ps.output.draw
original_animate = ps.output.animate


def _draw_testing(*args, **kwargs):
    kwargs["testing"] = True
    original_draw(*args, **kwargs)


def _animate_testing(*args, **kwargs):
    kwargs["testing"] = True
    original_animate(*args, **kwargs)


with (
    patch.object(ps.output, "draw", _draw_testing),
    patch.object(ps.output, "animate", _animate_testing),
    patch.object(plt, "show", lambda *_args, **_kwargs: None),
):
    runpy.run_path(sys.argv[1], run_name="__main__")
"""


def _run_example(script_path: Path, output_subdir: Path) -> bool:
    """Runs a single example script in a subprocess with its outputs directed to a
    subdirectory.

    Each example runs in its own process so memory is fully reclaimed between examples.
    The subprocess patches ps.output.draw and ps.output.animate to force testing=True so
    PyVista windows close automatically, and patches plt.show to prevent blocking. The
    working directory is set to the output subdirectory so all saved files land there.

    :param script_path: The path to the example script to run.
    :param output_subdir: The directory in which to collect the script's outputs.
    :return: True if the script ran successfully, False otherwise.
    """
    output_subdir.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [sys.executable, "-u", "-c", _SUBPROCESS_WRAPPER, str(script_path)],
        cwd=output_subdir,
    )

    if result.returncode != 0:
        print(f"  FAILED (exit code {result.returncode})")
        return False

    return True


def main() -> int:
    """Regenerates example outputs for all or a single example script.

    :return: An int representing the exit code, which is 0 if all examples succeeded and
        1 if any failed.
    """
    parser = argparse.ArgumentParser(
        description="Regenerate expected outputs for example scripts."
    )
    parser.add_argument(
        "example",
        nargs="?",
        default=None,
        help="The file name of a single example to regenerate (e.g., "
        "steady_horseshoe_vortex_lattice_method_solver.py). When omitted, "
        "all examples are regenerated.",
    )
    args = parser.parse_args()

    if args.example is not None:
        script_path = EXAMPLES_DIR / args.example
        if not script_path.exists():
            print(f"Example not found: {script_path}")
            return 1
        examples = [script_path]
    else:
        examples = _discover_examples()
        if not examples:
            print("No example scripts found.")
            return 1

        # Wipe and recreate the entire output directory.
        if OUTPUT_DIR.exists():
            shutil.rmtree(OUTPUT_DIR)
        OUTPUT_DIR.mkdir(parents=True)

    succeeded = []
    failed = []

    for script_path in examples:
        name = script_path.stem
        output_subdir = OUTPUT_DIR / name

        # For single-example mode, wipe only that subdirectory.
        if args.example is not None and output_subdir.exists():
            shutil.rmtree(output_subdir)

        print(f"Running {script_path.name}...", flush=True)
        if _run_example(script_path, output_subdir):
            n_files = sum(1 for p in output_subdir.rglob("*") if p.is_file())
            if n_files > 0:
                print(
                    f"  Saved {n_files} file(s) to {output_subdir.relative_to(PROJECT_ROOT)}"
                )
                _rerender_oversized_webps(output_subdir, script_path)
                solver_file = _find_solver_file(output_subdir)
                if solver_file is not None:
                    solver_file.unlink()
                succeeded.append(name)
            else:
                # An example that writes nothing can still have left the subdirectories
                # it meant to write into, so the whole tree goes rather than the one
                # directory.
                shutil.rmtree(output_subdir)
                print("  No output files produced, removed empty directory.")
                succeeded.append(name)
        else:
            failed.append(name)

    print()
    print(f"Done: {len(succeeded)} succeeded, {len(failed)} failed.")
    if failed:
        print("Failed examples:")
        for name in failed:
            print(f"  {name}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
