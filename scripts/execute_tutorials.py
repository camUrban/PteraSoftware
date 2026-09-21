"""Executes the tutorial notebooks and stores their outputs.

Runs each notebook in tutorials/ from top to bottom and writes the resulting outputs
back into the notebook file so the documentation site can render them without executing
anything. When run without arguments, executes every notebook. When given a notebook
file name, executes only that one.

The kernel runs with git pointed at an empty directory. Ptera Software's save and load
functions record and check git provenance when the package lives inside a repository,
and those checks warn whenever the working tree is dirty, which it always is while a
tutorial is being edited. Hiding the repository makes the kernel behave like a
pterasoftware installation from PyPI, so the stored outputs match what a user sees.

The kernel also renders PyVista scenes off screen with PyVista's Jupyter backend turned
off, so a notebook can call draw or animate the ordinary way and the script still runs
unattended. Run interactively, those calls display a static render inline in a notebook,
or open a window that waits for a key press outside one.

Each notebook's working directory is tutorials/, so files a notebook writes, such as a
saved render, land next to it. Images must be saved to files and embedded from markdown
cells rather than displayed inline, which keeps base64 image data out of the notebook
JSON. The script fails if any cell stores an image output.

After execution, consecutive stream outputs on the same stream are merged into one,
which is how a Jupyter frontend stores them, and all cell metadata and the notebook's
language_info block are stripped so the file holds no execution timestamps and no
details about the machine that ran it.
"""

import argparse
import importlib.util
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

import nbformat
import zmq
from nbclient import NotebookClient

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TUTORIALS_DIR = PROJECT_ROOT / "tutorials"

_KERNEL_NAME = "python3"
_CELL_TIMEOUT_S = 600


def _discover_notebooks() -> list[Path]:
    """Discovers all notebooks in the tutorials directory.

    :return: A sorted list of Paths to the notebooks.
    """
    return sorted(TUTORIALS_DIR.glob("*.ipynb"))


def _merge_stream_outputs(
    outputs: list[nbformat.NotebookNode],
) -> list[nbformat.NotebookNode]:
    """Merges consecutive stream outputs on the same stream into one output.

    A logging handler that flushes each line separately produces one stream output per
    line, whereas a Jupyter frontend stores consecutive text on the same stream as a
    single output.

    :param outputs: The list of a cell's outputs in order.
    :return: A new list of outputs with consecutive same-stream outputs merged.
    """
    merged: list[nbformat.NotebookNode] = []
    for output in outputs:
        previous = merged[-1] if merged else None
        if (
            previous is not None
            and output.output_type == "stream"
            and previous.output_type == "stream"
            and output.name == previous.name
        ):
            previous.text += output.text
        else:
            merged.append(output)
    return merged


def _execute_notebook(notebook_path: Path) -> bool:
    """Executes a single notebook in place and cleans up the stored result.

    :param notebook_path: The path to the notebook to execute.
    :return: True if the notebook executed successfully, False otherwise.
    """
    nb = nbformat.read(notebook_path, as_version=4)

    # Point git at an empty directory for the kernel so the package's provenance checks
    # take their no-repository path. Render PyVista scenes off screen, and turn off
    # PyVista's Jupyter backend so it doesn't display a static render inline as an image
    # output, so draw and animate just write their files. See the module docstring for
    # why.
    with tempfile.TemporaryDirectory() as git_free_dir:
        kernel_env = dict(os.environ)
        kernel_env["GIT_DIR"] = git_free_dir
        kernel_env["PYVISTA_OFF_SCREEN"] = "true"
        kernel_env["PYVISTA_JUPYTER_BACKEND"] = "none"

        client = NotebookClient(
            nb,
            timeout=_CELL_TIMEOUT_S,
            kernel_name=_KERNEL_NAME,
            record_timing=False,
            resources={"metadata": {"path": str(TUTORIALS_DIR)}},
        )
        # Have the kernel manager provision CurveZMQ keys when pyzmq supports them,
        # which encrypts the kernel's local sockets and keeps ipykernel from warning
        # that it is running over TCP without encryption. The "auto" policy falls back
        # to an unencrypted connection if the kernelspec doesn't declare support.
        start_kwargs: dict[str, Any] = {"env": kernel_env}
        if zmq.has("curve"):
            start_kwargs["transport_encryption"] = "auto"

        try:
            client.execute(**start_kwargs)
        except Exception as error:
            print(f"  FAILED: {error}")
            return False

    for cell in nb.cells:
        cell.metadata = nbformat.NotebookNode()
        if cell.cell_type != "code":
            continue
        cell.outputs = _merge_stream_outputs(cell.outputs)
        for output in cell.outputs:
            data = output.get("data", {})
            image_types = [mime for mime in data if mime.startswith("image/")]
            if image_types:
                print(
                    f"  FAILED: cell {cell.id} stored a {image_types[0]} output. Save "
                    "figures to files and embed them from markdown cells instead."
                )
                return False

    nb.metadata = nbformat.NotebookNode(kernelspec=nb.metadata.kernelspec)
    nbformat.write(nb, notebook_path)
    return True


def main() -> int:
    """Executes all or a single tutorial notebook.

    :return: An int representing the exit code, which is 0 if all notebooks succeeded
        and 1 if any failed.
    """
    parser = argparse.ArgumentParser(
        description="Execute the tutorial notebooks and store their outputs."
    )
    parser.add_argument(
        "notebook",
        nargs="?",
        default=None,
        help="The file name of a single notebook to execute (e.g., "
        "getting_started.ipynb). When omitted, all notebooks are executed.",
    )
    args = parser.parse_args()

    # Without ipywidgets, tqdm warns on import inside a Jupyter kernel, and that warning
    # would be stored as an output of the first cell that imports pterasoftware.
    if importlib.util.find_spec("ipywidgets") is None:
        print(
            "ipywidgets is not installed. Install the development requirements with "
            "pip install -r requirements.txt -r requirements_dev.txt and try again."
        )
        return 1

    if args.notebook is not None:
        notebook_path = TUTORIALS_DIR / args.notebook
        if not notebook_path.exists():
            print(f"Notebook not found: {notebook_path}")
            return 1
        notebooks = [notebook_path]
    else:
        notebooks = _discover_notebooks()
        if not notebooks:
            print("No notebooks found.")
            return 1

    succeeded = []
    failed = []

    for notebook_path in notebooks:
        print(f"Executing {notebook_path.name}...", flush=True)
        if _execute_notebook(notebook_path):
            print(f"  Stored outputs in {notebook_path.relative_to(PROJECT_ROOT)}")
            succeeded.append(notebook_path.name)
        else:
            failed.append(notebook_path.name)

    print()
    print(f"Done: {len(succeeded)} succeeded, {len(failed)} failed.")
    if failed:
        print("Failed notebooks:")
        for name in failed:
            print(f"  {name}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
