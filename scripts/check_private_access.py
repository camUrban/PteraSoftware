"""Verify that code outside the package reaches Ptera Software only through public
names.

This pre-commit hook parses every tracked Python file in examples/, scripts/, and
validation/, the code cells of every tracked notebook in tutorials/, and the Python code
blocks in README.md. It flags any reference to a Ptera Software name that has a leading
underscore, whether it appears in an import (import pterasoftware._core or from
pterasoftware import _core) or in an attribute chain rooted at a name bound to the
package or to something imported from it (ps._core or ps.geometry._meshing). Dunder
names are allowed.

Each violation is reported with its location and the dotted reference that contains the
internal name.
"""

import ast
import json
import re
import subprocess
import sys
from pathlib import Path

PACKAGE_NAME = "pterasoftware"

# These are the tracked paths outside the package whose Python code is checked.
PYTHON_PATHSPECS = ("examples/*.py", "scripts/*.py", "validation/*.py")
NOTEBOOK_PATHSPECS = ("tutorials/*.ipynb",)
MARKDOWN_PATHSPECS = ("README.md",)

# This matches a fenced Python code block in a Markdown file and captures its body.
PYTHON_FENCE_PATTERN = re.compile(
    r"^```python[ \t]*\n(.*?)^```", re.DOTALL | re.MULTILINE
)


def is_internal(name: str) -> bool:
    """Returns whether a name is internal, meaning it has a leading underscore and is
    not a dunder name.

    :param name: The name to classify.
    :return: True if the name is internal and False otherwise.
    """
    is_dunder = name.startswith("__") and name.endswith("__") and len(name) > 4
    return name.startswith("_") and not is_dunder


def get_attribute_chain(node: ast.Attribute) -> tuple[str, list[str]] | None:
    """Returns the root name and the attribute names of a dotted attribute chain.

    :param node: The outermost Attribute node of the chain.
    :return: A tuple of the root name and the list of attribute names from the root
        outward, or None if the chain is not rooted at a plain name (for example, if it
        is rooted at a call or a subscript).
    """
    attributes: list[str] = []
    current: ast.expr = node
    while isinstance(current, ast.Attribute):
        attributes.append(current.attr)
        current = current.value
    if not isinstance(current, ast.Name):
        return None
    return current.id, attributes[::-1]


def find_violations(source: str) -> list[tuple[int, str]]:
    """Returns (line, reference) tuples for each internal Ptera Software reference in a
    piece of Python source.

    :param source: The Python source to check.
    :return: A list of (line, reference) tuples, one per internal reference, where line
        is the one-based line number within the source.
    """
    tree = ast.parse(source)

    # Collect the names bound to the package or to anything imported from it, so the
    # attribute chains rooted at them can be checked below. Bindings are collected from
    # the whole tree, regardless of scope, which is conservative enough for scripts.
    package_bindings: set[str] = set()
    violations: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] != PACKAGE_NAME:
                    continue
                if any(is_internal(part) for part in parts):
                    violations.append((node.lineno, alias.name))
                package_bindings.add(alias.asname or parts[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level != 0 or node.module is None:
                continue
            parts = node.module.split(".")
            if parts[0] != PACKAGE_NAME:
                continue
            module_is_internal = any(is_internal(part) for part in parts)
            for alias in node.names:
                if module_is_internal or is_internal(alias.name):
                    violations.append((node.lineno, f"{node.module}.{alias.name}"))
                package_bindings.add(alias.asname or alias.name)

    # Check each attribute chain rooted at a package binding. Only the outermost
    # Attribute node of a chain is checked, so a chain is reported once.
    inner_attributes: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Attribute):
            inner_attributes.add(id(node.value))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Attribute) or id(node) in inner_attributes:
            continue
        chain = get_attribute_chain(node)
        if chain is None:
            continue
        root, attributes = chain
        if root not in package_bindings:
            continue
        if any(is_internal(attribute) for attribute in attributes):
            violations.append((node.lineno, ".".join([root, *attributes])))

    return sorted(violations)


def strip_notebook_magics(source: str) -> str:
    """Blanks the IPython magic and shell lines in a notebook cell's source so that it
    parses as Python, keeping the line numbering intact.

    :param source: The cell's source.
    :return: The cell's source with each magic and shell line replaced by an empty line.
    """
    lines = source.split("\n")
    return "\n".join(
        "" if line.lstrip().startswith(("%", "!")) else line for line in lines
    )


def list_tracked_files(pathspecs: tuple[str, ...]) -> list[Path]:
    """Returns the tracked files that match the given pathspecs.

    :param pathspecs: The git pathspecs to match.
    :return: A list of Paths to the matching tracked files, relative to the repository
        root.
    """
    result = subprocess.run(
        ["git", "ls-files", "--", *pathspecs],
        capture_output=True,
        check=True,
        text=True,
    )
    return [Path(line) for line in result.stdout.splitlines() if line]


def main() -> int:
    """Checks every in scope file and prints each violation.

    :return: An int representing the exit code, which is 0 if no violations were found
        and 1 otherwise.
    """
    # Each entry is a location prefix, a piece of source, and the line offset that maps
    # the source's line numbers to the file's.
    sources: list[tuple[str, str, int]] = []
    for path in list_tracked_files(PYTHON_PATHSPECS):
        sources.append((str(path), path.read_text(encoding="utf-8"), 0))
    for path in list_tracked_files(NOTEBOOK_PATHSPECS):
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for cell_index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] != "code":
                continue
            cell_source = "".join(cell["source"])
            sources.append(
                (
                    f"{path} (cell {cell_index})",
                    strip_notebook_magics(cell_source),
                    0,
                )
            )
    for path in list_tracked_files(MARKDOWN_PATHSPECS):
        text = path.read_text(encoding="utf-8")
        for match in PYTHON_FENCE_PATTERN.finditer(text):
            line_offset = text.count("\n", 0, match.start(1))
            sources.append((str(path), match.group(1), line_offset))

    exit_code = 0
    for location, source, line_offset in sources:
        try:
            violations = find_violations(source)
        except SyntaxError as exc:
            print(f"{location}: could not parse ({exc})", file=sys.stderr)
            exit_code = 1
            continue
        for line, reference in violations:
            exit_code = 1
            print(
                f"{location}:{line + line_offset}: {reference} uses an internal name. "
                f"Code outside the package must use only Ptera Software's public names."
            )
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
