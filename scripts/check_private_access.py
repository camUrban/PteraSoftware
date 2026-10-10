"""Verify that code outside the package reaches Ptera Software only through public
names, and that the package names its modules and members by its naming rules.

This pre-commit hook parses every tracked Python file in examples/, scripts/, and
validation/, the code cells of every tracked notebook in tutorials/, and the Python code
blocks in README.md. It flags any reference to a Ptera Software name that has a leading
underscore, whether it appears in an import (import pterasoftware._core or from
pterasoftware import _core) or in an attribute chain rooted at a name bound to the
package or to something imported from it (ps._core or ps.Movement._lcm_period). Dunder
names are allowed.

It also parses every tracked Python file in pterasoftware/ and flags any name that
breaks one of the package's naming rules, which follow.

Every module and package directly inside pterasoftware/, other than __init__.py, starts
with an underscore, and the modules and packages nested inside them do not. The
deprecated modules and packages at the old public paths, which are the keys of the
package's _LAZY_MODULES table, are exempt.

A name defined at module level in pterasoftware/__init__.py starts with an underscore if
and only if it is not in __all__. A name defined at module level in any other module
never starts with an underscore.

A member defined on an internal-only class never starts with an underscore, unless it
backs a property of the same name without the underscore. A class is internal-only if it
is not a public class or a base of one, where the public classes are the classes in the
package's _LAZY_CALLABLES table. A member is a name defined in the class body, an entry
in its __slots__, or an attribute assigned through self in one of its methods.

Imported names, function-local names, and dunder names are exempt from these rules. Base
classes defined outside the package are skipped when resolving class hierarchies.

Each violation is reported with its location and the name that breaks the rule.
"""

import ast
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

PACKAGE_NAME = "pterasoftware"

# These are the tracked paths outside the package whose Python code is checked.
PYTHON_PATHSPECS = ("examples/*.py", "scripts/*.py", "validation/*.py")
NOTEBOOK_PATHSPECS = ("tutorials/*.ipynb",)
MARKDOWN_PATHSPECS = ("README.md",)

# These are the tracked package files whose names are checked against the naming rules.
PACKAGE_PATHSPECS = ("pterasoftware/*.py",)

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


def get_module_name(path: Path) -> str:
    """Returns the dotted name of the module or package that a package file defines.

    :param path: The file's Path, relative to the repository root.
    :return: The dotted module name, which is the package's name for an __init__.py.
    """
    parts = list(path.with_suffix("").parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def get_defined_names(body: list[ast.stmt]) -> list[tuple[int, str]]:
    """Returns the names that a block of statements defines in its own scope.

    The names come from def and class statements, from assignment, for loop, and with
    statement targets, and from the same statements nested inside if, try, for, while,
    and with blocks. Imported names are not included.

    :param body: The statements to scan.
    :return: A list of (line, name) tuples, one per defined name.
    """
    names: list[tuple[int, str]] = []
    for node in body:
        targets: list[ast.expr] = []
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.append((node.lineno, node.name))
        elif isinstance(node, ast.Assign):
            targets.extend(node.targets)
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign, ast.For, ast.AsyncFor)):
            targets.append(node.target)
        elif isinstance(node, (ast.With, ast.AsyncWith)):
            targets.extend(
                item.optional_vars
                for item in node.items
                if item.optional_vars is not None
            )
        for target in targets:
            for target_node in ast.walk(target):
                if isinstance(target_node, ast.Name):
                    names.append((node.lineno, target_node.id))
        if isinstance(
            node, (ast.If, ast.Try, ast.For, ast.AsyncFor, ast.While, ast.With)
        ):
            for field in ("body", "orelse", "finalbody"):
                names.extend(get_defined_names(getattr(node, field, [])))
            for handler in getattr(node, "handlers", []):
                names.extend(get_defined_names(handler.body))
    return names


def get_literal_assignment(tree: ast.Module, name: str) -> Any:
    """Returns the literal value assigned to a module-level name.

    :param tree: The parsed module.
    :param name: The name whose assigned value to return.
    :return: The evaluated literal value.
    """
    for node in tree.body:
        if isinstance(node, ast.Assign):
            target_names = [
                target.id for target in node.targets if isinstance(target, ast.Name)
            ]
            if name in target_names:
                return ast.literal_eval(node.value)
    raise ValueError(f"pterasoftware/__init__.py does not assign {name}.")


def find_naming_violations() -> list[tuple[str, int, str]]:
    """Returns the names in the package that break its naming rules.

    :return: A list of (path, line, message) tuples, one per violation, where path is
        the file's path relative to the repository root.
    """
    paths = {
        get_module_name(path): path for path in list_tracked_files(PACKAGE_PATHSPECS)
    }
    trees = {
        module: ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for module, path in paths.items()
    }
    init_tree = trees[PACKAGE_NAME]
    public_names = set(get_literal_assignment(init_tree, "__all__"))
    deprecated_names = set(get_literal_assignment(init_tree, "_LAZY_MODULES"))
    lazy_callables = get_literal_assignment(init_tree, "_LAZY_CALLABLES")

    violations: list[tuple[str, int, str]] = []

    # Check the names of the modules and packages. Each module and package is checked by
    # its own name alone, since its enclosing packages are checked by theirs.
    for module, path in paths.items():
        parts = module.split(".")
        if len(parts) == 2:
            if not is_internal(parts[1]) and parts[1] not in deprecated_names:
                violations.append(
                    (
                        str(path),
                        1,
                        f"{parts[1]} is directly inside pterasoftware/, so its name "
                        f"must start with an underscore.",
                    )
                )
        elif len(parts) > 2 and is_internal(parts[-1]):
            violations.append(
                (
                    str(path),
                    1,
                    f"{parts[-1]} is nested inside a package in pterasoftware/, so its "
                    f"name must not start with an underscore.",
                )
            )

    # Check the module-level names. A global statement inside a function also defines a
    # module-level name.
    for module, tree in trees.items():
        path = paths[module]
        defined_names = get_defined_names(tree.body)
        for node in ast.walk(tree):
            if isinstance(node, ast.Global):
                defined_names.extend((node.lineno, name) for name in node.names)
        for line, name in defined_names:
            if module == PACKAGE_NAME:
                is_dunder = name.startswith("__") and name.endswith("__")
                if is_dunder or is_internal(name) != (name in public_names):
                    continue
                if name in public_names:
                    message = (
                        f"{name} is in __all__, so it must not start with an "
                        f"underscore."
                    )
                else:
                    message = (
                        f"{name} is not in __all__, so it must start with an "
                        f"underscore."
                    )
                violations.append((str(path), line, message))
            elif is_internal(name):
                violations.append(
                    (
                        str(path),
                        line,
                        f"{name} is a module-level name outside "
                        f"pterasoftware/__init__.py, so it must not start with an "
                        f"underscore.",
                    )
                )

    # Check the members of the internal-only classes. First, find each module's import
    # bindings and module-level classes, so that each class's bases can be resolved to
    # the classes they name.
    bindings: dict[str, dict[str, tuple[str, str]]] = {}
    classes: dict[tuple[str, str], ast.ClassDef] = {}
    for module, tree in trees.items():
        is_package = paths[module].name == "__init__.py"
        module_bindings: dict[str, tuple[str, str]] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                if node.level == 0:
                    source = str(node.module)
                else:
                    base_parts = module.split(".")
                    if not is_package:
                        base_parts = base_parts[:-1]
                    base_parts = base_parts[: len(base_parts) - node.level + 1]
                    if node.module is not None:
                        base_parts.append(node.module)
                    source = ".".join(base_parts)
                for alias in node.names:
                    module_bindings[alias.asname or alias.name] = (source, alias.name)
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.asname is None:
                        root = alias.name.split(".")[0]
                        module_bindings[root] = ("", root)
                    else:
                        module_bindings[alias.asname] = ("", alias.name)
        bindings[module] = module_bindings
        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                classes[(module, node.name)] = node

    def resolve(module: str, expression: ast.expr) -> str | tuple[str, str] | None:
        """Resolves an expression in a module to the package module or class it names.

        :param module: The dotted name of the module the expression appears in.
        :param expression: The expression to resolve.
        :return: The dotted name of a package module, a (module, name) tuple naming a
            package class, or None if the expression names neither.
        """
        if isinstance(expression, ast.Subscript):
            return resolve(module, expression.value)
        if isinstance(expression, ast.Attribute):
            owner = resolve(module, expression.value)
            if isinstance(owner, str):
                return resolve_attribute(owner, expression.attr, set())
            return None
        if isinstance(expression, ast.Name):
            return resolve_attribute(module, expression.id, set())
        return None

    def resolve_attribute(
        module: str, name: str, visited: set[tuple[str, str]]
    ) -> str | tuple[str, str] | None:
        """Resolves a name looked up on a package module.

        :param module: The dotted name of the module.
        :param name: The name to look up on it.
        :param visited: The (module, name) lookups already followed, which guards
            against import cycles.
        :return: The dotted name of a package module, a (module, name) tuple naming a
            package class, or None if the name is neither.
        """
        if (module, name) in visited:
            return None
        visited.add((module, name))
        if (module, name) in classes:
            return module, name
        if name in bindings.get(module, {}):
            source, imported_name = bindings[module][name]
            if not source:
                return imported_name if imported_name in trees else None
            if f"{source}.{imported_name}" in trees:
                return f"{source}.{imported_name}"
            if source in trees:
                return resolve_attribute(source, imported_name, visited)
            return None
        submodule = f"{module}.{name}" if module else name
        return submodule if submodule in trees else None

    bases: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for (module, name), node in classes.items():
        bases[(module, name)] = []
        for base in node.bases:
            resolved = resolve(module, base)
            if isinstance(resolved, tuple):
                bases[(module, name)].append(resolved)

    def get_ancestors(key: tuple[str, str]) -> set[tuple[str, str]]:
        """Returns a package class and every package class it inherits from.

        :param key: The (module, name) tuple naming the class.
        :return: A set of (module, name) tuples, including the given class's.
        """
        ancestors: set[tuple[str, str]] = set()
        pending = [key]
        while pending:
            current = pending.pop()
            if current not in ancestors:
                ancestors.add(current)
                pending.extend(bases[current])
        return ancestors

    reachable: set[tuple[str, str]] = set()
    for module, name in lazy_callables.values():
        if (module, name) in classes:
            reachable |= get_ancestors((module, name))

    for key, node in classes.items():
        if key in reachable:
            continue
        properties = {
            item.name
            for ancestor in get_ancestors(key)
            for item in classes[ancestor].body
            if isinstance(item, ast.FunctionDef)
            and any(
                isinstance(decorator, ast.Name) and decorator.id == "property"
                for decorator in item.decorator_list
            )
        }
        members: list[tuple[int, str]] = []
        for item in node.body:
            if isinstance(item, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == "__slots__"
                for target in item.targets
            ):
                members.extend(
                    (item.lineno, constant.value)
                    for constant in ast.walk(item.value)
                    if isinstance(constant, ast.Constant)
                    and isinstance(constant.value, str)
                )
        members.extend(get_defined_names(node.body))
        for subnode in ast.walk(node):
            if (
                isinstance(subnode, ast.Attribute)
                and isinstance(subnode.ctx, ast.Store)
                and isinstance(subnode.value, ast.Name)
                and subnode.value.id == "self"
            ):
                members.append((subnode.lineno, subnode.attr))
        for line, member in members:
            if is_internal(member) and member[1:] not in properties:
                violations.append(
                    (
                        str(paths[key[0]]),
                        line,
                        f"{key[1]}.{member} is a member of an internal-only class "
                        f"and does not back a property, so it must not start with an "
                        f"underscore.",
                    )
                )

    return sorted(set(violations))


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
    """Checks every in scope file and package file and prints each violation.

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
    for path_name, line, message in find_naming_violations():
        exit_code = 1
        print(f"{path_name}:{line}: {message}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
