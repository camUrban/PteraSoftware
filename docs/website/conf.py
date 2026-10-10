import ast
import os
import re
import shutil
import sys
from collections.abc import Iterator
from datetime import datetime
from pathlib import Path
from typing import Any

# Add project root to sys.path so sphinx.ext.autodoc can import pterasoftware.
sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))

REPO_ROOT = Path(__file__).resolve().parents[2]
_TUTORIALS_SOURCE = REPO_ROOT / "tutorials"
_TUTORIALS_TARGET = Path(__file__).resolve().parent / "tutorials"
_PACKAGE_DIR = REPO_ROOT / "pterasoftware"
_API_TARGET = Path(__file__).resolve().parent / "api"

# Parameter annotations to show in place of the source annotation, keyed by the fully
# qualified class (for constructor parameters) or method where it is defined, and then
# by parameter name. These are signatures whose source annotation names an internal
# class, which the API reference does not document, and is wider than what the
# implementation accepts: the hook methods are widened to the shared parent solver type
# because an override cannot narrow a parameter type, and the base unsteady solver's
# constructor is widened to the shared parent problem type so the derived solvers can
# pass their own problem types through it. The docs show the type that actually works,
# and contributors can read the source for the formal contract.
_ANNOTATION_OVERRIDES = {
    "pterasoftware._unsteady_ring_vortex_lattice_method.UnsteadyRingVortexLatticeMethodSolver": {
        "unsteady_problem": "~pterasoftware.UnsteadyProblem",
    },
    "pterasoftware._problems.FreeFlightUnsteadyProblem.initialize_next_problem": {
        "solver": "~pterasoftware.FreeFlightUnsteadyRingVortexLatticeMethodSolver",
    },
    "pterasoftware._problems.AeroelasticUnsteadyProblem.initialize_next_problem": {
        "solver": "~pterasoftware.AeroelasticUnsteadyRingVortexLatticeMethodSolver",
    },
}

# Mock all runtime dependencies so autodoc can import pterasoftware without them
# installed. Each documented object's annotations still render as written in the source
# (for example, np.ndarray), since autodoc falls back to the source text for the types
# it cannot resolve.
autodoc_mock_imports = [
    "fontTools",
    "matplotlib",
    "mujoco",
    "numba",
    "numpy",
    "pyvista",
    "scipy",
    "threadpoolctl",
    "tqdm",
    "webp",
]

# -- Project information -----------------------------------------------------

project = "PteraSoftware"
author = "Cameron Urban and contributors"
CURRENT_YEAR = datetime.now().year
# noinspection PyShadowingBuiltins
copyright = f"{CURRENT_YEAR}, {author}"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx_design",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

# Render the tutorial notebooks from their committed outputs instead of executing them.
# The docs build does not install pterasoftware or its runtime dependencies (see
# autodoc_mock_imports above), so it cannot run the notebooks.
nb_execution_mode = "off"


def _load_benchmark_host_info() -> dict[str, str]:
    """Read docs/_static/benchmarks/host.json into MyST substitutions.

    The benchmark publish workflow at PteraSoftwareBenchmarks drops host.json alongside
    the chart artifacts under docs/_static/benchmarks/. Fallback strings are returned
    when the file is absent (fresh checkout before the first benchmark publish has
    landed) so docs/website/performance.md still builds.
    """
    import json

    path = Path(__file__).parent.parent / "_static" / "benchmarks" / "host.json"
    fallback = {
        "host_os": "Pending first benchmark publish",
        "host_cpu": "Pending first benchmark publish",
        "host_cores": "n/a",
        "host_governor": "n/a",
        "host_memory_mb": "n/a",
        "host_swappiness": "n/a",
        "host_thp": "n/a",
        "host_gpu": "Pending first benchmark publish",
        "host_gpu_driver": "n/a",
        "host_cuda": "n/a",
        "host_storage": "Pending first benchmark publish",
    }
    if not path.exists():
        return fallback
    payload = json.loads(path.read_text())
    os_info = payload.get("os", {})
    cpu = payload.get("cpu", {})
    memory = payload.get("memory", {})
    gpu = payload.get("gpu", {})
    storage = payload.get("storage", {})
    os_label = (
        f"{os_info.get('name', '')} {os_info.get('version', '')}".strip() or "unknown"
    )
    return {
        "host_os": os_label,
        "host_cpu": cpu.get("model", "unknown"),
        "host_cores": str(cpu.get("logical_cores", "unknown")),
        "host_governor": cpu.get("governor", "unknown"),
        "host_memory_mb": str(memory.get("total_mb", "unknown")),
        "host_swappiness": str(memory.get("swappiness", "unknown")),
        "host_thp": memory.get("transparent_hugepages", "unknown"),
        "host_gpu": gpu.get("model", "unknown"),
        "host_gpu_driver": gpu.get("driver_version", "unknown"),
        "host_cuda": gpu.get("cuda_version", "unknown"),
        "host_storage": storage.get("model", "unknown"),
    }


myst_substitutions = _load_benchmark_host_info()

# Suppress warnings that are informational or unavoidable
suppress_warnings = [
    # Duplicate labels from {include} directive pulling in source files.
    "autosectionlabel.*",
]

autosectionlabel_prefix_document = True

# Render every Python signature with one parameter per line. Any signature longer than
# this threshold (in characters) wraps so that each parameter sits on its own indented
# line with a trailing comma and the closing parenthesis on its own line. A threshold of
# one forces this multi-line layout for every signature that takes at least one
# parameter, which keeps long class and function parameter lists readable instead of
# running together on one line.
python_maximum_signature_line_length = 1

# Keep straight quotes in the rendered site. Docstrings write str values and code
# examples in plain prose with double quotes (see
# docs/TYPE_HINT_AND_DOCSTRING_STYLE.md), and the default conversion to curly quotes
# would turn examples like set_up_logging(level="Info") into text that cannot be copied
# back into Python.
smartquotes = False

# Use README as the landing page (instead of index)
root_doc = "README"

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    # The notebook extension writes its executed notebooks into a folder next to the
    # build output folder, so a local build into docs/website/_build/ leaves this folder
    # in the source directory, where the next build would read it as stray documents.
    "jupyter_execute",
    "Thumbs.db",
    ".DS_Store",
    "venv",
    ".venv",
    # Exclude directories in parent docs/ folder
    "../katz_plotkin_13_12",
    "../lambert_2015_2_3__2_4",
    "../RUNNING_TESTS_AND_TYPE_CHECKS.md",
    "../examples_expected_output",
    # Exclude source markdown files to avoid duplicate autosectionlabel labels (they are
    # included into docs/website/ files, so we only want one copy processed).
    "../ANGLE_VECTORS_AND_TRANSFORMATIONS.md",
    "../AXES_POINTS_AND_FRAMES.md",
    "../CLASSES_AND_IMMUTABILITY.md",
    "../CODE_STYLE.md",
    "../MUJOCO_CONVENTIONS.md",
    "../STRONG_COUPLING.md",
    "../TYPE_HINT_AND_DOCSTRING_STYLE.md",
    "../WRITING_STYLE.md",
    # Exclude brand files directory
    "Ptera Software Logo and Brand Files",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_attr_annotations = True

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_title = "PteraSoftware"
html_favicon = "favicon/favicon.ico"
# Drop the "View this page" source link (and the _sources/*.txt dump it points to) so no
# page exposes a link to its underlying source.
html_show_sourcelink = False
html_copy_source = False
html_static_path = ["_static", "../_static", "Black_Text_Logo.png", "Logo.png"]
# Optionally also copy to site root (may be ignored by some builders)
html_extra_path = ["favicon"]

# Custom CSS with Ptera brand styling
html_css_files = [
    "custom.css",
]

# Custom JavaScript
html_js_files = [
    "custom.js",
]

html_theme_options = {
    # Use black text logo in light mode (better contrast), normal logo in dark mode
    "light_logo": "Black_Text_Logo.png",
    "dark_logo": "Logo.png",
    # Hide the "PteraSoftware" text below logo (logo already contains the name)
    "sidebar_hide_name": True,
}

# Sphinx only renders documents that live inside the source directory, so copy the
# tutorial notebooks (and the images they embed) from the repo root's tutorials/
# directory into docs/website/tutorials/. The copies are gitignored, and tutorials/
# stays the single source of truth.
_TUTORIALS_TARGET.mkdir(exist_ok=True)
for _tutorial_file in sorted(_TUTORIALS_SOURCE.iterdir()):
    if _tutorial_file.suffix in {".ipynb", ".png", ".webp"}:
        shutil.copy2(_tutorial_file, _TUTORIALS_TARGET / _tutorial_file.name)


def _read_package_tables() -> tuple[list[str], dict[str, tuple[str, str]]]:
    """Read __all__ and the lazy callable table from the package's __init__.py.

    The tables are read from the source rather than by importing the package, so this
    does not depend on the package's runtime dependencies. The lazy callable table maps
    each public name to the internal module that defines it and its name there.
    """
    tree = ast.parse((_PACKAGE_DIR / "__init__.py").read_text())
    tables: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in {
                "__all__",
                "_LAZY_CALLABLES",
            }:
                tables[target.id] = ast.literal_eval(node.value)
    return tables["__all__"], tables["_LAZY_CALLABLES"]


def _is_class(module: str, name: str) -> bool:
    """Report whether a module in the package defines the given name as a class."""
    path = REPO_ROOT.joinpath(*module.split("."))
    path = path / "__init__.py" if path.is_dir() else path.with_suffix(".py")
    return any(
        isinstance(node, ast.ClassDef) and node.name == name
        for node in ast.parse(path.read_text()).body
    )


_PUBLIC_NAMES, _PUBLIC_DEFINITIONS = _read_package_tables()

# The public name of each public object, keyed by the internal module that defines it
# and its name there.
_PUBLIC_NAMES_BY_DEFINITION = {
    definition: name for name, definition in _PUBLIC_DEFINITIONS.items()
}

# Write one API reference page per public name into docs/website/api/, which is
# gitignored, so that the reference always matches __all__. api.md lists the pages, and
# the build warns about any page it lists that does not exist and any page it does not
# list.
shutil.rmtree(_API_TARGET, ignore_errors=True)
_API_TARGET.mkdir()
for _name in _PUBLIC_NAMES:
    if _is_class(*_PUBLIC_DEFINITIONS[_name]):
        _title, _directive = f"`pterasoftware.{_name}`", "autoclass"
    else:
        _title, _directive = f"`pterasoftware.{_name}()`", "autofunction"
    (_API_TARGET / f"{_name}.md").write_text(
        f"# {_title}\n\n```{{eval-rst}}\n.. {_directive}:: pterasoftware.{_name}\n```\n"
    )


def _rewrite_repo_root_links(app: Any, docname: str, source: list[str]) -> None:
    """Rewrite relative links in files included from the repo root.

    Files like CONTRIBUTING.md live at the repo root and use paths like
    docs/CODE_STYLE.md which are correct on GitHub. When Sphinx includes them via
    {include}, those paths are resolved relative to docs/website/ where the wrapper
    lives, so docs/CODE_STYLE.md cannot be found. This handler replaces the wrapper's
    {include} directive with the actual file content, stripping the leading docs/ from
    each such path so they resolve correctly in the Sphinx build.
    """
    contributing_path = REPO_ROOT / "CONTRIBUTING.md"
    if docname == "CONTRIBUTING" and contributing_path.exists():
        text = contributing_path.read_text()
        text = re.sub(r"\(docs/([A-Z_]+\.md)\)", r"(\1)", text)
        source[0] = text


def _is_deprecation_warning(node: ast.AST) -> bool:
    """Reports whether a node is a statement calling warnings.warn with
    DeprecationWarning."""
    if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
        return False
    call = node.value
    if not (
        isinstance(call.func, ast.Attribute)
        and call.func.attr == "warn"
        and isinstance(call.func.value, ast.Name)
        and call.func.value.id == "warnings"
    ):
        return False
    category: ast.expr | None = call.args[1] if len(call.args) > 1 else None
    for keyword in call.keywords:
        if keyword.arg == "category":
            category = keyword.value
    return isinstance(category, ast.Name) and category.id == "DeprecationWarning"


def _iter_functions(package_dir: Path) -> Iterator[tuple[str, ast.FunctionDef]]:
    """Yield each module-level function and method in a package from its source.

    Each function is yielded with its fully qualified name, which is the module's dotted
    name followed by the function's qualified name, as in the function's own __module__
    and __qualname__. A constructor is named by its class instead, since the API
    reference renders its parameters on the class signature.
    """
    for path in sorted(package_dir.rglob("*.py")):
        parts = path.relative_to(package_dir).with_suffix("").parts
        if parts[-1] == "__init__":
            parts = parts[:-1]
        module = ".".join((package_dir.name, *parts))
        for node in ast.parse(path.read_text()).body:
            if isinstance(node, ast.FunctionDef):
                yield f"{module}.{node.name}", node
            elif isinstance(node, ast.ClassDef):
                class_id = f"{module}.{node.name}"
                for member in node.body:
                    if isinstance(member, ast.FunctionDef):
                        if member.name == "__init__":
                            yield class_id, member
                        else:
                            yield f"{class_id}.{member.name}", member


def _get_parameter_names(function: ast.FunctionDef) -> set[str]:
    """Return the names of a function's parameters, other than * and ** parameters."""
    return {
        argument.arg
        for argument in (
            *function.args.posonlyargs,
            *function.args.args,
            *function.args.kwonlyargs,
        )
    }


def _find_deprecations(package_dir: Path) -> tuple[set[str], dict[str, list[str]]]:
    """Find the deprecated members and parameters in a package from its source.

    A deprecation is whatever emits a DeprecationWarning, so this reads the package's
    source rather than any docstring marker. A function, method, or property is
    deprecated when its body calls warnings.warn with DeprecationWarning as a top-level
    statement. A parameter is deprecated when such a call sits inside a top-level if
    statement whose test names that parameter. The first return value holds the fully
    qualified names of the deprecated members. The second maps the fully qualified name
    of each function or method with deprecated parameters to those parameters' names,
    with a constructor's parameters keyed by its class, since the API reference renders
    them on the class signature.
    """
    members: set[str] = set()
    parameters: dict[str, list[str]] = {}
    for function_id, function in _iter_functions(package_dir):
        if any(_is_deprecation_warning(statement) for statement in function.body):
            members.add(function_id)
        parameter_names = _get_parameter_names(function)
        for statement in function.body:
            if not isinstance(statement, ast.If):
                continue
            if not any(_is_deprecation_warning(node) for node in ast.walk(statement)):
                continue
            for node in ast.walk(statement.test):
                if isinstance(node, ast.Name) and node.id in parameter_names:
                    parameters.setdefault(function_id, []).append(node.id)
    return members, parameters


_DEPRECATED_MEMBERS, _DEPRECATED_PARAMETERS = _find_deprecations(_PACKAGE_DIR)

# Fail the build on an annotation override whose function or parameter no longer exists,
# so a rename cannot silently drop an override.
_PARAMETER_NAMES = {
    function_id: _get_parameter_names(function)
    for function_id, function in _iter_functions(_PACKAGE_DIR)
}
for _target_id, _overrides in _ANNOTATION_OVERRIDES.items():
    for _parameter in _overrides:
        if _parameter not in _PARAMETER_NAMES.get(_target_id, set()):
            raise ValueError(
                f"Annotation override target {_target_id} has no parameter "
                f"{_parameter}."
            )


def _get_object_id(obj: Any) -> str:
    """Return the fully qualified name of an object that autodoc documents.

    The name is the object's module followed by its qualified name, which matches the
    names that _iter_functions gives. A property is named by its getter.
    """
    target = obj.fget if isinstance(obj, property) else obj
    return f"{getattr(target, '__module__', '')}.{getattr(target, '__qualname__', '')}"


def _skip_deprecated_members(
    app: Any, what: str, name: str, obj: Any, skip: bool, options: Any
) -> bool | None:
    """Keep deprecated members out of the API reference.

    Returns True to skip a deprecated member and None to leave autodoc's own decision in
    place for everything else.
    """
    if _get_object_id(obj) in _DEPRECATED_MEMBERS:
        return True
    return None


def _parameter_end(signature: str, start: int) -> int:
    """Find where the parameter starting at the given index ends in a signature.

    The end is the next top-level comma or closing parenthesis, ignoring the commas
    inside subscripted annotations such as list[dict[str, int]].
    """
    depth = 0
    end = start
    while end < len(signature):
        character = signature[end]
        if character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
        elif character in ",)" and depth == 0:
            break
        end += 1
    return end


def _override_annotation(signature: str, parameter: str, annotation: str) -> str:
    """Replace one parameter's annotation in a signature, keeping any default."""
    match = re.search(rf"(?<![\w]){parameter}: ", signature)
    if match is None:
        raise ValueError(f"Parameter {parameter} not found in signature {signature}.")
    start = match.end()
    end = _parameter_end(signature, start)
    default = signature[start:end].find(" = ")
    if default != -1:
        end = start + default
    return signature[:start] + annotation + signature[end:]


def _remove_parameter(signature: str, parameter: str) -> str:
    """Remove one parameter, with its annotation and default, from a signature."""
    match = re.search(rf"(?<![\w]){parameter}(?=[:=,)])", signature)
    if match is None:
        raise ValueError(f"Parameter {parameter} not found in signature {signature}.")
    start = match.start()
    end = _parameter_end(signature, start)
    if signature.startswith(", ", end):
        end += 2
    elif signature.startswith(", ", start - 2):
        start -= 2
    return signature[:start] + signature[end:]


# This matches a dotted path into the package as autodoc writes it into a signature,
# such as ~pterasoftware._geometry.wing.Wing.
_PACKAGE_PATH_PATTERN = re.compile(r"~?pterasoftware(?:\.\w+)+")


def _flatten_public_paths(text: str) -> str:
    """Replace each path to a public object's internal definition with its public name.

    This makes the reference link each such annotation to the object's page. Paths to
    internal objects are left as they are.
    """

    def flatten(match: re.Match[str]) -> str:
        module, _, name = match.group().lstrip("~").rpartition(".")
        public_name = _PUBLIC_NAMES_BY_DEFINITION.get((module, name))
        if public_name is None:
            return match.group()
        return f"~pterasoftware.{public_name}"

    return _PACKAGE_PATH_PATTERN.sub(flatten, text)


def _rewrite_signature(
    app: Any,
    what: str,
    name: str,
    obj: Any,
    options: Any,
    signature: str | None,
    return_annotation: str | None,
) -> tuple[str | None, str | None]:
    """Show public names in a signature, apply its annotation overrides, and remove its
    deprecated parameters."""
    if signature is not None:
        object_id = _get_object_id(obj)
        signature = _flatten_public_paths(signature)
        for parameter, annotation in _ANNOTATION_OVERRIDES.get(object_id, {}).items():
            signature = _override_annotation(signature, parameter, annotation)
        for parameter in _DEPRECATED_PARAMETERS.get(object_id, []):
            signature = _remove_parameter(signature, parameter)
    if return_annotation:
        return_annotation = _flatten_public_paths(return_annotation)
    return signature, return_annotation


def _process_docstring(
    app: Any, what: str, name: str, obj: Any, options: Any, lines: list[str]
) -> None:
    """Prepare a docstring for the API reference.

    This removes the "The initialization method." line that opens each constructor
    docstring, which autodoc shows below the class docstring. It also removes the
    :param: field of each deprecated parameter, together with its continuation lines,
        which are indented deeper than the field line.
    """
    kept_lines: list[str] = []
    field_indent: int | None = None
    deprecated_parameters = _DEPRECATED_PARAMETERS.get(_get_object_id(obj), [])
    for line in lines:
        indent = len(line) - len(line.lstrip())
        if field_indent is not None:
            if line.strip() and indent > field_indent:
                continue
            field_indent = None
        if any(
            line.lstrip().startswith(f":param {parameter}:")
            for parameter in deprecated_parameters
        ):
            field_indent = indent
            continue
        if line.strip() == "The initialization method.":
            continue
        kept_lines.append(line)
    lines[:] = kept_lines


def setup(app: Any) -> None:
    app.connect("source-read", _rewrite_repo_root_links)
    app.connect("autodoc-process-signature", _rewrite_signature)
    app.connect("autodoc-process-docstring", _process_docstring)
    app.connect("autodoc-skip-member", _skip_deprecated_members)

    # Copy extra assets to the site root after build
    # noinspection PyShadowingNames
    def _copy_extra_assets(app: Any, exception: Exception | None) -> None:
        if exception is not None:
            return
        outdir = Path(app.outdir)

        # Copy favicon assets so browsers can find them at root
        src = Path(__file__).parent / "favicon"
        names = [
            "favicon.ico",
            "favicon.svg",
            "apple-touch-icon.png",
            "favicon-96x96.png",
            "site.webmanifest",
            "web-app-manifest-192x192.png",
            "web-app-manifest-512x512.png",
        ]
        for n in names:
            p = src / n
            if p.exists():
                (outdir / n).write_bytes(p.read_bytes())

        # Create index.html redirect to README.html
        index_html = outdir / "index.html"
        index_html.write_text(
            "<!DOCTYPE html>\n"
            "<html>\n"
            "<head>\n"
            '    <meta charset="utf-8">\n'
            "    <title>Redirecting to Ptera Software</title>\n"
            '    <meta http-equiv="refresh" content="0; url=README.html">\n'
            '    <link rel="canonical" href="README.html">\n'
            "</head>\n"
            "<body>\n"
            '    <p>Redirecting to <a href="README.html">Ptera Software</a>...</p>\n'
            '    <script>window.location.href = "README.html";</script>\n'
            "</body>\n"
            "</html>\n"
        )

    app.connect("build-finished", _copy_extra_assets)


# -- Autodoc configuration ---------------------------------------------------

# Document each class's public members, including those it inherits from its internal
# parents, in the order the source defines them.
autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "member-order": "bysource",
}

# Include __init__ docstrings (which contain parameter descriptions) with class docs
autoclass_content = "both"
