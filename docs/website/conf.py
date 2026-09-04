import os
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

# Add project root to sys.path so sphinx.ext.autodoc can import pterasoftware.
sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))

# Mock all runtime dependencies so autodoc can import pterasoftware (via the
# autofunction directives for save, load, and set_up_logging) without them installed.
# This is safe because the documented functions only use stdlib types in their
# signatures.
autodoc_mock_imports = [
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
current_year = datetime.now().year
# noinspection PyShadowingBuiltins
copyright = f"{current_year}, {author}"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_nb",
    "autoapi.extension",
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
    "toc.not_included",  # Template files not in toctree (expected)
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

# AutoAPI renders docstrings as reStructuredText, where a vector magnitude written with
# bars (for example, "|g_E|") parses as a substitution reference. Define those tokens so
# the reference resolves to the literal barred text instead of emitting an "Undefined
# substitution referenced" build error, without annotating the docstrings themselves.
rst_prolog = r"""
.. |g_E| replace:: \|g_E\|
"""

# Use README as the landing page (instead of index)
root_doc = "README"

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "venv",
    ".venv",
    # Exclude autoapi templates from being parsed as RST documents
    "_autoapi_templates",
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

REPO_ROOT = Path(__file__).resolve().parents[2]

# Sphinx only renders documents that live inside the source directory, so copy the
# tutorial notebooks (and the images they embed) from the repo root's tutorials/
# directory into docs/website/tutorials/. The copies are gitignored, and tutorials/
# stays the single source of truth.
_tutorials_source = REPO_ROOT / "tutorials"
_tutorials_target = Path(__file__).resolve().parent / "tutorials"
_tutorials_target.mkdir(exist_ok=True)
for _tutorial_file in sorted(_tutorials_source.iterdir()):
    if _tutorial_file.suffix in {".ipynb", ".png", ".webp"}:
        shutil.copy2(_tutorial_file, _tutorials_target / _tutorial_file.name)


def _rewrite_repo_root_links(app: Any, docname: str, source: list[str]) -> None:
    """Rewrite relative links in files included from the repo root.

    Files like CONTRIBUTING.md live at the repo root and use paths like
    ``docs/CODE_STYLE.md`` which are correct on GitHub. When Sphinx includes them via
    ``{include}``, those paths are resolved relative to ``docs/website/`` where the
    wrapper lives, so ``docs/CODE_STYLE.md`` cannot be found. This handler replaces the
    wrapper's ``{include}`` directive with the actual file content, rewriting
    ``docs/*.md`` paths to ``*.md`` so they resolve correctly in the Sphinx build.
    """
    contributing_path = REPO_ROOT / "CONTRIBUTING.md"
    if docname == "CONTRIBUTING" and contributing_path.exists():
        text = contributing_path.read_text()
        text = re.sub(r"\(docs/([A-Z_]+\.md)\)", r"(\1)", text)
        source[0] = text


def setup(app: Any) -> None:
    app.connect("source-read", _rewrite_repo_root_links)

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


# -- AutoAPI configuration ---------------------------------------------------

# Parse the codebase directly (no imports) for API docs.
autoapi_type = "python"
autoapi_dirs = [os.path.abspath(os.path.join("..", "..", "pterasoftware"))]
autoapi_root = "api"
autoapi_add_toctree_entry = False
autoapi_keep_files = False
autoapi_member_order = "bysource"
autoapi_ignore = [
    "*/ui_resources/*",
    "*/airfoils/*",
    "*/models/*",
]
autoapi_options = [
    "members",
    "show-module-summary",
    "show-inheritance",
    "inherited-members",
]
autoapi_template_dir = "_autoapi_templates"

# Include __init__ docstrings (which contain parameter descriptions) with class docs
autoapi_python_class_content = "both"


def _first_paragraph(docstring: str) -> str:
    """Extract the first paragraph from a docstring (up to the first blank line)."""
    if not docstring:
        return ""
    lines = docstring.split("\n")
    result_lines = []
    for line in lines:
        if line.strip() == "":
            break
        result_lines.append(line)
    return "\n".join(result_lines)


def _strip_init_boilerplate(docstring: str) -> str:
    """Remove 'The initialization method.' line from docstrings."""
    if not docstring:
        return ""
    lines = docstring.split("\n")
    filtered_lines = []
    for line in lines:
        if line.strip() == "The initialization method.":
            continue
        filtered_lines.append(line)
    # Remove leading blank lines that may result from filtering
    while filtered_lines and filtered_lines[0].strip() == "":
        filtered_lines.pop(0)
    return "\n".join(filtered_lines)


def _first_sentence(docstring: str) -> str:
    """Extract the first sentence from a docstring, joining hard-wrapped lines."""
    if not docstring:
        return ""
    # Join all lines with spaces, then find the first sentence
    lines = docstring.split("\n")
    # Collect lines until we hit a blank line (paragraph break)
    paragraph_lines = []
    for line in lines:
        if line.strip() == "":
            break
        paragraph_lines.append(line.strip())
    # Join into one string
    text = " ".join(paragraph_lines)
    # Find the first sentence (ends with period followed by space or end)
    match = re.match(r"^(.*?\.)\s", text + " ")
    if match:
        return match.group(1)
    return text


def autoapi_prepare_jinja_env(jinja_env: Any) -> None:
    """Add custom Jinja filters for AutoAPI templates."""
    jinja_env.filters["first_paragraph"] = _first_paragraph
    jinja_env.filters["strip_init_boilerplate"] = _strip_init_boilerplate
    jinja_env.filters["first_sentence"] = _first_sentence
