import os
import re
import sys
from datetime import datetime
from importlib.machinery import ModuleSpec
from pathlib import Path
from unittest.mock import MagicMock

# Add project root to sys.path so sphinx.ext.autodoc can import pterasoftware.
sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))

# Mock all runtime dependencies so autodoc can import pterasoftware without them
# installed. This is safe because the documented functions (save, load,
# set_up_logging) only use stdlib types in their signatures. autodoc_mock_imports
# only takes effect during autodoc's build phase, not when conf.py executes, so
# we also install a meta path finder to mock them in sys.modules before the
# pterasoftware imports below trigger the full package import chain.
autodoc_mock_imports = [
    "cmocean",
    "matplotlib",
    "numba",
    "numpy",
    "PySide6",
    "pyvista",
    "scipy",
    "tqdm",
    "webp",
]


class _MockModuleFinder:
    """Meta path finder that returns MagicMock for mocked top-level packages."""

    def __init__(self, mock_names):
        self._mock_names = set(mock_names)

    # noinspection PyUnusedLocal
    def find_spec(self, fullname, path, target=None):
        if fullname.split(".")[0] in self._mock_names:
            return ModuleSpec(fullname, self, is_package=True)
        return None

    # noinspection PyUnusedLocal
    # noinspection PyMethodMayBeStatic
    def create_module(self, spec):
        return MagicMock()

    def exec_module(self, module):
        pass


sys.meta_path.insert(0, _MockModuleFinder(autodoc_mock_imports))

# noinspection PyProtectedMember
from pterasoftware._logging import set_up_logging as _set_up_logging

# noinspection PyProtectedMember
# noinspection PyProtectedMember
from pterasoftware._serialization import load as _load
from pterasoftware._serialization import save as _save

# Override __module__ for public functions that live in private modules so
# autodoc renders them with the public package path (e.g., pterasoftware.save
# instead of pterasoftware._serialization.save).
_save.__module__ = "pterasoftware"
_load.__module__ = "pterasoftware"
_set_up_logging.__module__ = "pterasoftware"

# -- Project information -----------------------------------------------------

project = "PteraSoftware"
author = "Cameron Urban and contributors"
current_year = datetime.now().year
# noinspection PyShadowingBuiltins
copyright = f"{current_year}, {author}"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "autoapi.extension",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3


def _load_benchmark_host_info() -> dict[str, str]:
    """Read docs/_static/benchmarks/host.json into MyST substitutions.

    The benchmark publish workflow at PteraSoftwareBenchmarks drops host.json
    alongside the chart artifacts under docs/_static/benchmarks/. Fallback
    strings are returned when the file is absent (fresh checkout before the
    first benchmark publish has landed) so docs/website/performance.md still
    builds.
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
    "autosectionlabel.*",  # Duplicate labels from {include} directive pulling in source files
]

autosectionlabel_prefix_document = True

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
    # Exclude source markdown files to avoid duplicate autosectionlabel labels
    # (they are included into docs/website/ files, so we only want one copy processed)
    "../ANGLE_VECTORS_AND_TRANSFORMATIONS.md",
    "../AXES_POINTS_AND_FRAMES.md",
    "../CLASSES_AND_IMMUTABILITY.md",
    "../CODE_STYLE.md",
    "../TYPE_HINT_AND_DOCSTRING_STYLE.md",
    "../WRITING_STYLE.md",
    # Exclude brand files directory
    "Ptera Software Logo and Brand Files",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_attr_annotations = True

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_title = "PteraSoftware"
html_favicon = "favicon/favicon.ico"
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

# Furo: enable "Edit this page" with GitHub
html_theme_options = {
    "source_repository": "https://github.com/camUrban/PteraSoftware/",
    "source_branch": "main",
    "source_directory": "docs/website/",
    # Use black text logo in light mode (better contrast), normal logo in dark mode
    "light_logo": "Black_Text_Logo.png",
    "dark_logo": "Logo.png",
    # Hide the "PteraSoftware" text below logo (logo already contains the name)
    "sidebar_hide_name": True,
}

# For AutoAPI-generated pages, the default "Edit this page" points to a
# generated .rst path that doesn't exist in the repo. Override the URL to
# point to the corresponding Python source file in GitHub.
REPO_ROOT = Path(__file__).resolve().parents[2]


def _repo_rel_for_autoapi_page(pagename: str) -> str | None:
    """Return a repo-relative path for an AutoAPI page's corresponding source file."""
    if not pagename.startswith("api/"):
        return None
    rel = pagename[len("api/") :]
    if rel.endswith("/index"):
        rel = rel[: -len("/index")]
    py_path = REPO_ROOT / (rel.replace("/", os.sep) + ".py")
    if py_path.exists():
        target = py_path
    else:
        init_path = REPO_ROOT / rel.replace("/", os.sep) / "__init__.py"
        if init_path.exists():
            target = init_path
        else:
            return None
    return target.relative_to(REPO_ROOT).as_posix()


# noinspection PyUnusedLocal
def _html_page_context(app, pagename, templatename, context, doctree):
    repo_rel = _repo_rel_for_autoapi_page(pagename)
    if repo_rel:
        context["theme_source_edit_link"] = (
            f"https://github.com/camUrban/PteraSoftware/edit/main/{repo_rel}"
        )
        context["theme_source_view_link"] = (
            f"https://github.com/camUrban/PteraSoftware/blob/main/{repo_rel}?plain=true"
        )


def _rewrite_repo_root_links(app, docname, source):
    """Rewrite relative links in files included from the repo root.

    Files like CONTRIBUTING.md live at the repo root and use paths like
    ``docs/CODE_STYLE.md`` which are correct on GitHub. When Sphinx
    includes them via ``{include}``, those paths are resolved relative to
    ``docs/website/`` where the wrapper lives, so ``docs/CODE_STYLE.md``
    cannot be found. This handler replaces the wrapper's ``{include}``
    directive with the actual file content, rewriting ``docs/*.md`` paths
    to ``*.md`` so they resolve correctly in the Sphinx build.
    """
    contributing_path = REPO_ROOT / "CONTRIBUTING.md"
    if docname == "CONTRIBUTING" and contributing_path.exists():
        text = contributing_path.read_text()
        text = re.sub(r"\(docs/([A-Z_]+\.md)\)", r"(\1)", text)
        source[0] = text


def setup(app):
    app.connect("source-read", _rewrite_repo_root_links)
    app.connect("html-page-context", _html_page_context)

    # Copy extra assets to the site root after build
    # noinspection PyShadowingNames
    def _copy_extra_assets(app, exception):
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
    "*/movements/free_flight_*",
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


def autoapi_prepare_jinja_env(jinja_env):
    """Add custom Jinja filters for AutoAPI templates."""
    jinja_env.filters["first_paragraph"] = _first_paragraph
    jinja_env.filters["strip_init_boilerplate"] = _strip_init_boilerplate
    jinja_env.filters["first_sentence"] = _first_sentence
