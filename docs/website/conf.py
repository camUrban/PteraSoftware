import os
from datetime import datetime

# -- Project information -----------------------------------------------------

project = "PteraSoftware"
author = "Cameron Urban and contributors"
current_year = datetime.now().year
copyright = f"{current_year}, {author}"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "autoapi.extension",
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
html_static_path = ["_static", "favicon", "Black_Text_Logo.png", "Logo.png"]
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
from pathlib import Path

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


def _html_page_context(app, pagename, templatename, context, doctree):
    repo_rel = _repo_rel_for_autoapi_page(pagename)
    if repo_rel:
        context["theme_source_edit_link"] = (
            f"https://github.com/camUrban/PteraSoftware/edit/main/{repo_rel}"
        )
        context["theme_source_view_link"] = (
            f"https://github.com/camUrban/PteraSoftware/blob/main/{repo_rel}?plain=true"
        )


def setup(app):
    app.connect("html-page-context", _html_page_context)

    # Copy extra assets to the site root after build
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
]
autoapi_options = [
    "members",
    "show-module-summary",
    "show-inheritance",
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
    import re

    match = re.match(r"^(.*?\.)\s", text + " ")
    if match:
        return match.group(1)
    return text


def autoapi_prepare_jinja_env(jinja_env):
    """Add custom Jinja filters for AutoAPI templates."""
    jinja_env.filters["first_paragraph"] = _first_paragraph
    jinja_env.filters["strip_init_boilerplate"] = _strip_init_boilerplate
    jinja_env.filters["first_sentence"] = _first_sentence
