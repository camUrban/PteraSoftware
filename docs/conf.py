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

autosectionlabel_prefix_document = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "venv", ".venv"]

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
html_logo = "Logo.png"
html_favicon = "favicon/favicon.ico"
html_static_path = ["_static", "favicon"]
# Optionally also copy to site root (may be ignored by some builders)
html_extra_path = ["favicon"]

# Furo: enable "Edit this page" with GitHub
html_theme_options = {
    "source_repository": "https://github.com/camUrban/PteraSoftware/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# For AutoAPI-generated pages, the default "Edit this page" points to a
# generated .rst path that doesn't exist in the repo. Override the URL to
# point to the corresponding Python source file in GitHub.
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _repo_rel_for_autoapi_page(pagename: str) -> str | None:
    """Return a repo-relative path for an AutoAPI page's corresponding source file."""
    if not pagename.startswith("api/"):
        return None
    rel = pagename[len("api/"):]
    if rel.endswith("/index"):
        rel = rel[:-len("/index")]
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

    # Copy favicon assets to the site root so browsers that look for
    # e.g. "/favicon.ico" or "/apple-touch-icon.png" can find them.
    def _copy_favicons(app, exception):
        if exception is not None:
            return
        outdir = Path(app.outdir)
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

    app.connect("build-finished", _copy_favicons)


# -- AutoAPI configuration ---------------------------------------------------

# Parse the codebase directly (no imports) for API docs.
autoapi_type = "python"
autoapi_dirs = [os.path.abspath(os.path.join("..", "pterasoftware"))]
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
