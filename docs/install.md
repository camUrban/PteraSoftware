# Installation

PteraSoftware targets Python 3.13 for runtime usage. For building this
documentation (including on Read the Docs), only Sphinx and doc build
dependencies are installed — the main package and heavy scientific
dependencies are not required.

## User installation (package)

```shell
pip install pterasoftware
```

## From source (developers)

```shell
git clone https://github.com/camUrban/PteraSoftware.git
cd PteraSoftware
python -m venv .venv
. .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
```

## Build docs locally (optional)

```shell
python -m venv .venv
. .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r docs/requirements_docs.txt
sphinx-build -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html` in a browser.
