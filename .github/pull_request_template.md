# Description

Provide a concise description of the changes in this pull request. Appropriate titles for PRs use sentence-case and no trailing period. Also, add applicable labels and assign the PR to yourself. If using the `bug` or `feature` label, prefix your title with [BUG] or [FEATURE]. If a first-time contributor, welcome! If you'd like, feel free to add yourself to the end of the list of contributors in the [README](https://github.com/camUrban/PteraSoftware/blob/main/README.md), following one of the styles shown.

## Motivation

Explain why you are making these changes and what problem they solve.

## Relevant Issues

Link any related issues using GitHub's syntax. For bugs, write, "Fixes #<issue number>." For features, "Closes #<issue number>." If there are no related issues, write "None."

## Changes

* List the changes you made in bullet points.

## Dependency Updates

List any new dependencies (including dev dependencies) added in this PR. Also list any updates to required dependency versions. If no new dependencies have been added or version requirements changed, write "None."

## Change Magnitude

Identify the option that best describes the impact of your change, then delete the other two lines and this sentence.

**Major**: Large change that adds significant new functionality, changes existing behavior, or may affect many parts of the codebase.

**Moderate**: Medium-sized change that adds or modifies a feature without large-scale impact.

**Minor**: Small change such as a bug fix, small enhancement, or documentation update.

## Checklist (check each item when completed or not applicable)

* [ ] I am familiar with the current [contribution guidelines](https://github.com/camUrban/PteraSoftware/blob/main/CONTRIBUTING.md).
* [ ] PR description links all relevant issues and follows this template.
* [ ] My branch is based on `main` and is up to date with the upstream `main` branch.
* [ ] All calculations use S.I. units.
* [ ] Code is formatted with [black](https://github.com/psf/black) (line length = 88).
* [ ] Code is well documented with block comments where appropriate.
* [ ] Any external code, algorithms, or equations used have been cited in comments or docstrings.
* [ ] All new modules, classes, functions, and methods have docstrings in [reStructuredText format](https://realpython.com/documenting-python-code/), and are formatted using [docformatter](https://github.com/PyCQA/docformatter) (`--in-place --black`). See the [style guide for type hints and docstrings](https://github.com/camUrban/PteraSoftware/blob/main/docs/TYPE_HINT_AND_DOCSTRING_STYLE.md) for more details.
* [ ] All new classes, functions, and methods in the `pterasoftware` package use type hints. See the [style guide for type hints and docstrings](https://github.com/camUrban/PteraSoftware/blob/main/docs/TYPE_HINT_AND_DOCSTRING_STYLE.md) for more details.
* [ ] If any major functionality was added or significantly changed, I have added or updated tests in the `tests` package.
* [ ] Code locally passes all tests in the `tests` package.
* [ ] This PR passes the ReadTheDocs build check (this runs automatically with the other workflows).
* [ ] This PR passes the `ascii-only`, `black`, `codespell`, `docformatter`, `isort`, and `pre-commit-hooks` GitHub actions.
* [ ] This PR passes the `mypy` GitHub action.
* [ ] This PR passes all the `tests` GitHub actions.
