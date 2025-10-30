# Description
Provide a concise description of the changes in this pull request.

## Motivation
Explain why you are making these changes and what problem they solve.

## Relevant Issues
Link any related issues using GitHub's syntax, e.g., `Fixes #123` or `Closes #456`.

## Changes
List the changes you made in bullet points.

## New Dependencies
List any new dependencies (including dev dependencies) added in this PR.

## Change Magnitude
Select the option that best describes the size and impact of your change:

- **Major**: Large change that adds significant new functionality, changes existing 
behavior, or may affect many parts of the codebase.
- **Moderate**: Medium-sized change that adds or modifies a feature without large-scale 
impact.
- **Minor**: Small change such as a bug fix, small enhancement, or documentation update.

---

# Checklist

- [ ] I have created or claimed an issue for this work as described in 
[Contributing Code](https://github.com/camUrban/PteraSoftware/blob/main/CONTRIBUTING.md#contributing-code).
- [ ] My branch is based on `main` and is up to date with the upstream `main` branch.
- [ ] All calculations use S.I. units.
- [ ] Code is formatted with [black](https://github.com/psf/black) (line length = 88).
- [ ] Code is well documented with block comments where appropriate.
- [ ] Any external code, algorithms, or equations used have been cited in comments or 
docstrings.
- [ ] All new modules, classes, functions, and methods have docstrings in 
[reStructuredText format](https://realpython.com/documenting-python-code/), and are 
formatted using [docformatter](https://github.com/PyCQA/docformatter) 
(`--in-place --black`).
- [ ] All new classes, functions, and methods use type hints.
- [ ] If any major functionality was added or significantly changed, I have added or 
updated tests in the `tests` package.
- [ ] Code locally passes all tests in the `tests` package.
- [ ] After pushing, PR passes all automated checks (`codespell`, `black`, and `tests`).
- [ ] PR description links all relevant issues and follows this template.

---