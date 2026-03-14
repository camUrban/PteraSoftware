# Contributing Guidelines

We are excited that you are interested in contributing to **Ptera Software**! This guide will help you get started. If you have any questions not answered here, please [open a new discussion](https://github.com/camUrban/PteraSoftware/discussions) and we will respond.

---

## Before Contributing

Please review the following documents before making contributions. These documents are also available on the [Ptera Software documentation website](https://pterasoftware.readthedocs.io/).

1. [README](README.md)
2. [Code of Conduct](CODE_OF_CONDUCT.md)
3. [Security Policy](SECURITY.md)
4. [License](LICENSE.md)

---

## Ways to Contribute

There are three main ways you can contribute:

1. [**Report a bug**](#reporting-a-bug)
    - Identify and document issues that prevent Ptera Software from working as expected.
    - This includes errors, crashes, incorrect outputs, or unexpected behavior.
2. [**Request a feature**](#requesting-a-feature)
    - Suggest new features, improvements to existing functionality, or usability changes.
    - This can be small quality-of-life improvements or larger feature proposals.
3. [**Contribute code**](#contributing-code)
    - Submit changes that add new features, fix bugs, improve performance, or enhance documentation.
    - Code contributions can address your own ideas or work on existing open issues.
    - If you find an issue labeled `good first issue` that you want to work on, **comment on the issue to claim it** before starting work. This prevents duplicate efforts.

---

### Reporting a Bug

**If the bug is a security vulnerability, do not post it as a public issue**. Follow the [security policy](SECURITY.md) instead.

For all other bugs:

1. Search the [issues page](https://github.com/camUrban/PteraSoftware/issues) to check if it's already been reported.
    - If it exists, you can comment or add an emoji reaction to indicate that you are also affected.
2. If it has not been reported, open a new issue using the [bug report template](https://github.com/camUrban/PteraSoftware/issues/new?template=bug_report.md).
3. Add the `bug` label plus any others that apply.

---

### Requesting a Feature

For feature requests:

1. Search the [issues page](https://github.com/camUrban/PteraSoftware/issues) to avoid duplicates.
    - If it exists, comment or react to indicate interest.
2. If it has not been requested, open a new issue using the [feature request template](https://github.com/camUrban/PteraSoftware/issues/new?template=feature_request.md).
3. Add the `feature` label plus any others that apply.

---

### Contributing Code

Ptera Software now uses GitHub Flow to manage code contributions. If this is new to you, it's a good idea to read through [this guide](https://docs.github.com/en/get-started/using-github/github-flow) first. Once you understand the process, here's how to implement it:

1. **Choose what to work on**
    - Look for issues labeled `good first issue`.
        - If you want to work on one, check that no one else has already commented claiming it. If unclaimed, comment on the issue to claim it and a maintainer will then "assign" the issue to you.
        - If you want to work on a claimed issue that hasn't been updated in a while, write a comment asking if the user who originally claimed it is still actively working on it.
    - If you have your own idea, search the issues to ensure it hasn't already been proposed.
        - If you can't find anything, you could open a new issue describing your change and state that you will be implementing it. This step isn't strictly required, but it is recommended for first-time contributors.

   #### Coordination on In-Progress and Draft Work

   Some issues and pull requests (PRs) may represent work that is actively under design or refinement.

   If a PR exists, or if an issue has been assigned to or claimed by someone else, please do **not** start parallel implementation work without first checking in via a comment. A short message like "Are you still actively working on this?" or "Is this a good time to help with X?” is usually sufficient.

   This helps avoid duplicated effort and ensures that contributions align with the current design direction.
2. **Set up your local environment**
    - Fork the repository on GitHub to your own account.
    - Clone your fork to your local machine:
    ```shell
    git clone https://github.com/<your-username>/PteraSoftware.git
    cd PteraSoftware
    ```
    - Add the main repository as a remote named `upstream` so you can keep your fork up to date:
    ```shell
    git remote add upstream https://github.com/camUrban/PteraSoftware.git
    ```
    - Create a virtual environment and install dependencies for development. Note that Ptera Software requires Python 3.11, 3.12, or 3.13 (3.13 is recommended):
    ```shell
    python3 -m venv .venv # On Windows use python instead of python3
    source .venv/bin/activate # On Windows use .venv\Scripts\activate
    python3 -m pip install --upgrade pip # On Windows use python instead of python3
    pip install -r requirements.txt # Install dependencies for running simulations
    pip install -r requirements_dev.txt # Install dependencies for development (e.g. black, codespell, etc.)
    pre-commit install # Install git hooks for automatic code formatting checks
    deactivate
    ```
3. **Create a new branch**
    - Branch from main for each change.
    - Use descriptive branch names, such as `feature/add-new-plot` or `bug/fix-units`.
    ```shell
    git switch main
    git switch -c <branch-name>
    ```
4. **Make your changes**
    - Follow the [code](docs/CODE_STYLE.md) and [writing](docs/WRITING_STYLE.md) guidelines.

   #### Other Guidelines

   Depending on what you're working on, the following documents may also be helpful:

    - [Type Hint and Docstring Style](docs/TYPE_HINT_AND_DOCSTRING_STYLE.md): if you're adding or modifying type hints or docstrings.
    - [Classes and Immutability](docs/CLASSES_AND_IMMUTABILITY.md): if you're adding or modifying classes or their attributes.
    - [Angle Vectors and Transformations](docs/ANGLE_VECTORS_AND_TRANSFORMATIONS.md): if you're working with rotation matrices, angle vectors, or coordinate transformations.
    - [Axes, Points, and Frames](docs/AXES_POINTS_AND_FRAMES.md): if you're working with vector-valued quantities such as positions, velocities, forces, or moments.

   #### Core Modeling Contributions

   Some parts of Ptera Software, such as the UVLM formulation, wake evolution, force and moment integration, and related aerodynamic models, are closely tied to underlying theoretical assumptions and numerical constraints.

   Contributions in these areas are very welcome, but they typically require coordination and shared understanding of the relevant aerodynamics. When opening a PR that touches core modeling logic, please:
    - Clearly describe the physical assumptions being made.
    - Note the regime in which the change is expected to be valid.
    - Include references, derivations, or reasoning where appropriate.
    - Add tests that validate the behavior against known physical expectations (e.g., symmetry, limiting cases, conservation behavior).

   If you’re unsure whether a change falls into this category, feel free to ask in the issue or PR thread before investing significant effort.
5. **Commit your work**
    - Run automated checks locally:
    ```shell
    source .venv/bin/activate # On Windows use .venv\Scripts\activate
    pre-commit run --all-files
    mypy pterasoftware
    python -m unittest discover -s tests
    ```
    - If the checks pass, check which files you've modified:
    ```shell
    git status
    ```
    - Stage the changes you'd like to commit:
    ```shell
    git add <file-1> <file=2> ...
    ```
    - Before committing, note that you need to know how to use Git's text editor, which defaults to Vim. If you've never used Vim before, I recommend telling Git to use Nano instead, as many people find it less confusing. To do so, run the following (you only need to run it once):
    ```shell
    git config --global core.editor "nano"
    ```
    - Now, commit your changes:
    ```shell
    git commit
    ```

   This will open Nano, where you enter your commit message above the commented out `#` lines. Please follow best practices for commit messages. Below is a lightly adapted version of [Tim Pope's famous sample commit message](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html):
    ```shell
    Capitalized, short (50 chars or less) summary
    
    More detailed explanatory text, if necessary. Wrap it to about 72
    characters or so. In some contexts, the first line is treated as the
    subject of an email and the rest of the text as the body. The blank
    line separating the summary from the body is critical (unless you omit
    the body entirely); tools like rebase can get confused if you run the
    two together.
    
    Write your commit message in the imperative: "Fix bug" and not "Fixed
    bug" or "Fixes bug." This convention matches up with commit messages
    generated by commands like git merge and git revert.
    
    Further paragraphs come after blank lines.
    
    - Bullet points are okay, too
    
    - Typically a hyphen or asterisk is used for the bullet, followed by a
      single space, with blank lines in between, but conventions vary here
    
    - Use a hanging indent
    ```
   Once you are satisfied with your commit message, enter Ctrl+O and Ctrl+X to save and exit, which will trigger the commit. Alternatively, with an empty commit message, enter Ctrl+X to exit and abort the commit.
6. **Push your changes and open a PR**
    - Push your branch to your fork:
    ```shell
    git push -u origin <branch-name> # After the first push, just use git push
    ```
    - Open a PR from your branch to the main branch of the upstream repository.
        - You can open the PR as a **Draft** to get feedback early before the work is complete. Draft PRs indicate that design details may still be changing.
        - In the PR description, follow the PR template, and link any related issues.

7. **Keeping your branch up to date**
    - If main changes before your PR is merged, sync your branch to avoid merge conflicts:
    ```shell
    git fetch upstream
    git merge upstream/main
    ```
8. **Review and approval process**
    - Only the repository owner (currently [camUrban](https://github.com/camUrban)) can approve merges to main.
    - Your PR will be reviewed, and changes may be requested.
    - Once approved, it will be merged into main and included in the next release.