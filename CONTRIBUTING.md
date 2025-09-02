# Contributing Guidelines

We are excited that you are interested in contributing to **Ptera Software**! This 
guide will help you get started. If you have any questions not answered here, please 
[open an issue](https://github.com/camUrban/PteraSoftware/issues) and we will respond.  

---

## Before Contributing

Please review the following documents before making contributions:  

1. [README](README.md)  
2. [Code of Conduct](CODE_OF_CONDUCT.md)  
3. [Security Policy](SECURITY.md)  
4. [License](LICENSE.txt)  

---

## Ways to Contribute

There are three main ways you can contribute:  

1. [**Report a bug**](#reporting-a-bug)  
   - Identify and document issues that prevent Ptera Software from working as 
   expected.  
   - This includes errors, crashes, incorrect outputs, or unexpected behavior.  

2. [**Request a feature**](#requesting-a-feature)  
   - Suggest new features, improvements to existing functionality, or usability 
   changes.  
   - This can be small quality-of-life improvements or larger feature proposals.  

3. [**Contribute code**](#contributing-code)  
   - Submit changes that add new features, fix bugs, improve performance, or enhance 
   documentation.  
   - Code contributions can address your own ideas or work on existing open issues.  
   - If you find an issue labeled `help wanted` or `good first issue` that you want to 
   work on, **comment on the issue to claim it** before starting work. This prevents 
   duplicate efforts.  

---

### Reporting a Bug

**If the bug is a security vulnerability, do not post it as a public issue**. Follow 
the [security policy](SECURITY.md) instead.  

For all other bugs:  
1. Search the [issues page](https://github.com/camUrban/PteraSoftware/issues) to check 
   if it's already been reported.  
   - If it exists, you can comment or add an emoji reaction to indicate that you are 
   also affected.  
2. If it has not been reported, open a new issue using the 
   [bug report template](.github/ISSUE_TEMPLATE/bug_report.md).  
3. Add the `bug` and `help wanted` labels and any other relevant labels.  

---

### Requesting a Feature

For feature requests or enhancements:  
1. Search the [issues page](https://github.com/camUrban/PteraSoftware/issues) to avoid 
   duplicates.  
   - If it exists, comment or react to indicate interest.  
2. If it has not been requested, open a new issue using the 
   [feature request template](.github/ISSUE_TEMPLATE/feature_request.md).  
3. Add the `feature` and `help wanted` labels plus any others that apply.  

---

### Contributing Code

Ptera Software now uses the GitHub Flow to manage code contributions. If this is new to 
you, it's a good idea to read through 
[this guide](https://docs.github.com/en/get-started/using-github/github-flow) first. 
Once you understand the process, here's how to implement it:

1. **Choose what to work on**
   - Look for issues labeled `help wanted` or `good first issue`.  
     - If you want to work on one, check that no one else has already commented 
     claiming it. If unclaimed, comment on the issue to claim it.  
     - If you want to work on a claimed issue that hasn't been updated in a while, 
     write a comment asking if the user who originally claimed it is still actively 
     working on it.
   - If you have your own idea, search the issues to ensure it hasn't already been 
   proposed.  
     - If you can't find anything, open a new issue describing your change and state 
     that you will be implementing it.  

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
   - Create a virtual environment and install dependencies for development:  
     ```shell
     python -m venv .venv
     .venv\Scripts\activate # On Mac or Linux use source .venv/bin/activate
     python -m pip install --upgrade pip setuptools wheel
     pip install -r requirements.txt # Install dependencies for running simulations
     pip install -r requirements_dev.txt # Install dependencies for development (e.g. black, codespell, etc.)
     deactivate
     ```  

3. **Create a new branch**
   - Branch from main for each change.  
   - Use descriptive branch names, such as feature/add-new-plot or 
   bugfix/fix-units.  
   ```shell
   git checkout main
   git branch <branch-name>
   git checkout <branch-name>
   git commit -m "Created a new branch for <description>"
   git push origin <branch-name>
   ```  

4. **Make your changes**  
   - Commit frequently with clear, descriptive messages.  
   - Follow the [code style and standards](#code-style-and-standards).  
   - Run automated checks locally before pushing:  
     ```shell
     .venv\Scripts\activate # On Mac or Linux use source .venv/bin/activate
     codespell --ignore-words=.codespell-ignore.txt --skip="*.dat"
     black .
     python -m unittest discover -s tests
     ```  

5. **Push your changes and open a pull request**  
   - Push your branch to your fork:  
     ```shell
     git push origin <branch-name>
     ```  
   - Open a pull request (PR) from your branch to the main branch of the upstream repository.  
   - You can open the PR as a **Draft** to get feedback early before the work is complete.  
   - In the PR description, follow the [pull request template](.github/pull_request_template.md) and link any related issues.

6. **Keeping your branch up to date**  
   - If main changes before your PR is merged, sync your branch to avoid merge conflicts:  
     ```shell
     git fetch upstream
     git merge upstream/main # Or git rebase upstream/main
     ```  

7. **Review and approval process**  
   - Only the repository owner (currently @camUrban) can approve merges to main.  
   - Your PR will be reviewed, and changes may be requested.  
   - Once approved, it will be merged into main and included in the next release.

---

## Code Style and Standards

- Always use S.I. units in calculations and results.  
- Format code using [black](https://github.com/psf/black) with a line length of 88.  
- Include docstrings for all new modules, classes, functions, and methods in 
[reStructuredText format](https://realpython.com/documenting-python-code/).  
- Use block comments where needed for clarity.  
- Tag comments with `TODO`, `BUG`, or similar where applicable.  
- Cite any external sources such as code, equations, or algorithms in comments or 
  docstrings.

---