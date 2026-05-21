# Running Tests and Type Checks

> This document describes how AI coding agents (Claude Code and similar) should invoke linters, formatters, type checkers, tests, and ad-hoc scripts when working on Ptera Software. **Human contributors should follow [`CONTRIBUTING.md`](../CONTRIBUTING.md) instead**, which documents the developer-facing workflow.

The guidance below assumes the project's virtual environment has already been created and activated, and that the package has been installed in editable mode via `pip install -e .` per the setup steps in `CONTRIBUTING.md`. These steps are part of normal project setup and are not repeated here. When the venv is active, `python`, `mypy`, `pre-commit`, and the other tools resolve to their venv copies on `PATH`.

## Linters, Formatters, and Spell-Checkers

Every tool configured as a hook in `.pre-commit-config.yaml` (currently `end-of-file-fixer`, `trailing-whitespace`, `mixed-line-ending`, `fix-byte-order-marker`, `requirements-txt-fixer`, `isort`, `black`, `docformatter`, `check-merge-conflict`, `check-case-conflict`, `check-illegal-windows-names`, `check-symlinks`, `destroyed-symlinks`, `check-yaml`, `check-toml`, `check-xml`, `check-json`, `detect-private-key`, `forbid-submodules`, `check-vcs-permalinks`, `check-executables-have-shebangs`, `check-shebang-scripts-are-executable`, `check-added-large-files`, `debug-statements`, `codespell`, `ascii-only`, `mypy`) must always be invoked through pre-commit, never directly:

```shell
pre-commit run --all-files                                 # every hook, every tracked file
pre-commit run --all-files codespell                       # one hook against the whole tree
pre-commit run --files pterasoftware/foo.py docformatter   # one hook against specific files
```

Note that `pre-commit run <hook>` without `--all-files` or `--files` only runs the hook against files currently staged in git, which is rarely what you want when verifying work in progress.

Invoking a bare entry point such as `codespell ...` or `docformatter ...` is forbidden. The bare version can disagree with the pre-commit version because pre-commit pins each hook to a specific release in its own isolated environment, while the venv copy may have drifted ahead or behind. The pre-commit run is the source of truth that CI uses.

## Type Checking

mypy is configured as a pre-commit hook and reads `mypy.ini` from the project root automatically. Run it with:

```shell
pre-commit run --all-files mypy
```

## Tests

Ptera Software uses the standard library's `unittest`, not pytest. Run the full suite from the activated venv:

```shell
python -u -m unittest discover -s tests
```

For a single module:

```shell
python -u -m unittest tests.unit.test_module_name -v
```

For a single test method:

```shell
python -u -m unittest tests.unit.test_module_name.TestClass.test_method -v
```

## Example and Experimental Scripts

The package is editable-installed, so any script that imports `pterasoftware` works from any working directory with the bare interpreter:

```shell
python -u examples/steady_horseshoe_vortex_lattice_method_solver.py
python -u experimental/<whatever>.py
```

Do not prefix with `.venv/bin/python` or `.venv/Scripts/python.exe`, do not set `PYTHONPATH`, and do not chain `cd <project-root> && ...`. The activated venv puts the right interpreter on `PATH` and the editable install handles import resolution.

## Output Discipline

The harness captures Bash stdout through a pipe, so Python defaults to block-buffered output. Without intervention, a ten minute simulation can sit silent for the entire run and then dump everything at exit. To keep the user able to see what is happening:

- **Always pass `python -u`** (or `python -u -m`) when invoking Python directly so stdout and stderr are line-buffered. For a tool that buffers regardless of how it is invoked, set `PYTHONUNBUFFERED=1` in the shell.
- **Never pipe** a running command's output through `tail`, `head`, `grep`, `sed`, `awk`, or any filter that batches the stream. These swallow output and emit only the trimmed remainder at exit, defeating real-time visibility.
- **Never redirect to `/dev/null`.** If output should be discarded, say so in the conversation first so the user knows nothing is being captured.
- **Short runs (under about five minutes): foreground, no redirect.** Let output stream directly into the conversation.
- **Long runs: redirect to `experimental/<descriptive>.log`** with the path written plainly in the Bash command, combined with `run_in_background: true` so the harness tracks completion. The `experimental/` directory is gitignored (`/experimental/` appears in `.gitignore`), so the log cannot accidentally be committed. Read excerpts from the log at meaningful checkpoints and surface them in the conversation. Never let a backgrounded job finish with no output visible to the user.

Example pattern for a long run:

```shell
python -u experimental/long_running_thing.py > experimental/long_running_thing.log 2>&1
```

Launched with `run_in_background: true`, the file grows line by line as the script runs and can be Read at any time.

The trailing `2>&1` is required: it duplicates stderr to wherever stdout is going, so Python's logging output, NumPy/PyVista warnings, and tracebacks (all of which default to stderr) land in the same visible log instead of splitting off into the harness's separate stderr capture. Order matters: `> <file> 2>&1` works; `2>&1 > <file>` leaves stderr on the terminal.
