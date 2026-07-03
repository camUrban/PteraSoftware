---
description: Find and delete unused fixtures and dead setUp attributes across the whole test suite with scripts/find_unused_fixtures.py
---

# Delete Unused Fixtures

Find and remove unused test fixtures and dead `setUp`/`setUpClass` attributes across the entire test suite, driving `scripts/find_unused_fixtures.py`. This command takes no argument: the script scans the whole suite at once (both `tests/unit/` and `tests/integration/`), which is necessary because a fixture defined in one module is routinely used by tests in another.

## What the script detects

The script parses the suite with the `ast` module instead of grepping for names, so it resolves each call to the specific fixture definition it refers to and does not raise false matches on coincidentally shared names. It reports three categories:

- **Directly unused fixtures**: `make_*` fixture functions never called anywhere under `tests/`.
- **Transitively unused fixtures**: fixtures called only by other unused fixtures. Removing the directly unused ones can expose these, so the analysis iterates to a fixed point.
- **Dead `setUp`/`setUpClass` attributes**: attributes assigned in a test class's `setUp` or `setUpClass` that no test method, and no other `setUp` attribute, ever reads.

## Why deletion is safe for `__init__.py`

Both `tests/unit/fixtures/__init__.py` and `tests/integration/fixtures/__init__.py` are bare package markers: a single docstring with no imports. Deleting a fixture function, or even an entire fixture module, therefore cannot break an `__init__.py`, and you do not edit any `__init__.py` in this command.

## Environment constraints

- The script's verification modes run the full test suite under `coverage` (pinned in `requirements_dev.txt`). This takes several minutes. Always invoke the script with `python -u` so its output streams, and never pipe the run through `tail`, `head`, `grep`, or any other filter, per `docs/RUNNING_TESTS_AND_TYPE_CHECKS.md`.
- The script never modifies a file until you pass `--delete-verified`. The plain run and `--verify` leave every file untouched.

## Steps

1. **Report.** Run the script with no flags to list candidates without changing anything:
   ```bash
   python -u scripts/find_unused_fixtures.py
   ```
   It prints the three categories above with each candidate's file and line, then a summary count. It exits non-zero whenever any candidate is found; that is expected and is not an error.
2. **Present the candidates** to the user: the directly unused fixtures, the transitively unused fixtures, and the dead `setUp` attributes, each with its location. If the report is empty, report that the suite is clean and stop.
3. **Optionally preview the verification.** If the user wants to see which candidates survive verification before authorizing any deletion, run the dry-run verification, which deletes each candidate in memory, runs the full suite under coverage, checks that all tests still pass, that the set of discoverable test IDs is unchanged, and that coverage does not drop, then restores every file:
   ```bash
   python -u scripts/find_unused_fixtures.py --verify
   ```
   Report which candidates it confirmed and which it flagged as false positives (a candidate whose name surfaces in a failing traceback after deletion).
4. **Delete, on confirmation.** Only after the user confirms, run the deletion pass. It re-runs the same verification, deletes only the candidates that pass, excludes any false positives, and restores files if verification fails outright:
   ```bash
   python -u scripts/find_unused_fixtures.py --delete-verified
   ```
5. **Reformat the touched files** through pre-commit, never a bare formatter, per `docs/RUNNING_TESTS_AND_TYPE_CHECKS.md`. Find the modified files with `git status -sb`, then:
   ```bash
   pre-commit run --files <each modified fixture or test file>
   ```
6. **Re-scan for newly exposed candidates.** Removing a fixture can leave another one unused. The deletion pass reports whether new candidates were exposed; if so, repeat from step 1 until a report comes back clean.

## Important Reminders

- This command operates on the whole suite and takes no argument; do not try to scope it to a single module.
- Never edit an `__init__.py` after deleting fixtures: the fixture packages are bare markers with no imports, so a deleted fixture or module cannot break one.
- Always show the user the report, and the verification result if you ran one, before any deletion. The script only mutates files under `--delete-verified`.
- Run the script with `python -u` and never pipe it through a filter; the verification suite can run for minutes and its output must stay visible.
- Reformat modified files with `pre-commit run --files`, never a bare formatter.
