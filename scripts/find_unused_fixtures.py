"""Finds unused fixture functions and dead setUp attributes in the test suite."""

import argparse
import ast
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = PROJECT_ROOT / "tests"

FIXTURE_DIRS = [
    TESTS_DIR / "unit" / "fixtures",
    TESTS_DIR / "integration" / "fixtures",
]

TEST_DIRS = [
    TESTS_DIR / "unit",
    TESTS_DIR / "integration",
]


def find_fixture_functions() -> dict[tuple[str, str], tuple[int, int]]:
    """Finds all fixture function definitions across all fixture modules.

    :return: A dictionary mapping (module_path, function_name) to a (lineno, end_lineno)
        tuple.
    """
    fixtures = {}
    for fixture_dir in FIXTURE_DIRS:
        for filepath in sorted(fixture_dir.glob("*.py")):
            if filepath.name == "__init__.py":
                continue
            tree = ast.parse(filepath.read_text(), filename=str(filepath))
            for node in ast.walk(tree):
                if isinstance(node, ast.FunctionDef) and node.name.startswith("make_"):
                    end_lineno = node.end_lineno or node.lineno
                    fixtures[(str(filepath), node.name)] = (node.lineno, end_lineno)
    return fixtures


def find_all_python_files() -> list[Path]:
    """Finds all Python files in the test directories.

    :return: A list of Path objects for all .py files.
    """
    files = []
    for test_dir in TEST_DIRS:
        files.extend(sorted(test_dir.rglob("*.py")))
    return files


def _resolve_import_to_filepath(
    module_str: str, level: int, caller_filepath: str
) -> str | None:
    """Resolves a Python import to a filesystem path.

    Handles both absolute and relative imports by converting the dotted module string to
    a path and checking if the corresponding .py file exists.

    :param module_str: The dotted module string (e.g., "tests.unit.fixtures.foo").
    :param level: The number of leading dots for relative imports (0 for absolute).
    :param caller_filepath: The path of the file containing the import.
    :return: The resolved file path as a string, or None if the file does not exist.
    """
    if level > 0:
        base = Path(caller_filepath).parent
        for _ in range(level - 1):
            base = base.parent
        if module_str:
            resolved = base / Path(*module_str.split("."))
        else:
            resolved = base
    else:
        resolved = PROJECT_ROOT / Path(*module_str.split("."))

    py_file = resolved.with_suffix(".py")
    if py_file.exists():
        return str(py_file)

    return None


def _parse_fixture_imports(
    tree: ast.Module, caller_filepath: str, fixture_file_set: set[str]
) -> dict[str, str]:
    """Parses a file's imports to map local names to their source fixture files.

    Handles qualified calls through imported fixture modules and unqualified calls
    within fixture files that reference functions defined in the same file.

    :param tree: The parsed AST of the file.
    :param caller_filepath: The path of the file being parsed.
    :param fixture_file_set: The set of all known fixture file paths (as strings).
    :return: A dictionary mapping local module alias names (e.g., "foo_fixtures") to
        their resolved fixture file paths.
    """
    module_aliases = {}

    for node in ast.iter_child_nodes(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.module is None and node.level == 0:
            continue

        for alias in node.names:
            combined = f"{node.module}.{alias.name}" if node.module else alias.name
            source = _resolve_import_to_filepath(combined, node.level, caller_filepath)
            if source and source in fixture_file_set:
                local_alias = alias.asname or alias.name
                module_aliases[local_alias] = source

    return module_aliases


def find_fixture_call_sites(
    all_files: list[Path], fixture_functions: dict[tuple[str, str], tuple[int, int]]
) -> dict[tuple[str, str], set[tuple[str, str, str]]]:
    """Finds where each fixture function is called.

    Uses import resolution to attribute calls to specific fixture definitions. When a
    call cannot be resolved to a specific fixture, it is conservatively attributed to
    all fixtures with that name to avoid false positives.

    :param all_files: The list of all Python files to search.
    :param fixture_functions: The dictionary of all fixture function definitions, keyed
        by (filepath, function_name).
    :return: A dictionary mapping (filepath, function_name) to a set of
        (caller_filepath, caller_context, caller_name) tuples.
    """
    fixture_file_set = {fp for fp, _ in fixture_functions}

    name_to_keys: dict[str, list[tuple[str, str]]] = {}
    for fp, name in fixture_functions:
        name_to_keys.setdefault(name, []).append((fp, name))

    call_sites: dict[tuple[str, str], set[tuple[str, str, str]]] = {
        key: set() for key in fixture_functions
    }

    for filepath in all_files:
        try:
            source = filepath.read_text()
            tree = ast.parse(source, filename=str(filepath))
        except SyntaxError:
            continue

        filepath_str = str(filepath)
        module_aliases = _parse_fixture_imports(tree, filepath_str, fixture_file_set)

        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue

            caller_context = "fixture" if node.name.startswith("make_") else "test"

            for child in ast.walk(node):
                if not isinstance(child, ast.Call):
                    continue

                func = child.func
                resolved_source = None
                func_name = None

                if isinstance(func, ast.Attribute) and func.attr.startswith("make_"):
                    # Qualified call: module.make_foo().
                    func_name = func.attr
                    if isinstance(func.value, ast.Name):
                        resolved_source = module_aliases.get(func.value.id)
                        # Same file call within a fixture module.
                        if resolved_source is None and filepath_str in fixture_file_set:
                            resolved_source = filepath_str

                elif isinstance(func, ast.Name) and func.id.startswith("make_"):
                    # Unqualified call: make_foo() (same file).
                    func_name = func.id
                    if filepath_str in fixture_file_set:
                        resolved_source = filepath_str

                if func_name is None:
                    continue

                caller_info = (filepath_str, caller_context, node.name)

                if resolved_source:
                    key = (resolved_source, func_name)
                    if key in call_sites:
                        call_sites[key].add(caller_info)
                else:
                    # Fallback: attribute to all definitions with this name.
                    for key in name_to_keys.get(func_name, []):
                        call_sites[key].add(caller_info)

    return call_sites


def find_directly_unused(
    fixture_functions: dict[tuple[str, str], tuple[int, int]],
    call_sites: dict[tuple[str, str], set[tuple[str, str, str]]],
) -> set[tuple[str, str]]:
    """Finds fixture functions that are never called anywhere.

    :param fixture_functions: The dictionary of all fixture function definitions.
    :param call_sites: The dictionary of call sites keyed by (filepath, name).
    :return: A set of (filepath, function_name) keys that are never called.
    """
    return {key for key in fixture_functions if not call_sites.get(key)}


def find_transitively_unused(
    fixture_functions: dict[tuple[str, str], tuple[int, int]],
    call_sites: dict[tuple[str, str], set[tuple[str, str, str]]],
) -> set[tuple[str, str]]:
    """Finds fixture functions that are only called by other unused fixtures.

    This iteratively marks fixtures as unused if all their callers are themselves
    unused, until no more changes occur. Caller identity is resolved using the caller's
    filepath and function name.

    :param fixture_functions: The dictionary of all fixture function definitions.
    :param call_sites: The dictionary of call sites keyed by (filepath, name).
    :return: A set of (filepath, function_name) keys that are transitively unused.
    """
    all_keys = set(fixture_functions.keys())
    unused: set[tuple[str, str]] = set()

    changed = True
    while changed:
        changed = False
        for key in all_keys - unused:
            sites = call_sites.get(key, set())
            if not sites:
                unused.add(key)
                changed = True
                continue

            all_callers_unused = all(
                context == "fixture"
                and (caller_filepath, caller_name) in all_keys
                and (caller_filepath, caller_name) in unused
                for caller_filepath, context, caller_name in sites
            )
            if all_callers_unused:
                unused.add(key)
                changed = True

    return unused


def find_dead_setup_attributes(
    all_files: list[Path],
) -> list[tuple[str, str, str, int, int]]:
    """Finds setUp/setUpClass attributes that are never used in the class.

    Parses each test file for classes that have a setUp or setUpClass method. For each
    attribute assigned in setUp/setUpClass (e.g., self.foo = ...), checks whether that
    attribute name is ever referenced in any other method of the same class or read (in
    a Load context) within setUp/setUpClass itself. This catches attributes that are
    only used to construct other setUp attributes.

    :param all_files: The list of all Python files to search.
    :return: A list of (filepath, class_name, attr_name, lineno, end_lineno) tuples.
    """
    dead_attrs = []

    for filepath in all_files:
        if not filepath.name.startswith("test_"):
            continue

        try:
            source = filepath.read_text()
            tree = ast.parse(source, filename=str(filepath))
        except SyntaxError:
            continue

        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue

            setup_attrs = {}
            setup_methods = []
            other_methods = []

            for item in node.body:
                if not isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    continue

                if item.name in ("setUp", "setUpClass"):
                    setup_methods.append(item)
                    for stmt in ast.walk(item):
                        if not isinstance(stmt, ast.Assign):
                            continue
                        for target in stmt.targets:
                            if (
                                isinstance(target, ast.Attribute)
                                and isinstance(target.value, ast.Name)
                                and target.value.id in ("self", "cls")
                            ):
                                end_lineno = stmt.end_lineno or stmt.lineno
                                setup_attrs[target.attr] = (stmt.lineno, end_lineno)

                else:
                    other_methods.append(item)

            has_test_methods = any(m.name.startswith("test_") for m in other_methods)
            if not setup_attrs or not has_test_methods:
                continue

            # Collect attribute references from non-setUp methods.
            used_attrs = set()
            for method in other_methods:
                for child in ast.walk(method):
                    if (
                        isinstance(child, ast.Attribute)
                        and isinstance(child.value, ast.Name)
                        and child.value.id in ("self", "cls")
                    ):
                        used_attrs.add(child.attr)

            # Also collect attribute READS within setUp/setUpClass itself. This catches
            # attributes only used to construct other setUp attributes (e.g.,
            # self.derived = make_derived(self.base)). Only Load context counts; Store
            # context is the assignment itself.
            for method in setup_methods:
                for child in ast.walk(method):
                    if (
                        isinstance(child, ast.Attribute)
                        and isinstance(child.value, ast.Name)
                        and child.value.id in ("self", "cls")
                        and isinstance(child.ctx, ast.Load)
                    ):
                        used_attrs.add(child.attr)

            for attr_name, (lineno, end_lineno) in sorted(setup_attrs.items()):
                if attr_name not in used_attrs:
                    dead_attrs.append(
                        (str(filepath), node.name, attr_name, lineno, end_lineno)
                    )

    return dead_attrs


def _batch_delete(deletions: list[tuple[str, int, int]]) -> dict[str, str]:
    """Deletes multiple line ranges across multiple files.

    Within each file, deletions are applied from bottom to top so that line numbers
    remain valid.

    :param deletions: A list of (filepath, start_line, end_line) tuples.
    :return: A dictionary mapping filepath to original content for restoration.
    """
    by_file: dict[str, list[tuple[int, int]]] = {}
    for filepath, start_line, end_line in deletions:
        by_file.setdefault(filepath, []).append((start_line, end_line))

    originals = {}
    for filepath, ranges in by_file.items():
        content = Path(filepath).read_text()
        originals[filepath] = content
        lines = content.splitlines(keepends=True)
        for start_line, end_line in sorted(ranges, reverse=True):
            lines = lines[: start_line - 1] + lines[end_line:]
        Path(filepath).write_text("".join(lines))

    return originals


def _restore_files(originals: dict[str, str]) -> None:
    """Restores files to their original content.

    :param originals: A dictionary mapping filepath to original content.
    """
    for filepath, content in originals.items():
        Path(filepath).write_text(content)


def _discover_test_ids() -> list[str]:
    """Discovers all test IDs without running them.

    Uses unittest's test loader to walk the test suite and collect fully qualified test
    IDs (e.g., ``tests.unit.test_foo.TestBar.test_baz``). This is fast, deterministic,
    and not affected by output interleaving from warnings or multi-line docstrings.

    :return: A sorted list of test ID strings.
    """
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import unittest, sys\n"
                "loader = unittest.TestLoader()\n"
                f"suite = loader.discover({str(TESTS_DIR)!r})\n"
                "def collect(s):\n"
                "    ids = []\n"
                "    for t in s:\n"
                "        if isinstance(t, unittest.TestSuite):\n"
                "            ids.extend(collect(t))\n"
                "        else:\n"
                "            ids.append(t.id())\n"
                "    return ids\n"
                "for tid in sorted(collect(suite)):\n"
                "    print(tid)\n"
            ),
        ],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
        timeout=60,
    )
    if result.returncode != 0:
        print("  [ERROR] Test discovery failed:")
        print(result.stderr)
        return []
    return sorted(result.stdout.strip().splitlines())


def _run_full_test_suite() -> tuple[bool, str, str, float | None]:
    """Runs the full test suite under coverage in a subprocess.

    :return: A (success, summary, stderr, coverage_pct) tuple where success is a bool,
        summary is a short string, stderr is the full stderr output, and coverage_pct is
        the overall coverage percentage (or None if coverage reporting failed).
    """
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "coverage",
            "run",
            "--source=pterasoftware",
            "-m",
            "unittest",
            "discover",
            "-s",
            str(TESTS_DIR),
        ],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
        timeout=600,
    )
    stderr = result.stderr

    # Extract coverage percentage.
    coverage_pct = None
    cov_result = subprocess.run(
        [sys.executable, "-m", "coverage", "report"],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
        timeout=30,
    )
    if cov_result.returncode == 0:
        for line in cov_result.stdout.splitlines():
            if line.startswith("TOTAL"):
                parts = line.split()
                for part in parts:
                    if part.endswith("%"):
                        try:
                            coverage_pct = float(part.rstrip("%"))
                        except ValueError:
                            pass
                        break
                break

    if result.returncode == 0:
        for line in stderr.splitlines():
            if line.startswith("Ran "):
                return True, line.strip().lower(), stderr, coverage_pct
        return True, "all tests passed", stderr, coverage_pct
    stderr_lines = stderr.strip().splitlines()
    last_line = stderr_lines[-1] if stderr_lines else "unknown error"
    return False, last_line, stderr, coverage_pct


def _find_false_positives_from_traceback(
    stderr: str, all_deletions: list[tuple[str, int, int, str]]
) -> set[str]:
    """Analyzes test suite stderr to identify which deletions caused failures.

    Scans the traceback for references to deleted fixture function names and deleted
    setUp attribute names to determine which candidates are false positives.

    :param stderr: The full stderr output from the test suite run.
    :param all_deletions: The list of (filepath, start_line, end_line, label) tuples for
        all candidates.
    :return: A set of labels that appear to be false positives based on the traceback.
    """
    false_positives = set()
    for _, _, _, label in all_deletions:
        # The label is like "path:line  name" or "path:line  Class.attr". Extract the
        # name or attr from the end.
        candidate_name = label.split()[-1]

        # For "Class.attr" form, extract just the attr.
        if "." in candidate_name:
            candidate_name = candidate_name.split(".")[-1]

        if candidate_name in stderr:
            false_positives.add(label)

    return false_positives


def _format_key(key: tuple[str, str]) -> tuple[Path, str]:
    """Formats a (filepath, name) fixture key for display.

    :param key: A (filepath, function_name) tuple.
    :return: A (relative_path, function_name) tuple for display.
    """
    filepath, name = key
    rel = Path(filepath).relative_to(PROJECT_ROOT)
    return rel, name


def main() -> int:
    """Runs the unused fixture analysis and optionally verifies each candidate.

    :return: An int representing the exit code, which is 0 if no issues were found and 1
        otherwise.
    """
    parser = argparse.ArgumentParser(
        description="Find unused fixture functions in the test suite."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--verify",
        action="store_true",
        help=(
            "Verify candidates by deleting them and running the full test "
            "suite with coverage. Reports results but does not modify files."
        ),
    )
    group.add_argument(
        "--delete-verified",
        action="store_true",
        help=(
            "Verify candidates (same as --verify), then permanently delete "
            "all confirmed candidates if verification passes."
        ),
    )
    args = parser.parse_args()

    print("Scanning fixture definitions...", flush=True)
    fixture_functions = find_fixture_functions()
    print(f"  Found {len(fixture_functions)} fixture functions.\n")

    print("Scanning for call sites...", flush=True)
    all_files = find_all_python_files()
    call_sites = find_fixture_call_sites(all_files, fixture_functions)

    directly_unused = find_directly_unused(fixture_functions, call_sites)
    transitively_unused = find_transitively_unused(fixture_functions, call_sites)
    only_transitively_unused = transitively_unused - directly_unused

    all_unused_fixtures = sorted(transitively_unused)

    print("=" * 80)
    print("DIRECTLY UNUSED FIXTURES (never called anywhere)")
    print("=" * 80)
    if directly_unused:
        for key in sorted(directly_unused):
            rel, name = _format_key(key)
            filepath, _ = key
            lineno, _ = fixture_functions[key]
            print(f"  {rel}:{lineno}  {name}")
    else:
        print("  None found.")

    print()
    print("=" * 80)
    print("TRANSITIVELY UNUSED FIXTURES (only called by other unused fixtures)")
    print("=" * 80)
    if only_transitively_unused:
        for key in sorted(only_transitively_unused):
            rel, name = _format_key(key)
            lineno, _ = fixture_functions[key]
            callers = call_sites[key]
            caller_names = sorted({c for _, _, c in callers})
            print(f"  {rel}:{lineno}  {name}")
            print(f"    Called only by: {', '.join(caller_names)}")
    else:
        print("  None found.")

    print()
    print("=" * 80)
    print("DEAD setUp/setUpClass ATTRIBUTES (assigned but never used in tests)")
    print("=" * 80)
    dead_attrs = find_dead_setup_attributes(all_files)
    if dead_attrs:
        for filepath, class_name, attr_name, lineno, end_lineno in dead_attrs:
            rel = Path(filepath).relative_to(PROJECT_ROOT)
            print(f"  {rel}:{lineno}  {class_name}.{attr_name}")
    else:
        print("  None found.")

    total_issues = len(transitively_unused) + len(dead_attrs)
    print()
    print(
        f"Total: {len(directly_unused)} directly unused, "
        f"{len(only_transitively_unused)} transitively unused, "
        f"{len(dead_attrs)} dead setUp attributes."
    )

    if not args.verify and not args.delete_verified:
        if total_issues > 0:
            print(
                "\nRe-run with --verify to confirm, or --delete-verified to "
                "confirm and delete."
            )
        return 1 if total_issues > 0 else 0

    # Build the list of all deletions: (filepath, start_line, end_line, label).
    all_deletions = []
    for key in all_unused_fixtures:
        filepath, name = key
        lineno, end_lineno = fixture_functions[key]
        rel = Path(filepath).relative_to(PROJECT_ROOT)
        all_deletions.append((filepath, lineno, end_lineno, f"{rel}:{lineno}  {name}"))

    for filepath, class_name, attr_name, lineno, end_lineno in dead_attrs:
        rel = Path(filepath).relative_to(PROJECT_ROOT)
        all_deletions.append(
            (filepath, lineno, end_lineno, f"{rel}:{lineno}  {class_name}.{attr_name}")
        )

    if not all_deletions:
        return 0

    print()
    print("=" * 80)
    print("VERIFICATION (full test suite + coverage)")
    print("=" * 80)

    # Baseline: discover test IDs (fast, no execution) and run the full suite under
    # coverage. Post-deletion runs must match on three dimensions:
    # 1. All tests pass.
    # 2. Same test IDs are discoverable (no test methods lost or broken).
    # 3. Coverage does not drop (no subtests or code paths silently lost).
    print()
    print(
        "Baseline: Discovering test IDs and running test suite with coverage...",
        flush=True,
    )
    baseline_ids = _discover_test_ids()
    if not baseline_ids:
        print("  [FAIL] Baseline discovery found no tests.")
        return 1

    baseline_ok, baseline_summary, _, baseline_cov = _run_full_test_suite()
    if not baseline_ok:
        print(f"  [FAIL] Baseline test suite failed ({baseline_summary})")
        print("Cannot verify candidates against a failing baseline.")
        return 1
    cov_str = f", {baseline_cov:.2f}% coverage" if baseline_cov is not None else ""
    print(f"  Baseline: {len(baseline_ids)} test IDs, " f"{baseline_summary}{cov_str}.")

    # Phase 1: batch-delete ALL candidates, then verify all three dimensions.
    print()
    print(
        f"Phase 1: Batch-deleting all {len(all_deletions)} candidates and "
        f"running the full test suite with coverage...",
        flush=True,
    )

    batch = [(fp, sl, el) for fp, sl, el, _ in all_deletions]
    originals = _batch_delete(batch)
    try:
        ok, summary, stderr, post_cov = _run_full_test_suite()
        post_ids = _discover_test_ids() if ok else []
    finally:
        _restore_files(originals)

    if ok:
        # Check test ID match.
        missing = sorted(set(baseline_ids) - set(post_ids))
        extra = sorted(set(post_ids) - set(baseline_ids))
        if missing or extra:
            print(f"  [FAIL] Tests passed but test ID mismatch ({summary})")
            if missing:
                print(f"  Missing {len(missing)} test(s) after deletion:")
                for tid in missing:
                    print(f"    - {tid}")
            if extra:
                print(f"  Unexpected {len(extra)} new test(s) after deletion:")
                for tid in extra:
                    print(f"    + {tid}")
            return 1

        # Check coverage match.
        if baseline_cov is not None and post_cov is not None:
            if post_cov < baseline_cov:
                print(
                    f"  [FAIL] Coverage dropped: {baseline_cov:.2f}% -> "
                    f"{post_cov:.2f}% ({summary})"
                )
                return 1
            cov_str = f", coverage {post_cov:.2f}%"
        else:
            cov_str = ""

        print(
            f"  [PASS] All {len(all_deletions)} candidates confirmed "
            f"({summary}, {len(post_ids)} test IDs match baseline{cov_str})"
        )
        confirmed_deletions = all_deletions
    else:
        # Some candidates are false positives. Identify them.
        print(f"  [FAIL] Batch test failed ({summary})")
        print()
        false_positives = _find_false_positives_from_traceback(stderr, all_deletions)
        if not false_positives:
            print("Could not identify specific false positives from the " "traceback.")
            print("Full stderr from the test run:")
            print()
            print(stderr)
            return 1

        n_confirmed = len(all_deletions) - len(false_positives)
        print(
            f"Traceback analysis identified {len(false_positives)} false "
            f"positive(s):"
        )
        print()
        for _, _, _, label in all_deletions:
            status = "FAIL" if label in false_positives else "PASS"
            print(f"  [{status}] {label}")

        confirmed_deletions = [d for d in all_deletions if d[3] not in false_positives]

        # Phase 2: re-batch with only confirmed candidates, verify.
        print()
        print(
            f"Phase 2: Batch-deleting {n_confirmed} confirmed candidates and "
            f"re-running the full test suite with coverage...",
            flush=True,
        )
        batch = [(fp, sl, el) for fp, sl, el, _ in confirmed_deletions]
        originals = _batch_delete(batch)
        try:
            ok, summary, stderr, post_cov = _run_full_test_suite()
            post_ids = _discover_test_ids() if ok else []
        finally:
            _restore_files(originals)

        if ok:
            missing = sorted(set(baseline_ids) - set(post_ids))
            extra = sorted(set(post_ids) - set(baseline_ids))
            if missing or extra:
                print(f"  [FAIL] Tests passed but test ID mismatch ({summary})")
                if missing:
                    print(f"  Missing {len(missing)} test(s) after deletion:")
                    for tid in missing:
                        print(f"    - {tid}")
                if extra:
                    print(f"  Unexpected {len(extra)} new test(s) after deletion:")
                    for tid in extra:
                        print(f"    + {tid}")
                return 1

            if baseline_cov is not None and post_cov is not None:
                if post_cov < baseline_cov:
                    print(
                        f"  [FAIL] Coverage dropped: {baseline_cov:.2f}% -> "
                        f"{post_cov:.2f}% ({summary})"
                    )
                    return 1
                cov_str = f", coverage {post_cov:.2f}%"
            else:
                cov_str = ""

            print(
                f"  [PASS] All {n_confirmed} confirmed candidates pass "
                f"({summary}, {len(post_ids)} test IDs match baseline"
                f"{cov_str})"
            )
        else:
            print(f"  [FAIL] Confirmed-only batch also failed ({summary})")
            print("There may be interaction effects. Manual investigation " "required.")
            return 1

    # Phase 3: with confirmed deletions applied, re-run static analysis to find newly
    # exposed candidates.
    print()
    print("=" * 80)
    print("RE-SCAN FOR NEWLY EXPOSED CANDIDATES")
    print("=" * 80)
    print(flush=True)

    batch = [(fp, sl, el) for fp, sl, el, _ in confirmed_deletions]
    originals = _batch_delete(batch)
    try:
        new_fixture_functions = find_fixture_functions()
        new_all_files = find_all_python_files()
        new_call_sites = find_fixture_call_sites(new_all_files, new_fixture_functions)
        new_transitively_unused = find_transitively_unused(
            new_fixture_functions, new_call_sites
        )
        new_dead_attrs = find_dead_setup_attributes(new_all_files)
    finally:
        _restore_files(originals)

    # Filter to only genuinely new candidates (not already in confirmed list).
    confirmed_keys = set()
    confirmed_attr_keys = set()
    for fp, sl, el, label in confirmed_deletions:
        parts = label.split()[-1]
        if "." in parts:
            confirmed_attr_keys.add(parts)
        else:
            confirmed_keys.add((fp, parts))

    newly_exposed_fixtures = sorted(new_transitively_unused - confirmed_keys)
    newly_exposed_attrs = [
        (fp, cls, attr, ln, eln)
        for fp, cls, attr, ln, eln in new_dead_attrs
        if f"{cls}.{attr}" not in confirmed_attr_keys
    ]

    has_newly_exposed = bool(newly_exposed_fixtures or newly_exposed_attrs)

    if has_newly_exposed:
        print("Removing confirmed candidates exposed new unused items:")
        print()
        if newly_exposed_fixtures:
            print("  New unused fixture functions:")
            for key in newly_exposed_fixtures:
                rel, name = _format_key(key)
                lineno, _ = new_fixture_functions[key]
                print(f"    {rel}:{lineno}  {name}")
        if newly_exposed_attrs:
            print("  New dead setUp/setUpClass attributes:")
            for filepath, class_name, attr_name, lineno, _ in newly_exposed_attrs:
                rel = Path(filepath).relative_to(PROJECT_ROOT)
                print(f"    {rel}:{lineno}  {class_name}.{attr_name}")
    else:
        print("No newly exposed candidates found. The confirmed list is complete.")

    # If --delete-verified, permanently apply the confirmed deletions.
    if args.delete_verified:
        print()
        print("=" * 80)
        print("DELETING CONFIRMED CANDIDATES")
        print("=" * 80)
        print()

        batch = [(fp, sl, el) for fp, sl, el, _ in confirmed_deletions]
        _batch_delete(batch)

        print(f"Deleted {len(confirmed_deletions)} confirmed candidate(s):")
        for _, _, _, label in confirmed_deletions:
            print(f"  {label}")

        if has_newly_exposed:
            n_newly = len(newly_exposed_fixtures) + len(newly_exposed_attrs)
            print()
            print(
                f"Re-run this script to verify and delete the {n_newly} "
                f"newly exposed candidate(s)."
            )

    print()
    total_confirmed = len(confirmed_deletions)
    total_failed = len(all_deletions) - total_confirmed
    n_newly = len(newly_exposed_fixtures) + len(newly_exposed_attrs)
    print(
        f"Final summary: {total_confirmed} confirmed, {total_failed} false "
        f"positive(s), {n_newly} newly exposed."
    )

    return 1 if total_issues > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
