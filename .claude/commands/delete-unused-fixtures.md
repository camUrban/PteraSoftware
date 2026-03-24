---
description: Find and delete unused fixtures from a specified fixture module
---

# Delete Unused Fixtures

Find and delete unused fixtures for `$ARGUMENTS`.

## Important: Cross-Module Usage

Fixtures may be used by tests in ANY module, not just the corresponding test file. For example, `wing_movement.py` has fixtures in `wing_movement_fixtures.py` and tests in `test_wing_movement.py`. However, fixtures from `wing_movement_fixtures.py` might also be used by `test_airplane_movement.py`. You MUST search the entire `tests/` directory for usage.

## Progress Tracking

Follow these steps carefully and track your progress:
- [ ] Read the fixture file for the module specified in `$ARGUMENTS`
- [ ] Extract all fixture function names (functions that don't start with `_`)
- [ ] For each fixture, search the entire `tests/` directory for usage
- [ ] Identify fixtures with zero usages (excluding their own definition)
- [ ] Generate an unused fixtures report
- [ ] Delete unused fixtures from the file
- [ ] Update imports in `tests/unit/fixtures/__init__.py` if needed

## Omissions

Please do NOT:
- Run or debug any tests you create
- Reformat any files with black

These steps are both handled by other slash commands.

## Detailed Steps

1. **Read the fixture file** for the module specified in `$ARGUMENTS`:
    - Note all public function names (fixture functions)
    - Ignore private functions (those starting with `_`)
2. **For each fixture function**, search for usage:
   Use the Grep tool to search the entire `tests/` directory:
   ```
   pattern: fixture_function_name
   path: tests/
   output_mode: files_with_matches
   ```

   A fixture is USED if it appears in any file OTHER than:
    - Its own fixture file (the definition)
    - `__init__.py` files (just imports)
3. **Classify each fixture**:
    - **Used**: Found in at least one test file
    - **Unused**: Only found in its definition and/or `__init__.py` imports
4. **Generate an unused fixtures report**:
   ```
   UNUSED FIXTURES REPORT for [fixture_file]

   UNUSED (will be deleted):
   - [fixture_name]: no usages found in test files

   USED (keeping):
   - [fixture_name]: used in [list of test files]

   SUMMARY:
   - Total fixtures: [N]
   - Used: [N]
   - Unused: [N]
   ```
5. **Delete unused fixtures**:
    - Remove the function definitions from the fixture file
    - Remove any imports that were only used by deleted fixtures
    - Keep the file structure clean (no orphaned blank lines)
6. **Update `tests/unit/fixtures/__init__.py`**:
    - Remove imports for deleted fixtures
    - Update the module docstring if needed

## Efficiency Tips

- Use Grep with `output_mode: files_with_matches` for fast searching
- Search for multiple fixtures in parallel when possible
- Don't read test files unless you need to verify ambiguous matches

## Error Handling

- If a fixture name is common (e.g., `valid_wing`), verify the grep matches are actual usages, not coincidental string matches
- If uncertain whether a fixture is used, keep it and flag for manual review
- If the fixture file doesn't exist, report the error and stop

## Quality Checklist

Before finalizing:
- [ ] Searched the ENTIRE `tests/` directory, not just corresponding test file
- [ ] Verified each "unused" fixture truly has no usages
- [ ] Generated the unused fixtures report
- [ ] Removed only unused fixtures
- [ ] Updated `__init__.py` imports
- [ ] No useful fixtures were accidentally deleted
