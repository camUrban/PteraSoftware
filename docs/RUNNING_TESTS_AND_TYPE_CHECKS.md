# Running Tests and Type Checks
Guidelines for Claude on running tests and type checks for Ptera Software. For developers, see the `CONTRIBUTING.md` file for more instructions.

## Running Tests

```
Bash(.venv/Scripts/python.exe:-m:unittest:discover:-s:${WORKSPACE}/tests:*)
```

If that fails, try:

```
Bash(.venv/bin/python:-m:unittest:discover:-s:${WORKSPACE}/tests:*)
```

To run a specific test module in `tests\unit\`:

```
".venv/Scripts/python.exe" -m unittest tests.unit.test_module_name -v
```

For example:

```
".venv/Scripts/python.exe" -m unittest tests.unit.test_wing_cross_section_movement -v
```

## Running MyPy Type Checking

```bash
cd ${WORKSPACE} && ".venv/Scripts/python.exe" -m mypy pterasoftware
```