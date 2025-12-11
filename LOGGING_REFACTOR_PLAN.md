# Logging Refactor Plan

## Summary

Refactor the logging system to follow Python best practices: libraries should not configure logging, only emit log records. Users configure logging themselves.

## Changes Required

### 1. Expose `setup_logging` as public API

Add to `pterasoftware/__init__.py`:

```python
from pterasoftware._logging import setup_logging
```

Update the module docstring to document this function.

### 2. Remove `logging_level` parameter from all solver `run()` methods

Files to modify:
- `pterasoftware/steady_horseshoe_vortex_lattice_method.py`
- `pterasoftware/steady_ring_vortex_lattice_method.py`
- `pterasoftware/unsteady_ring_vortex_lattice_method.py`

In each file:
- Remove `logging_level` parameter from `run()` signature
- Remove `logging_level` from docstring
- Remove the `ensure_logging_configured()` call

### 3. Always show progress bar in unsteady solver

In `unsteady_ring_vortex_lattice_method.py`, change the progress bar `disable` parameter:

```python
# Before
disable=logging_level_value == logging.DEBUG,

# After
disable=False,
```

Also remove the `logging_level_value` variable since it's no longer needed.

### 4. Clean up `_logging.py`

The `ensure_logging_configured()` function can be removed since it's no longer called by solvers. Keep:
- `setup_logging()`
- `get_logger()`
- `_TqdmLoggingHandler`
- `convert_logging_level_name_to_value()`

### 5. Update tests

- Remove/update any tests that use the `logging_level` parameter
- Remove tests for `ensure_logging_configured()` if that function is removed
- Keep the test we updated: `test_preserves_user_level_if_already_configured`

### 6. Update example file

In `examples/steady_horseshoe_vortex_lattice_method_solver.py`, change:

```python
# Before
from pterasoftware import _logging
_logging.setup_logging(level="Debug")

# After
import pterasoftware as ps
ps.setup_logging(level="Debug")
```

## Design Rationale

- **Libraries should not configure logging** - only emit records, let users configure
- **Progress bar shown by default** - if users configure their own logging (not using `ps.setup_logging()`) and get progress bar corruption, that's documented behavior
- **`ps.setup_logging()` provides TQDM-compatible handler** - users who want both logging AND progress bars should use this function
- **Breaking change** - removing `logging_level` parameter, but results in cleaner API

## Already Completed

- Fixed bug in `ensure_logging_configured()` that was overwriting user's log level
- Updated test `test_preserves_user_level_if_already_configured`