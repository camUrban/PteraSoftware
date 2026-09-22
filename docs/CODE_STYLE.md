# Code Style

## Code Formatting

- Follow existing code style (black) and conventions
- Maintain consistent indentation and spacing
- When initializing NumPy arrays (e.g., using `np.zeros`, `np.ones`, `np.empty`, etc.), always specify the `dtype` parameter (e.g., `dtype=float`, `dtype=int`, `dtype=object`, etc.)
- In runtime strings such as error and log messages, name str values in double quotes, as the docstring rules in [Type Hints and Docstrings](TYPE_HINT_AND_DOCSTRING_STYLE.md) do. Black then keeps the outer string single-quoted to avoid escapes (`f'force_method must be "joukowski" or "katz", got "{force_method}".'`). Write the value inline with double quotes rather than through `!r`, which always emits single quotes.

## Variable Naming

- Use descriptive variable names that clearly indicate their purpose
- Use underscores for variable names
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the [Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md) and [Angle Vectors and Transformations](ANGLE_VECTORS_AND_TRANSFORMATIONS.md) documents when naming vector-valued variables or things such as transformation and rotation matrices.
- Do not use `wcs` (or any other abbreviation) for "wing cross section" or `WingCrossSection` in variable names. Instead, always write it out in full (e.g., `wing_cross_section`, `wing_cross_section_movement`, etc.). Exception: a capitalization variant of `wcs` is allowed when the name is written in exactly one of the forms required by the [Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md) or [Angle Vectors and Transformations](ANGLE_VECTORS_AND_TRANSFORMATIONS.md) documents (e.g., the axes IDs in `angles_Wcsp_to_Wcs_ixyz` and the point IDs in `Lp_Wcsp_Lpp`).
- For non-vector-valued variables whose values are in radians-based units (rad, rad/s, rad/s^2, etc.), append `_rad` to the variable name. Noting the units in a docstring or comment alone is not enough. Vector-valued variables and their components instead mark radians with `Rad` in their camelCase portions, as described in the [Angle Vectors and Transformations](ANGLE_VECTORS_AND_TRANSFORMATIONS.md) document.
- When naming collections of vector-valued quantities, use the `stack`, `grid`, and `list` camelCase prefixes as described in the [Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md) document.
- Name module-level constants in UPPER_SNAKE_CASE, and keep the leading underscore on private names (e.g., `_SEED`, `_FOUR_PI`, `PACKAGE_LOGGER_NAME`). There are three exceptions. First, module-level names that hold state rather than a constant value stay lowercase (e.g., `_logger`, `_solve_loop_lock`, and `_indent_level`). Second, a constant whose name is written in one of the forms required by the [Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md) or [Angle Vectors and Transformations](ANGLE_VECTORS_AND_TRANSFORMATIONS.md) documents keeps that spelling (e.g., `_freeFlightViewDirection_E`). Third, a module-level name that an external tool reads by a fixed name keeps that name (e.g., the Sphinx configuration values in `docs/website/conf.py`, such as `project` and `html_theme`).
- In variable names, refer to a simulation time increment as a step (e.g., `step`, `num_steps`, `step_discards`). `time_step` is also acceptable, but never use `timestep`.

## Running Black

Black is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files black
```

## Running codespell

codespell is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files codespell
```

## Running mypy

mypy is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files mypy
```

## Imports

- Import Ptera Software using the following pattern: `import pterasoftware as ps`
- By default, place import statements at the top of the file and avoid imports inside functions or methods. The only exceptions are intentional lazy-import patterns (for example, using `importlib.import_module` inside `__getattr__` for lazy loading) and cases where there is no other way to avoid circular imports.

## Miscellaneous Guidelines

- Use `np.deg2rad` and `np.rad2deg` for angle conversions instead of `np.radians` and `np.degrees` or manual conversions.
