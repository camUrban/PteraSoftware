# Code Style

## Code Formatting

- Follow existing code style (black) and conventions
- Maintain consistent indentation and spacing
- When initializing numpy arrays (e.g., using `np.zeros`, `np.ones`, `np.empty`, etc.), always specify the `dtype` parameter (e.g., `dtype=float`, `dtype=int`, `dtype=object`, etc.)

## Variable Naming

- Use descriptive variable names that clearly indicate their purpose
- Use underscores for variable names
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the `AXES_POINTS_AND_FRAMES.md` and `ANGLE_VECTORS_AND_TRANSFORMATIONS.md` documents when naming vector-valued variables or things such as transformation and rotation matrices.
- Do not use `wcs` (or any other abbreviation) for "wing cross section" or "WingCrossSection" in variable names. Instead, always write it out in full (e.g., `wing_cross_section`, `wing_cross_section_movement`, etc.).
- For non-vector-valued variables whose values are in radians-based units (rad, rad/s, rad/s^2, etc.), append `_rad` to the variable name; noting the units in a docstring or comment alone is not enough. Vector-valued variables and their components instead mark radians with `Rad` in their camelCase portions, as described in the `ANGLE_VECTORS_AND_TRANSFORMATIONS.md` document.
- When naming collections of vector-valued quantities, use the `stack`, `grid`, and `list` camelCase prefixes as described in the `AXES_POINTS_AND_FRAMES.md` document.
- In variable names, refer to a simulation time increment as a step (e.g., `step`, `num_steps`, `step_discards`); `time_step` is also acceptable, but never use `timestep`.

## Running Black

Black is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files black
```

## Running CodeSpell

CodeSpell is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files codespell
```

## Running mypy

mypy is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files mypy
```

## Imports

- Import Ptera Software using the following pattern: ```import pterasoftware as ps```
- By default, place import statements at the top of the file and avoid imports inside functions or methods. The only exceptions are intentional lazy-import patterns (for example, using `importlib.import_module` inside `__getattr__` for lazy loading) and cases where there is no other way to avoid circular imports.

## Miscellaneous Guidelines

- Use `np.deg2rad` and `np.rad2deg` for angle conversions instead of `np.radians` and `np.degrees` or manual conversions.
