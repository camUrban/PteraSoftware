# Code Style

## Code Formatting

- Follow existing code style (black) and conventions
- Maintain consistent indentation and spacing
- When initializing numpy arrays (e.g., using `np.zeros`, `np.ones`, `np.empty`, etc.), always specify the `dtype` parameter (e.g., `dtype=float`, `dtype=int`, `dtype=object`, etc.)

## Variable Naming

- Use descriptive variable names that clearly indicate their purpose
- Use underscores for variable names
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the `AXES_AND_COORDINATE_SYSTEMS.md` and `AXES_POINTS_AND_FRAMES.md` documents when naming vector-valued variables or things such as transformation and rotation matrices.
- Do not use `wcs` (or any other abbreviation) for "wing cross section" or "WingCrossSection" in variable names. Instead, always write it out in full (e.g., `wing_cross_section`, `wing_cross_section_movement`, etc.).

## Formatting with Black

```bash
cd ${WORKSPACE} && ".venv/Scripts/python.exe" black .
```

## Running a CodeSpell Spell Check

```bash
cd ${WORKSPACE} && ".venv/Scripts/codespell.exe" --ignore-words=.codespell-ignore.txt --skip="*/_build/*,*.dat"
```

## Imports

- Import Ptera Software using the following pattern: ```import pterasoftware as ps```
- By default, place import statements at the top of the file and avoid imports inside functions or methods. The only exceptions are intentional lazy-import patterns (for example, using `importlib.import_module` inside `__getattr__` for lazy loading) and cases where there is no other way to avoid circular imports.

## Miscellaneous Guidelines

- Use `np.deg2rad` and `np.rad2deg` for angle conversions instead of `np.radians` and `np.degrees` or manual conversions.