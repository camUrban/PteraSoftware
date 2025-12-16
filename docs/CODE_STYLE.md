# Code Style Guidelines
Guidelines for developers and Claude when writing code.

## Code Formatting
- Follow existing code style (black) and conventions
- Maintain consistent indentation and spacing
- When initializing numpy arrays (e.g., using `np.zeros`, `np.ones`, `np.empty`, etc.), always specify the `dtype` parameter (e.g., `dtype=float`, `dtype=int`, `dtype=object`, etc.)

## Variable Naming
- Use descriptive variable names that clearly indicate their purpose  
- Use underscores for variable names
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the `AXES_AND_COORDINATE_SYSTEMS.md` and `AXES_POINTS_AND_FRAMES.md` documents when naming vector-valued variables or things such as transformation and rotation matrices.

## Formatting with Black
For Claude, use the following command. For developers, see the `CONTRIBUTING.md` file for instructions.

```bash
cd ${WORKSPACE} && ".venv/Scripts/python.exe" black .
```

## Miscellaneous Guidelines
- Import Ptera Software using the following pattern: ```import pterasoftware as ps```
- Use `np.deg2rad` and `np.rad2deg` for angle conversions instead of `np.radians` and `np.degrees` or manual conversions.