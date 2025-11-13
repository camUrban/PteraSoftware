# Type Hint and Docstring Style Guide

This document defines the conventions for type hints and docstrings in the Ptera Software codebase.

## Table of Contents
- [Type Hints](#type-hints)
- [Docstring Format](#docstring-format)
  - [Module-Level Docstrings](#module-level-docstrings)
  - [Class Docstrings](#class-docstrings)
  - [Function and Method Docstrings](#function-and-method-docstrings)
- [Examples](#examples)

---

## Type Hints

### General Principles

1. **Type hints should reflect what functions ACCEPT, not just what they store internally**
2. **Always include return type hints** (use `-> None` for functions that don't return values)
3. **Shape information belongs in docstrings**, not type hints (Python's type system can't express array shapes)

### Import Requirements

```python
from collections.abc import Sequence

import numpy as np
```

### Type Hint Patterns by Parameter Type

#### Basic Types

| Parameter Description                 | Type Hint           |
|---------------------------------------|---------------------|
| String                                | `str`               |
| Boolean                               | `bool`              |
| Boolean (accepting numpy bools)       | `bool \| np.bool_`  |
| Integer                               | `int`               |
| Number (int or float)                 | `float \| int`      |
| Float only                            | `float`             |

#### Array and Array-Like Types

| Parameter Description                       | Type Hint                                        | Notes                                                      |
|---------------------------------------------|--------------------------------------------------|------------------------------------------------------------|
| Array-like accepting numbers (int or float) | `np.ndarray \| Sequence[float \| int]`           | For user-facing parameters that accept tuples/lists/arrays |
| Array-like accepting floats only            | `np.ndarray \| Sequence[float]`                  | When only floats are valid                                 |
| Array of arrays (2D array-like)             | `np.ndarray \| Sequence[Sequence[float \| int]]` | For nested sequences like `[[1,2], [3,4]]`                 |
| Already a numpy array                       | `np.ndarray`                                     | For internal functions or returns                          |

#### Class Types

| Parameter Description              | Type Hint                | Notes                                      |
|------------------------------------|--------------------------|--------------------------------------------|
| Class from same package            | `ClassName`              | Direct reference                           |
| Class from imported module         | `module_alias.ClassName` | Use module alias to avoid circular imports |
| Self-reference (forward reference) | `"ClassName"`            | Use string quotes for return types         |

#### Optional and Union Types

| Parameter Description   | Type Hint                 |
|-------------------------|---------------------------|
| Optional parameter      | `Type \| None`            |
| Union of multiple types | `Type1 \| Type2 \| Type3` |

### Module Alias Pattern

Import modules with aliases:

```python
from . import wing as wing_mod
from . import airfoil as airfoil_mod
from . import wing_cross_section as wing_cross_section_mod

# In function signature
def mesh_wing(wing: wing_mod.Wing) -> None:
    ...

def _get_mcl_points(
    inner_airfoil: airfoil_mod.Airfoil,
    outer_airfoil: airfoil_mod.Airfoil,
    ...
) -> list[np.ndarray]:
    ...
```

### Avoiding Circular Imports with Type Hints

#### Preferred Method: `from __future__ import annotations`

To avoid circular import errors when type hinting, use `from __future__ import annotations` as the **first import** in your module. This defers evaluation of type hints, treating them as strings automatically:

```python
from __future__ import annotations

from collections.abc import Sequence
import numpy as np

from .. import geometry
from . import wing_movement as wing_movement_mod

# In function signature - no quotes needed!
def __init__(
    self,
    base_airplane: geometry.airplane.Airplane,
    wing_movements: list[wing_movement_mod.WingMovement],
) -> None:
    ...
```

This approach:
- Keeps all imports at the top of the file
- Prevents circular import errors
- Requires no string quotes around type hints
- Is the default behavior in Python 3.11+

---

## Docstring Format

### General Principles

1. **Use reStructuredText (rST) format**
2. **Short description starts on same line as triple quotes**
3. **Parameter descriptions are inline** (no separate type line)
4. **Begin descriptions with article + shape/type info for arrays** (e.g., "A (4,4) ndarray of floats...")
5. **Use present tense for descriptions** (e.g., "Returns..." not "Will return...")
6. **Avoid starting descriptions with "This..."**
7. **Never use em-dashes (—) or en-dashes (–); always use hyphens (-)**
8. **Never use multiplication sign (×); use lowercase x**
9. **Never use pi symbol (π); write "pi"**
10. **Never use approximately-equal sign (≈); use "~"**
11. **Place closing triple-quotes on their own line**

### Module-Level Docstrings

Module-level docstrings appear at the very top of each Python file and describe the module's contents. The style varies based on the type of module.

#### Public Package __init__.py Files

`__init__.py` files for a public package list subpackages, directories, and modules:

```python
"""Contains the geometry classes.

**Contains the following subpackages:**

None

**Contains the following directories:**

None

**Contains the following modules:**

airfoil.py: Contains the Airfoil class.

airplane.py: Contains the Airplane class.

wing.py: Contains the Wing class.

wing_cross_section.py: Contains the WingCrossSection class.
"""
```

**Pattern:**
- Brief description using "Contains" (present tense)
- List public subpackages (or "None")
- List public directories (or "None")
- List public modules with one-line descriptions
- Use blank lines between subpackages/directories/module entries for readability

#### Public Modules

Public modules (e.g., `airfoil.py`, `wing.py`) have structured docstrings listing their contents:

```python
"""Contains the Airfoil class.

**Contains the following classes:**

Airfoil: A class used to contain the Airfoil of a WingCrossSection.

**Contains the following functions:**

None
"""
```

**Pattern:**
- Brief description using "Contains" (present tense)
- List public classes with brief descriptions (use "A class used to..." or similar)
- List public functions (or "None")
- Use blank lines between class/function entries if there are multiple

#### Private Modules

Private modules (e.g., `_meshing.py`, `_functions.py`) have minimal docstrings:

```python
"""Contains the function for meshing Wings."""
```

**Pattern:**
- Single brief sentence
- No listing of functions or classes
- Keep it concise since these are internal implementation details

### Function and Method Docstrings

```python
def function_name(
    param1: Type1,
    param2: Type2,
    param3: Type3,
) -> ReturnType:
    """Short description of what the function does.

    Optional longer description providing more context. This provides detailed explanations of the function's behavior. It can be one or more paragraphs.
    
    Optional citation block:

    **Citation:**
    
    Adapted from (can be more specific if the whole function wasn't adapted): <source>
    
    Author (or "Authors"): <author>
    
    Date of retrieval (don't include if not known): <date>

    :param param1: A (shape) dtype description of param1. Additional details about
        what it represents, valid ranges, units, etc. Can wrap to multiple lines.
    :param param2: Description of param2.
    :param param3: Description of param3.
    :return: A (shape) dtype description of what is returned. Additional details
        about the return value.
    """
```

### Array Parameter Descriptions

For numpy arrays, always include:
1. **Shape**: e.g., "(4,4)", "(M,N,3)", "(N,)"
2. **Dtype**: e.g., "ndarray of floats", "ndarray of ints", "ndarray of bools"
3. **Coordinate system and reference point** (when applicable)
4. **Units** (when applicable)
5. **Default value** (when applicable)

#### Pattern for Array Parameters

```python
:param parameter_name: A (shape) ndarray of dtype representing <description>.
    Additional context about coordinate systems, valid ranges, units, default value, etc.
```

#### Pattern for Array-Like Parameters

```python
:param parameter_name: An array-like object of numbers (int or float) with shape
    (N,M) representing <description>. Can be a tuple, list, or ndarray. Values are
    converted to floats internally. The units are <units>. The default is <default>.
```

### Class Docstrings

```python
class ClassName:
    """Short description of the class.

    **Contains the following methods:**
        public_method_1: Short description (identical to method's docstring's short 
        description.

        public_method_2: Short description (identical to method's docstring's short 
        description.

    Optional notes block
    
    **Notes:**
    
    Detailed description of the class's purpose, behavior, or usage. Can be one or more paragraphs. Avoid numbered or bulleted lists.
    
    Optional citation block:
    
    **Citation:**
    
    Adapted from (can be more specific if the whole class wasn't adapted): <source>
    
    Author (or "Authors"): <author>
    
    Date of retrieval (don't include if not known): <date>
    """
```

### Property Docstring Template

```python
@property
def property_name(self) -> ReturnType:
    """Short description of what the property represents.

    :return: Description of what is returned, including type, shape, units.
    """
```

---

## Examples

### Example 1: Module-Level Docstrings

#### Public Package __init__.py
```python
"""Contains the geometry classes.

**Contains the following subpackages:**

None

**Contains the following directories:**

None

**Contains the following modules:**

airfoil.py: Contains the Airfoil class.

airplane.py: Contains the Airplane class.

wing.py: Contains the Wing class.

wing_cross_section.py: Contains the WingCrossSection class.
"""
```

#### Public Module
```python
"""Contains the Airfoil class.

**Contains the following classes:**

Airfoil: A class used to contain the Airfoil of a WingCrossSection.

**Contains the following functions:**

None
"""
```

#### Private Module
```python
"""Contains the function for meshing Wings."""
```

### Example 2: Function with Array Parameters (Internal)

```python
def _get_mcl_points(
    inner_airfoil: airfoil_mod.Airfoil,
    outer_airfoil: airfoil_mod.Airfoil,
    chordwise_coordinates: np.ndarray,
) -> list[np.ndarray]:
    """Takes in the inner and outer Airfoils of a wing section and its normalized 
    chordwise coordinates. It returns a list of four column vectors containing the 
    normalized components of the positions of points along the mean camber line (MCL) 
    (in each Airfoil's axes, relative to each Airfoil's leading point).

    :param inner_airfoil: The wing section's inner Airfoil.
    :param outer_airfoil: The wing section's outer Airfoil.
    :param chordwise_coordinates: A (N,) ndarray of floats for the normalized
        chordwise coordinates where we'd like to sample each Airfoil's MCL. The values
        are normalized from 0.0 to 1.0 and are unitless.
    :return: A list of four (N,1) ndarrays of floats, where N is the number of points
        at which we'd like to sample each Airfoil's MCL. The ndarrays contain components
        of the positions of points along each Airfoil's MCL. In order, the ndarrays
        returned are, (1) the inner Airfoil's MCL points' y-components, (2) the inner
        Airfoil's MCL points' x-components (3) the outer Airfoil's MCL points'
        y-components, and (4) the outer Airfoil's MCL points' x-components. The values
        are normalized from 0.0 to 1.0 and are unitless.
    """
```

### Example 2: Function with Transformation Matrices

```python
def _get_mcs_points(
    T_pas_Wcsi_Lpi_Wn_Ler: np.ndarray,
    T_pas_Wcso_Lpo_Wn_Ler: np.ndarray,
    inner_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    outer_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    inner_mcl_pointsY_Ai_lpAi: np.ndarray,
    inner_mcl_pointsX_Ai_lpAi: np.ndarray,
    outer_mcl_pointsY_Ao_lpAo: np.ndarray,
    outer_mcl_pointsX_Ao_lpAo: np.ndarray,
    spanwise_coordinates: np.ndarray,
) -> list[np.ndarray]:
    """Calculates the points on a wing section's mean camber surface (MCS) (in wing 
    axes, relative to the leading edge root point).

    :param T_pas_Wcsi_Lpi_Wn_Ler: A (4,4) ndarray of floats representing a passive
        transformation matrix which maps in homogeneous coordinates from the inner
        WingCrossSection's axes, relative to its leading point to wing axes relative to
        the leading edge root point.
    :param T_pas_Wcso_Lpo_Wn_Ler: A (4,4) ndarray of floats representing a passive
        transformation matrix which maps in homogeneous coordinates from the outer
        WingCrossSection's axes, relative to its leading point to wing axes relative to
        the leading edge root point.
    :param inner_wing_cross_section: The wing section's inner WingCrossSection.
    :param outer_wing_cross_section: The wing section's outer WingCrossSection.
    :param inner_mcl_pointsY_Ai_lpAi: A (M,1) ndarray of floats, where M is the
        number of chordwise points in the mesh. Each element represents the y-component
        of the inner Airfoil's MCL points (in the inner Airfoil's axes, relative to the
        inner Airfoil's leading point). The values are normalized from 0.0 to 1.0 and
        are unitless.
    :return: A list of four (M,N,3) ndarrays of floats, where M is the number of
        chordwise points and N is the number of spanwise points. The four ndarrays are,
        in order, this wing section's Panel's (1) forward inner, (2) forward outer,
        (3) backward inner, and (4) backward outer panel points (in wing axes, relative
        to the leading edge root point). The units are in meters.
    """
```

### Example 3: Public Method with Array-Like Parameters

```python
def __init__(
    self,
    name: str = "NACA0012",
    outline_A_lp: np.ndarray | Sequence[Sequence[float | int]] | None = None,
    resample: bool = True,
    n_points_per_side: int = 400,
) -> None:
    """The initialization method.

    :param name: The name of the Airfoil. It should correspond to the name of a file
        the airfoils directory, or to a valid NACA 4-series airfoil (once converted to
        lower-case and stripped of leading and trailing whitespace) unless you are
        passing in your own array of points using outline_A_lp. Note that NACA0000 isn't
        a valid NACA-series airfoil. The default is "NACA0012".
    :param outline_A_lp: An array-like object of numbers (int or float) with shape
        (N,2) representing the 2D points making up the Airfoil's outline (in airfoil
        axes, relative to the leading point). If you wish to load coordinates from the
        airfoils directory, leave this as None, which is the default. Can be a tuple,
        list, or ndarray. Values are converted to floats internally. Make sure all
        x-component values are in the range [0.0, 1.0]. The default value is None.
    :param resample: Determines whether to resample the points defining the Airfoil's
        outline. This applies to points passed in by the user or to those from the
        airfoils directory. I highly recommended setting this to True. Can be a bool or
        a numpy bool and will be converted internally to a bool. The default is True.
    :param n_points_per_side: The number of points to use when creating the Airfoil's
        MCL and when resampling the upper and lower parts of the Airfoil's outline. It
        must be a positive int greater than or equal to 3. The resampled outline will
        have a total number of points equal to (2 * n_points_per_side) - 1. I highly
        recommend setting this to at least 100. The default value is 400.
    :return: None
    """
```

### Example 4: Method Returning Self-Reference

```python
def add_control_surface(
    self, deflection: float | int, hinge_point: float | int
) -> "Airfoil":
    """Returns a version of the Airfoil with a control surface added at a given point. 
    It is called during meshing.

    :param deflection: The control deflection in degrees. Deflection downwards is
        positive. It must be a number (int or float) in the range [-5.0, 5.0] degrees.
        Values are converted to floats internally.
    :param hinge_point: The location of the hinge as a fraction of chord length. It
        must be a number (int or float) in the range (0.0, 1.0). Values are converted
        to floats internally.
    :return: The new Airfoil with the control surface added.
    """
```

### Example 5: Method with Optional Return

```python
def get_plottable_data(self, show: bool = False) -> list[np.ndarray] | None:
    """Returns plottable data for this Airfoil's outline and mean camber line.

    :param show: Determines whether to display the plot. Can be a bool or a numpy bool,
        and will be converted internally to a bool. If True, the method displays the
        plot and returns None. If False, the method returns the data without displaying.
        The default is False.
    :return: A list of two ndarrays containing the outline and MCL data, or None if
        show is True.
    """
```

### Example 6: Method Returning Array

```python
def get_resampled_mcl(
    self, mcl_fractions: np.ndarray | Sequence[float]
) -> np.ndarray:
    """Returns a ndarray of points along the mean camber line (MCL), resampled from the 
    mcl_A_outline attribute. It is used to discretize the MCL for meshing.

    :param mcl_fractions: A (N,) array-like object of floats representing normalized
        distances along the MCL (from the leading to the trailing edge) at which to
        return the resampled MCL points. Can be a tuple, list, or ndarray. The first
        value must be 0.0, the last must be 1.0, and the remaining must be in the range
        [0.0, 1.0]. All values must be non-duplicated and in ascending order.
    :return: A (N,2) ndarray of floats that contains the positions of the resampled
        MCL points (in airfoil axes, relative to the leading point).
    """
```

---

## Quick Reference

### Common Type Hint Patterns

```python
# Simple types
param: str
param: bool
param: bool | np.bool_  # Accepts both Python and NumPy booleans
param: int
param: float | int

# Array-like (user input)
param: np.ndarray | Sequence[float | int]
param: np.ndarray | Sequence[float]
param: np.ndarray | Sequence[Sequence[float | int]]

# Already numpy arrays
param: np.ndarray

# Classes
param: ClassName  # Same module
param: module_alias.ClassName  # Different module
-> "ClassName"  # Self-reference

# Optional/Union
param: Type | None
param: Type1 | Type2

# Collections
-> list[np.ndarray]
-> list[ClassName]
```

### Common Docstring Phrases

```python
# Module level
"Contains the <description>."
"Contains the following subpackages:"
"Contains the following directories:"
"Contains the following modules:"
"Contains the following classes:"
"Contains the following functions:"

# Array parameters
":param name: A (shape) ndarray of dtype representing..."
":param name: An array-like object of numbers (int or float) with shape..."

# Boolean parameters (accepting numpy bools)
":param name: A bool that... Can be a bool or a numpy bool and will be converted internally to a bool."

# Return values
":return: A (shape) ndarray of dtype that..."
":return: A list of N (shape) ndarrays..."
":return: None"

# Function descriptions (avoid "This function/method")
"Takes in... and returns..."
"Calculates..."
"Returns..."
"Validates..."

# Units and ranges
"The units are meters."
"The values are normalized from 0.0 to 1.0 and are unitless."
"It must be in the range [0.0, 1.0]."
"It must be a positive int."

# Coordinate systems
"(in geometry axes, relative to the CG)"
"(in wing axes, relative to the leading edge root point)"
"(in airfoil axes, relative to the leading point)"
```

---

## Notes

- This style guide should be updated as new patterns emerge
- All existing code should gradually be updated to match this style
- Use `docformatter` or similar tools to help maintain consistent formatting
- Shape information is critical and must always be included in docstrings for arrays