# Type Hints and Docstrings

This document defines the conventions for type hints and docstrings in the Ptera Software codebase.

## Table of Contents

- [Type Hints](#type-hints)
    - [Type Hints in Tests](#type-hints-in-tests)
- [Docstring Format](#docstring-format)
    - [Module-Level Docstrings](#module-level-docstrings)
    - [Class Docstrings](#class-docstrings)
    - [Public Subclasses of Private Parents](#public-subclasses-of-private-parents)
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

| Parameter Description           | Type Hint         |
|---------------------------------|-------------------|
| String                          | `str`             |
| Boolean                         | `bool`            |
| Boolean (accepting numpy bools) | `bool \| np.bool` |
| Integer                         | `int`             |
| Number (int or float)           | `float \| int`    |
| Float only                      | `float`           |

#### Array and Array-Like Types

| Parameter Description                       | Type Hint                                        | Notes                                                      |
|---------------------------------------------|--------------------------------------------------|------------------------------------------------------------|
| Array-like accepting numbers (int or float) | `np.ndarray \| Sequence[float \| int]`           | For user-facing parameters that accept tuples/lists/arrays |
| Array-like accepting floats only            | `np.ndarray \| Sequence[float]`                  | When only floats are valid                                 |
| Array of arrays (2D array-like)             | `np.ndarray \| Sequence[Sequence[float \| int]]` | For nested sequences like `[[1,2], [3,4]]`                 |
| Already a numpy array                       | `np.ndarray`                                     | For internal functions or returns                          |

#### Class Types

| Parameter Description      | Type Hint                | Notes                                      |
|----------------------------|--------------------------|--------------------------------------------|
| Class from same package    | `ClassName`              | Direct reference                           |
| Class from imported module | `module_alias.ClassName` | Use module alias to avoid circular imports |

#### Optional and Union Types

| Parameter Description   | Type Hint                 |
|-------------------------|---------------------------|
| Optional parameter      | `Type \| None`            |
| Union of multiple types | `Type1 \| Type2 \| Type3` |

### Type Narrowing Patterns

When working with attributes that may be `None`, use these patterns:

#### Optional Attributes

For class attributes initialized to `None` but populated later:

```python
class Solver:
    def __init__(self):
        # Attribute that will be populated before use
        self.current_airplanes: list[geometry.airplane.Airplane] | None = None
```

#### Narrowing with Assertions

Use `assert` when you have a programming invariant (the value should never be `None` at this point):

```python
def compute_forces(self):
    # At this point in the code, these should always be populated
    assert self.current_airplanes is not None
    for airplane in self.current_airplanes:
        # mypy now knows current_airplanes is not None
        ...
```

**Use `assert` when:**

- `None` represents a bug, not a valid state
- You want runtime safety during development
- The invariant should always hold

#### Narrowing with `cast()`

Use `cast()` sparingly, only when the type checker cannot infer what you know to be true:

```python
from typing import cast

# For dtype=object arrays where we know the element type
panel = cast(_panel.Panel, object_array[i, j])
```

**Use `cast()` when:**

- Working around type checker limitations (e.g., numpy `dtype=object` arrays)
- You're certain of the type but can't prove it to the type checker
- No runtime check is needed

**Avoid `cast()` for `Type | None` -> `Type` narrowing.** Use `assert` instead for runtime safety.

When the type being cast to lives across a circular import boundary, use the string form (`cast("OtherClass", value)`). See "Casting Across a Circular Dependency" below.

### Module Alias Pattern

Import modules with aliases:

```python
from . import airfoil as airfoil_mod
from . import wing as wing_mod
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

from . import wing_movement as wing_movement_mod

from .. import geometry


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

#### Casting Across a Circular Dependency

`from __future__ import annotations` defers type-hint evaluation but does not help with `cast()`, which is a runtime call whose first argument is evaluated. When narrowing a type that lives across a circular boundary, use the string form of `cast()` so the type name does not need to be importable at runtime:

```python
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from .other_module import OtherClass


def narrow(value):
    other = cast("OtherClass", value)
    ...
```

The `TYPE_CHECKING` import gives mypy the symbol for static resolution. The string argument keeps the runtime call free of any reference to `OtherClass`. Prefer this over importing `OtherClass` inside the function body: in-function imports are reserved for genuine lazy-load or circular cases, and the string-form `cast()` resolves the circularity without that escape hatch, keeping all imports at the top of the file.

### Type Hints in Tests

The test suite is type-checked with the same mypy configuration as the package, including `disallow_untyped_defs`. Test code has recurring situations that package code does not, and this section defines the pragma-free convention for each of them. The project contains no `type: ignore` pragmas anywhere, and none of these situations justifies adding one.

#### Attributes Assigned in `setUpClass`

mypy cannot see attributes assigned through `cls` inside `setUpClass`, so every use site reports an attribute error. Declare the attributes with class-level annotations directly below the class docstring, taking the types from the fixture factories' return annotations. Do not wrap the annotations in `ClassVar`, and do not convert `setUpClass` to `setUp` just to satisfy the type checker:

```python
class TestUnsteadyProblem(unittest.TestCase):
    """This is a class with functions to test UnsteadyProblems."""

    basic_unsteady_problem: ps.problems.UnsteadyProblem

    @classmethod
    def setUpClass(cls) -> None:
        """Set up the shared test fixtures."""
        cls.basic_unsteady_problem = (
            problem_fixtures.make_basic_unsteady_problem_fixture()
        )
```

#### Deliberately Invalid Arguments

A rejection test passes a value whose type is intentionally wrong. Route the value through a local annotated as `Any`, and pass the local. The call site then type-checks without a pragma, the invalid value itself stays unchanged, and the annotation marks the invalidity as intentional right at the assignment:

```python
def test_wings_validation(self) -> None:
    """Test that non-list wings inputs are rejected."""
    bad_wings: Any = "not a list"
    with self.assertRaises(TypeError):
        ps.geometry.airplane.Airplane(wings=bad_wings)
```

Lists of invalid values follow the same recipe: `invalid_values: list[Any] = [0, -5, 2.5, "three"]`. Lists of valid values never use `Any`. Annotate them precisely, as in `valid_positions: list[np.ndarray | Sequence[float | int]]` for an array-like acceptance test.

#### Read-Only Property Tests

A test that proves a property is read-only assigns through `setattr` inside the `assertRaises` block. A direct assignment to a read-only property is a mypy error, and `setattr` exercises the same descriptor protocol at runtime, so the test still proves that the property rejects assignment:

```python
def test_chord_is_read_only(self) -> None:
    """Test that chord cannot be reassigned."""
    with self.assertRaises(AttributeError):
        setattr(self.basic_wing_cross_section, "chord", 2.0)
```

#### Narrowing After `assertIsNotNone`

mypy does not narrow a type on `self.assertIsNotNone(x)`. Keep the unittest assertion, and add a bare `assert x is not None` after it before the first use that needs the narrowed type:

```python
self.assertIsNotNone(wing.panels)
assert wing.panels is not None
self.assertEqual(wing.panels.shape, (4, 8))
```

#### Returning Floats from numpy Expressions

With `warn_return_any` enabled, a def annotated `-> float` cannot return a bare numpy expression, because the numpy stubs type many such expressions as `Any`. Wrap the return expression in `float(...)`, which is the same pattern the package uses in `_oscillation.py`:

```python
def custom_spacing(x: float) -> float:
    return float(np.sin(x))
```

---

## Docstring Format

### General Principles

1. **Use reStructuredText (rST) format**
2. **Short description starts on same line as triple quotes**
3. **Parameter descriptions are inline** (no separate type line)
4. **Begin descriptions with article + shape/type info for arrays** (e.g., "A (4,4) ndarray of floats...")
5. **Use present tense for descriptions** (e.g., "Returns..." not "Will return...")
6. **Avoid starting descriptions with "This..."**
7. **Follow the ASCII Only rule in [WRITING_STYLE.md](WRITING_STYLE.md)**, which covers all character substitutions (dashes, math symbols, smart quotes, ellipsis, arrows, emojis, and other typographic Unicode) used across the project's prose, comments, and docstrings.
8. **If the docstring is multiple paragraphs or contains any blank lines, place closing triple-quotes on their own line (otherwise, place them directly behind the last sentence)**
9. **Summary line is a single sentence.** Any additional description goes in a new paragraph after a blank line. docformatter enforces this: if the first paragraph contains multiple sentences, it moves all but the first into a new paragraph.
10. **No blank line between the closing triple-quotes and the next line of code.** docformatter enforces this too: a blank gap after the docstring will be removed.
11. **No backticks in prose.** Write identifiers, expressions, calls, and keyword assignments bare (the free_wake parameter, passive=True, get_logger("trim")). Single backticks are not code markup in rST (they render as italics), and the identifier casing already sets names apart from prose. The one place double backticks belong is inside an rST line block (a line starting with `|`) holding a standalone code example, as in the use-case blocks of `_transformations.py`. Comments follow the same rule, as the Markup and Quoting section of [WRITING_STYLE.md](WRITING_STYLE.md) records, along with the code span rules for Markdown files.
12. **Double-quote string values, never paths.** A str value (an accepted parameter value such as "sine", a dict key such as "position_E_Eo", a default such as "draw.webp", or an extension that is itself the str being passed or checked such as ".webp") is written in double quotes, matching black's quoting in code. Paths and glob patterns named as references (docs/AXES_POINTS_AND_FRAMES.md, a .psz file) are never quoted. Comments follow the same quoting, as described in [WRITING_STYLE.md](WRITING_STYLE.md), and so do runtime strings such as error messages, as described in [CODE_STYLE.md](CODE_STYLE.md).
13. **Bare text is still rST.** Since prose carries no literal markup, avoid sequences rST parses as markup: `|x|` (a substitution reference), `*x*` (emphasis), and `word_` followed by whitespace (a hyperlink reference). Reword instead, for example abs(angleY) rather than a barred magnitude.

### Module-Level Docstrings

Module-level docstrings appear at the very top of each Python file and describe the module's contents. For a package's `__init__.py` module, it should instead describe the package's contents.

```python
"""Contains the <placeholder> classes/functions/subpackages/directories/modules.

<Optional longer description block.>

<Optional citation block.>
"""
```

**Pattern:**

- Brief description using "Contains" (present tense)

### Function and Method Docstrings

```python
def function_name(
    param1: Type1,
    param2: Type2,
    param3: Type3,
) -> ReturnType:
    """Short description of what the function/method does.

    <Optional longer description block.>

    <Optional citation block.>

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

    <Optional longer description block.>

    <Optional citation block.>
    """
```

### Subclass Docstrings

When a class inherits from another class, use a modified pattern that avoids duplicating documentation from the parent class.

#### Subclass Class Docstring Template

```python
class ChildClass(ParentClass):
    """A subclass of ParentClass used to <description>.

    Inherits all parameters and methods from ParentClass without modification.

    <Optional longer description block.>

    <Optional citation block.>
    """
```

**Key points:**

- Short description explicitly mentions "A subclass of ParentClass"
- States what is inherited from the parent

#### Subclass `__init__` Docstring Template

```python
def __init__(
    self,
    inherited_param1: Type1,
    inherited_param2: Type2,
    new_param: Type3,
) -> None:
    """The initialization method.

    See ParentClass's initialization method for descriptions of inherited
    parameters.

    <Optional longer description block.>

    <Optional citation block.>

    :param new_param: Description of the new parameter unique to this subclass.
    :return: None
    """
    super().__init__(inherited_param1, inherited_param2)
    self.new_param = new_param
```

**Key points:**

- Reference the parent class's `__init__` docstring for inherited parameters
- Only document parameters that are NEW to the subclass
- Call `super().__init__()` with inherited parameters

### Public Subclasses of Private Parents

When a public class inherits from a private parent (a class in an underscore prefixed module like `_core.py`), the conventions above are inverted. The public child keeps a self-contained docstring that documents inherited methods and parameters as its own. The private parent's class docstring and `__init__` docstring use minimal descriptions that reference the public child. However, public methods and properties on the private parent must have self-contained docstrings because they appear on the public child's RTD page via inheritance (see "Private Parent Method and Property Docstrings" below).

This is because:

1. The ReadTheDocs site is purely public API. Private parents are excluded from the generated documentation, but inherited public methods and properties DO appear on the public child's page.
2. Users should never need to navigate to a private module to understand the public API.
3. Contributors reading the private parent's source code can easily navigate to the public child for full documentation of the class and `__init__`.

#### Private Parent Class Docstring Template

```python
class _CoreClass:
    """A core class used to contain the shared foundation of PublicClass and its
    feature variant siblings.

    See PublicClass for full documentation of the shared interface.

    <Optional longer description block.>

    <Optional citation block.>
    """
```

**Key points:**

- Reference the public child for full documentation of the shared interface
- Include a brief architectural description for contributors

#### Private Parent `__init__` Docstring Template

```python
def __init__(
    self,
    param1: Type1,
    param2: Type2,
) -> None:
    """The initialization method.

    See PublicClass's initialization method for full parameter descriptions.

    <Optional longer description block.>

    <Optional citation block.>

    :param param1: Brief description.
    :param param2: Brief description.
    :return: None
    """
```

**Key points:**

- Reference the public child's `__init__` docstring for full parameter descriptions
- Include brief parameter descriptions (enough for contributors to understand the code without navigating away)

#### Public Child Class Docstring Template

```python
class PublicClass(_core.CoreClass):
    """A class used to <description>.

    <Optional longer description block.>

    <Optional citation block.>
    """
```

**Key points:**

- Do not mention the private parent in the short description
- The class reads as a standalone public API entry point

#### Public Child `__init__` Docstring Template

```python
def __init__(
    self,
    inherited_param1: Type1,
    inherited_param2: Type2,
    new_param: Type3,
) -> None:
    """The initialization method.

    <Optional longer description block.>

    <Optional citation block.>

    :param inherited_param1: Full description.
    :param inherited_param2: Full description.
    :param new_param: Full description.
    :return: None
    """
    super().__init__(inherited_param1, inherited_param2)
    self.new_param = new_param
```

**Key points:**

- Document all parameters fully (inherited and new)
- Do not reference the private parent

#### Private Parent Method and Property Docstrings

Public methods and properties defined on a private parent are inherited by all public children and appear on their RTD pages. Their docstrings must therefore be self-contained and written for a public audience, unlike the class and `__init__` docstrings which can defer to the public child.

**Rules:**

1. **No deferral language.** Do not write "see child class for full details" or similar, since the docstring IS the documentation the user sees on the child's page.
2. **No references to specific sibling types.** A `CoreWingMovement` method docstring must not mention `WingCrossSectionMovement` or `AeroelasticWingCrossSectionMovement`, because the docstring appears on all siblings' RTD pages. Instead, reference the universal geometry class that the movement class manages (e.g., `WingCrossSection`), since geometry classes have no feature subclasses and are always correct.
3. **Use "each X's movement class" framing** when referring to child movement objects. For example, write "each WingCrossSection's movement class" rather than "its WingCrossSections' movement classes". This avoids implying that a movement class owns geometry objects (movement classes own other movement classes, and geometry classes own geometry classes).

**Example (correct):**

```python
# On CoreWingMovement
@property
def wing_cross_section_movements(self) -> tuple:
    """The movement classes for each of this Wing's WingCrossSections.

    :return: A tuple of movement classes, one per WingCrossSection.
    """
```

**Example (incorrect):**

```python
# References a specific sibling type
@property
def wing_cross_section_movements(self) -> tuple:
    """The WingCrossSectionMovements for this WingMovement.
    ...
    """

# Uses deferral language
def generate_wing_at_time_step(self, ...) -> Wing:
    """Generates a Wing at a single time step.

    See WingMovement for full details.
    ...
    """
```

**Scope:** This rule applies only to public methods and properties on core classes (those that will be inherited and displayed on RTD). Private helper methods (underscore prefixed) on core classes are internal and can use any convenient wording.

#### Multiple Public Siblings

When multiple public classes share the same private parent (e.g., `Movement`, `FreeFlightMovement`, and `AeroelasticMovement` all extending `CoreMovement`), each sibling maintains its own self-contained docstring. The inherited method descriptions can be tailored to each sibling's context (e.g., "Movement's sub movement objects" vs "FreeFlightMovement's sub movement objects").

### Private Names in Public Docstrings and Signatures

The API reference documents only the public modules. Anything defined in a private module has no page, and that includes classes whose own names carry no underscore, such as `Panel` in `_panel.py` and `CoupledUnsteadyRingVortexLatticeMethodSolver` in `_coupled_unsteady_ring_vortex_lattice_method.py`. A public docstring or signature that names one of them therefore renders as dead text: a class name the reader cannot look up, or a fully qualified path such as `pterasoftware._core.CoreUnsteadyProblem` in a signature. The rules below cover every way a private name can reach a rendered page. "Rendered" means the docstring of a public module, class, function, method, or property, including methods and properties inherited from a private parent (see "Private Parent Method and Property Docstrings").

#### Prose

1. **Never name a private class in rendered prose.** Point at the referent instead of naming its type: "The list of Wings associated with this movement", not "associated with this CoreWingMovement", and "The solver driving this problem, which provides the aerodynamic data from the current time step", not "The CoupledUnsteadyRingVortexLatticeMethodSolver instance providing aerodynamic data". Where the private class is a parent, describe what the public class adds rather than what it extends: "A class used to solve AeroelasticUnsteadyProblems with the unsteady ring vortex lattice method" and "**Key additions over the unsteady ring vortex lattice method:**", not "A subclass of CoupledUnsteadyRingVortexLatticeMethodSolver".
2. **Do not substitute a specific public sibling when several would work.** A statement must not become incorrect by omission. "The AirplaneMovement that owns this Wing's movement" is wrong when an `AeroelasticAirplaneMovement` also fits, so write "the Airplane movement class that owns this Wing's movement". When only one public class fits, name it.
3. **Never name a private hook or helper method.** Describe when the work happens instead of which override does it: "resets them at the start of each time step, and computes the moments about the strip leading edge points once those loads are known", not "overrides _reinitialize_step_arrays_hook to reset the SLEP arrays and overrides _process_panel_loads_hook to compute the moments".
4. **Do not defer to a private parent.** "See _CoupledUnsteadyProblem's initialization method for descriptions of inherited parameters" points the reader at a page that does not exist. Document the inherited parameters in the public child, as "Public Subclasses of Private Parents" requires.
5. **Contributor detail that needs private names goes in a comment.** The justification for why `Airplane.deep_copy_with_Cg_GP1_CgP1` copies what it copies names `_T_pas_G_Cg_to_GP1_CgP1` and `Panel.__deepcopy__`, so it lives in a comment at the top of the method body while the docstring keeps a one-sentence public summary. The comment is the right home for anything a contributor needs and a user does not.
6. **`Panel` is the standing exception.** Public docstrings name `Panel` throughout because it is the vocabulary of the mesh, and whether it becomes a public class or is reworded is an open decision. Leave existing `Panel` mentions as they are and do not add new private names on the strength of this exception.

#### Signatures

1. **Annotate with public types wherever the implementation allows.** A private type in a parameter annotation renders as an unlinked fully qualified path.
2. **When an annotation must be a private type, add an override.** Two shapes force this: a hook method whose override cannot narrow the parameter type, so the hook is annotated with the shared parent solver, and a base solver constructor that accepts the shared parent problem type so the derived solvers can pass their own problems through it. For those, add an entry to `_ANNOTATION_OVERRIDES` in `docs/website/conf.py`, keyed by the fully qualified class (for constructor parameters) or method, then by parameter name, giving the one public type that actually works. The build resolves each key against AutoAPI's object tree and fails on a stale key, so a rename cannot silently drop an override. The parameter's docstring stays exact without naming the private type: "The UnsteadyProblem to be solved. The derived solvers pass their own problem types through this parameter."
3. **Private bases are hidden automatically.** The class template drops any base whose path contains a private module or underscore prefixed name from the Bases line, so a public class extending a private parent shows no Bases line at all. Its docstring must stand on its own for that reason.
4. **No private parameters or sentinels in public signatures.** A private parameter such as a `_trust` token, or a private sentinel such as `_UNSET` as a default, renders with the signature. For a construction path that must skip the constructor's validation, allocate with `object.__new__(Cls)` and set the slots directly inside the class's own module, as `Airfoil.__deepcopy__` and `Airfoil.add_control_surface` do. A sentinel default is acceptable only on a deprecated parameter, which the reference filters out (next section).

#### Deprecated API

Deprecated functions, methods, properties, and parameters are filtered out of the API reference, so they need no docstring marker. The build detects them from the source, and the detection only works when the deprecation takes this exact shape:

- A function, method, or property is deprecated when a top-level statement of its body is `warnings.warn(..., DeprecationWarning)` (positional or `category=` keyword). Deprecated members are omitted from the reference entirely.
- A parameter is deprecated when such a call sits inside a top-level `if` statement whose test names the parameter, as in `if outline_A_lp is not _UNSET:`. Deprecated parameters are removed from the rendered signature along with their `:param:` field, while the rest of the signature and docstring render unchanged.

Do not move the `warnings.warn` call into a helper function, since the detection reads the deprecated member's own body. For a member, the call must be a top-level statement of that body, not nested in a block. For a parameter, the call may sit anywhere inside the guarding `if`, but the `if` itself must be top-level and its test must name the parameter. The docstring still documents the deprecated member or parameter for source readers and `help()`, in the form "A deprecated alias for outline_A_Lp. Passing it emits a DeprecationWarning, and it will be removed in v6.0.0."

### Property Docstring Template

```python
@property
def property_name(self) -> ReturnType:
    """Short description of what the property represents.

    <Optional longer description block.>

    <Optional citation block.>

    :return: Description of what is returned, including type, shape, units.
    """
```

For simple getter properties that just return a stored attribute, the docstring can be omitted if the attribute is already thoroughly documented in `__init__`'s docstring. This applies broadly to all simple getter properties, not only to cached properties with invalidating setters.

### Cached Properties with Invalidating Setters

When implementing a caching pattern where computed properties are lazily evaluated and cached, with setters that invalidate dependent caches, follow these conventions:

#### Property Getters for Cached Attributes

For simple getters that return a cached value (e.g., corner point positions that were previously plain attributes), use a brief docstring or omit the docstring entirely if the attribute is already thoroughly documented in `__init__`:

```python
@property
def Frpp_G_Cg(self) -> np.ndarray:
    # No docstring as this attribute is documented in __init__()'s docstring
    return self._Frpp_G_Cg


@property
def Frpp_GP1_CgP1(self) -> np.ndarray:
    """The position of the Panel's front right vertex (in the first Airplane's geometry
    axes, relative to the first Airplane's CG).

    :return: A (3,) ndarray of floats representing the position of the Panel's front
        right vertex (in the first Airplane's geometry axes, relative to the first
        Airplane's CG). The units are in meters. Returns None if not yet set or if
        Frpp_G_Cg has been modified since last set.
    """
    return self._Frpp_GP1_CgP1
```

#### Property Setters that Invalidate Caches

Do not add docstrings to setters. The cache invalidation behavior is an implementation detail.

```python
@Frpp_G_Cg.setter
def Frpp_G_Cg(self, newFrpp_G_Cg: np.ndarray) -> None:
    # No docstring as this is a setter method
    self._rightLeg_G = None
    self._frontLeg_G = None
    self._Frbvp_G_Cg = None
    self._Cpp_G_Cg = None
    self._unitNormal_G = None
    self._area = None
    self._aspect_ratio = None

    self.Frpp_GP1_CgP1 = None

    self._Frpp_G_Cg = newFrpp_G_Cg
```

#### Cached Computed Properties

For computed properties that are now cached (e.g., `rightLeg_G`, `area`), the existing docstring remains unchanged. Caching is an implementation detail that does not affect the public interface.

#### Class Docstring for Classes with Caching

When a class uses this caching pattern, add a section to the class docstring explaining the caching behavior:

```python
class Panel:
    """A class used to contain the panels of a Wing.

    Computed geometric properties (leg vectors, bound vortex points, collocation points,
    unit normals, area, and aspect ratio) are lazily evaluated and cached. Setting any
    corner point position invalidates all dependent cached values, ensuring consistency
    while avoiding redundant computation. Setting a corner point's local position
    (one of the parameters with a _G_Cg suffix), sets the corresponding global position
    (_GP1_CgP1 suffix) to None. It also sets this Panel's bound vortices and the loads
    on the Panel to None.
    """
```

### Optional Longer Description Blocks

Provides detailed explanations of the function/method's behavior. It can be one or more paragraphs. It can also be broken up with sections separated by sentence-case headers, wrapped with double-asterisks and padded with a blank line above and below. Avoid numbered or bulleted lists.

### Optional Citation Blocks

```python
"""
**Citation(s):**

Adapted from (can be more specific if the whole function/method wasn't adapted): <source>

Author(s): <author>

Date of retrieval (don't include if not known): <date>
"""
```

---

## Examples

### Example 1: Module-Level Docstrings

#### Public Package `__init__.py`

```python
"""Contains the geometry classes."""
```

#### Public Module

```python
"""Contains the Airfoil class."""
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
    chordwise coordinates.

    It returns a list of four column vectors containing the normalized components of
    the positions of points along the mean camber line (MCL) (in each Airfoil's axes,
    relative to each Airfoil's leading point).

    :param inner_airfoil: The wing section's inner Airfoil.
    :param outer_airfoil: The wing section's outer Airfoil.
    :param chordwise_coordinates: A (N,) ndarray of floats for the normalized
        chordwise coordinates where we'd like to sample each Airfoil's MCL. The values
        are normalized from 0.0 to 1.0 and are unitless.
    :return: A list of four (N,1) ndarrays of floats, where N is the number of points
        at which we'd like to sample each Airfoil's MCL. The ndarrays contain components
        of the positions of points along each Airfoil's MCL. In order, the ndarrays
        returned are, (1) the inner Airfoil's MCL points' y components, (2) the inner
        Airfoil's MCL points' x components (3) the outer Airfoil's MCL points'
        y components, and (4) the outer Airfoil's MCL points' x components. The values
        are normalized from 0.0 to 1.0 and are unitless.
    """
```

### Example 3: Function with Transformation Matrices

```python
def _get_mcs_points(
    T_pas_Wcsi_Lpi_Wn_Ler: np.ndarray,
    T_pas_Wcso_Lpo_Wn_Ler: np.ndarray,
    inner_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    outer_wing_cross_section: wing_cross_section_mod.WingCrossSection,
    inner_mcl_pointsY_Ai_LpAi: np.ndarray,
    inner_mcl_pointsX_Ai_LpAi: np.ndarray,
    outer_mcl_pointsY_Ao_LpAo: np.ndarray,
    outer_mcl_pointsX_Ao_LpAo: np.ndarray,
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
    :param inner_mcl_pointsY_Ai_LpAi: A (M,1) ndarray of floats, where M is the
        number of chordwise points in the mesh. Each element represents the y component
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

### Example 4: Public Method with Array-Like Parameters

```python
def __init__(
    self,
    name: str = "NACA0012",
    outline_A_Lp: np.ndarray | Sequence[Sequence[float | int]] | None = None,
    resample: bool = True,
    n_points_per_side: int = 400,
) -> None:
    """The initialization method.

    :param name: The name of the Airfoil. It should correspond to the name of a file in
        the airfoils directory, or to a valid NACA 4-series airfoil (once converted to
        lower-case and stripped of leading and trailing whitespace) unless you are
        passing in your own array of points using outline_A_Lp. Note that NACA0000 isn't
        a valid NACA-series airfoil. The default is "NACA0012".
    :param outline_A_Lp: An array-like object of numbers (int or float) with shape
        (N,2) representing the 2D points making up the Airfoil's outline (in airfoil
        axes, relative to the leading point). If you wish to load coordinates from the
        airfoils directory, leave this as None, which is the default. Can be a tuple,
        list, or ndarray. Values are converted to floats internally. Make sure all
        x component values are in the range [0.0, 1.0]. The default value is None.
    :param resample: Determines whether to resample the points defining the Airfoil's
        outline. This applies to points passed in by the user or to those from the
        airfoils directory. I highly recommend setting this to True. Can be a bool or
        a numpy bool and will be converted internally to a bool. The default is True.
    :param n_points_per_side: The number of points to use when creating the Airfoil's
        MCL and when resampling the upper and lower parts of the Airfoil's outline. It
        must be a positive int greater than or equal to 3. The resampled outline will
        have a total number of points equal to (2 * n_points_per_side) - 1. I highly
        recommend setting this to at least 100. The default value is 400.
    :return: None
    """
```

### Example 5: Method Returning Self-Reference

```python
def add_control_surface(
    self, deflection: float | int, hinge_point: float | int
) -> Airfoil:
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

### Example 6: Method with Optional Return

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

### Example 7: Method Returning Array

```python
def get_resampled_mcl(
    self, mcl_fractions: np.ndarray | Sequence[float]
) -> np.ndarray:
    """Returns a ndarray of points along the mean camber line (MCL), resampled from the
    mcl_A_Lp attribute.

    It is used to discretize the MCL for meshing.

    :param mcl_fractions: A (N,) array-like object of floats representing normalized
        distances along the MCL (from the leading to the trailing edge) at which to
        return the resampled MCL points. Can be a tuple, list, or ndarray. The first
        value must be 0.0, the last must be 1.0, and the remaining must be in the range
        [0.0, 1.0]. All values must be non duplicated and in ascending order.
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
param: bool | np.bool  # Accepts both Python and numpy bools
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
