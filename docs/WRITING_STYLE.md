# Writing Style Guidelines
Guidelines for developers and Claude when writing comments, docstrings, and documentation for Ptera Software.

## Terminology
- **"Ptera Software"**: When writing as text, always write as two words without a hyphen, each being capitalized (never "ptera", "ptera software", or "PteraSoftware")
- **"cross section"**: Always write as two words, never hyphenated (not "cross-section")  
- **Object references**: When referring to code objects, use proper class naming convention. The capitalization indicates that we are talking about a code object, not an abstraction. You don't need to add "object" or "objects" after the class name since the capitalization already makes this clear (e.g. "update the Wings" instead of "update the Wing objects"). In summary, when talking about code objects:
  - GOOD: "the previous WingCrossSection"
  - BAD: "the previous cross section"
  - GOOD: "this Wing"
  - BAD: "this wing"
  - GOOD: "update the Wings"
  - BAD: "update the Wing objects" (unnecessary)
- **Abstract references**: When referring to abstractions, use lowercase and separate individual words with a space (e.g. "an airplane's wings are used to generate lift" and "the cross section of a wing typically has a streamlined shape known as an airfoil"). This is to distinguish them from code objects.
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the `AXES_AND_COORDINATE_SYSTEMS.md` and `AXES_POINTS_AND_FRAMES.md` documents when writing about or referencing in text vector-valued variables or things such as transformation and rotation matrices.

## Running a CodeSpell Spell Check
For Claude, use the following command. For developers, see the `CONTRIBUTING.md` file for instructions.

```bash
cd ${WORKSPACE} && ".venv/Scripts/codespell.exe" --ignore-words=.codespell-ignore.txt --skip="*/_build/*,*.dat"
```

## Formatting Docstrings with docformatter
For Claude, use the following command. For developers, see the `CONTRIBUTING.md` file for instructions.

```bash
cd ${WORKSPACE} && ".venv/Scripts/python.exe" docformatter --black --in-place pterasoftware -r
```

## Miscellaneous Guidelines
- Avoid abbreviations in text unless they are well-known in the context.
- Never hyphenate words in docstrings or comments. This is because they are often incorrectly wrapped by docformatter and incorrectly rendered in PyCharm's quick documentation. For example, even though not standard grammar, it's okay to write "non symmetric" instead of "non-symmetric".
- In documentation, docstrings, and comments, represent subtraction using hyphens surrounded by spaces ( - ); never use em-dashes (—) or en-dashes (–).
- In documentation, docstrings, and comments, never use a multiplication sign (×); always use a lowercase x or an asterisk surrounded by spaces (" x " or " * ").
- In documentation, docstrings, and comments, never use the pi symbol (π); always write "pi" instead (e.g., "2 * pi" not "2 π"). The same goes for other Greek letters (e.g., use "alpha" instead of "α").
- In documentation, docstrings, and comments, never use the approximately-equal sign (≈); always write " ~ " instead (e.g., "a ~ b" not "a≈b").
- When referring to axes, coordinates, or planes, use  lowercase letters without hyphens between coordinate letters and descriptors (e.g., "x axis", "y component", "xz plane", "z direction"). Never use uppercase letters for axis references in text.
- Never use emojis in code, comments, docstrings, or documentation.
- Always use straight single and double quotes, not curly ones.
- In markdown files, use "…" instead of "...". In all other files, use "...".
- Always end *.py file with an empty line.
- Preserve existing comment structure and detail level.
- Write comments as complete sentences ending with a period.
- Prefer comments on their own line above the code they describe.
- Use American English spelling (e.g., "color" not "colour", "center" not "centre").
- Write "time step" instead of "timestep", "time-step", or "step".