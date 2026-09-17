# Writing Style

Guidelines when writing comments, docstrings, and documentation for Ptera Software.

## Terminology

- **"Ptera Software"**: When referring to the project, package, or codebase by its proper name in prose, always write "Ptera Software" as two capitalized words without a hyphen. Never use "ptera", "ptera software", "PteraSoftware", or "Ptera" alone in prose. This includes possessives ("Ptera Software's", not "Ptera's") and shortenings ("the Ptera Software solver", not "the Ptera solver"). Identifier strings keep their canonical form and are not affected: the GitHub repo name `PteraSoftware`, the Python package `pterasoftware`, and paths like `~/Documents/GitHub/PteraSoftware/` stay as-is.
- **Object references**: When referring to code objects, use proper class naming convention. The capitalization indicates that we are talking about a code object, not an abstraction. You don't need to add "object" or "objects" after the class name since the capitalization already makes this clear (e.g. "update the Wings" instead of "update the Wing objects"). The examples below are written as they appear in comments and docstrings. Markdown files, issue bodies, and pull request descriptions also put the class name in a code span, as the [Markup and Quoting](#markup-and-quoting) section describes. In summary, when talking about code objects:
    - GOOD: "the previous WingCrossSection"
    - BAD: "the previous cross section"
    - GOOD: "this Wing"
    - BAD: "this wing"
    - GOOD: "update the Wings"
    - BAD: "update the Wing objects" (unnecessary)
- **`WingCrossSection` and "wing cross section"**: Do not use "WCS" (or any other abbreviation). Also, never hyphenate "cross section". Exception: a capitalization variant of "WCS" is allowed in a variable name when that name is written in exactly one of the forms required by [AXES_POINTS_AND_FRAMES.md](AXES_POINTS_AND_FRAMES.md) or [ANGLE_VECTORS_AND_TRANSFORMATIONS.md](ANGLE_VECTORS_AND_TRANSFORMATIONS.md) (e.g., the axes IDs in `angles_Wcsp_to_Wcs_ixyz` and the point IDs in `Lp_Wcsp_Lpp`). Such names keep their canonical form when referenced in text.
- **Abstract references**: When referring to abstractions, use lowercase and separate individual words with a space (e.g. "an airplane's wings are used to generate lift" and "the cross section of a wing typically has a streamlined shape known as an airfoil"). This is to distinguish them from code objects.
- **"strip leading edge point"**: When spelled out in prose, do not capitalize "strip leading edge point" unless it starts a sentence.
- **"V-tail"**: Always write "V-tail" (capital V, hyphenated) in prose. Never use "v-tail", "v tail", "Vtail", or other variants. Identifier strings keep their canonical form and are not affected: variable names like `v_tail_movement` and string literals like `"V-Tail"` stay as-is.
- **"time step"**: Write a simulation time increment as either "time step" or "step". The two terms are interchangeable. When using the two-word form, always write it as "time step", never "timestep" or "time-step".
- **Axis references**: When referring to axes, coordinates, or planes, use lowercase letters without hyphens between coordinate letters and descriptors (e.g., "x axis", "y component", "xz plane", "z direction"). Never use uppercase letters for axis references in text.
- **Abbreviations**: Avoid abbreviations in text unless they are well-known in the context.
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the [AXES_POINTS_AND_FRAMES.md](AXES_POINTS_AND_FRAMES.md) and [ANGLE_VECTORS_AND_TRANSFORMATIONS.md](ANGLE_VECTORS_AND_TRANSFORMATIONS.md) documents when writing about or referencing in text vector-valued variables or things such as transformation and rotation matrices.

## Sentence Structure

- Write comments as complete sentences starting with a capital letter and ending with a period. This means actual complete sentences with a subject and verb, not fragments with punctuation stapled on.
- Preserve existing comment structure and detail level.
- Prefer comments on their own line above the code they describe.
- Avoid semicolons in prose. Replace them with commas, colons, or sentence breaks.
- Capitalize the first word after a colon only if what follows is a complete sentence. Keep it lowercase for fragments, lists, or continuations. When the capitalized word reads as awkward or looks like a typo, split the two clauses into separate sentences instead. The colon is not required, and two plain sentences are usually clearer than a colon followed by a capital.
- Use American English spelling (e.g., "color" not "colour", "center" not "centre").
- Do not use asterisks for bold or italics, and do not capitalize or use all caps for emphasis in comments.
- Never use doubled hyphens (`--`) as an em dash substitute. Restructure the sentence instead: split it into two sentences, use a comma, use a colon, or place the clause in parentheses. See [How to Handle Each Forbidden Case](#how-to-handle-each-forbidden-case) for the full em dash and en dash rules.
- Use "and" instead of "&".

## Math and Numbers

- For subtraction, use a hyphen surrounded by spaces (e.g., "a - b").
- For multiplication, use a lowercase x or an asterisk, both surrounded by spaces (e.g., "8 x 8 panels" or "2 * pi"). Compound units are exempt: write them without spaces (e.g., "N\*m", "kg\*m^2", and "N\*m\*s/rad").
- For division, use a slash surrounded by spaces (e.g., "a / b"). Derivative notation is exempt: write it without spaces (e.g., "d(theta)/dt", "d^2(theta)/dt^2", and "dGamma/dt").
- For equals signs in inline equations and value descriptions, surround with spaces (e.g., "area = 0.5", not "area=0.5").
- For approximately-equal, use a tilde surrounded by spaces (e.g., "a ~ b").
- When writing tuples or vectors in prose, include a space after each comma (e.g., "(1.0, 2.0, 3.0)", not "(1.0,2.0,3.0)").
- When writing numerals that represent float values in prose, include a trailing ".0" (e.g., "1.0" not "1", "0.0" not "0") when the context describes float-valued quantities.
- Include a space between a numeral and its unit (e.g., "0.5 m" not "0.5m", "1.0 s" not "1.0s").

## Markup and Quoting

The same name is marked up differently depending on where it appears. Comments and docstrings write it bare, while Markdown files, issue bodies, and pull request descriptions put it in a code span. The rules below say which applies where.

### Comments and Docstrings

- No backticks. Write identifiers, expressions, calls, and keyword assignments bare (the free_wake parameter, passive=True, get_logger("trim")). Docstrings are rST, where single backticks render as italics rather than code, and comments follow the same rule so that the two read alike. The capitalization of a class name already sets it apart from prose, as the Terminology section describes.
- Write a str value in double quotes, matching black's quoting in code: an accepted parameter value such as "sine", a dict key such as "position_E_Eo", a default such as "draw.webp", or an extension that is itself the str being passed or checked such as ".webp". Never quote a path or glob pattern named as a reference (docs/AXES_POINTS_AND_FRAMES.md, a .psz file).
- Docstrings carry two further rules, for the one place double backticks belong and for sequences that rST parses as markup. See [TYPE_HINT_AND_DOCSTRING_STYLE.md](TYPE_HINT_AND_DOCSTRING_STYLE.md). Runtime strings such as error messages follow the same quoting, as described in [CODE_STYLE.md](CODE_STYLE.md).

### Markdown Files, Issue Bodies, and Pull Request Descriptions

- Put every code token in a code span (single backticks): shell commands (`pre-commit run --all-files codespell`), CLI flags (`--all-files`), keyword assignments (`dtype=float`), slash paths (`examples/`, `tests/unit/`), class names (`WingCrossSection`), variable names of every form (`wing_cross_section`, `omegasRad_BP1__E`, `_UNSET`), file names (`AXES_POINTS_AND_FRAMES.md`, `_meshing.py`), an extension used as a file-type reference (a `.psz` file), and any other code, including in headers.
- A str value is code too, so its double quotes go inside the code span: `"sine"`, `"edge_defined"`, and `".webp"` when the extension is itself the str being passed or checked.
- Attach plurals and possessives outside the code span: `WingMovement`s, `OperatingPoint`'s.
- Numerals, array shapes, and inline math stay bare and follow the Math and Numbers section: 1.0e-10, (M, 4), and r1 * r2.
- When pointing the reader at another document in this repository, use a Markdown link ([Axes, Points, and Frames](AXES_POINTS_AND_FRAMES.md)). When naming a file without sending the reader to it, use the backticked file name.
- Never escape underscores with backslashes. Any text containing an underscore belongs in a code span, where Markdown parses nothing, and that includes naming templates such as `[variable name]_[axes ID]`.
- Issue and pull request titles contain no backticks. Write identifiers bare there.

## File Formatting

- Always end `*.py` files with an empty line.
- In markdown files, always include a blank line after a header line. Also precede them with a blank line, except for header lines that happen to also be the first line in their file.
- In markdown files, do not use trailing whitespace for line breaks. Markdown already handles breaks between paragraphs, list items, and headings.
- In markdown files, do not use hard line wrapping.

## ASCII Only

In all written contributions to this project (code, comments, docstrings, documentation, commit messages, file contents), use only the 95 printable ASCII characters (0x20 through 0x7E: space, letters, digits, and standard punctuation).

The constraint is on what you author. Content quoted verbatim from external sources (third-party error messages, file contents you are citing, output captured from tools, web references) does not need to be transliterated.

### Forbidden Characters

- Em dashes and en dashes.
- Non-ASCII arrows of any kind (rightwards-arrow, leftwards-arrow, up-arrow, down-arrow, double-arrows, etc.).
- Emojis.
- Mathematical or scientific symbol characters (Greek pi, infinity sign, plus-or-minus sign, multiplication sign, division sign, micro sign, degree sign, superscripts, subscripts, set-theory symbols, etc.).
- Smart or curly quotes and apostrophes. Use straight ASCII `"` and `'`.
- Ellipsis character. Use three periods `...`.
- Bullets, middle dots, non-breaking spaces, or any other typographic Unicode.

### How to Handle Each Forbidden Case

- **Em dashes are special.** Do not substitute with doubled hyphens (`--`). Instead, restructure the sentence: split into two sentences, use a comma, use a colon, or place the clause in parentheses. The doubled-hyphen substitute is banned because it tends to be overused as a verbal tic.
- **En dashes used in ranges** (for example, between two numbers or two dates) may be replaced with a plain hyphen (`-`). This is technically incorrect typography, but the en dash is forbidden, so the hyphen is the accepted substitute in range contexts only.
- **En dashes used elsewhere** (for example, as a substitute for em dash in some style guides, or to join compound modifiers) must be restructured the same way as em dashes. Do not paper over them with a plain hyphen.
- **Arrows**: write `->`, `<-`, `=>`, `<=`, or use words.
- **Math symbols**: spell them out (`pi`, `infinity` or `inf`, `+/-`, `*`, `/`, `micro`, `deg`, `^2`, etc.).
- **Quotes**: straight ASCII `"` and `'`.
- **Ellipsis**: three periods `...`.

### When in Doubt

If a legitimate reason to violate this rule arises (a contribution must include a file with Unicode contents, or correctness requires a specific Unicode character), flag it in the PR description rather than assume.

## Running CodeSpell

CodeSpell is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files codespell
```

## Running docformatter

docformatter is configured as a pre-commit hook. Run it with:

```shell
pre-commit run --all-files docformatter
```

## Running `ascii-only`

`ascii-only` is configured as a pre-commit hook that enforces the [ASCII Only](#ascii-only) rule above. Run it with:

```shell
pre-commit run --all-files ascii-only
```

The hook reports each violation with its line, column, the offending character, its Unicode code point and name, and its UTF-8 byte sequence. The check is read-only. It never modifies files. Qt Designer `.ui` files are excluded because Qt regenerates them on every save.
