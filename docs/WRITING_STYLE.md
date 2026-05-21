# Writing Style

Guidelines when writing comments, docstrings, and documentation for Ptera Software.

## Terminology

- **"Ptera Software"**: When referring to the project, package, or codebase by its proper name in prose, always write "Ptera Software" as two capitalized words without a hyphen. Never use "ptera", "ptera software", "PteraSoftware", or "Ptera" alone in prose. This includes possessives ("Ptera Software's", not "Ptera's") and shortenings ("the Ptera Software solver", not "the Ptera solver"). Identifier strings keep their canonical form and are not affected: the GitHub repo name "PteraSoftware", the Python package "pterasoftware", and paths like "~/Documents/GitHub/PteraSoftware/" stay as-is.
- **Object references**: When referring to code objects, use proper class naming convention. The capitalization indicates that we are talking about a code object, not an abstraction. You don't need to add "object" or "objects" after the class name since the capitalization already makes this clear (e.g. "update the Wings" instead of "update the Wing objects"). In summary, when talking about code objects:
    - GOOD: "the previous WingCrossSection"
    - BAD: "the previous cross section"
    - GOOD: "this Wing"
    - BAD: "this wing"
    - GOOD: "update the Wings"
    - BAD: "update the Wing objects" (unnecessary)
- **"WingCrossSection" and "wing cross section"**: Do not use "WCS" (or any other abbreviation). Also, never hyphenate "cross section".
- **Abstract references**: When referring to abstractions, use lowercase and separate individual words with a space (e.g. "an airplane's wings are used to generate lift" and "the cross section of a wing typically has a streamlined shape known as an airfoil"). This is to distinguish them from code objects.
- **CRITICAL**: Follow the formalized coordinate system naming conventions exactly as described in the `AXES_AND_COORDINATE_SYSTEMS.md` and `AXES_POINTS_AND_FRAMES.md` documents when writing about or referencing in text vector-valued variables or things such as transformation and rotation matrices.

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

## Running ascii-only

ascii-only is configured as a pre-commit hook that enforces the [ASCII Only](#ascii-only) rule below. Run it with:

```shell
pre-commit run --all-files ascii-only
```

The hook reports each violation with its line, column, the offending character, its Unicode code point and name, and its UTF-8 byte sequence. The check is read-only; it never modifies files. Qt Designer `.ui` files are excluded because Qt regenerates them on every save.

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

## Miscellaneous Guidelines

- Avoid abbreviations in text unless they are well-known in the context.
- For subtraction, use a hyphen surrounded by spaces (e.g., "a - b").
- For multiplication, use a lowercase x or an asterisk, both surrounded by spaces (e.g., "8 x 8 panels" or "2 * pi").
- For approximately-equal, use a tilde surrounded by spaces (e.g., "a ~ b").
- When referring to axes, coordinates, or planes, use lowercase letters without hyphens between coordinate letters and descriptors (e.g., "x axis", "y component", "xz plane", "z direction"). Never use uppercase letters for axis references in text.
- Always end *.py file with an empty line.
- Preserve existing comment structure and detail level.
- Write comments as complete sentences ending with a period.
- Prefer comments on their own line above the code they describe.
- Use American English spelling (e.g., "color" not "colour", "center" not "centre").
- Write "time step" instead of "timestep", "time-step", or "step".
- In markdown files, always include a blank line after a header line. Also precede them with a blank line, except for header lines that happen to also be the first line in their file.
- In markdown files, do not use trailing whitespace for line breaks. Markdown already handles breaks between paragraphs, list items, and headings.
- In markdown files, do not use hard line wrapping.
