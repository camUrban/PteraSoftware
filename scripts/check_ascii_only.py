"""Verify that the given files contain only ASCII-safe characters.

This pre-commit hook flags any character that is not tab (U+0009), line
feed (U+000A), carriage return (U+000D), or printable ASCII (U+0020
through U+007E). It exists to enforce the project-wide convention that
committed text files stay within the 95 printable ASCII characters plus
the usual whitespace.

Each violation is reported with its line, column, the offending
character, its Unicode code point, its Unicode name (when available),
and its UTF-8 byte sequence.
"""

import sys
import unicodedata
from pathlib import Path

ALLOWED: frozenset[int] = frozenset({0x09, 0x0A, 0x0D} | set(range(0x20, 0x7F)))


def describe(char: str) -> str:
    """Format a one-line description of a disallowed character."""
    codepoint = ord(char)
    try:
        name = unicodedata.name(char)
    except ValueError:
        name = "no Unicode name"
    utf8 = " ".join(f"{b:02X}" for b in char.encode("utf-8"))
    return f"U+{codepoint:04X} {char!r} ({name}; UTF-8: {utf8})"


def find_violations(path: Path) -> list[tuple[int, int, str]]:
    """Return (line, col, description) tuples for each disallowed character."""
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        return [(0, 0, f"could not decode as UTF-8 ({exc})")]
    violations: list[tuple[int, int, str]] = []
    line = 1
    col = 0
    for char in text:
        if char == "\n":
            line += 1
            col = 0
            continue
        col += 1
        if ord(char) not in ALLOWED:
            violations.append((line, col, describe(char)))
    return violations


def main(argv: list[str]) -> int:
    exit_code = 0
    for arg in argv:
        path = Path(arg)
        try:
            violations = find_violations(path)
        except OSError as exc:
            print(f"{path}: could not read ({exc})", file=sys.stderr)
            exit_code = 1
            continue
        if violations:
            exit_code = 1
            for line, col, description in violations:
                print(f"{path}:{line}:{col}: {description}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
