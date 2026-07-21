"""Verify that the environment holds the versions pinned in requirements_dev.txt.

The octowrap and mypy hooks run with language: system, so pre-commit executes whatever
tool versions the active virtual environment holds instead of managing pinned copies
itself. This pre-commit hook closes that gap: it reads every exact pin (==) from
requirements_dev.txt, compares each against the version installed in the running
environment, and fails with a reinstall hint when any package is missing or has drifted
from its pin.
"""

import re
import sys
from importlib import metadata
from pathlib import Path

PIN_PATTERN = re.compile(r"^\s*([A-Za-z0-9][A-Za-z0-9._-]*)\s*==\s*([^\s;#]+)")


def read_pins(path: Path) -> dict[str, str]:
    """Returns the exact version pins from a requirements file.

    :param path: The path to the requirements file to read.
    :return: A dict mapping each exactly pinned (==) package name to its pinned version.
        Packages with range constraints are not included.
    """
    pins: dict[str, str] = {}
    for line in path.read_text().splitlines():
        match = PIN_PATTERN.match(line)
        if match:
            pins[match.group(1)] = match.group(2)
    return pins


def main() -> int:
    """Checks each pinned package's installed version and prints any mismatches.

    :return: An int representing the exit code, which is 0 if every pinned package is
        installed at its pinned version and 1 otherwise.
    """
    requirements_path = Path("requirements_dev.txt")
    problems: list[str] = []
    for name, pinned in read_pins(requirements_path).items():
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            problems.append(
                f"{requirements_path} pins {name} == {pinned}, but it is not"
                " installed in the active environment."
            )
            continue
        if installed != pinned:
            problems.append(
                f"{requirements_path} pins {name} == {pinned}, but the active"
                f" environment has {installed}."
            )
    if problems:
        for problem in problems:
            print(problem)
        print()
        print("Re-run: pip install -r requirements.txt -r requirements_dev.txt")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
