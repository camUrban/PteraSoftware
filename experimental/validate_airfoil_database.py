"""Script to validate all airfoils in the database.

This script attempts to load each airfoil from the database and reports any
that fail validation. This helps identify bad airfoil files that should be
fixed or removed from the database.
"""

import importlib.resources
import sys
from pathlib import Path

# Add parent directory to path so we can import pterasoftware
sys.path.insert(0, str(Path(__file__).parent.parent))

import pterasoftware as ps


def get_all_airfoil_names() -> list[str]:
    """Get all airfoil names from the database."""
    airfoils_dir = importlib.resources.files("pterasoftware.geometry").joinpath(
        "_airfoils"
    )
    names = []
    for item in airfoils_dir.iterdir():
        if item.name.endswith(".dat"):
            names.append(item.name[:-4])  # Remove .dat extension
    return sorted(names)


def validate_all_airfoils(verbose: bool = False) -> None:
    """Try to load each airfoil and report failures."""
    airfoil_names = get_all_airfoil_names()
    print(f"Found {len(airfoil_names)} airfoils in database.\n")

    failed: list[tuple[str, str]] = []
    passed: list[str] = []

    for name in airfoil_names:
        try:
            _ = ps.geometry.airfoil.Airfoil(name=name)
            passed.append(name)
        except Exception as e:
            failed.append((name, str(e)))

    print(f"PASSED: {len(passed)} airfoils loaded successfully.")
    print(f"FAILED: {len(failed)} airfoils failed validation.\n")

    if failed:
        # Group failures by error message
        errors_by_type: dict[str, list[str]] = {}
        for name, error in failed:
            if error not in errors_by_type:
                errors_by_type[error] = []
            errors_by_type[error].append(name)

        print("FAILURE SUMMARY BY ERROR TYPE:")
        print("=" * 80)
        for error, names in sorted(errors_by_type.items(), key=lambda x: -len(x[1])):
            print(f"\n{error}")
            print(f"  Count: {len(names)}")
            if verbose:
                print(f"  Airfoils: {', '.join(names)}")
            else:
                preview = names[:5]
                if len(names) > 5:
                    print(f"  Examples: {', '.join(preview)}, ...")
                else:
                    print(f"  Airfoils: {', '.join(preview)}")
    else:
        print("All airfoils passed validation!")


if __name__ == "__main__":
    validate_all_airfoils()
