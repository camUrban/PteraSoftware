"""Script to validate all valid NACA 4 series airfoils.

This script generates all valid NACA 4 series airfoil designations based on
the constraints defined in the Airfoil class, attempts to load each one, and
reports any failures.

NACA 4 series constraints:
1. Must be exactly 4 digits (0000-9999)
2. Cannot be "0000"
3. Thickness (last two digits) must be <= 30%
4. First two digits must either both be zero (symmetric) or both be non-zero
   (cambered)
5. For cambered airfoils: camber_loc >= max_camber + thickness/2
   (second digit * 10% >= first digit * 1% + last two digits * 0.5%)
"""

import sys
from pathlib import Path

# Add parent directory to path so we can import pterasoftware
sys.path.insert(0, str(Path(__file__).parent.parent))

import pterasoftware as ps


def generate_valid_naca4_names() -> list[str]:
    """Generate all valid NACA 4 series airfoil names."""
    valid_names = []

    for first_digit in range(10):  # 0-9: max camber (%)
        for second_digit in range(10):  # 0-9: camber location (x10%)
            for thickness in range(1, 31):  # 01-30: thickness (%)

                # Constraint: Cannot be 0000
                if first_digit == 0 and second_digit == 0 and thickness == 0:
                    continue

                # Constraint: First two digits must both be zero or both be non-zero
                max_camber = first_digit * 0.01
                camber_loc = second_digit * 0.1
                if (max_camber > 0) != (camber_loc > 0):
                    continue

                # Constraint: For cambered airfoils, camber_loc >= max_camber + thickness/2
                if max_camber > 0:
                    thickness_fraction = thickness * 0.01
                    if camber_loc < max_camber + thickness_fraction / 2:
                        continue

                # Build the NACA name
                name = f"NACA{first_digit}{second_digit}{thickness:02d}"
                valid_names.append(name)

    return valid_names


def validate_all_naca4_airfoils(verbose: bool = False) -> None:
    """Try to load each valid NACA 4 series airfoil and report failures."""
    naca_names = generate_valid_naca4_names()
    print(f"Generated {len(naca_names)} valid NACA 4 series designations.\n")

    failed: list[tuple[str, str]] = []
    passed: list[str] = []

    for name in naca_names:
        try:
            _ = ps.geometry.airfoil.Airfoil(name=name)
            passed.append(name)
        except Exception as e:
            failed.append((name, str(e)))

    print(f"PASSED: {len(passed)} airfoils loaded successfully.")
    print(f"FAILED: {len(failed)} airfoils failed.\n")

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
        print("All NACA 4 series airfoils passed validation!")


if __name__ == "__main__":
    # Print some examples of valid names first
    names = generate_valid_naca4_names()
    print("Examples of valid NACA 4 series names:")
    print(f"  Symmetric: {[n for n in names if n[4:6] == '00'][:5]}")
    print(f"  Cambered:  {[n for n in names if n[4:6] != '00'][:5]}")
    print()

    validate_all_naca4_airfoils()
