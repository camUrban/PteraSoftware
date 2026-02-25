"""Helper functions for convergence analysis of GammaBot simulations.

This module provides geometry resampler functions that can be used with the
analyze_unsteady_convergence_non_trapezoidal function in pterasoftware.convergence.
"""

from pathlib import Path
from typing import Callable

import dxf_to_csv
import numpy as np


def create_gammabot_geometry_resampler(
    dxf_filepath: Path | str,
) -> Callable[[int, int], np.ndarray]:
    """Create a geometry resampler function for GammaBot wings.

    GammaBot has two identical wings (left and right), so the wing_id parameter
    is ignored and the same geometry is returned for both wings.

    :param dxf_filepath: Path to the DXF file with wing outline.
    :return: Resampler function that takes (wing_id, num_sections) and returns
        an (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.

    Example usage::

        from pathlib import Path
        import pterasoftware as ps
        from convergence_helpers import create_gammabot_geometry_resampler

        # Create the resampler
        dxf_path = Path(__file__).parent / "gammabot_approximate_wing.dxf"
        resampler = create_gammabot_geometry_resampler(dxf_path)

        # Use with convergence analysis
        result = ps.convergence.analyze_unsteady_convergence_non_trapezoidal(
            ref_problem=my_problem,
            wing_geometry_resampler=resampler,
            ...
        )
    """
    dxf_filepath_str = str(dxf_filepath)

    def resampler(wing_id: int, num_sections: int) -> np.ndarray:
        """Resample wing geometry at the specified number of spanwise sections.

        :param wing_id: The wing ID (ignored for GammaBot as both wings are identical).
        :param num_sections: Number of spanwise sections to generate.
        :return: (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
        """
        # GammaBot has identical left/right wings, so wing_id is ignored
        return dxf_to_csv.process_dxf_to_wing_section_data(
            dxf_filepath_str, num_sections
        )

    return resampler


def create_rectangular_wing_resampler(
    span: float,
    chord: float,
) -> Callable[[int, int], np.ndarray]:
    """Create a geometry resampler for a simple rectangular wing.

    This is useful for testing and validation purposes.

    :param span: Total span of the wing.
    :param chord: Constant chord length.
    :return: Resampler function that takes (wing_id, num_sections) and returns
        an (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
    """

    def resampler(wing_id: int, num_sections: int) -> np.ndarray:
        """Resample rectangular wing geometry at the specified sections.

        :param wing_id: The wing ID (ignored for rectangular wing).
        :param num_sections: Number of spanwise sections to generate.
        :return: (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
        """
        dy = span / num_sections
        result = np.zeros((num_sections + 1, 4))
        result[0, :] = [0.0, 0.0, 0.0, chord]  # root
        for i in range(1, num_sections + 1):
            result[i, :] = [0.0, dy, 0.0, chord]
        return result

    return resampler


def create_tapered_wing_resampler(
    span: float,
    root_chord: float,
    tip_chord: float,
) -> Callable[[int, int], np.ndarray]:
    """Create a geometry resampler for a tapered wing.

    This is useful for testing and validation purposes.

    :param span: Total span of the wing.
    :param root_chord: Chord length at the root.
    :param tip_chord: Chord length at the tip.
    :return: Resampler function that takes (wing_id, num_sections) and returns
        an (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
    """

    def resampler(wing_id: int, num_sections: int) -> np.ndarray:
        """Resample tapered wing geometry at the specified sections.

        :param wing_id: The wing ID (ignored for tapered wing).
        :param num_sections: Number of spanwise sections to generate.
        :return: (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
        """
        dy = span / num_sections
        result = np.zeros((num_sections + 1, 4))

        for i in range(num_sections + 1):
            # Linear interpolation of chord from root to tip
            t = i / num_sections
            this_chord = root_chord + t * (tip_chord - root_chord)

            if i == 0:
                result[i, :] = [0.0, 0.0, 0.0, this_chord]
            else:
                result[i, :] = [0.0, dy, 0.0, this_chord]

        return result

    return resampler


def create_elliptical_wing_resampler(
    span: float,
    root_chord: float,
) -> Callable[[int, int], np.ndarray]:
    """Create a geometry resampler for an elliptical wing.

    The chord distribution follows the elliptical formula:
        c(y) = c_root * sqrt(1 - (2*y/b)^2)

    where b is the span and c_root is the root chord.

    :param span: Total span of the wing.
    :param root_chord: Chord length at the root (center).
    :return: Resampler function that takes (wing_id, num_sections) and returns
        an (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
    """

    def resampler(wing_id: int, num_sections: int) -> np.ndarray:
        """Resample elliptical wing geometry at the specified sections.

        :param wing_id: The wing ID (ignored for elliptical wing).
        :param num_sections: Number of spanwise sections to generate.
        :return: (N+1, 4) ndarray with [dx, dy, dz, chord] per cross section.
        """
        dy = span / num_sections
        result = np.zeros((num_sections + 1, 4))

        for i in range(num_sections + 1):
            # Calculate spanwise position (0 at root, span at tip)
            y = i * dy
            # Elliptical chord distribution
            # Using semi-span since formula uses 2*y/b
            eta = y / span  # normalized spanwise position (0 to 1)
            this_chord = root_chord * np.sqrt(max(0.0, 1.0 - (2 * eta - 1) ** 2))

            # Ensure minimum chord at tip
            this_chord = max(this_chord, root_chord * 0.01)

            if i == 0:
                result[i, :] = [0.0, 0.0, 0.0, this_chord]
            else:
                result[i, :] = [0.0, dy, 0.0, this_chord]

        return result

    return resampler
