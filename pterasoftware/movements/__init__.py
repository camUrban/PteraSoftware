"""Contains the deprecated module paths for the movements classes.

The movements classes are now available at the package top level. Accessing one through
these modules emits a DeprecationWarning, and these modules will be removed in v6.0.0.
"""

from . import (
    aeroelastic_airplane_movement,
    aeroelastic_movement,
    aeroelastic_wing_cross_section_movement,
    aeroelastic_wing_movement,
    airplane_movement,
    free_flight_movement,
    free_flight_operating_point_movement,
    movement,
    operating_point_movement,
    wing_cross_section_movement,
    wing_movement,
)
