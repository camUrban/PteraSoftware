"""Contains the deprecated module paths for the geometry classes.

The geometry classes are now available at the package top level. Accessing one through
these modules emits a DeprecationWarning, and these modules will be removed in v6.0.0.
"""

from . import airfoil, airplane, wing, wing_cross_section
