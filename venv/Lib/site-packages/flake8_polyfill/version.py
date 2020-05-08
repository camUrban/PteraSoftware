"""Version information for Flake8 2.x and 3.x."""

import flake8

version_info = getattr(flake8, '__version_info__', None)
if version_info is None:
    version_info = tuple(
        int(i) for i in flake8.__version__.split('.') if i.isdigit()
    )
