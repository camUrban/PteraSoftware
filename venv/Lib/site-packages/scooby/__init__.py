# coding=utf-8
"""
A Great Dane turned Python environment detective
================================================

A lightweight toolset to easily report your Python environment's package
versions and hardware resources.

History
-------
The scooby reporting is derived from the versioning-scripts created by Dieter
Werthmüller for ``empymod``, ``emg3d``, and the ``SimPEG`` framework
(https://empymod.github.io; https://simpeg.xyz). It was heavily inspired by
``ipynbtools.py`` from ``qutip`` (https://github.com/qutip) and
``watermark.py`` from https://github.com/rasbt/watermark.
"""

from scooby.report import Report, get_version
from scooby.knowledge import (get_standard_lib_modules, in_ipykernel,
                              in_ipython, meets_version, version_tuple)
from scooby.tracker import TrackedReport, track_imports, untrack_imports

__all__ = ['Report', 'TrackedReport', 'get_standard_lib_modules', 'in_ipython',
           'in_ipykernel', 'get_version', 'track_imports', 'untrack_imports']


__author__ = 'Dieter Werthmüller & Bane Sullivan'
__license__ = 'MIT'
__copyright__ = '2019, Dieter Werthmüller & Bane Sullivan'
__version__ = '0.5.2'
