"""Configures the test process to quiet known sources of test run noise.

Three categories of noise are addressed here:

1. Two pterasoftware._serialization warnings about uncommitted git state. These
   are correct in production but always fire during test runs because tests are
   commonly executed on a dirty working tree.

2. tqdm progress bars from solver runs invoked through paths that do not
   forward a show_progress=False parameter (such as the convergence and trim
   helpers). Setting TQDM_DISABLE=1 in the environment disables every tqdm
   instance in the process without requiring per call site changes.

3. The VTK "bad X server connection" warning that fires on Linux when no X
   display is reachable (typically over SSH without X11 forwarding). When this
   is detected, the process is switched to PyVista off-screen mode. On macOS,
   Windows, and Linux sessions with a working display (including CI runs that
   use pyvista/setup-headless-display-action to provide xvfb), no change is
   made and plotters render normally.

Importing this module installs the three suppressions as side effects. It is
imported as the first line of tests/__init__.py so the configuration is in
place before any pterasoftware, pyvista, or tqdm module loads.
"""

import logging
import os
import sys

_SUPPRESSED_SERIALIZATION_MESSAGES = frozenset(
    {
        "The current working tree has uncommitted changes.",
        (
            "The file was saved with uncommitted changes (_dirty=true). The "
            "_commit hash may not fully represent the code state at save time."
        ),
    }
)


class _SuppressDirtyProvenanceWarnings(logging.Filter):
    """Drops the two known pterasoftware._serialization dirty tree warnings."""

    def filter(self, record: logging.LogRecord) -> bool:
        return record.getMessage() not in _SUPPRESSED_SERIALIZATION_MESSAGES


def _install_serialization_log_filter() -> None:
    logger = logging.getLogger("pterasoftware._serialization")
    for existing in logger.filters:
        if isinstance(existing, _SuppressDirtyProvenanceWarnings):
            return
    logger.addFilter(_SuppressDirtyProvenanceWarnings())


def _disable_tqdm_progress_bars() -> None:
    """Forces every tqdm instance to be silent.

    Setting TQDM_DISABLE=1 only takes effect when the caller passes
    disable=None. pterasoftware passes disable=not show_progress, which is an
    explicit False whenever show_progress is True (the default for several
    helpers like convergence and trim). Monkey-patching tqdm.__init__ ensures
    disable=True regardless of caller intent.
    """
    os.environ.setdefault("TQDM_DISABLE", "1")
    try:
        import tqdm
    except ImportError:
        return

    original_init = tqdm.tqdm.__init__

    def _silent_init(self, *args, **kwargs):
        kwargs["disable"] = True
        return original_init(self, *args, **kwargs)

    tqdm.tqdm.__init__ = _silent_init


def _enable_pyvista_off_screen_if_no_display() -> None:
    """Switches PyVista to off-screen and silences VTK warnings on headless Linux.

    PyVista off-screen mode prevents window creation, but VTK still emits a
    "bad X server connection" warning from its render window constructor on
    Linux when DISPLAY is unset, even with off-screen enabled. Dropping the
    VTK stderr verbosity to ERROR suppresses that warning while keeping
    genuine VTK errors visible. This is scoped to headless Linux only;
    sessions with a display and other platforms are left untouched.
    """
    if not sys.platform.startswith("linux"):
        return
    if os.environ.get("DISPLAY"):
        return
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
    try:
        import pyvista as pv
    except ImportError:
        return
    pv.OFF_SCREEN = True
    try:
        from vtkmodules.vtkCommonCore import vtkLogger
    except ImportError:
        return
    vtkLogger.SetStderrVerbosity(vtkLogger.VERBOSITY_ERROR)


_install_serialization_log_filter()
_disable_tqdm_progress_bars()
_enable_pyvista_off_screen_if_no_display()
