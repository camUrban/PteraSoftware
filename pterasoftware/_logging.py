"""Contains the centralized logging configuration for the pterasoftware package."""

from __future__ import annotations

import contextvars
import logging
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

from tqdm import tqdm

# Package level logger. All module loggers should be children of this logger.
PACKAGE_LOGGER_NAME = "pterasoftware"

# The current log-message nesting level. Each level indents messages formatted with
# indent() by two spaces.
_indent_level: contextvars.ContextVar[int] = contextvars.ContextVar(
    "indent_level", default=0
)


def indent(levels: int = 0) -> str:
    """Returns the leading whitespace for a log message at the current nesting level.

    :param levels: The number of levels to indent the message beyond the current nesting
        level. It must be a non negative int. The default is 0.
    :return: The leading whitespace, two spaces per level.
    """
    return "  " * (_indent_level.get() + levels)


@contextmanager
def nested(levels: int = 1) -> Iterator[None]:
    """Deepens the log-message nesting level for the duration of a with block.

    Log messages formatted with indent() inside the with block are indented by the given
    number of additional levels. The nesting level is tracked with a context variable,
    so a function's log messages nest under its caller's without the function knowing
    its caller's nesting level.

    :param levels: The number of levels to deepen the nesting level by. It must be a
        positive int. The default is 1.
    :return: None
    """
    token = _indent_level.set(_indent_level.get() + levels)
    try:
        yield
    finally:
        _indent_level.reset(token)


class _TqdmLoggingHandler(logging.Handler):
    """A logging handler that writes messages through tqdm.write().

    **Contains the following methods:**

    emit: Emits a log record using tqdm.write().

    flush: Flushes the stream.

    **Notes:**

    This prevents log messages from breaking TQDM progress bars by using tqdm's write
    mechanism which properly handles terminal output.
    """

    def __init__(self, stream: TextIO | None = None) -> None:
        """The initialization method.

        :param stream: The stream to write to. Defaults to sys.stderr.
        :return: None
        """
        super().__init__()
        self.stream = stream or sys.stderr

    def emit(self, record: logging.LogRecord) -> None:
        """Emits a log record using tqdm.write().

        :param record: The log record to emit.
        :return: None
        """
        # Catch all exceptions to prevent logging failures from crashing the program.
        # noinspection PyBroadException
        try:
            msg = self.format(record)
            tqdm.write(msg, file=self.stream)
            self.flush()
        except Exception:
            self.handleError(record)

    def flush(self) -> None:
        """Flushes the stream.

        :return: None
        """
        if self.stream and hasattr(self.stream, "flush"):
            self.stream.flush()


def _max_module_logger_display_name_length() -> int:
    """Returns the length of the longest possible module logger display name in the
    package.

    A logger's display name is its name with the package prefix stripped (see
    _PackageLogFormatter). Loggers are named after their modules via get_logger, so the
    package's module file names determine every display name it can create. The modules
    are found by walking the package's source tree rather than by importing them, which
    keeps this compatible with the package's lazy imports and free of import side
    effects. Only modules that create a logger contribute to the width, and a module is
    detected as creating one by its source containing the "_logging.get_logger(" call
    form, which every module logger definition uses.

    :return: The length of the longest possible module logger display name.
    """
    package_dir = Path(__file__).parent
    name_lengths = [len(PACKAGE_LOGGER_NAME)]
    for module_path in package_dir.rglob("*.py"):
        relative_path = module_path.relative_to(package_dir)
        if relative_path.name == "__init__.py":
            continue
        if "_logging.get_logger(" not in module_path.read_text():
            continue
        display_name = ".".join((*relative_path.parts[:-1], relative_path.stem))
        name_lengths.append(len(display_name))
    return max(name_lengths)


class _PackageLogFormatter(logging.Formatter):
    """A logging formatter that adds each record's display name.

    **Contains the following methods:**

    format: Formats a log record after adding its display name.

    **Notes:**

    The display name is the logger's name with the leading "pterasoftware." stripped.
    Every line of a Ptera Software log would otherwise carry that constant prefix, so
    stripping it from the displayed name keeps the fixed part of each line compact. The
    display name is available to format strings as %(display_name)s.
    """

    def format(self, record: logging.LogRecord) -> str:
        """Formats a log record after adding its display name.

        :param record: The log record to format.
        :return: The formatted log message.
        """
        record.display_name = record.name.removeprefix(f"{PACKAGE_LOGGER_NAME}.")
        return super().format(record)


def get_logger(name: str) -> logging.Logger:
    """Gets a logger with proper hierarchical naming.

    This function ensures all loggers are children of the pterasoftware package logger,
    enabling unified configuration. For example, ``get_logger("trim")`` returns a logger
    named "pterasoftware.trim", and ``get_logger("movements.movement")`` returns a
    logger named "pterasoftware.movements.movement".

    :param name: The module name (e.g., "trim", "convergence", "movements.movement").
        Should use dots for hierarchy, not slashes.
    :return: A properly configured logger instance.
    """
    if name.startswith(PACKAGE_LOGGER_NAME):
        return logging.getLogger(name)
    return logging.getLogger(f"{PACKAGE_LOGGER_NAME}.{name}")


def set_up_logging(
    level: int | str = logging.WARNING,
    handler: logging.Handler | None = None,
    format_string: str | None = None,
) -> logging.Logger:
    """Configures logging for the pterasoftware package.

    Call this function when you want to see log messages from pterasoftware. Without it,
    log messages are silently discarded (the default Python behavior for libraries).

    If you want log messages printed to the console, call this function without passing
    a handler (e.g., ``ps.set_up_logging(level="Info")``). This installs a handler that
    is compatible with TQDM progress bars, so log messages and progress bars can coexist
    on the console without garbling each other.

    If you want log messages written to a file (while any progress bars continue to
    display on the console), pass a ``logging.FileHandler`` (e.g.,
    ``ps.set_up_logging(level="Info", handler=logging.FileHandler("simulation.log"))``).

    :param level: The logging level. Can be an int (e.g., logging.DEBUG) or a string
        (either "debug", "info", "warning", "error", or "critical", case insensitive).
        The default is logging.WARNING.
    :param handler: A custom logging handler. If None, a TQDM compatible handler is
        created that prints to the console without interfering with progress bars. Pass
        a handler (such as ``logging.FileHandler``) to route log messages elsewhere.
    :param format_string: Custom format string for log messages. If None, uses the
        format "%(levelname)-8s|%(display_name)-<width>s|%(message)s", where
        display_name is the logger's name with the leading "pterasoftware." stripped.
        The level name is right padded to the width of the longest standard level name
        (CRITICAL) and the display name to the width of the package's longest possible
        module logger display name, so the messages align across levels and modules.
        Custom format strings may also reference %(display_name)s.
    :return: The configured package level logger.
    """
    # Validate level.
    if isinstance(level, str):
        level = _convert_logging_level_name_to_value(level)
    elif not isinstance(level, int):
        raise TypeError("level must be an int or a str.")

    # Validate handler.
    if handler is not None and not isinstance(handler, logging.Handler):
        raise TypeError("handler must be a logging.Handler or None.")

    # Validate format_string.
    if format_string is not None and not isinstance(format_string, str):
        raise TypeError("format_string must be a str or None.")

    # Get the package level logger.
    logger = logging.getLogger(PACKAGE_LOGGER_NAME)

    # Clear any existing handlers to avoid duplicates.
    logger.handlers.clear()

    # Set the level.
    logger.setLevel(level)

    # Create handler if not provided.
    if handler is None:
        handler = _TqdmLoggingHandler()

    # Set up formatting. The level name is right padded to the width of the longest
    # standard level name (CRITICAL) and the display name to the width of the package's
    # longest possible module logger display name, so the messages align across levels
    # and modules. Aligned messages let the indentation of nested messages read as one
    # continuous tree.
    if format_string is None:
        name_width = _max_module_logger_display_name_length()
        format_string = f"%(levelname)-8s|%(display_name)-{name_width}s|%(message)s"

    formatter = _PackageLogFormatter(format_string)
    handler.setFormatter(formatter)
    handler.setLevel(level)

    # Add handler to logger.
    logger.addHandler(handler)

    # Prevent propagation to root logger to avoid duplicate messages.
    logger.propagate = False

    return logger


def _convert_logging_level_name_to_value(name: str) -> int:
    """Converts a logging level name string to its integer value.

    :param name: The string representation of the logging level. The options are
        "debug", "info", "warning", "error", and "critical" (case insensitive).
    :return: The int that can be used to set the appropriate logging level.
    """
    logging_levels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }
    try:
        return logging_levels[name.lower()]
    except KeyError:
        raise ValueError(f"{name} is not a valid value of name.")
