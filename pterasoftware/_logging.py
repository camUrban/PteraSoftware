"""Contains the centralized logging configuration for the pterasoftware package."""

from __future__ import annotations

import logging
import sys
from typing import TextIO

from tqdm import tqdm


# Package level logger. All module loggers should be children of this logger.
PACKAGE_LOGGER_NAME = "pterasoftware"


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
    use_tqdm_handler: bool = True,
) -> logging.Logger:
    """Configures logging for the pterasoftware package that is compatible with TQDM
    progress bars.

    This function sets up the package level logger with consistent formatting and
    optionally uses a TQDM compatible handler to prevent progress bar interference.

    :param level: The logging level. Can be an int (e.g., logging.DEBUG) or a string
        (e.g., "Debug", "Info", "Warning", "Error", "Critical"). The default is
        logging.WARNING.
    :param handler: A custom logging handler. If None, a _TqdmLoggingHandler or
        StreamHandler will be created based on use_tqdm_handler.
    :param format_string: Custom format string for log messages. If None, uses a
        sensible default.
    :param use_tqdm_handler: Determines whether to use _TqdmLoggingHandler to prevent
        interference with progress bars. If False, uses standard StreamHandler behavior.
        The default is True.
    :return: The configured package level logger.
    """
    # Convert string level to an int if needed.
    if isinstance(level, str):
        level = _convert_logging_level_name_to_value(level)

    # Get the package level logger.
    logger = logging.getLogger(PACKAGE_LOGGER_NAME)

    # Clear any existing handlers to avoid duplicates.
    logger.handlers.clear()

    # Set the level.
    logger.setLevel(level)

    # Create handler if not provided.
    if handler is None:
        if use_tqdm_handler:
            handler = _TqdmLoggingHandler()
        else:
            handler = logging.StreamHandler(sys.stderr)

    # Set up formatting.
    if format_string is None:
        format_string = "%(levelname)s - %(name)s - %(message)s"

    formatter = logging.Formatter(format_string)
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
        "Debug", "Info", "Warning", "Error", and "Critical".
    :return: The int that can be used to set the appropriate logging level.
    """
    logging_levels = {
        "Debug": logging.DEBUG,
        "Info": logging.INFO,
        "Warning": logging.WARNING,
        "Error": logging.ERROR,
        "Critical": logging.CRITICAL,
    }
    try:
        return logging_levels[name]
    except KeyError:
        raise ValueError(f"{name} is not a valid value of name.")
