"""Centralized logging configuration for pterasoftware.

This module provides a consistent logging infrastructure with:
- Hierarchical logger naming (pterasoftware.module.submodule)
- TQDM-compatible logging handler
- User-configurable logging setup
"""

from __future__ import annotations

import logging
import sys
from typing import TextIO

from tqdm import tqdm


# Package-level logger - all module loggers should be children of this
PACKAGE_LOGGER_NAME = "pterasoftware"

# Track whether logging has been configured
_logging_configured = False


class TqdmLoggingHandler(logging.Handler):
    """A logging handler that writes messages through tqdm.write().

    This prevents log messages from breaking TQDM progress bars by using
    tqdm's write mechanism which properly handles terminal output.
    """

    def __init__(self, stream: TextIO | None = None) -> None:
        """Initialize the handler.

        :param stream: The stream to write to. Defaults to sys.stderr.
        """
        super().__init__()
        self.stream = stream or sys.stderr

    def emit(self, record: logging.LogRecord) -> None:
        """Emit a log record using tqdm.write().

        :param record: The log record to emit.
        """
        try:
            msg = self.format(record)
            tqdm.write(msg, file=self.stream)
            self.flush()
        except Exception:
            self.handleError(record)

    def flush(self) -> None:
        """Flush the stream."""
        if self.stream and hasattr(self.stream, "flush"):
            self.stream.flush()


def get_logger(name: str) -> logging.Logger:
    """Get a logger with proper hierarchical naming.

    This function ensures all loggers are children of the pterasoftware
    package logger, enabling unified configuration.

    :param name: The module name (e.g., 'trim', 'convergence', 'movements.movement').
        Should use dots for hierarchy, not slashes.
    :return: A properly configured logger instance.

    Example:
        logger = get_logger("trim")  # Returns logger named "pterasoftware.trim"
        logger = get_logger("movements.movement")  # Returns "pterasoftware.movements.movement"
    """
    if name.startswith(PACKAGE_LOGGER_NAME):
        return logging.getLogger(name)
    return logging.getLogger(f"{PACKAGE_LOGGER_NAME}.{name}")


def setup_logging(
    level: int | str = logging.WARNING,
    handler: logging.Handler | None = None,
    format_string: str | None = None,
    use_tqdm_handler: bool = True,
) -> logging.Logger:
    """Configure logging for the pterasoftware package.

    This function sets up the package-level logger with consistent formatting
    and optionally uses a TQDM-compatible handler to prevent progress bar
    interference.

    :param level: The logging level. Can be an int (e.g., logging.DEBUG) or
        a string (e.g., "Debug", "Info", "Warning", "Error", "Critical").
        Default is logging.WARNING.
    :param handler: A custom logging handler. If None, a TqdmLoggingHandler
        or StreamHandler will be created based on use_tqdm_handler.
    :param format_string: Custom format string for log messages. If None,
        uses a sensible default.
    :param use_tqdm_handler: If True (default), uses TqdmLoggingHandler to
        prevent interference with progress bars. Set to False for standard
        StreamHandler behavior.
    :return: The configured package-level logger.

    Example:
        # Basic setup with default WARNING level
        setup_logging()

        # Debug level with TQDM-compatible output
        setup_logging(level=logging.DEBUG)

        # Using string level names (matches solver API)
        setup_logging(level="Info")
    """
    global _logging_configured

    # Convert string level to int if needed
    if isinstance(level, str):
        level = convert_logging_level_name_to_value(level)

    # Get the package-level logger
    logger = logging.getLogger(PACKAGE_LOGGER_NAME)

    # Clear any existing handlers to avoid duplicates
    logger.handlers.clear()

    # Set the level
    logger.setLevel(level)

    # Create handler if not provided
    if handler is None:
        if use_tqdm_handler:
            handler = TqdmLoggingHandler()
        else:
            handler = logging.StreamHandler(sys.stderr)

    # Set up formatting
    if format_string is None:
        format_string = "%(levelname)s - %(name)s - %(message)s"

    formatter = logging.Formatter(format_string)
    handler.setFormatter(formatter)
    handler.setLevel(level)

    # Add handler to logger
    logger.addHandler(handler)

    # Prevent propagation to root logger to avoid duplicate messages
    logger.propagate = False

    _logging_configured = True

    return logger


def convert_logging_level_name_to_value(name: str) -> int:
    """Convert a logging level name string to its integer value.

    :param name: The string representation of the logging level. The options are
        "Debug", "Info", "Warning", "Error", and "Critical".
    :return: The int that can be used to set the appropriate logging level.
    :raises ValueError: If the name is not a valid logging level.
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


def ensure_logging_configured(level: int = logging.WARNING) -> None:
    """Ensure logging is configured, setting up defaults if not.

    This is called internally by solvers to ensure logging works even if
    the user hasn't explicitly called setup_logging().

    :param level: The logging level to use if setting up defaults.
    """
    global _logging_configured
    if not _logging_configured:
        setup_logging(level=level)
    else:
        # Just update the level if already configured
        logger = logging.getLogger(PACKAGE_LOGGER_NAME)
        logger.setLevel(level)
        for handler in logger.handlers:
            handler.setLevel(level)


def is_logging_configured() -> bool:
    """Check if logging has been explicitly configured.

    :return: True if setup_logging() has been called, False otherwise.
    """
    return _logging_configured
