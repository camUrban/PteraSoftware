"""Tests for the centralized logging module."""

import io
import logging
import unittest

# noinspection PyProtectedMember
from pterasoftware import _logging


class TestGetLogger(unittest.TestCase):
    """Tests for the get_logger function."""

    def test_returns_logger_with_package_prefix(self):
        """get_logger should return a logger with the pterasoftware prefix."""
        logger = _logging.get_logger("trim")
        self.assertEqual(logger.name, "pterasoftware.trim")

    def test_hierarchical_naming_with_dots(self):
        """get_logger should handle hierarchical names with dots."""
        logger = _logging.get_logger("movements.movement")
        self.assertEqual(logger.name, "pterasoftware.movements.movement")

    def test_already_prefixed_name_not_duplicated(self):
        """get_logger should not duplicate prefix if already present."""
        logger = _logging.get_logger("pterasoftware.test")
        self.assertEqual(logger.name, "pterasoftware.test")


class TestConvertLoggingLevelNameToValue(unittest.TestCase):
    """Tests for the convert_logging_level_name_to_value function."""

    def test_debug_level(self):
        """Should convert 'Debug' to logging.DEBUG."""
        self.assertEqual(
            _logging.convert_logging_level_name_to_value("Debug"),
            logging.DEBUG
        )

    def test_info_level(self):
        """Should convert 'Info' to logging.INFO."""
        self.assertEqual(
            _logging.convert_logging_level_name_to_value("Info"),
            logging.INFO
        )

    def test_warning_level(self):
        """Should convert 'Warning' to logging.WARNING."""
        self.assertEqual(
            _logging.convert_logging_level_name_to_value("Warning"),
            logging.WARNING
        )

    def test_error_level(self):
        """Should convert 'Error' to logging.ERROR."""
        self.assertEqual(
            _logging.convert_logging_level_name_to_value("Error"),
            logging.ERROR
        )

    def test_critical_level(self):
        """Should convert 'Critical' to logging.CRITICAL."""
        self.assertEqual(
            _logging.convert_logging_level_name_to_value("Critical"),
            logging.CRITICAL
        )

    def test_invalid_level_raises_value_error(self):
        """Should raise ValueError for invalid level names."""
        with self.assertRaises(ValueError):
            _logging.convert_logging_level_name_to_value("InvalidLevel")


class TestTqdmLoggingHandler(unittest.TestCase):
    """Tests for the TqdmLoggingHandler class."""

    def test_handler_inherits_from_logging_handler(self):
        """TqdmLoggingHandler should inherit from logging.Handler."""
        handler = _logging.TqdmLoggingHandler()
        self.assertIsInstance(handler, logging.Handler)

    def test_handler_uses_provided_stream(self):
        """TqdmLoggingHandler should use the provided stream."""
        stream = io.StringIO()
        handler = _logging.TqdmLoggingHandler(stream=stream)
        self.assertEqual(handler.stream, stream)


class TestSetupLogging(unittest.TestCase):
    """Tests for the setup_logging function."""

    def setUp(self):
        """Reset logging configuration before each test."""
        # Reset the configured flag
        _logging._logging_configured = False
        # Clear any handlers from the package logger
        pkg_logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        pkg_logger.handlers.clear()

    def test_returns_package_logger(self):
        """setup_logging should return the package-level logger."""
        logger = _logging.setup_logging()
        self.assertEqual(logger.name, _logging.PACKAGE_LOGGER_NAME)

    def test_sets_configured_flag(self):
        """setup_logging should set the configured flag."""
        _logging.setup_logging()
        self.assertTrue(_logging.is_logging_configured())

    def test_accepts_int_level(self):
        """setup_logging should accept integer log levels."""
        logger = _logging.setup_logging(level=logging.DEBUG)
        self.assertEqual(logger.level, logging.DEBUG)

    def test_accepts_string_level(self):
        """setup_logging should accept string log levels."""
        logger = _logging.setup_logging(level="Info")
        self.assertEqual(logger.level, logging.INFO)

    def test_clears_existing_handlers(self):
        """setup_logging should clear existing handlers to avoid duplicates."""
        _logging.setup_logging()
        _logging.setup_logging()
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertEqual(len(logger.handlers), 1)

    def test_uses_tqdm_handler_by_default(self):
        """setup_logging should use TqdmLoggingHandler by default."""
        _logging.setup_logging()
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertIsInstance(logger.handlers[0], _logging.TqdmLoggingHandler)

    def test_uses_stream_handler_when_tqdm_disabled(self):
        """setup_logging should use StreamHandler when use_tqdm_handler=False."""
        _logging.setup_logging(use_tqdm_handler=False)
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertIsInstance(logger.handlers[0], logging.StreamHandler)
        self.assertNotIsInstance(logger.handlers[0], _logging.TqdmLoggingHandler)


class TestEnsureLoggingConfigured(unittest.TestCase):
    """Tests for the ensure_logging_configured function."""

    def setUp(self):
        """Reset logging configuration before each test."""
        _logging._logging_configured = False
        pkg_logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        pkg_logger.handlers.clear()

    def test_configures_logging_if_not_configured(self):
        """ensure_logging_configured should set up logging if not already done."""
        _logging.ensure_logging_configured()
        self.assertTrue(_logging.is_logging_configured())

    def test_uses_provided_level(self):
        """ensure_logging_configured should use the provided level."""
        _logging.ensure_logging_configured(level=logging.DEBUG)
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertEqual(logger.level, logging.DEBUG)

    def test_updates_level_if_already_configured(self):
        """ensure_logging_configured should update level if already configured."""
        _logging.setup_logging(level=logging.WARNING)
        _logging.ensure_logging_configured(level=logging.DEBUG)
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertEqual(logger.level, logging.DEBUG)


class TestLoggerHierarchy(unittest.TestCase):
    """Tests for proper logger hierarchy behavior."""

    def setUp(self):
        """Reset logging configuration before each test."""
        _logging._logging_configured = False
        pkg_logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        pkg_logger.handlers.clear()

    def test_child_logger_inherits_from_package_logger(self):
        """Child loggers should inherit settings from package logger."""
        _logging.setup_logging(level=logging.DEBUG)
        child_logger = _logging.get_logger("test_child")

        # Child should be able to log at DEBUG level
        self.assertTrue(child_logger.isEnabledFor(logging.DEBUG))

    def test_child_logger_messages_go_to_package_handler(self):
        """Child logger messages should go through package logger handlers."""
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        handler.setFormatter(logging.Formatter("%(name)s - %(message)s"))

        _logging.setup_logging(level=logging.INFO, handler=handler)
        child_logger = _logging.get_logger("test_module")

        child_logger.info("Test message")
        output = stream.getvalue()

        self.assertIn("pterasoftware.test_module", output)
        self.assertIn("Test message", output)


if __name__ == "__main__":
    unittest.main()
