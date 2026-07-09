"""Tests for the centralized logging module."""

import io
import logging
import unittest

# noinspection PyProtectedMember
from pterasoftware import _logging


class TestIndent(unittest.TestCase):
    """Tests for the indent function."""

    def test_returns_empty_string_at_default_nesting_level(self):
        """indent should return an empty string at the default nesting level."""
        self.assertEqual(_logging.indent(), "")

    def test_returns_two_spaces_per_level(self):
        """indent should return two spaces per level beyond the nesting level."""
        self.assertEqual(_logging.indent(3), "      ")


class TestNested(unittest.TestCase):
    """Tests for the nested context manager."""

    def test_deepens_nesting_level_by_one_by_default(self):
        """nested should deepen the nesting level by one level by default."""
        with _logging.nested():
            self.assertEqual(_logging.indent(), "  ")

    def test_deepens_nesting_level_by_given_levels(self):
        """nested should deepen the nesting level by the given number of levels."""
        with _logging.nested(3):
            self.assertEqual(_logging.indent(), "      ")

    def test_nested_blocks_accumulate(self):
        """nested blocks should accumulate their levels."""
        with _logging.nested():
            with _logging.nested(2):
                self.assertEqual(_logging.indent(), "      ")

    def test_indent_levels_add_to_nesting_level(self):
        """indent's levels should add to the current nesting level."""
        with _logging.nested():
            self.assertEqual(_logging.indent(2), "      ")

    def test_restores_nesting_level_after_block(self):
        """nested should restore the nesting level after the block."""
        with _logging.nested(2):
            pass
        self.assertEqual(_logging.indent(), "")

    def test_restores_nesting_level_after_exception(self):
        """nested should restore the nesting level when the block raises."""
        with self.assertRaises(ValueError):
            with _logging.nested():
                raise ValueError("Test error")
        self.assertEqual(_logging.indent(), "")


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
            _logging._convert_logging_level_name_to_value("Debug"), logging.DEBUG
        )

    def test_info_level(self):
        """Should convert 'Info' to logging.INFO."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("Info"), logging.INFO
        )

    def test_warning_level(self):
        """Should convert 'Warning' to logging.WARNING."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("Warning"), logging.WARNING
        )

    def test_error_level(self):
        """Should convert 'Error' to logging.ERROR."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("Error"), logging.ERROR
        )

    def test_critical_level(self):
        """Should convert 'Critical' to logging.CRITICAL."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("Critical"), logging.CRITICAL
        )

    def test_invalid_level_raises_value_error(self):
        """Should raise ValueError for invalid level names."""
        with self.assertRaises(ValueError):
            _logging._convert_logging_level_name_to_value("InvalidLevel")

    def test_case_insensitive_lowercase(self):
        """Should accept lowercase level names."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("debug"), logging.DEBUG
        )

    def test_case_insensitive_uppercase(self):
        """Should accept uppercase level names."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("WARNING"), logging.WARNING
        )

    def test_case_insensitive_mixed_case(self):
        """Should accept mixed case level names."""
        self.assertEqual(
            _logging._convert_logging_level_name_to_value("CrItIcAl"), logging.CRITICAL
        )


class TestTqdmLoggingHandler(unittest.TestCase):
    """Tests for the _TqdmLoggingHandler class."""

    def test_handler_inherits_from_logging_handler(self):
        """_TqdmLoggingHandler should inherit from logging.Handler."""
        handler = _logging._TqdmLoggingHandler()
        self.assertIsInstance(handler, logging.Handler)

    def test_handler_uses_provided_stream(self):
        """_TqdmLoggingHandler should use the provided stream."""
        stream = io.StringIO()
        handler = _logging._TqdmLoggingHandler(stream=stream)
        self.assertEqual(handler.stream, stream)

    def test_handler_emits_formatted_log_record(self):
        """_TqdmLoggingHandler should emit formatted log records."""
        stream = io.StringIO()
        handler = _logging._TqdmLoggingHandler(stream=stream)
        handler.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))

        # Create and emit a log record.
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None,
        )
        handler.emit(record)

        # Verify the message was written.
        output = stream.getvalue()
        self.assertIn("INFO - Test message", output)

    def test_handler_flush_works(self):
        """_TqdmLoggingHandler should flush the stream when flush() is called."""
        stream = io.StringIO()
        handler = _logging._TqdmLoggingHandler(stream=stream)
        handler.setFormatter(logging.Formatter("%(message)s"))

        # Emit a log record.
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test",
            args=(),
            exc_info=None,
        )
        handler.emit(record)

        # Call flush explicitly.
        handler.flush()

        # Verify stream was flushed (content should be available).
        output = stream.getvalue()
        self.assertIn("Test", output)


class TestMaxModuleLoggerDisplayNameLength(unittest.TestCase):
    """Tests for the _max_module_logger_display_name_length function."""

    def test_returns_at_least_the_package_logger_name_length(self):
        """Should return at least the length of the package logger name."""
        self.assertGreaterEqual(
            _logging._max_module_logger_display_name_length(),
            len(_logging.PACKAGE_LOGGER_NAME),
        )

    def test_returns_the_longest_display_name_length(self):
        """Should return the length of the longest module logger display name."""
        self.assertEqual(
            _logging._max_module_logger_display_name_length(),
            len("_coupled_unsteady_ring_vortex_lattice_method"),
        )


class TestPackageLogFormatter(unittest.TestCase):
    """Tests for the _PackageLogFormatter class."""

    def test_strips_package_prefix_from_display_name(self):
        """_PackageLogFormatter should strip the package prefix from the name."""
        formatter = _logging._PackageLogFormatter("%(display_name)s|%(message)s")
        record = logging.LogRecord(
            name="pterasoftware.trim",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None,
        )
        self.assertEqual(formatter.format(record), "trim|Test message")

    def test_keeps_hierarchy_below_package_prefix(self):
        """_PackageLogFormatter should keep the hierarchy below the prefix."""
        formatter = _logging._PackageLogFormatter("%(display_name)s|%(message)s")
        record = logging.LogRecord(
            name="pterasoftware.movements.movement",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None,
        )
        self.assertEqual(formatter.format(record), "movements.movement|Test message")

    def test_keeps_bare_package_logger_name(self):
        """_PackageLogFormatter should leave the bare package name unchanged."""
        formatter = _logging._PackageLogFormatter("%(display_name)s|%(message)s")
        record = logging.LogRecord(
            name="pterasoftware",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None,
        )
        self.assertEqual(formatter.format(record), "pterasoftware|Test message")


class TestSetupLogging(unittest.TestCase):
    """Tests for the set_up_logging function."""

    def setUp(self):
        """Reset logging configuration before each test."""
        # Clear any handlers from the package logger.
        pkg_logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        pkg_logger.handlers.clear()

    def test_returns_package_logger(self):
        """set_up_logging should return the package-level logger."""
        logger = _logging.set_up_logging()
        self.assertEqual(logger.name, _logging.PACKAGE_LOGGER_NAME)

    def test_accepts_int_level(self):
        """set_up_logging should accept integer log levels."""
        logger = _logging.set_up_logging(level=logging.DEBUG)
        self.assertEqual(logger.level, logging.DEBUG)

    def test_accepts_string_level(self):
        """set_up_logging should accept string log levels."""
        logger = _logging.set_up_logging(level="Info")
        self.assertEqual(logger.level, logging.INFO)

    def test_clears_existing_handlers(self):
        """set_up_logging should clear existing handlers to avoid duplicates."""
        _logging.set_up_logging()
        _logging.set_up_logging()
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertEqual(len(logger.handlers), 1)

    def test_uses_tqdm_handler_by_default(self):
        """set_up_logging should use _TqdmLoggingHandler when no handler is provided."""
        _logging.set_up_logging()
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertIsInstance(logger.handlers[0], _logging._TqdmLoggingHandler)

    def test_uses_custom_handler_when_provided(self):
        """set_up_logging should use the provided handler instead of the default."""
        custom_handler = logging.StreamHandler()
        _logging.set_up_logging(handler=custom_handler)
        logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        self.assertIs(logger.handlers[0], custom_handler)

    def test_default_format_right_pads_level_names(self):
        """set_up_logging's default format should right pad the level names."""
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        _logging.set_up_logging(level=logging.INFO, handler=handler)

        _logging.get_logger("trim").info("Test message")

        self.assertTrue(stream.getvalue().startswith("INFO    |"))

    def test_default_format_aligns_messages_across_modules(self):
        """set_up_logging's default format should align messages across modules."""
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        _logging.set_up_logging(level=logging.INFO, handler=handler)

        _logging.get_logger("trim").info("First message")
        _logging.get_logger("movements.movement").info("Second message")
        first_line, second_line = stream.getvalue().splitlines()

        self.assertEqual(
            first_line.index("First message"), second_line.index("Second message")
        )

    def test_invalid_level_type_raises_type_error(self):
        """set_up_logging should raise TypeError for invalid level types."""
        with self.assertRaises(TypeError) as context:
            _logging.set_up_logging(level=3.14)
        self.assertIn("level must be an int or a str", str(context.exception))

    def test_invalid_handler_type_raises_type_error(self):
        """set_up_logging should raise TypeError for invalid handler types."""
        with self.assertRaises(TypeError) as context:
            _logging.set_up_logging(handler="not a handler")
        self.assertIn(
            "handler must be a logging.Handler or None", str(context.exception)
        )

    def test_invalid_format_string_type_raises_type_error(self):
        """set_up_logging should raise TypeError for invalid format_string types."""
        with self.assertRaises(TypeError) as context:
            _logging.set_up_logging(format_string=123)
        self.assertIn("format_string must be a str or None", str(context.exception))


class TestLoggerHierarchy(unittest.TestCase):
    """Tests for proper logger hierarchy behavior."""

    def setUp(self):
        """Reset logging configuration before each test."""
        pkg_logger = logging.getLogger(_logging.PACKAGE_LOGGER_NAME)
        pkg_logger.handlers.clear()

    def test_child_logger_inherits_from_package_logger(self):
        """Child loggers should inherit settings from package logger."""
        _logging.set_up_logging(level=logging.DEBUG)
        child_logger = _logging.get_logger("test_child")

        # Child should be able to log at DEBUG level.
        self.assertTrue(child_logger.isEnabledFor(logging.DEBUG))

    def test_child_logger_messages_go_to_package_handler(self):
        """Child logger messages should go through package logger handlers."""
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)

        _logging.set_up_logging(level=logging.INFO, handler=handler)
        child_logger = _logging.get_logger("test_module")

        child_logger.info("Test message")
        output = stream.getvalue()

        self.assertIn("test_module", output)
        self.assertIn("Test message", output)
