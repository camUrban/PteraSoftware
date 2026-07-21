"""Tests for the test environment's serialization warning suppression."""

import unittest
from unittest.mock import patch

# noinspection PyProtectedMember
from pterasoftware import _serialization


class TestSuppressDirtyProvenanceWarnings(unittest.TestCase):
    """Tests for the dirty tree warning filter that tests/_test_environment.py installs.

    The filter is installed as an import side effect of the tests package, so these
    tests exercise it through the production warning paths in
    pterasoftware._serialization. Driving the real paths pins the suppression to the
    real message text, so a future rewording that defeats the filter fails these tests
    instead of quietly re-polluting the test output.
    """

    def test_dirty_file_load_warning_is_suppressed(self) -> None:
        """Loading provenance flagged _dirty should log no warning."""
        with self.assertNoLogs("pterasoftware._serialization", level="WARNING"):
            _serialization._log_load_warnings({"_dirty": True})

    def test_dirty_working_tree_load_warning_is_suppressed(self) -> None:
        """Loading with a dirty working tree should log no warning."""
        # The first check_output call is git rev-parse HEAD, and the second is git
        # status --porcelain. Matching the stored and current commits keeps the
        # unrelated commit-mismatch warning out of this test, and the non-empty status
        # triggers the dirty working tree warning.
        with patch.object(
            _serialization.subprocess,
            "check_output",
            side_effect=[b"abc123\n", b" M some_file.py\n"],
        ):
            with self.assertNoLogs("pterasoftware._serialization", level="WARNING"):
                _serialization._log_load_warnings(
                    {"_dirty": False, "_commit": "abc123"}
                )

    def test_unrelated_serialization_warnings_still_pass(self) -> None:
        """The filter should not suppress other serialization warnings."""
        # The commit-mismatch warning shares the provenance path but is deliberately
        # outside the suppression's scope.
        with patch.object(
            _serialization.subprocess,
            "check_output",
            side_effect=[b"def456\n", b""],
        ):
            with self.assertLogs(
                "pterasoftware._serialization", level="WARNING"
            ) as logs:
                _serialization._log_load_warnings(
                    {"_dirty": False, "_commit": "abc123"}
                )
        self.assertEqual(len(logs.output), 1)
        self.assertIn("saved at commit", logs.output[0])
