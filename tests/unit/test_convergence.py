"""This module contains classes to test the convergence analysis functions."""

import unittest

# noinspection PyProtectedMember
from pterasoftware import convergence


class TestConvergedParameterId(unittest.TestCase):
    """This class contains methods for testing
    convergence._converged_parameter_id.
    """

    def test_single_returns_this_id(self):
        """Test that a single tested value returns this iteration's own index."""
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=0, single=True, converged=False
            ),
            0,
        )

    def test_single_takes_precedence_over_converged(self):
        """Test that a single tested value returns this iteration's own index even
        when the converged flag is set.
        """
        self.assertEqual(
            convergence._converged_parameter_id(this_id=3, single=True, converged=True),
            3,
        )

    def test_converged_returns_coarser_id(self):
        """Test that a converged iteration returns the incrementally coarser
        index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=4, single=False, converged=True
            ),
            3,
        )

    def test_saturated_returns_this_id(self):
        """Test that an iteration that passed without converging returns this
        iteration's own index.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=5, single=False, converged=False
            ),
            5,
        )

    def test_converged_id_depends_on_this_id(self):
        """Test that the coarser index tracks this iteration's index rather than a
        fixed value.
        """
        self.assertEqual(
            convergence._converged_parameter_id(
                this_id=7, single=False, converged=True
            ),
            6,
        )

    def test_returns_int(self):
        """Test that the converged index is returned as an int."""
        result = convergence._converged_parameter_id(
            this_id=2, single=False, converged=True
        )
        self.assertIsInstance(result, int)
