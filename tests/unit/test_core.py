"""This module contains classes to test the core classes in _core.py."""

import unittest

import pterasoftware as ps


class TestLcmFunctions(unittest.TestCase):
    """Tests for _lcm and _lcm_multiple functions in _core.py."""

    def test_lcm_of_two_positive_numbers(self):
        """Test _lcm returns correct LCM for two positive numbers."""
        result = ps._core.lcm(2.0, 3.0)
        self.assertEqual(result, 6.0)

    def test_lcm_of_same_numbers(self):
        """Test _lcm returns the number when both inputs are the same."""
        result = ps._core.lcm(4.0, 4.0)
        self.assertEqual(result, 4.0)

    def test_lcm_with_first_zero(self):
        """Test _lcm returns 0.0 when first input is zero."""
        result = ps._core.lcm(0.0, 5.0)
        self.assertEqual(result, 0.0)

    def test_lcm_with_second_zero(self):
        """Test _lcm returns 0.0 when second input is zero."""
        result = ps._core.lcm(5.0, 0.0)
        self.assertEqual(result, 0.0)

    def test_lcm_with_both_zero(self):
        """Test _lcm returns 0.0 when both inputs are zero."""
        result = ps._core.lcm(0.0, 0.0)
        self.assertEqual(result, 0.0)

    def test_lcm_of_one_and_any_number(self):
        """Test _lcm of 1.0 and any number returns that number."""
        result = ps._core.lcm(1.0, 7.0)
        self.assertEqual(result, 7.0)

        result = ps._core.lcm(7.0, 1.0)
        self.assertEqual(result, 7.0)

    def test_lcm_of_multiples(self):
        """Test _lcm of a number and its multiple returns the larger number."""
        result = ps._core.lcm(3.0, 9.0)
        self.assertEqual(result, 9.0)

        result = ps._core.lcm(9.0, 3.0)
        self.assertEqual(result, 9.0)

    def test_lcm_multiple_empty_list(self):
        """Test _lcm_multiple returns 0.0 for empty list."""
        result = ps._core.lcm_multiple([])
        self.assertEqual(result, 0.0)

    def test_lcm_multiple_all_zeros(self):
        """Test _lcm_multiple returns 0.0 when all periods are zero."""
        result = ps._core.lcm_multiple([0.0, 0.0, 0.0])
        self.assertEqual(result, 0.0)

    def test_lcm_multiple_single_nonzero(self):
        """Test _lcm_multiple returns the value for a single non zero period."""
        result = ps._core.lcm_multiple([5.0])
        self.assertEqual(result, 5.0)

    def test_lcm_multiple_single_zero(self):
        """Test _lcm_multiple returns 0.0 for a single zero period."""
        result = ps._core.lcm_multiple([0.0])
        self.assertEqual(result, 0.0)

    def test_lcm_multiple_mixed_with_zeros(self):
        """Test _lcm_multiple correctly ignores zeros in the list."""
        result = ps._core.lcm_multiple([0.0, 2.0, 0.0, 3.0, 0.0])
        self.assertEqual(result, 6.0)

    def test_lcm_multiple_three_periods(self):
        """Test _lcm_multiple returns correct LCM for three periods."""
        result = ps._core.lcm_multiple([2.0, 3.0, 4.0])
        self.assertEqual(result, 12.0)

    def test_lcm_multiple_many_same_periods(self):
        """Test _lcm_multiple returns the period when all are the same."""
        result = ps._core.lcm_multiple([5.0, 5.0, 5.0, 5.0])
        self.assertEqual(result, 5.0)

    def test_lcm_multiple_coprime_periods(self):
        """Test _lcm_multiple of coprime numbers returns their product."""
        # 2, 3, and 5 are coprime, so LCM = 2 * 3 * 5 = 30.
        result = ps._core.lcm_multiple([2.0, 3.0, 5.0])
        self.assertEqual(result, 30.0)

    def test_lcm_non_integer_periods(self):
        """Test _lcm returns correct LCM for non integer periods."""
        # LCM(1.5, 2.5) = 7.5 (both divide 7.5 evenly: 7.5/1.5=5, 7.5/2.5=3).
        result = ps._core.lcm(1.5, 2.5)
        self.assertAlmostEqual(result, 7.5, places=6)

    def test_lcm_multiple_non_integer_periods(self):
        """Test _lcm_multiple returns correct LCM for non integer periods."""
        # LCM(1.5, 2.0, 2.5) = 30.0.
        # 30.0/1.5=20, 30.0/2.0=15, 30.0/2.5=12.
        result = ps._core.lcm_multiple([1.5, 2.0, 2.5])
        self.assertAlmostEqual(result, 30.0, places=6)

    def test_lcm_small_periods(self):
        """Test _lcm handles small periods correctly without precision issues."""
        # LCM(0.001, 0.002) = 0.002.
        result = ps._core.lcm(0.001, 0.002)
        self.assertAlmostEqual(result, 0.002, places=9)

    def test_lcm_multiple_small_periods(self):
        """Test _lcm_multiple handles small periods correctly."""
        # LCM(0.01, 0.02, 0.03) = 0.06.
        result = ps._core.lcm_multiple([0.01, 0.02, 0.03])
        self.assertAlmostEqual(result, 0.06, places=9)
