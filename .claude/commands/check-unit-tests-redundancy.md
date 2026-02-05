---
description: Check the unit tests of a specified module for redundancy before debugging
---

# Check Test Redundancy

Review unit tests for `$ARGUMENTS` to identify and remove redundant tests. Do NOT debug or run the tests after doing so.

## Progress Tracking

Follow these steps carefully and track your progress:

- [ ] Read the test file(s) for `$ARGUMENTS`
- [ ] Catalog all tests by what behavior they verify
- [ ] Identify redundant tests (same behavior tested multiple times)
- [ ] Identify tests that could be consolidated
- [ ] Document findings in a redundancy report
- [ ] Remove or consolidate redundant tests
- [ ] Update `__init__.py` files if any files were removed

## Omissions

Please do NOT:
- Run or debug any tests you create
- Reformat any files with black

These steps are both handled by other slash commands.

## Detailed Steps

1. **Read the test and fixture files** for `$ARGUMENTS`:
   - Find corresponding test files in `tests/unit/`
   - Find corresponding fixture files in `tests/unit/fixtures/`
   - Note all test methods and fixture functions

2. **Catalog all tests** by the behavior they verify:
   - Group tests by the method or function they test
   - Note what specific aspect each test verifies
   - Identify the input conditions and expected outcomes

3. **Identify redundant tests**:

   Tests are redundant if they:
   - Test the exact same behavior with the same inputs
   - Test the same edge case multiple times
   - Verify the same error condition in the same way
   - Differ only in naming but not in what they verify

   Tests are NOT redundant if they:
   - Test the same method but with different inputs
   - Test the same behavior but for different edge cases
   - Verify the same error but triggered by different conditions

4. **Document findings** in a redundancy report:
   ```
   REDUNDANCY REPORT for [module]

   REDUNDANT TESTS:
   - [test_name_1] duplicates [test_name_2]: [reason]
   - ...

   CONSOLIDATION OPPORTUNITIES:
   - [test_name_1] and [test_name_2] could be parameterized: [reason]
   - ...
   ```

5. **Remove or consolidate** redundant tests and fixtures:
   - Delete clearly redundant tests
   - Consolidate tests that verify the same thing
   - Keep the most comprehensive version when consolidating

6. **Update `__init__.py` files** if any modules were removed or renamed.

## Guidelines for Judgment Calls

When uncertain if something is redundant:

### Keep both tests if:
- They document different requirements or use cases
- They would fail for different reasons if the source broke
- Removing one would reduce confidence in the code

### Remove/consolidate if:
- One test is strictly a subset of another
- The tests would always pass or fail together
- Keeping both adds maintenance burden without value

## Error Handling

If any step fails or information is missing:
- Explicitly state what's missing
- Ask for clarification before removing tests
- Document any assumptions made
- When in doubt, keep the test and flag for review

## Quality Checklist

Before finalizing:
- [ ] All tests have been reviewed for redundancy
- [ ] Redundancy report has been generated
- [ ] No useful tests were accidentally removed
- [ ] Remaining tests still provide comprehensive coverage
- [ ] `__init__.py` files are updated if needed