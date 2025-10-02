---
description: Run and debug unit tests, carefully analyzing failures for real bugs
---

# Run and Debug Tests

Run and debug unit tests for the contents of $ARGUMENTS, with careful analysis of failures.

## Critical Rule

**NEVER automatically change expected test values to make tests pass.** When a test fails, you must determine if:
1. The test expectation is wrong (fix the test), OR
2. The source code has a bug (STOP and alert the user)

**BEFORE modifying any test**, you MUST explicitly document your reasoning in a TEST BUG ANALYSIS showing why the test (not the source) is incorrect.

If there's ANY possibility you've found a real bug based on the source code logic or documentation, **STOP IMMEDIATELY** and alert the user.

## Progress Tracking

Follow these steps carefully and track your progress:

- [ ] Identify which test files correspond to $ARGUMENTS
- [ ] Run the relevant tests and capture output
- [ ] For each failure, analyze the source code logic
- [ ] Check documentation/docstrings for intended behavior
- [ ] Determine if failure indicates test error or source bug
- [ ] Fix test errors ONLY when certain source is correct
- [ ] STOP and alert user if potential source bug found
- [ ] Verify all tests pass after fixes
- [ ] Document any assumptions or edge cases discovered

## Detailed Steps

1. **Identify test files** for $ARGUMENTS:
   - Find corresponding test files in tests/unit/
   - Note any fixture dependencies

2. **Run initial tests**:
   - Capture and analyze the full output
   - Note all failures and their error messages

3. **For EACH test failure**:
   
   a. **Read the failing test carefully**:
      - Understand what the test is trying to verify
      - Note the expected vs actual values
      - Identify which source method/function is being tested
   
   b. **Analyze the source code**:
      - Read the source implementation in $ARGUMENTS
      - Trace through the logic step by step
      - Check all relevant docstrings and comments
      - Look for any documentation about intended behavior
   
   c. **Determine the root cause**:
      - **If source code logic matches test expectation**: Likely a test implementation issue
      - **If source code logic contradicts test expectation**: Possible source bug - STOP
      - **If unclear or ambiguous**: STOP and ask for clarification
   
   d. **Take appropriate action**:
      - **Test is clearly wrong**: 
        - FIRST, articulate your reasoning:
          ```
          TEST BUG ANALYSIS for [test_name]:
          1. What the test expects: [expected_behavior]
          2. What the source code actually does: [step-by-step trace]
          3. Why the source is correct: [reference to documentation/logic]
          4. Why the test is wrong: [specific reason]
          5. Proposed fix: [what will be changed and why]
          ```
        - ONLY AFTER documenting this analysis, fix the test implementation or expected values
      - **Source might have bug**: 
        ```
        ⚠️ POTENTIAL BUG DETECTED
        Test: [test_name]
        Expected: [expected_behavior]
        Actual: [actual_behavior]
        Source logic suggests: [analysis]
        
        Please review the source code in [file:line] to confirm intended behavior.
        ```
      - **Ambiguous**: Ask user to clarify intended behavior
   
   e. **Document your reasoning**:
      - Explain why you made each change
      - Note any assumptions about intended behavior
      - Flag any edge cases discovered

4. **Re-run tests after each fix**:
   - Run tests incrementally to verify fixes
   - Ensure no regressions are introduced

5. **Review CLAUDE.md** if you encounter testing patterns questions

6. **Final verification**:
   - Run all tests for the module
   - Confirm all tests pass

## Analysis Guidelines

When analyzing test failures:

### Signs the TEST is wrong:
- Test logic doesn't align with method docstring
- Test makes incorrect assumptions about input validation
- Test has obvious typos or copy-paste errors

### Signs the SOURCE has a bug:
- Source code logic clearly produces different result than documented
- Mathematical calculations are incorrect
- Boundary conditions aren't handled as specified
- Missing validation that documentation claims exists
- Inconsistent behavior with similar methods in codebase

### When to STOP and alert:
- Source behavior differs from documentation
- Test comment says "should do X" but source does Y
- Mathematical/logical error in source is apparent
- You're tempted to change a test just to make it pass
- The test failure reveals an edge case not handled in source

## Error Handling

If you encounter:
- **Import errors**: Check fixture dependencies and module paths
- **Fixture errors**: Verify fixtures exist and are properly imported
- **Unexpected exceptions**: Analyze if source should handle these cases
- **Flaky tests**: Note timing issues or external dependencies

## Quality Checklist

Before considering debugging complete:
- [ ] All tests pass
- [ ] No test expectations were changed without confirming source is correct
- [ ] All potential bugs have been flagged to user
- [ ] Changes are documented with reasoning
- [ ] Edge cases discovered are documented