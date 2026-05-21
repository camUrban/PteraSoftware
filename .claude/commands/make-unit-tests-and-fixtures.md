---
description: Generate comprehensive unit tests and fixtures for a specified module
---

# Create Unit Tests and Fixtures

Generate, but do not run, comprehensive unit tests and fixtures for the contents of `$ARGUMENTS`.

## Progress Tracking

Follow these steps carefully and track your progress:
- [ ] Read `$ARGUMENTS` carefully, noting the public classes and functions it contains
- [ ] Understand the structure of each public class and function
- [ ] Search the codebase for other relevant files
- [ ] Read and understand any relevant files found
- [ ] Read `tests/unit/test_wing_cross_section_movement.py` for examples of test patterns
- [ ] Read `tests/unit/fixtures/wing_cross_section_movement_fixtures.py` for examples of fixture patterns
- [ ] Understand how `movements/wing_cross_section_movement.py` informs design of example tests and fixtures
- [ ] Create comprehensive tests and fixtures for each public class/function in `$ARGUMENTS`
- [ ] Review the tests and fixtures to verify they follow style guidelines
- [ ] Update `tests/unit/__init__.py` and `tests/unit/fixtures/__init__.py` with any new files or edited docstrings

## Omissions

Please do NOT:
- Run or debug any tests you create
- Reformat any files with black

These steps are both handled by other slash commands.

## Detailed Steps

1. **Read `$ARGUMENTS` carefully**, noting the classes and functions it contains.
2. **Understand the structure** of each class and function.
3. **Search the codebase** for other relevant files.
4. **If you find any relevant files**, read and understand those as well.
5. **Study existing test and fixture patterns**:
    - Read `tests/unit/test_wing_cross_section_movement.py`
    - Read `tests/unit/fixtures/wing_cross_section_movement_fixtures.py`
    - Understand the patterns for creating fixtures and unit tests
6. **Study the implementation**:
    - Read `movements/wing_cross_section_movement.py`
    - Consider how the implementation reveals test design patterns
7. **For each class and function within `$ARGUMENTS`**:
   a. **Determine file placement**:
      - Consider if you should make new testing and fixture files
      - Or if you should add to currently existing ones
   b. **Create a test plan**:
      - List all methods to test
      - Identify edge cases and error conditions
      - Plan fixtures needed for comprehensive testing
      - Note if these will be added to existing files or written in new files
   c. **Create the fixtures**:
      - Follow existing naming conventions
      - Include both valid and invalid test cases
      - Document fixture purpose and usage
   d. **Create the tests** (but don't run them - we will check redundancy and debug later):
   e. **Review for style compliance**:
      - Re-read the fixtures and tests you just wrote
      - Check if they match the guidance in your docs on code and writing styles
      - Verify they match the unit testing patterns
      - Polish any style issues found
   f. **Repeat these sub-steps** for all remaining classes and functions within `$ARGUMENTS`
8. **Read through `tests/unit/__init__.py` and `tests/unit/fixtures/__init__.py`** and update them with any new files or any changes to the docstrings of existing modules. Also add to or update the import list. Make sure both the import lists and the docstrings list modules in alphabetical order!

## Error Handling

If any step fails or information is missing:
- Explicitly state what's missing
- Ask for clarification before proceeding
- Document any assumptions made
- Note any dependencies or prerequisites

## Quality Checklist

Before finalizing:
- [ ] All public classes and functions have tests
- [ ] Edge cases are covered
- [ ] Error conditions are tested
- [ ] Fixtures are reusable
- [ ] Tests are independent
- [ ] Naming is consistent with project patterns
- [ ] Style guidelines are followed
- [ ] `__init__.py` files are updated
