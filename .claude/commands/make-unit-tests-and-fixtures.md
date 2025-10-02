---
description: Generate comprehensive unit tests and fixtures for specified module
---

# Create Unit Tests and Fixtures

Generate comprehensive unit tests and fixtures for the contents of $ARGUMENTS.

## Progress Tracking

Follow these steps carefully and track your progress:

- [ ] Read $ARGUMENTS carefully, noting the public classes and functions it contains
- [ ] Understand the structure of each public class and function
- [ ] Search the codebase for other relevant files
- [ ] Read and understand any relevant files found
- [ ] Read test pattern examples (tests/unit/test_wing_cross_section.py, tests/unit/test_wing_cross_section_movement.py)
- [ ] Read fixture pattern examples (tests/unit/fixtures/geometry_fixtures.py, tests/unit/fixtures/wing_cross_section_movement_fixtures.py)
- [ ] Understand how geometry/wing_cross_section.py and movements/wing_cross_section_movement.py inform design of example tests and fixtures
- [ ] Review **all** of CLAUDE.md for guidance
- [ ] Search existing tests/fixtures to avoid redundancy
- [ ] Create comprehensive tests and fixtures for each public class/function in $ARGUMENTS
- [ ] Review the tests and fixtures for each public class/function in $ARGUMENTS, both those you created any pre-existing ones
- [ ] Explain how you checked if the tests and fixtures are non-redundant
- [ ] Explain how you checked if the tests and fixtures follow your guidelines
- [ ] Edit or delete any tests and fixtures relevant to $ARGUMENTS that are redundant or don't follow your guidelines
- [ ] Update tests/unit/__init__.py and tests/unit/fixtures/__init__.py with any new files or edited docstrings

## Detailed Steps

1. **Read $ARGUMENTS carefully**, noting the classes and functions it contains.

2. **Understand the structure** of each class and function.

3. **Search the codebase** for other relevant files.

4. **If you find any relevant files**, read and understand those as well.

5. **Study existing test patterns**:
   - Read tests/unit/test_wing_cross_section.py
   - Read tests/unit/test_wing_cross_section_movement.py
   - Read tests/unit/fixtures/geometry_fixtures.py
   - Read tests/unit/fixtures/wing_cross_section_movement_fixtures.py
   - Understand the patterns for creating fixtures and unit tests

6. **Study the implementation**:
   - Read geometry/wing_cross_section.py
   - Read movements/wing_cross_section_movement.py
   - Consider how the implementation reveals test design patterns

7. **Re-read through CLAUDE.md**, ensuring you understand all its guidance. Don't just read the first 100 lines. Read the entire thing!

8. **Search the testing and fixture files** to understand what has already been implemented so you don't write unnecessary fixtures and tests.

9. **For each class and function within $ARGUMENTS**:
   
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
   
   d. **Create the tests** (but don't run them - we will run and debug later):
   
   e. **Review and polish**:
      - Re-read the fixtures and tests you just wrote (and any previously existing ones for this class/function)
      - Check if they match the guidance in CLAUDE.md
      - Verify they match the unit testing patterns
      - Ensure they aren't redundant with existing tests
      - Explain what you checked and your findings
      - Polish any identified areas of improvement, including potentially removing redundant tests/fixtures.
   
   f. **Repeat these sub-steps** for all remaining classes and functions within $ARGUMENTS

10. **Read through tests/unit/__init__.py and tests/unit/fixtures/__init__.py** and update them with any new files or any changes to the docstrings of existing modules. Also add to or update the import list. Make sure both the import lists and the docstrings list modules in alphabetical order!

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
- [ ] Documentation is clear
- [ ] Python files end with a blank newline