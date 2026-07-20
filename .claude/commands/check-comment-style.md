---
description: Check and fix block comments in specified files against the writing style guidelines
argument-hint: <file paths>
---

# Check Comment Style

Check and fix block comments in `$ARGUMENTS` against `docs/WRITING_STYLE.md`. This command runs three focused agent passes over each file, one per rule group, applying fixes and flagging uncertain cases with `# FIXME:` markers for manual review.

## Scope

This command checks comment prose style only. It does not check:

- Line wrapping (handled by octowrap)
- Non-ASCII character detection (handled by the ascii-only pre-commit hook, which is read-only; this command does fix doubled hyphens and restructures em dash sentences, since doubled hyphens are pure ASCII and invisible to the hook)
- Docstrings (handled by docformatter and manual review)

## Design Constraints

- Each pass edits the file, so the three passes for a given file must run strictly in sequence: a pass must not start until the previous pass's agent has finished. Different files are independent, so their pass chains may run in parallel.
- Each agent receives only the verbatim text of its own rule sections, never the full style doc. An agent cannot misapply rules it has never seen, and a short closed checklist is more reliable than a long filtered one.
- Never hard-code style rules in this command or in agent prompts. All rules come from reading `docs/WRITING_STYLE.md` at runtime and copying its sections verbatim. If the style doc changes, the agents' behavior changes automatically.
- Use `model: "sonnet"` for all agents to balance cost and quality.
- Agents never see this command file. Anything an agent must know has to appear in its prompt, so include the Common Agent Instructions below verbatim in every agent prompt.

## Steps

1. If `$ARGUMENTS` is empty, stop and ask the user which files to check. Otherwise, confirm each file exists and contains block comments by running `grep -c '^\s*#' <file>` on each. Do not read the files yourself; the agents will read them. Skip (and later report) any file with no block comments.

2. Read `docs/WRITING_STYLE.md` in full. Extract the following sections verbatim (heading through last bullet, no summarizing or paraphrasing) for use in agent prompts:
   - The "Sentence Structure" section and the "How to Handle Each Forbidden Case" subsection (for pass 1)
   - The "Terminology" section (for pass 2)
   - The "Math and Numbers" section (for pass 3)

   The "File Formatting" and remaining "ASCII Only" material is intentionally given to no pass: it is out of scope for block comments or covered by hooks.

3. Build the class name inventory for pass 2. Run:

   ```bash
   grep -rn '^class ' pterasoftware/ --include='*.py' | sed 's/.*class \([A-Za-z0-9_]*\).*/\1/' | sort -u
   ```

   This list tells the pass 2 agent which words in comments should use class-name capitalization (e.g., "Wing" not "wing" when referring to the code object). Always build it fresh.

4. For each file, run the three passes in this order. Sentence restructuring changes wording that the later passes inspect, so it goes first. Math and number formatting is nearly orthogonal to the other two, so it goes last.

   Pass 1, sentence structure: the agent receives the verbatim "Sentence Structure" section, the verbatim "How to Handle Each Forbidden Case" subsection, the file path, and the Common Agent Instructions.

   Pass 2, terminology and class references: the agent receives the verbatim "Terminology" section, the class name inventory from step 3, the file path, and the Common Agent Instructions, plus these two pass-specific instructions:
   - Use the class name inventory to decide whether a word in a comment refers to a code object (capitalize, no spaces) or an abstraction (lowercase, spaced words). If you are unsure which one a word refers to, leave it unchanged and insert `# FIXME: Manually review this comment.` on its own line above it.
   - Do not attempt to fix coordinate system naming (the conventions referenced in the CRITICAL bullet); those conventions live in documents you have not been given. If a comment looks like it violates them, insert a FIXME marker instead of editing.

   Pass 3, math and number formatting: the agent receives the verbatim "Math and Numbers" section, the file path, and the Common Agent Instructions.

5. Wait for every file's chain to complete and collect all agent reports.

6. Reformat. Run the modified files through pre-commit to handle formatting side effects (octowrap rewrapping in particular):

   ```bash
   pre-commit run --files <each modified file>
   ```

   This runs before the summary so that reported line numbers are not stale.

7. Summarize the results. Present a single consolidated report:
   - Changes made: list each change with file, rule, and reasoning (line numbers may have shifted during reformatting, so quote the comment text rather than relying on line numbers alone).
   - FIXME markers left: list each `# FIXME:` marker with file and why the agent was unsure.
   - Files with no issues: list any files that passed all three passes cleanly.
   - Files skipped: list any files skipped in step 1 and why.

## Common Agent Instructions

Include these verbatim in every agent prompt:

- Read the file before editing.
- Edit block comments only: lines whose first non-whitespace character is `#`. Never edit docstrings, string literals, or code.
- Do not change line wrapping. Octowrap handles that separately.
- Only fix comments with clear, obvious violations of the rules you were given. Do not go on a fishing expedition.
- If a comment looks structurally broken (mid-sentence truncation, dangling preposition, sentence that does not parse) but you are not confident in the fix, insert `# FIXME: Manually review this comment.` on its own line above it and leave the comment unchanged. Do not insert a marker if one is already present above that comment.
- In your final report, list each change with the comment's original text, the rule violated, and your reasoning. Also list each FIXME marker you inserted and why you were unsure.
