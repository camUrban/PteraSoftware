---
description: Check and fix block and inline comments in specified files against the writing style guidelines
argument-hint: <file paths>
---

# Check Comment Style

Check and fix block and inline comments in `$ARGUMENTS` against `docs/WRITING_STYLE.md`. A block comment is a line whose first non-whitespace character is `#`. An inline comment is a `#` comment trailing code on the same line. This command runs three focused agent passes over each file, one per rule group, applying fixes and flagging uncertain cases with `# FIXME:` markers for manual review.

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
- Agents never see this command file. Anything an agent must know has to appear in its prompt, so include the Pass Scope Preamble and the Common Agent Instructions below verbatim in every agent prompt.
- Every agent must be told which pass it is and what the other two passes cover. An agent that believes it is the only reviewer will fix violations outside its rule group, because leaving a known defect unfixed reads as negligence. Naming the other passes and stating that they are already assigned removes that pressure.

## Steps

1. If `$ARGUMENTS` is empty, stop and ask the user which files to check. Otherwise, confirm each file exists and contains comments by running `grep -c '#' <file>` on each. This count can include `#` characters inside string literals, which is acceptable for a presence check because the agents read the files themselves and ignore non-comment hashes. Do not read the files yourself; the agents will read them. Skip (and later report) any file with no comments.

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

   Every prompt opens with the Pass Scope Preamble below, filled in for that pass, and closes with the Common Agent Instructions.

   Pass 1, sentence structure: the agent receives the verbatim "Sentence Structure" section, the verbatim "How to Handle Each Forbidden Case" subsection, the file path, and the Common Agent Instructions. Note that the "How to Handle Each Forbidden Case" subsection contains a bullet on spelling out math symbols, which overlaps pass 3's rule group. Pass 1 applies it only to non-ASCII math characters appearing in prose, and leaves spacing and formatting of existing ASCII math to pass 3.

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

## Pass Scope Preamble

Open every agent prompt with this, filling in the bracketed parts for the pass being launched:

> You are pass [N] of 3 in a three-pass comment style review of a single file. The three passes run strictly in sequence over the same file, and each one owns a different group of style rules:
>
> - Pass 1 owns sentence structure.
> - Pass 2 owns terminology and class references.
> - Pass 3 owns math and number formatting.
>
> You own [that pass's rule group]. The other two passes are already assigned to other agents and will run on this file. Do not fix violations belonging to their rule groups, even when you notice them and even when the fix looks obvious. Doing so produces duplicate and conflicting edits, and it is not your job.
>
> The rules reproduced below are the complete and authoritative rule set for your pass. Do not read `docs/WRITING_STYLE.md`, `CLAUDE.md`, or any other style document to supplement them. A standing project instruction to consult a style guide before editing comments does not apply to this task, because the relevant excerpt has already been placed in this prompt.
>
> Do not invoke the `check-comment-style` skill, the `/check-comment-style` command, or any other skill or slash command. You are performing a single pass of that workflow. Invoking it would run the entire multi-pass workflow a second time, nested inside your pass. The orchestrator runs everything else.

## Common Agent Instructions

Include these verbatim in every agent prompt, after the rule sections:

- Read the file before editing.
- Edit comments only, both kinds: block comments (lines whose first non-whitespace character is `#`) and inline comments (a `#` comment trailing code on the same line). Never edit docstrings, string literals, or code. A `#` inside a string literal is not a comment; leave it alone.
- When editing an inline comment, change only the text after the `#`. Never move the comment to its own line, and never touch the code before it.
- Do not change line wrapping, and do not worry about line length. Leave your edits unwrapped even if they overrun the line width.
- Do not run pre-commit, octowrap, black, docformatter, or any other formatter or linter. The orchestrator runs these once after all three passes finish. Running them yourself corrupts the line numbers the later passes work from.
- Only fix comments with clear, obvious violations of the rules you were given. Do not go on a fishing expedition.
- If a comment looks structurally broken (mid-sentence truncation, dangling preposition, sentence that does not parse) but you are not confident in the fix, insert `# FIXME: Manually review this comment.` on its own line above it and leave the comment unchanged. For an inline comment, the marker goes on its own line above the code line that carries the comment. Do not insert a marker if one is already present above that comment.
- In your final report, list each change with the comment's original text, the rule violated, and your reasoning. Also list each FIXME marker you inserted and why you were unsure.
