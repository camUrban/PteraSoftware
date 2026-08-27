---
description: Draft a GitHub issue from the conversation following this repository's streamlined issue conventions, then open it with the GitHub CLI
---

# Open an Issue

Draft a GitHub issue for this repository from the current conversation and the command's argument, then open it with `gh issue create`. The issue follows this repository's lived, streamlined issue conventions, which are a maintainer-oriented adaptation of the templates in `.github/ISSUE_TEMPLATE/`, not the verbose template forms.

Arguments: $ARGUMENTS

The argument is a short description of what the issue is about. If its first token is the keyword `gfi`, mark the issue as a good first issue (add the `good_first_issue` label and a "Hints for New Contributors" section), and treat the rest of the argument as the description. The command also draws on the current conversation for context.

## Environment constraints

- `gh api` is denied. Use only the `gh issue` porcelain (`gh issue list`, `gh issue view`, `gh issue create`).
- The title and body must contain only printable ASCII characters. Do not add a footer or any "Generated with" or "Co-Authored-By" line. The only permitted trailing line is the AI-use policy's `Assisted-by:` disclosure described in step 4.
- No git state is needed: an issue is not tied to a branch or commit, so do not run push or diff commands and do not require a pushed branch.

## Steps

1. **Determine the issue type and gather its substance** from the conversation and the argument:
    - Decide whether this is a bug or a non-bug issue (feature, maintenance, performance, or a question). A bug reports incorrect behavior; everything else is a problem statement or proposal.
    - Gather the concrete substance: what the problem is, where it lives (specific files, functions, or classes), and the proposed fix or next step. Draw this from the conversation and the argument, supplementing only with light, read-only lookups (for example, confirming a file path or function name) when needed to make the issue precise.
    - If the available context is too thin to write a specific, accurate issue, stop and ask the user for the missing details (for a bug, the symptom and where it occurs; for a feature or task, the problem and the proposed approach). Never fabricate an issue.
2. **Check for duplicates** with the pre-allowed read-only commands (no permission prompts):
    - `gh issue list --search "<keywords>" --state all --limit 20 --json number,title,state,url` to find open or closed issues on the same topic.
    - If any look like plausible duplicates, list them (number, state, title, and URL) in your draft and ask the user whether to proceed before creating anything.
3. **Choose the title** following the repository's issue conventions:
    - Sentence case, no trailing period. The title may contain backticked identifiers, and there is no length limit for issue titles.
    - Set the prefix from the issue type and label: if it warrants the `bug` label, prefix with `[BUG] `; if it warrants the `feature` label, prefix with `[FEATURE] `; otherwise use no prefix. Issues have no branch, so the prefix is label-driven, unlike pull request titles.
4. **Draft the body** following the lived, streamlined convention. Do not copy the verbose `.github/ISSUE_TEMPLATE/` placeholder prose, the end-user-reporter sections (`Reproduction`, `Screenshots`, `Desktop`), or `Alternative Solutions`. Keep each section to a few concise sentences.

   Apply these writing conventions throughout the body:
    - Do not hard-wrap. Write each paragraph and bullet as a single continuous line; GitHub reflows Markdown for display.
    - Follow `docs/WRITING_STYLE.md` for tone, terminology, and mechanics.
    - Backtick every file, path, module, class, function, method, and inline code span.

   Use this structure:
    - For a **bug**:
        - `# Bug Description`: a concise description of the incorrect behavior, followed by a `**Location(s):**` line naming the affected files in backticks.
        - `## Expected Behavior`: what should happen instead.
        - `## Additional Context` (optional): extra notes, including any `TODO`, `TEST`, or `TEMP` comment that references the issue.
    - For a **non-bug** (feature, maintenance, performance, or question):
        - `# Problem Statement`: a concise description of the problem or opportunity, followed by a `**Location(s):**` line naming the affected files in backticks.
        - `## Proposed Solution`: the fix or next step, as a short paragraph or a numbered list.
        - `## Additional Context` (optional): extra notes, including any referencing `TODO`, `TEST`, or `TEMP` comment.
    - If the `gfi` keyword was given, append a `## Hints for New Contributors` section: remind the contributor to read `CONTRIBUTING.md` and set up the development environment, point to the specific files and any relevant docs (for example, `docs/RUNNING_TESTS_AND_TYPE_CHECKS.md` for running tests, or the convention docs when relevant), give a short numbered plan, and remind them to run the tests afterward.
    - Always end the body with the AI-use policy's disclosure (docs/AI_USE_POLICY.md): a final line of the form `Assisted-by: MODEL_OR_TOOL_NAME`, separated from the last section by a single blank line. Running this command is AI drafting assistance by definition, so the line is unconditional. Write the tag's casing exactly as shown, include no email address, and use your own model name.
5. **Choose labels and assignment**:
    - Read `.github/labels.yml`, the canonical label set, and pick the applicable labels from it. Include `bug` for a bug or `feature` for a feature so the label matches any title prefix, and add `maintenance`, `question`, `performance`, and so on as the topic warrants, matching the combinations seen in existing issues (for example, `maintenance` with `question` for an open-ended investigation).
    - If the `gfi` keyword was given, include the `good_first_issue` label.
    - Never assign the issue to anyone; do not pass `--assignee`. Leave it unassigned so it can be picked up, which matters especially for good first issues.
    - Do not attach a milestone or project. Never pass `--milestone` or `--project`.
6. **Present the draft** in your reply: the title, the full body, the selected labels, any possible duplicates found, and the fact that the issue will be unassigned. This is the user's review opportunity.
7. **Open the issue** with a single command, providing the body on stdin via a quoted heredoc so the permission prompt shows the exact body as the final review gate:
   ```bash
   gh issue create \
     --title "Your title here" \
     --label "label-one" --label "label-two" \
     --body-file - <<'EOF'
   # Problem Statement

   ...full body following the convention...
   EOF
   ```
    - Repeat `--label` once per label. Do not pass `--assignee`, `--milestone`, or `--project`.
    - The `gh issue create` permission prompt is the final gate. If the user denies it, treat any feedback as revision input, update the draft, and repeat from step 6. If they deny without feedback, stop and report that no issue was opened.
8. **Confirm** by reporting the new issue's URL, which `gh issue create` prints on success.

## Important Reminders

- Draw the issue's substance from the conversation and the argument; never fabricate. If context is too thin, stop and ask for specifics.
- Follow the lived, streamlined convention, not the verbose `.github/ISSUE_TEMPLATE/` forms: include the `**Location(s):**` line, and omit `Reproduction`, `Screenshots`, `Desktop`, and `Alternative Solutions`.
- Keep the title and body ASCII-only and unwrapped, with backticked identifiers and no footer beyond the policy's `Assisted-by:` disclosure line.
- Never use `gh api`; use only the `gh issue` porcelain.
- Never assign the issue, and never attach a milestone or project.
- Add the `good_first_issue` label and a `## Hints for New Contributors` section only when the `gfi` keyword is given.
