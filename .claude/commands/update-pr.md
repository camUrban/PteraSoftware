---
description: Refresh an existing pull request's title, body, and labels to reflect new commits, following the same template and conventions as create-pr
---

# Update a Pull Request

Refresh the open pull request for the current branch so its title, body, and labels reflect the commits added (or undone) since it was last written, then apply the changes with `gh pr edit`. This is an incremental refresh: it preserves human-authored content and changes only what the new or removed work requires.

This command takes no arguments. It operates on the single open PR whose head is the current branch, and it leaves that PR's base branch and draft status unchanged.

## Shared conventions

Read `.claude/commands/create-pr.md` and apply every convention it defines: the `.github/pull_request_template.md` section structure and per-section detail level, the body writing conventions (no hard-wrapping, `docs/WRITING_STYLE.md`, backticking every identifier, path, and inline code span), the title rules (branch-derived `[FEATURE]`/`[BUG]` prefix and the 42-character limit with its awk check), label selection from `.github/labels.yml`, the checklist rule (never check an item that asserts a GitHub action or the ReadTheDocs build passes), and the environment constraints (an already-pushed branch; only the `gh pr` porcelain and never `gh api`; an ASCII-only title and body with no footer beyond the AI-use policy's `Assisted-by:` disclosure line; and never attaching a milestone or project or requesting reviewers). Only the differences specific to updating an existing PR are spelled out below.

## Environment constraints

- `git push` is denied, so this command cannot push. The new commits must already be on the remote; if the branch has unpushed commits, stop and ask the user to push, then re-run.
- `gh api` is denied. Use only the `gh pr` porcelain (`gh pr view`, `gh pr edit`).

## Steps

1. **Locate the PR and gather context** with these read-only commands (all run without permission prompts):
    - `git status -sb` to confirm the branch is pushed and current, using the same first-line interpretation as create-pr. If the branch has no upstream, a `[gone]` upstream, or any `[ahead N]` marker, stop and tell the user to push first so the PR reflects the new commits.
    - `gh pr list --head <branch> --state open` to find the open PR. If there is none, stop and tell the user to run `/create-pr` first. If there is more than one, stop and report them rather than guessing.
    - `gh pr view <number> --json number,url,title,body,labels,baseRefName,isDraft` to read the PR's current state. Use `baseRefName` as the base for all diffs; do not assume `main`.
    - `git log --oneline <base>..HEAD`, `git diff --stat <base>...HEAD`, and `git diff <base>...HEAD` to see the full, current change set the PR should describe.
    - Read `.github/pull_request_template.md` and `.github/labels.yml` so the refreshed body matches the current template and the labels come from the canonical set.
2. **Diff the PR against reality.** Compare the existing body's `Description` and `Changes` against the current diff to identify (a) new changes not yet described, (b) described changes that have since been undone or reverted and no longer appear in the diff, and (c) prose that is now inaccurate or has formatting or template defects.
3. **Refresh the body incrementally**, preserving human-authored content and changing only what is required:
    - Keep the author's `Motivation`, `Relevant Issues`, manually written `Description` prose, and the current checklist marks, except where the edits below require a change.
    - Add `Changes` bullets for new work, and remove bullets describing work that has been undone. Extend the `Description` to cover genuinely new scope, and trim it where scope was removed.
    - Update `Dependency Updates` if dependencies were added, changed, or removed, and re-evaluate `Change Magnitude` against the new total change set, changing the selected line only if the magnitude genuinely changed.
    - Re-evaluate the checklist against the new state: preserve existing marks, add any template items that are missing, uncheck any item the new work makes no longer satisfied, reconsider items that were not applicable but now are (for example, package type hints once package code is touched), and uncheck any item asserting a GitHub action or the ReadTheDocs build passes.
    - Correct typos and formatting defects anywhere in the body, including removing hard wraps, fixing backticking, and repairing any drift from the template's structure, even in otherwise-preserved sections.
    - Apply the AI-use policy's disclosure (docs/AI_USE_POLICY.md): keep any existing `Assisted-by:` line as the body's last line. When none is present, add one of the form `Assisted-by: MODEL_OR_TOOL_NAME` (preceded by a single blank line) if the conversation context makes clear the user used AI on the PR, or if the requested updates to the description go beyond formatting and typo corrections. Write the tag's casing exactly as shown, include no email address, and use your own model name unless the context names a different assisting tool.
4. **Re-evaluate the title and labels** against the new change set:
    - Recompute the title prefix from the branch name and confirm the title still fits and still describes the change set. Change it only if the new work makes it inaccurate or it breaks the rules, and re-run the 42-character and ASCII awk check from create-pr on any new title.
    - Recompute the applicable labels from `.github/labels.yml`. Plan to add newly applicable labels and remove labels that no longer apply, keeping them consistent with any `[FEATURE]`/`[BUG]` title prefix.
5. **Present the planned update** in your reply: the PR number and URL, the old and new title (or "unchanged"), the labels to add and remove (or "unchanged"), and the full new body. This is the user's review opportunity. If nothing needs to change, say so and stop without editing.
6. **Apply the update** with `gh pr edit`. Provide the body on stdin via a quoted heredoc so the permission prompt shows the exact new body as the final review gate:
   ```bash
   gh pr edit <number> \
     --title "Updated title here" \
     --add-label "new-label" \
     --remove-label "stale-label" \
     --body-file - <<'EOF'
   # Description

   ...full refreshed body...
   EOF
   ```
    - Include `--title` only if the title changed. Repeat `--add-label` and `--remove-label` once per label, and omit them when labels are unchanged. Do not pass `--milestone`, `--project`, `--reviewer`, `--base`, or any draft flag; this command leaves the base branch and draft status as they are.
    - The `gh pr edit` permission prompt is the final gate. If the user denies it, treat any feedback as revision input, update the plan, and repeat from step 5. If they deny without feedback, stop and report that the PR was not changed.
7. **Confirm** by reporting the PR's URL.

## Important Reminders

- This command edits an existing PR; if none exists for the branch, stop and direct the user to `/create-pr`.
- Preserve human-authored content by default; change prose only to add new scope, remove undone scope, or fix typos, formatting (including hard wraps), and template drift.
- Never push; if the new commits are not on the remote, stop and ask the user to push.
- Never use `gh api`; use only the `gh pr` porcelain.
- Apply all create-pr conventions to the refreshed title, body, and labels: ASCII-only with no hard-wrapping, branch-based prefix and the 42-character title limit, backticked identifiers, no checked action or ReadTheDocs checklist items, the policy's `Assisted-by:` line as the only permitted footer, and no milestone, project, or reviewer changes.
- Leave the PR's base branch and draft status unchanged.
