---
description: Draft a pull request that strictly follows the repository template, then open it with the GitHub CLI
---

# Create a Pull Request

Draft a pull request for the current branch whose body strictly follows `.github/pull_request_template.md`, then open it with `gh pr create`.

This command takes two independent, optional inputs through its argument, in any order. They are unrelated: setting one has no effect on the other.

- A base branch name, which is any token other than `draft`. If omitted, the base defaults to `main`.
- The literal keyword `draft`, which opens the PR as a draft. If omitted, the PR is opened ready for review.

Examples: `/create-pr` (base `main`, ready), `/create-pr develop` (base `develop`, ready), `/create-pr draft` (base `main`, draft), and `/create-pr develop draft` (base `develop`, draft).

## Environment constraints

- `git push` is denied in this environment, so this command cannot push the branch. The branch must already exist on the remote and be up to date. If it is not, stop and ask the user to push, then re-run.
- `gh api` is denied. Use only the `gh pr` porcelain subcommands. Do not fall back to the REST or GraphQL API.
- The title and body must contain only printable ASCII characters. Do not add a footer or any "Generated with" or "Co-Authored-By" line. The only permitted trailing line is the AI-use policy's `Assisted-by:` disclosure described in step 3.

## Steps

1. **Gather context** by running these read-only commands (all run without permission prompts):
    - `git status -sb` to read the current branch and its upstream state from the first line:
        - `## <branch>` with no `...<remote>` segment means the branch has no upstream. Stop: tell the user to push the branch first.
        - `...<remote> [gone]` means the upstream was deleted. Stop: tell the user to push the branch first.
        - `...<remote> [ahead N]` or `[ahead N, behind M]` means there are unpushed local commits. Stop: tell the user to push first so the PR includes them.
        - `...<remote>` with no `ahead`/`gone` marker means the branch is pushed and current. Proceed.
        - If the branch is the base branch itself (e.g., `main`), stop: a PR cannot be opened from the base onto itself.
    - `gh pr list --head <branch> --state open` to check for an existing open PR for this branch. If one exists, stop and report its number and URL rather than opening a duplicate.
    - `git log --oneline <base>..HEAD` to see the commits the PR will contain. If this is empty, stop: there is nothing to propose.
    - `git diff --stat <base>...HEAD` and `git diff <base>...HEAD` to understand the actual changes. Base the description and the `Changes` bullets on this diff, not on the commit messages alone.
    - Read `.github/pull_request_template.md` so the body reproduces its exact, current section structure.
2. **Choose and verify the title** following the template's title rules and this repository's conventions:
    - Sentence case, no trailing period.
    - Set the prefix from the current branch name: if it begins with `feature/`, prefix the title with `[FEATURE] `; if it begins with `bug/`, prefix it with `[BUG] `; otherwise use no prefix. Each prefix applies if and only if the branch carries the matching prefix, and these are the only two prefixes.
    - A title may state an outcome (e.g., "[FEATURE] Add free flight simulations").
    - The entire title, including any `[FEATURE]` or `[BUG]` prefix, must be at most 42 characters. Verify its length and ASCII-only content with a single awk call (this runs without a permission prompt, so it costs the user nothing):
      ```bash
      title="Your title here"
      awk -v s="$title" 'BEGIN {
        fail = 0
        if (length(s) > 42) { print "FAIL: title is " length(s) " chars (max 42)"; fail = 1 }
        if (s ~ /[^ -~]/) { print "FAIL: title contains non-ASCII characters"; fail = 1 }
        if (fail == 0) print "OK: title is " length(s) " chars and ASCII-only"
      }'
      ```
      Avoid the `!` character anywhere in this call (write `fail == 0`, not `!fail`); the harness escapes `!` in Bash commands, which corrupts the awk program. If the check prints `FAIL`, shorten the title and re-run this check before proceeding.
3. **Draft the body** by reproducing every section heading from the template in order and replacing each placeholder with real content. Do not copy the template's instructional placeholder prose, and do not omit, rename, or reorder sections. Match the level of detail seen in this repository's substantive PRs.

   Apply these writing conventions throughout the body:
    - Do not hard-wrap. Write each paragraph and bullet as a single continuous line; this repository never hard-wraps Markdown and GitHub reflows it for display. This is the opposite of the commit-message convention, which wraps at 72 characters.
    - Consult `docs/WRITING_STYLE.md` and follow its conventions for tone, terminology, and mechanics.
    - Enclose every file, path, directory, module, package, class, object, function, method, attribute, and any other identifier or inline code in backticks (for example, `_transformations.py`, `OperatingPoint`, and `alpha_and_beta_from_vInf_BP1`).

   Fill each section as follows:
    - **Description**: One or two prose paragraphs describing what the change does, scaled to its magnitude. For a small change, a few sentences, noting backward compatibility when relevant (e.g., a new default that leaves existing callers unaffected). For a larger change, a fuller paragraph covering the mechanism and the scope, including anything deliberately left out of scope.
    - **Motivation**: A genuine paragraph explaining why the change is needed and the problem it solves, never a throwaway line. For a bug fix, characterize the symptom, distinguish a numerical or incidental artifact from intended behavior, and state what behavior is left unchanged.
    - **Relevant Issues**: Use GitHub's closing syntax: `Fixes #<n>` for a bug, `Closes #<n>` for a feature. Link only genuinely related issues. If there are none, write `None.`
    - **Changes**: The most detailed section. One bullet per logical change. Name the specific modules, classes, functions, and methods touched (in backticks, per the conventions above), and state precisely what changed in each. Include quantitative results where they apply. Aim for roughly 3 to 5 bullets for a minor change and 8 to 12 for a moderate or major one. Be concrete; avoid vague bullets like "fixed a bug".
    - **Dependency Updates**: List any new or version-changed runtime or dev dependencies with their constraints, derived from changes to `requirements*.txt`, `setup.cfg`, or `pyproject.toml`. If none changed, write `None.`
    - **Change Magnitude**: Choose exactly one of Major, Moderate, or Minor by the nature of the change, not its line count. Major is significant new functionality, a behavior change, or broad impact. Moderate is a medium feature or refactor without large-scale impact (even a large refactor that preserves behavior is Moderate). Minor is a bug fix, small enhancement, or documentation update. Keep only the chosen bolded line and its description; delete the other two lines and the instruction sentence.
    - **Checklist**: Reproduce every item. Mark `[x]` only for items that are genuinely satisfied or not applicable; never blanket-check. Verify each locally checkable item against the actual state (for example, formatting, docstrings, type hints, and, only if you actually ran them, the local test suite). Always leave unchecked (`[ ]`) every item that asserts a GitHub action or the ReadTheDocs build check passes (the `ascii-only`, `black`, `codespell`, `docformatter`, `isort`, and `pre-commit-hooks` actions; the `mypy` action; the `tests` actions; and the ReadTheDocs build check). Those workflows are not triggered until the PR is created, so however certain you are that they will pass, they have not actually run yet and the box stays empty. If the listed set of actions does not match the repository's actual workflows, edit the affected line to match reality, but still leave it unchecked.
   Always end the body with the AI-use policy's disclosure (docs/AI_USE_POLICY.md): a final line of the form `Assisted-by: MODEL_OR_TOOL_NAME`, separated from the last section by a single blank line. Running this command is AI drafting assistance by definition, so the line is unconditional. Write the tag's casing exactly as shown, include no email address, and use your own model name.
4. **Choose labels, assignee, and draft status**:
    - Read `.github/labels.yml`, the canonical label set, and pick the applicable labels from it. If the title carries a `[FEATURE]` or `[BUG]` prefix, include the matching `feature` or `bug` label accordingly.
    - Always assign the PR to the current user with `--assignee "@me"`. The `@me` token resolves to whoever `gh` is authenticated as.
    - Do not attach a milestone or a project, and do not request any reviewers. Never pass `--milestone`, `--project`, or `--reviewer`.
    - Open the PR as a draft only if the argument requested it; otherwise open it ready for review.
5. **Present the draft** in your reply: the chosen base branch, title, full body, selected labels, assignee, and whether it is a draft. This is the user's review opportunity.
6. **Open the PR** with a single command. Provide the body on stdin via a quoted heredoc so its backticks and markdown are preserved literally and the permission prompt shows the user the exact body as a final review gate:
   ```bash
   gh pr create \
     --base main \
     --title "Your title here" \
     --assignee "@me" \
     --label "label-one" --label "label-two" \
     --body-file - <<'EOF'
   # Description

   ...full body following the template...
   EOF
   ```
    - Repeat `--label` once per label. Add `--draft` if requested. The head branch defaults to the current branch, so `--head` is not needed. Do not add `--milestone`, `--project`, or `--reviewer`.
    - The `gh pr create` permission prompt is the final gate. If the user denies it, treat any feedback as revision input, update the draft, and repeat from step 5. If they deny without feedback, stop and report that no PR was opened.
7. **Confirm** by reporting the new PR's URL, which `gh pr create` prints on success.

## Important Reminders

- Never push, and never use `--no-verify` or any flag that bypasses checks. If the branch is not pushed and current, stop and ask the user to push.
- Never use `gh api`; stick to the `gh pr` porcelain.
- Follow `.github/pull_request_template.md` exactly: every section, in order, with real content.
- Set the title prefix from the branch name (`feature/` gives `[FEATURE]`, `bug/` gives `[BUG]`, otherwise none), and keep the whole title, prefix included, to at most 42 characters, verified with the awk check.
- Never check a checklist item that asserts a GitHub action or the ReadTheDocs build passes; those run only after the PR is created, so they have not run yet. Leave every such item `[ ]`, no matter how certain the outcome.
- Keep the title and body ASCII-only, with no footer beyond the policy's `Assisted-by:` disclosure line.
- The body is Markdown: do not hard-wrap it, follow `docs/WRITING_STYLE.md`, and backtick every identifier, path, and inline code span.
- Do not open a PR from the base branch, and do not open a duplicate when one already exists for the branch.
- Do not attach a milestone or project, and do not request reviewers. Always self-assign with `--assignee "@me"`, which targets whoever `gh` is authenticated as.
