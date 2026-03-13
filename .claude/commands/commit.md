---
description: Analyze the current staged changes, generate a commit message following proper line width guidelines, and create the commit
---

# Create a Commit

Generate a commit message for the current staged changes, then create the commit.

## Steps

1. **Gather context** by running these git commands in parallel:
   - `git status` to see all changed and untracked files
   - `git diff --staged` to see staged changes
   - Do not stage additional files.
   - Stop and warn the user if there are no staged changes.
   - Stop and warn the user if you see staged changes that likely contain secrets (e.g., `.env`, credentials, etc.).
2. **Draft the commit message** following these rules:

   ### Subject line

   - Use imperative mood (e.g., "Add feature" not "Added feature").
   - Keep it to 50 characters or fewer. This is a hard limit.
   - Do not end with a period.
   - Capitalize the first letter of the first word.

   ### Body (only if needed for non-trivial changes)

   - Separate from the subject with one blank line.
   - Wrap every line at 72 characters. This is a hard limit.
   - Explain *why* the change was made, not *what* was changed (the diff shows that).
   - Use complete sentences ending with periods, each of which written using an imperative mood.
   
   Do not include a footer or co-authorship line.
3. **Validate the subject line** by checking its length with `awk`:
   ```bash
   awk 'BEGIN { s = "Your subject line here"; if (length(s) > 50) print "FAIL: subject is " length(s) " chars (max 50)"; else print "OK: " length(s) " chars" }'
   ```
   If it exceeds 50 characters, shorten it and re-check before proceeding.
4. **Wrap the body** (if present) using Python's `textwrap` module:
   ```bash
   python -c "import textwrap; print(textwrap.fill('''Your body text here.''', width=72))"
   ```
   Use the wrapped output as the final body text.
5. **Present the draft** to the user for approval before committing.
6. **Create the commit** using a heredoc to preserve formatting:
   ```bash
   git commit -m "$(cat <<'EOF'
   Subject line here (50 characters or less)

   Optional body here (wrapped at 72 characters per line).
   EOF
   )"
   ```
7. **Verify** by running `git status` after the commit completes.

## Important Reminders

- Never use `--no-verify` or skip pre-commit hooks.
- If the pre-commit hook fails, fix the issue and create a NEW commit (do not amend).
- Never push unless the user explicitly asks.
- If there are no changes to commit, inform the user and stop.
