# AI-Use Policy

## Scope

* This policy applies to commits, issue descriptions, pull request descriptions (PRDs), review comments, and security reports.
* Disclosure is strongly encouraged for commits, issue descriptions, and PRDs. Disclosure should not be included for review comments or security reports.
* This policy covers generative AI and LLM tools only. Deterministic automation using tools such as Dependabot and GitHub Actions is out of scope. The disclosure section notes how their trailers are handled.
* These rules apply to all contributors, including the maintainer(s).
* These rules are effective for commits merged after [7e84ae05c](https://github.com/camUrban/PteraSoftware/commit/7e84ae05ca20e9f5eff75a03b71beef385a49124). Earlier history is unchanged.

## Principles

* AI and LLM tools are permitted.
* A contribution must be worth more to the project than it costs to review. Anything failing that test may be closed regardless of how it was produced.
* The contributor is always the author and is fully accountable.
* Read and review all generated content before requesting review, and be able to answer questions about it.
* PRDs and issue descriptions must reflect your own understanding of the change. Drafting them with a tool is fine as long as you review and correct the result.
* AI must not be the final arbiter on whether a contribution is accepted.
* AI agents must not act in project spaces without human approval and direct oversight. This includes bots that open PRs or post comments without a human reviewing the content first.

## Copyright

* Contributors are responsible for having the right to license their contribution under Ptera Software's MIT license.
* Regenerating copyrighted material through a tool does not clear its copyright.
* Violating contributions are removed like any other.

## Disclosure

* Disclosure is strongly encouraged but not required. Non-disclosure is undetectable, so a mandate would only penalize the honest.
* The purpose of disclosure is to help reviewers calibrate, not to track which code is generated.
* Commits:
  * Location: the commit message. Don't redundantly add disclosure to code comments or in committed documentation.
  * Form: `Co-authored-by: MODEL_OR_TOOL_NAME <MODEL_OR_TOOL_EMAIL>` (e.g., `Co-authored-by: Claude Fable 5 <noreply@anthropic.com>`).
  * `Assisted-by:` should not be used. Other disclosure trailers may be normalized on merge.
  * `Signed-off-by:` should not be used. Ptera Software has no DCO, so the trailer carries no meaning here. Tooling that emits it unconditionally, such as Dependabot, is exempt.
* PRDs and issue descriptions:
  * Location: as the last line of the text, separated from the body above by a single blank line.
  * Form: `Assisted-by: MODEL_OR_TOOL_NAME`.
  * Other disclosure messages may be removed or normalized by the maintainer(s).
* Disclosure isn't necessary for routine grammar or formatting assistance.
* The accepted forms for disclosure declaration tags (the string preceding the colon) must be written exactly as shown. Rejected forms are rejected in any casing.

## Good first issues

* Be deliberate and conservative with AI use here. These exist as learning opportunities, and automating them squanders that while adding little.

## Testing

* Existing unit and integration test requirements apply unchanged.
* Be able to explain what behavior each new test actually pins down.

## Amendments

Policy changes are recorded below with each change's commit SHA and rationale.
