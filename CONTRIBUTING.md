# Contributing Guidelines

We are excited to welcome new users, so thank you for considering contributing to Ptera
Software! Hopefully this help you get started. If you have any questions not answered
here, open an issue request. We'll respond and to help out!

## Before Contributing

Before contributing, you should check out Ptera Software's other reference documents (
listed in recommended reading order :wink:):

* [README](README.md)
* [Code of Conduct](CODE_OF_CONDUCT.md)
* [Security Policy](SECURITY.md)
* [License](LICENSE.txt)

## Contributing

There are three pathways to contributing: [reporting a bug](#reporting-a-bug),
[requesting an enhancement](#requesting-an-enhancement),
and [creating an enhancement](#creating-an-enhancement).

### Reporting a Bug

**If the bug you've found is a security vulnerability, please do not post it as an
issue. Instead, follow the guidelines in our [security policy](SECURITY.md).**

If you discover a bug in Ptera Software, head over to our
[issues page](https://github.com/camUrban/PteraSoftware/issues). Then, use to search bar
to check that it doesn't already have an open issue. If there is an open issue, feel
free to comment on this issue or use one of the emoji reactions. This will let us know
the problem is being experienced by multiple users, and make us move extra quickly to
fixing it!

If no one else has reported this bug, do so by creating a new issue! Please use our
[bug report template](.github/ISSUE_TEMPLATE/bug_report.md) so that we get all the
information we need to help you ASAP. Also, add the `bug` and `help wanted ` labels (and
any other applicable labels) to the issue.

### Requesting an Enhancement

As with bugs, the best way to request a feature or to suggest an enhancement is to open
issue. Again, make sure you search
the [issues page](https://github.com/camUrban/PteraSoftware/issues) to see if anyone
else has requested the same thing before creating a new issue. If there is an open
issue, feel free to comment on this issue or use one of the emoji reactions. This will
let us know that multiple users would benefit from this change, and we'll try extra hard
to implement it quickly!

If no one else has requested this change, do so by creating a new issue! Please use our
[feature request template](.github/ISSUE_TEMPLATE/feature_request.md) so that we get all
the information we need to help you ASAP. Also, add the `enhancement` and `help wanted`
labels (and any other applicable labels) to the issue.

### Creating an Enhancement

If you can't think of a particular enhancement you'd benefit from, search for issues
labeled `help wanted` or
`good first issue` for inspiration.

If you have your own idea for an enhancement, we recommend searching the
[issues page](https://github.com/camUrban/PteraSoftware/issues) to make sure there is no
other active development on this enhancement. If you find someone else working on the
same feature or change, add a comment letting them know you'd like to help out!

If you do not see anything on the issues page, we recommend opening an enhancement issue
as its possible the feature you're about to implement already exists, and we'll comment
letting you know! Also, be sure to indicate in your post that you'll take charge of
implementing this change (otherwise, someone else might start working on your request).

Next, we recommend you clone the Ptera Software repository, and then create new feature
branch from your develop branch. You should create a different feature branch for each
enhancement you are working on. Follow the
[GitFlow methodology](https://nvie.com/posts/a-successful-git-branching-model/?fbclid=IwAR3F9IwEXG1T6oMn5Bnk84_u0mv_RAEI5qTJQE7Puovj0hbGZcA8ly_KXYI)
for creating and managing feature branches.

In this feature branch, implement your enhancement. Be sure to commit frequently and
write detailed commit messages. Additionally, please try to use the follow style guide:

* Always use S.I. units for all calculations and results.
* Format your code using the [black](https://github.com/psf/black).
* Use a line length of 88 (which is black's default).
* All modules, scripts, classes, methods, and functions should have 
[docstrings](https://realpython.com/documenting-python-code/#docstring-types).
* Use the reStructuredText docstring format.
* Your code itself should be well documented using block comments.
* Inline comments should be avoided if possible.
* Use tagging (`TODO`, `BUG`, etc.) to mark areas of the code that need changing.
* Major new functionality should also include new tests, which you should add to the 
tests package.

Once you've finished implementing the enhancement, push your feature branch to your
clone of Ptera Software. Then, make a pull request to merge your feature branch with 
the original Ptera Software's develop branch using our 
[pull request template](.github/pull_request_template.md).
