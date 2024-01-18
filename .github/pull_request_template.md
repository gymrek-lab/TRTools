## Checklist

* [ ] I've checked to ensure there aren't already other open [pull requests](../../../pulls) for the same update/change
* [ ] I've prefixed the title of my PR according to [the conventional commits specification](https://www.conventionalcommits.org/). If your PR fixes a bug, please prefix the PR with `fix: `. Otherwise, if it introduces a new feature, please prefix it with `feat: `. If it introduces a breaking change, please add an exclamation before the colon, like `feat!: `. If the scope of the PR changes because of a revision to it, please update the PR title, since the title will be used in our CHANGELOG.
* [ ] At the top of the PR, I've [listed any open issues that this PR will resolve](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue#linking-a-pull-request-to-an-issue-using-a-keyword). For example, "resolves #0" if this PR resolves issue #0
- [ ] I've explained my changes in a manner that will make it possible for both users and maintainers of TRTools to understand them
* [ ] I've added tests for any new functionality. Or, if this PR fixes a bug, I've added test(s) that replicate it
* [ ] All directories with large test files are listed in [the "exclude" section](https://python-poetry.org/docs/pyproject/#include-and-exclude) of our pyproject.toml so that they do not appear in our PyPI distribution
* [ ] I've updated the relevant REAMDEs with any new usage information and checked that the newly built documentation is formatted properly
* [ ] All functions, modules, classes etc. still conform to [numpy docstring standards](https://numpydoc.readthedocs.io/en/latest/format.html)
* [ ] (if applicable) I've updated the pyproject.toml file with any changes I've made to TRTools's dependencies, and I've run `poetry lock --no-update` to ensure the lock file stays up to date and that our dependencies are locked to their minimum versions
* [ ] In the body of this PR, I've included a short address to the reviewer highlighting one or two items that might deserve their focus
