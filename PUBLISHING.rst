Publishing
----------

Only maintainers of the trtools repository may publish changes to the package.
If you are a community member and want to contribute new code, see the contributing section in the README.
If you are a community member and have already contributed new code and want us to publish it
now, please contact us (our contact info is in the README)

This document explains how trtools maintainers should publish new changes.
Maintainers should reach consensus before going ahead with publishing changes.
Note that the publishing to PyPI step below will require credentials
that are only made available to core maintainers of TRTools.

We use a simplified version of
`git flow <http://web.archive.org/web/20200520162709/https://nvie.com/posts/a-successful-git-branching-model/>`_
to maintain and publish trtools.
We use the master branch as the default branch with the latest stable codebase.
The builds from this branch are distributed to PyPI and conda.
The develop branch contains new features that have yet to make their way into master.

New Dependencies
----------------
If you've added dependencies to trtools or its tests, those dependencies should be listed in

  * pyproject.toml
  * the .readthedocs_conda_env.yml file in the root of the repository that's used for building
    TRTool's Read The Docs webpage.
  * the appropriate section of the bioconda recipe (see below)


Publishing Steps
----------------

Once changes have been made to develop that are ready to be published, first choose the new version number according to `semantic versioning <https://semver.org>`_. Then set up the environment you're going to publish TRTools from:

#. Create a clean environment.
#. Install setuptools with version >= 40.8.0
#. Additionally, install ``pytest``, ``wheel``, ``build``, and ``twine``
#. Clone the `trtools repo <https://github.com/gymrek-lab/TRTools>`_
#. Check out the develop branch
#. Run :code:`pip install --upgrade pip && pip install -e .`

Then go through the steps of merging the changes into the master branch:

#. Run :code:`pytest` and make sure all the tests pass. Then run :code:`./test/cmdline_tests.sh` and make sure those tests pass.
#. Change the 'Unreleased Changes' section of :code:`RELEASE_NOTES.rst` to the new version number.
#. Check if any changes have been made that have not yet been documented in the release notes. If so, document them.
#. Submit a pull request from develop into master on the github webiste.
#. If the code review checks pass, merge the pull request.
#. Tag the merge commit with the package version in vX.Y.Z format. (For more details on tagging, see `below`)

Then go through the steps of publishing the changed code to PyPI:

1. :code:`cd` into the root of your clone of the trtools repo, checkout master and pull the latest change. Note that the most recent commit *must* be tagged.
2. Run :code:`rm -rf build dist *.egg-info` to make sure all previous build artifacts are removed
3. Run :code:`python -m build` to build the package with the version number you just tagged. (Note: you might need to install ``build`` first.)
5. Run :code:`twine upload dist/*` to upload the distribution to PyPI

Lastly, the change needs to be published to bioconda.

A bioconda bot will automatically open a pull request (within a day?) updating the version number
and the PyPI reference. If there are no new dependencies, no changes to the build,
and no new tests that need to be integrated into the build, and all we need to do is mark that PR as okay.
(Jonathan Margoliash is currently the bioconda recipe maintainer and will get pinged by this. Please notify him to look out for that ping.
If you'd like to be a bioconda recipe maintainer, let's set that up.)

If there are new dependencies or build changes, then we'll need to close the automatic PR without accepting it and make our own.
To do that, go through the following steps to publish the changed code to bioconda: (see `here <http://bioconda.github.io/contributor/workflow.html>`_ for official documentation)

1. Run :code:`cd <project-root> && openssl sha256 dist/trtools-<version>.tar.gz` and save the generated hash code for later
2. Create a fork of the `bioconda recipes <https://github.com/bioconda/bioconda-recipes>`_ repo and clone it.
3. Run ``git remote add upstream https://github.com/bioconda/bioconda-recipes`` so that you can pull from not only your fork but also upstream.
4. Create a new branch from master with your package's name
5. In `the recipe for trtools <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/trtools/meta.yaml#L1-L2>`_, update the version number to the latest version, and change the sha256 hashcode to the code you recorded in step 1 above
6. Also update any other part of the yaml file as needed to account for the new dependencies/builds steps/etc.
7. Run :code:`git pull upstream master` to make sure you're up to date with the central repo's master branch
8. Check the recipe by cd'ing to the bioconda-recipes project root and running:

.. code-block:: bash

  # This will only need to be run occasionally. Other times it will fail because it has already created a temporary miniconda installation in this location. That's okay
  ./bootstrap.py /tmp/miniconda
  
  source ~/.config/bioconda/activate
  
  # Run bioconda's linter
  bioconda-utils lint --packages <conda package name>
  
  # build and test
  bioconda-utils build --docker --mulled-test --packages <conda package name>

#. Commit and push the changes to your fork repo
#. Create a pull request from your new branch in your fork into master in the central repository
#. Follow the instructions in the pull request to get the update published

Possible Issues:

* bioconda packages should not include large test data files. If the dist/trtools-<version>.tar.gz file contains such files, you'll need to modify the MANIFEST.in file to exclude them,
  fix the test_trtools.sh script to download them manually and point pytest to them, confirm the tests run in a :code:`conda build` and then restart the publishing process.
  (This should not happen if new test files are just put in :code:`trtools/testsupport/sample_vcfs` or :code:`trtools/testsupport/sample_regions`)

Git Tagging
-----------

Git tags are used to mark specific commits with certain names (i.e. v1.2.0).
Please note that tags are assigned to commits, not branches.
You can tag a commit in two different ways.

#. Command line:

.. code-block:: bash

  git tag -a vX.Y.Z -m vX.Y.Z
  git push --tags

2. Web interface: you can go to the releases page of the repository and create a new release.
