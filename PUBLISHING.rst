Publishing
----------

Only maintainers of the trtools repository (see @gymrek-lab/trtools) may publish changes to the package.
If you are a community member and want to contribute new code, see the contributing section in the README.

This document explains how trtools maintainers should publish new changes.
Maintainers should reach consensus before going ahead with publishing changes.
Note that the publishing to PyPI step below will require credentials
that are only made available to core maintainers of TRTools.

Note that we use the master branch as the default branch with the latest stable codebase.
The builds from this branch are distributed to PyPI and conda.
Other branches are used for development.

New Dependencies
----------------
If you've added dependencies to trtools or its tests, those dependencies should be listed in

  * pyproject.toml
  * the .readthedocs_conda_env.yml file in the root of the repository that's used for building
    TRTool's Read The Docs webpage.
  * the appropriate section of the bioconda recipe (see below)

Publishing Steps
----------------

To publish a new version of trtools:

1. First, locate the most recent PR prefixed "chore(main)" created by our Github actions bot
2. List a maintainer of our repository (@gymrek-lab/trtools) as a reviewer of the PR and ask them to merge it
3. The bot will automatically create a new version on PyPI and tag a release on Github

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
