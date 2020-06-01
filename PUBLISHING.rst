Publishing
----------

Only maintainers of the trtools repository may publish changes to the package.
If you are a community member and want to contribute new code, see the contributing section in the README.
If you are a community member and have already contributed new code and want us to publish it
now, please contact us (our contact info is in the README)

This document explains how trtools maintainers should publish new changes. 
Maintainers should reach consensus before going ahead with publishing changes.

We use a simplified version of 
`git flow <http://web.archive.org/web/20200520162709/https://nvie.com/posts/a-successful-git-branching-model/>`_
to maintain and publish trtools.
We use the master branch as the default branch with the latest stable codebase.
The builds from this branch are distributed to PyPI and conda.
The develop branch contains new features that have yet to make their way into master.

Once changes have been made to develop that are ready to be published, first set up the environment you're going to publish TRTools from:

#. Create a clean environment.
#. Install setuptools with version >= 40.8.0
#. Install all the requirements in requirements.txt
#. Additionally, install ``pytest``, ``wheel`` and ``twine``

Then go through the steps of merging the changes into the master branch:

#. Clone the `trtools repo <https://github.com/gymreklab/TRTools>`_
#. Check out the develop branch
#. Run pytest and make sure all the tests pass
#. Update the version number in setup.py
#. Run ``python setup.py sdist bdist_wheel`` (this ensures that trtools/version.py contains the updated version number)
#. Commit the changes to setup.py and trtools/version.py and push them>
#. Submit a pull request from develop into master on the github webiste.
#. If the code review and travis checks pass, merge the pull request.
#. Tag the merge commit with the package version in vX.Y.Z format. (For more details on tagging, see `below`)

Then go through the steps of publishing the changed code to PyPI

1. :code:`cd` into the root of your clone of the trtools repo, checkout master and pull the latest change.
2. Run :code:`rm -rf build dist *.egg-info` to make sure all previous build artifacts are removed
3. Run :code:`python setup.py sdist bdist_wheel` to build the package.

 This will create the warning::

   UserWarning: Unknown distribution option: 'license_file'  warnings.warn(msg)

 You can ignore this warning: the 'license_file' option is necessary for creating the build artifacts

4. Run :code:`twine upload dist/*` to upload the build to PyPI

Lastly go through the following stesp to publish the changed code to bioconda (see `here <http://bioconda.github.io/contributor/workflow.html>`_ for official documentation)

1. Run :code:`cd <project-root> && openssl sha256 dist/trtools-<version>.tar.gz` and save the generated hash code for later
2. Create a fork of the `bioconda recipes <https://github.com/bioconda/bioconda-recipes>`_ repo and clone it.
3. Run ``git remote add upstream https://github.com/bioconda/bioconda-recipes`` so that you can pull from not only your fork but also upstream.
4. Create a new branch from master with your package's name
5. In `the recipe for trtools <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/trtools/meta.yaml#L1-L2>`_, update the version number to the latest version, and change the sha256 hashcode to the code you recorded in step 1 above
6. Check the recipe by running the following

.. code-block:: bash

  # This will only need to be run occasionally, other times it will fail because the temporary miniconda installation it creates already exists, that's okay
  ./bootstrap.py /tmp/miniconda
  
  source ~/.config/bioconda/activate
  
  # Run bioconda's linter
  bioconda-utils lint --packages <conda package name>
  
  # build and test
  bioconda-utils build --docker --mulled-test --packages <conda package name>

7. Before committing, run `git pull upstream master` to make sure you're up to date with the central repo's master branch
#. Commit and push the changes to your fork repo
#. Create a pull request from your new branch in your fork into master in the central repository
#. Follow the instructions in the pull request to get the update published

Possible Issues:

* bioconda packages should not include large test data files. If the dist/trtools-<version>.tar.gz file contains such files, you'll need to modify the MANIFEST.in file to exclude them, fix the test_trtools.sh script to download them manually and point pytest to them, and then recommit everything. (This should not happen if new test files are just put in trtools/testsupport/sample_vcfs or trtools/testsupport/sample_regions)
* If you've added dependencies to trtools or its tests, those dependencies should be listed in

  * setup.py
  * requirements.txt (list a specific version of the dependency that is up to date and that we know will work)
  * the appropriate section of the bioconda recipe

Git Tagging
-----------

Git tags are used to mark specific commits with certain names (i.e. v1.2.0). 
Please note that tags are assigned to commits, not branches. 
You can tag a commit in two different ways.

#. Command line::

  git tag -a <tag-name> -m <tag-description>
  git push --tags

2. Web interface: you can go to the releases page of the repository and create a new release.
