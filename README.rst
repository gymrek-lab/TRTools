
.. a location that the doc/index.rst uses for including this file
.. before_header

.. image:: https://github.com/codespaces/badge.svg
  :width: 160
  :target: https://codespaces.new/gymrek-lab/TRTools

.. image:: https://github.com/gymrek-lab/trtools/workflows/Tests/badge.svg
    :target: https://github.com/gymrek-lab/trtools/workflows/Tests/badge.svg

.. image:: https://codecov.io/gh/gymrek-lab/TRTools/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/gymrek-lab/TRTools


.. a location that the doc/index.rst uses for including this file
.. after_header

TRTools
=======

.. a location that the doc/index.rst uses for including this file
.. after_title

TRTools includes a variety of utilities for filtering, quality control and analysis of tandem repeats downstream of genotyping them from next-generation sequencing. It supports multiple recent genotyping tools (see below).

See full documentation and examples at https://trtools.readthedocs.io/en/stable/.

If you use TRTools in your work, please cite: Nima Mousavi, Jonathan Margoliash, Neha Pusarla, Shubham Saini, Richard Yanicky, Melissa Gymrek. (2020) TRTools: a toolkit for genome-wide analysis of tandem repeats. Bioinformatics. (https://doi.org/10.1093/bioinformatics/btaa736)

Install
-------

Note: TRTools supports Python versions 3.8 and up. We do not officially support python version 3.7 as it is `end of life <https://devguide.python.org/versions/#status-of-python-versions>`_, but we believe TRTools likely works with it from previous testing results.

With conda
^^^^^^^^^^

::

        conda install -c conda-forge -c bioconda trtools

Optionally install :code:`bcftools` which is used to prepare input files for TRTools (and :code:`ART` which is used by simTR) by running:

::

        conda install -c conda-forge -c bioconda bcftools art

With pip
^^^^^^^^

First install :code:`htslib` (which contains :code:`tabix` and :code:`bgzip`). Optionally install :code:`bcftools`.
These are used to prepare input files for TRTools and aren't installed by pip.

Then run:

::

        pip install --upgrade pip
        pip install trtools

Note: TRTools installation may fail for pip version 10.0.1, hence the need to upgrade pip first

Note: if you will run or test :code:`simTR`, you will also need to install `ART <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`_. The simTR tests will only run if the executable :code:`art_illumina` is found on your :code:`PATH`. If it has been installed, :code:`which art_illumina` should return a path.

From source
^^^^^^^^^^^

To install from source (only recommended for development) clone the TRTools repository from `github <https://github.com/gymrek-lab/TRTools/>`_ and checkout the branch you're interested in::

        git clone -b master https://github.com/gymrek-lab/TRTools
        cd TRTools/

Now, create 1) a conda environment with our development tools and 2) a virtual environment with our dependencies and an editable install of TRTools::

        conda env create -n trtools -f dev-env.yml
        conda run -n trtools poetry install

Now, whenever you'd like to run/import pytest or TRTools, you will first need to activate both environments::

        conda activate trtools
        poetry shell

With Docker
^^^^^^^^^^^

Please refer to `the biocontainers registry for TRTools <https://biocontainers.pro/tools/trtools>`_ for all of our images. To use the most recent release, run the following command::

        docker pull quay.io/biocontainers/trtools:latest

Tools
-----
TRTools includes the following tools.

* `mergeSTR <https://trtools.readthedocs.io/en/stable/source/mergeSTR.html>`_: a tool to merge VCF files across multiple samples genotyped using the same tool
* `dumpSTR <https://trtools.readthedocs.io/en/stable/source/dumpSTR.html>`_: a tool for filtering VCF files with TR genotypes
* `qcSTR <https://trtools.readthedocs.io/en/stable/source/qcSTR.html>`_: a tool for generating various quality control plots for a TR callset
* `statSTR <https://trtools.readthedocs.io/en/stable/source/statSTR.html>`_: a tool for computing various statistics on VCF files
* `compareSTR <https://trtools.readthedocs.io/en/stable/source/compareSTR.html>`_: a tool for comparing TR callsets
* `associaTR <https://trtools.readthedocs.io/en/stable/source/associaTR.html>`_: a tool for testing TR length-phenotype associations (e.g., running a TR GWAS)
* `prancSTR <https://trtools.readthedocs.io/en/stable/source/prancSTR.html>`_: a tool for identifying somatic mosacisim at TRs. Currently only compatible with HipSTR VCF files. (*beta mode*)
* `simTR <https://trtools.readthedocs.io/en/stable/source/simTR.html>`_: a tool for simulating next-generation sequencing reads from TR regions. (*beta mode*)

Type :code:`<command> --help` to see a full set of options.

It additionally includes a python library, :code:`trtools`, which can be accessed from within Python scripts. e.g.::

        import trtools.utils.utils as stls
        allele_freqs = {5: 0.5, 6: 0.5} # 50% of alleles have 5 repeat copies, 50% have 6
        stls.GetHeterozygosity(allele_freqs) # should return 0.5

Usage
-----

We recommend new users start with the example commands described in the `command-line interface for each tool <https://trtools.readthedocs.io/en/stable/UTILITIES.html>`_.
We also suggest going through our `vignettes <https://trtools.readthedocs.io/en/stable/VIGNETTES.html>`_ that walk through some example workflows using TRTools.

Supported TR Callers
--------------------
TRTools supports VCFs from the following TR genotyping tools:

* AdVNTR_
* ExpansionHunter_
* GangSTR_ version 2.4 or higher
* HipSTR_
* PopSTR_ version 2 or higher

See our description of the `features and example use-cases <https://trtools.readthedocs.io/en/stable/CALLERS.html>`_ of each of these tools.

..
    please ensure this list of links remains the same as the one in the main README

.. _AdVNTR: https://advntr.readthedocs.io/en/latest/
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR

Testing
-------
After you've installed TRTools, we recommend running our tests to confirm that TRTools works properly on your system. Just execute the following::

        test_trtools.sh

Development Notes
-----------------

* TRTools only currently supports diploid genotypes. Haploid calls, such as those on male chrX or chrY, are not yet supported but should be coming soon.

Contact Us
----------
Please submit an issue on the `trtools github <https://github.com/gymrek-lab/TRTools>`_

.. _Contributing:

Contributing
------------
We appreciate contributions to TRTools. If you would like to contribute a fix or new feature, follow these guidelines:

1. Consider `discussing <https://github.com/gymrek-lab/TRTools/issues>`_ your solution with us first so we can provide help or feedback if necessary.
#. Install TRTools from source `as above <From source_>`_.
#. Fork the TRTools repository.
#. Create a branch off of :code:`master` titled with the name of your feature.
#. Make your changes. 
#. If you need to add a dependency or update the version of a dependency, you can use the :code:`poetry add` command.

    * You should specify a `version constraint <https://python-poetry.org/docs/master/dependency-specification#version-constraints>`_ when adding a dependency. Use the oldest version compatible with your code. Don't worry if you're not sure at first, since you can (and should!) always update it later. For example, to specify a version of :code:`numpy>=1.23.0`, you can run :code:`poetry add 'numpy>=1.23.0'`.
    * Afterwards, double-check that the :code:`poetry.lock` file contains 1.23.0 in it. **All of our dependencies should be locked to their minimum versions at all times.** To downgrade to a specific version of :code:`numpy` in our lock file, you can explicitly add the version via :code:`poetry add 'numpy==1.23.0'`, manually edit the pyproject.toml file to use a :code:`>=` sign in front of the version number, and then run :code:`poetry lock --no-update`.

#. Document your changes.

   * Ensure all functions, modules, classes etc. conform to `numpy docstring standards <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

    If applicable, update the REAMDEs in the directories of the files you changed with new usage information.

   * New doc pages for `the website <https://trtools.readthedocs.io/en/stable/>`_ can be created under :code:`<project-root>/doc` and linked to as appropriate.
   * If you have added significant amounts of documentation in any of these ways, build the documentation locally to ensure it looks good.

    :code:`cd` to the :code:`doc` directory and run :code:`make clean && make html`, then view :code:`doc/_build/html/index.html` and navigate from there

#. Add tests to test any new functionality. Add them to the :code:`tests/` folder in the directory of the code you modified.

   * :code:`cd` to the root of the project and run :code:`poetry run pytest --cov=. --cov-report term-missing` to make sure that (1) all tests pass and (2) any code you have added is covered by tests. (Code coverage may **not** go down).
   * :code:`cd` to the root of the project and run :code:`nox` to make sure that the tests pass on all versions of python that we support.

#. Submit a pull request (PR) **to the master branch** of the central repository with a description of what changes you have made. Prefix the title of the PR according to the `conventional commits spec <https://www.conventionalcommits.org>`_.
   A member of the TRTools team will reply and continue the contribution process from there, possibly asking for additional information/effort on your part.

   * If you are reviewing a pull request, please double-check that the PR addresses each item in `our PR checklist <https://github.com/gymrek-lab/TRTools/blob/master/.github/pull_request_template.md>`_

Publishing
----------
If you are a TRTools maintainer and wish to publish changes and distribute them to PyPI and bioconda, please see PUBLISHING.rst in the root of the git repo.
If you are a community member and would like that to happen, contact us (see above).
