.. TRTools documentation master file, created by
   sphinx-quickstart on Mon Feb 17 11:00:31 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TRTools: a toolkit for genome-wide tandem repeat analysis
=========================================================

.. image:: https://github.com/codespaces/badge.svg
  :width: 160
  :target: https://codespaces.new/gymrek-lab/TRTools

.. image:: https://github.com/gymrek-lab/trtools/workflows/Tests/badge.svg
    :target: https://github.com/gymrek-lab/trtools/workflows/Tests/badge.svg

.. image:: https://codecov.io/gh/gymrek-lab/TRTools/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/gymrek-lab/TRTools


TRTools
=======

TRTools includes a variety of utilities for filtering, quality control and analysis of tandem repeats downstream of genotyping them from next-generation sequencing. It supports multiple recent genotyping tools (see below).

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

If you would like to develop or edit the TRTools source code, you will need to perform a "dev install" directly from the source.

**Note:** Instead of performing the following steps, you can also just open `a GitHub codespace <https://docs.github.com/en/codespaces/overview>`_. Simply type a comma "," when viewing a branch on `our GitHub <https://github.com/gymrek-lab/TRTools>`_ to open an editor with our development setup pre-installed. You can then run :code:`conda activate trtools` and :code:`poetry shell` in the terminal there.

You can clone the TRTools repository from `github <https://github.com/gymrek-lab/TRTools/>`_ and checkout the branch you're interested in::

        git clone -b master https://github.com/gymrek-lab/TRTools
        cd TRTools/

Now, create 1) a conda environment with our development tools and 2) a virtual environment with our dependencies and an editable install of TRTools::

        conda env create -n trtools -f dev-env.yml
        conda run -n trtools poetry install

Now, whenever you'd like to run/import pytest or TRTools, you will first need to activate both environments::

        conda activate trtools
        poetry shell

.. note::
    There's no need to install TRTools this way if you aren't planning to develop or edit the source code! If you want the latest version from our master branch and just can't wait for us to release it, you only need to run::

        pip install --upgrade --force-reinstall git+https://github.com/gymrek-lab/trtools.git@master

With Docker
^^^^^^^^^^^

Please refer to `the biocontainers registry for TRTools <https://biocontainers.pro/tools/trtools>`_ for all of our images. To use the most recent release, run the following command::

        docker pull quay.io/biocontainers/trtools:latest

Tools
-----
TRTools includes the following tools.

* :doc:`mergeSTR </source/mergeSTR>`: a tool to merge VCF files across multiple samples genotyped using the same tool
* :doc:`dumpSTR </source/dumpSTR>`: a tool for filtering VCF files with TR genotypes
* :doc:`qcSTR </source/qcSTR>`: a tool for generating various quality control plots for a TR callset
* :doc:`statSTR </source/statSTR>`: a tool for computing various statistics on VCF files
* :doc:`compareSTR </source/compareSTR>`: a tool for comparing TR callsets
* :doc:`associaTR </source/associaTR>`: a tool for testing TR length-phenotype associations (e.g., running a TR GWAS)
* :doc:`prancSTR </source/prancSTR>`: a tool for identifying somatic mosacisim at TRs. Currently only compatible with HipSTR VCF files. (*beta mode*)
* :doc:`simTR </source/simTR>`: a tool for simulating next-generation sequencing reads from TR regions. (*beta mode*)
* :doc:`annotaTR </source/annotaTR>`: a tool for annotating TR VCF files with dosage or other metadata and optionally converting to PGEN output.

Type :code:`<command> --help` to see a full set of options.

It additionally includes a python library, :code:`trtools`, which can be accessed from within Python scripts. e.g.::

        import trtools.utils.utils as stls
        allele_freqs = {5: 0.5, 6: 0.5} # 50% of alleles have 5 repeat copies, 50% have 6
        stls.GetHeterozygosity(allele_freqs) # should return 0.5

Usage
-----

We recommend new users start with the example commands described in the :doc:`command-line interface for each tool </UTILITIES>`.
We also suggest going through our :doc:`vignettes </VIGNETTES>` that walk through some example workflows using TRTools.

Supported TR Callers
--------------------
TRTools supports VCFs from the following TR genotyping tools:

* AdVNTR_
* ExpansionHunter_
* GangSTR_ version 2.4 or higher
* HipSTR [`main repo <https://github.com/tfwillems/HipSTR>`_] [`Gymrek Lab repo <https://github.com/gymrek-lab/hipstr>`_]
* PopSTR_ version 2 or higher
* LongTR_

See our description of the :doc:`features and example use-cases </CALLERS>` of each of these tools.

..
    please ensure this list of links remains the same as the one in CALLERS.rst

.. _AdVNTR: https://advntr.readthedocs.io
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR
.. _LongTR: https://github.com/gymrek-lab/longtr

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

    * You should specify a `version constraint <https://python-poetry.org/docs/master/dependency-specification#version-constraints>`_ when adding a dependency. Use the oldest version compatible with your code. For example, to specify a version of :code:`numpy>=1.23.0`, you can run :code:`poetry add 'numpy>=1.23.0'`.
    * Afterward, double-check that the :code:`poetry.lock` file contains 1.23.0 in it. **All of our dependencies should be locked to their minimum versions at all times.** To downgrade to a specific version of a package in our lock file, you can `explicitly add the version <https://python-poetry.org/docs/dependency-specification/#projectdependencies-and-toolpoetrydependencies>`_ to the :code:`tool.poetry.dependencies` of our pyproject.toml file and then run :code:`poetry lock`.
    * Only PyPI packages can be added to our pyproject.toml file. So if a dependency is only available on conda, then you can add it to our :code:`dev-env.yml` file instead. Please note that anyone who installs TRTools from PyPI will not be guaranteed to have your dependency installed, so you should design your code accordingly. 
    * Any changes to our dependencies must also added to our bioconda recipe at the time of publication. See `PUBLISHING.rst <https://github.com/gymrek-lab/TRTools/blob/master/PUBLISHING.rst>`_ for more details.

#. Document your changes.

   * Ensure all functions, modules, classes etc. conform to `numpy docstring standards <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

    If applicable, update the REAMDEs in the directories of the files you changed with new usage information.

   * New doc pages for `the website <https://trtools.readthedocs.io>`_ can be created under :code:`<project-root>/doc` and linked to as appropriate.
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


Table of Contents
=================

.. toctree::
   :maxdepth: 2

   VIGNETTES
   UTILITIES
   CALLERS
   LIBRARY
   Release Notes <https://github.com/gymrek-lab/TRTools/blob/master/CHANGELOG.md>
   site_indices

