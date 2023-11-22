
.. a location that the doc/index.rst uses for including this file
.. before_header

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

See full documentation and examples at https://trtools.readthedocs.io/en/latest/.

If you use TRTools in your work, please cite: Nima Mousavi, Jonathan Margoliash, Neha Pusarla, Shubham Saini, Richard Yanicky, Melissa Gymrek. (2020) TRTools: a toolkit for genome-wide analysis of tandem repeats. Bioinformatics. (https://doi.org/10.1093/bioinformatics/btaa736)

Install
-------

With conda
^^^^^^^^^^

::

        conda install -c conda-forge -c bioconda trtools

Optionally install :code:`bcftools` which is used to prepare input files for TRTools by running:

::

        conda install -c conda-forge -c bioconda bcftools

Note: Bioconda only supports python versions 3.6-3.8 currently,
so that is all TRTools supports in conda.
If you are using a different version of python we support (3.5 or >= 3.9),
install TRTools using pip.

With pip
^^^^^^^^

First install :code:`htslib` (which contains :code:`tabix` and :code:`bgzip`). Optionally install :code:`bcftools`.
These are used to prepare input files for TRTools and aren't installed by pip.

Then run:

::

        pip install --upgrade pip
        pip install trtools

Note: TRTools installation may fail for pip version 10.0.1, hence the need to upgrade pip first

From source
^^^^^^^^^^^

To install from source (only recommended for development) download the TRTools repository from `github <https://github.com/gymrek-lab/TRTools/>`_,
checkout the branch you're interested in, and run the following command from the base directory of the repo. e.g.::

        git clone https://github.com/gymrek-lab/TRTools
        cd TRTools/
        pip install --upgrade pip
        pip install -e .

Note: required package :code:`pybedtools` requires zlib. If you receive an error about a missing file :code:`zlib.h`, you can install on Ubuntu using :code:`sudo apt-get install zlib1g-dev` or CentOS using :code:`sudo yum install zlib-devel`.

Note: make sure TRTools is not installed in the environment via a different method before installing from source. :code:`which dumpSTR` should return nothing.

Note: if you will run or test :code:`simTR`, you will also need to install `ART <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`_. The simTR tests will only run if the executable :code:`art_illumina` is found on your :code:`PATH`. If it has been installed, :code:`which art_illumina` should return a path.

Tools
-----
TRTools includes the following tools.

* `mergeSTR <https://trtools.readthedocs.io/en/latest/source/mergeSTR.html>`_: a tool to merge VCF files across multiple samples genotyped using the same tool
* `dumpSTR <https://trtools.readthedocs.io/en/latest/source/dumpSTR.html>`_: a tool for filtering VCF files with TR genotypes
* `qcSTR <https://trtools.readthedocs.io/en/latest/source/qcSTR.html>`_: a tool for generating various quality control plots for a TR callset
* `statSTR <https://trtools.readthedocs.io/en/latest/source/statSTR.html>`_: a tool for computing various statistics on VCF files
* `compareSTR <https://trtools.readthedocs.io/en/latest/source/compareSTR.html>`_: a tool for comparing TR callsets
* `associaTR <https://trtools.readthedocs.io/en/latest/source/associaTR.html>`_: a tool for testing TR length-phenotype associations (e.g., running a TR GWAS)
* `prancSTR <https://trtools.readthedocs.io/en/latest/source/prancSTR.html>`_: a tool for identifying somatic mosacisim at TRs. Currently only compatible with HipSTR VCF files. (*beta mode*)
* `simTR <https://trtools.readthedocs.io/en/latest/source/simTR.html>`_: a tool for simulating next-generation sequencing reads from TR regions. (*beta mode*)

Type :code:`<command> --help` to see a full set of options.

It additionally includes a python library, :code:`trtools`, which can be accessed from within Python scripts. e.g.::

        import trtools.utils.utils as stls
        allele_freqs = {5: 0.5, 6: 0.5} # 50% of alleles have 5 repeat copies, 50% have 6
        stls.GetHeterozygosity(allele_freqs) # should return 0.5

Usage
-----

We recommend new users start with the example commands described in the `command-line interface for each tool <https://trtools.readthedocs.io/en/latest/UTILITIES.html>`_.
We also suggest going through our `vignettes <https://trtools.readthedocs.io/en/latest/VIGNETTES.html>`_ that walk through some example workflows using TRTools.

Supported TR Callers
--------------------
TRTools supports VCFs from the following TR genotyping tools:

* AdVNTR_
* ExpansionHunter_
* GangSTR_ version 2.4 or higher
* HipSTR_
* PopSTR_ version 2 or higher

See our description of the `features and example use-cases <https://trtools.readthedocs.io/en/latest/CALLERS.html>`_ of each of these tools.

..
    please ensure this list of links remains the same as the one in the main README

.. _AdVNTR: https://advntr.readthedocs.io/en/latest/
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR

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
#. Additionally, install :code:`pytest`, `pytest-cov <https://anaconda.org/conda-forge/pytest-cov>`_, :code:`sphinx>=3` and :code:`sphinx_rtd_theme`, in your environment.
#. Fork the TRTools repository.
#. The :code:`develop` branch contains the latest pre-release codebase. Create a branch off of :code:`develop` titled with the name of your feature.
#. Make your changes. 
#. Document your changes.

   * Add bullet point(s) to the 'Unreleased Changes' section of :code:`RELEASE_NOTES.rst` describing all the user facing changes you've made (if that section doesn't exist, create it at the top of the file). See prior releases in that file for examples.
   * Ensure all functions, modules, classes etc. conform to `numpy docstring standards <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

    If applicable, update the REAMDEs in the directories of the files you changed with new usage information.

   * New doc pages for `the website <https://trtools.readthedocs.io/en/latest/>`_ can be created under :code:`<project-root>/doc` and linked to as appropriate.
   * If you have added significant amounts of documentation in any of these ways, build the documentation locally to ensure it looks good.

    :code:`cd` to the :code:`doc` directory and run :code:`make clean && make html`, then view :code:`doc/_build/html/index.html` and navigate from there

#. Add tests to test any new functionality. Add them to the :code:`tests/` folder in the directory of the code you modified.

   * :code:`cd` to the root of the project and run :code:`python -m pytest --cov=. --cov-report term-missing` to make sure that (1) all tests pass and (2) any code you have added is covered by tests. (Code coverage may **not** go down).

#. Submit a pull request **to the develop branch** of the central repository with a description of what changes you have made.
   A member of the TRTools team will reply and continue the contribution process from there, possibly asking for additional information/effort on your part.

Publishing
----------
If you are a TRTools maintainer and wish to publish changes from the develop branch into master and distribute them to PyPI and bioconda,
please see PUBLISHING.rst in the root of the git repo.
If you are a community member and would like that to happen, contact us (see above).


