
.. image:: https://travis-ci.org/gymreklab/TRTools.svg?branch=master
    :target: https://travis-ci.org/gymreklab/TRTools


.. image:: https://codecov.io/gh/gymreklab/TRTools/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/gymreklab/TRTools


TRTools
============

TRTools includes a variety of utilities for filtering, quality control and analysis of short tandem repeats (STRs) and variable number tandem repeats (VNTRs) downstream of genotyping them from next-generation sequencing. It supports multiple recent genotyping tools (see below).

See full documentation and examples at https://trtools.readthedocs.io/en/latest/.

Install
-------

You can obtain TRTools from pip::

        pip install --upgrade pip
	pip install trtools

(Note: trtools installation may fail for pip version 10.0.1, hence the need to upgrade pip first)

Or, to install from source, run the following command from the base directory of the TRTools repo::

	python setup.py install [--prefix=PREFIX]

to install locally, set :code:`--prefix=$HOME` and ensure :code:`$HOME` is on your :code:`PYTHONPATH`.

(Note, required package :code:`pybedtools` requires zlib. If you receive an error about a missing file :code:`zlib.h`, you can install on Ubuntu using :code:`sudo apt-get install zlib1g-dev` or CentOS using :code:`sudo yum install zlib-devel`.)

Tools
-----
TRTools includes the following tools.

* dumpSTR: a tool for filtering VCF files with STR/VNTR genotypes
* mergeSTR: a tool to merge VCF files across multiple samples genotyped using the same tool
* statSTR: a tool for computing various statistics on VCF files
* compareSTR: a tool for comparing TR callsets
* qcSTR: a tool for generating various quality control plots for a TR callset

Type :code:`<command> --help` to see a full set of options.

It additionally includes a python library, :code:`trtools`, which can be accessed from within Python scripts. e.g.::

	import trtools.utils.utils as stls
	allele_freqs = {5: 0.5, 6: 0.5} # 50% of alleles have 5 repeat copies, 50% have 6
	stls.GetHeterozygosity(allele_freqs) # should return 0.5

Usage
-----

See the README in each subdirectory of <repo-root>/trtools for usage details.

Supported Tools
---------------
TRTools supports VCFs from the following STR/VNTR genotyping tools:

* GangSTR_ version 2.4 or higher.
* HipSTR_ 
* PopSTR_
* ExpansionHunter_
* AdVNTR_

.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://github.com/tfwillems/HipSTR
.. _PopSTR: https://github.com/DecodeGenetics/popSTR
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _AdVNTR: https://github.com/mehrdadbakhtiari/adVNTR

.. _Contributing 

Contributing
------------
We appreciate contributions to TRTools. If you would like to contribute a fix or new feature, follow these guidelines:

* Consider `discussing <https://github.com/gymreklab/TRTools/issues>`_ your solution with us first so we can provide help or feedback if necessary.
* Create a clean environment with the dependencies in requirements.txt installed.
* Additionally, install `pytest` and `pytest-cov <https://anaconda.org/conda-forge/pytest-cov>`_ in your environment.
* Fork the trtools repository. 
* The `develop` branch contains the latest pre-release codebase. Create a branch off of `develop` titled with the name of your feature
* Make your changes. 
* Ensure all functions, modules, classes etc. conform to `numpy docstring standards <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
* Add tests to test any new functionality. Add them to the :code:`tests/` folder in the directory of the code you modified.
* :code:`cd` to the root of the project and run :code:`python -m pytest --cov=. --cov-report term-missing` to make sure that (1) all tests pass and (2) any code you have added is covered by tests. (Code coverage may **not** go down).
* If applicable, update the REAMDEs in the directories of the files you changed with new usage information.
* Submit a pull request **to `develop` branch** of the central repository with a description of what changes you have made.
* A member of the TRTools team will reply and continue the contribution process from there, possibly asking for additional information/effort on your part.

Publishing
----------
If you are a trtools maintainer and wish to publish changes from the develop branch into master and distribute them to PyPI and bioconda,
please see PUBLISHING.rst.
If you are a community member and would like that to happen, contact us (see below).


.. _`Contact Us`
Contact Us
----------
Please submit an issue on the `trtools github <https://github.com/gymreklab/TRTools>`_
# NOTE: I'm intentionally removing Nima's email b/c he will eventually graduate and it will go dead

