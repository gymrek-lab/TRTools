
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

	pip install trtools

Or, to install from source, run the following command from the base directory of the TRTools repo::

	python setup.py install [--prefix=PREFIX]

to install locally, set :code:`--prefix=$HOME` and ensure :code:`$HOME` is on your :code:`PYTHONPATH`.

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

See the README in each subdirectory for usage details.

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

Contributing
------------

If you would like to contribute a fix or new tool to TRTools, follow these guidelines:

* Fork the repository.
* Make your changes. 
* Ensure all functions, modules, classes etc. conform to Numpy docstring standards (https://numpydoc.readthedocs.io/en/latest/format.html).
* Add tests (see :code:`tests/` folder to find the appropriate location for new tests) to test any new functionality. To make sure pytest knows about them, you may need to edit :code:`pytest.ini`.
* Run :code:`pytest --cov=. --cov-report term-missing` to make sure that (1) all tests pass and (2) any code you have added is covered by tests.
* If applicable, update REAMDEs with new usage information.
* Submit a pull request, with a reasonably descriptive message of what changes you have made.

