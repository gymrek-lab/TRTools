[![Build Status](https://travis-ci.org/gymreklab/STRTools.svg?branch=master)](https://travis-ci.org/gymreklab/STRTools)
[![codecov](https://codecov.io/gh/gymreklab/STRTools/branch/master/graph/badge.svg)](https://codecov.io/gh/gymreklab/STRTools)

# STRTools: Toolkit for genome-wide analysis of STRs

<a href="#install">INSTALL</a> | <a href="#tools">TOOLS</a> | <a href="#usage">USAGE</a> | <a href="#supported">SUPPORTED PROGRAMS</a> | <a href="#contributing">CONTRIBUTING</a>

STRTools includes a variety of utilities for filtering, quality control and analysis of short tandem repeats (STRs) and variable number tandem repeats (VNTRs) downstream of genotyping them from next-generation sequencing. It supports multiple recent genotyping tools (see below).

<a name="install"></a>
## Install

You can obtain STRTools from pip:

```
pip install strtools
```

Or, to install from source, run the following command from the base directory of the STRTools repo:

```
python setup.py install [--prefix=PREFIX]
```

to install locally, set `--prefix=$HOME` and ensure `$HOME` is on your `PYTHONPATH`.

<a name="tools"></a>
## Tools
STRTools includes the following tools.

* dumpSTR: a tool for filtering VCF files with STR/VNTR genotypes
* mergeSTR: a tool to merge VCF files across multiple samples genotyped using the same tool
* statSTR: a tool for computing various statistics on VCF files
* compareSTR: a tool for comparing TR callsets

Type `<command> --help` to see a full set of options.

It additionally includes a python library, `strtools`, which can be accessed from within Python scripts. e.g.:

```
import strtools.utils.utils as stls
allele_freqs = {5: 0.5, 6: 0.5} # 50% of alleles have 5 repeat copies, 50% have 6
stls.GetHeterozygosity(allele_freqs) # should return 0.5
```

<a name="usage"></a>
## Usage
See the README in each subdirectory for usage details.

<a name="supported"></a>
## Supported Tools
STRTools supports VCFs from the following STR/VNTR genotyping tools:

* [GangSTR](https://github.com/gymreklab/gangstr) version 2.4 or higher.
* [HipSTR](https://github.com/tfwillems/HipSTR)
* [PopSTR](https://github.com/DecodeGenetics/popSTR)
* [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)
* [AdVNTR](https://github.com/mehrdadbakhtiari/adVNTR)

<a name="contributing"></a>
## Contributing

If you would like to contribute a fix or new tool to STRTools, follow these guidelines:

* Fork the repository.
* Make your changes. 
* Ensure all functions, modules, classes etc. conform to [Numpy docstring standards](https://numpydoc.readthedocs.io/en/latest/format.html).
* Add tests (see `tests/` folder to find the appropriate location for new tests) to test any new functionality. To make sure pytest knows about them, you may need to edit `pytest.ini`.
* Run `pytest --cov=. --cov-report term-missing` to make sure that (1) all tests pass and (2) any code you have added is covered by tests.
* If applicable, update REAMDEs with new usage information.
* Submit a pull request, with a reasonably descriptive message of what changes you have made.

