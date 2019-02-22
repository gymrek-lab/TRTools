# STRTools
Toolkit for genome-wide analysis of STRs

## Install

Run the following command to install:

```
python setup.py install [--prefix=PREFIX]
```
to install locally, set `--prefix=$HOME` and ensure `$HOME` is on your `PYTHONPATH`.

## Tools
STRTools includes the following tools.

* dumpSTR: a tool for filtering VCF files with STR genotypes (from [HipSTR](https://github.com/tfwillems/HipSTR) or [GangSTR](https://github.com/gymreklab/gangstr))
* plinkSTR: a tool for performing STR association studies

Type `<command> --help` to see a full set of options.

## Usage
See the README in each subidrectory for usage details.

## Docker
The Dockerfile in this container sets up a Docker with GangSTR and STRTools installed. See https://hub.docker.com/r/gymreklab/str-toolkit.
