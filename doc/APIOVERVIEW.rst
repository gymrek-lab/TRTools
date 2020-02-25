API Overview
============

In addition to the command-line utilities, TRTools has a Python library which can be accessed in custom scripts.

TR String Utilities
-------------------

The module :code:`trtools.utils.utils` contains various helpful functions, e.g. for inferring a repeat sequence motif from a string, converting repeat motifs to canonical representations, or computing functions like heterozygosity from allele frequency distributions::

  import trtools.utils.utils as strutils  
  strutils.InferRepeatSequence('ATATATATATATA', 2) # returns 'AT'
  strutils.GetCanonicalMotif('CAG') # returns 'ACG'
  afreqs = {10:0.25, 11:0.5, 12: 0.25} # frequency of each length
  strutils.GetHeterozygosity(afreqs) # returns 0.625

See https://trtools.readthedocs.io/en/latest/trtools.utils.utils.html for a complete list of utility functions.

TR Harmonization across tools
-----------------------------

TODO
