API Overview
============

In addition to the command-line utilities, STRTools has a Python library which can be accessed in custom scripts.

TR String Utilities
-------------------

The module :code:`strtools.utils.utils` contains various helpful functions, e.g. for inferring a repeat sequence motif from a string, converting repeat motifs to canonical representations, or computing functions like heterozygosity from allele frequency distributions::

  import strtools.utils.utils as strutils  
  strutils.InferRepeatSequence('ATATATATATATA', 2) # returns 'AT'
  strutils.GetCanonicalMotif('CAG') # returns 'ACG'
  afreqs = {10:0.25, 11:0.5, 12: 0.25} # frequency of each length
  strutils.GetHeterozygosity(afreqs) # returns 0.625

See https://strtools.readthedocs.io/en/latest/strtools.utils.utils.html for a complete list of utility functions.

TR Harmonization across tools
-----------------------------

The module :code:`strtools.utils.tr_harmonizer` is responsible for providing a caller agnostic view of a VCF containing TR records. The class you will work with that provides this functionality is :code:`strtools.utils.tr_harmonizer.TRRecord`. There are two coding paradigms for accessing this API. If you just want to iterate through the TRRecords in a vcf, use the TRRecordHarmonizer::

  import vcf
  import strtools.utils.tr_harmonizer as trh
  
  invcf = vcf.Reader(filename = ...)
  harmonizer = trh.TRRecordHarmonizer(invcf)
  for trrecord in harmonizer:
        # do something with the trrecord 
        for sample in trrecord:
                print(sample, trrecord.GetLengthGenotype(sample))

If you want to first work with the underlying record that pyvcf returns, and then later convert it to a TRRecord, use the module method HarmonizeRecord::

  import vcf
  import strtools.utils.tr_harmonizer as trh

  invcf = vcf.Reader(filename = ...)
  #make sure to grab the vcf's filetype
  vcftype = trh.InferVCFType(invcf)
  for record in invcf:
    # do some filtering on the record
    passesFilters = filter(record)

    # later on, convert the record to a trrecord
    if not passesFilters:
       continue

    trrecord = trh.HarmonizeRecord(vcftype, record)
    # do something with the tr record
    for sample in trrecord:
       print(sample, trrecord.GetLengthGenotype(sample))

See https://strtools.readthedocs.io/en/latest/strtools.utils.tr_harmonizer.html for a complete list of all functions used for inspecting vcfs and TRRecords.

