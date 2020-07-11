TRHarmonizer Library Details
============================

Here are details about how the harmonizer library (and the TRTools which are built on top of it)
handle VCFs produced by the different supported genotypers.

.. _Quality Scores:

Quality Scores
--------------

Quality scores are numbers associated with each call in a VCF where a higher number means a
better call. The reliability of these numbers is dependent on the genotyper that emitted them.
Here is how TRHarmonizer infers quality scores for each supported genotyper:

* AdVNTR quality scores are taken directly from the :code:`ML` (maximum likelihood)
  format field.
* ExpansionHunter does not output a quality score. It does output a confidence interval, 
  but this field is not currently used by TRTools.
* GangSTR quality scores are taken directly from the :code:`Q` format field.
* HipSTR quality scores are taken directly from the :code:`Q` format field. TRHarmonizer ignores the
  :code:`PQ` field, which represents the quality of a call and its phasing.
* PopSTR does not output a genotype quality score. It does output a :code:`PL` field 
  with Phred-scaled genotype likelihoods, but this field is
  not currently used by TRTools,

More details about these fields can be found in the documentation for each genotyper.
