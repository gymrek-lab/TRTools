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
* ExpansionHunter does not output anything that can be reasonably interpreted as a
  quality score
* GangSTR quality scores are taken directly from the :code:`Q` format field.
* HipSTR quality scores are taken directly from the :code:`Q` format field, which represents
  the quality of a call up to choice of phasing. TRHarmonizer ignores the
  :code:`PQ` field, which represents the quality of a call and its phasing.
* TRTools does not currently support PopSTR quality scores. We are working to figure out
  how best to interpret the :code:`PL` field as a quality score.

Read the relevant genotyper's documentation to understand the meaning of the format field
that TRHarmonzier uses.
