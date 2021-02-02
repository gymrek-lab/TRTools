TRHarmonizer Library Details
============================

Here are details about how the harmonizer library (and the TRTools which are built on top of it)
handle VCFs produced by the different supported genotypers.

Alleles
-------

Different callers will call alleles differently (sequences or just lengths, if sequences
with or without impurities such as partial repeats or internal SNPs, if sequences with or without
flanking SNPs)

* AdVNTR - Calls sequence genotypes with impurities.
* ExpansionHunter - Calls lengths genotypes.
* GangSTR - Calls sequence genotypes without impurities.
* HipSTR - Calls sequence genotypes with impurities. Can call flanking base pairs and SNPs,
  if so TRRecord will have :py:meth:`trtools.utils.tr_harmonizer.TRRecord.GetFullStringGenotypes`.
  HipSTR does not report the underlying motif, so it is inferred from the sequence, which
  is usually but not always correct.
* PopSTR - Calls the reference as a sequence genotype with impurities. Calls alt alleles 
  as length genotypes.

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
  not currently used by TRTools.

More details about these fields can be found in the documentation for each genotyper.
