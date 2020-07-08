TRHarmonizer Library Details
============================

Here are details about how the harmonizer handles VCFs produced by the
different supported genotypers.

Quality Scores
--------------
Quality scores are 0-1 probabilities representing the estimated likelihood that a 
given call is correct.

* AdVNTR quality scores are taken directly from the :code:`ML` (maximum likelihood) 
  format field.
* ExpansionHunter does not output anything that can be interpreted as a quality score
* GangSTR quality scores are taken directly from the :code:`Q` 
  format field.
* HipSTR quality scores are taken directly from the :code:`Q` format field, which represents
  the probability that the call is correct up to choice of phasing. TRHarmonizer ignores the 
  :code:`PQ` field, which represents the probability that the call and its phasing are both correct.
* PopSTR quality scores are taken from the :code:`PL` format field. The meaning of the :code:`PL` 
  field is specified by the `VCF spec <https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf>`_ 
  to mean "Phred-scaled genotype likelihoods rounded to the closest integer". 
  `This wikipedia article <https://en.wikipedia.org/wiki/Phred_quality_score>`_ describes
  Phred scaling. TRHarmonizer takes the score corresponding to the called genotype and 
  converts it to a 0-1 probability. Due to the integer rounding, the Phred scores 
  1, 2, 3, 4, etc. correspond to probabilities in the ranges 
  89%-71%-56%-45-35-28-22-18-14-11-... , as such TRHarmonizer when run on PopSTR records 
  produces quality scores binned by those ranges.   

  Note that the integer rounding is only appropriate for accurately specifying low probability 
  events (where the bin sizes are small), but not for differentiating between the highest 
  probability events (89%-100%). However, a common use-case of quality scores is to distinguish
  between call probabilities that are all in the 90th percentile. As a workaround, TRHarmonizer 
  converts PHRED scores of 0 to

  .. 
    \max (1 - \sum_{\text{other genotypes} Prob(\text{other genotype}), Prob(\text{PHRED score 1}))

  .. math::
     \max \big(
       Prob(\text{PHRED score 1}),
       1 - \sum_{\text{other genotypes}} Prob(\text{other genotype})
     \big)

  where the probabilities for other genotypes are taken from the :code:`PL` field for the same call.
  This provides more distinct binning for calls with quality near 1.

