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
* TRTools does not currently support PopSTR quality scores. We are working to figure out
  if the :code:`PL` field can be appropriately interpreted (after scaling) as a 0-1
  probability.

